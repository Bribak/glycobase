library(shiny)
library(httr)
library(rjson)
library(DT)
library(plotly)

# profvis::profvis({ runApp() })
DEBUG <- Sys.getenv('DEBUG') == 'TRUE'
OFFLINE <- FALSE
REQUIRE_LOGIN <- FALSE

# TODO
# Add GBID to the tSNE tables for easier mapping and filtering by species, etc

# ----------------------------------- App UI --------------------------------- #

# Check whether user has auth'd
has_auth_code <- function(params) {
  # params is a list object containing the parsed URL parameters. Return TRUE if
  # based on these parameters, it looks like auth code is present that we can
  # use to get an access token. If not, it means we need to go through the OAuth
  # flow.
  return(!is.null(params$code))
}

# UI will change depending on whether the user has logged in
uiFunc <- function(req) {
  if (OFFLINE | !REQUIRE_LOGIN){
    AuthenticatedUI
  } else {
    if (!has_auth_code(parseQueryString(req$QUERY_STRING))) {
      # Login button
      AnonymousUI
    } else {
      # App UI
      AuthenticatedUI
    }
  }
}

# Import UI to be shown after user before and after auth'd
source('app_ui.R')
if (!dir.exists('data')){
  dir.create('data')
}

# ------------------------ Virtualenv setup -------------------------- #
if (Sys.info()[['sysname']] != 'Darwin'){
  # When running on shinyapps.io, create a virtualenv 
  reticulate::virtualenv_create(envname = 'python35_gly_env', 
                                python = '/usr/bin/python3')
  reticulate::virtualenv_install('python35_gly_env', 
                                 packages = c('synapseclient', 'requests',
                                              'pandas'))
}
reticulate::use_virtualenv('python35_gly_env', required = T)

# ---------------------------- OAuth --------------------------------- #

if (!OFFLINE){
  reticulate::source_python('connect_to_synapse.py')
  # Initialize Synapse client
  login_to_synapse(username = Sys.getenv('SYN_USERNAME'),
                   api_key = Sys.getenv('SYN_API_KEY'))
  logged_in <- reactiveVal(FALSE)
  source('oauth.R')
}

# ----------------------------------- Server --------------------------------- #

server <- function(input, output, session) {
  
  if (REQUIRE_LOGIN){
    # Click on the 'Log in' button to kick off the OAuth round trip
    observeEvent(input$action, {
      session$sendCustomMessage("customredirect", oauth2.0_authorize_url(API, APP, scope = SCOPE))
      return()
    })
    
    params <- parseQueryString(isolate(session$clientData$url_search))
    if (!has_auth_code(params)) {
      return()
    }
    
    url <- paste0(API$access, '?', 'redirect_uri=', APP_URL, '&grant_type=', 
                 'authorization_code', '&code=', params$code)
    
    # Get the access_token and userinfo token
    token_request <- POST(url,
                          encode = 'form',
                          body = '',
                          authenticate(APP$key, APP$secret, type = 'basic'),
                          config = list()
    )
    
    stop_for_status(token_request, task = 'Get an access token')
    token_response <- httr::content(token_request, type = NULL)
    
    access_token <- token_response$access_token
    id_token <- token_response$id_token
    if (token_request$status_code == 201){
      logged_in(T)
    }
    
    # ------------------------------ App --------------------------------- #
    
    # Get information about the user
    user_response = get_synapse_userinfo(access_token)
    user_id = user_response$userid
    user_content_formatted = paste(lapply(names(user_response), 
                                          function(n) paste(n, user_response[n])), collapse="\n")
    
    # Get user profile
    profile_response <- get_synapse_user_profile()
    
    # Cache responses
    if (DEBUG){
      saveRDS(token_response, 'cache/token_response.rds')
      saveRDS(user_response, 'cache/user_response.rds')
      saveRDS(profile_response, 'cache/profile_response.rds')
    }
    
    output$userInfo <- renderText(user_content_formatted)
    output$teamInfo <- renderText(teams_content_formatted)
    # See in app_ui.R with verbatimTextOutput("userInfo")
  
    # ---------------------------- Menus --------------------------------- #
    
    # Logout modal
    observeEvent(input$user_account_modal, {
      showModal(
        modalDialog(title = "Synapse Account Information",
                    h4(paste0(profile_response$firstName, ' ', profile_response$lastName)),
                    p(profile_response$company),
                    p(user_response$email, style = 'color: #00B07D;'),
                    easyClose = T,
                    footer = tagList(
                      actionButton("button_view_syn_profile", "View Profile on Synapse",
                                   style = 'color: #ffffff; background-color:  #00B07D; border-color: #0f9971ff;',
                                   onclick = paste0("window.open('https://www.synapse.org/#!Profile:", profile_response$ownerId, "', '_blank')")),
                      modalButton("Back to Analysis")
                      #actionButton("button_logout", "Log Out")
                    )
        )
      )
    })
    
    output$logged_user <- renderText({
      if(logged_in()){
        return(paste0('Welcome, ', profile_response$firstName, '!'))
      }
    })
    
  } else {
    
    # Logout modal
    observeEvent(input$user_account_modal, {
      showModal(
        modalDialog(title = "Welcome, Guest!",
                    p('Thanks for using GlycoBase. Looking for raw data?'),
                    a('View the GlycoBase data on Synapse', 
                      href='https://www.synapse.org/#!Synapse:syn21568077/wiki/600880',
                      target='_blank'),
                    easyClose = T,
                    footer = tagList(
                      modalButton("Back to Analysis")
                    )
        )
      )
    })
    
    output$logged_user <- renderText({
      paste0('Welcome, Guest!')
    })
    
  }
  # end REQUIRE_LOGIN
  
  # Citation info modal
  observeEvent(input$citation_modal, {
    showModal(modalDialog(
      title = 'Citing GlycoBase',
      p('When using GlycoBase in your research, please cite the following:'),
      div(style = 'padding-left: 30px;',
        p('D. Bojar, D.M. Camacho, J.J. Collins. Using Natural Language Processing to Learn the Grammar of Glycans.'),
        a('Preprint available on bioRxiv', href = 'https://www.biorxiv.org/content/10.1101/2020.01.10.902114v1',
          target = '_blank',
          style = 'color: #00B07D;')
      ),
      easyClose = T,
      footer = NULL
    ))
  })

  
  # ----------------------------- TAB 1: OVERVIEW --------------------------- #
  
  glycobaseData <- reactiveValues(glycobase_df = NULL,
                                  monosaccharides = NULL,
                                  species = NULL,
                                  tsne_glycans_df = NULL,
                                  tsne_glycoletters_df = NULL,
                                  tsne_glycowords_df = NULL,
                                  num_glycans = 19299,
                                  num_glycoletters = 1027,
                                  num_glycowords = 19866)
  
  # Load the current version of GlycoBase to display
  DATASET = PROJECT_CONFIG$DATASETS$V2
  
  if (OFFLINE){
    glycobase_csv = 'data/v2_glycobase.csv'
    tsne_glycans_csv = 'data/v2_tsne_glycans_isomorph.csv'
    tsne_glycoletters_csv = 'data/v2_tsne_glycoletters.csv'
    tsne_glycowords_csv = 'data/v2_tsne_glycowords.csv'
    
  } else{
    glycobase_csv = fetch_synapse_filepath(DATASET$glycobase)
    tsne_glycans_csv = fetch_synapse_filepath(DATASET$tsne_glycans)
    tsne_glycoletters_csv = fetch_synapse_filepath(DATASET$tsne_glycoletters)
    tsne_glycowords_csv = fetch_synapse_filepath(DATASET$tsne_glycowords)
  }
  
  # Load glycobase
  glycobase_df <- read.csv(glycobase_csv, stringsAsFactors = F)
  glycobase_df$glycan_id = paste0('GBID', glycobase_df$glycan_id)
  glycobase_df$species = gsub("\\['|\\']|'", '', glycobase_df$species)
  glycobase_df$species = gsub("_", ' ', glycobase_df$species)
  humans = glycobase_df[grepl('Homo sapiens', glycobase_df$species), ]
  others = glycobase_df[!grepl('Homo sapiens', glycobase_df$species), ]
  glycobase_df = rbind(humans, others)
  glycobase_df$immunogenicity[glycobase_df$immunogenicity == 0] = 'No'
  glycobase_df$immunogenicity[glycobase_df$immunogenicity == 1] = 'Yes'
  glycobase_df$immunogenicity[is.na(glycobase_df$immunogenicity)] = 'Unknown'
  glycobase_df$link[glycobase_df$link == ''] = 'None'
  glycobase_df$link[glycobase_df$link == 'free'] = 'Free'
  names(glycobase_df)[c(1:4,6)] = c('GlycoBase_ID', 'Glycan', 'Species', 'Immunogenic', 'Link')
  glycobaseData$glycobase_df <- glycobase_df
  
  # Unique monosaccharides, bonds, species for filtering
  glycobaseData$monosaccharides <- readRDS('rdata/v2_monosaccharides.rds')
  glycobaseData$bonds <- readRDS('rdata/v2_bonds.rds')
  glycobaseData$species <- readRDS('rdata/v2_species.rds')
  glycobaseData$kingdoms <- readRDS('rdata/v2_kingdoms.rds')
  
  # All data loaded
  glyco_data <- reactive({ glycobaseData$glycobase_df })
  tsne_glycans_data <- reactive({ glycobaseData$tsne_glycans_df })
  tsne_glycowords_data <- reactive({ glycobaseData$tsne_glycowords_df })
  tsne_glycoletters_data <- reactive({ glycobaseData$tsne_glycoletters_df })
  n_glycans <- reactive({ glycobaseData$num_glycans })
  n_glycoletters <- reactive({ glycobaseData$num_glycoletters })
  n_glycowords <- reactive({ glycobaseData$num_glycowords })
  monos <- reactive({ glycobaseData$monosaccharides })
  bonds <- reactive({ glycobaseData$bonds })
  kingdoms <- reactive({ glycobaseData$kingdoms })
  specs <- reactive({ glycobaseData$species })
  
  # Three boxes on first row
  output$num_glycans <- renderText({ n_glycans() })
  output$num_glycowords <- renderText({ n_glycowords() })
  output$num_glycoletters <- renderText({ n_glycoletters() })
  
  # tSNE plot for glycans
  output$tsne_glycans <- renderPlotly({
    
    # Load glycans tSNE if needed
    if (is.null(glycobaseData$tsne_glycans_df)){
      tsne_glycans <- read.csv(tsne_glycans_csv,
                               stringsAsFactors = F)
      names(tsne_glycans) = c('Glycan', 'Dim1', 'Dim2')
      glycobaseData$tsne_glycans_df <- tsne_glycans
    }
    
    plot_df = tsne_glycans_data()
    p <- plot_ly(data = plot_df, x = ~Dim1, y = ~Dim2,
                 symbols = 21,
                 text = ~Glycan,
                 hovertemplate = '<b>Glycan:</b> %{text}<extra></extra>',
                 type = 'scatter', mode = 'markers',
                 showlegend = FALSE,
                 marker = list(size = 8,
                               # TODO color them by a meaningful grouping?
                               color = '#01B07D',
                               line = list(
                                 color = '#212D32',
                                 width = 1
                               )
                 )
    ) %>%
      layout(title = 'Glycans (tSNE plot)')
    
    # Display Plotly plot in UI
    return(p)
  })
  
  # tSNE plot for glycowords
  output$tsne_glycowords <- renderPlotly({
    
    # Load glycowords tSNE if needed
    if (is.null(glycobaseData$tsne_glycowords_df)){
      tsne_glycowords = read.csv(tsne_glycowords_csv,
                                 stringsAsFactors = F)
      tsne_glycowords = tsne_glycowords[,c(3,1,2)]
      names(tsne_glycowords) = c('Glycoword', 'Dim1', 'Dim2')
      tsne_glycowords$Glycoword = gsub("\\(|\\)|'|,", "", tsne_glycowords$Glycoword)
      glycobaseData$tsne_glycowords_df <- tsne_glycowords
    }
    
    # Load glycans tSNE if needed
    if (is.null(glycobaseData$tsne_glycocans_df)){
      tsne_glycans <- read.csv(tsne_glycans_csv,
                               stringsAsFactors = F)
      names(tsne_glycans) = c('Glycan', 'Dim1', 'Dim2')
      glycobaseData$tsne_glycocans_df <- tsne_glycans
    }
    
    plot_df = tsne_glycowords_data()
    p <- plot_ly(data = plot_df, x = ~Dim1, y = ~Dim2,
                 symbols = 21,
                 text = ~Glycoword,
                 hovertemplate = '<b>Glycoword:</b> %{text}<extra></extra>',
                 type = 'scatter', mode = 'markers',
                 showlegend = FALSE,
                 marker = list(size = 8,
                               # TODO color them by a meaningful grouping?
                               color = '#019DB0',
                               line = list(
                                 color = '#212D32',
                                 width = 1
                               )
                 )
    ) %>%
      layout(title = 'Glycowords (tSNE plot)')
    
    # Display Plotly plot in UI
    return(p)
  })
  
  # tSNE plot for glycowords
  output$tsne_glycowords <- renderPlotly({
    
    # Load glycowords tSNE if needed
    if (is.null(glycobaseData$tsne_glycowords_df)){
      tsne_glycowords = read.csv(tsne_glycowords_csv,
                                 stringsAsFactors = F)
      tsne_glycowords = tsne_glycowords[,c(3,1,2)]
      names(tsne_glycowords) = c('Glycoword', 'Dim1', 'Dim2')
      tsne_glycowords$Glycoword = gsub("\\(|\\)|'|,", "", tsne_glycowords$Glycoword)
      glycobaseData$tsne_glycowords_df <- tsne_glycowords
    }
    
    plot_df = tsne_glycowords_data()
    p <- plot_ly(data = plot_df, x = ~Dim1, y = ~Dim2,
                 symbols = 21,
                 text = ~Glycoword,
                 hovertemplate = '<b>Glycoword:</b> %{text}<extra></extra>',
                 type = 'scatter', mode = 'markers',
                 showlegend = FALSE,
                 marker = list(size = 8,
                               # TODO color them by a meaningful grouping?
                               color = '#019DB0',
                               line = list(
                                 color = '#212D32',
                                 width = 1
                               )
                 )
    ) %>%
      layout(title = 'Glycowords (tSNE plot)')
    
    # Display Plotly plot in UI
    return(p)
  })
  
  
  # tSNE plot for glycoletters
  output$tsne_glycoletters <- renderPlotly({
    
    # Load glycoletters tSNE if needed
    if (is.null(glycobaseData$tsne_glycoletters_df)){
      tsne_glycoletters <- read.csv(tsne_glycoletters_csv,
                                    stringsAsFactors = F)
      tsne_glycoletters = tsne_glycoletters[,c(3,1,2)]
      names(tsne_glycoletters) = c('Glycoletter', 'Dim1', 'Dim2')
      tsne_glycoletters = tsne_glycoletters[tsne_glycoletters$Glycoletter != '',]
      glycobaseData$tsne_glycoletters_df <- tsne_glycoletters
    }
    
    plot_df = tsne_glycoletters_data()
    p <- plot_ly(data = plot_df, x = ~Dim1, y = ~Dim2,
                 symbols = 21,
                 text = ~Glycoletter,
                 hovertemplate = '<b>Glycoletter:</b> %{text}<extra></extra>',
                 type = 'scatter', mode = 'markers',
                 showlegend = FALSE,
                 marker = list(size = 8,
                               # TODO color them by a meaningful grouping?
                               color = '#903C83',
                               line = list(
                                 color = '#212D32',
                                 width = 1
                               )
                 )
    ) %>%
      layout(title = 'Glycoletters (tSNE plot)')
    
    # Display Plotly plot in UI
    return(p)
  })
  
  # Glycans modal
  observeEvent(input$modal_glycans, {
    showModal(modalDialog(
      title = 'Unique glycans in GlycoBase',
      div(style='padding-left: 20px;',
          p(style='color: #D2D6DD;', 'Hover over a point to view glycan')
      ),
      div(withSpinner(plotlyOutput('tsne_glycans'), type = 4, color = '#00B07D'),
          style = "overflow-y: auto;"),
      selectInput('select_tsne_glycans', 'Highlight glycans containing monosaccharide:', 
                  choices = c('All', monos()), selected = 'All', selectize = T),
      easyClose = T,
      footer = NULL
    ))
  })
  
  # Glycowords modal
  observeEvent(input$modal_glycowords, {
    showModal(modalDialog(
      title = 'Unique glycowords in GlycoBase',
      div(style='padding-left: 20px;',
        p(style='color: #D2D6DD;', 'Hover over a point to view glycoword')
      ),
      div(withSpinner(plotlyOutput('tsne_glycowords'), type = 4, color = '#019DB0'),
          style = "overflow-y: auto;"),
      selectInput('select_tsne_glycowords', 'Highlight glycowords containing monosaccharide:', 
                  choices = c('All', monos()), selected = 'All', selectize = T),
      easyClose = T,
      footer = NULL
    ))
  })
  
  # Glycoletters modal
  observeEvent(input$modal_glycoletters, {
    showModal(modalDialog(
      title = 'Unique glycoletters in GlycoBase',
      div(style='padding-left: 20px;',
          p(style='color: #D2D6DD;', 'Hover over a point to view glycoletter')
      ),
      div(withSpinner(plotlyOutput('tsne_glycoletters'), type = 4, color = '#903C83'),
          style = "overflow-y: auto;"),
      selectInput('select_tsne_glycoletters', 'Highlight monosaccharide:', 
                  choices = c('All', monos()), selected = 'All', selectize = T),
      easyClose = T,
      footer = NULL
    ))
  })
  
  
  # Show the whole glycobase table
  output$table_glycobase <- DT::renderDT({
    df = glyco_data()
    df = df[,c('GlycoBase_ID', 'Glycan', 'Link', 'Species', 'Immunogenic')]
    df$Link = factor(df$Link, c('N', 'O', 'Free', 'None'))
    df$Immunogenic = factor(df$Immunogenic, c('Yes', 'No', 'Unknown'))
    df$Species = factor(df$Species, specs())
    
    if (!is.null(df)){
      
      # Display
      return(datatable(df, rownames = F, selection = 'none',
                       style = 'bootstrap', escape = F, 
                       filter = 'top',
                       options = list(
                         dom = 'tipr',
                         pageLength = 20,
                         autoWidth = TRUE,
                         columnDefs = list(list(width = '100px', targets = c(0,2,4)),
                                           list(width = '200px', targets = 3))
                       )))
    } else{
      return(NULL)
    }
  })
  
  # ---------------------- TAB 2: STRUCTURAL CONTEXT ------------------------ #
  reticulate::source_python('structural_context.py')
  
  context_query_criteria <- reactive({ input$select_context_criteria })
  selected_context_glycoletter <- reactive({ input$select_context_glycoletter })
  context_tax_level <- reactive({ input$select_context_taxonomy_level })
  selected_tax_value <- reactive({ input$select_context_taxonomy_value })
  
  observeEvent(context_query_criteria(), {
    criteria = context_query_criteria()
    if (criteria == 'Observed monosaccharides making bond:'){
      updateSelectInput(session, 'select_context_glycoletter', 
                        choices = bonds())
    } else{
      updateSelectInput(session, 'select_context_glycoletter', 
                        choices = monos())
    }
  })
  
  observeEvent(context_tax_level(), {
    tax_level = context_tax_level()
    if (tax_level == 'Kingdom'){
      updateSelectInput(session, 'select_context_taxonomy_value', 
                        label = 'Kingdom',
                        choices = c('All', kingdoms()))
    } else{
      updateSelectInput(session, 'select_context_taxonomy_value',
                        label = 'Species',
                        choices = c('All', specs()))
    }
  })
  
  output$context_observed_barplot <- renderPlotly({
    
    motif = selected_context_glycoletter()
    criteria = context_query_criteria()
    if (criteria == 'Observed monosaccharides making bond:'){
      mode = 'bond'
    } else if (criteria == 'Observed monosaccharides paired with:'){
      mode = 'sugar'
    } else{
      mode = 'sugarbond'
    }
    
    context_results = characterize_context(motif, mode = mode, 
                                           taxonomy_filter=context_tax_level(), taxonomy_value = selected_tax_value())
    
    plot_title = context_results[[1]]
    plot_x = factor(context_results[[2]], levels = context_results[[2]]) # keep the order that characterize_context() returns
    plot_y = context_results[[3]]
    plot_ly(x = plot_x, y = plot_y, type = 'bar', marker = list(color = '#00B07D')) %>%
      layout(title = plot_title,
             xaxis = list(title = 'Monosaccharide'),
             yaxis = list(title = 'Number of glycans'))
    
  })
  
  output$context_main_side_barplot <- renderPlotly({
    
    motif = selected_context_glycoletter()
    
    main_side = main_v_side_branch(motif, taxonomy_filter=context_tax_level(), taxonomy_value = selected_tax_value())
    
    plot_title = motif
    plot_x = c('Main branch', 'Side branch')
    plot_y = main_side
    plot_ly(x = plot_x, y = plot_y, type = 'bar', 
            marker = list(color = c('#019DB0', '#903C83'))) %>%
      layout(title = plot_title,
             xaxis = list(title = 'Position'),
             yaxis = list(title = 'Ocurrence'))
    
  })
  
}

# uiFunc instead of ui
shinyApp(uiFunc, server)
