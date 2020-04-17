from collections import Counter
from glycan_processing_helpers import df_species, link_find
import pandas as pd
import pickle
import re


def main_v_side_branch(glycoletter, taxonomy_filter = 'Kingdom', taxonomy_value = 'All'):
  """gets frequency of glycoletter in main versus side branch of glycan"""
  
  if taxonomy_value == 'All':
    glycans_df = df_species
  else:
    if taxonomy_filter == 'Kingdom':
      filter_list = pd.Series([taxonomy_value in str(k) for k in df_species['kingdom']],
                              index = df_species.index)
      glycans_df = df_species[filter_list]
    else:
      filter_list = pd.Series([taxonomy_value in str(k) for k in df_species['species']],
                              index = df_species.index)
      glycans_df = df_species[filter_list]
  glycan_list = glycans_df.glycan.values.tolist()
      
  main = 0
  side = 0
  for k in range(len(glycan_list)):
    starts = [m.start() for m in re.finditer(glycoletter, glycan_list[k])]
    init = 0
    for i in starts:
      gly = glycan_list[k][init:i]
      if '[' in gly and ']' not in gly:
        side += 1
      else:
        main += 1
      init = i
      
  return main, side

def characterize_context(glycoletter, mode = 'bond', taxonomy_filter = 'Kingdom', taxonomy_value = 'All'):
  """get characteristic microenvironment for glycoletter"""
  
  if taxonomy_value == 'All':
    glycans_df = df_species
  else:
    if taxonomy_filter == 'Kingdom':
      filter_list = pd.Series([taxonomy_value in str(k) for k in df_species['kingdom']],
                              index = df_species.index)
      glycans_df = df_species[filter_list]
    else:
      filter_list = pd.Series([taxonomy_value in str(k) for k in df_species['species']],
                              index = df_species.index)
      glycans_df = df_species[filter_list]
    
  pool = [link_find(k) for k in glycans_df.glycan.values.tolist()]
  pool = [item for sublist in pool for item in sublist]
  if mode == 'bond':
    pool = [k.split('*')[0] for k in pool if k.split('*')[1] == glycoletter]
    lab = 'Observed monosaccharides making bond %s' % glycoletter # input is a bond
  elif mode == 'sugar':
    pool = [k.split('*')[2] for k in pool if k.split('*')[0] == glycoletter]
    lab = 'Observed monosaccharides paired with %s' % glycoletter # input is a sugar
  elif mode == 'sugarbond':
    pool = [k.split('*')[1] for k in pool if k.split('*')[0] == glycoletter]
    lab = 'Observed bonds made by %s' % glycoletter # input is a sugar
  cou = Counter(pool).most_common()
  cou_k = [k[0] for k in cou if k[1]>10]
  cou_v = [k[1] for k in cou if k[1]>10]
  
  lab = lab + ' (' + taxonomy_filter + ' = ' + taxonomy_value + ')'
  
  return lab, cou_k, cou_v
