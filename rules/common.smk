import os
import pandas as pd

###### Config file and sample sheets #####


configfile: 
  "config/config.yaml"

# load sample info
samples_ALL = pd.read_csv("list/Da_Resteggs_4lakes.metadata.csv", sep=',', index_col=False)
samples_ALL_prefix = list(samples_ALL['prefix'])
# extract ref clones, except outgroup
samples_subset = samples_ALL[(samples_ALL['location']!='VAR') & (samples_ALL['location']!='ANN')]
samples_subset_prefix = list(samples_subset['prefix'])
# extract ref clones, except outgroup
#refclones = samples_ALL[(samples_ALL['period']=='REF') & (samples_ALL['species']!='curvirostris')]
#refclone_names = list(refclones['prefix'])
# load sample prefix from duplicate removed bam files, except outgroup 
Dlgc = samples_ALL[samples_ALL['species']!='curvirostris']
Dlgc_names = list(Dlgc['prefix'])
# ougroup info
#outgroup = samples_ALL[samples_ALL['species']=='curvirostris']
#outgroup_names = list(outgroup['prefix'])

# load depth filters
if os.path.isfile("depth/stats/DLGC_depthFilter.list"):
  print ("Depth filter file exists")
  depth_information = pd.read_csv("depth/stats/DLGC_depthFilter.list", sep='\t')
  #depth_structure = depth_information[(depth_information['pops']=='LC') | (depth_information['pops']=='LZ') | (depth_information['pops']=='LCwithoutREF') | (depth_information['pops']=='LZwithoutREF')]
  #depth_saf = depth_information[(depth_information['pops']!='LC') & (depth_information['pops']!='LZ') & (depth_information['pops']!='LCwithoutREF') & (depth_information['pops']!='LZwithoutREF')]
        

  # number of samples (individuals)
  N = list(depth_information['Number_of_samples'])
  #N_structure = list(depth_structure['Number_of_samples'])
  #N_saf = list(depth_saf['Number_of_samples'])
  # minimum depth over all samples: 1 x number of samples
  MinDepth = 1*N
  #MinDepth_structure = 1*N_structure
  #MinDepth_saf = 1*N_saf
  # maximum over all samples: (mean depth + 3 x standard deviation) x number of samples
  #MaxDepth = int(depth_information['HengLi_max']*N)
  MaxDepth = [int(max) for max in depth_information['HengLi_max']*depth_information['Number_of_samples']]
  #MaxDepth_structure = [int(max) for max in depth_structure['HengLi_max']*depth_structure['Number_of_samples']]
  #MaxDepth_saf = [int(max) for max in depth_saf['HengLi_max']*depth_saf['Number_of_samples']]

  #sets = list(depth_information['pops'])
  #sets_structure = list(depth_structure['pops'])
  #sets_saf = list(depth_saf['pops'])
else:
  print ("Depth file does not exist")

GL_structure = ['2', '2', '2', '2']
minMaf = ['0.05', '0.05', '0.05', '0.05']


## angsd parameters
# minor allele frequency for PCA and admixture plots for 1% and 5% of the data
#minMaf = ['0.01', '0.05']
# minor allele frequency for SFS, keeping singletons
# SFS
#minMaf_SFS = round(2/(2*N),3)

# load minimum number of individuals a read has to be present
#nInd = pd.read_csv("list/cutoff_nInd.txt", sep='\t')
# in 50 %, 75 % and 100 % of individuals
#IND = list(nInd['nInd'])


# number of Ks for admixture proportions
#admix_K = ['2', '2', '3', '3', '4', '4', '5', '5', '6', '6', '7', '7', '8', '8', '9', '9', '10', '10']
admix_K = ['2', '3', '4', '5', '6']

POP = ['longispina_March21', 'longispina_June21']

## Population pairs for 2Dsfs and fst
if os.path.isfile("list/pop_pairs.list"):
  print ("Population pair file for 2Dsfs and fst exists")
  pop_pair = pd.read_csv("list/pop_pairs.list", sep='\t', index_col=False)
  POP1 = list(pop_pair['POP1'])
  POP2 = list(pop_pair['POP2'])
else:
  print ("Population pair file for 2Dsfs and fst does not exist")


## 4-population combinations for ABBA-BABA (Dstats)
if os.path.isfile("list/abbababa/4pop_combinations.list"):
  print ("Population combination file for Dstats exists")
  pop_combi = pd.read_csv("list/abbababa/4pop_combinations.list", sep=',', index_col=False)
  pop_combi['combi'] = pop_combi[['P1', 'P2', 'P3', 'outgroup']].apply(lambda x: '_'.join(x), axis=1)
  Dstats_combi = list(pop_combi['combi'])
  pop_combi2 = pd.read_csv("list/abbababa/4pop_combinations2.list", sep=',', index_col=False)
  pop_combi2['combi'] = pop_combi2[['P1', 'P2', 'P3', 'outgroup']].apply(lambda x: '_'.join(x), axis=1)
  Dstats_combi2 = list(pop_combi2['combi'])
else:
  print ("Population combination file for Dstats does not exist")



