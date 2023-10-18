localrules: genoFile, posFile, LDdecay, LDpruned_SNPlist

rule genoFile:
  input:
    'angsd/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done'
  output:
    'LD_decay/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_geno.beagle.gz'
  log: 'log/{sets}/prepare_subgenoFile_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.log'
  threads: 12
  message:
    """ Prepare beagle formatted genotype likelihood file generated from ANGSD (-doGlf 2) by removing the header row and the first three columns (i.e. positions, major allele, minor allele) """
  shell:
    """
    zcat angsd/{wildcards.sets}/angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}.beagle.gz | cut -f 4- | awk 'NR != 1'| gzip  > {output} 2> {log}
    """


rule posFile:
  input:
    'angsd/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done'
  output:
    'LD_decay/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_pos.gz'
  log: 'log/{sets}/prepare_posFile_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.log'
  message:
    """ Prepare position file generated from ANGSD (-doGlf 2) by removing the header row and the first three columns (i.e. positions, major allele, minor allele) """
  shell:
    """
    zcat angsd/{wildcards.sets}/angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}.mafs.gz | cut -f 1,2 |  awk 'NR != 1' | sed 's/:/_/g'| gzip > {output}
    """


rule estimate_ngsLD:
  input:
    position = 'LD_decay/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_pos.gz',
    geno = 'LD_decay/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_geno.beagle.gz'
  output:
    'LD_decay/{sets}/ngsLD_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{kb}kb_dist.ld.gz'
  log: 'log/{sets}/ngsLD_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{kb}kb_dist.ld.log'
  threads: 2
  resources: mem_mb=120000, walltime="120:00:00"
  message:
    """ Estimate LD using ngsLD, which employs a maximum likelihood approach to account for genotype uncertainty in low-coverage whole genome sequencing data. """
  shell:
    """
    module load singularity/3.8.7-python-3.10.8-gcc-8.5.0-e6f6onc
    N_SITES=`zcat {input.position} | wc -l` &&
    singularity exec --home $PWD:$HOME /scratch/c7701178/bio/ngs+tools.sif ngsLD --geno {input.geno} --pos {input.position} --probs --n_ind {wildcards.IND} --n_sites $N_SITES --max_kb_dist {wildcards.kb} --max_snp_dist 0 --min_maf 0.05 --rnd_sample 0.01 --n_threads {threads} | gzip --best > {output} 2> {log}
    """


rule LDdecay:
  input:
    'LD_decay/{sets}/ngsLD_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{kb}kb_dist.ld.gz'
  output:
    'LD_decay/{sets}/ngsLD_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{kb}kb_dist.LDdecay.pdf'
  message:
    """ Estimate LDdecay """
  shell:
    """
    module load singularity/3.8.7-python-3.10.8-gcc-8.5.0-e6f6onc
    singularity exec --home $PWD:$HOME /scratch/c7701178/bio/ngs+tools.sif Rscript --vanilla --slave /opt/ngsTools/ngsLD/scripts/fit_LDdecay.R --ld_files ld_files.list --out {output} --fit_level 3 --plot_data --ld=r2 --n_ind=330 --max_kb_dist=10000 --plot_size=2,4
    """


rule run_LDpruning:
  input:
    ld = 'LD_decay/{sets}/ngsLD_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{kb}kb_dist.ld.gz',
    decay = 'LD_decay/{sets}/ngsLD_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{kb}kb_dist.LDdecay.pdf'
  output:
    touch('LD_decay/{sets}/run_LDpruning_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{kb}kb_dist_{Minweight}Minweight.done')
  log: 'log/{sets}/run_LDpruning_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{kb}kb_dist_{Minweight}Minweight.log'
  threads: 2
  resources: mem_mb=120000, walltime="120:00:00"
  message:
    """ Prune your data, remove SNPs with high LD """
  shell:
    """
     module load singularity/3.8.7-python-3.10.8-gcc-8.5.0-e6f6onc 
     zcat {input} | singularity exec --home $PWD:$HOME /scratch/c7701178/bio/ngs+tools.sif /opt/ngsLD/scripts/prune_graph.pl --max_kb_dist 2500 --min_weight {wildcards.Minweight} --out ngsLD/{wildcards.sets}/angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}_{wildcards.Minweight}Minweight_sub_unlinked.id --print_excl ngsLD/{wildcards.sets}/angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}_{wildcards.Minweight}Minweight_sub_excluded_nodes.csv
    """


rule LDpruned_SNPlist:
  input:
    touched = 'LD_decay/{sets}/run_LDpruning_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{kb}kb_dist_{Minweight}Minweight.done'
  output:
    arg3 = 'LD_decay/{sets}/LDpruned_snps_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{kb}kb_dist_{Minweight}Minweight.list'  
  log: 'log/{sets}/LDpruned_SNPlist_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{kb}kb_dist_{Minweight}Minweight.log'
  message:
    """ Make a list of LDpruned SNPs in R and index SNP list """
  shell:
    """
    Rscript scripts/getLDpruned_SNPlist.R ngsLD/{wildcards.sets}/angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}_{wildcards.Minweight}Minweight_sub_unlinked.id angsd/{wildcards.sets}/angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}.mafs.gz {output.arg3} 2> {log}
    """


#rule LDblocks:
#  input:
#   'ngsLD/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_sub.ld.gz'
#  output:
#   touch('ngsLD/{sets}/run_LDpruning_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done')
#  log: 'log/{sets}/LD_pruned_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.log'
#  threads: 12
#  message:
#    """ Prune your data, remove SNPs with high LD """
#  shell:
#    """
#    module load singularity/2.x
#    zcat {input} | singularity exec /apps/uibk/ngsld/1.1.1/ngsLD.sandbox /opt/ngsLD/scripts/prune_graph.pl --max_kb_dist 2000 --min_weight 0.5 --out ngsLD/{wildcards.sets}/angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}_sub_unlinked.id --print_excl ngsLD/{wildcards.sets}/angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}_sub_excluded_nodes.csv
#    """


