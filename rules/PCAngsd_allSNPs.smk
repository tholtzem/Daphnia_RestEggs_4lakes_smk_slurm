
rule PCAngsd:
  input: 
    touched = 'angsd/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done' 
  output:
    touch('pcangsd/{sets}/PCAngsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_cov_admix.done')
  log: 'log/{sets}/PCAngsd_covmat_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_cov_admix.log'
  threads: 2
  resources: mem_mb=120000, walltime="48:00:00"
  message:
    """ Estimate covariance matrix and admixture proportions from GLs using PCAngsd """
  shell:
    """
    module load singularity/3.8.7-python-3.10.8-gcc-8.5.0-e6f6onc
    for i in {{2..10}}; do
    	singularity exec --home $PWD:$HOME /scratch/c7701178/bio/ngs+tools.sif pcangsd -beagle angsd/{wildcards.sets}/angsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth}.beagle.gz -o pcangsd/{wildcards.sets}/PCAngsd_GL{wildcards.GL}_minInd{wildcards.IND}_maf{wildcards.minMaf}_minDepth{wildcards.MinDepth}_maxDepth{wildcards.MaxDepth} -minMaf {wildcards.minMaf} -admix -admix_K $i -e 10 -threads {threads} 2> {log}
    done
    """
