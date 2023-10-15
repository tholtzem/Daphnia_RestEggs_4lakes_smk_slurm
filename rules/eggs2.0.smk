localrules: ls_depth, plot_summary

rule sym_link_REFs:
  input:
    '/scratch/c7701178/mach2/DAPHNIA/Daphnia_RestEggs_snakemake_pbs_2.0_HiC/realigned/{sample}.realigned.bam'
  output:
    touch('realigned/{sample}.symlink.done')
  log: 'log/{sample}.symlink.log'
  message: """--- Create a symlink for clipped REF files ---"""
  shell:
    """
    ln -s {input} realigned/{wildcards.sample}.realigned.bam 2> {log}
    """

rule samtools_coverage:
  input:
    realigned = 'realigned/{sample}.realigned.bam'
  output:
    coverage = 'depth/{sample}.realigned.bam.coverage.hist'
  log: 'log/{sample}.realigned.bam.coverage.hist.log'
  threads: 12
  message: """ Produce an ASCII-art histogram of coverage per chromosome   """
  shell:
    """
    samtools coverage -A -o {output.coverage} {input.realigned} 2> {log}
    """


rule samtools_depth:
  input:
    realigned = 'realigned/{sample}.realigned.bam'
  output:
    depth = 'depth/{sample}.realigned.bam.depth.gz'
  log: 'log/{sample}.realigned.bam.depth.log'
  threads: 1
  resources: mem_mb=1000, walltime="2:00:00"
  message: """ Count per position depth per sample using samtools depth """
  shell:
    """
    samtools depth -aa {input.realigned} | cut -f3 | gzip > {output.depth} 2> {log}
    """


rule ls_depth:
  input:
    'depth/{sample}.realigned.bam.depth.gz'
  output:
    touch('depth/{sample}.added2DepthList.done')
  log: 'log/{sample}.depth.list.log'
  message: """ --- Creating a sample list of samtools-depth files for Rscript --- """
  shell:
    """
    ls {input} | cut -f2 -d '/' >> list/depth.list 2> {log}
    """


rule read_depth:
  input:
    'list/depth.list'
  output:
    'depth/stats/depth_statistics.tsv'
  log: 'log/depth_statistics.log'
  threads: 2
  resources: mem_mb=60000, walltime="6:00:00"
  message:
    """ --- Running Rscript to plot the genome-wide distribution of coverage --- """
  shell:
    """
    Rscript scripts/read_depth.R {input} {output} 2> {log} 
    """


rule plot_summary:
  input:
    'depth/stats/depth_statistics_{sets}.tsv'
  output:
    touch('depth/plots/plot_summary_{sets}.done')
  log: 'log/plot_summary_{sets}.log'
  message:
    """ --- Running Rscript to plot depth summary per sample and output depth filters --- """
  shell:
    """
    Rscript scripts/plot_summary.R {input} {wildcards.sets} 2> {log}
    """


#rule rbind_depthFilter:
#  output:
#    touch('depth/stats/rbind_depthFilter.done')
#  log: 'log/rbind_dfFilter.log'
#  message:
#    """ --- Running Rscript to plot depth summary per sample and output depth filters --- """
#  shell:
#    """
#    Rscript scripts/rbind_dfFilter.R depthFilter.list 2> {log}
#    """

#rule genome_coverage_bed:
#  input:
#    'realigned/{sample}.realigned.bam'
#  output:
#    'bedtools/{sample}.realigned.genomecov.bed'
#  log: 'log/{sample}.realigned.genomecov.log'
#  threads: 12
#  message:
#    """ Computes BED summaries using bedtools """
#  shell:
#    """
#    bedtools genomecov -ibam {input} > {output} 2> {log}
#    """
#
#
#rule plot_gencov:
#  input:
#    bed= 'bedtools/{sample}.realigned.genomecov.bed'
#  output:
#    pdf = 'bedtools/plots/{sample}.realigned.genomecov.pdf'
#  log: 'log/{sample}.realigned.genomecov_plot.log'
#  threads: 4
#  message:
#    """ Running Rscript to plot the genome-wide distribution of coverage """
#  shell:
#    """
#    Rscript scripts/plot_gene_covs.R {input.bed} {output.pdf} 2> {log}
#    """
