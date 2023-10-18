include: "rules/common.smk"



# -----------------------------------------------

rule all:
	input:
		expand('realigned/{sample}.symlink.done', sample=samples_subset_prefix),
		#expand('depth/{sample}.{ext}', sample=samples_ALL_prefix, ext=['realigned.bam.coverage.hist']),
		expand('depth/{sample}.{ext}', sample=samples_ALL_prefix, ext=['realigned.bam.depth.gz']),
		expand('depth/{sample}.added2DepthList.done', sample=Dlgc_names),
		expand('depth/stats/depth_statistics_{sets}.tsv', sets=['DLGC']),
		expand('depth/plots/plot_summary_{sets}.done', sets=['DLGC']),
		expand('angsd/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done', zip, sets=['DLGC'], GL=['2'], minMaf=['0.05'], IND=N, MinDepth=MinDepth, MaxDepth=MaxDepth),
		expand('angsd/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_globalSNP.list', zip, sets=['DLGC'], GL=['2'], minMaf=['0.05'], IND=N, MinDepth=MinDepth, MaxDepth=MaxDepth),
		expand('angsd/{sets}/index_globalSNP_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.done', zip, sets=['DLGC'], GL=['2'], minMaf=['0.05'], IND=N, MinDepth=MinDepth, MaxDepth=MaxDepth),
		expand('angsd/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_globalSNP.chr', zip, sets=['DLGC'], GL=['2'], minMaf=['0.05'], IND=N, MinDepth=MinDepth, MaxDepth=MaxDepth),
		expand('pcangsd/{sets}/PCAngsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_cov_admix.done', zip, sets=['DLGC'], GL=['2'], minMaf=['0.05'], IND=N, MinDepth=MinDepth, MaxDepth=MaxDepth),
		expand('LD_decay/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_geno.beagle.gz', zip, sets=['DLGC'], GL=['2'], minMaf=['0.05'], IND=N, MinDepth=MinDepth, MaxDepth=MaxDepth),
		expand('LD_decay/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_pos.gz', zip, sets=['DLGC'], GL=['2'], minMaf=['0.05'], IND=N, MinDepth=MinDepth, MaxDepth=MaxDepth),
		expand('LD_decay/{sets}/ngsLD_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{kb}kb_dist.ld.gz', zip, sets=['DLGC'], GL=['2'], minMaf=['0.05'], IND=N, MinDepth=MinDepth, MaxDepth=MaxDepth, kb=['0', '100']),
		expand('LD_decay/{sets}/ngsLD_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{kb}kb_dist.LDdecay.pdf', zip, sets=['DLGC'], GL=['2'], minMaf=['0.05'], IND=N, MinDepth=MinDepth, MaxDepth=MaxDepth, kb=['0', '100'])
		#expand('LD_deacy/{sets}/run_LDpruning_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{kb}kb_dist_{Minweight}Minweight.done', zip, sets=sets_structure, GL=GL_structure, minMaf=minMaf, IND=N_structure, MinDepth=MinDepth_structure, MaxDepth=MaxDepth_structure, Minweight=['0.1', '0.1', '0.1', '0.1']),
		#expand('ngsLD/{sets}/LDpruned_snps_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{Minweight}Minweight.list', zip, sets=['DLGC'], GL=['2'], minMaf=['0.05'], IND=N, MinDepth=MinDepth, MaxDepth=MaxDepth, kb=['0', '100'], Minweight=['0.1', '0.5']),
		#expand('ngsLD/{sets}/index_SNPs_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{Minweight}Minweight.done', zip, sets=sets_structure, GL=GL_structure, minMaf=minMaf, IND=N_structure, MinDepth=MinDepth_structure, MaxDepth=MaxDepth_structure, Minweight=['0.5', '0.5', '0.5', '0.5']),
		#expand('ngsLD/{sets}/LDpruned_snps_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{Minweight}Minweight.chr', zip, sets=sets_structure, GL=GL_structure, minMaf=minMaf, IND=N_structure, MinDepth=MinDepth_structure, MaxDepth=MaxDepth_structure, Minweight=['0.5', '0.5', '0.5', '0.5']),
		#expand('angsd/{sets}/LDpruned_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{Minweight}Minweight.done', zip, sets=sets_structure, GL=GL_structure, minMaf=minMaf, IND=N_structure, MinDepth=MinDepth_structure, MaxDepth=MaxDepth_structure, Minweight=['0.1', '0.1', '0.1', '0.1']),
		#expand('angsd/{sets}/LDpruned_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{Minweight}Minweight.done', zip, sets=sets_structure, GL=GL_structure, minMaf=minMaf, IND=N_structure, MinDepth=MinDepth_structure, MaxDepth=MaxDepth_structure, Minweight=['0.5', '0.5', '0.5', '0.5']),
		#expand('pcangsd/{sets}/PCAngsd_LDpruned_angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_{Minweight}Minweight.done', zip, sets=sets_structure, GL=GL_structure, minMaf=minMaf, IND=N_structure, MinDepth=MinDepth_structure, MaxDepth=MaxDepth_structure, Minweight=['0.5', '0.5', '0.5', '0.5']),
		#expand('LD_decay/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_geno.beagle.gz', sets=['LC'], GL=['2'], minMaf=['0.05'], IND=['202'], MinDepth=['202'], MaxDepth=['9051']),
		#expand('LD_decay/{sets}/angsd_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}_pos.gz', sets=['LC'], GL=['2'], minMaf=['0.05'], IND=['202'], MinDepth=['202'], MaxDepth=['9051']),
		#expand('LD_decay/{sets}/ngsLD_windows_GL{GL}_minInd{IND}_maf{minMaf}_minDepth{MinDepth}_maxDepth{MaxDepth}.ld.gz', sets=['LC'], GL=['2'], minMaf=['0.05'], IND=['202'], MinDepth=['202'], MaxDepth=['9051']),
                #expand('saf/saf_samples/{sample}.GL2.saf.idx.done', sample=ALL_names),
                #expand('saf/saf_samples/{sample}.GL2.est.ml', sample=ALL_names),
		#expand('list/saf/{POPS}.list', POPS=sets_saf),
		#expand('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.saf.idx.done', zip, POPS=sets_saf, IND=N_saf, MinDepth=MinDepth_saf, MaxDepth=MaxDepth_saf),
		#expand('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.1D.sfs', zip, POPS=sets_saf, IND=N_saf, MinDepth=MinDepth_saf, MaxDepth=MaxDepth_saf),
		#expand('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.sfs2theta.done', zip, POPS=sets_saf, IND=N_saf, MinDepth=MinDepth_saf, MaxDepth=MaxDepth_saf),
		#expand('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.theta_stats_chrom.done', zip, POPS=sets_saf, IND=N_saf, MinDepth=MinDepth_saf, MaxDepth=MaxDepth_saf),
		#expand('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.theta_stats_Window10kb.done', zip, POPS=sets_saf, IND=N_saf, MinDepth=MinDepth_saf, MaxDepth=MaxDepth_saf),
		#expand('saf/POPS/{POP1}_vs_{POP2}_2D.sfs', zip, POP1=POP1, POP2=POP2),
		#expand('saf/POPS/2pops/{POP1}_vs_{POP2}.fstIndex.done', zip, POP1=POP1, POP2=POP2),
		#expand('saf/POPS/3pops/{POP1}_vs_{POP2}_vs_{POP3}.fstIndex.done', zip, POP1='PRE_long_LC', POP2='ALL_gal_LC', POP3='PEL_long_LC'),
		#expand('saf/POPS/2pops/{POP1}_vs_{POP2}.fst.global.tsv', zip, POP1=POP1, POP2=POP2),
		#expand('saf/POPS/2pops/{POP1}_vs_{POP2}.fst.Window10kb.type{types}.tsv', zip, POP1=POP1, POP2=POP2, types=['1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1']),
		#expand('saf/POPS/3pops/{POP1}_vs_{POP2}_vs_{POP3}.fst.global.tsv', POP1='PRE_long_LC', POP2='ALL_gal_LC', POP3='PEL_long_LC'),
		#expand('saf/POPS/3pops/{POP1}_vs_{POP2}_vs_{POP3}.fst.Window10kb.type{types}.tsv', POP1='PRE_long_LC', POP2='ALL_gal_LC', POP3='PEL_long_LC', types=['1', '2', '3']),
		#expand('saf/POPS/2pops/{POP1}_vs_{POP2}.A_B.tsv.gz', zip, POP1=POP1, POP2=POP2),
		#expand('saf/POPS/3pops/{POP1}_vs_{POP2}_vs_{POP3}.A_B.tsv.gz', POP1='PRE_long_LC', POP2='ALL_gal_LC', POP3='PEL_long_LC'),
		#expand('ngsRelate/{POPS}.glf3.freq', POPS=['PEL_cuc_LC', 'PEL_long_LC', 'PRE_long_LC', 'PEL_F1_LC']),
		#expand('ngsRelate/{POPS}.glf3.freq', POPS=sets_saf),
		#expand('ngsRelate/{POPS}.glf3.stats.txt', POPS=sets_saf),
		#expand('mafs/{POPS}_1000random_SNPs.mafs.done', POPS=['PRE_long_LC', 'ALL_gal_LC', 'PEL_long_LC']),
		#expand('mafs/{POPS}_1000random_sites.mafs.done', POPS=['PRE_long_LC', 'ALL_gal_LC', 'PEL_long_LC']),
		#expand('mafs/{POPS}_1000random_fst_outliers.mafs.done', POPS=['PRE_long_LC', 'ALL_gal_LC', 'PEL_long_LC']),
		#expand('mafs/{POPS}_top10_NOTfixed_fst_outliers.mafs.done', POPS=['PRE_long_LC', 'ALL_gal_LC', 'PEL_long_LC']),
		#expand('mafs/{POPS}_10random_fixed_fst_outliers.mafs.done', POPS=['PRE_long_LC', 'ALL_gal_LC', 'PEL_long_LC']),
		#expand('list/abbababa/4pop_combinations.bam.list.done'),
		#expand('list/abbababa/4pop_combinations.sizefile.done'),
		#expand('abbababa/{prefix}/{combi}.abbababa.{blocksize}blocksize.done', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], combi=Dstats_combi[4:12], blocksize=['2500000']),
		#expand('list/abbababa/errorFile.txt'),
		#expand('list/abbababa/4pop_combinations.Popnames.list.done'),
		#expand('abbababa/{prefix}/{combi}.Dstat_results.{blocksize}blocksize.done', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], combi=Dstats_combi[4:12], blocksize=['2500000']),
		##expand('realigned/{sample}.{ext}', sample=samples_ALL_names, ext=['realigned.bai']),
		#expand('saf/POPS/sites/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.sites.saf.idx.done', zip, POPS=sets_saf, IND=N_saf, MinDepth=MinDepth_saf, MaxDepth=MaxDepth_saf),
		#expand('saf/POPS/sites/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.sites.1D.sfs', zip, POPS=sets_saf, IND=N_saf, MinDepth=MinDepth_saf, MaxDepth=MaxDepth_saf),
		#expand('saf/POPS/sites/{POP1}_vs_{POP2}.sites.2D.sfs', zip, POP1=POP1, POP2=POP2),
		#expand('saf/POPS/sites/2pops/{POP1}_vs_{POP2}.sites.fstIndex.done', zip, POP1=POP1, POP2=POP2),
		#expand('saf/POPS/sites/2pops/{POP1}_vs_{POP2}.sites.fst.global.tsv', zip, POP1=POP1, POP2=POP2),
		#expand('saf/POPS/sites/3pops/{POP1}_vs_{POP2}_vs_{POP3}.sites.fstIndex.done', zip, POP1='PRE_long_LC', POP2='ALL_gal_LC', POP3='PEL_long_LC'),
		#expand('saf/POPS/sites/3pops/{POP1}_vs_{POP2}_vs_{POP3}.sites.fst.global.tsv', POP1='PRE_long_LC', POP2='ALL_gal_LC', POP3='PEL_long_LC'),
		#expand('saf/POPS/sites/3pops/{POP1}_vs_{POP2}_vs_{POP3}.sites.fst.Window10kb.type{types}.tsv', POP1='PRE_long_LC', POP2='ALL_gal_LC', POP3='PEL_long_LC', types=['1', '2', '3']),
		#expand('saf/POPS/sites/3pops/{POP1}_vs_{POP2}_vs_{POP3}.sites.A_B.tsv.gz', POP1='PRE_long_LC', POP2='ALL_gal_LC', POP3='PEL_long_LC'),
		#expand('synthesis_daphnia/list/PRE_long_LC_ALL_gal_LC_PEL_long_LC.SNPs.chr'),
		#expand('synthesis_daphnia/list/PRE_long_LC_ALL_gal_LC_PEL_long_LC.SNPs_index.done'),
		#expand('synthesis_daphnia/mafs/{bams}.mafs.done', bams=bams),
		#expand('synthesis_daphnia/mafs/df{df}.mafs.done', df=range(1,19)),
		#expand('ancestry/LC/data/get_bcf_{prefix}.done', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051']),
		#expand('ancestry/LC/data/{prefix}.{ext}', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], ext=['vcf.gz']),
		#expand('ancestry/LC/data/{prefix}_{ext}', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], ext=['INFODP.txt', 'nbrSites.txt']),
		#expand('ancestry/LC/data/{prefix}_{ext}', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], ext=['INFODP.pdf', 'INFODP.stats.txt']),
		#expand('ancestry/LC/data/{prefix}_hf.{ext}', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], ext=['vcf.gz', 'nbrSites.txt']),
		#expand('ancestry/LC/data/{prefix}_hf.DP{genoDP}.{ext}', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], genoDP=['10'], ext=['vcf.gz', 'nbrSites.txt']),
		#expand('ancestry/LC/data/{prefix}_hf.DP{genoDP}.{ext}', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], genoDP=['10'], ext=['imiss']),
		#expand('ancestry/LC/data/{prefix}_hf.DP{genoDP}.imiss.{ext}', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], genoDP=['10'], ext=['pdf', 'tsv']),
		#expand('ancestry/LC/data/imissRM/{prefix}_hf.DP{genoDP}.imissRM.{ext}', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], genoDP=['10'], ext=['vcf.gz']),
		#expand('ancestry/LC/data/imissRM/{prefix}_hf.DP{genoDP}.imissRM.vmiss20.{ext}', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], genoDP=['10'], ext=['vcf', 'nbrSites.txt']),
		#expand('ancestry/LC/all/{pop1}_{pop2}.{ext}', pop1=['PRE_long_LC'], pop2=['galeata9'], ext=['fixed_sites.txt', 'fixed_sites.report.tsv']),
		#expand('ancestry/LC/all/{pop1}_{pop2}.{ext}', pop1=['PRE_long_LC'], pop2=['galeata9'], ext=['fixed_sites_100percent_thinned1000bp.svg', 'fixed_sites_100percent_thinned1000bp.report.tsv']),
		#expand('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.saf.idx.done', zip, POPS=['cucREF6'], IND=['6'], MinDepth=['6'], MaxDepth=['284']),
		#expand('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.1D.sfs', zip, POPS=['cucREF6'], IND=['6'], MinDepth=['6'], MaxDepth=['284']),
		#expand('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.sfs2theta.done', zip, POPS=['cucREF6'], IND=['6'], MinDepth=['6'], MaxDepth=['284']),
		#expand('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.theta_stats_chrom.done', zip, POPS=['cucREF6'], IND=['6'], MinDepth=['6'], MaxDepth=['284']),
		#expand('saf/POPS/{POPS}_{IND}_{MinDepth}_{MaxDepth}.GL2.theta_stats_Window10kb.done', zip, POPS=['cucREF6'], IND=['6'], MinDepth=['6'], MaxDepth=['284']),
		#expand('list/abbababa/4pop_combinations.bam.list.done'),
		#expand('list/abbababa/4pop_combinations.sizefile.done'),
		#expand('abbababa/{prefix}/{combi}.abbababa.{blocksize}blocksize.done', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], combi=Dstats_combi2, blocksize=['2500000']),
		#expand('list/abbababa/errorFile.txt'),
		#expand('list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_globalSNP.chr'),
		#expand('list/angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051_globalSNP_index.done'),
		#expand('mafs/{bams}.mafs.done', bams=bams),
		#expand('ancestry/LC/unfiltered/{pop1}_{pop2}.{ext}', pop1=['PRE_long_LC'], pop2=['galeata9'], ext=['fixed_sites.txt', 'fixed_sites.report.tsv']),
		#expand('ancestry/LC/notfiltered/{pop1}_{pop2}.{ext}', pop1=['PRE_long_LC'], pop2=['galeata9'], ext=['fixed_sites_100percent_thinned1000bp.svg', 'fixed_sites_100percent_thinned1000bp.report.tsv']),
		#expand('list/abbababa/4pop_combinations.Popnames.list.done'),
		#expand('abbababa/{prefix}/{combi}.Dstat_results.{blocksize}blocksize.done', prefix=['angsd_GL2_minInd202_maf0.05_minDepth202_maxDepth9051'], combi=Dstats_combi2, blocksize=['2500000']),


# -----------------------------------------------


include: "rules/eggs2.0.smk"
include: "rules/get_beagle_allSNPs.smk"
include: "rules/PCAngsd_allSNPs.smk"
include: "rules/ngsLD.smk"
#include: "rules/get_beagle_LDprunedSNPs.smk"
#include: "rules/PCAngsd_LDprunedSNPs.smk"
#include: "rules/estimate_ngsLD.smk"
#include: "rules/realSFS.smk"
#include: "rules/get_allelefreq_perpop.smk"
#include: "rules/ngsRelate.smk"
#include: "rules/realSFS_Fst_sites.smk"
#include: "rules/get_allelefreq_perpop_sites.smk"
#include: "rules/abbababa2.smk"
#include: "rules/get_vcf.smk"
#include: "rules/ancestry.smk"
