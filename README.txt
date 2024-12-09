### Set up
setup.R
setup_snps.R defines snp names, top SNPs from GALA II, etc
munge.R compiles phenotype/covariate/environment features into simple matrices, with helper functions in pull_data.R

### run GALA II SNPs + window
extract_candidate_loci.sh pulls 100kb windows around galaii SNPs
targeted_gxewas.sh calls gxewas.R (using negative indices for these special genotype chunks)
compile_targets.R compiles these candidate windows, and plot_galaII_cands.R plots them

### run GWAS/GxEWAS
gxewas.sh calls gwas.R and gxewas.R across chunks of the genome 
compile_gxewas.R concatenates output across chunks 
manhattan_figure.R plots Manhattans+QQplots, uses fxns/my_man.R and fxns/qqplot_fxns.R
gwas_forest.R reruns interactions on highlight SNPs to get per-context effects and plot them 

### run PRSxC
build_prs_eo.sh builds PRS using Price lab summstats for (whole genome + even/odd) PS
PRSxE.R fits interactions with each PRS and E, and plots results

### run GxCMM
gxemm_4e4.sh calls gxemm.R to fit different data subsets and contexts, which is plotted by plot_gxemm.R

### run CE
ce.sh calls ce.R to fit EO PRS*PRS interactions, including permutations, which plot_ce.R plots
