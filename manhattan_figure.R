rm( list=ls() )
source( 'fxns/qqplot_fxns.R' )
source( 'fxns/my_man.R' )
source( 'fxns/R2plink.R' )
load( 'Rdata/setup.Rdata' )
load( 'Rdata/setup_snps.Rdata' )

envs	<- c( 'disease_ALLERGY_ECZEMA_DIAGNOSED', 'eo_qn', 'neuroticism' )

sink( "Rout/gxe_manhattans.Rout" )
png(	"figs/gxe_manhattans.png", width=(5.1+1)*6, height=3*6, units= "in", res=200, pointsize=4 )
layout( matrix( 1:6, 3, 2, byrow=T ), width=c(5.1,1) ) 
par( cex.axis=2, cex.lab=2.3, cex=2.3 ) 

for( env in envs ){
	load( paste0( 'Rdata/compiled_gxewas/', env, '_pcxe.Rdata' ) ) 

	candfile	<- paste0( 'Rdata/compiled_gxewas/', env, '_pcxe_cands.Rdata' )
	if( ! file.exists( candfile ) ){
		cand_rs_raw	<- snpnames[which( allpv['hom0',]< 5e-8 )] 
		cand_rs			<- ldclump( cand_rs_raw, allpv['hom0',cand_rs_raw], r2=.1 ) 
		save(cand_rs_raw,cand_rs,file=candfile)
	} else {
		load(candfile)
	}
	cat( 'GWAS hits:', length( cand_rs ), '\n' )

	jointfile	<- paste0( 'Rdata/compiled_gxewas/', env, '_pcxe_joint.Rdata' )
	if( ! file.exists( jointfile ) ){
		joint_rs_raw	<- snpnames[which( allpv['glob',]< 5e-8 )] 
		joint_rs			<- ldclump( joint_rs_raw, allpv['glob',joint_rs_raw], r2=.1 ) 
		print( table( joint_rs %in% cand_rs_raw ) )
		print( quantile( ( allpv['hom0',joint_rs] )[! joint_rs %in% cand_rs_raw ] ) )

		### clump against additive GWAS
		joint_rs			<- setdiff( joint_rs, cand_rs_raw ) 
		joint_rs_new	<- ldclump( c(joint_rs, cand_rs_raw), rep(c(.5,1e-12),c(length(joint_rs),length(cand_rs_raw))), r2=.1 )
		joint_rs_new	<- intersect( joint_rs, joint_rs_new )
		save(joint_rs_raw,joint_rs,joint_rs_new,file=jointfile)
	} else {
		load(jointfile)
	}
	cat( 'Joint G+GXE GWAS hits:', joint_rs_new, '\n\n' )

	try({
	cat( 'Binom test for ', env, ':', 
		binom.test( sum( allpv['het' ,cand_rs] > -log10(.05) ), length( cand_rs ), alternative="greater", p=.05)$p.value, '\n' )
	})

	cat( 'GxEWAS hits for ', env, ':\n' )
	print( allpv[,which( allpv['het' ,] < 5e-8 )] )

	par( mar=c(8,12,1,1.2) )
	my_man( allpv['het' ,], chrnames, bpnames, snpnames, chr="CHR", bp="BP", snp="SNP", p="P", suggestiveline=FALSE, highlight=cand_rs, cex.highlight=3.0 ) 
	mtext( side=2, srt=90, niceEnvs[env], cex=10, line=8 )

	par( mar=c(8,8,1,1) )
	my_qq( y=-log10( allpv['het',] ), ylab='Interaction log10(p)', cand_rs=list( cand_rs ), cex.lab=2.3 ) 
	if( env == envs[1] )
		legend( 'topleft', bty='n', fill=c(2,5,4,3), leg=c( 'GWAS', 'Bonferroni', 'FDR<.1', 'N.S.' ), cex=2 )

	rm( cand_rs, allbeta, allpv ) 
}
dev.off() 

png(	 "figs/qqplots.png", width=5*(3+2/5), height=5, units= "in", res=400, pointsize=14 )
layout( matrix( 1:4, 1, 4, byrow=T ), width=c(5,5,5,2.3) ) 
par( cex.axis=1.3, cex.lab=1.5 )
par( mar=c(4.5,4.5,4,2) )

for( env in envs ){
load( paste0( 'Rdata/compiled_gxewas/', env, '_pcxe.Rdata' ) ) 
load( paste0( 'Rdata/compiled_gxewas/', env, '_pcxe_cands.Rdata' ) ) 

highlights	<- list( c( 'rs532965' ), c( 'rs1837253','rs34645399', 'rs146125856','rs72823646','rs12123821','rs12964116','rs1701704') )
my_qq( y=-log10( allpv['het',] ), ylab='Interaction log10(p)', cand_rs=list( cand_rs ), cex.main=1.9, main=paste0( 'Gx', niceEnvs[env] ), highlights=highlights )
rm( cand_rs, allbeta, allpv ) 
}

par( mar=c(4,0,4,0) )
plot.new()

legend( 'top'		, bty='n', title='All SNPs'								, col=c(2,'grey'), leg=c( 'GWAS Hit', 'N.S.' )											, pch=16, cex=1.8, pt.cex=c(2,1) )
legend( 'bottom', bty='n', title='Hom. Asthma\nGWAS hits'	, col=c(5,4,3), leg=c( 'Bonferroni', 'FDR<.1', 'N.S.' )	, pch=16, cex=1.8, pt.cex=c(2,2,1) )

dev.off() 


### clump all new hits:
allsnps <- NULL 
for( env in envs ){
	load( paste0( 'Rdata/compiled_gxewas/', env, '_pcxe_joint.Rdata' ) )
	allsnps			<- c( allsnps, joint_rs_new )
	rm( joint_rs_new )
}
out <- ldclump( allsnps, runif(length(allsnps),0,.1), r2=.1 )
length( allsnps )
length( out )
out

sink()
