rm( list=ls() )
load( 'Rdata/setup.Rdata' )
load( 'Rdata/setup_snps.Rdata' )

nperm	<- 100 
envs	<- c( 'disease_ALLERGY_ECZEMA_DIAGNOSED', 'eo_qn', 'neuroticism' )
outps	<- array( NA, dim=c(4,3,length(keysnps)), dimnames=list( c('snp','snp-perm','best-add','locus'), envs, keysnps ) )
mode	<- '_pcxe'

batches <- -(11:17)
for( batch in batches )
		for( env in envs )
{
	for( perm.i in 0:nperm ){
		suff	<- ifelse( perm.i > 0, paste0( '_', perm.i ), '' ) 
		savefile	<- paste0( 'Rdata/gxewas/batch_', env, '_', batch, mode, suff, '.Rdata' )
		load( savefile )

		if( perm.i == 0 ){
			pvs		<- array( NA, dim=c(nperm+1,ncol(pv)) )
			betas	<- beta
			colnames(betas)<-
			colnames(pvs)	<- sapply( strsplit( colnames( pv ), '_' ), function(x) x[[1]] )
			best_add <- colnames(pvs)[which.min( pv['hom',] )]
			best_het <- colnames(pvs)[which.min( pv['het',] )]
			print( cor( pv['het',], pv['het_chisq',] ) )
			print( cor( -log10(pv['het',]), -log10(pv['het_chisq',]) ) )
		}
		pvs[perm.i+1,]	<- pv['het',]
		rm( pv, beta ) 
	}

	minps	<- apply( pvs, 1, min ) 
	outps['snp'  ,env,keysnps[-batch]]	<- pvs[ 1,keysnps[-batch]]
	outps['locus',env,keysnps[-batch]]	<- mean( minps <= minps[1], na.rm=T )
	outps['snp-perm'  ,env,keysnps[-batch]]	<- mean( pvs[,keysnps[-batch]] <= pvs[ 1,keysnps[-batch]], na.rm=T )
	outps['best-add'  ,env,keysnps[-batch]]	<- pvs[ 1,best_add ]

	rm( pvs, minps ) 
}
save( outps, file='Rdata/cand_pvals_compiled.Rdata' )


pdf( 'figs/target/snp_vs_locus.pdf', width=12, height=12 ) 
par( mfrow=c(2,2) ) 

plot1	<- function(x,y,dx=.009,dy=.005,srt=65, lims=range(c(x,y),na.rm=T),ylab='100kb locus (perm test)',...){
	plot( x, y, xlim=lims, ylim=lims, xlab='Cand SNP from GALAII', ylab=ylab, ... )
	abline( col=2, v=.01 )
	abline( col=2, h=.01 )
	abline( col=2, v=.05 )
	abline( col=2, h=.05 ) 
}
plot1(outps['snp',,], outps['locus',,]   ,.03,.09,srt=70,main='', lims=0:1 )
plot1(outps['snp',,], outps['best-add',,],.03,.09,srt=70,main='', lims=0:1, ylab='Best Add SNP' )
plot1(outps['snp',,], outps['snp-perm',,],.03,.09,srt=70,main='', lims=0:1, ylab='Permtest' )
plot(-log10(outps['snp',,]), -log10(outps['snp-perm',,]),xlim=c(0,3),ylim=c(0,3)); abline( a=0, b=1 )

dev.off() 
