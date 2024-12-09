rm( list=ls() )
library(BEDMatrix)
load( 'Rdata/setup.Rdata' )
load( 'Rdata/setup_snps.Rdata' )

batch			<- as.numeric( commandArgs(TRUE)[[1]] )
if( batch > nbatch ) stop('batch too big' ) 

nperm	<- ifelse( batch > 0, 0, 1e2 )  ## permutation test is only for testing windows around GALA II loci -- and didn't make it to the paper as it was basically identical to just testing the GALA II SNPs directly
for( mode in c( '_pcxe','_juvenile', '_noage' ) )
	for( perm.i in 0:nperm )
		for( env in envs )
{
	if( perm.i > 0 ){
if(  mode != '_pcxe' ) next
		if( batch >= -3 ) next
		suff	<- paste0( '_', perm.i )
	} else {
		suff	<- ''
	}

	savefile	<- paste0( 'Rdata/gxewas/batch_', env, '_', batch, mode, suff, '.Rdata' )
	sinkfile	<- paste0( 'Rout/gxewas/batch_'	, env, '_', batch, mode, suff, '.Rout' )
	if( file.exists( savefile ) | file.exists( sinkfile ) )	next
	print(	sinkfile )
	sink(	sinkfile )

	datfile	<- paste0( 'Rdata/munged/', env, mode, '.Rdata' )
	if( !file.exists( datfile ) ) stop( datfile )
	load( datfile )

	if( perm.i > 0 ){
		set.seed(perm.i)
		stopifnot( ncol(Z) == 1 )
		Z[,1]	<- sample(Z[,1])
	}

	bedfile	<- ifelse( batch > 0, paste0( 'split_ukb_genos/block.', batch ), target_beds[ -1*batch ] )
	Gsnp	<- as.matrix( BEDMatrix( bedfile ) )
	iids	<- sapply( strsplit( rownames( Gsnp ), '_' ), function(x) x[[1]] )
	if( !all( IID %in% iids ) ){
		stopifnot( batch < 0 )
		names(y)		<- IID
		rownames(X)	<- IID 
		rownames(Z)	<- IID 
		IID	<- intersect( IID, iids )
		y		<- y[IID]
		X		<- X[IID,]
		Z		<- Z[IID,]
	}
	rownames(Gsnp)	<- iids 

	if( batch < 0 ){
		ncols	<- ncol(Gsnp)
		colnames_loc	<- colnames(Gsnp)
	} else {
		ncols	<- nperbatch
		colnames_loc	<- snpnames[(batch-1)*ncols+1:ncols]
	}

	beta	<- array( NA, dim=c(3,ncols), dimnames=list( c( 'hom', 'hom1', 'het' ), colnames_loc ) )
	se		<- array( NA, dim=c(3,ncols), dimnames=list( c( 'hom', 'hom1', 'het' ), colnames_loc ) )
	times	<- array(	NA, dim=c(3,ncols), dimnames=list( c( 'hom', 'het', 'null' ), colnames_loc ) )
	pv		<- array(	NA, dim=c(4,ncols), dimnames=list( c( 'hom', 'het', 'het_chisq', 'glob' ), colnames_loc ) )
	for( jj in 1:ncol(Gsnp) ){
		print( jj )
		g		<- as.numeric( Gsnp[IID,jj] )
		sub	<- which( !is.na(g) )

		times['hom',jj]	<- system.time({ fit_hom	<- glm( y ~ X + Z + g				, subset=sub, family='binomial') })[3]
		times['het',jj]	<- system.time({ fit_het	<- glm( y ~ X + Z + g + Z:g	, subset=sub, family='binomial') })[3]

		beta['hom'	,jj]	<- summary(fit_hom)$coef['g',1]
		se  ['hom'	,jj]	<- summary(fit_hom)$coef['g',2]
		pv	['hom'	,jj]	<- summary(fit_hom)$coef['g',4] 
		beta['hom1'	,jj]	<- summary(fit_het)$coef['g',1]
		se  ['hom1'	,jj]	<- summary(fit_het)$coef['g',2]
		beta['het'	,jj]	<- summary(fit_het)$coef['Z:g',1]
		se  ['het'	,jj]	<- summary(fit_het)$coef['Z:g',2]
		pv	['het'	,jj]	<- summary(fit_het)$coef['Z:g',4]

		pv	['het_chisq',jj]	<- anova( fit_hom, fit_het, test='Chisq')$Pr[2] ### almost identical to 'het'
		rm( fit_hom ); gc()

		times['null',jj]	<- system.time({ 
			pv	['glob'	,jj]	<- anova( glm( y ~ X + Z, subset=sub, family='binomial'	), fit_het, test='Chisq')$Pr[2]
		})[3] 

		rm( g, sub, fit_het ); gc()
	} 
	save( beta, se, pv, times, file=savefile )
	print(warnings()); print('Done')
	rm( y, beta, pv, X, Z, Gsnp ); gc() 
	sink()
}
