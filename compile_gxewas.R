rm( list=ls() )
load( 'Rdata/setup.Rdata' )
load( 'Rdata/setup_snps.Rdata' ) 

for( env in envs ){ 
	compfile	<- paste0( 'Rdata/compiled_gxewas/'	, env, '_pcxe.Rdata' )
	sinkfile	<- paste0( 'Rout/compiled_gxewas/'	, env, '_pcxe.Rdata' )
	if( file.exists( compfile ) )	next

	allbeta	<- array( NA, dim=c(4,nsnp), dimnames=list( c( 'hom0', 'hom', 'hom1', 'het' )							, snpnames ) )
	allpv		<- array(	NA, dim=c(5,nsnp), dimnames=list( c( 'hom0', 'hom', 'het', 'het_chisq', 'glob' ), snpnames ) )
	for( batch in 1:nbatch ){ 

		savefile	<- paste0( 'Rdata/gxewas/batch_', env, '_', batch, '_pcxe.Rdata' )
		if( ! file.exists( savefile ) )	next

		cols	<- (batch-1)*nperbatch+1:nperbatch 

		load( savefile )
		if( batch == nbatch ){
			cols	<- intersect( cols, 1:nsnp )
			beta <- beta[,1:length(cols)]
			pv   <- pv  [,1:length(cols)]
		} 
		stopifnot( all.equal( snpnames[cols], colnames(beta) ) )

		allbeta	[-1,cols]<- beta
		allpv		[-1,cols]<- pv
		rm( beta, pv ) 

		load( paste0( 'Rdata/gwas/batch_', batch, '', '.Rdata' ) ) 
		if( batch == nbatch ){
			beta <- beta[,1:length(cols)]
			pv   <- pv  [,1:length(cols)]
		} 
		stopifnot( all.equal( snpnames[cols], colnames(beta) ) )

		allbeta	[1,cols]<- beta	['hom',]
		allpv		[1,cols]<- pv		['hom',]
		rm( beta, pv ) 
	} 
	save( allbeta, allpv, file=compfile )
}
