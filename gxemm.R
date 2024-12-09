rm( list=ls() )
library(GxEMM)
load( 'Rdata/setup.Rdata' ) 

foldsize	<- as.numeric( commandArgs(TRUE)[[1]] ) 
for( env in sample(envs) )
	for( fold.i in sample( floor(maxN/foldsize) ) )
{
	savefile	<- paste0( 'Rdata/gxemm/'	, env, '_', fold.i, '_', foldsize, '_TRUE_1.Rdata')
	sinkfile	<- paste0( 'Rout/gxemm/'	, env, '_', fold.i, '_', foldsize, '_TRUE_1.Rout' ) 
	if( file.exists( savefile ) | file.exists( sinkfile ) )	next
	print( sinkfile )
	sink(	sinkfile )

	datfile	<- paste0( 'Rdata/munged/', env, '_pcxe.Rdata' )
	if( !file.exists( datfile ) ) stop( datfile )
	load( datfile ) 

	N			<- length(y)
	nfold	<- floor( N / foldsize )
	if( fold.i > nfold ){ print( c( fold.i, nfold ) ); sink(); next }

	set.seed(12346)
	sub		<- split( 1:N, sample( 1:nfold, N, rep=T ) )[[fold.i]]

	IID	<- IID[sub]
	y		<- y[sub] + 1
	Z		<- Z[sub,,drop=F]
	X		<- X[sub,] 
	X		<- scale(cbind( X, Z ))

	if( ncol(Z) > 1 ) stop(paste0( 'ncol(Z)=', ncol(Z) ))
	Z	<- cbind( 1-Z, Z )

	grmstem		<- paste0( 'kins/grm.'	, env, '_', fold.i, '_', foldsize, '_1' )
	keepfile	<- paste0( 'kins/keep.'	, env, '_', fold.i, '_', foldsize, '_1.txt' )
	write.table( cbind( IID, IID ), file=keepfile, col.names=F, row.names=F, quote=F )
	system( paste0( '/u/home/a/andyd/ldak5.linux --bfile /u/project/zaitlenlab/nadavrap/UKBB/GEN/ukbb --calc-kins-direct ', grmstem, ' --keep ', keepfile, ' --ignore-weights YES --power -1 --kinship-raw YES' ) )
	GRM	<- as.matrix( read.table( paste0( grmstem, '.grm.raw' ) ) )

	out	<- GxEMM:::full_GxEMM( y=y, X=X, K=GRM, Z=Z, binary=TRUE, ldak_loc='/u/home/a/andyd/ldak5.linux', tmpdir=paste0( grmstem, '_gxemm_tmp' ), hom_IID=FALSE, hom_free=FALSE, IID=FALSE )
	print( out )
	system( paste0( 'rm ', grmstem, '.{grm.adjust,grm.bin,grm.details,grm.id,grm.raw,progress}' ) )
	file.remove( keepfile )
	
	save( out, file=savefile )
	rm(kintime, out,y,X,Z,GRM); gc()
	warnings(); sink()
}
