rm( list=ls() )
library(BEDMatrix)
load( 'Rdata/setup.Rdata' )
load( 'Rdata/setup_snps.Rdata' )

batch			<- as.numeric( commandArgs(TRUE)[[1]] )
if( batch > nbatch ) stop('batch too big' ) 

for( mode in c( '', '_juvenile' ) ){ # , '_noage', '_lm', '', '_flip' 
if( mode == '_juvenile' & batch >= -3 ) next

savefile	<- paste0( 'Rdata/gwas/batch_', batch, mode, '.Rdata' )
sinkfile	<- paste0( 'Rout/gwas/batch_'	, batch, mode, '.Rout' )
if( file.exists( savefile ) | file.exists( sinkfile ) )	next
print( sinkfile )
sink(	sinkfile )

env			<- 'disease_ALLERGY_ECZEMA_DIAGNOSED'
datfile	<- paste0( 'Rdata/munged/', env, mode, '.Rdata' )
if( !file.exists( datfile ) ) stop( datfile )
load( datfile )

### prepare SNPs
if( batch > 0 ){
	bedfile	<- paste0( 'split_ukb_genos/block.', batch )
	if( !file.exists( paste0( bedfile, '.bed' ) ) ){
		snpfile	<- paste0( 'split_ukb_genos/block.', batch, '.mysnps.txt' )
		mysnps	<- snpnames[intersect( (batch-1)*nperbatch+1:nperbatch, 1:nsnp )]
		write.table( mysnps, file=snpfile, quote=F, row.names=F, col.names=F )
		system( paste0( 'plink --bfile /u/project/zaitlenlab/nadavrap/UKBB/GEN/ukbb --extract ', snpfile, ' --keep /u/project/sriram/nadavrap/UKBB/CSV/final_covars.txt --make-bed --out ', bedfile, ' --threads 1 --memory 4000' ) )
		rm( mysnps, snpfile )
	}
} else { 
	bedfile	<- target_beds[ -1*batch ]
} 

### load SNPs
Gsnp	<- as.matrix( BEDMatrix( bedfile ) )
if( batch < 0 ){
	ncols	<- ncol(Gsnp)
	colnames_loc	<- colnames(Gsnp)
} else {
	ncols	<- nperbatch
	colnames_loc	<- snpnames[(batch-1)*ncols+1:ncols]
} 

####
iids	<- sapply( strsplit( rownames( Gsnp ), '_' ), function(x) x[[1]] )
if( !all( IID %in% iids ) ){
	names(y)		<- IID
	rownames(X)	<- IID 
	IID	<- intersect( IID, iids )
	y		<- y[IID]
	X		<- X[IID,]
}
rownames(Gsnp)	<- iids 

beta	<- array( NA, dim=c(3,ncols), dimnames=list( c( 'hom', 'hom1', 'het' ), colnames_loc ) )
times	<- array(	NA, dim=c(3,ncols), dimnames=list( c( 'hom', 'het', 'null' ), colnames_loc ) )
pv		<- array(	NA, dim=c(4,ncols), dimnames=list( c( 'hom', 'het', 'het_chisq', 'glob' ), colnames_loc ) )
for( jj in 1:ncol(Gsnp) ){
	print( jj )
	g		<- as.numeric( Gsnp[IID,jj] )
	sub	<- which( !is.na(g) )

	if( mode == '_lm' ){
	times['hom',jj]	<- system.time({ fit_hom	<- lm( y ~ X + g, subset=sub) })[3]
	} else {
	times['hom',jj]	<- system.time({ fit_hom	<- glm( y ~ X + g, subset=sub, family='binomial') })[3]
	}

	beta['hom'	,jj]	<- summary(fit_hom)$coef ['g',1]
	pv	['hom'	,jj]	<- summary(fit_hom)$coef ['g',4]

	rm( g, sub, fit_hom ); gc()
} 
print( gc() )
print( sum( times )/60 )

save( beta, pv, times, file=savefile )
print(warnings()); print('Done')
rm( y, beta, pv, X, Gsnp ); gc()

sink()
}
