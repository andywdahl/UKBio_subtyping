ldclump <- function( rsids, pvs, r2=.1, maf=.01, bedfile='/u/project/zaitlenlab/nadavrap/UKBB/GEN/ukbb', plink='plink' ){ 
	if( length(rsids) %in% 0:1 ) return(rsids) 
	cat( '# unclumped rsids:', length( rsids ), '\n' ) 
	mat	<- cbind(rsids,pvs)
	colnames(mat)	<- c( 'SNP', 'P' )
	write.table( mat	, file='Rout/tmp.clump.input'	, quote=F, row.names=F, col.names=TRUE )
	write.table( rsids, file='Rout/tmp.extract'			, quote=F, row.names=F, col.names=FALSE ) 
	system( paste0( plink, ' --threads 1 --bfile ', bedfile, ' --extract Rout/tmp.extract --clump Rout/tmp.clump.input --clump-p1 .99 --clump-p2 .99 --clump-r2 ', r2, ' --clump-kb 1000000 --out Rout/tmp --maf ', maf, ' --chr 1-22 --geno .1' ) ) 
	out	<- as.character( read.table( 'Rout/tmp.clumped', head=T )[,3] )
	system( 'rm Rout/tmp.{clump.input,clumped,extract,log,nosex}' )
	cat( '# clumped rsids:', length(out), '\n' ) 
	out 
}

ldprune <- function( rsids, r2=.1, bedfile='/u/project/zaitlenlab/nadavrap/UKBB/GEN/ukbb', plink='plink' ){
	if( length(rsids) %in% 0:1 ) return(rsids)
	write.table( rsids, file='Rout/tmp.extract', quote=F, row.names=F, col.names=F )
	system( paste0( plink, ' --threads 1 --bfile ', bedfile, ' --extract Rout/tmp.extract --indep-pairwise 10000 5 ', r2, ' --out Rout/tmp --maf .01 --chr 1-22 --geno .1' ) )
	as.character( read.table( 'Rout/tmp.prune.in' )[,1] )
}
