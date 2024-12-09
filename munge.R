rm( list=ls() )
source( 'pull_data.R' )
load( 'Rdata/setup.Rdata' )

for( env in envs ){ 
	sinkfile	<- paste0( 'Rout/munge/'  , env, '_pcxe.Rout' )
	datfile		<- paste0( 'Rdata/munged/', env, '_pcxe.Rdata' )
	if( file.exists( datfile ) | file.exists( sinkfile ) )	next
	print( sinkfile )
	sink(	sinkfile )

	y	<- getPheno('disease_ASTHMA_DIAGNOSED',binary=TRUE)
	X	<- getCovars(nPC=nPC)
	Z	<- getPheno(loadnames[env], binary=( env %in% discrEnv ), qn=( env == 'eo_qn' ), scale01=TRUE )
	dat	<- merge( y, Z )
	dat	<- merge( dat, X ) 
	dat	<- dat[complete.cases(dat),] 
	nZ	<- ncol(Z)-2
	rm( y, X, Z ) 

	IID	<- as.character(dat[,1])
	y		<- as.numeric(	dat[,3] )
	Z		<- as.matrix(		dat[,3+1:nZ,drop=F] )
	X		<- as.matrix(		dat[,-(1:(3+nZ))] )
	X		<- scale(X)
	rm( dat ); gc()

	X	<- X[,c( "female", "age", "age2", paste0( 'PC_', 1:nPC ))]
	X	<- cbind(X, (Z[,1] %o% rep(1,nPC)) * X[,paste0( 'PC_', 1:nPC )] )
	X	<- scale(X)

	print( head( y )  )
	print( head( X )  )
	print( head( Z )  )

	save( y, X, Z, IID, file=datfile )
	rm( y, X, Z, IID ) 
	print(warnings()); print('Done')
	sink()
}
