rm( list=ls() )
library(Matrix)
library(rgwas)
library(gplots)
source( 'pull_data.R' )
load( 'Rdata/setup.Rdata' )

nperm	<- 1e3 
for( perm.i in 0:nperm )
	for( p1 in sample(prs) )
		for( p2 in prs )
{

	savefile	<- paste0( 'Rdata/ce/', p1, '_', p2, '_', perm.i, '.Rdata' )
	sinkfile	<- paste0( 'Rout/ce/'	, p1, '_', p2, '_', perm.i, '.Rout' ) 
	if( file.exists( savefile ) | file.exists( sinkfile ) )	next
	print(sinkfile )
	sink(	sinkfile )

	y	<- getPheno('disease_ASTHMA_DIAGNOSED',binary=TRUE)
	X	<- getCovars(nPC=nPC)
	X	<- X[,c( "IID", "FID", "female", "age", "age2", paste0( 'PC_', 1:nPC ))]

	Z1	<- read.table( paste0( 'parsed_data/', p1, '_', 'odd.1.0.profile'  ), head=T ) [,c(1,2,6)]
	Z2	<- read.table( paste0( 'parsed_data/', p2, '_', 'even.1.0.profile' ), head=T ) [,c(1,2,6)]
	colnames(Z1)[3]	<- paste0( 'prs_', p1, '_odd'  )
	colnames(Z2)[3]	<- paste0( 'prs_', p2, '_even' )

	if( perm.i > 0 ){
		set.seed( perm.i + 123 )
		myperm	<- sample( nrow(Z1) ) 
		Z1[,3]	<- Z1[myperm,3]
		Z2[,3]	<- Z2[myperm,3] 
	}

	dat	<- merge( y		, Z1 )
	dat	<- merge( dat	, Z2 )
	dat	<- merge( dat	, X )
	dat	<- dat[complete.cases(dat),]

	y		<- dat[,3]
	Z1	<- scale( dat[,4] )
	Z2	<- scale( dat[,5] )
	X		<- scale( dat[,-(1:5)] )
	rm( dat ); gc()

	PCs	<- X[,paste0( 'PC_', 1:nPC )]
	X	<- cbind( X, matrix(Z1,nrow(Z1),nPC) * PCs, matrix(Z2,nrow(Z2),nPC) * PCs )
	X	<- scale(X)

	homlm	<- glm( y ~ X + Z1 + Z2 , family='binomial' )
	hetlm <- glm( y ~ X + Z1 * Z2	, family='binomial' )

	hetp		<- anova( hetlm , homlm , test='Chisq')$Pr[2]
	coefs	<- summary( hetlm )$coef

	save( hetp, coefs, file=savefile )

	sink()
	rm( homlm, hetlm, y, X, Z1, Z2 ); gc()
}
