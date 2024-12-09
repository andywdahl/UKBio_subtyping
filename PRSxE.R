rm( list=ls() )
library(Matrix)
library(rgwas)
library(gplots)
load( 'Rdata/setup.Rdata' )

for( p1 in prs )
	for( env in envs )
{ 
	savefile	<- paste0( 'Rdata/PRSxE/'	, p1, '_', env, '.Rdata' )
	sinkfile	<- paste0( 'Rout/PRSxE/'	, p1, '_', env, '.Rout' ) 
	if( file.exists( savefile ) | file.exists( sinkfile ) )	next
	print(sinkfile )
	sink(	sinkfile )

	load( paste0( 'Rdata/munged/', env, '_pcxe.Rdata' ) ) 

	PS	<- read.table( paste0( 'parsed_data/', p1, '.1.0.profile' ), head=T ) [,c(1,2,6)]  # only use 1.0 bc in-sample
	IID <- data.frame( as.numeric(IID) )
	colnames(IID) <- 'IID'
	PS	<- scale( merge( IID, PS )[,3] )

	hetp		<- anova(
		glm( y ~ X + Z + PS				, family='binomial' ),
		glm( y ~ X + Z + PS + Z:PS, family='binomial' ),
		test='Chisq')$Pr[2]

	negZ	<- 1-Z
	coefs	<- summary( glm( y ~ X + Z + PS:Z + PS:negZ	, family='binomial' ) )$coef
	homcoefs	<- summary( glm( y ~ X + PS							, family='binomial' ) )$coef

	save( hetp, homcoefs, coefs, file=savefile )
	sink()
	rm( homlm, hetlm, y, X, PS, Z ); gc()
}


ps	  <- array( NA, dim=c(length(prs),1+  length(envs)), dimnames=list( prs, c( 'hom', envs ) ) )
betas	<- array( NA, dim=c(length(prs),1+2*length(envs)), dimnames=list( prs, c( 'hom', paste0( rep( envs, each=2 ), rep(c('-','+'),3) ) ) ) )
ses		<- array( NA, dim=c(length(prs),1+2*length(envs)), dimnames=list( prs, c( 'hom', paste0( rep( envs, each=2 ), rep(c('-','+'),3) ) ) ) )
for( p1 in prs )
	for( env in envs ) 
{
	load( paste0( 'Rdata/PRSxE/'	, p1, '_', env, '.Rdata' ) )

	betas [p1,'hom']	<- homcoefs['PS',1]
	ses   [p1,'hom']	<- homcoefs['PS',2]
	ps    [p1,'hom']	<- homcoefs['PS',4]

	betas [p1,paste0(env,c('-','+'))]	<- coefs[c('PS:negZ','Z:PS'),1]
	ses   [p1,paste0(env,c('-','+'))]	<- coefs[c('PS:negZ','Z:PS'),2]
	ps    [p1,env]	<- hetp

	rm( coefs, homcoefs, hetp )
}


xlabs	<- c( 'Additive', 'Atopy-', 'Atopy+', 'Low-Eos', 'High-Eos', 'Low-Neurot', 'High-Neurot' ) 
cols	<- c('grey',rep(c('pink','purple'),3)) 
at0s	<- c( .2, 1:3*.5 + .1 ) 
ats		<- c( .2, rep(1:3,each=2)*.5 + rep(0:1)*.2 )
xlim	<- range(ats,na.rm=T)+c(-2,1)*.05 

pdf( 'figs/PSxE.pdf', width=5, height=7  ) 
layout( cbind( 1, 1+matrix(1:5,5,1) ), width=c(1,16), height=c( rep(1,4), 1.1 ) )

par( mar=rep(0,4) )
plot.new()
mtext( side=2, srt=90, 'Odds Ratio for Asthma Risk', cex=1.5, line=-2.8, adj=.7 )


par( mar=c(.5,6.5,.5,.5) )
for( p1 in prs ){

	y	<- betas[p1,]
	se<- ses[p1,]

	ylim		<- range( c( y-2*se, y+2*se ), na.rm=T ) 
	ylim[2]	<- ylim[2]+.4*(ylim[2]-ylim[1])
	if( p1 != 'disease_ASTHMA_DIAGNOSED' )
		ylim		<- range( c( ylim, 0 ) )

	plot(  xlim, ylim, type='n', axes=F, ylab=paste0( prs0[p1], ' PRS' ), xlab='', cex.lab=1.3 )
	box()
	axis(2,cex.axis=.9)
	abline( h=0, col=2, lwd=1.1, lty=1 )

	points(	ats, y										, col=cols, cex=1.9, pch=16) 
	arrows(	ats, y-2*se, ats, y+2*se	, col=cols, lwd=1.5, length=0.08, angle=90, code=3)

	pvseq	<- ps[p1,]
	labs2	<- round( pvseq, 2 )
	small <- which( pvseq < .01 )
	if( length(small) > 0 )
		labs2[small]<- format( pvseq[small], digits=1 )
	labs2<- as.expression(sapply( labs2, function(pv) bquote( p==.(pv)  ) ) )
	ytxt	<- ylim[2]*.94 + ylim[1]*.06
	text(		at0s,rep(ytxt,3)		, col=1, cex=1.2, lab=labs2 ) 

} 
axis( 1, at=ats	, lab=xlabs			, las=2, cex.axis=1.3 ) 
mtext( side=1, srt=90, 'Observed Risk Factors', cex=1.5, line=9.5 )

dev.off()
