rm( list=ls() )
load( 'Rdata/setup.Rdata' )

sig2names	<- c('sig2g', 'v_1', 'v_2')
nsig2			<- length(sig2names)

foldsizes	<- c( 1e4, 3e4, 4e4 ) 
max_nfolds<- 100

sigmas		<- array( NA, dim=c( length(envs), nsig2, length(foldsizes), max_nfolds ), dimnames=list( envs, sig2names, foldsizes, 1:max_nfolds ) )
for( foldsize in foldsizes[3] )
	for( env in envs )
		for( fold.i in 1:floor(maxN/foldsize) )
{
	savefile	<- paste0( 'Rdata/gxemm/'	, env, '_', fold.i, '_', foldsize, '_TRUE_1.Rdata')
	if( ! file.exists(savefile) ){ print(savefile); next }
	load( savefile )

	sigmas[env,'sig2g',which(foldsizes==foldsize),fold.i]<- out$fits$hom$h2[1]
	sigmas[env,1+1:2  ,which(foldsizes==foldsize),fold.i]<- (out$fits$free$sig2s * out$fits$free$ws )[1+1:2] 
	rm( out, kintime )
}

for( e in envs )
print( apply( !is.na(sigmas[e,1:2,3,]), 1, sum, na.rm=T ) )

print( apply( sigmas[,,3,], 1:2, mean, na.rm=T ) )
print( apply( sigmas[,,3,], 1:2, function(x) sd(x,na.rm=T)/sqrt(sum(!is.na(x))) ) )

waldtest	<- function(ys){
	y		<- mean(ys,na.rm=T )
	se	<- sd(	ys,na.rm=T ) / sqrt(sum( !is.na(ys) ))
	pnorm( y / se, lower.tail=F )
}

for( foldsize in c( 3e4, 4e4 ) ){

	pdf( paste0( 'figs/gxemm_', foldsize, '.pdf' ), width=5.7, height=6.2 )
	par( mar=c(7,5,.5,.5) )

	ys	<- sigmas[1,'sig2g',which(foldsizes==foldsize),]
	pvseq	<- waldtest( ys )
	for( env in envs ){
		load_foldsize	<- foldsize
		ys_e	<- sigmas[env,c('v_1','v_2'),which(foldsizes==load_foldsize),]
		pvs		<- apply( ys_e, 1, function(x) waldtest(x) )
		ys		<- rbind( ys, ys_e )
		pvseq	<- c( pvseq, pvs )
	}

	labs2			<- rep(NA,length(pvseq))
	smalls		<- which( pvseq < .01 )
	labs2[-smalls]<- format( pvseq[-smalls], nsmall=2, digits=1 )
	labs2[ smalls]<- format( pvseq[ smalls], digits=1 )
	labs2<- as.expression(sapply( labs2, function(pv) bquote( p==.(pv)  ) ) )

	ny	<- nrow(ys)
	ats	<- 1+c( .3, rep(1:3,each=2)*.6 + rep(0:1)*.25 )
	xlim<- range(ats,na.rm=T)+c(-1,1)*.05

	xlabs	<- expression(
		h[g]^2,
		h['Atopy-']^2, h['Atopy+']^2,
		h['Low-Eos']^2, h['High-Eos']^2,
		h['Low-Neurot']^2, h['High-Neurot']^2
	)

	ylim	<- c(-.4,.7)

	plot(  xlim, ylim, type='n', axes=F, ylab='Asthma Variance Explained (%)' , xlab='', cex.lab=1.3 )
	text( ats[1]-.05, .65, cex=1.3, 'e' )
	legend( 'topright', bty='n', cex=0.9, leg=c( 'Meta-analysis', 'Each subsample' ), pch=c(16,16), pt.cex=c(1.4,.8) )

	axis(2,at=0:3/5,lab=0:3/5*100)
	abline( h=0, col=2, lwd=1.1, lty=1 )

	nfold	<- length(ys)/length(ats)
	axis( 1, line=0, at=ats			, lab=rep('',ny), las=2, cex.axis=1.5 )
	axis( 1, line=0, at=ats-.06	, lab=xlabs			, las=2, cex.axis=1.4, tick=F )
	axis( 1, line=0, at=ats+.06	, lab=labs2			, las=2, cex.axis=0.8, tick=F ) 

	cols	<- c('grey',rep(c('pink','purple'),3)) 
	xs		<- rep(ats,nfold) + runif(nfold*ny,.04,.09)+.01
	y			<- apply(ys,1,mean,na.rm=T)
	se		<- apply(ys,1,sd,na.rm=T)/sqrt(apply(!is.na(ys),1,sum))

	points(xs	, c(ys)								, col=rep(cols,nfold)	, cex=0.8, pch=16 )
	points(ats, y										, col=cols						, cex=2.4, pch=16 )
	arrows(ats, y-2*se, ats, y+2*se	, col=cols						, lwd=2.7, length=0.08, angle=90, code=3) 

	dev.off()
}
