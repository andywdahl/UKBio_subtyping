library(qqman)
my_qq <- function(y,cand_rs,
	xlab='Expected -log10(p)', ylab='Observed -log10(p)', main='', 
	maftol=.01,misstol=.10,
	adjust_gc=FALSE,
	highlights=NULL,
	mafs=rep(1,length(y)),miss=rep(0,length(y)), ... 
){

	y[mafs < maftol ] <- NA
	y[miss > misstol] <- NA

	y <- sort( y[!is.na(y)] )

	lam_gc<- lam_gc_fxn( y )
	if( adjust_gc )
	y     <- correct_gc( y  , lam_gc )

	lims  <- range( c( y , 8.1 ) )
	plot( lims, lims, type='n', xlab=xlab, ylab=ylab, main=main, ... )
	abline( a=0, b=1 )

	add_y_qq( y, highlights[[1]] )

	snps_gw	<- names(y)[ y > -log10( 5e-8 ) ]
	topsnps	<- list( snps_gw )

	if( !missing( cand_rs ) )
		for( i in 1:length(cand_rs) )
	try({
		y0  <- y[cand_rs[[i]]]
		y0  <- y0[!is.na(y0)]
		y1	<- y0 
		snps_cand	<- add_y_qq( y1, highlights[[i+1]] ) 
		if( length(snps_cand) > 0 ){
			snps_cand0	<- names(y0)[ y0 >= min(y1[snps_cand]) ]
			stopifnot( all( snps_cand %in% snps_cand0 ) )
		} else {
			snps_cand0	<- NULL
		} 
		topsnps	<- c( topsnps, list( snps_cand0 ) )
	})
	legend( 'bottomright', bty='n', leg=bquote( lambda['GC']==.(round( lam_gc, 3 )) ), cex=1.6 )

	list( topsnps=topsnps, lam_gc=lam_gc )
}

lam_gc_fxn  <- function( log10p ){
	chi2s   <- qchisq(1-10^-log10p,df=1)
	median(chi2s,na.rm=T)/qchisq(0.5,1)
}
correct_gc  <- function( log10p, lam_gc ){
	lam_gc  <- max( lam_gc, 1 )
	chi2s   <- qchisq(-log10p*log(10),df=1, lower.tail=F, log=T)
	-pchisq( chi2s/lam_gc, df=1, lower.tail=F, log=T )/log(10)
}

add_y_qq  <- function(y,highlights=NULL,col){
	y <- sort( y[!is.na(y)] )
	n <- length(y)
	x <- sort( -log10( 1:n/(n+1) ) )

	if( length(y) > 1e4 ){
		signif	<- 10^-y < 5e-8
		highlights <- intersect( highlights, names(y)[signif] )
		cexs	<- rep(.8,length(y))
		cols	<- rep( 'grey',length(y))
		if( any(signif) ){
			cols[signif]	<- 2
			cexs[signif]	<- 2.2
		}
		xadj	<- 0
		yadj	<- .5
		srt		<- 0
	} else {
		bhs		<- p.adjust( 10^-y, 'BH' )
		signif	<- bhs < .1
		bonfs   <- p.adjust( 10^-y, 'bonferroni' ) < .05 
		highlights <- intersect( highlights, names(y)[signif] )

		print( 'Bonferronis:' )
		print( names(y)[bonfs] )

		print( 'FDR <.1:' )
		print( names(y)[setdiff(bhs,bonfs)] )

		cols  <- c( 3, 4, 5 )[1 + ( bhs < .1 ) + bonfs]
		cexs  <- .8 + 1.4*( bhs < .1 ) 

		xadj	<- seq( -.7, .7, len=length(highlights) )
		yadj	<- c(.2,2.0,1.9)[1:length(highlights)]
		srt		<- 45
	}
	if(!missing(col))	cols	<- col
	points( x, y, pch=16, col=cols, cex=cexs )

	if( length( highlights ) > 0 ){
		names(x)	<- names(y) 
		text(		x[highlights]+xadj, y[highlights] + yadj, highlights, cex=1.2, srt=srt )
		arrows( x[highlights]			, y[highlights], 
						x[highlights]+xadj, y[highlights] + yadj-.3, code=3, angle=90, length=0, lwd=.9 )
	}
} 

xyplot	<- function(x,y,...){ 
	plot(		x, y, type='n', xlim=range(c(x,y),na.rm=T), ylim=range(c(x,y),na.rm=T), ... )
	points( x, y, pch=16 )
	abline( a=0, b=1, col=2, lwd=2 ) 
}
