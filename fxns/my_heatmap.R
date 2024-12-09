my_heatmap	<- function( be, pv, #pdfname, width=12.1, height=7, 
	xtext='', ytext='', 
	xline=7.5, yline=7.5,ncolor=7,
	symm, xlabs=colnames(be), ylabs=rownames(be), main=expression( hat(gamma)['EO'] ),
	signifs
){ 

	cols    <- redblue(ncolor)
	breaks  <- seq( -1, 1, len=length(cols)+1 )*max(abs(be),na.rm=T)
	Px			<- ncol(be)
	Py			<- nrow(be)

	if( missing( signifs ) )
	if( symm ){
		P		<- Px
		signifs				<- ( pv      < .05 ) + ( pv      < .05/(P*(P-1)/2) ) 
		diag(signifs)	<- ( diag(pv)< .05 ) + ( diag(pv)< .05/P )
	} else {
		signifs				<- ( pv      < .05 ) + ( pv      < .05/(Px*Py) )
	}

	image( x=1:Px, y=1:Py, t(be)[1:Px,Py:1], axes=F, breaks=breaks, col=cols, xlab='', ylab='', main=main, cex.main=1.7 )
	for( i in 1:Px )
		for( j in 1:Py )
			points( i, j, pch='*', cex=1.5*(t(signifs)[1:Px,Py:1])[i,j] )
	axis(1, at=1:Px, lab=xlabs, las=2 )
	axis(2, at=Py:1, lab=ylabs, las=2 )
	mtext( side=1, line=xline, xtext, cex=1.5 )
	mtext( side=2, line=yline, ytext, cex=1.5 )

	list(
		cols    = cols  ,
		breaks  = breaks 
	) 
}
