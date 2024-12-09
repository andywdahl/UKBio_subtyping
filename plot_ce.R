rm( list=ls() )
library(gplots)
load( 'Rdata/setup.Rdata' )
nperm	<- 1e3
cols  <- redblue(20)

metap_fxn	<- function( b1, b2, se1, se2 ){
	w1	<- 1/(se1)^2
	w2	<- 1/(se2)^2
	z		<- (w1*b1 + w2*b2)/sqrt(w1+w2)
	2*min( pnorm(z), pnorm(-z) )
} 

metps	<- array( NA, dim=c(length(prs),length(prs),nperm+1), dimnames=list( prs, prs, 0:nperm ) )
hetps	<- array( NA, dim=c(length(prs),length(prs),nperm+1), dimnames=list( prs, prs, 0:nperm ) )
betas	<- array( NA, dim=c(length(prs),length(prs),nperm+1), dimnames=list( prs, prs, 0:nperm ) )
ses		<- array( NA, dim=c(length(prs),length(prs),nperm+1), dimnames=list( prs, prs, 0:nperm ) )
for( p1 in prs )
	for( p2 in prs )
		for( perm.i in 0:nperm )
{
	load( paste0( 'Rdata/ce/', p1, '_', p2, '_', perm.i, '.Rdata' ) )
	hetps[p1,p2,perm.i+1]	<- hetp
	betas[p1,p2,perm.i+1]	<- coefs['Z1:Z2',1]
	ses	 [p1,p2,perm.i+1]	<- coefs['Z1:Z2',2]
	if( p1 == p2 ){
	metps[p1,p2,perm.i+1]	<- hetps[p1,p1,perm.i+1]
	} else {
	metps[p1,p2,perm.i+1]	<- metap_fxn( betas[p1,p2,perm.i+1], betas[p1,p2,perm.i+1], ses[p1,p2,perm.i+1], ses[p1,p2,perm.i+1] ) 
	}
}

emp_ps		<- apply( hetps, 1:2, function(x) mean( x[1] >=x ) )
emp_ps_m	<- apply( metps, 1:2, function(x) mean( x[1] >=x ) )

print(hetps[,,1])
print(metps[,,1])
print(round( emp_ps  , 4 ))
print(round( emp_ps_m, 4 ))

pdf( paste0( 'figs/ce/CE_permuted_pvalues.pdf' ), width=16, height=16 )
par( mfcol=c(np,np) )
par( mar=c(1,1,5,1) )
for( i in 1:np )
	for( j in 1:np )
{ 
	if( i > j ){ plot.new(); next } 
	hist( hetps[i,j,], breaks=21, xlim=0:1, main=paste0( prs0[i], ':', prs0[j] ), xlab='', ylab='' )
	abline( v=hetps[i,j,1], col=2, lwd=5 )
}
dev.off() 


be <- betas[,,1]
pv <- hetps[,,1]
pv_meta <- metps[,,1]

pdf( 'figs/ce/gammas_heatmap.pdf', width=6.05, height=7 )
par( mar=c(9,9,4,1) )

breaks  <- seq( -1, 1, len=length(cols)+1 )*max(abs(be))
P				<- ncol(be)
image( x=1:P, y=1:P, t(be)[1:P,P:1], axes=F, breaks=breaks, col=cols, xlab='', ylab='', main=expression( hat(gamma)['EO'] ), cex.main=1.7 )

tex	<- apply( pv, 1:2, function(p) ifelse( p>.05, '', paste0( 'p=', format(p,digits=2,scientific=T) ) ) )
for( i in 1:P )
	for( j in 1:P )
		text( i, j, lab=(t(tex)[1:P,P:1])[i,j], cex=1 )
axis(1, at=1:P, lab=prs0, las=2 )
axis(2, at=P:1, lab=prs0, las=2 )
mtext( side=1, line=7.5, 'Odd Chrom PRS'	, cex=1.5 )
mtext( side=2, line=7.5, 'Even Chrom PRS'	, cex=1.5 )
dev.off()


pdf( 'figs/ce/gammas.pdf', width=6, height=6 )
par(mar=c(5.7,5.7,1,1))
lims	<- max(abs(be))*c(-1,1) + c(-.005,.005)
plot( lims, lims, type='n', xlab=expression( hat(gamma)['OE'] ), ylab=expression( hat(gamma)['EO'] ), axes=F, cex.lab=1.6, main='' )
axis(1);axis(2)
abline( h=0, col=2, lwd=.5 ); abline( v=0, col=2, lwd=.5 ); abline( a=0, b=1  , lwd=.5 )

ut	<- function(x) x[upper.tri(x)]
tut	<- function(x) t(x)[upper.tri(x)]

signifs_m				<- pv_meta < .05
points( ut(be)	,tut(be)	, pch=c(NA,'*')	[1+ut		(signifs_m)], cex=1.8 )
points( ut(be)	,tut(be)	, pch=c(16,NA)	[1+ut		(signifs_m)], cex=0.6 )
points( diag(be),diag(be)	, pch=c(NA,'*')	[1+diag	(signifs_m)], cex=1.8 )
points( diag(be),diag(be)	, pch=c(16,NA)	[1+diag	(signifs_m)], cex=0.6 ) 

tex	<- matrix( paste0( rep( prs0, each=np ), ':', rep( prs0, np ) ), np, np )
tex[which(signifs_m==0)]	<- '' 
text(	ut	(be)	,tut	(be-.003)	, ut	(tex)	, cex=.78 )
text(	diag(be)	,diag	(be-.003)	, diag(tex)	, cex=.78 )

tex	<- apply( pv_meta, 1:2, function(p) ifelse( p>.05, '', paste0( 'p=', format(p,digits=2,scientific=T) ) ) )
#tex	<- matrix( tex, np, np )
text(	ut	(be)	,tut	(be-.005)	, ut	(tex)	, cex=.68 )
text(	diag(be)	,diag	(be-.005)	, diag(tex)	, cex=.68 )

dev.off()



eigvals <- apply( betas, 3, function(be) eigen( (be+t(be))/2, symm=T )$val )
#eigvals <- apply( betas, 3, function(be) svd(be)$d )

pdf( 'figs/ce/gamma-eigen.pdf', width=8.2, height=6 ) 
layout( matrix(2:1,1,2), width=c(3.9,7) )

par( mar=c(5,5,2,1) ) 
plot( c(1,np), range( eigvals ), type='n', axes=F, xlab='Index', ylab='Eigenvalues of cross-trait gamma matrix', main='' )
axis( 1, at=1:np )
axis( 2 )
for( j in 1+1:nperm )
lines( 1:np, eigvals[,j], col='lightgrey', lwd=.3 ) 
abline( h=0, lty=1, lwd=2, col=2 ) 
lines( 1:np, apply( eigvals, 1, quantile, .025 ), col=1, lwd=1.2 )
lines( 1:np, apply( eigvals, 1, quantile, .975 ), col=1, lwd=1.2 )
lines( 1:np, eigvals[,1], col=3, lwd=3 )
points(1:np, eigvals[,1], col=3, pch=16, cex=1.2 )
legend( 'topright', col=c( 3, 'grey', 1 ), leg=c( 'Real data', 'Permuted PRS', '95% Empirical CI' ), bty='n', lty=rep(1,3), lwd=2 )

par( mar=c(3,6.5,4,2) ) 
u		<- eigen( (be+t(be))/2, symm=T )$vec[,1:2,drop=F] 
u		<- apply( u, 2, function(x) x*sign(x[1]) )
breaks  <- seq( -1, 1, len=length(cols)+1 )*max(abs(u))
image( x=1:2, y=1:np, t(u[np:1,,drop=F]), axes=F, xlab='', ylab='', main=expression( PC[1] ~~~~ PC[2] ), cex.main=1.7, breaks=breaks, col=cols ) 
axis(2, at=np:1, lab=prs0, las=2 )

dev.off() 
