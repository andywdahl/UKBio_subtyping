rm( list=ls() )
#source( 'pull_data.R' )
load( 'Rdata/setup.Rdata' )
load( 'Rdata/setup_snps.Rdata' )
load( 'Rdata/cand_pvals_compiled.Rdata' )

envs	<- c( 'disease_ALLERGY_ECZEMA_DIAGNOSED', 'eo_qn', 'neuroticism' ) 
batches <- -(11:17) 
outses	<-
outbets	<- array( NA, dim=c(2,3,length(keysnps),2), dimnames=list( c('snp','locus'), envs, keysnps, c( 'hom', 'e=1' ) ) )
for( batch in batches )
	for( env in envs )
{

	load( paste0( 'Rdata/gwas/batch_', batch, '', '.Rdata' ) )
	pv0		<- pv
	beta0	<- beta
	rm(pv,beta)

	load( paste0( 'Rdata/gxewas/batch_', env, '_', batch, '_pcxe.Rdata' ) ) 
	colnames(beta0)	<-
	colnames(pv0)		<-
	colnames(beta)	<-
	colnames(se)		<-
	colnames(pv)		<- sapply( strsplit( colnames( pv ), '_' ), function(x) x[[1]] ) 

	recover_se	<- function( beta, pv ){
		z	<- qnorm( pv / 2 )
		beta / z
	}

	snp <- keysnps[-batch]
	top <- colnames(pv)[which.min( pv[2,] )]
	outbets['snp'  ,env,snp,]	<- c( beta0['hom',snp]															, beta['het',snp]) * sign(beta0['hom',snp]) 
	outbets['locus',env,snp,]	<- c( beta0['hom',top]															, beta['het',top]) * sign(beta0['hom',top])
	outses ['snp'  ,env,snp,]	<- c( recover_se( beta0['hom',snp], pv0['hom',snp] ), se  ['het',snp])
	outses ['locus',env,snp,]	<- c( recover_se( beta0['hom',top], pv0['hom',top] ), se  ['het',top]) 

	rm( beta0, pv0, beta, se, pv ) 
}

batches	<- -batches

ys		<- t(sapply( batches, function(i) i+c(-.2,-.4) ))
cols	<- rep( 1:2, each=length(batches) )
outbets	<- outbets
outses	<- outses 

for( testmode in c( 'locus', 'snp' ) ){

pdf( paste0( 'figs/galaIIsnps_forest_', testmode, '.pdf' ), width=5, height=6.5 ) 
layout( matrix( c(1:4,rep(5,4)), 2, 4, byrow=T ), width=c(.9,rep(6,3)), height=c(10,.5) )
par( mar=rep(0,4) )
plot.new()

par( mar=c(4.5,0.6,3.0,0.6) ) 
for( env in envs ){
	xlim	<- c(-.3,.5)
	plot(	xlim, range(ys), type='n', main=paste0( 'Gx', niceEnvs[env] ), axes=F, xlab='Log OR on Asthma', ylab='' )
	box()
	axis(1)
	if( env == envs[1] )
		axis(2, at=batches, lab=keysnps[batches], cex.axis=0.95 ) 

	abline( v=0, col='grey' )
	points(		outbets[testmode,env,batches,], ys, pch=16, col=cols, cex=1.4 )
	arrows(
		(outbets - outses*1.96)[testmode,env,batches,], ys,
		(outbets + outses*1.96)[testmode,env,batches,], ys, 
		length=0.05, angle=90, code=3, col=cols )

	nstars	<- sapply( outps[testmode,env,keysnps[batches]],function(pv) sum( pv < c( .05, .005 ) ) )
	points(	rep(0,length(batches)), ys[,1]+.20, pch=c(NA,8,8,15)[ nstars + 1 ], cex=c( 1, 1.5, 2.5, 4 )[ nstars + 1 ], col=c( 1, 3, 4, 5 )[ nstars + 1 ] ) 
}

par( mar=rep(0,4) )
plot.new()
legend( 'center', bty='n', horiz=T, leg=c( 'Additive', 'Interaction', 'p < .05', 'p < .005' ), col=c(1,2,3,4), cex=1.3, pch=c(16,16,8,8), pt.cex=c(1.5,1.5,1,2) )

dev.off() 
} 

pvec <- outps['snp',,keysnps[batches]]
cat( 'binomial test:', binom.test( sum( pvec < .05 ), length( pvec ), alternative="greater", p=.05 )$p.value, '\n')
