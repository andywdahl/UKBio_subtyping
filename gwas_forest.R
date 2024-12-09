rm( list=ls() )
library(BEDMatrix)
source( 'fxns/R2plink.R' )
load( 'Rdata/setup.Rdata' )
load( 'Rdata/setup_snps.Rdata' )

envs	<- c( 'disease_ALLERGY_ECZEMA_DIAGNOSED', 'eo_qn', 'neuroticism' )
snps	<- c( 'rs532965', 'rs1837253', 'rs12964116', 'rs146125856', 'rs34645399', 'rs72823646', 'rs12123821', 'rs1701704', 'rs6932730', 'rs241437', 'rs17500510', 'rs9268515' ) 

bedfile	<- 'split_ukb_genos/hits'
if( !file.exists( paste0( bedfile, '.bed' ) ) ){
	snpfile	<- 'split_ukb_genos/hits.mysnps.txt'
	write.table( snps, file=snpfile, quote=F, row.names=F, col.names=F )
	system( paste0( 'plink --bfile /u/project/zaitlenlab/nadavrap/UKBB/GEN/ukbb --extract ', snpfile, ' --keep /u/project/sriram/nadavrap/UKBB/CSV/final_covars.txt --make-bed --out ', bedfile, ' --threads 1 --memory 4000' ) )
}

for( env in envs ){ 
	outfile <- paste0( 'Rdata/hits_', env, '_pcxe.Rdata' )
	if( ! file.exists( outfile ) ){
		sink( paste0( 'Rout/hits_', env, '_pcxe.Rout' ) )
		load( paste0( 'Rdata/munged/', env, '_pcxe.Rdata' ) )
		Gsnp	<- as.matrix( BEDMatrix( bedfile ) )
		iids	<- sapply( strsplit( rownames( Gsnp ), '_' ), function(x) x[[1]] )
		stopifnot( all( IID %in% iids ) )
		rownames(Gsnp)	<- iids 
		colnames(Gsnp)	<- sapply( strsplit( colnames(Gsnp), '_' ), function(x) x[[1]] ) 
		stopifnot( all( IID %in% iids ) )

		Z0 <- 1-Z
		Z1 <- Z

		beta	<- array( NA, dim=c(3,length(snps)), dimnames=list( c( 'hom', 'g:Z0', 'g:Z1' ), snps ) )
		se		<- array( NA, dim=c(3,length(snps)), dimnames=list( c( 'hom', 'g:Z0', 'g:Z1' ), snps ) )
		times	<- array(	NA, dim=c(3,length(snps)), dimnames=list( c( 'hom', 'het', 'null' ) , snps ) )
		pv		<- array(	NA, dim=c(4,length(snps)), dimnames=list( c( 'hom', 'het', 'het_chisq', 'glob' ), snps ) )
		for( snp in snps ){
			g		<- as.numeric( Gsnp[IID,snp] )
			sub	<- which( !is.na(g) )

			times['hom',snp]	<- system.time({ fit_hom	<- glm( y ~ X + Z + g				   , subset=sub, family='binomial') })[3]
			times['het',snp]	<- system.time({ fit_het	<- glm( y ~ X + Z + g:Z0 + g:Z1, subset=sub, family='binomial') })[3]

			beta['hom'	,snp]	<- summary(fit_hom)$coef['g',1]
			se  ['hom'	,snp]	<- summary(fit_hom)$coef['g',2]
			pv	['hom'	,snp]	<- summary(fit_hom)$coef['g',4] 
			beta['g:Z0'	,snp]	<- summary(fit_het)$coef['g:Z0',1]
			se  ['g:Z0'	,snp]	<- summary(fit_het)$coef['g:Z0',2]
			beta['g:Z1'	,snp]	<- summary(fit_het)$coef['g:Z1',1]
			se  ['g:Z1'	,snp]	<- summary(fit_het)$coef['g:Z1',2]

			pv	['het_chisq'  ,snp]	<- anova(
				fit_hom,
				glm( y ~ X + Z + g + g:Z1, subset=sub, family='binomial'),
			test='Chisq')$Pr[2]
			rm( fit_hom ); gc()

			times['null',snp]	<- system.time({ 
				pv	['glob'	,snp]	<- anova( glm( y ~ X + Z, subset=sub, family='binomial'	), fit_het, test='Chisq')$Pr[2]
			})[3] 

			rm( g, sub, fit_het ); gc()
		}
		save( beta, se, pv, times, file=outfile )
		print(warnings()); print('Done')
		sink()
	}
}

dat		<- array( NA, dim=c(4,3,3,length(snps)), dimnames=list( c('gwas','g:Z0','g:Z1','glob'), c('beta','se','pv'), envs, snps) ) 
for( env in envs ){
	load( paste0( 'Rdata/hits_', env, '_pcxe.Rdata' ) )
	dat[-4,1,env,]	<- beta[c('hom','g:Z0','g:Z1'),snps]
	dat[-4,2,env,]	<- se  [c('hom','g:Z0','g:Z1'),snps]
	dat[-3,3,env,]	<- pv	 [c('hom','het_chisq','glob') ,snps]
	rm( beta, se, pv ) 
}

xlabs	<- c(
	'Additive',
	'Atopy-', 'Atopy+',
	'Low-Eos', 'High-Eos',
	'Low-Neurot', 'High-Neurot'
)
ats		<- 1+c( -.05, rep(1:3,each=2)*.6 + rep(0:1,3)*.20 )
at0s	<- 1+c(  .05, 1:3*.6 + .5*.25 )
xlim	<- range(ats,na.rm=T)+c(-2,1)*.05

cols	<- c('grey',rep(c('pink','purple'),3))

for( snp in snps ){

	pdf( paste0( 'figs/forest/snp_', snp, '.pdf' ), width=5.3, height=7 )
	par( mar=c(7,5,3,.5) )

	y		<- c( dat[1,1,1,snp], rbind( dat[2,1,,snp], dat[3,1,,snp] ) )
	y.se<- c( dat[1,2,1,snp], rbind( dat[2,2,,snp], dat[3,2,,snp] ) )
	y		<- y * sign(y[1])
	pv	<- c( dat[1,3,1,snp], rbind( dat[2,3,,snp], dat[2,3,,snp] ) )
	ylim<- range( c( 0, y-2*y.se, y+2*y.se ), na.rm=T ) 

	ytxt		<- ylim[2] + (ylim[2]-ylim[1])*.1
	ylim[2]	<- ytxt+.01

	plot(  xlim, ylim, type='n', axes=F, ylab=paste0( 'Asthma Odds Ratio' ), main=snp, xlab='', cex.lab=1.0 )
	box()
	axis( 1, cex.axis=1.1, at=ats, lab=xlabs, las=2 )
	axis( 2, cex.axis=.9 )
	abline( h=0, col=2, lwd=1.1, lty=1 )

	points(	ats, y										, col=cols, cex=1.9, pch=16) 
	arrows(	ats, y-2*y.se, ats, y+2*y.se	, col=cols, lwd=2.5, length=0.08, angle=90, code=3)

	pvs		<- c( dat[1,'pv',1,snp], dat[2,'pv',,snp] )
	smalls<- which( pvs < .01 )
	if( length(smalls) > 0 ){
	labs2	<- rep(NA,length(pvs))
	labs2[-smalls]<- format( pvs[-smalls], nsmall=2, digits=1 )
	labs2[ smalls]<- format( pvs[ smalls], digits=1 )
	} else {
	labs2	<- format( pvs, nsmall=2, digits=1 )
	}

	labs2	<- as.expression(sapply( labs2, function(pv) bquote( p==.(pv)  ) ))
	text(		at0s,rep(ytxt,3)		, col=1, cex=1.1, lab=labs2 )
	dev.off()
} 

#### heatmap of Z-scores
source( 'fxns/my_heatmap.R' )
dat		<- array( NA, dim=c(4,2,3,length(snps)), dimnames=list( c('gwas','hom','het','glob'), c('beta','pv'), envs, snps) ) 
for( env in envs ){ 
	load( paste0( 'Rdata/compiled_gxewas/', env, '_pcxe.Rdata' ) ) 
	dat[-4,1,env,]	<- allbeta[c(1,3,4),snps]
	dat[  ,2,env,]	<- allpv	[c(1:3,5),snps]
	rm( allbeta, allpv ) 
} 

zmat	<- (
	apply( dat['het','pv',,], 2, function(x) -qnorm(x/2) )  *
	apply( dat['het','beta',,], 2, function(x) sign(x) ) ) %*% diag( sign(dat['gwas','beta',1,]) )
rownames(zmat) <- c( 'GxAtopy', 'GxEos', 'GxNeurot' )
colnames(zmat) <- snps
zmat <- zmat[,sort.list(zmat[2,])]

library(gplots)
pdf( 'figs/topsnp_betamat.pdf', w=6, h=3  )
par( mar=c(7,6,1,1) )
my_heatmap( zmat, symm=F, pv=zmat*0+1, main='', ncolor=13 )
dev.off()


