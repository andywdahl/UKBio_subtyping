library(rgwas)
## below are simplified versions of code used to load data in Sheppard et al 2021

getCovars <- function(
	fname=paste0(ukbdir,'CSV/final_covars.txt'), cachefname='output/final_covars.rds',
	keep_bmi=FALSE,keep_batch=FALSE,keep_center=FALSE,nPC
){ 
  
  if(file.exists(cachefname)) return(readRDS(cachefname)) 

  covs <- read.table(fname, header = TRUE, sep='')

  covs$female				<- 1-as.numeric(covs$sex_f31_0_0)
  covs$sex_f31_0_0	<- NULL

	if( keep_batch )
		covs$batch	<- as.factor(covs$genotype_measurement_batch_f22000_0_0)
	covs$genotype_measurement_batch_f22000_0_0	<- NULL

	if( keep_center )
		covs$centre		<- as.factor(covs$uk_biobank_assessment_centre_f54_0_0)
	covs$uk_biobank_assessment_centre_f54_0_0		<- NULL

  covs$age <- as.numeric(covs$age_when_attended_assessment_centre_f21003_0_0)
  covs$age2<- ( covs$age - mean(covs$age,na.rm=T) )^2
	covs$age_when_attended_assessment_centre_f21003_0_0	<- NULL

	if(keep_bmi)
		covs$BMI <- apply(covs[,grep('body_mass_index_bmi_f21001_', colnames(covs)),drop=F], 1, function(x) median(x, na.rm = TRUE)) # median if multiple measurments
	covs[,grep('body_mass_index_bmi_f21001_', colnames(covs))] <- NULL 

	colnames( covs )[ colnames(covs) %in% paste0( 'genetic_principal_components_f22009_0_', 1:nPC ) ]	<- paste0( 'PC_', 1:nPC )

	bads	<- which( apply( covs[,grep('genotype_measurement_batch_f22000_0_0', colnames(covs))], 1, function(x) any(x==-9) ) )
	if( length(bads) > 0 ) stop( length(bads) )

  saveRDS(covs, file=cachefname)
  covs
} 

getPheno <- function(phen,binary,qn=FALSE,scale01=(!binary)){ 
	if( phen %in% new_phens ){
		prsfile		<- paste0(ukbdir, '/phenotypes/', phen, '.pheno' )
		y					<- read.table(prsfile,head=T)
		if( phen %in% 'asthma_onset' ){
			y[ which(y[,3]==-3), -(1:2) ]	<- NA
			y[ which(y[,3]==-1), -(1:2) ]	<- NA
		}
		y			<- as.matrix( y[,1:3] )
	} else {
		loadfile<- paste0( ukbdir, phen, '/', phen, '.pheno' )
		y				<- read.table(loadfile,head=T)[,1:3]
	}
	colnames(y)		<- c( 'FID', 'IID', phen )

	if( phen %in% c( 'disease_ALLERGY_ECZEMA_DIAGNOSED', 'disease_ASTHMA_DIAGNOSED' ) )
		y[ y[,3]==-9, 3 ]	<- NA 
	if( binary )
		y[,3]	<- y[,3] - 1 
	if( qn )
		y[,3] <- as.numeric(phenix::quantnorm(as.matrix(y[,3]))) 
	if( scale01 )
		y[,3] <- rgwas:::scale01( y[,3] ) 
	y
}
