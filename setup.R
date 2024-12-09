rm( list=ls() )
if( file.exists( 'Rdata/setup.Rdata' ) ) stop()

ukbdir	<- '/u/project/sriram/nadavrap/UKBB/' 
maxN		<- 342815 
nPC			<- 10 ## cf Pividori et al on asthma AoO

envs	<- c( 'disease_ALLERGY_ECZEMA_DIAGNOSED', 'eo_qn', 'neuroticism' )
loadnames	<- c( 'disease_ALLERGY_ECZEMA_DIAGNOSED'	,'blood_EOSINOPHIL_COUNT', 'neuroticism' )
niceEnvs	<- c( 'Atopy', 'Eosinophils', 'Neuroticism') 
names(loadnames)	<- names(niceEnvs)	<- envs

## misc flags to handle the way UKB is stored on our cluster
new_phens<- c( 'neuroticism', 'asthma_onset' )
discrEnv <- c( 'disease_ALLERGY_ECZEMA_DIAGNOSED' )
quantEnv <- c( 'eo_qn', 'neuroticism' ) 

prs	<- c( 'disease_ASTHMA_DIAGNOSED', 'disease_ALLERGY_ECZEMA_DIAGNOSED', 'blood_EOSINOPHIL_COUNT', 'mental_NEUROTICISM' )
prs0<- c( 'Asthma', 'Atopy', 'Eosinophils', 'Neuroticism' ) 
np	<- length(prs)
names(prs0)	<- prs 

save.image( file='Rdata/setup.Rdata' )
