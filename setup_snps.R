rm( list=ls() ) 
if( file.exists( 'Rdata/setup_snps.Rdata' ) ) stop() 
if( ! file.exists( 'output/pruned.snplist' ) )
	system('plink --threads 1 --bfile /u/project/zaitlenlab/nadavrap/UKBB/GEN/ukbb --write-snplist --keep /u/project/sriram/nadavrap/UKBB/CSV/final_covars.txt --maf .01 --hwe 1e-12 --chr 1-22 --geno .1 --out output/pruned')

snpnames0	<- as.character(read.table( 'output/pruned.snplist' )[,1])

chrnames	<- as.numeric(as.character( read.table( '/u/project/zaitlenlab/nadavrap/UKBB/GEN/ukbb.bim' )[,1] ))
snpnames	<-						as.character( read.table( '/u/project/zaitlenlab/nadavrap/UKBB/GEN/ukbb.bim' )[,2] )

sub	<- which( snpnames %in% snpnames0 )
chrnames	<- chrnames [sub]
snpnames	<- snpnames [sub]
names(chrnames)	<- snpnames

nsnp	<- length( snpnames )

nperbatch	<- 200
nbatch		<- ceiling( nsnp/nperbatch ) 

target_beds	<- c( 
	rep( NA, 10 ),
	'parsed_data/locus_rs13194312_100',
	'parsed_data/locus_rs2033784_100',
	'parsed_data/locus_rs2396255_100',
	'parsed_data/locus_rs2544523_100',
	'parsed_data/locus_rs4129267_100',
	'parsed_data/locus_rs7593948_100',
	'parsed_data/locus_rs992969_100'
) 

keysnps	<- c( 
	rep( NA, 10 ),
	'rs13194312',
	'rs2033784',
	'rs2396255',
	'rs2544523',
	'rs4129267',
	'rs7593948',
	'rs992969'
) 

save.image( file='Rdata/setup_snps.Rdata' ) 
