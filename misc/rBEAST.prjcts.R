rbeast.mixing<- function()
{
	require(coda)
	
	indir	<- file	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_ReductionIncidence/150417_files'
	#	read log file
	infiles	<- list.files(indir, pattern='log$')
	#	get options
	infiles	<- data.table(FILE=infiles)
	infiles[, SUBSTM:= NA_character_]
	set(infiles, infiles[, which(grepl('GTR',FILE))], 'SUBSTM', 'GTR')
	set(infiles, infiles[, which(grepl('HKY',FILE))], 'SUBSTM', 'HKY')
	infiles[, CODON:= 'No']
	set(infiles, infiles[, which(grepl('CODON',FILE))], 'CODON', 'Yes')
	infiles[, GENE:= NA_character_]
	set(infiles, infiles[, which(grepl('pol',FILE))], 'GENE', 'pol')
	infiles[, SEQ:= infiles[, substring(regmatches(FILE, regexpr('seq[0-9]+',FILE)),4)]]	
	set(infiles, NULL, 'SEQ', infiles[, as.numeric(SEQ)])
	set(infiles, infiles[, which(grepl('cseq3',FILE))], 'SEQ', 1350)
	infiles[, CLU_TYPE:= NA_character_]
	set(infiles, infiles[, which(grepl('mseq',FILE))], 'CLU_TYPE', 'largest')
	set(infiles, infiles[, which(grepl('clrnd',FILE))], 'CLU_TYPE', 'random')
	set(infiles, infiles[, which(grepl('cseq',FILE))], 'CLU_TYPE', 'all')
	set(infiles, infiles[, which(grepl('clsm',FILE))], 'CLU_TYPE', 'smallest')
	
	
	dfl	<- do.call('rbind',lapply(infiles[, unique(FILE)], function(f)
			{
				cat(paste('\nprocess file', f))
				file	<- paste(indir, '/', f, sep='')				
				z		<- beast.read.log(file, select=c('state','likelihood','posterior','branchRates','gtr','hky','site','logPopSize'), verbose=0)
				z		<- melt(z, id.vars=c('state'))
				z[, FILE:=f]
				z
			}))	
	dfm		<- dfl[, {				
						z	<- effectiveSize(matrix(data=value,ncol=1))
						list(ESS=z)
					}, by=c('FILE','variable')]
	dfm		<- merge(dfm, infiles, by='FILE')
	
	#	GTR SEQ 1000
	dfsp	<- subset(dfm, SUBSTM=='GTR' & SEQ==1000)
	setkey(dfsp, ESS)	
	ggplot(dfsp, aes(x=variable, y=ESS, fill=FILE)) + geom_bar(stat='identity') + 
			coord_flip() +
			scale_y_continuous(breaks=c(0,100,300,500,700), minor_breaks=NULL, expand=c(0,0)) +
			labs(y='effective sample size\n(20e6 iterations, burn-in 10%)', x='') +
			theme_bw() +
			theme(legend.position='bottom',panel.grid.minor = element_line(colour='grey70', size=0.2), panel.grid.major = element_line(colour='grey70', size=0.3)) 	+
			facet_grid(SUBSTM~SEQ) +
			guides(fill=guide_legend(ncol=1))
	file	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_ReductionIncidence/150417_analysis/150413_Mixing_fixedtree_pol_GTR_mseq1000.pdf'	
	ggsave(file=file, h=15, w=10)
	
	#	HKY SEQ 1000
	dfsp	<- subset(dfm, SUBSTM=='HKY' & SEQ==1000 & CLU_TYPE=='largest')
	setkey(dfsp, ESS)	
	ggplot(dfsp, aes(x=variable, y=ESS, fill=FILE)) + geom_bar(stat='identity') + 
			coord_flip() +
			scale_y_continuous(breaks=c(0,100,300,500,700), minor_breaks=NULL, expand=c(0,0)) +
			labs(y='effective sample size\n(20e6 iterations, burn-in 10%)', x='') +
			theme_bw() +
			theme(legend.position='bottom',panel.grid.minor = element_line(colour='grey70', size=0.2), panel.grid.major = element_line(colour='grey70', size=0.3)) 	+
			facet_grid(SUBSTM~SEQ) +
			guides(fill=guide_legend(ncol=1))
	file	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_ReductionIncidence/150417_analysis/150413_Mixing_fixedtree_pol_HKY_mseq1000.pdf'	
	ggsave(file=file, h=15, w=10)
}

rbeast.mixing.uweight140420<- function()
{
	require(rBEAST)
	require(coda)
	
	indir	<- file	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_ReductionIncidence/150419_files'
	#	read log file
	infiles	<- list.files(indir, pattern='log$')
	#	get options
	infiles	<- data.table(FILE=infiles)
	infiles[, SUBSTM:= NA_character_]
	set(infiles, infiles[, which(grepl('GTR',FILE))], 'SUBSTM', 'GTR')
	set(infiles, infiles[, which(grepl('HKY',FILE))], 'SUBSTM', 'HKY')
	infiles[, CODON:= 'No']
	set(infiles, infiles[, which(grepl('CODON',FILE))], 'CODON', 'Yes')
	infiles[, GENE:= NA_character_]
	set(infiles, infiles[, which(grepl('pol',FILE))], 'GENE', 'pol')
	infiles[, SEQ:= infiles[, substring(regmatches(FILE, regexpr('seq[0-9]+',FILE)),4)]]	
	set(infiles, NULL, 'SEQ', infiles[, as.numeric(SEQ)])
	set(infiles, infiles[, which(grepl('cseq3',FILE))], 'SEQ', 1350)
	infiles[, CLU_TYPE:= NA_character_]
	set(infiles, infiles[, which(grepl('mseq',FILE))], 'CLU_TYPE', 'largest')
	set(infiles, infiles[, which(grepl('clrnd',FILE))], 'CLU_TYPE', 'random')
	set(infiles, infiles[, which(grepl('cseq',FILE))], 'CLU_TYPE', 'all')
	set(infiles, infiles[, which(grepl('clsm',FILE))], 'CLU_TYPE', 'smallest')
	infiles[, WGHT_TYPE:= NA_character_]
	set(infiles, infiles[, which(grepl('uweight',FILE))], 'WGHT_TYPE', 'same')
	set(infiles, infiles[, which(grepl('uwnoSkgrdBlck',FILE))], 'WGHT_TYPE', 'same ex skygrid')
	set(infiles, infiles[, which(grepl('uwnoSkgrdBlckBrRtCtgrs',FILE))], 'WGHT_TYPE', 'same ex skygrid and branch rates')	
	infiles[, WGHT:= infiles[, regmatches(FILE, regexpr('[0-9]+\\.log',FILE))]]
	set(infiles, NULL, 'WGHT', infiles[, substr(WGHT, 1, nchar(WGHT)-4)])
		
	dfl	<- do.call('rbind',lapply(infiles[, unique(FILE)], function(f)
					{
						cat(paste('\nprocess file', f))
						file	<- paste(indir, '/', f, sep='')				
						z		<- beast.read.log(file, select=c('state','likelihood','posterior','ucld','branchRates','gtr','hky','site','logPopSize'), verbose=0)
						z		<- melt(z, id.vars=c('state'))
						z[, FILE:=f]
						z
					}))		
	dfl	<- subset(dfl, state>2e6)
	dfm		<- dfl[, {				
				z	<- effectiveSize(matrix(data=value,ncol=1))
				list(ESS=z)
			}, by=c('FILE','variable')]
	dfm		<- merge(dfm, infiles, by='FILE')
	dfm[, WGHT_LGND:= paste(WGHT_TYPE, WGHT, sep='\n')]
	#	SEQ 400
	dfsp	<- subset(dfm, SEQ==400)
	setkey(dfsp, ESS)	
	ggplot(dfsp, aes(x=variable, y=ESS, fill=FILE)) + geom_bar(stat='identity') + 
			coord_flip() +
			scale_y_continuous(breaks=c(0,100,300,500,700), minor_breaks=NULL, expand=c(0,0)) +
			labs(y='effective sample size\n(20e6 iterations, burn-in 10%)', x='') +
			theme_bw() +
			theme(legend.position='bottom',panel.grid.minor = element_line(colour='grey70', size=0.2), panel.grid.major = element_line(colour='grey70', size=0.3)) 	+
			facet_grid(SEQ~WGHT_LGND) +
			guides(fill=guide_legend(ncol=1))
	file	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_ReductionIncidence/150419_files/150420_Weights_fixedtree_pol_GTR_mseq400.pdf'	
	ggsave(file=file, h=20, w=30)
	
	#	SEQ 1000
	dfsp	<- subset(dfm, SEQ==1000)
	setkey(dfsp, ESS)	
	ggplot(dfsp, aes(x=variable, y=ESS, fill=FILE)) + geom_bar(stat='identity') + 
			coord_flip() +
			scale_y_continuous(breaks=c(0,100,300,500,700), minor_breaks=NULL, expand=c(0,0)) +
			labs(y='effective sample size\n(20e6 iterations, burn-in 10%)', x='') +
			theme_bw() +
			theme(legend.position='bottom',panel.grid.minor = element_line(colour='grey70', size=0.2), panel.grid.major = element_line(colour='grey70', size=0.3)) 	+
			facet_grid(SEQ~WGHT_LGND) +
			guides(fill=guide_legend(ncol=1))
	file	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_ReductionIncidence/150419_files/150420_Weights_fixedtree_pol_GTR_mseq1000.pdf'	
	ggsave(file=file, h=20, w=30)
}

rbeast.mixing.for.commonparameters<- function()
{
	require(coda)
	
	indir	<- '/Users/Oliver/duke/2014_Gates/methods_comparison_pipeline/150414'
	#	read log file
	infiles	<- list.files(indir, pattern='log$')
	#	get options
	infiles	<- data.table(FILE=infiles)
	infiles[, SUBSTM:= NA_character_]
	set(infiles, infiles[, which(grepl('GTR',FILE))], 'SUBSTM', 'GTR')
	set(infiles, infiles[, which(grepl('HKY',FILE))], 'SUBSTM', 'HKY')
	infiles[, CODON:= 'No']
	set(infiles, infiles[, which(grepl('CODON',FILE))], 'CODON', 'Yes')
	infiles[, GENE:= NA_character_]
	set(infiles, infiles[, which(grepl('pol',FILE))], 'GENE', 'pol')
	infiles[, SEQ:= infiles[, substring(regmatches(FILE, regexpr('seq[0-9]+',FILE)),4)]]	
	set(infiles, NULL, 'SEQ', infiles[, as.numeric(SEQ)])
	set(infiles, infiles[, which(grepl('cseq3',FILE))], 'SEQ', 1350)
	infiles[, CLU_TYPE:= NA_character_]
	set(infiles, infiles[, which(grepl('mseq',FILE))], 'CLU_TYPE', 'largest')
	set(infiles, infiles[, which(grepl('clrnd',FILE))], 'CLU_TYPE', 'random')
	set(infiles, infiles[, which(grepl('cseq',FILE))], 'CLU_TYPE', 'all')
	set(infiles, infiles[, which(grepl('clsm',FILE))], 'CLU_TYPE', 'smallest')
	
	dfl	<- infiles[, {
				cat(paste('\nprocess file', FILE))
				file	<- paste(indir, '/', FILE, sep='')				
				beast.read.log(file, select=c('state','likelihood','posterior','branchRates'), verbose=0)				
			}, by='FILE']
	dfl		<- subset(dfl, state>2e6)
	
	dfm		<- do.call('rbind',lapply( dfl[, unique(FILE)], function(f)
			{
				tmp	<- subset(dfl, FILE==f)
				z	<- effectiveSize(subset(tmp, select=3:ncol(tmp)))
				data.table(FILE=f, VAR=names(z), ESS=z)
			}))
	dfm		<- merge(dfm, infiles, by='FILE')
	
	#	GTR model vs HKY model
	dfsp	<- subset(dfm, CLU_TYPE=='largest')	
	setkey(dfsp, SEQ)
	ggplot(dfsp, aes(x=VAR, y=ESS, fill=FILE)) + geom_bar(stat='identity') + 
			coord_flip() +
			#scale_y_log10() +
			scale_y_continuous(breaks=c(0,100,300,500,700), minor_breaks=NULL, expand=c(0,0)) +
			labs(y='effective sample size\n(20e6 iterations, burn-in 10%)', x='') +
			theme_bw() +
			theme(legend.position='bottom',panel.grid.minor = element_line(colour='grey70', size=0.2), panel.grid.major = element_line(colour='grey70', size=0.3)) 	+
			facet_grid(SUBSTM~SEQ) +
			guides(fill=guide_legend(ncol=1))
	file	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_ReductionIncidence/150413_Mixing_fixedtree_pol_BY_seq_substm.pdf'	
	ggsave(file=file, h=6, w=12)
}

rbeast.skygrid<- function()
{
	indir	<- '/Users/Oliver/duke/2014_Gates/methods_comparison_pipeline/150414'
	#	read log file
	infiles	<- list.files(indir, pattern='log$')
	#	get options
	infiles	<- data.table(FILE=infiles)
	infiles[, SUBSTM:= NA_character_]
	set(infiles, infiles[, which(grepl('GTR',FILE))], 'SUBSTM', 'GTR')
	set(infiles, infiles[, which(grepl('HKY',FILE))], 'SUBSTM', 'HKY')
	infiles[, CODON:= 'No']
	set(infiles, infiles[, which(grepl('CODON',FILE))], 'CODON', 'Yes')
	infiles[, GENE:= NA_character_]
	set(infiles, infiles[, which(grepl('pol',FILE))], 'GENE', 'pol')
	infiles[, SEQ:= infiles[, substring(regmatches(FILE, regexpr('seq[0-9]+',FILE)),4)]]	
	set(infiles, NULL, 'SEQ', infiles[, as.numeric(SEQ)])
	set(infiles, infiles[, which(grepl('cseq3',FILE))], 'SEQ', 1350)
	infiles[, CLU_TYPE:= NA_character_]
	set(infiles, infiles[, which(grepl('mseq',FILE))], 'CLU_TYPE', 'largest')
	set(infiles, infiles[, which(grepl('clrnd',FILE))], 'CLU_TYPE', 'random')
	set(infiles, infiles[, which(grepl('cseq',FILE))], 'CLU_TYPE', 'all')
	set(infiles, infiles[, which(grepl('clsm',FILE))], 'CLU_TYPE', 'smallest')
	
	dfl	<- infiles[, {
				cat(paste('\nprocess file', FILE))
				file	<- paste(indir, '/', FILE, sep='')				
				beast.read.log(file, select=c('state','skygrid'), verbose=0)				
			}, by='FILE']
	dfl	<- melt(dfl, id.vars=c('FILE','state','skygrid.cutOff','skygrid.numGridPoints','skygrid.precision'))
	setnames(dfl, 'value', 'LPS')
	setnames(dfl, 'variable', 'X')
	set(dfl, NULL, 'X', dfl[,gsub('skygrid.logPopSize','',X)])
	set(dfl, NULL, 'X', dfl[,as.numeric(X)])
	
	dfl		<- subset(dfl, FILE!="150129_HPTN071_scA-mFP25_TEST_pol_grid-cseq3.log")
	
	max.samplingtime	<- 2020
	#dfs		<- subset(dfl, FILE=='150129_HPTN071_scA-mFP25_TEST_pol_CODON-GTR-grid-cseq3.log' & state>2e6)
	dfs		<- subset(dfl, state>2e6)
	dfs		<- dfs[, list(	skygrid.cutOff=skygrid.cutOff[1], skygrid.numGridPoints=skygrid.numGridPoints[1], 
							LPSm= median(exp(LPS)), LPSlq= quantile(exp(LPS), p=0.025), LPSuq= quantile(exp(LPS), p=0.975),
							YR=max.samplingtime-(X-1)/(skygrid.numGridPoints[1])*skygrid.cutOff[1]), by=c('FILE','X')]
	dfs		<- merge(dfs, infiles, by='FILE')
	
	#	GTR model change # sequences used
	dfsp	<- subset(dfs, SUBSTM=='GTR') 
	setkey(dfsp, SEQ, YR)
	ggplot(dfsp, aes(x=YR, fill=FILE)) + geom_ribbon(aes(ymin=LPSlq, ymax=LPSuq), alpha=0.5, colour=NA) + geom_line(aes(y=LPSm)) + 
			scale_y_log10() +
			scale_x_continuous(expand=c(0,0)) +
			labs(x='', y='effective population size\n(Skygrid)') +
			theme_bw() +
			theme(legend.position='bottom',panel.grid.minor = element_line(colour='grey70', size=0.2), panel.grid.major = element_line(colour='grey70', size=0.3)) 	+
			facet_grid(~SEQ) +
			guides(fill=guide_legend(ncol=1))
	file	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_ReductionIncidence/150413_Skygrid_CODONGTR_fixedtree_pol_BY_seq.pdf'	
	ggsave(file=file, h=6, w=10)
	
	#	GTR model vs HKY model
	dfsp	<- subset(dfs, CLU_TYPE=='largest')
	setkey(dfsp, SEQ, YR)
	ggplot(dfsp, aes(x=YR, fill=FILE)) + geom_ribbon(aes(ymin=LPSlq, ymax=LPSuq), alpha=0.5, colour=NA) + geom_line(aes(y=LPSm)) + 
			scale_y_log10() +
			scale_x_continuous(expand=c(0,0)) +
			labs(x='', y='effective population size\n(Skygrid)') +
			theme_bw() +
			theme(legend.position='bottom',panel.grid.minor = element_line(colour='grey70', size=0.2), panel.grid.major = element_line(colour='grey70', size=0.3)) 	+
			facet_grid(SUBSTM~SEQ) +
			guides(fill=guide_legend(ncol=1))
	file	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_ReductionIncidence/150413_Skygrid_fixedtree_pol_BY_seq_substm.pdf'	
	ggsave(file=file, h=10, w=10)
	
	#	HKY model clustering type
	dfsp	<- subset(dfs, SUBSTM=='HKY' & SEQ>500 & SEQ<1300)
	setkey(dfsp, SEQ, YR)
	ggplot(dfsp, aes(x=YR, fill=FILE)) + geom_ribbon(aes(ymin=LPSlq, ymax=LPSuq), alpha=0.5, colour=NA) + geom_line(aes(y=LPSm)) + 
			scale_y_log10() +
			scale_x_continuous(expand=c(0,0)) +
			labs(x='', y='effective population size\n(Skygrid)') +
			theme_bw() +
			theme(legend.position='bottom',panel.grid.minor = element_line(colour='grey70', size=0.2), panel.grid.major = element_line(colour='grey70', size=0.3)) 	+
			facet_grid(CLU_TYPE~SEQ) +
			guides(fill=guide_legend(ncol=1))
	file	<- '/Users/Oliver/Dropbox (Infectious Disease)/OR_Work/2015/2015_PANGEA_ReductionIncidence/150413_Skygrid_HKY_fixedtree_pol_BY_clutype.pdf'	
	ggsave(file=file, h=10, w=7)
	
}

rbeast.skygrid.hky.weight.allexceptlogpopsize.140421<- function()
{
	require(rBEAST)
	require(data.table)
	require(ape)
	require(XML)
	tree.id.labelsep		<- '|'
	tree.id.label.idx.ctime	<- 4 	
	#select		<- 'grid-mseq500'
	#select		<- 'grid-cseq3'
	select		<- 'grid-mseq1000'
	indir		<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_ReductionIncidence/150419_files'
	#indir		<- '/Users/Oliver/git/HPTN071sim/tmp140914/140716_RUN001_INTERNAL'  
	outdir		<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_ReductionIncidence/150421_files'
	infiles		<- list.files(indir, '.*INTERNAL.R$', full.names=FALSE)
	#stopifnot(length(infiles)==1)
	
	#	run: no UpDown operator
	selects		<- paste('grid-mseq',c(400, 1000), sep='')
	for(w in c(5, 10, 20, 50))
	{			
		for(select in selects)
		{
			for(i in seq_along(infiles))
			{
				infile			<- infiles[i]
				#	load simulated data
				file			<- paste(indir, '/', infile, sep='')
				file.name		<- paste(outdir, gsub('_SIMULATED_INTERNAL.R',paste('_TEST_pol_HKY_fixedtree_',select,"_w4UCLD",w,'.xml',sep=''),infile), sep='/')
				cat(paste('\nLoading file', file))
				load(file)		#expect "df.epi"    "df.trms"   "df.inds"   "df.sample" "df.seq"
				set( df.seq, NULL, 'IDCLU', df.seq[, as.integer(IDCLU)] )
				setnames(df.seq, c("LABEL", "IDCLU", "IDPOP"), c("TAXON_NAME", "CLU_ID", "TAXON_ID"))
				df.seq[, SAMPLINGTIME:=df.seq[, as.numeric(sapply(strsplit(TAXON_NAME,'|',fixed=1),'[[',4))]]	
				seq.select		<- beast.choose.seq.by.clusters(df.seq, select, verbose=1)			
				#
				tmp				<- list.files(indir, '_DATEDTREE.newick$', full.names=FALSE)
				tmp				<- tmp[ grepl(substr(infile, 1, regexpr('_SIMULATED',infile)), tmp) ]
				phd				<- read.tree(paste(indir, tmp, sep='/'))
				tmp2			<- data.table(TAXON_ID= sapply(phd, function(x) x$tip.label[1]), IDX=seq_along(phd))
				set( tmp2, NULL, 'TAXON_ID', tmp2[, as.integer(substring(sapply(strsplit(TAXON_ID, tree.id.labelsep, fixed=TRUE),'[[',1),7)) ] )
				tmp2			<- merge(subset(seq.select, select=c(TAXON_ID, CLU_ID)), tmp2, by='TAXON_ID')			
				phd				<- lapply(tmp2[,IDX], function(i) phd[[i]] )
				names(phd)		<- tmp2[, CLU_ID]
				#
				#	create BEAST XML POL
				#
				seq.select.pol	<- subset(seq.select, select=c("CLU_ID", "TAXON_ID", "TAXON_NAME", "SAMPLINGTIME", "POL" ))
				setnames(seq.select.pol, 'POL', 'SEQ')									
				bxml			<- beastxml.multilocus.hky( file.name, seq.select.pol, phd, verbose=1 )
				#
				#	adjust parameter weights: all 1 but higher for UCLD
				#
				bxmlo		<- getNodeSet(bxml, "//operators")[[1]]			
				tmp			<- unlist(lapply( c("hky","site","skygrid.precision","skygrid","branchRates.categories_CLU"), function(x)
								{															
									getNodeSet(bxmlo, paste('*[descendant::*[contains(@idref,"',x,'")]]',sep=''))
								}))
				for(x in tmp)
					xmlAttrs(x)["weight"]	<- 1
				tmp			<- unlist(lapply( c("ucld"), function(x)
								{															
									getNodeSet(bxmlo, paste('*[descendant::*[contains(@idref,"',x,'")]]',sep=''))
								}))
				for(x in tmp)
					xmlAttrs(x)["weight"]	<- w
				
				cat(paste("\nwrite xml file to",file.name))
				saveXML(bxml, file=file.name)				
			}
		}
	}
	#	run: with UpDown operator
	selects		<- paste('grid-mseq',c(400, 1000), sep='')
	for(w in c(5, 10, 20, 50))
	{			
		for(select in selects)
		{
			for(i in seq_along(infiles))
			{
				infile			<- infiles[i]
				#	load simulated data
				file			<- paste(indir, '/', infile, sep='')
				file.name		<- paste(outdir, gsub('_SIMULATED_INTERNAL.R',paste('_TEST_pol_HKY_fixedtree_',select,"_w4UCLDUpDown",w,'.xml',sep=''),infile), sep='/')
				cat(paste('\nLoading file', file))
				load(file)		#expect "df.epi"    "df.trms"   "df.inds"   "df.sample" "df.seq"
				set( df.seq, NULL, 'IDCLU', df.seq[, as.integer(IDCLU)] )
				setnames(df.seq, c("LABEL", "IDCLU", "IDPOP"), c("TAXON_NAME", "CLU_ID", "TAXON_ID"))
				df.seq[, SAMPLINGTIME:=df.seq[, as.numeric(sapply(strsplit(TAXON_NAME,'|',fixed=1),'[[',4))]]	
				seq.select		<- beast.choose.seq.by.clusters(df.seq, select, verbose=1)			
				#
				tmp				<- list.files(indir, '_DATEDTREE.newick$', full.names=FALSE)
				tmp				<- tmp[ grepl(substr(infile, 1, regexpr('_SIMULATED',infile)), tmp) ]
				phd				<- read.tree(paste(indir, tmp, sep='/'))
				tmp2			<- data.table(TAXON_ID= sapply(phd, function(x) x$tip.label[1]), IDX=seq_along(phd))
				set( tmp2, NULL, 'TAXON_ID', tmp2[, as.integer(substring(sapply(strsplit(TAXON_ID, tree.id.labelsep, fixed=TRUE),'[[',1),7)) ] )
				tmp2			<- merge(subset(seq.select, select=c(TAXON_ID, CLU_ID)), tmp2, by='TAXON_ID')			
				phd				<- lapply(tmp2[,IDX], function(i) phd[[i]] )
				names(phd)		<- tmp2[, CLU_ID]
				#
				#	create BEAST XML POL
				#
				seq.select.pol	<- subset(seq.select, select=c("CLU_ID", "TAXON_ID", "TAXON_NAME", "SAMPLINGTIME", "POL" ))
				setnames(seq.select.pol, 'POL', 'SEQ')									
				bxml			<- beastxml.multilocus.hky( file.name, seq.select.pol, phd, verbose=1 )
				#
				#	adjust parameter weights: all 1 but higher for UCLD and add UpDown operator
				#
				beast.add.upDownOperator(bxml, "ucld.mean", "ucld.stdev", scaleFactor=0.75, weight=1)
				bxmlo		<- getNodeSet(bxml, "//operators")[[1]]		
				#
				tmp			<- unlist(lapply( c("hky","site","skygrid.precision","skygrid","branchRates.categories_CLU"), function(x)
								{															
									getNodeSet(bxmlo, paste('*[descendant::*[contains(@idref,"',x,'")]]',sep=''))
								}))
				for(x in tmp)
					xmlAttrs(x)["weight"]	<- 1
				tmp			<- unlist(lapply( c("ucld"), function(x)
								{															
									getNodeSet(bxmlo, paste('*[descendant::*[contains(@idref,"',x,'")]]',sep=''))
								}))
				for(x in tmp)
					xmlAttrs(x)["weight"]	<- w				
				cat(paste("\nwrite xml file to",file.name))
				saveXML(bxml, file=file.name)				
			}
		}
	}
}

rbeast.skygrid.hky.weight.allexceptlogpopsize.140420<- function()
{
	require(rBEAST)
	require(data.table)
	require(ape)
	require(XML)
	tree.id.labelsep		<- '|'
	tree.id.label.idx.ctime	<- 4 	
	#select		<- 'grid-mseq500'
	#select		<- 'grid-cseq3'
	select		<- 'grid-mseq1000'
	indir		<- '/Users/Oliver/Dropbox\ (Infectious Disease)/OR_Work/2015/2015_PANGEA_ReductionIncidence/150419_files'
	#indir		<- '/Users/Oliver/git/HPTN071sim/tmp140914/140716_RUN001_INTERNAL'  
	outdir		<- indir
	infiles		<- list.files(indir, '.*INTERNAL.R$', full.names=FALSE)
	#stopifnot(length(infiles)==1)
	#
	#	run  
	#
	selects		<- paste('grid-mseq',c(400, 1000), sep='')
	for(w in c(1, 5))
	{			
		for(select in selects)
		{
			for(i in seq_along(infiles))
			{
				infile			<- infiles[i]
				#	load simulated data
				file			<- paste(indir, '/', infile, sep='')
				file.name		<- paste(indir, gsub('_SIMULATED_INTERNAL.R',paste('_TEST_pol_HKY_fixedtree_',select,"_uweight",w,'.xml',sep=''),infile), sep='/')
				cat(paste('\nLoading file', file))
				load(file)		#expect "df.epi"    "df.trms"   "df.inds"   "df.sample" "df.seq"
				set( df.seq, NULL, 'IDCLU', df.seq[, as.integer(IDCLU)] )
				setnames(df.seq, c("LABEL", "IDCLU", "IDPOP"), c("TAXON_NAME", "CLU_ID", "TAXON_ID"))
				df.seq[, SAMPLINGTIME:=df.seq[, as.numeric(sapply(strsplit(TAXON_NAME,'|',fixed=1),'[[',4))]]	
				seq.select		<- beast.choose.seq.by.clusters(df.seq, select, verbose=1)			
				#
				tmp				<- list.files(indir, '_DATEDTREE.newick$', full.names=FALSE)
				tmp				<- tmp[ grepl(substr(infile, 1, regexpr('_SIMULATED',infile)), tmp) ]
				phd				<- read.tree(paste(indir, tmp, sep='/'))
				tmp2			<- data.table(TAXON_ID= sapply(phd, function(x) x$tip.label[1]), IDX=seq_along(phd))
				set( tmp2, NULL, 'TAXON_ID', tmp2[, as.integer(substring(sapply(strsplit(TAXON_ID, tree.id.labelsep, fixed=TRUE),'[[',1),7)) ] )
				tmp2			<- merge(subset(seq.select, select=c(TAXON_ID, CLU_ID)), tmp2, by='TAXON_ID')			
				phd				<- lapply(tmp2[,IDX], function(i) phd[[i]] )
				names(phd)		<- tmp2[, CLU_ID]
				#
				#	create BEAST XML POL
				#
				seq.select.pol	<- subset(seq.select, select=c("CLU_ID", "TAXON_ID", "TAXON_NAME", "SAMPLINGTIME", "POL" ))
				setnames(seq.select.pol, 'POL', 'SEQ')									
				bxml			<- beastxml.multilocus.hky( file.name, seq.select.pol, phd, verbose=1 )
				#
				#	adjust parameter weights: all the same
				#
				bxmlo		<- getNodeSet(bxml, "//operators")[[1]]			
				tmp			<- unlist(lapply( c("hky","site","ucld","skygrid.precision","skygrid","branchRates.categories_CLU"), function(x)
						{															
							getNodeSet(bxmlo, paste('*[descendant::*[contains(@idref,"',x,'")]]',sep=''))
						}))
				for(x in tmp)
					xmlAttrs(x)["weight"]	<- w
				cat(paste("\nwrite xml file to",file.name))
				saveXML(bxml, file=file.name)				
			}
		}
	}
	#	weights all the same except skygrid block updater
	for(w in c(1, 5, 20, 50))
	{	
		for(select in selects)
		{
			for(i in seq_along(infiles))
			{
				infile			<- infiles[i]
				#	load simulated data
				file			<- paste(indir, '/', infile, sep='')
				file.name		<- paste(indir, gsub('_SIMULATED_INTERNAL.R',paste('_TEST_pol_HKY_fixedtree_',select,"_uwnoSkgrdBlck",w,'.xml',sep=''),infile), sep='/')
				cat(paste('\nLoading file', file))
				load(file)		#expect "df.epi"    "df.trms"   "df.inds"   "df.sample" "df.seq"
				set( df.seq, NULL, 'IDCLU', df.seq[, as.integer(IDCLU)] )
				setnames(df.seq, c("LABEL", "IDCLU", "IDPOP"), c("TAXON_NAME", "CLU_ID", "TAXON_ID"))
				df.seq[, SAMPLINGTIME:=df.seq[, as.numeric(sapply(strsplit(TAXON_NAME,'|',fixed=1),'[[',4))]]	
				seq.select		<- beast.choose.seq.by.clusters(df.seq, select, verbose=1)			
				#
				tmp				<- list.files(indir, '_DATEDTREE.newick$', full.names=FALSE)
				tmp				<- tmp[ grepl(substr(infile, 1, regexpr('_SIMULATED',infile)), tmp) ]
				phd				<- read.tree(paste(indir, tmp, sep='/'))
				tmp2			<- data.table(TAXON_ID= sapply(phd, function(x) x$tip.label[1]), IDX=seq_along(phd))
				set( tmp2, NULL, 'TAXON_ID', tmp2[, as.integer(substring(sapply(strsplit(TAXON_ID, tree.id.labelsep, fixed=TRUE),'[[',1),7)) ] )
				tmp2			<- merge(subset(seq.select, select=c(TAXON_ID, CLU_ID)), tmp2, by='TAXON_ID')			
				phd				<- lapply(tmp2[,IDX], function(i) phd[[i]] )
				names(phd)		<- tmp2[, CLU_ID]
				#
				#	create BEAST XML POL
				#
				seq.select.pol	<- subset(seq.select, select=c("CLU_ID", "TAXON_ID", "TAXON_NAME", "SAMPLINGTIME", "POL" ))
				setnames(seq.select.pol, 'POL', 'SEQ')									
				bxml			<- beastxml.multilocus.hky( file.name, seq.select.pol, phd, verbose=1 )
				#
				#	adjust parameter weights: all the same
				#
				bxmlo		<- getNodeSet(bxml, "//operators")[[1]]			
				tmp			<- unlist(lapply( c("hky","site","ucld","skygrid.precision","branchRates.categories_CLU"), function(x)
								{															
									getNodeSet(bxmlo, paste('*[descendant::*[contains(@idref,"',x,'")]]',sep=''))
								}))
				for(x in tmp)
					xmlAttrs(x)["weight"]	<- w
				cat(paste("\nwrite xml file to",file.name))
				saveXML(bxml, file=file.name)				
			}
		}
	}
	#	weights all the same except skygrid block updater
	for(w in c(1, 5, 20, 50))
	{
		for(select in selects)
		{
			for(i in seq_along(infiles))
			{
				infile			<- infiles[i]
				#	load simulated data
				file			<- paste(indir, '/', infile, sep='')
				file.name		<- paste(indir, gsub('_SIMULATED_INTERNAL.R',paste('_TEST_pol_HKY_fixedtree_',select,"_uwnoSkgrdBlckBrRtCtgrs",w,'.xml',sep=''),infile), sep='/')
				cat(paste('\nLoading file', file))
				load(file)		#expect "df.epi"    "df.trms"   "df.inds"   "df.sample" "df.seq"
				set( df.seq, NULL, 'IDCLU', df.seq[, as.integer(IDCLU)] )
				setnames(df.seq, c("LABEL", "IDCLU", "IDPOP"), c("TAXON_NAME", "CLU_ID", "TAXON_ID"))
				df.seq[, SAMPLINGTIME:=df.seq[, as.numeric(sapply(strsplit(TAXON_NAME,'|',fixed=1),'[[',4))]]	
				seq.select		<- beast.choose.seq.by.clusters(df.seq, select, verbose=1)			
				#
				tmp				<- list.files(indir, '_DATEDTREE.newick$', full.names=FALSE)
				tmp				<- tmp[ grepl(substr(infile, 1, regexpr('_SIMULATED',infile)), tmp) ]
				phd				<- read.tree(paste(indir, tmp, sep='/'))
				tmp2			<- data.table(TAXON_ID= sapply(phd, function(x) x$tip.label[1]), IDX=seq_along(phd))
				set( tmp2, NULL, 'TAXON_ID', tmp2[, as.integer(substring(sapply(strsplit(TAXON_ID, tree.id.labelsep, fixed=TRUE),'[[',1),7)) ] )
				tmp2			<- merge(subset(seq.select, select=c(TAXON_ID, CLU_ID)), tmp2, by='TAXON_ID')			
				phd				<- lapply(tmp2[,IDX], function(i) phd[[i]] )
				names(phd)		<- tmp2[, CLU_ID]
				#
				#	create BEAST XML POL
				#
				seq.select.pol	<- subset(seq.select, select=c("CLU_ID", "TAXON_ID", "TAXON_NAME", "SAMPLINGTIME", "POL" ))
				setnames(seq.select.pol, 'POL', 'SEQ')									
				bxml			<- beastxml.multilocus.hky( file.name, seq.select.pol, phd, verbose=1 )
				#
				#	adjust parameter weights: all the same
				#
				bxmlo		<- getNodeSet(bxml, "//operators")[[1]]			
				tmp			<- unlist(lapply( c("hky","site","ucld","skygrid.precision"), function(x)
								{															
									getNodeSet(bxmlo, paste('*[descendant::*[contains(@idref,"',x,'")]]',sep=''))
								}))
				for(x in tmp)
					xmlAttrs(x)["weight"]	<- w
				cat(paste("\nwrite xml file to",file.name))
				saveXML(bxml, file=file.name)				
			}
		}
	}
}