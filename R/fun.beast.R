
######################################################################################
#	get BEAST taxon labels:		cluster	FASTASampleCode	NegT	AnyPosT	SeqT  -> turn date into numerical format
#' @export
beast.addBEASTLabel<- function( df, df.resetTipDate=NA )
{
	if(is.na(df.resetTipDate))
		df[,dummy:=NA]
	else if(df.resetTipDate=="LdTd")
		df	<- merge(df, df[,list(dummy=max(AnyPos_T1)),by="cluster"])
	else if(df.resetTipDate=="LsTd")
		df	<- merge(df, df[,list(dummy=max(PosSeqT)),by="cluster"])
	else if(df.resetTipDate=="UmTd")
		df[,dummy:=max(df[,PosSeqT])]
	else
		stop("Unexpected df.resetTipDate")
	df	<- merge(df, df[,list(PosSeqTadj= max(PosSeqT, dummy, na.rm=T)),by="FASTASampleCode"], by="FASTASampleCode")
	
	tmp	<- df[,		{
				z	<- as.POSIXlt(c(NegT, AnyPos_T1, PosSeqT, PosSeqTadj))
				tmp	<- z$year + 1900
				z	<- tmp + round( z$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )												
				list(BEASTlabel= paste(c(cluster, z, FASTASampleCode), collapse='_', sep=''))
			}, by="FASTASampleCode"]
	df	<- merge(df, tmp, by="FASTASampleCode")
	subset(df,select=-which(colnames(df)=="dummy"))	
}
######################################################################################
#	pool clusters into sets containing roughly 'pool.ntip' sequences
#' @export
beast.poolclusters<- function(cluphy.df, pool.ntip= 130, pool.includealwaysbeforeyear=NA, verbose=1)
{	
	df			<- cluphy.df[, list(clu.ntip=clu.ntip[1], clu.AnyPos_T1=clu.AnyPos_T1[1]), by="cluster"]	
	if(!is.na(pool.includealwaysbeforeyear))
	{
		if(verbose) cat(paste("\nalways include clusters starting before ",pool.includealwaysbeforeyear,"and then pool evenly across clu.AnyPos_T1"))
		df.always	<- subset(df,as.POSIXlt(clu.AnyPos_T1)$year<(pool.includealwaysbeforeyear-1900))
		df			<- subset(df,as.POSIXlt(clu.AnyPos_T1)$year>=(pool.includealwaysbeforeyear-1900))
		pool.ntip	<- pool.ntip - sum(df.always[,clu.ntip])		
		setkey(df, clu.AnyPos_T1)
		pool.n		<- ceiling( sum( df[,clu.ntip] ) / pool.ntip )
		tmp			<- lapply( seq_len(pool.n), function(x)	seq.int(x,nrow(df),by=pool.n) )		
		pool.df		<- lapply(seq_along(tmp), function(i) merge(subset(rbind(df.always, df[tmp[[i]],]), select=cluster), cluphy.df, by="cluster") )		
	}
	else
	{
		if(verbose) cat(paste("\npool evenly across clu.AnyPos_T1"))
		setkey(df, clu.AnyPos_T1)
		pool.n	<- ceiling( sum( df[,clu.ntip] ) / pool.ntip )
		tmp		<- lapply( seq_len(pool.n), function(x)	seq.int(x,nrow(df),by=pool.n) )
		pool.df	<- lapply(seq_along(tmp), function(i) merge(subset(df[tmp[[i]],], select=cluster), cluphy.df, by="cluster") )		
	}
	if(verbose) cat(paste("\nnumber of pools is n=",pool.n))		
	if(verbose) cat(paste("\nnumber of seq in pools is n=",paste( sapply(pool.df, nrow), sep='', collapse=', ' )))
	list(pool.df=pool.df, pool.ntip=pool.ntip)
}
######################################################################################
#	add taxa and alignment in bxml from BEASTlabels in df and alignment in seq.PROT.RT
#	beast.label.datepos= 4; beast.label.sep= '_'; beast.date.direction= "forwards"; beast.date.units= "years"; beast.alignment.dataType= "nucleotide"; xml.resetTipDate2LastDiag=1
#' @export
beast.add.taxa<- function(bxml, df, beast.label.datepos= 4, beast.label.sep= '_', beast.date.direction= "forwards", beast.date.units= "years", verbose=1)
{			
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	#	check if any taxa to be added
	bxml.beast.existing.taxa.id		<- unlist( xpathApply(bxml.beast, "taxa/taxon", xmlGetAttr, "id" ) )
	tmp								<- subset( df, !BEASTlabel%in%bxml.beast.existing.taxa.id )
	#	add taxa if needed
	if(nrow(tmp))
	{		
		tmp.label	<- tmp[,BEASTlabel]	
		if(verbose)	cat(paste("\nsetting tip date to time at pos x in label, x=", beast.label.datepos))
		tmp.date	<- sapply( strsplit(tmp.label, beast.label.sep, fixed=1), function(x) x[beast.label.datepos] )		
		dummy		<- newXMLCommentNode(text="The list of taxa to be analysed (can also include dates/ages).", parent=bxml.beast, doc=bxml, addFinalizer=T)
		dummy		<- newXMLCommentNode(text=paste("ntax=",nrow(tmp),sep=''), parent=bxml.beast, doc=bxml, addFinalizer=T)	
		seqtaxa		<- newXMLNode("taxa", attrs= list(id="taxa"), parent=bxml.beast, doc=bxml, addFinalizer=T)
		dummy		<- lapply(seq_along(tmp.label), function(i)
				{
					taxon	<- newXMLNode("taxon", attrs= list(id=tmp.label[i]), parent=seqtaxa, doc=bxml, addFinalizer=T )
					dummy	<- newXMLNode("date", attrs= list(value=tmp.date[i], direction=beast.date.direction, units=beast.date.units), parent=taxon, doc=bxml, addFinalizer=T )
					taxon
				})	
		if(verbose)	cat(paste("\nadded new seq taxa, n=", xmlSize(seqtaxa)))		
	}
	invisible(bxml)
}	
######################################################################################
#	add taxa and alignment in bxml from BEASTlabels in df and alignment in seq.PROT.RT
#	beast.label.datepos= 4; beast.label.sep= '_'; beast.date.direction= "forwards"; beast.date.units= "years"; beast.alignment.dataType= "nucleotide"; xml.resetTipDate2LastDiag=1
#' @export
beast.add.alignment<- function(bxml, seq.PROT.RT, beast.alignment.id="alignment", beast.alignment.dataType= "nucleotide", verbose=1)
{			
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	df			<- data.table(FASTASampleCode=rownames(seq.PROT.RT), BEASTlabel=rownames(seq.PROT.RT))	
	#	add alignment	
	dummy		<- newXMLCommentNode(text="The sequence alignment (each sequence refers to a taxon above).", parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text=paste("ntax=",nrow(df)," nchar=",ncol(seq.PROT.RT),sep=''), parent=bxml.beast, doc=bxml, addFinalizer=T)	
	seqalign	<- newXMLNode("alignment", attrs= list(id=beast.alignment.id, dataType=beast.alignment.dataType), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- lapply( seq_len(nrow(df)), function(i)
			{
				seq		<- newXMLNode("sequence", parent=seqalign, doc=bxml, addFinalizer=T)
				dummy	<- newXMLNode("taxon", attrs= list(idref= df[i, BEASTlabel]), parent=seq, doc=bxml, addFinalizer=T)
				tmp		<- which( rownames(seq.PROT.RT)==df[i, FASTASampleCode] )
				if(length(tmp)!=1)	stop("unexpected row in seq.PROT.RT selected")
				tmp		<- paste(as.character(seq.PROT.RT[tmp,])[1,],collapse='',sep='')						
				dummy	<- newXMLTextNode(text=tmp, parent=seq, doc=bxml,addFinalizer=T)
				seq
			})	
	if(verbose)	cat(paste("\nadded new alignments, n=", xmlSize(seqalign)))	
	invisible(bxml)
}	
######################################################################################
#	add taxa and alignment in bxml from BEASTlabels in df and alignment in seq.PROT.RT
#	beast.label.datepos= 4; beast.label.sep= '_'; beast.date.direction= "forwards"; beast.date.units= "years"; beast.alignment.dataType= "nucleotide"; xml.resetTipDate2LastDiag=1
#' @export
beast.add.seq<- function(bxml, seq.PROT.RT, df=NULL, beast.label.datepos= 4, beast.label.sep= '_', beast.date.direction= "forwards", beast.date.units= "years", beast.alignment.id="alignment", beast.alignment.dataType= "nucleotide", verbose=1)
{			
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	if(is.null(df))
			df	<- data.table(FASTASampleCode=rownames(seq.PROT.RT), BEASTlabel=rownames(seq.PROT.RT))
	#	check if any taxa to be added
	bxml.beast.existing.taxa.id		<- unlist( xpathApply(bxml.beast, "taxa/taxon", xmlGetAttr, "id" ) )
	tmp								<- subset( df, !BEASTlabel%in%bxml.beast.existing.taxa.id )
	#	add taxa if needed
	if(nrow(tmp))
	{		
		tmp.label	<- tmp[,BEASTlabel]	
		if(verbose)	cat(paste("\nsetting tip date to time at pos x in label, x=", beast.label.datepos))
		tmp.date	<- sapply( strsplit(tmp.label, beast.label.sep, fixed=1), function(x) x[beast.label.datepos] )		
		dummy		<- newXMLCommentNode(text="The list of taxa to be analysed (can also include dates/ages).", parent=bxml.beast, doc=bxml, addFinalizer=T)
		dummy		<- newXMLCommentNode(text=paste("ntax=",nrow(tmp),sep=''), parent=bxml.beast, doc=bxml, addFinalizer=T)	
		seqtaxa		<- newXMLNode("taxa", attrs= list(id="taxa"), parent=bxml.beast, doc=bxml, addFinalizer=T)
		dummy		<- lapply(seq_along(tmp.label), function(i)
				{
					taxon	<- newXMLNode("taxon", attrs= list(id=tmp.label[i]), parent=seqtaxa, doc=bxml, addFinalizer=T )
					dummy	<- newXMLNode("date", attrs= list(value=tmp.date[i], direction=beast.date.direction, units=beast.date.units), parent=taxon, doc=bxml, addFinalizer=T )
					taxon
				})	
		if(verbose)	cat(paste("\nadded new seq taxa, n=", xmlSize(seqtaxa)))		
	}
	#	add alignment	
	dummy		<- newXMLCommentNode(text="The sequence alignment (each sequence refers to a taxon above).", parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text=paste("ntax=",nrow(df)," nchar=",ncol(seq.PROT.RT),sep=''), parent=bxml.beast, doc=bxml, addFinalizer=T)	
	seqalign	<- newXMLNode("alignment", attrs= list(id=beast.alignment.id, dataType=beast.alignment.dataType), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- lapply( seq_len(nrow(df)), function(i)
			{
				seq		<- newXMLNode("sequence", parent=seqalign, doc=bxml, addFinalizer=T)
				dummy	<- newXMLNode("taxon", attrs= list(idref= df[i, BEASTlabel]), parent=seq, doc=bxml, addFinalizer=T)
				tmp		<- which( rownames(seq.PROT.RT)==df[i, FASTASampleCode] )
				if(length(tmp)!=1)	stop("unexpected row in seq.PROT.RT selected")
				tmp		<- paste(as.character(seq.PROT.RT[tmp,])[1,],collapse='',sep='')						
				dummy	<- newXMLTextNode(text=tmp, parent=seq, doc=bxml,addFinalizer=T)
				seq
			})	
	if(verbose)	cat(paste("\nadded new alignments, n=", xmlSize(seqalign)))	
	invisible(bxml)
}	
######################################################################################
#	write nexus file for all sequences specified in df. assumes df has BEASTlabel. assumes seq.DNAbin.matrix and ph contain FASTASampleCode in df.
#' @export
beast.writeNexus4Beauti<- function( seq.DNAbin.matrix, df, ph=NULL, file=NULL )
{
	# 	select sequences and relabel					
	seq.DNAbin.matrix			<- seq.DNAbin.matrix[ df[,FASTASampleCode], ]	
	rownames(seq.DNAbin.matrix)	<- df[,BEASTlabel]
	#	generate starting tree and relabel
	if(!is.null(ph))
	{			
		tmp					<- setdiff( ph$tip.label, df[,FASTASampleCode] )
		tmp					<- match( tmp, ph$tip.label)
		ph.start			<- drop.tip(ph, tmp)
		ph.start$node.label	<- NULL
		setkey(df, FASTASampleCode)
		ph.start$tip.label	<- df[ph.start$tip.label,][,BEASTlabel]		
	}
	#	produce nexus text			
	ans<- seq.write.dna.nexus(seq.DNAbin.matrix, ph=ph.start, file=file)
	ans
}
######################################################################################
#	read tip stem samples after burn in from log file and return upper left points of a histogram of the Monte Carlo sample
#' @export
beast.read.log2tstem<- function(file.log, file.xml, beastlabel.idx.samplecode=6, burn.in= 5e6, breaks.n= 30, verbose=0)
{
	require(data.table)
	require(XML)	
	if(verbose)	cat(paste("\nReading file ",file.log))
	df.log	<- read.delim2(file.log, header = TRUE, sep = "\t", quote="\"", dec=".", fill = TRUE, comment.char="#")
	df.log	<- as.data.table(df.log)		
	if(verbose)	cat(paste("\nRead log file with ncol=", ncol(df.log)))
	tmp		<- c("state", colnames(df.log)[ grepl("tstem", colnames(df.log)) ] )
	if(length(tmp)==1)
	{
		if(verbose)	cat(paste("\nNo tstem found, return NA"))
		return( data.table( FASTASampleCode=NA, tstem= NA, density= NA ) )
	}
	df.log	<- subset( df.log, state>burn.in, select=tmp )
	if(verbose)	cat(paste("\nFound tstem data for n=", ncol(df.log)-1))
	#	translate tips back into FASTASampleCode
	if(verbose)	cat(paste("\nReading file ",file.xml))
	bxml				<- xmlTreeParse(file.xml, useInternalNodes=TRUE, addFinalizer = TRUE)
	tmp					<- regexpr("tip[0-9]+",colnames(df.log))
	if(!length(which(tmp>0)))
	{
		if(verbose)	cat(paste("\nNo tip tstem found, return NA"))
		return( data.table( FASTASampleCode=NA, tstem= NA, density= NA ) )
	}	
	log.tips			<- sapply(seq_along(tmp)[tmp>0], function(i)	substr(colnames(df.log)[i],tmp[i],tmp[i]+attr(tmp,"match.length")[i]-1)		)
	log.FASTASampleCode	<- sapply(log.tips, function(x) unlist( xpathApply(bxml, paste("//taxa[@id='",x,"']/taxon",sep=''), xmlGetAttr, "idref") )	)
	log.FASTASampleCode	<- sapply( strsplit(log.FASTASampleCode, '_'), function(x) x[beastlabel.idx.samplecode])		
	setnames(df.log, colnames(df.log)[tmp>0], log.FASTASampleCode)
	#	compute histograms for each tip stem
	ans		<- lapply(log.FASTASampleCode, function(x)
			{
				#x<- "R11-11357"
				tmp	<- hist( as.numeric(unlist( df.log[, x, with=0] )), breaks=breaks.n, plot=0 )
				data.table( FASTASampleCode=x, tstem= tmp$breaks, density= c(tmp$density,0) )										
			})
	if(verbose)	cat(paste("\nProcessed tstem data for tips n=", length(log.FASTASampleCode)))
	ans	<- rbindlist(ans)
}		
######################################################################################
#' @export
beast.add.taxonsets4clusters<- function(bxml, df, xml.monophyly4clusters=1, verbose=1)
{
	bxml.treeModel.id			<- unlist(xpathApply(bxml, "//treeModel[@id]", xmlGetAttr, "id"))
	
	#get taxon sets for each cluster
	btaxonsets.clusters			<- beast.get.taxonsets4clusters(bxml, df)	
	#get tmrcaStatistic for each cluster
	btmrcaStatistics.clusters	<- beast.get.tmrcaStatistic(bxml, btaxonsets.clusters, bxml.treeModel.id, includeStem="false") 		
	
	#	modify from template
	bxml.beast					<- getNodeSet(bxml, "//beast")[[1]]
	#	add 'btaxonsets.clusters' after last taxa
	bxml.idx					<- which(xmlSApply(bxml.beast, xmlName)=="taxa")		
	addChildren(bxml.beast, btaxonsets.clusters, at=bxml.idx[length(bxml.idx)] )
	if(verbose) cat(paste("\nadded taxon sets comprising each cluster, n=",length(btaxonsets.clusters)))
	#	add 'btmrcaStatistics.clusters' after last treeModel
	bxml.idx					<- which(xmlSApply(bxml.beast, xmlName)=="treeModel")
	addChildren(bxml.beast, btmrcaStatistics.clusters, at=bxml.idx[length(bxml.idx)] )
	if(verbose) cat(paste("\nadded tmrcaStatistics for mrca of each cluster, n=",length(btmrcaStatistics.clusters)))
	#add tmrcaStatistics to log 
	bxml.fileLog				<- getNodeSet(bxml, "//log[@id='fileLog']")
	tmrcaStatistics.id			<- sapply(btmrcaStatistics.clusters, function(x)	xmlGetAttr(x,"id")	)
	tmp							<- lapply(tmrcaStatistics.id, function(x)	newXMLNode("tmrcaStatistic", attrs= list(idref=x), parent=bxml.fileLog, doc=bxml, addFinalizer=T )	)
	if(verbose) cat(paste("\nadded tmrcaStatistics to fileLog"))
	#
	#	get monophylyStatistic and add after last 'tmrcaStatistic'
	#	populate booleanLikelihood with monophylyStatistic constraints
	#
	if(xml.monophyly4clusters)
	{
		bmStatistics.clusters	<- beast.get.monophylyStatistic(bxml, btaxonsets.clusters, bxml.treeModel.id)
		bxml.idx				<- which(xmlSApply(bxml.beast, xmlName)=="tmrcaStatistic")
		addChildren(bxml.beast, bmStatistics.clusters, at=bxml.idx[length(bxml.idx)] )
		if(verbose) cat(paste("\nadded monophylyStatistics for mrca of each cluster, n=",length(bmStatistics.clusters)))
		beast.add.monophylylkl(bxml, bmStatistics.clusters)		
		if(verbose) cat(paste("\nadded monophylyStatistics to booleanLikelihood, n=",length(bmStatistics.clusters)))
	}
	
	bxml
}

######################################################################################
#	get a list of monopylyStatistics. assumes each taxonset in 'btaxonsets' is monophyletic
#' @export
beast.get.monophylyStatistic<- function(bxml, btaxonsets, treeModel.id) 
{		
	btaxonset.id	<- sapply(btaxonsets, function(x)	xmlGetAttr(x,"id")	)
	ans				<- lapply(btaxonset.id, function(x)
			{
				monophylyStatistic	<- newXMLNode("monophylyStatistic", attrs= list(id=paste("monophyly(",x,")",sep='')), doc=bxml, addFinalizer=T )
				mrca			<- newXMLNode("mrca", parent=monophylyStatistic, doc=bxml, addFinalizer=T)
				newXMLNode("taxa", attrs= list(idref=x), parent=mrca, doc=bxml, addFinalizer=T )
				newXMLNode("treeModel", attrs= list(idref=treeModel.id), parent=monophylyStatistic, doc=bxml, addFinalizer=T )
				monophylyStatistic					
			})
	ans
}	
######################################################################################
#	add a list of monopylyStatistics to the BEAST prior. assumes all monophylyStatistics are referenced in the list 'monophylyStatistics'
#' @export
beast.add.monophylylkl<- function(bxml, monophylyStatistics)
{
	bxml.prior				<- getNodeSet(bxml, "//*[@id='prior']")
	if(length(bxml.prior)!=1)	stop("unexpected length of bxml.prior")
	bxml.prior				<- bxml.prior[[1]]
	#see if there is a booleanLikelihood, and if yes use it, otherwise add a new XML node
	bxml.bool.lkl			<- getNodeSet(bxml.prior,"//booleanLikelihood")
	if(length(bxml.bool.lkl)>1)	stop("unexpected length of bxml.bool.lkl")
	if(length(bxml.bool.lkl)<1)
		bxml.bool.lkl		<- newXMLNode("booleanLikelihood", parent=bxml.prior, doc=bxml, addFinalizer=T )
	else
		bxml.bool.lkl		<- bxml.bool.lkl[[1]]
	
	monophylyStatistics.id	<- sapply(monophylyStatistics, function(x)	xmlGetAttr(x,"id")	)
	dummy					<- lapply(monophylyStatistics.id, function(x)
			{
				dummy		<- newXMLNode("monophylyStatistic", attrs= list(idref=x), parent=bxml.bool.lkl, doc=bxml, addFinalizer=T )
			})
	bxml
}
######################################################################################
#	For each taxonset, get a list of tmrcaStatistics. Assumes all tmrcaStatistics share a treeModel.id and the same includeStem attribute
#' @export
beast.get.tmrcaStatistic<- function(bxml, btaxonsets, treeModel.id, includeStem="false") 
{		
	btaxonset.id	<- sapply(btaxonsets, function(x)	xmlGetAttr(x,"id")	)
	prefix.id		<- ifelse(includeStem=="false","tmrca","tstem")
	ans				<- lapply(btaxonset.id, function(x)
			{
				tmrcaStatistic	<- newXMLNode("tmrcaStatistic", attrs= list(id=paste(prefix.id,"(",x,")",sep=''), includeStem=includeStem), doc=bxml, addFinalizer=T )
				mrca			<- newXMLNode("mrca", parent=tmrcaStatistic, doc=bxml, addFinalizer=T)
				newXMLNode("taxa", attrs= list(idref=x), parent=mrca, doc=bxml, addFinalizer=T )
				newXMLNode("treeModel", attrs= list(idref=treeModel.id), parent=tmrcaStatistic, doc=bxml, addFinalizer=T )
				tmrcaStatistic					
			})
	ans
}	
######################################################################################
#' @export
beast.get.sequences<- function(bxml, verbose=1)
{			
	bxml.ali		<- getNodeSet(bxml, "//alignment[@id]")
	if(verbose)	cat(paste('\nFound alignments, n=',length(bxml.ali)))	
	ans				<- lapply(bxml.ali, function(ali)
			{
				bxml.seq		<- xpathApply(ali, "sequence", function(z){ 	xmlValue(z[['text']])		})
				bxml.seq		<- sapply(bxml.seq, function(z) gsub('\\s', '',z))	
				tmp				<- unlist( xpathApply(ali, "sequence/taxon", xmlGetAttr, "idref" ) )
				data.table(SEQ=bxml.seq, TAXON_ID=tmp, ALIGNMENT_ID=xmlGetAttr(ali,'id'))				
			})
	ans				<- do.call('rbind',ans)	
	if(verbose)	cat(paste('\nFound sequences, n=',nrow(ans)))
	ans
}
######################################################################################
#	For each tip: construct a prior for the corresponding tmrcaStatistics
#' @export
beast.get.tipPrior<- function(bxml, df, btmrcaStatistics.tips, xml.prior4stem="uniform", beast.label.negpos=2, beast.label.diagpos=3, beast.label.datepos=4, verbose=1)
{
	if(xml.prior4stem!="uniform")	stop("unexpected xml.tipprior")
	#
	if(verbose)	cat(paste("\nuse tip date found at pos x in label, x=", beast.label.datepos))
	df.height	<- t( sapply( strsplit(df[,BEASTlabel],'_',fixed=1), function(x)  as.numeric( x[c(beast.label.negpos, beast.label.diagpos, beast.label.datepos)]) ) )
	df.height	<- data.table(BEASTlabel=df[,BEASTlabel], NegT=df.height[,1], AnyPos_T1=df.height[,2], PosSeqT=df.height[,3])
	tmp			<- max( df.height[, PosSeqT])		#TODO should this be height or length ?
	df.height	<- df.height[, list(BEASTlabel=BEASTlabel, NegT=tmp-NegT, AnyPos_T1=tmp-AnyPos_T1, PosSeqT=tmp-PosSeqT)]
	#add uniform prior according to last NegT and first PosT
	ans			<- lapply(btmrcaStatistics.tips, function(x)
			{						
				tmrcaStatistics.id	<- xmlGetAttr(x,"id")
				tmp					<- xpathApply(x, "mrca/taxa", xmlGetAttr, "idref" )
				if(length(tmp)!=1)	stop("unexpected length of idref for mrca/taxa") 
				tip					<- as.numeric( substr(tmp[[1]],4,nchar(tmp[[1]])) )												
				bxml.tipprior		<- newXMLNode("uniformPrior", attrs= list(lower=df.height[tip,AnyPos_T1], upper=df.height[tip,NegT]), doc=bxml, addFinalizer=T )
				dummy				<- newXMLNode("statistic", attrs= list(idref=tmrcaStatistics.id), parent=bxml.tipprior, doc=bxml, addFinalizer=T )
				bxml.tipprior
			})
	ans			
}
######################################################################################
#	For each tip: add a taxonset, tmrcaStatistic, and potentially a reference to fileLog to 'bxml'
#	The aim of this as standalone is to log the length of tip stems. 
#' @export
beast.add.taxonsets4tips<- function(bxml, df, log=1, verbose=1)
{
	bxml.treeModel.id			<- unlist(xpathApply(bxml, "//treeModel[@id]", xmlGetAttr, "id"))
	
	#get taxon sets for each tip
	btaxonsets.tips				<- beast.get.taxonsets4tips(bxml, df)	
	#get tmrcaStatistic for each tip
	btmrcaStatistics.tips		<- beast.get.tmrcaStatistic(bxml, btaxonsets.tips, bxml.treeModel.id, includeStem="true") 			
	#	modify from template
	bxml.beast					<- getNodeSet(bxml, "//beast")[[1]]
	#	add 'btaxonsets.tips' after last taxa
	bxml.idx					<- which(xmlSApply(bxml.beast, xmlName)=="taxa")		
	addChildren(bxml.beast, btaxonsets.tips, at=bxml.idx[length(bxml.idx)] )
	if(verbose) cat(paste("\nadded taxon sets for each tip, n=",length(btaxonsets.tips)))
	#	add 'btmrcaStatistics.tips' after last treeModel or tmrcaStatistic
	bxml.idx					<- which(xmlSApply(bxml.beast, xmlName)%in%c("treeModel","tmrcaStatistic"))
	addChildren(bxml.beast, btmrcaStatistics.tips, at=tail(bxml.idx,1) )
	if(verbose) cat(paste("\nadded tmrcaStatistics for stem of each tip, n=",length(btmrcaStatistics.tips)))
	#add tmrcaStatistics to log 
	if(log)
	{
		bxml.fileLog				<- getNodeSet(bxml, "//log[@id='fileLog']")
		tmrcaStatistics.id			<- sapply(btmrcaStatistics.tips, function(x)	xmlGetAttr(x,"id")	)
		tmp							<- lapply(tmrcaStatistics.id, function(x)	newXMLNode("tmrcaStatistic", attrs= list(idref=x), parent=bxml.fileLog, doc=bxml, addFinalizer=T )	)
		if(verbose) cat(paste("\nadded tmrcaStatistics for tip stems to fileLog, n=",length(tmrcaStatistics.id)))
	}
	
	bxml
}
######################################################################################
#	For each tip: add a taxonset, tmrcaStatistic, prior for the tmrcaStatistics and potentially a reference to fileLog to 'bxml'
#' @export
beast.add.prior4tips<- function(bxml, df, xml.prior4stem="uniform", beast.label.datepos=4, verbose=1)
{
	#find list of tmrcaStatistics with id containing 'tip'
	btmrcaStatistics.tips		<- getNodeSet(bxml, "//tmrcaStatistic[starts-with(@id,'tstem(tip')]")
	print(btmrcaStatistics.tips)
	#get prior for each tip stem
	bprior.tips					<- beast.get.tipPrior(bxml, df, btmrcaStatistics.tips, xml.prior4stem=xml.prior4stem, beast.label.datepos=beast.label.datepos, verbose=verbose)		
	#	add 'bprior.tips' to prior
	bxml.prior					<- getNodeSet(bxml, "//prior[@id='prior']")
	if(length(bxml.prior)!=1)	stop("unexpected number of //prior[@id='prior']")
	bxml.prior					<- bxml.prior[[1]]
	addChildren(bxml.prior, bprior.tips)
	if(verbose) cat(paste("\nadded priors for stem of each tip, n=",length(bprior.tips)))	
	bxml
}
######################################################################################
#	create xml file from btemplate and seq.PROT.RT, using seq in df 
# 	beast.label.datepos= 4; beast.label.sep= '_'; beast.date.direction= "forwards"; beast.date.units= "years"; verbose=1; xml.prior4tipstem="uniform"; xml.resetTipDate2LastDiag=1
#' @export
beast.get.xml<- function(	btemplate, seq.PROT.RT, df, file, ph=NULL, xml.monophyly4clusters=0, xml.taxon4tipstem=0, xml.prior4tipstem=NA, 
		beast.label.datepos= 4, beast.label.sep= '_', beast.date.direction= "forwards", beast.usingDates="true", beast.date.units= "years", beast.mcmc.chainLength=50000000, verbose=1)
{
	
	bxml		<- newXMLDoc(addFinalizer=T)
	bxml.beast	<- newXMLNode("beast", doc=bxml, addFinalizer=T)
	newXMLCommentNode(text=paste("Generated by HIVCLUST from template",file), parent=bxml.beast, doc=bxml, addFinalizer=T)
	#	add new set of sequences
	dummy		<- beast.add.seq(bxml, seq.PROT.RT, df, beast.label.datepos=beast.label.datepos, beast.label.sep=beast.label.sep, beast.date.direction=beast.date.direction, beast.date.units=beast.date.units, verbose=verbose)
	#	copy everything after alignment up to but not including constantSize from template
	bt.beast	<- getNodeSet(btemplate, "//beast")[[1]]
	dummy		<- sapply(seq.int( which( xmlSApply(bt.beast, xmlName)=="alignment" )+1, which( xmlSApply(bt.beast, xmlName)=="patterns" ) ), function(i)
			{
				if( class(bt.beast[[i]])[1]=="XMLInternalCommentNode" )
					dummy<- newXMLCommentNode(text=xmlValue(bt.beast[[i]]), parent=bxml.beast, doc=bxml, addFinalizer=T)
				else
					dummy<- addChildren( bxml.beast, xmlClone( bt.beast[[i]], addFinalizer=T, doc=bxml ) )
			})
	#	add startingTree
	if(!is.null(ph))		#	add startingTree if provided
		dummy	<- beast.add.startingtree(bxml, ph, df, beast.usingDates=beast.usingDates)		
	else					#	otherwise copy upgmaTree
	{
		dummy		<- sapply(seq.int( which( xmlSApply(bt.beast, xmlName)=="constantSize" )-2, which( xmlSApply(bt.beast, xmlName)=="upgmaTree" ) ), function(i)
				{
					if( class(bt.beast[[i]])[1]=="XMLInternalCommentNode" )
						dummy<- newXMLCommentNode(text=xmlValue(bt.beast[[i]]), parent=bxml.beast, doc=bxml, addFinalizer=T)
					else
						dummy<- addChildren( bxml.beast, xmlClone( bt.beast[[i]], addFinalizer=T, doc=bxml ) )
				})
	}
	#	copy everything from 'treeModel'-1 from template
	dummy		<- sapply(seq.int( which( xmlSApply(bt.beast, xmlName)=="treeModel" )-1, xmlSize(bt.beast) ), function(i)
			{
				if( class(bt.beast[[i]])[1]=="XMLInternalCommentNode" )
					dummy<- newXMLCommentNode(text=xmlValue(bt.beast[[i]]), parent=bxml.beast, doc=bxml, addFinalizer=T)
				else
					dummy<- addChildren( bxml.beast, xmlClone( bt.beast[[i]], addFinalizer=T, doc=bxml ) )
			})
	#	if user-provided startingTree, reset <upgmaTree idref="startingTree"/> to <newick idref="startingTree"/>
	if(!is.null(ph))
	{
		tmp					<- getNodeSet(bxml, "//*[@idref='startingTree']")
		if(length(tmp)!=1) stop("unexpected number of //*[@idref='startingTree']")
		xmlName(tmp[[1]])	<- "newick"
	}
	# 	reset dimension of GMRF likelihood
	tmp			<- getNodeSet(bxml, "//*[@id='skyride.logPopSize']")
	if(length(tmp)!=1)	stop("unexpected number of *[@id='skyride.logPopSize'")
	tmp			<- tmp[[1]]
	xmlAttrs(tmp)["dimension"]	<-	nrow(df)-1  
	tmp			<- getNodeSet(bxml, "//*[@id='skyride.groupSize']")
	if(length(tmp)!=1)	stop("unexpected number of *[@id='skyride.groupSize'")
	tmp			<- tmp[[1]]
	xmlAttrs(tmp)["dimension"]	<-	nrow(df)-1
	#	if uniform prior for rootheight, set minimum to earliest sample time in data set
	dummy		<- beast.adjust.rootheightprior(bxml, df, verbose=verbose)
	dummy		<- beast.adjust.mcmc(bxml, beast.mcmc.chainLength=beast.mcmc.chainLength, verbose=verbose)
	#	for tips, add taxon sets and tmrcaStatistics
	if(xml.taxon4tipstem)
		dummy	<- beast.add.taxonsets4tips(bxml, df, log=1, verbose=1)
	if(!is.na(xml.prior4tipstem))
		dummy	<- beast.add.prior4tips(bxml, df, xml.prior4stem=xml.prior4tipstem, beast.label.datepos=beast.label.datepos, verbose=verbose)
	
	#	for clusters, add taxon sets and tmrcaStatistics 
	dummy		<- beast.add.taxonsets4clusters(bxml, df, xml.monophyly4clusters=xml.monophyly4clusters, verbose=verbose)		
	#	reset output fileNames
	bxml.onodes	<- getNodeSet(bxml, "//*[@fileName]")
	tmp			<- sapply(bxml.onodes, function(x) xmlGetAttr(x,"fileName"))
	tmp			<- gsub("(time).","time",tmp,fixed=1)
	tmp			<- gsub("(subst).","subst",tmp,fixed=1)
	tmp			<- sapply(strsplit(tmp,'.',fixed=1), function(x)	paste(file, '.', x[2], sep=''))		
	dummy		<- sapply(seq_along(bxml.onodes), function(i){		xmlAttrs(bxml.onodes[[i]])["fileName"]<- tmp[i]		})
	#
	bxml
}	
######################################################################################
#	adjust mcmc BEAST XML element
#' @export
beast.adjust.mcmc<- function(bxml, beast.mcmc.chainLength=NA, beast.mcmc.logEvery=NA, verbose=1)
{
	if(!is.na(beast.mcmc.chainLength))
	{
		tmp			<- getNodeSet(bxml, "//*[@id='mcmc']")
		if(length(tmp)!=1)	stop("unexpected number of *[@id='mcmc']")
		tmp			<- tmp[[1]]
		if(verbose)	cat(paste("\nSet MCMC chainLength to",beast.mcmc.chainLength))
		xmlAttrs(tmp)["chainLength"]	<-	sprintf("%d",beast.mcmc.chainLength)		
	}
	if(!is.na(beast.mcmc.logEvery))
	{
		bxml.logs	<- getNodeSet(bxml, "//*[@logEvery]")		
		for(i in seq_along(bxml.logs))
		{
			if(verbose)	cat(paste("\nSet logEvery to",beast.mcmc.logEvery))
			xmlAttrs(bxml.logs[[i]])["logEvery"]	<- sprintf("%d",beast.mcmc.logEvery)
		}			
	}
	bxml
}
######################################################################################
#	if rootheight prior uniform, sets lower bound to earliest sampling time in data set
#' @export
beast.adjust.rootheightprior<- function(bxml, df, verbose=1)
{
	bxml.prior	<- getNodeSet(bxml, "//uniformPrior[descendant::parameter[@idref='treeModel.rootHeight']]")
	if(length(bxml.prior)>1)	stop("unexpected length of //uniformPrior[descendant::parameter[@idref='treeModel.rootHeight']]")
	if(length(bxml.prior)==1)
	{
		tmp								<- df[,range(AnyPos_T1)]
		if(verbose)	cat(paste("\nfound uniformPrior for treeModel.rootHeight. Range of tip dates is",tmp[1],tmp[2]))
		tmp								<- difftime(tmp[2],tmp[1], units="days")	
		bxml.prior						<- bxml.prior[[1]]
		if(verbose)	cat(paste("\nset lower bound of uniformPrior for treeModel.rootHeight to", floor( tmp/365 )))
		xmlAttrs(bxml.prior)["lower"]	<- floor( tmp/365 )	
	}
	bxml
}
######################################################################################
#	For each cluster, create a taxonset. Assumes df has BEASTlabel and cluster
#' @export
beast.get.taxonsets4clusters	<- function(bxml, df)
{	
	ans	<- lapply( unique( df[,cluster] ), function(clu)
			{							
				taxonset	<- newXMLNode("taxa", attrs= list(id=paste("c",clu,sep='')), doc=bxml, addFinalizer=T )
				tmp<- df[cluster==clu,][,BEASTlabel]
				tmp<- lapply(tmp, function(x)		newXMLNode("taxon", attrs= list(idref=x), parent=taxonset, doc=bxml, addFinalizer=T )	)
				taxonset
			})
	ans		
}
######################################################################################
#	For each tip, create a taxonset. Assumes df has BEASTlabel
#' @export
beast.get.taxonsets4tips	<- function(bxml, df)
{	
	ans	<- lapply( seq_len(nrow(df)), function(i)
			{							
				taxonset	<- newXMLNode("taxa", attrs= list(id=paste("tip",i,sep='')), doc=bxml, addFinalizer=T )
				newXMLNode("taxon", attrs= list(idref=df[i,BEASTlabel]), parent=taxonset, doc=bxml, addFinalizer=T )
				taxonset
			})
	ans		
}
######################################################################################
#' @export
beast.add.variableDemographic<- function(bxml, demographic.id, treeModel.ids, coalescent.id, type='stepwise', useMidpoints='true', 
												popSize.id='demographic.popSize', popSize.value='0.091', 
												demographic.indicators.id='demographic.indicators', demographic.indicators.value='0.0',
												demographic.popMeanDist.id='demographic.populationMeanDist', demographic.popMean.id='demographic.populationMean', demographic.popMean.value='1.0', 
												sumStatistic.id='demographic.populationSizeChanges', sumStatistic.elementwise='true', 
												verbose=1)
{
	bxml.beast		<- getNodeSet(bxml, "//beast")[[1]]
	if(verbose)
		cat(paste('\nadd variableDemographic',demographic.id))
	bxml.dem	<- newXMLNode("variableDemographic", attrs= list(id=demographic.id, type=type, useMidpoints=useMidpoints), parent= bxml.beast, doc=bxml, addFinalizer=T)
	tmp			<- newXMLNode("populationSizes", parent=bxml.dem, doc=bxml, addFinalizer=T)
	dummy		<- newXMLNode("parameter", attrs= list(id=popSize.id, value=popSize.value ), parent=tmp, doc=bxml, addFinalizer=T)
	tmp			<- newXMLNode("indicators", parent=bxml.dem, doc=bxml, addFinalizer=T)
	dummy		<- newXMLNode("parameter", attrs= list(id=demographic.indicators.id, value=demographic.indicators.value ), parent=tmp, doc=bxml, addFinalizer=T)
	bxml.trees	<- newXMLNode("trees", parent=bxml.dem, doc=bxml, addFinalizer=T)
	for(x in treeModel.ids)
	{
		if(verbose)
			cat(paste('\nadd treeModel',x))
		tmp		<- newXMLNode("ptree", attrs= list(ploidy="2.0"), parent=bxml.trees, doc=bxml, addFinalizer=T)
		dummy	<- newXMLNode("treeModel", attrs= list(idref=x ), parent=tmp, doc=bxml, addFinalizer=T)	
	}
	#	coalescent
	bxml.coal	<- newXMLNode("coalescentLikelihood", attrs= list(id=coalescent.id), parent= bxml.beast, doc=bxml, addFinalizer=T)
	tmp			<- newXMLNode("model", parent=bxml.coal, doc=bxml, addFinalizer=T)
	dummy		<- newXMLNode("variableDemographic", attrs= list(idref=demographic.id), parent=tmp, doc=bxml, addFinalizer=T)
	#	sum statistic
	tmp			<- newXMLNode("sumStatistic", attrs= list(id=sumStatistic.id, elementwise=sumStatistic.elementwise), parent= bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLNode("parameter", attrs= list(idref=demographic.indicators.id), parent=tmp, doc=bxml, addFinalizer=T)
	#	demographic.populationMeanDist
	dummy		<- beast.add.exponentialDistributionModel(bxml, bxml.beast, demographic.popMeanDist.id, demographic.popMean.id, demographic.popMean.value)
	bxml 
}
######################################################################################
#' @export
beast.set.fileNameTrunk<- function(bxml, fileNameTrunk)
{
	bxml.onodes	<- getNodeSet(bxml, "//*[@fileName]")
	tmp			<- sapply(bxml.onodes, function(x) xmlGetAttr(x,"fileName"))
	tmp			<- gsub("(time).","time",tmp,fixed=1)
	tmp			<- gsub("(subst).","subst",tmp,fixed=1)	
	tmp			<- sapply(strsplit(tmp,'.',fixed=1), function(x)	paste(fileNameTrunk, '.', tail(x,1), sep=''))
	dummy		<- sapply(seq_along(bxml.onodes), function(i){		xmlAttrs(bxml.onodes[[i]])["fileName"]<- tmp[i]		})	
}
######################################################################################
#' @export
beast.add.VDAnalysis<- function(bxml, VDAnalysis.id, logFileName, treeFileNames, csvfileName, populationModelType='stepwise', populationFirstColumn='demographic.popSize1', indicatorsFirstColumn='demographic.indicators1', burnIn="0.1", useMidpoints="true", verbose=1)
{
	if(verbose)
		cat(paste('\nadd VDAnalysis',VDAnalysis.id))
	bxml.beast		<- getNodeSet(bxml, "//beast")[[1]]
	bxml.vda		<- newXMLNode("VDAnalysis", attrs= list(id=VDAnalysis.id, burnIn=burnIn, useMidpoints=useMidpoints), parent=bxml.beast, doc=bxml, addFinalizer=T)	
	tmp				<- newXMLNode("logFileName", parent=bxml.vda, doc=bxml, addFinalizer=T)	
	dummy			<- newXMLTextNode(text=logFileName, parent=tmp, doc=bxml,addFinalizer=T)
	tmp				<- newXMLNode("treeFileNames", parent=bxml.vda, doc=bxml, addFinalizer=T)
	for(x in treeFileNames)
	{
		tmp2		<- newXMLNode("treeOfLoci", parent=tmp, doc=bxml, addFinalizer=T)
		dummy		<- newXMLTextNode(text=x, parent=tmp2, doc=bxml, addFinalizer=T)
	}
	tmp				<- newXMLNode("populationModelType", parent=bxml.vda, doc=bxml, addFinalizer=T)
	dummy			<- newXMLTextNode(text=populationModelType, parent=tmp, doc=bxml, addFinalizer=T)	
	tmp				<- newXMLNode("populationFirstColumn", parent=bxml.vda, doc=bxml, addFinalizer=T)
	dummy			<- newXMLTextNode(text=populationFirstColumn, parent=tmp, doc=bxml, addFinalizer=T)	
	tmp				<- newXMLNode("indicatorsFirstColumn", parent=bxml.vda, doc=bxml, addFinalizer=T)
	dummy			<- newXMLTextNode(text=indicatorsFirstColumn, parent=tmp, doc=bxml, addFinalizer=T)
	#	CSVexport
	bxml.csv		<- newXMLNode("CSVexport", attrs= list(fileName=csvfileName, separator=","), parent=bxml.beast, doc=bxml, addFinalizer=T)
	tmp				<- newXMLNode("columns", parent=bxml.csv, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("VDAnalysis", attrs= list(idref=VDAnalysis.id), parent=tmp, doc=bxml, addFinalizer=T)
	bxml 
}
######################################################################################
#' @export
beast.add.logParameter<- function(bxml, logParameter.ids, verbose=1)
{
	tmp				<- getNodeSet(bxml, "//*[@id='fileLog']")[[1]]
	if(verbose)
		cat(paste('\nAdd log parameters, n=',length(logParameter.ids)))
	for(x in logParameter.ids)
		dummy		<- newXMLNode("parameter", attrs= list(idref=x), parent=tmp, doc=bxml, addFinalizer=T)
	bxml 
}
######################################################################################
#' @export
beast.add.logTree<- function(bxml, logTree.id, treeModel.id, discretizedBranchRates.id, posterior.id, logEvery=1e3, fileName='trees.trees', nexusFormat='true', sortTranslationTable='true', verbose=1)
{
	if(verbose)
		cat(paste('\nadd logTree',logTree.id))
	bxml.mcmc		<- getNodeSet(bxml, "//mcmc")[[1]]
	bxml.lt			<- newXMLNode("logTree", attrs= list(id=logTree.id, logEvery=sprintf("%d",logEvery), fileName=fileName, nexusFormat=nexusFormat, sortTranslationTable=sortTranslationTable), parent=bxml.mcmc, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("treeModel", attrs= list(idref=treeModel.id), parent=bxml.lt, doc=bxml, addFinalizer=T)
	tmp				<- newXMLNode("trait", attrs= list(name='rate', tag='rate'), parent=bxml.lt, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("discretizedBranchRates", attrs= list(idref=discretizedBranchRates.id), parent=tmp, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("posterior", attrs= list(idref=posterior.id), parent=bxml.lt, doc=bxml, addFinalizer=T)	
	invisible(bxml) 
}
######################################################################################
#' @export
beast.add.deltaExchangeOperator<- function(bxml, parameter.id, delta=0.75, parameterWeights=c(948, 948, 948), weight=9)
{
	bxml.o			<- getNodeSet(bxml, "//operators")[[1]]
	tmp				<- newXMLNode("deltaExchange", attrs= list(delta=as.character(delta), parameterWeights=paste(parameterWeights, collapse=' '), weight=as.character(weight)), parent=bxml.o, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(idref=parameter.id), parent=tmp, doc=bxml, addFinalizer=T)
	invisible(bxml) 
}
######################################################################################
#' @export
beast.add.gmrfGridBlockUpdateOperator<- function(bxml, parameter.id, scaleFactor=2, weight=2)
{
	bxml.o			<- getNodeSet(bxml, "//operators")[[1]]
	tmp				<- newXMLNode("gmrfGridBlockUpdateOperator", attrs= list(scaleFactor=as.character(scaleFactor), weight=as.character(weight)), parent=bxml.o, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("gmrfSkyrideLikelihood", attrs= list(idref=parameter.id), parent=tmp, doc=bxml, addFinalizer=T)
	invisible(bxml) 
}
######################################################################################
#' @export
beast.add.scaleOperator<- function(bxml, parameter.id, scaleFactor=0.75, weight=0.1)
{
	bxml.o			<- getNodeSet(bxml, "//operators")[[1]]
	tmp				<- newXMLNode("scaleOperator", attrs= list(scaleFactor=as.character(scaleFactor), weight=as.character(weight)), parent=bxml.o, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(idref=parameter.id), parent=tmp, doc=bxml, addFinalizer=T)
	invisible(bxml) 
}
######################################################################################
#' @export
beast.add.uniformIntegerOperator<- function(bxml, parameter.id, weight=0.1)
{
	bxml.beast		<- getNodeSet(bxml, "//operators")[[1]]
	tmp				<- newXMLNode("uniformIntegerOperator", attrs= list(weight=as.character(weight)), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(idref=parameter.id), parent=tmp, doc=bxml, addFinalizer=T)
	invisible(bxml) 
}
######################################################################################
#' @export
beast.add.swapOperator<- function(bxml, parameter.id, size=1, weight=0.1, autoOptimize='false')
{
	bxml.beast		<- getNodeSet(bxml, "//operators")[[1]]
	tmp				<- newXMLNode("swapOperator", attrs= list(size=as.character(size), weight=as.character(weight), autoOptimize=autoOptimize), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(idref=parameter.id), parent=tmp, doc=bxml, addFinalizer=T)
	invisible(bxml) 
}
######################################################################################
#' @export
beast.add.likelihood<- function(bxml, id, likelihood.term.ids)
{
	bxml.beast		<- getNodeSet(bxml, "//beast")[[1]]
	tmp				<- newXMLNode("likelihood", attrs= list(id=id), parent=bxml.beast, doc=bxml, addFinalizer=T)
	for(x in likelihood.term.ids)
		dummy		<- newXMLNode("treeLikelihood", attrs= list(idref=x), parent=tmp, doc=bxml, addFinalizer=T)	
	invisible(bxml) 
}
######################################################################################
#' @export
beast.add.priors<- function(bxml, prior.id)
{
	bxml.beast		<- getNodeSet(bxml, "//beast")[[1]]
	tmp				<- newXMLNode("prior", attrs= list(id=prior.id), parent=bxml.beast, doc=bxml, addFinalizer=T)	
	invisible(bxml)
}
######################################################################################
#' @export
beast.add.exponentialPrior<- function(bxml, idref, mean, offset)
{
	bxml.beast		<- getNodeSet(bxml, "//prior")[[1]]
	tmp				<- newXMLNode("exponentialPrior", attrs= list(mean=as.character(mean), offset=as.character(offset)), parent=bxml.beast, doc=bxml, addFinalizer=T)
	tmp2			<- newXMLNode("parameter", attrs= list(idref=idref), parent=tmp, doc=bxml, addFinalizer=T)
	invisible(bxml)
}
######################################################################################
#' @export
beast.add.gammaPrior<- function(bxml, idref, shape, scale, offset)
{
	bxml.beast		<- getNodeSet(bxml, "//prior")[[1]]
	tmp				<- newXMLNode("gammaPrior", attrs= list(shape=as.character(shape), scale=as.character(scale), offset=as.character(offset)), parent=bxml.beast, doc=bxml, addFinalizer=T)
	tmp2			<- newXMLNode("parameter", attrs= list(idref=idref), parent=tmp, doc=bxml, addFinalizer=T)
	invisible(bxml)
}
######################################################################################
#' @export
beast.add.uniformPrior<- function(bxml, idref, lower, upper)
{
	bxml.beast		<- getNodeSet(bxml, "//prior")[[1]]
	tmp				<- newXMLNode("uniformPrior", attrs= list(lower=as.character(lower), upper=as.character(upper)), parent=bxml.beast, doc=bxml, addFinalizer=T)
	tmp2			<- newXMLNode("parameter", attrs= list(idref=idref), parent=tmp, doc=bxml, addFinalizer=T)
	invisible(bxml)
}
######################################################################################
#' @export
beast.add.gmrfSkyGridPrior<- function(bxml, idref)
{
	bxml.beast		<- getNodeSet(bxml, "//prior")[[1]]
	tmp				<- newXMLNode("gmrfSkyGridLikelihood", attrs= list(idref=idref), parent=bxml.beast, doc=bxml, addFinalizer=T)	
	invisible(bxml)
}
######################################################################################
#' @export
beast.add.mcmc<- function(bxml, mcmc.id, posterior.id='posterior', prior.idref='prior', likelihood.idref='likelihood', mcmc.chainLength=2e6, mcmc.autoOptimize='true', mcmc.operatorAnalysis=NA)
{
	bxml.beast		<- getNodeSet(bxml, "//beast")[[1]]
	bxml.m			<- newXMLNode("mcmc", attrs= list(id=mcmc.id, chainLength=sprintf("%d",mcmc.chainLength), autoOptimize=mcmc.autoOptimize, operatorAnalysis=mcmc.operatorAnalysis), parent=bxml.beast, doc=bxml, addFinalizer=T)
	tmp				<- newXMLNode("posterior", attrs= list(id=posterior.id), parent=bxml.m, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("prior", attrs= list(idref=prior.idref), parent=tmp, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("likelihood", attrs= list(idref=likelihood.idref), parent=tmp, doc=bxml, addFinalizer=T)
	invisible(bxml)
}
######################################################################################
#' @export
beast.add.operators<- function(bxml, operators.id, operators.optimizationSchedule="log")
{
	bxml.beast		<- getNodeSet(bxml, "//mcmc")[[1]]
	tmp				<- newXMLNode("operators", attrs= list(id=operators.id, optimizationSchedule=operators.optimizationSchedule), parent=bxml.beast, doc=bxml, addFinalizer=T)
	invisible(bxml)
}
######################################################################################
#' @export
beast.add.screenLog<- function(bxml, log.id='screenLog', posterior.idref=NA, likelihood.idref=NA, prior.idref=NA, meanrate.idref=NA, logEvery=1e4, dp=4, sf=6, width=12)
{
	bxml.beast		<- getNodeSet(bxml, "//mcmc")[[1]]	
	bxml.log		<- newXMLNode("log", attrs= list(id=log.id, logEvery=sprintf("%d",logEvery)), parent=bxml.beast, doc=bxml, addFinalizer=T)
	if(!is.na(prior.idref))
	{
		tmp				<- newXMLNode("column", attrs= list(label=prior.idref, dp=as.character(dp), width=as.character(width)), parent=bxml.log, doc=bxml, addFinalizer=T)	
		dummy			<- newXMLNode("prior", attrs= list(idref=prior.idref), parent=tmp, doc=bxml, addFinalizer=T)		
	}
	if(!is.na(likelihood.idref))
	{		
		tmp				<- newXMLNode("column", attrs= list(label=likelihood.idref, dp=as.character(dp), width=as.character(width)), parent=bxml.log, doc=bxml, addFinalizer=T)	
		dummy			<- newXMLNode("likelihood", attrs= list(idref=likelihood.idref), parent=tmp, doc=bxml, addFinalizer=T)
	}
	if(!is.na(posterior.idref))
	{		
		tmp				<- newXMLNode("column", attrs= list(label=posterior.idref, dp=as.character(dp), width=as.character(width)), parent=bxml.log, doc=bxml, addFinalizer=T)	
		dummy			<- newXMLNode("posterior", attrs= list(idref=posterior.idref), parent=tmp, doc=bxml, addFinalizer=T)
	}
	if(!is.na(meanrate.idref))
	{
		tmp				<- newXMLNode("column", attrs= list(label=meanrate.idref, sf=as.character(sf), width=as.character(width)), parent=bxml.log, doc=bxml, addFinalizer=T)	
		dummy			<- newXMLNode("parameter", attrs= list(idref=meanrate.idref), parent=tmp, doc=bxml, addFinalizer=T)
	}
	invisible(bxml) 
}
######################################################################################
#' @export
beast.add.fileLog<- function(bxml, log.id='fileLog', select.id=c('rootHeight', 'ucld', 'hky', 'skygrid','branchRates'), deselect.id=c('frequencies','categories'), posterior.idref=NA, likelihood.idref=NA, prior.idref=NA, logEvery=1e4, log.fileName=NA, log.overwrite='false')
{
	bxml.beast		<- getNodeSet(bxml, "//mcmc")[[1]]	
	bxml.log		<- newXMLNode("log", attrs= list(id=log.id, logEvery=sprintf("%d",logEvery), fileName=log.fileName, overwrite=log.overwrite), parent=bxml.beast, doc=bxml, addFinalizer=T)
	if(!is.na(prior.idref))
	{		
		dummy		<- newXMLNode("prior", attrs= list(idref=prior.idref), parent=bxml.log, doc=bxml, addFinalizer=T)		
	}
	if(!is.na(likelihood.idref))
	{					
		dummy		<- newXMLNode("likelihood", attrs= list(idref=likelihood.idref), parent=bxml.log, doc=bxml, addFinalizer=T)
	}
	if(!is.na(posterior.idref))
	{					
		dummy		<- newXMLNode("posterior", attrs= list(idref=posterior.idref), parent=bxml.log, doc=bxml, addFinalizer=T)
	}
	#	add selected parameter ids
	tmp				<- unique(unlist(xpathApply(bxml, "//parameter[@id]", xmlGetAttr, "id")))
	tmp				<- tmp[ as.logical(apply(sapply(select.id, function(x) as.numeric(grepl(x, tmp)) ), 1, sum)) ]
	tmp				<- tmp[ as.logical(apply(sapply(deselect.id, function(x) as.numeric(!grepl(x, tmp)) ), 1, prod)) ]
	for(x in tmp)
		dummy		<- newXMLNode("parameter", attrs= list(idref=x), parent=bxml.log, doc=bxml, addFinalizer=T)
	#	add selected rateStatistic ids
	tmp				<- unique(unlist(xpathApply(bxml, "//rateStatistic[@id]", xmlGetAttr, "id")))
	if(!is.null(tmp))
	{
		tmp				<- tmp[ as.logical(apply(sapply(select.id, function(x) as.numeric(grepl(x, tmp)) ), 1, sum)) ]
		tmp				<- tmp[ as.logical(apply(sapply(deselect.id, function(x) as.numeric(!grepl(x, tmp)) ), 1, prod)) ]
		for(x in tmp)
			dummy		<- newXMLNode("rateStatistic", attrs= list(idref=x), parent=bxml.log, doc=bxml, addFinalizer=T)			
	}
	#	add selected rateCovarianceStatistic ids
	tmp				<- unique(unlist(xpathApply(bxml, "//rateCovarianceStatistic[@id]", xmlGetAttr, "id")))
	if(!is.null(tmp))
	{
		tmp				<- tmp[ as.logical(apply(sapply(select.id, function(x) as.numeric(grepl(x, tmp)) ), 1, sum)) ]
		tmp				<- tmp[ as.logical(apply(sapply(deselect.id, function(x) as.numeric(!grepl(x, tmp)) ), 1, prod)) ]
		for(x in tmp)
			dummy		<- newXMLNode("rateCovarianceStatistic", attrs= list(idref=x), parent=bxml.log, doc=bxml, addFinalizer=T)
		
	}	
	invisible(bxml) 
}
######################################################################################
#' @export
beast.add.rateStatistics<- function(bxml, prefix.id, treeModel.id, discretizedBranchRates.id)
{
	bxml.beast		<- getNodeSet(bxml, "//beast")[[1]]
	tmp				<- newXMLNode("rateStatistic", attrs= list(id=paste(prefix.id,'_meanRate',sep=''), name=paste(prefix.id,'_meanRate',sep=''), mode='mean', internal='true', external='true'), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("treeModel", attrs= list(idref=treeModel.id), parent=tmp, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("discretizedBranchRates", attrs= list(idref=discretizedBranchRates.id), parent=tmp, doc=bxml, addFinalizer=T)	
	tmp				<- newXMLNode("rateStatistic", attrs= list(id=paste(prefix.id,'_coefficientOfVariation',sep=''), name=paste(prefix.id,'_coefficientOfVariation',sep=''), mode='coefficientOfVariation', internal='true', external='true'), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("treeModel", attrs= list(idref=treeModel.id), parent=tmp, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("discretizedBranchRates", attrs= list(idref=discretizedBranchRates.id), parent=tmp, doc=bxml, addFinalizer=T)	
	tmp				<- newXMLNode("rateStatistic", attrs= list(id=paste(prefix.id,'_covariance',sep=''), name=paste(prefix.id,'_covariance',sep=''), mode='variance', internal='true', external='true' ), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("treeModel", attrs= list(idref=treeModel.id), parent=tmp, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("discretizedBranchRates", attrs= list(idref=discretizedBranchRates.id), parent=tmp, doc=bxml, addFinalizer=T)	
	invisible(bxml) 
}
######################################################################################
#' @export
beast.add.siteModel<- function(bxml, siteModel.id, substitutionModel.id, substitutionModel.class='gtrModel', 
										relativeRate.id=NA, relativeRate.idref=NA, relativeRate.value=NA, relativeRate.lower=NA,
										gamma.id=NA, gamma.idref=NA, gamma.value=NA, gamma.lower=NA,
										gammaCategories=NA, verbose=1)
{
	bxml.beast		<- getNodeSet(bxml, "//beast")[[1]]
	if(verbose)
		cat(paste('\nadd siteModel',siteModel.id))		
	bxml.site		<- newXMLNode("siteModel", attrs= list(id=siteModel.id), parent=bxml.beast, doc=bxml, addFinalizer=T)
	tmp				<- newXMLNode("substitutionModel", parent=bxml.site, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode(substitutionModel.class, attrs= list(idref=substitutionModel.id), parent=tmp, doc=bxml, addFinalizer=T)
	if(!is.na(relativeRate.idref) | !is.na(relativeRate.id))
	{
		tmp				<- newXMLNode("relativeRate", parent=bxml.site, doc=bxml, addFinalizer=T)
		if(!is.na(relativeRate.id))
			dummy		<- newXMLNode("parameter", attrs= list(id=relativeRate.id, value=relativeRate.value, lower=relativeRate.lower), parent=tmp, doc=bxml, addFinalizer=T)
		if(!is.na(relativeRate.idref))
			dummy		<- newXMLNode("parameter", attrs= list(idref=relativeRate.idref), parent=tmp, doc=bxml, addFinalizer=T)		
	}
	tmp				<- newXMLNode("gammaShape", attrs= list(gammaCategories=gammaCategories), parent=bxml.site, doc=bxml, addFinalizer=T)
	if(is.na(gamma.idref))
		dummy		<- newXMLNode("parameter", attrs= list(id=gamma.id, value=gamma.value, lower=gamma.lower), parent=tmp, doc=bxml, addFinalizer=T)
	if(!is.na(gamma.idref))
		dummy		<- newXMLNode("parameter", attrs= list(idref=gamma.idref), parent=tmp, doc=bxml, addFinalizer=T)
	invisible(bxml) 
}
######################################################################################
#' @export
beast.add.hkyModel<- function(bxml, hkyModel.id, patterns.id, frequencies.id,
		rate.id=NA, rate.idref=NA, rate.value=NA, rate.lower=NA, frequencies.dimension='4', frequencyModel.dataType='nucleotide', verbose=1)
{
	bxml.beast		<- getNodeSet(bxml, "//beast")[[1]]
	if(verbose)
		cat(paste('\nadd hkyModel',hkyModel.id))	
	bxml.hky		<- newXMLNode("HKYModel", attrs= list(id=hkyModel.id), parent=bxml.beast, doc=bxml, addFinalizer=T)
	tmp				<- newXMLNode("frequencies", parent=bxml.hky, doc=bxml, addFinalizer=T)
	tmp2			<- newXMLNode("frequencyModel", attrs= list(dataType=frequencyModel.dataType), parent=tmp, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("patterns", attrs= list(idref=patterns.id), parent=tmp2, doc=bxml, addFinalizer=T)
	tmp3			<- newXMLNode("frequencies", parent=tmp2, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(id=frequencies.id, dimension=frequencies.dimension), parent=tmp3, doc=bxml, addFinalizer=T)
	
	tmp				<- newXMLNode("kappa", parent=bxml.hky, doc=bxml, addFinalizer=T)
	if(is.na(rate.idref))
		dummy		<- newXMLNode("parameter", attrs= list(id=rate.id, value=rate.value, lower=rate.lower), parent=tmp, doc=bxml, addFinalizer=T)
	if(!is.na(rate.idref))
		dummy		<- newXMLNode("parameter", attrs= list(idref=rate.idref), parent=tmp, doc=bxml, addFinalizer=T)	
	invisible(bxml) 
}
######################################################################################
#' @export
beast.add.gtrModel<- function(bxml, gtrModel.id, patterns.id, frequencies.id,
									rateAC.id=NA, rateAC.idref=NA, rateAC.value=NA, rateAC.lower=NA,
									rateAG.id=NA, rateAG.idref=NA, rateAG.value=NA, rateAG.lower=NA,
									rateAT.id=NA, rateAT.idref=NA, rateAT.value=NA, rateAT.lower=NA,
									rateCG.id=NA, rateCG.idref=NA, rateCG.value=NA, rateCG.lower=NA,
									rateGT.id=NA, rateGT.idref=NA, rateGT.value=NA, rateGT.lower=NA,
									frequencies.dimension='4', frequencyModel.dataType='nucleotide', verbose=1)
{
	bxml.beast		<- getNodeSet(bxml, "//beast")[[1]]
	if(verbose)
		cat(paste('\nadd gtrModel',gtrModel.id))	
	bxml.gtr		<- newXMLNode("gtrModel", attrs= list(id=gtrModel.id), parent=bxml.beast, doc=bxml, addFinalizer=T)
	tmp				<- newXMLNode("frequencies", parent=bxml.gtr, doc=bxml, addFinalizer=T)
	tmp2			<- newXMLNode("frequencyModel", attrs= list(dataType=frequencyModel.dataType), parent=tmp, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("patterns", attrs= list(idref=patterns.id), parent=tmp2, doc=bxml, addFinalizer=T)
	tmp3			<- newXMLNode("frequencies", parent=tmp2, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(id=frequencies.id, dimension=frequencies.dimension), parent=tmp3, doc=bxml, addFinalizer=T)
	
	tmp				<- newXMLNode("rateAC", parent=bxml.gtr, doc=bxml, addFinalizer=T)
	if(is.na(rateAC.idref))
		dummy		<- newXMLNode("parameter", attrs= list(id=rateAC.id, value=rateAC.value, lower=rateAC.lower), parent=tmp, doc=bxml, addFinalizer=T)
	if(!is.na(rateAC.idref))
		dummy		<- newXMLNode("parameter", attrs= list(idref=rateAC.idref), parent=tmp, doc=bxml, addFinalizer=T)	
	tmp				<- newXMLNode("rateAG", parent=bxml.gtr, doc=bxml, addFinalizer=T)
	if(is.na(rateAG.idref))
		dummy		<- newXMLNode("parameter", attrs= list(id=rateAG.id, value=rateAG.value, lower=rateAG.lower), parent=tmp, doc=bxml, addFinalizer=T)
	if(!is.na(rateAG.idref))
		dummy		<- newXMLNode("parameter", attrs= list(idref=rateAG.idref), parent=tmp, doc=bxml, addFinalizer=T)
	tmp				<- newXMLNode("rateAT", parent=bxml.gtr, doc=bxml, addFinalizer=T)
	if(is.na(rateAT.idref))
		dummy		<- newXMLNode("parameter", attrs= list(id=rateAT.id, value=rateAT.value, lower=rateAT.lower), parent=tmp, doc=bxml, addFinalizer=T)
	if(!is.na(rateAT.idref))
		dummy		<- newXMLNode("parameter", attrs= list(idref=rateAT.idref), parent=tmp, doc=bxml, addFinalizer=T)
	tmp				<- newXMLNode("rateCG", parent=bxml.gtr, doc=bxml, addFinalizer=T)
	if(is.na(rateCG.idref))
		dummy		<- newXMLNode("parameter", attrs= list(id=rateCG.id, value=rateCG.value, lower=rateCG.lower), parent=tmp, doc=bxml, addFinalizer=T)
	if(!is.na(rateCG.idref))
		dummy		<- newXMLNode("parameter", attrs= list(idref=rateCG.idref), parent=tmp, doc=bxml, addFinalizer=T)
	tmp				<- newXMLNode("rateGT", parent=bxml.gtr, doc=bxml, addFinalizer=T)
	if(is.na(rateGT.idref))
		dummy		<- newXMLNode("parameter", attrs= list(id=rateGT.id, value=rateGT.value, lower=rateGT.lower), parent=tmp, doc=bxml, addFinalizer=T)
	if(!is.na(rateGT.idref))
		dummy		<- newXMLNode("parameter", attrs= list(idref=rateGT.idref), parent=tmp, doc=bxml, addFinalizer=T)	
	invisible(bxml) 
}
######################################################################################
#' @export
beast.add.exponentialDistributionModel<- function(bxml, parent, id, mean.id, mean.value)
{
	bxml.beast		<- getNodeSet(bxml, "//beast")[[1]]
	bxml.ln			<- newXMLNode("exponentialDistributionModel", attrs= list(id=id),parent=parent, doc=bxml, addFinalizer=T)
	tmp				<- newXMLNode("mean", parent=bxml.ln, doc=bxml, addFinalizer=T)	
	dummy			<- newXMLNode("parameter", attrs= list(id=mean.id, value=mean.value), parent=tmp, doc=bxml, addFinalizer=T)	
	bxml 
}
######################################################################################
#' @export
beast.add.logNormalDistributionModel<- function(bxml, parent, mean.id, mean.idref, mean.value, mean.lower, sd.id, sd.idref, sd.value, sd.lower, meanInRealSpace="true")
{
	bxml.beast		<- getNodeSet(bxml, "//beast")[[1]]
	bxml.ln			<- newXMLNode("logNormalDistributionModel", attrs= list(meanInRealSpace=meanInRealSpace), parent=parent, doc=bxml, addFinalizer=T)
	tmp				<- newXMLNode("mean", parent=bxml.ln, doc=bxml, addFinalizer=T)
	if(is.na(mean.idref))
		dummy		<- newXMLNode("parameter", attrs= list(id=mean.id, value=mean.value, lower=mean.lower ), parent=tmp, doc=bxml, addFinalizer=T)	
	if(!is.na(mean.idref))
		dummy		<- newXMLNode("parameter", attrs= list(idref=mean.idref ), parent=tmp, doc=bxml, addFinalizer=T)		
	tmp				<- newXMLNode("stdev", parent=bxml.ln, doc=bxml, addFinalizer=T)
	if(is.na(sd.idref))
		dummy		<- newXMLNode("parameter", attrs= list(id=sd.id, value=sd.value, lower=sd.lower ), parent=tmp, doc=bxml, addFinalizer=T)	
	if(!is.na(sd.idref))
		dummy		<- newXMLNode("parameter", attrs= list(idref=sd.idref ), parent=tmp, doc=bxml, addFinalizer=T)		
	invisible(bxml) 
}
######################################################################################
#' @export
beast.add.gmrfskygrid<- function(bxml, gmrfSkyGridLikelihood.id, populationSizes.id, precisionParameter.id, numGridPoints.id, cutOff.id, treeModel.ids,  
										populationSizes.dimension, populationSizes.value,
										precisionParameter.value, precisionParameter.lower,
										cutOff.value, numGridPoints.value=populationSizes.dimension-1, verbose=1)
{
	bxml.beast		<- getNodeSet(bxml, "//beast")[[1]]
	if(verbose)
		cat(paste('\nadd gmrfSkyGridLikelihood',gmrfSkyGridLikelihood.id))
	bxml.grid		<- newXMLNode("gmrfSkyGridLikelihood", attrs= list(id=gmrfSkyGridLikelihood.id), parent= bxml.beast, doc=bxml, addFinalizer=T)
	tmp				<- newXMLNode("populationSizes", parent=bxml.grid, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(id=populationSizes.id, dimension=as.character(populationSizes.dimension), value=as.character(populationSizes.value) ), parent=tmp, doc=bxml, addFinalizer=T)
	tmp				<- newXMLNode("precisionParameter", parent=bxml.grid, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(id=precisionParameter.id, value=as.character(precisionParameter.value), lower=as.character(precisionParameter.lower) ), parent=tmp, doc=bxml, addFinalizer=T)
	tmp				<- newXMLNode("numGridPoints", parent=bxml.grid, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(id=numGridPoints.id, value=as.character(numGridPoints.value) ), parent=tmp, doc=bxml, addFinalizer=T)
	tmp				<- newXMLNode("cutOff", parent=bxml.grid, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(id=cutOff.id, value=as.character(cutOff.value) ), parent=tmp, doc=bxml, addFinalizer=T)
	tmp				<- newXMLNode("populationTree", parent=bxml.grid, doc=bxml, addFinalizer=T)
	for(x in treeModel.ids)
	{
		if(verbose)
			cat(paste('\nadd treeModel to GMRF skygrid',x))		
		dummy	<- newXMLNode("treeModel", attrs= list(idref=x ), parent=tmp, doc=bxml, addFinalizer=T)	
	}
	invisible(bxml)	 
}
######################################################################################
#' @export
beast.add.discretizedBranchRates<- function(bxml, discretizedBranchRates.id, treeModel.id, rateCategories.id,
													mean.id=NA, mean.idref=NA, mean.value=NA, mean.lower=NA, sd.id=NA, sd.idref=NA, sd.value=NA, sd.lower=NA, meanInRealSpace='true', rateCategories.dimension="NA", verbose=1)
{
	bxml.beast		<- getNodeSet(bxml, "//beast")[[1]]
	if(verbose)
		cat(paste('\nadd discretizedBranchRates',discretizedBranchRates.id))
	bxml.brr		<- newXMLNode("discretizedBranchRates", attrs= list(id=discretizedBranchRates.id), parent= bxml.beast, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("treeModel", attrs= list(idref=treeModel.id ), parent=bxml.brr, doc=bxml, addFinalizer=T)
	tmp				<- newXMLNode("distribution", parent=bxml.brr, doc=bxml, addFinalizer=T)
	dummy			<- beast.add.logNormalDistributionModel(bxml, tmp, mean.id, mean.idref, mean.value, mean.lower, sd.id, sd.idref, sd.value, sd.lower, meanInRealSpace=meanInRealSpace)
	tmp				<- newXMLNode("rateCategories", parent=bxml.brr, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(id=rateCategories.id, dimension=as.character(rateCategories.dimension) ), parent=tmp, doc=bxml, addFinalizer=T)
	invisible(bxml)	 
}
######################################################################################
#' @export
beast.add.treeLikelihood<- function(bxml, treeLikelihood.id, patterns.id, treeModel.id, siteModel.id, discretizedBranchRates.id, useAmbiguities='false', verbose=1)
{
	bxml.beast		<- getNodeSet(bxml, "//beast")[[1]]
	if(verbose)
		cat(paste('\nadd treeLikelihood',treeLikelihood.id))
	bxml.treelkl	<- newXMLNode("treeLikelihood", attrs= list(id=treeLikelihood.id, useAmbiguities=useAmbiguities), parent= bxml.beast, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("patterns", attrs= list(idref=patterns.id ), parent=bxml.treelkl, doc=bxml, addFinalizer=T)	
	dummy			<- newXMLNode("treeModel", attrs= list(idref=treeModel.id ), parent=bxml.treelkl, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("siteModel", attrs= list(idref=siteModel.id ), parent=bxml.treelkl, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("discretizedBranchRates", attrs= list(idref=discretizedBranchRates.id ), parent=bxml.treelkl, doc=bxml, addFinalizer=T)	
	invisible(bxml) 
}
######################################################################################
#' @export
beast.add.treemodel<- function(bxml, treemodel.id, rootHeight.id, internalNodeHeights.id, allInternalNodeHeights.id, newick.id=NA, internalNodes='true', rootNode='true', verbose=1)
{
	bxml.beast		<- getNodeSet(bxml, "//beast")[[1]]
	if(verbose)
		cat(paste('\nadd treeModel',treemodel.id))
	bxml.treemodel	<- newXMLNode("treeModel", attrs= list(id=treemodel.id), parent= bxml.beast, doc=bxml, addFinalizer=T)
	if(!is.na(newick.id))
		dummy		<- newXMLNode("newick", attrs= list(idref=newick.id ), parent=bxml.treemodel, doc=bxml, addFinalizer=T)
	rootHeight	<- newXMLNode("rootHeight", parent=bxml.treemodel, doc=bxml, addFinalizer=T)
	dummy		<- newXMLNode("parameter", attrs= list(id=rootHeight.id ), parent=rootHeight, doc=bxml, addFinalizer=T)	
	nodeHeights	<- newXMLNode("nodeHeights", attrs= list(internalNodes=internalNodes ),parent=bxml.treemodel, doc=bxml, addFinalizer=T)
	dummy		<- newXMLNode("parameter", attrs= list(id=internalNodeHeights.id ), parent=nodeHeights, doc=bxml, addFinalizer=T)	
	nodeHeights2<- newXMLNode("nodeHeights", attrs= list(internalNodes=internalNodes, rootNode=rootNode),parent=bxml.treemodel, doc=bxml, addFinalizer=T)
	dummy		<- newXMLNode("parameter", attrs= list(id=allInternalNodeHeights.id ), parent=nodeHeights2, doc=bxml, addFinalizer=T)	
	invisible(bxml) 
}
######################################################################################
#' @export
beast.add.compoundParameter<- function(bxml, compoundParameter.id, select.id='site.mu', verbose=1)
{
	bxml.beast		<- getNodeSet(bxml, "//beast")[[1]]
	tmp				<- unique(unlist(xpathApply(bxml, "//parameter[@id]", xmlGetAttr, "id")))
	tmp				<- tmp[ grepl(select.id, tmp) ]
	if(length(tmp))
	{
		bxml.p		<- newXMLNode("compoundParameter", attrs= list(id=compoundParameter.id), parent= bxml.beast, doc=bxml, addFinalizer=T)
		for(x in tmp)
		{
			if(verbose)
				cat(paste('\nadd parameter',x,'to compound parameter',compoundParameter.id))			
			dummy	<- newXMLNode("parameter", attrs= list(idref=x), parent= bxml.p, doc=bxml, addFinalizer=T)
		}
	}
	invisible(bxml) 
}
######################################################################################
#' @export
beast.add.patterns<- function(bxml, beast.patterns.id, alignment.id, beast.patterns.from, beast.patterns.every=3, beast.patterns.strip='false', verbose=1)
{
	bxml.beast		<- getNodeSet(bxml, "//beast")[[1]]
	if(verbose)
		cat(paste('\nadd pattern',beast.patterns.id,'for alignment',alignment.id))
	if(is.na(beast.patterns.every))
		bxml.pattern	<- newXMLNode("patterns", attrs= list(id=beast.patterns.id, from=as.character(beast.patterns.from), strip=beast.patterns.strip), parent= bxml.beast, doc=bxml, addFinalizer=T)	
	if(!is.na(beast.patterns.every))
		bxml.pattern	<- newXMLNode("patterns", attrs= list(id=beast.patterns.id, from=as.character(beast.patterns.from), every=as.character(beast.patterns.every), strip=beast.patterns.strip), parent= bxml.beast, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("alignment", attrs= list(idref=alignment.id ), parent=bxml.pattern, doc=bxml, addFinalizer=T)
	invisible(bxml) 
}
######################################################################################
# 	extract starting tree from 'ph' by tips in 'df'. Only keeps tree topology and resets branch lengths so that the maximum root distance is 'beast.rootHeight'
#	beast.rootHeight= 35; beast.usingDates= "false"; beast.newickid= "startingTree"
#' @export
beast.add.startingtree<- function(bxml, ph.start, df=NULL, beast.rootHeight= NA, beast.usingDates="true", beast.newickid= "startingTree", beast.brlunits="years", verbose=1)
{
	require(adephylo)
	if(verbose) cat(paste("\ncreate startingTree"))
	if(!is.null(df))
	{
		tmp					<- setdiff( ph.start$tip.label, df[,FASTASampleCode] )
		tmp					<- match( tmp, ph.start$tip.label)
		ph.start			<- drop.tip(ph.start, tmp)		
		ph.start$node.label	<- NULL
		setkey(df, FASTASampleCode)
		ph.start$tip.label	<- df[ph.start$tip.label,][,BEASTlabel]
		if(verbose) cat(paste("\nselected tips for startingTree, n=",Ntip(ph.start)))		
	}
	if(is.null(df) & verbose)
		cat(paste("\nuse tips for startingTree, n=",Ntip(ph.start)))
	#	adjust rootHeight to 'beast.rootHeight'
	if(!is.na(beast.rootHeight))
	{
		if(verbose) cat(paste("\nSet root height of starting tree=",beast.rootHeight))
		tmp					<- beast.rootHeight / max(distRoot(ph.start))
		ph.start$edge.length<- ph.start$edge.length*tmp		
	}
	#	write ph.start as newick tree to bxml
	tmp					<- write.tree( ph.start )	
	bxml.beast			<- getNodeSet(bxml, "//beast")[[1]]
	dummy				<- newXMLCommentNode(text="The user-specified starting tree in a newick tree format", parent=bxml.beast, doc=bxml, addFinalizer=T)
	bxml.startingTree	<- newXMLNode("newick", attrs= list(id=beast.newickid, usingDates=beast.usingDates, units=beast.brlunits), parent= bxml.beast, doc=bxml, addFinalizer=T)
	dummy				<- newXMLTextNode(text=tmp, parent=bxml.startingTree, doc=bxml, addFinalizer=T) 
	invisible(bxml)
}