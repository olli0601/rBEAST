
######################################################################################
#' @title Set TAXON_NAME field from columns in the taxon data.table.
#' @return Taxon data.table with additional column TAXON_NAME.
#' @export
beast.set.TAXON_NAME<- function( df, select=c('TAXON_ID'), select.sep='_')
{	
	stopifnot('TAXON_ID'%in%names(df), select%in%names(df))
	tmp		<- df[, list(TAXON_NAME=paste(.SD ,sep='', collapse=select.sep)), by='TAXON_ID', .SDcols=select]
	merge(df, tmp, by="TAXON_ID")	
}
######################################################################################
#' @title Set SEQ field from sequences in sequences in \code{DNAbin} matrix format
#' @description Sequences whose name matches the TAXON_ID field are added to the taxon data.table. 
#' @return Taxon data.table with additional column SEQ that contains sequences as a character string.
#' @export
beast.set.SEQ	<- function(df, seq, verbose=1)
{
	stopifnot( class(seq)=='DNAbin', is.matrix(seq) )
	stopifnot( 'TAXON_ID'%in%names(df))
	setkey(df, TAXON_ID)
	tmp	<- setdiff( df[, TAXON_ID], rownames(seq) )
	if(length(tmp))
	{
		if(verbose)
			cat(paste('\nFound TAXON_IDs with no sequence, removing data.table entries, n=', length(tmp)))
		df	<- df[rownames(seq),]
	}
	tmp	<- setdiff( rownames(seq), df[, TAXON_ID] )
	if(length(tmp))
	{
		if(verbose)
			cat(paste('\nFound sequences with no TAXON_ID in data.table, removing sequences, n=', length(tmp)))
		seq	<- seq[df[, TAXON_ID],]
	}
	tmp		<- as.character(seq)
	tmp		<- data.table( SEQ=apply(tmp,1,function(x) paste(x, collapse='')), TAXON_ID=rownames(tmp) )	
	merge(df, tmp, by='TAXON_ID')		
}
######################################################################################
#	pool clusters into sets containing roughly 'pool.ntip' sequences
#' @title Pick sequence data set by type of phylogenetic cluster  
#' @export
beast.choose.seq.by.clusters<- function(df.seq, select, verbose=1)
{
	thresh.NSEQ		<- as.numeric(substring(regmatches(select, regexpr('seq[0-9]+',select)), 4))
	stopifnot(!is.na(thresh.NSEQ))
	#	'how to select' options
	#	pick phylogenetic clusters randomly
	if(grepl('clrndseq',select))
	{		
		seq.select		<- df.seq[, list(CLU_N=length(TAXON_NAME)), by='CLU_ID']			
		setkey(seq.select, CLU_N, CLU_ID)
		seq.select		<- subset(unique(seq.select), CLU_N>2)
		seq.select[, DUMMY:= sample(nrow(seq.select), nrow(seq.select))]
		setkey(seq.select, DUMMY)
		seq.select[, CLU_CN:= seq.select[, cumsum(CLU_N)]]
		seq.select		<- seq.select[seq_len( seq.select[, which(CLU_CN>=thresh.NSEQ)[1]] ), ]			
		seq.select		<- merge( df.seq, subset(seq.select, select=CLU_ID), by='CLU_ID' )				
	}
	#	or pick sequences from largest phylogenetic clusters
	if(grepl('mseq',select))
	{	
		seq.select		<- df.seq[, list(CLU_N=-length(TAXON_NAME)), by='CLU_ID']			
		setkey(seq.select, CLU_N, CLU_ID)
		seq.select		<- subset(unique(seq.select), CLU_N< -2)	
		seq.select[, CLU_CN:= seq.select[, cumsum(-CLU_N)]]
		seq.select		<- seq.select[seq_len( seq.select[, which(CLU_CN>=thresh.NSEQ)[1]] ), ]			
		seq.select		<- merge( df.seq, subset(seq.select, select=CLU_ID), by='CLU_ID' )									
	}
	#	or pick sequences from smallest phylogenetic clusters
	if(grepl('clsmseq',select))
	{
		seq.select		<- df.seq[, list(CLU_N=length(TAXON_NAME)), by='CLU_ID']			
		setkey(seq.select, CLU_N, CLU_ID)
		seq.select		<- subset(unique(seq.select), CLU_N>2)	
		seq.select[, CLU_CN:= seq.select[, cumsum(CLU_N)]]
		seq.select		<- seq.select[seq_len( seq.select[, which(CLU_CN>=thresh.NSEQ)[1]] ), ]			
		seq.select		<- merge( df.seq, subset(seq.select, select=CLU_ID), by='CLU_ID' )					
	}
	#	or pick sequences of all phylogenetic clusters of size >= x
	if(grepl('cseq',select))
	{				
		seq.select		<- df.seq[, list(CLU_N=length(TAXON_NAME)), by='CLU_ID']			
		setkey(seq.select, CLU_N, CLU_ID)
		seq.select		<- subset(unique(seq.select), CLU_N>=thresh.NSEQ)	
		seq.select		<- merge( df.seq, subset(seq.select, select=CLU_ID), by='CLU_ID' )					
	}
	if(verbose)
	{
		cat(paste('\nFound clusters, n=', seq.select[, length(unique(CLU_ID))])) 
		cat(paste('\nFound sequences, n=', seq.select[, length(unique(TAXON_ID))]))					
	}
	seq.select
}
######################################################################################
#	pool clusters into sets containing roughly 'pool.ntip' sequences
#' @title Pool phylogenetic clusters for separate BEAST runs  
#' @description This function returns a list of data tables. Each data table contains the pooled taxa for a single BEAST run. See Examples.
#' @import XML phytools ape data.table reshape2 ggplot2 
#' @export
#' @example example/ex.beast.pool.cluster.R
#' @example example/ex.beast.from.template.R
#' @seealso \code{\link{beastxml.from.template}}
beast.pool.clusters<- function(cluphy.df, how=NA, verbose=1, ...)
{	
	stopifnot(!is.na(how))
	stopifnot(how%in%c('by.samplingtime', 'by.samplingtime.alwaysincludebeforeyear', 'by.required.seq.per.samplingperiod'))
	if(how=='by.samplingtime')
		return( beast.pool.clusters.by.samplingtime(cluphy.df, verbose=verbose, ...) )
	if(how=='by.samplingtime.alwaysincludebeforeyear')
		return( beast.pool.clusters.by.samplingtime.alwaysincludebeforeyear(cluphy.df, verbose=verbose, ...) )
	if(how=='by.required.seq.per.samplingperiod')
		return( beast.pool.clusters.by.requiredseq.per.samplingperiod(cluphy.df, verbose=verbose, ...) )
}
######################################################################################
beast.pool.clusters.by.samplingtime.alwaysincludebeforeyear<- function(cluphy.df, pool.ntip= 130, pool.includealwaysbeforeyear=1996, verbose=1)
{
	stopifnot( c('CLU_ID','SAMPLINGTIME','TAXON_ID')%in%names(cluphy.df) )
	df			<- subset(cluphy.df, select=c(TAXON_ID, CLU_ID, SAMPLINGTIME))	
	df			<- df[, list( clu.ntip=length(TAXON_ID), clu.T1=min(SAMPLINGTIME)), by="CLU_ID"]
	
	df.always	<- subset(df, clu.T1<pool.includealwaysbeforeyear)
	df			<- subset(df, clu.T1>=pool.includealwaysbeforeyear)	
	pool.ntip	<- pool.ntip - sum(df.always[,clu.ntip])		
	if(pool.ntip<50)
		warning('number of non-ancestral taxas small')
	stopifnot(pool.ntip>0)
	setkey(df, clu.T1)
	pool.n		<- ceiling( sum( df[,clu.ntip] ) / pool.ntip )
	tmp			<- lapply( seq_len(pool.n), function(x)	seq.int(x,nrow(df),by=pool.n) )			
	pool.df		<- lapply(seq_along(tmp), function(i) merge(subset(rbind(df.always, df[tmp[[i]],]), select=CLU_ID), cluphy.df, by="CLU_ID") )
	
	if(verbose)
	{		
		cat(paste("\ncall to= pool.clusters.by.samplingtime.alwaysincludebeforeyear"))
		cat(paste("\nrequested number of seqs per pool is n=", pool.ntip+sum(df.always[,clu.ntip])))
		cat(paste("\nalways include clusters starting before ",pool.includealwaysbeforeyear))
		cat(paste("\nnumber of pools is n=",pool.n))
		cat(paste("\nnumber of seq in pools is n=",paste( sapply(pool.df, nrow), sep='', collapse=', ' )))
	} 
	pool.df
}
######################################################################################
beast.pool.clusters.by.samplingtime<- function(cluphy.df, pool.ntip= 130, verbose=1)
{
	stopifnot( c('CLU_ID','SAMPLINGTIME','TAXON_ID')%in%names(cluphy.df) )
	df			<- subset(cluphy.df, select=c(TAXON_ID, CLU_ID, SAMPLINGTIME))	
	df			<- df[, list( clu.ntip=length(TAXON_ID), clu.T1=min(SAMPLINGTIME)), by="CLU_ID"]		
	setkey(df, clu.T1)
	pool.n		<- ceiling( sum( df[,clu.ntip] ) / pool.ntip )
	tmp			<- lapply( seq_len(pool.n), function(x)	seq.int(x, nrow(df),by=pool.n) )
	pool.df		<- lapply( seq_along(tmp), function(i) merge(subset(df[tmp[[i]],], select=CLU_ID), cluphy.df, by="CLU_ID") )		
	
	if(verbose)
	{
		cat(paste("\ncall to= pool.clusters.by.samplingtime"))
		cat(paste("\nrequested number of seqs per pool is n=", pool.ntip))
		cat(paste("\nnumber of pools is n=",pool.n))
		cat(paste("\nnumber of seq in pools is n=",paste( sapply(pool.df, nrow), sep='', collapse=', ' )))
	} 
	pool.df
}
######################################################################################
beast.pool.clusters.by.requiredseq.per.samplingperiod<- function(cluphy.df, pool.ntip.guide=150, pool.ntip.min=c(50, 70, 70, NA, NA), pool.breaks=c(0, 1.596, 3.596, 5.596, 9.596, Inf), verbose=1)
{	
	stopifnot( c('CLU_ID','SAMPLINGTIME','TAXON_ID')%in%names(cluphy.df) )
	df				<- copy(cluphy.df)
	#	add tip heights and tip periods 
	set(df, NULL, 'TipHeight', df[, max(SAMPLINGTIME)-SAMPLINGTIME] )
	set(df, NULL, 'TipPeriod', df[, cut(TipHeight, pool.breaks, right=FALSE)] )	
	#	first attempt of pooling
	clu.df			<- df[, list( clu.ntip=length(TAXON_ID), clu.T1=min(SAMPLINGTIME)), by="CLU_ID"]
	pool.n			<- ceiling( sum( clu.df[,clu.ntip] ) / pool.ntip.guide )
	setkey(clu.df, clu.T1)
	tmp				<- lapply( seq_len(pool.n), function(x)	seq.int(x,nrow(clu.df),by=pool.n) )
	#	always add the earliest and latest cluster to each pool
	df.mxclu		<- clu.df[, { tmp<- c(which.min(clu.T1),which.max(clu.T1)); list(CLU_ID= CLU_ID[tmp])}]
	pool.df			<- lapply( seq_along(tmp), function(i) unique(rbind(subset(clu.df[tmp[[i]],], select=CLU_ID), df.mxclu)) )			
	pool.df			<- lapply(seq_along(tmp), function(i) merge(pool.df[[i]], df, by="CLU_ID") )
	#	compute the intial numbers of selected sequences per time period
	cnts			<- sapply(seq_along(pool.df), function(i)	table( pool.df[[i]][,TipPeriod])		)
	if(verbose)
	{
		cat(paste("\ncall to= pool.clusters.by.requiredseq.per.samplingperiod"))
		cat(paste("\nrequested number of seqs per pool is n=", pool.ntip.guide))		
		cat('\ninitial number of sequences by tip height (row) and pool (col)\n')
		print(cnts)
	}		
	if( any( na.omit(rowSums(cnts)<pool.ntip.min) ) )	
		stop('total number of available sequences in one of the periods is smaller than required pool.ntip.min')	
	#	compute the number of sequences to be added so that pool.ntip.min is satisfied
	setkey(df, TAXON_ID)	
	cnts							<- pool.ntip.min - cnts
	cnts[is.na(cnts)|cnts<=0]		<- 0
	#	add sequences to pool	
	for(i in seq_along(pool.df))
		for(j in seq_len(nrow(cnts)))
			if(cnts[j,i]>0)
			{
				#cat(paste('\nadding sequences to pool',i,'for period',rownames(cnts)[j]))
				#	only if cnts>0
				#	determine clusters that can be added
				pool.df.notin			<- setdiff( subset(df,TipPeriod==rownames(cnts)[j])[,TAXON_ID], pool.df[[i]][,TAXON_ID] )
				clu.df.notin			<- unique( subset( df[pool.df.notin,], select=CLU_ID) )
				clu.df.notin			<- merge(df, clu.df.notin, by='CLU_ID')
				clu.df.notin			<- clu.df.notin[, list(clu.ntipperiod=length(which(TipPeriod==rownames(cnts)[j]))), by=CLU_ID]	
				#	determine clusters that will be added
				tmp						<- clu.df.notin[,tail(which(cumsum(clu.ntipperiod)<cnts[j,i]),1),]
				tmp						<- ifelse(length(tmp), tmp[1]+1, 1)
				clu.df.notin			<- clu.df.notin[ seq_len( min( nrow(clu.df.notin), tmp ) ), ]
				#	add clusters
				pool.df[[i]]			<- rbind( pool.df[[i]], merge( df, subset(clu.df.notin, select=CLU_ID), by='CLU_ID') )
				tmp						<- pool.ntip.min - as.numeric( table( pool.df[[i]][,TipPeriod]) )
				tmp[is.na(tmp)|tmp<0]	<- 0 
				cnts[,i]				<- tmp
				#cat(paste('\nnew number of sequences in pool',i,'is',nrow(pool.df[[i]])))
			}
	#	compute the final numbers of selected sequences per time period
	cnts	<- sapply(seq_along(pool.df), function(i)	table( pool.df[[i]][,TipPeriod])		)		
	if(verbose)
	{		
		cat('\nfinal number of sequences by tip height (row) and pool (col)\n')
		print(cnts)
	}		
	pool.df	
}
######################################################################################
#	add taxa and alignment in bxml from BEASTlabels in df and alignment in seq.PROT.RT
#	beast.label.datepos= 4; beast.label.sep= '_'; beast.date.direction= "forwards"; beast.date.units= "years"; beast.alignment.dataType= "nucleotide"; xml.resetTipDate2LastDiag=1
#' @title	Add taxa to XML file
#' @export
beast.add.taxa<- function(bxml, df, beast.date.direction= "forwards", beast.date.units= "years", verbose=1)
{			
	stopifnot(c('TAXON_NAME','SAMPLINGTIME')%in%names(df))
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	#	check if any taxa to be added
	bxml.beast.existing.taxa.id		<- unlist( xpathApply(bxml.beast, "taxa/taxon", xmlGetAttr, "id" ) )
	tmp								<- subset( df, !TAXON_NAME%in%bxml.beast.existing.taxa.id )
	#	add taxa if needed
	if(nrow(tmp))
	{		
		invisible(newXMLCommentNode(text="The list of taxa to be analysed (can also include dates/ages).", parent=bxml.beast, doc=bxml, addFinalizer=T))
		invisible(newXMLCommentNode(text=paste("ntax=",nrow(tmp),sep=''), parent=bxml.beast, doc=bxml, addFinalizer=T))	
		seqtaxa		<- newXMLNode("taxa", attrs= list(id="taxa"), parent=bxml.beast, doc=bxml, addFinalizer=T)
		invisible(tmp[, {
					taxon	<- newXMLNode("taxon", attrs= list(id=TAXON_NAME), parent=seqtaxa, doc=bxml, addFinalizer=T )
					dummy	<- newXMLNode("date", attrs= list(value=SAMPLINGTIME, direction=beast.date.direction, units=beast.date.units), parent=taxon, doc=bxml, addFinalizer=T )
					NULL
				}, by='TAXON_NAME'])
		if(verbose)	cat(paste("\nadded new seq taxa, n=", xmlSize(seqtaxa)))		
	}
	invisible(bxml)
}	
######################################################################################
#	add taxa and alignment in bxml from BEASTlabels in df and alignment in seq.PROT.RT
#	beast.label.datepos= 4; beast.label.sep= '_'; beast.date.direction= "forwards"; beast.date.units= "years"; beast.alignment.dataType= "nucleotide"; xml.resetTipDate2LastDiag=1
#' @export
beast.add.alignment<- function(bxml, df, df.seqfield='SEQ', beast.alignment.id="alignment", beast.alignment.dataType= "nucleotide", verbose=1)
{			
	stopifnot(df.seqfield%in%names(df))
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	tmp			<- unlist( xpathApply(bxml.beast, "taxa/taxon", xmlGetAttr, "id" ) )
	tmp			<- length(setdiff( df[, TAXON_NAME], tmp ))
	if( tmp )
		stop(paste('\nFound TAXON_NAME that are not an xml taxa. Call first beast.add.taxa. n=', tmp))
	if( length(getNodeSet(bxml, paste("//alignment[@id='",beast.alignment.id,"']", sep='')))>0 )
		stop(paste('\nAlready exists: alignment with id',beast.alignment.id )) 
	setnames(df, df.seqfield, 'x_THISSEQ')
	#	add alignment	
	invisible(newXMLCommentNode(text="The sequence alignment (each sequence refers to a taxon above).", parent=bxml.beast, doc=bxml, addFinalizer=T))
	invisible(newXMLCommentNode(text=paste("ntax=",nrow(df)," nchar=",nchar(df[[df.seqfield]][1]),sep=''), parent=bxml.beast, doc=bxml, addFinalizer=T))	
	seqalign	<- newXMLNode("alignment", attrs= list(id=beast.alignment.id, dataType=beast.alignment.dataType), parent=bxml.beast, doc=bxml, addFinalizer=T)
	invisible(df[, {
				seq		<- newXMLNode("sequence", parent=seqalign, doc=bxml, addFinalizer=T)
				invisible(newXMLNode("taxon", attrs= list(idref=TAXON_NAME), parent=seq, doc=bxml, addFinalizer=T))
				invisible(newXMLTextNode(text=x_THISSEQ, parent=seq, doc=bxml,addFinalizer=T))
				NULL
			}, by='TAXON_NAME'])	
	setnames(df, 'x_THISSEQ', df.seqfield)
	if(verbose)	cat(paste("\nadded new alignment with seqs, n=", xmlSize(seqalign)))	
	invisible(bxml)
}
######################################################################################
beast.add.alignment.fromDNAbin<- function(bxml, seq.PROT.RT, beast.alignment.id="alignment", beast.alignment.dataType= "nucleotide", verbose=1)
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
#' @title	Read BEAST ops file
#' @export
beast.read.ops<- function(file)
{
	z		<- readLines(file)
	#only data
	z		<- z[seq(4,length(z)-1)]
	#rm whitespace
	z		<- sapply(z, function(x) gsub('\\s+',' ',gsub("^\\s+|\\s+$", "", x) ) )
	names(z)<- NULL
	#
	z		<- strsplit(z,' ')
	ans		<- do.call('rbind',lapply(z, function(x)
			{
				tmp	<- suppressWarnings(as.numeric(x))
				op	<- paste(x[which(is.na(tmp))], collapse='_')
				tmp	<- tmp[!is.na(tmp)]
				if(length(tmp)<5)
					tmp<- c(NA, tmp)
				data.table(Operator=op, Tuning=tmp[1], Count=tmp[2], Time=tmp[3], TimeOp=tmp[4], PrAccept=tmp[5])				
			}))
	ans
}
######################################################################################
#' @title	Read BEAST log file
#' @export
beast.read.log<- function(file, select=c('state','likelihood'), verbose=1)
{
	df		<- as.data.table(read.delim(file, comment.char='#'))
	if(verbose)
		cat(paste('\nFound columns', paste(names(df), collapse=',')))
	tmp		<- as.logical(apply(sapply(select, function(x)  as.numeric(grepl(x, names(df)))), 1, sum))
	if(verbose)
		cat(paste('\nSelect columns', paste(names(df)[tmp], collapse=',')))		
	subset(df, select=names(df)[tmp])		
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
#' @title Add monophylyStatistics for each phylogenetic cluster to XML file
#' @export
beast.add.monophylyStatistic.for.taxonsets<- function(bxml, treeModel.id, taxonset.prefix='c') 
{		
	btaxonset.id	<- xpathSApply(bxml, paste("//taxa[starts-with(@id,'",taxonset.prefix,"')]", sep=''), xmlGetAttr, "id")
	if(length(btaxonset.id)==0)
		stop(paste('Cannot find taxonsets with prefix that identifies cluster taxonsets, =',taxonset.prefix))
	#	create monophyly statistics
	tmp				<- lapply(btaxonset.id, function(x)
			{
				monophylyStatistic	<- newXMLNode("monophylyStatistic", attrs= list(id=paste("monophyly_",x,sep='')), doc=bxml, addFinalizer=T )
				mrca			<- newXMLNode("mrca", parent=monophylyStatistic, doc=bxml, addFinalizer=T)
				invisible(newXMLNode("taxa", attrs= list(idref=x), parent=mrca, doc=bxml, addFinalizer=T ))
				invisible(newXMLNode("treeModel", attrs= list(idref=treeModel.id), parent=monophylyStatistic, doc=bxml, addFinalizer=T ))
				monophylyStatistic					
			})
	# 	link monophyly statistics
	bxml.beast		<- getNodeSet(bxml, "//beast")[[1]]
	bxml.idx		<- which(xmlSApply(bxml.beast, xmlName)=="tmrcaStatistic")
	if(length(bxml.idx)==0)
		stop('Cannot find tmrcaStatistic.')
	invisible(addChildren(bxml.beast, tmp, at=bxml.idx[length(bxml.idx)] ))
	if(verbose) cat(paste("\nadded monophylyStatistics for mrca of each cluster, n=",length(tmp)))
	# 	add monophyly statistics to prior
	bxml.prior				<- getNodeSet(bxml, "//*[@id='prior']")
	if(length(bxml.prior)>1)	
		stop("Found more than one prior")
	if(length(bxml.prior)==0)	
		stop("Cannot find prior")	
	bxml.prior				<- bxml.prior[[1]]
	#	see if there is a booleanLikelihood, and if yes use it, otherwise add a new XML node
	bxml.bool.lkl			<- getNodeSet(bxml.prior,"//booleanLikelihood")
	if(length(bxml.bool.lkl)>1)	
		stop("Found more than one booleanLikelihood")
	if(length(bxml.bool.lkl)==1)
		bxml.bool.lkl		<- bxml.bool.lkl[[1]]
	if(length(bxml.bool.lkl)<1)
		bxml.bool.lkl		<- newXMLNode("booleanLikelihood", parent=bxml.prior, doc=bxml, addFinalizer=T )
	tmp						<- xpathSApply(bxml, paste("//monophylyStatistic[starts-with(@id,'","monophyly_",taxonset.prefix,"')]", sep=''), xmlGetAttr, "id")	
	invisible(lapply(tmp, function(x)
			{
				newXMLNode("monophylyStatistic", attrs= list(idref=x), parent=bxml.bool.lkl, doc=bxml, addFinalizer=T )
			}))
	if(verbose) cat(paste("\nadded monophylyStatistics to booleanLikelihood, n=",length(tmp)))
	invisible(bxml) 
}	
######################################################################################
#	For each taxonset, get a list of tmrcaStatistics. Assumes all tmrcaStatistics share a treeModel.id and the same includeStem attribute
#' @title	Add tmrcaStatistic for phylogenetic clusters to XML file 
#' @export
beast.add.tmrcaStatistic.for.taxonsets<- function(bxml,  treeModel.id='treeModel', taxonset.prefix='c', includeStem="false", verbose=1) 
{			
	btaxonset.id	<- xpathSApply(bxml, paste("//taxa[starts-with(@id,'",taxonset.prefix,"')]", sep=''), xmlGetAttr, "id")
	if(length(btaxonset.id)==0)
		stop(paste('Cannot find taxonsets with prefix that identifies cluster taxonsets, =',taxonset.prefix))
	prefix.id		<- ifelse(includeStem=="false","tmrca","tstem")
	#	create trmca statistics
	tmp				<- lapply(btaxonset.id, function(x)
			{
				tmrcaStatistic	<- newXMLNode("tmrcaStatistic", attrs= list(id=paste(prefix.id,"_",x,sep=''), includeStem=includeStem), doc=bxml, addFinalizer=T )
				mrca			<- newXMLNode("mrca", parent=tmrcaStatistic, doc=bxml, addFinalizer=T)
				invisible(newXMLNode("taxa", attrs= list(idref=x), parent=mrca, doc=bxml, addFinalizer=T ))
				invisible(newXMLNode("treeModel", attrs= list(idref=treeModel.id), parent=tmrcaStatistic, doc=bxml, addFinalizer=T ))
				tmrcaStatistic					
			})	
	# 	link trmca statistics to XML	
	bxml.beast					<- getNodeSet(bxml, "//beast")[[1]]
	bxml.idx					<- which(xmlSApply(bxml.beast, xmlName)=="treeModel")
	if(length(bxml.idx)==0)
		stop('Cannot find treeModel')
	invisible(addChildren(bxml.beast, tmp, at=bxml.idx[length(bxml.idx)] ))
	if(verbose) cat(paste("\nadded tmrcaStatistics for each cluster, n=",length(tmp)))
	#	add tmrcaStatistics to log 
	bxml.fileLog				<- getNodeSet(bxml, "//log[contains(@id,'fileLog')]")
	if(length(bxml.fileLog)>1)
		stop('Found more than one fileLog')
	if(length(bxml.fileLog)==0)
		warning("Could not find fileLog whose id contains 'fileLog'. Not adding trmcaStatistics to fileLog")
	if(length(bxml.fileLog)==1)
	{
		tmrcaStatistics.id		<- xpathSApply(bxml, paste("//tmrcaStatistic[starts-with(@id,'",prefix.id,"_",taxonset.prefix,"')]", sep=''), xmlGetAttr, "id")
		stopifnot(length(tmrcaStatistics.id)>0)
		invisible(lapply(tmrcaStatistics.id, function(x)	newXMLNode("tmrcaStatistic", attrs= list(idref=x), parent=bxml.fileLog, doc=bxml, addFinalizer=T )	))
		if(verbose) cat(paste("\nadded tmrcaStatistics to fileLog, n=", length(tmrcaStatistics.id)))
	}
	invisible(bxml) 
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
#' @title Copy XML elements from template
#' @export
beast.add.from.template<- function(bxml, btemplate, select.after.elt=NA, select.to.elt=NA, verbose=1)
{
	bxml.beast			<- getNodeSet(bxml, "//beast")[[1]]
	bt.beast			<- getNodeSet(btemplate, "//beast")[[1]]
	bt.beast.elt.nms	<- xmlSApply(bt.beast, xmlName)
	select.after		<- 0
	if(!is.na(select.after.elt))
	{
		stopifnot( select.after.elt%in%bt.beast.elt.nms )
		select.after	<- which(grepl(select.after.elt,bt.beast.elt.nms))			
	}	
	select.up.to		<- length(bt.beast.elt.nms)
	if(!is.na(select.to.elt))
	{
		stopifnot( select.to.elt%in%bt.beast.elt.nms )
		select.up.to	<- which(grepl(select.to.elt,bt.beast.elt.nms))			
	}			
	stopifnot(select.after<select.up.to)
	tmp					<- which( !grepl('comment', bt.beast.elt.nms[seq.int(select.after+1,select.up.to)] ) )
	#	see if starting tree has been added already	
	bxml.nwck.id		<- xpathSApply(bxml.beast, "//newick[count(@id)]", xmlGetAttr, "id")
	#	if yes, make sure same id is not added through template 
	if(length(bxml.nwck.id))
	{
		tmp2			<- sapply(tmp, function(i) length(xpathSApply(bt.beast[[i]], paste("descendant-or-self::*[@id='",bxml.nwck.id,"']",sep=''), xmlGetAttr, "id" ))==0 )
		tmp				<- tmp[tmp2]		
	}
	if(verbose)
		cat(paste('\nCopy elements from template=', paste(bt.beast.elt.nms[tmp], collapse=', ', sep='')))
	invisible(sapply(tmp, function(i)
					{
						if( class(bt.beast[[i]])[1]=="XMLInternalCommentNode" )
							dummy<- newXMLCommentNode(text=xmlValue(bt.beast[[i]]), parent=bxml.beast, doc=bxml, addFinalizer=T)
						else
							dummy<- addChildren( bxml.beast, xmlClone( bt.beast[[i]], addFinalizer=T, doc=bxml ) )
					}))
	invisible(bxml)			
}
######################################################################################
#' @title	Adjust dimension parameter of gmrfSkyrideLikelihood
#' @export
beast.adjust.gmrfSkyrideLikelihood<- function(bxml, df, verbose=1)
{
	bxml.tmp	<- getNodeSet(bxml, "//gmrfSkyrideLikelihood[count(@id)]")
	if(length(bxml.tmp)>1)
		stop("unexpected number of gmrfSkyrideLikelihood")
	if(length(bxml.tmp)==1)
	{
		tmp			<- getNodeSet(bxml.tmp[[1]], "//*[count(@dimension)]")
		for(x in tmp)
		{
			if(verbose)
				cat(paste('\nAdjusting dimension of ',xmlName(x),'with id',xmlGetAttr(x,'id'),'to',nrow(df)-1))
			xmlAttrs(x)["dimension"]	<-	nrow(df)-1  
		}
	}		
	invisible(bxml)
}
######################################################################################
#' @title	Adjust MCMC element
#' @export
beast.adjust.mcmc<- function(bxml, mcmc.chainLength=NA, mcmc.logEvery=NA, mcmc.outfile=NA, verbose=1)
{
	bxml.tmp		<- getNodeSet(bxml, "//mcmc[count(@id)]")
	if(length(bxml.tmp)>1)
		stop("unexpected number of mcmc")	
	if(!is.na(mcmc.chainLength))
	{
		if(verbose)	
			cat(paste("\nAdjust MCMC chainLength to",mcmc.chainLength))
		xmlAttrs(bxml.tmp[[1]])["chainLength"]	<-	sprintf("%d",mcmc.chainLength)		
	}
	if(!is.na(mcmc.logEvery))
	{
		tmp			<- getNodeSet(bxml.tmp[[1]], "//*[@logEvery]")		
		for(x in tmp)
		{
			if(verbose)	
				cat(paste("\nAdjust logEvery of",xmlName(x),'with id',xmlGetAttr(x,'id'),"to", mcmc.logEvery))
			xmlAttrs(x)["logEvery"]	<- sprintf("%d", mcmc.logEvery)
		}			
	}
	if(!is.na(mcmc.outfile))
	{
		if(regexpr('.xml',mcmc.outfile)<0)
			stop('Expected file name ending in .xml')
		#
		xmlAttrs(bxml.tmp[[1]])["operatorAnalysis"]	<- gsub('.xml','.ops',mcmc.outfile)
		#
		tmp		<- getNodeSet(bxml.tmp[[1]], "//log[contains(@id,'fileLog')]")
		if(length(tmp)>1)
			stop('Found more than one fileLog')
		x		<- tmp[[1]]
		xmlAttrs(x)["fileName"]	<- gsub('.xml','.log',mcmc.outfile)
		if(verbose)	
			cat(paste("\nAdjust fileName of",xmlName(x),'with id',xmlGetAttr(x,'id'),"to", xmlAttrs(x)["fileName"]))
		#
		tmp		<- getNodeSet(bxml.tmp[[1]], "//logTree[count(@branchLengths)=0]")
		if(length(tmp)>1)
			stop('Found more than one subst logTree')
		x		<- tmp[[1]]
		xmlAttrs(x)["fileName"]	<- gsub('.xml','.tree',mcmc.outfile)
		if(verbose)	
			cat(paste("\nAdjust fileName of",xmlName(x),'with id',xmlGetAttr(x,'id'),"to", xmlAttrs(x)["fileName"]))
		#
		tmp		<- getNodeSet(bxml.tmp[[1]], "//logTree[@branchLengths='substitutions']")
		if(length(tmp)>1)
			stop('Found more than one subst logTree')
		x		<- tmp[[1]]
		xmlAttrs(x)["fileName"]	<- gsub('.xml','_subst.tree',mcmc.outfile)
		if(verbose)	
			cat(paste("\nAdjust fileName of",xmlName(x),'with id',xmlGetAttr(x,'id'),"to", xmlAttrs(x)["fileName"]))		
	}
	invisible(bxml)
}
######################################################################################
#	if rootheight prior uniform, sets lower bound to earliest sampling time in data set
#' @title Ajust lower bound of uniform prior on root height
#' @export
beast.adjust.prior.rootheight<- function(bxml, df, rootHeight.idref='treeModel.rootHeight', verbose=1)
{
	bxml.prior	<- getNodeSet(bxml, paste("//uniformPrior[descendant::parameter[@idref=",rootHeight.idref,"]]",sep=''))
	if(length(bxml.prior)>1)	
		stop('Found more than one rootHeight with same id.')
	if(length(bxml.prior)==1)
	{
		tmp									<- df[, range(SAMPLINGTIME)]
		tmp									<- floor(diff(tmp))
		xmlAttrs(bxml.prior[[1]])["lower"]	<- tmp		
		if(verbose)	
			cat(paste("\nset lower bound of uniformPrior for rootHeight to", tmp))			
	}
	invisible(bxml)
}
######################################################################################
#	For each cluster, create a taxonset. Assumes df has BEASTlabel and cluster
#' @export
beast.add.taxonsets.for.clusters	<- function(bxml, df, taxonset.prefix='c')
{	
	stopifnot(c('CLU_ID','TAXON_NAME')%in%names(df))
	#	create taxonsets
	tmp	<- lapply( unique( df[,CLU_ID] ), function(clu)
			{							
				taxonset	<- newXMLNode("taxa", attrs= list(id=paste(taxonset.prefix,clu,sep='')), doc=bxml, addFinalizer=T )				
				invisible( lapply(subset(df, CLU_ID==clu)[, TAXON_NAME], function(x)
								{
									newXMLNode("taxon", attrs= list(idref=x), parent=taxonset, doc=bxml, addFinalizer=T )	
								})	)
				taxonset
			})
	tmp	
	# 	link taxonsets to XML	
	bxml.beast					<- getNodeSet(bxml, "//beast")[[1]]
	bxml.idx					<- which(xmlSApply(bxml.beast, xmlName)=="taxa")
	if(length(bxml.idx)==0)
		stop('Cannot find taxa')
	invisible(addChildren(bxml.beast, tmp, at=bxml.idx[length(bxml.idx)] ))
	if(verbose) cat(paste("\nadded taxon sets for each cluster, n=",length(tmp)))
	invisible(bxml)
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
beast.add.upDownOperator<- function(bxml, up.parameter.id, down.parameter.id, scaleFactor=0.75, weight=1)
{
	bxml.o			<- getNodeSet(bxml, "//operators")[[1]]
	tmp				<- newXMLNode("upDownOperator", attrs= list(scaleFactor=as.character(scaleFactor), weight=as.character(weight)), parent=bxml.o, doc=bxml, addFinalizer=T)
	tmp2			<- newXMLNode("up", parent=tmp, doc=bxml, addFinalizer=T)
	for(x in up.parameter.id)
		invisible(newXMLNode("parameter", attrs= list(idref=x), parent=tmp2, doc=bxml, addFinalizer=T))
	tmp2			<- newXMLNode("down", parent=tmp, doc=bxml, addFinalizer=T)
	for(x in down.parameter.id)
		invisible(newXMLNode("parameter", attrs= list(idref=x), parent=tmp2, doc=bxml, addFinalizer=T))	
	invisible(bxml) 
}
######################################################################################
#' @export
beast.add.uniformOperator<- function(bxml, parameter.id, weight=1)
{
	bxml.o			<- getNodeSet(bxml, "//operators")[[1]]
	tmp				<- newXMLNode("uniformOperator", attrs= list(weight=as.character(weight)), parent=bxml.o, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(idref=parameter.id), parent=tmp, doc=bxml, addFinalizer=T)
	invisible(bxml) 
}
######################################################################################
#' @export
beast.add.deltaExchangeOperator<- function(bxml, parameter.id, delta=0.75, parameterWeights=c(948, 948, 948), weight=1)
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
beast.add.scaleOperator<- function(bxml, parameter.id, scaleFactor=0.75, weight=0.1, autoOptimize="true")
{
	bxml.o			<- getNodeSet(bxml, "//operators")[[1]]
	tmp				<- newXMLNode("scaleOperator", attrs= list(scaleFactor=as.character(scaleFactor), weight=as.character(weight), autoOptimize=autoOptimize), parent=bxml.o, doc=bxml, addFinalizer=T)
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
beast.add.logNormalPrior<- function(bxml, idref, mean, stdev, offset=0, meanInRealSpace='false')
{
	bxml.beast		<- getNodeSet(bxml, "//prior")[[1]]
	tmp				<- newXMLNode("logNormalPrior", attrs= list(mean=as.character(mean), stdev=as.character(stdev), offset=as.character(offset), meanInRealSpace=meanInRealSpace), parent=bxml.beast, doc=bxml, addFinalizer=T)
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
	tmp				<- newXMLNode("rateStatistic", attrs= list(id=paste(prefix.id,'_mean',sep=''), name=paste(prefix.id,'_mean',sep=''), mode='mean', internal='true', external='true'), parent=bxml.beast, doc=bxml, addFinalizer=T)
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
#' @title Generate starting tree   
#' @description Generate starting tree with tips that correspond to the taxa in \code{df}. See Examples.
#' @return Starting tree in newick format
#' @import data.table adephylo 
#' @example example/ex.beast.from.template.R
#' @seealso \code{\link{beastxml.from.template}}
#' @export 
beast.get.startingtree<- function(ph, df, starttree.rootHeight=NA, origin.value=NA, verbose=1)
{
	stopifnot(c('TAXON_ID','TAXON_NAME')%in%names(df))
	tmp		<- c( 	length(intersect( ph$tip.label, df[, TAXON_ID] )),
					length(intersect( ph$tip.label, df[, TAXON_NAME] ))	)	#tree must contain all taxa in df
	stopifnot(any(tmp))	#	tip labels must either match ID or NAME
	if(tmp[1])
	{
		if(verbose)
			cat('Found tip labels that match TAXON_ID\nset tip labels to TAXON_NAME')
		stopifnot( tmp[1]==nrow(df) )		#tree must contain all taxa in df
		tmp2				<- setdiff( ph$tip.label, df[, TAXON_ID] )
		tmp2				<- match( tmp2, ph$tip.label)
		ph.start			<- drop.tip(ph, tmp2)		
		ph.start$node.label	<- NULL
		setkey(df, TAXON_ID)
		ph.start$tip.label	<- df[ph.start$tip.label,][, TAXON_NAME]		
	}
	if(tmp[2])
	{
		if(verbose)
			cat('Found tip labels that match TAXON_NAME\nkeep tip labels as is')
		if(tmp[2]!=nrow(df))
			stop('Found taxa in df that are not in tree'  )		#tree must contain all taxa in df
		tmp2				<- setdiff( ph$tip.label, df[, TAXON_NAME] )
		tmp2				<- match( tmp2, ph$tip.label)
		ph.start			<- drop.tip(ph, tmp2)		
		ph.start$node.label	<- NULL
	}
	#	adjust rootHeight to 'beast.rootHeight'	
	if(!is.na(starttree.rootHeight))
	{
		if(verbose) cat(paste("\ncreate startingTree with root height=",starttree.rootHeight))
		ph.start$edge.length<- ph.start$edge.length	*	starttree.rootHeight / max(distRoot(ph.start))	
	}
	if(!is.na(origin.value) & starttree.rootHeight>=origin.value)	
		stop('Detected invalid rootHeight. Must be smaller than origin.value.')
	if(verbose) 
	{
		cat(paste("\nselected tips for startingTree, n=",Ntip(ph.start)))
		cat(paste("\nroot height=",max(distRoot(ph.start))))
	}
	write.tree( ph.start )		
}
######################################################################################
# 	extract starting tree from 'ph' by tips in 'df'. Only keeps tree topology and resets branch lengths so that the maximum root distance is 'beast.rootHeight'
#	beast.rootHeight= 35; beast.usingDates= "false"; beast.newickid= "startingTree"
#' @title Add starting tree to XML
#' @export
beast.add.startingtree<- function(bxml, start.tree.newick, beast.newickid= "startingTree", beast.usingDates="true", beast.brlunits="years", verbose=1)
{
	if( class(start.tree.newick)=='phylo' )
		start.tree.newick	<- write.tree( start.tree.newick )	
	stopifnot(class(start.tree.newick)=='character')
	bxml.beast				<- getNodeSet(bxml, "//beast")[[1]]
	if( length(getNodeSet(bxml, paste("//newick[@id='",beast.newickid,"']", sep='')))>0 )
		stop(paste('\nAlready exists: newick with id',beast.newickid )) 
	invisible(newXMLCommentNode(text="The user-specified starting tree in a newick tree format", parent=bxml.beast, doc=bxml, addFinalizer=T))
	bxml.startingTree		<- newXMLNode("newick", attrs= list(id=beast.newickid, usingDates=beast.usingDates, units=beast.brlunits), parent= bxml.beast, doc=bxml, addFinalizer=T)
	invisible(newXMLTextNode(text=start.tree.newick, parent=bxml.startingTree, doc=bxml, addFinalizer=T)) 
	invisible(bxml)
}