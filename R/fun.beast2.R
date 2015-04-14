
######################################################################################
beast2.add.data<- function(bxml, seq.PROT.RT, df, beast2.spec, verbose=1)
{
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	dummy		<- newXMLCommentNode(text="The sequences.", parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text=paste("ntax=",nrow(df)," nchar=",ncol(seq.PROT.RT),sep=''), parent=bxml.beast, doc=bxml, addFinalizer=T)	
	seqalign	<- newXMLNode("data", attrs= list(id=beast2.spec$data.id, dataType=beast2.spec$data.dataType, missing=beast2.spec$data.missing), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- lapply( seq_len(nrow(df)), function(i)
			{
				tmp		<- which( rownames(seq.PROT.RT)==df[i, FASTASampleCode] )
				if(length(tmp)!=1)	stop("unexpected row in seq.PROT.RT selected")
				tmp		<- paste(as.character(seq.PROT.RT[tmp,])[1,],collapse='',sep='')										
				seq		<- newXMLNode("sequence", attrs=list(id= paste('seq',df[i, BEASTlabel],sep=''), taxon=df[i, BEASTlabel], value=tmp), parent=seqalign, doc=bxml, addFinalizer=T)
				seq
			})	
	if(verbose)	cat(paste("\nadded new sequences, n=", xmlSize(seqalign)))
	bxml
}
######################################################################################
#	df<- cluphy.df
beast2.poolclusters.mincnts<- function(df, beast2.spec, verbose=1)
{
	if(verbose) cat(paste("\npool evenly across clu.AnyPos_T1 and fill sequences so that min cnts requested per time period is met"))
	pool.ntip		<- beast2.spec$pool.ntip
	cnts.requested	<- beast2.spec$pool.cnts.requested
	breaks			<- c(Inf, beast2.spec$bdsky.sprop.changepoint.value)
	#	add tip heights and tip periods to cluphy.df
	tmp				<- df[, max(PosSeqT)]
	tmp				<- as.numeric( df[, difftime(tmp, PosSeqT, units='days') / 365] )	
	df				<- merge( df, data.table(TipHeight=tmp, TipPeriod=as.character(cut(tmp, breaks, right=FALSE)), FASTASampleCode=df[,FASTASampleCode]), by='FASTASampleCode')	
	#	first attempt of pooling
	clu.df			<- df[, list(clu.ntip=clu.ntip[1], clu.AnyPos_T1=clu.AnyPos_T1[1]), by="cluster"]
	df.mxclu		<- clu.df[, { tmp<- c(which.min(clu.AnyPos_T1),which.max(clu.AnyPos_T1)); list(cluster= cluster[tmp])}]		
	setkey(df, clu.AnyPos_T1)
	pool.n			<- ceiling( sum( clu.df[,clu.ntip] ) / pool.ntip )
	tmp				<- lapply( seq_len(pool.n), function(x)	seq.int(x,nrow(clu.df),by=pool.n) )
	#	select the clusters evenly across clu.AnyPos_T1 and always add the earliest and latest cluster
	pool.df			<- lapply(seq_along(tmp), function(i) unique(rbind(subset(clu.df[tmp[[i]],], select=cluster), df.mxclu)) )		
	pool.df			<- lapply(seq_along(tmp), function(i) merge(pool.df[[i]], df, by="cluster") )
	#	compute the intial numbers of selected sequences per time period
	cnts			<- sapply(seq_along(pool.df), function(i)	table( pool.df[[i]][,TipPeriod])		)
	if(verbose)		print(cnts)
	if( any( na.omit(rowSums(cnts)<cnts.requested) ) )	
		stop('total number of sequences is smaller than cnts.requested')
	#	compute the number of sequences to be added
	setkey(df, FASTASampleCode)	
	cnts							<- cnts.requested - cnts
	cnts[is.na(cnts)|cnts<=0]		<- 0
	#	add sequences to pool	
	for(i in seq_along(pool.df))
		for(j in seq_len(nrow(cnts)))
			if(cnts[j,i]>0)
			{
				cat(paste('\nadding sequences to pool',i,'for period',rownames(cnts)[j]))
				#	only if cnts>0
				#	determine clusters that can be added
				pool.df.notin			<- setdiff( subset(df,TipPeriod==rownames(cnts)[j])[,FASTASampleCode], pool.df[[i]][,FASTASampleCode] )
				clu.df.notin			<- unique( subset( df[pool.df.notin,], select=cluster) )
				clu.df.notin			<- merge(df, clu.df.notin, by='cluster')
				clu.df.notin			<- clu.df.notin[, list(clu.ntipperiod=length(which(TipPeriod==rownames(cnts)[j]))), by=cluster]	
				#	determine clusters that will be added
				tmp						<- clu.df.notin[,tail(which(cumsum(clu.ntipperiod)<cnts[j,i]),1),]
				tmp						<- ifelse(length(tmp), tmp[1]+1, 1)
				clu.df.notin			<- clu.df.notin[ seq_len( min( nrow(clu.df.notin), tmp ) ), ]
				#	add clusters
				pool.df[[i]]			<- rbind( pool.df[[i]], merge( df, subset(clu.df.notin, select=cluster), by='cluster') )
				tmp						<- cnts.requested - as.numeric( table( pool.df[[i]][,TipPeriod]) )
				tmp[is.na(tmp)|tmp<0]	<- 0 
				cnts[,i]				<- tmp
				cat(paste('\nnew number of sequences in pool',i,'is',nrow(pool.df[[i]])))
			}
	#	compute the final numbers of selected sequences per time period
	cnts	<- sapply(seq_along(pool.df), function(i)	table( pool.df[[i]][,TipPeriod])		)
	if(verbose)		print(cnts)
	list(pool.df=pool.df, pool.ntip=pool.ntip)
}
######################################################################################
beast2.extract.distinct.topologies<- function(mph.clu)					
{
	#compute if topology between retained mph.clu's is identical (branch lengths may differ)
	tmp				<- sapply( seq_along(mph.clu)[-1], function(mph.i) all.equal(mph.clu[[mph.i-1]], mph.clu[[mph.i]], use.edge.length=FALSE, use.tip.label=TRUE, index.return=FALSE ) )
	tmp				<- data.table( mph.i= seq_along(mph.clu)[-1], equal.to.previous= tmp )				
	#	dtopo -> distinct topologies
	mph.clu.dtopo	<- data.table( mph.i= c(1, subset(tmp, !equal.to.previous)[, mph.i]), freq= as.vector(table(cumsum(!c(TRUE,tmp[,equal.to.previous])))), collapsed=FALSE	)					
	#	itopo -> identical topologies 
	mph.clu.itopo	<- cbind(tmp, equal.to= 1+cumsum(!tmp[,equal.to.previous]))
	mph.clu.itopo	<- rbind( data.table(mph.i=1, equal.to=1), subset(mph.clu.itopo, select=c(mph.i, equal.to)) )
	#	select by index
	set(mph.clu.itopo, NULL, 'equal.to', mph.clu.dtopo[mph.clu.itopo[,equal.to],][,mph.i])
	#	now select by key				
	setkey(mph.clu.dtopo, mph.i)
	#for each sequentially different topology, check if subsequent ones are identical and if yes collapse					
	while(mph.clu.dtopo[,any(!collapsed)])
	{
		cat(paste('\nprogress: number of seq distinct topologies =',nrow(mph.clu.dtopo)))
		mph.indexrow	<- mph.clu.dtopo[, which(!collapsed)][1]
		if(mph.indexrow<nrow(mph.clu.dtopo))
		{
			mph.index		<- mph.clu.dtopo[mph.indexrow, mph.i]
			mph.is			<- mph.clu.dtopo[seq.int(mph.indexrow+1,nrow(mph.clu.dtopo)), ][,mph.i]						
			tmp				<- sapply( mph.is, function(mph.i) all.equal(mph.clu[[mph.index]], mph.clu[[mph.i]], use.edge.length=FALSE, use.tip.label=TRUE, index.return=FALSE ) )
			topo.duplicates	<- data.table( mph.index= mph.index, mph.i=mph.is, equal.to.index= tmp )
			topo.duplicates	<- merge( mph.clu.dtopo, subset(topo.duplicates, equal.to.index), by='mph.i' )
			if(nrow(topo.duplicates))
			{
				tmp				<- topo.duplicates[, list(index=seq_len(freq)), by='mph.i']
				set(mph.clu.itopo, as.integer( tmp[, mph.i+index-1] ), 'equal.to', mph.index)
				set(mph.clu.dtopo, mph.indexrow, 'freq', mph.clu.dtopo[mph.indexrow, freq] + topo.duplicates[, sum(freq)])							
				collapsed.i		<- setdiff(mph.clu.dtopo[,mph.i], topo.duplicates[,mph.i])
				mph.clu.dtopo	<- mph.clu.dtopo[J(collapsed.i),]	
			}							
		}
		set(mph.clu.dtopo, mph.indexrow, 'collapsed', TRUE)
	}
	if( mph.clu.dtopo[, sum(freq)]!=length(mph.clu) )	stop('unexpected freq in mph.clu.dtopo')
	if( !setequal(mph.clu.dtopo[, mph.i], mph.clu.itopo[, unique(equal.to)]) )  stop('unexpected mph.i difference between mph.clu.dtopo and mph.clu.itopo')
	cat(paste('\nnumber of distinct topologies in mph.clu.dtopo=',nrow(mph.clu.dtopo)))
	cat(paste('\nnumber of distinct topologies in mph.clu.itopo=',length(unique(mph.clu.itopo[, equal.to]))))
	
	mph.clu.dtopo	<- subset(mph.clu.dtopo, select=c(mph.i, freq))	
	mph.clu.dtopo[, dens:= round(mph.clu.dtopo[, freq]/length(mph.clu),d=3) ]	
	mph.clu.dtopo[, topo.n:= nrow(mph.clu.dtopo) ]
	mph.clu.dtopo[, tip.n:= Ntip(mph.clu[[1]]) ]			
	#	creates warnings because list(DT) creates a copy of the data.table
	list(dtopo= mph.clu.dtopo, itopo=mph.clu.itopo)
}
######################################################################################
beast2.add.alignment<- function(bxml, beast2.spec, verbose=1)
{
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	dummy		<- newXMLCommentNode(text="The sequence alignments.", parent=bxml.beast, doc=bxml, addFinalizer=T)
	if(length(beast2.spec$alignment.filter)==1 && is.na(beast2.spec$alignment.filter))
	{		
		dummy	<- newXMLCommentNode(text="The <alignment> is the <data> block if there is no alignment filter.", parent=bxml.beast, doc=bxml, addFinalizer=T)
		dummy	<- newXMLNode("taxa", attrs= list(id=paste("TaxonSet.t",beast2.spec$alignment.id,sep=':'), spec=beast2.spec$alignment.taxa.spec, alignment=paste('@',beast2.spec$alignment.id,sep='')), parent=bxml.beast, doc=bxml, addFinalizer=T)
	}
	else	
		dummy	<- lapply(seq_along(beast2.spec$alignment.filter), function(i)
				{
					dummy	<- newXMLNode("alignment", attrs= list(id=beast2.spec$alignment.id[i], filter=beast2.spec$alignment.filter[i], data=paste('@',beast2.spec$data.id,sep=''), spec=beast2.spec$alignment.spec), parent=bxml.beast, doc=bxml, addFinalizer=T)
					dummy	<- newXMLNode("taxa", attrs= list(id=paste("TaxonSet.t",beast2.spec$alignment.id[i],sep=':'), spec=beast2.spec$alignment.taxa.spec, alignment=paste('@',beast2.spec$alignment.id[i],sep='')), parent=bxml.beast, doc=bxml, addFinalizer=T)
				})		
	if(verbose)	cat(paste("\nadded new alignments and taxonsets, n=", length(beast2.spec$alignment.id)))	
	bxml
}
######################################################################################
beast2.add.datetrait<- function(bxml, df, beast2.spec, verbose=1)	
{			
	tmp			<- df[,BEASTlabel]	
	if(verbose)	cat(paste("\nsetting tip date to time at pos x in label, x=", beast2.spec$beast.label.datepos))
	tmp			<- sapply( strsplit(tmp, beast2.spec$beast.label.sep, fixed=1), function(x) paste(paste(x,collapse=beast2.spec$beast.label.sep,sep=''),'=',x[beast2.spec$beast.label.datepos],sep='') )
	
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	dummy		<- newXMLCommentNode(text="The tip dates.", parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text=paste("ntax=",nrow(df),sep=''), parent=bxml.beast, doc=bxml, addFinalizer=T)	
	dummy		<- newXMLNode("trait", attrs= list(id=beast2.spec$datetrait.id, 
					spec= beast2.spec$datetrait.spec,
					units= beast2.spec$datetrait.units,
					taxa= paste('@',beast2.spec$datetrait.taxa,sep=''),
					traitname= beast2.spec$datetrait.traitname, 
					value=paste(tmp, collapse=", ",sep='')), parent=bxml.beast, doc=bxml, addFinalizer=T)
	if(verbose)	cat(paste("\nadded new date trait for sequences, n=", length(tmp)))
	#define taxonset in here
	bxml
}
######################################################################################
beast2.add.satree<- function(bxml, beast2.spec, verbose=verbose)
{
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	if(1 || is.na(beast2.spec$starttree.newick))
	{
		if(verbose)	cat(paste('\nadd SA tree, initialized with UPGMA clustering'))
		dummy		<- newXMLNode("tree", attrs= list(	id=beast2.spec$tree.id, 													 
						trait=paste('@',beast2.spec$datetrait.id,sep=''),
						nodetype=beast2.spec$sasky.tree.nodetype,
						clusterType=beast2.spec$sasky.tree.cluster.type,
						spec=beast2.spec$sasky.tree.cluster.spec,
						taxa=paste('@',beast2.spec$tree.taxonset,sep='')), parent=bxml.beast, doc=bxml, addFinalizer=T)			
	}
	else
	{
		if(verbose)	cat(paste('\nadd SA tree, initialized with newick tree'))
		bxml.tree	<- newXMLNode("tree", attrs= list(	id=beast2.spec$tree.id, 													 
						trait=paste('@',beast2.spec$datetrait.id,sep=''),
						nodetype=beast2.spec$sasky.tree.nodetype,						
						spec=beast2.spec$sasky.tree.parser.spec,
						taxa=paste('@',beast2.spec$tree.taxonset,sep='')), parent=bxml.beast, doc=bxml, addFinalizer=T)
		dummy		<- newXMLCommentNode(text=paste("start: Newick starting tree"), parent=bxml.tree, doc=bxml, addFinalizer=T)				
		tmp			<- newXMLNode("input", attrs= list(name='newick'), parent=bxml.tree, doc=bxml, addFinalizer=T)
		dummy		<- newXMLTextNode(text=beast2.spec$starttree.newick, parent=tmp, doc=bxml, addFinalizer=T)
		dummy		<- newXMLCommentNode(text=paste("end: Newick starting tree"), parent=bxml.tree, doc=bxml, addFinalizer=T)		
	}
	if(verbose)	cat(paste("\nadded SA trees for taxonsets, n=", length(beast2.spec$tree.taxonset)))
	bxml
}
######################################################################################
beast2.add.tree<- function(bxml, beast2.spec, verbose=verbose)
{
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	dummy		<- newXMLNode("tree", attrs= list(	id=beast2.spec$tree.id, 													 
					trait=paste('@',beast2.spec$datetrait.id,sep=''),
					taxonset=paste('@TaxonSet.t:',beast2.spec$tree.taxonset,sep='')), 
			parent=bxml.beast, doc=bxml, addFinalizer=T)	
	if(verbose)	cat(paste("\nadded trees for taxonsets, n=", length(beast2.spec$tree.taxonset)))
	bxml
}
######################################################################################
beast2.add.treemodel.bdsky<- function(bxml, beast2.spec, verbose=1)
{
	bxml.beast		<- getNodeSet(bxml, "//beast")[[1]]
	beast.treemodel	<- newXMLNode("BirthDeathSkylineModel", attrs= list(	id=beast2.spec$treemodel.id, 
					name=beast2.spec$treemodel.id, 
					tree=paste('@',beast2.spec$tree.id,sep=''),
					spec=beast2.spec$bdsky.spec,
					intervalNumber=as.character(beast2.spec$bdsky.intervalNumber)), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.sprop.id, 
					name="samplingProportion",
					dimension=beast2.spec$bdsky.intervalNumber,
					lower=as.character(beast2.spec$bdsky.sprop.lower),
					upper=as.character(beast2.spec$bdsky.sprop.upper),
					value=paste(beast2.spec$bdsky.sprop.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.R0.id, 
					name="R0",
					dimension=beast2.spec$bdsky.intervalNumber,
					lower=as.character(beast2.spec$bdsky.R0.lower),
					upper=as.character(beast2.spec$bdsky.R0.upper),
					value=paste(beast2.spec$bdsky.R0.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)											
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.notInf.id, 
					name="becomeUninfectiousRate",
					dimension=beast2.spec$bdsky.intervalNumber,
					lower=as.character(beast2.spec$bdsky.notInf.lower),
					upper=as.character(beast2.spec$bdsky.notInf.upper),
					value=paste(beast2.spec$bdsky.notInf.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.origin.id, 
					name="origin",
					lower=as.character(beast2.spec$bdsky.origin.lower),
					upper=as.character(beast2.spec$bdsky.origin.upper),
					value=paste(beast2.spec$bdsky.origin.value)), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	tmp				<- rep("false",4)
	tmp[which(!sapply(c(beast2.spec$bdsky.R0.changepoint.id[1],beast2.spec$bdsky.notInf.changepoint.id[1],beast2.spec$bdsky.sprop.changepoint.id[1]),is.null))]		<- "true"
	dummy			<- newXMLNode("reverseTimeArrays", attrs= list(	id=beast2.spec$bdsky.reverseTimeArrays.id, 
					spec=beast2.spec$bdsky.reverseTimeArrays.spec,
					value=paste(tmp, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	if(!is.null(beast2.spec$bdsky.R0.changepoint.id[1]))												
		dummy		<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.R0.changepoint.id, 
						name="birthRateChangeTimes",
						value=paste(beast2.spec$bdsky.R0.changepoint.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	
	if(!is.null(beast2.spec$bdsky.notInf.changepoint.id[1]))												
		dummy		<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.notInf.changepoint.id, 
						name="deathRateChangeTimes",
						value=paste(beast2.spec$bdsky.notInf.changepoint.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	if(!is.null(beast2.spec$bdsky.sprop.changepoint.id[1]))												
		dummy		<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.sprop.changepoint.id, 
						name="samplingRateChangeTimes",
						value=paste(beast2.spec$bdsky.sprop.changepoint.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	if(verbose)	cat(paste("\nadded BDSKY tree models for taxonsets, n=", length(beast2.spec$tree.taxonset)))
	bxml	
}
######################################################################################
beast2.add.treemodel.sasky<- function(bxml, beast2.spec, verbose=1)
{
	bxml.beast		<- getNodeSet(bxml, "//beast")[[1]]
	beast.treemodel	<- newXMLNode("BirthDeathSkylineModel", attrs= list(	id=beast2.spec$treemodel.id, 
					name=beast2.spec$treemodel.id, 
					tree=paste('@',beast2.spec$tree.id,sep=''),
					spec=beast2.spec$sasky.spec,
					intervalNumber=as.character(beast2.spec$bdsky.intervalNumber)), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.sprop.id, 
					name="samplingProportion",
					dimension=beast2.spec$bdsky.intervalNumber,
					lower=as.character(beast2.spec$bdsky.sprop.lower),
					upper=as.character(beast2.spec$bdsky.sprop.upper),
					value=paste(beast2.spec$bdsky.sprop.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.R0.id, 
					name="R0",
					dimension=beast2.spec$bdsky.intervalNumber,
					lower=as.character(beast2.spec$bdsky.R0.lower),
					upper=as.character(beast2.spec$bdsky.R0.upper),
					value=paste(beast2.spec$bdsky.R0.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)											
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.notInf.id, 
					name="becomeUninfectiousRate",
					dimension=beast2.spec$bdsky.intervalNumber,
					lower=as.character(beast2.spec$bdsky.notInf.lower),
					upper=as.character(beast2.spec$bdsky.notInf.upper),
					value=paste(beast2.spec$bdsky.notInf.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$sasky.r.id, 
					name="removalProbability",
					dimension=beast2.spec$bdsky.intervalNumber,
					lower=as.character(beast2.spec$sasky.r.lower),
					upper=as.character(beast2.spec$sasky.r.upper),
					value=paste(beast2.spec$sasky.r.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)	
	dummy			<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.origin.id, 
					name="origin",
					lower=as.character(beast2.spec$bdsky.origin.lower),
					upper=as.character(beast2.spec$bdsky.origin.upper),
					value=paste(beast2.spec$bdsky.origin.value)), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	if(!is.null(beast2.spec$bdsky.reverseTimeArrays.id[1]))
		dummy		<- newXMLNode("reverseTimeArrays", attrs= list(	id=beast2.spec$bdsky.reverseTimeArrays.id, 
					spec=beast2.spec$bdsky.reverseTimeArrays.spec,
					value=paste(beast2.spec$bdsky.reverseTimeArrays.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	if(!is.null(beast2.spec$bdsky.R0.changepoint.id[1]))												
		dummy		<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.R0.changepoint.id, 
						name="birthRateChangeTimes",
						value=paste(beast2.spec$bdsky.R0.changepoint.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)	
	if(!is.null(beast2.spec$bdsky.notInf.changepoint.id[1]))												
		dummy		<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.notInf.changepoint.id, 
						name="deathRateChangeTimes",
						value=paste(beast2.spec$bdsky.notInf.changepoint.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	if(!is.null(beast2.spec$bdsky.sprop.changepoint.id[1]))												
		dummy		<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.sprop.changepoint.id, 
						name="samplingRateChangeTimes",
						value=paste(beast2.spec$bdsky.sprop.changepoint.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)
	if(!is.null(beast2.spec$bdsky.r.changepoint.id[1]))												
		dummy		<- newXMLNode("parameter", attrs= list(	id=beast2.spec$bdsky.r.changepoint.id, 
						name="removalProbabilityChangeTimes",
						value=paste(beast2.spec$bdsky.r.changepoint.value, collapse=' ')), parent=beast.treemodel, doc=bxml, addFinalizer=T)	
	if(verbose)	cat(paste("\nadded SASKY tree models for taxonsets, n=", length(beast2.spec$tree.taxonset)))
	bxml	
}
######################################################################################
beast2.init.xml<- function( beast2.spec=NULL, verbose=1)
{
	bxml		<- newXMLDoc(addFinalizer=T)
	bxml.beast	<- newXMLNode("beast", attrs=list(beautitemplate='HIVCLUST', beautistatus='', version="2.0", namespace=paste(beast2.spec$namespace,sep='',collapse=':')), doc=bxml, addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="Beta"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.Beta, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="ExcludablePrior"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.ExcludablePrior, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="Exponential"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.Exponential, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="InverseGamma"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.InverseGamma, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="LogNormal"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.LogNormal, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="Gamma"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.Gamma, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="Uniform"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.Uniform, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="prior"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.prior, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="LaplaceDistribution"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.Laplace, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="OneOnX"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.OneOnX, parent=dummy, doc=bxml,addFinalizer=T)
	dummy		<- newXMLNode("map", attrs= list(name="Normal"), parent=bxml.beast, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(beast2.spec$map.Normal, parent=dummy, doc=bxml,addFinalizer=T)
	bxml
}	
######################################################################################
beast2.add.tiplogstem<- function(bxml, beast2.spec, verbose=0)
{
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	taxons		<- xpathApply(bxml.beast, "//sequence", xmlGetAttr, "taxon")
	if(verbose)	cat(paste("\nadd taxonsets for tips, n=", length(taxons)))
	#	add taxonsets	
	tmp			<- newXMLCommentNode(text="Tip Taxonsets, used to log the tip stem height", doc=bxml, addFinalizer=T)
	tmp			<- c(tmp, lapply(seq_along(taxons), function(i)
					{				
						taxonset	<- newXMLNode("taxonset", attrs= list(	id=paste(beast2.spec$tip.taxonset.id.prefix,i,sep=''), spec=beast2.spec$taxonset.spec ), doc=bxml, addFinalizer=T)
						if(length(getNodeSet(bxml, paste("//taxon[@id='",taxons[[i]],"']",sep=''))))
							dummy	<- newXMLNode("taxon", attrs= list(	idref=taxons[[i]], spec=beast2.spec$taxon.spec ), parent=taxonset, doc=bxml, addFinalizer=T)
						else
							dummy	<- newXMLNode("taxon", attrs= list(	id=taxons[[i]], spec=beast2.spec$taxon.spec ), parent=taxonset, doc=bxml, addFinalizer=T)
						taxonset
					}))	
	dummy		<- addChildren(bxml.beast, tmp, at= tail(which(xmlSApply(bxml.beast, xmlName)=="data"),1) )
	#	add MRCA priors	
	if(verbose)	cat(paste("\nadd MRCA priors for tips, n=", length(taxons)))
	bxml.prior	<- getNodeSet(bxml, "//*[@id='prior']")
	if(length(bxml.prior)!=1)	stop("unexpected length of bxml.prior")
	bxml.prior	<- bxml.prior[[1]]	
	dummy		<- newXMLCommentNode(text="start: Tip Taxonset priors, used to log the tip stem height", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_along(taxons), function(i)
			{
				newXMLNode("distribution", attrs= list(	id=paste(beast2.spec$tip.prior.id.prefix,i,sep=''), taxonset=paste('@',beast2.spec$tip.taxonset.id.prefix,i,sep=''), useOriginate='true', monophyletic='false', tree=paste('@',beast2.spec$tree.id,sep=''), spec=beast2.spec$mrca.prior.spec ), parent=bxml.prior, doc=bxml, addFinalizer=T)
			})
	dummy		<- newXMLCommentNode(text="end: Tip Taxonset priors, used to log the tip stem height", parent=bxml.prior, doc=bxml, addFinalizer=T)
	if(verbose)	cat(paste("\nadd tracelog for stem of tips, n=", length(taxons)))
	#	add log entry to tracelog
	bxml.log	<- getNodeSet(bxml, "//logger[@id='tracelog']")
	if(length(bxml.log)!=1)	stop("unexpected length of bxml.log")
	bxml.log	<- bxml.log[[1]]
	dummy		<- newXMLCommentNode(text="start: Tip Taxonset log entries, used to log the tip stem height", parent=bxml.log, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_along(taxons), function(i)
			{
				newXMLNode("log", attrs= list(	idref=paste(beast2.spec$tip.prior.id.prefix,i,sep='')	), parent=bxml.log, doc=bxml, addFinalizer=T)
			})
	dummy		<- newXMLCommentNode(text="end: Tip Taxonset log entries, used to log the tip stem height", parent=bxml.log, doc=bxml, addFinalizer=T)
	bxml.beast
}
######################################################################################
beast2.add.startingtree.random<- function(bxml, beast2.spec, verbose=0)
{
	bxml.run	<- getNodeSet(bxml, "//run")[[1]]
	if(length(bxml.run)!=1)	stop("unexpected length of bxml.run")
	if(verbose)	cat(paste('\nadd random starting tree'))
	dummy		<- newXMLCommentNode(text=paste("Random starting tree"), parent=bxml.run, doc=bxml, addFinalizer=T)
	tmp			<- newXMLNode("init", attrs= list(	id=beast2.spec$starttree.id, estimate='false', initial=beast2.spec$tree.id, 
					taxa=paste('@',beast2.spec$rndstarttree.taxonset,sep=''), taxonset=paste('@TaxonSet.t:',beast2.spec$rndstarttree.taxonset,sep='')), parent=bxml.run, doc=bxml, addFinalizer=T)		
	tmp			<- newXMLNode("populationModel", attrs= list(	id=paste('ConstantPopulation',beast2.spec$tree.id,sep=''), spec='ConstantPopulation'	), parent=tmp, doc=bxml, addFinalizer=T)
	tmp			<- newXMLNode("parameter", attrs= list(	id=paste('PopSize',beast2.spec$tree.id,sep=''), name='popSize', value='1.0'), parent=tmp, doc=bxml, addFinalizer=T)
	bxml
}
######################################################################################
beast2.add.startingtree.newick<- function(bxml, beast2.spec, verbose=0)
{
	bxml.run	<- getNodeSet(bxml, "//run")[[1]]
	if(length(bxml.run)!=1)	stop("unexpected length of bxml.run")
	if(verbose)	cat(paste('\nadd newick starting tree'))
	dummy		<- newXMLCommentNode(text=paste("start: Newick starting tree"), parent=bxml.run, doc=bxml, addFinalizer=T)
	tmp			<- newXMLNode("init", attrs= list(	id=beast2.spec$starttree.id, initial=beast2.spec$tree.id, IsLabelledNewick=beast2.spec$starttree.islabelledtree, spec=beast2.spec$starttree.spec), parent=bxml.run, doc=bxml, addFinalizer=T)		
	tmp			<- newXMLNode("input", attrs= list(name='newick'), parent=tmp, doc=bxml, addFinalizer=T)
	dummy		<- newXMLTextNode(text=beast2.spec$starttree.newick, parent=tmp, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text=paste("end: Newick starting tree"), parent=bxml.run, doc=bxml, addFinalizer=T)
	bxml
}
######################################################################################
beast2.add.sasky.serialpriors<- function(bxml, beast2.spec, verbose=1)
{
	bxml.prior	<- getNodeSet(bxml, "//*[@id='prior']")
	if(length(bxml.prior)!=1)	stop("unexpected length of bxml.prior")
	bxml.prior	<- bxml.prior[[1]]	
	dummy		<- newXMLCommentNode(text="start: serial SASKY priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text="The SASKY model prior.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	tmp			<- newXMLNode("distribution", attrs= list(	id=paste("prior",beast2.spec$treemodel.id,sep='-'),spec=beast2.spec$compoundprior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)	
	dummy		<- newXMLNode("distribution", attrs= list(	idref=beast2.spec$treemodel.id), parent=tmp, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text="The origin prior.", parent=bxml.prior, doc=bxml, addFinalizer=T)	
	tmp			<- newXMLNode("prior", attrs= list(	id=paste("sprior",beast2.spec$bdsky.origin.id,sep='-'),
					name="distribution",
					x= paste('@',beast2.spec$bdsky.origin.id,sep='')), parent=bxml.prior, doc=bxml, addFinalizer=T)	
	dummy		<- beast2.get.prior(beast2.spec$bdsky.origin.prior, tmp, xmlAttrs(tmp)["id"], bxml)
	dummy		<- newXMLCommentNode(text="The serial samplingProb priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_len(beast2.spec$bdsky.intervalNumber), function(serial.id)
			{
				tmp				<- rep(0,beast2.spec$bdsky.intervalNumber)
				tmp[serial.id]	<- 1
				serialwrapper	<- newXMLNode("distribution", attrs= list(	id=paste("sprior",serial.id,beast2.spec$bdsky.sprop.id,sep='-'),
								name="distribution",
								x= paste('@',beast2.spec$bdsky.sprop.id,sep=''),
								xInclude=paste(tmp, collapse=' ',sep=''), 					
								spec= beast2.spec$bdsky.prior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)
				dummy			<- beast2.get.prior(beast2.spec$bdsky.sprop.prior[serial.id], serialwrapper, xmlAttrs(serialwrapper)["id"], bxml)				
			})
	if(verbose)	cat(paste("\nadded serial samplingProb priors, n=", length(dummy)))
	dummy		<- newXMLCommentNode(text="The serial R0 priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_len(beast2.spec$bdsky.intervalNumber), function(serial.id)
			{
				tmp				<- rep(0,beast2.spec$bdsky.intervalNumber)
				tmp[serial.id]	<- 1
				serialwrapper	<- newXMLNode("distribution", attrs= list(	id=paste("sprior",serial.id,beast2.spec$bdsky.R0.id,sep='-'),
								name="distribution",
								x= paste('@',beast2.spec$bdsky.R0.id,sep=''),
								xInclude=paste(tmp, collapse=' ',sep=''), 					
								spec= beast2.spec$bdsky.prior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)
				dummy			<- beast2.get.prior(beast2.spec$bdsky.R0.prior[serial.id], serialwrapper, xmlAttrs(serialwrapper)["id"], bxml)				
			})												
	if(verbose)	cat(paste("\nadded serial R0 priors, n=", length(dummy)))
	dummy		<- newXMLCommentNode(text="The serial becomingUninfectious priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_len(beast2.spec$bdsky.intervalNumber), function(serial.id)
			{
				tmp				<- rep(0,beast2.spec$bdsky.intervalNumber)
				tmp[serial.id]	<- 1
				serialwrapper	<- newXMLNode("distribution", attrs= list(	id=paste("sprior",serial.id,beast2.spec$bdsky.notInf.id,sep='-'),
								name="distribution",
								x= paste('@',beast2.spec$bdsky.notInf.id,sep=''),
								xInclude=paste(tmp, collapse=' ',sep=''), 					
								spec= beast2.spec$bdsky.prior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)
				dummy			<- beast2.get.prior(beast2.spec$bdsky.notInf.prior[serial.id], serialwrapper, xmlAttrs(serialwrapper)["id"], bxml)				
			})	
	if(verbose)	cat(paste("\nadded serial becomingUninfectious priors, n=", length(dummy)))
	dummy		<- newXMLCommentNode(text="The serial becomeNoninfectiousAfterSamplingProbability priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_len(beast2.spec$bdsky.intervalNumber), function(serial.id)
			{
				tmp				<- rep(0,beast2.spec$bdsky.intervalNumber)
				tmp[serial.id]	<- 1
				serialwrapper	<- newXMLNode("distribution", attrs= list(	id=paste("sprior",serial.id,beast2.spec$sasky.r.id,sep='-'),
								name="distribution",
								x= paste('@',beast2.spec$sasky.r.id,sep=''),
								xInclude=paste(tmp, collapse=' ',sep=''), 					
								spec= beast2.spec$bdsky.prior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)
				dummy			<- beast2.get.prior(beast2.spec$sasky.r.prior[serial.id], serialwrapper, xmlAttrs(serialwrapper)["id"], bxml)				
			})	
	if(verbose)	cat(paste("\nadded serial becomeNoninfectiousAfterSamplingProbability priors, n=", length(dummy)))	
	dummy		<- newXMLCommentNode(text="end: serial SASKY priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	bxml.prior
}
######################################################################################
beast2.get.startingtree<- function(ph, df, beast2.spec, verbose=1)
{
	require(adephylo)
	if(verbose) cat(paste("\ncreate startingTree with root height=",beast2.spec$starttree.rootHeight))
	tmp					<- match( setdiff( ph$tip.label, df[,FASTASampleCode] ), ph$tip.label)
	ph.start			<- drop.tip(ph, tmp)		
	ph.start$node.label	<- NULL
	setkey(df, FASTASampleCode)
	ph.start$tip.label	<- df[ph.start$tip.label,][,BEASTlabel]
	if(verbose) cat(paste("\nselected tips for startingTree, n=",Ntip(ph.start)))
	#	adjust rootHeight to 'beast.rootHeight'	
	ph.start$edge.length<- ph.start$edge.length	*	beast2.spec$starttree.rootHeight / max(distRoot(ph.start))
	if(verbose) cat(paste("\nadjusted root height=",max(distRoot(ph.start))))	
	if(beast2.spec$starttree.rootHeight>=beast2.spec$bdsky.origin.value)	stop('Detected invalid rootHeight. Must be smaller than origin.value.')
	write.tree( ph.start )		
}
######################################################################################
beast2.get.specifications	<- function(xml.dir=NA, xml.filename=NA, mcmc.length=20e6, bdsky.intervalNumber=4, alignment.filter=NA, tip.log.stem=FALSE, cluster.log=FALSE, cluster.monophyletic=FALSE)
{
	beast2.spec<- list()	
	beast2.spec$namespace						<- c("beast.core","beast.evolution.alignment","beast.evolution.tree.coalescent","beast.core.util","beast.evolution.nuc","beast.evolution.operators","beast.evolution.sitemodel","beast.evolution.substitutionmodel","beast.evolution.likelihood","beast.evolution.speciation","beast.core.parameter")
	beast2.spec$xml.dir							<- xml.dir		
	beast2.spec$xml.filename					<- xml.filename
	beast2.spec$pool.cnts.requested				<- rep(NA, bdsky.intervalNumber)
	beast2.spec$pool.ntip						<- 130
	beast2.spec$pool.fNegT						<- 0.8
	beast2.spec$map.Beta						<- "beast.math.distributions.Beta"
	beast2.spec$map.Exponential					<- "beast.math.distributions.Exponential"
	beast2.spec$map.ExcludablePrior				<- "beast.math.distributions.ExcludablePrior"
	beast2.spec$map.InverseGamma				<- "beast.math.distributions.InverseGamma"
	beast2.spec$map.LogNormal					<- "beast.math.distributions.LogNormalDistributionModel"
	beast2.spec$map.Gamma						<- "beast.math.distributions.Gamma"
	beast2.spec$map.Uniform						<- "beast.math.distributions.Uniform"
	beast2.spec$map.prior						<- "beast.math.distributions.Prior"
	beast2.spec$map.Laplace						<- "beast.math.distributions.LaplaceDistribution"
	beast2.spec$map.OneOnX						<- "beast.math.distributions.OneOnX"
	beast2.spec$map.Normal						<- "beast.math.distributions.Normal"
	beast2.spec$mcmc.length						<- mcmc.length
	beast2.spec$beast.label.datepos				<- 4
	beast2.spec$beast.label.sep					<- '_'		
	beast2.spec$data.missing					<- "-?"
	beast2.spec$data.dataType					<- 'nucleotide'
	beast2.spec$alignment.spec					<- "FilteredAlignment"
	beast2.spec$alignment.taxa.spec				<- "TaxonSet"
	beast2.spec$alignment.filter				<- alignment.filter
	beast2.spec$taxon.spec						<- "Taxon"
	beast2.spec$taxonset.spec					<- "TaxonSet"	
	beast2.spec$mrca.prior.spec					<- "beast.math.distributions.MRCAPrior"
	if(length(alignment.filter)==1 && is.na(alignment.filter))
	{
		beast2.spec$data.id						<- 'ds'
		beast2.spec$alignment.id				<- 'ds'
	}
	else
	{
		beast2.spec$data.id						<- 'data'
		beast2.spec$alignment.id				<- paste('ds', seq_along(beast2.spec$alignment.filter), sep='_')
	}
	beast2.spec$tree.taxonset					<- beast2.spec$alignment.id[1]
	beast2.spec$tree.id							<- paste('Tree',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$sequence.totalcount				<- 4
	beast2.spec$datetrait.spec					<- "beast.evolution.tree.TraitSet"
	beast2.spec$datetrait.taxa.spec				<- "TaxonSet"
	beast2.spec$datetrait.taxa					<- paste("TaxonSet.t",beast2.spec$alignment.id[1],sep=':')
	beast2.spec$datetrait.id					<- paste("dateTrait.t",beast2.spec$alignment.id[1],sep=':')
	beast2.spec$datetrait.units					<- "year"
	beast2.spec$datetrait.traitname				<- "date"
	beast2.spec$treemodel						<- "BirthDeathSkylineModel"
	beast2.spec$treemodel.id					<- paste("birthDeath",beast2.spec$tree.taxonset,sep='.t:')	
	beast2.spec$bdsky.spec						<- "beast.evolution.speciation.BirthDeathSkylineModel"
	beast2.spec$bdsky.prior.spec				<- "beast.math.distributions.ExcludablePrior"
	beast2.spec$bdsky.intervalNumber			<- bdsky.intervalNumber
	beast2.spec$bdsky.origin.id					<- paste('originS',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$bdsky.origin.value				<- rep(38, length(beast2.spec$bdsky.origin.id))
	beast2.spec$bdsky.origin.lower				<- 0.0
	beast2.spec$bdsky.origin.upper				<- 1000.0
	beast2.spec$bdsky.origin.prior				<- "Uniform/20.0/40.0"
	beast2.spec$bdsky.sprop.id					<- paste('samplingProportionS',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$bdsky.sprop.value				<- rep(0.4, beast2.spec$bdsky.intervalNumber)
	beast2.spec$bdsky.sprop.lower				<- 0.0
	beast2.spec$bdsky.sprop.upper				<- 1.0
	beast2.spec$bdsky.sprop.prior				<- rep("Uniform/0.2/1.0",beast2.spec$bdsky.intervalNumber)
	beast2.spec$bdsky.sprop.changepoint.id		<- paste('samplingRateChangeTimes',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$bdsky.sprop.changepoint.value	<- c(1.596, 5.596, 9.596, 0.)	
	beast2.spec$bdsky.R0.id						<- paste('R0S',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$bdsky.R0.value					<- rep(1.2, beast2.spec$bdsky.intervalNumber)
	beast2.spec$bdsky.R0.lower					<- 0.0
	beast2.spec$bdsky.R0.upper					<- 10.0
	beast2.spec$bdsky.R0.prior					<- rep("Gamma/1.5/1.5/0",beast2.spec$bdsky.intervalNumber)
	beast2.spec$bdsky.R0.changepoint.id			<- paste('birthRateChangeTimes',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$bdsky.R0.changepoint.value		<- c(1.596, 5.596, 9.596, 0.)	
	beast2.spec$bdsky.notInf.id					<- paste('becomeUninfectiousRateS',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$bdsky.notInf.value				<- rep(0.1, beast2.spec$bdsky.intervalNumber)
	beast2.spec$bdsky.notInf.lower				<- 0.0
	beast2.spec$bdsky.notInf.upper				<- 10.0
	beast2.spec$bdsky.notInf.prior				<- rep("OneOnX/0",beast2.spec$bdsky.intervalNumber)
	beast2.spec$bdsky.notInf.changepoint.id		<- paste('deathRateChangeTimes',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$bdsky.notInf.changepoint.value	<- c(1.596, 5.596, 9.596, 0.)
	beast2.spec$bdsky.reverseTimeArrays.id		<- paste('reverseTimeArrays',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$bdsky.reverseTimeArrays.value	<- c('true', 'true', 'false', 'true')
	beast2.spec$bdsky.reverseTimeArrays.spec	<- "parameter.BooleanParameter"	
	beast2.spec$sasky.spec						<- "beast.evolution.speciation.SABDSkylineModel"		
	beast2.spec$sasky.tree.nodetype				<- "beast.evolution.tree.ZeroBranchSANode"
	beast2.spec$sasky.tree.parser.spec			<- "beast.util.ZeroBranchSATreeParser"
	beast2.spec$sasky.tree.cluster.spec			<- "beast.util.ClusterZBSATree"
	beast2.spec$sasky.tree.cluster.type			<- "upgma"
	beast2.spec$sasky.r.id						<- paste('r',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$sasky.r.value					<- rep(0.5, beast2.spec$bdsky.intervalNumber)
	beast2.spec$sasky.r.lower					<- 0.0
	beast2.spec$sasky.r.upper					<- 1.0
	beast2.spec$sasky.r.prior					<- rep("Uniform/0.0/1.0",beast2.spec$bdsky.intervalNumber)
	beast2.spec$sasky.r.changepoint.id			<- paste('removalProbabilityChangeTimes',beast2.spec$tree.taxonset,sep='.t:')
	beast2.spec$sasky.r.changepoint.value		<- c(1.596, 5.596, 9.596, 0.)		
	beast2.spec$compoundprior.spec				<- "util.CompoundDistribution"
	beast2.spec$tip.log.stem					<- tip.log.stem
	beast2.spec$tip.taxonset.id.prefix			<- 'tip'
	beast2.spec$tip.prior.id.prefix				<- 'tipprior'
	beast2.spec$cluster.taxonset.id.prefix		<- 'c'
	beast2.spec$cluster.prior.id.prefix			<- 'cprior'
	beast2.spec$cluster.monophyletic			<- ifelse(cluster.monophyletic, 'true', 'false')
	beast2.spec$cluster.log						<- cluster.log	
	beast2.spec$starttree.rootHeight			<- 35
	beast2.spec$starttree.usingDates			<- 'true'
	beast2.spec$starttree.brlunits				<- 'years'
	beast2.spec$starttree.islabelledtree		<- 'true'
	beast2.spec$starttree.id					<- paste('Starting',beast2.spec$tree.id,sep='')	
	beast2.spec$starttree.newick				<- NA
	beast2.spec$starttree.spec					<- 'beast.util.TreeParser'	
	beast2.spec$rndstarttree.taxonset			<- beast2.spec$tree.taxonset
	beast2.spec$rndstarttree.spec				<- 'beast.evolution.tree.RandomTree'
	beast2.spec
} 
######################################################################################
beast2.get.xml<- function(	bxml.template, seq.PROT.RT, df, beast2.spec, ph=NULL, verbose=1)
{	
	require(XML)
	#	init XML
	bxml		<- beast2.init.xml( beast2.spec, verbose=verbose)
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]	
	#	add data	
	dummy		<- beast2.add.data(bxml, seq.PROT.RT, df, beast2.spec, verbose=verbose)
	#	add alignment, alignment filters, alignment taxonsets
	dummy		<- beast2.add.alignment(bxml, beast2.spec, verbose=verbose)
	#	add tip dates
	dummy		<- beast2.add.datetrait(bxml, df, beast2.spec, verbose=verbose)
	#	add tree and starting tree for alignment
	if(beast2.spec$treemodel=="BirthDeathSkylineModel")
		dummy	<- beast2.add.tree(bxml, beast2.spec, verbose=verbose)
	if(beast2.spec$treemodel=="SampledAncestorSkylineModel")
		dummy	<- beast2.add.satree(bxml, beast2.spec, verbose=verbose)
	#	add tree model -- TODO move starting tree into here
	if(beast2.spec$treemodel=="BirthDeathSkylineModel")
		dummy	<- beast2.add.treemodel.bdsky(bxml, beast2.spec, verbose=verbose)
	#	add tree model and starting tree
	if(beast2.spec$treemodel=="SampledAncestorSkylineModel")
		dummy	<- beast2.add.treemodel.sasky(bxml, beast2.spec, verbose=verbose)	
	#	copy branchRateModel from template
	tmp			<- getNodeSet(bxml.template, "//branchRateModel")
	if(length(tmp)!=1) stop("unexpected number of //branchRateModel")
	dummy<- addChildren( bxml.beast, xmlClone( tmp[[1]], addFinalizer=T, doc=bxml ) )
	if(verbose)	cat(paste("\nadded branchRateModel from template, size=", xmlSize(tmp[[1]])))	
	#	copy siteModel from template
	tmp			<- getNodeSet(bxml.template, "//siteModel")
	if(length(tmp)<1) stop("not found //siteModel")
	dummy		<- lapply(seq_along(tmp),function(i)	addChildren( bxml.beast, xmlClone( tmp[[i]], addFinalizer=T, doc=bxml ) )		) 		 
	if(verbose)	cat(paste("\nadded siteModel from template, numer=", length(tmp)))
	#	copy run from template
	tmp			<- getNodeSet(bxml.template, "//run")
	if(length(tmp)!=1) stop("unexpected number of //run")
	dummy		<- addChildren( bxml.beast, xmlClone( tmp[[1]], addFinalizer=T, doc=bxml ) )
	if(verbose)	cat(paste("\nadded run from template, size=", xmlSize(tmp[[1]])))
	#	add initial tree now if bdsky
	if(beast2.spec$treemodel=="BirthDeathSkylineModel" && is.na(beast2.spec$starttree.newick))
		dummy	<- beast2.add.startingtree.random(bxml, beast2.spec, verbose=verbose)
	if(beast2.spec$treemodel=="BirthDeathSkylineModel" && !is.na(beast2.spec$starttree.newick))
		dummy	<- beast2.add.startingtree.newick(bxml, beast2.spec, verbose=verbose)		
	#	add tree model prior
	if(beast2.spec$treemodel=="BirthDeathSkylineModel")
		dummy	<- beast2.add.bdsky.serialpriors(bxml, beast2.spec, verbose=verbose)
	if(beast2.spec$treemodel=="SampledAncestorSkylineModel")
		dummy	<- beast2.add.sasky.serialpriors(bxml, beast2.spec, verbose=verbose)
	
	#	reset output fileNames
	bxml.onodes	<- getNodeSet(bxml, "//*[@fileName]")
	tmp			<- sapply(bxml.onodes, function(x) xmlGetAttr(x,"fileName"))	
	tmp			<- sapply(strsplit(tmp,'.',fixed=1), function(x)	paste(beast2.spec$xml.filename, '.', rev(x)[1], sep=''))		
	dummy		<- sapply(seq_along(bxml.onodes), function(i){		xmlAttrs(bxml.onodes[[i]])["fileName"]<- tmp[i]		})
	if(verbose)	cat(paste("\nchanged trunk filename to", beast2.spec$xml.filename))
	#	reset chain length and logEvery
	bxml.onodes	<- getNodeSet(bxml, "//*[@chainLength]")
	dummy		<- sapply(seq_along(bxml.onodes), function(i){		xmlAttrs(bxml.onodes[[i]])["chainLength"]<- sprintf("%d",beast2.spec$mcmc.length)		})
	if(verbose)	cat(paste("\nchanged chain length to", beast2.spec$mcmc.length))
	bxml.onodes	<- getNodeSet(bxml, "//*[@logEvery]")
	dummy		<- sapply(seq_along(bxml.onodes), function(i){		xmlAttrs(bxml.onodes[[i]])["logEvery"]<- sprintf("%d",round(beast2.spec$mcmc.length/1e4))		})
	if(verbose)	cat(paste("\nchanged logEvery to", round(beast2.spec$mcmc.length/1e4)))	
	#	add logs for height of tip stem
	if(verbose)	cat(paste("\nadd log for height of tip stems=", beast2.spec$tip.log.stem))
	if(beast2.spec$tip.log.stem)
		dummy	<- beast2.add.tiplogstem(bxml, beast2.spec, verbose=verbose)
	#	add cluster taxonsets
	if(verbose)	cat(paste("\nadd taxonsets for clusters=", beast2.spec$cluster.monophyletic))
	if(beast2.spec$cluster.monophyletic)
		dummy	<- beast2.add.cluster.taxonsets(bxml, df, beast2.spec, verbose=verbose)
	
	bxml	
}
######################################################################################
beast2.add.cluster.taxonsets<- function(bxml, df, beast2.spec, verbose=0)
{
	bxml.beast	<- getNodeSet(bxml, "//beast")[[1]]
	clusters	<- df[,unique(cluster)]	
	if(verbose)	cat(paste("\nadd taxonsets for clusters, n=", length( clusters )))
	#	add taxonsets	
	tmp			<- newXMLCommentNode(text="start: The cluster taxonsets, used to log their TMRCA or to enforce monophyly", doc=bxml, addFinalizer=T)
	tmp			<- c(tmp, lapply(clusters, function(clu)
					{				
						taxonset	<- newXMLNode("taxonset", attrs= list(	id=paste(beast2.spec$cluster.taxonset.id.prefix,clu,sep=''), spec=beast2.spec$taxonset.spec ), doc=bxml, addFinalizer=T)
						dummy		<- sapply( subset(df, cluster==clu)[, BEASTlabel], function(x)
								{
									if(length(getNodeSet(bxml, paste("//taxon[@id='",x,"']",sep=''))))									
										newXMLNode("taxon", attrs= list(	idref=x, spec=beast2.spec$taxon.spec ), parent=taxonset, doc=bxml, addFinalizer=T)
									else
										newXMLNode("taxon", attrs= list(	id=x, spec=beast2.spec$taxon.spec ), parent=taxonset, doc=bxml, addFinalizer=T)
								})										
						taxonset
					}))	
	tmp			<- c(tmp, newXMLCommentNode(text="end: The cluster taxonsets, used to log their TMRCA or to enforce monophyly", doc=bxml, addFinalizer=T))
	dummy		<- addChildren(bxml.beast, tmp, at= tail(which(xmlSApply(bxml.beast, xmlName)=="data"),1) )
	#	add MRCA priors	
	if(verbose)	cat(paste("\nadd MRCA priors for clusters, n=", length(clusters)))
	bxml.prior	<- getNodeSet(bxml, "//*[@id='prior']")
	if(length(bxml.prior)!=1)	stop("unexpected length of bxml.prior")
	bxml.prior	<- bxml.prior[[1]]	
	dummy		<- newXMLCommentNode(text="start: Cluster taxonset priors, used to log their TMRCA or to enforce monophyly", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(clusters, function(clu)
			{
				newXMLNode("distribution", attrs= list(	id=paste(beast2.spec$cluster.prior.id.prefix,clu,sep=''), taxonset=paste('@',beast2.spec$cluster.taxonset.id.prefix,clu,sep=''), useOriginate='false', monophyletic=beast2.spec$cluster.monophyletic, tree=paste('@',beast2.spec$tree.id,sep=''), spec=beast2.spec$mrca.prior.spec ), parent=bxml.prior, doc=bxml, addFinalizer=T)
			})
	dummy		<- newXMLCommentNode(text="end: Cluster taxonset priors, used to log their TMRCA or to enforce monophyly", parent=bxml.prior, doc=bxml, addFinalizer=T)
	if(verbose)	cat(paste("\nadd tracelog for cluster tmrca=", beast2.spec$cluster.log))
	if(beast2.spec$cluster.log)
	{
		if(verbose)	cat(paste("\nadd tracelog for cluster tmrca, n=", length(clusters)))				
		#	add log entry to tracelog
		bxml.log	<- getNodeSet(bxml, "//logger[@id='tracelog']")
		if(length(bxml.log)!=1)	stop("unexpected length of bxml.log")
		bxml.log	<- bxml.log[[1]]
		dummy		<- newXMLCommentNode(text="start: Cluster taxonset log entries, used to log their TMRCA", parent=bxml.log, doc=bxml, addFinalizer=T)
		dummy		<- lapply(clusters, function(clu)
				{
					newXMLNode("log", attrs= list(	idref=paste(beast2.spec$cluster.prior.id.prefix,clu,sep='')	), parent=bxml.log, doc=bxml, addFinalizer=T)
				})
		dummy		<- newXMLCommentNode(text="end: Cluster taxonset log entries, used to log their TMRCA", parent=bxml.log, doc=bxml, addFinalizer=T)		
	}
	bxml.beast	
}
######################################################################################
beast2.get.prior<- function(args, parent, parentid, bxml)
{
	args	<- strsplit(args,'/')[[1]]
	if(args[1]=="Uniform")
		prior	<- newXMLNode("Uniform", attrs= list( id= paste('U',parentid,sep='-'), name="distr", lower=args[2], upper=args[3]), parent=parent, doc=bxml, addFinalizer=T)											
	else if(args[1]=="OneOnX")
		prior	<- newXMLNode("OneOnX", attrs= list( id= paste('OneOnX.',parentid,sep=''), name="distr"), parent=parent, doc=bxml, addFinalizer=T)
	else if(args[1]=="Exponential")
	{
		prior	<- newXMLNode("Exponential", attrs= list( id= paste('Exponential',parentid,sep='-'), name="distr", offset=args[3]), parent=parent, doc=bxml, addFinalizer=T)
		dummy	<- newXMLNode("parameter", attrs= list( id= paste('pExponential',parentid,sep='-'), name="mean", value=args[2], estimate="false"), parent=prior, doc=bxml, addFinalizer=T)
	}
	else if(args[1]=="Beta")
	{
		prior	<- newXMLNode("Beta", attrs= list( id= paste('Beta',parentid,sep='-'), name="distr", offset=args[4]), parent=parent, doc=bxml, addFinalizer=T)
		dummy	<- newXMLNode("parameter", attrs= list( id= paste('pBeta1',parentid,sep='-'), lower="0.0", upper="10.0", name="alpha", value=args[2], estimate="false"), parent=prior, doc=bxml, addFinalizer=T)
		dummy	<- newXMLNode("parameter", attrs= list( id= paste('pBeta2',parentid,sep='-'), lower="0.0", upper="10.0", name="beta", value=args[3], estimate="false"), parent=prior, doc=bxml, addFinalizer=T)
	}
	else if(args[1]=="Gamma")
	{
		prior	<- newXMLNode("Gamma", attrs= list( id= paste('Gamma',parentid,sep='-'), name="distr", offset=args[4]), parent=parent, doc=bxml, addFinalizer=T)
		dummy	<- newXMLNode("parameter", attrs= list( id= paste('pGamma1',parentid,sep='-'), name="alpha",  value=args[2], estimate="false"), parent=prior, doc=bxml, addFinalizer=T)
		dummy	<- newXMLNode("parameter", attrs= list( id= paste('pGamma2',parentid,sep='-'), name="beta",  value=args[3], estimate="false"), parent=prior, doc=bxml, addFinalizer=T)
	}
	else if(args[1]=="LogNormal")
	{
		prior	<- newXMLNode("LogNormal", attrs= list( id= paste('LogNormal',parentid,sep='-'), name="distr", offset=args[4], meanInRealSpace=args[5]), parent=parent, doc=bxml, addFinalizer=T)
		dummy	<- newXMLNode("parameter", attrs= list( id= paste('pLogNormal1',parentid,sep='-'), name="M",  value=args[2], estimate="false"), parent=prior, doc=bxml, addFinalizer=T)
		dummy	<- newXMLNode("parameter", attrs= list( id= paste('pLogNormal2',parentid,sep='-'), name="S",  value=args[3], estimate="false"), parent=prior, doc=bxml, addFinalizer=T)
	}
	else	stop("prior not implemented")	
	prior
}
######################################################################################
beast2.add.bdsky.serialpriors<- function(bxml, beast2.spec, verbose=1)
{
	bxml.prior	<- getNodeSet(bxml, "//*[@id='prior']")
	if(length(bxml.prior)!=1)	stop("unexpected length of bxml.prior")
	bxml.prior	<- bxml.prior[[1]]	
	dummy		<- newXMLCommentNode(text="start: serial BDSKY priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text="The BDSKY model prior.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	tmp			<- newXMLNode("distribution", attrs= list(	id=paste("prior",beast2.spec$treemodel.id,sep='-'),spec=beast2.spec$compoundprior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)	
	dummy		<- newXMLNode("distribution", attrs= list(	idref=beast2.spec$treemodel.id), parent=tmp, doc=bxml, addFinalizer=T)
	dummy		<- newXMLCommentNode(text="The origin prior.", parent=bxml.prior, doc=bxml, addFinalizer=T)	
	tmp			<- newXMLNode("prior", attrs= list(	id=paste("sprior",beast2.spec$bdsky.origin.id,sep='-'),
					name="distribution",
					x= paste('@',beast2.spec$bdsky.origin.id,sep='')), parent=bxml.prior, doc=bxml, addFinalizer=T)	
	dummy		<- beast2.get.prior(beast2.spec$bdsky.origin.prior, tmp, xmlAttrs(tmp)["id"], bxml)
	dummy		<- newXMLCommentNode(text="The serial samplingProb priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_len(beast2.spec$bdsky.intervalNumber), function(serial.id)
			{
				tmp				<- rep(0,beast2.spec$bdsky.intervalNumber)
				tmp[serial.id]	<- 1
				serialwrapper	<- newXMLNode("distribution", attrs= list(	id=paste("sprior",serial.id,beast2.spec$bdsky.sprop.id,sep='-'),
								name="distribution",
								x= paste('@',beast2.spec$bdsky.sprop.id,sep=''),
								xInclude=paste(tmp, collapse=' ',sep=''), 					
								spec= beast2.spec$bdsky.prior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)
				dummy			<- beast2.get.prior(beast2.spec$bdsky.sprop.prior[serial.id], serialwrapper, xmlAttrs(serialwrapper)["id"], bxml)				
			})
	if(verbose)	cat(paste("\nadded serial samplingProb priors, n=", length(dummy)))
	dummy		<- newXMLCommentNode(text="The serial R0 priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_len(beast2.spec$bdsky.intervalNumber), function(serial.id)
			{
				tmp				<- rep(0,beast2.spec$bdsky.intervalNumber)
				tmp[serial.id]	<- 1
				serialwrapper	<- newXMLNode("distribution", attrs= list(	id=paste("sprior",serial.id,beast2.spec$bdsky.R0.id,sep='-'),
								name="distribution",
								x= paste('@',beast2.spec$bdsky.R0.id,sep=''),
								xInclude=paste(tmp, collapse=' ',sep=''), 					
								spec= beast2.spec$bdsky.prior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)
				dummy			<- beast2.get.prior(beast2.spec$bdsky.R0.prior[serial.id], serialwrapper, xmlAttrs(serialwrapper)["id"], bxml)				
			})												
	if(verbose)	cat(paste("\nadded serial R0 priors, n=", length(dummy)))
	dummy		<- newXMLCommentNode(text="The serial becomingUninfectious priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	dummy		<- lapply(seq_len(beast2.spec$bdsky.intervalNumber), function(serial.id)
			{
				tmp				<- rep(0,beast2.spec$bdsky.intervalNumber)
				tmp[serial.id]	<- 1
				serialwrapper	<- newXMLNode("distribution", attrs= list(	id=paste("sprior",serial.id,beast2.spec$bdsky.notInf.id,sep='-'),
								name="distribution",
								x= paste('@',beast2.spec$bdsky.notInf.id,sep=''),
								xInclude=paste(tmp, collapse=' ',sep=''), 					
								spec= beast2.spec$bdsky.prior.spec), parent=bxml.prior, doc=bxml, addFinalizer=T)
				dummy			<- beast2.get.prior(beast2.spec$bdsky.notInf.prior[serial.id], serialwrapper, xmlAttrs(serialwrapper)["id"], bxml)				
			})	
	if(verbose)	cat(paste("\nadded serial becomingUninfectious priors, n=", length(dummy)))
	dummy		<- newXMLCommentNode(text="end: serial BDSKY priors.", parent=bxml.prior, doc=bxml, addFinalizer=T)
	bxml.prior
}
