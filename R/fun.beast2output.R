######################################################################################
#' No description yet
#' @author Oliver Ratmann
beast2out.combine.clu.trees<- function(indir, file.info, beastlabel.idx.clu=1, beastlabel.idx.hivn=2, beastlabel.idx.hivd=3, beastlabel.idx.hivs=4, beastlabel.idx.samplecode= 6, beastlabel.idx.rate= NA, method.nodectime='any', verbose=FALSE)
{
	#	collect consensus tree and further info for plotting
	tmp			<- lapply(seq_len(nrow(file.info)), function(i)
			{
				#	load dated cluster phylogenies
				file				<- paste(indir, file.info[i,file], sep='/')
				if(verbose) cat(paste('\nload file=',file,'i=',i))
				tmp					<- load(file)
				topo.map			<- mph.clu.dtopo[which.max(freq),]
				tmp					<- which( grepl( paste('mph.i=',topo.map[,mph.i],'_',sep=''), names(ph.consensus) ) )
				if(length(tmp)!=1)	stop('unexpected ph.consensus index')
				topo.map.ph			<- ph.consensus[[ tmp ]]
				topo.map.SA			<- subset(mph.SA.cnt, equal.to==topo.map[,mph.i])
				set(topo.map.SA, NULL, 'SA.freq', topo.map.SA[,SA.freq/n])
				topo.map.nodectime	<- subset(mph.node.ctime, equal.to==topo.map[,mph.i])
				topo.map			<- merge(subset(topo.map, select=c(cluster, dens)), subset(topo.map.SA, select=c(cluster, tip, SA.freq)), by='cluster')
				topo.map.nodectime	<- subset(topo.map.nodectime, select=c(cluster, node, q, cdf, pdf))
				topo.any.nodectime	<- subset(mph.mapnode.pctime, select=c(cluster, node, q, cdf, pdf))
				if(method.nodectime=='any')
					node.ctime		<- topo.any.nodectime
				else if(method.nodectime=='map')
					node.ctime		<- topo.map.nodectime
				else	stop('unknown method.nodectime')
				list(map= topo.map, nodectime=node.ctime, ph=topo.map.ph)
			})
	cluphy.map					<- do.call('rbind',lapply(tmp, function(x) x$map))
	cluphy.map.nodectime		<- do.call('rbind',lapply(tmp, function(x) x$nodectime))
	cluphy.subtrees				<- lapply(tmp, function(x) x$ph)
	names(cluphy.subtrees)		<- file.info[,cluster]
	#	determine subtree root times and set the root time of the combined phylogeny to the mean
	cluphy.subtrees.root.ctime	<- sapply(seq_along(cluphy.subtrees), function(i)
			{
				tmp		<- treeannotator.tiplabel2df(cluphy.subtrees[[i]], beastlabel.idx.clu=beastlabel.idx.clu, beastlabel.idx.hivn=beastlabel.idx.hivn, beastlabel.idx.hivd=beastlabel.idx.hivd, beastlabel.idx.hivs=beastlabel.idx.hivs, beastlabel.idx.samplecode=beastlabel.idx.samplecode, beastlabel.idx.rate=beastlabel.idx.rate)
				tmp[1,TipT]-node.depth.edgelength(cluphy.subtrees[[i]])[1]-cluphy.subtrees[[i]]$root.edge
			})
	cluphy.root.ctime			<- mean(cluphy.subtrees.root.ctime)
	for(i in seq_along(cluphy.subtrees))
		cluphy.subtrees[[i]]$root.edge	<- cluphy.subtrees[[i]]$root.edge+cluphy.subtrees.root.ctime[i]-cluphy.root.ctime
	brl.error	<- max( sapply(seq_along(cluphy.subtrees), function(i)
					{
						tmp		<- treeannotator.tiplabel2df(cluphy.subtrees[[i]], beastlabel.idx.clu=beastlabel.idx.clu, beastlabel.idx.hivn=beastlabel.idx.hivn, beastlabel.idx.hivd=beastlabel.idx.hivd, beastlabel.idx.hivs=beastlabel.idx.hivs, beastlabel.idx.samplecode=beastlabel.idx.samplecode, beastlabel.idx.rate=beastlabel.idx.rate)
						tmp		<- node.depth.edgelength(cluphy.subtrees[[i]])[seq_len(Ntip(cluphy.subtrees[[i]]))]+cluphy.subtrees[[i]]$root.edge+cluphy.root.ctime-tmp[, TipT]
						max(abs(tmp))
					}) )
	if(brl.error>2*EPS)	warning(paste('found brl error, max error=',brl.error))
	#	combine subtree phylogenies
	cluphy						<- clu.polyphyletic.clusters(cluphy.subtrees=cluphy.subtrees)$cluphy
	cluphy.info					<- treeannotator.tiplabel2df(cluphy, beastlabel.idx.clu=beastlabel.idx.clu, beastlabel.idx.hivn=beastlabel.idx.hivn, beastlabel.idx.hivd=beastlabel.idx.hivd, beastlabel.idx.hivs=beastlabel.idx.hivs, beastlabel.idx.samplecode=beastlabel.idx.samplecode, beastlabel.idx.rate=beastlabel.idx.rate)
	cluphy$tip.label			<- cluphy.info[, FASTASampleCode]
	cluphy.tip.ctime			<- cluphy.info[, TipT]
	if(cluphy.root.ctime!=cluphy.info[1,TipT]-node.depth.edgelength(cluphy)[1]-cluphy$root.edge)	stop('unexpected root time of the combined phylogeny')
	#	need to map nodes from subtrees to nodes in combined phylogeny
	cluphy.info					<- merge(cluphy.info, cluphy.info[,	list(mrca= clu.mrca(cluphy, FASTASampleCode)$mrca),	by='cluster'], by='cluster')
	cluphy.vertexmap			<- cluphy.info[,	{
				tmp		<- cluphy.subtrees[[as.character(cluster)]]
				list( 	ph.vertex=c(mrca[1], Descendants(cluphy, mrca[1], type="all")), clu.vertex=c(Ntip(tmp)+1, Descendants(tmp, Ntip(tmp)+1, type="all")) )
			},by='cluster']
	#	cluphy.vertexmap does not contain time of root to mrca (node 0) and we don t need it.
	#	cluphy.vertexmap contains tip indices and we don t need those.
	setnames(cluphy.vertexmap,'clu.vertex','node')
	cluphy.map.nodectime		<- merge(cluphy.map.nodectime, cluphy.vertexmap, by=c('cluster','node'))
	set(cluphy.map.nodectime, NULL, 'node', cluphy.map.nodectime[,ph.vertex])
	cluphy.map.nodectime[, ph.vertex:=NULL]
	#	set node.label to MAP prob
	tmp											<- merge( unique(subset(cluphy.info, select=c(cluster, mrca))), unique(subset(cluphy.map, select=c(cluster, dens))), by='cluster' )
	cluphy$node.label							<- rep('',Nnode(cluphy))
	cluphy$node.label[tmp[,mrca-Ntip(cluphy)]]	<- as.character(tmp[,dens])
	#set(cluphy.info, NULL, 'Patient', cluphy.info[, as.character(Patient)])
	list(cluphy=cluphy, cluphy.subtrees=cluphy.subtrees, cluphy.info=cluphy.info, cluphy.map.nodectime=cluphy.map.nodectime)
}
######################################################################################
#' No description yet
#' @author Oliver Ratmann
beast2out.plot.cluster.trees<- function(df.all, df.immu, df.viro, df.treatment, ph, ph.root.ctime, ph.tip.ctime, ph.prob=NA, df.node.ctime=NULL, df.rates=NULL, df.tips=NULL, end.ctime=2013.3,  cex.nodelabel=0.5,  cex.tiplabel=0.5,  file=NULL,  pdf.width=7, pdf.height=20, pdf.xlim=NULL, label.select=c("cluster","PatientA","Trm","isAcute","RegionHospital"))
{
	#df.all, df.immu, df.viro, df.treatment,
	#df.tips	<- df.tpairs.plot; ph<- cluphy; ph.root.ctime=cluphy.root.ctime; ph.tip.ctime=cluphy.tip.ctime; df.node.ctime=cluphy.map.nodectime;  cex.nodelabel=0.5;  cex.tiplabel=0.5;  file=NULL;  pdf.width=7; pdf.height=120; pdf.xlim=pdf.xlim
	require(RColorBrewer)
	if(class(file)=="character")
		pdf(file, width=pdf.width, height=pdf.height)
	par(mar=rep(0,4))
	youngest.tip.ctime	<- max(ph.tip.ctime)
	cols				<- brewer.pal(12,"Paired")
	cols[1]				<- cols[8]
	#	get tip labels
	ph.tiplabel			<- clu.get.tiplabels(ph, 	df.all, col.notmsm="#4EB3D3", col.Early="#EF9708", col.highVL="#FEE391", col.AfterTreat="#D4B9DA", col.green="#D9F0A3", col.latePres="#FA9FB5", select=label.select )
	if(is.null(pdf.xlim))
	{
		tmp				<- max( apply(ph.tiplabel$text, 2, function(x)  sum(xinch(strwidth(x, units="inches", cex=cex.tiplabel)))  ) )
		pdf.xlim		<- c(-3, ceiling(end.ctime)-ph.root.ctime+max(2.5,tmp))
	}
	pdf.ylim			<- c(1,Ntip(ph)) + c(-1,1)
	plot(ph, x.lim=pdf.xlim, y.lim= pdf.ylim, show.tip.label=0, edge.color = 1, tip.color = 0)
	# add calendar timeline
	treeannotator.plot.ctimeline(ph, youngest.tip.ctime, end.ctime, add.yinch= 0.5)
	# add NegT and AnyPos_T1
	ph.seronodeheight	<- treeannotator.sero.getnodeheight.range(ph, df.all, youngest.tip.ctime)
	treeannotator.plot.seronodeheightrange(ph, ph.seronodeheight, add.yinch= -0.03, width.yinch= 0.03, width.yinch.past.AnyPos_T1= 0, col=cols[2])
	# add lRNA timeline
	ph.viro.timeline	<- treeannotator.get.viro.timeline(ph, df.all, df.viro, youngest.tip.ctime, df.treatment=df.treatment)
	treeannotator.plot.viro.timeline(ph, ph.viro.timeline, viro.min= log10(300), width.yinch= 0.15, add.yinch= 0.005, col.bg= cols[c(5,10,12)], col.legend= cols[6], cex.txt= 0.2, lines.lwd=0.1)
	# add CD4 timeline
	ph.immu.timeline	<- treeannotator.get.immu.timeline(ph, df.all, df.immu, youngest.tip.ctime, end.ctime=2013.3)
	treeannotator.plot.immu.timeline(ph, ph.immu.timeline, immu.min= 150, immu.max= 800, width.yinch= 0.15, add.yinch= -0.005, col.bg= cols[3], col.legend= cols[4], cex.txt= 0.2, lines.lwd=0.1)
	# add BEAST posterior density of nodes where available
	if(!is.null(df.node.ctime))
		treeannotator.plot.node.ctime(copy(df.node.ctime), ph.root.ctime, width.yinch=0.1, add.yinch=0.005, col.bg=cols[1] )
	# re-plot phylogeny
	if(!is.null(ph$node.label))
		ph$node.label	<- as.numeric(sapply( strsplit( ph$node.label, '_' ), function(x)	x[1] ))
	edge.width			<- treeannotator.get.edgewidth(ph, df.rates, scale.edgewidth= 8)
	phy.plotupon(ph, show.tip.label=0, show.node.label=ifelse(is.null(ph$node.label),0,1), cex=cex.nodelabel, edge.width=edge.width[,width])
	# add root edge
	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
	lines(c(-ph$root.edge,0),rep(lastPP$yy[Ntip(ph)+1],2), lwd=edge.width[1,width])
	# plot tick marks
	if(!is.null(df.tips))
		treeannotator.plot.tipmarks(ph, df.tips, add.xinch=0, add.yinch=0)
	# add rate labels
	if(!is.null(df.rates))
		treeannotator.plot.rates(ph, edge.width, add.xinch=-0.1, cex.rate=0.3)
	# add tip labels
	tmp					<- rep( max(node.depth.edgelength(ph)) - (youngest.tip.ctime-ceiling(end.ctime)), Ntip(ph))
	clu.plot.tiplabels(seq_len(Ntip(ph)), ph.tiplabel$text, ph.tiplabel$col, xx=tmp, adj = c(-0.05, 0.5), cex=cex.tiplabel, add.xinch= 0.03, add.yinch= 0.02)
	# add legend
	legend("topright", fill= cols[c(1,2,3,5,10,12)], legend=c("SA-BEAST2 TMRCA pdf", "interval [last HIV-, diagnosis]", "CD4 timeline", "VL timeline", "VL timeline under treatment", "VL timeline under treatment"), bty='n', border=NA, cex=cex.tiplabel)
	if(!is.na(ph.prob))
		legend("bottomright", legend=paste('prob=',ph.prob),bty='n', border=NA, cex=cex.tiplabel*2)
	if(class(file)=="character")
		dev.off()
}
######################################################################################
#' No description yet
#' @author Oliver Ratmann
beast2out.pool.cluster.trees<- function(files)
{
	tmp				<- regmatches( files, regexpr('_pool_[0-9]+',files))
	run				<- as.numeric( regmatches(tmp, regexpr('[0-9]+',tmp))	)
	mph.clu.pool	<- NULL
	for(i in seq_along(run))
	{
		options(show.error.messages = FALSE)
		readAttempt		<- try(suppressWarnings(load(files[i])))
		if(!inherits(readAttempt, "try-error"))	cat(paste("\nresumed file",files[i]))
		if(inherits(readAttempt, "try-error"))	stop(paste("\ncould not read file",files[i]))
		names(mph.clu)	<- paste('RUN_',run[i],'_',names(mph.clu),sep='')
		mph.clu.pool	<- c(mph.clu.pool, mph.clu)
	}
	mph.clu.pool
}
######################################################################################
#' No description yet
#' @author Oliver Ratmann
beast2out.extract.cluster.trees<- function(mph, mph.info)
{
	mph.clusters		<- unique(mph.info[,cluster])
	mph.by.clu			<- lapply(mph.clusters, function(clu)
			{
				cat(paste('\nprocess cluster',clu))
				mph.clu			<- lapply( seq_along(mph), function(mph.i)
						{
							clu.seq			<- subset(mph.info, cluster==clu)[,BEASTlabel]
							clu.mrca		<- clu.mrca(mph[[mph.i]], clu.seq)
							extract.clade( mph[[mph.i]], clu.mrca$mrca, root.edge=clu.mrca$mrca.height, interactive=FALSE )
						})
				names(mph.clu)	<- names(mph)
				mph.clu
			})
	names(mph.by.clu)	<- mph.clusters
	mph.by.clu
}
######################################################################################
#' No description yet
#' @author Oliver Ratmann
beast2out.tip.date.check<- function(ph, fun, ...)
{
	ph.info	<- fun(ph, ...)
	ph.info[, tip:= match(ph.info[,BEASTlabel], ph$tip.label)]
	setkey(ph.info, tip)
	tmp		<- node.depth.edgelength(ph)[seq_len(Ntip(ph))]
	max( abs( ph.info[1,TipT]-tmp[1]+tmp  -  ph.info[,TipT] ) )
}
######################################################################################
#' Extract a list of phylogenies from a BEAST2 posterior file
#' @param file name of the BEAST2 posterior filename, usually ends with '.trees'
#' @param opt.rescale.edge.length Rescaling factor of the phylogenies,
#'   a value of 1.0 denotes that the original edge lengths are maintained
#' @param opt.burnin how many phylogenies to discard, an value of will keep all trees
#' @return a list of phylogenies of type 'phylo'
#' @examples
#'   trees_file <- "vignettes/example.trees"
#'   testit::assert(file.exists(trees_file))
#'   posterior <- beast2out.read.trees(trees_file)
#'   testit::assert(length(posterior) == 10)
#'   testit::assert(class(posterior[[1]]) == "phylo")
#' @export
#' @author Oliver Ratmann
beast2out.read.trees<- function(
  file,
  opt.rescale.edge.length= 1.,
  opt.burnin=0
)
{
	tmp			<- readLines(file, n=2e3, warn = FALSE)
	tmp			<- which( grepl('#NEXUS', tmp) )
	if(length(tmp)>1)
	{
		cat(paste('\nFound #NEXUS headers, n=',length(tmp),'.\nDiscard all lines before last entry on line', tail(tmp,1)))
		cmd		<- paste('sed -i".bak" 1,',tail(tmp,1)-1,'d ', file, sep='')
		system(cmd)
		cmd		<- paste('sed -i".bak2" 1s/\\;// ', file, sep='')
		system(cmd)
		cmd		<- list.files(paste(rev(rev(strsplit(file, '/')[[1]])[-1]),collapse='/'), pattern='*bak*', full.names=TRUE)
		cat(paste('\nrm files\n', paste(cmd, collapse='\n')))
		file.remove(cmd)
	}
	mph			<- ape::read.nexus(file)
	#	remove burn in
	tmp			<- regexpr('[0-9]+',names(mph))
	if(any(tmp<0))	stop('unexpected nexus file without STATE iteration numbers')
	mph.it		<- as.numeric( regmatches( names(mph), tmp) )
	mph			<- lapply( which( mph.it>opt.burnin), function(j)	mph[[j]]	)
	mph.it		<- mph.it[ mph.it > opt.burnin ]
	names(mph)	<- paste('STATE_',mph.it,sep='')
	#	rescale edge lengths
	if(opt.rescale.edge.length!=1.)
		for(j in seq_along(mph))
			mph[[j]]$edge.length	<- mph[[j]]$edge.length * opt.rescale.edge.length
	mph
}
######################################################################################
#' No description yet
#' @author Oliver Ratmann
beast2out.get.log.evolrate<- function(files, opt.burnin=opt.burnin)
{
	if(!opt.burnin)	warning('burn in equals zero')
	#	collect posterior TreeHeight samples
	df.log		<- lapply(seq_along(files), function(i)
			{
				file		<- files[i]
				cat(paste("\nReading file ",file))
				df.log		<- read.delim2(file, header = TRUE, sep = "\t", quote="\"", dec=".", fill = TRUE, comment.char="#")
				df.log		<- as.data.table(df.log)
				is.BEAST2	<- ifelse('Sample'%in%colnames(df.log), 1, 0)
				cat(paste("\ndetected is.BEAST2=",is.BEAST2))
				if(is.BEAST2)
				{
					if(!'Sample'%in%colnames(df.log)) stop('expect column Sample in log file')
					if(!'ucldMean'%in%colnames(df.log)) stop('expect column ucldMean in log file')
					df.log	<- subset(df.log, Sample>opt.burnin, select=c(Sample, ucldMean))
				}
				else
				{
					if(!'state'%in%colnames(df.log)) stop('expect column state in log file')
					if(!'ucld.mean'%in%colnames(df.log)) stop('expect column ucld.mean in log file')
					df.log	<- subset(df.log, state>opt.burnin, select=c(state, ucld.mean))
					setnames(df.log, 'ucld.mean','ucldMean')
				}
				df.log[, file.i:= i]
				df.log
			})
	df.log		<- do.call('rbind', df.log)
	#	compute overall mean and file specific mean
	ans			<- df.log[, list(ucldMean.mean.i=mean(ucldMean)), by='file.i']
	ans[,ucldMean.mean:= mean(ucldMean.mean.i)]
	ans
}
######################################################################################
#' No description yet
#' @export
#' @author Oliver Ratmann
beast2out.read.nexus.and.stats<- function(file, tree.id=NA, method.node.stat='any.node')
{
	stopifnot(method.node.stat%in%c('any.node','inner.node'))

	X				<- scan(file = file, what = "", sep = "\n", quiet = TRUE)
	#	read TRANSLATE chunk
	X.endblock		<- grep("END;|ENDBLOCK;|End;", X, ignore.case = TRUE)
	X.semico 		<- grep(";", X)
	X.i1 			<- grep("BEGIN TREES;|Begin trees;", X, ignore.case = TRUE)
	X.i2 			<- grep("TRANSLATE|Translate", X, ignore.case = TRUE)
	tmp 			<- X.semico[X.semico > X.i2][1]
	tmp 			<- X[(X.i2 + 1):tmp]
	tmp				<- gsub('[,;]$','',gsub('^\\s+','',tmp))
	tmp				<- tmp[nzchar(tmp)]
	tmp				<- strsplit(tmp, ' ')
	df.translate	<- data.table(NEXUS_ID= sapply(tmp, '[[', 1), NEXUS_LABEL=sapply(tmp, '[[', 2) )
	set(df.translate, NULL, 'NEXUS_LABEL', df.translate[, gsub("\'","",NEXUS_LABEL)])
	cat(paste('\nFound taxa, n=', nrow(df.translate)))

	if(!is.na(tree.id))
	{
		#	read one newick tree with id 'tree.id'
		bstr		<- X[grep(paste(tree.id,"[[:space:]]+",sep=''), X)]
		node.stat	<- beast2out.read.nodestats(bstr)
		cat(paste('\nFound node statistics, n=', nrow(node.stat)))
		set(node.stat, NULL, 'tree.id', tree.id[i] )
		btree		<- beast2out.read.nodeidtree(bstr, method.node.stat=method.node.stat)
		#
		# link node.stats with tree nodes (tip + inner node)
		# NODE_ID is index of node in 'btree' phylo object
		#
		tmp			<- strsplit( btree$tip.label, 'NODE_PARSE_ID' )
		df.link		<- data.table(NODE_ID=seq_along(btree$tip.label), NEXUS_ID=sapply(tmp,'[[',1), NODE_PARSE_ID=sapply(tmp,'[[',2))
		df.link		<- merge(df.link, df.translate, by='NEXUS_ID')
		cat(paste('\nFound tree tips with taxon name, n=', nrow(df.link)))
		tmp			<- strsplit( btree$node.label, 'NODE_PARSE_ID' )
		tmp			<- data.table(NODE_ID=Ntip(btree)+seq_along(btree$node.label), NODE_PARSE_ID=sapply(tmp,'[[',2), NEXUS_LABEL=NA_character_)
		df.link		<- rbind(subset(df.link,select=c(NODE_ID, NODE_PARSE_ID, NEXUS_LABEL)), tmp)
		set(df.link,NULL,'NODE_PARSE_ID',df.link[, as.integer(NODE_PARSE_ID)])
		set(df.link,NULL,'NODE_ID',df.link[, as.integer(NODE_ID)])
		set(df.link,NULL,'TREE_ID',tree.id)
		node.stat	<- merge( node.stat, subset(df.link, select=c(NODE_PARSE_ID, NODE_ID, TREE_ID)), by='NODE_PARSE_ID' )
		set(node.stat,NULL,'NODE_PARSE_ID',NULL)
		cat(paste('\nLinked node statistics to tree nodes, n=', nrow(node.stat)))
		#
		# set tip.labels and rm node.labels
		#
		setkey(df.link, NODE_ID)
		btree$tip.label		<- df.link[seq_len(Ntip(btree)),][,NEXUS_LABEL]
		btree$node.label	<- NULL
	}
	if(is.na(tree.id))
	{
		#	read all newick trees in nexus file
		tmp			<- regexpr('^tree\\s\\S+',X)
		tree.id		<- sapply( regmatches(X,tmp), function(x) substr(x, 5, nchar(x)))
		tree.id		<- gsub('\\s','',tree.id)
		cat(paste('\nFound tree id=', paste(tree.id, collapse=' ')))
		X			<- X[ which(tmp>0) ]
		cat(paste('\nFound trees, n=',length(tree.id)))
		node.stat	<- lapply(seq_along(tree.id), function(i)
				{

					bstr	<- X[grep(paste(tree.id[i],"[[:space:]]+",sep=''), X)]
					cat(paste('\nGet node statistics for tree id=',tree.id[i]))
					tmp		<- beast2out.read.nodestats(bstr)
					set(tmp, NULL, 'TREE_ID', tree.id[i] )
					tmp
				})
		if(length(node.stat)>1)
			node.stat	<- do.call('rbind',node.stat)
		if(length(node.stat)==1)
			node.stat	<- node.stat[[1]]
		node.stat[, NODE_ID:=NA_integer_]
		setkey(node.stat, TREE_ID, NODE_PARSE_ID)
		btree		<- vector('list',length(tree.id))
		for(i in seq_along(tree.id))
		{
			bstr		<- X[grep(paste(tree.id[i],"[[:space:]]+",sep=''), X)]
			cat(paste('\nRead tree for tree id=',tree.id[i]))
			btree.i		<- beast2out.read.nodeidtree(bstr, method.node.stat=method.node.stat)
			#
			# link node.stats with tree nodes (tip + inner node)
			# NODE_ID is index of node in 'btree.i' phylo object
			#
			tmp			<- strsplit( btree.i$tip.label, 'NODE_PARSE_ID' )
			df.link		<- data.table(NODE_ID=seq_along(btree.i$tip.label), NEXUS_ID=sapply(tmp,'[[',1), NODE_PARSE_ID=sapply(tmp,'[[',2))
			df.link		<- merge(df.link, df.translate, by='NEXUS_ID')
			cat(paste('\nFound tree tips with taxon name, n=', nrow(df.link)))
			tmp			<- strsplit( btree.i$node.label, 'NODE_PARSE_ID' )
			tmp			<- data.table(NODE_ID=Ntip(btree.i)+seq_along(btree.i$node.label), NODE_PARSE_ID=sapply(tmp,'[[',2), NEXUS_LABEL=NA_character_)
			df.link		<- rbind(subset(df.link,select=c(NODE_ID, NODE_PARSE_ID, NEXUS_LABEL)), tmp)
			set(df.link,NULL,'NODE_PARSE_ID',df.link[, as.integer(NODE_PARSE_ID)])
			set(df.link,NULL,'NODE_ID',df.link[, as.integer(NODE_ID)])
			for(j in seq_len(nrow(df.link)))
				set(node.stat, node.stat[, which(TREE_ID==tree.id[i] & NODE_PARSE_ID==df.link[j,NODE_PARSE_ID])], 'NODE_ID', df.link[j,NODE_ID])
			tmp			<- node.stat[, length(which(!is.na(NODE_ID)))]
			cat(paste('\nTotal linked node statistics to tree nodes, n=', tmp  ))
			#
			# set tip.labels and rm node.labels
			#
			setkey(df.link, NODE_ID)
			btree.i$tip.label	<- df.link[seq_len(Ntip(btree.i)),][,NEXUS_LABEL]
			btree.i$node.label	<- NULL
			btree[[i]]			<- btree.i
		}
		if(length(btree)>=2)
		{
			names(btree)	<- tree.id
			class(btree)	<- "multiPhylo"
		}
		if(length(btree)<2)
			btree	<- btree[[1]]
		tmp				<- node.stat[, length(which(is.na(NODE_ID)))]
		cat(paste('\nTotal unlinked node statistics [should be zero], n=', tmp  ))
		set(node.stat,NULL,'NODE_PARSE_ID',NULL)
	}
	list(tree=btree, node.stat=node.stat)
}
######################################################################################
#' Private funtion to read tree such that tip.label and nodel.label include
#' information on the index on when a tip/node occurs in bstr
#' @author Oliver Ratmann
beast2out.read.nodeidtree <- function(bstr, method.node.stat='any.node')
{
	# strip all meta variables and ; at end
	bstr		<- gsub("\\[[^]]*\\]", "", bstr)
	bstr		<- gsub(';','',bstr)
	# for each node, add a dummy node label NODE_PARSE_IDxx
	dummy.tree	<- unlist(strsplit(bstr, ":"))
	if(method.node.stat=='inner.node')
	{
		#	interior branch length: 	previous index ends in ). so tmp is the index of the dummy.tree chunks that gives the start of a branch length of an inner node
		tmp			<- which( c(FALSE, grepl(')$',dummy.tree)[-length(dummy.tree)]) )
		#	prepend NODE_PARSE_IDxx before the branch length of an inner node
		tmp			<- tmp-1
	}
	if(method.node.stat=='any.node')
		tmp			<- seq_along(dummy.tree)
	dummy.tree	<- sapply(seq_along(dummy.tree), function(i)
			{
				z<- which(i==tmp)
				ifelse(length(z),	paste(dummy.tree[i],'NODE_PARSE_ID',z,sep=''),	dummy.tree[i] )
			})
	dummy.tree	<- paste(dummy.tree, collapse=':',sep='')
	dummy.tree	<- regmatches(dummy.tree, regexpr('\\(.*',dummy.tree))
	dummy.tree	<- paste(dummy.tree, ';', sep='')
	ph<-  tryCatch(
			{
				read.tree(text=dummy.tree)
			}, error=function(e)
			{
				cat(paste('\nerror in read.tree\n',e$message,'\ntry seq.read.newick'))
				return( seq.read.newick(text=dummy.tree) )
			}, warning=function(e)
			{
				cat(paste('\nwarning in read.tree\n',e$message,'\ntry seq.read.newick'))
				return( seq.read.newick(text=dummy.tree) )
			})
	ph
}
######################################################################################
#' No description yet
#' @author Oliver Ratmann
beast2out.read.nodestats <- function(bstr)
{
	#	remove anything before first '('
	bstr	<- regmatches(bstr, regexpr('\\(.*',bstr))
	# 	store meta info for inner nodes that is given in [], and not in :[] which is meta info for edges
	tmp		<- unlist(regmatches(bstr,gregexpr('[^:]\\[[^]]+',bstr)))
	tmp		<- sapply( tmp, function(x) substr(x, 4, nchar(x)) )
	#	for each inner node, extract stats
	tmp		<- strsplit(tmp, ',')
	tmp		<- lapply(seq_along(tmp), function(i)
			{
				z<- strsplit(tmp[[i]],'=')
				data.table(NODE_PARSE_ID=i, STAT=sapply(z,'[',1), VALUE=sapply(z,'[',2))
			})
	node.stat	<- do.call('rbind', tmp)
	tmp			<- node.stat[, unique(STAT)]
	cat(paste('\nFound node statistics=',paste(tmp,collapse=' ')))
	tmp			<- node.stat[, list(has.all.stats= !length(setdiff(tmp, STAT))  ) , by='NODE_PARSE_ID']
	tmp			<- subset(tmp, !has.all.stats)[, NODE_PARSE_ID]
	cat(paste('\nSome statistics missing for nodes=',paste(tmp,collapse=' ')))
	node.stat
}
