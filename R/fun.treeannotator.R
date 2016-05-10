######################################################################################
#' No description yet
#' @author Oliver Ratmann
treeannotator.tiplabel2df<- function(x, beastlabel.idx.clu=1, beastlabel.idx.hivn=2, beastlabel.idx.hivd=3, beastlabel.idx.hivs=4, beastlabel.idx.samplecode=5, beastlabel.idx.rate=6)
{
	tmp		<- t( sapply(strsplit(x$tip.label,'_'), function(z)	z[c(beastlabel.idx.clu,beastlabel.idx.hivn, beastlabel.idx.hivd, beastlabel.idx.hivs, beastlabel.idx.samplecode, beastlabel.idx.rate)] ) )
	data.table(cluster=as.numeric(tmp[,1]), NegT= suppressWarnings(as.numeric(tmp[,2])), AnyPos_T1= as.numeric(tmp[,3]), TipT= as.numeric(tmp[,4]), FASTASampleCode=tmp[,5], rate=as.numeric(tmp[,6]), BEASTlabel=x$tip.label )
}
######################################################################################
#' No description yet
#' @author Oliver Ratmann
treeannotator.get.clusterprob<- function(ph.beast, beastlabel.idx.clu=1, beastlabel.idx.samplecode=5, verbose=1)
{
	#	for each of the clusters, compute the posterior probability of a common MRCA
	clu.df	<- lapply(ph.beast, function(x)
			{
				tmp		<- data.table(	tip=seq_len(Ntip(x)),
						cluster=as.numeric( sapply( strsplit(x$tip.label,'_'),function(z)  z[beastlabel.idx.clu] ) ),
						FASTASampleCode=sapply( strsplit(x$tip.label,'_'),function(z)  z[beastlabel.idx.samplecode] )
				)
				ans		<- merge( tmp[, list(node=clu.mrca(x, x.tip=tip)$mrca, FASTASampleCode=FASTASampleCode), by=cluster], subset(x$node.label, select=c(node, posterior)), by="node" )
				subset(ans, select=c(cluster, FASTASampleCode, posterior))
			})
	clu.df	<- rbindlist(clu.df)
	if(verbose)	cat(paste("\nRange of posterior probabilities that each of the putative clusters each has a common MRCA, min=",min(clu.df[,posterior])," max=",max(clu.df[,posterior]) ))
	clu.df
}
######################################################################################
#' No description yet
#' @author Oliver Ratmann
treeannotator.get.edgewidth<- function(ph, rates.df, scale.edgewidth= 12)
{
	if(is.null(rates.df))
		rates.df<- data.table(node=Nnode(ph,internal.only=0), rate=NA)

	edge.width	<- merge(rates.df, data.table(node=ph$edge[,2], edge=seq_len(nrow(ph$edge))), all.y=1, by="node")
	edge.width[, width:=edge.width[,rate]/mean(edge.width[,rate], na.rm=1)]
	set(edge.width, NULL, "width", (edge.width[,width]-1)*scale.edgewidth + 1)
	set(edge.width, which(edge.width[,is.na(width)]), "width", 1)
	set(edge.width, which(edge.width[,width<0.1]), "width", 0.1)
	setkey(edge.width, edge)
	edge.width
}
######################################################################################
#' No description yet
#' @author Oliver Ratmann
#	ph<- cluphy; end.ctime=2013.3; cex.nodelabel=0.5; cex.tiplabel=0.5; file=NULL; pdf.width=7; pdf.height=20
treeannotator.plot<- function(ph, ph.root.ctime, youngest.tip.ctime, df.all, df.viro, df.immu, df.treatment=NULL, df.tstem=NULL, df.rates=NULL, end.ctime=2013.3, cex.nodelabel=0.5, cex.tiplabel=0.5, file=NULL, pdf.width=7, pdf.height=20)
{
	require(RColorBrewer)
	if(class(file)=="character")
		pdf(file, width=pdf.width, height=pdf.height)
	par(mar=c(0,0,0,0))

	cols			<- brewer.pal(12,"Paired")
	ph.xlim			<- end.ctime-ph.root.ctime+ c(-22,6)
	ph.ylim			<- c(1,Ntip(ph)) + c(-1,1)

	plot(ph, x.lim=ph.xlim, y.lim= ph.ylim, show.tip.label=0, edge.color = 0, tip.color = 0)
	# add calendar timeline
	treeannotator.plot.ctimeline(ph, youngest.tip.ctime, end.ctime, add.yinch= 0.5)
	# add BEAST TMRCA 95% credibility interval
	ph.nodeheighthpds	<- treeannotator.nodelabels.getnodeheightHPD(ph, youngest.tip.ctime)
	#do not plot the BEAST TMRCA 95% credibility interval for those nodes for which more data than the 95% interval is available
	if(!is.null(df.tstem))
		ph.nodeheighthpds	<- subset( ph.nodeheighthpds, !node%in%unique(df.tstem[,mrca]) )
	treeannotator.plot.hpdbars(ph, ph.nodeheighthpds, col=cols[1], lwd=4)
	# add NegT and AnyPos_T1
	ph.seronodeheight	<- treeannotator.sero.getnodeheight.range(ph, df.all, youngest.tip.ctime)
	treeannotator.plot.seronodeheightrange(ph, ph.seronodeheight, add.yinch= -0.03, width.yinch= 0.03, width.yinch.past.AnyPos_T1= 0, col=cols[2])
	# add lRNA timeline
	ph.viro.timeline	<- treeannotator.get.viro.timeline(ph, df.all, df.viro, youngest.tip.ctime, df.treatment=df.treatment)
	treeannotator.plot.viro.timeline(ph, ph.viro.timeline, viro.min= log10(300), width.yinch= 0.15, add.yinch= 0.005, col.bg= cols[c(5,10,12)], col.legend= cols[6], cex.txt= 0.2)
	# add BEAST posterior density of TMRCAs where available
	if(!is.null(df.tstem))
		treeannotator.plot.tipstem.timeline(ph, youngest.tip.ctime, df.tstem, width.yinch=0.1, add.yinch=0.005, col.bg=cols[1] )
	# add CD4 timeline
	ph.immu.timeline	<- treeannotator.get.immu.timeline(ph, df.all, df.immu, youngest.tip.ctime, end.ctime=2013.3)
	treeannotator.plot.immu.timeline(ph, ph.immu.timeline, immu.min= 150, immu.max= 800, width.yinch= 0.15, add.yinch= -0.005, col.bg= cols[3], col.legend= cols[4], cex.txt= 0.2)
	# re-plot phylogeny
	ph$node.label		<- as.numeric(sapply( strsplit( ph$node.label, '_' ), function(x)	x[1] ))
	edge.width			<- treeannotator.get.edgewidth(ph, df.rates, scale.edgewidth= 8)
	phy.plotupon(ph, show.tip.label=0, show.node.label=1, cex=cex.nodelabel, edge.width=edge.width[,width],)
	# add rate labels
	if(!is.null(df.rates))
		treeannotator.plot.rates(ph, edge.width, add.xinch=-0.1, cex.rate=0.3)
	# add tip labels
	ph.tiplabel			<- clu.get.tiplabels(ph, 	df.all, col.notmsm="#4EB3D3", col.Early="#EF9708", col.highVL="#FEE391", col.AfterTreat="#D4B9DA", col.green="#D9F0A3", col.latePres="#FA9FB5", select=c("CountryInfection","Trm","Sex","isAcute","lRNA.early","Patient","RegionHospital") )
	tmp					<- rep( max(node.depth.edgelength(ph)) - (youngest.tip.ctime-ceiling(end.ctime)), Ntip(ph))
	clu.plot.tiplabels(seq_len(Ntip(ph)), ph.tiplabel$text, ph.tiplabel$col, xx=tmp, adj = c(-0.05, 0.5), cex=cex.tiplabel, add.xinch= 0.03, add.yinch= 0.02)
	# add legend
	legend("topright", fill= cols[c(1,2,3,5,10,12)], legend=c("BEAST 95% TMRCA", "interval [last HIV-, diagnosis]", "CD4 timeline", "VL timeline", "VL timeline under treatment", "VL timeline under treatment"), bty='n', border=NA, cex=cex.tiplabel)

	if(class(file)=="character")
		dev.off()
}
######################################################################################
#' No description yet
#' @author Oliver Ratmann
treeannotator.plot.tipmarks<- function(ph, df.tips, add.xinch=0, add.yinch=0)
{
	tmp	<- match(df.tips[, FASTASampleCode], ph$tip.label)
	if(any(is.na(tmp)))		stop('unexpected missing tip names in ph for df.tips')
	lastPP 	<- get("last_plot.phylo", envir = .PlotPhyloEnv)
	tmp		<- data.table(xx= lastPP$xx[ tmp ]+xinch(add.xinch), yy= lastPP$yy[ tmp ], tips=tmp)
	tmp		<- cbind(df.tips, tmp)
	points(tmp[,xx], tmp[,yy], col=tmp[,col], pch=tmp[,pch], cex=tmp[,cex] )
}
######################################################################################
#' No description yet
#' @author Oliver Ratmann
treeannotator.plot.rates<- function(ph, edge.width, add.xinch=-0.1, cex.rate=0.5)
{
	lastPP 	<- get("last_plot.phylo", envir = .PlotPhyloEnv)
	tmp	<- edge.width[, list(xx= lastPP$xx[ ph$edge[edge,1] ]+xinch(add.xinch), yy= mean( lastPP$yy[ ph$edge[edge,] ] ), rate=rate), by="edge"]
	text(tmp[,xx], tmp[,yy], tmp[,rate], cex=cex.rate)
}
######################################################################################
#' No description yet
#' @author Oliver Ratmann
treeannotator.plot.node.ctime<- function(df.node.ctime, ph.root.ctime, width.yinch=0.1, add.yinch=0.001, col.bg="black", density=30, lwd=0.5)
{
	lastPP 	<- get("last_plot.phylo", envir = .PlotPhyloEnv)
	df.node.ctime[, yyi:=pdf]
	df.node.ctime[, xx:=q-ph.root.ctime]
	setkey(df.node.ctime, node)
	df.node.ctime	<- merge(df.node.ctime, df.node.ctime[,	list(scale=yinch(width.yinch)/max(yyi))	,by="node"], by="node")
	set(df.node.ctime, NULL, "yyi",  df.node.ctime[,yyi*scale])
	df.node.ctime[, yy:= yinch(add.yinch)+lastPP$yy[df.node.ctime[,node]]]
	df.node.ctime[, col:=col.bg]
	setkey(df.node.ctime, node)
	dummy<- sapply( unique(df.node.ctime[,node]), function(x)
			{
				z		<- df.node.ctime[J(x)]
				setkey(z, xx)
				polygon( c( z[,xx],z[nrow(z),xx] ), z[1,yy]+c( z[,yyi],0 ), border=z[1,col], col=z[1,col], density=density, lwd=lwd )
			})
}
######################################################################################
#' No description yet
#' @author Oliver Ratmann
treeannotator.plot.tipstem.timeline<- function(ph, youngest.tip.ctime, df.tstem, width.yinch=0.15, add.yinch=0, col.bg="grey75", density=30, lwd=0.5)
{
	lastPP 	<- get("last_plot.phylo", envir = .PlotPhyloEnv)
	set(df.tstem, NULL, "tstem", youngest.tip.ctime - df.tstem[,tstem])
	set(df.tstem, NULL, "tstem", max(node.depth.edgelength(ph)) - df.tstem[,tstem])
	df.tstem[, yyi:= density]
	setkey(df.tstem, tip)
	df.tstem<- merge(df.tstem, df.tstem[,	list(scale=yinch(width.yinch)/max(yyi))	,by="tip"], by="tip")
	set(df.tstem, NULL, "yyi",  df.tstem[,yyi*scale])
	df.tstem[, yy:= yinch(add.yinch)+lastPP$yy[df.tstem[,tip]]]
	df.tstem[, col:=col.bg]
	setkey(df.tstem, tip)
	dummy<- sapply( unique(df.tstem[,tip]), function(x)
			{
				z		<- df.tstem[J(x)]
				setkey(z, tstem)
				polygon( c( z[,tstem],z[nrow(z),tstem] ), z[1,yy]+c( z[,yyi],0 ), border=z[1,col], col=z[1,col], density=density, lwd=lwd )
			})
}
######################################################################################
#' No description yet
#' @author Oliver Ratmann
treeannotator.plot.ctimeline<- function(ph, youngest.tip.ctime, end.ctime, col.bg= c(my.fade.col("black",0.15),my.fade.col("black",0.05)), col.txt= c(my.fade.col("black",1),"transparent"), cex.txt= 0.5, add.yinch= 0.5)
{
	lastPP 	<- get("last_plot.phylo", envir = .PlotPhyloEnv)
	tmp 	<- seq( max(lastPP$xx) - ( youngest.tip.ctime-floor(youngest.tip.ctime) )+1, lastPP$x.lim[1], -1 )
	df.time	<- data.table(ctime= seq( floor(youngest.tip.ctime), by=-1, len= length(tmp) ), xx.u=tmp)
	#	if tips don t reach to end.ctime, add more bars
	if(1+floor(youngest.tip.ctime)<floor(end.ctime))
	{
		tmp		<- seq( 1+floor(youngest.tip.ctime),floor(end.ctime),by=1 )
		tmp		<- data.table( ctime=rev(tmp), xx.u=rev(seq(df.time[1,xx.u]+1, len=length(tmp), by=1)))
		df.time	<- rbind(tmp, df.time)
	}
	df.time	<- cbind( 	df.time[-nrow(df.time), ],
			data.table(	xx.l	= df.time[-1,xx.u],
					col.bg	= rep(col.bg,ceiling(nrow(df.time)/length(col.bg)))[seq_len(nrow(df.time)-1)],
					col.txt = rep(col.txt,ceiling(nrow(df.time)/length(col.txt)))[seq_len(nrow(df.time)-1)]		) )

	rect(df.time[,xx.l], lastPP$y.lim[1]-yinch(add.yinch), df.time[,xx.u], lastPP$y.lim[2]+yinch(add.yinch), col=df.time[,col.bg], border=NA)
	text(df.time[,xx.l+(xx.u-xx.l)/2.1], lastPP$y.lim[1]-yinch(add.yinch)/4, df.time[,ctime], cex=cex.txt, col=df.time[,col.txt], offset=0)
	text(df.time[,xx.l+(xx.u-xx.l)/2.1], lastPP$y.lim[2]+yinch(add.yinch)/4, df.time[,ctime], cex=cex.txt, col=df.time[,col.txt], offset=0)
}
######################################################################################
#' No description yet
#' @author Oliver Ratmann
treeannotator.plot.immu.timeline<- function(ph, ph.immu.timeline, immu.min= 150, immu.max= 800, immu.legend= c(200, 350, 500, immu.max), width.yinch= 0.15, add.yinch= -0.005, col.bg= cols[3], col.legend= cols[4], cex.txt= 0.2, lines.lwd=0.2)
{
	lastPP 	<- get("last_plot.phylo", envir = .PlotPhyloEnv)
	set(ph.immu.timeline, NULL, "PosCD4", max(node.depth.edgelength(ph)) - ph.immu.timeline[,PosCD4])
	ph.immu.timeline[, yyCD4:= CD4]
	set(ph.immu.timeline, which(ph.immu.timeline[,yyCD4]>immu.max), "yyCD4", immu.max)
	set(ph.immu.timeline, NULL, "yyCD4", ph.immu.timeline[,yyCD4]-immu.min)
	set(ph.immu.timeline, which(ph.immu.timeline[,yyCD4]<0), "yyCD4", 0.)
	scale	<- yinch(width.yinch) / max( ph.immu.timeline[,yyCD4])
	set(ph.immu.timeline, NULL, "yyCD4", ph.immu.timeline[,yyCD4] * scale )
	ph.immu.timeline[, yy:= yinch(add.yinch)+lastPP$yy[ph.immu.timeline[,tip]]]

	dummy<- sapply( unique(ph.immu.timeline[,tip]), function(x)
			{
				z<- ph.immu.timeline[J(x)]
				polygon( c( z[,PosCD4], z[nrow(z),PosCD4], z[1,PosCD4] ), c( z[,yy-yyCD4], z[nrow(z),yy], z[1,yy] ), border=NA, col=col.bg	)
				sapply(z[1,yy]-(immu.legend-immu.min)*scale,function(i)		lines(z[c(1,nrow(z)),PosCD4], rep(i,2), col=col.legend, lty=3, lwd=lines.lwd)		)
				text(rep(z[nrow(z),PosCD4],3),z[1,yy]-(immu.legend-immu.min)*scale,immu.legend,cex=cex.txt, col=col.legend)
			})
}
######################################################################################
#' No description yet
#' @author Oliver Ratmann
treeannotator.get.immu.timeline<- function(ph, df, df.immu, youngest.tip.ctime, end.ctime=2013.3)
{
	setkey(df, FASTASampleCode)
	tmp				<- df[J(ph$tip.label)][, list(tip=seq_along(ph$tip.label), FASTASampleCode=FASTASampleCode, Patient=Patient, DateDied=DateDied)]
	ans				<- merge(subset(df.immu, select=c(Patient, PosCD4, CD4)), tmp, by="Patient")
	set(ans,NULL,"PosCD4",	youngest.tip.ctime - db.Date2numeric(ans[,PosCD4]))
	set(ans,NULL,"DateDied",	youngest.tip.ctime - db.Date2numeric(ans[,DateDied]))
	set(ans,which(is.na(ans[,DateDied])), "DateDied", youngest.tip.ctime - end.ctime)
	ans				<- ans[,list(PosCD4=c(PosCD4,DateDied[1]), CD4=c(CD4,tail(CD4,1))),by="tip"]
	setkey(ans,tip)
	ans
}
######################################################################################
#' No description yet
#' @author Oliver Ratmann
#	viro.min= log10(300); width.yinch= 0.15; add.yinch= 0.005; col.bg= cols[c(5,9,10)]; col.legend= cols[6]; cex.txt= 0.2
treeannotator.plot.viro.timeline<- function(ph, ph.viro.timeline, viro.min= log10(300), viro.legend=c(3,4,5,6), width.yinch= 0.2, add.yinch= 0.005, col.bg= "red", col.legend="red", cex.txt= 0.2, lines.lwd=0.2)
{
	lastPP 	<- get("last_plot.phylo", envir = .PlotPhyloEnv)
	set(ph.viro.timeline, NULL, "PosRNA", max(node.depth.edgelength(ph)) - ph.viro.timeline[,PosRNA])
	ph.viro.timeline[, yylRNA:= lRNA]
	set(ph.viro.timeline, NULL, "yylRNA", ph.viro.timeline[,yylRNA]-viro.min)
	set(ph.viro.timeline, ph.viro.timeline[,which(yylRNA<0)], "yylRNA", 0.)
	scale	<- yinch(width.yinch) / max( ph.viro.timeline[,yylRNA])
	set(ph.viro.timeline, NULL, "yylRNA", ph.viro.timeline[,yylRNA] * scale )
	ph.viro.timeline[, yy:= yinch(add.yinch)+lastPP$yy[ph.viro.timeline[,tip]]]
	#define treatment periods if there
	if(!any(ph.viro.timeline[, NoDrug!=0]))
		ph.viro.timeline[, col:=col.bg[1]]
	else
	{
		if(length(col.bg)==1)	stop("for TPeriod != NA, expect more than one 'col.bg'")
		#print(max(ph.viro.timeline[, TPeriod], na.rm=TRUE))
		#print(ph.viro.timeline)
		#print(subset(ph.viro.timeline, is.na(TPeriod)))
		col.bg.nNA	<- rep(col.bg[-1], ceiling(max(ph.viro.timeline[, TPeriod],na.rm=TRUE)/(length(col.bg)-1)))
		ph.viro.timeline[, col:=""]
		set(ph.viro.timeline, which(ph.viro.timeline[,NoDrug==0]), "col", col.bg[1])
		set(ph.viro.timeline, which(ph.viro.timeline[,NoDrug!=0]), "col", col.bg.nNA[ subset(ph.viro.timeline, NoDrug!=0)[,TPeriod] ])
	}
	setkey(ph.viro.timeline, tip)
	dummy<- sapply( unique(ph.viro.timeline[,tip]), function(x)
			{
				z		<- ph.viro.timeline[J(x)]
				#reset TPeriod because there can be multiple off treatment periods (ie TPeriod==0) and we cannot lump them together in the next line
				#NOT NEEDED ANY LONGER	set(z, NULL, "TPeriod",cumsum(c(0,as.numeric(abs(diff(z[,TPeriod]))>0))))
				dummy	<- z[,	{
							polygon( c( PosRNA, PosRNA[length(PosRNA)], PosRNA[1] ), c( yylRNA+yy, yy[length(yy)], yy[1] ), border=NA, col=col[1]	)
						}, by="TPeriod"]
				sapply(z[1,yy]+(viro.legend-viro.min)*scale,function(i)		lines(z[c(1,nrow(z)),PosRNA], rep(i,2), col=col.legend, lty=3, lwd=lines.lwd)		)
				text(rep(z[nrow(z),PosRNA],3),z[1,yy]+(viro.legend-viro.min)*scale,paste("1e",viro.legend,sep=''),cex=cex.txt, col=col.legend)
				#stop()
			})
}
######################################################################################
#' return lRNA measurements for treatment periods (going to be different
#' colors) for tips for which at least one RNA value is available.
#' PosRNA is translated relative to youngest tip calendar time
#' @author Oliver Ratmann
treeannotator.get.viro.timeline<- function(ph, df, df.viro, youngest.tip.ctime, df.treatment=NULL, end.ctime=2013.3)
{
	#df<- df.all; end.ctime=2013.3
	setkey(df, FASTASampleCode)

	if(is.null(df.treatment))		#prepare a single viral load timeline without treatment periods
	{
		tmp				<- df[J(ph$tip.label)][, list(tip=seq_along(ph$tip.label), FASTASampleCode=FASTASampleCode, Patient=Patient, DateDied=DateDied)]
		ans				<- merge(subset(df.viro, select=c(Patient, PosRNA, lRNA)), tmp, by="Patient")	#all.y=1,
		set(ans,NULL,"PosRNA",	youngest.tip.ctime - db.Date2numeric(ans[,PosRNA]))
		set(ans,NULL,"DateDied",	youngest.tip.ctime - db.Date2numeric(ans[,DateDied]))
		set(ans,NULL,"DateDied",	youngest.tip.ctime - ans[,DateDied])
		set(ans,which(is.na(ans[,DateDied])), "DateDied", youngest.tip.ctime - end.ctime)
		ans				<- ans[,list(PosRNA=c(PosRNA,DateDied[1]), lRNA=c(lRNA,tail(lRNA,1), TPeriod=0, NoDrug=0)),by="tip"]
		setkey(ans,tip)
	}
	else							#prepare a single viral load timeline with treatment periods
	{
		#as above except TPeriod=NA
		tmp				<- df[J(ph$tip.label)][, list(tip=seq_along(ph$tip.label), FASTASampleCode=FASTASampleCode, Patient=Patient, DateDied=DateDied)]
		ans				<- merge(subset(df.viro, select=c(Patient, PosRNA, lRNA)), tmp, by="Patient")	#all.y=1,
		set(ans,NULL,"PosRNA",	youngest.tip.ctime - db.Date2numeric(ans[,PosRNA]))
		set(ans,NULL,"DateDied",	youngest.tip.ctime - db.Date2numeric(ans[,DateDied]))
		set(ans,which(is.na(ans[,DateDied])), "DateDied", youngest.tip.ctime - end.ctime)
		ans				<- ans[,list(Patient= rep(Patient[1], length(Patient)+1), PosRNA=c(PosRNA,DateDied[1]), lRNA=c(lRNA,tail(lRNA,1))),by="tip"]
		#prepare subset of df.treatment
		tmp				<- merge(subset(df.treatment, select=c(Patient, StartTime, StopTime, NoDrug)), unique(subset(ans,select=Patient)), all.y=1, by="Patient")
		#modify ans if there is at least one patient with a treatment period
		if( tmp[,any(!is.na(StopTime))] )
		{
			set(tmp,NULL,"StartTime",	youngest.tip.ctime - db.Date2numeric(tmp[,StartTime]))
			set(tmp,NULL,"StopTime",	youngest.tip.ctime - db.Date2numeric(tmp[,StopTime]))
			set(tmp, which( tmp[,StopTime==min(StopTime, na.rm=1)] ), "StopTime", youngest.tip.ctime - 2013.3)
			ans				<- merge(tmp, ans, by="Patient", allow.cartesian=1)
			tmp				<- which(ans[,is.na(StartTime)])
			set(ans,tmp,"StartTime",youngest.tip.ctime - end.ctime)
			set(ans,tmp,"StopTime",youngest.tip.ctime - end.ctime)
			set(ans,tmp,"NoDrug",0L)
			#now have all treatments and viral loads measures together
			ans					<- ans[,	{
						x		<- data.table(StartTime, StopTime, NoDrug, tip, PosRNA,  lRNA)
						#select outside any treatment period
						tmp		<- subset(x, PosRNA>max(StartTime))								#handle no drug before first treatment
						tmp		<- rbind(tmp, subset(x, PosRNA<min(StopTime)) )					#handle no drug after stop treatment
						if(nrow(tmp))
						{
							setkey(tmp, PosRNA)
							tmp		<- unique(tmp)
							set(tmp,NULL,"NoDrug",0L)
							#select only those viral loads within a particular treatment period
							tmp		<- rbind(tmp, subset(x, StartTime>=PosRNA & PosRNA>=StopTime))
						}
						else
							tmp		<- subset(x, StartTime>=PosRNA & PosRNA>=StopTime)
						#assign integer value to different treatment periods
						setkey(tmp, NoDrug, StartTime)
						tmp2	<- subset(unique(tmp), select=c(StartTime, NoDrug))[order(-StartTime, NoDrug)]
						tmp2	<- tmp2[, list(StartTime=StartTime, NoDrug=NoDrug, TPeriod=seq_along(NoDrug))]
						tmp		<- merge(tmp, tmp2, all.x=1, by=c("StartTime", "NoDrug"))
						tmp		<- subset(tmp[order(-PosRNA)], select=c(PosRNA,  lRNA, NoDrug, TPeriod))
						#add endpoints for each treatment period
						tmp2	<- subset(tmp, TPeriod>min(TPeriod))[,  list(PosRNA=PosRNA[1], lRNA=lRNA[1] ),by=TPeriod]
						set(tmp2, NULL, "TPeriod", tmp2[,TPeriod]-1)
						tmp2	<- merge(tmp2, unique(subset(tmp, select=c(TPeriod, NoDrug))), by="TPeriod")
						tmp		<- rbind(tmp, subset(tmp2, select=c(PosRNA,  lRNA, NoDrug, TPeriod)))
						tmp
					},by="tip"]
		}
		else
		{
			ans[,TPeriod:=1]
			ans[,NoDrug:=0]
		}
		ans					<- subset(ans, !is.na(TPeriod))
		ans					<- ans[order(tip, -PosRNA, TPeriod)]
	}
	ans
}
######################################################################################
#' No description yet
#' @author Oliver Ratmann
treeannotator.sero.getnodeheight.range<- function(ph, df, youngest.tip.ctime)
{
	setkey(df, FASTASampleCode)
	ans					<- cbind( data.table(tip=seq_along(ph$tip.label)), subset(df[J(ph$tip.label)], select=c(NegT, AnyPos_T1, DateDied)) )
	set(ans,NULL,"NegT",		youngest.tip.ctime - db.Date2numeric(ans[,NegT]))
	set(ans,NULL,"AnyPos_T1",	youngest.tip.ctime - db.Date2numeric(ans[,AnyPos_T1]))
	set(ans,NULL,"DateDied",	youngest.tip.ctime - db.Date2numeric(ans[,DateDied]))
	set(ans,which(is.na(ans[,DateDied])), "DateDied", 0.)
	ans
}
######################################################################################
#' No description yet
#' @author Oliver Ratmann
treeannotator.plot.seronodeheightrange<- function(ph, ph.seronodeheight, add.yinch= -0.05, width.yinch= 0.1, width.yinch.past.AnyPos_T1= 0.02, col="red")
{
	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
	if (!lastPP$use.edge.length)
		stop("function needs edge length information")
	if (lastPP$type != "phylogram")
		stop("currently only 'type == phylogram' supported")
	if (lastPP$dir!="rightwards")
		stop("currently only rightwards supported")


	tmp				<- max(node.depth.edgelength(ph)) - subset(ph.seronodeheight, select=c(NegT, AnyPos_T1, DateDied))		#from node heights to root heights, assuming root is plotted at 0
	yy.l 			<- lastPP$yy[ph.seronodeheight[,tip]] + yinch(add.yinch)
	rect(tmp[,NegT], yy.l, tmp[,AnyPos_T1], yy.l+yinch(width.yinch), col = col, border=NA)
	rect(tmp[,AnyPos_T1], yy.l, tmp[,DateDied], yy.l+yinch(width.yinch.past.AnyPos_T1), col=col, border=NA)
}
######################################################################################
#' No description yet
#' @author Oliver Ratmann
treeannotator.nodelabels.getnodeheightHPD<- function(ph, youngest.tip.ctime, node.label.hpd.l= 3, node.label.hpd.u= 4)
{
	nodes			<- which(!is.na(ph$node.label))
	node.hpd		<- t( sapply(strsplit(ph$node.label[nodes],'_'), function(x)	as.numeric(x[c(node.label.hpd.l, node.label.hpd.u)])) )
	node.hpd		<- youngest.tip.ctime - node.hpd				#from calendar time to node heights
	data.table(node=nodes, hpd.l=node.hpd[,1], hpd.u=node.hpd[,2])
}
######################################################################################
#' No description yet
#' @author Oliver Ratmann
treeannotator.plot.hpdbars<- function(ph, ph.nodeheighthpds, col="grey75", lwd=4 )
{
	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
	if (!lastPP$use.edge.length)
		stop("function needs edge length information")
	if (lastPP$type != "phylogram")
		stop("currently only 'type == phylogram' supported")
	if (lastPP$dir!="rightwards")
		stop("currently only rightwards supported")
	tmp				<- max(node.depth.edgelength(ph)) - subset(ph.nodeheighthpds, select=c(hpd.l,hpd.u))		#from node heights to root heights, assuming root is plotted at 0
	yy 				<- lastPP$yy[Ntip(ph)+ph.nodeheighthpds[,node]]
	segments(tmp[,hpd.l], yy, tmp[,hpd.u], yy, col = col, lwd = lwd)
}
######################################################################################
#' No description yet
#' @author Oliver Ratmann
treeannotator.get.rates<- function(ph, tip.df, nodelabel.idx.edgewidth=5)
{
	tmp				<- as.numeric( sapply( strsplit(ph$node.label,'_'),function(x)	x[nodelabel.idx.edgewidth] ) )
	rates.df		<- data.table(node=seq_len(Nnode(ph))+Ntip(ph), rate=tmp)
	ans				<- rbind(rates.df, subset(tip.df, select=c(node, rate)))
	#set(ans, which(is.na(ans[,rate])),"rate",mean(ans[,rate],na.rm=1))
	setkey(ans,node)
	ans
}
######################################################################################
#' No description yet
#' @author Oliver Ratmann
treeannotator.get.phy<- function(ph.beast, beastlabel.idx.clu=1, beastlabel.idx.hivs=4, beastlabel.idx.samplecode=5, beastlabel.idx.rate=6, verbose=1, debug=0)
{
	#	get root height for final tree in calendar time
	ph.tip.ctime	<- sapply(ph.beast, function(x) max( as.numeric( sapply(strsplit(x$tip.label,'_'), function(x)	x[beastlabel.idx.hivs] ) ) ))
	ph.root.ctime	<- min( sapply(seq_along(ph.beast), function(i)	ph.tip.ctime[i]-max(node.depth.edgelength(ph.beast[[i]]))	) )
	#	extract tip label information
	tip.df			<- lapply(ph.beast, function(x) treeannotator.tiplabel2df(x, beastlabel.idx.clu=beastlabel.idx.clu, beastlabel.idx.hivs=beastlabel.idx.hivs, beastlabel.idx.samplecode=beastlabel.idx.samplecode, beastlabel.idx.rate=beastlabel.idx.rate) )
	#	prepare cluster subtrees
	clu.subtrees	<- lapply( seq_along(ph.beast), function(i)
			{
				x<- ph.beast[[i]]
				#	convert heights into calendar time and collapse node.label
				tmp					<- ph.tip.ctime[i]
				#subset(x$node.label, node==147, select=c(node, height_median, height_95_HPD_MIN, height_95_HPD_MAX, posterior))
				tmp					<- x$node.label[, list(	node=node,
								height_median=ph.tip.ctime[i]-height_median,
								height_95_HPD_MIN=ph.tip.ctime[i]-height_95_HPD_MAX,
								height_95_HPD_MAX=ph.tip.ctime[i]-height_95_HPD_MIN,
								rate_median=rate_median,
								posterior=posterior)]
				tmp					<- tmp[, list(node.label= paste(posterior,height_median, height_95_HPD_MIN, height_95_HPD_MAX, rate_median, sep='_')),by="node"]
				x$node.label		<- tmp[, node.label]
				x$node.label.format	<- "posterior height_median height_95_HPD_MIN height_95_HPD_MAX rate_median"
				#
				#	extract rooted ExaML clusters
				#
				x$tip.label	<- tip.df[[i]][, FASTASampleCode]
				tmp			<- tip.df[[i]][, list(node=getMRCA(x,FASTASampleCode)),by=cluster]
				clu.subtrees<- lapply(tmp[,node], function(z)
						{
							ans						<- extract.clade(x, z, root.edge= 1, interactive = FALSE)
							ans$root.edge			<- as.numeric(strsplit(ans$node.label[1],'_')[[1]][2])-ph.root.ctime		#reset root edge against root of all runs combined
							ans$node.label.format	<- x$node.label.format
							ans
						})
				names(clu.subtrees)<- tmp[,cluster]
				clu.subtrees
			})
	clu.subtrees	<- eval(parse(text= paste("c(",paste('clu.subtrees[[',seq_along(clu.subtrees),']]', sep='',collapse=','),")",sep='') ))
	if(verbose)	cat(paste("\nFound ExaML clusters in treeannotator files, number of clusters is n=", length(clu.subtrees) ))
	if(debug)
		clu.subtrees	<- lapply(1:3, function(i) clu.subtrees[[i]] )
	#	join all clusters
	cluphy				<- eval(parse(text=paste('clu.subtrees[[',seq_along(clu.subtrees),']]', sep='',collapse='+')))
	if(verbose)	cat(paste("\nFound ExaML clusters in treeannotator files, number of sequences is n=", Ntip(cluphy) ))
	#	retain tip info for those tip labels in cluphy
	tip.df			<- rbindlist(tip.df)
	tip.df			<- merge(data.table(FASTASampleCode=cluphy$tip.label), tip.df, by="FASTASampleCode")
	tip.df			<- cbind(tip.df, node=seq_len(Ntip(cluphy)))

	list(cluphy=cluphy, ph.tip.df=tip.df, ph.tip.ctime=ph.tip.ctime, ph.root.ctime=ph.root.ctime)
}
######################################################################################
#' No description yet
#' @author Oliver Ratmann
treeannotator.get.tmrcas<- function(ph.beast, beastlabel.idx.hivs=4)
{
	#	for each tree, return 	height_median height_95_HPD_MIN height_95_HPD_MAX for the parent of each tip
	ph.trmca	<- lapply(ph.beast, function(x)
			{
				#select heights and convert heights into calendar time
				tip.latest			<- max( as.numeric( sapply(strsplit(x$tip.label,'_'), function(x)	x[beastlabel.idx.hivs] ) ) )
				tip.parents			<- unique( x$edge[x$edge[,2]<=Ntip(x),1] )
				ans					<- x$node.label[J(tip.parents)][, list(node=node, height_median=tip.latest-height_median, height_95_HPD_MIN=tip.latest-height_95_HPD_MAX, height_95_HPD_MAX=tip.latest-height_95_HPD_MIN)]
				tmp					<- sapply(ans[,node], function(z)	x$edge[ x$edge[,1]==z, 2 ] )
				tmp[tmp>Ntip(x)]	<- NA
				tmp					<- t( apply(tmp,2,function(z) x$tip.label[ sort(z,na.last=1) ]) )
				ans[,height_95_diff:= height_95_HPD_MAX-height_95_HPD_MIN]
				ans[,tip1:= tmp[,1]]
				ans[,tip2:= tmp[,2]]
				ans
			})
	ph.trmca	<- rbindlist( ph.trmca )
}
######################################################################################
#' Read a BEAST treeannotator file. The node.label is set to a data.table
#' that contains the SIMMAP annotation for the interior nodes in
#' the newick tree.
#' @author Oliver Ratmann
treeannotator.read<- function(file, add.to.tiplabel=NA, rate.multiplier=NA, round.digit=NA, verbose=1)
{
	require(data.table)
	require(ape)
	ph 			<- read.nexus(file)
	if(verbose)	cat(paste("\nReading BEAST file",file,sep=''))
	X 			<- scan(file = file, what = "", sep = "\n", quiet = TRUE)
	#	read annotations of node in order as they appear in 'file'
	tab			<- X[grep("tree TREE1[[:space:]]+=", X)]
	tab 		<- gsub("tree TREE1[[:space:]]+= \\[&R\\] ", "", tab)
	tab 		<- unlist(strsplit(tab, "\\["))[-1]
	tab 		<- gsub("&|;|\\]", "", tab)
	tab 		<- gsub(":.+$", "", tab)
	tab 		<- lapply(tab, function(x) unlist(strsplit(x, ","))	)
	tab			<- lapply(tab, function(x)
			{
				x			<- gsub('%','',x,fixed=1)
				ind 		<- grep("[{]", x)
				names 		<- gsub("=.+$", "", x[ind])
				x[ind] 		<- gsub("[{]", "", x[ind])
				x[ind] 		<- gsub("=", "_MIN=", x[ind])
				x[ind + 1] 	<- gsub("[}]", "", x[ind + 1])
				x[ind + 1] 	<- paste(paste(names, "MAX=", sep = "_"), x[ind + 1],sep='')
				x
			})
	colnames	<- unique(gsub("=.+$", "", unlist(tab)))
	if(verbose)	cat(paste("\nFound BEAST variables ",paste(colnames,collapse=', '),sep=''))
	tab			<- c(list(paste(colnames,-1,sep='=')), tab)		#rbindlist bug fix
	tab			<- lapply(tab, function(x)
			{
				ans									<- rep(NA, length(colnames))
				names(ans)							<- colnames
				x									<- strsplit(x,'=',fixed=1)
				ans[ sapply(x, function(z) z[1]) ]	<- sapply(x, function(z) z[2])
				ans									<- paste(apply( rbind( names(ans), ans ), 2, function(z) paste(z,collapse='=',sep='')),collapse=',')
				eval(parse(text=paste("data.table(",ans,")",sep='')))
			})
	df.beast	<- rbindlist(tab)[-1,]
	if(!is.na(rate.multiplier))
	{
		tmp	<- colnames(df.beast)[grepl("rate",colnames(df.beast))]
		sapply(tmp, function(x)		set(df.beast, NULL, x, as.numeric(unlist(df.beast[,x,with=F]))*rate.multiplier)			)
	}
	if(!any(is.na(round.digit)))
	{
		if(length(round.digit)!=ncol(df.beast))
			round.digit<- rep(round.digit[1], ncol(df.beast))
		tmp	<- colnames(df.beast)
		sapply(seq_along(tmp), function(i)
				{
					if(class(df.beast[[i]])=="numeric")
						set(df.beast, NULL, tmp[i], round(as.numeric(unlist(df.beast[,tmp[i],with=F])), d=round.digit[i]))
				})
	}

	tmp			<- length(which(df.beast[,!is.na(posterior)]))
	if(verbose)	cat(paste("\nFound annotated nodes, n=", tmp))
	if(verbose)	cat(paste("\nFound annotated tips, n=", nrow(df.beast)-tmp))
	#	determine node index for 'df.beast':
	#
	#	- delete SIMMAP information from 'X'
	LEFT 		<- grep("\\[", X)
	RIGHT 		<- grep("\\]", X)
	if (length(LEFT))
	{
		w 	<- LEFT == RIGHT
		if (any(w))
		{
			s 		<- LEFT[w]
			X[s] 	<- gsub("\\[[^]]*\\]", "", X[s])
		}
		w <- !w
		if(any(w))
		{
			s 		<- LEFT[w]
			X[s] 	<- gsub("\\[.*", "", X[s])
			sb	 	<- RIGHT[w]
			X[sb] 	<- gsub(".*\\]", "", X[sb])
			if(any(s < sb - 1))
				X 	<- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
		}
	}
	#	- read tree block
	endblock 			<- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
	semico 				<- grep(";", X)
	i1 					<- grep("BEGIN TREES;", X, ignore.case = TRUE)
	i2 					<- grep("TRANSLATE", X, ignore.case = TRUE)
	tree 				<- X[(semico[semico > i2][1] + 1):(endblock[endblock > i1][1] - 1)]
	tree 				<- gsub("^.*= *", "", tree)
	tree				<- substr(tree, 1, nchar(tree)-1)
	# 	- For each node, add a dummy node label that is the index in 'df.beast'
	tmp					<- unlist(strsplit(tree, ":"))
	interiorm1			<- which( df.beast[,!is.na(posterior)] )
	tmp					<- sapply(seq_along(tmp), function(i) 			ifelse(i %in% interiorm1, 	paste(tmp[i],i,sep=''),	tmp[i])				)
	tmp					<- paste(paste(tmp, collapse=':'),';',sep='')
	#	- read this newick string and determine the node index in 'df.beast'
	tmp					<- read.tree(text=tmp)
	ph[["node.label"]]	<- cbind(data.table(node=Ntip(ph) + seq_len(Nnode(ph))), df.beast[as.numeric( tmp$node.label ),])
	setkey(ph[["node.label"]], node)

	if(!any(is.na(add.to.tiplabel)))
	{
		if(length(intersect(add.to.tiplabel,colnames(df.beast)))!=length(add.to.tiplabel))
			stop("Cannot find add.to.tiplabel")
		tmp				<- as.matrix(subset( df.beast[as.numeric( tmp$tip.label ),], select=add.to.tiplabel, with=F))
		tmp				<- cbind(ph$tip.label, tmp)
		ph$tip.label	<- apply(tmp, 1, function(x) paste(x,collapse='_'))
	}

	ph
}
