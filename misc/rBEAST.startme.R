#! /Library/Frameworks/R.framework/Versions/3.1/Resources/bin/Rscript
##	first line in shell script starts with #! and start Rscript
##	CHANGE  as needed
##! /Library/Frameworks/R.framework/Versions/2.11/Resources/bin/Rscript
##! /Library/Frameworks/R.framework/Versions/2.15/Resources/bin/Rscript
##! /apps/R/2.15/lib64/R/bin/Rscript
##! /usr/bin/Rscript
##! /Library/Frameworks/R.framework/Versions/2.11/Resources/bin/Rscript
###############################################################################
#
#
# Author: oliver ratmann
# file: hivclu.startme.R
#
# usage from R:
#> setwd("/Users/Oliver/git/hivclust/pkg"); source("misc/hivclu.startme.R")
# usage from bash:
#> cd /Users/Oliver/git/hivclust/pkg
#> misc/hivclu.startme.R --help
#
#
# Installation instructions:
#	1 make sure the first line in this file points to your Rscript
#	2.1	create a directory CODE.HOME and set the CODE.HOME variable below to this path
#	2.2	create CODE.HOME/src_tipclust and copy all the R files into this directory
#	2.3 create a directory HOME and set the HOME variable below to this path
#
#
###############################################################################
args <- commandArgs()
if(!any(args=='--args'))
	args<- vector("numeric",0)
if(any(args=='--args'))
	args<- args[-(1:match("--args", args)) ]

require(data.table)
require(grid)
require(reshape2)
require(ggplot2)
require(hivclust)


CODE.HOME	<<- "/Users/Oliver/git/hivclust/pkg"
#CODE.HOME	<<- "/Users/Stephane/Phylogenetics/github/hivclust/pkg"
#CODE.HOME	<<- "/work/or105/libs/hivclust/pkg"
INST		<<- paste(CODE.HOME,"inst",sep='/')
HOME		<<- "/Users/Oliver/duke/2013_HIV_NL/ATHENA_2013"
#HOME		<<- "/Users/Oliver/duke/2013_HIV_Hue/UKCA_1309"
#HOME		<<- "/work/or105/UKCA_1309"
#HOME		<<- "/Users/Stephane/Desktop/CASCADE_phylo/hivclust"
#HOME		<<- "/home/koelle/or7/phylody"
#HOME		<<- "/work/or105/ATHENA_2013"
DATA		<<- paste(HOME,"data",sep='/')
HIVC.DEBUG	<<- 0
LIB.LOC		<<- NULL
#LIB.LOC		<<- paste(CODE.HOME,"../",sep='')
EPS			<<- 1e-12

#default.fun		<- "my.make.documentation"
#default.fun		<- "project.hivc.clustalo"
#default.fun		<- "project.hivc.examl"
#default.fun	<- "hivc.prog.get.geneticdist"
#default.fun 	<- "hivc.prog.get.allseq"
#default.fun	<- "project.hivc.test"
#default.fun		<- "project.hivc.Excel2dataframe"
#default.fun	<- "project.hivc.check"
#default.fun		<- "project.hivc.collectpatientdata"
#default.fun	<- "project.hivc.get.geneticdist.from.sdc"
#default.fun		<- "project.bezemer2013a.figs"
#default.fun		<-"project.bezemer2013b.rates"
#default.fun		<- "project.hivc.clustering"
#default.fun	<- "project.hivc.beast"
#default.fun	<- "hivc.prog.get.clustering.precompute"
#default.fun	<- "project.gccontent"
#default.fun 	<- "hivc.proj.pipeline"
#default.fun 	<- "hivc.prog.remove.resistancemut"
#default.fun 	<- "hivc.prog.BEAST.poolrunxml"
#default.fun 	<- "hivc.prog.BEAST.evalpoolrun"
#default.fun 	<- "hivc.pipeline.recombination"
#default.fun 	<- "hivc.pipeline.ExaML"
#default.fun		<- "hivc.pipeline.clustering"
#default.fun 	<- "hivc.pipeline.BEAST"
#default.fun		<- "hivc.pipeline.BEASTout"
default.fun		<- "hivc.pipeline.props_univariate"
#default.fun		<- "project.athena.Fisheretal.numbers"
#default.fun		<- 'project.bezemer2013a.figs.v131023_DB'
#default.fun		<- "project.athena.Fisheretal.exact.repro"
#default.fun		<- "hivc.pipeline.various"
###############################################################################
#if(length(args) && !is.loaded("tipc_tabulate_after_sample"))
#{
#	file<- paste(CODE.HOME,"src_tipclust",paste("libtipcr",.Platform$dynlib.ext,sep=''),sep='/')
#	cat(paste("\nloading",file,'\n',sep=' '))
#	dyn.load(file)
#}
#cat(paste("is.loaded('tipcr')->",is.loaded("tipc_tabulate_after_sample"),'\n'))
###############################################################################
function.list<-c(list.files(path= paste(CODE.HOME,"R",sep='/'), pattern = ".R$", all.files = FALSE,
		full.names = TRUE, recursive = FALSE),paste(CODE.HOME,"misc","hivclu.prjcts.R",sep='/'),paste(CODE.HOME,"misc","gcproject2013.R",sep='/'))
sapply(function.list,function(x){ print(x); source(x,echo=FALSE,print.eval=FALSE, verbose=FALSE) })
###############################################################################	
if(0)
{
	PR.CLUSTALO<<- paste(INST,PR.CLUSTALO,sep='/')	
}
###############################################################################
.ls.objects <- function (pos = 1, pattern, order.by,
		decreasing=FALSE, head=FALSE, n=5) {
	napply <- function(names, fn) sapply(names, function(x)
					fn(get(x, pos = pos)))
	names <- ls(pos = pos, pattern = pattern)
	obj.class <- napply(names, function(x) as.character(class(x))[1])
	obj.mode <- napply(names, mode)
	obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
	obj.prettysize <- napply(names, function(x) {
				capture.output(print(object.size(x), units = "auto")) })
	obj.size <- napply(names, object.size)
	obj.dim <- t(napply(names, function(x)
						as.numeric(dim(x))[1:2]))
	vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
	obj.dim[vec, 1] <- napply(names, length)[vec]
	out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
	names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
	if (!missing(order.by))
		out <- out[order(out[[order.by]], decreasing=decreasing), ]
	if (head)
		out <- head(out, n)
	out
}

# from: http://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session
lsos <- function(..., n=10) {
	.ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

my.mkdir<-function(root,data.name)
{
	if(length(dir(root,pattern=paste('^',data.name,'$',sep='')))==0)
		system(paste("mkdir ",paste(root,data.name,sep='/'),sep=''))
}

my.make.documentation<- function()
{
	require(roxygen2)		
	roxygenize("~/git/rBEAST")
}

#index of row i (i=1..n) and col j (j=1..n) of nxn matrix in a distance object (as.vector of lower triangular matrix)
my.lower.tri.index<- function(n, i, j)
{
	n*(j-1) - j*(j-1)/2 + i-j
}

my.fade.col<-function(col,alpha=0.5)
{
	return(rgb(col2rgb(col)[1]/255,col2rgb(col)[2]/255,col2rgb(col)[3]/255,alpha))
}

my.sample <- function(x, ...)
{
	x[sample.int(length(x), ...)]	
} 

my.dumpframes<- function()
{
	geterrmessage()
	dump.frames()
	cat(paste("\nmy.dumpframes dump 'last.dump' to file",paste(DATA,paste("debug_",paste(strsplit(date(),' ')[[1]],collapse='_'),".rda\n",sep=''),sep='/')))
	save(last.dump, file=paste(DATA,paste("debug_",paste(strsplit(date(),' ')[[1]],collapse='_'),".rda",sep=''),sep='/'))
	q()
}

dziexp	<-function(x, z, rate, log=FALSE)
{
	ans					<- numeric(length(x))
	tmp					<- x==0.
	ans[tmp]			<- z+(1-z)*dexp(x[tmp], rate=rate)
	ans[!tmp]			<- (1-z)*dexp(x[!tmp], rate=rate)
	if(log)
		ans				<- log(ans)
	ans
}

pziexp	<-function(x, z, rate, log=FALSE)
{
	ans					<- numeric(length(x))
	tmp					<- x==0.
	ans[tmp]			<- z
	ans[!tmp]			<- z+(1-z)*pexp(x[!tmp], rate=rate)
	if(log)
		ans				<- log(ans)
	ans
}

my.aggregate<- function(x, bins)
{
	bins	<- cbind(bins[-length(bins)],bins[-1])			
	ans		<- numeric(length(x))
	for(i in seq_len(nrow(bins)))
		ans[which(x<bins[i,2] & x>=bins[i,1])]<- i
	factor(ans, levels=seq_len(nrow(bins)), labels=apply(bins,1,function(z) paste(z,collapse=',')))		
}

my.intersect.n<- function( list.to.intersect )
{
	if(length(list.to.intersect)==2)
		ans<- intersect(list.to.intersect[[1]],list.to.intersect[[2]])
	else
	{
		ans	<- table(unlist(list.to.intersect))														
		ans <- as.numeric(attr(ans,"dimnames")[[1]])[ans==length(list.to.intersect)]	
		
	}	
	ans
}

print.v<- function(x,cut=3,digits=4,prefix= "simu_",print.char= TRUE, as.R= TRUE)
{
	if(as.R)
	{
		tmp<- paste("c(",paste(c(x,recursive=T),collapse=',',sep=''),')',sep='')
		if(!is.null(names(x)))
			tmp<- paste("{tmp<-", tmp, "; names(tmp)<- ", paste('c("',paste(c(names(x),recursive=T),collapse='", "',sep=''),'")',sep=''), "; tmp}", sep= '', collapse= '')
	}
	else
	{
		if(!is.null(names(x)))
		{
			m<- matrix(NA,nrow=2,ncol=length(x))
			m[1,]<- substr(names(x),1,cut)
			m[2,]<- round( x, digits=digits )
			if(cut==0)		m<- m[2,]
			tmp<- gsub('.',',',paste(prefix,paste(as.vector(m), collapse='_',sep=''),sep=''),fixed=T)
		}
		else
			tmp<- gsub('.',',',paste(prefix,paste(round( x, digits=digits ), collapse='_',sep=''),sep=''),fixed=T)
	}
	if(print.char) print(tmp)
	tmp
}

my.barplot.table<- function(x, xh, xlab, ylab, cols, x.baseline=0, xax=as.numeric(colnames(x)), legend.loc=NULL)
{
	z		<- rbind( rep(x.baseline, ncol(x)), apply(x, 2, cumsum) )
	xlim	<- range(xax)
	ylim	<- range(z)
	
	plot(1,1,type='n',bty='n',xlab=xlab,ylab=ylab, xlim=xlim, ylim=ylim)
	lapply(2:nrow(z), function(i)
			{
				lapply(seq_len(ncol(z)), function(j)
						{
							polygon(	c(xax[j]-xh,xax[j]+xh,xax[j]+xh,xax[j]-xh),		c(rep(z[i-1,j],2),rep(z[i,j],2)), border=NA, col=cols[i-1] 	)
						})
			} )
	if(!is.null(legend.loc))
		legend(legend.loc, fill= cols, legend=rownames(x), border=NA, bty='n')
}

plot.2D.dens<- function(x,y,xlab,ylab,xlim=NA,ylim=NA,nbin=NA,width.infl=2,n.hists=5,method="gauss", palette= "topo", persp.theta= -30, persp.phi= 30, zero.abline=TRUE, ...)
{
	if(!method%in%c("gauss","ash","persp"))	stop("plot.2D.dens: exception 1a")
	if(!palette%in%c("topo","heat","gray"))	stop("plot.2D.dens: exception 1b")
	switch(method,
			gauss=
					{
						require(KernSmooth)
						require(fields)
						x.bw<- width.infl*diff(summary(x)[c(2,5)])
						y.bw<- width.infl*diff(summary(y)[c(2,5)])
						if(!x.bw) x.bw<- EPS
						if(!y.bw) y.bw<- EPS
						if(any(is.na(xlim)))	xlim<- range(x)+c(-1.5,1.5)*x.bw
						if(any(is.na(ylim)))	ylim<- range(y)+c(-1.5,1.5)*y.bw
						dens <- bkde2D(cbind(x, y), range.x=list(xlim,ylim),bandwidth=c(x.bw,y.bw))
						contour(dens$x1, dens$x2, dens$fhat,xlab=xlab,ylab=ylab)
						if(zero.abline) abline(v=0,col="black",lty=3,lwd=1.5)
						if(zero.abline) abline(h=0,col="black",lty=3,lwd=1.5)
					},
			ash={
				require(ash)
				if(any(is.na(xlim))) xlim<- range(x)*1.05
				if(any(is.na(ylim))) ylim<- range(y)*1.05
				if(any(is.na(nbin))) nbin<- 2*c(nclass.Sturges(x),nclass.Sturges(y))
				bins<- bin2(cbind(x, y), ab=rbind(xlim,ylim),nbin=nbin)
				f <- ash2(bins,rep(n.hists,2))
				#image(f$x,f$y,f$z, col=rainbow(50,start= 4/6,end=0),xlab=xlab,ylab=ylab)
				if(palette=="topo")					image(f$x,f$y,f$z, col=tail( topo.colors(trunc(50*1.4)), 50 ),xlab=xlab,ylab=ylab,...)
				else if(palette=="gray")			image(f$x,f$y,f$z, col=head( rev(gray(seq(0,.95,len=trunc(50*1.4)))), 50),xlab=xlab,ylab=ylab,...)					
				else								image(f$x,f$y,f$z, col=heat.colors( 50 ),xlab=xlab,ylab=ylab,...)
				contour(f$x,f$y,f$z,add=TRUE, nlevels= 5)
				if(zero.abline) abline(v=0,col="black",lty=2,lwd=2)
				if(zero.abline) abline(h=0,col="black",lty=2,lwd=2)
			},
			persp={
				require(ash)
				require(MASS)
				if(any(is.na(xlim))) xlim<- range(x)*1.05
				if(any(is.na(ylim))) ylim<- range(y)*1.05
				bins<- bin2(cbind(x, y), ab=rbind(xlim,ylim),nbin=2*c(nclass.Sturges(x),nclass.Sturges(y)))
				f <- ash2(bins,rep(n.hists,2))
				
				nrz <- nrow(f$z)
				ncz <- ncol(f$z)
				col<- tail(topo.colors(trunc(1.4 * 50)),50)
				fcol      <- col[trunc(f$z / max(f$z)*(50-1))+1]
				dim(fcol) <- c(nrz,ncz)
				fcol      <- fcol[-nrz,-ncz]
				par(mar=c(1/2,1/2,1/2,1/2))
				persp(x=f$x,y=f$y,z=f$z,col= fcol,theta=persp.theta,phi=persp.phi,xlab=xlab,ylab=ylab,zlab='', ticktype= "detailed" )
			})
}

plot.persplocfit<- function(x, pv, theta= 30, phi= 20, palette= "gray",tcl=-0.05,...)
{	
	d <- x$mi["d"]
	ev <- x$mi["ev"]
	where <- "grid"
	
	pv <- match(pv, x$vnames)
	tv <- (1:d)[-pv]
	vrs <- c(pv, tv)
	if (any(duplicated(vrs))) 
		warning("Duplicated variables in pv, tv")
	if (any((vrs <= 0) | (vrs > d))) 
		stop("Invalid variable numbers in pv, tv")
	m <- ifelse(d == 1, 100, 40) 
	m <- rep(m, d)
	m[tv] <- mtv<- 6		
	xl <- x$box
	marg <- lfmarg(xl, m)
	pred <- locfit:::preplot.locfit(x, marg, band = "none", tr = NULL, what = "coef", get.data = 0, f3d = 0)
	z<- matrix(pred$fit, nrow = length(marg[[1]]))
	nbcol <- 100
	if(palette=="gray")	
		color<- head( rev(gray(seq(0,0.95,len=trunc(nbcol*1.4)))), nbcol)
	else
	{
		jet.colors <- colorRampPalette( c("blue", "green") )
		color <- jet.colors(nbcol)
	}
	# Compute the z-value at the facet centres
	nrz <- nrow(z)
	ncz <- ncol(z)
	zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
	# Recode facet z-values into color indices
	facetcol <- cut(zfacet, nbcol)		
	par(mar=c(0,2.5,0,0), tcl=tcl)
	pmat<- persp(marg[[1]], marg[[2]], z, zlim= range(z)*1.1, col = color[facetcol], shade= 0.1, border=NA, ticktype = "detailed", ltheta = 120,theta = theta, phi = phi, expand = 0.75, box=1, ... )
	
	list(pmat=pmat, x= marg[[1]], y= marg[[2]], z= z)
}

###############################################################################
my.mkdir(HOME,"data")
my.mkdir(HOME,"pdf")
my.mkdir(HOME,"script")
argv<- list()
if(length(args))
{
	tmp<- na.omit(sapply(args,function(arg)
					{
						switch(substr(arg,2,4),
								exe= return(substr(arg,6,nchar(arg))),
								NA)
					}))
	if(length(tmp)!=0)
	{
		if(length(tmp)>1) stop("hivclu.startme.R: duplicate -exe")
		else default.fun<- switch(tmp[1],
					MAKE.DOCUMENTATION		= "my.make.documentation",
					MALIGN					= "project.hivc.clustalo",
					FIRSTSEQ				= "hivc.prog.get.firstseq",
					BOOTSTRAPSEQ			= "hivc.prog.get.bootstrapseq",
					GENDISTMAT				= "hivc.prog.get.geneticdist",
					PRECLUST				= "hivc.prog.get.clustering.precompute",
					CLUST					= "hivc.prog.get.clustering",
					CLUSTTPTN				= "hivc.prog.get.clustering.TPTN",
					CLUSTMSM				= "hivc.prog.get.clustering.MSM",
					BEASTPOOLRUN			= "hivc.prog.BEAST.generate.xml",
					BEASTEVALRUN			= "hivc.prog.BEAST.evalpoolrun",
					BEAST.READNEXUS			= "hivc.prog.BEAST.read.nexus.and.stats",
					BEAST2.PIPE.CLUTREES	= "hivc.pipeline.BEASTout.get.cluster.trees",
					BEAST2.CLUTREES			= "hivc.prog.BEAST2.get.cluster.trees",
					BEAST2.CLUPOSTERIOR		= "hivc.prog.BEAST2.process.cluster.trees",					
					BEAST2.PLOTCLUTREES		= "hivc.prog.BEAST2.plot.cluster.trees",
					PH.DISTTIPS				= "hivc.prog.get.dist.tips",
					RECOMB.PROCESS3SEQOUT	= "hivc.prog.recombination.process.3SEQ.output",
					RECOMB.CHECKCANDIDATES	= "hivc.prog.recombination.check.candidates",
					RECOMB.PLOTINCONGRUENCE	= "hivc.prog.recombination.plot.incongruence",
					PROPS.ESTIMATE			= "hivc.prog.props_univariate",
					BETAREG.NUMBERS			= "project.athena.Fisheretal.numbers",
					VARIOUS					= "project.hivc.examl.median.brl"
					)
	}
	tmp<- na.omit(sapply(args,function(arg)
					{
						switch(substr(arg,2,6),
								debug= 1,
								NA)
					}))		
	if(length(tmp)!=0)	HIVC.DEBUG<<- tmp[1]	
	argv<<- args
}
###############################################################################
stop()
if(HIVC.DEBUG)	options(error= my.dumpframes)	
cat(paste("\nhivclu.startme.R: ",ifelse(HIVC.DEBUG,"debug",""),"call",default.fun,"\n"))
do.call(default.fun,list()) 	
cat("\nhivclu: ",ifelse(HIVC.DEBUG,"debug","")," end\n")
