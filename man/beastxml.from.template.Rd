% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/script.beast.R
\name{beastxml.from.template}
\alias{beastxml.from.template}
\title{Create BEAST XML file from template}
\usage{
beastxml.from.template(btemplate, df, start.tree = NULL,
  xml.taxonsets4clusters = 1, xml.monophyly4clusters = 1, verbose = 1)
}
\arguments{
\item{btemplate}{BEAST XML template or file name of BEAST XML template}

\item{df}{Taxon data.table}

\item{start.tree}{Starting tree containing the taxa in the taxon data.table (optional)}

\item{xml.taxonsets4clusters}{If true, add taxon sets for each phylogenetic cluster in the taxon data.table}

\item{xml.monophyly4clusters}{If true, add monophyly contraints to each phylogenetic cluster in the taxon data.table}

\item{verbose}{If true, run in verbose mode}
}
\description{
See Examples.
}
\examples{
require(rBEAST)
infile.seq				<- system.file(package="rBEAST", 'extra','sim_150416_seq.R')
infile.tree				<- system.file(package="rBEAST", 'extra','sim_150416_njtree.R')
infile.xml.template		<- system.file(package="rBEAST", 'xml_templates','BEAST_template_um192rhU2080.xml')
outdir.xml				<- getwd()
outfile.xml				<- 'sim_150416_um192rhU2080.xml'

#	load data
load(infile.seq)		
#	expect data.table "df.seq" with columns
#	"TAXON_ID"		taxon ID
#	"TAXON_NAME" 	taxon name
#	"SAMPLINGTIME"	taxon sampling time
#	"GAG"   		gag sequence as character string
#	"POL"   		pol sequence as character string
#	"ENV"   		env sequence as character string
#	"CLU_ID" 		phylogenetic cluster ID of sequence

#	use POL sequences
setnames(df.seq,'POL','SEQ')
#	load tree
load(infile.tree)
#	expect ape tree 'ph'

#	illustrate automatic generation of multiple XML files, each for one pool of sequences
#	with monophyly constraints on phylogenetic clusters.
#	with this option, an appropriate starting tree should be specified
#
#	group clusters into separate pools
df.clupool	<- beast.pool.clusters(df.seq, how='by.samplingtime', pool.ntip=300, verbose=1)
#	generate XML file from template
for( pid in seq_along(df.clupool) )
{
	df			<- df.clupool[[pid]]	
	#	get newick start tree for all sequences in pool
	start.tree	<- beast.get.startingtree(ph, df, verbose=1)				
	#	get xml file from template			
	bxml		<- beastxml.from.template(infile.xml.template, df, start.tree=start.tree, xml.taxonsets4clusters=1, xml.monophyly4clusters=1,  verbose=1)
	#	adjust xml file depending on the template that is used
	beast.adjust.gmrfSkyrideLikelihood(bxml, df, verbose=1)				
	beast.adjust.prior.rootheight(bxml, df, rootHeight.idref='treeModel.rootHeight', verbose=1)
	beast.adjust.mcmc(bxml, 
			mcmc.chainLength=beast.mcmc.chainLength, 
			mcmc.logEvery=beast.mcmc.chainLength/1e4, 
			mcmc.outfile=gsub('.xml', paste('_pool',pid,'.xml',sep=''),outfile.xml), 
			verbose=1)
	#	write to file
	saveXML(bxml, file=paste(outdir.xml,gsub('.xml', paste('_pool',pid,'.xml',sep=''),outfile.xml),sep='/'))	
}
}
\author{
Oliver Ratmann
}

