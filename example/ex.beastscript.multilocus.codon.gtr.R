require(rBEAST)
tree.label.sep			<- '|'
tree.label.idx.idpop	<- 1
tree.label.idx.ctime	<- 4 	
#select					<- 'grid-cseq3'		#generates large file to illustrate rBEAST functionality
select					<- 'grid-mseq400'	#generates relatively small file for testing
infile.seq				<- '~/git/rBEAST/inst/extra/sim_150414a_seq.R'
infile.starttrees		<- '~/git/rBEAST/inst/extra/sim_150414a_starttrees.newick'
#
#	load data
#
file			<- infile.seq	
cat(paste('\nLoading file', file))
load(file)		
#	expect data.table "df.seq" with columns
#	"LABEL" 	taxon name
#	"GAG"   	gag sequence as character string
#	"POL"   	pol sequence as character string
#	"ENV"   	env sequence as character string
#	"IDCLU" 	phylogenetic cluster ID of sequence, NA if not in any cluster
#	"IDPOP"		individual ID
set( df.seq, NULL, 'IDCLU', df.seq[, as.integer(IDCLU)] )

#	'size' option: eg how many sequences to select
thresh.NSEQ		<- as.numeric(substring(regmatches(select, regexpr('seq[0-9]+',select)), 4))
stopifnot(!is.na(thresh.NSEQ))
#	'how to select' options
#	pick phylogenetic clusters randomly
if(grepl('clrndseq',select))
{		
	seq.select		<- df.seq[, list(CLU_N=length(LABEL)), by='IDCLU']			
	setkey(seq.select, CLU_N, IDCLU)
	seq.select		<- subset(unique(seq.select), CLU_N>1)
	seq.select[, DUMMY:= sample(nrow(seq.select), nrow(seq.select))]
	setkey(seq.select, DUMMY)
	seq.select[, CLU_CN:= seq.select[, cumsum(CLU_N)]]
	seq.select		<- seq.select[seq_len( seq.select[, which(CLU_CN>=thresh.NSEQ)[1]] ), ]			
	seq.select		<- merge( df.seq, subset(seq.select, select=IDCLU), by='IDCLU' )				
}
#	or pick sequences from largest phylogenetic clusters
if(grepl('mseq',select))
{	
	seq.select		<- df.seq[, list(CLU_N=-length(LABEL)), by='IDCLU']			
	setkey(seq.select, CLU_N, IDCLU)
	seq.select		<- subset(unique(seq.select), CLU_N< -1)	
	seq.select[, CLU_CN:= seq.select[, cumsum(-CLU_N)]]
	seq.select		<- seq.select[seq_len( seq.select[, which(CLU_CN>=thresh.NSEQ)[1]] ), ]			
	seq.select		<- merge( df.seq, subset(seq.select, select=IDCLU), by='IDCLU' )									
}
#	or pick sequences from smallest phylogenetic clusters
if(grepl('clsmseq',select))
{
	seq.select		<- df.seq[, list(CLU_N=length(LABEL)), by='IDCLU']			
	setkey(seq.select, CLU_N, IDCLU)
	seq.select		<- subset(unique(seq.select), CLU_N>1)	
	seq.select[, CLU_CN:= seq.select[, cumsum(CLU_N)]]
	seq.select		<- seq.select[seq_len( seq.select[, which(CLU_CN>=thresh.NSEQ)[1]] ), ]			
	seq.select		<- merge( df.seq, subset(seq.select, select=IDCLU), by='IDCLU' )					
}
#	or pick sequences of all phylogenetic clusters of size >= x
if(grepl('cseq',select))
{				
	seq.select		<- df.seq[, list(CLU_N=length(LABEL)), by='IDCLU']			
	setkey(seq.select, CLU_N, IDCLU)
	seq.select		<- subset(unique(seq.select), CLU_N>=thresh.NSEQ)	
	seq.select		<- merge( df.seq, subset(seq.select, select=IDCLU), by='IDCLU' )					
}
cat(paste('\nFound clusters, n=', seq.select[, length(unique(IDCLU))])) 
cat(paste('\nFound sequences, n=', seq.select[, length(unique(IDPOP))]))
#
#	read starting trees for each cluster phylogeny
#		
phd					<- read.tree(infile.starttrees)
tmp2				<- data.table(IDPOP= sapply(phd, function(x) x$tip.label[1]), IDX=seq_along(phd))
set( tmp2, NULL, 'IDPOP', tmp2[, as.integer(substring(sapply(strsplit(IDPOP, tree.label.sep, fixed=TRUE),'[[',tree.label.idx.idpop),7)) ] )
tmp2				<- merge(subset(seq.select, select=c(IDPOP, IDCLU)), tmp2, by='IDPOP')			
phd					<- lapply(tmp2[,IDX], function(i) phd[[i]] )
names(phd)			<- tmp2[, IDCLU]
# 	plot starting trees
\dontrun{
phd.plot			<- eval(parse(text=paste('phd[[',seq_along(phd),']]', sep='',collapse='+')))			
phd.plot			<- ladderize(phd.plot)		
pdf(file=gsub('.newick',paste('_',select,'.pdf',sep=''),infile.starttrees), w=10, h=Ntip(phd.plot)*0.1)
plot(phd.plot, show.tip=TRUE, cex=0.5)
dev.off()	
}
#
#	create XML file for POL
#			
seq.select.pol	<- subset(seq.select, select=c("IDCLU", "IDPOP", "LABEL", "POL" ))
setnames(seq.select.pol, 'POL', 'SEQ')
file.name		<- gsub('_seq.R',paste('_HKY_fixedtree_',select,'.xml',sep=''),infile.seq)
bxml			<- beastscript.multilocus.codon.gtr( file.name, seq.select.pol, phd, verbose=1 )
#	write to file
\dontrun{
cat(paste("\nwrite xml file to",file.name))
saveXML(bxml, file=file.name)	
}
