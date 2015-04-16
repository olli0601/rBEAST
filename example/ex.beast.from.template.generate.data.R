require(rBEAST)

infile.seq				<- '~/git/rBEAST/inst/extra/sim_150414a_seq.R'
outfile.seq				<- '~/git/rBEAST/inst/extra/sim_150416_seq.R'
outfile.tree			<- '~/git/rBEAST/inst/extra/sim_150416_njtree.R'

#	load data
load(infile.seq)		
#	remove singletons
df.seq		<- merge(df.seq, subset( df.seq[, list(CLU_N= length(TAXON_ID)), by='CLU_ID'], CLU_N>1, CLU_ID ), by='CLU_ID')	
#	create NJ tree
file			<- system.file(package="rPANGEAHIVsim", "misc",'PANGEA_SSAfg_HXB2outgroup.R')
cat(paste('\nLoading outgroup seq from file', file))
load(file)		#expect "outgroup.seq.gag" "outgroup.seq.pol" "outgroup.seq.env"
df.seq.nj		<- copy(df.seq)
#	concatenate sequences	
tmp				<- tolower(do.call('rbind',strsplit(df.seq.nj[, paste(GAG,POL,ENV,sep='')],'')))		
rownames(tmp)	<- df.seq.nj[, paste(CLU_ID,'|',TAXON_NAME,sep='')]	
seq				<- as.DNAbin(tmp)
tmp				<- cbind(outgroup.seq.gag[,1:nchar(df.seq.nj[1, GAG])], outgroup.seq.pol, outgroup.seq.env)
seq				<- rbind(seq,tmp)	
seq.ph			<- nj(dist.dna(seq, model='raw'))		
tmp				<- which(seq.ph$tip.label=="HXB2")
seq.ph			<- reroot(seq.ph, tmp, seq.ph$edge.length[which(seq.ph$edge[,2]==tmp)])
seq.ph			<- ladderize(seq.ph)
seq.ph			<- drop.tip(seq.ph,'HXB2')
ph				<- seq.ph
#	cluster based on distance 	
dist.brl		<- hivc.clu.brdist.stats(seq.ph, eval.dist.btw="leaf", stat.fun=hivc.clu.min.transmission.cascade)	 
thresh.brl		<- 0.07
clustering		<- hivc.clu.clusterbythresh(seq.ph, thresh.brl=thresh.brl, dist.brl=dist.brl, retval="all")	
tmp				<- data.table(TAXON_NAME=seq.ph$tip.label, NJ_CLU_ID=clustering$clu.mem[ seq_len(Ntip(seq.ph))] )
#	merge
set(tmp, NULL, 'TAXON_NAME', tmp[, gsub('^[0-9]+\\|','',TAXON_NAME)])
df.seq			<- merge(df.seq, tmp, by='TAXON_NAME')
df.seq			<- subset(df.seq, !is.na(NJ_CLU_ID))
set(df.seq, NULL, 'CLU_ID', NULL)
setnames(df.seq, 'NJ_CLU_ID', 'CLU_ID')
#	save
save(df.seq, file=outfile.seq)
ph$tip.label	<- gsub('^[0-9]+\\|','',ph$tip.label)
save(ph, file=outfile.tree)
