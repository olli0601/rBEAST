require(rBEAST)
file				<- '~/git/rBEAST/inst/extra/exclupool_150414.R'
load(file)	
#expect data.table 'df' with (at least) columns
#	"TAXON_ID" 			taxon id
#	"CLU_ID" 			phylogenetic cluster ID of sequence, NA if not in any cluster
#	"SAMPLINGTIME"		sampling time of taxon

#	option 'by.samplingtime'
#	defines cluster sampling time as the earliest taxon sampling time
#	ranks clusters by cluster sampling time
#	pools clusters by evenly spaced rank, ie for 3 clusters, ranks 1, 4, 7, ... are grouped into one pool
#
#	pool.ntip: guide to the number of sequences in each pool. This is typically not met exactly, because sequences contain differing numbers of sequences.
df.clupool			<- beast.pool.clusters(df, how='by.samplingtime', pool.ntip=150, verbose=1)

#	option 'by.samplingtime.alwaysincludebeforeyear'
#	as option 'by.samplingtime'
#	and in addition includes all early cluster with cluster sampling time before 'pool.includealwaysbeforeyear' into each pool
#
#	pool.ntip: guide to the number of sequences in each pool. This is typically not met exactly, because sequences contain differing numbers of sequences.
#	pool.includealwaysbeforeyear: year before which clusters are always included.
df.clupool			<- beast.pool.clusters(df, how='by.samplingtime.alwaysincludebeforeyear', pool.ntip=150, pool.includealwaysbeforeyear=1996, verbose=1)

#	option 'by.required.seq.per.samplingperiod'

#	pool.ntip.guide: guide to the number of sequences in each pool. This is typically not met exactly, because sequences contain differing numbers of sequences.
#	pool.breaks: sampling periods defined in terms of height (max sampling time - sampling time)
#	pool.ntip.min: minimum number of sequences required per sampling period; NA if none required
df.clupool			<- beast.pool.clusters(df, how='by.required.seq.per.samplingperiod', pool.ntip.guide=150, pool.ntip.min=c(50, 70, 70, NA, NA), pool.breaks=c(0, 1.596, 3.596, 5.596, 9.596, Inf), verbose=1)
