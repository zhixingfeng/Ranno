library(genomeIntervals)
gff3torefFlat <- function(gff.file, outfile)
{
	gff <- 	readGff3(gff.file)
	gff.attr <- parseGffAttributes(gff)
	ref.names <- sapply(gff.attr, function(x) {if(length(x[names(x)=='Name'])>0) x[names(x)=='Name'] else 'unknown'} )
	ref.names <- as.character(ref.names)
	idx.names <- which(ref.names!='unknown')	

	refFlat <- cbind(ref.names[idx.names], ref.names[idx.names], as.character(gff$seq_name)[idx.names],
			as.character(gff$strand)[idx.names], gff[idx.names,1], gff[idx.names,2],
			gff[idx.names,1], gff[idx.names,2], rep(1,length(idx.names)), 
			paste(gff[idx.names,1],',',sep=''), paste(gff[idx.names,2],',',sep=''))
	refFlat <- as.data.frame(refFlat)
	names(refFlat) <- c('locus_tag', 'name', 'seq', 'strand', 'start', 'end', 
			'start_cds', 'end_cds', 'exon_num', 'start_exons', 'end_exons')	
	write.table(refFlat, file=outfile, quote=FALSE, sep='\t', row.names=FALSE)
	
}


