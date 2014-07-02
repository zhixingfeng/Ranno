getGFF <- function(file.name)
{
	gff <- readGff3(file.name)
	gff.attr <- parseGffAttributes(gff)
	
	idx.genes <- which(gff$type == 'gene')
	rl.out <- as.data.frame(matrix(nrow = length(idx.genes), ncol = 8))
	names(rl.out) <- c('ID', 'locus_tag', 'name', 'start', 'end', 'strand', 'product', 'children')
	for (i in 1:length(idx.genes)){
		cur.start <- gff[idx.genes[i],1]
		cur.end <- gff[idx.genes[i],2]				
		idx.children <- which(gff[,1]==cur.start & gff[,2]==cur.end)
		idx.children <- setdiff(idx.children, idx.genes[i])
		
		cur.attr <- gff.attr[[idx.genes[i]]]
		rl.out$ID[i] <- cur.attr['ID']		
		rl.out$locus_tag[i] <- cur.attr['locus_tag']
		rl.out$name[i] <- cur.attr['Name']
		rl.out$start[i] <- cur.start
		rl.out$end[i] <- cur.end
		rl.out$strand[i] <- as.character(gff$strand[idx.genes[i]])
		
		if (length(idx.children)>0){
                	#cur.attr.children <- gff.attr[[idx.children[1]]]
			cur.ID <- sapply(idx.children, function(x, d) d[[x]]['ID']  ,d=gff.attr)
			cur.product <- sapply(idx.children, function(x, d) d[[x]]['product']  ,d=gff.attr)
			rl.out$children[i] <- paste(cur.ID, collapse=';')
			rl.out$product[i] <- paste(cur.product, collapse=';')
		}
		if (i==100*floor(i/100))
			cat(i,'\r')
	}		
	cat(length(idx.genes), '\n')
	
	rl.out	
}

