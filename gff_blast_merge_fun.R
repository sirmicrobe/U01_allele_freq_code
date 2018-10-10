gff <-"/Users/chrismarshall/Documents/Pitt/Cooper_Lab/Misc_projects/Pig_wound/nanopore/first_nanopore_run/05_prokka/05_npore.gff"
blast <- "/Users/chrismarshall/Documents/Pitt/Cooper_Lab/Misc_projects/Pig_wound/nanopore/first_nanopore_run/05_blastresult.txt"

contig_match<- function(gff,blast,name){
  library(tidyr)
  gff_df <- read.delim(file=gff,sep="\t",header=F,comment.char="#")
  blast_df <- read.delim(file = blast,sep="\t",header=F)
  colnames(blast_df) <- c("ID","pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle", "staxids", "sseqid")
  gff_df <- gff_df[!grepl("##", gff_df[,1],fixed=T),]
  gff_df <- subset(gff_df, gff_df[,2] == "Prodigal:2.6")
  gff_df <- separate(gff_df,V9, c("ID", "EC","Name","gene","inference","locus_tag","product","protein_id"), ";")
  gff_df$ID <- gsub("ID=","",gff_df$ID)
  contig_tag_df <- merge(blast_df,gff_df,by="ID")
  write.csv(contig_tag_df,file = paste(name,"contig_tag.csv",sep="_"))
  return(contig_tag_df)
}

contig_tag_df <- contig_match(gff,blast)
write.csv(contig_tag_df,file = "contig_tag.csv")
