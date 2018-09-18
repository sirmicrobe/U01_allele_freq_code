df_melt_cast_fun <- function(df_snp,sample_key) {
  df_snp$chr_pos <- paste(df_snp$Seq.ID,df_snp$Position,sep="_")
  df_snp <- (merge(df_snp,sample_key,by="Sample"))
  #Mutation data frames found in common reference = false positives
  df_1_10ref_snps <- read.csv("/Users/chrismarshall/Documents/Pitt/Cooper_Lab/Abaum_UO1/experimental/references_breseq_out/final_ref_NZ_CP012004_pAb123/snp_final_ref.csv", header=TRUE)
  df_1_10ref_snps$chr_pos <- paste(df_1_10ref_snps$Seq.ID,df_1_10ref_snps$Position,sep="_")
  #now do same with gene level instead of position level (more aggressive, could miss things)
  df_snp_gene_noref <- df_snp[ !(df_snp$Gene %in% df_1_10ref_snps$Gene), ]
  #see what actually matches - this is what has been removed
  df_snp_refmatch <- df_snp[ (df_snp$chr_pos %in% df_1_10ref_snps$chr_pos), ]
  #create a column for gene_annot_desc
  df_snp_gene_noref$gene_description <- paste(df_snp_gene_noref$Gene,df_snp_gene_noref$Description, sep="::") 
  #combine annotation::gene::description in one column
  df_snp_gene_noref$desc_gene_annot <-  paste(df_snp_gene_noref$Annotation,df_snp_gene_noref$gene_description, sep="::") 
  df_snp_gene_noref$desc_gene_annot <- as.factor(df_snp_gene_noref$desc_gene_annot)
  #remove '%' symbol
  df_snp_gene_noref$Frequency <- gsub( "%", "", as.character(df_snp_gene_noref$Frequency), n)
  df_snp_gene_noref$Frequency <- as.numeric(as.character(df_snp_gene_noref$Frequency))
  #convert Position to numeric
  df_snp_gene_noref$Position <- gsub(":.*","",df_snp_gene_noref$Position)
  df_snp_gene_noref$Position <- as.numeric(as.character(gsub(",","",df_snp_gene_noref$Position)))
  #melted data frame for future casting
  m_df_snp_gene_noref <- melt(df_snp_gene_noref, id=c("Seq.ID","Position","sample_name","replicate","day","biofilm_plank","xmic","Annotation","Gene","Description","Mutation","desc_gene_annot"),measure.vars = c("Frequency"))
  #cast data frame
  c_df_snp_gene_noref <- t(dcast(m_df_snp_gene_noref,sample_name+biofilm_plank+day+xmic+replicate~desc_gene_annot,mean, value.var = "value",fill=0))
  c_df_snp_gene_noref <- as.data.frame(c_df_snp_gene_noref,header=TRUE)
  colnames(c_df_snp_gene_noref) <- as.character(unlist(c_df_snp_gene_noref[1,]))
  c_df_snp_gene_noref  <- c_df_snp_gene_noref[-1,]
  #return(m_df_snp_gene_noref)
  cast_df <- c_df_snp_gene_noref
  #set column order
  setcolorder(cast_df,intersect(as.character(sample_key$sample_name),colnames(cast_df)))
  #convert all values from factor to numeric
  #cast_df_bp <- cast_df[1,]
  cast_df_rows <- rownames(cast_df)
  cast_df <-(sapply(cast_df[2:nrow(cast_df),], function(x) as.numeric(as.character(x)))) 
  #cast_df <- (rbind(cast_df,cast_df_bp))
  #rownames(cast_df) <- c(cast_df_rows[2:nrow(cast_df)],"biofilm_plank")
  rownames(cast_df) <- cast_df_rows[2:length(cast_df_rows)]
  cast_df <- as.data.frame(cast_df)
  cast_df_sum <- transform(cast_df, count=rowSums(cast_df!=0.0), sum=rowSums(cast_df))
  return(cast_df_sum)
}
