# This is R code for processing data from the Abaum ciprofloxacin resistance evolution experiments
library("vegan")
library("plyr")
library("RColorBrewer")
library("ggplot2")
library("data.table")
library("dplyr")
library("reshape2")

theme_set(theme_bw())

getwd()
df_snp_cip <- read.csv("./cipro/snp_cipro_Breseq_Output_noutf.csv",header=TRUE) # read in concatenated breseq output

#get the data frame in shape for melting
#this function removes reference errors (provide path to reference below)
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

#cast data using function
cast_cip <- df_melt_cast_fun(df_snp_cip,sample_key_cip) 

############# reduce dataset to mutations >10% #########

cast_cip_orig <- cast_cip
#### Filtering ####
#select only rows with greater than 10% total frequency
cast_cip_cutoff <- (subset(cast_cip, cast_cip$sum >= 10)) #greater than 10% = 70 rows
cast_cip_cutoff <- (subset(cast_cip_cutoff, cast_cip_cutoff$sum < 4600)) #68 rows
nrow(subset(cast_cip_cutoff, cast_cip$count > 1)) #51 rows
write.csv(cast_cip_cutoff,file="snp_cip_noref_cutoff10percent.csv")

t_cast_cip_cutoff <- as.data.frame(t(cast_cip_cutoff[,-c(48,49)]))
cip_b_p_string<- c(rep("P",times=23),rep("B",times=24))
t_cast_cip_cutoff$biofilm_planktonic <-cip_b_p_string
t_cast_cip_cutoff$population <- rownames(t_cast_cip_cutoff)

melt_cip_snp <- (melt(t_cast_cip_cutoff,id.vars = c("population","biofilm_planktonic","day","xmic","replicate"),value.name = "frequency"))
melt_cip_snp$xmic <- as.factor(melt_cip_snp$xmic)
melt_cip_snp$replicate <- as.factor(melt_cip_snp$replicate)

#subset planktonic or biofilm
melt_cip_snp_P <- subset(melt_cip_snp, biofilm_planktonic == "P")
melt_cip_snp_Pr <- subset(melt_cip_snp_P,xmic !=0) 
melt_cip_snp_P0 <- subset(melt_cip_snp_P, xmic ==0)
melt_cip_snp_B <- subset(melt_cip_snp, biofilm_planktonic == "B")
#drop unused factor levles
melt_cip_snp_P[] <- lapply(melt_cip_snp_P, function(x) if(is.factor(x)) factor(x) else x)
melt_cip_snp_Pr[] <- lapply(melt_cip_snp_Pr, function(x) if(is.factor(x)) factor(x) else x)
melt_cip_snp_B[] <- lapply(melt_cip_snp_B, function(x) if(is.factor(x)) factor(x) else x)

################## setting colors for plots ##########

levels(melt_cip_snp$variable)#65

colScale <- c("#4B2012","royalblue", "chartreuse3","darkorange","#0089a0", "darkgreen","blue4", "yellow1","deepskyblue", "mediumorchid3","orangered3","violetred", "#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030",   "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861","cyan2","darkcyan","orchid1","orange1","mediumpurple1","skyblue4","salmon1","steelblue2","yellowgreen", "lightslategrey","peachpuff","turquoise4","grey26", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7","#ABC6DF","#272617","#D6D676","#8C378E","#56133D","#C9734D","#6E7424","#A0276E","#6D6944","#DC6025","#4B2A30","#AB7187","#4E66B8","#E6B554","#6670EA","#A62E4C","#5EA568","#DE3E36","#7A5A1E","#6890AA","#BEE82C","#8E3631","#9294C6","#6D4C67","#412457","#B2E967","#E0937A","#A9895A","#406B3C","red","#E4C738","#605453","#91A787","#D7BFE7","#D35BEA","#326D5B","#2A412B","#274146","#DDE3DE","#392039","#A69695","#005465")


colscale_cip <- data.frame(levels(as.factor(melt_cip_snp$variable)),colScale[1:65])
colnames(colscale_cip) <- c("variable","color")
names(colScale) <- as.character(colscale_cip$variable)
colScale <- colScale[1:65]

#if you want to manipulate color scale in excel and reimport
#color_key <- read.csv("./plots/color_key.csv",header=T)
#color_key <- color_key[1:73,]
#colnames(color_key) <- c("variable","color")
#colScale <- as.character(color_key$color)
#names(colScale) <- as.character(color_key$variable)
#class(colScale)


#add colors for each variable
melt_cip_snp_Prc <- merge(melt_cip_snp_Pr,colscale_cip,by="variable")
melt_cip_snp_Bc <- merge(melt_cip_snp_B,colscale_cip,by="variable")

########### allele plotting ###########
allele_print_fun <- function(melt_snp,repl,b_p,color_key,abx) {
  melt_snp_rep <- subset(melt_snp,replicate==repl)
  melt_snp_rep_legend <- subset(melt_snp_rep,frequency>0)
  melt_snp_rep_legend <- as.character(melt_snp_rep_legend$variable)
  melt_snp_rep_droplevel <- subset(melt_snp_rep,variable %in% melt_snp_rep_legend)
  melt_snp_rep_droplevel[] <- lapply(melt_snp_rep_droplevel, function(x) if(is.factor(x)) factor(x) else x)
  if(abx==1) {
    #subset sub-mic level to plot
    melt_snp_rep_half <- subset(melt_snp_rep,xmic==0.5)
    melt_snp_rep_half[] <- lapply(melt_snp_rep_half, function(x) if(is.factor(x)) factor(x) else x)
    #add days 1 and 3 data points to ramp plot
    melt_snp_rep_ramp <- subset(melt_snp_rep_half, day %in% c(1,3))
    melt_snp_rep_ramp[] <- lapply(melt_snp_rep_ramp, function(x) if(is.factor(x)) factor(x) else x)
    levels(melt_snp_rep_ramp$xmic) <- "ramp"
    melt_snp_rep <- rbind(melt_snp_rep,melt_snp_rep_ramp)
    #subset just ramp to plot
    melt_snp_rep_ramp <- subset(melt_snp_rep,xmic %in% c( "1", "2", "4", "ramp"))
    melt_snp_rep_ramp[] <- lapply(melt_snp_rep_ramp, function(x) if(is.factor(x)) factor(x) else x)
    #remove unused factor levels from 0.5 and ramp
    melt_snp_rep_ramp_legend <- subset(melt_snp_rep_ramp,frequency>0)
    melt_snp_rep_ramp_legend <- as.character(melt_snp_rep_ramp_legend$variable)
    melt_snp_rep_ramp_droplevel <- subset(melt_snp_rep_ramp,variable %in% melt_snp_rep_ramp_legend)
    melt_snp_rep_ramp_droplevel[] <- lapply(melt_snp_rep_ramp_droplevel, function(x) if(is.factor(x)) factor(x) else x)
    levels(melt_snp_rep_ramp_droplevel$xmic) <- c("ramp","ramp","ramp","ramp")
    melt_snp_rep_half_legend <- subset(melt_snp_rep_half,frequency>0)
    melt_snp_rep_half_legend <- as.character(melt_snp_rep_half_legend$variable)
    melt_snp_rep_half_droplevel <- subset(melt_snp_rep_half,variable %in% melt_snp_rep_half_legend)
    melt_snp_rep_half_droplevel[] <- lapply(melt_snp_rep_half_droplevel, function(x) if(is.factor(x)) factor(x) else x)
    #add day 0 time point where everything starts at 0 frequency
    df_half_from0 <- data.frame(melt_snp_rep_half_droplevel$variable,rep(levels(as.factor(melt_snp_rep_half_droplevel$population)),length(levels(melt_snp_rep_half_droplevel$variable))),rep("P",length(levels(melt_snp_rep_half_droplevel$variable))),rep(0,length(levels(melt_snp_rep_half_droplevel$variable))),rep(0.5,length(levels(melt_snp_rep_half_droplevel$variable))),rep(5,length(levels(melt_snp_rep_half_droplevel$variable))),rep(0,length(levels(melt_snp_rep_half_droplevel$variable))))
    colnames(df_half_from0) <- colnames(melt_snp_rep_half_droplevel[1:7])
    df_half_from0 <- merge(df_half_from0,color_key,by="variable") 
    melt_snp_rep_half_droplevel <- rbind(df_half_from0,melt_snp_rep_half_droplevel)
    
    df_ramp_from0 <- data.frame(melt_snp_rep_ramp_droplevel$variable,rep(levels(as.factor(melt_snp_rep_ramp_droplevel$population)),length(levels(melt_snp_rep_ramp_droplevel$variable))),rep("P",length(levels(melt_snp_rep_ramp_droplevel$variable))),rep(0,length(levels(melt_snp_rep_ramp_droplevel$variable))),rep("ramp",length(levels(melt_snp_rep_ramp_droplevel$variable))),rep(5,length(levels(melt_snp_rep_ramp_droplevel$variable))),rep(0,length(levels(melt_snp_rep_ramp_droplevel$variable))))
    colnames(df_ramp_from0) <- colnames(melt_snp_rep_ramp_droplevel[1:length(colnames(melt_snp_rep_ramp_droplevel))-1])
    df_ramp_from0 <- merge(df_ramp_from0,color_key,by="variable") 
    melt_snp_rep_ramp_droplevel <- rbind(df_ramp_from0,melt_snp_rep_ramp_droplevel)
    
    
    #plot 0.5 mic
    p_p_0.5_l <- ggplot(melt_snp_rep_half_droplevel,aes(x=day,y=frequency,color=variable,shape=xmic))+
      theme(axis.title.x = element_text(face="bold",size=22),axis.text.x = element_text(angle=0,face="bold", colour = "black", vjust=1, hjust = 1, size=22), axis.text.y = element_text(colour = "black", face="bold",size=22),axis.title.y = element_text(face="bold", size=22), plot.title = element_text(face="bold",size = 20), legend.title = element_text(size=14), legend.text = element_text(size = 13), legend.position="bottom",strip.text.x = element_text(size=12, face="bold"), strip.text.y = element_text(size=12, face="bold"),strip.background = element_rect(colour="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      geom_point(size=4)+ 
      ylab("Allele Frequency") +
      scale_x_continuous(breaks = c(0,2,4,6,8,10,12), limits = c(0,12))+
      guides(colour=guide_legend(nrow=10,title.position="top",title="Mutation:Gene:Description")) + 
      ggtitle(paste(b_p,"Replicate",repl,sep=" ")) +  
      geom_line(aes(linetype=xmic),size=2) + 
      scale_colour_manual(values=colScale,drop=T) + 
      facet_grid(~xmic)
    pdf(paste(b_p,"0.5_Replicate",repl,"l.pdf",sep="_"), width= 24, height = 10, useDingbats=F) 
    print(p_p_0.5_l) #print our plot
    dev.off() #stop making pdfs
    
    #plot ramp
    p_p_ramp_l <- ggplot(melt_snp_rep_ramp_droplevel,aes(x=day,y=frequency,color=variable,shape=xmic))+
      theme(axis.title.x = element_text(face="bold",size=22),axis.text.x = element_text(angle=0,face="bold", colour = "black", vjust=1, hjust = 1, size=22), axis.text.y = element_text(colour = "black", face="bold",size=22),axis.title.y = element_text(face="bold", size=22), plot.title = element_text(face="bold",size = 20), legend.title = element_text(size=14), legend.text = element_text(size = 13), legend.position="bottom",strip.text.x = element_text(size=12, face="bold"), strip.text.y = element_text(size=12, face="bold"),strip.background = element_rect(colour="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      geom_point(size=4)+ 
      ylab("Allele Frequency") +
      scale_x_continuous(breaks = c(0,2,4,6,8,10,12), limits = c(0,12))+
      guides(colour=guide_legend(nrow=10,title.position="top",title="Mutation:Gene:Description")) + 
      ggtitle(paste(b_p,"Replicate",repl,sep=" ")) +    
      geom_line(aes(linetype=xmic),size=2) + 
      scale_colour_manual(values=colScale,drop=T) + 
      facet_grid(~xmic)
    pdf(paste(b_p,"Ramp_Replicate",repl,"l.pdf",sep="_"), width= 24, height = 10, useDingbats=F) 
    print(p_p_ramp_l) #print our plot
    dev.off() #stop making pdfs
    
    #No legends on plot
    p_p_0.5 <- ggplot(melt_snp_rep_half_droplevel,aes(x=day,y=frequency,color=variable,shape=xmic))+
      theme(axis.title.x = element_text(face="bold",size=22),axis.text.x = element_text(angle=0,face="bold", colour = "black", vjust=1, hjust = 1, size=22), axis.text.y = element_text(colour = "black", face="bold",size=22),axis.title.y = element_text(face="bold", size=22), plot.title = element_text(face="bold",size = 20), legend.title = element_text(size=14), legend.text = element_text(size = 13), legend.position="none",strip.text.x = element_text(size=12, face="bold"), strip.text.y = element_text(size=12, face="bold"),strip.background = element_rect(colour="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      geom_point(size=4)+ 
      ylab("Allele Frequency") +
      scale_x_continuous(breaks = c(0,2,4,6,8,10,12), limits = c(0,12))+
      guides(colour=guide_legend(nrow=10,title.position="top",title="Mutation:Gene:Description")) + 
      ggtitle(paste(b_p,"Replicate",repl,sep=" ")) +  
      geom_line(aes(linetype=xmic),size=2) + 
      scale_colour_manual(values=colScale,drop=T) + 
      facet_grid(~xmic)
    pdf(paste(b_p,"0.5_Replicate",repl,".pdf",sep="_"), width= 6, height = 6, useDingbats=F) 
    print(p_p_0.5) #print our plot
    dev.off() #stop making pdfs
    
    #plot ramp
    p_p_ramp <- ggplot(melt_snp_rep_ramp_droplevel,aes(x=day,y=frequency,color=variable,shape=xmic))+
      theme(axis.title.x = element_text(face="bold",size=22),axis.text.x = element_text(angle=0,face="bold", colour = "black", vjust=1, hjust = 1, size=22), axis.text.y = element_text(colour = "black", face="bold",size=22),axis.title.y = element_text(face="bold", size=22), plot.title = element_text(face="bold",size = 20), legend.title = element_text(size=14), legend.text = element_text(size = 13), legend.position="none",strip.text.x = element_text(size=12, face="bold"), strip.text.y = element_text(size=12, face="bold"),strip.background = element_rect(colour="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      geom_point(size=4)+ 
      ylab("Allele Frequency") +
      scale_x_continuous(breaks = c(0,2,4,6,8,10,12), limits = c(0,12))+
      guides(colour=guide_legend(nrow=10,title.position="top",title="Mutation:Gene:Description")) + 
      ggtitle(paste(b_p,"Replicate",repl,sep=" ")) +  
      geom_line(aes(linetype=xmic),size=2) + 
      scale_colour_manual(values=colScale,drop=T) + 
      facet_grid(~xmic)
    pdf(paste(b_p,"Ramp_Replicate",repl,".pdf",sep="_"), width= 6, height = 6, useDingbats=F) 
    print(p_p_ramp) #print our plot
    dev.off() #stop making pdfs
    } else {
      df_0_from0 <- data.frame(melt_snp_rep_droplevel$variable,rep(levels(as.factor(melt_snp_rep_droplevel$population)),length(levels(melt_snp_rep_droplevel$variable))),rep("P",length(levels(melt_snp_rep_droplevel$variable))),rep(0,length(levels(melt_snp_rep_droplevel$variable))),rep(0,length(levels(melt_snp_rep_droplevel$variable))),rep(repl,length(levels(melt_snp_rep_droplevel$variable))),rep(0,length(levels(melt_snp_rep_droplevel$variable))))
      colnames(df_0_from0) <- c("variable", "population","biofilm_planktonic","day","xmic","replicate","frequency")
      df_0_from0 <- merge(df_0_from0,color_key,by="variable") 
      melt_snp_rep_droplevel <- rbind(df_0_from0,melt_snp_rep_droplevel)
    
    #plot 0 mic
      p_p_0l <- ggplot(melt_snp_rep_droplevel,aes(x=day,y=frequency,color=variable))+
        theme(axis.title.x = element_text(face="bold",size=22),axis.text.x = element_text(angle=0,face="bold", colour = "black", vjust=1, hjust = 1, size=22), axis.text.y = element_text(colour = "black", face="bold",size=22),axis.title.y = element_text(face="bold", size=22), plot.title = element_text(face="bold",size = 20), legend.title = element_text(size=14), legend.text = element_text(size = 13), legend.position="bottom",strip.text.x = element_text(size=12, face="bold"), strip.text.y = element_text(size=12, face="bold"),strip.background = element_rect(colour="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        geom_point(size=4)+ 
        ylab("Allele Frequency") +
        scale_x_continuous(breaks = c(0,2,4,6,8,10,12), limits = c(0,12))+
        guides(colour=guide_legend(nrow=10,title.position="top",title="Mutation:Gene:Description")) + 
        ggtitle(paste(b_p,"Replicate",repl,sep=" ")) +  
        geom_line(aes(linetype=xmic),size=2) + 
        scale_colour_manual(values=colScale,drop=T)
      pdf(paste(b_p,"No_Abx_Replicate",repl,"l.pdf",sep="_"), width= 24, height = 10, useDingbats=F) 
      print(p_p_0l) #print our plot
      dev.off() #stop making pdfs
    
    #No legends on plot
      p_p_0 <- ggplot(melt_snp_rep_droplevel,aes(x=day,y=frequency,color=variable))+
        theme(axis.title.x = element_text(face="bold",size=22),axis.text.x = element_text(angle=0,face="bold", colour = "black", vjust=1, hjust = 1, size=22), axis.text.y = element_text(colour = "black", face="bold",size=22),axis.title.y = element_text(face="bold", size=22), plot.title = element_text(face="bold",size = 20), legend.title = element_text(size=14), legend.text = element_text(size = 13), legend.position="none",strip.text.x = element_text(size=12, face="bold"), strip.text.y = element_text(size=12, face="bold"),strip.background = element_rect(colour="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        geom_point(size=4)+ 
        ylab("Allele Frequency") +
        scale_x_continuous(breaks = c(0,2,4,6,8,10,12), limits = c(0,12))+
        guides(colour=guide_legend(nrow=10,title.position="top",title="Mutation:Gene:Description")) + 
        ggtitle(paste(b_p,"Replicate",repl,sep=" ")) +  
        geom_line(aes(linetype=xmic),size=2) + 
        scale_colour_manual(values=colScale,drop=T)
      pdf(paste(b_p,"No_Abx_Replicate",repl,".pdf",sep="_"), width= 24, height = 10, useDingbats=F)
      print(p_p_0) #print our plot
      dev.off() #stop making pdfs
  }
} 


allele_print_fun(melt_cip_snp_Bc ,3,"Biofilm",colscale_cip,1) 

allele_print_fun(melt_cip_snp_Prc ,3,"Planktonic",colscale_cip,1) # 1 for antibiotics




######## GD tools output ###########
#using GD tools output for mutation type counting

setwd("~/Documents/Chris_Argonne_Work_Mac/Pitt/Cooper_Lab/Abaum_UO1/experimental/alfonso_cipro_evol/populations_new_ref/anvio/gd_files")

temp = list.files(pattern="*.txt") # no mutations in sample 58 (look at NJE)
all_snp <- lapply(temp,read.delim) #create a list with each object in list a data frame of mutation table
# headers we want to keep
headers <- c("aa_new_seq","aa_position","aa_ref_seq","codon_new_seq","codon_number","codon_position","codon_ref_seq","gene_name","gene_position","gene_product","frequency","mutation_category","new_seq","position","seq_id","snp_type","type","title")
#take a look at the first data frame in the list
df_snp_01 <- all_snp[[1]]
df_snp_01 <- df_snp_01[,headers] #select only headers we want to keep
head(df_snp_01)
# convert all data frames in the list to one large data frame
all_snp_df <- as.data.frame(data.table::rbindlist(all_snp,fill=T),colClasses = c("character"))
all_snp_df <- all_snp_df[,headers] #keep headers we want
all_snp_df$mutation <- paste(all_snp_df$aa_ref_seq,all_snp_df$aa_position,all_snp_df$aa_new_seq, sep=":") #create mutation column
View(all_snp_df)
nrow(all_snp_df) #753
all_snp_df$title <- gsub("sample", "breseq",all_snp_df$title)
all_snp_df$sample <- all_snp_df$title


sample_key$sample <- sample_key$Sample
sample_key$sample <- gsub("_v_vc17978","",sample_key$sample)
sample_key$Sample <- NULL


all_snp_df <- (merge(all_snp_df,sample_key,by="sample"))

write.csv(all_snp_df,file="df_snp_noref.csv")


long_snp_df <- merge(all_snp_df,df_snp_gene_noref,by="sample_seq_pos")

########## tansitions vs transversions ratio function #######
mtype <- function(subset_df,repday){
  transitions = c('AG', 'CT', 'GA', 'TC')
  transversion_vs_transition <- data.frame('mtype' = character(), 'freq' = integer(),'count'=integer(), stringsAsFactors=FALSE) #create data frame with 3 columns
  transversion_vs_transition[1, ] <- c('transversion', 0,0) #assign values to df
  transversion_vs_transition[2, ] <- c('transition', 0,0)
  transversion_vs_transition[transversion_vs_transition$mtype == "transition", ]$freq = sum(subset_df[subset_df$Mutation %in% transitions, ]$frequency) # sum the mutation frequencies of all transitions
  transversion_vs_transition[transversion_vs_transition$mtype == "transversion", ]$freq = sum(subset_df[!(subset_df$Mutation %in% transitions), ]$frequency) # sum frequencies of transversions
  transversion_vs_transition[transversion_vs_transition$mtype == "transition", ]$count <- as.numeric(nrow(subset_df[subset_df$Mutation %in% transitions, ])) # count number of transitions
  transversion_vs_transition[transversion_vs_transition$mtype == "transversion", ]$count <-(nrow(subset_df[!(subset_df$Mutation %in% transitions), ]))#count number of transversions
  transversion_vs_transition$freq <- as.numeric(transversion_vs_transition$freq) 
  transversion_vs_transition$x <- 'all'
  transversion_vs_transition_ratio <- transversion_vs_transition$freq[2] / transversion_vs_transition$freq[1]
  
  transversion_vs_transition_colors <- c('transversion' = "#451910",
                                         'transition' = '#8C8983')

  transversion_vs_transition$count <- as.numeric(transversion_vs_transition$count)

  
  #plots
  p <- ggplot(transversion_vs_transition, aes(x=factor(x), y=freq)) + geom_bar(aes(fill=mtype), stat="identity")
  p <- p + scale_color_manual(values = transversion_vs_transition_colors, guide = guide_legend(override.aes=aes(fill=NA)))
  p <- p + scale_fill_manual(values = transversion_vs_transition_colors)
  p <- p + coord_flip() +ggtitle(repday)
  p <- p + xlab(paste('k:', round(transversion_vs_transition_ratio, digits=2), sep=" "))
  p <- p + theme_bw() + theme(axis.text.y=element_blank(),
                              axis.ticks.y=element_blank(),
                              legend.position="bottom",
                              axis.line.y = element_line(colour = "black"),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.border = element_blank(),
                              panel.background = element_blank(),
                              axis.title.x=element_blank(),
                              axis.ticks.x=element_blank(),
                              axis.text.x = element_blank(),
                              plot.margin=unit(c(0.20,0,0.20,0), "cm"))
  p2 <- ggplot(transversion_vs_transition, aes(x=x, y=count,fill=mtype)) + geom_bar( stat="identity",position="stack")
  p2 <- p2 + scale_color_manual(values = transversion_vs_transition_colors, guide = guide_legend(override.aes=aes(fill=NA)))
  p2 <- p2 + scale_fill_manual(values = transversion_vs_transition_colors)
  p2 <- p2 + coord_flip() + ggtitle(paste(repday," n=",length(unique(long_snp_df_SNP_P_4x_d12$replicate.y)))) + ylab(paste("Number of Mutations, total=",sum(transversion_vs_transition$count)))
  p2 <- p2 + theme_bw() + theme(axis.text.y=element_blank(),
                                legend.position="bottom",
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                plot.margin=unit(c(0.20,0,0.20,0), "cm")) 
  #p2 #+ scale_y_continuous(breaks=c(0,2,4,6,8,10))
  return(list(p,p2,transversion_vs_transition_ratio))

  #bonus - synonymous vs nonsynonymous
  
} #if you save the output into an object you can access the different data in the list: mtype_list[[1]] = ratio plot, [[2]] = count plot, [[3]] = transition / transversion ratio

mtype(long_snp_df_P_05x,"all") # 1.5716, 23 mutations
mtype(long_snp_df_B_05x,"all") # 0.247, 58 mutations

##### other genome dynamics - dN/dS, indel rate ####
uniq_SNP_P <- long_snp_df_SNP_P[!duplicated(long_snp_df_SNP_P[,"desc_gene_annot"]),]
uniq_SNP_B <- long_snp_df_SNP_B[!duplicated(long_snp_df_SNP_B[,"desc_gene_annot"]),]

NvS_fun <- function(df,pop) {
  NvS <- (length(df$snp_type) - length(which(df$snp_type == "synonymous"))) / (length(which(df$snp_type == "synonymous")))
  print(paste0(pop,": Nonsynonymous over synonymous = ", (length(df$snp_type) - length(which(df$snp_type == "synonymous"))), "/",(length(which(df$snp_type == "synonymous")))))
  print(NvS)
}
  
NvS_fun(uniq_SNP_B,"all_b")

#count intergenic sequences
(length(which(uniq_Pr$snp_type=="intergenic"))+3) / (length(uniq_Pr$snp_type)-3) #.32 (2/5 = 0.2 in P0 no abx)
(length(which(uniq_Br$snp_type=="intergenic"))+1 )/ (length(uniq_Br$snp_type)-1) #0.297 (2/12=0.167 in B0 no abx)



############## Diversity #################

snp_gene_noref_cutoff <- (apply(t_cast_cip_cutoff,2,as.character))
snp_gene_noref_cutoff <- (apply(snp_gene_noref_cutoff,2,as.numeric))
rownames(snp_gene_noref_cutoff) <- rownames(tt_cast_cip_cutoff)
snp_gene_noref_cutoff <- as.data.frame(snp_gene_noref_cutoff)


colnames(snp_gene_noref_cutoff) 
rownames(snp_gene_noref_cutoff)
b_p_string<- c(rep("P",times=42),rep("B",times=42))
snp_gene_noref_cutoff$biofilm_planktonic <-b_p_string
snp_gene_noref_cutoff$population <- rownames(snp_gene_noref_cutoff)

gene_divers_shan <- diversity(snp_gene_noref_cutoff[,4:76],index="shannon")
gene_divers_invsimps <- diversity(snp_gene_noref_cutoff[,4:76],index="invsimpson")

gene_divers_df <- (cbind(snp_gene_noref_cutoff,gene_divers_shan,gene_divers_invsimps,rowSums(snp_gene_noref_cutoff[,4:76]>0)))
colnames(gene_divers_df)[81] <- "richness"
gene_divers_df$bp_day <- paste(gene_divers_df$biofilm_planktonic,gene_divers_df$day,sep="_")

kruskal.test(gene_divers_df$gene_divers_shan ~ as.factor(gene_divers_df$biofilm_planktonic),data=gene_divers_df)
#overall significance (of no abx,ramp,+subinhibitory) is that biofilm is significantly more diverse than planktonic - p=0.005, chi-squared=7.723

#shannon and richness
m_gene_divers_df <- melt(gene_divers_df,id.vars=c("day","xmic","biofilm_planktonic","replicate","population"),measure.vars=c("gene_divers_shan",'gene_divers_invsimps',"richness"))
m_gene_divers_df <- subset(m_gene_divers_df, xmic !=0)
mr_gene_divers_df <- subset(m_gene_divers_df, xmic !=0.5)
m0.5_gene_divers_df <- subset(m_gene_divers_df, xmic ==0.5)
m0.5d1_gene_divers_df <- subset(m0.5_gene_divers_df,day %in% c(1,3))
mr_gene_divers_df <- rbind(mr_gene_divers_df,m0.5d1_gene_divers_df)

#plot shannon diversity of evolutionary rescue
p_gene_shan <- ggplot(subset(mr_gene_divers_df,variable=="gene_divers_shan") , aes(as.factor(day), value, color=biofilm_planktonic))
p_gene_shan  + theme(text = element_text(size=16, face="bold"), axis.text.x = element_text(vjust = 0.5),legend.title=element_blank()) +  xlab("Time (days)") + ylab("Shannon Index") + geom_boxplot()+ geom_point(aes(fill = biofilm_planktonic), size = 3, shape = 21, position = position_jitterdodge())+ scale_color_manual(values=c("royalblue","red")) +scale_fill_manual(values=c("royalblue","red"))

#plot shannon of subinhibitory
p_gene0.5_shan <- ggplot(subset(m0.5_gene_divers_df,variable=="gene_divers_shan") , aes(as.factor(day), value, color=biofilm_planktonic))
p_gene0.5_shan   + theme(text = element_text(size=16, face="bold"), axis.text.x = element_text(vjust = 0.5),legend.title=element_blank()) +  xlab("Time (days)") + ylab("Shannon Index") + geom_boxplot()+ geom_point(aes(fill = biofilm_planktonic), size = 3, shape = 21, position = position_jitterdodge())+ scale_color_manual(values=c("royalblue","red")) +scale_fill_manual(values=c("royalblue","red"))


#Dunn (1964) Kruskal-Wallis multiple comparison p-values adjusted with the Benjamini-Hochberg method.

library("FSA")

dunnTest(gene_divers_df_0.5$gene_divers_shan ~ as.factor(gene_divers_df_0.5$bp_day),
         data=gene_divers_df_0.5,
         method="bh") #multiple comparison test 
#            z             p
#B_1 - P_1  1.25996814 0.207680856
#B_3 - P_3  2.05328141 0.040045293
#B_4 - P_4  2.79992920 0.005111382
#B_9 - P_9  0.93330973 0.350660071
#B_12 - P_12  0.46665487 0.640746826

p0.5 <- data.frame("day"=c(1,3,4,9,12),"z"=c(1.25996814,2.05328141,2.79992920,0.93330973,0.46665487),"Raw.p"=c(0.207680856,0.040045293,0.005111382,0.350660071,0.640746826))
p0.5 = p0.5[order(p0.5$Raw.p),]
p0.5$Bonferroni = 
  p.adjust(p0.5$Raw.p, 
           method = "bonferroni")
p0.5$BH = 
  p.adjust(p0.5$Raw.p, 
           method = "BH")

dunnTest(gene_divers_df_ramp$gene_divers_shan ~ as.factor(gene_divers_df_ramp$bp_day),
         data=gene_divers_df_ramp,
         method="bh") #multiple comparison test 
#                   Z       punadj
#4  B_4 - P_4  0.19404598 0.846139878
#6 B_6 - P_6  0.00000000 1.000000000
#7 B_7 - P_7  0.15523679 0.876634650
#9 B_9 - P_9  0.46571036 0.641422839
#10  B_10 - P_10 -0.50451955 0.613896327
#12 B_12 - P_12 -0.07761839 0.938131611

                                            
############ Plotting MIC over time ################   
mic <- read.csv("stats_MIC.csv",header=F)
colnames(mic) <- mic[1,]
colnames(mic) <- c("pop",0,1,3,4,6,7,9,10,12)
mic <- mic[-1,]

m_mic<- (melt(mic,id.vars = c("pop")))

sum_mic <- summarySE(m_mic,measurevar=c("value"),groupvars=c("day","variable"))
colnames(sum_mic) <- c("pop","day","N", "mic","sd","se","ci")
sum_mic$b_p <- c(rep("b",27),rep("p",27)) 
sum_mic[45,4] <- 32
sum2_mic <- summarySE(sum_mic,measurevar=c("mic"),groupvars=c("b_p","day"))
sum_mic$day <- as.numeric(as.character(sum_mic$day))

m2_mic <- sum2_mic[,c(1,2,4,7)]
m2_mic$day <- as.numeric(as.character(m2_mic$day))
m2_mic$cilend <- m2_mic$mic-0.5*m2_mic$ci
m2_mic$ciuend <- m2_mic$mic+0.5*m2_mic$ci

p  + geom_vline(xintercept=c(3.1,6.1,9.1),linetype="dashed",alpha=0.4) + # add vertical line
  geom_hline(yintercept=c(.5),linetype="dashed") +  #add horizontal line
  scale_color_manual("",values=c("blue","red")) + # color points red and blue
  geom_ribbon(data=m2_mic, aes(ymin=mic,ymax=ciuend,linetype=NA),alpha=0.2) + #add shaded CI using m2_mic
  scale_fill_manual("",values="grey70") + # fill the CI with grey color
  scale_x_continuous(breaks=c(0,2,4,6,8,10,12)) + #add the tick marks to x-axis
  xlab("Time (days)") + #label x-axis
  ylab(expression(Minimum~Inhibitory~Concentration~(mg~L^-1))) + #label y-axis (use expression() for exponent)
  theme(text = element_text(size=18, face="bold")) +#bold text and increase font
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #remove gridlines
  scale_y_continuous(trans="log2") #convert y axis to log2 scale

                                            
  #END
