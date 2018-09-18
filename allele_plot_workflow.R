# This is an example R code for plotting allele frequencies from breseq outputs
getwd()
df_snp_tob <- read.csv("./tobra/snp_Tobra_Breseq_Output_noutf.csv",header=TRUE) # read in concatenated breseq output

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
cast_tob <- df_melt_cast_fun(df_snp_tob,sample_key_tob) 

############# reduce dataset to mutations >10% #########

cast_tob_orig <- cast_tob
#### Filtering ####
#select only rows with greater than 10% total frequency
cast_tob_cutoff <- (subset(cast_tob, cast_tob$sum >= 10)) #greater than 10% = 70 rows
cast_tob_cutoff <- (subset(cast_tob_cutoff, cast_tob_cutoff$sum < 4600)) #68 rows
nrow(subset(cast_tob_cutoff, cast_tob$count > 1)) #51 rows
write.csv(cast_tob_cutoff,file="snp_tob_noref_cutoff10percent.csv")

t_cast_tob_cutoff <- as.data.frame(t(cast_tob_cutoff[,-c(48,49)]))
tob_b_p_string<- c(rep("P",times=23),rep("B",times=24))
t_cast_tob_cutoff$biofilm_planktonic <-tob_b_p_string
t_cast_tob_cutoff$population <- rownames(t_cast_tob_cutoff)

melt_tob_snp <- (melt(t_cast_tob_cutoff,id.vars = c("population","biofilm_planktonic","day","xmic","replicate"),value.name = "frequency"))
melt_tob_snp$xmic <- as.factor(melt_tob_snp$xmic)
melt_tob_snp$replicate <- as.factor(melt_tob_snp$replicate)

#subset planktonic or biofilm
melt_tob_snp_P <- subset(melt_tob_snp, biofilm_planktonic == "P")
melt_tob_snp_Pr <- subset(melt_tob_snp_P,xmic !=0) 
melt_tob_snp_P0 <- subset(melt_tob_snp_P, xmic ==0)
melt_tob_snp_B <- subset(melt_tob_snp, biofilm_planktonic == "B")
#drop unused factor levles
melt_tob_snp_P[] <- lapply(melt_tob_snp_P, function(x) if(is.factor(x)) factor(x) else x)
melt_tob_snp_Pr[] <- lapply(melt_tob_snp_Pr, function(x) if(is.factor(x)) factor(x) else x)
melt_tob_snp_B[] <- lapply(melt_tob_snp_B, function(x) if(is.factor(x)) factor(x) else x)

################## setting colors for plots ##########

levels(melt_tob_snp$variable)#65

colScale <- c("#4B2012","royalblue", "chartreuse3","darkorange","#0089a0", "darkgreen","blue4", "yellow1","deepskyblue", "mediumorchid3","orangered3","violetred", "#89C5DA", "#DA5724", "#74D944", "#C84248", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030",   "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861","cyan2","darkcyan","orchid1","orange1","mediumpurple1","skyblue4","salmon1","steelblue2","yellowgreen", "lightslategrey","peachpuff","turquoise4","grey26", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7","#ABC6DF","#272617","#D6D676","#8C378E","#56133D","#C9734D","#6E7424","#A0276E","#6D6944","#DC6025","#4B2A30","#AB7187","#4E66B8","#E6B554","#6670EA","#A62E4C","#5EA568","#DE3E36","#7A5A1E","#6890AA","#BEE82C","#8E3631","#9294C6","#6D4C67","#412457","#B2E967","#E0937A","#A9895A","#406B3C","red","#E4C738","#605453","#91A787","#D7BFE7","#D35BEA","#326D5B","#2A412B","#274146","#DDE3DE","#392039","#A69695","#005465")


colscale_tob <- data.frame(levels(as.factor(melt_tob_snp$variable)),colScale[1:65])
colnames(colscale_tob) <- c("variable","color")
names(colScale) <- as.character(colscale_tob$variable)
colScale <- colScale[1:65]

#if you want to manipulate color scale in excel and reimport
#color_key <- read.csv("./plots/color_key.csv",header=T)
#color_key <- color_key[1:73,]
#colnames(color_key) <- c("variable","color")
#colScale <- as.character(color_key$color)
#names(colScale) <- as.character(color_key$variable)
#class(colScale)


#add colors for each variable
melt_tob_snp_Prc <- merge(melt_tob_snp_Pr,colscale_tob,by="variable")
melt_tob_snp_Bc <- merge(melt_tob_snp_B,colscale_tob,by="variable")

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


allele_print_fun(melt_tob_snp_Bc ,3,"Biofilm",colscale_tob,1) 

allele_print_fun(melt_tob_snp_Prc ,3,"Planktonic",colscale_tob,1) # 1 for antibiotics

#END
