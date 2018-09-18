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
