require(limma)
require(dplyr)
require(tidyr)

import_data <- function( dirin, targetfile ){
  
  targets <- read.csv( targetfile, sep = "\t")
  
  alldata = data.frame() # empty data frame
  
  for( fname in targets$FileName  ){
    
    filename <- strsplit(fname, ".txt")[[1]]
    
    # print( filename)
    
    RG <- read.maimages(  paste(dirin,fname,sep=""),source="agilent.mean")  
    
    # Normalize within arrays
    MAw <- normalizeWithinArrays(RG, method="loess");
    
    #Normalize between arrays # NOT USED NOW
    #MAw <- normalizeBetweenArrays(MAw, method="quantile");
    
    #Convert normalized MA values back to RG values
    RGnorm <- RG.MA(MAw);  
    
    #Target infos
    tinfo <- targets %>% filter( FileName == fname  )
    
    # Green cy3
    G <- as.data.frame(RGnorm$G)
    colnames(G) <- "intensity"
    gdf <- select( merge( RGnorm$gene, G, by = 0 ) ,
                   ControlType, ProbeName, 
                   GeneName, intensity, Row.names) %>%
            mutate(rownumber = as.numeric(Row.names)) %>%
            select(-Row.names) %>%  arrange(rownumber)
    gdf$dye <- rep("cy3",length(RG$G))
    

    numexp <- as.character( tinfo$Cy3_numexp )
    batch <- as.character( tinfo$original.set )
    line <- as.character( tinfo$Cy3_line )
    genotype <- as.character( tinfo$Cy3_genotype )
    treatmenttype <- as.character( tinfo$Cy3_treatment )
    treatmentduration <- as.character(tinfo$Cy3_treatment.duration )
    agebeginshock <- as.character( tinfo$Cy3_agebeginshock  )
    abstime <- as.character( tinfo$Cy3_time )
    include <- as.character( tinfo$Cy3_excl )

    
    shockduration <- 0
    if( treatmenttype == "HS" ) shockduration <- as.numeric(treatmentduration)
    if( treatmenttype == "Rec" ) shockduration <- (as.numeric(abstime)-as.numeric(agebeginshock) ) - as.numeric(treatmentduration)  

    gdf$codeexp <- rep(numexp,length(RG$G))    
    gdf$batch <- rep(batch,length(RG$G))  
    gdf$line <- rep(line,length(RG$G))
    gdf$genotype <- rep(genotype,length(RG$G))
    gdf$experimentype <- rep(treatmenttype,length(RG$G))  
    gdf$agebeginshock <- rep(as.numeric(agebeginshock),length(RG$G))  
    gdf$shockduration <- rep(as.numeric(shockduration),length(RG$G))  
    gdf$abstime <- rep(as.numeric(abstime),length(RG$G))  
    gdf$include <- rep(include,length(RG$G))  
    #gdf$rownumber <- 1:length(RG$G) ### SUBSTITUTE
    
    #add data to the data frame    
    alldata <- rbind(alldata, gdf)
    
    
    # Red cy5
    R <- as.data.frame(RGnorm$R)
    colnames(R) <- "intensity"
    rdf <- select( merge( RGnorm$gene, R, by = 0 ) ,
                   ControlType, ProbeName, 
                   GeneName, intensity, Row.names) %>%
      mutate(rownumber = as.numeric(Row.names)) %>%
      select(-Row.names) %>%  arrange(rownumber)
    rdf$dye <- rep("cy5",length(RG$R))
    
    numexp <- as.character( tinfo$Cy5_numexp )
    batch <- as.character( tinfo$original.set )
    line <- as.character( tinfo$Cy5_line )
    genotype <- as.character( tinfo$Cy5_genotype )
    treatmenttype <- as.character( tinfo$Cy5_treatment )
    treatmentduration <- as.character(tinfo$Cy5_treatment.duration )
    agebeginshock <- as.character( tinfo$Cy5_agebeginshock  )
    abstime <- as.character( tinfo$Cy5_time )
    include <- as.character( tinfo$Cy5_excl )
    
    shockduration <- 0
    if( treatmenttype == "HS" ) shockduration <- as.numeric(treatmentduration)
    if( treatmenttype == "Rec" ) shockduration <- (as.numeric(abstime)-as.numeric(agebeginshock) ) - as.numeric(treatmentduration)  
    
    
    rdf$codeexp <- rep(numexp,length(RG$R))    
    rdf$batch <- rep(batch,length(RG$R))  
    rdf$line <- rep(line,length(RG$R))
    rdf$genotype <- rep(genotype,length(RG$R))
    rdf$experimentype <- rep(treatmenttype,length(RG$R))  
    rdf$agebeginshock <- rep(as.numeric(agebeginshock),length(RG$R))  
    rdf$shockduration <- rep(as.numeric(shockduration),length(RG$R))  
    rdf$abstime <- rep(as.numeric(abstime),length(RG$R))  
    rdf$include <- rep(include,length(RG$R))  
    #rdf$rownumber <- 1:length(RG$R) ### SUBSTITUTE
    
    #add data to the data frame    
    alldata <- rbind(alldata, rdf)
    
    
  }
  
  return(alldata)
  
}

dirin = "../../../Data/MicroArrays/"

targetfile <- paste(dirin,"Tagets_overall.csv",sep="")

targets <- read.csv( targetfile, sep = "\t")

alldata_expr <- import_data( dirin, targetfile )

save(alldata_expr,file=paste(dirin,"allexperiments.Rdata",sep=""))
