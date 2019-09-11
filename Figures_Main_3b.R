library(tidyverse)
require(RColorBrewer)

#load wbid
wbid <- read.csv("../../R_zone/R_Resources/MA_Info_files/agi_id_2014.txt",
                 sep = "\t") %>% select( FeatureNum, WBID ) %>%
  rename( rownumber = FeatureNum )

mytheme <- theme_bw() + theme(
  legend.title  = element_text( size=17),
  #  legend.position = "bottom",
  #	legend.direction = "horizontal",
  legend.key = element_blank(),
  legend.text  = element_text( size=17),
  panel.background = element_blank(),
  panel.grid = element_blank(),
  text = element_text( family="Helvetica", size=19),
  panel.border = element_rect( colour = "black", size=2.5),
  axis.ticks = element_line(size = 2.),
  strip.background = element_rect( fill = "transparent", size = 2.5, colour = "black"  ),
  strip.text = element_text(size = 19)
)

scalecols <- scale_colour_manual(values = c("#009900", "#FF0000", "#0000CC", "#FF8000", "#9933FF", "#00CCCC"))
### very useful: http://www.rapidtables.com/web/color/RGB_Color.htm


hsqtl <- read.csv("../../../Katharina_Trial/TansQTL/eQTL-Table_hs_eQTL.txt", sep = "\t") 
cislist <- (hsqtl %>%  filter( Chr == "IV", Type == "Cis", Eff > 0 ))$WBID
translist <- (hsqtl %>%  filter( Chr == "IV", Type == "Trans", Eff > 0 ))$WBID

load("../../../Data/MicroArrays/projection_axis.Rdata")

p <- ggplot(  ) + theme_bw() +
  geom_freqpoly( data = axis_genes_wbid, aes(
    x = hs_RILs,
    y = ..density..
  ) ) +
  geom_freqpoly( data = axis_genes_wbid %>% 
                   filter( WBID %in% translist ), aes(
                     x = hs_RILs,
                     y = ..density..
                   ) , colour = "red") +
  scale_y_log10()
p

ESfromTag <- function(tag.indicator, norm.tag, norm.no.tag){
  
  no.tag.indicator <- 1 - tag.indicator 
  tag.indicator <- tag.indicator * norm.tag
  no.tag.indicator <- no.tag.indicator * norm.no.tag
  diffvect <- tag.indicator - no.tag.indicator
  
  RES <- cumsum(diffvect)      
  max.ES <- max(RES)
  min.ES <- min(RES)
  ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)
  
  return(ES)  
  
}
GSEA.EnrichmentScore <- function(gene.list, gene.set, cutoff = 10, Nrnd = 1000) {  
  
  N <- length(gene.list) 
  Nh <- length(gene.set) 
  Nm <-  N - Nh 
  
  
  if(Nh < cutoff ) return("NA")
  else{
    
    norm.tag    <- 1./Nh
    norm.no.tag <- 1./Nm
    tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
    
    ES <- ESfromTag(tag.indicator,norm.tag, norm.no.tag)
    
    p <- 0
    mESr <- 0
    vESr <- 0
    for (i in 1:Nrnd) {
      ESr <- ESfromTag(sample(tag.indicator),norm.tag, norm.no.tag)
      mESr <- mESr + ESr
      vESr <- vESr + ESr^2
      p <- p + 1*(ESr > ES)
      if( i %% 1000 == 0) print(i/Nrnd)
    }
    p <- p / Nrnd
    mESr <- mESr / Nrnd
    vESr <- vESr / Nrnd - mESr^2
    zsc <- (ES - mESr) / sqrt(vESr)
    
    return( paste(ES,p, zsc, sep = "_") )
    
  }
  
}

sortgenelist <- (axis_genes_wbid %>%arrange( abs(hs_RILs) ))$WBID
enrichment_chrom <- data.frame()
for( c in as.vector(hsqtl$Chr  %>% unique()  %>% sort) ){
  
  translist <- (hsqtl %>%  filter( Chr == c, Type == "Trans", Eff > 0 ))$WBID
  es <- GSEA.EnrichmentScore(sortgenelist, translist , Nrnd = 10000)
  ess <- as.numeric(strsplit(es, "_")[[1]][1])
  pval <- as.numeric(strsplit(es, "_")[[1]][2])
  zsc <- as.numeric(strsplit(es, "_")[[1]][3])
  d <- as.data.frame(list(es = ess,
                          pval = pval,
                          zsc = zsc,
                          eff = "N2", chr = c))
  enrichment_chrom <- rbind(d,enrichment_chrom)
  
  translist <- (hsqtl %>%  filter( Chr == c, Type == "Trans", Eff < 0 ))$WBID
  es <- GSEA.EnrichmentScore(sortgenelist, translist , Nrnd = 10000)
  ess <- as.numeric(strsplit(es, "_")[[1]][1])
  pval <- as.numeric(strsplit(es, "_")[[1]][2])
  zsc <- as.numeric(strsplit(es, "_")[[1]][3])
  d <- as.data.frame(list(es = ess,
                          pval = pval,
                          zsc = zsc,
                          eff = "CB", chr = c))
  enrichment_chrom <- rbind(d,enrichment_chrom)
  
}

enrichment_chrom <- enrichment_chrom %>%  mutate( 
  zsc = as.numeric(zsc) , pval = as.numeric(pval), es = as.numeric(es) )

sign_thre <- 0.05 / 6
chr_levels <- enrichment_chrom$chr[order(desc(enrichment_chrom$chr))]
enrichment_chrom$chr2 <- factor(enrichment_chrom$chr, levels = chr_levels)
p <- ggplot(enrichment_chrom %>% arrange(chr) ) + mytheme + 
  scale_fill_manual(values = c("#FF8000", "#990099",  "#4C9900") ) +
  aes(x = chr2, y = zsc,  alpha = 0.75 + 0.25*(pval < sign_thre  |  pval > (1-sign_thre)  ) ) +
  geom_col( fill = "#FF8000" ) + theme( legend.position = "none" ) +
  scale_x_discrete( "Chromosome" ) + scale_y_continuous("Z score"  ) +
  facet_wrap( ~ eff) 
p

#ggsave(p, filename = "../../../Figures/Main/Fig3/enrichment_chr.pdf", width = 8.24, height = 3.91)
ggsave(p, filename = "../../../Figures/Main/Fig3/enrichmentabs_chr.pdf", width = 8.24, height = 3.91)

