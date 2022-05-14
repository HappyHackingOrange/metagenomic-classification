
### DESCRIPTION ################################################################

# Plots a heat map with Hamming dissimilarity percent values between each 
# assembly sequence and each reference sequence.

# This applies to the assembly sequences assembled by PhyloFlash for the mock 
# community dataset and the 26 reference 16s sequences provided to me.



### HEADER CODE ################################################################

# Clean out environment
rm(list=ls())
gc()



### LIBRARIES ##################################################################

# Enhanced data.frame
library(data.table)

# Elegant data visualizations using grammar of graphics
library(ggplot2)



### CODE ######################################################################

# Binds together all CSV files containing Hamming distance info that was
# generated from other script, for each assembly sequence and each reference
# sequence.
# DT <- rbindlist(lapply(dir("results/amanda/hd",pattern=".csv",full.names = TRUE),fread))
DT <- fread("results/amanda/hd/hd_2022-04-07.csv")

# Calculate Hamming dissimilarity percentage.
# DT[,.(`HD for closest hit`,`Target sequence length`,`Reference sequence length`)]
DT[,Percent := round(100*`HD for closest hit`/pmin(`Target sequence length`,`Reference sequence length`),1)]

# Get assembly sequence header names.
levels <- unique(DT$`PhyloFlash header`)[order(sapply(strsplit(unique(DT$`PhyloFlash header`),"_"),function (x) as.numeric(x[2])))]
DT[,`PhyloFlash header` := factor(`PhyloFlash header`,levels=levels)]

# Order Phyloflash headers
order.pf <- DT[,order(as.numeric(sapply(as.character(unique(`PhyloFlash header`)),function(x) strsplit(x,"_")[[1]][1])))]
levels.pf <- DT[,as.character(unique(`PhyloFlash header`))[order.pf]]
DT[,`PhyloFlash header` := factor(`PhyloFlash header`,levels=levels.pf)]
setorder(DT,"PhyloFlash header",Percent)

# Print out all header names into vector format which can be copied and pasted
# into this script.
# DT[,cat(paste0("c(\n",paste0(paste0('  "',sort(unique(`Closest hit to reference 16S`)),'"'),collapse = ",\n"),"\n)"),sep = "\n")]

# For renaming the header names into just first two words and ordering them
# in a way that the best matched reference sequence to each assembly sequence
# are neatly lined up.
# levels2 <- c(
#   "Caldimicrobium thiodismutans strain TF1",
#   "Candidatus Solibacter usitatus Ellin6076",
#   "Candidatus Thermochlorobacteriaceae bacterium GBChlB",
#   "Chloracidobacterium thermophilum B",
#   "Chloroflexus aggregans DSM 9485",
#   "Chloroflexus aurantiacus J-10-fl",
#   "Chloroflexus aurantiacus Y-400-fl",
#   "Chloroflexus islandicus strain isl-2",
#   "Chloroflexus sp. MS-G",
#   "Chloroflexus sp. Y-396-1",
#   "Elioraea tepidiphila DSM 17972",
#   "Gloeomargarita lithophora Alchichica-D10",
#   "Ignavibacterium album JCM 16511",
#   "Meiothermus ruber DSM 1279",
#   "Meiothermus silvanus DSM 9946",
#   "Meiothermus taiwanensis WR-220",
#   "Oscillochloris trichoides DG-6",
#   "Pseudanabaena sp. ABRG5-3",
#   "Pseudanabaena sp. PCC 7367",
#   "Pseudothermotoga hypogea DSM 11164 = NBRC 106472",
#   "Pseudothermotoga lettingae TMO",
#   "Roseiflexus castenholzii DSM 13941",
#   "Roseiflexus sp. RS-1",
#   "Synechococcus sp. JA-2-3B'a(2-13)",
#   "Synechococcus sp. JA-3-3Ab",
#   "Thermocrinis ruber strain DSM 23557",
#   "Thermodesulfovibrio sp. N1",
#   "Thermodesulfovibrio thiophilus DSM 17215",
#   "Thermus aquaticus Y51MC23",
#   "Thermus brockianus strain GE-1"
# )

# Order reference names so that the percents line up in the matrix
DT.PF <- DT[`Closest hit to reference 16S`%in%levels.pf,.(`PhyloFlash header`=`PhyloFlash header`[which.min(Percent)],Percent=Percent[which.min(Percent)]),.(`Closest hit to reference 16S`)]
setorder(DT.PF,`PhyloFlash header`,Percent)
levels3 <- DT.PF[,as.character(`Closest hit to reference 16S`)]
DT[,`Closest hit to reference 16S` := factor(
  `Closest hit to reference 16S`,
  levels=levels3
)]
levels.pf <- DT[,levels(`PhyloFlash header`)]

# Plots the heat map.
{
  gg <- ggplot(
    data = DT,
    mapping = aes(
      # x = `Closest hit to reference 16S`,
      # x = ordered(`Closest hit to reference 16S`,levels=levels2),
      x = ordered(`Closest hit to reference 16S`,levels=levels3),
      y = ordered(`PhyloFlash header`,levels=rev(levels.pf)),
      fill = Percent
    )
  )
  gg <- gg + geom_tile()
  gg <- gg + geom_text(
    mapping = aes(
      # x = `Closest hit to reference 16S`,
      x = ordered(`Closest hit to reference 16S`,levels=levels3),
      y = ordered(`PhyloFlash header`,levels=rev(levels.pf)),
      label = Percent
    )
  )
  gg <- gg + scale_fill_gradient2(
    low = "green",
    mid = "yellow",
    high = "red",
    midpoint = (min(DT$Percent)+max(DT$Percent))/2,
    limit = c(min(DT$Percent),max(DT$Percent))
  )
  gg <- gg + theme_minimal()
  gg <- gg + theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(hjust = 0.5)
  )
  gg <- gg + xlab("Closest hit to reference 16S")
  gg <- gg + ylab("Assembly")
  gg <- gg + ggtitle("Hamming Dissimilarity Percentage")
  print(gg)
}


