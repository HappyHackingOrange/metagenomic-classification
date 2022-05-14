
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
DT <- rbindlist(lapply(dir("results/hd",pattern=".csv",full.names = TRUE),fread))

# Calculate Hamming dissimilarity percentage.
DT[,Percent := round(100*`HD for closest hit`/pmin(`Target sequence length`,`Reference sequence length`),1)]

# Get assembly sequence header names.
levels <- unique(DT$`PhyloFlash header`)[order(sapply(strsplit(unique(DT$`PhyloFlash header`),"_"),function (x) as.numeric(x[2])))]
DT[,`PhyloFlash header` := factor(`PhyloFlash header`,levels=levels)]

# Print out all header names into vector format which can be copied and pasted
# into this script.
# DT[,cat(paste0("c(\n",paste0(paste0('  "',sort(unique(`Closest hit to reference 16S`)),'"'),collapse = ",\n"),"\n)"),sep = "\n")]

# For renaming the header names into just first two words and ordering them
# in a way that the best matched reference sequence to each assembly sequence
# are neatly lined up.
levels2 <- c(
  "Clostridium perfringensATCC 13124",
  "Clostridium thermocellumVPI 7372 ATCC 27405",
  "Coraliomargarita akajimensis DSM 45221",
  "Corynebacterium glutamicum ATCC 13032",
  "Desulfosporosinus acidiphilus SJ4 DSM 22704",
  "Desulfosporosinus meridiei DSM 13257",
  "Desulfotomaculum gibsoniae DSM 7213",
  "E.coli K12 ATCC 700926",
  "Echinicola vietnamensis DSM 17526",
  "Fervidobacterium pennivorans DSM 9078",
  "Frateuria aurantia DSM 6220",
  "Halovivax ruber XH-70",
  "Hirschia baltica ATCC 49814",
  "Meiothermus silvanus DSM 9946",
  "Natronobacterium gregoryi SP2",
  "Natronococcus occultus DSM 3396",
  "Nocardiopsis dassonvillei DSM 43111",
  "Olsenella uli DSM 7084",
  "Pseudomonas stutzeri RCH2",
  "Salmonella bongori NCTC 12419",
  "Salmonella enterica subsp. arizonae serovar 62 z4 z23 - strain RSK2980",
  "Segniliparus rotundus DSM 44985",
  "Spirochaeta smaragdinae DSM 11293",
  "Streptococcus pyogenes M1 GAS",
  "Terriglobus roseus DSM 18391",
  "Thermobacillus composti KWC4"
)
labels2 <- c(
  "Clostridium perfringens",
  "Clostridium thermocellum",
  "Coraliomargarita akajimensis",
  "Corynebacterium glutamicum",
  "Desulfosporosinus acidiphilus",
  "Desulfosporosinus meridiei",
  "Desulfotomaculum gibsoniae",
  "Escherichia coli",
  "Echinicola vietnamensis",
  "Fervidobacterium pennivorans",
  "Frateuria aurantia",
  "Halovivax ruber",
  "Hirschia baltica",
  "Meiothermus silvanus",
  "Natronobacterium gregoryi",
  "Natronococcus occultus",
  "Nocardiopsis dassonvillei",
  "Olsenella uli",
  "Pseudomonas stutzeri",
  "Salmonella bongori",
  "Salmonella enterica",
  "Segniliparus rotundus",
  "Spirochaeta smaragdinae",
  "Streptococcus pyogenes",
  "Terriglobus roseus",
  "Thermobacillus composti"
)
levels3 <- c(
  "Streptococcus pyogenes", #1
  "Hirschia baltica", #2
  "Frateuria aurantia", #3
  "Thermobacillus composti", #4
  "Coraliomargarita akajimensis", #5
  "Pseudomonas stutzeri", #6
  "Natronobacterium gregoryi", #7
  "Spirochaeta smaragdinae", #8
  "Clostridium perfringens", #9
  "Salmonella enterica", #10
  "Segniliparus rotundus", #11
  "Olsenella uli", #12
  "Natronococcus occultus", #13
  "Terriglobus roseus", #14
  "Halovivax ruber", #15
  "Corynebacterium glutamicum", #16
  "Desulfosporosinus acidiphilus", #18
  "Clostridium thermocellum", #19
  "Echinicola vietnamensis", #20
  "Meiothermus silvanus", #22
  "Desulfosporosinus meridiei", #23
  "Desulfotomaculum gibsoniae", #24
  "Escherichia coli",
  "Salmonella bongori",
  "Fervidobacterium pennivorans",
  "Nocardiopsis dassonvillei"
)
DT[,`Closest hit to reference 16S` := factor(
  `Closest hit to reference 16S`,levels=levels2,
  labels = labels2
)]
setorder(DT,"PhyloFlash header",Percent)

# Plots the heat map.
{
  gg <- ggplot(
    data = DT,
    mapping = aes(
      x = ordered(`Closest hit to reference 16S`,levels=levels3),
      y = ordered(`PhyloFlash header`,levels=rev(levels)),
      fill = Percent
    )
  )
  gg <- gg + geom_tile()
  gg <- gg + geom_text(
    mapping = aes(
      x = ordered(`Closest hit to reference 16S`,levels=levels3),
      y = ordered(`PhyloFlash header`,levels=rev(levels)),
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


