
### DESCRIPTION ################################################################

# Calculates pairwise BLAST identity percent values, all-against-all between 
# assembled 16Ss and real 16Ss, then plots a heat map with percent values.



### PREREQUISTES ###############################################################

# Need to install blast first locally (command example for Ubuntu machine):
# sudo apt-get install ncbi-blast+

# Packages to install first manually
# install.packages("devtools") # Already installed by default
# devtools::install_bioc("BiocGenerics")
# devtools::install_bioc("Biostrings")
# devtools::install_github("mhahsler/rBLAST")



### HEADER CODE ################################################################

# Clean out environment
rm(list=ls())
gc()



### LIBRARIES ##################################################################

# Enhanced data.frame
library(data.table)

# Memory efficient string containers, string matching algorithms, and other
# utilities, for fast manipulation of large biological sequences or sets of
# sequences
library(Biostrings)

# Interfaces the Basic Local Alignment Search Tool (BLAST) to search genetic
# sequence data bases with the Bioconductor infrastructure
library(rBLAST)

# Exploratory data analysis and data visualization for biological sequence (DNA
# and protein) data. Include also utilities for sequence data management under
# the ACNUC system.
library(seqinr)

# Elegant data visualizations using grammar of graphics
library(ggplot2)



### DEFINITIONS ################################################################

# Names of 26 organisms
species <- c(
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



### CODE ######################################################################

# Set path assembly FASTA files
dir.ass <- "data/mbarc26/assembly"
file.ass <- "LIB.all.final.fasta"

# Set path to original Barrnap files (reference FASTA files)
dir.ref <- "data/mbarc26/blast"
file.ref <- "mockcommunity.barrnap.fasta"

# Load assembly data
seqs <- readDNAStringSet(filepath = file.path(dir.ass,file.ass))
names(seqs)

# Load reference data as "blast database
bl <- blast(db=file.path(dir.ref,file.ref))

# Order the assembly names and species in a way that the best matched reference
# sequence to each assembly sequence are neatly lined up.
# cat(paste0("levels.ass <- c(\n",paste0(paste0('  "',names(seqs),'"'),collapse=",\n"),'\n)'))
levels.ass <- c(
  "LIB.PFspades_1_46.766972",
  "LIB.PFspades_2_170.172388",
  "LIB.PFspades_3_223.359275",
  "LIB.PFspades_4_439.483061",
  "LIB.PFspades_5_82.268872",
  "LIB.PFspades_6_78.150487",
  "LIB.PFspades_7_116.642510",
  "LIB.PFspades_8_166.545678",
  "LIB.PFspades_9_40.886155",
  "LIB.PFspades_10_48.414023",
  "LIB.PFspades_11_26.585160",
  "LIB.PFspades_12_69.894048",
  "LIB.PFspades_13_157.423986",
  "LIB.PFspades_14_54.010634",
  "LIB.PFspades_15_51.667115",
  "LIB.PFspades_16_26.310202",
  "LIB.PFspades_18_1021.366635",
  "LIB.PFspades_19_15.066084",
  "LIB.PFspades_20_16.678230",
  "LIB.PFspades_21_3.386923",
  "LIB.PFspades_22_355.590116",
  "LIB.PFspades_23_419.524158",
  "LIB.PFspades_24_594.813644"
)
levels.ref <- c(
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

# Query sequences using BLAST
{ # Blast with e-value of 10
  cl <- data.table(predict(bl,seqs,BLAST_args = '-task dc-megablast -num_threads 3 -evalue 10'))
  cl[, `Query ID` := factor(QueryID,levels=levels.ass)]
  cl[, `Subject ID` := sapply(strsplit(SubjectID,":"), function(x) sub("_"," ",x[1]))]
  cl[, `Subject ID` := factor(`Subject ID`,levels=levels.ref)]
  DT10 <- cl[, .(
    `Percent identity` = round(max(Perc.Ident),1)
  ),by=.(`Query ID`,`Subject ID`)]
}
{ # Blast with e-value of 1e-30 (this gives NA values for some pairs)
  cl <- data.table(predict(bl,seqs,BLAST_args = '-task dc-megablast -num_threads 3 -evalue 1e-30'))
  cl[, `Query ID` := factor(QueryID,levels=levels.ass)]
  cl[, `Subject ID` := sapply(strsplit(SubjectID,":"), function(x) sub("_"," ",x[1]))]
  cl[, `Subject ID` := factor(`Subject ID`,levels=levels.ref)]
  DT <- cl[, .(
    `Percent identity` = round(max(Perc.Ident),1)
  ),by=.(`Query ID`,`Subject ID`)]
}

# Replace NA values calculated with e-value of 1e-30 with ones calculated with
# e-value of 10
DT <- rbindlist(list(DT,DT10[!paste(`Query ID`,`Subject ID`) %in% DT[,paste(`Query ID`,`Subject ID`)]]))
# DT[,min(`Percent identity`)]

# Plot blast ID percentage matrix
{
  gg <- ggplot(
    data = DT,
    mapping = aes(
      x = `Subject ID`,
      y = ordered(`Query ID`,levels=rev(levels.ass)),
      fill = `Percent identity`
    )
  )
  gg <- gg + geom_tile(color="white")
  gg <- gg + geom_text(
    mapping = aes(
      x =`Subject ID`,
      y = ordered(`Query ID`,levels=rev(levels.ass)),
      label = `Percent identity`
    )
  )
  gg <- gg + scale_fill_gradient2(
    low = "red",
    mid = "yellow",
    high = "green",
    # midpoint = (min(DT$`Percent identity`)+max(DT$`Percent identity`))/2,
    # limit = c(min(DT$`Percent identity`),max(DT$`Percent identity`))
    midpoint = (65+100)/2,
    limit = c(65,100)
  )
  gg <- gg + ylab("Subject ID")
  gg <- gg + theme_minimal()
  gg <- gg + theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(hjust = 0.5)
  )
  gg <- gg + ggtitle("BLAST ID Percentage")
  # print(gg)
}

ggsave(
  filename = paste0("blast.id.perc.heatmap_mbarc26_assembled16s.vs.real16s_",Sys.Date(),".png"),
  path = "visualizations",
  width = 13,
  height = 7
)


