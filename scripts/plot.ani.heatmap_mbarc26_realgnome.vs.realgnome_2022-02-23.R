
### DESCRIPTION ################################################################

# Estimates pairwise ANI percentages, all-against-all, for the 26 reference
# genomes then plots a heat map with the ANI percent values.



### PREREQUISTES ###############################################################

# Need to install FastANI first to be able to calculate ANI. More info:
# https://github.com/ParBLiSS/FastANI

# To download and install:
# git clone https://github.com/ParBLiSS
# cd FastANI
# ./bootstrap.sh
# ./configure
# make



### HEADER CODE ################################################################

# Clean out environment
rm(list=ls())
gc()



### LIBRARIES ##################################################################

# Enhanced data.frame
library(data.table)

# Exploratory data analysis and data visualization for biological sequence (DNA
# and protein) data. Include also utilities for sequence data management under
# the ACNUC system.
library(seqinr)



### ANI CALCULATION ############################################################

# Set path to the 26 reference genomes and get the file names
dir.ref <- "data/mbarc26/ref/"
files.ref <- dir(dir.ref)

# Make a list of genus and species names for each organism
# cat(paste0("species -> c(\n",paste0('  "',gsub("_"," ",sub("\\.barrnap\\.fasta","",files.ref)),'"',collapse = ",\n"),"\n)"))
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

# Set path to new directory where updated FASTA files will be stored
dir.ref2 <- "data/mbarc26/ref2"
dir.create(dir.ref2,recursive = TRUE)

# Rename genome names to contain only genus and species, and save updated FASTA
# files to new path
idx <- 1
for (idx in 1:length(files.ref)) {
  file.in <- files.ref[idx]
  seqs <- read.fasta(file.path(dir.ref,file.in),as.string = TRUE,forceDNAtolower = FALSE)
  write.fasta(seqs,paste0(sub(" ","_",species[idx]),":",1:length(seqs)),file.path(dir.ref2,file.in))
}

# Get the names of updated FASTA files
files.ref2 <- dir(dir.ref2)

# Directory to save results
dir.res <- "results/fastANI"
dir.create(dir.res,recursive=TRUE)

# Calculate ANI all-against-all for all reference gnomes
# idx1 <- 1 # For debugging purposes
# idx2 <- 2 # For debugging purposes
DT <- data.table()
for (idx1 in 1:length(files.ref2)) {
  for (idx2 in 1:length(files.ref2)) {
    file.out <- paste0(
      dir.res,
      "/mbarc26/",
      gsub(" ","_",species[idx1]),
      ".",
      gsub(" ","_",species[idx2]),".fastani"
    )
    system(paste0("FastANI/fastANI -q data/mbarc26/ref2/",files.ref2[idx1]," -r data/mbarc26/ref2/",files.ref2[idx2]," -o ",file.out," --minFraction 0 -t 3 --fragLen 1000"))
    ani.est <- data.table(read.delim(file.out,header = FALSE,sep="\t"))[,round(V3,1)]
    DT <- rbindlist(list(DT,data.table(
      `Query ID` = species[idx1],
      `Subject ID` = species[idx2],
      `ANI estimate` = ani.est
    )))
  }
}
print(DT)

# To specify how the values in each variable should be sorted
DT[, `Query ID` := factor(`Query ID`,levels=species)]
DT[, `Subject ID` := factor(`Subject ID`,levels=species)]

# Save the results
saveRDS(DT,file.path(dir.res,"ani.mbarc26.rds"))



### PLOT #######################################################################

# Plot a heat map of an ANI matrix
{
  gg <- ggplot(
    data = DT,
    mapping = aes(
      x = `Query ID`,
      y = ordered(`Subject ID`,levels=rev(species)),
      fill = `ANI estimate`
    )
  )
  gg <- gg + geom_tile(color="white")
  gg <- gg + geom_text(
    mapping = aes(
      x = `Query ID`,
      y = ordered(`Subject ID`,levels=rev(species)),
      label = `ANI estimate`
    )
  )
  gg <- gg + scale_fill_gradient2(
    low = "red",
    mid = "yellow",
    high = "green",
    midpoint = (min(DT$`ANI estimate`)+max(DT$`ANI estimate`))/2,
    limit = c(min(DT$`ANI estimate`),max(DT$`ANI estimate`))
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
  gg <- gg + ggtitle("ANI estimate of reference genomes")
  print(gg)
}

