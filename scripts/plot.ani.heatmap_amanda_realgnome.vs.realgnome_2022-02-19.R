
### DESCRIPTION ################################################################

# Estimates pairwise ANI percentages, all-against-all, of the reference whole
# genomes for Amanda's dataset, then plots a heat map with the ANI percent
# values.



### PREREQUISTES ###############################################################

# Need to install FastANI first to be able to calculate ANI. More info:
# https://github.com/ParBLiSS/FastANI

# To download and install:
# git clone https://github.com/ParBLiSS/FastANI
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
dir.ref <- "data/amanda/ref/whole"
files.ref <- dir(dir.ref)

# Make a list of genus and species names for each organism
# cat(paste0("species -> c(\n",paste0('  "',gsub("_"," ",sub("\\.barrnap\\.fasta","",files.ref)),'"',collapse = ",\n"),"\n)"))
scientific.names <- c(
  "Candidatus Caldatribacterium saccharofermentans isolate SpSt-82",
  "Candidatus Thermochlorobacteriaceae bacterium GBChlB",
  "Synechococcus sp. JA-3-3Ab",
  "Synechococcus sp. JA-2-3B'a(2-13)",
  "Candidatus Solibacter usitatus Ellin6076",
  "Roseiflexus sp. RS-1",
  "Roseiflexus castenholzii DSM 13941",
  "Pseudothermotoga lettingae TMO",
  "Chloroflexus aurantiacus J-10-fl",
  "Chloroflexus aggregans DSM 9485",
  "Chloroflexus aurantiacus Y-400-fl",
  "Meiothermus ruber DSM 1279",
  "Meiothermus silvanus DSM 9946",
  "Chloracidobacterium thermophilum B",
  "Ignavibacterium album JCM 16511",
  "Pseudanabaena sp. PCC 7367",
  "Pseudothermotoga hypogea DSM 11164 = NBRC 106472",
  "Oscillochloris trichoides DG-6",
  "Caldimicrobium thiodismutans strain TF1",
  "Pseudanabaena sp. ABRG5-3",
  "Elioraea tepidiphila DSM 17972",
  "Thermodesulfovibrio thiophilus DSM 17215",
  "Chloroflexus sp. Y-396-1",
  "Thermocrinis ruber strain DSM 23557",
  "Thermus aquaticus Y51MC23",
  "Thermus brockianus strain GE-1",
  "Gloeomargarita lithophora Alchichica-D10",
  "Meiothermus taiwanensis WR-220",
  "Chloroflexus sp. MS-G",
  "Chloroflexus islandicus strain isl-2",
  "Thermodesulfovibrio sp. N1"
)
names(scientific.names) <- files.ref

# Set path to new directory where updated FASTA files will be stored
dir.ref2 <- "data/amanda/ref2/whole"
dir.create(dir.ref2,recursive = TRUE)

# Rename genome names to contain only genus and species, and save updated FASTA
# files to new path
# idx <- 1 # For debugging purposes
for (idx in 1:length(files.ref)) {
  file.in <- files.ref[idx]
  seqs <- read.fasta(file.path(dir.ref,file.in),as.string = TRUE,forceDNAtolower = FALSE)
  write.fasta(seqs,paste0(sub(" ","_",scientific.names[file.in]),":",1:length(seqs)),file.path(dir.ref2,file.in))
}

# Get the names of updated FASTA files
files.ref2 <- dir(dir.ref2)

# Directory to save results
dir.res <- "results/fastANI"
dir.create(file.path(dir.res,"amanda"),recursive=TRUE)

# Calculate ANI all-against-all for all reference gnomes
# idx1 <- 1 # For debugging purposes
# idx2 <- 2 # For debugging purposes
DT <- data.table()
for (idx1 in 1:length(files.ref2)) {
  for (idx2 in 1:length(files.ref2)) {
    file.out <- paste0(
      dir.res,
      "/amanda/",
      gsub("[()']","",gsub(" ","_",scientific.names[idx1])),
      ",",
      gsub("[()']","",gsub(" ","_",scientific.names[idx2])),
      ".fastani"
    )
    system(paste0("FastANI/fastANI -q ",file.path(dir.ref2,files.ref2[idx1])," -r ",file.path(dir.ref2,files.ref2[idx2])," -o ",file.out," --minFraction 0 -t 3"))
    ani.est <- data.table(read.delim(file.out,header = FALSE,sep="\t"))[,round(V3,1)]
    DT <- rbindlist(list(DT,data.table(
      `Query ID` = scientific.names[idx1],
      `Subject ID` = scientific.names[idx2],
      `ANI estimate` = ani.est
    )))
  }
}
print(DT)

# To specify how the values in each variable should be sorted
DT[, `Query ID` := factor(`Query ID`,levels=sort(scientific.names))]
DT[, `Subject ID` := factor(`Subject ID`,levels=sort(scientific.names))]

# Save ANI results
saveRDS(DT,file.path(dir.res,"ani.amanda.rds"))



### PLOT #######################################################################

# Plot a heat map of an ANI matrix
{
  gg <- ggplot(
    data = DT,
    mapping = aes(
      x = `Query ID`,
      y = ordered(`Subject ID`,levels=rev(sort(scientific.names))),
      fill = `ANI estimate`
    )
  )
  gg <- gg + geom_tile(color="white")
  gg <- gg + geom_text(
    mapping = aes(
      x = `Query ID`,
      y = ordered(`Subject ID`,levels=rev(scientific.names)),
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


