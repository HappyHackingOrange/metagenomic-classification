
### DESCRIPTION ################################################################

# Calculates pairwise BLAST identity percent values, all-against-all between 
# each of the Amanda's reference 16s sequences, then plots a heat map with
# percent values.



### PREREQUISTES ###############################################################

# Need to install blast first locally (command example for Ubuntu machine):
# sudo apt-get install ncbi-blast+

# Packages to install first manually:
# install.packages("devtools") # Already installed by default
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

# Species and strain names of Amanda's 31 reference organisms
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

# Files not to use
files.dont.use <- c(
  "NZ_JAAFGM000000000.fna.fna",
  "NZ_NKXO00000000.fna.fna",
  "NZ_QMDI00000000.fna.fna",
  "NZ_SIJK00000000.fna.fna",
  "NZ_VFAC00000000.fna.fna",
  "NZ_VIKA00000000.fna.fna",
  "PHFL00000000.fna.fna"
)



### DATA PRE-PROCESSING ########################################################

# Set path to original Barrnap files (reference FASTA files)
dir.ref <- "data/amanda/ref/16s"
files.ref <- dir(dir.ref)
files.ref <- files.ref[!files.ref %in% files.dont.use]

# Set new path to directory where new reference FASTA files will be generated
dir.ref16s <- "data/amanda/ref2/16s"
dir.create(dir.ref16s,recursive = TRUE)

# Select only 16s sequences, rename the sequence using the names above, and
# save the updated FASTA files to new directory path
# idx <- 2 # For debugging purposes
for (idx in 1:length(files.ref)) {
  file.in <- files.ref[idx]
  if (file.in %in% files.dont.use) next
  seqs <- read.fasta(file.path(dir.ref,file.in),as.string = TRUE,forceDNAtolower = FALSE)
  seqs <- seqs[grepl("16S",names(seqs))]
  # names(seqs)
  if (length(seqs)==0) next
  new.names <- paste0(scientific.names[idx],":",sapply(strsplit(names(seqs),":"),function (x) x[length(x)]))
  new.names <- gsub(" ","_",new.names)
  # print(new.names)
  file.out <- sub("\\.fna\\.fna",".16s.fasta",basename(file.in))
  write.fasta(seqs,new.names,file.path(dir.ref16s,file.out))
}

# Get the updated file names
files.ref16s <- dir(dir.ref16s)

# Combine all sequences into one FASTA file
seqs <- do.call(c,lapply(file.path(dir.ref16s,files.ref16s),function(seq) {
  read.fasta(seq,as.string = TRUE,forceDNAtolower = FALSE)
}))
dir.blast <- "data/amanda/blast"
file.blast <- "amanda.16s.fasta"
write.fasta(seqs,names(seqs),file.path(dir.blast,file.blast))

# Create BLAST database from the combined FASTA file
makeblastdb(file.path(dir.blast,file.blast))



### CODE ######################################################################

# Load reference data
seqs <- readDNAStringSet(
  filepath = file.path(dir.blast,file.blast)
)

# Load reference data as "blast database
bl <- blast(db=file.path(dir.blast,file.blast))

# Query sequences using BLAST
{ # Blast with default e-value of 10
  cl <- data.table(predict(bl,seqs,BLAST_args = '-task dc-megablast -num_threads 3'))
  cl[, `Query ID` := sapply(strsplit(QueryID,":"), function(x) gsub("_"," ",x[1]))]
  cl[, `Query ID` := factor(`Query ID`,levels=sort(scientific.names))]
  cl[, `Subject ID` := sapply(strsplit(SubjectID,":"), function(x) gsub("_"," ",x[1]))]
  cl[, `Subject ID` := factor(`Subject ID`,levels=sort(scientific.names))]
  DT <- cl[, .(
    `Percent identity` = round(max(Perc.Ident),1)
  ),by=.(`Query ID`,`Subject ID`)]
  setorder(DT,"Query ID","Subject ID")
}

# Plot blast ID percentage matrix
{
  gg <- ggplot(
    data = DT,
    mapping = aes(
      x = `Query ID`,
      y = ordered(`Subject ID`,levels=rev(sort(scientific.names))),
      fill = `Percent identity`
    )
  )
  gg <- gg + geom_tile(color="white")
  gg <- gg + geom_text(
    mapping = aes(
      x = `Query ID`,
      y = ordered(`Subject ID`,levels=rev(sort(scientific.names))),
      label = `Percent identity`
    )
  )
  gg <- gg + scale_fill_gradient2(
    low = "red",
    mid = "yellow",
    high = "green",
    midpoint = (min(DT$`Percent identity`)+max(DT$`Percent identity`))/2,
    limit = c(min(DT$`Percent identity`),max(DT$`Percent identity`))
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
  print(gg)
}


