
### DESCRIPTION ################################################################

# Calculates pairwise BLAST identity percent values, all-against-all between 
# each of the 26 reference 16s sequences, then plots a heat map  with percent
# values.



### PREREQUISTES ###############################################################

# Need to install blast first locally (command example for Ubuntu machine):
# sudo apt-get install ncbi-blast+

# Packages to install first manually
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



### DATA PRE-PROCESSING ########################################################

# Set path to original Barrnap files (reference FASTA files)
dir.ref <- "data/barrnap"
files.ref <- dir(dir.ref)

# Set new path to directory where new reference FASTA files will be generated
dir.ref16s <- "data/blast/barrnap_16s"
files.ref16s <- dir(dir.ref16s)
dir.create(dir.ref,recursive = TRUE)

# Select only 16s sequences, rename the sequence using the names above, and
# save as individual reference FASTA files
idx <- 1
for (idx in 1:length(files.ref)) {
  file.in <- files.ref[idx]
  seqs <- read.fasta(file.path(dir.ref,file.in),as.string = TRUE,forceDNAtolower = FALSE)
  seqs <- seqs[grepl("16S",names(seqs))]
  new.names <- paste0(species[idx],":",sapply(strsplit(names(seqs),":"),function (x) x[length(x)]))
  file.out <- sub("\\.fasta",".16s.fasta",basename(file.in))
  write.fasta(seqs,sub(" ","_",new.names),file.path(dir.ref16s,file.out))
}

# Combine all sequences into one FASTA file
seqs <- do.call(c,lapply(file.path(dir.ref16s,files.ref16s),function(seq) {
  read.fasta(seq,as.string = TRUE,forceDNAtolower = FALSE)
}))
write.fasta(seqs,names(seqs),"data/blast/mockcommunity.barrnap.fasta")

# Create BLAST database from the combined FASTA file
makeblastdb("data/blast/mockcommunity.barrnap.fasta")



### CODE ######################################################################

# Load reference data
seqs <- readDNAStringSet(
  filepath = "data/blast/mockcommunity.barrnap.fasta"
)

# Load reference data as "blast database
bl <- blast(db="./data/blast/mockcommunity.barrnap.fasta")


tmp <- predict(bl,seqs,BLAST_args = '-task dc-megablast -num_threads 3 -evalue 10 -qcov_hsp_perc 90')

# Query sequences using BLAST
{ # Blast with e-value of 10
  cl <- data.table(predict(bl,seqs,BLAST_args = '-task dc-megablast -num_threads 3 -evalue 10'))
  cl[, `Query ID` := sapply(strsplit(QueryID,":"), function(x) sub("_"," ",x[1]))]
  cl[, `Query ID` := factor(`Query ID`,levels=species)]
  cl[, `Subject ID` := sapply(strsplit(SubjectID,":"), function(x) sub("_"," ",x[1]))]
  cl[, `Subject ID` := factor(`Subject ID`,levels=species)]
  DT10 <- cl[, .(
    `Percent identity` = round(max(Perc.Ident),1)
  ),by=.(`Query ID`,`Subject ID`)]
}
{ # Blast with e-value of 1e-30 (this gives NA values for some pairs)
  cl <- data.table(predict(bl,seqs,BLAST_args = '-task dc-megablast -num_threads 3 -evalue 1e-30'))
  cl[, `Query ID` := sapply(strsplit(QueryID,":"), function(x) sub("_"," ",x[1]))]
  cl[, `Query ID` := factor(`Query ID`,levels=species)]
  cl[, `Subject ID` := sapply(strsplit(SubjectID,":"), function(x) sub("_"," ",x[1]))]
  cl[, `Subject ID` := factor(`Subject ID`,levels=species)]
  DT <- cl[, .(
    `Percent identity` = round(max(Perc.Ident),1)
  ),by=.(`Query ID`,`Subject ID`)]
}

# Replace NA values calculated with e-value of 1e-30 with ones calculated with
# e-value of 10
DT <- rbindlist(list(DT,DT10[!paste(`Query ID`,`Subject ID`) %in% DT[,paste(`Query ID`,`Subject ID`)]]))


# Plot blast ID percentage matrix
{
  gg <- ggplot(
    data = DT,
    mapping = aes(
      x = `Query ID`,
      y = ordered(`Subject ID`,levels=rev(species)),
      fill = `Percent identity`
    )
  )
  gg <- gg + geom_tile(color="white")
  gg <- gg + geom_text(
    mapping = aes(
      x = `Query ID`,
      y = ordered(`Subject ID`,levels=rev(species)),
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


