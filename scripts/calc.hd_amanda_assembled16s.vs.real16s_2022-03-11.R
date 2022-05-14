
### DESCRIPTION ################################################################

# Basically, calculates Hamming distances between each assembly sequence to each
# reference sequence and stores the results in CSV format, for creating 
# heat maps later in other script.

# However, each reference FASTA file can contain  more than one 16s sequence.
# So first, for each sequence in the reference FASTA file, the smaller sequence
# "slides" across the larger sequence and calculates Hamming distances along
# the way. The minimum distance out from this sliding procedure is selected, 
# then the minimum out of all minimum distances from the assembly sequence to 
# all the sequences in that reference FAST file is selected.

# This applies to the assembly sequences assembled by PhyloFlash for the mock 
# community dataset and the 26 reference 16s sequences provided to me.



### HEADER CODE ################################################################

# Clean out environment
rm(list=ls())
gc()



### LIBRARIES ##################################################################

# Enhanced data.frame
library(data.table)

# Biological sequences retrieval and analysis
# install.packages('seqinr', repos='http://cran.us.r-project.org')
library(seqinr)

# Compute distance metrics between strings
# install.packages('stringdist', repos='http://cran.us.r-project.org')
library(stringdist)



### DEFINITIONS ################################################################

# Distance methods to use. We think Hamming distance is good choice to measure
# or ocunt nucleotide differences in DNA/RNA.
methods <- c(
  "osa", "lv", "dl", "hamming", "lcs", "qgram",
  "cosine", "jaccard", "jw", "soundex"
)
method <- methods[4]

# Path where assembled sequences are
path.in <- "data/amanda/assembly"

# Choose a path to store results
path.out <- "results/amanda/hd/"



### CODE #######################################################################

# Load assembly sequence file in FASTA format.
file.fasta.phyloflash <- "cami2_match_v4.all.final.fasta"
fasta.phyloflash <- read.fasta(file.path(path.in,file.fasta.phyloflash))
# length(fasta.phyloflash) # To check how many sequences there are

# Set directory which reference 16s rRNA sequence files in FASTA format are
# stored
files.fasta.ref <- "data/amanda/ref2/16s"

# Set up data table template to store results after Hamming distances are
# measured.
DT0 <- data.table(
  `PhyloFlash header` = character(),
  `Target sequence length` = integer(),
  `Closest hit to reference 16S` = character(),
  `Reference sequence length` = integer(),
  `HD for closest hit` = integer()
)

# Go through each sequence in assembly file
# idx1 <- 1 # For debugging purposes
for(idx1 in 1:length(fasta.phyloflash)) {
  
  # Select assembly file to focus n
  seq.pf <- fasta.phyloflash[[idx1]]
  
  # Create other data table template to use to store temporary info
  # of min Hamming distances of all reference 16s sequences to this 
  # assembly sequence
  DT <- data.table(
    `PhyloFlash header` = character(),
    `Target sequence length` = integer(),
    `Closest hit to reference 16S` = character(),
    `Reference sequence length` = integer(),
    `HD for closest hit` = integer()
  )
  
  # Go through each reference 16s sequences. Each reference FASTA file
  # can contain more than one sequence so we calculate Hamming distances
  # of each sequence in the reference file to the assembly sequence and take
  # only one that has min Hamming distance and store that to DT table.
  
  # idx2 <- 1 # For debugging purposes
  for (idx2 in 1:length(dir(files.fasta.ref,full.names = TRUE))) {
    
    # Read reference FASTA file
    file.fasta.ref <- dir(files.fasta.ref,full.names = TRUE)[idx2]
    fasta.ref <- read.fasta(file.fasta.ref)
    
    # The reference FASTA file may have 23s, 16s, or 5s sequences. We want only
    # 1ht 16S rRNA sequences .
    # fasta.ref <- fasta.ref[grepl("16S_rRNA",names(fasta.ref))]
    # length(fasta.ref) # To show how many 16s sequences that was stored in 
                        # the reference FASTA file.
    
    # The assembly sequence length might not be equal to the reference sequence
    # length. So, the smaller sequence will "slide" across the larger sequence
    # and the Hamming distance values are calculated along the way. Then only
    # the minimum Hamming distance are appended to the DT table.
    
    min.ref <- "" # Ref. sequence that currently has minimum Hamming distance.
    min.dist <- Inf # The current Hamming distance value.
    min.pos <- NA # The position of where the smaller sequence lies on the 
                  # larger sequence, which currently has min Hamming distance
    
    # Go through each sequence in the current reference FASTA file
    # idx3 <- 1
    for (idx3 in 1:length(fasta.ref)) {
      
      # Print info to monitor the whole progress.
      cat(idx1,"/",length(fasta.phyloflash)," : ",sep="")
      cat(idx2,"/",length(dir(files.fasta.ref,full.names = TRUE))," : ",sep="")
      cat(idx3,"/",length(fasta.ref),"\n",sep="")
      
      # Get the current name and sequence from the list to focus on 
      seq.ref <- fasta.ref[[idx3]]
      name.ref <- attr(seq.ref,"name")
      name.ref <- gsub(":.*","",gsub("_"," ",name.ref))
      
      # Delegate the assembly sequence as the smaller or larger sequence, and
      # vice versa for the reference sequence.
      if (length(seq.ref) < length(seq.pf)) {
        short <- seq.ref
        long <- seq.pf
      } else {
        short <- seq.pf
        long <- seq.ref
      }
      n <- length(short)
      N <- length(long)
      
      # Store the position of the smaller sequence on the larger sequence, and
      # calculate the Hamming distance.
      current.pos <- which.min(sapply(1:(N-n+1),function(i) sum(stringdist(short, long[i:(n+i-1)], method = method))))
      current.dist <- min(sapply(1:(N-n+1),function(i) sum(stringdist(short, long[i:(n+i-1)], method = method))))
      
      # Compare the current Hamming distance to the current minimum Hamming
      # distance and replace if current distance value is smaller than the
      # current min value.
      if (current.dist < min.dist) {
        min.pos <- current.pos
        min.dist <- current.dist
        # min.ref <- seq.ref
        min.ref <- name.ref
      }
      
    }
    
    # Then append the info that corresponds to the min Hamming distance to the
    # DT table.
    DT2 <- data.table(
      `PhyloFlash header` = gsub("cami2_match_v4\\.PFspades_","",attr(seq.pf,"name")),
      `Target sequence length` = length(seq.pf),
      # `Closest hit to reference 16S` = gsub("_"," ",sub("\\.barrnap\\.fasta","",basename(file.fasta.ref))),
      `Closest hit to reference 16S` = name.ref,
      `Reference sequence length` = length(seq.ref),
      `HD for closest hit` = min.dist
    )
    DT <- rbindlist(list(DT,DT2))
    
  }
  
  # Save current work in CSV format in case if R session crashes.
  fwrite(DT,file.path(path.out,paste0("hd",idx1,".csv")))
  
  # Also, periodically back up the CSV files on cloud in case if this machine's
  # memory gets overloaded and you have no other choice to restart this machine,
  # losing all work that was stored in the memory :(
  # (This only works if you have Google Cloud Storage and has 
  # gsutil installed)
  # system(paste0("gsutil cp hd",idx1,".csv gs://BACKUP_",Sys.Date(),"/"))
  
  # Append to global datatable
  DT0 <- rbindlist(list(DT,DT0))
  
}

# Save table containing all minimum Hamming distance values and the
# corresponding names of assembly and reference sequences
fwrite(DT0,file.path(path.out,paste0("hd_",Sys.Date(),".csv")))
print(DT0)





