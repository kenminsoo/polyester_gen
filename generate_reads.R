#check if needed packages exist
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!require("polyester", quietly = TRUE))
    BiocManager::install("polyester")

if (!require("Biostrings", quietly = TRUE))
    BiocManager::install("Biostrings")

library(polyester)
library(Biostrings)

#this initial experiment creates simulated sequencing reads on the illumina 5 system with 
#enter location of the fasta transcript file
fasta_loc <- '/out/subset_transcripts.fasta'

fasta_file <- "/Users/kenminsoo/Desktop/Projects/JABSOM/Nakatsu/projects/polyester/out/subset_transcripts.fasta"

#our current subset is at 1400, but change as please
num_transcripts <- count_transcripts(fasta_file)

set.seed(482)
fold_change_values <- sample(c(0.5, 1, 2, 4), size = num_transcripts, prob=c(0.025, 0.9, 0.05, 0.025), replace=TRUE)
fold_changes <- matrix(fold_change_values, nrow=num_transcripts)

transcripts <- gsub("'", "", names(readDNAStringSet(fasta_file)))

simulate_experiment(fasta_file, transcriptid=transcripts, fold_changes = fold_changes, outdir = 'out/sim_reads',seed = 313, readlen = 35, error_rate = 0, paired = FALSE, num_reps=c(5,5))
