# polyester_gen

# This is a script to generate reads with package polyester which simulates DEG reads

Requires an annotation file and a fasta file with the transcript sequences. 

# Usage
0) Clean annotation with the python script. 
1) Source the R file in the terminal after choosing between interactive and default. 
Interactive will allow more customization of parameters while default will generate reads with the following parameters:
read length: 35
pair ended: false
error model: uniform, 0%
stranded: true
2) Test favorite tools. 

# Features
## Customizable config with the main parameters from the package. 
## Choose where you want your reads to be put

## Cohort Reads 1
Generation features

- No merged features
- No duplicated features (i.e. large repeats) 
