source("Code/Tidyverse_function.R")


Variantes_5_prd <- read_lines("/home/pedro.paixao/Code/Previsoes/Variantes_5_porcento_prd.fasta")
Variantes_5_prd <-str_subset(Variantes_5_prd, "^CCCCCCCCCCCCCCCCC")
Variantes_5_prd <- strsplit(Variantes_5_prd, "")

EGFR_fasta_prd <- readLines("/home/pedro.paixao/Code/Previsoes/EGFR_prd.fasta")
EGFR_fasta_prd <- str_subset(EGFR_fasta_prd, "^CCCCCCCCCCCCCCCCC")
EGFR_fasta_prd <- strsplit(EGFR_fasta_prd, "")

Variantes_100_prd <- str_subset(Teste, "^CCCCCCCCCCCCCCCCC")
Variantes_100_prd <- strsplit(Variantes_100_prd, "")





Comparacao <- mean(Variantes_100_prd == EGFR_fasta_prd) * 100
class(Variantes_100_prd)
class(EGFR_fasta_prd)
Similaridade <- vector()
Similaridade_5 <- vector()

i <- 1

for (i in 1:length(Variantes_5_prd)){
  Similaridade_5[[i]] <- mean(Variantes_5_prd[[i]] == EGFR_fasta_prd[[1]])
  i <- i + 1
}

hist(Similaridade_5)
