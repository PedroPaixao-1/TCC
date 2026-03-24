source("/home/pedro.paixao/Code/FunctionsV2.R")

Variantes_100 <- readLines("/home/pedro.paixao/Code/Fastas/Variantes/Variantes_100_porcento_SED.fasta")


Variantes_5_prd <- read_lines("/home/pedro.paixao/Code/Previsoes/Variantes_5_porcento_prd.fasta")
Variantes_5_prd <-str_subset(Variantes_5_prd, "^CCCCCCCCCCCCCCCCC")
Variantes_5_prd <- strsplit(Variantes_5_prd, "")

EGFR_fasta_prd <- readLines("/home/pedro.paixao/Code/Previsoes/EGFR_prd.fasta")
EGFR_fasta_prd <- str_subset(EGFR_fasta_prd, "^CCCCCCCCCCCCCCCCC")
EGFR_fasta_prd <- strsplit(EGFR_fasta_prd, "")

Variantes_100_prd <- readLines("/home/pedro.paixao/Code/Previsoes/Variantes_100_porcento_prd.fasta")
Variantes_100_prd <- str_subset(Variantes_100_prd, "^CCCCCCCCCCCCCCCCC")
Variantes_100_prd <- strsplit(Variantes_100_prd, "")

Similaridade_100 <- vector()
Similaridade_5 <- vector()

i <- 1

for (i in 1:length(Variantes_100_prd)){
  Similaridade_100[[i]] <- mean(Variantes_100_prd[[i]] == EGFR_fasta_prd[[1]])
  i <- i + 1
}

hist(Similaridade_100)
min(Similaridade_100)
