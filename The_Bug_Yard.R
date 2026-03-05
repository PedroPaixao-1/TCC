source("Code/Tidyverse_function.R")

script_py <- "/home/pedro.paixao/Code/s4pred-main/run_model.py"
arquivo_entrada <- "/home/pedro.paixao/Code/Fastas/Variantes/Variantes_100_porcento_SED.fasta"
arquivo_saida <- "/home/pedro.paixao/Code/Variantes_100_porcento_prd.fasta"

Texto <- paste("~/anaconda3/envs/esm2-env/bin/python3", script_py, "-t fas -z", arquivo_entrada, ">", arquivo_saida)

# Execute
system(Texto, wait = TRUE)

EGFR_fasta_prd <- read.fasta("/home/pedro.paixao/Code/Previsoes/EGFR_prd.fasta")

Variantes_5_porcento_fasta_pred <- read.fasta("/home/pedro.paixao/Code/Previsoes/Variantes_5_porcento_prd.fasta")

Variantes_100_porcento_fasta_pred$`EGFR|CHAIN` <- read.fasta("/home/pedro.paixao/Code/Previsoes/Variantes_100_porcento_prd.fasta")
Variantes_100_porcento_fasta_pred$`EGFR|CHAIN`[1][1]
attr(Variantes_100_porcento_fasta_pred,"name")
Comparacao <- list()

Comparacao <- Variantes_5_porcento_fasta_pred[[i]][614:length(EGFR_fasta_prd)] %in% 
  EGFR_fasta_prd[614:length(EGFR_fasta_prd[[1]])]

Comparacao_diferencas <- sum(unlist(Variantes_5_porcento_fasta_pred[[1]]) != EGFR_fasta_prd)
attr(Variantes_100_porcento_fasta_pred,"class")

Vetor_1 <- list(1,2,3,4,5,6,7,8)
a <- Variantes_5_porcento_fasta_pred %in% EGFR_fasta_prd 


# for (i in 1:length(Variantes_5_porcento_fasta_pred)){
#   if (Variantes_5_porcento_fasta_pred %in% EGFR_fasta_prd == TRUE){
#     Comparacao[[i]] <- 0
#   } else {
    sum(ulist(Variantes_5_porcento_fasta_pred) != EGFR_fasta_prd)
#   } return(Comparacao)
  
  
  # Comparacao[[i]] <- sum(Variantes_5_porcento_fasta_pred[[i]] != EGFR_fasta_prd[[1]])
  # return(Comparacao)
  # # i <- i +1
# }

Comparacao <- Variantes_5_porcento_fasta_pred %in% EGFR_fasta_prd

