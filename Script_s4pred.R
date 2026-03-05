source("/home/pedro.paixao/Code//FunctionsV2.R")

Sequencia <- pdbseq(Arquivo_pdb)
Sequencia <- Sequencia[1:613]
Sequencia <- paste(Sequencia, collapse = "")

Header_EGFR <- paste(">EGFR|CHAIN A|Controle")

EGFR_fasta <- paste(Header_EGFR, Sequencia, sep = "\n")
writeLines(unlist(EGFR_fasta), "EGFR.fasta")

Texto <- paste(" python /home/pedro.paixao/Code/s4pred-main/run_model.py -t fas -z", Variantes.fasta)

system(Texto, wait = TRUE)

Previsoes <- list.files("/home/pedro.paixao/Code/s4pred-main/preds", pattern = "\\.fas$", full.names = T) 

Previsoes <- lapply(Previsoes, read.fasta)

# #lê as preds
# Dados <- lapply(Arquivos_pred, read.fasta)
# 
# return(Dados)
# }


# length(seq) + 1



#-----------------------------------------------Teste--------------------------------------#

Fasta_teste <- read.fasta("/home/pedro.paixao/Code/s4pred-main/preds/s4_out_0.fas")

Header <- attr(Fasta_teste,"name")

Similaridade_df <-data_frame(Header)


script_py <- "/home/pedro.paixao/Code/s4pred-main/run_model.py"
arquivo_entrada <- "/home/pedro.paixao/Code/Fastas/Variantes/Variantes_5_porcento.fasta"
arquivo_saida <- "/home/pedro.paixao/Code/Variantes_5_porcento_prd.fasta"

Texto <- paste("~/anaconda3/envs/esm2-env/bin/python3", script_py, "-t fas -z", arquivo_entrada, ">", arquivo_saida)

# Execute
system(Texto, wait = TRUE)
