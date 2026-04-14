source("/home/pedro.paixao/Code//FunctionsV2.R")

Arquivo_pdb <- read.pdb("/home/pedro.paixao/Code/egfr_padronizado2.pdb")
Arquivo_pdb_ca <- atom.select(Arquivo_pdb, "calpha")
Arquivo_pdb_ca <- trim.pdb(Arquivo_pdb, inds = Arquivo_pdb_ca)

Contatos_EGFR <- Calcular_Contato(Arquivo_pdb_ca, 8)

Amostras <- 200

Sequencia <- pdbseq(Arquivo_pdb)
Sequencia <- Sequencia[1:613]
Sequencia <- paste(Sequencia, collapse = "")

Temperatura <- 1.5

Porcentagem <- c(20,30,40,50,60,70,80,90)
Variantes <- list()
Multi_fasta <- list()

j <- 1
k <- 1


length(Porcentagem)

for (k in 1:length(Porcentagem)){
  j <- 1
  while (j <= 10) {
  
  
    Posicao <- sample(Contatos_EGFR, size = ceiling(Porcentagem[[k]]/100*length(Contatos_EGFR)))
    unlist(Posicao)
  
    Comando <- paste("/home/conda/condabin/conda run -n esm2","python /home/pedro.paixao/Code/generate_sequence_esm2.py", "--sequence", Sequencia, "--position", Posicao,"--temperature", Temperatura)

    system(Comando, wait = TRUE)
  
    Saida <- readLines("completed_sequence.txt")

    if (Saida %in% Variantes){
    next
  }
    Header <- paste(">EGFR|CHAIN A|Variante",j)
    Multi_fasta[[j]] <- paste(Header,Saida, sep = "\n")
    Arquivo_fasta <- paste("Variante",Porcentagem[[k]],".fasta", sep = "_")
    writeLines(unlist(Multi_fasta), Arquivo_fasta)
    j = j + 1
  }
  k = k + 1
}
# -----------------------------------------------------------------------------------
library(stringr)

base_path <- "C:/Users/Pichau/OneDrive/Documentos/Code/TCC-main/Fastas/Variantes/"

# Leitura
Var_10  <- readLines(paste0(base_path, "Variante_10_.fasta"))
Var_20  <- readLines(paste0(base_path, "Variante_20_.fasta"))
Var_30  <- readLines(paste0(base_path, "Variante_30_.fasta"))
Var_40  <- readLines(paste0(base_path, "Variante_40_.fasta"))
Var_50  <- readLines(paste0(base_path, "Variante_50_.fasta"))
Var_60  <- readLines(paste0(base_path, "Variante_60_.fasta"))
Var_70  <- readLines(paste0(base_path, "Variante_70_.fasta"))
Var_80  <- readLines(paste0(base_path, "Variante_80_.fasta"))
Var_90  <- readLines(paste0(base_path, "Variante_90_.fasta"))
Var_100 <- readLines(paste0(base_path, "Variante_100_.fasta"))

# Extração do domínio (posições 306-511) + collapse
extrair_dominio <- function(fasta) {
  seqs <- str_split(str_subset(fasta, "^EEKKVCQ"), "")
  dom  <- lapply(seqs, `[`, 306:511)
  sapply(dom, paste, collapse = "")
}

Var_10_dom_seq  <- extrair_dominio(Var_10)
Var_20_dom_seq  <- extrair_dominio(Var_20)
Var_30_dom_seq  <- extrair_dominio(Var_30)
Var_40_dom_seq  <- extrair_dominio(Var_40)
Var_50_dom_seq  <- extrair_dominio(Var_50)
Var_60_dom_seq  <- extrair_dominio(Var_60)
Var_70_dom_seq  <- extrair_dominio(Var_70)
Var_80_dom_seq  <- extrair_dominio(Var_80)
Var_90_dom_seq  <- extrair_dominio(Var_90)
Var_100_dom_seq <- extrair_dominio(Var_100)

substituir_seq <- function(fasta_original, seq_recortadas) {
  headers    <- str_subset(fasta_original, "^>")
  novo_fasta <- as.character(c(rbind(headers, seq_recortadas)))
  return(novo_fasta)
}

Var_10_novo  <- substituir_seq(Var_10,  Var_10_dom_seq)
Var_20_novo  <- substituir_seq(Var_20,  Var_20_dom_seq)
Var_30_novo  <- substituir_seq(Var_30,  Var_30_dom_seq)
Var_40_novo  <- substituir_seq(Var_40,  Var_40_dom_seq)
Var_50_novo  <- substituir_seq(Var_50,  Var_50_dom_seq)
Var_60_novo  <- substituir_seq(Var_60,  Var_60_dom_seq)
Var_70_novo  <- substituir_seq(Var_70,  Var_70_dom_seq)
Var_80_novo  <- substituir_seq(Var_80,  Var_80_dom_seq)
Var_90_novo  <- substituir_seq(Var_90,  Var_90_dom_seq)
Var_100_novo <- substituir_seq(Var_100, Var_100_dom_seq)

# Salvamento
writeLines(Var_10_novo,  paste0(base_path, "Variante_10_recortada.fasta"))
writeLines(Var_20_novo,  paste0(base_path, "Variante_20_recortada.fasta"))
writeLines(Var_30_novo,  paste0(base_path, "Variante_30_recortada.fasta"))
writeLines(Var_40_novo,  paste0(base_path, "Variante_40_recortada.fasta"))
writeLines(Var_50_novo,  paste0(base_path, "Variante_50_recortada.fasta"))
writeLines(Var_60_novo,  paste0(base_path, "Variante_60_recortada.fasta"))
writeLines(Var_70_novo,  paste0(base_path, "Variante_70_recortada.fasta"))
writeLines(Var_80_novo,  paste0(base_path, "Variante_80_recortada.fasta"))
writeLines(Var_90_novo,  paste0(base_path, "Variante_90_recortada.fasta"))
writeLines(Var_100_novo, paste0(base_path, "Variante_100_recortada.fasta"))
