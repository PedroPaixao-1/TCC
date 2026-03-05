#scp -r Code pedro.paixao@157.86.114.54:/home/pedro.paixao
source("/home/pedro.paixao/Code//FunctionsV2.R")

Aminoacidos <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 
                 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')

Arquivo_pdb <- read.pdb("/home/pedro.paixao/Code/egfr_padronizado2.pdb")
Arquivo_pdb_ca <- atom.select(Arquivo_pdb, "calpha")
Arquivo_pdb_ca <- trim.pdb(Arquivo_pdb, inds = Arquivo_pdb_ca)

Contatos_EGFR <- Calcular_Contato(Arquivo_pdb_ca, 8)

Sequencia <- pdbseq(Arquivo_pdb)
Sequencia <- Sequencia[1:613]
Sequencia <- paste(Sequencia, collapse = "")

Temperatura <- 1.5

Variantes <- list()
Multi_fasta <- list()

Contador <- 1

for (i in 1:length(Contatos_EGFR)) {
  for (j in 1:length(Aminoacidos)) {
    Variantes <- Sequencia
    substring(Variantes,Contatos_EGFR[i],Contatos_EGFR[i]) <- Aminoacidos[j]
    
    Header <- paste(">EGFR|CHAIN A|Variante",Contador)
    Multi_fasta[[Contador]] <- paste(Header,Variantes, sep = "\n")
    Contador = Contador + 1
  }
}

writeLines(unlist(Multi_fasta), "Variantes_5_porcento.fasta")

