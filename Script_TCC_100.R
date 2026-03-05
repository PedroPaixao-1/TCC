#scp -r Code pedro.paixao@157.86.114.54:/home/pedro.paixao
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

Variantes <- list()
Multi_fasta <- list()

j <- 1

while (j <= 200) {
  
    
    Posicao <- sample(Contatos_EGFR, size = ceiling(100/100*length(Contatos_EGFR)))
    unlist(Posicao)
    
    Comando <- paste("/home/conda/condabin/conda run -n esm2","python /home/pedro.paixao/Code/generate_sequence_esm2.py", "--sequence", Sequencia, "--position", Posicao,"--temperature", Temperatura)
    #Isso aqui vai dar merda por causa do diretório na rtx
    system(Comando, wait = TRUE)
    
    Saida <- readLines("completed_sequence.txt")
    #Variantes[[i]] <- Saida
    if (Saida %in% Variantes){
      next
    }
    Header <- paste(">EGFR|CHAIN A|Variante",j)
    Multi_fasta[[j]] <- paste(Header,Saida, sep = "\n")
    writeLines(unlist(Multi_fasta), "Variantes.fasta")
    j = j+1
  }





