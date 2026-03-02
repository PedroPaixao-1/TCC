library("reticulate")
library("bio3d")
library("dplyr")
library("readr")
library("seqinr")

#Funções de fundação

Renomear_cadeias <- function(Arquivo_pdb){
  pdb_atom <- Arquivo_pdb$atom
  
  #Salva as cadeias como um vetor e renomeia elas usando LETTERS
  
  Cadeias <- unique(pdb_atom$chain)
  pdb_atom$chain <- LETTERS[match(pdb_atom$chain, Cadeias)] #match retorna um vetor das POSIÇÕES então independentemente de qual seja o nome da cadeia vai virar A,B,C... nessa ordem e tu sabe o pq ent n se prolonga!
  
  return(pdb_atom)
  
}

Renumerar_resid <- function(pdb_atom){

 # Cria pares de coordenadas (cadeia,resno)
  
  Id_resd <- paste(pdb_atom$chain,pdb_atom$resno)
  
  # Muda os resnos seguindo a ordem dos pares, ou seja, vai resolver TODOS os da cadeia A antes e depois a da B.
  
  Resno_nv <- match(Id_resd, unique(Id_resd)) # (A;Na) e (B;(A;Na)+Nb)
  pdb_atom$resno <- Resno_nv
  
  return(pdb_atom)
}

Renumerar_bfactor <- function(pdb_atom,Limpar_b_factor){
  #Muda o b facotr de todos os átomos para 0
  if (Limpar_b_factor == TRUE){pdb_atom <- pdb_atom %>%
    mutate(b = 0)} 
  
  return(pdb_atom)
}

Calcular_Distância <- function(Interesse){
 
  x <- Interesse[, "x"]
  y <- Interesse[, "y"]
  z <- Interesse[, "z"]
  
  n <- length(x)
  Matriz_distancias <- matrix(0, nrow = n, ncol = n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      Distância <- sqrt((x[i] - x[j])**2 + (y[i] - y[j])**2 + (z[i] - z[j])**2)
      Matriz_distancias[i, j] <- Distância
    }
  }
  
  return(Matriz_distancias)
}

Indices_ct <- function(Arquivo_pdb,Matriz_ct,Distancia_contato){
  pdb_atom <- Arquivo_pdb$atom
  
  chains <- tapply(pdb_atom$chain, pdb_atom$resno, unique)
  chains <- as.character(chains)
  
  contatos <- Matriz_ct <= Distancia_contato
  
  i <- row(contatos)[contatos]
  j <- col(contatos)[contatos]
  
  vetor_log <- i < j
    
  validos <-
    (chains[i] == "A" & chains[j] == "B") |
    (chains[i] == "B" & chains[j] == "A")
  
  Ids <- cbind(i[validos & vetor_log], j[validos & vetor_log])
  
  return(Ids)
}

Sequenciar_interesse <- function(pdb){
  seq <- pdbseq(pdb)
  paste(seq, collapse = "")
}

#Funções de 1° Nível

Padronizar_pdb <- function(Arquivo_pdb,Limpar_b_factor){
  
  pdb_atom <- Renomear_cadeias(Arquivo_pdb)
  
  pdb_atom <- Renumerar_resid(pdb_atom)
  
  if (Limpar_b_factor == TRUE){
    pdb_atom <- Renumerar_bfactor(pdb_atom,Limpar_b_factor)
  }
  
  return(pdb_atom)
}

Calcular_Contato <- function(Arquivo_pdb,Distancia_contato){
  Atom <- Arquivo_pdb$atom
  
  Matriz_ct <- Calcular_Distância(Interesse = Atom)
  
  Indices <- Indices_ct(Arquivo_pdb,Matriz_ct,Distancia_contato)
  
  Contatos_A <- paste(unique(cbind(Indices)[,1]))

  return(Contatos_A)
  
}


#Funções de 2° Nível



#--------------------------------Bug yard---------------------------------#

Rodar_ESM2 <- function(Arquivo_pdb,Temperatura) {
  Porcentagens <- c(5,100)
  
  Sequencia <- Sequenciar_interesse(Arquivo_pdb)
  
  Contatos <- Calcular_Contato(Arquivo_pdb, Distancia_contato = 8)
  for (i in Amostras) {
    for (i in 1:length(Porcentagens)) {
      Posicao <- Randomizar_posicao(Contatos,Porcentagens)  
    }
  }
  
  
  
  
  Comando <- paste("/home/pedro.paixao/anaconda3/condabin/conda run -n esm2-env","python /home/pedro.paixao/Code/generate_sequence_esm2.py", "--sequence", Sequencia, "--position", Posicao,"--temperature", Temperatura)
  
  system(Comando, wait = TRUE)
  
  saida <- readLines("completed_sequence.txt")
  
}

Randomizar_posicao <- function(Contatos_unicos,Porcentagem){
  Posicao <- sample(Contatos_unicos, size = ceiling(Porcentagem/100*length(Contatos_unicos)))
  return(unlist(Posicao))
}

#--------------------------------TESTE-----------------------------------#


