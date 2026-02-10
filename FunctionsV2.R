library("reticulate")
library("bio3d")
library("dplyr")
library("readr")
library("seqinr")

#Funções de fundação

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

Indices_ct <- function(pdb_atom,Matriz_ct,Distancia_contato){
  
  chains <- tapply(pdb_atom$chain, pdb_atom$resno, unique)
  chains <- as.character(chains)
  
  contatos <- Matriz_ct <= Distancia_contato
  
  i <- row(contatos)[contatos]
  j <- col(contatos)[contatos]
  
  vetor_log <- i < j
    
  validos <-
    (chains[i] == "A" & chains[j] == "B") |
    (chains[i] == "B" & chains[j] == "A")
  
  cbind(i[validos & vetor_log], j[validos & vetor_log])
  
  
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
  
  Matriz_ct <- Calcular_Distância(Interesse = Arquivo_pdb)
  
  Indices <- Indices_ct(Arquivo_pdb,Matriz_ct,Distancia_contato)

  
return(Indices)
  
}


#Funções de 2° Nível







#--------------------------------TESTE-----------------------------------#




