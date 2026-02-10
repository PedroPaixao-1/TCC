source("Code/Tidyverse_function.R")

Arquivo_pdb <- read.pdb("/home/pedro.paixao/Code/egfr_padronizado2.pdb")
EGFR_Atom <- Arquivo_pdb$atom
Ca_a <- atom.select(Arquivo_pdb, "calpha", chain = "A")

Arquivo_pdb_ca_a <- trim.pdb(Arquivo_pdb, inds = Ca_a)

Ca_a_atom <- Arquivo_pdb_ca_a$atom

EGFR_atom_padronizado <- Padronizar_pdb(Arquivo_pdb, TRUE)

EGFR_ca <- EGFR_atom_padronizado %>%
  filter(elety == "CA")



Seq_EGFR <- Sequenciar_interesse(Arquivo_pdb,1:613)

Ct_EGFR <- Calcular_Contato(EGFR_ca, 8)

Teste <- Rodar_ESM2(Arquivo_pdb,1)
