#  Jesus Leon / LANGEBIO
# January 18th, 2021
library(tidyverse)
library(dplyr)
library(readr)

# Remember to move output blast files to Workdir
setwd("Workdir/")

# Loading and processing raw datasests ####
lxcxe_dataframes <- list.files(pattern = "*.fasta.csv", 
                               full.names = TRUE)
# adding columns and filtering low quality hits by bitscore/evalue
for (lxcxe_dataframe in lxcxe_dataframes) {
                  data <- read_csv(lxcxe_dataframe, 
                   col_names =  c("query_id",
                                  "prot_id",
                                  "evalue", 
                                  "bitscore")) %>%
    group_by(query_id) %>%
    filter(evalue == min(evalue))  %>%
    filter(bitscore == max(bitscore))  %>%
    write_csv(lxcxe_dataframe,
              ".csv")
}


annotation_dataframes <- list.files(pattern = "*Athaliana.csv",
                                    full.names = TRUE)
# BOECHERA annotation ####
boechera_tophits <- read_csv("Bstricta_LXCXE.fasta.csv")
boechera_arabidopsis <- read_csv("Bstricta_v_Athaliana.csv",
                               col_names = c("gene_id",
                                             "prot_id",
                                             "uniprotkb",
                                             "defline",
                                             "symbol")) %>%
  inner_join(boechera_tophits,
             boechera_arabidopsis,
             by = "prot_id") %>%
  write_csv(boechera_arabidopsis,
            col_names = TRUE,
            path = "Boechera-strictaLXCXE_BLASTp.csv")

# SOLANUM annotation ####
solanum_tophits <- read_csv("Slycopersicum_LXCXE.fasta.csv")
solanum_arabidopsis <- read_csv("Slycopersicum_v_Athaliana.csv",
                                 col_names = c("gene_id",
                                               "prot_id",
                                               "uniprotkb",
                                               "defline",
                                               "symbol")) %>%
  inner_join(solanum_tophits,
             solanum_arabidopsis,
             by = "prot_id") %>%
  write_csv(solanum_arabidopsis,
            col_names = TRUE,
            path = "Solanum-lycopersicumLXCXE_BLASTp.csv")

# ZEA annotation ####
zea_tophits <- read_csv("Zmays_LXCXE.fasta.csv")
zea_arabidopsis <- read_csv("Zmays_v_Athaliana.csv",
                                col_names = c("gene_id",
                                              "prot_id",
                                              "uniprotkb",
                                              "defline",
                                              "symbol")) %>%
  inner_join(zea_tophits,
             zea_arabidopsis,
             by = "prot_id") %>%
  write_csv(zea_arabidopsis,
            col_names = TRUE,
            path = "Zea-maysLXCXE_BLASTp.csv")


# ORYZA annotation ####
oryza_tophits <- read_csv("Osativa_LXCXE.fasta.csv")
oryza_arabidopsis <- read_csv("Osativa_v_Athaliana.csv",
                            col_names = c("gene_id",
                                          "prot_id",
                                          "uniprotkb",
                                          "defline",
                                          "symbol")) %>%
  inner_join(oryza_tophits,
             oryza_arabidopsis,
             by = "prot_id") %>%
  write_csv(oryza_arabidopsis,
            col_names = TRUE,
            path = "Oryza-sativaLXCXE_BLASTp.csv")

# SETARIA annotation ####
setaria_tophits <- read_csv("Sviridis_LXCXE.fasta.csv")
setaria_arabidopsis <- read_csv("Sviridis_v_Athaliana.csv",
                                 col_names = c("gene_id",
                                               "prot_id",
                                               "uniprotkb",
                                               "defline",
                                               "symbol")) %>%
  inner_join(setaria_tophits,
             setaria_arabidopsis,
             by = "prot_id") %>%
  write_csv(setaria_arabidopsis,
            col_names = TRUE,
            path = "Setaria-viridisLXCXE_BLASTp.csv")

# ZOOSTERA annotation ####
zoostera_tophits <- read_csv("Zmarina_LXCXE.fasta.csv")
zoostera_arabidopsis <- read_csv("Zmarina_v_Athaliana.csv",
                              col_names = c("gene_id",
                                            "prot_id",
                                            "uniprotkb",
                                            "defline",
                                            "symbol")) %>%
  inner_join(zoostera_tophits,
             zoostera_arabidopsis,
             by = "prot_id") %>%
  write_csv(zoostera_arabidopsis,
            col_names = TRUE,
            path = "Zoostera-marinaLXCXE_BLASTp.csv")

# AMBORELLA annotation ####
amborella_tophits <- read_csv("Atrichopoda_LXCXE.fasta.csv")
amborella_arabidopsis <- read_csv("Atrichopoda_v_Athaliana.csv",
                                 col_names = c("gene_id",
                                               "prot_id",
                                               "uniprotkb",
                                               "defline",
                                               "symbol")) %>%
  inner_join(amborella_tophits,
             amborella_arabidopsis,
             by = "prot_id") %>%
  write_csv(amborella_arabidopsis,
            col_names = TRUE,
            path = "Amborella-trichopodaLXCXE_BLASTp.csv")

# PICEA annotation ####
picea_tophits <- read_csv("Pabies_LXCXE.fasta.csv")
picea_arabidopsis <- read_csv("Pabies_v_Athaliana.csv",
                                  col_names = c("gene_id",
                                                "prot_id",
                                                "uniprotkb",
                                                "defline",
                                                "symbol")) %>%
  inner_join(picea_tophits,
             picea_arabidopsis,
             by = "prot_id") %>%
  write_csv(picea_arabidopsis,
            col_names = TRUE,
            path = "Picea-abiesLXCXE_BLASTp.csv")

# PINUS annotation ####
pinus_tophits <- read_csv("Ptaeda_LXCXE.fasta.csv")
pinus_arabidopsis <- read_csv("Ptaeda_v_Athaliana.csv",
                              col_names = c("gene_id",
                                            "prot_id",
                                            "uniprotkb",
                                            "defline",
                                            "symbol")) %>%
  inner_join(pinus_tophits,
             pinus_arabidopsis,
             by = "prot_id") %>%
  write_csv(pinus_arabidopsis,
            col_names = TRUE,
            path = "Pinus-taedaLXCXE_BLASTp.csv")

# GNETUM annotation ####
gnetum_tophits <- read_csv("Gmontanum_LXCXE.fasta.csv",
                           col_names =  c("query_id",
                                          "prot_id",
                                          "evalue", 
                                          "bitscore")) %>%
group_by(query_id) %>%
  filter(evalue == min(evalue))  %>%
  filter(bitscore == max(bitscore)) %>%
  write_csv(path = "GmontanumTopHits.csv")
gnetum_arabidopsis <- read_csv("Gmontanum_v_Athaliana.csv",
                              col_names = c("gene_id",
                                            "prot_id",
                                            "uniprotkb",
                                            "defline",
                                            "symbol")) %>%
  inner_join(gnetum_tophits,
             gnetum_arabidopsis,
             by = "prot_id") %>%
  write_csv(gnetum_arabidopsis,
            col_names = TRUE,
            path = "Gnetum-montanumLXCXE_BLASTp.csv")

# AZOLLA annotation ####
azolla_tophits <- read_csv("Afilliculoides_LXCXE.fasta.csv")
azolla_arabidopsis <- read_csv("Afilliculoides_v_Athaliana.csv",
                              col_names = c("gene_id",
                                            "prot_id",
                                            "uniprotkb",
                                            "defline",
                                            "symbol")) %>%
  inner_join(azolla_tophits,
             azolla_arabidopsis,
             by = "prot_id") %>%
  write_csv(azolla_arabidopsis,
            col_names = TRUE,
            path = "Azolla-filliculoidesLXCXE_BLASTp.csv")

# SALVINIA annotation ####
salvinia_tophits <- read_csv("Scucullata_LXCXE.fasta.csv")
salvinia_arabidopsis <- read_csv("Scucullata_v_Athaliana.csv",
                               col_names = c("gene_id",
                                             "prot_id",
                                             "uniprotkb",
                                             "defline",
                                             "symbol")) %>%
  inner_join(salvinia_tophits,
             salvinia_arabidopsis,
             by = "prot_id") %>%
  write_csv(salvinia_arabidopsis,
            col_names = TRUE,
            path = "Salvinia-cucullataLXCXE_BLASTp.csv")

# CERATOPTERIS annotation ####
ceratopteris_tophits <- read_csv("Crichardii_LXCXE.fasta.csv")
ceratopteris_arabidopsis <- read_csv("Crichardii_v_Athaliana.csv",
                                     col_names = c("gene_id",
                                                   "prot_id",
                                                   "uniprotkb",
                                                   "defline",
                                                   "symbol")) %>%
  inner_join(ceratopteris_tophits,
             ceratopteris_arabidopsis,
             by = "prot_id") %>%
  write_csv(ceratopteris_arabidopsis_arabidopsis,
            col_names = TRUE,
            path = "Ceratopteris-richardiiLXCXE_BLASTp.csv")

# SELAGINELLA annotation ####
selaginella_tophits <- read_csv("Smoellendorffii_LXCXE.fasta.csv")
seglaginella_arabidopsis <- read_csv("Smoellendorffii_v_Athaliana.csv",
                                 col_names = c("gene_id",
                                               "prot_id",
                                               "uniprotkb",
                                               "defline",
                                               "symbol")) %>%
  inner_join(selaginella_tophits,
             seglaginella_arabidopsis,
             by = "prot_id") %>%
  write_csv(seglaginella_arabidopsis,
            col_names = TRUE,
            path = "Selaginella-moellendorffiiLXCXE_BLASTp.csv")

# PHYSCOMITRELLA annotation ####
physcomitrella_tophits <- read_csv("Ppatens_LXCXE.fasta.csv")
physcomitrella_arabidopsis <- read_csv("Ppatens_v_Athaliana.csv",
                                     col_names = c("gene_id",
                                                   "prot_id",
                                                   "uniprotkb",
                                                   "defline",
                                                   "symbol")) %>%
  inner_join(physcomitrella_tophits,
             physcomitrella_arabidopsis,
             by = "prot_id") %>%
  write_csv(physcomitrella_arabidopsis,
            col_names = TRUE,
            path = "Physcomitrella-patensLXCXE_BLASTp.csv")

# SPHAGNUM annotation ####
sphagnum_tophits <- read_csv("Sfallax_LXCXE.fasta.csv")
sphagnum_arabidopsis <- read_csv("Sfallax_v_Athaliana.csv",
                                       col_names = c("gene_id",
                                                     "prot_id",
                                                     "uniprotkb",
                                                     "defline",
                                                     "symbol")) %>%
  inner_join(sphagnum_tophits,
             sphagnum_arabidopsis,
             by = "prot_id") %>%
  write_csv(sphagnum_arabidopsis,
            col_names = TRUE,
            path = "Sphagnum-fallaxLXCXE_BLASTp.csv")

# MARCHANTIA annotation ####
marchantia_tophits <- read_csv("Mpolymorpha_LXCXE.fasta.csv")
marchantia_arabidopsis <- read_csv("Mpolymorpha_v_Athaliana.csv",
                                 col_names = c("gene_id",
                                               "prot_id",
                                               "uniprotkb",
                                               "defline",
                                               "symbol")) %>%
  inner_join(marchantia_tophits,
             marchantia_arabidopsis,
             by = "prot_id") %>%
  write_csv(marchantia_arabidopsis,
            col_names = TRUE,
            path = "Marchantia-polymorphaLXCXE_BLASTp.csv")

# ANTHOCEROS annotation ####
anthoceros_tophits <- read_csv("Aagrestis_LXCXE.fasta.csv")
anthoceros_arabidopsis <- read_csv("Aagrestis_v_Athaliana.csv",
                                   col_names = c("gene_id",
                                                 "prot_id",
                                                 "uniprotkb",
                                                 "defline",
                                                 "symbol")) %>%
  inner_join(anthoceros_tophits,
             anthoceros_arabidopsis,
             by = "prot_id") %>%
  write_csv(anthoceros_arabidopsis,
            col_names = TRUE,
            path = "Anthoceros-agrestisLXCXE_BLASTp.csv")
# PENIUM annotation ####
penium_tophits <- read_csv("Pmargaritaceum_LXCXE.fasta.csv")
penium_arabidopsis <- read_csv("Pmagaritaceum_v_Athaliana.csv",
                                 col_names = c("gene_id",
                                               "prot_id",
                                               "uniprotkb",
                                               "defline",
                                               "symbol")) %>%
  inner_join(penium_tophits,
             penium_arabidopsis,
             by = "prot_id") %>%
write_csv(penium_arabidopsis,
          col_names = TRUE,
          path = "Penium-margaritaceumLXCXE_BLASTp.csv")

# Spirogloea annotation ####
spirogloea_tophits <- read_csv("Smuscicola_LXCXE.fasta.csv")
spirogloea_arabidopsis <- read_csv("Smuscicola_v_Athaliana.csv",
                               col_names = c("gene_id",
                                             "prot_id",
                                             "uniprotkb",
                                             "defline",
                                             "symbol")) %>%
  inner_join(spirogloea_tophits,
             spirogloea_arabidopsis,
             by = "prot_id") %>%
  write_csv(spirogloea_arabidopsis,
            col_names = TRUE,
            path = "Spirogloea-muscicolaLXCXE_BLASTp.csv")

# Mesotaenium annotation ####
mesotaenium_tophits <- read_csv("Mendlicheranium_LXCXE.fasta.csv")
mesotaenium_arabidopsis <- read_csv("Mendlicheranium_v_Athaliana.csv",
                                   col_names = c("gene_id",
                                                 "prot_id",
                                                 "uniprotkb",
                                                 "defline",
                                                 "symbol")) %>%
  inner_join(mesotaenium_tophits,
             mesotaenium_arabidopsis,
             by = "prot_id") %>%
  write_csv(mesotaenium_arabidopsis,
            col_names = TRUE,
            path = "Mesotaenium-endlicheraniumLXCXE_BLASTp.csv")

# Chara annotation ####
chara_tophits <- read_csv("Cbraunii_LXCXE.fasta.csv")
chara_arabidopsis <- read_csv("Cbraunii_v_Athaliana.csv",
                                    col_names = c("gene_id",
                                                  "prot_id",
                                                  "uniprotkb",
                                                  "defline",
                                                  "symbol")) %>%
  inner_join(chara_tophits,
             chara_arabidopsis,
             by = "prot_id") %>%
  write_csv(chara_arabidopsis,
            col_names = TRUE,
            path = "Chara-brauniiLXCXE_BLASTp.csv")

# Klebsormidium annotation ####
klebsormidium_tophits <- read_csv("Knitens_LXCXE.fasta.csv")
klebsormidium_arabidopsis <- read_csv("Knitens_v_Athaliana.csv",
                              col_names = c("gene_id",
                                            "prot_id",
                                            "uniprotkb",
                                            "defline",
                                            "symbol")) %>%
  inner_join(klebsormidium_tophits,
             klebsormidium_arabidopsis,
             by = "prot_id") %>%
  write_csv(klebsormidium_arabidopsis,
            col_names = TRUE,
            path = "Klebsormidium-nitensLXCXE_BLASTp.csv")

# Chlorokybus annotation ####
chlorokybus_tophits <- read_csv("Catmophyticus_LXCXE.fasta.csv")
chlorokybus_arabidopsis <- read_csv("Catmophyticus_v_Athaliana.csv",
                                      col_names = c("gene_id",
                                                    "prot_id",
                                                    "uniprotkb",
                                                    "defline",
                                                    "symbol")) %>%
  inner_join(chlorokybus_tophits,
             chlorokybus_arabidopsis,
             by = "prot_id") %>%
  write_csv(chlorokybus_arabidopsis,
            col_names = TRUE,
            path = "Chlorokybus-atmophyticusLXCXE_BLASTp.csv")

# Mesostigma annotation ####
mesostigma_tophits <- read_csv("MvirideLIANG_LXCXE.fasta.csv")
mesostigma_arabidopsis <- read_csv("Mviride_v_Athaliana.csv",
                                    col_names = c("gene_id",
                                                  "prot_id",
                                                  "uniprotkb",
                                                  "defline",
                                                  "symbol")) %>%
  inner_join(mesostigma_tophits,
             mesostigma_arabidopsis,
             by = "prot_id") %>%
  write_csv(mesostigma_arabidopsis,
            col_names = TRUE,
            path = "Mesostigma-virideLXCXE_BLASTp.csv")

# Volvox annotation ####
volvox_tophits <- read_csv("Vcarteri_LXCXE.fasta.csv")
volvox_arabidopsis <- read_csv("Vcarteri_v_Athaliana.csv",
                                   col_names = c("gene_id",
                                                 "prot_id",
                                                 "uniprotkb",
                                                 "defline",
                                                 "symbol")) %>%
  inner_join(volvox_tophits,
             volvox_arabidopsis,
             by = "prot_id") %>%
  write_csv(volvox_arabidopsis,
            col_names = TRUE,
            path = "Volvox-carteriLXCXE_BLASTp.csv")

# CHLAMY annotation ####
chlamydomonas_tophits <- read_csv("Creinhardtii_LXCXE.fasta.csv")
chlamydomonas_arabidopsis <- read_csv("Creinhardtii_v_Athaliana.csv",
                               col_names = c("gene_id",
                                             "prot_id",
                                             "uniprotkb",
                                             "defline",
                                             "symbol")) %>%
  inner_join(chlamydomonas_tophits,
             chlamydomonas_arabidopsis,
             by = "prot_id") %>%
  write_csv(volvox_arabidopsis,
            col_names = TRUE,
            path = "Chlamydomonas-reinhardtiiLXCXE_BLASTp.csv")


# OSTREOCOCCUS annotation ####
ostreococcus_tophits <- read_csv("Olucimarinus_LXCXE.fasta.csv")
ostreococcus_arabidopsis <- read_csv("Olucimarinus_v_Athaliana.csv",
                                     col_names = c("gene_id",
                                                   "prot_id",
                                                   "uniprotkb",
                                                   "defline",
                                                   "symbol")) %>%
  inner_join(ostreococcus_tophits,
             ostreococcus_arabidopsis,
             by = "prot_id") %>%
  write_csv(ostreococcus_arabidopsis,
            col_names = TRUE,
            path = "Ostreococcus-lucimarinusLXCXE_BLASTp.csv")

