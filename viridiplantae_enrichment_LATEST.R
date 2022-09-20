# This script is for performing functional profiling
# of LxCxE motif containing  gene lists across the Viridiplantae phylogeny
# Jesus A. Leon Ruiz
# alberto.leon@cinvestav.mx
# Summer 2020

# Note: to detach org.At.tair.db package to use select more easily, run
# detach("package:org.At.tair.db", unload = T)

# Loading packages
library(clusterProfiler)
library(tidyverse)
library(dplyr)
library(pheatmap)
library(viridis)

#### REMEMBER to move *_BLASTp.csv files to Workdir ####

#Setting working directory
setwd("Workdir/")
# Importing all data ####
# arabidopsis 
at <-read_csv("Arabidopsis-thalianaLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
at_ids <- at$gene_id

# Boechera stricta
bs <- read_csv("Boechera-strictaLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
bs_ids <- bs$gene_id

# Solanum lycopercisum
sl <- read_csv("Solanum-lycopersicumLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
sl_ids <- sl$gene_id

# Zea mays
zm <- read_csv("Zea-maysLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
zm_ids <- zm$gene_id

# Oryza sativa
os <- read_csv("Oryza-sativaLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
os_ids <- os$gene_id

# Setaria viridis
sv <- read_csv("Setaria-viridisLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
sv_ids <- sv$gene_id

# Amborella trichopoda
atr <- read_csv("Amborella-trichopodaLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
atr_ids <- atr$gene_id

# Picea abies
pa <- read_csv("Picea-abiesLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
pa_ids <- pa$gene_id

# Pinus taeda 
pt <- read_csv("Pinus-taedaLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
pt_ids <- pt$gene_id

# Gnetum montanum
gm <- read_csv("Gnetum-montanumLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
gm_ids <- gm$gene_id

# Azolla filliculoides
af <- read_csv("Azolla-filliculoidesLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
af_ids <- af$gene_id

# Salvinia cucullata
sc <- read_csv("Salvinia-cucullataLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
sc_ids <- sc$gene_id

# Ceratopteris richardii
cr <- read_csv("Ceratopteris-richardiiLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
cr_ids <- cr$gene_id

# Selaginella moellendorffii
sm <- read_csv("Selaginella-moellendorffiiLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
sm_ids <- sm$gene_id

# Isoetes taiwanensis
it <- read_csv("Isoetes-taiwanensisLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
it_ids <- it$gene_id

# Physcomitrium patens
pp <- read_csv("Physcomitrella-patensLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
pp_ids <- pp$gene_id

# Sphagnum fallax
sf <- read_csv("Sphagnum-fallaxLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
sf_ids <- sf$gene_id

# Marchantia polymorpha
mp <- read_csv("Marchantia-polymorphaLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
mp_ids <- mp$gene_id

# Anthoceros agrestis
aa <- read_csv("Anthoceros-agrestisLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
aa_ids <- aa$gene_id

# Spirogloea muscicola
smu <- read_csv("Spirogloea-muscicolaLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
smu_ids <- smu$gene_id

# Mesotaenium endlicheranium
me <- read_csv("Mesotaenium-endlicheraniumLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
me_ids <- me$gene_id

# Penium margaritaceum
pma <- read_csv("Penium-margaritaceumLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
pma_ids <- pma$gene_id

# Chara braunii
cb <- read_csv("Chara-brauniiLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
cb_ids <- cb$gene_id

# Klebsormidium nitens
kn <- read_csv("Klebsormidium-nitensLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
kn_ids <- kn$gene_id

# Chlorokybus atmophyticus
ca <- read_csv("Chlorokybus-atmophyticusLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
ca_ids <- ca$gene_id

# Mesostigma viride
mv <- read_csv("Mesostigma-virideLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
mv_ids <- mv$gene_id

# Chlamydomonas reinhardtii
cre <- read_csv("Chlamydomonas-reinhardtiiLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
cre_ids <- cr$gene_id

# Volvox carteri
vc <- read_csv("Volvox-carteriLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
vc_ids <- vc$gene_id

# Ostreococcus lucimarinus
ol <- read_csv("Ostreococcus-lucimarinusLXCXE_BLASTp.csv") %>%
  dplyr::select(gene_id)
ol_ids <- ol$gene_id

# Retrieving all TAIR keys from org.At.tair.db to create universe of genes
library(org.At.tair.db)
at_universe<- AnnotationDbi::keys(org.At.tair.db,
                                  keytype = 'TAIR')
detach("package:dplyr", unload = TRUE)

# Nuclear enrichments ####
# Arabidopsis thaliana - Nuclear ####
ath <- AnnotationDbi::select(org.At.tair.db,
                             keys=at_ids,
                             columns = c("ONTOLOGY",
                                         "GO",
                                         "SYMBOL"),
                             keytype = "TAIR")
ath_nuclear <- ath %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634")) %>%  # CC GO for nucleus
write_csv("AthalianaLXCXE_nuclear.csv")
ath_nuclear$TAIR 

# Boechera stricta - Nuclear ####
bst <- AnnotationDbi::select(org.At.tair.db,
                             keys=bs_ids,
                             columns = c("ONTOLOGY",
                                         "GO",
                                         "SYMBOL"),
                             keytype = "TAIR")
bst_nuclear <- bst %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634")) %>%  # CC GO for nucleus
  write_csv("BstrictaLXCXE_nuclear.csv") 
bst_nuclear$TAIR

# Solanum lycopersicum - Nuclear ####
sly <- AnnotationDbi::select(org.At.tair.db,
                             keys=sl_ids,
                             columns = c("ONTOLOGY",
                                         "GO",
                                         "SYMBOL"),
                             keytype = "TAIR")
sly_nuclear <- sly %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634")) %>%  # CC GO for nucleus
  write_csv("SlycopersicumLXCXE_nuclear.csv")  # CC GO for nucleus
sly_nuclear$TAIR

# Zea mays - Nuclear ####
zma <- AnnotationDbi::select(org.At.tair.db,
                             keys=zm_ids,
                             columns = c("ONTOLOGY",
                                         "GO",
                                         "SYMBOL"),
                             keytype = "TAIR")
zma_nuclear <- zma %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634")) %>%  # CC GO for nucleus
  write_csv("ZmaysLXCXE_nuclear.csv") 
zma_nuclear$TAIR

# Oryza sativa - Nuclear ####
osa <- AnnotationDbi::select(org.At.tair.db,
                             keys=os_ids,
                             columns = c("ONTOLOGY",
                                         "GO",
                                         "SYMBOL"),
                             keytype = "TAIR")
osa_nuclear <- osa %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634")) %>%  # CC GO for nucleus
  write_csv("OsativaLXCXE_nuclear.csv") 
osa_nuclear$TAIR

# Setaria viridis - Nuclear ####
svi <- AnnotationDbi::select(org.At.tair.db,
                             keys=sv_ids,
                             columns = c("ONTOLOGY",
                                         "GO",
                                         "SYMBOL"),
                             keytype = "TAIR")
svi_nuclear <- svi %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634")) %>%  # CC GO for nucleus
  write_csv("SviridisLXCXE_nuclear.csv")  
svi_nuclear$TAIR


# Amborella trichopoda - Nuclear ####
atr <- AnnotationDbi::select(org.At.tair.db,
                             keys=atr_ids,
                             columns = c("ONTOLOGY",
                                         "GO",
                                         "SYMBOL"),
                             keytype = "TAIR")
atr_nuclear <- atr %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634")) %>%  # CC GO for nucleus
  write_csv("AtrichopodaLXCXE_nuclear.csv")
atr_nuclear$TAIR

# Angiosperms compareCluster ####
angiosperms_list <- list("Ath" = ath_nuclear$TAIR,
                         "Bst" = bst_nuclear$TAIR,
                         "Sly" = sly_nuclear$TAIR,
                         "Zma" = zma_nuclear$TAIR,
                         "Osa" = osa_nuclear$TAIR,
                         "Svi" = svi_nuclear$TAIR,
                         "Atr" = atr_nuclear$TAIR) 

angio_compare <- compareCluster(angiosperms_list, fun="enrichGO",
                                pvalueCutoff=0.05,
                                pAdjustMethod="BH",
                                keyType = "TAIR",
                                OrgDb=org.At.tair.db,
                                ont="BP")
dotplot(angio_compare, showCategory = 20) + scale_color_viridis()

# Converting to tibble
angio_tibble <- as_tibble(angio_compare) %>%
  dplyr::select(Cluster,
                Description,
                p.adjust, geneID)

# Extracting each cluster and then combining them with full_angio_tibble
# A_thaliana
arabidopsis_tibble <- angio_tibble %>%
  dplyr::filter(Cluster == "Ath") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Ath = p.adjust)
# B_stricta
boechera_tibble <- angio_tibble %>%
  dplyr::filter(Cluster == "Bst") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Bst = p.adjust)
# S_lycopersicum
solanum_tibble <- angio_tibble %>%
  dplyr::filter(Cluster == "Sly") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Sly = p.adjust)
# Z_mays
zea_tibble <- angio_tibble %>%
  dplyr::filter(Cluster == "Zma") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Zma = p.adjust)
# O_sativa
oryza_tibble <- angio_tibble %>%
  dplyr::filter(Cluster == "Osa") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Osa = p.adjust)
# S_viridis
setaria_tibble <- angio_tibble %>%
  dplyr::filter(Cluster == "Svi") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Svi = p.adjust)
# A_trichopoda
amborella_tibble <- angio_tibble %>%
  dplyr::filter(Cluster == "Atr") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Atr = p.adjust)

# Writing files
write_csv(arabidopsis_tibble, path = "Athaliana-GO_LXCXE.csv")
write_csv(boechera_tibble, path = "Bstricta-GO_LXCXE.csv")
write_csv(solanum_tibble, path = "Slycopersicum-GO_LXCXE.csv")
write_csv(zea_tibble, path = "Zmays-GO_LXCXE.csv")
write_csv(oryza_tibble, path = "Osativa-GO_LXCXE.csv")
write_csv(setaria_tibble, path = "Sviridis-GO_LXCXE.csv")
write_csv(amborella_tibble, path = "Atrichopoda-GO_LXCXE.csv")




# Performing joins to get final tibble
angiosperms <- dplyr::full_join(arabidopsis_tibble,
                                boechera_tibble,
                                by = "Description") %>%
  dplyr::full_join(solanum_tibble,
                   by = "Description") %>%
  dplyr::full_join(zea_tibble,
                   by = "Description") %>%
  dplyr::full_join(oryza_tibble,
                   by = "Description") %>%
  dplyr::full_join(setaria_tibble,
                   by = "Description") %>%
  dplyr::full_join(amborella_tibble,
                   by = "Description") %>%
  dplyr::select(Description, Ath, Bst, Sly, Zma, Osa, Svi, Atr) %>%
  dplyr::slice(1:40)

# Angiosperms matrix/Heatmap #####
write_csv(angiosperms, "angiosperms_mat.csv", col_names = T)
angiosperms_mat <- as.matrix(read.csv("angiosperms_mat.csv",
                                      sep = ",",
                                      row.names = "Description"))
angiosperms_matrix <- -log10(angiosperms_mat)
angiosperm_heatmap <- pheatmap(
  angiosperms_matrix,
  color = viridis(500, option = "D", begin = 0, end = 0.8),
  bias = 1, 
  border_color = F,
  cellwidth = 7,
  cellheight = 4,
  fontsize = 5,
  na_col = "light grey",
  clustering_distance_rows = "euclidean", 
  cluster_rows = F,
  cluster_cols = F,
  show_colnames = T,
  show_rownames = TRUE,
  main = "Angiosperms_LxCxE")
angiosperm_heatmap #### Angiosperm heatmap
pdf("angiosperm_heatmap.pdf")
print(angiosperm_heatmap) # Plot 1 --> in the first page of PDF
dev.off()

# Gymnosperms 
# Picea abies - Nuclear ####
pab <- AnnotationDbi::select(org.At.tair.db,
              keys=pa_ids,
              columns = c("ONTOLOGY",
                          "GO",
                          "SYMBOL"),
              keytype = "TAIR")
pab_nuclear <- pab %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634")) %>%  # CC GO for nucleus
  write_csv("PabiesLXCXE_nuclear.csv")
pab_nuclear$TAIR

# Pinus taeda - Nuclear ####
pta <- AnnotationDbi::select(org.At.tair.db,
              keys=pt_ids,
              columns = c("ONTOLOGY",
                          "GO",
                          "SYMBOL"),
              keytype = "TAIR")
pta_nuclear <- pta %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634"))  %>%  # CC GO for nucleus
  write_csv("PtaedaLXCXE_nuclear.csv")
pta_nuclear$TAIR

# Gnetum montanum - Nuclear ####
gmo <- AnnotationDbi::select(org.At.tair.db,
              keys=gm_ids,
              columns = c("ONTOLOGY",
                          "GO",
                          "SYMBOL"),
              keytype = "TAIR")
gmo_nuclear <- gmo %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634")) %>%  # CC GO for nucleus
  write_csv("GmontanumLXCXE_nuclear.csv")
gmo_nuclear$TAIR

# Gymnosperms compareCluster ####
gymnosperms_list <- list("Pab" = pab_nuclear$TAIR,
                         "Pta" = pta_nuclear$TAIR,
                         "Gmo" = gmo_nuclear$TAIR)

gymno_compare <- compareCluster(gymnosperms_list, 
                                fun="enrichGO",
                                pvalueCutoff=0.05,
                                pAdjustMethod="BH",
                                keyType = "TAIR",
                                OrgDb=org.At.tair.db,
                                ont="BP")
dotplot(gymno_compare, showCategory = 35) + scale_color_viridis()

# Converting to tibble!
gymno_tibble <- as_tibble(gymno_compare) %>%
  dplyr::select(Cluster,
                Description,
                p.adjust,
                geneID)
# Extracting each cluster and then combining them with full_angio_tibble
# P_abies
picea_tibble <- gymno_tibble %>%
  dplyr::filter(Cluster == "Pab") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Pab = p.adjust)
# P_taeda
pinus_tibble <- gymno_tibble %>%
  dplyr::filter(Cluster == "Pta") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Pta = p.adjust)
# G_montanum
gnetum_tibble <- gymno_tibble %>%
  dplyr::filter(Cluster == "Gmo") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Gmo = p.adjust)

# Writing tibbles
write_csv(picea_tibble, path = "Picea_abiesGO.csv")
write_csv(pinus_tibble, path = "Pinus_taedaGO.csv")
write_csv(gnetum_tibble, path = "Gnetum_montanumGO.csv")

# Performing joints to get final tibble
gymnosperms <- dplyr::full_join(picea_tibble,
                                pinus_tibble,
                                by = "Description") %>%
  dplyr::full_join(gnetum_tibble,
                   by = "Description") %>%
  dplyr::select(Description, Pab, Pta, Gmo) %>%
  dplyr::slice(1:40)

# Gymnosperms matrix/Heatmap #####
write_csv(gymnosperms, "gymnosperms_mat.csv", col_names = T)
gymnosperms_mat <- as.matrix(read.csv("gymnosperms_mat.csv",
                                      sep = ",",
                                      row.names = "Description"))
gymnosperms_matrix <- -log10(gymnosperms_mat)
gymnosperm_heatmap <- pheatmap(
  gymnosperms_matrix,
  color = viridis(500, option = "D", begin = 0, end = 0.8),
  bias = 1, 
  border_color = F,
  cellwidth = 7,
  cellheight = 4,
  fontsize = 5,
  na_col = "light grey",
  clustering_distance_rows = "euclidean", 
  cluster_rows = F,
  cluster_cols = F,
  show_colnames = T,
  show_rownames = TRUE,
  main = "Gymnosperms")
gymnosperm_heatmap
pdf("gymnosperm_heatmap.pdf")
print(gymnosperm_heatmap) # Plot 1 --> in the first page of PDF
dev.off()

# Ferns ####
# Azolla filliculoides - Nuclear ####
af <- AnnotationDbi::select(org.At.tair.db,
             keys=af_ids,
             columns = c("ONTOLOGY",
                         "GO",
                         "SYMBOL"),
             keytype = "TAIR")
af_nuclear <- af %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634")) %>%  # CC GO for nucleus
  write_csv("AfilliculoidesLXCXE_nuclear.csv")
af_nuclear$TAIR

# Salvinia cucullata Nuclear ####
sc <- AnnotationDbi::select(org.At.tair.db,
             keys=sc_ids,
             columns = c("ONTOLOGY",
                         "GO",
                         "SYMBOL"),
             keytype = "TAIR")
sc_nuclear <- sc %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634")) %>%  # CC GO for nucleus
  write_csv("ScucullataLXCXE_nuclear.csv")
sc_nuclear$TAIR

# Ceratopteris richardii Nuclear ####
cr <- AnnotationDbi::select(org.At.tair.db,
             keys=cr_ids,
             columns = c("ONTOLOGY",
                         "GO",
                         "SYMBOL"),
             keytype = "TAIR")
cr_nuclear <- cr %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634")) %>%  # CC GO for nucleus
  write_csv("CrichardiiLXCXE_nuclear.csv")
cr_nuclear$TAIR

# Ferns compareCluster ####
ferns_list <- list("Afi" = af_nuclear$TAIR,
                   "Scu" = sc_nuclear$TAIR,
                   "Cri" = cr_nuclear$TAIR)

ferns_compare <- compareCluster(ferns_list, 
                                fun="enrichGO",
                                pvalueCutoff=0.05,
                                pAdjustMethod="BH",
                                keyType = "TAIR",
                                OrgDb=org.At.tair.db,
                                ont="BP")
dotplot(ferns_compare, showCategory = 35) + scale_color_viridis()

# Converting to tibble!
fern_tibble <- as_tibble(ferns_compare) %>%
  dplyr::select(Cluster,
                Description,
                p.adjust,
                geneID)
# A_filliculoides
azolla_tibble <- fern_tibble %>%
  dplyr::filter(Cluster == "Afi") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Afi = p.adjust)
# S_cucullata
salvinia_tibble <- fern_tibble %>%
  dplyr::filter(Cluster == "Scu") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Scu = p.adjust)
# C_richardii
ceratopteris_tibble <- fern_tibble %>%
  dplyr::filter(Cluster == "Cri") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Cri = p.adjust)
# Writing tibbles
write_csv(azolla_tibble, path = "Azolla_filliculoidesGO.csv")
write_csv(salvinia_tibble, path = "Salvinia_cucullataGO.csv")
write_csv(ceratopteris_tibble, path = "Ceratopteris_richardiiGO.csv")


# Performing joins to get final tibble
ferns <- dplyr::full_join(azolla_tibble,
                          salvinia_tibble,
                          by = "Description") %>%
  dplyr::full_join(ceratopteris_tibble,
                   by = "Description") %>%
  dplyr::select(Description, Afi, Scu, Cri) %>%
  dplyr::slice(1:40)

# Ferns matrix/Heatmap #####
write_csv(ferns, "ferns_mat.csv", col_names = T)
ferns_mat <- as.matrix(read.csv("ferns_mat.csv",
                                sep = ",",
                                row.names = "Description"))
ferns_matrix <- -log10(ferns_mat)
ferns_heatmap <- pheatmap(
  ferns_matrix,
  color = viridis(500, option = "D", begin = 0, end = 0.8),
  bias = 1, 
  border_color = F,
  cellwidth = 7,
  cellheight = 4,
  fontsize = 5,
  na_col = "light grey",
  clustering_distance_rows = "euclidean", 
  cluster_rows = F,
  cluster_cols = F,
  show_colnames = T,
  show_rownames = TRUE,
  main = "Ferns")
ferns_heatmap
pdf("ferns_heatmap.pdf")
print(ferns_heatmap) # Plot 1 --> in the first page of PDF
dev.off()


# Lycophytes ####
# Selaginella moellendorffii - Nuclear ####
sm <- AnnotationDbi::select(org.At.tair.db,
             keys=sm_ids,
             columns = c("ONTOLOGY",
                         "GO",
                         "SYMBOL"),
             keytype = "TAIR")
sm_nuclear <- sm %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634")) %>%  # CC GO for nucleus
  write_csv("SmoellendorffiiLXCXE_nuclear.csv")
sm_nuclear$TAIR

# Isoetes taiwanensis - Nuclear ####
it <- AnnotationDbi::select(org.At.tair.db,
             keys=it_ids,
             columns = c("ONTOLOGY",
                         "GO",
                         "SYMBOL"),
             keytype = "TAIR")
it_nuclear <- sm %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634")) %>%  # CC GO for nucleus
  write_csv("ItaiwanensisLXCXE_nuclear.csv")
it_nuclear$TAIR

#Lycophytes compareCluster ####
lycophytes_list <- list("Smo" = sm_nuclear$TAIR,
                   "Ita" = it_nuclear$TAIR)

lycophytes_compare <- compareCluster(lycophytes_list, 
                                fun="enrichGO",
                                pvalueCutoff=0.05,
                                pAdjustMethod="BH",
                                keyType = "TAIR",
                                OrgDb=org.At.tair.db,
                                ont="BP")
dotplot(lycophytes_compare, showCategory = 35) + scale_color_viridis()

# Converting to tibble!
lycophyte_tibble <- as_tibble(lycophytes_compare) %>%
  dplyr::select(Cluster,
                Description,
                p.adjust,
                geneID)
# S_moellendorffii
selaginella_tibble <- lycophyte_tibble %>%
  dplyr::filter(Cluster == "Smo") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Smo = p.adjust)
# I_taiwanensis
isoetes_tibble <- lycophyte_tibble %>%
  dplyr::filter(Cluster == "Ita") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Ita = p.adjust)

#Writing tibbles
write_csv(selaginella_tibble, path = "Selaginella_moellendorffiiGO.csv")
write_csv(isoetes_tibble, path = "Isoetes_taiwanensisGO.csv")


# Performing joins to get final tibble
lycophytes <- dplyr::full_join(selaginella_tibble,
                                      isoetes_tibble,
                                      by = "Description") %>%
  dplyr::select(Description, Smo, Ita) %>%
  dplyr::slice(1:40)

# Lycophytes matrix/Heatmap #####
write_csv(lycophytes, "lycophytes_mat.csv", col_names = T)
lycophytes_mat <- as.matrix(read.csv("lycophytes_mat.csv",
                                sep = ",",
                                row.names = "Description"))
lycophytes_matrix <- -log10(lycophytes_mat)
lycophytes_heatmap <- pheatmap(
  lycophytes_matrix,
  color = viridis(500, option = "D", begin = 0, end = 0.8),
  bias = 1, 
  border_color = F,
  cellwidth = 7,
  cellheight = 4,
  fontsize = 5,
  na_col = "light grey",
  clustering_distance_rows = "euclidean", 
  cluster_rows = F,
  cluster_cols = F,
  show_colnames = T,
  show_rownames = TRUE,
  main = "Lycophytes")
lycophytes_heatmap
pdf("lycophytes_heatmap.pdf")
print(lycophytes_heatmap) # Plot 1 --> in the first page of PDF
dev.off()

# Bryophytes ####
# Physcomitrella patens - Nuclear ####
pp <-  AnnotationDbi::select(org.At.tair.db,
                             keys=pp_ids,
                             columns = c("ONTOLOGY",
                                         "GO",
                                         "SYMBOL"),
                             keytype = "TAIR")
pp_nuclear <- pp %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634")) %>%  # CC GO for nucleus
  write_csv("PpatensLXCXE_nuclear.csv")
pp_nuclear$TAIR

# Sphagnum fallax - Nuclear ####
sf <-  AnnotationDbi::select(org.At.tair.db,
                             keys=sf_ids,
                             columns = c("ONTOLOGY",
                                         "GO",
                                         "SYMBOL"),
                             keytype = "TAIR")
sf_nuclear <- sf %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634")) %>%  # CC GO for nucleus
  write_csv("SfallaxLXCXE_nuclear.csv")
sf_nuclear$TAIR

# Marchantia polymorpha - Nuclear ####
mp <-  AnnotationDbi::select(org.At.tair.db,
                             keys=mp_ids,
                             columns = c("ONTOLOGY",
                                         "GO",
                                         "SYMBOL"),
                             keytype = "TAIR")
mp_nuclear <- mp %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634")) %>%  # CC GO for nucleus
  write_csv("MpolymorphaLXCXE_nuclear.csv")
mp_nuclear$TAIR

# Anthoceros agrestis - Nuclear ####
aa <-  AnnotationDbi::select(org.At.tair.db,
                             keys=aa_ids,
                             columns = c("ONTOLOGY",
                                         "GO",
                                         "SYMBOL"),
                             keytype = "TAIR")
aa_nuclear <- aa %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634")) %>%  # CC GO for nucleus
  write_csv("AagrestisLXCXE_nuclear.csv")
aa_nuclear$TAIR

# Bryophytes compareCluster ####
bryophytes_list <- list("Ppa" = pp_nuclear$TAIR,
                        "Sfa" = sf_nuclear$TAIR,
                        "Mpo" = mp_nuclear$TAIR,
                        "Aag" = aa_nuclear$TAIR)

bryophytes_compare <- compareCluster(bryophytes_list, 
                                     fun="enrichGO",
                                     pvalueCutoff=0.05,
                                     pAdjustMethod="BH",
                                     keyType = "TAIR",
                                     OrgDb=org.At.tair.db,
                                     ont="BP")
dotplot(bryophytes_compare, showCategory = 35) + scale_color_viridis()

# Converting to tibble!
bryophytes_tibble <- as_tibble(bryophytes_compare) %>%
  dplyr::select(Cluster,
                Description,
                p.adjust,
                geneID)

# Extracting each cluster and then combining them with full_angio_tibble
# P_patens
physcomitrella_tibble <- bryophytes_tibble %>%
  dplyr::filter(Cluster == "Ppa") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Ppa = p.adjust)
# S_fallax
sphagnum_tibble <- bryophytes_tibble %>%
  dplyr::filter(Cluster == "Sfa") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Sfa = p.adjust)
# M_polymorpha
marchantia_tibble <- bryophytes_tibble %>%
  dplyr::filter(Cluster == "Mpo") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Mpo = p.adjust)
# A_agrestis
anthoceros_tibble <- bryophytes_tibble %>%
  dplyr::filter(Cluster == "Aag") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Aag = p.adjust)

# Writing tibbles
write_csv(physcomitrella_tibble, path = "Physcomitrella_patensGO.csv")
write_csv(marchantia_tibble, path = "Marchantia_polymorphaGO.csv")
write_csv(sphagnum_tibble, path = "Sphagnum_fallaxGO.csv")
write_csv(anthoceros_tibble, path = "Anthoceros_agrestisGO.csv")

# Performing joints to get final tibble
bryophytes <- dplyr::full_join(physcomitrella_tibble,
                               sphagnum_tibble,
                               by = "Description") %>%
  dplyr::full_join(marchantia_tibble,
                   by = "Description") %>%
  dplyr::full_join(anthoceros_tibble,
                   by = "Description") %>%
  dplyr::select(Description, Ppa, Sfa, Mpo, Aag) %>%
  dplyr::slice(1:40)

# Building matrices - bryophytes #####
write_csv(bryophytes, "bryophytes_mat.csv", col_names = T)
bryophytes_mat <- as.matrix(read.csv("bryophytes_mat.csv",
                                     sep = ",",
                                     row.names = "Description"))
bryophytes_matrix <- -log10(bryophytes_mat)
bryophytes_heatmap<- pheatmap(
  bryophytes_matrix,
  color = viridis(500, option = "D", begin = 0, end = 0.8),
  bias = 1, 
  border_color = F,
  cellwidth = 7,
  cellheight = 4,
  fontsize = 5,
  na_col = "light grey",
  clustering_distance_rows = "euclidean", 
  cluster_rows = F,
  cluster_cols = F,
  show_colnames = T,
  show_rownames = TRUE,
  main = "Bryophytes")
bryophytes_heatmap # Bryophytes heatmap
pdf("bryophytes_heatmap.pdf")
print(bryophytes_heatmap) # Plot 1 --> in the first page of PDF
dev.off()

# Charophytes ####
# Spirogloea muscicola - Nuclear ####
smu <-  AnnotationDbi::select(org.At.tair.db,
                              keys=smu_ids,
                              columns = c("ONTOLOGY",
                                          "GO",
                                          "SYMBOL"),
                              keytype = "TAIR")
smu_nuclear <- smu %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634")) %>%  # CC GO for nucleus
  write_csv("SmuscicolaLXCXE_nuclear.csv")
smu_nuclear$TAIR

# Mesotaenium endlicheranium - Nuclear ####
me <-  AnnotationDbi::select(org.At.tair.db,
                             keys=me_ids,
                             columns = c("ONTOLOGY",
                                         "GO",
                                         "SYMBOL"),
                             keytype = "TAIR")
me_nuclear <- me %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634")) %>%  # CC GO for nucleus
  write_csv("MendlicheraniumLXCXE_nuclear.csv")
me_nuclear$TAIR

# Peniums margaritaceum - Nuclear ####
pma <-  AnnotationDbi::select(org.At.tair.db,
                              keys=pma_ids,
                              columns = c("ONTOLOGY",
                                          "GO",
                                          "SYMBOL"),
                              keytype = "TAIR")
pma_nuclear <- pma %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634"))  %>%  # CC GO for nucleus
  write_csv("PmargaritaceumLXCXE_nuclear.csv")
pma_nuclear$TAIR

# Chara braunii - Nuclear ####
cb <-  AnnotationDbi::select(org.At.tair.db,
                             keys=cb_ids,
                             columns = c("ONTOLOGY",
                                         "GO",
                                         "SYMBOL"),
                             keytype = "TAIR")
cb_nuclear <- cb %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634"))  %>%  # CC GO for nucleus
  write_csv("CbrauniiLXCXE_nuclear.csv")
cb_nuclear$TAIR

# Klebsormidium nitens - Nuclear ####
kn <-  AnnotationDbi::select(org.At.tair.db,
                             keys=kn_ids,
                             columns = c("ONTOLOGY",
                                         "GO",
                                         "SYMBOL"),
                             keytype = "TAIR")
kn_nuclear <- kn %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634"))  %>%  # CC GO for nucleus
  write_csv("KnitensLXCXE_nuclear.csv")
kn_nuclear$TAIR

# Chlorokybus atmophyticus - Nuclear ####
ca <-  AnnotationDbi::select(org.At.tair.db,
                             keys=ca_ids,
                             columns = c("ONTOLOGY",
                                         "GO",
                                         "SYMBOL"),
                             keytype = "TAIR")
ca_nuclear <- ca %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634"))  %>%  # CC GO for nucleus
  write_csv("CatmophyticusLXCXE_nuclear.csv")
ca_nuclear$TAIR

# Mesostigma viride - Nuclear ####
mv <-  AnnotationDbi::select(org.At.tair.db,
                             keys=mv_ids,
                             columns = c("ONTOLOGY",
                                         "GO",
                                         "SYMBOL"),
                             keytype = "TAIR")
mv_nuclear <- mv %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634"))  %>%  # CC GO for nucleus
  write_csv("MvirideLXCXE_nuclear.csv")
mv_nuclear$TAIR

# Charophytes compareCluster ####
charophytes_list <- list("Smu" = smu_nuclear$TAIR,
                         "Men" = me_nuclear$TAIR,
                         "Pma" = pma_nuclear$TAIR,
                         "Cbr" = cb_nuclear$TAIR,
                         "Kni" = kn_nuclear$TAIR,
                         "Cat" = ca_nuclear$TAIR,
                         "Mvi" = mv_nuclear$TAIR)

charophytes_compare <- compareCluster(charophytes_list, 
                                      fun="enrichGO",
                                      pvalueCutoff=0.05,
                                      pAdjustMethod="BH",
                                      keyType = "TAIR",
                                      OrgDb=org.At.tair.db,
                                      ont="BP")
dotplot(charophytes_compare, showCategory = 35) + scale_color_viridis()

# Converting to tibble!
charophytes_tibble <- as_tibble(charophytes_compare) %>%
  dplyr::select(Cluster,
                Description,
                p.adjust, geneID)

# Extracting each cluster and then combining them with full_angio_tibble
# S_muscicola
spirogloea_tibble <- charophytes_tibble %>%
  dplyr::filter(Cluster == "Smu") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Smu = p.adjust)
# M_endlicheranium
mesotaenium_tibble <- charophytes_tibble %>%
  dplyr::filter(Cluster == "Men") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Men = p.adjust)
# P_margaritaceum
penium_tibble <- charophytes_tibble %>%
  dplyr::filter(Cluster == "Pma") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Pma = p.adjust)
# C_braunii
chara_tibble <- charophytes_tibble %>%
  dplyr::filter(Cluster == "Cbr") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Cbr = p.adjust)
# K_nitens
klebsormidium_tibble <- charophytes_tibble %>%
  dplyr::filter(Cluster == "Kni") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Kni = p.adjust)
# C_atmophyticus
chlorokybus_tibble <- charophytes_tibble %>%
  dplyr::filter(Cluster == "Cat") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Cat = p.adjust)
# M_viride
mesostigma_tibble <- charophytes_tibble %>%
  dplyr::filter(Cluster == "Mvi") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Mvi = p.adjust)

# Writing tibbles
write_csv(penium_tibble, path = "Penium_margaritaceumGO.csv")
write_csv(spirogloea_tibble, path = "Spirogloea_muscicolaGO.csv")
write_csv(mesotaenium_tibble, path = "Mesotaenium_endlicheraniumGO.csv")
write_csv(chara_tibble, path = "Chara_brauniGO.csv")
write_csv(klebsormidium_tibble, path = "Klebsormidium_nitensGO.csv")
write_csv(chlorokybus_tibble, path = "Chlorokybus_atmophyticusGO.csv")
write_csv(mesostigma_tibble, path = "Mesostigma_virideGO.csv")

# Performing joints to get final tibble
charophytes <- dplyr::full_join(spirogloea_tibble,
                                mesotaenium_tibble,
                                by = "Description") %>%
  dplyr::full_join(penium_tibble,
                   by = "Description") %>%
  dplyr::full_join(chara_tibble,
                   by = "Description") %>%
  dplyr::full_join(klebsormidium_tibble,
                   by = "Description") %>%
  dplyr::full_join(chlorokybus_tibble,
                   by = "Description") %>%
  dplyr::full_join(mesostigma_tibble,
                   by = "Description") %>%
  dplyr::select(Description, Smu, Men, Pma, Cbr, Kni, Cat, Mvi) %>%
  dplyr::slice(1:40)

# Charophytes matrix/Heatmap ####
write_csv(charophytes, "charophytes_mat.csv", col_names = T)
charophytes_mat <- as.matrix(read.csv("charophytes_mat.csv",
                                      sep = ",",
                                      row.names = "Description"))
charophytes_matrix <- -log10(charophytes_mat)
charophytes_heatmap<- pheatmap(
  charophytes_matrix,
  color = viridis(500, option = "D", begin = 0, end = 0.8),
  bias = 1, 
  border_color = F,
  cellwidth = 7,
  cellheight = 4,
  fontsize = 5,
  na_col = "light grey",
  clustering_distance_rows = "euclidean", 
  cluster_rows = F,
  cluster_cols = F,
  show_colnames = T,
  show_rownames = TRUE,
  main = "Charophytes")
charophytes_heatmap
pdf("charophytes_heatmap.pdf")
print(charophytes_heatmap) # Plot 1 --> in the first page of PDF
dev.off()

# Chlorophytes ####
# Volvox carteri - Nuclear ####
vc <-  AnnotationDbi::select(org.At.tair.db,
                             keys=vc_ids,
                             columns = c("ONTOLOGY",
                                         "GO",
                                         "SYMBOL"),
                             keytype = "TAIR")
vc_nuclear <- vc %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634"))  %>%  # CC GO for nucleus
  write_csv("VcarteriXCXE_nuclear.csv")
vc_nuclear$TAIR

# Chlamydomonas reinhardtii - Nuclear ####
cre <-  AnnotationDbi::select(org.At.tair.db,
                             keys=cre_ids,
                             columns = c("ONTOLOGY",
                                         "GO",
                                         "SYMBOL"),
                             keytype = "TAIR")
cre_nuclear <- cr %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634"))  %>%  # CC GO for nucleus
  write_csv("CreinhardtiiLXCXE_nuclear.csv")
cre_nuclear$TAIR

# Ostreococcus lucimarinus - Nuclear ####
ol <-  AnnotationDbi::select(org.At.tair.db,
                             keys=ol_ids,
                             columns = c("ONTOLOGY",
                                         "GO",
                                         "SYMBOL"),
                             keytype = "TAIR")
ol_nuclear <- ol %>%
  dplyr::select(TAIR,
                ONTOLOGY,
                GO,
                SYMBOL) %>%
  dplyr::filter(GO %in% as.vector("GO:0005634"))  %>%  # CC GO for nucleus
  write_csv("OlucimarinusLXCXE_nuclear.csv")
ol_nuclear$TAIR

# Chlorophytes compareCluster ####
chlorophytes_list <- list("Vca" = vc_nuclear$TAIR,
                          "Cre" = cr_nuclear$TAIR,
                          "Olu" = ol_nuclear$TAIR)

chlorophytes_compare <- compareCluster(chlorophytes_list, 
                                       fun="enrichGO",
                                       pvalueCutoff=0.05,
                                       pAdjustMethod="BH",
                                       keyType = "TAIR",
                                       OrgDb=org.At.tair.db,
                                       ont="BP")
dotplot(chlorophytes_compare, showCategory = 35) + scale_color_viridis()

# Converting to tibble!
chlorophyte_tibble <- as_tibble(chlorophytes_compare) %>%
  dplyr::select(Cluster,
                Description,
                p.adjust,
                geneID)
# Extracting each cluster and then combining them with full_angio_tibble
# V_carteri
volvox_tibble <- chlorophyte_tibble %>%
  dplyr::filter(Cluster == "Vca") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Vca = p.adjust)
# C_reinhardtii
chlamydomonas_tibble <- chlorophyte_tibble %>%
  dplyr::filter(Cluster == "Cre") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Cre = p.adjust)
# O_lucimarinus
ostreococcus_tibble <- chlorophyte_tibble %>%
  dplyr::filter(Cluster == "Olu") %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::rename(Olu= p.adjust)
# Writing tibbles
write_csv(volvox_tibble, path = "Volvox_carteriGO.csv")
write_csv(chlamydomonas_tibble, path = "Chlamydomonas_reinhardtiiGO.csv")
write_csv(ostreococcus_tibble, path = "Ostreococcus_lucimarinusGO.csv")

# Performing joints to get final tibble
chlorophytes_complete <- dplyr::full_join(volvox_tibble,
                                          chlamydomonas_tibble,
                                          by = "Description") %>%
  dplyr::full_join(ostreococcus_tibble,
                   by = "Description") %>%
  dplyr::select(Description, Vca, Cre, Olu) %>%
  dplyr::slice(1:40)

# Building matrices - chlorophytes #####
write_csv(chlorophytes_complete, "chlorophytes_mat.csv", col_names = T)
chlorophytes_mat <- as.matrix(read.csv("chlorophytes_mat.csv",
                                       sep = ",",
                                       row.names = "Description"))
chlorophytes_matrix <- -log10(chlorophytes_mat)
chlorophyte_heatmap<- pheatmap(
  chlorophytes_matrix,
  color = viridis(500, option = "D", begin = 0, end = 0.8),
  bias = 1, 
  border_color = F,
  cellwidth = 7,
  cellheight = 4,
  fontsize = 5,
  na_col = "light grey",
  clustering_distance_rows = "euclidean", 
  cluster_rows = F,
  cluster_cols = F,
  show_colnames = T,
  show_rownames = TRUE,
  main = "Chlorophytes")
chlorophyte_heatmap
pdf("chlorophyte_heatmap.pdf")
print(chlorophyte_heatmap) # Plot 1 --> in the first page of PDF
dev.off()


#  Tibbles - Viridiplantae Enrichment ####
# Complete Angiosperms tibble
angiosperms_complete <- dplyr::full_join(arabidopsis_tibble,
                                         boechera_tibble,
                                         by = "Description") %>%
  dplyr::full_join(solanum_tibble,
                   by = "Description") %>%
  dplyr::full_join(zea_tibble,
                   by = "Description") %>%
  dplyr::full_join(oryza_tibble,
                   by = "Description") %>%
  dplyr::full_join(setaria_tibble,
                   by = "Description") %>% 
  dplyr::full_join(amborella_tibble,
                   by = "Description") %>%
  dplyr::select(Description, Ath, Bst, Sly, Zma, Osa, Svi, Atr)

# Complete Gymnosperms tibble
gymnosperms_complete<- dplyr::full_join(picea_tibble,
                                        pinus_tibble,
                                        by = "Description") %>%
  dplyr::full_join(gnetum_tibble,
                   by = "Description") %>%
  dplyr::select(Description, Pab, Pta, Gmo)

# Complete Ferns tibble
ferns_complete <- dplyr::full_join(azolla_tibble,
                                   salvinia_tibble,
                                   by = "Description") %>%
  dplyr::full_join(ceratopteris_tibble,
                   by = "Description") %>%
  dplyr::select(Description, Afi, Scu, Cri)

# Complete lycophyte tibble
lycophytes_complete <- dplyr::full_join(selaginella_tibble,
                                   isoetes_tibble,
                                   by = "Description") %>%
  dplyr::select(Description, Smo, Ita) 

# Complete Bryophyte tibble
bryophytes_complete <- dplyr::full_join(physcomitrella_tibble,
                                        sphagnum_tibble,
                                        by = "Description") %>%
  dplyr::full_join(marchantia_tibble,
                   by = "Description") %>%
  dplyr::full_join(anthoceros_tibble,
                   by = "Description") %>%
  dplyr::select(Description, Ppa, Sfa, Mpo, Aag) 

# Complete Charophyte tibble
charophytes_complete <- dplyr::full_join(spirogloea_tibble,
                                         mesotaenium_tibble,
                                         by = "Description") %>%
  dplyr::full_join(penium_tibble,
                   by = "Description") %>%
  dplyr::full_join(chara_tibble,
                   by = "Description") %>%
  dplyr::full_join(klebsormidium_tibble,
                   by = "Description") %>%
  dplyr::full_join(chlorokybus_tibble,
                   by = "Description") %>%
  dplyr::full_join(mesostigma_tibble,
                   by = "Description") %>%
  dplyr::select(Description, Smu, Men, Pma, Cbr, Kni, Cat, Mvi)

# Complete Chlorophyte tibble
chlorophytes_complete <- dplyr::full_join(volvox_tibble,
                                          chlamydomonas_tibble,
                                          by = "Description") %>%
  dplyr::full_join(ostreococcus_tibble,
                   by = "Description") %>%
  dplyr::select(Description, Vca, Cre, Olu)

# Viridiplantae enrichment - inner join
viridiplantae <- dplyr::inner_join(angiosperms_complete,
                                   gymnosperms_complete,
                                   by = "Description") %>%
  dplyr::inner_join(ferns_complete,
                    by = "Description") %>%
  dplyr::inner_join(lycophytes_complete,
                    by = "Description") %>%
  dplyr::inner_join(bryophytes_complete,
                    by = "Description") %>%
  dplyr::inner_join(charophytes_complete,
                    by = "Description") %>%
  dplyr::inner_join(chlorophytes_complete,
                    by = "Description") %>%
  dplyr::select(Description, Ath, Bst, Sly, Zma, Osa, Svi, Atr,
                Pab, Pta, Gmo, Afi, Scu, Cri, Smo, Ita, Ppa, Sfa, Mpo, Aag,
                Smu, Men, Pma, Cbr, Kni, Cat, Mvi, Vca, Cre, Olu) %>%
  dplyr::slice(1:50) 
# Building matrices - viridiplantae #####
write_csv(viridiplantae, "viridiplantae_mat.csv", col_names = T)
viridiplantae_mat <- as.matrix(read.csv("viridiplantae_mat.csv",
                                        sep = ",",
                                        row.names = "Description"))
viridiplantae_matrix <- -log10(viridiplantae_mat)
viridiplantae_heatmap_inner<- pheatmap(
  viridiplantae_matrix,
  color = viridis(500, option = "D", begin = 0, end = 0.8),
  bias = 1, 
  border_color = F,
  cellwidth = 7,
  cellheight = 4,
  fontsize = 5,
  na_col = "light grey",
  clustering_distance_rows = "euclidean", 
  cluster_rows = F,
  cluster_cols = F,
  show_colnames = T,
  show_rownames = TRUE,
  main = "LxCxE Viridiplantae enriched GO BP")
viridiplantae_heatmap_inner
pdf("viridiplantae_heatmap_inner-joinBP.pdf")
print(viridiplantae_heatmap_inner) # Plot 1 --> in the first page of PDF
dev.off()

# Viridiplantae enrichment - full join
viridiplantae_full <- dplyr::full_join(angiosperms_complete,
                                       gymnosperms_complete,
                                       by = "Description") %>%
  dplyr::full_join(ferns_complete,
                   by = "Description") %>%
  dplyr::full_join(lycophytes_complete,
                   by = "Description") %>%
  dplyr::full_join(bryophytes_complete,
                   by = "Description") %>%
  dplyr::full_join(charophytes_complete,
                   by = "Description") %>%
  dplyr::full_join(chlorophytes_complete,
                   by = "Description") %>%
  dplyr::select(Description, Ath, Bst, Sly, Zma, Osa, Svi, Atr,
                Pab, Pta, Gmo, Afi, Scu, Cri, Smo, Ita, Ppa, Sfa, Mpo, Aag,
                Smu, Men, Pma, Cbr, Kni, Cat, Mvi, Vca, Cre, Olu) %>%
  dplyr::slice(1:50) 
# Building matrices - viridiplantae #####
write_csv(viridiplantae_full, "viridiplantae_mat_full.csv", col_names = T)
viridiplantae_mat_full <- as.matrix(read.csv("viridiplantae_mat_full.csv",
                                             sep = ",",
                                             row.names = "Description"))
viridiplantae_matrix <- -log10(viridiplantae_mat_full)
viridiplantae_heatmap_full<- pheatmap(
  viridiplantae_matrix,
  color = viridis(500, option = "D", begin = 0, end = 0.8),
  bias = 1, 
  border_color = F,
  cellwidth = 7,
  cellheight = 4,
  fontsize = 5,
  na_col = "light grey",
  clustering_distance_rows = "euclidean", 
  cluster_rows = F,
  cluster_cols = F,
  show_colnames = T,
  show_rownames = TRUE,
  main = "LxCxE Viridiplantae enriched GO BP")
viridiplantae_heatmap_full

pdf("viridiplantae_heatmap_full-joinBP.pdf")
print(viridiplantae_heatmap_full) # Plot 1 --> in the first page of PDF
dev.off()

# ggplot objects ####
library(ggplot2)
library(ggplotify)
angiosperm_ggplot <- as.ggplot(angiosperm_heatmap)
gymnosperm_ggplot <- as.ggplot(gymnosperm_heatmap)
fern_ggplot <- as.ggplot(ferns_heatmap)
bryophyte_ggplot <- as.ggplot(bryophytes_heatmap)
charophyte_ggplot <- as.ggplot(charophytes_heatmap)
chlorophyte_ggplot <- as.ggplot(chlorophyte_heatmap)

library(cowplot)
pdf(file="viridiplantae_GO-heatmap-species_ORTHOFINDER.pdf")
viridiplantae<-plot_grid(angiosperm_ggplot,
                         gymnosperm_ggplot,
                         fern_ggplot,
                         bryophyte_ggplot,
                         charophyte_ggplot,
                         chlorophyte_ggplot,
                         nrow = 2,
                         labels = "AUTO",
                         label_size = 12,
                         align = "h",
                         scale = 1)
save_plot("viridiplantae.pdf",
          viridiplantae,
          ncol = 3,
          nrow = 2,
          base_height = 3,
          base_width = 3.7)
save_plot("viridiplantae.png",
          viridiplantae,
          ncol = 3,
          nrow = 2,
          base_height = 3,
          base_width = 3.7)
# Done
