# Package for analyzing LxCxE-Transcription Associated Proteins
library(readr)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(gapminder)
library(viridis)
library(dplyr)
# Setting working directory to Workdir
setwd("Workdir/")

# Loading TAPdb family-domain data
TAPdb <- read_csv("taps_simplified.csv") %>%
  dplyr::select(Family, type, domain)

# Loading and processing raw datasests ####
lxcxe_dataframes <- list.files(pattern = "*.fasta.csv", 
                               full.names = TRUE)

# adding column names and processing data FOR FIRST OUTPUT!!
for (lxcxe_dataframe in lxcxe_dataframes) {
  data <- read_table2(lxcxe_dataframe, 
                   col_names = c("domain",
                                  "accession",
                                  "query",
                                  "accession2",
                                  "evalue",
                                  "score",
                                  "bias",
                                  "exp",
                                  "reg",
                                  "clu",
                                  "ov",
                                  "env",
                                  "dom",
                                  "rep",
                                  "inc",
                                  "description")) %>%
    dplyr::select(domain,
                  accession,
                  query,
                  evalue,
                  score)  %>%
    dplyr::inner_join(TAPdb,
                      by = "domain") %>%
    dplyr::distinct(query,
                    Family,
                    .keep_all = T) %>%
    dplyr::mutate_at(vars(evalue, score),
                     funs(as.numeric)) %>%
    dplyr::filter(score > 10,
                  evalue < 0.00005) %>% 
    dplyr::arrange(domain, .by_group = TRUE) %>%
    dplyr::group_by(Family) %>%
    dplyr::mutate(count = n()) %>%
    drop_na() %>%
    write_csv(lxcxe_dataframe, "*.csv")
}



# Loading processed data
ath_hmmer <- read_csv("Athaliana_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na()
bst_hmmer <- read_csv("Bstricta_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na()
sly_hmmer <- read_csv("Slycopersicum_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na()
zma_hmmer <- read_csv("Zmays_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na() 
osa_hmmer <- read_csv("Osativa_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na()
svi_hmmer <- read_csv("Sviridis_LXCXE.fasta.csv") %>% group_by(Family )%>% drop_na()
atr_hmmer <- read_csv("Atrichopoda_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na()
pta_hmmer <- read_csv("Ptaeda_LXCXE.csv") %>% group_by(Family) %>% drop_na()
pab_hmmer <- read_csv("Pabies_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na()
gmo_hmmer <- read_csv("Gmontanum_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na()
afi_hmmer <- read_csv("Afilliculoides_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na()
scu_hmmer <- read_csv("Scucullata_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na()
cri_hmmer <- read_csv("Crichardii_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na()
smo_hmmer <- read_csv("Smoellendorffii_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na()
ppa_hmmer <- read_csv("Ppatens_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na()
sfa_hmmer <- read_csv("Sfallax_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na()
mpo_hmmer <- read_csv("Mpolymorpha_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na()
aag_hmmer <- read_csv("Aagrestis_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na()
pma_hmmer <- read_csv("Pmargaritaceum_LXCXE.fasta.csv") %>% group_by(Family)  %>% drop_na()
men_hmmer <- read_csv("Mendlicheranium_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na()
smu_hmmer <- read_csv("Smuscicola_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na()
cbr_hmmer <- read_csv("Cbraunii_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na()
kni_hmmer <- read_csv("Knitens_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na()
cat_hmmer <- read_csv("Catmophyticus_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na()
mvi_hmmer <- read_csv("MvirideLIANG_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na()
vca_hmmer <- read_csv("Vcarteri_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na()
cre_hmmer <- read_csv("Creinhardtii_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na()
olu_hmmer <- read_csv("Olucimarinus_LXCXE.fasta.csv") %>% group_by(Family) %>% drop_na()

# Writing files to manually curate some data #
#write_csv(ath_hmmer, file = "AthalianaLXCXE_TAPs.csv")
#write_csv(bst_hmmer, file = "BstrictaLXCXE_TAPs.csv")
#write_csv(sly_hmmer, file = "SlycopersicumLXCXE_TAPs.csv")
#write_csv(zma_hmmer, file = "ZmaysLXCXE_TAPs.csv")
#write_csv(osa_hmmer, file = "OsativaLXCXE_TAPs.csv")
#write_csv(svi_hmmer, file = "SviridisLXCXE_TAPs.csv")
#write_csv(atr_hmmer, file = "AtrichopodaLXCXE_TAPs.csv")
#write_csv(pta_hmmer, file = "PtaedaLXCXE_TAPs.csv")
#write_csv(pab_hmmer, file = "PabiesLXCXE_TAPs.csv")
#write_csv(gmo_hmmer, file = "GmontanumLXCXE_TAPs.csv")
#write_csv(afi_hmmer, file = "AfilliculoidesLXCXE_TAPs.csv")
#write_csv(scu_hmmer, file = "ScucullataLXCXE_TAPs.csv")
#write_csv(cri_hmmer, file = "CrichardiiLXCXE_TAPs.csv")
#write_csv(smo_hmmer, file = "SmoellendorffiiLXCXE_TAPs.csv")
#write_csv(ppa_hmmer, file = "PpatensLXCXE_TAPs.csv")
#write_csv(sfa_hmmer, file = "SfallaxLXCXE_TAPs.csv")
#write_csv(mpo_hmmer, file = "MpolymorphaLXCXE_TAPs.csv")
#write_csv(aag_hmmer, file = "AagrestisLXCXE_TAPs.csv")
#write_csv(pma_hmmer, file = "PmargaritaceumLXCXE_TAPs.csv")
#write_csv(smu_hmmer, file = "SmuscicolaLXCXE_TAPs.csv")
#write_csv(men_hmmer, file = "MendlicheraniumLXCXE_TAPs.csv")
#write_csv(cbr_hmmer, file = "CbrauniiLXCXE_TAPs.csv")
#write_csv(kni_hmmer, file = "KnitensLXCXE_TAPs.csv")
#write_csv(cat_hmmer, file = "CatmophyticusLXCXE_TAPs.csv")
#write_csv(mvi_hmmer, file = "MvirideLXCXE_TAPs.csv")
#write_csv(vca_hmmer, file = "VcarteriLXCXE_TAPs.csv")
#write_csv(cre_hmmer, file = "CreinhardtiiLXCXE_TAPs.csv")
#write_csv(olu_hmmer, file = "OlucimarinusLXCXE_TAPs.csv")

# Manual curation of TAPs for some families ...









# Viridiplantae dataframe ##### 
test  <- dplyr::full_join(ath_hmmer, sly_hmmer, by = "Family") %>%
  select(Family, count.x, count.y, type.x, type.y) %>%
  distinct(Family, .keep_all = T) %>% 
  rename(Ath = count.x,
         Sly = count.y) # worked


viridiplant_taps_full <- dplyr::full_join(ath_hmmer, #Viridiplantae TAPs Full Join ####
                                          bst_hmmer,
                                          by = "Family") %>%
  dplyr::select(Family, count.x, count.y) %>%
  rename(Ath = count.x,
         Bst = count.y) %>%
  distinct(Family, .keep_all = T) %>%
  full_join(sly_hmmer, # Solanum lycopersicum
            by = "Family") %>%
  select(Family,
         Ath,
         Bst,
         count) %>%
  rename(Sly = count) %>%
  full_join(zma_hmmer, # Zea mays
            by = "Family") %>%
  select(Family,
         Ath, 
         Bst,
         Sly,
         count) %>%
  rename(Zma = count) %>%
  full_join(osa_hmmer, # Oryza sativa
            by = "Family") %>%
  select(Family,
         Ath,
         Bst,
         Sly,
         Zma,
         count) %>%
  rename(Osa = count) %>%
  full_join(svi_hmmer, # Setaria viridis
            by = "Family") %>%
  select(Family,
         Ath,
         Bst,
         Sly,
         Zma,
         Osa,
         count) %>%
  rename(Svi = count) %>% distinct(Family,
                                   .keep_all = T) %>%
  full_join(atr_hmmer, # Amborella trichopoda
            by  = "Family") %>% # Gymnosperms
  select(Family,
         Ath,
         Bst,
         Sly,
         Zma,
         Osa,
         Svi,
         count) %>% 
  rename(Atr = count) %>%
  full_join(pab_hmmer, # Picea abies
            by = "Family") %>%
  select(Family,
         Ath,
         Bst,
         Sly,
         Zma,
         Osa,
         Svi,
         Atr,
         count) %>% 
  rename(Pab = count) %>% distinct(Family,
                                   .keep_all = T) %>%
  full_join(pta_hmmer, # Pinus taeda
            by = "Family") %>%
  select(Family,
         Ath,
         Bst,
         Sly,
         Zma,
         Osa,
         Svi,
         Atr,
         Pab,
         count) %>%  
  rename(Pta = count) %>%
  full_join(gmo_hmmer, # Gnetum montanum
            by = "Family") %>%
  select(Family,
         Ath,
         Bst,
         Sly,
         Zma,
         Osa,
         Svi,
         Atr,
         Pab,
         Pta,
         count) %>%
  rename(Gmo = count) %>%
  full_join(afi_hmmer, # Azolla filliculoides
            by = "Family") %>%
  select(Family,
         Ath,
         Bst,
         Sly,
         Zma,
         Osa,
         Svi,
         Atr,
         Pab,
         Pta,
         Gmo,
         count) %>%
  rename(Afi = count) %>%
  full_join(scu_hmmer, # Salvinia cucullata
            by = "Family") %>%
  select(Family,
         Ath,
         Bst,
         Sly,
         Zma,
         Osa,
         Svi,
         Atr,
         Pab,
         Pta,
         Gmo,
         Afi,
         count) %>%
  rename(Scu = count) %>% distinct(Family, .keep_all = T) %>%
  full_join(cri_hmmer,
            by = "Family") %>%  # Ceratopteris richardii 
  select(Family, 
         Ath,
         Bst,
         Sly,
         Zma,
         Osa,
         Svi,
         Atr,
         Pab,
         Pta,
         Gmo,
         Afi,
         Scu,
         count) %>%
  rename(Cri = count) %>% distinct(Family, .keep_all = T)%>%
  full_join(smo_hmmer, # Selaginella moellendorffii
            by = "Family") %>%
  select(Family,
         Ath,
         Bst,
         Sly,
         Zma,
         Osa,
         Svi,
         Atr,
         Pab,
         Pta,
         Gmo,
         Afi,
         Scu,
         Cri,
         count) %>%
  rename(Smo = count) %>% 
  full_join(ppa_hmmer, # Physcomitrella patens
            by = "Family") %>%
  select(Family,
         Ath,
         Bst,
         Sly,
         Zma,
         Osa,
         Svi,
         Atr,
         Pab,
         Pta,
         Gmo,
         Afi,
         Scu,
         Cri,
         Smo,
         count) %>%
  rename(Ppa = count) %>%
  full_join(sfa_hmmer, # Sphagnum fallax
            by = "Family") %>%
  select(Family,
         Ath,
         Bst,
         Sly,
         Zma,
         Osa,
         Svi,
         Atr,
         Pab,
         Pta,
         Gmo,
         Afi,
         Scu,
         Cri,
         Smo,
         Ppa,
         count) %>%
  rename(Sfa = count) %>% distinct(Family, .keep_all = T) %>%
  full_join(mpo_hmmer, # Marchantia polymorpha
            by = "Family") %>%
  select(Family,
         Ath,
         Bst,
         Sly,
         Zma,
         Osa,
         Svi,
         Atr,
         Pab,
         Pta,
         Gmo,
         Afi,
         Scu,
         Cri,
         Smo,
         Ppa,
         Sfa,
         count) %>%
  rename(Mpo = count) %>%
  full_join(aag_hmmer, # Anthoceros agrestis
            by = "Family") %>%
  select(Family,
         Ath,
         Bst,
         Sly,
         Zma,
         Osa,
         Svi,
         Atr,
         Pab,
         Pta,
         Gmo,
         Afi,
         Scu,
         Cri,
         Smo,
         Ppa,
         Sfa,
         Mpo,
         count) %>%
  rename(Aag = count) %>% distinct(Family, .keep_all = T) %>%
  full_join(pma_hmmer, # Penium margaritaceum
            by = "Family") %>%
  select(Family,
         Ath,
         Bst,
         Sly,
         Zma,
         Osa,
         Svi,
         Atr,
         Pab,
         Pta,
         Gmo,
         Afi,
         Scu,
         Cri,
         Smo,
         Ppa,
         Sfa,
         Mpo,
         Aag,
         count) %>%
  rename(Pma = count) %>%
  full_join(smu_hmmer, # Spirogloea muscicola
            by = "Family") %>%
  select(Family,
         Ath,
         Bst,
         Sly,
         Zma,
         Osa,
         Svi,
         Atr,
         Pab,
         Pta,
         Gmo,
         Afi,
         Scu,
         Cri,
         Smo,
         Ppa,
         Sfa,
         Mpo,
         Aag,
         Pma,
         count) %>% 
  rename(Smu = count) %>% distinct(Family, .keep_all = T) %>%
  full_join(men_hmmer, # Mesotaenium endlicheranium
            by = "Family") %>%
  select(Family,
         Ath,
         Bst,
         Sly,
         Zma,
         Osa,
         Svi,
         Atr,
         Pab,
         Pta,
         Gmo,
         Afi,
         Scu,
         Cri,
         Smo,
         Ppa,
         Sfa,
         Mpo,
         Aag,
         Pma,
         Smu,
         count) %>%
  rename(Men = count) %>%
  full_join(cbr_hmmer, # Chara braunii
            by = "Family") %>%
  select(Family,
         Ath,
         Bst,
         Sly,
         Zma,
         Osa,
         Svi,
         Atr,
         Pab,
         Pta,
         Gmo,
         Afi,
         Scu,
         Cri,
         Smo,
         Ppa,
         Sfa,
         Mpo,
         Aag,
         Pma,
         Smu,
         Men,
         count) %>%
  rename(Cbr = count) %>%
  full_join(kni_hmmer, # Klebsormidium nitens
            by = "Family") %>%
  select(Family,
         Ath,
         Bst,
         Sly,
         Zma,
         Osa,
         Svi,
         Atr,
         Pab,
         Pta,
         Gmo,
         Afi,
         Scu,
         Cri,
         Smo,
         Ppa,
         Sfa,
         Mpo,
         Aag,
         Pma,
         Smu,
         Men,
         Cbr,
         count) %>%
  rename(Kni = count) %>% distinct(Family, .keep_all = T) %>%
  full_join(cat_hmmer, # Chlorokybus atmophyticus
            by = "Family") %>%
  select(Family,
         Ath,
         Bst,
         Sly,
         Zma,
         Osa,
         Svi,
         Atr,
         Pab,
         Pta,
         Gmo,
         Afi,
         Scu,
         Cri,
         Smo,
         Ppa,
         Sfa,
         Mpo,
         Aag,
         Pma,
         Smu,
         Men,
         Cbr,
         Kni,
         count) %>%
  rename(Cat = count) %>%
  full_join(mvi_hmmer, # Mesostigma viride
            by = "Family") %>%
  select(Family,
         Ath,
         Bst,
         Sly,
         Zma,
         Osa,
         Svi,
         Atr,
         Pab,
         Pta,
         Gmo,
         Afi,
         Scu,
         Cri,
         Smo,
         Ppa,
         Sfa,
         Mpo,
         Aag,
         Pma,
         Smu,
         Men,
         Cbr,
         Kni,
         Cat,
         count) %>%
  rename(Mvi = count) %>%
  full_join(vca_hmmer, # Volvox carteri
            by = "Family") %>%
  select(Family,
         Ath,
         Bst,
         Sly,
         Zma,
         Osa,
         Svi,
         Atr,
         Pab,
         Pta,
         Gmo,
         Afi,
         Scu,
         Cri,
         Smo,
         Ppa,
         Sfa,
         Mpo,
         Aag,
         Pma,
         Smu,
         Men,
         Cbr,
         Kni,
         Cat,
         Mvi,
         count) %>%
  rename(Vca= count) %>% distinct(Family, .keep_all = T) %>%
  full_join(cre_hmmer,
            by = "Family") %>%
  select(Family,
         Ath,
         Bst,
         Sly,
         Zma,
         Osa,
         Svi,
         Atr,
         Pab,
         Pta,
         Gmo,
         Afi,
         Scu,
         Cri,
         Smo,
         Ppa,
         Sfa,
         Mpo,
         Aag,
         Pma,
         Smu,
         Men,
         Cbr,
         Kni,
         Cat,
         Mvi,
         Vca,
         count) %>%
  rename(Cre = count) %>% distinct(Family, .keep_all = T) %>%
  full_join(olu_hmmer, by = "Family") %>%
  select(Family,
         Ath,
         Bst,
         Sly,
         Zma,
         Osa,
         Svi,
         Atr,
         Pab,
         Pta,
         Gmo,
         Afi,
         Scu,
         Cri,
         Smo,
         Ppa,
         Sfa,
         Mpo,
         Aag,
         Pma,
         Smu,
         Men,
         Cbr,
         Kni,
         Cat,
         Mvi,
         Vca,
         Cre,
         count) %>%
  rename(Olu = count) %>% distinct(Family, .keep_all = T) %>%
  group_by("Family") %>%
  slice(1:100, .preserve = T) %>%
  ungroup() %>%
  select(Family, Ath, Bst, Sly, Zma, Osa, Svi, Atr, Pta, Pab, Gmo, Afi, Scu, Cri, Smo, Ppa, Sfa, Mpo, Aag, Pma, Smu,
         Men, Cbr, Kni, Cat, Mvi, Vca, Cre, Olu)
dim(viridiplant_taps_full)


# Building matrices - viridiplantae taps #####
write_csv(viridiplant_taps_full, "viridiplantae_taps_full_mat.csv", col_names = T)
viridiplantae_taps_full_mat <- as.matrix(read.csv("viridiplantae_taps_full_mat.csv",
                                                  sep = ",",
                                                  row.names = "Family"))
viridiplantae_taps_matrix <- log2(viridiplantae_taps_full_mat)
viridiplantae_taps_full_heatmap<- pheatmap(
  viridiplantae_taps_matrix,
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
  main = "LxCxE TAPs family distribution")
viridiplantae_taps_full_heatmap
# Saving plot
pdf("viridiplantae_taps_full_heatmap_HMMSCAN.pdf")
print(viridiplantae_taps_full_heatmap) # Plot 1 --> in the first page of PDF
dev.off()


##### DONE
# Angiosperm TAP averages
angiosperm_taps_averages <-full_join(at, bs, by = "Family") %>%
  select(Family, count.x, count.y) %>%
  rename(Ath = count.x,
         Bst = count.y) %>%
  full_join(sl, by = "Family") %>%
  select(Family, Ath, Bst, count) %>%
  rename(Sly = count) %>%
  full_join(zm, by = "Family") %>%
  select(Family, Ath, Bst, Sly, count) %>%
  rename(Zma = count) %>%
  full_join(os, by = "Family") %>%
  select(Family, Ath, Bst, Sly, Zma, count) %>%
  rename(Osa = count) %>%
  full_join(atr, by = "Family") %>%
  select(Family, Ath, Bst, Sly, Zma, Osa, count) %>%
  rename(Atr = count) %>% distinct(Family, .keep_all = T) %>%
  replace(is.na(.), 0) %>% # Angiosperms, debugging
  mutate(Average = mean(Ath,Bst, Sly, Zma, Osa,Atr)) %>%
  rename(Angiosperms = Average) %>%
  select(Family, Angiosperms)

# Gymnosperm TAP averages
gymnosperm_taps_average <- full_join(pt, pa, by = "Family") %>%
  select(Family, count.x, count.y) %>%
  rename(Pta = count.x,
         Pab = count.y) %>% 
  distinct(Family, .keep_all = T) %>%
  replace(is.na(.), 0) %>% # Angiosperms, debugging
  mutate(Average = mean(Pta + Pab)) %>%
  rename(Gymnosperms = Average) %>%
  select(Family, Gymnosperms)

# Fern TAP averages
fern_taps_average <- full_join(af, sc, by = "Family") %>%
  select(Family, count.x, count.y) %>%
  rename(Afi = count.x,
         Scu = count.y) %>% 
  distinct(Family, .keep_all = T) %>%
  replace(is.na(.), 0) %>% # Angiosperms, debugging
  mutate(Average = mean(Afi, Scu)) %>%
  rename(Ferns = Average) %>%
  select(Family, Ferns)

# Lycophyte TAP averages
lyco_taps_average <- sm %>%
  rename(Ppa = count) %>%
  group_by(Family) %>%
  select(Family, Ppa) %>%
  mutate(Lycophyte = mean(Ppa)) %>%
  select(Family, Lycophyte)


# Bryophytes TAP averages
bryo_taps_average <- full_join(ppa, mpo, by = "Family") %>%
  select(Family, count.x, count.y) %>%
  rename(Ppa = count.x,
         Mpo = count.y) %>% 
  full_join(aag, by = "Family") %>% 
  select (Family, Ppa, Mpo, count) %>%
  rename(Aag = count) %>%
  distinct(Family, .keep_all = T) %>%
  replace(is.na(.), 0) %>% # Angiosperms, debugging
  mutate(Bryophytes = mean(Ppa, Mpo, Aag)) %>%
  select(Family, Bryophytes)

# Charophytes TAP averages
charophyte_taps_average <- full_join(smu, men, by = "Family") %>%
  select(Family, count.x, count.y) %>%
  rename(Smu = count.x,
         Men = count.y) %>% 
  full_join(cbr, by = "Family") %>% 
  select (Family, Smu, Men, count) %>%
  rename(Cbr = count) %>%
  distinct(Family, .keep_all = T) %>%
  replace(is.na(.), 0) %>% #Angiosperms, debugging
  mutate(Charophytes = mean(Smu, Men, Cbr)) %>%
  select(Family, Charophytes)

# Chlorophytes TAP averages
chlorophyte_taps_average <- full_join(vca, cre, by = "Family") %>%
  select(Family, count.x, count.y) %>%
  rename(Vca = count.x,
         Cre = count.y) %>% 
  select(Family, Vca, Cre) %>%
  distinct(Family, .keep_all = T) %>%
  replace(is.na(.), 0) %>% 
  mutate(Chlorophytes = mean(Vca, Cre)) %>%
  select(Family, Chlorophytes)

# Creating matrix for averages
viridiplant_taps_means <-full_join(angiosperm_taps_averages,
                                   gymnosperm_taps_average,
                                   by = "Family") %>%
  full_join(fern_taps_average, by = "Family") %>%
  full_join(lyco_taps_average, by = "Family") %>%
  full_join(bryo_taps_average, by = "Family") %>%
  full_join(charophyte_taps_average, by = "Family") %>%
  full_join(chlorophyte_taps_average, by = "Family") %>%
  distinct(Family, .keep_all = T) %>%
  na_if(0) %>%
  distinct(Family, .keep_all = T) %>%
  group_by("Family") %>%
  slice(1:60, .preserve = T) %>%
  ungroup() 

# Building matrices - viridiplantae taps #####
write_csv(viridiplant_taps_means, "viridiplant_taps_means.csv", col_names = T)
viridiplantae_taps_means_mat <- as.matrix(read.csv("viridiplant_taps_means.csv",
                                                   sep = ",",
                                                   row.names = "Family"))
viridiplantae_mean_matrix <- log10(viridiplantae_taps_means_mat)
viridiplantae_taps_means_heatmap<- pheatmap(
  viridiplantae_mVean_matrix,
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
  main = "LxCxE TAPs average family distribution by lineage")
viridiplantae_taps_means_heatmap
pdf("viridiplantae_taps_means_heatmap.pdf")
print(viridiplantae_taps_means_heatmap) # Plot 1 --> in the first page of PDF
dev.off()


