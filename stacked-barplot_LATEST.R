library(ggplot2)
library(tidyverse)
library(viridis)
# Creating data frame 
## LxCxE proteins and TAPs detected on each lineage
## TAPS = Transcription associated proteins

setwd("Workdir/")

df_lxcxe_taps <- data.frame(x= c("AthLCE",
                                 "BstLCE",
                                 "SlyLCE", 
                                 "ZmaLCE", 
                                 "OsaLCE",
                                 "AtrLCE",
                                 "PtaLCE", 
                                 "PabLCE", 
                                 "AfiLCE", 
                                 "ScuLCE", 
                                 "SmoLCE",
                                 "PpaLCE", 
                                 "MpoLCE", 
                                 "AagLCE", 
                                 "SmuLCE", 
                                 "MenLCE", 
                                 "ChoLCE", 
                                 "CbrLCE", 
                                 "KniLCE",
                                 "CatLCE", 
                                 "MviLCE", 
                                 "VcaLCE", 
                                 "CreLCE"),
                            y = c(1308,  # Ath
                                  1102,  # Bst
                                  1053,  # Sly
                                  2012,  # Zma
                                  1130,  # Osa
                                  794,    # Atr
                                  2020, # Pta
                                  1478,  # Pab
                                  1658,  # Afi
                                  1067,   # Scu
                                  788,   # Smo
                                  3407,   # Ppa
                                  871,   # Mpo
                                  2163,   # Aag
                                  1014,  # Smu
                                  345,    # Men
                                  0,     # Cho
                                  1393,  # Cbr
                                  682,    # Kni
                                  309,    # Cat
                                  271,    # Mvi
                                  653,   # Vca
                                  701)) # Cre

write_csv(df_lxcxe_taps, path = "df_lxcxe_taps.csv")
lxcxe_taps <- read_csv("bioinfo")
select(x, y, group)
# Consider grouping "LxCxE"/"TAPs" in a different column (groups)
lxcxe_taps$x <- factor(lxcxe_taps$x,
                       levels = lxcxe_taps$x)
colnames(lxcxe_taps) <- c("Species", "Count", "Group")
# Horizontal version
p <- ggplot(lxcxe_taps, aes(x=Species, y=Count,
                            fill = Group,
                            colour = Group,
                            label = paste0(round(Count, 0)))) +
  geom_segment(aes(x=Species, xend=Species, y=0, yend=Count), color="skyblue") +
  geom_point(size=3, alpha=0.6) +
  theme_light() + 
  coord_flip() +  
  theme_vi
theme_vir
ggplot(OH_top10, aes(percollege, county, label = round(percollege, 1))) +
  geom_segment(aes(x = 0, y = county, xend = percollege, yend = county), color = "grey50") +
  geom_point() +
  geom_text(nudge_x = 1.5)

pdf("ggplot.pdf")
print(p) # Plot 1 --> in the first page of PDF
dev.off()

# STACKED BARPLOTS ####
library(Biostrings)
readAAStringSet("Proteomes/Athaliana_LXCXE.fasta")
readAAStringSet("Proteomes/Bstricta_LXCXE.fasta")
readAAStringSet("Proteomes/Slycopersicum_LXCXE.fasta")
readAAStringSet("Proteomes/Zmays_LXCXE.fasta")
readAAStringSet("Proteomes/Sviridis_LXCXE.fasta")
readAAStringSet("Proteomes/Osativa_LXCXE.fasta")
readAAStringSet("Proteomes/Atrichopoda_LXCXE.fasta")
readAAStringSet("Proteomes/Ptaeda_LXCXE.fasta")
readAAStringSet("Proteomes/Pabies_LXCXE.fasta")
readAAStringSet("Proteomes/Gmontanum_LXCXE.fasta")
readAAStringSet("Proteomes/Afilliculoides_LXCXE.fasta")
readAAStringSet("Proteomes/Scucullata_LXCXE.fasta")
readAAStringSet("Proteomes/Crichardii_LXCXE.fasta")
readAAStringSet("Proteomes/Smoellendorffii_LXCXE.fasta")
readAAStringSet("Proteomes/Itaiwanensis_LxCxE.fasta")
readAAStringSet("Proteomes/Ppatens_LXCXE.fasta")
readAAStringSet("Proteomes/Sfallax_LXCXE.fasta")
readAAStringSet("Proteomes/Mpolymorpha_LXCXE.fasta")
readAAStringSet("Proteomes/Aagrestis_LXCXE.fasta")
readAAStringSet("Proteomes/Pmargaritaceum_LXCXE.fasta")
readAAStringSet("Proteomes/Mendlicheranium_LXCXE.fasta")
readAAStringSet("Proteomes/Smuscicola_LXCXE.fasta")
readAAStringSet("Proteomes/Cbraunii_LXCXE.fasta")
readAAStringSet("Proteomes/Knitens_LXCXE.fasta")
readAAStringSet("Proteomes/Catmophyticus_LXCXE.fasta")
readAAStringSet("Proteomes/MvirideLIANG_LXCXE.fasta")
readAAStringSet("Proteomes/Vcarteri_LXCXE.fasta")
readAAStringSet("Proteomes/Creinhardtii_LXCXE.fasta")
readAAStringSet("Olucimarinus_LXCXE.fasta")

# Quantifying TAPs with unique(species$query) examples

# angiosperms
ath <- read_csv("AthalianaLXCXE_TAPs_final.csv")
bst <- read_csv("BstrictaLXCXE_TAPs_final.csv")
sly <- read_csv("SlycopersicumLXCXE_TAPs_final.csv")
zma <- read_csv("ZmaysLXCXE_TAPs_final.csv")
svi <- read_csv("SviridisLXCXE_TAPs_final.csv")
osa <- read_csv("OsativaLXCXE_TAPs_final.csv")
atr <- read_csv("AtrichopodaLXCXE_TAPs_final.csv")
unique(ath$query) #179
unique(bst$query) #145
unique(sly$query) #107
unique(zma$query) #153
unique(svi$query) #212
unique(osa$query) #137
unique(atr$query) #70

#gymnosperms
pta <- read_csv("PtaedaLXCXE_TAPs_final.csv")
pab <- read_csv("PabiesLXCXE_TAPs_final.csv")
gmo <- read_csv("GmontanumLXCXE_TAPs_final.csv")
unique(pta$query) #172
unique(pab$query) #104
unique(gmo$query) #70

#ferns
afi <- read_csv("AfiliculoidesLXCXE_TAPs_final.csv")
scu <- read_csv("ScucullataLXCXE_TAPs_final.csv")
cri <- read_csv("CrichardiiLXCXE_TAPs_final.csv")
unique(afi$query) #110
unique(scu$query) #81
unique(cri$query) #170

#lycophyte
smo <- read_csv("SmoellendorffiiLXCXE_TAPs_final.csv")
ita <- read_csv("Itaiwanensis_LxCxE.fasta.csv")
unique(smo$query) #49
unique(ita$query) #101

#bryophytes
ppa <- read_csv("PpatensLXCXE_TAPs_final.csv")
sfa <- read_csv("SfallaxLXCXE_TAPs_final.csv")
mpo <- read_csv("MpolymorphaLXCXE_TAPs_final.csv")
aag <- read_csv("AagrestisLXCXE_TAPs_final.csv")
unique(ppa$query) #232
unique (sfa$query) #112
unique(mpo$query) #68
unique(aag$query) #94

#zygnematales
pma <- read_csv("PmargaritaceumLXCXE_TAPs_final.csv")
men <- read_csv("MendlicheraniumLXCXE_TAPs_final.csv")
smu <- read_csv("SmuscicolaLXCXE_TAPs_final.csv")
unique(pma$query) #63
unique(men$query) #23
unique(smu$query) #85

#charales
cbr <- read_csv("CbrauniiLXCXE_TAPs_final.csv")
unique(cbr$query) #73

#klebsormidiales
kni <- read_csv("KnitensLXCXE_TAPs_final.csv")
unique(kni$query) #45

#chlorokybales
cat <- read_csv("CatmophyticusLXCXE_TAPs_final.csv")
unique(cat$query) #17

#mesostigmatales
mvi <- read_csv("MvirideLXCXE_TAPs_final.csv")
unique(mvi$query) #50

#chlorophytees
cre<- read_csv("CreinhardtiiLXCXE_TAPs_final.csv")
vca <- read_csv("VcarteriLXCXE_TAPs_final.csv")
olu <- read_csv("OlucimarinusLXCXE_TAPs_final.csv")
unique(cre$query) #30
unique(vca$query) #32
unique(olu$query) #11

# create a dataset
species <- c(rep("Ath", 2),
             rep("Bst", 2),
             rep("Sly", 2),
             rep("Zma", 2),
             rep("Svi", 2),
             rep("Osa", 2),
             rep("Atr", 2),
             rep("Pta", 2),
             rep("Pab", 2),
             rep("Gmo", 2),
             rep("Afi", 2),
             rep("Scu", 2),
             rep("Cri", 2),
             rep("Smo", 2),
             rep("Ita", 2),
             rep("Ppa", 2),
             rep("Sfa", 2),
             rep("Mpo", 2),
             rep("Aag", 2),
             rep("Pma", 2),
             rep("Men", 2),
             rep("Smu", 2),
             rep("Cbr", 2),
             rep("Kni", 2),
             rep("Cat", 2),
             rep("Mvi", 2),
             rep("Vca", 2),
             rep("Cre", 2),
             rep("Olu", 2))
condition <- rep(c("LxCxE-protein" , "LxCxE-TAP") , 29)
lxcxe_taps_values <- c(1308, 179, # Ath
                       1102, 145, # Bst
                       1053, 107, # Sly
                       2012, 153, # Zma
                       2084, 212, #Svi
                       1566, 137,  # Osa
                       794, 75,   # Atr
                       2020, 172, # Pta
                       1478, 104, # Pab
                       825, 70, # Gmo
                       1658, 110, # Afi
                       1067, 81,  # Scu
                       4242, 170, # Cri
                       788, 49,   # Smo
                       2024, 101, # Ita
                       2163, 94,  # Aag
                       3407, 232,  # Ppa
                       1023, 112, # Sfa
                       852, 68,   # Mpo
                       1129, 63, # Pma
                       345, 23,   # Men
                       1014, 85,  # Smu
                       1393, 73,  # Cbr
                       682, 45,   # Kni
                       309, 17,   # Cat
                       1108, 50,   # Mvi
                       653, 32,   # Vca
                       701, 30, # Cre
                       216, 11) # Olu

# Creating tibble to plot stacked bar plots
tap_data <- tibble(Species = species,
                   Type = condition,
                   Frequency = lxcxe_taps_values) %>%
  mutate(name = fct_relevel(Species, 
                            "Olu", "Cre", "Vca", "Mvi", 
                            "Cat", "Kni", "Cbr", 
                            "Smu","Men", "Pma",
                            "Mpo", "Sfa", "Ppa",
                            "Aag", "Ita", "Smo", 
                            "Cri", "Scu", "Afi", 
                            "Gmo", "Pab", "Pta", 
                            "Atr", "Osa", "Svi", 
                            "Zma", "Sly", "Bst", "Ath"))

# plotting stacked bar plot
p1 <- ggplot(tap_data, aes(fill=Type, y=Frequency, x=name)) + 
  geom_bar(position="stack", stat="identity") + coord_flip()

p2<-p1 + scale_fill_viridis(discrete = T) +
  ggtitle("LxCxE/TAP") + xlab("Frequency") + ylab ("Species") + 
  theme_minimal() + theme(
    axis.text.y = element_text(face = "italic", size = 11.5),
    axis.text.x = element_text(face = "bold")) + 
  scale_fill_manual(values = c("gray77", "#1088aa")) + 
  geom_text(size = 2.5, position = position_stack(vjust = 1),
            colour = "009193", label = lxcxe_taps_values, fontface= "bold")
pdf("p2.pdf")
print(p2) # Plot 1 --> in the first page of PDF
dev.off()



#ggplot(df,aes(x,y))+geom_bar(stat="identity")+theme(axis.text.x=element_text(face=c("italic","italic","italic","italic")))

#p2 <- p1 + scale_color_manual(values = c("#588033", "#1088AA"))



