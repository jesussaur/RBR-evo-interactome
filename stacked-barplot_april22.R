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
readAAStringSet("Crichardii_LXCXE.fasta")
readAAStringSet("Proteomes/Smoellendorffii_LXCXE.fasta")
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
ath <- read_csv("Athaliana_LXCXE.fasta.csv")
bst <- read_csv("Bstricta_LXCXE.fasta.csv")
sly <- read_csv("Slycopersicum_LXCXE.fasta.csv")
zma <- read_csv("Zmays_LXCXE.fasta.csv")
svi <- read_csv("Sviridis_LXCXE.fasta.csv")
osa <- read_csv("Osativa_LXCXE.fasta.csv")
atr <- read_csv("Atrichopoda_LXCXE.fasta.csv")
unique(ath$query) #197
unique(bst$query) #162
unique(sly$query) #115
unique(zma$query) #162
unique(svi$query) #221
unique(osa$query) #150
unique(atr$query) #75

#gymnosperms
pta <- read_csv("Ptaeda_LXCXE.csv")
pab <- read_csv("Pabies_LXCXE.fasta.csv")
gmo <- read_csv("Gmontanum_LXCXE.fasta.csv")
unique(pta$query) #178
unique(pab$query) #109
unique(gmo$query) #75

#ferns
afi <- read_csv("Afilliculoides_LXCXE.fasta.csv")
scu <- read_csv("Scucullata_LXCXE.fasta.csv")
cri <- read_csv("Crichardii_LXCXE.fasta.csv")
unique(afi$query) #112
unique(scu$query) #83
unique(cri$query) #182

#lycophyte
smo <- read_csv("Smoellendorffii_LXCXE.fasta.csv")
unique(smo$query) #55

#bryophytes
ppa <- read_csv("Ppatens_LXCXE.fasta.csv")
sfa <- read_csv("Sfallax_LXCXE.fasta.csv")
mpo <- read_csv("Mpolymorpha_LXCXE.fasta.csv")
aag <- read_csv("Aagrestis_LXCXE.fasta.csv")
unique(ppa$query) #245
unique (sfa$query) #117
unique(mpo$query) #71
unique(aag$query) #99

#zygnematales
pma <- read_csv("Pmargaritaceum_LXCXE.fasta.csv")
men <- read_csv("Mendlicheranium_LXCXE.fasta.csv")
smu <- read_csv("Smuscicola_LXCXE.fasta.csv")
unique(pma$query) #65
unique(men$query) #24
unique(smu$query) #87

#charales
cbr <- read_csv("Cbraunii_LXCXE.fasta.csv")
unique(cbr$query) #80

#klebsormidiales
kni <- read_csv("Knitens_LXCXE.fasta.csv")
unique(kni$query) #48

#chlorokybales
cat <- read_csv("Catmophyticus_LXCXE.fasta.csv")
unique(cat$query) #20

#mesostigmatales
mvi <- read_csv("MvirideLIANG_LXCXE.fasta.csv")
unique(mvi$query) #56

#chlorophytees
cre<- read_csv("Creinhardtii_LXCXE.fasta.csv")
vca <- read_csv("Vcarteri_LXCXE.fasta.csv")
olu <- read_csv("Olucimarinus_LXCXE.fasta.csv")
unique(cre$query) #32
unique(vca$query) #35
unique(olu$query) #14

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
condition <- rep(c("LxCxE-protein" , "LxCxE-TAP") , 28)
lxcxe_taps_values <- c(1308, 197, # Ath
                       1102, 162, # Bst
                       1053, 115, # Sly
                       2012, 162, # Zma
                       2084, 221, #Svi
                       1566, 150,  # Osa
                       794, 75,   # Atr
                       2020, 178, # Pta
                       1478, 109, # Pab
                       825, 75, # Gmo
                       1658, 112, # Afi
                       1067, 83,  # Scu
                       4242, 182, # Cri
                       788, 55,   # Smo
                       3407, 245,  # Ppa
                       1023, 117, # Sfa
                       852, 71,   # Mpo
                       2163, 99,  # Aag
                       1129, 65, # Pma
                       345, 24,   # Men
                       1014, 87,  # Smu
                       1393, 80,  # Cbr
                       682, 48,   # Kni
                       309, 20,   # Cat
                       1108, 56,   # Mvi
                       653, 35,   # Vca
                       701, 32, # Cre
                       216, 14) # Olu

# Creating tibble to plot stacked bar plots
tap_data <- tibble(Species = species,
                   Type = condition,
                   Frequency = lxcxe_taps_values) %>%
  mutate(name = fct_relevel(Species, 
                            "Olu", "Cre", "Vca", "Mvi", 
                            "Cat", "Kni", "Cbr", 
                            "Smu","Men", "Pma",
                            "Aag", "Mpo", "Sfa",
                            "Ppa", "Smo", "Cri", "Scu",
                            "Afi", "Gmo", "Pab",
                            "Pta", "Atr", "Osa",
                            "Svi", "Zma", "Sly", "Bst", "Ath"))

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

ggplot(df,aes(x,y))+geom_bar(stat="identity")+theme(axis.text.x=element_text(face=c("italic","italic","italic","italic")))

p2 <- p1 + scale_color_manual(values = c("#588033", "#1088AA"))



