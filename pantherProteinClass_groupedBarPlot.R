# April 14th, 2021. Jesus Leon 
# This script is for plotting PANTHER protein class data per groups
# Angiosperms, 

library(tidyverse)
library(viridis)
setwd("Workdir/")
# Importing and processing data per group
# Angiosperm Protein Class ####
angiosperm_lxcxe_pc <- read_csv("angiospermsLXCXE_ProteinClass.csv") %>%
  dplyr::select(gene_ID,
                Protein_class) %>%
  group_by(Protein_class) %>%
  mutate(Clade = "Angiosperms", count = n()) %>%
  dplyr::select(Protein_class, Clade, count) %>%
  distinct() %>%
  drop_na() %>%
  ungroup() %>%
  top_n(25) %>%
  arrange(desc(count))

# Tracheophyta Protein Class ####
tracheophyte_lxcxe_pc <- read_csv("tracheophyta-LXCXE_Athaliana.csv") %>%
  dplyr::select(gene_ID,
                Protein_class) %>%
  group_by(Protein_class) %>%
  mutate(Clade = "Tracheophyta", count = n()) %>%
  dplyr::select(Protein_class, Clade, count) %>%
  distinct() %>%
  drop_na() %>%
  ungroup() %>%
  top_n(25) %>%
  arrange(desc(count))

# Embryophyta Protein Class ####
embryophyte_lxcxe_pc <- read_csv("embryophytaLXCXE_ProteinClass.csv") %>%
  dplyr::select(gene_ID,
                Protein_class) %>%
  group_by(Protein_class) %>%
  mutate(Clade = "Embryophyta", count = n()) %>%
  dplyr::select(Protein_class, Clade, count) %>%
  distinct() %>%
  drop_na() %>%
  ungroup() %>%
  top_n(25) %>%
  arrange(desc(count))

# Streptophyta Protein Class ####
streptophyte_lxcxe_pc <- read_csv("streptophytaLXCXE_ProteinClass.csv") %>%
  dplyr::select(gene_ID,
                Protein_class) %>%
  group_by(Protein_class) %>%
  mutate(Clade = "Streptophyta", count = n()) %>%
  dplyr::select(Protein_class, Clade, count) %>%
  distinct() %>%
  drop_na() %>%
ungroup() %>%
  top_n(25) %>%
  arrange(desc(count))

# Charophyta Protein Class ####
charophyta_lxcxe_pc <- read_csv("charophytaLXCXE_ProteinClass.csv") %>%
  dplyr::select(gene_ID,
                Protein_class) %>%
  group_by(Protein_class) %>%
  mutate(Clade = "Charophyta", count = n()) %>%
  dplyr::select(Protein_class, Clade, count) %>%
  distinct() %>%
  drop_na() %>%
  ungroup() %>%
  top_n(25) %>%
  arrange(desc(count))

# Chlorophyta Protein Class ####
chlorophyte_lxcxe_pc <- read_csv("chlorophytaLXCXE_ProteinClass.csv") %>%
  dplyr::select(gene_ID,
                Protein_class) %>%
  group_by(Protein_class) %>%
  mutate(Clade = "Chlorophyta", count = n()) %>%
  dplyr::select(Protein_class, Clade, count) %>%
  distinct() %>%
  drop_na() %>%
  ungroup() %>%
  top_n(25) %>%
  arrange(desc(count))


viridiplantae_proteinclass_tibble <- tibble(clade = c(angiosperm_lxcxe_pc$Clade,
                                                      tracheophyte_lxcxe_pc$Clade,
                                                      embryophyte_lxcxe_pc$Clade,
                                                      streptophyte_lxcxe_pc$Clade,
                                                      charophyta_lxcxe_pc$Clade,
                                                      chlorophyte_lxcxe_pc$Clade),
                                            protein_class = c(angiosperm_lxcxe_pc$Protein_class,
                                                              tracheophyte_lxcxe_pc$Protein_class,
                                                              embryophyte_lxcxe_pc$Protein_class,
                                                              streptophyte_lxcxe_pc$Protein_class,
                                                              charophyta_lxcxe_pc$Protein_class,
                                                              chlorophyte_lxcxe_pc$Protein_class),
                                            value = c(angiosperm_lxcxe_pc$count,
                                                      tracheophyte_lxcxe_pc$count,
                                                      embryophyte_lxcxe_pc$count,
                                                      streptophyte_lxcxe_pc$count,
                                                      charophyta_lxcxe_pc$count,
                                                      chlorophyte_lxcxe_pc$count)) %>%
  arrange(desc(value))
  
test <-transform(viridiplantae_proteinclass_tibble,
          protein_class = reorder(protein_class, value))

# Grouped barplot test
ggplot(test, aes(fill=clade,
                 y=value,
                 x=protein_class,
                 width = 0.75)) + 
  geom_bar(position="dodge", stat="identity", width = 5) + 
  scale_fill_viridis(discrete = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7)) +
  coord_flip() + theme_minimal()

##### Plotting data for NUCLEAr LxCxE proteins #####

# Importing Protein Class data for Nuclear LxCxE proteins
# Angiosperm Nuclear LxCxE Protein Class ####
angiosperm_nuclear_pc <- read_csv("angiosperm-LXCXE_NuclearProteinClass.csv") %>%
  dplyr::select(gene_ID,
                Protein_class) %>%
  group_by(Protein_class) %>%
  mutate(Clade = "Angiosperms", count = n()) %>%
  dplyr::select(Protein_class, Clade, count) %>%
  distinct() %>%
  drop_na() %>%
  ungroup() %>%
  top_n(25) %>%
  arrange(desc(count))

# Tracheophyta Nuclear LxCxE Protein Class ####
tracheophyte_nuclear_pc <- read_csv("tracheophyta-LXCXE_NuclearProteinClass.csv") %>%
  dplyr::select(gene_ID,
                Protein_class) %>%
  group_by(Protein_class) %>%
  mutate(Clade = "Tracheophyta", count = n()) %>%
  dplyr::select(Protein_class, Clade, count) %>%
  distinct() %>%
  drop_na() %>%
  ungroup() %>%
  top_n(25) %>%
  arrange(desc(count))

# Spermatophyta Protein Class ####
spermatophyte_nuclear_pc <- read_csv("spermatophyta-LXCXE_NuclearProteinClass.csv") %>%
  dplyr::select(gene_ID,
                Protein_class) %>%
  group_by(Protein_class) %>%
  mutate(Clade = "Spermatophyta", count = n()) %>%
  dplyr::select(Protein_class, Clade, count) %>%
  distinct() %>%
  drop_na() %>%
  ungroup() %>%
  top_n(25) %>%
  arrange(desc(count))

# Embryophyta Protein Class ####
embryophyte_nuclear_pc <- read_csv("embryophyta-LXCXE_NuclearProteinClass.csv") %>%
  dplyr::select(gene_ID,
                Protein_class) %>%
  group_by(Protein_class) %>%
  mutate(Clade = "Embryophyta", count = n()) %>%
  dplyr::select(Protein_class, Clade, count) %>%
  distinct() %>%
  drop_na() %>%
  ungroup() %>%
  top_n(25) %>%
  arrange(desc(count))

# Streptophyta Protein Class ####
streptophyte_nuclear_pc <- read_csv("streptophyta-LXCXE_NuclearProteinClass.csv") %>%
  dplyr::select(gene_ID,
                Protein_class) %>%
  group_by(Protein_class) %>%
  mutate(Clade = "Streptophyta", count = n()) %>%
  dplyr::select(Protein_class, Clade, count) %>%
  distinct() %>%
  drop_na() %>%
  ungroup() %>%
  top_n(25) %>%
  arrange(desc(count))

# Charophyta Protein Class ####
charophyte_nuclear_pc <- read_csv("charophyta-LXCXE_NuclearProteinClass.csv") %>%
  dplyr::select(gene_ID,
                Protein_class) %>%
  group_by(Protein_class) %>%
  mutate(Clade = "Charophyta", count = n()) %>%
  dplyr::select(Protein_class, Clade, count) %>%
  distinct() %>%
  drop_na() %>%
  ungroup() %>%
  top_n(25) %>%
  arrange(desc(count))

# Chlorophyta Protein Class ####
chlorophyte_nuclear_pc <- read_csv("chlorophyta-LXCXE_NuclearProteinClass.csv") %>%
  dplyr::select(gene_ID,
                Protein_class) %>%
  group_by(Protein_class) %>%
  mutate(Clade = "Chlorophyta", count = n()) %>%
  dplyr::select(Protein_class, Clade, count) %>%
  distinct() %>%
  drop_na() %>%
  ungroup() %>%
  top_n(25) %>%
  arrange(desc(count))


viridiplantae_nuclear_tibble <- tibble(clade = c(angiosperm_nuclear_pc$Clade,
                                                      tracheophyte_nuclear_pc$Clade,
                                                      spermatophyte_nuclear_pc$Clade,
                                                      embryophyte_nuclear_pc$Clade,
                                                      streptophyte_nuclear_pc$Clade,
                                                      charophyte_nuclear_pc$Clade,
                                                      chlorophyte_nuclear_pc$Clade),
                                            protein_class = c(angiosperm_nuclear_pc$Protein_class,
                                                              tracheophyte_nuclear_pc$Protein_class,
                                                              spermatophyte_nuclear_pc$Protein_class,
                                                              embryophyte_nuclear_pc$Protein_class,
                                                              streptophyte_nuclear_pc$Protein_class,
                                                              charophyte_nuclear_pc$Protein_class,
                                                              chlorophyte_nuclear_pc$Protein_class),
                                            value = c(angiosperm_nuclear_pc$count,
                                                      tracheophyte_nuclear_pc$count,
                                                      spermatophyte_nuclear_pc$count,
                                                      embryophyte_nuclear_pc$count,
                                                      streptophyte_nuclear_pc$count,
                                                      charophyte_nuclear_pc$count,
                                                      chlorophyte_nuclear_pc$count)) %>%
  arrange(desc(value))

test_nuclear <-transform(viridiplantae_nuclear_tibble,
                 protein_class = reorder(protein_class, value))

p2 <- ggplot(df, aes(x = reorder(Category, -Count), y = Count)) +
  geom_bar(stat = "identity")
# Reordering legend
test_nuclear$clade <- factor(test_nuclear$clade, levels = c("Angiosperms",
                                                            "Spermatophyta",
                                                            "Tracheophyta",
                                                            "Embryophyta",
                                                            "Streptophyta",
                                                            "Charophyta",
                                                            "Chlorophyta"))
test_nuclear$clade <- factor(viridiplantae_nuclear_tibble$clade, levels = c("Angiosperms",
                                                            "Spermatophyta",
                                                            "Tracheophyta",
                                                            "Embryophyta",
                                                            "Streptophyta",
                                                            "Charophyta",
                                                            "Chlorophyta"))

# Grouped barplot test for nuclear proteins
nuclear_lxcxe_proteinClass <- ggplot(test_nuclear, aes(x = reorder(protein_class, -value),
                                                       y = value,
                                                       fill=clade,
                                                       width = 0.8)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust = 0.95)) + 
  geom_bar(position="dodge", stat="identity", width = 2) + 
  scale_fill_viridis(discrete = TRUE) + 
  labs(y = "Number of Protein Class hits") +
  theme(plot.margin = unit(c(3, 3, 3, 3), "cm"))
nuclear_lxcxe_proteinClass
pdf("PANTHERProteinClass_NuclearLxCxE.pdf")
print(nuclear_lxcxe_proteinClass) 
dev.off()


#________






# Grouped barplot test for nuclear proteins
nuclear_lxcxe_proteinClass <- ggplot(test_nuclear, aes(fill=clade,
                 y=value,
                 x=protein_class,
                 width = 0.8)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.95)) + 
  geom_bar(position="dodge", stat="identity", width = 2) + 
  scale_fill_viridis(discrete = TRUE) + 
  labs(y = "Number of Protein Class hits",
       x = "Protein Class")
nuclear_lxcxe_proteinClass
pdf("PANTHERProteinClass_NuclearLxCxE.pdf")
print(nuclear_lxcxe_proteinClass) 
dev.off()

ggplot(data, aes(x, y, fill = y)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90))


#### EXAMPLE GGPLOT GROUP BARPLOT
# create a dataset
specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
value <- abs(rnorm(12 , 0 , 15))
data <- data.frame(specie,condition,value)

# Grouped
ggplot(data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="dodge", stat="identity")
