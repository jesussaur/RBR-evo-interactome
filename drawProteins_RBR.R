# Loading drawProteins
library(drawProteins)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(viridis)
# Getting Arabidopsis RBR data from uniprot
drawProteins::get_features("Q9LKZ3") -> RBR_json
# Converting to json file to dataframe
drawProteins::feature_to_dataframe(RBR_json) -> RBR_data
rbr_tibble <- as_tibble(RBR_data) # Tibble with RBR features data

# Can be managed as a multidataframe for several proteins in one file
prot_data <- drawProteins::get_features("Q04206 Q01201 Q04864 P19838 Q00653")
drawProteins::feature_to_dataframe(prot_data) -> prot_df
prot_tibble <- as_tibble(prot_df)

p <- draw_canvas(prot_tibble)
p <- draw_chains(p, prot_tibble)
p <- draw_domains(p, prot_tibble, fill)




# Plotting several features 
p <- draw_canvas(rbr_tibble)
p <- draw_chains(p, rbr_tibble)
p <- draw_regions(p, rbr_tibble)

test <- read.csv("RBR_drawProteins/Mpo_RBR_features.csv")
t1 <- draw_canvas(test)
t1 <- draw_chains(t1, test)
t1 <- draw_domains(t1, test)



# Plotting RBR orthologues
RBR_df <- read_csv("Workdir/viridiplantae_RBR-features.csv") # Getting RBR data
RBR_prots <- draw_canvas(RBR_df) # Drawing canvas
RBR_prots <- draw_chains(RBR_prots, RBR_df, fill = "#287D8E", outline = "black", size = 0.5)
RBR_prots <- draw_domains(RBR_prots, RBR_df, label_domains = FALSE)

draw_phosphosites <- function (p, data = data, size = 2, fill = "yellow", show.legend = FALSE) 
{
  begin = end = description = NULL
  p <- p + ggplot2::geom_point(data = drawProteins::phospho_site_info(data), 
                               ggplot2::aes(x = begin, y = order + 0.35), shape = 21, 
                               colour = "black", fill = fill, size = size, show.legend = show.legend)
  return(p)
}

RBR_prots1 <- draw_phosphosites(RBR_prots, RBR_df, fill = "#B8DE29",size = 1.3) + theme_void()
RBR_prots2 <- draw_regions(RBR_prots1, RBR_df)

pdf("Viridplantae_RBR-features.pdf")
print(RBR_prots2) # Plot 1 --> in the first page of PDF
dev.off()




