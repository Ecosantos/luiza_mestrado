library(vegan)
library(ggplot2)
terra::rast(system.file("ex/elev.tif", package="terra"))
library(dplyr)

# SPECIES ACCUMULATION

mat_uni2 <- mat_uni1$PAM 
mat_uni2_clean <- mat_uni2[, sapply(mat_uni2, is.numeric)]
accumula <- specaccum(mat_uni2_clean)
plot(accumula)
save(accumula, file = "Data/accumula.RData")
save(mat_uni2_clean,file = "Data/mat_uni2_clean.RData")

# Plot Sites
plot_data <- data.frame("Locais" = c(0, accumula$sites),
                        "Riqueza" = c(0, accumula$richness),
                        "lower" = c(0, accumula$richness - accumula$sd),
                        "upper" = c(0, accumula$richness + accumula$sd))
g <- ggplot(plot_data, aes(x = Locais, y = Riqueza)) +
  geom_point(color = "blue", size = 2) +
  geom_line(color = "blue", lwd = 2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              linetype=2, alpha=0.3, fill = "orange") +
  ylab("Riqueza acumulada") +
  theme_classic() +
  theme(text = element_text(size = 16))
g

ggsave("Figures/Rarefac.png")

# RAREFACTION CURVE
library(forcats)
library(iNEXT)
library(vegan)
library(tidyverse)


site_totals <- rowSums(mat_uni2_clean)
mat_uni2_final <- mat_uni2_clean[site_totals > 0, ]

individual_plots <- list()

area_names <- rownames(mat_uni2_final)
area_names[grepl("[^A-Za-z0-9_-]", area_names)]

# Loop iNEXT for each area
for(i in 1:length(area_names)) {
  single_area_data <- mat_uni2_final[i, , drop = FALSE]
  out_single <- iNEXT(t(single_area_data), q = 0, 
                      datatype = "abundance")
  save(out_single, file = "Data/out_single.RData")
  # Plot
  individual_plots[[i]] <- ggiNEXT(out_single, type = 1) +
    theme_bw() +
    labs(title = UCs$nome_uc[i]) +
    xlab("Número de indivíduos") + 
    ylab("Riqueza de espécies") +
    theme(legend.position = "none")
  
  # Save
  clean_name <- gsub("[^A-Za-z0-9_-]", "_", UCs$nome_uc[i])
  filename <- paste0("Figures/", clean_name, "_rarefaction.png")
  ggsave(filename, individual_plots[[i]], 
         width = 8, height = 6, dpi = 300)
  print(paste("Created and saved:", filename))
}

# Name the plots
names(individual_plots) <- area_names
individual_plots[i]


