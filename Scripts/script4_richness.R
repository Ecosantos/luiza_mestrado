library(forcats)
library(iNEXT)
library(vegan)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggpubr)
library(tidyr)
library(sf)
library(viridis)
library(ggridges)

# ESTIMATED AND OBSERVED RICHNESS

# Open shapefiles
UCs <- terra::vect("Data/shp_cnuc_2024_02/cnuc_2024_02.shp")
UCs <- UCs[UCs$esfera == c("Federal", "Estadual"), ]
UCs <- UCs[UCs$categoria != "Reserva Particular do Patrimônio Natural",]
UCs <- UCs[is.na(UCs$marinho) | UCs$marinho == "", ]
BR <- terra::vect("Data/BR_UF_2024/BR_UF_2024.shp")

# Organize data

sample_sizes <- c(30, 50, 100, 500, 1000, 2500, 5000)

site_totals <- rowSums(mat_uni2_clean)
area_names <- rownames(mat_uni2_clean)  

# Results table
results_table <- data.frame(
  UC = area_names,
  Pontos_Totais = rowSums(mat_uni2_clean),
  Riqueza_Observada = rowSums(mat_uni2_clean > 0),
  stringsAsFactors = FALSE
)

# Add columns for each sample size
for(size in sample_sizes) {
  results_table[[paste0("Riqueza_", size, "_pontos")]] <- NA
}

# Process each UC 
for(i in 1:nrow(mat_uni2_clean))  {
  cat("Processing:", area_names[i], "\n")
  
  single_area_data <- mat_uni2_clean[i, , drop = FALSE]
  total_points <- sum(single_area_data)
  
  if(total_points < 10 || sum(single_area_data > 0) < 2) {
    next
  }

    # Get valid sample sizes
    valid_sizes <- sample_sizes
    
    if(length(valid_sizes) > 0) {
      
      #  iNEXT with all valid sizes at once 
      out_single <- iNEXT(t(single_area_data), q = 0, 
                          datatype = "abundance",
                          size = valid_sizes,
                          endpoint = total_points)
      
    
      estimates <- out_single$iNextEst$size_based
      
  }
}


# Completeness calculation
calc_completude <- function(observada, estimada) {
  if(is.na(estimada) || estimada == "N/A") {
    return(NA)
  }
  estimada_num <- as.numeric(estimada)
  return(round((observada / estimada_num) * 100, 1))
}

results_table$Completude_30 <- mapply(calc_completude, 
                                      results_table$Riqueza_Observada,
                                      results_table$Riqueza_30_pontos)

results_table$Completude_100 <- mapply(calc_completude, 
                                       results_table$Riqueza_Observada,
                                       results_table$Riqueza_100_pontos)

results_table$Completude_500 <- mapply(calc_completude, 
                                       results_table$Riqueza_Observada,
                                       results_table$Riqueza_500_pontos)



results_table <- results_table[order(-results_table$Pontos_Totais), ]

write.csv(results_table, "Results/riqueza_UCs.csv", row.names = FALSE)


# Species richness data
r <- ggplot(results_table, aes(x = Riqueza_Observada)) +
  geom_histogram(bins = 5, fill = "lightblue", color = "black", alpha = 0.7) +
  geom_vline(aes(xintercept = mean(Riqueza_Observada)), 
             color = "red", linetype = "dashed", linewidth = 1.2) +
  geom_vline(aes(xintercept = max(Riqueza_Observada)), 
             color = "darkgreen", linetype = "solid", linewidth = 1.2) +
  geom_vline(aes(xintercept = min(Riqueza_Observada)), 
             color = "blue", linetype = "solid", linewidth = 1.2) +
  labs(title = "Species Richness Distribution",
       x = "Species Richness", 
       y = "Frequency") +
  theme_minimal() 
 r



# PLOT PAIRED DATA


# For a 500 points

paired_data500_complete <- results_table %>%
  dplyr::select(UC, Riqueza_Observada, Riqueza_500_pontos, Pontos_Totais) %>%
  mutate(Riqueza_500_final = case_when(
      is.na(Riqueza_500_pontos) | 
        Riqueza_500_pontos %in% c("NA", "N/A", "") |
        is.na(suppressWarnings(as.numeric(Riqueza_500_pontos))) ~ Riqueza_Observada,
      TRUE ~ as.numeric(Riqueza_500_pontos)
    )
  ) %>%
  dplyr::select(UC, Riqueza_Observada, Riqueza_500_final) %>%
  pivot_longer(cols = c(Riqueza_Observada, Riqueza_500_final),
               names_to = "Type", 
               values_to = "Richness") %>%
  mutate(Type = case_when(
    Type == "Riqueza_Observada" ~ "Observed",
    Type == "Riqueza_500_final" ~ "Estimated"
  ))


gp500_complete <- ggpaired(paired_data500_complete,
                           x = "Type",
                           y = "Richness", 
                           id = "UC",
                           color = "Type",
                           palette = c("#00AFBB", "#E7B800"),
                           line.color = "gray60",
                           line.size = 0.3,
                           point.size = 1.5,
                           title = "Observed vs Estimated Richness (500 points)",
                           xlab = "Richness",
                           ylab = "Species",
                           ) +
  scale_y_log10()+
  theme_bw() +
  labs(color = NULL)

gp500_complete

ggsave("Results/gp500.png", gp500_complete)



# For a 1000 points

paired_data1000_complete <-  results_table %>%
  dplyr::select(UC, Riqueza_Observada, Riqueza_1000_pontos, Pontos_Totais) %>%
  mutate(Riqueza_1000_final = case_when(
    is.na(Riqueza_1000_pontos) | 
      Riqueza_1000_pontos %in% c("NA", "N/A", "") |
      is.na(suppressWarnings(as.numeric(Riqueza_1000_pontos))) ~ Riqueza_Observada,
    TRUE ~ as.numeric(Riqueza_1000_pontos)
  )
  ) %>%
  dplyr::select(UC, Riqueza_Observada, Riqueza_1000_final) %>%
  pivot_longer(cols = c(Riqueza_Observada, Riqueza_1000_final),
               names_to = "Type", 
               values_to = "Richness") %>%
  mutate(Type = case_when(
    Type == "Riqueza_Observada" ~ "Observed",
    Type == "Riqueza_1000_final" ~ "Estimated"
  ))

gp1000_complete <- ggpaired(paired_data1000_complete,
                           x = "Type",
                           y = "Richness", 
                           id = "UC",
                           color = "Type",
                           palette = c("#00AFBB", "#E7B800"),
                           line.color = "gray60",
                           line.size = 0.3,
                           point.size = 1.5,
                           title = "Observed vs Estimated Richness (1000 points)",
                           xlab = "Richness",
                           ylab = "Species") +
  scale_y_log10()+
  theme_bw() +
  labs(color = NULL)

gp1000_complete
ggsave("Results/gp1000.png", gp1000_complete)

# For a 100 points 

paired_data100_complete <- results_table %>%
  dplyr::select(UC, Riqueza_Observada, Riqueza_100_pontos, Pontos_Totais) %>%
  mutate(Riqueza_100_final = case_when(
    is.na(Riqueza_100_pontos) | 
      Riqueza_100_pontos %in% c("NA", "N/A", "") |
      is.na(suppressWarnings(as.numeric(Riqueza_100_pontos))) ~ Riqueza_Observada,
    TRUE ~ as.numeric(Riqueza_100_pontos)
  )
  ) %>%
  dplyr::select(UC, Riqueza_Observada, Riqueza_100_final) %>%
  pivot_longer(cols = c(Riqueza_Observada, Riqueza_100_final),
               names_to = "Type", 
               values_to = "Richness") %>%
  mutate(Type = case_when(
    Type == "Riqueza_Observada" ~ "Observed",
    Type == "Riqueza_100_final" ~ "Estimated"
  ))


gp100_complete <- ggpaired(paired_data100_complete,
                           x = "Type",
                           y = "Richness", 
                           id = "UC",
                           color = "Type",
                           palette = c("#00AFBB", "#E7B800"),
                           line.color = "gray60",
                           line.size = 0.3,
                           point.size = 1.5,
                           title = "Observed vs Estimated Richness (100 points)",
                           xlab = "Richness",
                           ylab = "Species") +
  scale_y_log10()+
  theme_bw() +
  labs(color = NULL)

gp100_complete
ggsave("Results/gp100.png", gp100_complete)


# For 5000 points

paired_data5000_complete <- results_table %>%
  dplyr::select(UC, Riqueza_Observada, Riqueza_5000_pontos, Pontos_Totais) %>%
  mutate(Riqueza_5000_final = case_when(
    is.na(Riqueza_5000_pontos) | 
      Riqueza_5000_pontos %in% c("NA", "N/A", "") |
      is.na(suppressWarnings(as.numeric(Riqueza_5000_pontos))) ~ Riqueza_Observada,
    TRUE ~ as.numeric(Riqueza_5000_pontos)
  )
  ) %>%
  dplyr::select(UC, Riqueza_Observada, Riqueza_5000_final) %>%
  pivot_longer(cols = c(Riqueza_Observada, Riqueza_5000_final),
               names_to = "Type", 
               values_to = "Richness") %>%
  mutate(Type = case_when(
    Type == "Riqueza_Observada" ~ "Observed",
    Type == "Riqueza_5000_final" ~ "Estimated"
  ))


gp5000_complete <- ggpaired(paired_data5000_complete,
                           x = "Type",
                           y = "Richness", 
                           id = "UC",
                           color = "Type",
                           palette = c("#00AFBB", "#E7B800"),
                           line.color = "gray60",
                           line.size = 0.3,
                           point.size = 1.5,
                           title = "Observed vs Estimated Richness (5000 points)",
                           xlab = "Richness",
                           ylab = "Species",
) +
  scale_y_log10()+
  theme_bw() +
  labs(color = NULL)

gp5000_complete

ggsave("Results/gp5000.png", gp5000_complete)

# MAPPING RICHNESS

terra::rast(system.file("ex/elev.tif", package="terra"))

UCs_sf <- st_as_sf(UCs)
UCs_data <- cbind(riqueza_UCs, UCs_sf)
BR_sf <- st_as_sf(BR)
UCs_data <- st_as_sf(UCs_data)

UCs_caat <- filter(UCs_data, caatinga != "NA")
UCs_ama <- filter(UCs_data, amazonia != "NA")
UCs_ma <- filter(UCs_data, matlantica != "NA")
UCs_cer <- filter(UCs_data, cerrado != "NA")
UCs_pampa <- filter(UCs_data, pampa != "NA")
UCs_pant <- filter(UCs_data, pantanal != "NA")


ama <- UCs_ama[,c("UC", "Riqueza_Observada")]  
caat <- UCs_caat[,c("UC", "Riqueza_Observada")]  
cer <- UCs_cer[,c("UC", "Riqueza_Observada")]  
ma <- UCs_ma[,c("UC", "Riqueza_Observada")]  
pampa <- UCs_pampa[,c("UC", "Riqueza_Observada")]  
pant <- UCs_pant[,c("UC", "Riqueza_Observada")]  



combined_biomes <- bind_rows(
  ama %>% mutate(Biome = "Amazonia"),
  caat %>% mutate(Biome = "Caatinga"),
  cer %>% mutate(Biome = "Cerrado"),
  ma %>% mutate(Biome = "Mata Atlantica"),
  pampa %>% mutate(Biome = "Pampa"),
  pant %>% mutate(Biome = "Pantanal")
) %>%
  rename(Richness = Riqueza_Observada)


save(combined_biomes, file = "Data/combined_biomes.RData")

biome_totals <- combined_biomes %>%
  group_by(Biome) %>%
  summarise(Richness = sum(Richness, na.rm = TRUE),
            .groups = 'drop')

# Plot


p <- ggplot(biome_totals, aes(x = Biome,
                              y = Richness,
                              fill = Biome)) +
  geom_col(alpha = 0.8, color = "black") +
  scale_fill_viridis_d(name = "Biome",
                       option = "plasma") +
  labs(title = "Total Species Richness in UCs by Biome",
       x = "Biome", 
       y = "Species Richness") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

p




# Plot
p2 <- ggplot(combined_biomes, aes(x = `Richness`, 
                               y = `Biome`, 
                               fill = Richness)) +
  geom_density_ridges_gradient(aes(fill = after_stat(x)),
                               rel_min_height = 0.01) +
  scale_fill_viridis( 
                     option = "plasma") +
  labs(
       x = "Species Richness", 
       y = "Biome") +
  theme_ridges() +
  theme(
    legend.position="right",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )
p2

gr3 <- ggplot() +
  geom_sf(data = BR_sf, 
          color = "black", 
          fill = "white") +
  geom_sf(data = UCs_data, 
          color = 0, 
          fill = "lightblue", 
          size = 0.1) +
  geom_sf(data = UCs_data, 
          aes(fill = Riqueza_Observada), 
          color = 0,
          size = 0.1) +
  scale_fill_viridis_c(name = "Species\nRichness", 
                       option = "plasma",
                       trans = "log10",
                       na.value = "blue") +
  theme_void() 

gr3

ggsave("Results/Invasive_Species_Richness.tiff",gr3)



p2 <- p2 + 
  theme(
    plot.margin = margin(5, 5, 5, 5),        
    axis.title.x = element_text(size = 8, hjust = 0.5),    
    axis.title.y = element_text(size = 8, hjust = 0.2, vjust = 0),
    axis.text = element_text(size = 7),       
    legend.title = element_text(size = 8),    
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.4, "cm")         
  )

combined_plot <- ggdraw() +
  draw_plot(gr3) +
  draw_plot(p2, 
            x = 0.60, y = 0,              
            width = 0.4, height = 0.35) +
  draw_plot_label(c("(a)", "(b)"), 
                  c(0.16, 0.60), c(0.95, 0.35),  
                  size = 12) +
  draw_label("Invasive Species Richness in Brazil by Conservations Units and Biome",
             x = 0.5, y = 0.99,           
             hjust = 0.5, vjust = 1,
             size = 12, fontface = "bold")
combined_plot
ggsave("Figures/combined_plot.pdf", combined_plot, 
       width = 14, height = 6, device = "pdf")

