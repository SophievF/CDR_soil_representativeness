### PCA analysis: global representative ###
### Sophie von Fromm ###
### 2023-10-30 ###

library(tidyverse)
library(FactoMineR)
library(factoextra)
library(ggpubr)
library(missMDA)

## Load data
data_global <- read_csv("./data/global_climate_data_2023-11-09.csv") %>% 
  rename(long_dec_deg = x, 
         lat_dec_deg = y) %>% 
  mutate(sampled = "no")

head(data_global)
skimr::skim_without_charts(data_global)

data_sampled <- read_csv("./data/all_data_combined_2024-09-29.csv") %>% 
  mutate(sampled = "yes")

head(data_sampled)
skimr::skim_without_charts(data_sampled)

## Combine both datasets and add labels for colors and grouping
data_all <- data_global %>% 
  full_join(data_sampled) 

skimr::skim_without_charts(data_all)

## Plot data distribution (w/o climate data)
world <- map_data("world") %>% 
  filter(region != "Antarctica")

# Create function for plotting
plot_map_fun <- function(variable, color_points){
  data_all %>% 
    filter(sampled == "yes") %>% 
    filter(label == variable) %>% 
    ggplot() +
    geom_map(
      data = world, map = world,
      aes(long, lat, map_id = region),
      color = "white", fill = "lightgrey", linewidth = 0.01)  +
    geom_point(aes(x = long_dec_deg, y = lat_dec_deg), color = color_points, 
               shape = 21, size = 0.5) +
    theme_bw(base_size = 12) +
    theme(axis.text = element_text(color = "black"),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_blank()) +
    coord_sf() +
    scale_color_manual(values = CDR_color) +
    scale_y_continuous(expand = c(0,0), breaks = seq(-40,80,40)) +
    scale_x_continuous(expand = c(0,0), breaks = seq(-150,150,75))
}

## Plot data by groups and arrange figures
# Good representation
c_stocks <- plot_map_fun(variable = "C stock", color_points = "#1E5664") +
  ggtitle("a) C stocks")

soil_properties <- plot_map_fun(variable = "soil properties", color_points = "#6A269C") +
  ggtitle("b) Soil properties")

soil_spectral <- plot_map_fun(variable = "spectral", color_points = "#6A269C") +
  ggtitle("c) Spectral")

plot_1 <- annotate_figure(
  ggarrange(c_stocks, soil_properties, soil_spectral, nrow = 1),
  left = text_grob("Good representation", face = "bold",
                   rot = 90, size = 14)
)

ggsave(file = paste0("./figures/map_group_1_", Sys.Date(), ".jpeg"),
       width = 11, height = 6)

# Limited representation
respiration <- plot_map_fun(variable = "respiration", color_points = "#1E5664") +
  ggtitle("d) Soil respiration")

microbial <- plot_map_fun(variable = "microbial biomass", color_points = "#6A269C") +
  ggtitle("e) Microbial biomass")

radiocarbon <- plot_map_fun(variable = "radiocarbon", color_points = "#CC1599") +
  ggtitle("f) Radiocarbon")

plot_2 <- annotate_figure(
  ggarrange(respiration, microbial, radiocarbon,  nrow = 1),
  left = text_grob("Limited representation", face = "bold",
                   rot = 90, size = 14)
)

ggsave(file = paste0("./figures/map_group_2_", Sys.Date(), ".svg"),
       width = 11, height = 6)

# Poor representation
maom <- plot_map_fun(variable = "MAOM", color_points = "#1E5664") +
  ggtitle("g) MAOM")

necro <- plot_map_fun(variable = "microbial necromass", color_points = "#6A269C") +
  ggtitle("h) Microbial necromass")

time <- plot_map_fun(variable = "time series", color_points = "#CC1599") +
  ggtitle("i) Timeseries")

plot_3 <- annotate_figure(
  ggarrange(maom, necro, time, nrow = 1),
  left = text_grob("Poor representation", face = "bold",
                   rot = 90, size = 14)
)

ggsave(file = paste0("./figures/map_group_3_", Sys.Date(), ".svg"),
       width = 11, height = 6)

## Prepare and perform PCA analysis
data_all_pca <- data_all %>% 
  dplyr::select(MAT:GPP)

skimr::skim_without_charts(data_all_pca)

# use imputePCA to gap-fill missing data for MAT, PET, GPP, MAP
# takes quite some time (~1hr on my windows machine)
nb <- estim_ncpPCA(data_all_pca, ncp.min = 0, ncp.max = 5, method.cv = "Kfold")

data_pca_impute <- imputePCA(data_all_pca, ncp = 3)

pca_all <- data_pca_impute$completeObs %>% 
  PCA(graph = FALSE, ind.sup = c((nrow(data_global) + 1):nrow(data_all)))

factoextra::fviz_pca_biplot(pca_all, col.var = "black", col.ind  = "grey", 
                            col.ind.sup = "red", geom.ind = "point") +
  theme_bw(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "none")

## Extract PCA data for easier plotting
ind <- data.frame(pca_all$ind$coord)
ind.sup <- data.frame(pca_all$ind.sup$coord)
pca_all_df <- rbind(ind, ind.sup)

var <- facto_summarize(pca_all, element = "var", 
                       result = c("coord", "contrib", "cos2"))

var %>% 
  dplyr::select(name, Dim.1, Dim.2)

#factor for scaling variables to space of individuals (for Dim.1 and Dim.2)
r <- min((max(ind[, "Dim.1"]) - min(ind[, "Dim.1"])/(max(var[, "Dim.1"]) - 
                                                       min(var[, "Dim.1"]))), 
         (max(ind[, "Dim.2"]) - min(ind[, "Dim.2"])/(max(var[, "Dim.2"]) - 
                                                       min(var[, "Dim.2"]))))
# r <- 9.015744

pca_data <- cbind(data_all, pca_all_df) %>% 
  tibble()

# Save data for later 
write_csv(pca_data, file = paste0("./data/PCA_output_all_data_", 
                                  Sys.Date(), ".csv"))

write_csv(var, file = paste0("./data/var_output_all_data_", 
                             Sys.Date(), ".csv"))

# pca_data <- read_csv("./data/PCA_output_all_data_2024-09-29.csv")
# 
# var <- read_csv("./data/var_output_all_data_2024-09-29.csv")

## PCA plotting
# Avoid overlapping of labels
var$jit <- with(var, ifelse(name == "GPP" | name == "MAP", 1, 2))

# function to plot data
plot_pca_fun <- function(variable, color_points){
  ggplot() +
    # PCA coordinate lines
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    # global PCA space
    geom_point(data = pca_data %>% 
                 filter(sampled == "no"),
               aes(x = Dim.1, y = Dim.2), color = "#bdbdbd", shape = 21) +
    # variable specific PCA space
    geom_point(data = pca_data %>% 
                 filter(sampled == "yes") %>% 
                 filter(label == variable),
               aes(x = Dim.1, y = Dim.2), color = color_points, shape = 21) +
    # PCA arrows (climate data)
    geom_segment(data = var, aes(x = 0, xend = Dim.1*r*0.7, y = 0, yend = Dim.2*r*0.7),
                 linewidth = 0.75, arrow = arrow(length = unit(0.02, "npc"))) +
    # PCA arrow labels (climate data)
    ggrepel::geom_text_repel(data = var, aes(x = Dim.1*r*0.75, y = Dim.2*r*0.75, label = name),
                             size = 4, color = "#B1615F", fontface = "bold", nudge_x = 0.4) +
    theme_bw(base_size = 12) +
    theme(axis.text = element_text(color = "black"),
          panel.grid = element_blank()) +
    scale_x_continuous(limits = c(-4,8), expand = c(0,0)) +
    scale_y_continuous(limits = c(-7,6), expand = c(0,0))
}

## Plot data by groups and arrange figures
# Good representation
c_stocks <- plot_pca_fun(variable = "C stock", color_points = "#1E5664") +
  ggtitle("a) C stocks")

soil_properties <- plot_pca_fun(variable = "soil properties", color_points = "#6A269C") +
  ggtitle("b) Soil properties")

soil_spectral <- plot_pca_fun(variable = "spectral", color_points = "#6A269C") +
  ggtitle("c) Spectral")

plot_1 <- annotate_figure(
  ggarrange(c_stocks, soil_properties, soil_spectral, nrow = 1),
  left = text_grob("Good representation", face = "bold",
                   rot = 90, size = 16)
)

# Limited representation
respiration <- plot_pca_fun(variable = "respiration", color_points = "#1E5664") +
  ggtitle("d) Soil respiration") 

microbial <- plot_pca_fun(variable = "microbial biomass", color_points = "#6A269C") +
  ggtitle("e) Microbial biomass")

radiocarbon <- plot_pca_fun(variable = "radiocarbon", color_points = "#CC1599") +
  ggtitle("f) Radiocarbon")

plot_2 <- annotate_figure(
  ggarrange(respiration, microbial, radiocarbon,  nrow = 1),
  left = text_grob("Limited representation", face = "bold",
                   rot = 90, size = 16)
)

# Poor representation
maom <- plot_pca_fun(variable = "MAOM", color_points = "#1E5664") +
  ggtitle("g) MAOM")

necro <- plot_pca_fun(variable = "microbial necromass", color_points = "#6A269C") +
  ggtitle("h) Microbial necromass")

time <- plot_pca_fun(variable = "time series", color_points = "#CC1599") +
  ggtitle("i) Timeseries")

plot_3 <- annotate_figure(
  ggarrange(maom, necro, time, nrow = 1),
  left = text_grob("Poor representation", face = "bold",
                   rot = 90, size = 16)
)

# Plot all groups together and save final figure
ggarrange(plot_1, plot_2, plot_3, nrow = 3)
ggsave(file = paste0("./figures/PCA_all_groups_", Sys.Date(), ".emf"),
       device = {function(filename, ...) devEMF::emf(file = filename, ...)},
       width = 9, height = 10)

## Density plots
plot_den_fun <- function(dimension, variable, color_line){
  ggplot() +
    geom_density(data = pca_data %>% 
                   filter(sampled == "no"),
                 aes(x = .data[[dimension]]), color = "#bdbdbd", linewidth = 1) +
    geom_density(data = pca_data %>% 
                   filter(sampled == "yes") %>% 
                   filter(label == variable), 
                 aes(x = .data[[dimension]]), color = color_line, linewidth = 1) +
    theme_bw(base_size = 12) +
    theme(axis.text = element_text(color = "black"),
          panel.grid = element_blank()) +
    scale_y_continuous(limits = c(0,1.2), expand = c(0,0))
}

## Plot data by groups and for dimension 1 and arrange figures
# Good representation
c_stocks <- plot_den_fun(variable = "C stock", color_line = "#1E5664",
                         dimension = "Dim.1") +
  scale_x_continuous(limits = c(-4,8), expand = c(0,0)) +
  ggtitle("a) C stocks")

soil_properties <- plot_den_fun(variable = "soil properties", color_line = "#6A269C",
                                dimension = "Dim.1") +
  scale_x_continuous(limits = c(-4,8), expand = c(0,0)) +
  ggtitle("b) Soil properties")

soil_spectral <- plot_den_fun(variable = "spectral", color_line = "#6A269C",
                              dimension = "Dim.1") +
  scale_x_continuous(limits = c(-4,8), expand = c(0,0)) +
  ggtitle("c) Spectral")

plot_1 <- annotate_figure(
  ggarrange(c_stocks, soil_properties, soil_spectral, nrow = 1),
  left = text_grob("Good representation", face = "bold",
                   rot = 90, size = 16)
)

# Limited representation
respiration <- plot_den_fun(variable = "respiration", color_line = "#1E5664",
                            dimension = "Dim.1") +
  scale_x_continuous(limits = c(-4,8), expand = c(0,0)) +
  ggtitle("d) Soil respiration")

microbial <- plot_den_fun(variable = "microbial biomass", color_line = "#6A269C",
                          dimension = "Dim.1") +
  scale_x_continuous(limits = c(-4,8), expand = c(0,0)) +
  ggtitle("e) Microbial biomass")

radiocarbon <- plot_den_fun(variable = "radiocarbon", color_line = "#CC1599",
                            dimension = "Dim.1") +
  scale_x_continuous(limits = c(-4,8), expand = c(0,0)) +
  ggtitle("f) Radiocarbon")

plot_2 <- annotate_figure(
  ggarrange(respiration, microbial, radiocarbon,  nrow = 1),
  left = text_grob("Limited representation", face = "bold",
                   rot = 90, size = 16)
)

# Poor representation
maom <- plot_den_fun(variable = "MAOM", color_line = "#1E5664",
                     dimension = "Dim.1") +
  scale_x_continuous(limits = c(-4,8), expand = c(0,0)) +
  ggtitle("g) MAOM")

necro <- plot_den_fun(variable = "microbial necromass", color_line = "#6A269C",
                      dimension = "Dim.1") +
  scale_x_continuous(limits = c(-4,8), expand = c(0,0)) +
  ggtitle("h) Microbial necromass")

time <- plot_den_fun(variable = "time series", color_line = "#CC1599",
                     dimension = "Dim.1") +
  scale_x_continuous(limits = c(-4,8), expand = c(0,0)) +
  ggtitle("i) Timeseries")

plot_3 <- annotate_figure(
  ggarrange(maom, necro, time, nrow = 1),
  left = text_grob("Poor representation", face = "bold",
                   rot = 90, size = 16)
)

# Plot all groups together and save final figure
ggarrange(plot_1, plot_2, plot_3, nrow = 3)
ggsave(file = paste0("./figures/Density_Dim.1_all_groups_", Sys.Date(), ".emf"),
       device = {function(filename, ...) devEMF::emf(file = filename, ...)},
       width = 9, height = 10)

## Plot data by groups and for dimension 2 and arrange figures
# Good representation
c_stocks <- plot_den_fun(variable = "C stock", color_line = "#1E5664",
                         dimension = "Dim.2") +
  scale_x_continuous(limits = c(-7,6), expand = c(0,0)) +
  ggtitle("a) C stocks")

soil_properties <- plot_den_fun(variable = "soil properties", color_line = "#6A269C",
                                dimension = "Dim.2") +
  scale_x_continuous(limits = c(-7,6), expand = c(0,0)) +
  ggtitle("b) Soil properties")

soil_spectral <- plot_den_fun(variable = "spectral", color_line = "#6A269C",
                              dimension = "Dim.2") +
  scale_x_continuous(limits = c(-7,6), expand = c(0,0)) +
  ggtitle("c) Spectral")

plot_1 <- annotate_figure(
  ggarrange(c_stocks, soil_properties, soil_spectral, nrow = 1),
  left = text_grob("Good representation", face = "bold",
                   rot = 90, size = 16)
)

# Limited representation
respiration <- plot_den_fun(variable = "respiration", color_line = "#1E5664",
                            dimension = "Dim.2") +
  scale_x_continuous(limits = c(-7,6), expand = c(0,0)) +
  ggtitle("d) Soil respiration")

microbial <- plot_den_fun(variable = "microbial biomass", color_line = "#6A269C",
                          dimension = "Dim.2") +
  scale_x_continuous(limits = c(-7,6), expand = c(0,0)) +
  ggtitle("e) Microbial biomass")

radiocarbon <- plot_den_fun(variable = "radiocarbon", color_line = "#CC1599",
                            dimension = "Dim.2") +
  scale_x_continuous(limits = c(-7,6), expand = c(0,0)) +
  ggtitle("f) Radiocarbon")

plot_2 <- annotate_figure(
  ggarrange(respiration, microbial, radiocarbon,  nrow = 1),
  left = text_grob("Limited representation", face = "bold",
                   rot = 90, size = 16)
)

# Poor representation
maom <- plot_den_fun(variable = "MAOM", color_line = "#1E5664",
                     dimension = "Dim.2") +
  scale_x_continuous(limits = c(-7,6), expand = c(0,0)) +
  ggtitle("g) MAOM")

necro <- plot_den_fun(variable = "microbial necromass", color_line = "#6A269C",
                      dimension = "Dim.2") +
  scale_x_continuous(limits = c(-7,6), expand = c(0,0)) +
  ggtitle("h) Microbial ecromass")

time <- plot_den_fun(variable = "time series", color_line = "#CC1599",
                     dimension = "Dim.2") +
  scale_x_continuous(limits = c(-7,6), expand = c(0,0)) +
  ggtitle("i) Timeseries")

plot_3 <- annotate_figure(
  ggarrange(maom, necro, time, nrow = 1),
  left = text_grob("Poor representation", face = "bold",
                   rot = 90, size = 16)
)

# Plot all groups together and save final figure
ggarrange(plot_1, plot_2, plot_3, nrow = 3)
ggsave(file = paste0("./figures/Density_Dim.2_all_groups_", Sys.Date(), ".tiff"),
       width = 9, height = 10)

ggsave(file = paste0("./figures/Density_Dim.2_all_groups_", Sys.Date(), ".emf"),
       device = {function(filename, ...) devEMF::emf(file = filename, ...)},
       width = 9, height = 10)

