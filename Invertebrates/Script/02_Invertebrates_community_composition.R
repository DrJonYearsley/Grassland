####################################################################
## Plotting Invertebrates data from the Grassland Resilience      ##
## Experiment in 2017                                             ##
## Data provided by Matthew Magilton                              ##
##                                                                ##
## Figures produced:                                              ##
## Community composition per site (each figure depicts one day)   ##
## Community composition per day (each figure depicts one site)   ##
## Relative abundance of each order per site and day              ##
##                                                                ##
## Author of the script:                                          ##
## Maja Ilic (QUB)                                                ##
##                                                                ##
## First modified: 11 Feb 2020                                    ##
## Last modified: 11 Feb 2020                                     ##
####################################################################

#===========================================
#### Clear objects from the workspace ####

rm(list = ls())

#===========================================
#### Set working directory ####

# Maja's desktop PC

setwd("C:/Users/3054311/Documents/My Documents/Grassland project/07_Invertebrates")

data.dir <- paste0(getwd(),"/Data/")
figures.dir <- paste0(getwd(),"/Figures/")
script.dir <- paste0(getwd(),"/Script/")

#======================================
#### Packages ####

library(ggplot2)      ## for plots
library(dplyr)        ## for easier data management
library(tidyr)        ## same as above
library(RColorBrewer) ## for managing colors
library(cowplot)  
library(gridExtra)
library(grid)

#======================================
# Import data ####

df <- read.csv(paste0(data.dir,"Grassland Resilience 2017 - Invertebrates raw data.csv"), sep = ",", header = T)  ## read in the dataset

#======================================
# Extract only the relevant columns and put them in a logical order ####

df.new <- df[,c(1:6,13:39)]

head(df.new)
str(df.new)

#======================================
# Change from wide to long format ####

df.long <- gather(df.new, key = "Order", value = "Abundance",
                  -Region, -Farm, -site_ID, -treatment, -Plot, -day)

#======================================
# Calculate relative abundance of each order in the community (per each plot!) ####

df.long <- df.long %>% 
  group_by(Region,Farm,site_ID,treatment,Plot,day) %>% 
  mutate(Relative.Abundance = Abundance/sum(Abundance, na.rm = T)*100)

#======================================
# Add another column "habitat" with the levels "soil" and "other" (talk to Matt about this) ####

soil.orders <- names(df.new)[7:17]
other.orders <- names(df.new)[18:33]

df.long[which(df.long$Order %in% soil.orders), "habitat"] <- "soil"
df.long[which(df.long$Order %in% other.orders), "habitat"] <- "other"

#======================================
# Plot the community composition per site and plot ####
# First, subset the dataset per day

day0 <- df.long %>% filter(day == 0)
day32 <- df.long %>% filter(day == 32)
day64 <- df.long %>% filter(day == 64)

# Get colors

cols <- colorRampPalette(brewer.pal(length(unique(df.long$Order)), "Paired"))
mypalette <- cols(length(unique(df.long$Order)))

# Day 0

g0 <- ggplot(day0, aes(x = Plot, y = Relative.Abundance)) +
  geom_bar(stat = "identity", width = 0.8, aes(color = Order, fill = Order)) +
  facet_grid(Region ~ Farm) +
  geom_vline(xintercept = 4.5, color = "gray50") +
  theme_minimal() +
  theme(plot.caption = element_text(face = "italic"),
        panel.spacing = unit(0.5, "line"),
        strip.background = element_rect()) +
  scale_y_continuous(limits = c(-0.5,100.5), breaks = c(0,25,50,75,100), exp = c(0,0)) +
  scale_color_manual(values = mypalette) +
  scale_fill_manual(values = mypalette) +
  labs(x = "Plot",
       y = "Relative Abundance (%)",
       title = "Grassland Resilience - Invertebrates",
       subtitle = "Day 0",
       caption = "Maja Ilic (QUB), 2020")

ggsave(paste0(figures.dir,"Grassland Resilience Invertebrates Day 0.png"), g0, width = 12, height = 9)

# Day 32

g32 <- ggplot(day32, aes(x = Plot, y = Relative.Abundance)) +
  geom_bar(stat = "identity", width = 0.8, aes(color = Order, fill = Order)) +
  facet_grid(Region ~ Farm) +
  geom_vline(xintercept = 4.5, color = "gray50") +
  theme_minimal() +
  theme(plot.caption = element_text(face = "italic"),
        panel.spacing = unit(0.5, "line"),
        strip.background = element_rect()) +
  scale_y_continuous(limits = c(-0.5,100.5), breaks = c(0,25,50,75,100), exp = c(0,0)) +
  scale_color_manual(values = mypalette) +
  scale_fill_manual(values = mypalette) +
  labs(x = "Plot",
       y = "Relative Abundance (%)",
       title = "Grassland Resilience - Invertebrates",
       subtitle = "Day 32",
       caption = "Maja Ilic (QUB), 2020")

ggsave(paste0(figures.dir,"Grassland Resilience Invertebrates Day 32.png"), g32, width = 12, height = 9)

# Day 64

g64 <- ggplot(day64, aes(x = Plot, y = Relative.Abundance)) +
  geom_bar(stat = "identity", width = 0.8, aes(color = Order, fill = Order)) +
  facet_grid(Region ~ Farm) +
  geom_vline(xintercept = 4.5, color = "gray50") +
  theme_minimal() +
  theme(plot.caption = element_text(face = "italic"),
        panel.spacing = unit(0.5, "line"),
        strip.background = element_rect()) +
  scale_y_continuous(limits = c(-0.5,100.5), breaks = c(0,25,50,75,100), exp = c(0,0)) +
  scale_color_manual(values = mypalette) +
  scale_fill_manual(values = mypalette) +
  labs(x = "Plot",
       y = "Relative Abundance (%)",
       title = "Grassland Resilience - Invertebrates",
       subtitle = "Day 64",
       caption = "Maja Ilic (QUB), 2020")

ggsave(paste0(figures.dir,"Grassland Resilience Invertebrates Day 64.png"), g64, width = 12, height = 9)

#======================================
#### Plot average community composition per treatment and site over time ####
# Calculate the mean (relative) abundance of each order per treatment, site and day

mean.comm <- df.long %>% 
  group_by(Region,Farm,site_ID,treatment,day,Order,habitat) %>% 
  summarize(Mean.Abu = mean(Abundance, na.rm = T),
            SD.Abu = sd(Abundance, na.rm = T),
            Mean.Rel.Abu = mean(Relative.Abundance, na.rm = T),
            SD.Rel.Abu = sd(Relative.Abundance, na.rm = T))

# Relative Abundance

regions <- unique(mean.comm$Region)

for(k in 1:length(regions)) {
  
  g.avg <- ggplot(mean.comm[which(mean.comm$Region == regions[k]),], 
                  aes(x = as.factor(day), y = Mean.Rel.Abu)) +
    geom_bar(stat = "identity", width = 0.8, aes(color = Order, fill = Order)) +
    facet_grid(treatment ~ Farm) +
    geom_vline(xintercept = 4.5, color = "gray50") +
    theme_bw() +
    theme(plot.caption = element_text(face = "italic"),
          panel.spacing = unit(0.5, "line")) +
    scale_y_continuous(limits = c(-0.5,100.5), breaks = c(0,25,50,75,100), exp = c(0,0)) +
    scale_color_manual(values = mypalette) +
    scale_fill_manual(values = mypalette) +
    labs(x = "Days after the end of the drought",
         y = "Relative Abundance (%)",
         title = paste("Invertebrates -", regions[k]),
         subtitle = "Average relative abundance of each order over time",
         caption = "Maja Ilic (QUB), 2020")
  
  ggsave(paste0(figures.dir,"Grassland Resilience Invertebrates in ", regions[k], " Rel Abu.png"), 
         g.avg, width = 10, height = 6)

}

#======================================
#### Plot average abundance of each order per treatment and site over time 

orders <- unique(mean.comm$Order)

for(m in 1:length(orders)) {
  
  g.order <- ggplot(mean.comm[which(mean.comm$Order == orders[m]),], 
                    aes(x = day, y = Mean.Abu)) +
    geom_point(size = 3, shape = 21, alpha = 0.5,
               aes(color = Order, fill = Order)) +
    geom_smooth(method = "lm", se = F, color = mypalette[m], size = 1) +
    facet_grid(treatment ~ Region) +
    theme_minimal() +
    theme(plot.caption = element_text(face = "italic"),
          panel.spacing = unit(0.5, "line"),
          strip.background = element_rect()) +
    scale_color_manual(values = mypalette[m]) +
    scale_fill_manual(values = mypalette[m]) +
    scale_x_continuous(limits = c(-1,65), breaks = c(0,32,64)) +
    labs(x = "Days after the end of the drought",
         y = "Abundance",
         title = paste("Invertebrates -", orders[m]),
         subtitle = "Average abundance of each order within one treatment plotted over time\nEach data point represents one site within the given region",
         caption = "Maja Ilic (QUB), 2020")
  
  ggsave(paste0(figures.dir,"Average abundance of ", orders[m],".png"), 
         g.order, width = 10, height = 6.25)
  
}

#======================================
#### Soil invertebrates plotted per Region ####

mean.soil <- mean.comm %>% filter(habitat == "soil")

mycols.soil <- mypalette[which(orders %in% soil.orders)]

col.C <- rgb(248, 118, 109, max = 255)
col.D <- rgb(0, 191, 196, max = 255)

for(k in 1:length(regions)) {
  
  mean.region <- mean.soil %>% filter(Region == regions[k])
  mean.region.C <- mean.region %>%  filter(treatment == "Control")
  mean.region.D <- mean.region %>% filter(treatment == "Drought")
  
  # Control
  
  g.C <- ggplot(mean.region.C, 
                aes(x = day, y = Mean.Abu)) +
    geom_point(size = 3, shape = 21, alpha = 0.5,
               aes(color = Order, fill = Order)) +
    geom_smooth(method = "lm", se = F, aes(color = Order), size = 1) +
    facet_wrap(~ Order, scales = "free_y", nrow = 1) +
    theme_bw() +
    theme(plot.caption = element_text(face = "italic"),
          panel.spacing = unit(0.5, "line"),
          legend.position = "none",
          axis.title = element_blank(),
          strip.background = element_rect(fill = col.C)) +
    scale_x_continuous(limits = c(-1,65), breaks = c(0,32,64)) +
    scale_color_manual(values = mycols.soil) +
    scale_fill_manual(values = mycols.soil)
  
  # Drought
  
  g.D <- ggplot(mean.region.D, 
                aes(x = day, y = Mean.Abu)) +
    geom_point(size = 3, shape = 21, alpha = 0.5,
               aes(color = Order, fill = Order)) +
    geom_smooth(method = "lm", se = F, aes(color = Order), size = 1) +
    facet_wrap(~ Order, scales = "free_y", nrow = 1) +
    theme_bw() +
    theme(plot.caption = element_text(face = "italic"),
          panel.spacing = unit(0.5, "line"),
          legend.position = "none",
          axis.title = element_blank(),
          strip.background = element_rect(fill = col.D)) +
    scale_x_continuous(limits = c(-1,65), breaks = c(0,32,64)) +
    scale_color_manual(values = mycols.soil) +
    scale_fill_manual(values = mycols.soil)
  
  # Combine plots 
  
  plot.region <- plot_grid(g.C,g.D, 
                           nrow = 2, align = "v")
  
  y.grob <- textGrob("Abundance", 
                     gp = gpar(fontsize = 12), rot = 90)
  
  x.grob <- textGrob("Days since the end of the drought", 
                     gp = gpar(fontsize = 12))
  
  title <- paste("Soil Invertebrates -", regions[k])
  sub.title.1 <- "Average abundance of each order within one treatment plotted over time"
  sub.title.2 <- "Each data point represents one site within the given region"
  sub.title.3 <- c("Control "," vs "," Drought")
  subtitle.colors <- c(col.C,"black",col.D)
  
  plot.final <- grid.arrange(arrangeGrob(plot.region, 
                                         left = y.grob, 
                                         bottom = x.grob))
  
  plot.final <- grid.arrange(plot.final,
                             top = tableGrob(t(sub.title.3),
                                            theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                   base_colour = subtitle.colors,
                                                                   base_size = 12)))
  
  plot.final <- grid.arrange(plot.final,
                             top = tableGrob(t(sub.title.2),
                                             theme = ttheme_minimal(padding = unit(c(0,4),'mm'),
                                                                    base_colour = "black",
                                                                    base_size = 12)))
  
  plot.final <- grid.arrange(plot.final,
                             top = tableGrob(t(sub.title.1),
                                             theme = ttheme_minimal(padding = unit(c(0,4),'mm'),
                                                                    base_colour = "black",
                                                                    base_size = 12)))
  
  plot.final <- grid.arrange(plot.final,
                             top = tableGrob(t(title),
                                             theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                    base_colour = "black",
                                                                    base_size = 16)))
  
  ggsave(paste0(figures.dir,"Average abundance of Soil Invertebrates in ", regions[k],".png"), 
         plot.final, width = 18, height = 6)
  
}

#======================================
#### Other invertebrates plotted per Region ####

mean.other <- mean.comm %>% filter(habitat == "other")

mycols.other <- mypalette[which(orders %in% other.orders)]

for(k in 1:length(regions)) {
  
  mean.region <- mean.other %>% filter(Region == regions[k])
  mean.region.C <- mean.region %>%  filter(treatment == "Control")
  mean.region.D <- mean.region %>% filter(treatment == "Drought")
  
  # Control
  
  g.C <- ggplot(mean.region.C, 
                aes(x = day, y = Mean.Abu)) +
    geom_point(size = 3, shape = 21, alpha = 0.5,
               aes(color = Order, fill = Order)) +
    geom_smooth(method = "lm", se = F, aes(color = Order), size = 1) +
    facet_wrap(~ Order, scales = "free_y", nrow = 1) +
    theme_bw() +
    theme(plot.caption = element_text(face = "italic"),
          panel.spacing = unit(0.5, "line"),
          legend.position = "none",
          axis.title = element_blank(),
          strip.background = element_rect(fill = col.C)) +
    scale_x_continuous(limits = c(-1,65), breaks = c(0,32,64)) +
    scale_color_manual(values = mycols.other) +
    scale_fill_manual(values = mycols.other)
  
  # Drought
  
  g.D <- ggplot(mean.region.D, 
                aes(x = day, y = Mean.Abu)) +
    geom_point(size = 3, shape = 21, alpha = 0.5,
               aes(color = Order, fill = Order)) +
    geom_smooth(method = "lm", se = F, aes(color = Order), size = 1) +
    facet_wrap(~ Order, scales = "free_y", nrow = 1) +
    theme_bw() +
    theme(plot.caption = element_text(face = "italic"),
          panel.spacing = unit(0.5, "line"),
          legend.position = "none",
          axis.title = element_blank(),
          strip.background = element_rect(fill = col.D)) +
    scale_x_continuous(limits = c(-1,65), breaks = c(0,32,64)) +
    scale_color_manual(values = mycols.other) +
    scale_fill_manual(values = mycols.other)
  
  # Combine plots 
  
  plot.region <- plot_grid(g.C,g.D, 
                           nrow = 2, align = "v")
  
  y.grob <- textGrob("Abundance", 
                     gp = gpar(fontsize = 12), rot = 90)
  
  x.grob <- textGrob("Days since the end of the drought", 
                     gp = gpar(fontsize = 12))
  
  title <- paste("Other Invertebrates -", regions[k])
  sub.title.1 <- "Average abundance of each order within one treatment plotted over time"
  sub.title.2 <- "Each data point represents one site within the given region"
  sub.title.3 <- c("Control "," vs "," Drought")
  subtitle.colors <- c(col.C,"black",col.D)
  
  plot.final <- grid.arrange(arrangeGrob(plot.region, 
                                         left = y.grob, 
                                         bottom = x.grob))
  
  plot.final <- grid.arrange(plot.final,
                             top = tableGrob(t(sub.title.3),
                                             theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                    base_colour = subtitle.colors,
                                                                    base_size = 12)))
  
  plot.final <- grid.arrange(plot.final,
                             top = tableGrob(t(sub.title.2),
                                             theme = ttheme_minimal(padding = unit(c(0,4),'mm'),
                                                                    base_colour = "black",
                                                                    base_size = 12)))
  
  plot.final <- grid.arrange(plot.final,
                             top = tableGrob(t(sub.title.1),
                                             theme = ttheme_minimal(padding = unit(c(0,4),'mm'),
                                                                    base_colour = "black",
                                                                    base_size = 12)))
  
  plot.final <- grid.arrange(plot.final,
                             top = tableGrob(t(title),
                                             theme = ttheme_minimal(padding = unit(c(0,8),'mm'),
                                                                    base_colour = "black",
                                                                    base_size = 16)))
  
  ggsave(paste0(figures.dir,"Average abundance of Other Invertebrates in ", regions[k],".png"), 
         plot.final, width = 22, height = 6)
  
}

