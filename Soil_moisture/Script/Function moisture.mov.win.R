moisture.mov.win <- function(data = df, dates = events, 
                             daytime.hr = 12, duration.hr = 24, 
                             ylim = c(0,3.5)) {
  
  Sys.setenv(tz='UTC')  ## change timezone to UTC
  
  sites <- unique(dates$Site)  ## extract all sites
  
  for (i in 1:length(sites)) {
    
    mysite <- data %>% 
      filter(Site == sites[i])
    
    ## Get all dates of sampling days starting at specific day time
    
    day0 <- as.POSIXct(mdy(dates$EndDrought[dates$Site == sites[i]])) + hours(daytime.hr)
    day8 <- as.POSIXct(mdy(dates$Day8[dates$Site == sites[i]])) + hours(daytime.hr)
    day16 <- as.POSIXct(mdy(dates$Day16[dates$Site == sites[i]])) + hours(daytime.hr)
    day32 <- as.POSIXct(mdy(dates$Day32[dates$Site == sites[i]])) + hours(daytime.hr)
    day64 <- as.POSIXct(mdy(dates$Day64[dates$Site == sites[i]])) + hours(daytime.hr)
    
    ## Create intervals for each sampling day of specific duration
    
    if (duration.hr == 0) {
      
      day0.interval <- day0
      day8.interval <- day8
      day16.interval <- day16 
      day32.interval <- day32 
      day64.interval <- day64
      
    }
    
    if (duration.hr > 0) {
      
      day0.interval <- day0 + hours(0:(duration.hr-1))
      day8.interval <- day8 + hours(0:(duration.hr-1))
      day16.interval <- day16 + hours(0:(duration.hr-1))
      day32.interval <- day32 + hours(0:(duration.hr-1))
      day64.interval <- day64 + hours(0:(duration.hr-1))
      
    }

    date.intervals <- c(day0.interval,
                        day8.interval,
                        day16.interval,
                        day32.interval,
                        day64.interval)
    
    moisture.site <- mysite %>% 
      filter(date.time %in% date.intervals) %>% 
      mutate(Day = case_when(date.time %in% day0.interval ~ 0,
                             date.time %in% day8.interval ~ 8,
                             date.time %in% day16.interval ~ 16,
                             date.time %in% day32.interval ~ 32,
                             date.time %in% day64.interval ~ 64))
    
    if (i == 1) {
      df.window <- moisture.site
    }
    
    if (i > 1) {
      df.window <- rbind(df.window,moisture.site)
    }
    
    rm(day0,day8,day16,day32,day64,
       day0.interval,day8.interval,day16.interval,day32.interval,day64.interval,
       date.intervals,mysite,moisture.site)
  }
  
  ## Save data
  
  write.table(df.window, paste0(mydir.data,"Soil RH ",duration.hr," h, ",daytime.hr,"-00.csv"), sep = ",", row.names = F)
  save(df.window, file = paste0(mydir.data,"Soil_RH_",duration.hr,"_h_",daytime.hr,"_00.RData"))

  #================================================
  # Calculate the average soil moisture over the period
  
  avg.window <- df.window %>% 
    group_by(Day, Site, Region, Farm, Plot, site_ID, treatment) %>% 
    summarize(Avg.Moisture = mean(Moisture, na.rm = T))
  
  ## Calculate ratio for each Drought replicate as Ratio = Di / Cmean for every sampling day and every site
  
  # Calculate Cmean (baseline) for each site
  
  baseline.window <- avg.window %>% 
    filter(treatment == "C") %>% 
    group_by(Day, Site, Region, Farm, site_ID) %>% 
    summarize(control_mean = mean(Avg.Moisture, na.rm = T))
  
  # Join datasets 
  
  moisture.window <- left_join(avg.window, baseline.window)
  
  # Calculate ratio
  
  moisture.window$Ratio <- moisture.window$Avg.Moisture/moisture.window$control_mean
  
  # Save data
  
  write.table(moisture.window, paste0(mydir.data,"Avg soil RH ",duration.hr," h, ",daytime.hr,"-00.csv"), sep = ",", row.names = F)
  save(moisture.window, file = paste0(mydir.data,"Avg_soil_RH_",duration.hr,"_h_",daytime.hr,"_00.RData"))

  # Print min and max
  
  minR <- min(moisture.window$Ratio, na.rm = T)
  maxR <- max(moisture.window$Ratio, na.rm = T)
  
  cat("Min: ",minR," in site ",moisture.window$Site[which(moisture.window$Ratio == minR)],"\n")
  cat("Max: ",maxR," in site ",moisture.window$Site[which(moisture.window$Ratio == maxR)],"\n")
  
  #================================================
  # Plot ratios per sampling day
  
  moisture.window$Day.Title <- paste("Day",moisture.window$Day)
  moisture.window$Day.Title <- as.factor(moisture.window$Day.Title)
  moisture.window$Day.Title <- factor(moisture.window$Day.Title, 
                                      levels = c("Day 0", "Day 8", "Day 16", "Day 32", "Day 64"))
  
  # Only drought treatment
  
  drought.window <- moisture.window %>% 
    filter(treatment == "D")
  
  ylab <- expression(paste(bold("Ratio"~D["i"]~":"~C["mean"])))
  
  if (daytime.hr < 10) {
    title <- paste0("Soil moisture ratio ",duration.hr, " h period, starting at 0",daytime.hr,":00")
    figure.title1 <- paste0("Soil moisture ratio ",duration.hr, " h, 0",daytime.hr,"-00 V1.png")
    figure.title2 <- paste0("Soil moisture ratio ",duration.hr, " h, 0",daytime.hr,"-00 V2.png")
    figure.title3 <- paste0("Soil resilience ",duration.hr, " h, 0",daytime.hr,"-00 V1.png")
    figure.title4 <- paste0("Soil resilience ",duration.hr, " h, 0",daytime.hr,"-00 V2.png")
  }
  
  if (daytime.hr >= 10) {
    title <- paste0("Soil moisture ratio ",duration.hr, " h period, starting at ",daytime.hr,":00")
    figure.title1 <- paste0("Soil moisture ratio ",duration.hr, " h, ",daytime.hr,"-00 V1.png")
    figure.title2 <- paste0("Soil moisture ratio ",duration.hr, " h, ",daytime.hr,"-00 V2.png")
    figure.title3 <- paste0("Soil resilience ",duration.hr, " h, ",daytime.hr,"-00 V1.png")
    figure.title4 <- paste0("Soil resilience ",duration.hr, " h, ",daytime.hr,"-00 V2.png")
  }
  
  g1 <- ggplot(drought.window, aes(x = Site, y = Ratio)) +
    geom_boxplot(aes(color = Region, fill = Region),
                 alpha = 0.5) +
    facet_wrap(~ Day.Title, nrow = 5) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 12),
          axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
          axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 12),
          legend.title = element_text(size = 14, face = "bold"),
          panel.spacing = unit(0.4, "cm")) +
    labs(x = "Day", y = ylab, title = title) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50")
  
  if (!is.null(ylim)) {
    g1 <- g1 + scale_y_continuous(limits = ylim)
  }
  
  ggsave(paste0(mydir,figure.title1),
         g1, width = 18, height = 24, units = "cm")
  
  ####
  
  g2 <- ggplot(drought.window, aes(x = as.factor(Day), y = Ratio)) +
    geom_boxplot(aes(color = Region, fill = Region),
                 alpha = 0.5) +
    facet_grid(Region ~ Farm) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 12),
          axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
          axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 12),
          legend.position = "none",
          panel.spacing = unit(0.4, "cm")) +
    labs(x = "Day", y = ylab, title = title) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50")
  
  if (!is.null(ylim)) {
    g2 <- g2 + scale_y_continuous(limits = ylim)
  }
  
  ggsave(paste0(mydir,figure.title2),
         g2, width = 24, height = 18, units = "cm")
  
  #================================================
  # Plot ratios per sampling day regardless of site (see Willson's script)
  
  ylab2 <- expression(paste(bold("Ratio"~"Plot"[italic("i")]~":"~"Plot"[italic("control mean")])))
  
  g3 <- ggplot(data = moisture.window, 
                aes(x = Day, y = Ratio, 
                    color = treatment, fill = treatment)) + 
    geom_point() + 
    geom_smooth(method = "loess", alpha = 0.2) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 12),
          axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
          axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
    labs(y = ylab2, title = title)
  
  ggsave(paste0(mydir,figure.title3),
         g3, width = 8, height = 5)
  
  #================================================
  # Plot ratios per sampling day for each region and farm
  
  g4 <- ggplot(data = moisture.window, 
                aes(x = Day, y = Ratio, 
                    color = treatment, fill = treatment)) + 
    geom_point() + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 12),
          axis.title.x = element_text(size = 14, face = "bold", margin = margin(15,0,0,0)),
          axis.title.y = element_text(size = 14, face = "bold", margin = margin(0,15,0,0)),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          strip.text = element_text(size = 12),
          panel.spacing = unit(0.4, "cm")) +
    geom_smooth(method = "loess", alpha = 0.2) +
    facet_grid(Region ~ Farm) +
    labs(y = ylab2, title = title)
  
  g4
  
  ggsave(paste0(mydir,figure.title4),
         g4, width = 24, height = 18, units = "cm")
  
  
  # Return objects
  
  return(list(df.window,moisture.window,g1,g2,g3,g4))
}