#######################################################################
### Analyses and graphics for stability metrics at 10 km resolution ###
#######################################################################

### Hannah White 11.03.2019
### Edited 28.03.2019

library(ggplot2)
library(cowplot)
library(raster)
library(viridis)
library(scales)
library(rasterVis)


## load stability measures

stability10km <- read.csv('J:\\Postdoc Grassland Resilience\\MODIS6\\Resilience\\stability10km.csv', header = TRUE)

## read in location of terrestrial hectads

corine10.terr <- read.csv('J:/Postdoc Grassland Resilience/LandCoverData/corine10.terr.csv', header = TRUE)

coords <- corine10.terr[,2:3]
coords[,1] <- coords[,1] + 5000
coords[,2] <- coords[,2] + 5000

stability10km <- merge(stability10km, coords, by.x = c('eastings', 'northings'), by.y = c('east', 'north'), all.x = FALSE, all.y = TRUE)

## Change recovery time to days rather than units

stability10km$evi.recmed <- stability10km$evi.recmed*8
stability.red <- stability10km[, c(1:3, 5, 9, 11, 13:14)]



### Maps

dfr <- rasterFromXYZ(stability.red)  #Convert first two columns as lon-lat and third as value    

map_colors<-colorRampPalette(viridis(12))

tiff('J:\\Postdoc Grassland Resilience\\Writing\\MODIS_resilience_metrics\\Graphics\\maps10km.tiff', height = 1000, width = 1000)
par(oma = c(0,0,0,3))
plot(dfr, box = FALSE, axes = FALSE, 
     main = c('Variability', 'Resistance', 'Recovery Time', 'Recovery Rate', 'Longest Return', 'Slowest Return'),
     cex.main = 2, legend.width = 2, 
     axis.args=list(cex.axis=1.5))
dev.off()

# variability

p.var <- ggplot(stability.red, aes(x = eastings, y = northings, col = evi.var)) + geom_point(pch = 15) + coord_equal()

p.var <- p.var + xlab('') + ylab('')


# magnitude

p.mag <- ggplot(stability.red, aes(x = eastings, y = northings, col = evi.magmed)) + geom_point(pch = 15) + coord_equal()

p.mag <- p.mag + xlab('') + ylab('')



# recovery time

p.rec <- ggplot(stability.red, aes(x = eastings, y = northings, col = evi.recmed)) + geom_point(pch = 15) + coord_equal()

p.rec <- p.rec + xlab('') + ylab('')


# recovery rate

p.rate <- ggplot(stability.red, aes(x = eastings, y = northings, col = evi.recrate.med)) + geom_point(pch = 15) + coord_equal()

p.rate <- p.rate + xlab('') + ylab('')


# length

p.long <- ggplot(stability.red, aes(x = eastings, y = northings, col = evi.long)) + geom_point(pch = 15) + coord_equal()

p.long <- p.long + xlab('') + ylab('')


# slow

p.slow <- ggplot(stability.red, aes(x = eastings, y = northings, col = evi.slow)) + geom_point(pch = 15) + coord_equal()


p.slow <- p.slow + xlab('') + ylab('')


plot_grid(p.var, p.mag, p.rec, p.rate, p.long, p.slow, labels = c('a)', 'b)', 'c)', 'd)', 'e)', 'f)'),
          nrow = 3, ncol = 2)





### Cross-correlations



pairs(stability.red[,3:8])


cor(stability.red[,3:8], use = 'na.or.complete', method = 'spearman') # use Spearman's as many correlations appear non linear (see pair graph)

pairs(stability.red[,3:8]) # low correlations, so seem to capturing different aspects of stability

panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1))
  txt <- as.character(format(cor(x, y, use = 'na.or.complete', method = 'spearman'), digits=2))
  text(0.5, 0.5, txt, cex = 4)
}

tiff('J:\\Postdoc Grassland Resilience\\Writing\\MODIS_resilience_metrics\\Graphics\\correlations10km.tiff', height = 1000, width = 1000)
pairs(stability.red[,3:8], labels = c('variability', 'resistance', 'recovery time', 'recovery rate', 'longest return', 'slowest return'), 
      upper.panel = panel.cor) ## need to change colour of some of them as non-sig
dev.off()

library(Hmisc)

rcorr(cbind(stability.red[,3], stability.red[,4], stability.red[,5], stability.red[,6], stability.red[,7], stability.red[,8]),  type = 'spearman') ## all significant



### Timing of events


load('D:\\MODIS6\\MODIS6_18.2.19\\modisv6_10km_all.RData')

rm(evi.mean, evi.median, evi.ncell, evi.sd, evi.mad,
   ndvi.mean, ndvi.median, ndvi.ncell, ndvi.sd, ndvi.mad)


load('J:\\Postdoc Grassland Resilience\\MODIS6\\Resilience\\resilience_measures_10km.RData')



### ggplot histogram
library(ggplot2)
library(scales)
library(lubridate)

#### All events

## stack values
dd <- evi.dates[which(!is.na(evi.dates))]
ddnum <- as.numeric(dd)

dd.df <- data.frame(dd, ddnum)


## eastings and northings
deast <- eastings[which(!is.na(evi.dates[,,1]))]

for (i in 2:dim(evi.dates)[3]){
  east.next <- eastings[which(!is.na(evi.dates[,,i]))]
  deast <- c(deast, east.next)
}


dnorth <- northings[which(!is.na(evi.dates[,,1]))]

for (j in 2:dim(evi.dates)[3]){
  north.next <- northings[which(!is.na(evi.dates[,,j]))]
  dnorth <- c(dnorth, north.next)
}

dd.df$eastings <- deast
dd.df$northings <- dnorth

coords$comb <- paste(coords$east, coords$north, sep = '.')

dd.df$comb <- paste(dd.df$eastings, dd.df$northings, sep = '.')


dd.df <- merge(dd.df, coords, by = 'comb', all.x = FALSE, all.y = TRUE) # extracts those with > 50% land


bin = 30

p <- ggplot(dd.df, aes(dd, ..count..))
p <- p + geom_histogram(data = dd.df, binwidth = bin)

# The numeric data is treated as a date,
# breaks are set to an interval equal to the binwidth,
# and a set of labels is generated and adjusted in order to align with bars
p <- p + scale_x_date(breaks = waiver(),
                      date_breaks = '1 year',
                      labels = date_format("%Y-%m"),
                      limits = c(as.Date("2003-01-01"), 
                                 as.Date("2019-01-25")))
p <- p + xlab('')
p <- p + scale_fill_manual(values = c(rep(c('grey20', 'grey80'), 7), 'grey20'))

# from here, format at ease
p <- p + theme(axis.text.x  = element_text(angle=45, hjust = 1, vjust = 1, size = 14), axis.text.y = element_text(size = 11), legend.position ='none',
               panel.background = element_blank(), panel.grid = element_blank(), axis.line = element_line(colour = 'black')) + ylab('')

p







#### Largest events

evi.dates.big <- evi.dates[,,1]
ddbig <- evi.dates.big[which(!is.na(evi.dates.big))]
ddnumbig <- as.numeric(ddbig)

ddbig.df <- data.frame(ddbig, ddnumbig)

deast <- eastings[which(!is.na(evi.dates[,,1]))]
dnorth <- northings[which(!is.na(evi.dates[,,1]))]

ddbig.df$east <- deast
ddbig.df$north <- dnorth

ddbig.df <- merge(ddbig.df, coords, by = c('east', 'north'), all.x = FALSE, all.y = TRUE) # extracts those with >50% land


p2 <- ggplot(ddbig.df, aes(x = ddbig, y = ..count..))
p2 <- p2 + geom_histogram(data = ddbig.df, binwidth = 30)

# The numeric data is treated as a date,
# breaks are set to an interval equal to the binwidth,
# and a set of labels is generated and adjusted in order to align with bars
p2 <- p2 + scale_x_date(breaks = waiver(),
                        date_breaks = '1 year',
                        labels = date_format("%Y-%m"),
                        limits = c(as.Date("2003-01-01"), 
                                   as.Date("2019-01-25")))
p2 <- p2 + xlab('')
p2 <- p2 + scale_fill_manual(values = c(rep(c('grey20', 'grey80'), 7), 'grey20'))

# from here, format at ease
p2 <- p2 + theme(axis.text.x  = element_text(angle=45, hjust = 1, vjust = 1, size = 14), axis.text.y = element_text(size = 11), legend.position ='none',
                 panel.background = element_blank(), panel.grid = element_blank(), axis.line = element_line(colour = 'black')) + ylab('')

p2

library(cowplot)

#tiff('J:\\Postdoc Grassland Resilience\\Writing\\MODIS_resilience_metrics\\Graphics\\events10km.tiff', width = 14, height = 7, units = 'in', res = 360)
#plot_grid(p, p2, nrow = 1, labels = c('a)', 'b)'))
#dev.off()


### Longest return rates



ddlong <- evi.longest.date[which(!is.na(evi.longest.date))]
ddnumlong <- as.numeric(ddlong)

ddlong.df <- data.frame(ddlong, ddnumlong)

deast <- eastings[which(!is.na(evi.longest.date))]
dnorth <- northings[which(!is.na(evi.longest.date))]

ddlong.df$east <- deast
ddlong.df$north <- dnorth

ddlong.df <- merge(ddlong.df, coords, by = c('east', 'north'), all.x = FALSE, all.y = TRUE) # extracts those with >50% land



p3 <- ggplot(ddlong.df, aes(x = ddlong, y = ..count..))
p3 <- p3 + geom_histogram(data = ddlong.df, binwidth = 30)

# The numeric data is treated as a date,
# breaks are set to an interval equal to the binwidth,
# and a set of labels is generated and adjusted in order to align with bars
p3 <- p3 + scale_x_date(breaks = waiver(),
                        date_breaks = '1 year',
                        labels = date_format("%Y-%m"),
                        limits = c(as.Date("2003-01-01"), 
                                   as.Date("2019-01-25")))
p3 <- p3 + xlab('')
p3 <- p3 + scale_fill_manual(values = c(rep(c('grey20', 'grey80'), 7), 'grey20'))

# from here, format at ease
p3 <- p3 + theme(axis.text.x  = element_text(angle=45, hjust = 1, vjust = 1, size = 14), axis.text.y = element_text(size = 11), legend.position ='none',
                 panel.background = element_blank(), panel.grid = element_blank(), axis.line = element_line(colour = 'black')) + ylab('')

p3

## Slowest events

ddslow <- evi.slowest.date[which(!is.na(evi.slowest.date))]
ddnumslow <- as.numeric(ddslow)

ddslow.df <- data.frame(ddslow, ddnumslow)

deast <- eastings[which(!is.na(evi.slowest.date))]
dnorth <- northings[which(!is.na(evi.slowest.date))]

ddslow.df$east <- deast
ddslow.df$north <- dnorth

ddslow.df <- merge(ddslow.df, coords, by = c('east', 'north'), all.x = FALSE, all.y = TRUE)



p4 <- ggplot(ddslow.df, aes(x = ddslow, y = ..count..))
p4 <- p4 + geom_histogram(data = ddslow.df, binwidth = 30)

# The numeric data is treated as a date,
# breaks are set to an interval equal to the binwidth,
# and a set of labels is generated and adjusted in order to align with bars
p4 <- p4 + scale_x_date(breaks = waiver(),
                        date_breaks = '1 year',
                        labels = date_format("%Y-%m"),
                        limits = c(as.Date("2003-01-01"), 
                                   as.Date("2019-01-25")))
p4 <- p4 + xlab('')
p4 <- p4 + scale_fill_manual(values = c(rep(c('grey20', 'grey80'), 7), 'grey20'))

# from here, format at ease
p4 <- p4 + theme(axis.text.x  = element_text(angle=45, hjust = 1, vjust = 1, size = 14), axis.text.y = element_text(size = 11), legend.position ='none',
                 panel.background = element_blank(), panel.grid = element_blank(), axis.line = element_line(colour = 'black')) + ylab('')

p4



tiff('J:\\Postdoc Grassland Resilience\\Writing\\MODIS_resilience_metrics\\Graphics\\events10km.tiff', width = 14, height = 14, units = 'in', res = 360)
plot_grid(p, p2, p3, p4, nrow = 2, ncol = 2, labels = c('a)', 'b)', 'c)', 'd)'))
dev.off()







### Maps


library(ggplot2)
library(raster)
library(sp)


### Maps 1km 


dfr <- rasterFromXYZ(stability10km)  #Convert first two columns as lon-lat and third as value                
plot(dfr)

tiff('J:\\Postdoc Grassland Resilience\\Writing\\MODIS_resilience_metrics\\Graphics\\maps10km.tiff', width = 14, height = 7, units = 'in', res = 360)
plot(dfr)
dev.off()

