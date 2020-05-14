### Make basic plots for DEST samples
### Final data has columsn: sampleId, country, city, collectionDate, lat, long, season, nFlies, locality, type (inbred/pooled), continent
### Adapted from Alan Bergland, Oct 3, 2018; updated Feb 2020


### ijob -c1 -p standard -A berglandlab
### module load gcc/7.1.0  openmpi/3.1.4 R/3.6.0; R


### libraries
library(data.table)
library(gdata)
library(cowplot)
library(data.table)
library(foreach)
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)


### set working directory
setwd("/scratch/cat7ep/simCline/biosampleresults")
# setwd("~/Downloads/GitHub/simCline/biosampleresults")

### load data
samps <- fread("./concatenated.csv")

### time plot
### find sites with multiple time points
samps.ag <- samps[,list(nSamps=length(row),
                        nTime=length( unique(samps[,c('year','month','day')]) ),
                        lat=mean(samps$lat[which(samps$lat!="NA")]),
                        long=mean(samps$long[which(samps$long!="NA")]) ),
                  list(state, year, country) ]

setkey(samps.ag, country, year)
setkey(samps, country, year)


### plot multi-sample populations

# multi_sample <- ggplot() +
#   geom_line(data= samps[J(samps.ag[maxDelta>10])], aes(x=as.Date(yday, origin = as.Date("2018-01-01")), y=lat, group=locality, linetype=continent)) +
#   geom_point(data=samps[J(samps.ag[maxDelta>10])], aes(x=as.Date(yday, origin = as.Date("2018-01-01")), y=lat, group=locality, color=season)) +
#   facet_grid(.~year) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="bottom", legend.direction="vertical") +
#   scale_x_date(date_labels = "%b", limits = as.Date(c(110,355), origin = as.Date("2018-01-01"))) +
#   xlab("Collection Date") + ylab("Latitude")
# 
# ggsave(multi_sample, file="./DEST/populationInfo/multiSample.pdf")


### world map plot

worldData <- as.data.table(map_data("world"))

samps.ag.ag <- samps.ag[,list(n=sum(nTime), lat=samps$lat, 
                                            long=samps$long ), 
                          list(country) ]
# samps.ag.ag <- samps.ag[,list(n=sum(nTime), lat=mean(samps$lat[which(samps$lat!="NA")]), 
#                               long=mean(samps$long[which(samps$long!="NA")]) ), 
#                         list(state) ]

### make maps

min.lat.eu <- 35
max.lat.eu <- 55
min.long.eu <- -10
max.long.eu <- 37
# [long>=min.long.eu & long<= max.long.eu][lat>=min.lat.eu & lat<=max.lat.eu]
#[longitude>=min.long.eu & longitude<= max.long.eu][latitude>=min.lat.eu & latitude<=max.lat.eu]

# size=I((n-1)/2 + 4))
world <- 	ggplot() +
  geom_polygon(data = worldData,
               aes(x=long, y = lat, group = group), fill="lightgrey") +
  geom_point(data = samps.ag.ag,
             aes(x=long, y=lat), size=3, alpha=.5) +
  xlab("Longitude") + ylab("Latitude") + scale_fill_manual(values="black")

ggsave(world, file="./worldPlot.pdf", height=4, width=6)


## north america
min.lat.na <- 25
max.lat.na <- 50
min.long.na <- -130
max.long.na <- -65


### num. flies plot
# 
# samps[,nFlies:=as.numeric(as.character(nFlies))]
# samps[,sampleId:=as.character(sampleId)]
# samps[,sampleId:=factor(sampleId,
#                         levels=c(samps[set=="DrosEU"]$sampleId,
#                                  samps[set=="DrosRTEC"][order(nFlies)]$sampleId,
#                                  samps[set=="dpgp"]$sampleId))]
# 
# 
# nFlies.plot <- ggplot(data=samps, aes(x=sampleId, y=nFlies, color=set)) + geom_point() +
#   theme(axis.text.x=element_blank(), axis.title.x=element_blank()) +
#   ylab("Num. flies sampled")
# 
# ggsave(nFlies.plot, file="~/numFlies.pdf")