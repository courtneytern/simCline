### Make basic plots for DEST samples
### Final data has columsn: sampleId, country, city, collectionDate, lat, long, season, nFlies, locality, type (inbred/pooled), continent
### Adapted from Alan Bergland, Oct 3, 2018; updated Feb 2020


### ijob -c1 -p standard -A berglandlab
### module load gcc/7.1.0  openmpi/3.1.4 R/3.6.0; R


### libraries
library(data.table)
library(ggplot2)
library(maps)


### set working directory
# setwd("/scratch/cat7ep/simCline/biosampleresults")
setwd("~/Downloads/GitHub/simCline/biosampleresults")

### load data
samps <- fread("./concatenated.csv")
# samps <- rbind(samps, fread("./concatenatedEuro.csv"))

### time plot
### find sites with multiple time points
samps.ag <- samps[,list(nSamps=length(row),
                        nTime=length( unique(samps[,c('year','month','day')]) ),
                        lat=mean(samps$lat[which(samps$lat!="NA")]),
                        long=mean(samps$long[which(samps$long!="NA")]) ),
                  list(state, year, country) ]

setkey(samps.ag, country, year)
setkey(samps, country, year)


### world map plot

worldData <- as.data.table(map_data("world"))

samps.ag.ag <- samps.ag[,list(n=sum(nTime), lat=samps$lat, 
                                            long=samps$long,
                                            pool.indiv=samps$`p/i`), 
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
             aes(x=long, y=lat, color=pool.indiv), size=2, alpha=.5) +
  xlab("Longitude") + ylab("Latitude") + scale_fill_manual(values="black")
world

ggsave(world, file="./metadata/worldPlot.pdf", height=4, width=6)


## north america
# min.lat.na <- 25
# max.lat.na <- 50
# min.long.na <- -130
# max.long.na <- -65


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