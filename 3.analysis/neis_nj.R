##This script will calculate nei's genetic distance 

# module load gcc/7.1.0  openmpi/3.1.4  R
# R

library(foreach)
library(ape)
library(data.table)

# summTable includes pop,chr,pos,ad,rd
summTable<- fread("/scratch/cat7ep/simCline/data/allPops_adrd.txt")
summTable[,freq:=ad/(ad+rd)]

pops <- unique(summTable$population)
setkey(summTable, population)

nei.dt <- foreach(p1=pops, .combine="rbind")%do%{
  foreach(p2=pops, .combine="rbind")%do%{
    # p1 <- pops[,1]; p2 <- pops[,2]
    
    print(paste("p1=",p1,"p2=",p2))
    if(p1!=p2){
      tmp <- summTable[J(c(p1, p2))]
      tmp.wide <- dcast(tmp, chr+pos~population, value.var="freq")
      
      #print("setnames")
      setnames(tmp.wide, c("chr","pos","p1", "p2"))
      tmp.wide <- tmp.wide[!is.na(p1) & !is.na(p2)]
      print("Nei")
      Nei <- -log(sum(tmp.wide$p1 * tmp.wide$p2) / (sqrt(sum(tmp.wide$p1^2)) * sqrt(sum(tmp.wide$p2^2))))
    
    } else{
      print("Nei=0")
      Nei <- 0
    }
       
    print("make data table")
    data.table(pop1=p1, pop2=p2, nei=Nei)
  }
}
fwrite(nei.dt,"/scratch/cat7ep/simCline/data/nei.dt.txt")

nei.dt<- fread("/scratch/cat7ep/simCline/data/nei.dt.txt")
#nei.dt<- fread("~/Downloads/nei.dt.txt")
Da_mat <- matrix(nei.dt$nei, nrow=length(pops), ncol=length(pops))
  colnames(Da_mat)<- unique(nei.dt$pop1)
  rownames(Da_mat)<- unique(nei.dt$pop1)
tree <- njs(Da_mat)
plot.phylo(tree,cex=0.4,no.margin=T,node.pos=2)
