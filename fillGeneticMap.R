## keep R from writing numbers as 3e06 for example
options(scipen = 999)

mapFile <-commandArgs(trailingOnly = TRUE)[1]

map<-read.table(mapFile,as.is=T,h=T)

linearPoints<-c(0,map[map$Genetic_Map.cM.!=0,"Genetic_Map.cM."])
linearPositions<-c(0,map[map$Genetic_Map.cM.!=0,"position"])

lowerPos<-c(0,map[map$Genetic_Map.cM.!=0,"position"])[1:(length(linearPositions)-1)]
upperPos<-c(0,map[map$Genetic_Map.cM.!=0,"position"])[2:length(linearPositions)]

lowerMorgan<-c(0,map[map$Genetic_Map.cM.!=0,"Genetic_Map.cM."])[1:(length(linearPositions)-1)]
upperMorgan<-c(0,map[map$Genetic_Map.cM.!=0,"Genetic_Map.cM."])[2:length(linearPositions)]

lowerUpper<-data.frame(lower=lowerPos,upper=upperPos,lowerMorgan=lowerMorgan,upperMorgan=upperMorgan,stringsAsFactors = F)

getValue<-function(pos){

  ## the interval that it is - meaning the one it is least above and least below (perhaps do sum of duff)
  ##which.min(rowSums(cbind(abs(pos-lowerUpper$lower),abs(pos-lowerUpper$upper))))
  if(pos>lowerUpper[which.max(lowerUpper$upper),"upper"]){

    return(lowerUpper[which.max(lowerUpper$upper),"upperMorgan"])
  }

  closeUpper<-lowerUpper[lowerUpper$upper>=pos,"upper"][which.min(lowerUpper[lowerUpper$upper>=pos,"upper"]-pos)]
  index<-which(lowerUpper$upper==closeUpper)

  morgan<-(pos-lowerUpper[index,"lower"])*((lowerUpper[index,"upperMorgan"]-lowerUpper[index,"lowerMorgan"])/(lowerUpper[index,"upper"]-lowerUpper[index,"lower"])) + lowerUpper[index,"lowerMorgan"]
  ##print(morgan)
  return(morgan)

}


centiMorgans<-unlist(sapply(map$position,getValue))

map[,3]<-centiMorgans

## differnece in position
map[,4]<-c(0,diff(map[,1]))
## difference in gpos
map[,5]<-c(0,diff(map[,3]))

map[,2]<-(map[,5]/map[,4])*1000000

map[1,2]<-map[2,2]
colnames(map)<-c("position", "COMBINED_rate(cM/Mb)", "Genetic_Map(cM)")

write.table(map[,1:3],paste("filledOut2_",mapFile,sep=""),col.names = T,row.names=F,quote=F)
