x=read.table("scale_Ic.txt",head=FALSE)

library(fpc)
ds <- dbscan(x,0.5^(1/2),MinPts=40)
write.table(ds[1],file ="data_temp_Ic.txt",row.names = F, col.names=F, quote = F)

