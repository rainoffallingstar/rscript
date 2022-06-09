names <- installed.packages()
nnames <- names[,1]
write.table(nnames,"packages_backup/installed.csv",row.names=FALSE,col.names=TRUE,sep="/n")
