Haplarithmisis_po2 <- function(Father,Mother,Gtypes,Family,EmbryoID,outPath){
	
	print("Parent of origin analysis...")
	PatGtype <- Gtypes[,Father]
	MatGtype <- Gtypes[,Mother]

	ChildGtype <- Gtypes[,EmbryoID]
	P <- matrix(0,nrow = nrow(Gtypes), ncol = 1)
	
	P[MatGtype== "AA" & PatGtype == "BB" & ChildGtype == "AA"] <- 1
	P[MatGtype== "AA" & PatGtype == "BB" & ChildGtype == "BB"] <- -1
	P[MatGtype== "BB" & PatGtype == "AA" & ChildGtype == "AA"] <- -1
	P[MatGtype== "BB" & PatGtype == "AA" & ChildGtype == "BB"] <- 1
	P[MatGtype== "AB" & PatGtype == "AA" & ChildGtype == "BB"] <- 0.5
	P[MatGtype== "AA" & PatGtype == "AB" & ChildGtype == "BB"] <- -0.5
	P[MatGtype== "BB" & PatGtype == "AB" & ChildGtype == "AA"] <- -0.5
	P[MatGtype== "AB" & PatGtype == "BB" & ChildGtype == "AA"] <- 0.5

        colnames(P) <-EmbryoID
        Ps <- P
        print(EmbryoID)

	dataPo <- data.frame(Gtypes[,c("Names","Chr","Position")],Ps,stringsAsFactors=FALSE,check.names=F)
	dataOut <- dataPo
	colnames(dataOut) <- colnames(dataPo)
	write.table(dataOut,paste(outPath,paste(Family,".poo",sep=""),sep="/"),quote=F,sep="\t",col.names=T,row.names=F)
	print(paste("Parent of origin origin file (",Family,".poo) is written on: ",outPath,sep=""))
	return(dataPo)
} #end PO function