Haplarithmisis_testmedfilt <- function(Hap,Window,Chroms){
	for(chr in Chroms){
		HapChr <- Hap[Hap$Chr==chr,4]
		#Defining chr-specific window size
		WindowChr <- round((Window*length(Hap$Chr==chr))/length(Hap$Chr==1))
		if((WindowChr%%2)==0){WindowChr<-WindowChr+1}
		#HapIntChr <- medfilt1(HapChr,WindowChr)
        HapIntChr <- runmed(HapChr,WindowChr,"median")
		if(chr==Chroms[1]){HapIntGenome <- HapIntChr}else{HapIntGenome <- c(HapIntGenome,HapIntChr)}
	   #print(chr)
   }#end chr loop
   return(HapIntGenome)
}#end function
