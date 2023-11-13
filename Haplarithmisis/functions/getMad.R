getMad <- function(x,k=25){
  
  #Remove observations that are equal to zero; are likely to be imputed, should not contribute to sd:
  x <- x[x!=0]
  
  #Calculate runMedian  
  runMedian <- medianFilter(x,k)
  
  dif <- x-runMedian
  SD <- mad(dif)
 
	return(SD)
}
medianFilter <- function(x,k){
  n <- length(x)
  filtWidth <- 2*k + 1

  #Make sure filtWidth does not exceed n
  if(filtWidth > n){
    if(n==0){
      filtWidth <- 1
    }else if(n%%2 == 0){
      #runmed requires filtWidth to be odd, ensure this:
      filtWidth <- n - 1
    }else{
      filtWidth <- n
    }
  }

  runMedian <- runmed(x,k=filtWidth,endrule="median")

  return(runMedian)

}
