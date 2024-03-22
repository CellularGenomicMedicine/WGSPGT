# Functions to calculate running median for a given window size (k = 25) obtained from Nik-Zainal et al. 2012
# Reference:
# Nik-Zainal S., Van Loo P., Wedge D.C., Alexandrov L.B., Greenman C.D., Lau K.W., Raine K., Jones D., Marshall J., Ramakrishna M.et al... The life history of 21 breast cancers. Cell. 2012; 149:994–1007.


getMad <- function(x, k = 25){
  
  #Remove observations that are equal to zero; are likely to be imputed, should not contribute to sd:
  x <- x[x!=0]
  
  #Calculate runMedian  
  runMedian <- medianFilter(x, k)
  
  dif <- x - runMedian
  SD <- mad(dif)
 
	return(SD)
}

# medianFilter - Function to calculate the running median for a given window size (k)
medianFilter <- function(x, k){
  n <- length(x)
  filtWidth <- 2 * k + 1

  #Make sure filtWidth does not exceed n
  if(filtWidth > n){
    if(n == 0){
      filtWidth <- 1
    }else if(n%%2 == 0){
      #runmed requires filtWidth to be odd, ensure this:
      filtWidth <- n - 1
    }else{
      filtWidth <- n
    }
  }

  runMedian <- runmed(x, k = filtWidth, endrule = "median")

  return(runMedian)

}
