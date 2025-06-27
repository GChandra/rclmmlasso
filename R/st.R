# soft-thresholding algorithm
st <- function(x, alpha){
  st_x <- rep(0, length(x))
  
  for (i in 1:length(x)){
    if (x[i] > alpha){ st_x[i] <- x[i] - alpha }
    else if (x[i] < -alpha){ st_x[i] <- x[i] + alpha }
    else { st_x[i] <- 0 }
  }
  
  return ( st_x )
}
