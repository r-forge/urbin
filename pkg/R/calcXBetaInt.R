calcXBetaInt <- function( allCoef, allXVal, xPos ){
  # number of coefficients
  nCoef <- length( allCoef )
  # number of intervals
  nInt <- length( xPos ) 
  # vector of probabilities of y=1 for each interval
  xBeta <- rep( NA, nInt )
  for( i in 1:nInt ){
    allXValTemp <- replace( allXVal, xPos, 0 )
    if( xPos[i] != 0 ) {
      allXValTemp[ xPos[i] ] <- 1
    }
    xBeta[ i ] <- sum( allCoef * allXValTemp )
  }
  return( xBeta )
}
