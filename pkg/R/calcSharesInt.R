calcSharesInt <- function( allXVal, xPos ){
  # number of intervals
  nInt <- length( xPos ) 
  # vector of shares of observations in each interval
  shareVec <- rep( NA, nInt )
  for( i in 1:nInt ){
    if( xPos[i] != 0 ) {
      shareVec[ i ] <- allXVal[ xPos[ i ] ] 
    }
  }
  shareVec[ xPos == 0 ] <- 1 - sum( shareVec[ xPos != 0 ] )
  return( shareVec )
}
