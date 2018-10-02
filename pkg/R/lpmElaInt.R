lpmElaInt <- function( allCoef, allXVal, xBound, xPos,
  allCoefVcov = NULL ){
  
  # number of coefficients
  nCoef <- length( allCoef )
  
  # Check position vector
  checkXPos( xPos, minLength = 2, 
    maxLength = length( allCoef ) + any( xPos == 0 ), 
    minVal = 0, maxVal = length( allCoef ) )
  xCoef <- rep( 0, length( xPos ) )
  for( i in 1:length( xPos ) ) {
    if( xPos[i] != 0 ) {
      xCoef[i] <- allCoef[ xPos[i] ]
    }
  }
  shareVec <- calcSharesInt( allXVal, xPos )
  if( "derivOnly" %in% names( attributes( allXVal ) ) ) {
    attr( shareVec, "derivOnly" ) <- 1 
  }
  # number of intervals
  nInt <- length( xPos ) 
  # Check x values
  if( length( allXVal ) != nCoef ) {
    stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
  }
  # check 'xBound' and replace infinite values
  xBound <- elaIntBounds( xBound, nInt )
  # check and prepare allCoefVcov
  allCoefVcov <- prepareVcov( allCoefVcov, nCoef, NA, xMeanSd = NULL )
  # weights
  weights <- elaIntWeights( shareVec )
  # semi-elasticities 'around' each inner boundary and their weights
  semElas <- rep( NA, nInt - 1 )
  for( m in 1:(nInt-1) ){
    semElas[m] <- 2 * ( xCoef[ m+1 ] - xCoef[ m ] ) * xBound[ m+1 ] /
      ( xBound[m+2] - xBound[m] )
  }
  # (average) semi-elasticity
  semElaAvg <- sum( semElas * weights )
  # derivatives of the (average) semi-elasticity wrt the coefficients
  derivCoef <- rep( 0, nCoef )
  derivCoef[ xPos[1] ] <- 
    -2 * weights[1] * xBound[2] / ( xBound[3] - xBound[1] )
  derivCoef[ xPos[nInt] ] <- 
    2 * weights[nInt-1] * xBound[nInt] / ( xBound[nInt+1] - xBound[nInt-1] ) 
  if( nInt > 2 ) {
    for( n in 2:( nInt-1 ) ) {
      derivCoef[ xPos[n] ] <- 
        2 * weights[n-1] * xBound[n] / ( xBound[n+1] - xBound[n-1] ) -
        2 * weights[n]   * xBound[n+1] / ( xBound[n+2] - xBound[n] )
    }
  }
  # if argument allXVal has attribute 'derivOnly',
  # return partial derivatives only (for testing partial derivatives)
  if( "derivOnly" %in% names( attributes( shareVec ) ) ) {
    return( derivCoef )
  }
  # standard error of the (average) semi-elasticity
  semElaSE <- drop( sqrt( t( derivCoef ) %*% allCoefVcov %*% derivCoef ) )
  # prepare object that will be returned
  result <- c( semEla = semElaAvg, stdEr = semElaSE )
  return( result )
}
