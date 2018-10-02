lpmElaInt <- function( xCoef, xShares, xBound, 
  allCoefVcov = NULL ){
  nInt <- length( xCoef )
  if( nInt < 2 || !is.vector( xCoef ) ) {
    stop( "argument 'xCoef' must be a vector with at least two elements" )
  }
  if( length( xShares ) != nInt ) {
    stop( "arguments 'xCoef' and 'xShares' must be vectors of the same length" )
  }
  if( any( xShares < 0 ) ) {
    stop( "all shares in argument 'xShares' must be non-negative" )
  }
  if( abs( sum( xShares ) - 1 ) > 0.015 ) {
    stop( "the shares in argument 'xShares' must sum to one" )
  }
  # check 'xBound' and replace infinite values
  xBound <- elaIntBounds( xBound, nInt )
  # check and prepare allCoefVcov
  allCoefVcov <- prepareVcov( allCoefVcov, length( xCoef ), NA, xMeanSd = NULL )
  # weights
  weights <- elaIntWeights( xShares )
  # semi-elasticities 'around' each inner boundary and their weights
  semElas <- rep( NA, nInt - 1 )
  for( m in 1:(nInt-1) ){
    semElas[m] <- 2 * ( xCoef[ m+1 ] - xCoef[ m ] ) * xBound[ m+1 ] /
      ( xBound[m+2] - xBound[m] )
  }
  # (average) semi-elasticity
  semElaAvg <- sum( semElas * weights )
  # derivatives of the (average) semi-elasticity wrt the coefficients
  derivCoef <- rep( NA, nInt )
  derivCoef[1] <- 
    -2 * weights[1] * xBound[2] / ( xBound[3] - xBound[1] )
  derivCoef[nInt] <- 
    2 * weights[nInt-1] * xBound[nInt] / ( xBound[nInt+1] - xBound[nInt-1] ) 
  if( nInt > 2 ) {
    for( n in 2:( nInt-1 ) ) {
      derivCoef[n] <- 
        2 * weights[n-1] * xBound[n] / ( xBound[n+1] - xBound[n-1] ) -
        2 * weights[n]   * xBound[n+1] / ( xBound[n+2] - xBound[n] )
    }
  }
  # standard error of the (average) semi-elasticity
  semElaSE <- drop( sqrt( t( derivCoef ) %*% allCoefVcov %*% derivCoef ) )
  # prepare object that will be returned
  result <- c( semEla = semElaAvg, stdEr = semElaSE )
  return( result )
}
