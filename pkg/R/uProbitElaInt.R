uProbitElaInt <- function( allCoef, allXVal, xPos, xBound, 
  allCoefVcov = NULL ){
  # number of coefficients
  nCoef <- length( allCoef )
  # number of intervals
  nInt <- length( xPos ) 
  # checking arguments
  if( length( allXVal ) != nCoef ) {
    stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
  }
  checkXPos( xPos, minLength = 2, maxLength = nCoef, 
    minVal = 0, maxVal = nCoef, requiredVal = 0 )
  if( any( allXVal[ xPos ] < 0 ) ) {
    stop( "all elements of argument 'allXVal'",
      " that are indicated by argument 'xPos'",
      " (i.e., the shares of observations in each interval)",
      " must be non-negative" )
  }
  if( sum( allXVal[ xPos ] > 1 ) ) {
    stop( "the sum of the elements of argument 'allXVal'",
      " that are indicated by argument 'xPos'",
      " (i.e., the shares of observations in each interval)",
      " must not be larger than one" )
  }
  # check 'xBound' and replace infinite values
  xBound <- elaIntBounds( xBound, nInt )
  # check and prepare allCoefVcov
  allCoefVcov <- prepareVcov( allCoefVcov, nCoef, xPos, xMeanSd = NULL )
  # vector of probabilities of y=1 for each interval and
  # vector of shares of observations in each interval
  xBeta <- shareVec <- rep( NA, nInt )
  for( i in 1:nInt ){
    allXValTemp <- replace( allXVal, xPos, 0 )
    if( xPos[i] != 0 ) {
      allXValTemp[ xPos[i] ] <- 1
      shareVec[ i ] <- allXVal[ xPos[ i ] ] 
    }
    xBeta[ i ] <- sum( allCoef * allXValTemp )
  }
  shareVec[ xPos == 0 ] <- 1 - sum( shareVec[ xPos != 0 ] )
  checkXBeta( xBeta )
  phiVec <- pnorm( xBeta )
  # weights
  weights <- elaIntWeights( shareVec )
  # calculation of the semi-elasticity
  semEla <- uLinElaInt( phiVec, shareVec, xBound )
  ### calculation of its standard error
  # partial derivatives of each semi-elasticity around each boundary
  # w.r.t. all estimated coefficients
  gradM <- matrix( 0, nCoef, nInt - 1 )
  gradPhiVec <- dnorm( xBeta )
  for( m in 1:( nInt - 1 ) ) {
    gradM[ -xPos, m ] <- 2 * ( gradPhiVec[m+1] - gradPhiVec[m] ) * 
      allXVal[ -xPos ] * xBound[m+1] / ( xBound[m+2] - xBound[m] )
    gradM[ xPos[m], m ] <- - 2 * gradPhiVec[m] * xBound[m+1] / 
      ( xBound[m+2] - xBound[m] )
    gradM[ xPos[m+1], m ] <- 2 * gradPhiVec[m+1] * xBound[m+1] / 
      ( xBound[m+2] - xBound[m] )
  }
  # partial derivative of the semi-elasticity 
  # w.r.t. all estimated coefficients
  derivCoef <- rep( 0, nCoef )
  for( m in 1:( nInt - 1 ) ){
    derivCoef <- derivCoef + weights[m] * gradM[ , m ]
  }
  # standard error of the (average) semi-elasticity
  semElaSE <- drop( sqrt( t( derivCoef ) %*% allCoefVcov %*% derivCoef ) )
  # prepare object that will be returned
  result <- c( semEla[1], stdEr = semElaSE )
  return( result )
}
