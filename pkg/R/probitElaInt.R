probitElaInt <- function( allCoef, allXVal, xPos, xBound, 
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
  # vector of probabilities of y=1 for each interval
  xBeta <- calcXBetaInt( allCoef, allXVal, xPos )
  checkXBeta( xBeta )
  phiVec <- pnorm( xBeta )
  # vector of shares of observations in each interval
  shareVec <- calcSharesInt( allXVal, xPos )
  # weights
  weights <- elaIntWeights( shareVec )
  # calculation of the semi-elasticity
  semEla <- lpmElaInt( phiVec, shareVec, xBound )
  ### calculation of its standard error
  # partial derivative of the semi-elasticity 
  # w.r.t. all estimated coefficients
  derivCoef <- probitElaIntDeriv( allCoef, allXVal, xPos, xBound )
  # standard error of the (average) semi-elasticity
  semElaSE <- drop( sqrt( t( derivCoef ) %*% allCoefVcov %*% derivCoef ) )
  # prepare object that will be returned
  result <- c( semEla[1], stdEr = semElaSE )
  return( result )
}
