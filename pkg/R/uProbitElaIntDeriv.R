uProbitElaIntDeriv <- function( allCoef, allXVal, xPos, xBound ){
  # number of coefficients
  nCoef <- length( allCoef )
  # number of intervals
  nInt <- length( xPos ) 
  # vector of probabilities of y=1 for each interval
  xBeta <- calcXBetaInt( allCoef, allXVal, xPos )
  # vector of shares of observations in each interval
  shareVec <- calcSharesInt( allXVal, xPos )
  # weights
  weights <- elaIntWeights( shareVec )
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
  d <- rep( 0, nCoef )
  for( m in 1:( nInt - 1 ) ){
    d <- d + weights[m] * gradM[ , m ]
  }
  return( d )
}
