urbinElaInt <- function( allCoef, allXVal, xPos, xBound, model,
  allCoefVcov = NULL ){
  
  if( model != "probit" ) {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  
  # number of coefficients
  nCoef <- length( allCoef )
  # Check position vector
  checkXPos( xPos, minLength = 2, maxLength = nCoef, 
    minVal = 0, maxVal = nCoef, requiredVal = 0 )
  # number of intervals
  nInt <- length( xPos ) 
  # Check x values
  if( length( allXVal ) != nCoef ) {
    stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
  }
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
  # prepare calculation of semi-elasticity 
  # vector of probabilities of y=1 for each interval
  xBeta <- calcXBetaInt( allCoef, allXVal, xPos )
  checkXBeta( xBeta )
  phiVec <- pnorm( xBeta )
  # vector of shares of observations in each interval
  shareVec <- calcSharesInt( allXVal, xPos )
  # weights
  weights <- elaIntWeights( shareVec )
  # calculation of the semi-elasticity
  semEla <- lpmElaInt( phiVec, shareVec, xBound )[1]
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
  # partial derivatives of semi-elasticities wrt coefficients
  derivCoef <- rep( 0, nCoef )
  for( m in 1:( nInt - 1 ) ){
    derivCoef <- derivCoef + weights[m] * gradM[ , m ]
  }
  # if argument allXVal has attribute 'derivOnly',
  # return partial derivatives only (for testing partial derivatives)
  if( "derivOnly" %in% names( attributes( allXVal ) ) ) {
    return( derivCoef )
  }
  # approximate standard error of the semi-elasticity
  semElaSE <- drop( sqrt( t( derivCoef ) %*% allCoefVcov %*% derivCoef ) )
  # create object that will be returned
  result <- c( semEla = unname( semEla ), stdEr = semElaSE )
  return( result )
}
