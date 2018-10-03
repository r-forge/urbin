urbinElaInt <- function( allCoef, allXVal, xPos, xBound, model,
  allCoefVcov = NULL ){
  
  # number of coefficients
  nCoef <- length( allCoef )
  # Check position vector
  checkXPos( xPos, minLength = 2, maxLength = length( allCoef ) + 1, 
    minVal = 0, maxVal = nCoef )#, requiredVal = 0 )
  # number of intervals
  nInt <- length( xPos ) 
  # Check x values
  if( length( allXVal ) != nCoef ) {
    stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
  }
  # check 'xBound' and replace infinite values
  xBound <- elaIntBounds( xBound, nInt )
  # check and prepare allCoefVcov
  allCoefVcov <- prepareVcov( allCoefVcov, nCoef, xPos, xMeanSd = NULL )
  # vector of shares of observations in each interval
  shareVec <- calcSharesInt( allXVal, xPos )
  # weights
  weights <- elaIntWeights( shareVec )
  # prepare calculation of semi-elasticity 
  if( model == "lpm" ) {
    xCoef <- rep( 0, length( xPos ) )
    for( i in 1:length( xPos ) ) {
      if( xPos[i] != 0 ) {
        xCoef[i] <- allCoef[ xPos[i] ]
      }
    }
  } else if( model == "probit" ) {
    # vector of probabilities of y=1 for each interval
    xBeta <- calcXBetaInt( allCoef, allXVal, xPos )
    checkXBeta( xBeta )
    xCoef <- pnorm( xBeta )
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  # semi-elasticities 'around' each inner boundary and their weights
  semElaBound <- rep( NA, nInt - 1 )
  for( m in 1:(nInt-1) ){
    semElaBound[m] <- 2 * ( xCoef[ m+1 ] - xCoef[ m ] ) * xBound[ m+1 ] /
      ( xBound[m+2] - xBound[m] )
  }
  # (average) semi-elasticity
  semEla <- sum( semElaBound * weights )
  # partial derivatives of semi-elasticities wrt coefficients
  if( model == "lpm" ) {
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
  } else if( model == "probit" ) {
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
    derivCoef <- rep( 0, nCoef )
    for( m in 1:( nInt - 1 ) ){
      derivCoef <- derivCoef + weights[m] * gradM[ , m ]
    }
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
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
