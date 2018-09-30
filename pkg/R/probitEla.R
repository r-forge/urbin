probitEla <- function( allCoef, allXVal, xPos, allCoefVcov = NULL,
    seSimplify = !is.matrix( allCoefVcov ), xMeanSd = NULL ){

  # check argument seSimplify
  if( length( seSimplify ) != 1 || !is.logical( seSimplify ) ) {
    stop( "argument 'seSimplify' must be TRUE or FALSE" )
  }
  # number of coefficients
  nCoef <- length( allCoef )
  # Check position vector
  checkXPos( xPos, minLength = 1, maxLength = 2, minVal = 1, maxVal = nCoef )
  # Check x values
  if( length( allXVal ) != nCoef ) {
    stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
  }
  # check and prepare allCoefVcov
  allCoefVcov <- prepareVcov( allCoefVcov, nCoef, xPos, xMeanSd )
  # Identify coefficients of interest (kth/tth covariate)
  if( length( xPos ) == 2 ){
    xCoef <- allCoef[ xPos ]
    if( !isTRUE( all.equal( allXVal[ xPos[2] ], allXVal[ xPos[1] ]^2 ) ) ) {
      stop( "the value of 'allXVal[ xPos[2] ]' must be equal",
        "to the squared value of 'allXVal[ xPos[1] ]' " )
    }
  } else if( length( xPos ) == 1 ) {
    xCoef <- c( allCoef[ xPos ], 0 )
  } else {
    stop( "argument 'xPos' must be a scalar or a vector with two elements" )
  }
  # prepare calculation of semi-elasticity 
  xVal <- allXVal[ xPos[ 1 ] ]
  xBeta <- sum( allCoef * allXVal )
  checkXBeta( xBeta )
  dfun <- dnorm( xBeta )
  semEla <- ( xCoef[ 1 ] + 2 * xCoef[ 2 ] * xVal ) * xVal * dfun
  # partial derivatives of semi-elasticities wrt coefficients
  if( seSimplify ) {
    derivCoef <- rep( 0, length( allCoef ) )
  } else {
    derivCoef <- ddnorm( xBeta ) * allXVal * 
      ( xCoef[ 1 ] + 2 * xCoef[ 2 ] * xVal ) * xVal
  }
  derivCoef[ xPos[1] ] <- derivCoef[ xPos[1] ] + dnorm( xBeta ) * xVal
  if( length( xPos ) == 2 ) {
    derivCoef[ xPos[2] ] <- derivCoef[ xPos[2] ] + dnorm( xBeta ) * 2 * xVal^2
  }
  # if argument allXVal has attribute 'derivOnly',
  # return partial derivatives only (for testing partial derivatives)
  if( "derivOnly" %in% names( attributes( allXVal ) ) ) {
    return( c( derivCoef ) )
  }
  # approximate standard error of the semi-elasticity
  semElaSE <- drop( sqrt( t( derivCoef ) %*% allCoefVcov %*% derivCoef ) )
  result <- c( semEla = semEla, stdEr = semElaSE )
  return( result )
} 
