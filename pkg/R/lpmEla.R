lpmEla <- function( allCoef, allXVal, allCoefVcov = NULL, xPos = NULL, xMeanSd = NULL ){
  
  # number of coefficients
  nCoef <- length( allCoef )
  # Obtain position vector if necessary
  if( is.null( xPos ) ) {
    xPos <- 1:nCoef
  }
  # Check position vector
  checkXPos( xPos, minLength = 1, maxLength = 2, minVal = 1, maxVal = nCoef )
  # Check x values
  if( nCoef == 2 && length( allXVal ) == 1 ) {
    temp <- c( allXVal, allXVal^2 )
    if( "derivOnly" %in% names( attributes( allXVal ) ) ) {
      attr( temp, "derivOnly" ) <- 1
    }
    allXVal <- temp
  }
  if( nCoef != length( allXVal ) ) {
    stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
  }
  # check and prepare allCoefVcov
  allCoefVcov <- prepareVcov( allCoefVcov, nCoef, xPos, xMeanSd )
  # Identify coefficients of interest (kth/tth covariate)
  if( length( xPos ) == 2 ){
    xCoef <- allCoef[ xPos ]
    if( !isTRUE( all.equal( allXVal[xPos[2]], allXVal[xPos[1]]^2 ) ) ) {
      stop( "the value of 'allXVal[ xPos[2] ]' must be equal",
        "to the squared value of 'allXVal[ xPos[1] ]' " )
    }
  } else if( length( xPos ) == 1 ) {
    xCoef <- c( allCoef[ xPos ], 0 )
  }
  # prepare calculation of semi-elasticity 
  xVal <- allXVal[ xPos[ 1 ] ]
  semEla <- ( xCoef[1] + 2 * xCoef[2] * xVal ) * xVal
  # partial derivatives of semi-elasticities wrt coefficients
  derivCoef <- rep( 0, nCoef ) 
  derivCoef[ xPos[1] ] <- xVal
  if( length( xPos ) == 2 ) {
    derivCoef[ xPos[2] ] <- 2 * xVal^2
  }
  # if argument allXVal has attribute 'derivOnly',
  # return partial derivatives only (for testing partial derivatives)
  if( "derivOnly" %in% names( attributes( allXVal ) ) ) {
    return( derivCoef )
  }
  # approximate standard error of the semi-elasticity
  semElaSE <- drop( sqrt( t( derivCoef ) %*% allCoefVcov %*% derivCoef ) )
  result <- c( semEla = semEla, stdEr = semElaSE )
  return( result )
}
