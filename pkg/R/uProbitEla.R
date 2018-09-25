uProbitEla <- function( allCoef, allXVal, allCoefVcov = NULL, xPos,
    seSimplify = !is.matrix( allCoefVcov ) ){

  nCoef <- length( allCoef )
  if( length( seSimplify ) != 1 || !is.logical( seSimplify ) ) {
    stop( "argument 'seSimplify' must be TRUE or FALSE" )
  }
  if( length( allCoef ) != length( allXVal ) ) {
    stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
  }
  errMsgVcov <- paste( "argument 'allCoefVcov' must be a vector of length",
    nCoef, "or a symmetric matrix with dimension", nCoef )
  if( is.null( allCoefVcov ) ) {
    allCoefVcov <- matrix( NA, nrow = nCoef, ncol = nCoef )
  } else if( is.matrix( allCoefVcov ) ) {
    if( nrow( allCoefVcov ) != nCoef || ncol( allCoefVcov ) != nCoef ) {
      stop( errMsgVcov )
    }
  } else if( is.vector( allCoefVcov ) ) {
    if( length( allCoefVcov ) != nCoef ) {
      stop( errMsgVcov )
    } else {
      allCoefVcov <- diag( allCoefVcov^2 )
    }
  } else {
    stop( errMsgVcov )
  }
  checkXPos( xPos, minLength = 1, maxLength = 2, minVal = 1, maxVal = nCoef )
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
  xVal <- allXVal[ xPos[ 1 ] ]
  xBeta <- sum( allCoef * allXVal )
  checkXBeta( xBeta )
  dfun <- dnorm( xBeta )
  semEla <- ( xCoef[ 1 ] + 2 * xCoef[ 2 ] * xVal ) * xVal * dfun
  derivCoef <- uProbitElaDeriv( allCoef, allXVal, xPos, seSimplify )
  semElaSE <- drop( sqrt( t( derivCoef ) %*% allCoefVcov %*% derivCoef ) )
  result <- c( semEla = semEla, stdEr = semElaSE )
  return( result )
} 
