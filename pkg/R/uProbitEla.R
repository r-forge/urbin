uProbitEla <- function( allCoef, allXVal, allCoefVcov = NULL, xPos,
    seSimplify = !is.matrix( allCoefVcov ), xMeanSd = NULL ){

  nCoef <- length( allCoef )
  if( length( seSimplify ) != 1 || !is.logical( seSimplify ) ) {
    stop( "argument 'seSimplify' must be TRUE or FALSE" )
  }
  if( length( allXVal ) != nCoef ) {
    stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
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
  errMsgVcov <- paste( "argument 'allCoefVcov' must be a vector of length",
    nCoef, "or a symmetric matrix with dimension", nCoef )
  if( is.null( allCoefVcov ) ) {
    allCoefVcov <- matrix( NA, nrow = nCoef, ncol = nCoef )
    if( !is.null( xMeanSd ) ) {
      warning( "argument 'xMeanSd' is ignored,",
        " because argument 'allCoefVcov' has not been specified" )
    }
  } else if( is.matrix( allCoefVcov ) ) {
    if( nrow( allCoefVcov ) != nCoef || ncol( allCoefVcov ) != nCoef ) {
      stop( errMsgVcov )
    }
    if( !is.null( xMeanSd ) ) {
      warning( "argument 'xMeanSd' is ignored,",
        " the full variance-covariance matrix has been specified",
        " by argument 'allCoefVcov'" )
    }
  } else if( is.vector( allCoefVcov ) ) {
    if( length( allCoefVcov ) != nCoef ) {
      stop( errMsgVcov )
    } else {
      allCoefVcov <- diag( allCoefVcov^2 )
    }
    if( !is.null( xMeanSd ) ) {
      if( length( xPos ) != 2 ) {
        warning( "argument 'xMeanSd' is ignored,",
          " because the model does not include a quadratic term",
          " of the explanatory variable of interest" )
      } else {
        if( length( xMeanSd ) != 2 || !is.numeric( xMeanSd ) ) {
          stop( "argument 'xMeanSd' must be a vector of two numeric values")
        }
        set.seed( 123 )
        x <- rnorm( 1000, xMeanSd[ 1 ], xMeanSd[ 2 ] )
        X <- cbind( 1, x, x^2 )
        XXinv <- solve( t(X) %*% X )
        sigmaSq <- sqrt( ( allCoefVcov[ xPos[1], xPos[1] ] / XXinv[2,2] ) * 
            ( allCoefVcov[ xPos[2], xPos[2] ] / XXinv[3,3] ) )
        allCoefVcov[ xPos[1], xPos[2] ] <- allCoefVcov[ xPos[2], xPos[1] ] <- 
          sigmaSq * XXinv[2,3]
      }
    }
  } else {
    stop( errMsgVcov )
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
