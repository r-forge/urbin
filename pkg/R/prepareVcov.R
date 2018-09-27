prepareVcov <- function( allCoefVcov, nCoef, xPos, xMeanSd ){

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
  return( allCoefVcov )
} 
