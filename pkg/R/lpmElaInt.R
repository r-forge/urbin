lpmElaInt <- function( allCoef, allXVal, xBound, xPos = NULL,
  allCoefVcov = NULL ){
  
  if( is.null( xPos ) ) {
    xCoef <- allCoef
    xShares <- allXVal
    xCoefVcov <- allCoefVcov
  } else {
    # Check position vector
    checkXPos( xPos, minLength = 2, maxLength = length( allCoef ) + 1, 
      minVal = 0, maxVal = length( allCoef ), requiredVal = 0 )
    xCoef <- rep( 0, length( xPos ) )
    for( i in 1:length( xPos ) ) {
      if( xPos[i] != 0 ) {
        xCoef[i] <- allCoef[ xPos[i] ]
      }
    }
    xShares <- calcSharesInt( allXVal, xPos )
    if( "derivOnly" %in% names( attributes( allXVal ) ) ) {
      attr( xShares, "derivOnly" ) <- 1 
    }
    if( is.null( allCoefVcov ) ) {
      xCoefVcov <- NULL
    } else {
      if( is.matrix( allCoefVcov ) ) {
        temp <- matrix( 0, length( xPos ), ncol( allCoefVcov ) )
        for( i in 1:length( xPos ) ) {
          if( xPos[i] != 0 ) {
            temp[i,] <- allCoefVcov[ xPos[i], ]
          }
        }
        xCoefVcov <- matrix( 0, length( xPos ), length( xPos ) )
        for( i in 1:length( xPos ) ) {
          if( xPos[i] != 0 ) {
            xCoefVcov[,i] <- temp[ , xPos[i] ]
          }
        }
      } else {
        xCoefVcov <- rep( 0, length( xPos ) )
        for( i in 1:length( xPos ) ) {
          if( xPos[i] != 0 ) {
            xCoefVcov[i] <- allCoefVcov[ xPos[i] ]
          }
        }
      }
    }
  }
  # number of intervals
  nInt <- length( xCoef )
  if( nInt < 2 || !is.vector( xCoef ) ) {
    stop( "argument 'xCoef' must be a vector with at least two elements" )
  }
  if( length( xShares ) != nInt ) {
    stop( "arguments 'xCoef' and 'xShares' must be vectors of the same length" )
  }
  if( any( xShares < 0 ) ) {
    stop( "all shares in argument 'xShares' must be non-negative" )
  }
  if( abs( sum( xShares ) - 1 ) > 0.015 ) {
    stop( "the shares in argument 'xShares' must sum to one" )
  }
  # check 'xBound' and replace infinite values
  xBound <- elaIntBounds( xBound, nInt )
  # check and prepare allCoefVcov
  xCoefVcov <- prepareVcov( xCoefVcov, length( xCoef ), NA, xMeanSd = NULL )
  # weights
  weights <- elaIntWeights( xShares )
  # semi-elasticities 'around' each inner boundary and their weights
  semElas <- rep( NA, nInt - 1 )
  for( m in 1:(nInt-1) ){
    semElas[m] <- 2 * ( xCoef[ m+1 ] - xCoef[ m ] ) * xBound[ m+1 ] /
      ( xBound[m+2] - xBound[m] )
  }
  # (average) semi-elasticity
  semElaAvg <- sum( semElas * weights )
  # derivatives of the (average) semi-elasticity wrt the coefficients
  derivCoef <- rep( NA, nInt )
  derivCoef[1] <- 
    -2 * weights[1] * xBound[2] / ( xBound[3] - xBound[1] )
  derivCoef[nInt] <- 
    2 * weights[nInt-1] * xBound[nInt] / ( xBound[nInt+1] - xBound[nInt-1] ) 
  if( nInt > 2 ) {
    for( n in 2:( nInt-1 ) ) {
      derivCoef[n] <- 
        2 * weights[n-1] * xBound[n] / ( xBound[n+1] - xBound[n-1] ) -
        2 * weights[n]   * xBound[n+1] / ( xBound[n+2] - xBound[n] )
    }
  }
  # if argument allXVal has attribute 'derivOnly',
  # return partial derivatives only (for testing partial derivatives)
  if( "derivOnly" %in% names( attributes( xShares ) ) ) {
    return( derivCoef )
  }
  # standard error of the (average) semi-elasticity
  semElaSE <- drop( sqrt( t( derivCoef ) %*% xCoefVcov %*% derivCoef ) )
  # prepare object that will be returned
  result <- c( semEla = semElaAvg, stdEr = semElaSE )
  return( result )
}
