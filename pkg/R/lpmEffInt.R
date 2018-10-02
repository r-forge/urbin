lpmEffInt <- function( allCoef, refBound, intBound, xPos,
  allCoefVcov = NULL ){

  # number of coefficients
  nCoef <- length( allCoef )
  # Check position vector
  checkXPos( xPos, minLength = 1, maxLength = 2, minVal = 1, maxVal = nCoef )
  # check and prepare allCoefVcov
  allCoefVcov <- prepareVcov( allCoefVcov, nCoef, xPos, xMeanSd = NULL )
  # check the boundaries of the intervals
  refBound <- elaIntBounds( refBound, 1, argName = "refBound" )
  intBound <- elaIntBounds( intBound, 1, argName = "intBound" )
  xCoef <- allCoef[xPos]
  if( length( xCoef ) == 1 ) {
    xCoef <- c( xCoef, 0 )
  }
  # difference between the xBars of the two intervals
  xDiff <- mean( intBound ) - mean( refBound )
  # difference between the xSquareBars of the two intervals
  xSquaredDiff <- 
    EXSquared( intBound[1], intBound[2] ) -
    EXSquared( refBound[1], refBound[2] )
  # effect E_{k,ml}
  eff <-  xCoef[1] * xDiff + xCoef[2] * xSquaredDiff
  # partial derivative of E_{k,ml} w.r.t. the beta_k and beta_{k+1}
  derivCoef <- rep( 0, nCoef ) 
  derivCoef[ xPos[1] ] <- xDiff
  if( length( xPos ) == 2 ) {
    derivCoef[ xPos[2] ] <- xSquaredDiff
  }
  # approximate standard error of the effect
  effSE <- drop( sqrt( t( derivCoef ) %*% allCoefVcov %*% derivCoef ) )
  # object to be returned
  result <- c( effect = eff, stdEr = effSE )
  return( result )
} 
