lpmEffInt <- function( allCoef, allXVal = NA, xPos, refBound, intBound,
  allCoefVcov = NULL, xMeanSd = NULL ){

  # number of coefficients
  nCoef <- length( allCoef )
  # Check position vector
  checkXPos( xPos, minLength = 1, maxLength = 2, minVal = 1, maxVal = nCoef )
  # check and prepare allCoefVcov
  allCoefVcov <- prepareVcov( allCoefVcov, nCoef, xPos, xMeanSd )
  # check the boundaries of the intervals
  refBound <- elaIntBounds( refBound, 1, argName = "refBound" )
  intBound <- elaIntBounds( intBound, 1, argName = "intBound" )
  # Check x values
  if( any( !is.na( allXVal ) ) ) {
    warning( "argument allXVal is ignored for lpm models",
      " (set this argument to 'NULL' or 'NA' to avoid this warning)" )
  }

  # calculate xBars
  intX <- mean( intBound )
  refX <- mean( refBound ) 
  if( length( xPos ) == 2 ) {
    intX <- c( intX, EXSquared( intBound[1], intBound[2] ) )
    refX <- c( refX, EXSquared( refBound[1], refBound[2] ) )
  }
  if( length( intX ) != length( xPos ) || 
      length( refX ) != length( xPos ) ) {
    stop( "internal error: 'intX' or 'refX' does not have the expected length" )
  }
  # define X' * beta 
  intXbeta <- sum( allCoef[ xPos ] * intX )
  refXbeta <- sum( allCoef[ xPos ] * refX )

  # effect E_{k,ml}
  eff <- intXbeta - refXbeta
  
  # partial derivative of the effect w.r.t. all estimated coefficients
  derivCoef <- rep( 0, nCoef ) 
  derivCoef[ xPos ] <- intX - refX
  # if argument allXVal has attribute 'derivOnly',
  # return partial derivatives only (for testing partial derivatives)
  if( "derivOnly" %in% names( attributes( allXVal ) ) ) {
    return( derivCoef )
  }
  
  # approximate standard error of the effect
  effSE <- drop( sqrt( t( derivCoef ) %*% allCoefVcov %*% derivCoef ) )
  
  # object to be returned
  result <- c( effect = eff, stdEr = effSE )
  return( result )
} 
