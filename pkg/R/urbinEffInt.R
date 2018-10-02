urbinEffInt <- function( allCoef, allXVal, xPos, refBound, intBound, model,
  allCoefVcov = NULL, xMeanSd = NULL ){
  
  if( model != "probit" ) {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  
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
  if( length( allXVal ) != nCoef ){
    stop( "argument 'allCoef' and 'allXVal' must have the same length" )
  }  
  if( any( !is.na( allXVal[ xPos ] ) ) ) {
    allXVal[ xPos ] <- NA
    warning( "values of argument 'allXVal[ xPos ]' are ignored",
      " (set these values to 'NA' to avoid this warning)" )
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
  intXbeta <- sum( allCoef * replace( allXVal, xPos, intX ) )
  refXbeta <- sum( allCoef * replace( allXVal, xPos, refX ) )
  checkXBeta( c( intXbeta, refXbeta ) )
  
  # effect E_{k,ml}
  eff <- pnorm( intXbeta ) - pnorm( refXbeta )
  
  # partial derivative of the effect w.r.t. all estimated coefficients
  derivCoef <- rep( NA, nCoef )
  derivCoef[ -xPos ] = ( dnorm( intXbeta ) - dnorm( refXbeta ) ) * 
    allXVal[ -xPos ] 
  derivCoef[ xPos ] = dnorm( intXbeta ) * intX - 
    dnorm( refXbeta ) * refX
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
