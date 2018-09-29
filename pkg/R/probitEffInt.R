probitEffInt <- function( allCoef, allXVal, xPos, refBound, intBound, 
  allCoefVcov = NULL ){
  # number of coefficients
  nCoef <- length( allCoef )
  # check arguments
  if( length( allXVal ) != nCoef ){
    stop( "argument 'allCoef' and 'allXVal' must have the same length" )
  }  
  checkXPos( xPos, minLength = 1, maxLength = 2, minVal = 1, maxVal = nCoef )
  # check and prepare allCoefVcov
  allCoefVcov <- prepareVcov( allCoefVcov, nCoef, xPos, xMeanSd = NULL )
  # check the boundaries of the intervals
  refBound <- elaIntBounds( refBound, 1, argName = "refBound" )
  intBound <- elaIntBounds( intBound, 1, argName = "intBound" )
  if( any( !is.na( allXVal[ xPos ] ) ) ) {
    allXVal[ xPos ] <- NA
    warning( "values of argument 'allXVal[ xPos ]' are ignored",
      " (set these values to 'NA' to avoid this warning)" )
  }
  
  model <- prepareEffInt( allCoef, allXVal, xPos, refBound, intBound )

  # effect E_{k,ml}
  eff <- pnorm( model$intXbeta ) - pnorm( model$refXbeta )
  
  # partial derivative of E_{k,ml} w.r.t. all estimated coefficients
  derivCoef <- probitEffIntDeriv( allCoef, allXVal, xPos, refBound, intBound )
  
  # approximate standard error of the effect
  effSE <- drop( sqrt( t( derivCoef ) %*% allCoefVcov %*% derivCoef ) )
  
  # object to be returned
  result <- c( effect = eff, stdEr = effSE )
  return( result )
}
