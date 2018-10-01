urbinEffInt <- function( allCoef, allXVal, xPos, refBound, intBound, model,
  allCoefVcov = NULL, xMeanSd = NULL ){
  
  if( model != "probit" ) {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  
  # number of coefficients
  nCoef <- length( allCoef )
  # check arguments
  if( length( allXVal ) != nCoef ){
    stop( "argument 'allCoef' and 'allXVal' must have the same length" )
  }  
  checkXPos( xPos, minLength = 1, maxLength = 2, minVal = 1, maxVal = nCoef )
  # check and prepare allCoefVcov
  allCoefVcov <- prepareVcov( allCoefVcov, nCoef, xPos, xMeanSd )
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
  
  # partial derivative of the effect w.r.t. all estimated coefficients
  derivCoef <- rep( NA, nCoef )
  derivCoef[ -xPos ] = ( dnorm( model$intXbeta ) - dnorm( model$refXbeta ) ) * 
    allXVal[ -xPos ] 
  derivCoef[ xPos ] = dnorm( model$intXbeta ) * model$intX - 
    dnorm( model$refXbeta ) * model$refX
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
