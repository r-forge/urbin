probitEffIntDeriv <- function( allCoef, allXVal, xPos, refBound, intBound ){
  # number of coefficients
  nCoef <- length( allCoef )

  # calculate linear predictors, etc
  model <- prepareEffInt( allCoef, allXVal, xPos, refBound, intBound )

  # partial derivative of E_{k,ml} w.r.t. all estimated coefficients
  derivCoef <- rep( NA, nCoef )
  derivCoef[ -xPos ] = ( dnorm( model$intXbeta ) - dnorm( model$refXbeta ) ) * 
    allXVal[ -xPos ] 
  derivCoef[ xPos ] = dnorm( model$intXbeta ) * model$intX - 
    dnorm( model$refXbeta ) * model$refX

  return( derivCoef )
}
