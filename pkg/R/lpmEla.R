lpmEla <- function( xCoef, xVal, xCoefVcov = NULL ){
  # number of coefficients
  nCoef<- length( xCoef )
  if( ! nCoef %in% c( 1, 2 ) ) {
    stop( "argument 'xCoef' must be a scalar or vector with 2 elements" )
  }
  # check and prepare allCoefVcov
  xCoefVcov <- prepareVcov( xCoefVcov, nCoef, xPos = 1:nCoef, xMeanSd = NULL )
  # check argument xVal
  if( length( xVal ) != 1 || !is.numeric( xVal ) ) {
    stop( "argument 'xVal' must be a single numeric value" )
  }
  if( nCoef == 1 ) {
    xCoef <- c( xCoef, 0 )
    xCoefVcov <- matrix( c( xCoefVcov, 0, 0, 0 ), nrow = 2, ncol = 2 )
  }
  semEla <- ( xCoef[1] + 2 * xCoef[2] * xVal ) * xVal
  derivCoef <- c( xVal, ifelse( xCoef[2] == 0, 0, 2 * xVal^2 ) )
  # approximate standard error of the semi-elasticity
  semElaSE <- drop( sqrt( t( derivCoef ) %*% xCoefVcov %*% derivCoef ) )
  result <- c( semEla = semEla, stdEr = semElaSE )
  return( result )
}
