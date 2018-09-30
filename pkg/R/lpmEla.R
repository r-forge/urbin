lpmEla <- function( xCoef, xVal, xCoefSE = rep( NA, length( xCoef ) ) ){
  if( ! length( xCoef ) %in% c( 1, 2 ) ) {
    stop( "argument 'xCoef' must be a scalar or vector with 2 elements" )
  }
  if( length( xCoef ) != length( xCoefSE ) ) {
    stop( "arguments 'xCoef' and 'xCoefSE' must have the same length" )
  }
  if( length( xVal ) != 1 || !is.numeric( xVal ) ) {
    stop( "argument 'xVal' must be a single numeric value" )
  }
  if( length( xCoef ) == 1 ) {
    xCoef <- c( xCoef, 0 )
    xCoefSE <- c( xCoefSE, 0 )
  }
  semEla <- ( xCoef[1] + 2 * xCoef[2] * xVal ) * xVal
  derivCoef <- c( xVal, ifelse( xCoef[2] == 0, 0, 2 * xVal^2 ) )
  vcovCoef <- diag( xCoefSE^2 )
  semElaSE <- drop( sqrt( t( derivCoef ) %*% vcovCoef %*% derivCoef ) )
  result <- c( semEla = semEla, stdEr = semElaSE )
  return( result )
}
