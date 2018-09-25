uProbitElaDeriv <- function( allCoef, allXVal, xPos, simplified = FALSE ){
  if( length( xPos ) == 2 ){
    xCoef <- allCoef[ xPos ]
  } else if( length( xPos ) == 1 ) {
    xCoef <- c( allCoef[ xPos ], 0 )
  } else {
    stop( "argument 'xPos' must be a scalar or a vector with two elements" )
  }
  xVal <- allXVal[ xPos[ 1 ] ]
  xBeta <- sum( allCoef * allXVal )
  if( simplified ) {
    d <- rep( 0, length( allCoef ) )
  } else {
    d <- ddnorm( xBeta ) * allXVal * ( xCoef[ 1 ] + 2 * xCoef[ 2 ] * xVal ) * xVal
  }
  d[ xPos[1] ] <- d[ xPos[1] ] + dnorm( xBeta ) * xVal
  if( length( xPos ) == 2 ) {
    d[ xPos[2] ] <- d[ xPos[2] ] + dnorm( xBeta ) * 2 * xVal^2
  }
  return( d )
}
