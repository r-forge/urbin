checkIPos <- function( iPos, xPos, allXVal ) {
  if( length( iPos ) != 1 ) {
    stop( "argument 'iPos' must be a single integer value" )
  }
  if( iPos != round( iPos ) ) {
    stop( "argument 'iPos' must be an integer" )
  }
  if( iPos < 0 ) {
    stop( "argument 'iPos' must be a non-negative value" )
  }
  if( iPos > length( allXVal ) ) {
    stop( "argument 'iPos' must be smaller than or equal to ",
      length( allXVal ) )
  }
  if( iPos > 0 & !is.null( allXVal ) & !is.list( allXVal ) ) {
    if( ! allXVal[ iPos ] %in% c( 1, -1 ) ) {
      stop( "the value of argument 'allXVal' indicated by argument 'iPos',",
        " i.e., allXVal[iPos], must be one",
        " (or minus one for ordered probit models)" )
    }
  }
}
