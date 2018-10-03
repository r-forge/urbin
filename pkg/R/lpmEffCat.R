lpmEffCat <- function( xCoef, xShares, Group, 
  xCoefSE = rep( NA, length( xCoef ) ) ){
  if( sum( xShares ) > 1 ){
    stop( "the shares in argument 'xShares' sum up to a value larger than 1" )
  }
  if( length( xCoef ) != length( xShares ) ){
    stop( "arguments 'xCoef' and 'xShares' must have the same length" )
  }
  if( length( xCoef ) != length( Group ) ){
    stop( "arguments 'xCoef' and 'Group' must have the same length" )
  }
  if( ! all( Group %in% c( -1, 0, 1 ) ) ){
    stop( "all elements of argument 'Group' must be -1, 0, or 1" )
  }
  # D_mr
  DRef <- xShares * ( Group == -1 ) / sum( xShares[ Group == -1 ] )
  # D_ml
  DEffect <- xShares * ( Group == 1 ) / sum( xShares[ Group == 1 ] )
  # effect: sum of delta_m * ( D_ml - D_mr )
  effeG <- sum( xCoef * ( DEffect - DRef ) )
  # approximate standard error
  effeGSE <- sqrt( sum( ( xCoefSE * ( DEffect - DRef ) )^2 ) )
  result <- c( effect = effeG, stdEr = effeGSE )
  return( result )
}
