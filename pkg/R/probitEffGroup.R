probitEffGroup <- function( allCoef, allXVal, xPos, xGroups, 
    allCoefVcov = NULL ){
  nCoef <- length( allCoef )
  xShares <- allXVal[ xPos ]
  xCoef <- allCoef[ xPos ]
  if( sum( xShares ) > 1 ){
    stop( "the shares in argument 'xShares' sum up to a value larger than 1" )
  }
  if( length( xCoef ) != length( xShares ) ){
    stop( "arguments 'xCoef' and 'xShares' must have the same length" )
  }
  if( length( xCoef ) != length( xGroups ) ){
    stop( "arguments 'xCoef' and 'xGroups' must have the same length" )
  }
  if( ! all( xGroups %in% c( -1, 0, 1 ) ) ){
    stop( "all elements of argument 'xGroups' must be -1, 0, or 1" )
  }
  # check and prepare allCoefVcov
  allCoefVcov <- prepareVcov( allCoefVcov, nCoef, xPos = NA, xMeanSd = NULL )
  # D_mr
  DRef <- sum( xCoef[ xGroups == -1 ] * xShares[ xGroups == -1 ]) / 
    sum( xShares[ xGroups == -1 ] )
  XBetaRef <- sum( allCoef[ -xPos ] * allXVal[ -xPos ]) + DRef
  # D_ml
  DEffect <- sum( xCoef[ xGroups == 1 ] * xShares[ xGroups == 1 ]) / 
    sum( xShares[ xGroups == 1 ] )
  XBetaEffect <- sum( allCoef[ -xPos ] * allXVal[ -xPos ]) + DEffect
  # effect
  effeG <- pnorm( XBetaEffect ) - pnorm( XBetaRef )
  # partial derivative of E_{k,ml} w.r.t. all estimated coefficients
  derivCoef <- rep( NA, nCoef )
  derivCoef[ -xPos ] = ( dnorm( XBetaEffect ) - dnorm( XBetaRef ) ) * 
    allXVal[ -xPos ] 
  derivCoef[ xPos ] = dnorm( XBetaEffect ) * DEffect - dnorm( XBetaRef ) * DRef
  # approximate standard error of the effect
  effeGSE <- drop( sqrt( t( derivCoef ) %*% allCoefVcov %*% derivCoef ) )
  result <- c( effect = effeG, stdEr = effeGSE )
  return( result )
}
