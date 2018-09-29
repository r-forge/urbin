probitEffGroup <- function( allCoef, allXVal, xPos, xGroups, 
    allCoefVcov = NULL ){
  nCoef <- length( allCoef )
  # check argument xPos
  checkXPos( xPos, minLength = 1, maxLength = nCoef, minVal = 1, 
    maxVal = nCoef )
  # shares in each category
  xShares <- allXVal[ xPos ]
  # coefficients of the dummy variables for the categories
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
  DRef <- ifelse( xGroups == -1, xShares, 0 ) / 
    sum( xShares[ xGroups == -1 ] )
  XBetaRef <- sum( allCoef[ -xPos ] * allXVal[ -xPos ]) + 
    sum( DRef * xCoef )
  # D_ml
  DEffect <- ifelse( xGroups == 1, xShares, 0 ) / 
    sum( xShares[ xGroups == 1 ] )
  XBetaEffect <- sum( allCoef[ -xPos ] * allXVal[ -xPos ]) + 
    sum( DEffect * xCoef )
  # effect
  effeG <- pnorm( XBetaEffect ) - pnorm( XBetaRef )
  # partial derivative of E_{k,ml} w.r.t. all estimated coefficients
  derivCoef <- rep( NA, nCoef )
  derivCoef[ -xPos ] = ( dnorm( XBetaEffect ) - dnorm( XBetaRef ) ) * 
    allXVal[ -xPos ] 
  derivCoef[ xPos ] = dnorm( XBetaEffect ) * DEffect - dnorm( XBetaRef ) * DRef
  # if argument allXVal has attribute 'derivOnly',
  # return partial derivatives only (for testing partial derivatives)
  if( "derivOnly" %in% names( attributes( allXVal ) ) ) {
    return( derivCoef )
  }
  # approximate standard error of the effect
  effeGSE <- drop( sqrt( t( derivCoef ) %*% allCoefVcov %*% derivCoef ) )
  result <- c( effect = effeG, stdEr = effeGSE )
  return( result )
}
