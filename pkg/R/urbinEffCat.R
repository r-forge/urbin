urbinEffCat <- function( allCoef, allXVal, xPos, xGroups, model,
    allCoefVcov = NULL ){
  
  if( model != "probit" ) {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  
  # number of coefficients
  nCoef <- length( allCoef )
  # Check position vector
  checkXPos( xPos, minLength = 1, maxLength = nCoef, minVal = 1, 
    maxVal = nCoef )
  # number of categories
  nCat <- length( xPos ) + 1
  # shares in each category
  xShares <- allXVal[ xPos ]
  if( any( xShares < 0 ) ){
    stop( "the share of the observations in at least one category",
      " is negative" )
  }
  if( sum( xShares ) > 1 ){
    stop( "the shares of the observations in the individual categories",
      " (without the reference category) sum up to a value larger than 1" )
  }
  xShares <- c( xShares, 1 - sum( xShares ) )
  if( length( xShares ) != nCat ) {
    stop( "internal error: length of shares not equal to number of categories" )
  }
  if( xShares[ nCat ] < 0.05  ) {
    warning( "there are only ", 100 * xShares[ length( xShares ) ],
      "% of the observations in the reference category --",
      " please check whether this is indeed the case" )
  }
  if( xShares[ nCat ] > 0.95  ) {
    warning( "there are ", 100 * xShares[ length( xShares ) ],
      "% of the observations in the reference category --",
      " please check whether this is indeed the case" )
  }
  # coefficients of the dummy variables for the categories
  xCoef <- c( allCoef[ xPos ], 0 )
  # check argument xGroups
  if( length( xGroups ) != ( length( xPos ) + 1 ) ){
    stop( "the vector specified by argument 'xGroups' must have",
      " one more element that the vector specified by argument 'xPos'" )
  }
  if( ! all( xGroups %in% c( -1, 0, 1 ) ) ){
    stop( "all elements of argument 'xGroups' must be -1, 0, or 1" )
  }
  # check and prepare allCoefVcov
  allCoefVcov <- prepareVcov( allCoefVcov, nCoef, xPos = NA, xMeanSd = NULL )
  # D_mr
  DRef <- ifelse( xGroups == -1, xShares, 0 ) / 
    sum( xShares[ xGroups == -1 ] )
  # D_ml
  DEffect <- ifelse( xGroups == 1, xShares, 0 ) / 
    sum( xShares[ xGroups == 1 ] )
  # linear predictors
  XBetaRef <- sum( allCoef[ -xPos ] * allXVal[ -xPos ]) + 
    sum( DRef * xCoef )
  XBetaEffect <- sum( allCoef[ -xPos ] * allXVal[ -xPos ]) + 
    sum( DEffect * xCoef )
  # effect
  effeG <- pnorm( XBetaEffect ) - pnorm( XBetaRef )
  # partial derivative of the effect w.r.t. all estimated coefficients
  derivCoef <- rep( NA, nCoef )
  derivCoef[ -xPos ] = ( dnorm( XBetaEffect ) - dnorm( XBetaRef ) ) * 
    allXVal[ -xPos ] 
  derivCoef[ xPos ] = dnorm( XBetaEffect ) * DEffect[ -nCat ] - 
    dnorm( XBetaRef ) * DRef[ -nCat ]
  # if argument allXVal has attribute 'derivOnly',
  # return partial derivatives only (for testing partial derivatives)
  if( "derivOnly" %in% names( attributes( allXVal ) ) ) {
    return( derivCoef )
  }
  # approximate standard error of the effect
  effeGSE <- drop( sqrt( t( derivCoef ) %*% allCoefVcov %*% derivCoef ) )
  # object to be returned
  result <- c( effect = effeG, stdEr = effeGSE )
  return( result )
}
