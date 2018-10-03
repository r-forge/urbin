lpmEffCat <- function( allCoef, allXVal, xPos, Group, 
  allCoefVcov = NULL ){
  
  # number of coefficients
  nCoef <- length( allCoef )
  # Check position vector
  checkXPos( xPos, minLength = 1, maxLength = nCoef, minVal = 1, 
    maxVal = nCoef )

  xCoef <- allCoef[ xPos ]
  xShares <- allXVal[ xPos ]

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
  # check and prepare allCoefVcov
  allCoefVcov <- prepareVcov( allCoefVcov, nCoef, xPos = NA, xMeanSd = NULL )
  # D_mr
  DRef <- xShares * ( Group == -1 ) / sum( xShares[ Group == -1 ] )
  # D_ml
  DEffect <- xShares * ( Group == 1 ) / sum( xShares[ Group == 1 ] )
  # effect: sum of delta_m * ( D_ml - D_mr )
  effeG <- sum( xCoef * ( DEffect - DRef ) )
  # partial derivative of the effect w.r.t. all estimated coefficients
  derivCoef <- rep( 0, nCoef )
  derivCoef[ xPos ] <- DEffect - DRef
  # approximate standard error of the effect
  effeGSE <- drop( sqrt( t( derivCoef ) %*% allCoefVcov %*% derivCoef ) )
  # object to be returned
  result <- c( effect = effeG, stdEr = effeGSE )
  return( result )
}
