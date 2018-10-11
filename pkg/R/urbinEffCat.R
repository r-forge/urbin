urbinEffCat <- function( allCoef, allXVal, xPos, xGroups, model,
  allCoefVcov = NULL, iPos = ifelse( all( xPos != 1 ), 1, 0 ), 
  yCat = NULL ){
  
  # number of coefficients
  nCoef <- length( allCoef )
  # number of explanatory variables
  nXVal <- length( allXVal )
  # Check position vector
  checkXPos( xPos, minLength = 1, maxLength = nCoef, minVal = 1, 
    maxVal = ifelse( model == "MNL", nXVal, nCoef ) )
  # check position of the intercept
  checkIPos( iPos, xPos, allXVal, model ) 
  # number of categories
  nCat <- length( xPos ) + 1
  # check allXVal and allCoef
  if( model %in% c( "lpm", "probit", "oprobit", "logit" ) ){
    xCoef <- c( allCoef[ xPos ], 0 )
    if( nXVal != nCoef ){
      stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
    }  
  } else if( model == "MNL" ){
    # number of alternative categories of the dependent variable
    nYCat <- round( nCoef / nXVal )
    if( nCoef != nXVal * nYCat ) {
      stop( "length of argument 'allCoef' must be a multiple",
        " of the length of argument 'allXVal'" )
    } 
    # create matrix of coefficients
    mCoef <- matrix( allCoef, nrow = nXVal, ncol = nYCat )
    xCoef <- rbind( mCoef[ xPos, ], 0 )
  } else if( model == "CondL" ){
    # number of categories of the dependent variable
    nYCat <- round( nXVal / nCoef )
    if( nXVal != nCoef * nYCat ) {
      stop( "length of argument 'allXVal' must be a multiple",
        " of the length of argument 'allCoef'" )
    } 
    # create matrix of explanatory variables
    mXVal <- matrix( allXVal, nrow = nCoef )
    xCoef <- c( allCoef[ xPos ], 0 )
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  # check argument yCat
  if( model %in% c( "MNL", "CondL" ) ) {
    checkYCat( yCat, nYCat ) 
  } else if( model != "NestedL" && !is.null( yCat ) ) {
    warning( "argument 'yCat' is ignored" )
  }
  # shares in each category
  if( model %in% c( "lpm", "probit", "oprobit", "logit", "MNL" ) ){
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
  } else if( model == "CondL" ){
    xShares <- mXVal[ xPos, ]
    if( any( xShares < 0 ) ){
      stop( "the share of the observations in at least one category",
        " is negative" )
    }
    for( p in 1:nYCat ){
      if( sum( xShares[ , p ] ) > 1 ){
        stop( "the shares in argument 'xShares' sum up to a value larger than 1" ) 
      }
    }
    xShares <- rbind( xShares, 1 - colSums( xShares ) )
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }  
  # check argument xGroups
  if( length( xGroups ) != ( length( xPos ) + 1 ) ){
    stop( "the vector specified by argument 'xGroups' must have",
      " one more element that the vector specified by argument 'xPos'" )
  }
  if( ! all( xGroups %in% c( -1, 0, 1 ) ) ){
    stop( "all elements of argument 'xGroups' must be -1, 0, or 1" )
  }
  # check and prepare allCoefVcov
  allCoefVcov <- prepareVcov( allCoefVcov, nCoef, xPos = NA, 
    pCall = match.call() )
  if( model %in% c( "lpm", "probit", "oprobit", "logit" ) ){
    # D_mr
    DRef <- ifelse( xGroups == -1, xShares, 0 ) / 
      sum( xShares[ xGroups == -1 ] )
    # D_ml
    DEffect <- ifelse( xGroups == 1, xShares, 0 ) / 
      sum( xShares[ xGroups == 1 ] )
    if( model == "lpm" ) {
      # effect: sum of delta_m * ( D_ml - D_mr )
      effeG <- sum( xCoef * ( DEffect - DRef ) )
    } else {
      # linear predictors
      XBetaRef <- sum( allCoef[ -xPos ] * allXVal[ -xPos ]) + 
        sum( DRef * xCoef )
      XBetaEffect <- sum( allCoef[ -xPos ] * allXVal[ -xPos ]) + 
        sum( DEffect * xCoef )
      checkXBeta( c( XBetaRef, XBetaEffect ) )
      if( model %in% c( "probit", "oprobit" ) ) {
        # effect
        effeG <- pnorm( XBetaEffect ) - pnorm( XBetaRef )
      } else if( model == "logit" ) {
        # effect
        effeG <- exp( XBetaEffect )/( 1 + exp( XBetaEffect ) ) - 
          exp( XBetaRef )/( 1 + exp( XBetaRef ) )
      }
    }
  } else if( model == "MNL" ){
    # D_mr  
    DRef <- ifelse( xGroups == -1, xShares, 0 ) / 
      sum( xShares[ xGroups == -1 ] )
    XBetaRef <- allXVal[ -xPos ] %*% mCoef[ -xPos, , drop = FALSE ] + 
      DRef %*% xCoef
    # D_ml
    DEffect <- ifelse( xGroups == 1, xShares, 0 ) / 
      sum( xShares[ xGroups == 1 ] )
    XBetaEffect <- allXVal[ -xPos ] %*% mCoef[ -xPos, , drop = FALSE ] + 
      DEffect %*% xCoef
    checkXBeta( c( XBetaRef, XBetaEffect ) )
    # effect
    effeG <- exp( XBetaEffect[ yCat ] )/( 1 + sum( exp( XBetaEffect ) ) ) -
      exp( XBetaRef[ yCat ] )/( 1 + sum( exp( XBetaRef ) ) )
  } else if( model == "CondL" ){
    # D_mr
    DRef <- matrix( NA, nrow = nCat, ncol = nYCat )
    for( p in 1:nYCat ) {
      DRef[ , p ] <- ifelse( xGroups == -1, xShares[ , p ], 0 ) /
        sum( xShares[ xGroups == -1, p ] )
    }
    XBetaRef <- allCoef[ -xPos ] %*% mXVal[ -xPos, , drop = FALSE ] + 
      xCoef %*% DRef
    # D_ml
    DEffect <- matrix( NA, nrow = nCat, ncol = nYCat )
    for( p in 1:nYCat ) {
      DEffect[ , p ] <- ifelse( xGroups == 1, xShares[ , p ], 0 ) /
        sum( xShares[ xGroups == 1, p ] )
    }
    XBetaEffect <- allCoef[ -xPos ] %*% mXVal[ -xPos, , drop = FALSE ] +
      xCoef %*% DEffect
    checkXBeta( c( XBetaRef, XBetaEffect ) )
    # effect
    effeG <- exp( XBetaEffect[ yCat ] )/( sum( exp( XBetaEffect ) ) ) -
      exp( XBetaRef[ yCat ] )/( sum( exp( XBetaRef ) ) )
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  # partial derivative of the effect w.r.t. all estimated coefficients
  if( model == "lpm" ) {
    derivCoef <- rep( 0, nCoef )
    derivCoef[ xPos ] <- DEffect[ -nCat ] - DRef[ -nCat ]
  } else if( model %in% c( "probit", "oprobit" ) ) {
    derivCoef <- rep( NA, nCoef )
    derivCoef[ -xPos ] = ( dnorm( XBetaEffect ) - dnorm( XBetaRef ) ) * 
      allXVal[ -xPos ] 
    derivCoef[ xPos ] = dnorm( XBetaEffect ) * DEffect[ -nCat ] - 
      dnorm( XBetaRef ) * DRef[ -nCat ]
  } else if( model == "logit" ){
    derivCoef <- rep( NA, nCoef )
    derivCoef[ -xPos ] <- ( exp( XBetaEffect )/( 1 + exp( XBetaEffect ))^2 - 
        exp( XBetaRef )/( 1 + exp( XBetaRef ))^2 ) * 
      allXVal[ -xPos ] 
    derivCoef[ xPos ] <- exp( XBetaEffect )/( 1 + exp( XBetaEffect))^2 * DEffect[ -nCat ] - 
      exp( XBetaRef )/( 1 + exp( XBetaRef ))^2 * DRef[ -nCat ]
  } else if( model == "MNL" ){
    derivCoef <- matrix( NA, nrow = nXVal, ncol = nYCat )
    for( p in 1:nYCat ){
      if( p == yCat ){
        derivCoef[ -xPos, p ] <- 
          ( exp( XBetaEffect[ p ] ) * 
              ( 1 + sum( exp( XBetaEffect[ -yCat ] ) ) )/
              ( 1 + sum( exp( XBetaEffect ) ) )^2 -
              exp( XBetaRef[ p ] ) * 
              ( 1 + sum( exp( XBetaRef[ -yCat ] ) ) )/
              ( 1 + sum( exp( XBetaRef ) ) )^2 ) * allXVal[ -xPos ]
        derivCoef[ xPos, p ] <- 
          ( exp( XBetaEffect[ p ] ) * 
              ( 1 + sum( exp( XBetaEffect[ -yCat ] ) ) )/
              ( 1 + sum( exp( XBetaEffect ) ) )^2 ) * DEffect[-nCat] -
          ( exp( XBetaRef[ p ] ) * 
              ( 1 + sum( exp( XBetaRef[ -yCat ] ) ) )/
              ( 1 + sum( exp( XBetaRef ) ) )^2 ) * DRef[-nCat]
      } else{  
        derivCoef[ -xPos, p ] <- 
          ( ( exp( XBetaRef[ yCat ] ) * exp( XBetaRef[ p ] ) )/
              ( 1 + sum( exp( XBetaRef ) ) )^2 -
              ( exp( XBetaEffect[ yCat ] ) * exp( XBetaEffect[ p ] ) )/
              ( 1 + sum( exp( XBetaEffect ) ) )^2 ) * allXVal[ -xPos ]
        derivCoef[ xPos, p ] <- 
          ( ( exp( XBetaRef[ yCat ] ) * exp( XBetaRef[ p ] ) )/
              ( 1 + sum( exp( XBetaRef ) ) )^2 ) * DRef[-nCat] -
          ( ( exp( XBetaEffect[ yCat ] ) * exp( XBetaEffect[ p ] ) )/
              ( 1 + sum( exp( XBetaEffect ) ) )^2 ) * DEffect[-nCat]
      }     
    }
    derivCoef <- c( derivCoef )
  } else if( model == "CondL" ){
    derivCoef <- rep( NA, nCoef )
    derivCoef[ -xPos ] = ( exp( XBetaEffect[ yCat ] ) * mXVal[ -xPos, yCat ] * 
        sum( exp( XBetaEffect ) ) -
        exp( XBetaEffect[ yCat ] ) * sum( exp( XBetaEffect ) * 
            mXVal[ -xPos, ] ) )/
      ( sum( exp( XBetaEffect ) ) )^2 - 
      ( exp( XBetaRef[ yCat ] ) * mXVal[ -xPos, yCat ] * 
          sum( exp( XBetaRef ) ) -
          exp( XBetaRef[ yCat ] ) * sum( exp( XBetaRef ) * 
              mXVal[ -xPos, ] ) )/
      ( sum( exp( XBetaRef ) ) )^2
    derivCoef[ xPos ] =  ( exp( XBetaEffect[ yCat ] ) * DEffect[ yCat ] * 
        sum( exp( XBetaEffect ) ) -
        exp( XBetaEffect[ yCat ] ) * 
        sum( exp( XBetaEffect ) * DEffect[ yCat ] ) )/
      ( sum( exp( XBetaEffect ) ) )^2 - 
      ( exp( XBetaRef[ yCat ] ) * DRef[ yCat ] * 
          sum( exp( XBetaRef ) ) -
          exp( XBetaRef[ yCat ] ) * 
          sum( exp( XBetaRef ) * DRef[ yCat ] ) )/
      ( sum( exp( XBetaRef ) ) )^2
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  # if argument allXVal has attribute 'derivOnly',
  # return partial derivatives only (for testing partial derivatives)
  if( "derivOnly" %in% names( attributes( allXVal ) ) ) {
    return( derivCoef )
  }
  # approximate standard error of the effect
  effeGSE <- drop( sqrt( t( derivCoef ) %*% allCoefVcov %*% derivCoef ) )
  # object to be returned
  result <- list()
  result$call <- match.call()
  result$allCoefVcov <- allCoefVcov
  result$derivCoef <- derivCoef
  result$effect <- effeG
  result$stdEr <- effeGSE
  class( result ) <- "urbin"
  return( result )
}
