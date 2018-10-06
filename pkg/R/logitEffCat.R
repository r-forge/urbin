logitEffCat <- function( allCoef, allXVal, xPos, xGroups, model,
  yCat = NA, allCoefVcov = NULL ){
  
  # number of coefficients
  nCoef <- length( allCoef )
  # number of explanatory variables
  nXVal <- length( allXVal )
  # Check position vector
  checkXPos( xPos, minLength = 1, maxLength = nCoef, minVal = 1, 
    maxVal = nCoef )
  # check allXVal and allCoef
  if( model == "logit" ){
    xCoef <- allCoef[ xPos ]
    if( nXVal != nCoef ){
      stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
    }  
  } else if( model == "MNL" ){
    # number of ???
    pCoef <- round( nCoef / nXVal )
    if( nCoef != nXVal * pCoef ) {
      stop( "length of argument 'allCoef' must be a multiple",
        " of the length of argument 'allXVal'" )
    } 
    # create matrix of coefficients
    mCoef <- matrix( allCoef, nrow = nXVal, ncol = pCoef )
    xCoef <- mCoef[ xPos, ]
  } else if( model == "CondL" ){
    # number of ???
    pCoef <- round( nXVal / nCoef )
    if( nXVal != nCoef * pCoef ) {
      stop( "length of argument 'allXVal' must be a multiple",
        " of the length of argument 'allCoef'" )
    } 
    # create matrix of explanatory variables
    mXVal <- matrix( allXVal, nrow = nCoef )
    xCoef <- allCoef[ xPos ]
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  # shares in each category
  if( model == "logit" || model == "MNL" ){
    xShares <- allXVal[ xPos ]
    if( sum( xShares ) > 1 ){
      stop( "the shares in argument 'xShares' sum up to a value larger than 1" )
    }
  } else if( model == "CondL" ){
    xShares <- mXVal[ xPos, ]
    for( p in 1:pCoef ){
      if( sum( xShares[ , p ] ) > 1 ){
        stop( "the shares in argument 'xShares' sum up to a value larger than 1" ) 
      }
    }
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }  
  if( model == "logit" ){
    if( length( xCoef ) != length( xShares ) ){
      stop( "arguments 'xCoef' and 'xShares' must have the same length" )
    }
  } else if( model == "MNL" ){
    if( dim( xCoef )[1] != length( xShares ) ){
      stop( "arguments 'xCoef' and 'xShares' must have the same length" )
    }
  } else if( model == "CondL" ){
    if( length( xCoef ) != dim( xShares )[1] ){
      stop( "arguments 'xCoef' and 'xShares' must have the same length" )
    }
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }  
  # check argument xGroups
  if( length( xGroups ) != length( xPos ) ){
    stop( "the vector specified by argument 'xGroups' must have",
      " one more element that the vector specified by argument 'xPos'" )
  }
  if( !all( xGroups %in% c( -1, 0, 1 ) ) ){
    stop( "all elements of argument 'xGroups' must be -1, 0, or 1" )
  }
  # check and prepare allCoefVcov
  allCoefVcov <- prepareVcov( allCoefVcov, nCoef, xPos = NA, xMeanSd = NULL )
  if( model == "logit" ){
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
    effeG <- exp( XBetaEffect )/( 1 + exp( XBetaEffect ) ) - 
      exp( XBetaRef )/( 1 + exp( XBetaRef ) )
  } else if( model == "MNL" ){
    # D_mr  
    DRef <- colSums( xCoef[ xGroups == -1, , drop = FALSE ] * 
        xShares[ xGroups == -1 ] )/ 
      sum( xShares[ xGroups == -1 ] )
    XBetaRef <- colSums( mCoef[ -xPos, , drop = FALSE ] * 
        allXVal[ -xPos ]) + DRef
    # D_ml
    DEffect <- colSums( xCoef[ xGroups == 1, , drop = FALSE ] * 
        xShares[ xGroups == 1 ] )/ 
      sum( xShares[ xGroups == 1 ] )
    XBetaEffect <- colSums( mCoef[ -xPos, , drop = FALSE ] * 
        allXVal[ -xPos ]) + DEffect  
    # effect
    effeG <- exp( XBetaEffect[ yCat ] )/( 1 + sum( exp( XBetaEffect ) ) ) -
      exp( XBetaRef[ yCat ] )/( 1 + sum( exp( XBetaRef ) ) )
  } else if( model == "CondL" ){
    # D_mr
    DRef <- colSums( xCoef[ xGroups == -1 ] * 
        xShares[ xGroups == -1, , drop = FALSE ] )/ 
      sum( xShares[ xGroups == -1, , drop = FALSE ] )
    XBetaRef <- colSums( allCoef[ -xPos ] * 
        mXVal[ -xPos, , drop = FALSE ] ) + DRef
    # D_ml
    DEffect <- colSums( xCoef[ xGroups == 1 ] * 
        xShares[ xGroups == 1, , drop = FALSE ] )/ 
      sum( xShares[ xGroups == 1, , drop = FALSE ] )
    XBetaEffect <- colSums( allCoef[ -xPos ] * 
        mXVal[ -xPos, , drop = FALSE ] ) + DEffect
    # effect
    effeG <- exp( XBetaEffect[ yCat ] )/( sum( exp( XBetaEffect ) ) ) -
      exp( XBetaRef[ yCat ] )/( sum( exp( XBetaRef ) ) )
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  # partial derivative of E_{k,ml} w.r.t. all estimated coefficients
  if( model == "logit" ){
    derivCoef <- rep( NA, nCoef )
    derivCoef[ -xPos ] <- ( exp( XBetaEffect )/( 1 + exp( XBetaEffect ))^2 - 
        exp( XBetaRef )/( 1 + exp( XBetaRef ))^2 ) * 
      allXVal[ -xPos ] 
    derivCoef[ xPos ] <- exp( XBetaEffect )/( 1 + exp( XBetaEffect))^2 * DEffect - 
      exp( XBetaRef )/( 1 + exp( XBetaRef ))^2 * DRef
  } else if( model == "MNL" ){
    derivCoef <- matrix( NA, nrow = nXVal, ncol = pCoef )
    for( p in 1:pCoef ){
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
              ( 1 + sum( exp( XBetaEffect ) ) )^2 ) * DEffect -
          ( exp( XBetaRef[ p ] ) * 
              ( 1 + sum( exp( XBetaRef[ -yCat ] ) ) )/
              ( 1 + sum( exp( XBetaRef ) ) )^2 ) * DRef
      } else{  
        derivCoef[ -xPos, p ] <- 
          ( ( exp( XBetaRef[ yCat ] ) * exp( XBetaRef[ p ] ) )/
              ( 1 + sum( exp( XBetaRef ) ) )^2 -
              ( exp( XBetaEffect[ yCat ] ) * exp( XBetaEffect[ p ] ) )/
              ( 1 + sum( exp( XBetaEffect ) ) )^2 ) * allXVal[ -xPos ]
        derivCoef[ xPos, p ] <- 
          ( ( exp( XBetaRef[ yCat ] ) * exp( XBetaRef[ p ] ) )/
              ( 1 + sum( exp( XBetaRef ) ) )^2 ) * DRef -
          ( ( exp( XBetaEffect[ yCat ] ) * exp( XBetaEffect[ p ] ) )/
              ( 1 + sum( exp( XBetaEffect ) ) )^2 ) * DEffect
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
  result <- c( effect = effeG, stdEr = effeGSE )
  return( result )
}
