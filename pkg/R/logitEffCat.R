logitEffCat <- function( allCoef, allXVal, xPos, xGroups, model = "logit",
  yCat = NA, allCoefSE = rep( NA, length( allCoef ) ) ){
  
  # number of coefficients
  nCoef <- length( allCoef )
  # number of explanatory variables
  nXVal <- length( allXVal )
  # Check position vector
  checkXPos( xPos, minLength = 1, maxLength = nCoef, minVal = 1, 
    maxVal = nCoef )
  if( model == "logit" ){
    xCoef <- allCoef[ xPos ]
    xShares <- allXVal[ xPos ]
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
    xShares <- allXVal[ xPos ]
  } else{
    # number of ???
    pCoef <- round( nXVal / nCoef )
    if( nXVal != nCoef * pCoef ) {
      stop( "length of argument 'allXVal' must be a multiple",
        " of the length of argument 'allCoef'" )
    } 
    # create matrix of explanatory variables
    mXVal <- matrix( allXVal, nrow = nCoef )
    xCoef <- allCoef[ xPos ]
    xShares <- mXVal[ xPos, ]
  }
  if( model == "logit" || model == "MNL" ){
    if( sum( xShares ) > 1 ){
      stop( "the shares in argument 'xShares' sum up to a value larger than 1" )
    }
  } else{
    for( p in 1:pCoef ){
      if( sum( xShares[ , p ] ) > 1 ){
        stop( "the shares in argument 'xShares' sum up to a value larger than 1" ) 
      }
    }
  }  
  if( model == "logit" ){
    if( length( xCoef ) != length( xShares ) ){
      stop( "arguments 'xCoef' and 'xShares' must have the same length" )
    }
    if( length( xCoef ) != length( xGroups ) ){
      stop( "arguments 'xCoef' and 'xGroups' must have the same length" )
    }
  } else if( model == "MNL" ){
    if( dim( xCoef )[1] != length( xShares ) ){
      stop( "arguments 'xCoef' and 'xShares' must have the same length" )
    }
    if( dim( xCoef )[1] != length( xGroups ) ){
      stop( "arguments 'xCoef' and 'xGroups' must have the same length" )
    }
  } else{
    if( length( xCoef ) != dim( xShares )[1] ){
      stop( "arguments 'xCoef' and 'xShares' must have the same length" )
    }
    if( length( xCoef ) != length( xGroups ) ){
      stop( "arguments 'xCoef' and 'xGroups' must have the same length" )
    }
  }  
  if( !all( xGroups %in% c( -1, 0, 1 ) ) ){
    stop( "all elements of argument 'xGroups' must be -1, 0, or 1" )
  }
  if( model == "logit" ){
    # D_mr  
    DRef <- sum( xCoef[ xGroups == -1 ] * xShares[ xGroups == -1 ]) / 
      sum( xShares[ xGroups == -1 ] )
    XBetaRef <- sum( allCoef[ -xPos ] * allXVal[ -xPos ]) + DRef
    # D_ml
    DEffect <- sum( xCoef[ xGroups == 1 ] * xShares[ xGroups == 1 ]) / 
      sum( xShares[ xGroups == 1 ] )
    XBetaEffect <- sum( allCoef[ -xPos ] * allXVal[ -xPos ]) + DEffect
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
  } else{
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
  }
  # partial derivative of E_{k,ml} w.r.t. all estimated coefficients
  if( model == "logit" ){
    derivCoef <- rep( NA, nCoef )
    derivCoef[ -xPos ] = ( exp( XBetaEffect )/( 1 + exp( XBetaEffect ))^2 - 
        exp( XBetaRef )/( 1 + exp( XBetaRef ))^2 ) * 
      allXVal[ -xPos ] 
    derivCoef[ xPos ] = exp( XBetaEffect )/( 1 + exp( XBetaEffect))^2 * DEffect - 
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
  } else{
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
  }
  # variance covariance of the coefficients (covariances set to zero)
  vcovCoef <- diag( allCoefSE^2 )
  # approximate standard error of the effect
  effeGSE <- drop( sqrt( t( derivCoef ) %*% vcovCoef %*% derivCoef ) )
  # object to be returned
  result <- c( effect = effeG, stdEr = effeGSE )
  return( result )
}
