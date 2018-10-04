logitEffCat <- function( allCoef, allXVal, xPos, Group, yCat = NA,
  allCoefSE = rep( NA, length( allCoef ) ), 
  method = "binary" ){
  if( method == "binary" ){
    nCoef <- length( allCoef )
    xCoef <- allCoef[ xPos ]
    xShares <- allXVal[ xPos ]
  } else if( method == "MNL" ){
    nCoef <- length( allCoef )
    mCoef <- matrix( allCoef, nrow = length( allXVal ) )
    NCoef <- dim( mCoef )[2]
    pCoef <- dim( mCoef )[1]
    xCoef <- mCoef[ xPos, ]
    xShares <- allXVal[ xPos ]
  } else{
    nCoef <- length( allCoef )
    xCoef <- allCoef[ xPos ]
    mXVal <- matrix( allXVal, nrow = nCoef )
    pCoef <- dim( mXVal )[2]
    xShares <- mXVal[ xPos, ]
  }
  if( method == "binary" || method == "MNL" ){
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
  if( method == "binary" ){
    if( length( xCoef ) != length( xShares ) ){
      stop( "arguments 'xCoef' and 'xShares' must have the same length" )
    }
    if( length( xCoef ) != length( Group ) ){
      stop( "arguments 'xCoef' and 'Group' must have the same length" )
    }
  } else if( method == "MNL" ){
    if( dim( xCoef )[1] != length( xShares ) ){
      stop( "arguments 'xCoef' and 'xShares' must have the same length" )
    }
    if( dim( xCoef )[1] != length( Group ) ){
      stop( "arguments 'xCoef' and 'Group' must have the same length" )
    }
  } else{
    if( length( xCoef ) != dim( xShares )[1] ){
      stop( "arguments 'xCoef' and 'xShares' must have the same length" )
    }
    if( length( xCoef ) != length( Group ) ){
      stop( "arguments 'xCoef' and 'Group' must have the same length" )
    }
  }  
  if( !all( Group %in% c( -1, 0, 1 ) ) ){
    stop( "all elements of argument 'Group' must be -1, 0, or 1" )
  }
  if( method == "binary" ){
    # D_mr  
    DRef <- sum( xCoef[ Group == -1 ] * xShares[ Group == -1 ]) / 
      sum( xShares[ Group == -1 ] )
    XBetaRef <- sum( allCoef[ -xPos ] * allXVal[ -xPos ]) + DRef
    # D_ml
    DEffect <- sum( xCoef[ Group == 1 ] * xShares[ Group == 1 ]) / 
      sum( xShares[ Group == 1 ] )
    XBetaEffect <- sum( allCoef[ -xPos ] * allXVal[ -xPos ]) + DEffect
    # effect
    effeG <- exp( XBetaEffect )/( 1 + exp( XBetaEffect ) ) - 
      exp( XBetaRef )/( 1 + exp( XBetaEffect ) )
  } else if( method == "MNL" ){
    # D_mr  
    DRef <- colSums( xCoef[ Group == -1, , drop = FALSE ] * 
        xShares[ Group == -1 ] )/ 
      sum( xShares[ Group == -1 ] )
    XBetaRef <- colSums( mCoef[ -xPos, , drop = FALSE ] * 
        allXVal[ -xPos ]) + DRef
    # D_ml
    DEffect <- colSums( xCoef[ Group == 1, , drop = FALSE ] * 
        xShares[ Group == 1 ] )/ 
      sum( xShares[ Group == 1 ] )
    XBetaEffect <- colSums( mCoef[ -xPos, , drop = FALSE ] * 
        allXVal[ -xPos ]) + DEffect  
    # effect
    effeG <- exp( XBetaEffect[ yCat ] )/( 1 + sum( exp( XBetaEffect ) ) ) -
      exp( XBetaRef[ yCat ] )/( 1 + sum( exp( XBetaRef ) ) )
  } else{
    # D_mr
    DRef <- colSums( xCoef[ Group == -1 ] * 
        xShares[ Group == -1, , drop = FALSE ] )/ 
      sum( xShares[ Group == -1, , drop = FALSE ] )
    XBetaRef <- colSums( allCoef[ -xPos ] * 
        mXVal[ -xPos, , drop = FALSE ] ) + DRef
    # D_ml
    DEffect <- colSums( xCoef[ Group == 1 ] * 
        xShares[ Group == 1, , drop = FALSE ] )/ 
      sum( xShares[ Group == 1, , drop = FALSE ] )
    XBetaEffect <- colSums( allCoef[ -xPos ] * 
        mXVal[ -xPos, , drop = FALSE ] ) + DEffect
    # effect
    effeG <- exp( XBetaEffect[ yCat ] )/( sum( exp( XBetaEffect ) ) ) -
      exp( XBetaRef[ yCat ] )/( sum( exp( XBetaRef ) ) )
  }
  # partial derivative of E_{k,ml} w.r.t. all estimated coefficients
  if( method == "binary" ){
    derivCoef <- rep( NA, nCoef )
    derivCoef[ -xPos ] = ( exp( XBetaEffect )/( 1 + exp( XBetaEffect ))^2 - 
        exp( XBetaRef )/( 1 + exp( XBetaRef ))^2 ) * 
      allXVal[ -xPos ] 
    derivCoef[ xPos ] = exp( XBetaEffect )/( 1 + exp( XBetaEffect))^2 * DEffect - 
      exp( XBetaRef )/( 1 + exp( XBetaRef ))^2 * DRef
  } else if( method == "MNL" ){
    derivCoef <- matrix( NA, nrow=pCoef, ncol=NCoef )
    for( p in 1:NCoef ){
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
  result <- c( effect = effeG, stdEr = effeGSE )
  return( result )
}
