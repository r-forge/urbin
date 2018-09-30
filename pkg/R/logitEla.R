logitEla <- function( allCoef, allXVal, xPos, 
  allCoefSE = rep( NA, length( allCoef ) ), 
  allCoefBra = NULL, allXValBra = NULL, yCat = NULL, yCatBra = NULL, 
  lambda = NULL, method =  "binary" ){
  
  if( method == "binary" || method == "CondL" ){
    nCoef <- length( allCoef )
    # Checking standard errors  
    if( nCoef != length( allCoefSE ) ) {
      stop( "arguments 'allCoef' and 'allCoefSE' must have the same length" )
    }  
  } else if( method == "MNL" ){
    NCoef <- length( allCoef )
    mCoef <- matrix( allCoef, nrow = length( allXVal ) )
    nCoef <- dim( mCoef )[1]
    pCoef <- dim( mCoef )[2]
    # Checking standard errors
    if( NCoef != length( allCoefSE ) ) {
      stop( "arguments 'allCoef' and 'allCoefSE' must have the same length" )
    }
  } else{
    nCoef <- length( allCoef )
    NCoef <- length( allCoefBra )
    # Checking standard errors
    if( nCoef != length( allCoefSE ) ){
      stop( "arguments 'allCoef' and 'allCoefSE' must have the same length" )
    }
  } 
  # Check position vector
  checkXPos( xPos, minLength = 1, maxLength = 2, minVal = 1, maxVal = nCoef )
  # Check x values
  if( method == "binary" || method == "MNL" ){
    if( nCoef != length( allXVal ) ) {
      stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
    }
  } else if( method == "CondL" ){
    mXVal <- matrix( allXVal, nrow = nCoef )
    nXVal <- dim( mXVal )[1]
    pXVal <- dim( mXVal )[2]
    if( nCoef != dim( mXVal )[1] ) {
      stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
    }
  } else{
    mXValBra <- matrix( allXValBra, nrow = NCoef )
    nXValBra <- dim( mXValBra )[1]
    pXValBra <- dim( mXValBra )[2]
    if( NCoef != nXValBra ) {
      stop( "arguments 'allCoefBra' and 'allXValBra' must have the same length" )
    }
    O <- length( allXVal ) 
    mXVal <- matrix( unlist( allXVal[ yCatBra ] ), nrow = nCoef ) 
    nXVal <- dim( mXVal )[1]
    pXVal <- dim( mXVal )[2]
    if( nCoef != nXVal ) {
      stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
    }
  } 
  # Identify coefficients of interest (kth/tth covariate)
  if( length( xPos ) == 2 ){
    if( method == "binary" ){
      xCoefSE <- allCoefSE[ xPos ]
      xCoef <- allCoef[ xPos ]
      if( !isTRUE( all.equal( allXVal[xPos[2]], allXVal[xPos[1]]^2 ) ) ) {
        stop( "the value of 'allXVal[ xPos[2] ]' must be equal",
          "to the squared value of 'allXVal[ xPos[1] ]' " )
      }
    } else if( method == "MNL" ){
      xCoefSE <- allCoefSE[ xPos ]
      xCoef <- mCoef[ xPos, ]
      if( !isTRUE( all.equal( allXVal[xPos[2]], allXVal[xPos[1]]^2 ) ) ) {
        stop( "the value of 'allXVal[ xPos[2] ]' must be equal",
          "to the squared value of 'allXVal[ xPos[1] ]' " ) 
      }    
    } else if( method == "CondL" ){
      xCoefSE <- allCoefSE[ xPos ]
      xCoef <- allCoef[ xPos ]
      for( p in 1:pXVal ){
        if( !isTRUE( all.equal( mXVal[xPos[2], p], mXVal[xPos[1], p]^2 ) ) ) {
          stop( "the value of 'allXVal[ xPos[2] ]' must be equal",
            "to the squared value of 'allXVal[ xPos[1] ]' " ) 
        }  
      }
    } else{
      xCoefSE <- allCoefSE[ xPos ] 
      xCoef <- allCoef[ xPos ]
      for( p in 1:pXVal ){
        if( !isTRUE( all.equal( mXVal[xPos[2], p], mXVal[xPos[1], p]^2 ) ) ) {
          stop( "the value of 'allXVal[ xPos[2] ]' must be equal",
            "to the squared value of 'allXVal[ xPos[1] ]' " ) 
        }
      }
    }
  } else if( length( xPos ) == 1 ) {
    if( method == "binary" || method == "CondL" ){
      xCoefSE <- c( allCoefSE[ xPos ], 0 )
      xCoef <- c( allCoef[ xPos ], 0 )
    } else if( method == "MNL" ){
      xCoefSE <- c( allCoefSE[ xPos ], 0 )
      xCoef <- matrix( c( mCoef[ xPos, ], rep( 0, dim( mCoef )[ 2 ] ) ), 
        nrow = 2, byrow = TRUE  ) 
    } else{
      xCoefSE <- c( allCoefSE[ xPos ], 0 ) 
      xCoef <- c( allCoef[ xPos ], 0 )
    }    
  } else {
    stop( "argument 'xPos' must be a scalar or a vector with two elements" )
  }
  if( method == "binary" ){
    xVal <- allXVal[ xPos[1] ]
    xBeta <- sum( allCoef * allXVal )
    checkXBeta( xBeta )
    dfun <- exp( xBeta )/( 1 + exp( xBeta ) )^2
    semEla <- ( xCoef[ 1 ] + 2 * xCoef[ 2 ] * xVal ) * xVal * dfun
  } else if( method == "MNL" ){     #checkXBeta missing
    xVal <- allXVal[ xPos[1] ]
    xBeta <- colSums( mCoef * allXVal )
    pfun <- rep( NA, length( xBeta ))
    term <- 0
    for( i in 1:length( xBeta )){
      pfun[i] <- exp( xBeta[i] )/( 1 + sum( exp( xBeta ) ) )
      term <- term + ( ( xCoef[ 1, yCat ] + 2 * xCoef[ 2, yCat ] * xVal ) -
          ( xCoef[ 1, i ] + 2 * xCoef[ 2, i ] * xVal ) * pfun[i] )
    }
    semEla <- xVal * pfun[ yCat ] * 
      ( ( xCoef[ 1, yCat ] + 2 * xCoef[ 2, yCat ] * xVal )/
          ( 1 + sum( exp( xBeta ))) + term )
    dfun <- pfun[ yCat ] * ( 1/( 1 + sum( exp( xBeta ) ) ) + term )
  } else if( method == "CondL" ){    #checkXBeta missing
    xVal <- rep( NA, pXVal )
    for( p in 1:pXVal ){
      xVal[p] <- mXVal[ xPos[ 1 ], p ]
    }
    xBeta <- colSums( allCoef * mXVal )
    pfun <- exp( xBeta[ yCat ] )/( sum( exp( xBeta ) ) )
    semEla <- ( xCoef[1] + 2 * xCoef[2] * xVal[ yCat ] ) * 
      xVal[ yCat ] * ( pfun - pfun^2 )
  } else{                            #checkXBeta missing
    xVal <- rep( NA, pXVal )
    for( p in 1:pXVal ){
      xVal[p] <- mXVal[ xPos[ 1 ], p ]
    }
    coef <- matrix( NA, nrow = O, ncol = nCoef )
    for( o in 1:O ){
      coef[o, ] <- allCoef/lambda[o] 
    }
    xBeta <- lapply( 1:O, function( i, m, v ){ colSums( m[[i]] * v[[i]] ) }, 
      m=allXVal, v=coef )
    IV <- unlist( lapply( 1:O, function( i, m ){ log( sum( exp( m[[i]] ) ) ) }, 
      m=xBeta ) ) 
    pfun <- exp( xBeta[[ yCatBra ]][ yCat ] )/
      ( sum( exp( xBeta[[ yCatBra ]] ) ) )
    xBetaBra <- colSums( allCoefBra * mXValBra )
    pfunBra <- exp( xBetaBra[ yCatBra ] + lambda[ yCatBra ] * IV[ yCatBra ] )/
      ( sum( exp( xBetaBra + lambda * IV ) ) )
    semEla <- ( xCoef[1] + 2 * xCoef[2] * xVal[ yCat ] ) * xVal[ yCat ] * 
      ( pfunBra * ( pfun - pfun^2 ) * 1/lambda[ yCatBra ] +
          pfun^2 * ( pfunBra - pfunBra^2 ) * lambda[ yCatBra ] * IV[ yCatBra ] )
  } 
  if( method == "binary" || method == "MNL" ){
    derivCoef <- c( dfun * xVal, 
      ifelse( length( xPos ) == 1, 0, dfun * 2 * xVal^2 ) )
  } else if( method == "CondL" ){
    derivCoef <- c( ( pfun - pfun^2 ) * xVal[ yCat ], 
      ifelse( length( xPos ) == 1, 0, 
        ( pfun - pfun^2 ) * 2 * xVal[ yCat ]^2 ) )
  } else{
    derivCoef <- c( ( pfunBra * ( pfun - pfun^2 ) * 1/lambda[ yCatBra ] + 
        pfun^2 * ( pfunBra - pfunBra^2 ) * 
        lambda[ yCatBra ] * IV[ yCatBra ] ) *
        xVal[ yCat ], 
      ifelse( length( xPos ) == 1, 0, 
        ( pfunBra * ( pfun - pfun^2 ) * 1/lambda[ yCatBra ] + 
            pfun^2 * ( pfunBra - pfunBra^2 ) * 
            lambda[ yCatBra ] * IV[ yCatBra ] ) *
          2 * xVal[ yCat ]^2 ) )
  }   
  vcovCoef <- diag( xCoefSE^2 )
  semElaSE <- drop( sqrt( t( derivCoef ) %*% vcovCoef %*% derivCoef ) )
  result <- c( semEla = semEla, stdEr = semElaSE )
  return( result )
} 
