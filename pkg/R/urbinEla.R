urbinEla <- function( allCoef, allXVal, xPos, model, allCoefVcov = NULL, 
  allCoefBra = NULL, allXValBra = NULL, yCat = NULL, yCatBra = NULL, 
  lambda = NULL, 
  seSimplify = !is.matrix( allCoefVcov ), xMeanSd = NULL ){
  
  # check argument seSimplify
  if( length( seSimplify ) != 1 || !is.logical( seSimplify ) ) {
    stop( "argument 'seSimplify' must be TRUE or FALSE" )
  }
  # number(s) of coefficients
  if( model == "binary" || model == "CondL" ){
    nCoef <- length( allCoef )
    # Checking standard errors  
  } else if( model == "MNL" ){
    NCoef <- length( allCoef )
    mCoef <- matrix( allCoef, nrow = length( allXVal ) )
    nCoef <- dim( mCoef )[1]
    pCoef <- dim( mCoef )[2]
  } else{
    nCoef <- length( allCoef )
    NCoef <- length( allCoefBra )
  } 
  # Check position vector
  checkXPos( xPos, minLength = 1, maxLength = 2, minVal = 1, maxVal = nCoef )
  # Check x values
  if( model == "binary" || model == "MNL" ){
    if( nCoef != length( allXVal ) ) {
      stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
    }
  } else if( model == "CondL" ){
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
  # check and prepare allCoefVcov
  allCoefVcov <- prepareVcov( allCoefVcov, length( allCoef ), xPos, xMeanSd )
  # Identify coefficients of interest (kth/tth covariate)
  if( length( xPos ) == 2 ){
    if( model == "binary" ){
      xCoef <- allCoef[ xPos ]
      if( !isTRUE( all.equal( allXVal[xPos[2]], allXVal[xPos[1]]^2 ) ) ) {
        stop( "the value of 'allXVal[ xPos[2] ]' must be equal",
          "to the squared value of 'allXVal[ xPos[1] ]' " )
      }
    } else if( model == "MNL" ){
      xCoef <- mCoef[ xPos, ]
      if( !isTRUE( all.equal( allXVal[xPos[2]], allXVal[xPos[1]]^2 ) ) ) {
        stop( "the value of 'allXVal[ xPos[2] ]' must be equal",
          "to the squared value of 'allXVal[ xPos[1] ]' " ) 
      }    
    } else if( model == "CondL" ){
      xCoef <- allCoef[ xPos ]
      for( p in 1:pXVal ){
        if( !isTRUE( all.equal( mXVal[xPos[2], p], mXVal[xPos[1], p]^2 ) ) ) {
          stop( "the value of 'allXVal[ xPos[2] ]' must be equal",
            "to the squared value of 'allXVal[ xPos[1] ]' " ) 
        }  
      }
    } else{
      xCoef <- allCoef[ xPos ]
      for( p in 1:pXVal ){
        if( !isTRUE( all.equal( mXVal[xPos[2], p], mXVal[xPos[1], p]^2 ) ) ) {
          stop( "the value of 'allXVal[ xPos[2] ]' must be equal",
            "to the squared value of 'allXVal[ xPos[1] ]' " ) 
        }
      }
    }
  } else if( length( xPos ) == 1 ) {
    if( model == "binary" || model == "CondL" ){
      xCoef <- c( allCoef[ xPos ], 0 )
    } else if( model == "MNL" ){
      xCoef <- matrix( c( mCoef[ xPos, ], rep( 0, dim( mCoef )[ 2 ] ) ), 
        nrow = 2, byrow = TRUE  ) 
    } else{
      xCoef <- c( allCoef[ xPos ], 0 )
    }    
  } else {
    stop( "argument 'xPos' must be a scalar or a vector with two elements" )
  }
  # prepare calculation of semi-elasticity 
  if( model == "binary" ){
    xVal <- allXVal[ xPos[1] ]
    xBeta <- sum( allCoef * allXVal )
    checkXBeta( xBeta )
    dfun <- exp( xBeta )/( 1 + exp( xBeta ) )^2
    semEla <- ( xCoef[ 1 ] + 2 * xCoef[ 2 ] * xVal ) * xVal * dfun
  } else if( model == "MNL" ){     #checkXBeta missing
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
  } else if( model == "CondL" ){    #checkXBeta missing
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
  # partial derivatives of semi-elasticities wrt coefficients
  if( model == "binary" ){
    if( seSimplify ) {
      derivCoef <- rep( 0, length( allCoef ) )
    } else {
      derivCoef <- allXVal * semEla *
        ( 1 - 2 * exp( xBeta ) / ( 1 + exp( xBeta ) ) )
    }
    derivCoef[ xPos[1] ] <- derivCoef[ xPos[1] ] + dfun * xVal
    if( length( xPos ) == 2 ) {
      derivCoef[ xPos[2] ] <- derivCoef[ xPos[2] ] + dfun * 2 * xVal^2
    }
  } else if( model == "MNL" ){
    if( !seSimplify ) {
      warning( "exact (non-simplified) calculation of derivatives of",
        " semi-elasticities wrt coefficients has not yet been implemented",
        " for MNL models" )
    }
    derivCoef <- rep( 0, length( allCoef ) )
    derivCoef[ xPos[1] ] <- dfun * xVal
    if( length( xPos ) == 2 ) {
      derivCoef[ xPos[2] ] <- dfun * 2 * xVal^2
    }
  } else if( model == "CondL" ){
    if( !seSimplify ) {
      warning( "exact (non-simplified) calculation of derivatives of",
        " semi-elasticities wrt coefficients has not yet been implemented",
        " for CondL models" )
    }
    derivCoef <- rep( 0, length( allCoef ) )
    derivCoef[ xPos[1] ] <- ( pfun - pfun^2 ) * xVal[ yCat ]
    if( length( xPos ) == 2 ) {
      derivCoef[ xPos[2] ] <- ( pfun - pfun^2 ) * 2 * xVal[ yCat ]^2
    }
  } else{
    if( !seSimplify ) {
      warning( "exact (non-simplified) calculation of derivatives of",
        " semi-elasticities wrt coefficients has not yet been implemented",
        " for NestL models" )
    }
    derivCoef <- rep( 0, length( allCoef ) )
    derivCoef[ xPos[1] ] <- ( 
      pfunBra * ( pfun - pfun^2 ) / lambda[ yCatBra ] +
        pfun^2 * ( pfunBra - pfunBra^2 ) * lambda[ yCatBra ] * 
        IV[ yCatBra ] ) *
      xVal[ yCat ]
    if( length( xPos ) == 2 ) {
      derivCoef[ xPos[2] ] <- ( 
        pfunBra * ( pfun - pfun^2 ) / lambda[ yCatBra ] + 
          pfun^2 * ( pfunBra - pfunBra^2 ) * lambda[ yCatBra ] * 
          IV[ yCatBra ] ) *
        2 * xVal[ yCat ]^2
    }
  }   
  # if argument allXVal has attribute 'derivOnly',
  # return partial derivatives only (for testing partial derivatives)
  if( "derivOnly" %in% names( attributes( allXVal ) ) ) {
    return( c( derivCoef ) )
  }
  # approximate standard error of the semi-elasticity
  semElaSE <- drop( sqrt( t( derivCoef ) %*% allCoefVcov %*% derivCoef ) )
  result <- c( semEla = semEla, stdEr = semElaSE )
  return( result )
} 
