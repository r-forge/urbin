urbinEla <- function( allCoef, allXVal, xPos, model, 
  allCoefVcov = NULL, seSimplify = !is.matrix( allCoefVcov ), 
  xMeanSd = NULL, yCat = NULL, 
  allCoefBra = NULL, allXValBra = NULL, yCatBra = NULL, lambda = NULL ){
  
  # check argument seSimplify
  if( length( seSimplify ) != 1 || !is.logical( seSimplify ) ) {
    stop( "argument 'seSimplify' must be TRUE or FALSE" )
  }
  if( model == "lpm" && !is.null( allCoefVcov ) && 
      ( seSimplify == is.matrix( allCoefVcov ) ) ) {
    warning( "argument 'seSimplify' is ignored in linear probability models" )    
  } else if( seSimplify && is.matrix( allCoefVcov ) ) {
    warning( "if the (full) covariance matrix is provided",
      " via argument 'allCoefCov'", 
      " it is NOT recommended to set argument 'seSimplify' to TRUE" )    
  } else if( !seSimplify && !is.null( allCoefVcov) && 
      !is.matrix( allCoefVcov ) ) {
    warning( "the returned standard error is likely upward biased;",
      " you can provide the full covariance matrix", 
      " via argument 'allCoefVcov' to avoid this bias",
      " or do NOT set argument 'seSimplify' to FALSE" )
  }
  # number of coefficients
  nCoef <- length( allCoef )
  # number of explanatory variables
  nXVal <- length( allXVal )
  # check allXVal and allCoef
  if( model %in% c( "lpm", "probit", "logit" ) ){
    # LPM model: allXVal can be a scalar even if there is a quadratic term
    if( model == "lpm" && length( xPos ) == 2 && length( allXVal ) == 1 ){
      temp <- c( allXVal, allXVal^2 )
      if( "derivOnly" %in% names( attributes( allXVal ) ) ) {
        attr( temp, "derivOnly" ) <- 1
      }
      allXVal <- temp
      nXVal <- length( allXVal )
    }
    if( nXVal != nCoef ) {
      stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
    } 
  } else if( model == "MNL" ){
    # number of alternative categories of the dependent variable
    pCoef <- round( nCoef / nXVal )
    if( nCoef != nXVal * pCoef ) {
      stop( "length of argument 'allCoef' must be a multiple",
        " of the length of argument 'allXVal'" )
    } 
    # create matrix of coefficients
    mCoef <- matrix( allCoef, nrow = nXVal, ncol = pCoef )
  } else if( model == "CondL" ){
    mXVal <- matrix( allXVal, nrow = nCoef )
    nXVal <- dim( mXVal )[1]
    pXVal <- dim( mXVal )[2]
    if( nCoef != dim( mXVal )[1] ) {
      stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
    }
  } else if( model == "NestedL" ){
    NCoef <- length( allCoefBra )
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
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  # Check position vector
  checkXPos( xPos, minLength = 1, maxLength = 2, minVal = 1, 
    maxVal = ifelse( model == "MNL", nXVal, nCoef ) )
  # check and prepare allCoefVcov
  allCoefVcov <- prepareVcov( allCoefVcov, length( allCoef ), xPos, xMeanSd )
  # Identify coefficients of interest (kth/tth covariate)
  if( length( xPos ) == 2 ){
    if( model %in% c( "lpm", "probit", "logit" ) ){
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
    } else if( model %in% c( "CondL", "NestedL" ) ){
      xCoef <- allCoef[ xPos ]
      for( p in 1:pXVal ){
        if( !isTRUE( all.equal( mXVal[xPos[2], p], mXVal[xPos[1], p]^2 ) ) ) {
          stop( "the value of 'allXVal[ xPos[2] ]' must be equal",
            " to the squared value of 'allXVal[ xPos[1] ]' " ) 
        }  
      }
    } else {
      stop( "argument 'model' specifies an unknown type of model" )
    }
  } else if( length( xPos ) == 1 ) {
    if( model %in% c( "lpm", "probit", "logit", "CondL", "NestedL" ) ){
      xCoef <- c( allCoef[ xPos ], 0 )
    } else if( model == "MNL" ){
      xCoef <- matrix( c( mCoef[ xPos, ], rep( 0, dim( mCoef )[ 2 ] ) ), 
        nrow = 2, byrow = TRUE  ) 
    } else {
      stop( "argument 'model' specifies an unknown type of model" )
    }    
  } else {
    stop( "argument 'xPos' must be a scalar or a vector with two elements" )
  }
  # prepare calculation of semi-elasticity 
  if( model == "lpm" ) {
    xVal <- allXVal[ xPos[ 1 ] ]
    semEla <- ( xCoef[1] + 2 * xCoef[2] * xVal ) * xVal
  } else if( model == "probit" ) {
    xVal <- allXVal[ xPos[ 1 ] ]
    xBeta <- sum( allCoef * allXVal )
    checkXBeta( xBeta )
    dfun <- dnorm( xBeta )
    semEla <- ( xCoef[ 1 ] + 2 * xCoef[ 2 ] * xVal ) * xVal * dfun
  } else if( model == "logit" ){
    xVal <- allXVal[ xPos[1] ]
    xBeta <- sum( allCoef * allXVal )
    checkXBeta( xBeta )
    dfun <- exp( xBeta )/( 1 + exp( xBeta ) )^2
    semEla <- ( xCoef[ 1 ] + 2 * xCoef[ 2 ] * xVal ) * xVal * dfun
  } else if( model == "MNL" ){     #checkXBeta missing
    xVal <- allXVal[ xPos[1] ]
    xCoefLinQuad <- xCoef[ 1, ] + 2 * xCoef[ 2, ] * xVal
    xBeta <- allXVal %*% mCoef
    pfun <- exp( xBeta ) / ( 1 + sum( exp( xBeta ) ) )
    term <- 0
    for( i in 1:length( xBeta )){
      term <- term + ( ( xCoef[ 1, yCat ] + 2 * xCoef[ 2, yCat ] * xVal ) -
          ( xCoef[ 1, i ] + 2 * xCoef[ 2, i ] * xVal ) * pfun[i] )
    }
    semEla <- xVal * pfun[ yCat ] * 
      ( xCoefLinQuad[ yCat ] - sum( xCoefLinQuad * pfun ) )
    dfun <- pfun[ yCat ] * ( 1/( 1 + sum( exp( xBeta ) ) ) + term )
  } else if( model == "CondL" ){    #checkXBeta missing
    xVal <- rep( NA, pXVal )
    for( p in 1:pXVal ){
      xVal[p] <- mXVal[ xPos[ 1 ], p ]
    }
    xBeta <- allCoef %*% mXVal
    pfun <- exp( xBeta[ yCat ] )/( sum( exp( xBeta ) ) )
    semEla <- ( xCoef[1] + 2 * xCoef[2] * xVal[ yCat ] ) * 
      xVal[ yCat ] * ( pfun - pfun^2 )
  } else if( model == "NestedL" ){                            #checkXBeta missing
    xVal <- rep( NA, pXVal )
    for( p in 1:pXVal ){
      xVal[p] <- mXVal[ xPos[ 1 ], p ]
    }
    coef <- matrix( NA, nrow = O, ncol = nCoef )
    for( o in 1:O ){
      coef[o, ] <- allCoef/lambda[o] 
    }
    xBeta <- lapply( 1:O, function( i, m, v ){ colSums( m[[i]] * v[[i]] ) }, 
      m=allXVal, v=coef )   #### v[[i]] is probably incorrect, because v=coef is not a list
    IV <- unlist( lapply( 1:O, function( i, m ){ log( sum( exp( m[[i]] ) ) ) }, 
      m=xBeta ) ) 
    pfun <- exp( xBeta[[ yCatBra ]][ yCat ] )/
      ( sum( exp( xBeta[[ yCatBra ]] ) ) )
    xBetaBra <- allCoefBra %*% mXValBra
    pfunBra <- exp( xBetaBra[ yCatBra ] + lambda[ yCatBra ] * IV[ yCatBra ] )/
      ( sum( exp( xBetaBra + lambda * IV ) ) )
    semEla <- ( xCoef[1] + 2 * xCoef[2] * xVal[ yCat ] ) * xVal[ yCat ] * 
      ( pfunBra * ( pfun - pfun^2 ) * 1/lambda[ yCatBra ] +
          pfun^2 * ( pfunBra - pfunBra^2 ) * lambda[ yCatBra ] * IV[ yCatBra ] )
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  } 
  # partial derivatives of semi-elasticities wrt coefficients
  if( model == "lpm" ){
    derivCoef <- rep( 0, nCoef ) 
    derivCoef[ xPos[1] ] <- xVal
    if( length( xPos ) == 2 ) {
      derivCoef[ xPos[2] ] <- 2 * xVal^2
    }
  } else if( model == "probit" ){
    if( seSimplify ) {
      derivCoef <- rep( 0, length( allCoef ) )
    } else {
      derivCoef <- ddnorm( xBeta ) * allXVal * 
        ( xCoef[ 1 ] + 2 * xCoef[ 2 ] * xVal ) * xVal
    }
    derivCoef[ xPos[1] ] <- derivCoef[ xPos[1] ] + dnorm( xBeta ) * xVal
    if( length( xPos ) == 2 ) {
      derivCoef[ xPos[2] ] <- derivCoef[ xPos[2] ] + dnorm( xBeta ) * 2 * xVal^2
    }
  } else if( model == "logit" ){
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
  } else if( model == "NestedL" ){
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
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }   
  # if argument allXVal has attribute 'derivOnly',
  # return partial derivatives only (for testing partial derivatives)
  if( "derivOnly" %in% names( attributes( allXVal ) ) ) {
    return( c( derivCoef ) )
  }
  # approximate standard error of the semi-elasticity
  semElaSE <- drop( sqrt( t( derivCoef ) %*% allCoefVcov %*% derivCoef ) )
  # create object that will be returned
  result <- c( semEla = semEla, stdEr = semElaSE )
  return( result )
} 
