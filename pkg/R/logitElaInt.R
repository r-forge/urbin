logitElaInt <- function( allCoef, allXVal, xPos, xBound, yCat = NA,
  allCoefVcov = NULL, 
  model = "logit" ){
  # number of coefficients
  nCoef <- length( allCoef )
  # number of explanatory variables
  nXVal <- length( allXVal )
  # check allXVal and allCoef
  if( model == "logit" ){
    if( nXVal != nCoef ) {
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
  } else if( model == "CondL" ) {
    # number of ???
    pCoef <- round( nXVal / nCoef )
    if( nXVal != nCoef * pCoef ) {
      stop( "length of argument 'allXVal' must be a multiple",
        " of the length of argument 'allCoef'" )
    } 
    # create matrix of explanatory variables
    mXVal <- matrix( allXVal, nrow = nCoef )
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  # Check position vector
  checkXPos( xPos, minLength = 2, maxLength = nCoef, 
    minVal = 0, maxVal = nCoef, requiredVal = 0 )
  # number of intervals
  nInt <- length( xPos ) 
  # check 'xBound' and replace infinite values
  xBound <- elaIntBounds( xBound, nInt )
  # check and prepare allCoefVcov
  allCoefVcov <- prepareVcov( allCoefVcov, length( allCoef ), xPos, xMeanSd = NULL )
  # vector of shares of observations in each interval
  shareVec <- rep( NA, nInt )
  for( i in 1:nInt ){
    if( xPos[i] != 0 ) {
      if( model == "logit" || model == "MNL" ){
        shareVec[i] <- allXVal[ xPos[i] ] 
      } else if( model == "CondL" ) {
        shareVec[i] <- mXVal[ xPos[i], yCat ]
      } else {
        stop( "argument 'model' specifies an unknown type of model" )
      }
    }
  }
  if( any( shareVec[ xPos != 0 ] < 0 ) ) {
    stop( "all elements of argument 'allXVal'",
      " that are indicated by argument 'xPos'",
      " (i.e., the shares of observations in each interval)",
      " must be non-negative" )
  }
  if( sum( shareVec[ xPos != 0 ] ) > 1 + sqrt(.Machine$double.eps) ) {
    stop( "the sum of the elements of argument 'allXVal'",
      " that are indicated by argument 'xPos'",
      " (i.e., the shares of observations in each interval)",
      " must not be larger than one" )
  }
  shareVec[ xPos == 0 ] <- 1 - sum( shareVec[ xPos != 0 ] )  
  # weights
  weights <- elaIntWeights( shareVec )
  # prepare calculation of semi-elasticity 
  if( model == "logit" ){
    xBeta <- rep( NA, nInt ) 
    for( i in 1:nInt ){
      allXValTemp <- replace( allXVal, xPos, 0 )
      if( xPos[i] != 0 ) {
        allXValTemp[ xPos[i] ] <- 1
      }
      xBeta[i] <- sum( allCoef * allXValTemp )
    }
    checkXBeta( xBeta )
  } else if( model == "MNL" ){
    xBeta <- matrix( rep( rep( NA, nInt ), pCoef ), ncol = pCoef ) 
    for( p in 1:pCoef ){
      for( i in 1:nInt ){
        allXValTemp <- replace( allXVal, xPos, 0 )
        if( xPos[i] != 0 ) {
          allXValTemp[ xPos[i] ] <- 1
        }
        xBeta[i,p] <- sum( mCoef[ ,p] * allXValTemp )
      }
      checkXBeta( xBeta[,p] )
    }
  } else if( model == "CondL" ) {
    xBeta <- matrix( rep( rep( NA, nInt ), pCoef ), ncol = pCoef ) 
    for( p in 1:pCoef ){
      for( i in 1:nInt ){
        allXValTemp <- replace( mXVal[ ,p], xPos, 0 )
        if( xPos[i] != 0 ) {
          allXValTemp[ xPos[i] ] <- 1 
        }
        xBeta[i,p] <- sum( allCoef * allXValTemp )    
      }
      checkXBeta( xBeta[,p] )
    }
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  # vector of probabilities of y=1 for each interval
  if( model == "logit" ){
    xCoef <- as.vector( exp( xBeta )/( 1 + exp( xBeta ) ) )  
  } else if( model == "MNL" ){
    xCoef <- as.vector( exp( xBeta[ , yCat ])/( 1 + rowSums( exp( xBeta ) ) ) )
  } else if( model == "CondL" ){
    xCoef <- as.vector( exp( xBeta[ , yCat ])/( rowSums( exp( xBeta ) ) ) )
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  # semi-elasticities 'around' each inner boundary and their weights
  semElaBound <- rep( NA, nInt - 1 )
  for( m in 1:(nInt-1) ){
    semElaBound[m] <- 2 * ( xCoef[ m+1 ] - xCoef[ m ] ) * xBound[ m+1 ] /
      ( xBound[m+2] - xBound[m] )
  }
  # (average) semi-elasticity
  semEla <- sum( semElaBound * weights )
  # partial derivatives of semi-elasticities wrt coefficients
  if( model == "logit" ){
    gradM <- matrix( 0, nCoef, nInt - 1 )
    gradExpVec <- exp( xBeta )/( 1 + exp( xBeta ) )^2
    for( m in 1:( nInt - 1 ) ) {
      gradM[ -xPos, m ] <- 2 * ( gradExpVec[m+1] - gradExpVec[m] ) * 
        allXVal[ -xPos ] * xBound[m+1] / ( xBound[m+2] - xBound[m] )
      gradM[ xPos[m], m ] <- - 2 * gradExpVec[m] * xBound[m+1] / 
        ( xBound[m+2] - xBound[m] )
      gradM[ xPos[m+1], m ] <- 2 * gradExpVec[m+1] * xBound[m+1] / 
        ( xBound[m+2] - xBound[m] )
    } 
  } else if( model == "MNL" ){
    gradM <- array( 0, c( nXVal, nInt - 1, pCoef ) )
    gradExpVecP <- ( exp( xBeta[ , yCat ] ) * 
        ( 1 + rowSums( exp( xBeta[ , -yCat, drop = FALSE ] ) ) ) )/
      ( 1 + rowSums( exp( xBeta ) ) )^2 
    for( p in 1:pCoef ){
      gradExpVecO <- ( exp( xBeta[ , yCat ] ) * exp( xBeta[ , p] ) )/
        ( 1 + rowSums( exp( xBeta ) ) )^2
      for( m in 1:( nInt - 1 ) ) {
        if( p == yCat ){
          gradM[ -xPos, m, p ] <- 2 * ( gradExpVecP[m+1] - gradExpVecP[m] ) *
            allXVal[ -xPos ] * xBound[m+1] / ( xBound[m+2] - xBound[m] )
          gradM[ xPos[ m ], m, p ] <- - 2 * gradExpVecP[m] * xBound[m+1] / 
            ( xBound[m+2] - xBound[m] )
          gradM[ xPos[ m + 1 ], m, p ] <- 2 * gradExpVecP[m+1] * xBound[m+1] / 
            ( xBound[m+2] - xBound[m] )
        } else {
          gradM[ -xPos, m, p ] <- 2 * ( gradExpVecO[m] - gradExpVecO[m+1] ) *
            allXVal[ -xPos ] * xBound[m+1] / ( xBound[m+2] - xBound[m] )
          gradM[ xPos[ m ], m, p ] <- 2 * gradExpVecO[m] * xBound[m+1] / 
            ( xBound[m+2] - xBound[m] )
          gradM[ xPos[ m + 1 ], m, p ] <- - 2 * gradExpVecO[m+1] * xBound[m+1] / 
            ( xBound[m+2] - xBound[m] )
        }  
      }
    }
    gradM <- apply( gradM, 2, function( x ) x )
  } else if( model == "CondL" ){
    gradM <- matrix( 0, nCoef, nInt - 1 )
    for( m in 1:( nInt - 1 ) ) {
      gradM[ -xPos, m ] <- 2 * 
        ( ( exp( xBeta[ m+1, yCat ] ) * mXVal[ -xPos, yCat ] * 
            sum( exp( xBeta[ m+1, ] ) ) - 
            exp( xBeta[ m+1, yCat ] ) * 
            rowSums( exp( xBeta[ m+1, ] ) * mXVal[ -xPos, , drop = FALSE ] ) )/
            ( sum( exp( xBeta[ m+1, ] ) ) )^2 -
            ( exp( xBeta[ m, yCat ] ) * mXVal[ -xPos, yCat ] * 
                sum( exp( xBeta[ m, ] ) ) - 
                exp( xBeta[ m, yCat ] ) * 
                rowSums( exp( xBeta[ m,  ] ) * mXVal[ -xPos, , drop = FALSE ] ) )/
            ( sum( exp( xBeta[ m, ] ) ) )^2 ) * 
        xBound[m+1] / ( xBound[m+2] - xBound[m] )
      gradM[ xPos[m], m ] <- 0
      gradM[ xPos[m+1], m ] <- 0
    } 
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  # partial derivative of the semi-elasticity 
  # w.r.t. all estimated coefficients
  derivCoef <- rep( 0, length( allCoef ) )
  for( m in 1:( nInt - 1 ) ){
    derivCoef <- derivCoef + weights[m] * gradM[,m]
  }
  # if argument allXVal has attribute 'derivOnly',
  # return partial derivatives only (for testing partial derivatives)
  if( "derivOnly" %in% names( attributes( allXVal ) ) ) {
    return( derivCoef )
  }
  # approximate standard error of the semi-elasticity
  semElaSE <- drop( sqrt( t( derivCoef ) %*% allCoefVcov %*% derivCoef ) )
  # create object that will be returned
  result <- c( semEla = semEla, stdEr = semElaSE )
  return( result )
}
