logitElaInt <- function( allCoef, allXVal, xPos, xBound, yCat = NA,
  allCoefSE = rep( NA, length( allCoef ) ), 
  method = "binary" ){
  # number of coefficients
  if( method == "binary" || method == "MNL" ){
    mCoef <- matrix( allCoef, nrow = length( allXVal ))
    nCoef <- dim( mCoef )[1]
    pCoef <- dim( mCoef )[2]
    # checking arguments
    if( length( allXVal ) != nCoef ) {
      stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
    } 
  } else{
    nCoef <- length( allCoef )
    mXVal <- matrix( allXVal, nrow = nCoef )
    pCoef <- dim( mXVal )[2]
    # checking arguments
    if( dim( mXVal )[1] != nCoef ) {
      stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
    } 
  }
  # number of intervals
  nInt <- length( xPos ) 
  checkXPos( xPos, minLength = 2, maxLength = nCoef, 
    minVal = 0, maxVal = nCoef, requiredVal = 0 )
  if( method == "binary" || method == "MNL" ){
    if( any( allXVal[ xPos ] < 0 ) ) {
      stop( "all elements of argument 'allXVal'",
        " that are indicated by argument 'xPos'",
        " (i.e., the shares of observations in each interval)",
        " must be non-negative" )
    }
    if( sum( allXVal[ xPos ] > 1 ) ) {
      stop( "the sum of the elements of argument 'allXVal'",
        " that are indicated by argument 'xPos'",
        " (i.e., the shares of observations in each interval)",
        " must not be larger than one" )
    }
  } else{
    for( p in 1:pCoef ){
      if( any( mXVal[ xPos, p ] < 0 ) ) {
        stop( "all elements of argument 'allXVal'",
          " that are indicated by argument 'xPos'",
          " (i.e., the shares of observations in each interval)",
          " must be non-negative" )
      }
      if( sum( mXVal[ xPos, p ] > 1 ) ) {
        stop( "the sum of the elements of argument 'allXVal'",
          " that are indicated by argument 'xPos'",
          " (i.e., the shares of observations in each interval)",
          " must not be larger than one" )
      }
    }  
  }  
  # check 'xBound' and replace infinite values
  xBound <- elaIntBounds( xBound, nInt )
  # vector of probabilities of y=1 for each interval and
  # vector of shares of observations in each interval
  xBeta <- matrix( rep( rep( NA, nInt ), pCoef ), ncol = pCoef ) 
  if( method == "binary" || method == "MNL" ){
    shareVec <- rep( NA, nInt )
    for( p in 1:pCoef ){
      for( i in 1:nInt ){
        allXValTemp <- replace( allXVal, xPos, 0 )
        if( xPos[i] != 0 ) {
          allXValTemp[ xPos[i] ] <- 1
          shareVec[i] <- allXVal[ xPos[i] ] 
        }
        xBeta[i,p] <- sum( mCoef[ ,p] * allXValTemp )
      }
    }
    shareVec[ xPos == 0 ] <- 1 - sum( shareVec[ xPos != 0 ] )  
  } else{
    shareVec <- matrix( rep( rep( NA, nInt ), pCoef ), ncol = pCoef )
    for( p in 1:pCoef ){
      for( i in 1:nInt ){
        allXValTemp <- replace( mXVal[ ,p], xPos, 0 )
        if( xPos[i] != 0 ) {
          allXValTemp[ xPos[i] ] <- 1 
          shareVec[i,p] <- mXVal[ xPos[i], p ]
        }
        xBeta[i,p] <- sum( allCoef * allXValTemp )    
      }
      shareVec[ xPos == 0, p ] <- 1 - sum( shareVec[ xPos != 0, p ] )
    }
    shareVec <- shareVec[ , yCat ]  
  }  
  #checkXBeta( xBeta )  #Please check this one with a matrix
  if( method == "binary" ){
    expVec <- as.vector( exp( xBeta )/( 1 + exp( xBeta ) ) )  
  } else if( method == "MNL" ){
    expVec <- as.vector( exp( xBeta[ , yCat ])/( 1 + rowSums( exp( xBeta ) ) ) )
  } else{
    expVec <- as.vector( exp( xBeta[ , yCat ])/( rowSums( exp( xBeta ) ) ) )
  }
  # weights
  weights <- elaIntWeights( shareVec )
  # calculation of the semi-elasticity
  semEla <- urbinElaInt( expVec, shareVec, 1:nInt, xBound, model = "lpm" )
  ### calculation of its standard error
  # partial derivatives of each semi-elasticity around each boundary
  # w.r.t. all estimated coefficients
  if( method == "binary" ){
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
  } else if( method == "MNL" ){
    gradM <- array( 0, c( nCoef, nInt - 1, pCoef ) )
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
  } else{
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
  }
  # partial derivative of the semi-elasticity 
  # w.r.t. all estimated coefficients
  derivCoef <- rep( 0, length( allCoef ) )
  for( m in 1:( nInt - 1 ) ){
    derivCoef <- derivCoef + weights[m] * gradM[,m]
  }
  # variance-covariance matrix of the coefficiencts
  vcovCoef <- diag( allCoefSE^2 )
  # standard error of the (average) semi-elasticity
  semElaSE <- drop( sqrt( t( derivCoef ) %*% vcovCoef %*% derivCoef ) )
  # prepare object that will be returned
  result <- c( semEla[1], stdEr = semElaSE )
  return( result )
}
