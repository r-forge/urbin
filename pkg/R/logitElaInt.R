logitElaInt <- function( allCoef, allXVal, xPos, xBound, yCat = NA,
  allCoefVcov = NULL, 
  model = "binary" ){
  # number of coefficients
  if( model == "binary" || model == "MNL" ){
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
  # Check position vector
  checkXPos( xPos, minLength = 2, maxLength = nCoef, 
    minVal = 0, maxVal = nCoef, requiredVal = 0 )
  # number of intervals
  nInt <- length( xPos ) 
  # check 'xBound' and replace infinite values
  xBound <- elaIntBounds( xBound, nInt )
  # check and prepare allCoefVcov
  allCoefVcov <- prepareVcov( allCoefVcov, length( allCoef ), xPos, xMeanSd = NULL )
  # vector of probabilities of y=1 for each interval and
  # vector of shares of observations in each interval
  xBeta <- matrix( rep( rep( NA, nInt ), pCoef ), ncol = pCoef ) 
  if( model == "binary" || model == "MNL" ){
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
  if( model == "binary" ){
    expVec <- as.vector( exp( xBeta )/( 1 + exp( xBeta ) ) )  
  } else if( model == "MNL" ){
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
  if( model == "binary" ){
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
  # if argument allXVal has attribute 'derivOnly',
  # return partial derivatives only (for testing partial derivatives)
  if( "derivOnly" %in% names( attributes( allXVal ) ) ) {
    return( derivCoef )
  }
  # approximate standard error of the semi-elasticity
  semElaSE <- drop( sqrt( t( derivCoef ) %*% allCoefVcov %*% derivCoef ) )
  # create object that will be returned
  result <- c( semEla[1], stdEr = semElaSE )
  return( result )
}
