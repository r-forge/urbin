logitEffInt <- function( allCoef, allXVal = NA, xPos, refBound, intBound, 
  model = "binary", allCoefBra = NULL, allXValBra = NULL, 
  yCat = NULL, yCatBra = NULL, lambda = NULL, 
  allCoefSE = rep( NA, length( allCoef ) ) ){
  
  # number of coefficients
  nCoef <- length( allCoef )
  # number of explanatory variables
  nXVal <- length( allXVal )
  if( model == "binary" ){  
    # check arguments
    if( length( allXVal ) != nCoef ){
      stop( "argument 'allCoef' and 'allXVal' must have the same length" )
    }  
    if( length( allCoefSE ) != nCoef ){
      stop( "argument 'allCoef' and 'allCoefSE' must have the same length" )
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
    if( length( allCoefSE ) != nCoef ){
      stop( "argument 'allCoef' and 'allCoefSE' must have the same length" )
    }
  } else if( model == "CondL"){
    # number of coefficients
    nCoef <- length( allCoef )
    mXVal <- matrix( allXVal, nrow = nCoef )
    pCoef <- dim( mXVal )[2]
    # check arguments
    if( dim( mXVal )[1] != nCoef ){
      stop( "argument 'allCoef' and 'allXVal' must have the same length" )
    }
    if( length( allCoefSE ) != nCoef ){
      stop( "argument 'allCoef' and 'allCoefSE' must have the same length" )
    }  
  } else if( model == "NestedL" ){
    nCoef <- length( allCoef )
    NCoef <- length( allCoefBra )
    mXValBra <- matrix( allXValBra, nrow = NCoef )
    nXValBra <- dim( mXValBra )[1]
    pXValBra <- dim( mXValBra )[2]
    # check arguments
    if( NCoef != nXValBra ){
      stop( "arguments 'allCoefBra' and 'allXValBra' must have the same length")
    }
    O <- length( allXVal )
    nXVal <- unlist( lapply( allXVal, function(x) dim( x )[1] ) )
    pCoef <- unlist( lapply( allXVal, function(x) dim( x )[2] ) )
    if( nCoef != nXVal[ yCatBra ] ){
      stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
    }
    if( nCoef != length( allCoefSE) ){
      stop( "arguments 'allCoef' and 'allCoefSE' must have the same length" )
    }
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  checkXPos( xPos, minLength = 1, maxLength = 2, minVal = 1, maxVal = nCoef )
  refBound <- elaIntBounds( refBound, 1, argName = "refBound" )
  intBound <- elaIntBounds( intBound, 1, argName = "intBound" )
  if( model == "binary" || model == "MNL" ){  
    if( any( !is.na( allXVal[ xPos ] ) ) ) {
      allXVal[ xPos ] <- NA
      warning( "values of argument 'allXVal[ xPos ]' are ignored",
        " (set these values to 'NA' to avoid this warning)" )
    }
  } else if( model == "CondL" ){
    for( p in 1:pCoef ){
      if( any( !is.na( mXVal[ xPos, p ] ) ) ){
        mXVal[ xPos, p ] <- NA
        warning( "values of argument 'allXVal[ xPos ]' are ignored",
          " (set these values to 'NA' to avoid this warning)" )
      }
    }
  } else if( model == "NestedL" ){
    for( p in 1:pCoef[ yCatBra ] ){
      if( any( !is.na( allXVal[[ yCatBra ]][ xPos, p ] ) ) ){
        mXVal[ xPos, p ] <- NA
        warning( "values of argument 'allXVal[ xPos ]' are ignored",
          " (set these values to 'NA' to avoid this warning)" )
      }
    }
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  } 
  # calculate xBars
  intX <- mean( intBound )
  refX <- mean( refBound ) 
  if( length( xPos ) == 2 ) {
    intX <- c( intX, EXSquared( intBound[1], intBound[2] ) )
    refX <- c( refX, EXSquared( refBound[1], refBound[2] ) )
  }
  if( length( intX ) != length( xPos ) || length( refX ) != length( xPos ) ) {
    stop( "internal error: 'intX' or 'refX' does not have the expected length" )
  }
  # define X' * beta 
  if( model == "binary" ){
    intXbeta <- sum( allCoef * replace( allXVal, xPos, intX ) )
    refXbeta <- sum( allCoef * replace( allXVal, xPos, refX ) )
    checkXBeta( c( intXbeta, refXbeta ) )
  } else if( model == "MNL" ){
    intXbeta <- colSums( mCoef * replace( allXVal, xPos, intX ) )
    refXbeta <- colSums( mCoef * replace( allXVal, xPos, refX ) )
  } else if( model == "CondL" ){
    mXValint <- mXValref <- mXVal
    for( p in 1:pCoef ){
      mXValint[ ,p] <- replace( mXValint[ ,p], xPos, intX )
      mXValref[ ,p] <- replace( mXValref[ ,p], xPos, refX )
    }
    intXbeta <- colSums( allCoef * mXValint )
    refXbeta <- colSums( allCoef * mXValref )
  } else if( model == "NestedL" ){
    mCoef <- matrix( rep( allCoef, O ), nrow = nCoef, O ) %*% diag( 1/ lambda )
    mXValint <- mXValref <- allXVal
    for( i in 1:O ){
      for( p in 1:pCoef[i] ){
        mXValint[[i]][ ,p] <- replace( mXValint[[i]][ ,p], xPos, intX )
        mXValref[[i]][ ,p] <- replace( mXValref[[i]][ ,p], xPos, refX )
      }
    }  
    refXbeta <- intXbeta <- rep( list( NA ), O )
    for( l in 1:O ){  
      intXbeta[[ l ]] <- colSums( mCoef[ ,l ] * mXValint[[ l ]] )
      refXbeta[[ l ]] <- colSums( mCoef[ ,l ] * mXValref[[ l ]] )
    }
    XbetaBra <- colSums( allCoefBra * mXValBra )
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  # effect E_{k,ml}
  if( model == "binary" ){  
    eff <- exp( intXbeta )/( 1 + exp( intXbeta ) ) - 
      exp( refXbeta )/( 1 + exp( refXbeta ) )
  } else if( model == "MNL" ){
    eff <- exp( intXbeta[ yCat ] )/( 1 + sum( exp( intXbeta ) ) ) - 
      exp( refXbeta[ yCat ] )/( 1 + sum( exp( refXbeta ) ) )
  } else if( model == "CondL"){
    eff <- exp( intXbeta[ yCat ] )/( sum( exp( intXbeta ) ) ) -
      exp( refXbeta[ yCat ] )/( sum( exp( refXbeta ) ) )    
  } else if( model == "NestedL" ){
    intBranch <- refBranch <- rep( list( NA ), O )
    for( l in 1:O ){
      intBranch[[ l ]] <- exp( XbetaBra[ l ] + lambda[ l ] * 
          log( sum( exp( intXbeta[[ l ]] ) ) ) ) 
      refBranch[[ l ]] <- exp( XbetaBra[ l ] + lambda[ l ] * 
          log( sum( exp( refXbeta[[ l ]] ) ) ) )
    }
    intBranch <- unlist( intBranch )
    refBranch <- unlist( refBranch )
    eff <- exp( intXbeta[[ yCatBra ]][ yCat ] )/( sum( exp( intXbeta[[ yCatBra ]] ) ) ) *
      intBranch[ yCatBra ]/ sum( intBranch ) - 
      exp( refXbeta[[ yCatBra ]][ yCat ] )/( sum( exp( refXbeta[[ yCatBra ]] ) ) ) *
      refBranch[ yCatBra ]/ sum( refBranch )
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  # calculating approximate standard error
  # partial derivative of E_{k,ml} w.r.t. all estimated coefficients
  if( model == "binary" ){
    derivCoef <- rep( NA, nCoef )
    derivCoef[ -xPos ] <- ( exp( intXbeta )/( 1 + exp( intXbeta ) )^2 - 
        exp( refXbeta )/( 1 + exp( refXbeta ) )^2 ) * 
      allXVal[ -xPos ] 
    derivCoef[ xPos ] <- exp( intXbeta )/( 1 + exp( intXbeta ) )^2 * intX - 
      exp( refXbeta )/( 1 + exp( refXbeta ) )^2 * refX
  } else if( model == "MNL" ){
    derivCoef <- matrix( NA, nrow = nXVal, ncol = pCoef )
    for( p in 1:pCoef ){
      if( p == yCat ){
        derivCoef[ -xPos, p ] <- 
          ( exp( intXbeta[ p ] ) * 
              ( 1 + sum( exp( intXbeta[ -yCat ] ) ) )/
              ( 1 + sum( exp( intXbeta ) ) )^2 -
              exp( refXbeta[ p ] ) * 
              ( 1 + sum( exp( refXbeta[ -yCat ] ) ) )/
              ( 1 + sum( exp( refXbeta ) ) )^2 ) * allXVal[ - xPos ]
        derivCoef[ xPos, p ] <- 
          ( exp( intXbeta[ p ] ) * 
              ( 1 + sum( exp( intXbeta[ -yCat ] ) ) )/
              ( 1 + sum( exp( intXbeta ) ) )^2 ) * intX -
          ( exp( refXbeta[ p ] ) * 
              ( 1 + sum( exp( refXbeta[ -yCat ] ) ) )/
              ( 1 + sum( exp( refXbeta ) ) )^2 ) * refX
      } else{  
        derivCoef[ -xPos, p ] <- 
          ( ( exp( refXbeta[ yCat ] ) * exp( refXbeta[ p ] ) )/
              ( 1 + sum( exp( refXbeta ) ) )^2 -
              ( exp( intXbeta[ yCat ] ) * exp( intXbeta[ p ] ) )/
              ( 1 + sum( exp( intXbeta ) ) )^2 ) * allXVal[ -xPos ]
        derivCoef[ xPos, p ] <- 
          ( ( exp( refXbeta[ yCat ] ) * exp( refXbeta[ p ] ) )/
              ( 1 + sum( exp( refXbeta ) ) )^2 ) * intX -
          ( ( exp( intXbeta[ yCat ] ) * exp( intXbeta[ p ] ) )/
              ( 1 + sum( exp( intXbeta ) ) )^2 ) * refX
      }     
    }
    derivCoef <- c( derivCoef )
  } else if( model == "CondL" ){
    derivCoef <- rep( NA, nCoef )
    derivCoef[ -xPos ] <- ( exp( intXbeta[ yCat] ) * mXVal[ -xPos, yCat] * 
        sum( exp( intXbeta ) ) -
        exp( intXbeta[ yCat] ) * rowSums( exp( intXbeta ) * 
            mXVal[ -xPos, ] ) )/
      ( sum( exp( intXbeta ) ) )^2 - 
      ( exp( refXbeta[ yCat] ) * mXVal[ -xPos, yCat] * 
          sum( exp( refXbeta ) ) -
          exp( refXbeta[ yCat] ) * rowSums( exp( refXbeta ) * 
              mXVal[ -xPos, ] ) )/
      ( sum( exp( refXbeta ) ) )^2 
    derivCoef[ xPos ] <-  ( exp( intXbeta[ yCat] ) * intX * 
        sum( exp( intXbeta ) ) -
        exp( intXbeta[ yCat] ) * sum( exp( intXbeta ) * intX ) )/
      ( sum( exp( intXbeta ) ) )^2 - 
      ( exp( refXbeta[ yCat] ) * refX * 
          sum( exp( refXbeta ) ) -
          exp( refXbeta[ yCat] ) * sum( exp( refXbeta ) * refX ) )/
      ( sum( exp( refXbeta ) ) )^2 
  } else if( model == "NestedL" ){
    derivCoef <- rep( NA, nCoef ) 
    PImp <- exp( intXbeta[[ yCatBra ]][ yCat ])/( sum( exp( intXbeta[[ yCatBra ]] ) ) )
    PIlp <- exp( refXbeta[[ yCatBra ]][ yCat ])/( sum( exp( refXbeta[[ yCatBra ]] ) ) )
    PImo <- intBranch[ yCatBra ]/ sum( intBranch )
    PIlo <- refBranch[ yCatBra ]/ sum( refBranch )
    Om <- matrix( 
      unlist( lapply( allXVal, function(x) rowSums( x[ -xPos, , drop = FALSE ] ) ) ), 
      ncol = O ) 
    derivCoef[ -xPos ] <- ( ( allXVal[[ yCatBra ]][ -xPos, yCat ]/lambda[ yCatBra ] -
        ( rowSums( 
          ( allXVal[[ yCatBra ]][ -xPos, ]/lambda[ yCatBra ] ) %*%
            diag( exp( intXbeta[[ yCatBra ]] ) ) ) )/
        ( sum( exp( intXbeta[[ yCatBra ]] ) ) ) ) + 
        ( rowSums( allXVal[[ yCatBra ]][ -xPos, ] ) -
            ( rowSums( Om %*% diag( exp( intBranch ) ) )/
                ( sum( intBranch ) ) ) ) ) * PImp * PImo -
      ( ( allXVal[[ yCatBra ]][ -xPos, yCat ]/lambda[ yCatBra ] -
          ( rowSums( 
            ( allXVal[[ yCatBra ]][ -xPos, ]/lambda[ yCatBra ] ) %*%
              diag( exp( refXbeta[[ yCatBra ]] ) ) ) )/
          ( sum( exp( refXbeta[[ yCatBra ]] ) ) ) ) + 
          ( rowSums( allXVal[[ yCatBra ]][ -xPos, ] ) -
              ( rowSums( Om %*% diag( exp( refBranch ) ) )/
                  ( sum( refBranch ) ) ) ) ) * PIlp * PIlo
    derivCoef[ xPos ] <-  ( ( intX/lambda[ yCatBra ] -
        ( sum( intX/lambda[ yCatBra ]  *
            exp( intXbeta[[ yCatBra ]] ) ) )/
        ( sum( exp( intXbeta[[ yCatBra ]] ) ) ) ) + 
        ( intX * pCoef[ yCatBra ] -
            ( sum( intX * exp( intBranch ) )/
                ( sum( intBranch ) ) ) ) ) * PImp * PImo -
      ( ( refX/lambda[ yCatBra ] -
          ( sum( refX/lambda[ yCatBra ]  *
              exp( refXbeta[[ yCatBra ]] ) ) )/
          ( sum( exp( refXbeta[[ yCatBra ]] ) ) ) ) + 
          ( refX * pCoef[ yCatBra ] -
              ( sum( refX * exp( refBranch ) )/
                  ( sum( refBranch ) ) ) ) ) * PImp * PImo  
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  # variance covariance of the coefficients (covariances set to zero)
  vcovCoef <- diag( allCoefSE^2 )
  # approximate standard error of the effect
  effSE <- drop( sqrt( t( derivCoef ) %*% vcovCoef %*% derivCoef ) )
  # object to be returned
  result <- c( effect = eff, stdEr = effSE )
  return( result )
}
