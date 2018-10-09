urbinEffInt <- function( allCoef, allXVal = NULL, xPos, refBound, intBound, model,
  allCoefVcov = NULL, xMeanSd = NULL, yCat = NULL, 
  allCoefBra = NULL, allXValBra = NULL, yCatBra = NULL, lambda = NULL ){
  
  # number of coefficients
  nCoef <- length( allCoef )
  # number of explanatory variables
  nXVal <- length( allXVal )
  # check allXVal and allCoef
  if( model %in% c( "lpm", "probit", "logit" ) ){
    if( model == "lpm" ) {
      if( any( !is.na( allXVal ) ) ) {
        warning( "argument allXVal is ignored for lpm models",
          " (set this argument to 'NULL' or 'NA' to avoid this warning)" )
      }
      temp <- rep( 0, nCoef )
      temp[ xPos ] <- NA
      if( "derivOnly" %in% names( attributes( allXVal ) ) ) {
        attr( temp, "derivOnly" ) <-1
      }
      allXVal <- temp
      nXVal <- length( allXVal )
    }  
    if( nXVal != nCoef ){
      stop( "arguments 'allCoef' and 'allXVal' must have the same length" )
    }  
    # if( any( !is.na( allXVal[ xPos ] ) ) ) {
    #   allXVal[ xPos ] <- NA
    #   warning( "values of argument 'allXVal[ xPos ]' are ignored",
    #     " (set these values to 'NA' to avoid this warning)" )
    # }
  } else if( model == "MNL" ){
    # number of alternative categories of the dependent variable
    nYCat <- round( nCoef / nXVal )
    if( nCoef != nXVal * nYCat ) {
      stop( "length of argument 'allCoef' must be a multiple",
        " of the length of argument 'allXVal'" )
    } 
    # check argument yCat
    checkYCat( yCat, nYCat ) 
    # create matrix of coefficients
    mCoef <- matrix( allCoef, nrow = nXVal, ncol = nYCat )
  } else if( model == "CondL"){
    # number of categories of the dependent variable
    nYCat <- round( nXVal / nCoef )
    if( nXVal != nCoef * nYCat ) {
      stop( "length of argument 'allXVal' must be a multiple",
        " of the length of argument 'allCoef'" )
    } 
    # create matrix of explanatory variables
    mXVal <- matrix( allXVal, nrow = nCoef )
  } else if( model == "NestedL" ){
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
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  
  
  # Check position vector
  checkXPos( xPos, minLength = 1, maxLength = 2, minVal = 1, 
    maxVal = ifelse( model == "MNL", nXVal, nCoef ) )
  # check and prepare allCoefVcov
  allCoefVcov <- prepareVcov( allCoefVcov, nCoef, xPos, xMeanSd,
    nXVal = nXVal, pCall = match.call() )
  # check the boundaries of the intervals
  refBound <- elaIntBounds( refBound, 1, argName = "refBound" )
  intBound <- elaIntBounds( intBound, 1, argName = "intBound" )

  # calculate xBars
  intX <- mean( intBound )
  refX <- mean( refBound ) 
  if( length( xPos ) == 2 ) {
    intX <- c( intX, EXSquared( intBound[1], intBound[2] ) )
    refX <- c( refX, EXSquared( refBound[1], refBound[2] ) )
  }
  if( length( intX ) != length( xPos ) || 
      length( refX ) != length( xPos ) ) {
    stop( "internal error: 'intX' or 'refX' does not have the expected length" )
  }
  # define X' * beta 
  if( model %in% c( "lpm", "probit", "logit" ) ){
    intXbeta <- sum( allCoef * replace( allXVal, xPos, intX ) )
    refXbeta <- sum( allCoef * replace( allXVal, xPos, refX ) )
    if( model != "lpm" ) {
      checkXBeta( c( intXbeta, refXbeta ) )
    }
  } else if( model == "MNL" ){
    intXbeta <- replace( allXVal, xPos, intX ) %*% mCoef
    refXbeta <- replace( allXVal, xPos, refX ) %*% mCoef
    checkXBeta( c( intXbeta, refXbeta ) )
  } else if( model == "CondL" ){
    mXValint <- mXValref <- mXVal
    for( p in 1:nYCat ){
      mXValint[ ,p] <- replace( mXValint[ ,p], xPos, intX )
      mXValref[ ,p] <- replace( mXValref[ ,p], xPos, refX )
    }
    intXbeta <- drop( allCoef %*% mXValint )
    refXbeta <- drop( allCoef %*% mXValref )
    checkXBeta( c( intXbeta, refXbeta ) )
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
      intXbeta[[ l ]] <- drop( mCoef[ ,l ] %*% mXValint[[ l ]] )
      refXbeta[[ l ]] <- drop( mCoef[ ,l ] %*% mXValref[[ l ]] )
    }
    XbetaBra <- allCoefBra %*% mXValBra
    checkXBeta( c( unlist(refXbeta), unlist(intXbeta), XbetaBra ) )
  } else {
    stop( "argument 'model' specifies an unknown type of model" )
  }
  
  # calculate the effect
  if( model == "lpm" ) {
    eff <- intXbeta - refXbeta
  } else if( model == "probit" ) {
    eff <- pnorm( intXbeta ) - pnorm( refXbeta )
  } else   if( model == "logit" ){  
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

  # partial derivative of the effect w.r.t. all estimated coefficients
  if( model == "lpm" ) {
    derivCoef <- rep( 0, nCoef ) 
    derivCoef[ xPos ] <- intX - refX
  } else if( model == "probit" ) {
    derivCoef <- rep( NA, nCoef )
    derivCoef[ -xPos ] = ( dnorm( intXbeta ) - dnorm( refXbeta ) ) * 
      allXVal[ -xPos ] 
    derivCoef[ xPos ] = dnorm( intXbeta ) * intX - 
      dnorm( refXbeta ) * refX
  } else if( model == "logit" ){
    derivCoef <- rep( NA, nCoef )
    derivCoef[ -xPos ] <- ( exp( intXbeta )/( 1 + exp( intXbeta ) )^2 - 
        exp( refXbeta )/( 1 + exp( refXbeta ) )^2 ) * 
      allXVal[ -xPos ] 
    derivCoef[ xPos ] <- exp( intXbeta )/( 1 + exp( intXbeta ) )^2 * intX - 
      exp( refXbeta )/( 1 + exp( refXbeta ) )^2 * refX
  } else if( model == "MNL" ){
    derivCoef <- matrix( NA, nrow = nXVal, ncol = nYCat )
    for( p in 1:nYCat ){
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
  
  # if argument allXVal has attribute 'derivOnly',
  # return partial derivatives only (for testing partial derivatives)
  if( "derivOnly" %in% names( attributes( allXVal ) ) ) {
    return( derivCoef )
  }
  
  # approximate standard error of the effect
  effSE <- drop( sqrt( t( derivCoef ) %*% allCoefVcov %*% derivCoef ) )
  
  # object to be returned
  result <- c( effect = eff, stdEr = effSE )
  return( result )
}
