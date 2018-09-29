prepareEffInt <- function( allCoef, allXVal, xPos, refBound, intBound ){
  # list for 'results' to be returned 
  result <- list()
  # calculate xBars
  result$intX <- mean( intBound )
  result$refX <- mean( refBound ) 
  if( length( xPos ) == 2 ) {
    result$intX <- c( result$intX, EXSquared( intBound[1], intBound[2] ) )
    result$refX <- c( result$refX, EXSquared( refBound[1], refBound[2] ) )
  }
  if( length( result$intX ) != length( xPos ) || 
      length( result$refX ) != length( xPos ) ) {
    stop( "internal error: 'intX' or 'refX' does not have the expected length" )
  }
  # define X' * beta 
  result$intXbeta <- sum( allCoef * replace( allXVal, xPos, result$intX ) )
  result$refXbeta <- sum( allCoef * replace( allXVal, xPos, result$refX ) )
  checkXBeta( c( result$intXbeta, result$refXbeta ) )
  
  return( result )
}
