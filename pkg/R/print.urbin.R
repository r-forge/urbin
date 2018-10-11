print.urbin <- function( x, ... ) {
  if( "semEla" %in% names( x ) ) {
    printVec <- c( x$semEla, x$stdEr )
  } else if( "effect" %in% names( x ) ) {
    printVec <- c( x$effect, x$stdEr )
  } else {
    stop( "internal error: object of class 'urbin' does not include",
      " expected components" )
  }
  print( unlist( x ) )
  invisible( x )
}
