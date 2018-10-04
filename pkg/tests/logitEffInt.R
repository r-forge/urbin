library( "urbin" )

# Example
eff6a <- urbin:::logitEffInt( allCoef = c( 0.33, 0.22, 0.05, 0.6 ),
  allXVal = c( 1, NA, 0.16, 0.13 ),
  xPos = 2, 
  refBound = c( 8, 12 ), intBound = c( 13, 15 ) )
eff6a

eff6b <- urbin:::logitEffInt( allCoef = c( 0.33, 0.22, 0.05, 0.6 ),
  allXVal = c( 1, NA, 0.16, 0.13 ),
  xPos = 2, 
  refBound = c( 8, 12 ), intBound = c( 13, 15 ), 
  allCoefSE = c( 0.003, 0.045, 0.007, 0.009 ) )
eff6b

# Example
eff7a <- urbin:::logitEffInt( allCoef = c( 0.33, 0.022, 0.005, 0.6 ),
  allXVal = c( 1, NA, NA, 0.0004 ),
  xPos = c( 2, 3 ), 
  refBound = c( 8, 12 ), intBound = c( 13, 15 ))
eff7a

eff7b <- urbin:::logitEffInt( allCoef = c( 0.33, 0.022, 0.005, 0.6 ),
  allXVal = c( 1, NA, NA, 0.13 ),
  xPos = c( 2, 3 ), 
  refBound = c( 8, 12 ), intBound = c( 13, 15 ), 
  allCoefSE = c( 0.003, 0.011, 0.0025, 0.009 ) )
eff7b

#Example
eff8a <- urbin:::logitEffInt( allCoef = c( 0.2, 0.3, 0.5, -0.2, 0.03, 0.6 ), 
  allXVal = c( 1, NA, 0.12 ), 
  xPos = 2, 
  refBound = c( 8, 12 ), intBound = c( 13, 15 ),
  yCat = 2, model = "MNL" )
eff8a

eff8b <- urbin:::logitEffInt( allCoef = c( 0.2, 0.3, 0.5, -0.2, 0.03, 0.6 ), 
  allXVal = c( 1, NA, 0.12 ), 
  xPos = 2, 
  refBound = c( 8, 12 ), intBound = c( 13, 15 ),
  yCat = 2, 
  allCoefSE = c( 0.003, 0.045, 0.007, 0.009, 0.0008, 0.9 ),
  model = "MNL" )
eff8b

#Example
eff9a <- urbin:::logitEffInt( allCoef = c( 0.2, 0.3, 0.5, -0.2, 0.03, 0.6 ), 
  allXVal = c( 1, NA, NA ), 
  xPos = c( 2, 3 ), 
  refBound = c( 1, 12 ), intBound = c( 13, 25 ),
  yCat = 2, model = "MNL" )
eff9a

eff9b <- urbin:::logitEffInt( allCoef = c( 0.2, 0.3, 0.5, -0.2, 0.03, 0.6 ), 
  allXVal = c( 1, NA, NA ), 
  xPos = c( 2, 3 ), 
  refBound = c( 1, 12 ), intBound = c( 13, 25 ),
  yCat = 2, 
  allCoefSE = c( 0.003, 0.045, 0.007, 0.009, 0.0008, 0.009 ),
  model = "MNL" )
eff9b

#Example
eff10a <- urbin:::logitEffInt( allCoef = c( 0.2, 0.3, 0.5, 0.091 ), 
  allXVal = c( 1, NA, NA, 2.45, 1, NA, NA, 0.79 ), 
  xPos = c( 2, 3 ), 
  refBound = c( 1, 12 ), intBound = c( 13, 15 ),
  yCat = 2, model = "CondL" )
eff10a

eff10b <- urbin:::logitEffInt( allCoef = c( 0.2, 0.3, 0.5, 0.091 ), 
  allXVal = c( 1, NA, NA, 2.45, 1, NA, NA, 0.79 ), 
  xPos = c( 2, 3 ), 
  refBound = c( 8, 12 ), intBound = c( 13, 15 ),
  allCoefSE = c( 0.003, 0.045, 0.007, 0.009 ),
  yCat = 2, model = "CondL" )
eff10b

# Example 
matrix1 <- matrix( c( 1, NA, 0.3, 0.09, 1, NA, 0.9, 1.8 ), nrow = 4 )
matrix2 <- matrix( c( 1, NA, 0.099, 0.211 ), nrow = 4 )
eff12a <- urbin:::logitEffInt( allCoefBra = c( 0.445, 0.03, -0.002 ), 
  allCoef = c( 1.8, 0.005, -0.12, 0.8 ), 
  allXValBra = c( 1, 3.3, 4.5, 1, 0.1, 0.987 ), 
  allXVal = list( matrix1, matrix2 ), 
  refBound = c( 0.5, 1.5 ), intBound = c( 2, 3.5 ),
  xPos = 2, yCatBra = 1, yCat = 2, lambda = c( 0.8, 1 ), 
  model = "NestedL" )
eff12a

matrix1 <- matrix( c( 1, NA, 0.3, 0.09, 1, NA, 0.9, 1.8 ), nrow = 4 )
matrix2 <- matrix( c( 1, NA, 0.099, 0.211 ), nrow = 4 )
eff12b <- urbin:::logitEffInt( allCoefBra = c( 0.445, 0.03, -0.002 ), 
  allCoef = c( 1.8, 0.005, -0.12, 0.8 ), 
  allXValBra = c( 1, 3.3, 4.5, 1, 0.1, 0.987 ), 
  allXVal = list( matrix1, matrix2 ), 
  allCoefSE = c( 0.003, 0.045, 0.007, 0.0032 ),
  refBound = c( 0.5, 1.5 ), intBound = c( 2, 3.5 ),
  xPos = 2, yCatBra = 1, yCat = 2, lambda = c( 0.8, 1 ), 
  model = "NestedL" )
eff12b
