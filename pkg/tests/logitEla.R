library( "urbin" )

ela6a <- urbinEla( allCoef = c( 0.445, 0.03, 0.00002, 0.067, 0.89, 0.124 ),
  allXVal = c( 1, 3.3, 4.5, 2.34, 0.1, 0.987 ), xPos = 2 )
ela6a

ela6b <- urbinEla( allCoef = c( 0.445, 0.03, 0.00002, 0.067, 0.89, 0.124 ),
  allXVal = c( 1, 3.3, 4.5, 2.24, 0.1, 0.987 ), 
  allCoefVcov = c( 0.001, 0.02, 0.000002, 0.05, 1.2, 0.03 ), 
  xPos = 2 )
ela6b  

# Example 
ela7a <- urbinEla( allCoef = c( 0.445, 0.03, 0.00002, 0.067, 0.89, 0.124 ),
  allXVal = c( 1, 3.3, 3.3^2, 2.34, 0.1, 0.987 ), 
  xPos = c( 2, 3 ) )
ela7a

ela7b <- urbinEla( allCoef = c( 0.445, 0.03, 0.00002, 0.067, 0.89, 0.124 ),
  allXVal = c( 1, 3.3, 3.3^2, 2.34, 0.1, 0.987 ), 
  allCoefVcov = c( 0.001, 0.02, 0.000002, 0.05, 1.2, 0.03 ),
  xPos = c( 2, 3 ) )
ela7b

# Example
ela8a <- urbinEla( allCoef = c( 0.2, 0.3, 0.5, -0.2, 0.03, 0.6 ), 
  allXVal = c( 1, 8.4, 0.06 ), xPos = 3,
  method = "MNL", yCat = 2 )
ela8a

ela8b <- urbinEla( allCoef = c( 0.2, 0.3, 0.5, -0.2, 0.03, 0.6 ), 
  allXVal = c( 1, 8.4, 0.06 ), 
  allCoefVcov = c( 0.002, 0.003, 0.004, 0.006, 0.00001, 0.08 ), 
  xPos = 3, 
  method = "MNL", yCat = 2 )
ela8b

# Example
ela9a <- urbinEla( allCoef = c( 0.2, 0.3, 0.5, -0.2, 0.03, 0.6 ), 
  allXVal = c( 1, 0.04, 0.0016 ), xPos = c( 2, 3 ),
  method = "MNL", yCat = 2 )
ela9a

ela9b <- urbinEla( allCoef = c( 0.2, 0.3, 0.5, -0.2, 0.03, 0.6 ), 
  allXVal = c( 1, 0.04, 0.0016 ), 
  allCoefVcov = c( 0.002, 0.003, 0.004, 0.006, 0.00001, 0.08 ), 
  xPos = c( 2, 3 ), 
  method = "MNL", yCat = 2 )
ela9b

# Example
ela10a <- urbinEla( allCoef = c( 0.445, 0.03, 0.00002 ),
  allXVal = c( 1, 3.3, 4.5, 1, 0.1, 0.987 ),
  xPos = 2, 
  method = "CondL", yCat = 2 )
ela10a

ela10b <- urbinEla( allCoef = c( 0.445, 0.03, -0.002 ),
  allXVal = c( 1, 0.3, 0.09, 1, 0.1, 0.01 ),
  xPos = c( 2, 3 ), 
  method = "CondL", yCat = 2 )
ela10b

# Example
ela11a <- urbinEla( allCoef = c( 0.445, 0.03, 0.00002 ),
  allXVal = c( 1, 3.3, 4.5, 1, 0.1, 0.987 ),
  allCoefVcov = c( 0.002, 0.003, 0.004 ),
  xPos = 2, 
  method = "CondL", yCat = 2 )
ela11a

ela11b <- urbinEla( allCoef = c( 0.445, 0.03, -0.002 ),
  allXVal = c( 1, 0.3, 0.09, 1, 0.1, 0.01 ),
  allCoefVcov = c( 0.002, 0.003, 0.004 ),
  xPos = c( 2, 3 ), 
  method = "CondL", yCat = 2 )
ela11b

# Example 
matrix1 <- matrix( c( 1, 2.5, 0.3, 0.09, 1, 0.33, 0.9, 1.8 ), nrow = 4 )
matrix2 <- matrix( c( 1, 2.8, 0.099, 0.211 ), nrow = 4 )
ela12a <- urbinEla( allCoefBra = c( 0.445, 0.03, -0.002 ), 
  allCoef = c( 1.8, 0.005, -0.12, 0.8 ), 
  allXValBra = c( 1, 3.3, 4.5, 1, 0.1, 0.987 ), 
  allXVal = list( matrix1, matrix2 ), 
  xPos = 2, yCatBra = 1, yCat = 2, lambda = c( 0.8, 1 ), 
  method = "NestedL" )
ela12a

matrix1 <- matrix( c( 1, 0.3, 0.09, 0.09, 1, 0.33, 0.1089, 1.8 ), nrow = 4 )
matrix2 <- matrix( c( 1, 0.31, 0.099, 0.211 ), nrow = 4 )
ela12b <- urbinEla( allCoefBra = c( 0.445, 0.03, -0.002 ), 
  allCoef = c( 1.8, 0.005, -0.12, 0.8 ), 
  allXValBra = c( 1, 3.3, 4.5, 1, 0.1, 0.987 ), 
  allXVal = list( matrix1, matrix2 ), 
  xPos = 2, yCatBra = 1, yCat = 2, lambda = c( 0.8, 1 ), 
  method = "NestedL" )
ela12b

# Example 
matrix1 <- matrix( c( 1, 2.5, 0.3, 0.09, 1, 0.33, 0.9, 1.8 ), nrow = 4 )
matrix2 <- matrix( c( 1, 2.8, 0.099, 0.211 ), nrow = 4 )
ela13a <- urbinEla( allCoefBra = c( 0.445, 0.03, -0.002 ), 
  allCoef = c( 1.8, 0.005, -0.12, 0.8 ), 
  allXValBra = c( 1, 3.3, 4.5, 1, 0.1, 0.987 ), 
  allXVal = list( matrix1, matrix2 ), 
  allCoefVcov = c( 0.001, 0.089, 0.0003, 0.12 ),
  xPos = 2, yCatBra = 1, yCat = 2, lambda = c( 0.8, 1 ), 
  method = "NestedL" )
ela13a

matrix1 <- matrix( c( 1, 0.3, 0.09, 0.09, 1, 0.33, 0.1089, 1.8 ), nrow = 4 )
matrix2 <- matrix( c( 1, 0.31, 0.099, 0.211 ), nrow = 4 )
ela13b <- urbinEla( allCoefBra = c( 0.445, 0.03, -0.002 ), 
  allCoef = c( 1.8, 0.005, -0.12, 0.8 ), 
  allXValBra = c( 1, 3.3, 4.5, 1, 0.1, 0.987 ), 
  allXVal = list( matrix1, matrix2 ), 
  allCoefVcov = c( 0.001, 0.089, 0.0003, 0.12 ),
  xPos = 2, yCatBra = 1, yCat = 2, lambda = c( 0.8, 1 ), 
  method = "NestedL" )
ela13b
