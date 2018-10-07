library( "urbin" )

# Example
eff10a <- urbin:::logitEffCat( allCoef = c( 0.28, 0.003, 0.0075, -0.034, 
  -0.005, 0.89, -1.2 ), 
  allXVal = c( 1, 0.1, 0.3, 0.15, 0.2, 2.34, 10.8 ), model = "logit", 
  xPos = c( 2:5 ), xGroups = c( 0, -1, 1, 1, -1 ) )
eff10a

eff10b <- urbin:::logitEffCat( allCoef = c( 0.28, 0.003, 0.0075, -0.034, 
  -0.005, 0.89, -1.2 ), 
  allXVal = c( 1, 0.1, 0.3, 0.15, 0.2, 2.34, 10.8 ), 
  xPos = c( 2:5 ), xGroups = c( 0, -1, 1, 1, -1 ), model = "logit",
  allCoefVcov = c( 0.03, 0.0001, 0.005, 0.01, 
    0.004, 0.05, 0.8 ) )
eff10b

# Example
eff11a <- urbin:::logitEffCat( allCoef = c( 0.28, 0.003, 0.0075, -0.034, 
  -0.005, 0.89, 0.005, 0.06, 1.7 ),
  allXVal = c( 1, 0.5, 0.3 ), xPos = c( 2:3 ), 
  xGroups = c( -1, -1, 1 ), yCat = 2, model = "MNL" )
eff11a

eff11b <- urbin:::logitEffCat( allCoef = c( 0.28, 0.003, 0.0075, -0.034, 
  -0.005, 0.89, 0.005, 0.06, 1.7 ),
  allXVal = c( 1, 0.5, 0.3 ), xPos = c( 2:3 ), 
  xGroups = c( -1, -1, 1 ), yCat = 2, model = "MNL", 
  allCoefVcov = c( 0.03, 0.0001, 0.005, 0.01, 0.004, 
    0.05, 0.004, 0.5, 0.0078 ) )
eff11b

# Example
eff12a <- urbin:::logitEffCat( allCoef = c( 0.28, 0.003, 0.0075 ),
  allXVal = c( 1, 0.5, 0.3, 1, 0.4, 0.4 ), 
  xPos = c( 2:3 ), 
  xGroups = c( -1, -1, 1 ), yCat = 2, model = "CondL" )
eff12a

eff12b <- urbin:::logitEffCat( allCoef = c( 0.28, 0.003, 0.0075 ),
  allXVal = c( 1, 0.5, 0.3, 1, 0.4, 0.4 ), 
  xPos = c( 2:3 ), 
  allCoefVcov = c( 0.03, 0.0001, 0.005 ),
  xGroups = c( -1, -1, 1 ), yCat = 2, model = "CondL" )
eff12b
