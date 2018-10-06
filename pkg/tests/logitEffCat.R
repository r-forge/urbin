library( "urbin" )

# Example
eff10a <- urbin:::logitEffCat( allCoef = c( 0.28, 0.003, 0.0075, 0, -0.034, 
  -0.005, 0.89, -1.2 ), 
  allXVal = c( 1, 0.1, 0.3, 0.25, 0.15, 0.2, 2.34, 10.8 ), model = "logit", 
  xPos = c( 2:6 ), xGroups = c( 0, -1, -1, 1, 1 ) )
eff10a

eff10b <- urbin:::logitEffCat( allCoef = c( 0.28, 0.003, 0.0075, 0, -0.034, 
  -0.005, 0.89, -1.2 ), 
  allXVal = c( 1, 0.1, 0.3, 0.25, 0.15, 0.2, 2.34, 10.8 ), 
  xPos = c( 2:6 ), xGroups = c( 0, -1, -1, 1, 1 ), model = "logit",
  allCoefSE = c( 0.03, 0.0001, 0.005, 0, 0.01, 
    0.004, 0.05, 0.8 ) )
eff10b

# Example
eff11a <- urbin:::logitEffCat( allCoef = c( 0.28, 0.003, 0.0075, 0, -0.034, 
  -0.005, 0.89, 0, 0.005, 0.06, 1.7, 0 ),
  allXVal = c( 1, 0.5, 0.3, 0.2 ), xPos = c( 2:4 ), 
  xGroups = c( -1, -1, 1 ), yCat = 2, model = "MNL" )
eff11a

eff11b <- urbin:::logitEffCat( allCoef = c( 0.28, 0.003, 0.0075, 0, -0.034, 
  -0.005, 0.89, 0, 0.005, 0.06, 1.7, 0 ),
  allXVal = c( 1, 0.5, 0.3, 0.2 ), xPos = c( 2:4 ), 
  xGroups = c( -1, -1, 1 ), yCat = 2, model = "MNL", 
  allCoefSE = c( 0.03, 0.0001, 0.005, 0, 0.01, 0.004, 
    0.05, 0, 0.004, 0.5, 0.0078, 0 ) )
eff11b

# Example
eff12a <- urbin:::logitEffCat( allCoef = c( 0.28, 0.003, 0.0075, 0 ),
  allXVal = c( 1, 0.5, 0.3, 0.2, 1, 0.4, 0.4, 0.1 ), 
  xPos = c( 2:4 ), 
  xGroups = c( -1, -1, 1 ), yCat = 2, model = "CondL" )
eff12a

eff12b <- urbin:::logitEffCat( allCoef = c( 0.28, 0.003, 0.0075, 0 ),
  allXVal = c( 1, 0.5, 0.3, 0.2, 1, 0.4, 0.4, 0.1 ), 
  xPos = c( 2:4 ), 
  allCoefSE = c( 0.03, 0.0001, 0.005, 0 ),
  xGroups = c( -1, -1, 1 ), yCat = 2, model = "CondL" )
eff12b
