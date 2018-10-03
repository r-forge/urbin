library( "urbin" )

# Example
ela8a <- logitElaInt( allCoef = c( 0.33, 0.22, 0.05, 0.6 ), 
  allXVal = c( 1, 0.4, 0.12, 0.13 ), 
  xPos = c( 2, 0, 3, 4 ),
  xBound = c( 0, 500, 1000, 1500, Inf ) )
ela8a 

ela8b <- logitElaInt( allCoef = c( 0.33, 0.22, 0.05, 0.6 ), 
  allXVal = c( 1, 0.4, 0.12, 0.13 ), 
  xPos = c( 2, 0, 3, 4 ),
  xBound = c( 0, 500, 1000, 1500, Inf ),
  allCoefSE = c( 0.003, 0.045, 0.007, 0.009 ) )
ela8b 

# Example
ela9a <- logitElaInt( allCoef = c( 0.2, 0.3, 0.5, -0.2, 0.03, 0.6 ), 
  allXVal = c( 1, 0.4, 0.12 ), 
  xPos = c( 2, 0, 3 ), 
  xBound = c( 0, 500, 1000, Inf ), yCat = 2, 
  method = "MNL" )
ela9a

ela9b <- logitElaInt( allCoef = c( 0.2, 0.3, 0.5, -0.2, 0.03, 0.6 ), 
  allXVal = c( 1, 0.4, 0.12 ), 
  xPos = c( 2, 0, 3 ), 
  xBound = c( 0, 500, 1000, Inf ), yCat = 2, 
  allCoefSE = c( 0.003, 0.045, 0.007, 0.009, 0.0008, 0.9 ),
  method = "MNL" )
ela9b

# Example
ela10a <- logitElaInt( allCoef = c( 1.33, 0.022, 0.58, 1.6 ), 
  allXVal = c( 1, 0.4, 0.12, 0.0002, 
    1, 0.28, 0.01, 0.000013 ), 
  xPos = c( 2, 0, 3 ),
  xBound = c( 0, 1000, 1500, Inf ), yCat = 2,
  method = "CondL" )
ela10a 

ela10b <- logitElaInt( allCoef = c( 1.33, 0.022, 0.58, 1.6 ), 
  allXVal = c( 1, 0.4, 0.12, 0.0002, 
    1, 0.28, 0.01, 0.000013 ), 
  xPos = c( 2, 0, 3 ),
  xBound = c( 0, 1000, 1500, Inf ), 
  allCoefSE = c( 0.003, 0.045, 0.007, 0.009 ), yCat = 2,
  method = "CondL" )
ela10b 
