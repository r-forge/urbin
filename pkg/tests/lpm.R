library( "urbin" )
library( "maxLik" )

# load data set
data( "Mroz87", package = "sampleSelection" )

# create dummy variable for kids
Mroz87$kids <- as.numeric( Mroz87$kids5 > 0 | Mroz87$kids618 > 0 )

### linear in age
estLpmLin <- lm( lfp ~ kids + age + educ, 
  data = Mroz87 )
summary( estLpmLin )
# mean values of the explanatory variables
xMeanLin <- c( 1, colMeans( Mroz87[ , c( "kids", "age", "educ" ) ] ) )
# semi-elasticity of age without standard errors
urbinEla( coef( estLpmLin )[ "age" ], xMeanLin[ "age" ], xPos = 1,
  model = "lpm" )
urbinEla( coef( estLpmLin ), xMeanLin, xPos = 3, model = "lpm" )
# semi-elasticity of age based on numerical derivation
100 * ( predict( estLpmLin, 
  newdata = as.data.frame( t( xMeanLin * c( 1, 1, 1.005, 1 ) ) ) ) -
    predict( estLpmLin, 
      newdata = as.data.frame( t( xMeanLin * c( 1, 1, 0.995, 1 ) ) ) ) )
# partial derivatives of the semi-elasticity wrt the coefficients
xMeanAgeAttr <- xMeanLin["age"]
attr( xMeanAgeAttr, "derivOnly" ) <- 1 
urbinEla( coef( estLpmLin )["age"], xMeanAgeAttr, xPos = 1, model = "lpm" )
xMeanLinAttr <- xMeanLin
attr( xMeanLinAttr, "derivOnly" ) <- 1 
urbinEla( coef( estLpmLin ), xMeanLinAttr, xPos = 3, model = "lpm" )
# numerically computed partial derivatives of the semi-elasticity wrt the coefficients
numericGradient( urbinEla, t0 = coef( estLpmLin )["age"], 
  allXVal = xMeanLin["age"], xPos = 1, model = "lpm" )
numericGradient( urbinEla, t0 = coef( estLpmLin ), 
  allXVal = xMeanLin, xPos = 3, model = "lpm" )
# semi-elasticity of age with standard errors (only standard errors)
urbinEla( coef( estLpmLin )["age"], xMeanLin["age"], xPos = 1, model = "lpm",
  sqrt( diag( vcov( estLpmLin ) ) )["age"] )
urbinEla( coef( estLpmLin ), xMeanLin, xPos = 3, model = "lpm",
  sqrt( diag( vcov( estLpmLin ) ) ) )

### quadratic in age
estLpmQuad <- lm( lfp ~ kids + age + I(age^2) + educ, 
  data = Mroz87 )
summary( estLpmQuad )
# mean values of the explanatory variables
xMeanQuad <- c( xMeanLin[ 1:3 ], xMeanLin[3]^2, xMeanLin[4] )
# semi-elasticity of age without standard errors
urbinEla( coef( estLpmQuad )[ c( "age", "I(age^2)" ) ], xMeanQuad[ "age" ], 
  xPos = c( 1, 2 ), model = "lpm" )
urbinEla( coef( estLpmQuad ), xMeanQuad, xPos = c( 3, 4 ), model = "lpm" )
# semi-elasticity of age based on numerical derivation
100 * ( predict( estLpmQuad, 
  newdata = as.data.frame( t( xMeanQuad * c( 1, 1, 1.005, 1.005^2, 1 ) ) ) ) -
    predict( estLpmQuad, 
      newdata = as.data.frame( t( xMeanQuad * c( 1, 1, 0.995, 0.995^2, 1 ) ) ) ) )
# partial derivatives of the semi-elasticity wrt the coefficients
urbinEla( coef( estLpmQuad )[ c( "age", "I(age^2)" ) ], xMeanAgeAttr, 
  xPos = c( 1, 2 ), model = "lpm" )
xMeanQuadAttr <- xMeanQuad
attr( xMeanQuadAttr, "derivOnly" ) <- 1 
urbinEla( coef( estLpmQuad ), xMeanQuadAttr, xPos = c( 3, 4 ), model = "lpm" )
# numerically computed partial derivatives of the semi-elasticity wrt the coefficients
numericGradient( urbinEla, t0 = coef( estLpmQuad )[ c( "age", "I(age^2)" ) ], 
  allXVal = xMeanQuad[ "age" ], xPos = c( 1, 2 ), model = "lpm" )
numericGradient( urbinEla, t0 = coef( estLpmQuad ), 
  allXVal = xMeanQuad, xPos = c( 3, 4 ), model = "lpm" )
# semi-elasticity of age with standard errors (full covariance matrix)
urbinEla( coef( estLpmQuad )[ c( "age", "I(age^2)" ) ], xMeanQuad["age"], 
  xPos = c( 1, 2 ), model = "lpm", 
  vcov( estLpmQuad )[ c( "age", "I(age^2)" ), c( "age", "I(age^2)" ) ] )
urbinEla( coef( estLpmQuad ), xMeanQuad, xPos = c( 3, 4 ), model = "lpm", 
  vcov( estLpmQuad ) )
# semi-elasticity of age with standard errors (only standard errors)
urbinEla( coef( estLpmQuad )[ c( "age", "I(age^2)" ) ], xMeanQuad[ "age" ], 
  xPos = c( 1, 2 ), model = "lpm", 
  sqrt( diag( vcov( estLpmQuad ) ) )[ c( "age", "I(age^2)" ) ] )
urbinEla( coef( estLpmQuad ), xMeanQuad, xPos = c( 3, 4 ), model = "lpm",
  sqrt( diag( vcov( estLpmQuad ) ) ) )
# approximate covariance between the coefficient of the linear term and 
# the coefficient of the quadratic term based on the original data
se <- sqrt( diag( vcov( estLpmQuad ) ) )
X <- cbind( Mroz87$age, Mroz87$age^2, 1 )
XXinv <- solve( t(X) %*% X )
sigmaSq <- sqrt( ( se["age"]^2 / XXinv[1,1] ) * ( se["I(age^2)"]^2 / XXinv[2,2] ) )
vcovApp <- diag( se^2 )
rownames( vcovApp ) <- colnames( vcovApp ) <- names( se )
vcovApp[ "age", "I(age^2)" ] <- vcovApp[ "I(age^2)", "age" ] <- 
  sigmaSq * XXinv[1,2]
urbinEla( coef( estLpmQuad )[ c( "age", "I(age^2)" ) ], xMeanQuad["age"], 
  xPos = c( 1, 2 ), model = "lpm", 
  vcovApp[ c( "age", "I(age^2)" ), c( "age", "I(age^2)" ) ] )
urbinEla( coef( estLpmQuad ), xMeanQuad, xPos = c( 3, 4 ), model = "lpm", 
  vcovApp )
# approximate covariance between the coefficient of the linear term and 
# the coefficient of the quadratic term based on simulated data
se <- sqrt( diag( vcov( estLpmQuad ) ) )
set.seed( 123 )
x <- rnorm( 1000, xMeanQuad[ "age" ], sd( Mroz87$age ) )
X <- cbind( x, x^2, 1 )
XXinv <- solve( t(X) %*% X )
sigmaSq <- sqrt( ( se["age"]^2 / XXinv[1,1] ) * ( se["I(age^2)"]^2 / XXinv[2,2] ) )
vcovApp <- diag( se^2 )
rownames( vcovApp ) <- colnames( vcovApp ) <- names( se )
vcovApp[ "age", "I(age^2)" ] <- vcovApp[ "I(age^2)", "age" ] <- 
  sigmaSq * XXinv[1,2]
urbinEla( coef( estLpmQuad )[ c( "age", "I(age^2)" ) ], xMeanQuad["age"], 
  xPos = c( 1, 2 ), model = "lpm", 
  vcovApp[ c( "age", "I(age^2)" ), c( "age", "I(age^2)" ) ] )
urbinEla( coef( estLpmQuad ), xMeanQuad, xPos = c( 3, 4 ), model = "lpm", 
  vcovApp )
urbinEla( coef( estLpmQuad )[ c( "age", "I(age^2)" ) ], xMeanQuad["age"], 
  xPos = c( 1, 2 ), model = "lpm",
  sqrt( diag( vcov( estLpmQuad ) ) )[ c( "age", "I(age^2)" ) ],
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ) )
urbinEla( coef( estLpmQuad ), xMeanQuad, xPos = c( 3, 4 ), model = "lpm",
  sqrt( diag( vcov( estLpmQuad ) ) ),
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ) )

### age is interval-coded (age is in the range 30-60)
# create dummy variables for age intervals
Mroz87$age30.37 <- Mroz87$age >= 30 & Mroz87$age <= 37
Mroz87$age38.44 <- Mroz87$age >= 38 & Mroz87$age <= 44
Mroz87$age45.52 <- Mroz87$age >= 45 & Mroz87$age <= 52
Mroz87$age53.60 <- Mroz87$age >= 53 & Mroz87$age <= 60
all.equal( 
  Mroz87$age30.37 + Mroz87$age38.44 + Mroz87$age45.52 + Mroz87$age53.60,
  rep( 1, nrow( Mroz87 ) ) )
# estimation
estLpmInt <- lm( lfp ~ kids + age30.37 + age38.44 + age53.60 + educ, 
  data = Mroz87 )
summary( estLpmInt )
# coefficients of the 'intervals'
coefLpmInt <- c( coef( estLpmInt )[3:4], 0, coef( estLpmInt )[5] )
# mean values of the explanatory variables
xMeanInt <- c( xMeanLin[1:2], mean( Mroz87$age30.37 ), 
  mean( Mroz87$age38.44 ), mean( Mroz87$age53.60 ), xMeanLin[4] )
# mean shares of the 'intervals'
xMeanIntShares <- c( xMeanInt[3:4], 1 - sum( xMeanInt[3:5] ), xMeanInt[5] )
# semi-elasticity of age without standard errors
urbinElaInt( coef( estLpmInt )[3:5], xMeanInt[3:5],
  c( 30, 37.5, 44.5, 52.5, 60 ), xPos = c( 1, 2, 0, 3 ), model = "lpm" )
urbinElaInt( coef( estLpmInt ), xMeanInt,
  c( 30, 37.5, 44.5, 52.5, 60 ), xPos = c( 3, 4, 0, 5 ), model = "lpm" )
# partial derivatives of the semi-elasticity wrt the coefficients
xMeanIntShares3Attr <- xMeanInt[3:5]
attr( xMeanIntShares3Attr, "derivOnly" ) <- 1 
urbinElaInt( coef( estLpmInt )[3:5], xMeanIntShares3Attr, 
  c( 30, 37.5, 44.5, 52.5, 60 ), xPos = c( 1, 2, 0, 3 ), model = "lpm" )
xMeanIntAttr <- xMeanInt
attr( xMeanIntAttr, "derivOnly" ) <- 1 
urbinElaInt( coef( estLpmInt ), xMeanIntAttr, 
  c( 30, 37.5, 44.5, 52.5, 60 ), xPos = c( 3, 4, 0, 5 ), model = "lpm" )
# numerically computed partial derivatives of the semi-elasticity wrt the coefficients
numericGradient( urbinElaInt, t0 = coef( estLpmInt )[3:5], allXVal = xMeanInt[3:5], 
  xBound = c( 30, 37.5, 44.5, 52.5, 60 ), xPos = c( 1, 2, 0, 3 ), 
  model = "lpm" )
numericGradient( urbinElaInt, t0 = coef( estLpmInt ), allXVal = xMeanInt, 
  xBound = c( 30, 37.5, 44.5, 52.5, 60 ), xPos = c( 3, 4, 0, 5 ), 
  model = "lpm" )
# semi-elasticity of age with standard errors (full covariance matrix)
vcovLpmInt <- vcov( estLpmInt )
vcovLpmInt <- rbind( vcovLpmInt[ 3:4, ], 0, vcovLpmInt[ 5, ] )
vcovLpmInt <- cbind( vcovLpmInt[ , 3:4 ], 0, vcovLpmInt[ , 5 ] )
urbinElaInt( coef( estLpmInt )[3:5], xMeanInt[3:5],
  c( 30, 37.5, 44.5, 52.5, 60 ), xPos = c( 1, 2, 0, 3 ), model = "lpm",
  allCoefVcov = vcov( estLpmInt )[ 3:5, 3:5 ] )
urbinElaInt( coef( estLpmInt ), xMeanInt,
  c( 30, 37.5, 44.5, 52.5, 60 ), xPos = c( 3, 4, 0, 5 ), model = "lpm",
  allCoefVcov = vcov( estLpmInt ) )
# semi-elasticity of age with standard errors (only standard errors)
urbinElaInt( coef( estLpmInt )[3:5], xMeanInt[3:5],
  c( 30, 37.5, 44.5, 52.5, 60 ), xPos = c( 1, 2, 0, 3 ), model = "lpm",
  allCoefVcov = sqrt( diag( vcov( estLpmInt ) ) )[3:5] )
urbinElaInt( coef( estLpmInt ), xMeanInt,
  c( 30, 37.5, 44.5, 52.5, 60 ), xPos = c( 3, 4, 0, 5 ), model = "lpm",
  allCoefVcov = sqrt( diag( vcov( estLpmInt ) ) ) )


### effect of age changing between discrete intervals 
### if age is used as linear explanatory variable 
# mean values of the 'other' explanatory variables
xMeanLinInt <- c( xMeanLin[ 1:2 ], NA, xMeanLin[4] )
# effects of age changing from the 30-40 interval to the 50-60 interval
# without standard errors
urbinEffInt( coef( estLpmLin )[3], NA, c( 30, 40 ), c( 50, 60 ), xPos = 1,
  model = "lpm" )
urbinEffInt( coef( estLpmLin ), NA, c( 30, 40 ), c( 50, 60 ), xPos = 3,
  model = "lpm" )
# effects of age changing from the 30-40 interval to the 50-60 interval
# based on predicted values
predict( estLpmLin, 
  newdata = as.data.frame( t( replace( xMeanLin, 3, 55 ) ) ) ) -
  predict( estLpmLin, 
    newdata = as.data.frame( t( replace( xMeanLin, 3, 35 ) ) ) )
# partial derivatives of the semi-elasticity wrt the coefficients
naAttr <- NA
attr( naAttr, "derivOnly" ) <- 1 
urbinEffInt( coef( estLpmLin ), naAttr, xPos = 3,
  c( 30, 40 ), c( 50, 60 ), model = "lpm" )
# numerically computed partial derivatives of the semi-elasticity wrt the coefficients
numericGradient( urbinEffInt, t0 = coef( estLpmLin )[3], allXVal = NA,
  refBound = c( 30, 40 ), intBound = c( 50, 60 ), xPos = 1, model = "lpm" )
numericGradient( urbinEffInt, t0 = coef( estLpmLin ), allXVal = NA,
  refBound = c( 30, 40 ), intBound = c( 50, 60 ), xPos = 3, model = "lpm" )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (full covariance matrix) 
urbinEffInt( coef( estLpmLin ), NA,
  c( 30, 40 ), c( 50, 60 ), xPos = 3, model = "lpm",
  allCoefVcov = vcov( estLpmLin ) )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (only standard errors) 
urbinEffInt( coef( estLpmLin )[3], NA, c( 30, 40 ), c( 50, 60 ), xPos = 1, 
  model = "lpm", allCoefVcov = sqrt( diag( vcov( estLpmLin ) ) )[3] )
urbinEffInt( coef( estLpmLin ), NA, c( 30, 40 ), c( 50, 60 ), xPos = 3, 
  model = "lpm", allCoefVcov = sqrt( diag( vcov( estLpmLin ) ) ) )


### effect of age changing between discrete intervals 
### if age is used as linear and quadratic explanatory variable 
# mean values of the 'other' explanatory variables
xMeanQuadInt <- c( xMeanLin[ 1:2 ], NA, NA, xMeanLin[4] )
# effects of age changing from the 30-40 interval to the 50-60 interval
# without standard errors
urbinEffInt( coef( estLpmQuad )[3:4], NA,
  c( 30, 40 ), c( 50, 60 ), xPos = 1:2, model = "lpm" )
urbinEffInt( coef( estLpmQuad ), NA,
  c( 30, 40 ), c( 50, 60 ), xPos = 3:4, model = "lpm" )
# effects of age changing from the 30-40 interval to the 50-60 interval
# based on predicted values
predict( estLpmQuad, 
  newdata = as.data.frame( t( replace( xMeanQuad, 3:4, c( 55, 55^2 ) ) ) ), 
  type = "response" ) -
  predict( estLpmQuad, 
    newdata = as.data.frame( t( replace( xMeanQuad, 3:4, c( 35, 35^2 ) ) ) ), 
    type = "response" )
# partial derivatives of the effect wrt the coefficients
urbinEffInt( coef( estLpmQuad ), naAttr, xPos = c( 3, 4 ),
  c( 30, 40 ), c( 50, 60 ), model = "lpm" )
# numerically computed partial derivatives of the effect wrt the coefficients
numericGradient( urbinEffInt, t0 = coef( estLpmQuad )[3:4], allXVal = NA,
  refBound = c( 30, 40 ), intBound = c( 50, 60 ), xPos = 1:2, model = "lpm" )
numericGradient( urbinEffInt, t0 = coef( estLpmQuad ),
  refBound = c( 30, 40 ), intBound = c( 50, 60 ), xPos = 3:4, model = "lpm" )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (full covariance matrix) 
urbinEffInt( coef( estLpmQuad )[3:4], NA, c( 30, 40 ), c( 50, 60 ), 
  xPos = 1:2, model = "lpm", allCoefVcov = vcov( estLpmQuad )[3:4,3:4] )
urbinEffInt( coef( estLpmQuad ), NA, c( 30, 40 ), c( 50, 60 ), 
  xPos = 3:4, model = "lpm", allCoefVcov = vcov( estLpmQuad ) )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (only standard errors) 
urbinEffInt( coef( estLpmQuad )[3:4], NA, c( 30, 40 ), c( 50, 60 ), 
  xPos = 1:2, model = "lpm", 
  allCoefVcov = sqrt( diag( vcov( estLpmQuad ) ) )[3:4] )
urbinEffInt( coef( estLpmQuad ), NA, c( 30, 40 ), c( 50, 60 ), 
  xPos = 3:4, model = "lpm", 
  allCoefVcov = sqrt( diag( vcov( estLpmQuad ) ) ) )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (standard errors + mean value and standard deviation of age)
urbinEffInt( coef( estLpmQuad ), NA, xPos = c( 3, 4 ),
  c( 30, 40 ), c( 50, 60 ), model = "lpm", 
  allCoefVcov = sqrt( diag( vcov( estLpmQuad ) ) ),
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ) )

### grouping and re-basing categorical variables
### effects of age changing from the 30-44 category to the 53-60 category
# without standard errors
urbinEffCat( coef( estLpmInt )[3:5], xMeanInt[3:5], 1:3, c( -1, -1, 1, 0 ), 
  model = "lpm" )
urbinEffCat( coef( estLpmInt ), xMeanInt, 3:5, c( -1, -1, 1, 0 ), 
  model = "lpm" )
# effects calculated based on predicted values
names( xMeanInt ) <- sub( "TRUE", "", names( coef( estLpmInt ) ) )
df30.37 <- df38.44 <- df45.52 <- df53.60 <- as.data.frame( t( xMeanInt ) ) 
df30.37[ , 3:5 ] <- c( TRUE, FALSE, FALSE )
df38.44[ , 3:5 ] <- c( FALSE, TRUE, FALSE )
df45.52[ , 3:5 ] <- c( FALSE, FALSE, FALSE )
df53.60[ , 3:5 ] <- c( FALSE, FALSE, TRUE )
predict( estLpmInt, newdata = df53.60 ) -
  sum( Mroz87$age30.37 ) / sum( Mroz87$age30.37 + Mroz87$age38.44 ) *
  predict( estLpmInt, newdata = df30.37 ) -
  sum( Mroz87$age38.44 ) / sum( Mroz87$age30.37 + Mroz87$age38.44 ) *
  predict( estLpmInt, newdata = df38.44 )
# partial derivatives of the effect wrt the coefficients
urbinEffCat( coef( estLpmInt )[3:5], xMeanIntShares3Attr, 1:3, c( -1, -1, 1, 0 ), 
  model = "lpm" )
urbinEffCat( coef( estLpmInt ), xMeanIntAttr, 3:5, c( -1, -1, 1, 0 ), 
  model = "lpm" )
# numerically computed partial derivatives of the effect wrt the coefficients
numericGradient( urbinEffCat, t0 = coef( estLpmInt )[3:5], xPos = 1:3,
  allXVal = xMeanInt[3:5], xGroups = c( -1, -1, 1, 0 ), model = "lpm" )
numericGradient( urbinEffCat, t0 = coef( estLpmInt ), xPos = 3:5,
  allXVal = xMeanInt, xGroups = c( -1, -1, 1, 0 ), model = "lpm" )
# with full covariance matrix
urbinEffCat( coef( estLpmInt )[3:5], xMeanInt[3:5], 1:3, c( -1, -1, 1, 0 ), 
  model = "lpm", vcov( estLpmInt )[3:5, 3:5] )
urbinEffCat( coef( estLpmInt ), xMeanInt, 3:5, c( -1, -1, 1, 0 ), 
  model = "lpm", vcov( estLpmInt ) )
# with standard errors only
urbinEffCat( coef( estLpmInt )[3:5], xMeanInt[3:5], 1:3, c( -1, -1, 1, 0 ), 
  model = "lpm", sqrt( diag( vcov( estLpmInt ) ) )[3:5] )
urbinEffCat( coef( estLpmInt ), xMeanInt, 3:5, c( -1, -1, 1, 0 ), 
  model = "lpm", sqrt( diag( vcov( estLpmInt ) ) ) )
### effects of age changing from the 53-60 category to the 38-52 category
# without standard errors
urbinEffCat( coef( estLpmInt )[3:5], xMeanInt[3:5], 1:3, c( 0, 1, -1, 1 ), 
  model = "lpm" )
urbinEffCat( coef( estLpmInt ), xMeanInt, 3:5, c( 0, 1, -1, 1 ), 
  model = "lpm" )
# effects calculated based on predicted values
sum( Mroz87$age38.44 ) / sum( Mroz87$age38.44 + Mroz87$age45.52 ) *
  predict( estLpmInt, newdata = df38.44 ) +
  sum( Mroz87$age45.52 ) / sum( Mroz87$age38.44 + Mroz87$age45.52 ) *
  predict( estLpmInt, newdata = df45.52 ) -
  predict( estLpmInt, newdata = df53.60 )
# partial derivatives of the effect wrt the coefficients
urbinEffCat( coef( estLpmInt )[3:5], xMeanIntShares3Attr, 1:3, c( 0, 1, -1, 1 ), 
  model = "lpm" )
urbinEffCat( coef( estLpmInt ), xMeanIntAttr, 3:5, c( 0, 1, -1, 1 ), 
  model = "lpm" )
# numerically computed partial derivatives of the effect wrt the coefficients
numericGradient( urbinEffCat, t0 = coef( estLpmInt )[3:5], xPos = 1:3,
  allXVal = xMeanInt[3:5], xGroups = c( 0, 1, -1, 1 ), model = "lpm" )
numericGradient( urbinEffCat, t0 = coef( estLpmInt ), xPos = 3:5,
  allXVal = xMeanInt, xGroups = c( 0, 1, -1, 1 ), model = "lpm" )
# with full covariance matrix
urbinEffCat( coef( estLpmInt )[3:5], xMeanInt[3:5], 1:3, c( 0, 1, -1, 1 ), 
  model = "lpm", vcov( estLpmInt )[3:5,3:5] )
urbinEffCat( coef( estLpmInt ), xMeanInt, 3:5, c( 0, 1, -1, 1 ), 
  model = "lpm", vcov( estLpmInt ) )
# with standard errors only
urbinEffCat( coef( estLpmInt )[3:5], xMeanInt[3:5], 1:3, c( 0, 1, -1, 1 ), 
  model = "lpm", sqrt( diag( vcov( estLpmInt ) ) )[3:5] )
urbinEffCat( coef( estLpmInt ), xMeanInt, 3:5, c( 0, 1, -1, 1 ), 
  model = "lpm", sqrt( diag( vcov( estLpmInt ) ) ) )

