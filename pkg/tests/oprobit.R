library( "urbin" )
library( "maxLik" )
library( "MASS" )

# load data set
data( "Mroz87", package = "sampleSelection" )

# create dummy variable for kids
Mroz87$kids <- as.numeric( Mroz87$kids5 > 0 | Mroz87$kids618 > 0 )

### create categorical variable
Mroz87$lfp3 <- factor( ifelse( Mroz87$hours == 0, "no",
  ifelse( Mroz87$hours <= 1300, "part", "full" ) ),
  levels = c( "no", "part", "full" ), ordered = TRUE )
table( Mroz87$lfp3 )
all.equal( Mroz87$lfp3 == "no", Mroz87$lfp == 0 )

### linear in age
estOProbitLin <- polr( lfp3 ~ kids + age + educ, data = Mroz87,
  method = "probit", Hess = TRUE )
summary( estOProbitLin )
# vector of coefficients and their variance covariance matrix (as if it were 
# a binary probit model with 'no' = 0 and 'part' = 'full' = 1)
coefOProbitLinNV <- coef( summary( estOProbitLin ) )[
  c( "no|part", "kids", "age", "educ" ), 1 ]
vcovOProbitLinNV <- vcov( estOProbitLin )[
  c( "no|part", "kids", "age", "educ" ),
  c( "no|part", "kids", "age", "educ" ) ]
# the same as above but using the negative threshold as intercept
coefOProbitLinNC <- coefOProbitLinNV * c( -1, 1, 1, 1 )
vcovOProbitLinNC <- diag( c( -1, 1, 1, 1 ) ) %*% vcovOProbitLinNV %*% 
  diag( c( -1, 1, 1, 1 ) )
# mean values of the explanatory variables
xMeanLinNC <- c( 1, colMeans( Mroz87[ , c( "kids", "age", "educ" ) ] ) )
# the same as above but setting the intercept to minus one
xMeanLinNV <- xMeanLinNC * c( -1, 1, 1, 1 )
# semi-elasticity of age without standard errors
urbinEla( coefOProbitLinNV, xMeanLinNV, xPos = 3, model = "probit" )
urbinEla( coefOProbitLinNC, xMeanLinNC, xPos = 3, model = "probit" )
# semi-elasticity of age based on numerical derivation
Mroz87Lower <- as.data.frame( t( xMeanLinNC * c( 1, 1, 0.995, 1 ) ) )
Mroz87Upper <- as.data.frame( t( xMeanLinNC * c( 1, 1, 1.005, 1 ) ) )
elaLinNum <- 100 * ( 
  predict( estOProbitLin, newdata = Mroz87Upper, type = "probs" ) -
    predict( estOProbitLin, newdata = Mroz87Lower, type = "probs" ) )
print( elaLinNum )
print( sum( elaLinNum[ c( "part", "full" ) ] ) )
# partial derivatives of the semi-elasticity wrt the coefficients
xMeanLinNVAttr <- xMeanLinNV
attr( xMeanLinNVAttr, "derivOnly" ) <- 1 
urbinEla( coefOProbitLinNV, xMeanLinNVAttr, 3, 
  seSimplify = FALSE, model = "probit" )
xMeanLinNCAttr <- xMeanLinNC
attr( xMeanLinNCAttr, "derivOnly" ) <- 1 
urbinEla( coefOProbitLinNC, xMeanLinNCAttr, 3, 
  seSimplify = FALSE, model = "probit" )
# numerically computed partial derivatives of the semi-elasticity wrt the coefficients
numericGradient( urbinEla, t0 = coefOProbitLinNV, 
  allXVal = xMeanLinNV, xPos = 3, model = "probit" )
numericGradient( urbinEla, t0 = coefOProbitLinNC, 
  allXVal = xMeanLinNC, xPos = 3, model = "probit" )
# simplified partial derivatives of the semi-elasticity wrt the coefficients
urbinEla( coefOProbitLinNV, xMeanLinNVAttr, 3, 
  model = "probit", seSimplify = TRUE )
urbinEla( coefOProbitLinNC, xMeanLinNCAttr, 3, 
  model = "probit", seSimplify = TRUE )
# semi-elasticity of age with standard errors (full covariance matrix)
urbinEla( coefOProbitLinNV, xMeanLinNV, 3, model = "probit", 
  vcovOProbitLinNV )
urbinEla( coefOProbitLinNC, xMeanLinNC, 3, model = "probit", 
  vcovOProbitLinNC )
# semi-elasticity of age with standard errors (only standard errors)
urbinEla( coefOProbitLinNV, xMeanLinNV, 3, model = "probit",
  sqrt( diag( vcovOProbitLinNV ) ), seSimplify = FALSE )
urbinEla( coefOProbitLinNC, xMeanLinNC, 3, model = "probit",
  sqrt( diag( vcovOProbitLinNC ) ), seSimplify = FALSE )
# semi-elasticity of age with standard errors (only standard errors, simplified)
urbinEla( coefOProbitLinNV, xMeanLinNV, 3, model = "probit", 
  sqrt( diag( vcovOProbitLinNV ) ) )
urbinEla( coefOProbitLinNC, xMeanLinNC, 3, model = "probit", 
  sqrt( diag( vcovOProbitLinNC ) ) )


### quadratic in age
estOProbitQuad <- polr( lfp3 ~ kids + age + I(age^2) + educ, 
  data = Mroz87, method = "probit", Hess = TRUE )
summary( estOProbitQuad )
# vector of coefficients and their variance covariance matrix (as if it were 
# a binary probit model with 'no' = 0 and 'part' = 'full' = 1)
coefOProbitQuadNV <- coef( summary( estOProbitQuad ) )[
  c( "no|part", "kids", "age", "I(age^2)", "educ" ), 1 ]
vcovOProbitQuadNV <- vcov( estOProbitQuad )[
  c( "no|part", "kids", "age", "I(age^2)", "educ" ),
  c( "no|part", "kids", "age", "I(age^2)", "educ" ) ]
# the same as above but using the negative threshold as intercept
coefOProbitQuadNC <- coefOProbitQuadNV * c( -1, 1, 1, 1, 1 )
vcovOProbitQuadNC <- diag( c( -1, 1, 1, 1, 1 ) ) %*% vcovOProbitQuadNV %*% 
  diag( c( -1, 1, 1, 1, 1 ) )
# mean values of the explanatory variables
xMeanQuadNC <- c( xMeanLinNC[ 1:3 ], xMeanLinNC[3]^2, xMeanLinNC[4] )
# the same as above but setting the intercept to minus one
xMeanQuadNV <- xMeanQuadNC * c( -1, 1, 1, 1, 1 )
# semi-elasticity of age without standard errors
urbinEla( coefOProbitQuadNV, xMeanQuadNV, c( 3, 4 ), 
  model = "probit" )
urbinEla( coefOProbitQuadNC, xMeanQuadNC, c( 3, 4 ), 
  model = "probit" )
# semi-elasticity of age based on numerical derivation
Mroz87Lower <- as.data.frame( 
  t( xMeanQuadNC * c( 1, 1, 0.995, 0.995^2, 1 ) ) )
Mroz87Upper <- as.data.frame( 
  t( xMeanQuadNC * c( 1, 1, 1.005, 1.005^2, 1 ) ) )
elaQuadNum <- 100 * ( 
  predict( estOProbitQuad, newdata = Mroz87Upper, type = "probs" ) -
    predict( estOProbitQuad, newdata = Mroz87Lower, type = "probs" ) )
print( elaQuadNum )
print( sum( elaQuadNum[ c( "part", "full" ) ] ) )
# partial derivatives of the semi-elasticity wrt the coefficients
xMeanQuadNVAttr <- xMeanQuadNV
attr( xMeanQuadNVAttr, "derivOnly" ) <- 1 
urbinEla( coefOProbitQuadNV, xMeanQuadNVAttr, c( 3, 4 ), 
  model = "probit", seSimplify = FALSE )
xMeanQuadNCAttr <- xMeanQuadNC
attr( xMeanQuadNCAttr, "derivOnly" ) <- 1 
urbinEla( coefOProbitQuadNC, xMeanQuadNCAttr, c( 3, 4 ), 
  model = "probit", seSimplify = FALSE )
# numerically computed partial derivatives of the semi-elasticity wrt the coefficients
numericGradient( urbinEla, t0 = coefOProbitQuadNV, 
  allXVal = xMeanQuadNV, xPos = c( 3, 4 ), model = "probit" )
numericGradient( urbinEla, t0 = coefOProbitQuadNC, 
  allXVal = xMeanQuadNC, xPos = c( 3, 4 ), model = "probit" )
# simplified partial derivatives of the semi-elasticity wrt the coefficients
urbinEla( coefOProbitQuadNV, xMeanQuadNVAttr, c( 3, 4 ), 
  model = "probit", seSimplify = TRUE )
urbinEla( coefOProbitQuadNC, xMeanQuadNCAttr, c( 3, 4 ), 
  model = "probit", seSimplify = TRUE )
# semi-elasticity of age with standard errors (full covariance matrix)
urbinEla( coefOProbitQuadNV, xMeanQuadNV, c( 3, 4 ), 
  model = "probit", vcovOProbitQuadNV )
urbinEla( coefOProbitQuadNC, xMeanQuadNC, c( 3, 4 ), 
  model = "probit", vcovOProbitQuadNC )
# semi-elasticity of age with standard errors (only standard errors)
urbinEla( coefOProbitQuadNV, xMeanQuadNV, c( 3, 4 ), 
  model = "probit", sqrt( diag( vcovOProbitQuadNV ) ), 
  seSimplify = FALSE )
urbinEla( coefOProbitQuadNC, xMeanQuadNC, c( 3, 4 ), 
  model = "probit", sqrt( diag( vcovOProbitQuadNC ) ), 
  seSimplify = FALSE )
# semi-elasticity of age with standard errors (only standard errors, simplified)
urbinEla( coefOProbitQuadNV, xMeanQuadNV, c( 3, 4 ), 
  model = "probit", sqrt( diag( vcovOProbitQuadNV ) ) )
urbinEla( coefOProbitQuadNC, xMeanQuadNC, c( 3, 4 ), 
  model = "probit", sqrt( diag( vcovOProbitQuadNC ) ) )
# semi-elasticity of age with standard errors (only standard errors, xMeanSd)
urbinEla( coefOProbitQuadNV, xMeanQuadNV, c( 3, 4 ), 
  model = "probit", sqrt( diag( vcovOProbitQuadNV ) ),
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ),
  seSimplify = FALSE )
urbinEla( coefOProbitQuadNC, xMeanQuadNC, c( 3, 4 ), 
  model = "probit", sqrt( diag( vcovOProbitQuadNC ) ),
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ),
  seSimplify = FALSE )
# semi-elasticity of age with standard errors (only standard errors, xMeanSd, simplified)
urbinEla( coefOProbitQuadNV, xMeanQuadNV, c( 3, 4 ), 
  model = "probit", sqrt( diag( vcovOProbitQuadNV ) ),
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ) )
urbinEla( coefOProbitQuadNC, xMeanQuadNC, c( 3, 4 ), 
  model = "probit", sqrt( diag( vcovOProbitQuadNC ) ),
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
estOProbitInt <- polr( lfp3 ~ kids + age30.37 + age38.44 + age53.60 + educ, 
  data = Mroz87, method = "probit", Hess = TRUE )
summary( estOProbitInt )
# vector of coefficients and their variance covariance matrix (as if it were 
# a binary probit model with 'no' = 0 and 'part' = 'full' = 1)
coefOProbitIntNV <- coef( summary( estOProbitInt ) )[
  c( "no|part", "kids", "age30.37TRUE", "age38.44TRUE", "age53.60TRUE", "educ" ), 
  1 ]
vcovOProbitIntNV <- vcov( estOProbitInt )[
  c( "no|part", "kids", "age30.37TRUE", "age38.44TRUE", "age53.60TRUE", "educ" ),
  c( "no|part", "kids", "age30.37TRUE", "age38.44TRUE", "age53.60TRUE", "educ" ) ]
# the same as above but using the negative threshold as intercept
coefOProbitIntNC <- coefOProbitIntNV * c( -1, 1, 1, 1, 1, 1 )
vcovOProbitIntNC <- diag( c( -1, 1, 1, 1, 1, 1 ) ) %*% vcovOProbitIntNV %*% 
  diag( c( -1, 1, 1, 1, 1, 1 ) )
# mean values of the explanatory variables
xMeanIntNC <- c( xMeanLinNC[1:2], mean( Mroz87$age30.37 ), 
  mean( Mroz87$age38.44 ), mean( Mroz87$age53.60 ), xMeanLinNC[4] )
# the same as above but setting the intercept to minus one
xMeanIntNV <- xMeanIntNC * c( -1, 1, 1, 1, 1, 1 )
# semi-elasticity of age without standard errors
urbinElaInt( coefOProbitIntNV, xMeanIntNV,
  c( 3, 4, 0, 5 ), c( 30, 37.5, 44.5, 52.5, 60 ), model = "probit" )
urbinElaInt( coefOProbitIntNC, xMeanIntNC,
  c( 3, 4, 0, 5 ), c( 30, 37.5, 44.5, 52.5, 60 ), model = "probit" )
# semi-elasticities based on numerical derivation
Mroz87Lower <- Mroz87
Mroz87Lower$age <- Mroz87$age * 0.95
Mroz87Lower$age30.37 <- Mroz87Lower$age <= 37.5
Mroz87Lower$age38.44 <- Mroz87Lower$age > 37.5 & Mroz87Lower$age <= 44.5
Mroz87Lower$age45.52 <- Mroz87Lower$age > 44.5 & Mroz87Lower$age <= 52.5
Mroz87Lower$age53.60 <- Mroz87Lower$age > 52.5 
all.equal( 
  Mroz87Lower$age30.37 + Mroz87Lower$age38.44 + Mroz87Lower$age45.52 + 
    Mroz87Lower$age53.60, rep( 1, nrow( Mroz87 ) ) )
Mroz87Upper <- Mroz87
Mroz87Upper$age <- Mroz87$age * 1.05
Mroz87Upper$age30.37 <- Mroz87Upper$age <= 37.5
Mroz87Upper$age38.44 <- Mroz87Upper$age > 37.5 & Mroz87Upper$age <= 44.5
Mroz87Upper$age45.52 <- Mroz87Upper$age > 44.5 & Mroz87Upper$age <= 52.5
Mroz87Upper$age53.60 <- Mroz87Upper$age > 52.5 
all.equal( 
  Mroz87Upper$age30.37 + Mroz87Upper$age38.44 + Mroz87Upper$age45.52 + 
    Mroz87Upper$age53.60, rep( 1, nrow( Mroz87 ) ) )
elaIntNum <- 10 * ( colMeans( 
  predict( estOProbitInt, newdata = Mroz87Upper, type = "probs" ) ) -
    colMeans(
      predict( estOProbitInt, newdata = Mroz87Lower, type = "probs" ) ) )
print( elaIntNum )
print( sum( elaIntNum[ c( "part", "full" ) ] ) )
# partial derivatives of the semi-elasticity wrt the coefficients
xMeanIntNVAttr <- xMeanIntNV
attr( xMeanIntNVAttr, "derivOnly" ) <- 1 
urbinElaInt( coefOProbitIntNV, xMeanIntNVAttr,
  c( 3, 4, 0, 5 ), c( 30, 37.5, 44.5, 52.5, 60 ), model = "probit" )
xMeanIntNCAttr <- xMeanIntNC
attr( xMeanIntNCAttr, "derivOnly" ) <- 1 
urbinElaInt( coefOProbitIntNC, xMeanIntNCAttr,
  c( 3, 4, 0, 5 ), c( 30, 37.5, 44.5, 52.5, 60 ), model = "probit" )
# numerically computed partial derivatives of the semi-elasticity wrt the coefficients
numericGradient( urbinElaInt, t0 = coefOProbitIntNV, 
  allXVal = xMeanIntNV, xPos = c( 3, 4, 0, 5 ), 
  xBound = c( 30, 37.5, 44.5, 52.5, 60 ), model = "probit" )
numericGradient( urbinElaInt, t0 = coefOProbitIntNC, 
  allXVal = xMeanIntNC, xPos = c( 3, 4, 0, 5 ), 
  xBound = c( 30, 37.5, 44.5, 52.5, 60 ), model = "probit" )
# semi-elasticity of age with standard errors (full covariance matrix)
urbinElaInt( coefOProbitIntNV, xMeanIntNV,
  c( 3, 4, 0, 5 ), c( 30, 37.5, 44.5, 52.5, 60 ), model = "probit",
  allCoefVcov = vcovOProbitIntNV )
urbinElaInt( coefOProbitIntNC, xMeanIntNC,
  c( 3, 4, 0, 5 ), c( 30, 37.5, 44.5, 52.5, 60 ), model = "probit",
  allCoefVcov = vcovOProbitIntNC )
# semi-elasticity of age with standard errors (only standard errors)
urbinElaInt( coefOProbitIntNV, xMeanIntNV,
  c( 3, 4, 0, 5 ), c( 30, 37.5, 44.5, 52.5, 60 ), model = "probit",
  allCoefVcov = sqrt( diag( vcovOProbitIntNV ) ) )
urbinElaInt( coefOProbitIntNC, xMeanIntNC,
  c( 3, 4, 0, 5 ), c( 30, 37.5, 44.5, 52.5, 60 ), model = "probit",
  allCoefVcov = sqrt( diag( vcovOProbitIntNC ) ) )


### effect of age changing between discrete intervals 
### if age is used as linear explanatory variable 
# mean values of the 'other' explanatory variables
xMeanLinIntNV <- c( xMeanLinNV[ 1:2 ], NA, xMeanLinNV[4] )
xMeanLinIntNC <- c( xMeanLinNC[ 1:2 ], NA, xMeanLinNC[4] )
# effects of age changing from the 30-40 interval to the 50-60 interval
# without standard errors
urbinEffInt( coefOProbitLinNV, allXVal = xMeanLinIntNV, 
  xPos = 3, refBound = c( 30, 40 ), intBound = c( 50, 60 ), model = "probit" )
urbinEffInt( coefOProbitLinNC, allXVal = xMeanLinIntNC, 
  xPos = 3, refBound = c( 30, 40 ), intBound = c( 50, 60 ), model = "probit" )
# effects of age changing from the 30-40 interval to the 50-60 interval
# based on predicted values
Mroz87Ref <- as.data.frame( t( replace( xMeanLinNC, 3, 35 ) ) )
Mroz87Int <- as.data.frame( t( replace( xMeanLinNC, 3, 55 ) ) )
effIntNum <- predict( estOProbitLin, newdata = Mroz87Int, type = "probs" ) -
  predict( estOProbitLin, newdata = Mroz87Ref, type = "probs" )
print( effIntNum )
print( sum( effIntNum[ c( "part", "full" ) ] ) )
# partial derivatives of the semi-elasticity wrt the coefficients
xMeanLinIntNVAttr <- xMeanLinIntNV
attr( xMeanLinIntNVAttr, "derivOnly" ) <- 1 
urbinEffInt( coefOProbitLinNV, xMeanLinIntNVAttr, 3,
  c( 30, 40 ), c( 50, 60 ), model = "probit" )
xMeanLinIntNCAttr <- xMeanLinIntNC
attr( xMeanLinIntNCAttr, "derivOnly" ) <- 1 
urbinEffInt( coefOProbitLinNC, xMeanLinIntNCAttr, 3,
  c( 30, 40 ), c( 50, 60 ), model = "probit" )
# numerically computed partial derivatives of the semi-elasticity wrt the coefficients
numericGradient( urbinEffInt, t0 = coefOProbitLinNV,
  allXVal = xMeanLinIntNV, xPos = 3,
  refBound = c( 30, 40 ), intBound = c( 50, 60 ), model = "probit" )
numericGradient( urbinEffInt, t0 = coefOProbitLinNC,
  allXVal = xMeanLinIntNC, xPos = 3,
  refBound = c( 30, 40 ), intBound = c( 50, 60 ), model = "probit" )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (full covariance matrix) 
urbinEffInt( coefOProbitLinNV, xMeanLinIntNV, 3,
  c( 30, 40 ), c( 50, 60 ), model = "probit", 
  allCoefVcov = vcovOProbitLinNV )
urbinEffInt( coefOProbitLinNC, xMeanLinIntNC, 3,
  c( 30, 40 ), c( 50, 60 ), model = "probit", 
  allCoefVcov = vcovOProbitLinNC )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (only standard errors) 
urbinEffInt( coefOProbitLinNV, allXVal = xMeanLinIntNV, 
  xPos = 3, refBound = c( 30, 40 ), intBound = c( 50, 60 ), model = "probit", 
  allCoefVcov = sqrt( diag( vcovOProbitLinNV ) ) )
urbinEffInt( coefOProbitLinNC, allXVal = xMeanLinIntNC, 
  xPos = 3, refBound = c( 30, 40 ), intBound = c( 50, 60 ), model = "probit", 
  allCoefVcov = sqrt( diag( vcovOProbitLinNC ) ) )


### effect of age changing between discrete intervals 
### if age is used as linear and quadratic explanatory variable 
# mean values of the 'other' explanatory variables
xMeanQuadIntNV <- c( xMeanLinNV[ 1:2 ], NA, NA, xMeanLinNV[4] )
xMeanQuadIntNC <- c( xMeanLinNC[ 1:2 ], NA, NA, xMeanLinNC[4] )
# effects of age changing from the 30-40 interval to the 50-60 interval
# without standard errors
urbinEffInt( coefOProbitQuadNV, allXVal = xMeanQuadIntNV, 
  xPos = c( 3, 4 ), refBound = c( 30, 40 ), intBound = c( 50, 60 ), 
  model = "probit" )
urbinEffInt( coefOProbitQuadNC, allXVal = xMeanQuadIntNC, 
  xPos = c( 3, 4 ), refBound = c( 30, 40 ), intBound = c( 50, 60 ), 
  model = "probit" )
# effects of age changing from the 30-40 interval to the 50-60 interval
# based on predicted values
Mroz87Ref <- as.data.frame( t( replace( xMeanQuadNC, 3:4, c( 35, 35^2 ) ) ) )
Mroz87Int <- as.data.frame( t( replace( xMeanQuadNC, 3:4, c( 55, 55^2 ) ) ) )
effIntQuadNum <- predict( estOProbitQuad, newdata = Mroz87Int, type = "probs" ) -
  predict( estOProbitQuad, newdata = Mroz87Ref, type = "probs" )
print( effIntQuadNum )
print( sum( effIntQuadNum[ c( "part", "full" ) ] ) )
# partial derivatives of the effect wrt the coefficients
xMeanQuadIntNVAttr <- xMeanQuadIntNV
attr( xMeanQuadIntNVAttr, "derivOnly" ) <- 1 
urbinEffInt( coefOProbitQuadNV, xMeanQuadIntNVAttr, 
  c( 3, 4 ), c( 30, 40 ), c( 50, 60 ), model = "probit" )
xMeanQuadIntNCAttr <- xMeanQuadIntNC
attr( xMeanQuadIntNCAttr, "derivOnly" ) <- 1 
urbinEffInt( coefOProbitQuadNC, xMeanQuadIntNCAttr, 
  c( 3, 4 ), c( 30, 40 ), c( 50, 60 ), model = "probit" )
# numerically computed partial derivatives of the effect wrt the coefficients
numericGradient( urbinEffInt, t0 = coefOProbitQuadNV,
  allXVal = xMeanQuadIntNV, xPos = c( 3, 4 ),
  refBound = c( 30, 40 ), intBound = c( 50, 60 ), model = "probit" )
numericGradient( urbinEffInt, t0 = coefOProbitQuadNC,
  allXVal = xMeanQuadIntNC, xPos = c( 3, 4 ),
  refBound = c( 30, 40 ), intBound = c( 50, 60 ), model = "probit" )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (full covariance matrix) 
urbinEffInt( coefOProbitQuadNV, xMeanQuadIntNV, 
  c( 3, 4 ), c( 30, 40 ), c( 50, 60 ), model = "probit", 
  allCoefVcov = vcovOProbitQuadNV )
urbinEffInt( coefOProbitQuadNC, xMeanQuadIntNC, 
  c( 3, 4 ), c( 30, 40 ), c( 50, 60 ), model = "probit", 
  allCoefVcov = vcovOProbitQuadNC )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (only standard errors) 
urbinEffInt( coefOProbitQuadNV, allXVal = xMeanQuadIntNV, 
  xPos = c( 3, 4 ), refBound = c( 30, 40 ), intBound = c( 50, 60 ), 
  model = "probit", sqrt( diag( vcovOProbitQuadNV ) ) )
urbinEffInt( coefOProbitQuadNC, allXVal = xMeanQuadIntNC, 
  xPos = c( 3, 4 ), refBound = c( 30, 40 ), intBound = c( 50, 60 ), 
  model = "probit", sqrt( diag( vcovOProbitQuadNC ) ) )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (standard errors + mean value and standard deviation of age)
urbinEffInt( coefOProbitQuadNV, xMeanQuadIntNV, c( 3, 4 ),
  c( 30, 40 ), c( 50, 60 ), model = "probit", 
  allCoefVcov = sqrt( diag( vcovOProbitQuadNV ) ),
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ) )
urbinEffInt( coefOProbitQuadNC, xMeanQuadIntNC, c( 3, 4 ),
  c( 30, 40 ), c( 50, 60 ), model = "probit", 
  allCoefVcov = sqrt( diag( vcovOProbitQuadNC ) ),
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ) )


### grouping and re-basing categorical variables
### effects of age changing from the 30-44 category to the 53-60 category
# without standard errors
urbinEffCat( coefOProbitIntNV, xMeanIntNV, 
  xPos = c( 3:5 ), xGroups = c( -1, -1, 1, 0 ), model = "probit" )
urbinEffCat( coefOProbitIntNC, xMeanIntNC, 
  xPos = c( 3:5 ), xGroups = c( -1, -1, 1, 0 ), model = "probit" )
# effects calculated based on predicted values
names( xMeanIntNC ) <- 
  gsub( "TRUE|full:", "", names( coefOProbitIntNC ) )
df30.37 <- df38.44 <- df45.52 <- df53.60 <- as.data.frame( t( xMeanIntNC ) ) 
df30.37[ , 3:5 ] <- c( TRUE, FALSE, FALSE )
df38.44[ , 3:5 ] <- c( FALSE, TRUE, FALSE )
df45.52[ , 3:5 ] <- c( FALSE, FALSE, FALSE )
df53.60[ , 3:5 ] <- c( FALSE, FALSE, TRUE )
effCatNum <- predict( estOProbitInt, newdata = df53.60, type = "probs" ) -
  sum( Mroz87$age30.37 ) / sum( Mroz87$age30.37 + Mroz87$age38.44 ) *
  predict( estOProbitInt, newdata = df30.37, type = "probs" ) -
  sum( Mroz87$age38.44 ) / sum( Mroz87$age30.37 + Mroz87$age38.44 ) *
  predict( estOProbitInt, newdata = df38.44, type = "probs" )
print( effCatNum )
print( sum( effCatNum[ c( "part", "full" ) ] ) )
# partial derivatives of the effect wrt the coefficients
urbinEffCat( coefOProbitIntNV, xMeanIntNVAttr, 
  c( 3:5 ), c( -1, -1, 1, 0 ), model = "probit" )
urbinEffCat( coefOProbitIntNC, xMeanIntNCAttr, 
  c( 3:5 ), c( -1, -1, 1, 0 ), model = "probit" )
# numerically computed partial derivatives of the effect wrt the coefficients
numericGradient( urbinEffCat, t0 = coefOProbitIntNV,
  allXVal = xMeanIntNV, xPos = c( 3:5 ), xGroups = c( -1, -1, 1, 0 ),
  model = "probit" )
numericGradient( urbinEffCat, t0 = coefOProbitIntNC,
  allXVal = xMeanIntNC, xPos = c( 3:5 ), xGroups = c( -1, -1, 1, 0 ),
  model = "probit" )
# with full covariance matrix
urbinEffCat( coefOProbitIntNV, xMeanIntNV, c( 3:5 ), 
  c( -1, -1, 1, 0 ), vcovOProbitIntNV, 
  model = "probit" )
urbinEffCat( coefOProbitIntNC, xMeanIntNC, c( 3:5 ), 
  c( -1, -1, 1, 0 ), vcovOProbitIntNC, 
  model = "probit" )
# with standard errors only
urbinEffCat( coefOProbitIntNV, xMeanIntNV, c( 3:5 ), 
  c( -1, -1, 1, 0 ), sqrt( diag( vcovOProbitIntNV ) ), 
  model = "probit" )
urbinEffCat( coefOProbitIntNC, xMeanIntNC, c( 3:5 ), 
  c( -1, -1, 1, 0 ), sqrt( diag( vcovOProbitIntNC ) ), 
  model = "probit" )

### effects of age changing from the 53-60 category to the 38-52 category
# without standard errors
urbinEffCat( coefOProbitIntNV, xMeanIntNV, c( 3:5 ), 
  c( 0, 1, -1, 1 ), model = "probit" )
urbinEffCat( coefOProbitIntNC, xMeanIntNC, c( 3:5 ), 
  c( 0, 1, -1, 1 ), model = "probit" )
# effects calculated based on predicted values
effCat2Num <- sum( Mroz87$age38.44 ) / sum( Mroz87$age38.44 + Mroz87$age45.52 ) *
  predict( estOProbitInt, newdata = df38.44, type = "probs" ) +
  sum( Mroz87$age45.52 ) / sum( Mroz87$age38.44 + Mroz87$age45.52 ) *
  predict( estOProbitInt, newdata = df45.52, type = "probs" ) -
  predict( estOProbitInt, newdata = df53.60, type = "probs" )
print( effCat2Num )
print( sum( effCat2Num[ c( "part", "full" ) ] ) )
# partial derivatives of the effect wrt the coefficients
urbinEffCat( coefOProbitIntNV, xMeanIntNVAttr, 
  c( 3:5 ), c( 0, 1, -1, 1 ), model = "probit" )
urbinEffCat( coefOProbitIntNC, xMeanIntNCAttr, 
  c( 3:5 ), c( 0, 1, -1, 1 ), model = "probit" )
# numerically computed partial derivatives of the effect wrt the coefficients
numericGradient( urbinEffCat, t0 = coefOProbitIntNV,
  allXVal = xMeanIntNV, xPos = c( 3:5 ), xGroups = c( 0, 1, -1, 1 ), 
  model = "probit" )
numericGradient( urbinEffCat, t0 = coefOProbitIntNC,
  allXVal = xMeanIntNC, xPos = c( 3:5 ), xGroups = c( 0, 1, -1, 1 ), 
  model = "probit" )
# with full covariance matrix
urbinEffCat( coefOProbitIntNV, xMeanIntNV, c( 3:5 ), 
  c( 0, 1, -1, 1 ), vcovOProbitIntNV, 
  model = "probit" )
urbinEffCat( coefOProbitIntNC, xMeanIntNC, c( 3:5 ), 
  c( 0, 1, -1, 1 ), vcovOProbitIntNC, 
  model = "probit" )
# with standard errors only
urbinEffCat( coefOProbitIntNV, xMeanIntNV, c( 3:5 ), 
  c( 0, 1, -1, 1 ), sqrt( diag( vcovOProbitIntNV ) ), 
  model = "probit" )
urbinEffCat( coefOProbitIntNC, xMeanIntNC, c( 3:5 ), 
  c( 0, 1, -1, 1 ), sqrt( diag( vcovOProbitIntNC ) ), 
  model = "probit" )
