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
# mean values of the explanatory variables and specification of the threshold
xMeanLin <- c( colMeans( Mroz87[ , c( "kids", "age", "educ" ) ] ), -1, 0 )
# semi-elasticity of age without standard errors
urbinEla( coef( summary( estOProbitLin ) )[,1], xMeanLin, xPos = 2, 
  iPos = 4, model = "probit" )
urbinEla( coef( summary( estOProbitLin ) )[-5,1], xMeanLin[-5], xPos = 2, 
  iPos = 4, model = "probit" )
# semi-elasticity of age based on numerical derivation
Mroz87Lower <- as.data.frame( t( xMeanLin * c( 1, 0.995, 1, 1, 1 ) ) )
Mroz87Upper <- as.data.frame( t( xMeanLin * c( 1, 1.005, 1, 1, 1 ) ) )
elaLinNum <- 100 * ( 
  predict( estOProbitLin, newdata = Mroz87Upper, type = "probs" ) -
    predict( estOProbitLin, newdata = Mroz87Lower, type = "probs" ) )
print( elaLinNum )
print( sum( elaLinNum[ c( "part", "full" ) ] ) )
# partial derivatives of the semi-elasticity wrt the coefficients
xMeanLinAttr <- xMeanLin
attr( xMeanLinAttr, "derivOnly" ) <- 1 
urbinEla( coef( summary( estOProbitLin ) )[,1], xMeanLinAttr, xPos = 2, 
  iPos = 4, seSimplify = FALSE, model = "probit" )
xMeanLin4Attr <- xMeanLin[-5]
attr( xMeanLin4Attr, "derivOnly" ) <- 1 
urbinEla( coef( summary( estOProbitLin ) )[-5,1], xMeanLin4Attr, xPos = 2, 
  iPos = 4, seSimplify = FALSE, model = "probit" )
# numerically computed partial derivatives of the semi-elasticity wrt the coefficients
numericGradient( urbinEla, t0 = coef( summary( estOProbitLin ) )[,1], 
  allXVal = xMeanLin, xPos = 2, iPos = 4, model = "probit" )
numericGradient( urbinEla, t0 = coef( summary( estOProbitLin ) )[-5,1], 
  allXVal = xMeanLin[-5], xPos = 2, iPos = 4, model = "probit" )
# simplified partial derivatives of the semi-elasticity wrt the coefficients
urbinEla( coef( summary( estOProbitLin ) )[,1], xMeanLinAttr, xPos = 2, 
  iPos = 4, model = "probit", seSimplify = TRUE )
urbinEla( coef( summary( estOProbitLin ) )[-5,1], xMeanLin4Attr, xPos = 2, 
  iPos = 4, model = "probit", seSimplify = TRUE )
# semi-elasticity of age with standard errors (full covariance matrix)
urbinEla( coef( summary( estOProbitLin ) )[,1], xMeanLin, xPos = 2, 
  iPos = 4, model = "probit", vcov( estOProbitLin ) )
urbinEla( coef( summary( estOProbitLin ) )[-5,1], xMeanLin[-5], xPos = 2, 
  iPos = 4, model = "probit", vcov( estOProbitLin )[-5,-5] )
# semi-elasticity of age with standard errors (only standard errors)
urbinEla( coef( summary( estOProbitLin ) )[,1], xMeanLin, 
  xPos = 2, iPos = 4, model = "probit", 
  sqrt( diag( vcov( estOProbitLin ) ) ), seSimplify = FALSE )
urbinEla( coef( summary( estOProbitLin ) )[-5,1], xMeanLin[-5], 
  xPos = 2, iPos = 4, model = "probit", 
  sqrt( diag( vcov( estOProbitLin ) ) )[-5], seSimplify = FALSE )
# semi-elasticity of age with standard errors (only standard errors, simplified)
urbinEla( coef( summary( estOProbitLin ) )[,1], xMeanLin, 
  xPos = 2, iPos = 4, model = "probit", 
  sqrt( diag( vcov( estOProbitLin ) ) ) )
urbinEla( coef( summary( estOProbitLin ) )[-5,1], xMeanLin[-5], 
  xPos = 2, iPos = 4, model = "probit", 
  sqrt( diag( vcov( estOProbitLin ) ) )[-5] )


### quadratic in age
estOProbitQuad <- polr( lfp3 ~ kids + age + I(age^2) + educ, 
  data = Mroz87, method = "probit", Hess = TRUE )
summary( estOProbitQuad )
# mean values of the explanatory variables and specification of the threshold
xMeanQuad <- c( xMeanLin[ 1:2 ], xMeanLin[2]^2, xMeanLin[3:5] )
# semi-elasticity of age without standard errors
urbinEla( coef( summary( estOProbitQuad ) )[,1], xMeanQuad, 
  xPos = c( 2, 3 ), iPos = 5, model = "probit" )
urbinEla( coef( summary( estOProbitQuad ) )[-6,1], xMeanQuad[-6], 
  xPos = c( 2, 3 ), iPos = 5, model = "probit" )
# semi-elasticity of age based on numerical derivation
Mroz87Lower <- as.data.frame( 
  t( xMeanQuad * c( 1, 0.995, 0.995^2, 1, 1, 1 ) ) )
Mroz87Upper <- as.data.frame( 
  t( xMeanQuad * c( 1, 1.005, 1.005^2, 1, 1, 1 ) ) )
elaQuadNum <- 100 * ( 
  predict( estOProbitQuad, newdata = Mroz87Upper, type = "probs" ) -
    predict( estOProbitQuad, newdata = Mroz87Lower, type = "probs" ) )
print( elaQuadNum )
print( sum( elaQuadNum[ c( "part", "full" ) ] ) )
# partial derivatives of the semi-elasticity wrt the coefficients
xMeanQuadAttr <- xMeanQuad
attr( xMeanQuadAttr, "derivOnly" ) <- 1 
urbinEla( coef( summary( estOProbitQuad ) )[,1], xMeanQuadAttr, 
  xPos = c( 2, 3 ), iPos = 5, model = "probit", seSimplify = FALSE )
xMeanQuad5Attr <- xMeanQuad[-6]
attr( xMeanQuad5Attr, "derivOnly" ) <- 1 
urbinEla( coef( summary( estOProbitQuad ) )[-6,1], xMeanQuad5Attr, 
  xPos = c( 2, 3 ), iPos = 5, model = "probit", seSimplify = FALSE )
# numerically computed partial derivatives of the semi-elasticity wrt the coefficients
numericGradient( urbinEla, t0 = coef( summary( estOProbitQuad ) )[,1], 
  allXVal = xMeanQuad, xPos = c( 2, 3 ), iPos = 5, model = "probit" )
numericGradient( urbinEla, t0 = coef( summary( estOProbitQuad ) )[-6,1], 
  allXVal = xMeanQuad[-6], xPos = c( 2, 3 ), iPos = 5, model = "probit" )
# simplified partial derivatives of the semi-elasticity wrt the coefficients
urbinEla( coef( summary( estOProbitQuad ) )[,1], xMeanQuadAttr, 
  xPos = c( 2, 3 ), iPos = 5, model = "probit", seSimplify = TRUE )
urbinEla( coef( summary( estOProbitQuad ) )[-6,1], xMeanQuad5Attr, 
  xPos = c( 2, 3 ), iPos = 5, model = "probit", seSimplify = TRUE )
# semi-elasticity of age with standard errors (full covariance matrix)
urbinEla( coef( summary( estOProbitQuad ) )[,1], xMeanQuad, 
  xPos = c( 2, 3 ), iPos = 5, model = "probit", 
  vcov( estOProbitQuad ) )
urbinEla( coef( summary( estOProbitQuad ) )[-6,1], xMeanQuad[-6], 
  xPos = c( 2, 3 ), iPos = 5, model = "probit", 
  vcov( estOProbitQuad )[-6,-6] )
# semi-elasticity of age with standard errors (only standard errors)
urbinEla( coef( summary( estOProbitQuad ) )[,1], xMeanQuad, 
  xPos = c( 2, 3 ), iPos = 5, model = "probit", 
  sqrt( diag( vcov( estOProbitQuad ) ) ), seSimplify = FALSE )
urbinEla( coef( summary( estOProbitQuad ) )[-6,1], xMeanQuad[-6], 
  xPos = c( 2, 3 ), iPos = 5, model = "probit", 
  sqrt( diag( vcov( estOProbitQuad ) ) )[-6], seSimplify = FALSE )
# semi-elasticity of age with standard errors (only standard errors, simplified)
urbinEla( coef( summary( estOProbitQuad ) )[,1], xMeanQuad, 
  xPos = c( 2, 3 ),iPos = 5, model = "probit", 
  sqrt( diag( vcov( estOProbitQuad ) ) ) )
urbinEla( coef( summary( estOProbitQuad ) )[-6,1], xMeanQuad[-6], 
  xPos = c( 2, 3 ),iPos = 5, model = "probit", 
  sqrt( diag( vcov( estOProbitQuad ) ) )[-6] )
# semi-elasticity of age with standard errors (only standard errors, xMeanSd)
urbinEla( coef( summary( estOProbitQuad ) )[,1], xMeanQuad, 
  xPos = c( 2, 3 ), iPos = 5, model = "probit", 
  sqrt( diag( vcov( estOProbitQuad ) ) ),
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ),
  seSimplify = FALSE )
urbinEla( coef( summary( estOProbitQuad ) )[-6,1], xMeanQuad[-6], 
  xPos = c( 2, 3 ), iPos = 5, model = "probit", 
  sqrt( diag( vcov( estOProbitQuad ) ) )[-6],
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ),
  seSimplify = FALSE )
# semi-elasticity of age with standard errors (only standard errors, xMeanSd, simplified)
urbinEla( coef( summary( estOProbitQuad ) )[,1], xMeanQuad, 
  xPos = c( 2, 3 ), iPos = 5, model = "probit", 
  sqrt( diag( vcov( estOProbitQuad ) ) ),
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ) )
urbinEla( coef( summary( estOProbitQuad ) )[-6,1], xMeanQuad[-6], 
  xPos = c( 2, 3 ), iPos = 5, model = "probit", 
  sqrt( diag( vcov( estOProbitQuad ) ) )[-6],
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
# mean values of the explanatory variables and specification of the threshold
xMeanInt <- c( xMeanLin[1], mean( Mroz87$age30.37 ), 
  mean( Mroz87$age38.44 ), mean( Mroz87$age53.60 ), xMeanLin[3:5] )
# semi-elasticity of age without standard errors
urbinElaInt( coef( summary( estOProbitInt ) )[,1], xMeanInt,
  c( 2, 3, 0, 4 ), iPos = 6, c( 30, 37.5, 44.5, 52.5, 60 ), model = "probit" )
urbinElaInt( coef( summary( estOProbitInt ) )[-7,1], xMeanInt[-7],
  c( 2, 3, 0, 4 ), iPos = 6, c( 30, 37.5, 44.5, 52.5, 60 ), model = "probit" )
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
xMeanIntAttr <- xMeanInt
attr( xMeanIntAttr, "derivOnly" ) <- 1 
urbinElaInt( coef( summary( estOProbitInt ) )[,1], xMeanIntAttr,
  c( 2, 3, 0, 4 ), iPos = 6, c( 30, 37.5, 44.5, 52.5, 60 ), model = "probit" )
xMeanInt6Attr <- xMeanInt[-7]
attr( xMeanInt6Attr, "derivOnly" ) <- 1 
urbinElaInt( coef( summary( estOProbitInt ) )[-7,1], xMeanInt6Attr,
  c( 2, 3, 0, 4 ), iPos = 6, c( 30, 37.5, 44.5, 52.5, 60 ), model = "probit" )
# numerically computed partial derivatives of the semi-elasticity wrt the coefficients
numericGradient( urbinElaInt, t0 = coef( summary( estOProbitInt ) )[,1], 
  allXVal = xMeanInt, xPos = c( 2, 3, 0, 4 ), iPos = 6, 
  xBound = c( 30, 37.5, 44.5, 52.5, 60 ), model = "probit" )
numericGradient( urbinElaInt, t0 = coef( summary( estOProbitInt ) )[-7,1], 
  allXVal = xMeanInt[-7], xPos = c( 2, 3, 0, 4 ), iPos = 6, 
  xBound = c( 30, 37.5, 44.5, 52.5, 60 ), model = "probit" )
# semi-elasticity of age with standard errors (full covariance matrix)
urbinElaInt( coef( summary( estOProbitInt ) )[,1], xMeanInt,
  c( 2, 3, 0, 4 ), iPos = 6, c( 30, 37.5, 44.5, 52.5, 60 ), model = "probit",
  allCoefVcov = vcov( estOProbitInt ) )
urbinElaInt( coef( summary( estOProbitInt ) )[-7,1], xMeanInt[-7],
  c( 2, 3, 0, 4 ), iPos = 6, c( 30, 37.5, 44.5, 52.5, 60 ), model = "probit",
  allCoefVcov = vcov( estOProbitInt )[-7,-7] )
# semi-elasticity of age with standard errors (only standard errors)
urbinElaInt( coef( summary( estOProbitInt ) )[,1], xMeanInt,
  c( 2, 3, 0, 4 ), iPos = 6, c( 30, 37.5, 44.5, 52.5, 60 ), model = "probit",
  allCoefVcov = sqrt( diag( vcov( estOProbitInt ) ) ) )
urbinElaInt( coef( summary( estOProbitInt ) )[-7,1], xMeanInt[-7],
  c( 2, 3, 0, 4 ), iPos = 6, c( 30, 37.5, 44.5, 52.5, 60 ), model = "probit",
  allCoefVcov = sqrt( diag( vcov( estOProbitInt ) ) )[-7] )


### effect of age changing between discrete intervals 
### if age is used as linear explanatory variable 
# mean values of the 'other' explanatory variables
xMeanLinInt <- c( xMeanLin[ 1 ], NA, xMeanLin[3:5] )
# effects of age changing from the 30-40 interval to the 50-60 interval
# without standard errors
urbinEffInt( coef( summary( estOProbitLin ) )[,1], allXVal = xMeanLinInt, 
  xPos = 2, iPos = 4, refBound = c( 30, 40 ), intBound = c( 50, 60 ), 
  model = "probit" )
urbinEffInt( coef( summary( estOProbitLin ) )[-5,1], allXVal = xMeanLinInt[-5], 
  xPos = 2, iPos = 4, refBound = c( 30, 40 ), intBound = c( 50, 60 ), 
  model = "probit" )
# effects of age changing from the 30-40 interval to the 50-60 interval
# based on predicted values
Mroz87Ref <- as.data.frame( t( replace( xMeanLin, 2, 35 ) ) )
Mroz87Int <- as.data.frame( t( replace( xMeanLin, 2, 55 ) ) )
effIntNum <- predict( estOProbitLin, newdata = Mroz87Int, type = "probs" ) -
  predict( estOProbitLin, newdata = Mroz87Ref, type = "probs" )
print( effIntNum )
print( sum( effIntNum[ c( "part", "full" ) ] ) )
# partial derivatives of the semi-elasticity wrt the coefficients
xMeanLinIntAttr <- xMeanLinInt
attr( xMeanLinIntAttr, "derivOnly" ) <- 1 
urbinEffInt( coef( summary( estOProbitLin ) )[,1], xMeanLinIntAttr, 
  xPos = 2, iPos = 4, c( 30, 40 ), c( 50, 60 ), model = "probit" )
xMeanLinInt4Attr <- xMeanLinInt[-5]
attr( xMeanLinInt4Attr, "derivOnly" ) <- 1 
urbinEffInt( coef( summary( estOProbitLin ) )[-5,1], xMeanLinInt4Attr, 
  xPos = 2, iPos = 4, c( 30, 40 ), c( 50, 60 ), model = "probit" )
# numerically computed partial derivatives of the semi-elasticity wrt the coefficients
numericGradient( urbinEffInt, t0 = coef( summary( estOProbitLin ) )[,1],
  allXVal = xMeanLinInt, xPos = 2, iPos = 4,
  refBound = c( 30, 40 ), intBound = c( 50, 60 ), model = "probit" )
numericGradient( urbinEffInt, t0 = coef( summary( estOProbitLin ) )[-5,1],
  allXVal = xMeanLinInt[-5], xPos = 2, iPos = 4,
  refBound = c( 30, 40 ), intBound = c( 50, 60 ), model = "probit" )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (full covariance matrix) 
urbinEffInt( coef( summary( estOProbitLin ) )[,1], xMeanLinInt, 
  xPos = 2, iPos = 4, c( 30, 40 ), c( 50, 60 ), model = "probit", 
  allCoefVcov = vcov( estOProbitLin ) )
urbinEffInt( coef( summary( estOProbitLin ) )[-5,1], xMeanLinInt[-5], 
  xPos = 2, iPos = 4, c( 30, 40 ), c( 50, 60 ), model = "probit", 
  allCoefVcov = vcov( estOProbitLin )[-5,-5] )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (only standard errors) 
urbinEffInt( coef( summary( estOProbitLin ) )[,1], allXVal = xMeanLinInt, 
  xPos = 2, iPos = 4, refBound = c( 30, 40 ), intBound = c( 50, 60 ), 
  model = "probit", allCoefVcov = sqrt( diag( vcov( estOProbitLin ) ) ) )
urbinEffInt( coef( summary( estOProbitLin ) )[-5,1], allXVal = xMeanLinInt[-5], 
  xPos = 2, iPos = 4, refBound = c( 30, 40 ), intBound = c( 50, 60 ), 
  model = "probit", allCoefVcov = sqrt( diag( vcov( estOProbitLin ) ) )[-5] )


### effect of age changing between discrete intervals 
### if age is used as linear and quadratic explanatory variable 
# mean values of the 'other' explanatory variables
xMeanQuadInt <- c( xMeanLin[1], NA, NA, xMeanLin[3:5] )
# effects of age changing from the 30-40 interval to the 50-60 interval
# without standard errors
urbinEffInt( coef( summary( estOProbitQuad ) )[,1], allXVal = xMeanQuadInt, 
  xPos = c( 2, 3 ), iPos = 5, refBound = c( 30, 40 ), intBound = c( 50, 60 ), 
  model = "probit" )
urbinEffInt( coef( summary( estOProbitQuad ) )[-6,1], allXVal = xMeanQuadInt[-6], 
  xPos = c( 2, 3 ), iPos = 5, refBound = c( 30, 40 ), intBound = c( 50, 60 ), 
  model = "probit" )
# effects of age changing from the 30-40 interval to the 50-60 interval
# based on predicted values
Mroz87Ref <- as.data.frame( t( replace( xMeanQuad, 2:3, c( 35, 35^2 ) ) ) )
Mroz87Int <- as.data.frame( t( replace( xMeanQuad, 2:3, c( 55, 55^2 ) ) ) )
effIntQuadNum <- predict( estOProbitQuad, newdata = Mroz87Int, type = "probs" ) -
  predict( estOProbitQuad, newdata = Mroz87Ref, type = "probs" )
print( effIntQuadNum )
print( sum( effIntQuadNum[ c( "part", "full" ) ] ) )
# partial derivatives of the effect wrt the coefficients
xMeanQuadIntAttr <- xMeanQuadInt
attr( xMeanQuadIntAttr, "derivOnly" ) <- 1 
urbinEffInt( coef( summary( estOProbitQuad ) )[,1], xMeanQuadIntAttr, 
  xPos = c( 2, 3 ), iPos = 5, c( 30, 40 ), c( 50, 60 ), model = "probit" )
xMeanQuadInt5Attr <- xMeanQuadInt[-6]
attr( xMeanQuadInt5Attr, "derivOnly" ) <- 1 
urbinEffInt( coef( summary( estOProbitQuad ) )[-6,1], xMeanQuadInt5Attr, 
  xPos = c( 2, 3 ), iPos = 5, c( 30, 40 ), c( 50, 60 ), model = "probit" )
# numerically computed partial derivatives of the effect wrt the coefficients
numericGradient( urbinEffInt, t0 = coef( summary( estOProbitQuad ) )[,1],
  allXVal = xMeanQuadInt, xPos = c( 2, 3 ), iPos = 5,
  refBound = c( 30, 40 ), intBound = c( 50, 60 ), model = "probit" )
numericGradient( urbinEffInt, t0 = coef( summary( estOProbitQuad ) )[-6,1],
  allXVal = xMeanQuadInt[-6], xPos = c( 2, 3 ), iPos = 5,
  refBound = c( 30, 40 ), intBound = c( 50, 60 ), model = "probit" )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (full covariance matrix) 
urbinEffInt( coef( summary( estOProbitQuad ) )[,1], xMeanQuadInt, 
  xPos = c( 2, 3 ), iPos = 5, c( 30, 40 ), c( 50, 60 ), model = "probit", 
  allCoefVcov = vcov( estOProbitQuad ) )
urbinEffInt( coef( summary( estOProbitQuad ) )[-6,1], xMeanQuadInt[-6], 
  xPos = c( 2, 3 ), iPos = 5, c( 30, 40 ), c( 50, 60 ), model = "probit", 
  allCoefVcov = vcov( estOProbitQuad )[-6,-6] )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (only standard errors) 
urbinEffInt( coef( summary( estOProbitQuad ) )[,1], allXVal = xMeanQuadInt, 
  xPos = c( 2, 3 ), iPos = 5, refBound = c( 30, 40 ), intBound = c( 50, 60 ), 
  model = "probit", sqrt( diag( vcov( estOProbitQuad ) ) ) )
urbinEffInt( coef( summary( estOProbitQuad ) )[-6,1], allXVal = xMeanQuadInt[-6], 
  xPos = c( 2, 3 ), iPos = 5, refBound = c( 30, 40 ), intBound = c( 50, 60 ), 
  model = "probit", sqrt( diag( vcov( estOProbitQuad ) ) )[-6] )
# effects of age changing from the 30-40 interval to the 50-60 interval
# (standard errors + mean value and standard deviation of age)
urbinEffInt( coef( summary( estOProbitQuad ) )[,1], xMeanQuadInt, 
  xPos = c( 2, 3 ), iPos = 5, c( 30, 40 ), c( 50, 60 ), model = "probit", 
  allCoefVcov = sqrt( diag( vcov( estOProbitQuad ) ) ),
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ) )
urbinEffInt( coef( summary( estOProbitQuad ) )[-6,1], xMeanQuadInt[-6], 
  xPos = c( 2, 3 ), iPos = 5, c( 30, 40 ), c( 50, 60 ), model = "probit", 
  allCoefVcov = sqrt( diag( vcov( estOProbitQuad ) ) )[-6],
  xMeanSd = c( mean( Mroz87$age ), sd( Mroz87$age ) ) )


### grouping and re-basing categorical variables
### effects of age changing from the 30-44 category to the 53-60 category
# without standard errors
urbinEffCat( coef( summary( estOProbitInt ) )[,1], xMeanInt, 
  xPos = c( 2:4 ), iPos = 6, xGroups = c( -1, -1, 1, 0 ), model = "probit" )
urbinEffCat( coef( summary( estOProbitInt ) )[-7,1], xMeanInt[-7], 
  xPos = c( 2:4 ), iPos = 6, xGroups = c( -1, -1, 1, 0 ), model = "probit" )
# effects calculated based on predicted values
names( xMeanInt ) <- 
  gsub( "TRUE|full:", "", rownames( coef( summary( estOProbitInt ) ) ) )
df30.37 <- df38.44 <- df45.52 <- df53.60 <- as.data.frame( t( xMeanInt ) ) 
df30.37[ , 2:4 ] <- c( TRUE, FALSE, FALSE )
df38.44[ , 2:4 ] <- c( FALSE, TRUE, FALSE )
df45.52[ , 2:4 ] <- c( FALSE, FALSE, FALSE )
df53.60[ , 2:4 ] <- c( FALSE, FALSE, TRUE )
effCatNum <- predict( estOProbitInt, newdata = df53.60, type = "probs" ) -
  sum( Mroz87$age30.37 ) / sum( Mroz87$age30.37 + Mroz87$age38.44 ) *
  predict( estOProbitInt, newdata = df30.37, type = "probs" ) -
  sum( Mroz87$age38.44 ) / sum( Mroz87$age30.37 + Mroz87$age38.44 ) *
  predict( estOProbitInt, newdata = df38.44, type = "probs" )
print( effCatNum )
print( sum( effCatNum[ c( "part", "full" ) ] ) )
# partial derivatives of the effect wrt the coefficients
urbinEffCat( coef( summary( estOProbitInt ) )[,1], xMeanIntAttr, 
  c( 2:4 ), iPos = 6, c( -1, -1, 1, 0 ), model = "probit" )
urbinEffCat( coef( summary( estOProbitInt ) )[-7,1], xMeanInt6Attr, 
  c( 2:4 ), iPos = 6, c( -1, -1, 1, 0 ), model = "probit" )
# numerically computed partial derivatives of the effect wrt the coefficients
numericGradient( urbinEffCat, t0 = coef( summary( estOProbitInt ) )[,1],
  allXVal = xMeanInt, xPos = c( 2:4 ), iPos = 6, xGroups = c( -1, -1, 1, 0 ),
  model = "probit" )
numericGradient( urbinEffCat, t0 = coef( summary( estOProbitInt ) )[-7,1],
  allXVal = xMeanInt[-7], xPos = c( 2:4 ), iPos = 6, xGroups = c( -1, -1, 1, 0 ),
  model = "probit" )
# with full covariance matrix
urbinEffCat( coef( summary( estOProbitInt ) )[,1], xMeanInt, c( 2:4 ), 
  iPos = 6, c( -1, -1, 1, 0 ), vcov( estOProbitInt ), 
  model = "probit" )
urbinEffCat( coef( summary( estOProbitInt ) )[-7,1], xMeanInt[-7], c( 2:4 ), 
  iPos = 6, c( -1, -1, 1, 0 ), vcov( estOProbitInt )[-7,-7], 
  model = "probit" )
# with standard errors only
urbinEffCat( coef( summary( estOProbitInt ) )[,1], xMeanInt, c( 2:4 ), 
  iPos = 6, c( -1, -1, 1, 0 ), sqrt( diag( vcov( estOProbitInt ) ) ), 
  model = "probit" )
urbinEffCat( coef( summary( estOProbitInt ) )[-7,1], xMeanInt[-7], c( 2:4 ), 
  iPos = 6, c( -1, -1, 1, 0 ), sqrt( diag( vcov( estOProbitInt ) ) )[-7], 
  model = "probit" )

### effects of age changing from the 53-60 category to the 38-52 category
# without standard errors
urbinEffCat( coef( summary( estOProbitInt ) )[,1], xMeanInt, c( 2:4 ), 
  iPos = 6, c( 0, 1, -1, 1 ), model = "probit" )
urbinEffCat( coef( summary( estOProbitInt ) )[-7,1], xMeanInt[-7], c( 2:4 ), 
  iPos = 6, c( 0, 1, -1, 1 ), model = "probit" )
# effects calculated based on predicted values
effCat2Num <- sum( Mroz87$age38.44 ) / sum( Mroz87$age38.44 + Mroz87$age45.52 ) *
  predict( estOProbitInt, newdata = df38.44, type = "probs" ) +
  sum( Mroz87$age45.52 ) / sum( Mroz87$age38.44 + Mroz87$age45.52 ) *
  predict( estOProbitInt, newdata = df45.52, type = "probs" ) -
  predict( estOProbitInt, newdata = df53.60, type = "probs" )
print( effCat2Num )
print( sum( effCat2Num[ c( "part", "full" ) ] ) )
# partial derivatives of the effect wrt the coefficients
urbinEffCat( coef( summary( estOProbitInt ) )[,1], xMeanIntAttr, 
  c( 2:4 ), iPos = 6, c( 0, 1, -1, 1 ), model = "probit" )
urbinEffCat( coef( summary( estOProbitInt ) )[-7,1], xMeanInt6Attr, 
  c( 2:4 ), iPos = 6, c( 0, 1, -1, 1 ), model = "probit" )
# numerically computed partial derivatives of the effect wrt the coefficients
numericGradient( urbinEffCat, t0 = coef( summary( estOProbitInt ) )[,1],
  allXVal = xMeanInt, xPos = c( 2:4 ), iPos = 6, xGroups = c( 0, 1, -1, 1 ), 
  model = "probit" )
numericGradient( urbinEffCat, t0 = coef( summary( estOProbitInt ) )[-7,1],
  allXVal = xMeanInt[-7], xPos = c( 2:4 ), iPos = 6, xGroups = c( 0, 1, -1, 1 ), 
  model = "probit" )
# with full covariance matrix
urbinEffCat( coef( summary( estOProbitInt ) )[,1], xMeanInt, c( 2:4 ), 
  iPos = 6, c( 0, 1, -1, 1 ), vcov( estOProbitInt ), 
  model = "probit" )
urbinEffCat( coef( summary( estOProbitInt ) )[-7,1], xMeanInt[-7], c( 2:4 ), 
  iPos = 6, c( 0, 1, -1, 1 ), vcov( estOProbitInt )[-7,-7], 
  model = "probit" )
# with standard errors only
urbinEffCat( coef( summary( estOProbitInt ) )[,1], xMeanInt, c( 2:4 ), 
  iPos = 6, c( 0, 1, -1, 1 ), sqrt( diag( vcov( estOProbitInt ) ) ), 
  model = "probit" )
urbinEffCat( coef( summary( estOProbitInt ) )[-7,1], xMeanInt[-7], c( 2:4 ), 
  iPos = 6, c( 0, 1, -1, 1 ), sqrt( diag( vcov( estOProbitInt ) ) )[-7], 
  model = "probit" )