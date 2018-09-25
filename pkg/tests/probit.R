library( "urbin" )
library( "maxLik" )

# load data set
data( "Mroz87", package = "sampleSelection" )

# create dummy variable for kids
Mroz87$kids <- as.numeric( Mroz87$kids5 > 0 | Mroz87$kids618 > 0 )

### linear in age
estProbitLin <- glm( lfp ~ kids + age + educ, 
  family = binomial(link = "probit"), 
  data = Mroz87 )
summary( estProbitLin )
# mean values of the explanatory variables
xMeanLin <- c( 1, colMeans( Mroz87[ , c( "kids", "age", "educ" ) ] ) )
# semi-elasticity of age without standard errors
uProbitEla( coef( estProbitLin ), xMeanLin, xPos = 3 )
# semi-elasticity of age based on numerical derivation
100 * ( predict( estProbitLin, 
  newdata = as.data.frame( t( xMeanLin * c( 1, 1, 1.005, 1 ) ) ), 
  type = "response" ) -
    predict( estProbitLin, 
      newdata = as.data.frame( t( xMeanLin * c( 1, 1, 0.995, 1 ) ) ), 
      type = "response" ) )
# semi-elasticity of age with standard errors
uProbitEla( coef( estProbitLin ), xMeanLin, 
  sqrt( diag( vcov( estProbitLin ) ) ), 3 )
# partial derivatives of the semi-elasticity wrt the coefficients
urbin:::uProbitElaDeriv( coef( estProbitLin ), xMeanLin, xPos = 3 )
# numerically computed partial derivatives of the semi-elasticity wrt the coefficients
numericGradient( uProbitEla, t0 = coef( estProbitLin ), 
  allXVal = xMeanLin, xPos = 3 )
# simplified partial derivatives of the semi-elasticity wrt the coefficients
urbin:::uProbitElaDeriv( coef( estProbitLin ), xMeanLin, xPos = 3,
  simplified = TRUE )

### quadratic in age
estProbitQuad <- glm( lfp ~ kids + age + I(age^2) + educ, 
  family = binomial(link = "probit"), 
  data = Mroz87 )
summary( estProbitQuad )
# mean values of the explanatory variables
xMeanQuad <- c( xMeanLin[ 1:3], xMeanLin[3]^2, xMeanLin[4] )
# semi-elasticity of age without standard errors
uProbitEla( coef( estProbitQuad ), xMeanQuad, xPos = c( 3, 4 ) )
# semi-elasticity of age based on numerical derivation
100 * ( predict( estProbitQuad, 
  newdata = as.data.frame( t( xMeanQuad * c( 1, 1, 1.005, 1.005^2, 1 ) ) ), 
  type = "response" ) -
    predict( estProbitQuad, 
      newdata = as.data.frame( t( xMeanQuad * c( 1, 1, 0.995, 0.995^2, 1 ) ) ), 
      type = "response" ) )
# semi-elasticity of age with standard errors
uProbitEla( coef( estProbitQuad ), xMeanQuad, 
  sqrt( diag( vcov( estProbitQuad ) ) ), c( 3, 4 ) )
# partial derivatives of the semi-elasticity wrt the coefficients
urbin:::uProbitElaDeriv( coef( estProbitQuad ), xMeanQuad, xPos = c( 3, 4 ) )
# numerically computed partial derivatives of the semi-elasticity wrt the coefficients
numericGradient( uProbitEla, t0 = coef( estProbitQuad ), 
  allXVal = xMeanQuad, xPos = c( 3, 4 ) )
# simplified partial derivatives of the semi-elasticity wrt the coefficients
urbin:::uProbitElaDeriv( coef( estProbitQuad ), xMeanQuad, xPos = c( 3, 4 ),
  simplified = TRUE )
