# GSOS
Graphical selector operator with shrinkage (GSOS) simulation codes

In this folder, one can find the R codes we need in the simulation study:


ROPE.R: simulation studies based on 5-fold cross-validation for ROPE


GLASSO.R: simulation studies based on 5-fold cross-validation for glasso


GELNET.R: simulation studies based on 5-fold cross-validation for gelnet


GSOS.R: simulation studies based on 5-fold cross-validation for GSOS


############################################################################


The Figure 1-4 folder includes all the results presented in Figure 1-4 in our main document.



The Sigma&Theta folder includes all Sigma and Theta matrices of different networks in our simulation study.



The gsos folder includes R and Fortran codes for our proposed estimator; we used the "GLassoElnetFast" package to write them.


#############################################################################


Require R packages:

MASS

mvtnorm

glasso

GLassoElnetFast
