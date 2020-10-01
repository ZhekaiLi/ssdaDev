% Netlab Toolbox
% Version unknown 	 10-Jan-2005
%
% c_class_bin-  Online update coeffs. for binary classification. 
% c_reg_exp-  Online update coeffs. for regression with positive exponential noise. 
% c_reg_gauss-  Online update coeffs. for regression with Gaussian noise. 
% c_reg_lapl-  Online update coeffs. for regression with Laplace noise. 
% cl_data  -  Generates two-dimensional dataset for classification. 
% cov_matern-  returns the Matern covariance function. 
% covgrad_matern-  Returns the gradient of the Matern kernel. 
% defoptions-  sets the ``gpopt'' structure to default values. 
% demogp_class-  Two-dimensional binary classification example. 
% demogp_class_gui-  Graphical frontend for Onlinde Gaussian Process classification. 
% demogp_fixed-  Gaussian Process regression using a fixed set of basis vectors. 
% demogp_matern-  Example to use an external covariance function. 
% demogp_reg-  One-dimensional regression example with different noise models. 
% demogp_reg_gui-  Graphical frontend to Online Gaussian Process regression. 
% demogp_simp-  Gaussian process regression using user inputs. 
% em_exp   -  Recomputes the likelihood parameter for pos.-exp. noise. 
% em_gauss -  Recomputes the (Gaussian) noise variance. 
% em_lapl  -  Recomputes the likelihood parameter for exponential noise. 
% err_2class-  Computes the binary classification error. 
% err_2logp-  Computes the log-predictive probability of the labels. 
% err_abs  -  Computes the absolute error. 
% err_mse  -  Computes the mean-square error. 
% g_l_gauss-  Recomputes the variance of the Gaussian noise using gradients 
% laplace  -  Sampling from a Laplace distribution. 
% logp_exp -  Computes the log-predictive probability for positive exponential noise 
% logp_g   -  Computes the log-predictive probability for Gaussian likelihood 
% logp_l   -  Computes the log-predictive probability for Gaussian noise 
% matern   -  Computes the Matern kernel 
% ogp      -  Initialises the global net structure for the OGP toolbox. 
% ogpadjgp -  Computes the GP coefficients from the TAP/EP ones. 
% ogpbvmin -  Finds the BV that contributes the least to the GP 
% ogpcovarf-  Calculate the covariance function for the OGP. 
% ogpcovarp-  Calculate the prior covariance for the Sparse GP. 
% ogpcovdiag-  Calculates the diagonal of the covariance function for the OGP. 
% ogpcovgrad-  Evaluate gradient for kernel parameters (except bias) for the Sparse GP 
% ogpdelbv -  Deletes the specified BVs from the BV set of the GP. 
% ogpemptybv-  Adds input elements to the BV set without altering the GP. 
% ogpevid  -  Evaluates the evidence for Sparse OGP. 
% ogpevidgrad-  Computes the gradient for the Sparse OGP. 
% ogpfwd   -  Forward propagation through a Sparse OGP. 
% ogphypcovpar-  Initialises the Sparse Gaussian Process hyperparameters. 
% ogpinit  -  Likelihood initialisation for the Online Gaussian Process structure. 
% ogpkl    -  Computes the KL-distance of the GP marginals. 
% ogppak   -  Puts the Sparse OGP hyperparametrs into a vector. 
% ogpparadj-  Adjusts GP such that the TAP/EP it. is not ill-conditioned. 
% ogppost  -  Calculation of the sparse posterior 
% ogpreset -  Resets the Gaussian Process. 
% ogpsample-  Generates samples from a Gaussian process. 
% ogpstep_ep-  Updates the TAP/EP parameters after an online sweep 
% ogpstep_full-  Performs a full online update step of the Gaussian Process 
% ogpstep_sp-  Performs a sparse online update step of the Gaussian Process. 
% ogptrain -  Inference for Sparse Gaussian Processes. 
% ogpunpak -  Puts hyperparameters back into the Sparse OGP. 
% sinc     -  Sin(pi*x)/(pi*x) function. (from MATLAB) 
% sinc2data-  Generates two-dimensional sinc test data. 
% sincdata -  Generates one-dimensional sinc test data. 
%
%	Copyright (c) Lehel Csato (2001-2004)
%
