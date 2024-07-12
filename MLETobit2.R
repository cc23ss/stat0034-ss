# Maximum likelihood estimation in Tobit-2 models
# sigma is reparameterised using a log link for unrestricted optimisation
# rho is reparameterised using a tanh link for unrestricted optimisation
#NOTE: This is the modified MLETobit2.R function from Example 1, which itself
#is based on Rubio (2020).  
MLETobit2 = function(init, out, sel, deso, dess, beta_par = 0,
                     method = "nlminb", maxit = 100){
  #Converting the inputs out and sel into vectors.
  out = as.vector(out); sel = as.vector(sel)
  #Defining x1 and w1 as the design matrices of the observed and selection 
  #processes respectively by adding back the missing column of 1's.
  x1 = as.matrix(cbind(1, deso)); w1 = as.matrix(cbind(1, dess))
  #Defining p and q as the number of columns of x1 and w1 respectively. 
  p = ncol(x1); q = ncol(w1)
  #Defining a function which evaluates the log-likelihood of the Tobit-2 model,
  #loglik. This takes in as inputs a parameter vector, par. 
  loglik = function(par){
    #Reparametrising sigma and rho in our model, using an exponential and a 
    #tanh link. 
    sigma = exp(par[1]); rho = tanh(par[2]) 
    #Defining our linear regression coefficients for the observed model, beta, 
    #and that for the selection model, alpha. 
    #TODO: Reparametrise beta = beta.tilde
    beta = par[c(3:(p+2))]; alpha = par[(p+3):(p+q+2)]
    #Reparametrisation if beta_par = 1
    if(beta_par){
      beta = beta/sigma
    }
    x_beta = x1%*%beta; w_alpha = w1%*%alpha
    #Defining inds0 and inda1 as indicators corresponding to the selection 
    #criterion, which can be seen in chapter 2. 
    inds0 = (sel == 0); inds1 = (sel == 1)
    #Evaluating the log-likelihood of the Tobit-2 model as given in Chapter 2. 
    ll = sum(pnorm(-w_alpha[inds0], log = T)) + 
      sum(pnorm( (w_alpha[inds1] + rho*(out[inds1] - x_beta[inds1])/sigma)/
                   sqrt(1 - rho^2), log = T)) + 
      #MODIFICATION: compressed the explicit evaluation of the normal part of
      #the observed processes into the dnorm function, as the distribution is a
      #N(x_beta, sigma^2) distributed random variable from proof (INSERT). 
      sum(dnorm(out[inds1], mean = x_beta[inds1], sd = sigma, log = T))
    #Returning the log-loss, which is the objective we seek to minimise. 
    return(-ll)
  }
  #Conditional to evaluate the optimal value using nlminb if specified
  #(by default).
  if(method == "nlminb"){ 
    OPT = nlminb(init, loglik, control = list(iter.max = maxit))
  } 
  #Conditional to evaluate the optimal value using the optim() function
  else{
    #Using optim to 
    OPT = optim(init, loglik, control = list(maxit = maxit), method = method)
  }
  #Storing the MLES in a vector, remembering to undo the links defined above. 
  MLE = c(exp(OPT$par[1]), tanh(OPT$par[2]), OPT$par[-c(1,2)])
  #Creating a vector of estimate names...
  names(MLE) = c("sigma", "rho", "intercept", paste(paste("beta_hat[,", 1:(p-1), 
    sep = "" ),"]", sep = ""), "intercept", paste(paste("alpha_hat[,", 1:(q-1), 
                                                                                                              sep = "" ),"]", sep = "")) 
  #Storing all our outputs in the list outf.
  outf = list(loglik = loglik, OPT = OPT, MLE = MLE)
  #Returning outf as the function output. 
  return(outf)
  
}