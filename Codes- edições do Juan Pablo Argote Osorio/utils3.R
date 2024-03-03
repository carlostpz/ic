transform <- function(data_mat
){
  D = max(data_mat[, 2])
  w = rep( list(NULL), D)
  for (d in 1:D){
    doc_d = data_mat[which(data_mat$doc == d), c(1,3)]
    print(paste('d = ', d))
    for (i in 1:nrow(doc_d) ){
      w[[d]] = c(w[[d]], rep( doc_d$word[i] - 1, doc_d$freq[i]) )
    } 
  }
  return(w)
}


DIC_LDA = function(data_mat, chain, mcmc_ind, output_file = FALSE){
    
  # extracts marginal chains
  beta_chain = chain$beta
  #chain$beta = 0
  theta_chain = chain$theta
  #chain$theta = 0
  
  w = transform(data)
  M = length(mcmc_ind)
  D = max(data_mat[,"doc"])
  
  for (d in 1:D){
    w[[d]] = w[[d]] + 1
  }
  
  print("Calculating mean D(theta) ...")
  
  mean_D_theta = 0
  for (m in mcmc_ind){
    print(m)
    for (d in 1:D){
	w_aux=na.omit(w[[d]])
      Nd = length(w_aux)
      for( n in 1:Nd){
        mean_D_theta = mean_D_theta + log(beta_chain[,w_aux[n],m] %*% theta_chain[d,,m])
      }
    }
  }
  mean_D_theta = -2*mean_D_theta/M
  
  beta_barra = apply(X = beta_chain,FUN = mean, MARGIN = c(1,2))
  theta_barra = apply(X = theta_chain,FUN = mean, MARGIN = c(1,2))
  
  print("Calculating D( mean(theta) ) ...")
  D_theta_mean = 0
  for (d in 1:D){
	w_aux=na.omit(w[[d]])
    Nd = length(w_aux)
    for (n in 1:Nd){
      D_theta_mean = D_theta_mean + log(beta_barra[,w_aux[n]]%*%theta_barra[d,])
    }
  }
  D_theta_mean = -2 * D_theta_mean
  
  p_D = mean_D_theta - D_theta_mean
  
  DIC = p_D + mean_D_theta
  
  print("Done.")
  if (output_file != FALSE){
    write.table(DIC, file = output_file)
  }
  return( DIC )
}

###Watanabe–Akaike Information Criterion (WAIC):

WAIC_LDA = function(data_mat, chain, mcmc_ind, output_file = FALSE){
    
  # extracts marginal chains
  beta_chain = chain$beta
  #chain$beta = 0
  theta_chain = chain$theta
  #chain$theta = 0
  
  w = transform(data)
  M = length(mcmc_ind)
  D = max(data_mat[,"doc"])
  
  for (d in 1:D){
    w[[d]] = w[[d]] + 1
  }

##Cálculo da log-verossimilhanca:

nw=0
for (d in 1:D){
w_aux=na.omit(w[[d]])
Nd = length(w_aux)
nw=nw+Nd
}

LL=matrix(0,nw,M)

  for (m in mcmc_ind){
	i=0
    print(m)
    for (d in 1:D){
	w_aux=na.omit(w[[d]])
      Nd = length(w_aux)
      for( n in 1:Nd){
	i=i+1
        LL[i,M] = log(beta_chain[,w_aux[n],m] %*% theta_chain[d,,m])
      }
    }
  }

##Cálculo do WAIC:

WAIC(LL)

}

################################################

WAIC2_LDA = function(data_mat, chain, mcmc_ind, output_file = FALSE){
    
  # extracts marginal chains
  beta_chain = chain$beta
  #chain$beta = 0
  theta_chain = chain$theta
  #chain$theta = 0
  
  w = transform(data)
  M = length(mcmc_ind)
  D = max(data_mat[,"doc"])
  
  for (d in 1:D){
    w[[d]] = w[[d]] + 1
  }

##Cálculo da log-verossimilhanca:

nw=0
for (d in 1:D){
w_aux=na.omit(w[[d]])
Nd = length(w_aux)
nw=nw+Nd
}

LL=matrix(0,nw,M)

  for (m in mcmc_ind){
	i=0
    print(m)
    for (d in 1:D){
	w_aux=na.omit(w[[d]])
      Nd = length(w_aux)
      for( n in 1:Nd){
	i=i+1
        LL[i,M] = log(beta_chain[,w_aux[n],m] %*% theta_chain[d,,m])
      }
    }
  }

LL = apply(X = LL,FUN = sum, MARGIN = 2)

##Cálculo do WAIC:

return(list(WAIC(LL),LL))

}


###AIC

AIC_LDA = function(data_mat, chain, mcmc_ind){
    
  # extracts marginal chains
  beta_chain = chain$beta
  #chain$beta = 0
  theta_chain = chain$theta
  #chain$theta = 0
  
  w = transform(data)
  M = length(mcmc_ind)
  D = max(data_mat[,"doc"])
  
  for (d in 1:D){
    w[[d]] = w[[d]] + 1
  }
  
  beta_barra = apply(X = beta_chain,FUN = mean, MARGIN = c(1,2))
  theta_barra = apply(X = theta_chain,FUN = mean, MARGIN = c(1,2))
  
LL=0

  for (d in 1:D){
	w_aux=na.omit(w[[d]])
    Nd = length(w_aux)
    for (n in 1:Nd){
      LL = LL + log(beta_barra[,w_aux[n]]%*%theta_barra[d,])
    }
  }

p = dim(beta_chain)[1]*dim(beta_chain)[2] + dim(theta_chain)[1]*dim(theta_chain)[2] #Número de parâmetros
aic_calc = 2*p-2*LL

return(list( AIC = aic_calc , beta_hat = beta_barra , theta_hat = theta_barra , Log_Verossimilhança = LL , Número_de_Parâmetros = p ))

}

###BIC

BIC_LDA = function(data_mat, chain, mcmc_ind){
    
  # extracts marginal chains
  beta_chain = chain$beta
  #chain$beta = 0
  theta_chain = chain$theta
  #chain$theta = 0
  
  w = transform(data)
  M = length(mcmc_ind)
  D = max(data_mat[,"doc"])
  
  for (d in 1:D){
    w[[d]] = w[[d]] + 1
  }

nw=0
for (d in 1:D){
w_aux=na.omit(w[[d]])
Nd = length(w_aux)
nw=nw+Nd
}

  beta_barra = apply(X = beta_chain,FUN = mean, MARGIN = c(1,2))
  theta_barra = apply(X = theta_chain,FUN = mean, MARGIN = c(1,2))
  
LL=0

  for (d in 1:D){
    w_aux=na.omit(w[[d]])
    Nd = length(w_aux)
    for (n in 1:Nd){
      LL = LL + log(beta_barra[,w_aux[n]]%*%theta_barra[d,])
    }
  }

p = dim(beta_chain)[1]*dim(beta_chain)[2] + dim(theta_chain)[1]*dim(theta_chain)[2] #Número de parâmetros
bic_calc = p*log(nw)-2*LL

return(list( BIC = bic_calc , beta_hat = beta_barra , theta_hat = theta_barra , Log_Verossimilhança = LL , Número_de_Parâmetros = p , Número_de_Observações = nw ))

}