###CDGHMM EM code with missing data + drop out

library(MASS)
library(mvtnorm)
library(ramify)
library(cluster)


#CDGHMM Models

#EEA model with lists -----------------------------------------------------
EEA_model = function(p, N, m, u, S, it, D.hat){
  #m is the number of states
  #S is the covariance matrix
  T.hat = diag(1,p)

  if (it == 2){
    #d = diag(1,p)
    d=replicate(1, diag(1,p), simplify=F)
    D.hat=replicate(m, matrix(ramify::eye(p),nrow=p,ncol=p), simplify=F)
  }
  else{
    d = D.hat
    D.hat = replicate(m, matrix(0,nrow=p,ncol=p), simplify=F)
  }

  S.hat = replicate(m, matrix(0,nrow=p,ncol=p), simplify=F)
  pi_g = rowSums(u)/N

  kappa.list = list()
  for (r in 2:p){
    if (r == 2){
      phi.hat = 0
      phi.num=0
      phi.doom=0
      ##g loop here
      for(g in 1:m){
        phi.num=phi.num+pi_g[g]*S[[g]][2,1]/d[[1]][2,2]
        phi.doom=phi.doom+pi_g[g]*S[[g]][1,1]/d[[1]][2,2]
      }
      phi.hat=-phi.num/phi.doom
      if (is.nan(phi.hat) == T){
        phi.hat = 0
      }
      T.hat[2,1] = phi.hat
    }
    else{
      kappa = matrix(0, r, r)
      for (i in 1:r){
        for (j in 1:r){
          for(g in 1:m){
            kappa[i,j]=kappa[i,j]+pi_g[g]*S[[g]][i,j]/d[[1]][r,r]
            if (is.nan(kappa[i,j]) == T){
              kappa[i,j] = 0
            }
          }
        }
      }
      kappa=-1*kappa
      phi.hat = matrix(NA, nrow = 1, ncol = r-1)
      phi.hat = t((-1)*ginv(kappa[1:r-1,1:r-1])%*%kappa[r,1:(r-1)])
      #T update
      T.hat[r,1:(r-1)] = phi.hat
    }

  }
  D.vec=rep(0,p) #1 x p
  for (g in 1:m){
    #D update
    D.hat[[g]]= pi_g[g]*(T.hat%*%S[[g]]%*%t(T.hat))
    D.hat[[g]]= diag(diag(D.hat[[g]]), p, p)
    D.vec[1:p]=D.vec[1:p]+rowSums(D.hat[[g]])
  }
  D.hat1 = replicate(1,diag(D.vec), simplify=F)

  for (g in 1:m){
    #S update
    S.hat[[g]] = t(T.hat)%*%ginv(D.hat1[[1]])%*%T.hat
    S.hat[[g]] = ginv(S.hat[[g]])
  }
  return(list(S.hat=S.hat, D.hat=D.hat1))
}


#VVA model with lists -----------------------------------------------------
VVA_model = function(p, m, S){
  #m is the number of states
  #S is the covariance matrix

  T.hat=replicate(m, matrix(0,nrow=p,ncol=p), simplify=F)
  D.hat=replicate(m, matrix(0,nrow=p,ncol=p), simplify=F)
  S.hat=replicate(m, matrix(0,nrow=p,ncol=p), simplify=F)

  for (g in 1:m){
    T.hat[[g]] = diag(1,p)
    for (r in 2:p){
      if (r == 2){
        phi.hat = c()
        phi.hat[g] = -S[[g]][2,1]/S[[g]][1,1]
        T.hat[[g]][2,1] = phi.hat[g]
      }
      else{
        phi.hat = matrix(NA, nrow = 1, ncol = r-1)
        phi.hat = t((-1)*ginv(S[[g]][1:r-1,1:r-1])%*%S[[g]][r,1:(r-1)])

        #T update
        T.hat[[g]][r,1:(r-1)] = phi.hat
      }
    }
    #D update
    D.hat[[g]]= T.hat[[g]]%*%S[[g]]%*%t(T.hat[[g]])
    D.hat[[g]] = diag(diag(D.hat[[g]]), p, p)
    #S update
    S.hat[[g]] = t(T.hat[[g]])%*%ginv(D.hat[[g]])%*%T.hat[[g]]
    S.hat[[g]] = ginv(S.hat[[g]])
  }
  S.hat
}


#VVI model with lists -----------------------------------------------------
VVI_model = function(p, m, S){
  #m is the number of states
  #S is the covariance matrix

  T.hat=replicate(m, matrix(0,nrow=p,ncol=p), simplify=F)
  D.hat=replicate(m, matrix(0,nrow=p,ncol=p), simplify=F)
  S.hat=replicate(m, matrix(0,nrow=p,ncol=p), simplify=F)
  iso_delta = c() #matrix(0.0, 1, m)

  for (g in 1:m){
    T.hat[[g]] = diag(1,p)
    for (r in 2:p){
      if (r == 2){
        phi.hat = c()
        phi.hat[g] = -S[[g]][2,1]/S[[g]][1,1]
        T.hat[[g]][2,1] = phi.hat[g]
      }
      else{
        phi.hat = matrix(NA, nrow = 1, ncol = r-1)
        phi.hat = t((-1)*ginv(S[[g]][1:r-1,1:r-1])%*%S[[g]][r,1:(r-1)])
        #T update
        T.hat[[g]][r,1:(r-1)] = phi.hat
      }
    }

    #D update
    iso_delta[g] = 1/p*tr(T.hat[[g]]%*%S[[g]]%*%t(T.hat[[g]]))
    D.hat[[g]] = iso_delta[g]*ramify::eye(p,p)

    #S update
    S.hat[[g]] = t(T.hat[[g]])%*%ginv(D.hat[[g]])%*%T.hat[[g]]
    S.hat[[g]] = ginv(S.hat[[g]])
  }
  S.hat
}


#VEA model lists ----------------------------------------------------------
VEA_model = function(p, N, m, u, S){
  # m is the number of states
  # S is the covariance matrix

  T.hat=replicate(m, matrix(0,nrow=p,ncol=p), simplify=F)
  D.hat=replicate(m, matrix(0,nrow=p,ncol=p), simplify=F)
  S.hat=replicate(m, matrix(0,nrow=p,ncol=p), simplify=F)
  pi_g = rowSums(u)/N

  for (g in 1:m){
    #T update
    T.hat[[g]] = diag(1,p)
    for (r in 2:p){
      if (r == 2){
        phi.hat = c()
        phi.hat[g] = -S[[g]][2,1]/S[[g]][1,1]
        T.hat[[g]][2,1] = phi.hat[g]
      }
      else{
        phi.hat = matrix(NA, nrow = 1, ncol = r-1)
        phi.hat = t((-1)*ginv(S[[g]][1:r-1,1:r-1])%*%S[[g]][r,1:(r-1)])
        #T update
        T.hat[[g]][r,1:(r-1)] = phi.hat
      }
    }
    #D update
    D.hat[[g]]= T.hat[[g]]%*%S[[g]]%*%t(T.hat[[g]])
    D.hat[[g]]= diag(diag(D.hat[[g]]), p, p)
    D.hat[[g]]= pi_g[g] * D.hat[[g]]
  }

  D.vec=rep(0,p)
  for (g in 1:m){
    D.vec[1:p]=D.vec[1:p]+rowSums(D.hat[[g]])
  }
  D.hat1 = diag(D.vec)

  for (g in 1:m){
    #S update
    S.hat[[g]] = t(T.hat[[g]])%*%ginv(D.hat1)%*%T.hat[[g]]
    S.hat[[g]] = ginv(S.hat[[g]])
  }
  S.hat
}


#VEI model with lists -----------------------------------------------------
VEI_model = function(p, N, m, u, S){
  # m is the number of states
  # S is the covariance matrix

  T.hat=replicate(m, matrix(0,nrow=p,ncol=p), simplify=F)
  D.hat=replicate(m, matrix(0,nrow=p,ncol=p), simplify=F)
  S.hat=replicate(m, matrix(0,nrow=p,ncol=p), simplify=F)
  pi_g = rowSums(u)/N
  iso_delta = c() #matrix(0.0, 1, m)

  for (g in 1:m){
    #T update
    T.hat[[g]] = diag(1,p)
    for (r in 2:p){
      if (r == 2){
        phi.hat = c()
        phi.hat[g] = -S[[g]][2,1]/S[[g]][1,1]
        T.hat[[g]][2,1] = phi.hat[g]
      }
      else{
        phi.hat = matrix(NA, nrow = 1, ncol = r-1)
        phi.hat = t((-1)*ginv(S[[g]][1:r-1,1:r-1])%*%S[[g]][r,1:(r-1)])
        #T update
        T.hat[[g]][r,1:(r-1)] = phi.hat
      }
    }
    #isotropic constraint update
    iso_delta[g] = pi_g[g]/p*tr(T.hat[[g]]%*%S[[g]]%*%t(T.hat[[g]]))
  }

  iso_delta1 = sum(iso_delta)
  D.hat1 = iso_delta1*ramify::eye(p,p)

  for (g in 1:m){
    #S update
    S.hat[[g]] = t(T.hat[[g]])%*%ginv(D.hat1)%*%T.hat[[g]]
    S.hat[[g]] = ginv(S.hat[[g]])
  }
  S.hat
}


#EEI model with lists -----------------------------------------------------
EEI_model = function(p, N, m, u, S, it, iso_delta){
  #m is the number of states
  #S is the covariance matrix

  T.hat=replicate(1, matrix(0,nrow=p,ncol=p), simplify=F)
  D.hat=replicate(m, matrix(0,nrow=p,ncol=p), simplify=F)
  S.hat=replicate(m, matrix(0,nrow=p,ncol=p), simplify=F)
  pi_g = rowSums(u)/N

  if (it == 2){
    iso_delta = 1
  }
  else{
    iso_delta = iso_delta
  }

  T.hat = diag(1,p)
  kappa.list = list()
  for (r in 2:p){
    if (r == 2){
      phi.hat = 0
      phi.num=0
      phi.doom=0
      for (g in 1:m){
        phi.num=phi.num+pi_g[g]*S[[g]][2,1]/iso_delta
        phi.doom=phi.doom+pi_g[g]*S[[g]][1,1]/iso_delta
      }
      phi.hat=-phi.num/phi.doom
      if (is.nan(phi.hat) == T){
        phi.hat = 0
      }
      T.hat[2,1] = phi.hat
    }
    else{
      kappa = matrix(0, r, r)
      for (i in 1:r){
        for (j in 1:r){
          for (g in 1:m){
            kappa[i,j]=kappa[i,j]+pi_g[g]*S[[g]][i,j]/iso_delta
            if (is.nan(kappa[i,j]) == T){
              kappa[i,j] = 0
            }
          }
        }
      }
      kappa=-1*kappa
      phi.hat = matrix(NA, nrow = 1, ncol = r-1)
      phi.hat = t((-1)*ginv(kappa[1:r-1,1:r-1])%*%kappa[r,1:(r-1)])
      #T update
      T.hat[r,1:(r-1)] = phi.hat
    }
  }

  iso_delta1 = 0
  for (g in 1:m){
    #isotropic constraint update
    iso_delta1 = iso_delta1 + pi_g[g]/p*tr(T.hat%*%S[[g]]%*%t(T.hat))
  }
  D.hat1 = iso_delta1*ramify::eye(p,p)

  for (g in 1:m){
    #S update
    S.hat[[g]] = t(T.hat)%*%ginv(D.hat1)%*%T.hat
    S.hat[[g]] = ginv(S.hat[[g]])
  }
  return(list(S.hat=S.hat, iso_delta=iso_delta1))
}


#EVA model with lists -----------------------------------------------------
EVA_model = function(p, N, m, u, S, it, D.hat){
  #m is the number of states
  #S is the covariance matrix

  T.hat=replicate(1, matrix(0,nrow=p,ncol=p), simplify=F)
  S.hat=replicate(m, matrix(0,nrow=p,ncol=p), simplify=F)
  pi_g = rowSums(u)/N

  if (it == 2){
    d=replicate(m, matrix(ramify::eye(p),nrow=p,ncol=p), simplify=F)
    D.hat=replicate(m, matrix(0,nrow=p,ncol=p), simplify=F)
  }
  else{
    d = D.hat
    D.hat=replicate(m, matrix(0,nrow=p,ncol=p), simplify=F)
  }

  T.hat = diag(1,p)
  kappa.list = list()
  for (r in 2:p){
    if (r == 2){
      phi.hat = c()
      phi.num=0
      phi.doom=0
      for (g in 1:m){
        phi.num=phi.num+pi_g[g]*S[[g]][2,1]/d[[g]][2,2]
        phi.doom=phi.doom+pi_g[g]*S[[g]][1,1]/d[[g]][2,2]
      }
      phi.hat=-phi.num/phi.doom
      if (is.nan(phi.hat) == T){
        phi.hat = 0
      }
      T.hat[2,1] = phi.hat
    }
    else{
      kappa = matrix(0, r, r)
      for (i in 1:r){
        for (j in 1:r){
          for (g in 1:m){
            kappa[i,j]=kappa[i,j]+pi_g[g]*S[[g]][i,j]/d[[g]][r,r]
            if (is.nan(kappa[i,j]) == T){
              kappa[i,j] = 0
            }
          }
        }
      }
      kappa=-1*kappa
      phi.hat = matrix(NA, nrow = 1, ncol = r-1)
      phi.hat = t((-1)*ginv(kappa[1:r-1,1:r-1])%*%kappa[r,1:(r-1)])
      #T update
      T.hat[r,1:(r-1)] = phi.hat
    }
  }

  for (g in 1:m){
    #D update
    D.hat[[g]]= (T.hat%*%S[[g]]%*%t(T.hat))
    D.hat[[g]]= diag(diag(D.hat[[g]]), p, p)
  }

  for (g in 1:m){
    #S update
    S.hat[[g]] = t(T.hat)%*%ginv(D.hat[[g]])%*%T.hat
    S.hat[[g]] = ginv(S.hat[[g]])
  }
  return(list(S.hat=S.hat, D.hat = D.hat))
}

#EVI model with lists -----------------------------------------------------
EVI_model = function(p, N, m, u, S, it, iso_delta){
  #m is the number of states
  #S is the covariance matrix

  T.hat=replicate(1, matrix(0,nrow=p,ncol=p), simplify=F)
  D.hat=replicate(m, matrix(0,nrow=p,ncol=p), simplify=F)
  S.hat=replicate(m, matrix(0,nrow=p,ncol=p), simplify=F)
  pi_g = rowSums(u)/N

  if (it == 2){
    iso_delta = rep(1,m)
  }
  else{
    iso_delta = iso_delta
  }

  T.hat = diag(1,p)
  kappa.list = list()

  for (r in 2:p){
    if (r == 2){
      phi.hat = c()
      phi.num=0
      phi.doom=0
      for (g in 1:m){
        phi.num=phi.num+S[[g]][2,1]/iso_delta[g]
        phi.doom=phi.doom+S[[g]][1,1]/iso_delta[g]
      }
      if (is.null(phi.hat) == T){
        phi.hat = 0
      }
      T.hat[2,1] = phi.hat
    }
    else{
      kappa = matrix(0, r, r)
      for (i in 1:r){
        for (j in 1:r){
          for (g in 1:m){
            kappa[i,j] = kappa[i,j]+S[[g]][i,j]/iso_delta[g]
            if (is.nan(kappa[i,j]) == T){
              kappa[i,j] = 0
            }
          }
        }
      }
      kappa=-1*kappa
      phi.hat = matrix(NA, nrow = 1, ncol = r-1)
      phi.hat = t((-1)*ginv(kappa[1:r-1,1:r-1])%*%kappa[r,1:(r-1)])
      #T update
      T.hat[r,1:(r-1)] = phi.hat
    }
  }

  for (g in 1:m){
    #delta update
    iso_delta[g] = 1/p*tr(T.hat%*%S[[g]]%*%t(T.hat))
  }

  D.hat1 = diag(iso_delta, p, p)

  for (g in 1:m){
    #S update
    S.hat[[g]] = t(T.hat)%*%ginv(D.hat1)%*%T.hat
    S.hat[[g]] = ginv(S.hat[[g]])
  }
  return(list(S.hat=S.hat, iso_delta = iso_delta))
}

#EM Support functions
#---------------------------------- Function to solve missingness parameters
solvealpha=function(missingMat,uvec,type,j,m,p,t,v1,uvec_orig){
  alpha.next2=array(NA,dim=c(1,p,t))
  beta.next2=NULL
  n=dim(missingMat)[3]
  if(type=="s"){ #state-dependent
    GLMfit = glm(c(missingMat[missingMat!=2])~1, family=quasibinomial(link="probit"), weights = rep(uvec[j,],p),control = list(maxit = 100))
    alpha.next2[1,,] = rep(c(GLMfit$coefficients),p)
    miss_p=m
  }
  if(type=="sv"){ #state & variable dependent
    for(h in 1:p){
      m2=missingMat[,h,]
      GLMfit = glm(c(m2)[m2!=2]~1, family=quasibinomial(link="probit"), weights = uvec[j,],control = list(maxit = 50))
      alpha.next2[1,h,] =GLMfit$coefficients
    }
    miss_p=m*p
  }
  if(type=="st"){ #state & time dependent
    a = 1:(n*t)
    for(ti in 1:t){
      m2=t(missingMat[ti,,])
      b = a[seq(ti, length(a), t)]
      GLMfit = glm(c(m2)[m2!=2]~1, family=quasibinomial(link="probit"), weights = rep(c(na.omit(uvec_orig[j,b][m2!=2])),p),control = list(maxit = 50))
      alpha.next2[1,,ti] =rep(c(GLMfit$coefficients),p)
    }
    miss_p=m*t
  }
  if(type=="st2"){ #state & time dependent with beta
    beta.next2=array(NA,dim=c(1,p,t))
    v1=rep(v1,p)
    GLMfit = glm(c(missingMat[missingMat!=2])~v1, family=quasibinomial(link="probit"), weights = rep(uvec[j,],p),control = list(maxit = 50))
    alpha.next2[1,,] = rep(c(GLMfit$coefficients[1]),p)
    beta.next2[1,,]= rep(c(GLMfit$coefficients[2]),p)
    miss_p=2*m
  }
  if(type=="svt"){ #state, variable, and time dependent
    a = 1:(n*t)
    for(ti in 1:t){
      b = a[seq(ti, length(a), t)]
      for(h in 1:p){
        m2=missingMat[ti,h,]
        GLMfit = glm(c(m2)[m2!=2]~1, family=quasibinomial(link="probit"), weights = c(na.omit(uvec_orig[j,b][m2!=2])),control = list(maxit = 50))
        alpha.next2[1,h,ti] =GLMfit$coefficients
      }
    }
    miss_p=m*p*t
  }
  if(type=="svt2"){ #state, variable, and time dependent with beta
    beta.next2=array(NA,dim=c(1,p,t))
    for(h in 1:p){
      m2=missingMat[,h,]
      GLMfit = glm(c(m2)[m2!=2]~v1, family=quasibinomial(link="probit"), weights = uvec[j,],control = list(maxit = 50))
      alpha.next2[1,h,] = GLMfit$coefficients[1]
      beta.next2[1,h,]= GLMfit$coefficients[2]
    }
    miss_p=2*m*p
  }
  return(list(alpha=alpha.next2,beta=beta.next2,miss_p=miss_p))
}

#-------------------------------------------- This function calculates the probability of missingness
probOfM=function(alpha_group,missingMat,p,beta_group,type){
  #input is alpha and/or beta, missingness indicators, and type of missingness

  if(type=="s" | type=="sv" | type=="st" | type=="svt"){
    beta_group=NULL
  }
  if(type=="mar"){
    prob=1
  }else{
    if(is.null(beta_group)){
      prob=1
      for(j in 1:p){
        prob = prob*(pnorm(alpha_group[j])^(missingMat[j]))*(1-pnorm(alpha_group[j]))^(1-(missingMat[j]))

      }
    } else{
      prob=1
      t1=1
      for(j in 1:p){
        prob = prob*(pnorm(alpha_group[j]+beta_group[j]*t1)^(missingMat[j]))*(1-pnorm(alpha_group[j]+beta_group[j]*t1))^(1-(missingMat[j]))
        t1=t1+1
      }
    }
  }

  return(prob)
}

#---------------------------------This function calculates the forward and backward probabilites, and log likelihood
alpha_beta=function(data,m,gamma,mu,sigma,alpha,id,delta=NULL,missingMat,beta,type){
  n=length(unique(id))
  t          = nrow(data)/n
  lalpha     =lbeta=matrix(NA,m,n*t)

  drop=any(data==999,na.rm=TRUE)
  miss=any(is.na(data))
  if(drop){lalpha     =lbeta=matrix(NA,m+1,n*t)}
  p=dim(data)[2]
  miss_mat = matrix(as.numeric(is.na(data)),ncol=ncol(data))
  x2=matrix(NA,nrow=nrow(data),ncol=1)
  for(i in 1:nrow(data)){
    if(all(complete.cases(data[i,]))){
      x2[i]=1
    } else{
      x2[i]=NA
    }
  }


  allprobs = array(NA, dim = c(n*t, 1, m))
  if(drop){allprobs = array(1, dim = c(n*t, 1, m+1))}
  if(miss){
    for(i in 1:m){
      for(j in 1:n){
        for(ti in 1:t){
          ind=(j-1)*t+ti
          if(!is.na(x2[ind])){
            allprobs[ind,,i]=pmax(dmvnorm(data[ind,],mean=mu[[i]],sigma=sigma[[i]])*probOfM(alpha[i,,ti],missingMat[ti,,j],p,beta[i,,ti],type),0.1^100)

          } else{
            var_obs = which(miss_mat[ind,]==0)
            if(length(var_obs)==0){
              allprobs[ind,,i]=1*probOfM(alpha[i,,ti],missingMat[ti,,j],p,beta[i,,ti],type) #added
            }else{
              Y_obs = as.matrix(data[ind,var_obs])
              mu_obs = lapply(mu, function(x) x[var_obs])
              sigma_obs = lapply(sigma , function(x) matrix(x[var_obs,var_obs],ncol=length(var_obs),nrow=length(var_obs)))
              allprobs[ind,,i]=pmax(dmvnorm(Y_obs,mean=mu_obs[[i]],sigma=sigma_obs[[i]])*probOfM(alpha[i,,ti],missingMat[ti,,j],p,beta[i,,ti],type),0.1^100)
            }
          }
        }
      }
    }
  }else{
    for(i in 1:m){
      allprobs[,,i]=pmax(dmvnorm(data,mean=mu[[i]],sigma=sigma[[i]]),0.1^100)

    }
  }

  if(drop){
    a = 1:(n*t)
    b = a[seq(1, length(a), t)]
    allprobs[b,1,m+1]=0
    d2=data
    d2[is.na(data)]=1
    tmp=which(d2[,1]!=999)
    allprobs[tmp,1,m+1]=0
    m=m+1

  }


  Li=array(NA,dim=c(m,n,t))
  for(j in 1:n){
    foo        = delta*allprobs[j*t-t+1,1,]
    sumfoo     = sum(foo)
    lscale     = log(sumfoo)
    Li[,j,1]=lscale
    foo        = foo/sumfoo
    lalpha[,j*t-t+1] = foo
    for (i in 2:t)
    {
      foo        = foo%*%gamma*allprobs[i+(j-1)*t,1,]
      sumfoo     = sum(foo)
      lscale     = lscale+log(sumfoo)
      Li[,j,i]=lscale
      foo        = foo/sumfoo
      lalpha[,i+(j-1)*t] = foo
    }

  }
  obsllk=matrix(NA,nrow=1,ncol=n)
  for(i in 1:n){
    c   =  max(Li[,i,t])
    obsllk[1,i]=c+log(sum(exp(Li[,i,t]-c)))
  }
  obs_llk_final=sum(obsllk[1,])

  for(j in 1:n){

    lbeta[,j*t]  = rep(1,m)
    foo        = rep(1,m)
    for (i in (t-1):1)
    {
      foo        = gamma%*%(allprobs[i+1+(j-1)*t,1,]*foo)
      sumfoo     = sum(foo)
      foo        = foo/sumfoo
      lbeta[,i+(j-1)*t]  = foo
    }
  }
  if(drop){
    lalpha[1:(m-1),-tmp]=0
  }
  list(la=lalpha,lb=lbeta,logLi=obsllk, allprob=allprobs,logLfinal=obs_llk_final)
}

state_probs = function(x,m,mu,sigma,gamma,id,alpha,delta,missingMat,beta,type)
{
  dropped=ifelse(x[,1]!=999 | is.na(x[,1]),0,1)
  n=length(unique(id))
  t          = nrow(x)/n
  fb         = alpha_beta(x,m,gamma,mu,sigma,alpha,id,
                          delta,missingMat,beta,type)
  la         = fb$la
  lb         = fb$lb
  stateprobs = matrix(NA,ncol=n*t,nrow=m)
  for(j in 1:m){
    stateprobs[j,]=matrix(c(la[j,]*lb[j,]/colSums(la*lb)),nrow=1,ncol=n*t)
  }
  stateprobs=rbind(stateprobs,rep(0,n*t))
  for(i in 1:(n*t)){
    if(dropped[i]==1){
      stateprobs[m+1,i]=1
      stateprobs[1:m,i]=0
    }
  }
  stateprobs
}

#------------------Local decoding
local_decoding =function(x,m,mu,sigma,gamma,id,alpha,delta,missingMat,beta,type)
{
  estimated=list()
  estimated$probs=state_probs(x,m,mu,sigma,gamma,id,alpha,delta,missingMat,beta,type)
  n=nrow(x)
  est_s = rep(NA,n)
  for (i in 1:n){
    est_s[i]=which.max(estimated$probs[,i])
  }
  estimated$id=est_s
  estimated
}

#---------------------------- EM Function
cdghmm = function(x,m,id,mu=NULL,sigma=NULL,gamma=NULL,delta=NULL,alpha=NULL,beta=NULL,
                   maxiter=10000,tol=1e-6,type="s",covtype="VVA")
{

  warning("Must use 999 to indicate a dropped out observation and NA to indicate a missing observation.")
  warning("If data (x) is not sorted by id and time/visit/response number, the results will be incorrect.")
  if(type!="mar" & type!="s" & type!="sv" & type!="st" & type!="svt" & type!="st2" & type!="svt2") {
    stop("No valid missing mechanism has been provided.")
  }
  if(covtype!="VVA" & covtype!="EEA" & covtype!="VEA" & covtype!="EVA" & covtype!="VVI" & covtype!="VEI" & covtype!="EVI" & covtype!="EEI") {
    stop("No valid CDGHMM family member has been selected")
  }

  drop=any(x==999,na.rm=TRUE)
  miss=any(is.na(x))
  n=length(unique(id))
  p=ncol(x)
  t          = nrow(x)/n
  missingMat=matrix(as.numeric(is.na(x)),nrow=n*t,ncol=p)
  xmiss=x
  xmiss[is.na(xmiss)] = 0
  meanV=colMeans(x,na.rm=T)
  xmiss=xmiss+sweep(missingMat, MARGIN=2,meanV, `*`)

  missingMat2=array(NA, dim = c(t, p,n))

  v1=c()
  v1=rep(1:t,n)
  for(i in 1:n){
    for(ti in 1:t){
      missingMat2[ti,,i]=ifelse(x[(i-1)*t+ti,]==999 & !is.na(x[(i-1)*t+ti,]),2,as.numeric(is.na(x[(i-1)*t+ti,])))
    }
  }


  #initializing delta (stationary distribution)
  if (is.null(delta) == T){
    delta = c()
    for (i in 1:m){
      delta[i] = 1/m
    }
    if(drop){
      delta=c(delta,0)
    }
  }

  ##initializing alpha
  if(miss){
  if (is.null(alpha) == T){
    alpha = array(-0.6755898,dim=c(m,p,t))
  }
  }
  #initializing beta
  if(miss){
  if(type =="st2" | type=="svt2" & is.null(beta)==T){
    beta = array(-0.6755898,dim=c(m,p,t))
  }
  }

  #initializing gamma (transition matrix)
  if (is.null(gamma) == T){
    if(!drop){
      gamma = matrix(1,m,m)+5*diag(m)
      gamma = diag(1/rowSums(gamma))%*%gamma
    }
    if(drop){
      gamma = matrix(1,m+1,m+1)+5*diag(m+1); gamma = diag(1/rowSums(gamma))%*%gamma;
      gamma[m+1,1:m] = 0
      gamma[m+1,m+1]=1
    }
  }

  #dimensions of mu is (p, m)
  if (is.null(mu) == T){
    mu = list()
    group = kmeans(na.omit(x[x[,1]!=999,1:p]),m)
    for (j in 1:m){
      mu[[j]] = group$centers[j,]
    }
  }

  # covariance matrices (p,p,m)
  if (is.null(sigma) == T){
    sigma =list()
    for (j in 1:m){
      sigma[[j]] = cov(na.omit(x[x[,1]!=999,1:p])) + rnorm(1, p)
    }
  }


  llk=c()
  x_imp=list()
  mu_imp=list()
  mu.next     = mu
  sigma.next=sigma
  gamma.next      = gamma
  delta.next      = delta
  alpha.next = alpha
  beta.next = beta
  it = 1

  a = 1:(n*t)
  b = a[seq(1, length(a), t)]

  expect2 = replicate(n*t, vector("list", m), simplify = FALSE)

  for(z in 1:m){
    x_imp[[z]]=xmiss
  }
  x_expected=list()
  if(drop){
    tmp=which(x[,1]==999)
    for(mi in 1:m){
      x_expected[[mi]]=x_imp[[mi]][-tmp,]
    }
    lengthND=dim(x)[1]-length(tmp)
  }else{
    tmp=n*t+1
    for(mi in 1:m){
      x_expected[[mi]]=x_imp[[mi]]
    }
    lengthND=n*t
  }
  missingMat3=missingMat[-tmp,]
  llk_i=matrix(NA,maxiter,n*t)
  for (iter in 1:maxiter)
  {
    it = it + 1
    fb  =  alpha_beta(x,m,gamma,mu,sigma,alpha,id,delta=delta,missingMat2,beta,type)

    if(!drop){
      lallprobs = array(NA, dim = c(n*t, 1, m))
      for(i in 1:m){
        lallprobs[,,i]=fb$allprob[,1,i]
      }
      la = array(NA, dim = c(n*t, m, m))
      lb = array(NA, dim = c(n*t, m, m))
    }
    if(drop){
      lallprobs = array(NA, dim = c(n*t, 1, m+1))
      for(i in 1:(m+1)){
        lallprobs[,,i]=fb$allprob[,1,i]
      }
      la = array(NA, dim = c(n*t, m+1, m+1))
      lb = array(NA, dim = c(n*t, m+1, m+1))
    }

    la=fb$la
    lb=fb$lb



    for(i in 1:n){
      llk_i[iter,(i*t-(t-1)):(i*t)]=fb$logLi[1,i]

    }
    sigma_imp= list()
    for(i in 1:lengthND){
      sigma_imp[[i]] = list()
      for (k in 1:m){
        sigma_imp[[i]][[k]] = matrix(0,nrow=p,ncol=p)
      }
    }

    uvec=matrix(NA,m,n*t)
    if(drop){
      v_vec=array(NA,dim=c(m,m+1,n*t))
      v_vec2=array(NA,dim=c(m,m+1,n*t))
    }else{
      v_vec=array(NA,dim=c(m,m,n*t))
      v_vec2=array(NA,dim=c(m,m,n*t))
    }
    for (j in 1:m)
    {

      if(drop){
        for (k in 1:(m+1))
        {
          for(i in 1:n){
            v_vec[j,k,(i*t-(t-2)):(i*t)]=gamma[j,k]*c(la[j,(i*t-(t-1)):(i*t-1)]*lallprobs[(i*t-(t-2)):(i*t),1,k]*lb[k,(i*t-(t-2)):(i*t)])
          }
        }
      }else{
        for (k in 1:(m))
        {
          for(i in 1:n){
            v_vec[j,k,(i*t-(t-2)):(i*t)]=gamma[j,k]*c(la[j,(i*t-(t-1)):(i*t-1)]*lallprobs[(i*t-(t-2)):(i*t),1,k]*lb[k,(i*t-(t-2)):(i*t)])


          }
        }

      }

      uvec[j,]=matrix(c(la[j,]*lb[j,]/colSums(la*lb)),nrow=1,ncol=n*t)

      for (i in 1:lengthND){
        indO = which(missingMat3[i,]==0)
        indM = which(missingMat3[i,]==1)
        if (length(indM)>0){

          if(length(indO)==0){
            mu_imp[[j]] = mu[[j]][indM]
            sigma_imp[[i]][[j]][indM,indM] = sigma[[j]][indM,indM]

          }else{
            mu_imp[[j]] = tryCatch(
              {
                mu[[j]][indM] + (sigma[[j]][indM,indO]%*%ginv(sigma[[j]][indO,indO]))%*%t(t(x_expected[[j]][i,indO]-mu[[j]][indO]))
              },
              error = function(e) {
                mu[[j]][indM] + (sigma[[j]][indM,indO]%*%ginv(sigma[[j]][indO,indO]))%*%t(x_expected[[j]][i,indO]-mu[[j]][indO])
              }
            )

            sigma_imp[[i]][[j]][indM,indM] = sigma[[j]][indM,indM]-sigma[[j]][indM,indO]%*%ginv(sigma[[j]][indO,indO])%*%sigma[[j]][indO,indM]
          }
          x_expected[[j]][i,which(missingMat3[i,]==1)] = mu_imp[[j]]

        }
        expect2[[i]][[j]] = tryCatch(
          {
            sigma_imp[[i]][[j]] + t(x_expected[[j]][i,]-mu[[j]])%*%t(t(x_expected[[j]][i,]-mu[[j]]))
          },
          error = function(e) {
            sigma_imp[[i]][[j]] + t(t(x_expected[[j]][i,]-mu[[j]]))%*%t(x_expected[[j]][i,]-mu[[j]])
          }
        )

      }

      sigma.next[[j]] = matrix(0,nrow=p,ncol=p)
      mu.next[[j]]=(uvec[j,-tmp]%*%t(t(x_expected[[j]])))/sum(uvec[j,-tmp])
      uvec2=uvec[,-tmp]
      for(i in 1:lengthND){
        sigma.next[[j]]=sigma.next[[j]] + expect2[[i]][[j]]*uvec2[j,i]
      }
      sigma.next[[j]]=sigma.next[[j]]/sum(uvec[j,])

      if(miss==T){
        if(type=="mar"){
          alpha.next=NULL
          beta.next=NULL
          miss_p=0
        }else{
          if(type=="st2" | type=="svt2"){
            t1=solvealpha(missingMat2,uvec[,-tmp],type,j,m,p,t,v1[-tmp],uvec) #v1[-c(tmp)]
            alpha.next[j,,]=t1$alpha
            beta.next[j,,]=t1$beta
            miss_p=t1$miss_p
          }else{
            t1=solvealpha(missingMat2,uvec[,-tmp],type,j,m,p,t,v1[-c(tmp)],uvec)
            beta.next=NULL
            alpha.next[j,,]=t1$alpha
            miss_p=t1$miss_p
          }
        }
      }else{
        miss_p=0
      }

      delta.next[j]=sum(uvec[j,b])

    }

    #update T,D,S
    if(covtype=="EEA"){
      if (it == 2){
        D.hat =replicate(1, matrix(ramify::eye(p),nrow=p,ncol=p), simplify=F)
      }
      else{
        D.hat = EEA_output$D.hat
      }
      EEA_output = EEA_model(p,lengthND, m, uvec[1:m,-tmp], sigma.next, it, D.hat) #should this be u without 999?
      sigma = EEA_output$S.hat
      np = p*(p-1)/2 + p + m*p + m*(m-1) + (m-1)+miss_p

    } else if(covtype=="VVA"){

      sigma = VVA_model(p, m, sigma.next)
      np= m*(p + (p*p - p)/2 + p) + m*(m-1) + (m-1)+miss_p

    } else if(covtype=="VEA"){

      sigma = VEA_model(p,lengthND, m, uvec[1:m,-tmp], sigma.next)
      np = m*(p*(p-1)/2) + p + m*p + m*(m-1) + (m-1)+miss_p

    } else if(covtype=="EVA"){

      if (it == 2){
        iso_delta = NULL
      }
      EVA_output = EVA_model(p,lengthND, m, uvec[1:m,-tmp], sigma.next, it, D.hat)
      sigma = EVA_output$S.hat
      D.hat = EVA_output$D.hat
      np = p*(p-1)/2 + m*p + m*p + m*(m-1) + (m-1)+miss_p

    } else if(covtype=="VVI"){

      sigma = VVI_model(p, m, sigma.next)
      np = m*(p*(p-1)/2) + m + m*p + m*(m-1) + (m-1)+miss_p

    } else if(covtype=="VEI"){
      sigma = VEI_model(p,lengthND, m, uvec[1:m,-tmp], sigma.next)
      np = m*(p*(p-1)/2) + 1 + m*p + m*(m-1) + (m-1)+miss_p

    } else if(covtype=="EVI"){
      if (it == 2){
        iso_delta = NULL
      }
      EVI_output = EVI_model(p,lengthND, m, uvec[1:m,-tmp], sigma.next, it, iso_delta)
      sigma = EVI_output$S.hat
      iso_delta = EVI_output$iso_delta
      np = p*(p-1)/2 + m + m*p + m*(m-1) + (m-1)+miss_p

    } else if(covtype=="EEI"){
      if (it == 2){
        iso_delta = NULL
      }
      EEI_output = EEI_model(p,lengthND, m, uvec[1:m,-tmp], sigma.next, it, iso_delta)
      sigma = EEI_output$S.hat
      iso_delta = EEI_output$iso_delta
      np = p*(p-1)/2 + 1 + m*p + m*(m-1) + (m-1)+miss_p

    } else{
      sigma=sigma.next
      np=(m*p*(p+1))/2+ m*p + m*(m-1) + (m-1) +miss_p
    }



    for(i in 1:(n*t)){
      v_vec2[,,i]=v_vec[,,i]/sum(v_vec[,,i])
    }

    if(drop){
      for(g1 in 1:m){
        for(g2 in 1:(m+1)){
          gamma.next[g1,g2]=sum(na.omit(v_vec2[g1,g2,]))
        }
      }
      gamma.next = gamma.next/apply(gamma.next,1,sum)

      gamma.next[m+1,1:m]=rep(0,m)
      gamma.next[m+1,m+1]=1
      delta.next = delta.next/n
      delta.next[m+1]=0
    }else{
      for(g1 in 1:m){
        for(g2 in 1:(m)){
          gamma.next[g1,g2]=sum(na.omit(v_vec2[g1,g2,]))
        }
      }
      gamma.next = gamma.next/apply(gamma.next,1,sum)
      delta.next = delta.next/n
    }

    llk[iter]=fb$logLfinal
    #cat("---------------------------------------\n")

    mu     = mu.next
    gamma      = gamma.next
    delta      = delta.next
    alpha = alpha.next
    beta = beta.next

    if(iter>1){
      if(abs(llk[iter]-llk[iter-1])<tol){
        AIC    = -2*(llk[iter]-np)
        BIC    = 2*llk[iter]-np*log(n*t)
        states=local_decoding(x,m,mu,sigma,gamma,id,alpha,delta,missingMat2,beta,type)
        penalty=matrix(NA,1,n*t)
        sil=silhouette(states$id, dist(x))
        avg_sil=summary(sil)$avg.width
        for(i in 1: (n*t)){
          if(all(uvec[1:m,i]!=0)){
            penalty[,i]=max(uvec[1:m,i])*sum(log(uvec[1:m,i]))
          }
        }
        ICL=BIC + 2*sum((penalty))

        return(list(mu=mu,sigma=sigma,gamma=gamma,delta=delta,alpha=alpha,beta=beta,llk=llk[iter],AIC=AIC,
                    BIC=BIC,ICL=ICL,Avg_Silhouette=avg_sil,probs=t(states$prob),states=states$id,beta=beta,mod=covtype))
      }
    }

  }
  message(paste("No convergence after",maxiter,"iterations"))
  NA
}
