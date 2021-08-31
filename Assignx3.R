##################################################Stage 1##################################################
###########################################################################################################
#Assigning Treatment Rules at Stage 1
A.assign.1 <- function(X1, X2,L, a_list = c(0L,1L,2L)){
  n_sample <- length(X1)
  if(n_sample != length(X2)) stop("Sample Size Don't Match")
  A_sample <-  vector(mode = "numeric", length = n_sample)
  for(k in 1:n_sample){
    x1 <- X1[k]
    x2 <- X2[k]
    l <- L[k]
    p <- c(exp(-0.2*x1+0.3*x2-0.2*(l==1)+0.5), exp(0.3*x2+1.5*(l==1)+0.5), 1)
    p <- p/sum(p)
    A_sample[k] <- sample(a_list, size = 1, prob = p)
  }
  return(A_sample) 
}

#True Optimal Treatment at Stage 1
A.opt.1 <- function(X1, X2, a_list = c(0L,1L,2L)){
  n_sample <- length(X1)
  if(n_sample != length(X2)) stop("Sample Size Don't Match")
  A_opt <-  vector(mode = "numeric", length = n_sample)
  for(k in 1:n_sample){
    x1 <- X1[k]
    x2 <- X2[k]
    if(x2>=0) A_opt[k] <- 1
    else{
      if(x1>=0) A_opt[k] <- 2
      else A_opt[k] <- 0
    }
  }
  return(A_opt)
}

#Outcome at Stage 1
R.1 <- function(X1, X2,X3, A1, rand = 1){
  n_sample <- length(X1)
  
  if(n_sample != length(X2)) stop("Sample Size Don't Match")
  if(n_sample != length(A1)) stop("Sample Size Don't Match")
  
  A_opt <- A.opt.1(X1, X2)
  R <-  vector(mode = "numeric", length = n_sample)
  R <- exp(1 + 0.05*X2 - 2*abs(A1-A_opt)) + rand*rnorm(n_sample,0,0.5)
  #R <- exp(1 + 0.5*X3 - 2*abs(A1-A_opt)) + rand*rnorm(n_sample,0,0.5)
}


##################################################Stage 2##################################################
###########################################################################################################
#Assign Health Indicator (Label)
Label <- function(X3){
  return(2 - as.numeric(X3>-0.5))#Pick people with higher R1
}

#Assigning Treatment Rules at Stage 2(Related to R1,X1,A1)
A.assign.2 <- function(X1, X2, R1, A1, L, a_list = c(0L,1L,2L)){#what is L doing? NOTHING
  n_sample <- length(X1)
  
  if(n_sample != length(X2)) stop("Sample Size Don't Match")
  if(n_sample != length(R1)) stop("Sample Size Don't Match")
  if(n_sample != length(R1)) stop("Sample Size Don't Match")
  
  A_sample <- vector(mode = "numeric", length = n_sample)
  for(k in 1:n_sample){
    if(A1[k]==1L) p <- c(exp(0.05*(R1[k])-0.05*X1[k]-0.2*(L[k]==1)-1),exp(0.08*(R1[k])-0.2*X1[k]-0.1*(L[k]==1)-1.6),1)
    else p <- c(exp(-0.2*(R1[k]-20)-0.05*X1[k]),exp(-0.08*(R1[k]-20)+0.6*X1[k]),2)
    p <- c(exp(-0.2*(R1[k])-0.05*X1[k]+4),exp(-0.08*(R1[k])+0.6*X1[k]+3),10)
    p <- p/sum(p)
    p[is.na(p)] <- 1
    A_sample[k] <- sample(a_list, size = 1, prob = p)
  }
  return(A_sample) 
}

#True Optimal Treatment at Stage 2
A.opt.2 <- function(X1, X2, R1, A1, a_list = c(0L,1L,2L)){
  n_sample <- length(X1)
  
  if(n_sample != length(X2)) stop("Sample Size Don't Match")
  if(n_sample != length(R1)) stop("Sample Size Don't Match")
  if(n_sample != length(A1)) stop("Sample Size Don't Match")
  
  A_opt <- vector(mode = "numeric", length = n_sample)#This is much too complicated.
  for(k in 1:n_sample){
    if(A1[k]==1L){
        if(R1[k]>1.5) A_opt[k] <- 2L
        else{
        if(R1[k]>0.5) A_opt[k] <- 1L
        else A_opt[k] <- 0L
      }
    }else{
      if(X2[k]>-0.5) A_opt[k] <- 0L
      else A_opt[k] <- 1L
    }
  }
  return(A_opt)
}

R.2 <- function(X1, X2, R1, A1, A2, L, a_list = c(0L,1L,2L), rand = 1){#extra R2 for L=1
  n_sample <- length(X1)
  
  if(n_sample != length(X2)) stop("Sample Size Don't Match")
  if(n_sample != length(R1)) stop("Sample Size Don't Match")
  if(n_sample != length(A1)) stop("Sample Size Don't Match")
  if(n_sample != length(A2)) stop("Sample Size Don't Match")
  if(n_sample != length(L)) stop("Sample Size Don't Match")
  
  R <- vector(mode = "numeric", length = n_sample)
  A_opt <- A.opt.2(X1,X2,R1,A1)
  
  R <- exp(1 +as.numeric(L==1)*(A1==1L) - 3*exp(abs(A2-A_opt))) + rand*rnorm(n_sample,0,0.5)
  
  return(R)
}


Opt <- function(X1,X2,X3) {
  n_sample <- length(X1)
  if(n_sample != length(X2)) stop("Sample Size Don't Match")
  A1.opt <- A.opt.1(X1,X2)
  R1.opt <- R.1(X1,X2,X3,A1.opt,rand=0)
  A2.opt <- A.opt.2(X1, X2, R1.opt, A1.opt)
  A.opt <- cbind(A1.opt,A2.opt)
  for(k in 1:n_sample){
        x1 <- X1[k];x2 <- X2[k];x3 <- X3[k]
        temp_a_list <- c()
        temp_a_value <- c()
        if(A.opt[k,1]==1L){
        for(a1 in c(0L,1L,2L)){
          for(a2 in c(0,1L,2L)){
            if(a1==1L&a2==2L) next
            temp_a_list <- rbind(temp_a_list,c(a1,a2))
            r1 <- R.1(x1,x2,x3,a1,rand=0) #Looking at A1=1 and A2!=2 or A1 can also be changed
            l <- ifelse((x3>-0.5),1,2)
            r2 <- R.2(x1,x2,r1,a1,a2,l,rand=0)
            temp_a_value <- c(temp_a_value,r1+r2)
          }
          }
        a.index.max <- which.max(temp_a_value)
        A.opt[k,] <- temp_a_list[a.index.max,]
        }
  }
  return(A.opt)
}






