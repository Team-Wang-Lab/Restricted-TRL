# functions
######################LOOK AT LINE 91
# function to summarize simulation results
summary2<-function(x){#x is a numerical vector or matrix
  s1<-summary(x)
  sd2<-function(y){
   return(sd(y,na.rm=T))
  }
  if(is.matrix(x)){
    SD<-apply(x,2,sd2)
  } else{
    SD<-sd2(x)
  }
  s2<-list(s1,SD)
  names(s2)<-c("summary","SD")#return summary and sd to each column
  return(s2)
}

# function randomly split the data into K folds, return the group indicators
split.sample<-function(K,N){ #K=number of fold N=number of samples
  rns<-runif(N,0,1)
  group<-rep(1,N)
  i=1
  while(i < K){
    group<-group+(rns>=quantile(rns,i/K))
    i=i+1
  }
  return(group)
}


# function to sample treatment A
# input matrix.pi as a matrix of sampling probabilities, which could be non-normalized
A.sim<-function(matrix.pi){
  N<-nrow(matrix.pi) # sample size
  K<-ncol(matrix.pi) # treatment options
  if(N<=1 | K<=1) stop("Sample size or treatment options are insufficient!")
  if(min(matrix.pi)<0) stop("Treatment probabilities should not be negative!")

  # normalize probabilities to add up to 1 and simulate treatment A for each row
  #pis<-apply(matrix.pi,1,sum)
  #probs<-matrix(NA,N,K)
  probs<-t(apply(matrix.pi,1,function(x){x/sum(x,na.rm = TRUE)}))
  A<-apply(probs,1,function(x) sample(0:(K-1),1,prob = x))
  #A<-rep(NA,N)
  #for(i in 1:N){
  #  probs[i,]<-matrix.pi[i,]/pis[i]
  #  A[i]<-sample(0:(K-1),1,prob=probs[i,]) #A is the actual treatment
  #}
  A
}

# function to estimate propensity score(by a multinomial model)
# input treatment vector=A and covariate matrix=Xs
M.propen<-function(A,Xs,covars = NULL){
  if(ncol(as.matrix(A))!=1) stop("Cannot handle multiple stages of treatments together!")
  if(length(A)!= nrow(as.matrix(Xs))) stop("A and Xs do not match in dimension!")
  if(length(unique(A))<=1) stop("Treament options are insufficient!")
  class.A<-sort(unique(A))#class.A=Unique treatments

  require(nnet)
  s.data<-data.frame(A,Xs)
  # multinomial regression with output suppressed
  if(is.null(covars)){
  model<-capture.output(mlogit<-multinom(A ~., data=s.data))
  }else{model<-capture.output(mlogit<-multinom(as.formula(paste("A ~",paste(covars,collapse = "+"))), data=s.data))}
  s.p<-predict(mlogit,s.data,"probs")#using model to predic the traning data and get probabilities of each outcome
  
  if(length(class.A)==2){
    s.p<-cbind(1-s.p,s.p)
  }
  colnames(s.p)<-paste("pi=",class.A,sep="")

  s.p#matrix of N*length(unique(A))
}

# function to estimate conditional means for multiple stages
# input Y as a continous outcome of interest
# input As = (A1, A2, ...) as a matrix of treatments at multiple stages; Stage t has treatment K_t options labeled as 0, 1, ..., K_t-1.
# input H as a matrix of covariates before assigning final treatment, excluding previous treatment variables
Reg.mu<-function(Y,As,H){
  if(nrow(as.matrix(As))!=nrow(as.matrix(H))) stop("Treatment and Covariates do not match in dimension!")
  Ts<-ncol(as.matrix(As)) # number of stages
  N<-nrow(as.matrix(As))
  if(Ts<0 | Ts>3) stop("Only support 1 to 3 stages!")#As as maxixmum 3 colmns
  H<-as.matrix(H)

  if(Ts==1L){ #one stage
    A1<-as.matrix(As)[,1]
    A1<-as.factor(A1)
    KT<-length(unique(A1)) # treatment options at last stage
    if(KT<2) stop("No multiple treatment options!")

    RegModel<-lm(Y ~ H*A1)
    #print(RegModel)
    #Build linear model between outcome Y~historyX+trt+historyX*trt ??????!!!!!! I think it should be H+A1 %in% H
    #Eg C1<-c(2,4,5,2,4);A1<-c(1,2,1,1,2);B1<-matrix(c(1,2,3,10,2,5,3,9,8,2),ncol = 2);lm(C1~B1*A1);lm(C1~B1+A1%in%B1)
    mus.reg<-matrix(NA,N,KT)
    for(k in 1L:KT) mus.reg[,k]<-predict(RegModel,newdata=data.frame(H,A1=factor(rep(sort(unique(A1))[k],N))))
    # use model to predict when treatment =each treatment what is the predicted outcomeY
    # Thus, mus.reg is a N*KT matrix where each column is the countradictory outcome value Y* given the input H and aT.
    }
  if(Ts==2L){
    A1<-as.matrix(As)[,1];A2<-as.matrix(As)[,2]
    A1<-as.factor(A1);A2<-as.factor(A2)#require A1 A2 to be ordered by the same way
    KT<-length(unique(A2))
    if(KT<2) stop("No multiple treatment options!")

    RegModel<-lm(Y ~ (H + A1)*A2)

    mus.reg<-matrix(NA,N,KT)
    for(k in 1L:KT) mus.reg[,k]<-predict(RegModel,newdata=data.frame(H,A1,A2=factor(rep(sort(unique(A2))[k],N))))#Has only change the last trt
  }
  if(Ts==3L){
    A1<-as.matrix(As)[,1];A2<-as.matrix(As)[,2];A3<-as.matrix(As)[,3]
    A1<-as.factor(A1);A2<-as.factor(A2);A3<-as.factor(A3)
    KT<-length(unique(A3))
    if(KT<2) stop("No multiple treatment options!")

    RegModel<-lm(Y ~ (H + A1 + A2)*A3)

    mus.reg<-matrix(NA,N,KT)
    for(k in 1L:KT) mus.reg[,k]<-predict(RegModel,newdata=data.frame(H,A1,A2,A3=factor(rep(sort(unique(A3))[k],N))))
  }

  output<-list(mus.reg, RegModel)
  names(output)<-c("mus.reg","RegModel")
  output#this is the Y*'s for all the possible treatments in stage K_T(last stage) given history and fixed previous trt
}

# function to calculate AIPW mus
# input outcome Y, treatment vector A, estimated propensity matrix pis.hat and regression-based conditional means mus.reg
mus.AIPW<-function(Y,A,pis.hat,mus.reg){
  class.A<-sort(unique(A))
  K<-length(class.A)
  N<-length(A)
  if(K<2 | N<2) stop("No multiple treatments or samples!")
  if(ncol(pis.hat)!=K | ncol(mus.reg)!=K | nrow(pis.hat)!=N | nrow(mus.reg)!=N) stop("Treatment, propensity or conditional means do not match!")

  #AIPW estimates
  mus.a<-matrix(NA,N,K)
  for(k in 1L:K){
    mus.a[,k]<-(A==class.A[k])*Y/pis.hat[,k]+(1-(A==class.A[k])/pis.hat[,k])*mus.reg[,k]
  }
  mus.a#AIPW Y* that weights the original Y and the Y* depend on model
}

# function to calculate AIPW adaptive contrasts and working orders
# input outcome Y, treatment vector A, estimated propensity matrix pis.hat and regression-based conditional means mus.reg
CL.AIPW<-function(Y,A,pis.hat,mus.reg){
  class.A<-sort(unique(A))
  K<-length(class.A)
  N<-length(A)
  if(K<2 | N<2) stop("No multiple treatments or samples!")
  if(ncol(pis.hat)!=K | ncol(mus.reg)!=K | nrow(pis.hat)!=N | nrow(mus.reg)!=N){
    stop("Treatment, propensity or conditional means do not match!")
  }

  #AIPW estimates
  mus.a<-matrix(NA,N,K)
  for(k in 1:K){
    mus.a[,k]<-(A==class.A[k])*Y/pis.hat[,k]+(1-(A==class.A[k])/pis.hat[,k])*mus.reg[,k]
  }
  #Just use one line
  #mus.a<-mus.AIPW(Y,A,pis.hat,mus.reg)

  # C.a1 and C.a2 are AIPW contrasts; l.a is AIPW working order
  C.a1<-C.a2<-l.a<-rep(NA,N)
  for(i in 1:N){
    # largest vs. second largest
    C.a1[i]<-max(mus.a[i,])-sort(mus.a[i,],decreasing=T)[2]
    # largest vs. smallest
    C.a2[i]<-max(mus.a[i,])-min(mus.a[i,])
    # minus 1 to match A's range of 0,...,K-1
    l.a[i]<-which(mus.a[i,]==max(mus.a[i,]))-1#l.a=the index of maximum outcome-1
  }
  output<-data.frame(C.a1, C.a2, l.a)
  output
}

##################################################################

# choose the optimal treatment with given mu's matrix in one given stage
Opt.A<-function(A,mus.hat){ #A are the treatment options given in this stage, mu.hat are the final counterfactural outcome according to the treatment.
  class.A<-sort(unique(A))
  if(length(class.A)==1){
    trt.opt1<-class.A
    Ey.opt1<-mean(mus.hat)
  } else{
    if(length(A)!= nrow(mus.hat) || length(class.A)!= ncol(mus.hat)){
      stop("Treatment options and mean matrix dimension do not match!")
    }

    # pick a single best treatment for all patients
    c.means<-apply(mus.hat,2,mean)
    Ey.opt1<-max(c.means)
    trt.opt1<-class.A[which(c.means==Ey.opt1)]###Problem: if the order of class.A is different from mus.hat becuase of the sort() this will cause problems
  }
  outs<-list(Ey.opt1,trt.opt1)
  names(outs)<-c("Ey.opt1","trt.opt1")
  outs
}

# combine matrices with different dimensions, add NA for additional rows/columns
combine.mat<-function(m1,m2,by="column"){
  nrow1<-nrow(m1);ncol1<-ncol(m1)
  nrow2<-nrow(m2);ncol2<-ncol(m2)
  if(by=="column"){
    combine<-matrix(NA,max(nrow1,nrow2),ncol1+ncol2)
    combine[1:nrow1,1:ncol1]<-m1
    combine[1:nrow2,(ncol1+1):(ncol1+ncol2)]<-m2
  }
  if(by=="row"){
    combine<-matrix(NA,nrow1+nrow2,max(ncol1,ncol2))
    combine[1:nrow1,1:ncol1]<-m1
    combine[(nrow1+1):(nrow1+nrow2),1:ncol2]<-m2
  }
  combine
}

# split data by X to fit child nodes, calculate new means for each child node
Split.X<-function(X,A,mus.hat,minsplit=20){#X is a covariate for patients; minsplit is the minimum cases in each node
  n<-length(X)
  X.val<-unique(X)
  n.X.val<-length(X.val)
  class.A<-sort(unique(A))

  if(n < 2*minsplit || n.X.val<2L || length(class.A)<2L) return(NULL)

  if(is.numeric(X)==TRUE || is.ordered(X)==TRUE || n.X.val==2L){# is.ordered=ordered factorial feature
    X.val<-sort(X.val)
    # reduce computation by using quantiles
    if(n.X.val>100L){
      X.val<-quantile(X,1:100/100)#change the unique x to quantiles of x(only test 100 possible x as candidates)
      n.X.val<-100L
    }
    Ey.opt1.X<-trt.left.X<-trt.right.X<-rep(NA,n.X.val-1)#initalize E(Y|optX)
    for(i in 1L:(n.X.val-1)){
      left<-which(X<=X.val[i])#left<- index of X's that is less than the evaluated optX candidate
      if(length(left)>=minsplit && length(left)<=n-minsplit){#Make sure after split the resulting two nodes has cases more than minsplit

        left.est<-Opt.A(A[left],mus.hat[left,which(class.A %in% unique(A[left]))])
        right.est<-Opt.A(A[-left],mus.hat[-left,which(class.A %in% unique(A[-left]))])

        trt.left.X[i]<-left.est$trt.opt1
        trt.right.X[i]<-right.est$trt.opt1
        Ey.opt1.X[i]<-length(left)/n*left.est$Ey.opt1+(1-length(left)/n)*right.est$Ey.opt1 #Average the optimum Ey for each cutoff point
      }
    }
    # pick the best split of X
    if(sum(!is.na(Ey.opt1.X))>0L){#check if one of the EY for candidate cutoffs is non NA
      mEy.opt1<-max(Ey.opt1.X, na.rm=T)
      cutoff1<-which(Ey.opt1.X==mEy.opt1)[1]#take the minimum cutoff
      X.cutoff<-X.val[cutoff1]
      trt.L<-trt.left.X[cutoff1]
      trt.R<-trt.right.X[cutoff1]

      output<-data.frame(X.cutoff, mEy.opt1, trt.L, trt.R)
      names(output)<-c("X.subset","mEy.opt1","trt.L","trt.R")
    } else{
      return(NULL)
    }

  }

  if(is.numeric(X)==F && is.ordered(X)==F && n.X.val>2L){
    #this is the case when X is not numerical or ordered catagorical or catagorical with 2 classes== catagorical with more than 2 classes
    n.X.combo<-2^(n.X.val-1)-1#Assume there are c classes for X, each class we can either pick or not pick. -1 for not pick any class. -1 for power for the ????
    X.combo<-combn(X.val,1)#???????unique(X.val)?
    if(n.X.val>3L && n.X.val%%2==1L){#EVEN AND ODD CASES
      for(k in 2L:(n.X.val-1)/2) X.combo<-combine.mat(X.combo,combn(X.val,k))
    }
    if(n.X.val>3L && n.X.val%%2==0L){
      for(k in 2L:(n.X.val/2)){
        if(k<(n.X.val/2)) X.combo<-combine.mat(X.combo,combn(X.val,k))
        if(k==(n.X.val/2)){
          temp.mat<-combn(X.val[-1],k-1)
          first.row<-rep(X.val[1],ncol(temp.mat))
          X.combo<-combine.mat(X.combo,rbind(first.row,temp.mat))
        }
      }
    }

    Ey.opt1.X<-trt.left.X<-trt.right.X<-rep(NA,n.X.combo)
    for(i in 1L:n.X.combo){
      left<-which(X %in% X.combo[,i])
      if(length(left)>=minsplit && length(left)<=n-minsplit){

        left.est<-Opt.A(A[left],mus.hat[left,which(class.A %in% unique(A[left]))])
        right.est<-Opt.A(A[-left],mus.hat[-left,which(class.A %in% unique(A[-left]))])

        trt.left.X[i]<-left.est$trt.opt1
        trt.right.X[i]<-right.est$trt.opt1
        Ey.opt1.X[i]<-length(left)/n*left.est$Ey.opt1+(1-length(left)/n)*right.est$Ey.opt1
      }
    }
    # pick the best split of X
    if(sum(!is.na(Ey.opt1.X))>0L){
      mEy.opt1<-max(Ey.opt1.X, na.rm=T)
      cutoff1<-which(Ey.opt1.X==mEy.opt1)[1]
      X.subset<-X.combo[,cutoff1]
      # change a vector into a single string while removing NA's
      X.subset<-paste(X.subset[!is.na(X.subset)], collapse=" ")
      trt.L<-trt.left.X[cutoff1]
      trt.R<-trt.right.X[cutoff1]

      output<-data.frame(X.subset, mEy.opt1, trt.L, trt.R)#data.frame only has one row
      names(output)<-c("X.subset","mEy.opt1","trt.L","trt.R")
    } else{
      return(NULL)
    }
  }
  return(output)#RETURN all the avg Y*|cutoff=each X
}


### pick the best X for split
best.H<-function(H,A,mus.hat,minsplit=20){
  p<-ncol(H)
  output<-as.data.frame(matrix(NA,p,5))
  output[,1]<-1:p
  colnames(output)<-c("X","X.subset","mEy.opt1","trt.L","trt.R")#output is a p*5 matrix p=number of features.

  for(i in 1:p){
    split.i<-Split.X(X=H[,i],A=A,mus.hat=mus.hat,minsplit=minsplit)

    if(!is.null(split.i)) output[i,-1]<-split.i#output column i other than the first cell is replaced by the best split under X
  }
  if(sum(!is.na(output$mEy.opt1))>0L){
    max.p<-which(output$mEy.opt1==max(output$mEy.opt1,na.rm=T))[1]#indicator for best split feature
    opt.output<-output[max.p,]
    if(opt.output$trt.L==opt.output$trt.R){
      return(NULL)
    } else{
      return(opt.output)#c("X","X.subset","mEy.opt1","trt.L","trt.R")
    }
  } else{
    return(NULL)
  }
}

#################################################
# lookahead
# split X with looking head for 1 more step
Split.X.lh<-function(X,A,H,mus.hat,minsplit=20){
  n<-length(X)
  X.val<-unique(X)
  n.X.val<-length(X.val)
  class.A<-sort(unique(A))

  if(n < 2*minsplit || n.X.val<2L || length(class.A)<2L) return(NULL)

  if(is.numeric(X)==T || is.ordered(X)==T || n.X.val==2L){
    X.val<-sort(X.val)
    # reduce computation by using quantiles
    if(n.X.val>100L){
      X.val<-quantile(X,1:100/100)
      n.X.val<-100L
    }
    Ey.opt1.X<-Ey.opt2.X<-trt.left.X<-trt.right.X<-rep(NA,n.X.val-1)
    for(i in 1L:(n.X.val-1)){
      left<-which(X<=X.val[i])
      if(length(left) >= minsplit && length(left) <= n-minsplit){
        # further split left and right
        # lookahead one step
        best.H.L<-best.H.R<-NULL
        if(length(left)>= 2*minsplit){#After the split left still have room for the next split
          best.H.L<-best.H(H=H[left,],A=A[left],
                           mus.hat=mus.hat[left,which(class.A %in% unique(A[left]))],minsplit=minsplit)#output the best variable to split output cutoff, mean,left, right optimal Y
        }
        if(n-length(left)>= 2*minsplit){
          best.H.R<-best.H(H=H[-left,],A=A[-left],
                           mus.hat=mus.hat[-left,which(class.A %in% unique(A[-left]))],minsplit=minsplit)
        }
        # if unable to further split, then pick 1 treatment
        left.est<-Opt.A(A[left],mus.hat[left,which(class.A %in% unique(A[left]))])#optimal Y in all left area
        right.est<-Opt.A(A[-left],mus.hat[-left,which(class.A %in% unique(A[-left]))])

        trt.left.X[i]<-left.est$trt.opt1 #option
        trt.right.X[i]<-right.est$trt.opt1

        if(!is.null(best.H.L)){
          left.Ey<-best.H.L$mEy.opt1#!!!! if can split always split???? what if the everage after split further is less than split just once?
        } else{
          left.Ey<-left.est$Ey.opt1
        }
        if(!is.null(best.H.R)){
          right.Ey<-best.H.R$mEy.opt1
        } else{
          right.Ey<-right.est$Ey.opt1
        }
        # the max Ey after assigning 2 treatments to two child nodes is used for later splits
        Ey.opt1.X[i]<-length(left)/n*left.est$Ey.opt1+(1-length(left)/n)*right.est$Ey.opt1

        # the max Ey if lookahead is used to choose the optimal X and its threshould
        Ey.opt2.X[i]<-length(left)/n*left.Ey+(1-length(left)/n)*right.Ey
      }
    }
    # pick the best split of X
    if(sum(!is.na(Ey.opt1.X))>0L & sum(!is.na(Ey.opt2.X))>0L){
      mEy.opt1<-max(Ey.opt1.X, na.rm=T)
      mEy.opt2<-max(Ey.opt2.X, na.rm=T)
      # based on lookahead
      cutoff<-which(Ey.opt2.X==mEy.opt2)[1]#!!!! if can split always split???? what if the everage after split further is less than split just once?
      X.cutoff<-X.val[cutoff]
      trt.L<-trt.left.X[cutoff]
      trt.R<-trt.right.X[cutoff]

      output<-data.frame(X.cutoff, mEy.opt1, trt.L, trt.R)
      names(output)<-c("X.subset","mEy.opt1","trt.L","trt.R")
    } else{
      return(NULL)
    }
  }else{
    stop("Lookahead currently only supports numerical or ordinal covariates!")
  }

  return(output)
}

#Best X based on split ahead
best.H.lh<-function(H,A,mus.hat,minsplit=20){
  p<-ncol(H)
  output<-as.data.frame(matrix(NA,p,5))
  output[,1]<-1:p
  colnames(output)<-c("X","X.subset","mEy.opt1","trt.L","trt.R")

  for(i in 1:p){
    split.i<-Split.X.lh(X=H[,i],A=A,H=H,mus.hat=mus.hat,minsplit=minsplit)

    if(!is.null(split.i)) output[i,-1]<-split.i
  }
  if(sum(!is.na(output$mEy.opt1))>0L){
    max.p<-which(output$mEy.opt1==max(output$mEy.opt1,na.rm=T))[1]
    opt.output<-output[max.p,]
    if(opt.output$trt.L==opt.output$trt.R){
      return(NULL)
    } else{
      return(opt.output)
    }
  } else{
    return(NULL)
  }
}

######################################################
### building tree

### limited to K depths of tree branches
### a single trt to each node

# propensity and conditional means are only estimated once using full data

DTRtree<-function(Y,A,H,pis.hat=NULL,m.method=c("AIPW","randomForest"),mus.reg=NULL,depth=5,lambda.pct=0.05,minsplit=20,lookahead=F){
  # initialization
  # indicator for subset data
  n<-length(Y)#number of people
  I.node<-rep(1,n)#indicator of nodes
  class.A<-sort(unique(A))
  output<-matrix(NA,1,5)
  colnames(output)<-c("node","X","cutoff","mEy","trt")

  # estimate mus.hat if not given
  if(m.method[1]=="AIPW"){
    # estimate propenstiy matrix if not given, using all data
    # same propensity for all subset data
    if(is.null(pis.hat)) pis.hat<-M.propen(A=A,Xs=H)
    if(is.null(mus.reg)) mus.reg<-Reg.mu(Y=Y,As=A,H=H)$mus.reg
    mus.hat<-mus.AIPW(Y=Y,A=A,pis.hat=pis.hat,mus.reg=mus.reg)
  } else if(m.method[1]=="randomForest"){
    require(randomForest)
    RF<-randomForest(Y~., data=data.frame(A,H))
    mus.hat<-matrix(NA,n,length(class.A))
    for(i in 1L:length(class.A)) mus.hat[,i]<-predict(RF,newdata=data.frame(A=rep(class.A[i],n),H))
  } else{
    stop("The method for estimating conditional means is not available!")
  }

  # expected outcome at root
  root<-Opt.A(A,mus.hat)
  Ey0<-root$Ey.opt1

  # split if improved at least lambda, as a percent of Ey0
  lambda<-abs(Ey0)*lambda.pct

  for(k in 1L:depth){#depth is the most number of split to reach one terminal node
    output<-rbind(output,matrix(NA,2^k,5))#2^k??????? originally output=1*5
    output[,1]<-1L:(2^(k+1)-1)# this does not equal to the number of rows(2^k+1)
    if(k==1L){
      # apply lookahead to the first split, the most important split
      # only to first split so as to save computation time
      # use a larger minsplit for the first split
      if(lookahead){
        best.H.1<-best.H.lh(H=H,A=A,mus.hat=mus.hat,minsplit=0.15*n)
      } else{
        best.H.1<-best.H(H=H,A=A,mus.hat=mus.hat,minsplit=minsplit)
      }
      if(is.null(best.H.1)==F && best.H.1$mEy.opt1>Ey0+lambda){#meet the split criteria
        output[k,-1]<-c(best.H.1$X, best.H.1$X.subset, best.H.1$mEy.opt1, NA)
        I.node[I.node==k & H[,best.H.1$X] <= best.H.1$X.subset]<-2*k
        output[2*k,-1]<-c(NA,NA,NA,best.H.1$trt.L)
        I.node[I.node==k & H[,best.H.1$X] > best.H.1$X.subset]<-2*k+1
        output[2*k+1,-1]<-c(NA,NA,NA,best.H.1$trt.R)
      } else{
        output[k,4:5]<-c(root$Ey.opt1,root$trt.opt1)
        break
      }
    } else{
      for(j in (2^(k-1)):(2^k-1)){
        if(!is.na(output[trunc(j/2),2])){
          best.H.j<-best.H(H=H[I.node==j,],A=A[I.node==j],mus.hat=mus.hat[I.node==j,],minsplit=minsplit)
          if(is.null(best.H.j)==F && best.H.j$mEy.opt1>output[trunc(j/2),4]+lambda){
            output[j,-1]<-c(best.H.j$X, best.H.j$X.subset, best.H.j$mEy.opt1, NA)
            I.node[I.node==j & H[,best.H.j$X] <= best.H.j$X.subset]<-2*j
            output[2*j,-1]<-c(NA,NA,NA,best.H.j$trt.L)
            I.node[I.node==j & H[,best.H.j$X] > best.H.j$X.subset]<-2*j+1
            output[2*j+1,-1]<-c(NA,NA,NA,best.H.j$trt.R)
          }
        }
      }
      if(sum(is.na(output[(2^(k-1)):(2^k-1),2]))==2^(k-1)) break
    }
  }
  output<-output[!is.na(output[,2]) | !is.na(output[,5]),]
  return(output)
}



#######################################################
# predit optimal treatment using the output from DTRtree

predict.DTR<-function(treeout,newdata){
  n<-nrow(newdata)
  predicts<-rep(NA,n)

  # treeout is supposed to be a matrix
  # if there is no split
  if(length(treeout)==5){
    predicts<-rep(treeout[5],n)
  } else{ # if there are splits
    treeout<-as.data.frame(treeout)
    newdata<-as.data.frame(newdata)

    for(i in 1:n){
      nd<-1
      while(is.na(treeout$trt[treeout$node==nd])){
        if(newdata[i,treeout$X[treeout$node==nd]] <= treeout$cutoff[treeout$node==nd]){#if the node<= cutoff
          nd=2*nd #yes proceed first
        } else{
          nd=2*nd+1#then no
        }
      }
      predicts[i]<-treeout$trt[treeout$node==nd]
    }
  }
  return(predicts)
}

###################################################
# function of cross-validation to decide lambda
CV.lambda<-function(Y,A,H,g.opt,lambda.pct,K=10,m.method=c("AIPW","randomForest"),depth=5,minsplit=20){
  N<-length(Y)
  I.group<-split.sample(K,N)
  ppower<-0

  for(i in 1L:K){
    tree1<-DTRtree(Y[I.group != i],A[I.group != i],H[I.group != i,],m.method=m.method,
                   depth=depth,lambda.pct=lambda.pct,minsplit=minsplit)
    g.tree<-predict.DTR(tree1,H[I.group==i,])
    ppower<-ppower + mean(g.tree==g.opt[I.group==i])
  }
  return(ppower/K)
}

##########################################################################
##########################################################################
# tree of Laber & Zhao

# choose the optimal treatment with given estimated max mu and propenstiy for observed A
Opt.A.LZ<-function(Y,A,mus.max,pis.A){
  N<-length(Y)
  class.A<-sort(unique(A))

  if(length(class.A)==1){
    trt.opt1<-class.A
    Ey.opt1<-mean(Y)
    weights<-N
  } else{
    # pick a single best treatment for all patients
    # IPW
    Resids<-rep(NA,length(class.A))
    for(a in 1L:length(class.A)) Resids[a]<-sum((Y-mus.max)*I(A==class.A[a])/pis.A)/sum(I(A==class.A[a])/pis.A)

    trt.opt1<-class.A[which(Resids==max(Resids))]

    Ey.opt1<-sum(Y*I(A==trt.opt1)/pis.A)/sum(I(A==trt.opt1)/pis.A)

    weights<-sum(I(A==trt.opt1)/pis.A)
  }

  outs<-list(Ey.opt1,trt.opt1,weights)
  names(outs)<-c("Ey.opt1","trt.opt1","weights")
  outs
}

Split.X.LZ<-function(X,Y,A,mus.max,pis.A,minsplit=20){
  n<-length(X)
  X.val<-unique(X)
  n.X.val<-length(X.val)
  class.A<-sort(unique(A))

  if(n < 2*minsplit || n.X.val<2L || length(class.A)<2L) return(NULL)

  if(is.numeric(X)==T || is.ordered(X)==T || n.X.val==2L){
    X.val<-sort(X.val)
    # reduce computation by using quantiles
    if(n.X.val>100L){
      X.val<-quantile(X,1:100/100)
      n.X.val<-100L
    }
    Ey.opt1.X<-trt.left.X<-trt.right.X<-rep(NA,n.X.val-1)
    for(i in 1L:(n.X.val-1)){
      left<-which(X<=X.val[i])
      if(length(left)>=minsplit && length(left)<=n-minsplit){

        left.est<-Opt.A.LZ(Y=Y[left],A=A[left],mus.max=mus.max[left],pis.A=pis.A[left])
        right.est<-Opt.A.LZ(Y=Y[-left],A=A[-left],mus.max=mus.max[-left],pis.A=pis.A[-left])

        trt.left.X[i]<-left.est$trt.opt1
        trt.right.X[i]<-right.est$trt.opt1
        Ey.opt1.X[i]<-(left.est$Ey.opt1*left.est$weights + right.est$Ey.opt1*right.est$weights)/(left.est$weights + right.est$weights)
      }
    }
    # pick the best split of X
    if(sum(!is.na(Ey.opt1.X))>0L){
      mEy.opt1<-max(Ey.opt1.X, na.rm=T)
      cutoff1<-which(Ey.opt1.X==mEy.opt1)[1]
      X.cutoff<-X.val[cutoff1]
      trt.L<-trt.left.X[cutoff1]
      trt.R<-trt.right.X[cutoff1]

      output<-data.frame(X.cutoff, mEy.opt1, trt.L, trt.R)
      names(output)<-c("X.subset","mEy.opt1","trt.L","trt.R")
    } else{
      return(NULL)
    }

  }

  if(is.numeric(X)==F && is.ordered(X)==F && n.X.val>2L){
    n.X.combo<-2^(n.X.val-1)-1
    X.combo<-combn(X.val,1)
    if(n.X.val>3L && n.X.val%%2==1L){
      for(k in 2L:(n.X.val-1)/2) X.combo<-combine.mat(X.combo,combn(X.val,k))
    }
    if(n.X.val>3L && n.X.val%%2==0L){
      for(k in 2L:(n.X.val/2)){
        if(k<(n.X.val/2)) X.combo<-combine.mat(X.combo,combn(X.val,k))
        if(k==(n.X.val/2)){
          temp.mat<-combn(X.val[-1],k-1)
          first.row<-rep(X.val[1],ncol(temp.mat))
          X.combo<-combine.mat(X.combo,rbind(first.row,temp.mat))
        }
      }
    }

    Ey.opt1.X<-trt.left.X<-trt.right.X<-rep(NA,n.X.combo)
    for(i in 1L:n.X.combo){
      left<-which(X %in% X.combo[,i])
      if(length(left)>=minsplit && length(left)<=n-minsplit){

        left.est<-Opt.A.LZ(Y=Y[left],A=A[left],mus.max=mus.max[left],pis.A=pis.A[left])
        right.est<-Opt.A.LZ(Y=Y[-left],A=A[-left],mus.max=mus.max[-left],pis.A=pis.A[-left])

        trt.left.X[i]<-left.est$trt.opt1
        trt.right.X[i]<-right.est$trt.opt1
        Ey.opt1.X[i]<-(left.est$Ey.opt1*left.est$weights + right.est$Ey.opt1*right.est$weights)/(left.est$weights + right.est$weights)
      }
    }
    # pick the best split of X
    if(sum(!is.na(Ey.opt1.X))>0L){
      mEy.opt1<-max(Ey.opt1.X, na.rm=T)
      cutoff1<-which(Ey.opt1.X==mEy.opt1)[1]
      X.subset<-X.combo[,cutoff1]
      # change a vector into a single string while removing NA's
      X.subset<-paste(X.subset[!is.na(X.subset)], collapse=" ")
      trt.L<-trt.left.X[cutoff1]
      trt.R<-trt.right.X[cutoff1]

      output<-data.frame(X.subset, mEy.opt1, trt.L, trt.R)
      names(output)<-c("X.subset","mEy.opt1","trt.L","trt.R")
    } else{
      return(NULL)
    }
  }
  return(output)
}


### pick the best X for split
best.H.LZ<-function(H,Y,A,mus.max,pis.A,minsplit=20){
  p<-ncol(H)
  output<-as.data.frame(matrix(NA,p,5))
  output[,1]<-1:p
  colnames(output)<-c("X","X.subset","mEy.opt1","trt.L","trt.R")

  for(i in 1:p){
    split.i<-Split.X.LZ(X=H[,i],Y=Y,A=A,mus.max=mus.max,pis.A=pis.A,minsplit=minsplit)

    if(!is.null(split.i)) output[i,-1]<-split.i
  }
  if(sum(!is.na(output$mEy.opt1))>0L){
    max.p<-which(output$mEy.opt1==max(output$mEy.opt1,na.rm=T))[1]
    opt.output<-output[max.p,]
    if(opt.output$trt.L==opt.output$trt.R){
      return(NULL)
    } else{
      return(opt.output)
    }
  } else{
    return(NULL)
  }
}

LZtree<-function(Y,A,H,mus.hat=NULL,pis.hat=NULL,depth=5,lambda.pct=0.05,minsplit=20){
  # initialization
  # indicator for subset data
  n<-length(Y)
  I.node<-rep(1,n)
  class.A<-sort(unique(A))
  output<-matrix(NA,1,5)
  colnames(output)<-c("node","X","cutoff","mEy","trt")

  # estimate propenstiy matrix if not given, using all data
  # propensity for observed A
  if(is.null(pis.hat)==T){
    pis.hat<-M.propen(A=A,Xs=H)
  }
  pis.A<-rep(NA,N)
  for(i in 1:N){
    pis.A[i]<-pis.hat[i,which(class.A==A[i])]
  }

  # estimate mus.hat if not given, using all data and random forest
  if(is.null(mus.hat)==T){
    require(randomForest)
    RF<-randomForest(Y~., data=data.frame(A,H))
    mus.hat<-matrix(NA,n,length(class.A))
    for(i in 1L:length(class.A)) mus.hat[,i]<-predict(RF,newdata=data.frame(A=rep(class.A[i],n),H))
  }
  mus.max<-apply(mus.hat,1,max)

  # expected outcome at root
  root<-Opt.A.LZ(Y,A,mus.max,pis.A)
  Ey0<-root$Ey.opt1

  # split if improved at least lambda, as a percent of Ey0
  lambda<-abs(Ey0)*lambda.pct

  for(k in 1L:depth){
    output<-rbind(output,matrix(NA,2^k,5))
    output[,1]<-1L:(2^(k+1)-1)
    if(k==1L){
      best.H.1<-best.H.LZ(H=H,Y=Y,A=A,mus.max=mus.max,pis.A=pis.A,minsplit=minsplit)
      if(is.null(best.H.1)==F && best.H.1$mEy.opt1>Ey0+lambda){
        output[k,-1]<-c(best.H.1$X, best.H.1$X.subset, best.H.1$mEy.opt1, NA)
        I.node[I.node==k & H[,best.H.1$X] <= best.H.1$X.subset]<-2*k
        output[2*k,-1]<-c(NA,NA,NA,best.H.1$trt.L)
        I.node[I.node==k & H[,best.H.1$X] > best.H.1$X.subset]<-2*k+1
        output[2*k+1,-1]<-c(NA,NA,NA,best.H.1$trt.R)
      } else{
        output[k,4:5]<-c(root$Ey.opt1,root$trt.opt1)
        break
      }
    } else{
      for(j in (2^(k-1)):(2^k-1)){
        if(!is.na(output[trunc(j/2),2])){
          best.H.j<-best.H.LZ(H=H[I.node==j,],Y=Y[I.node==j],A=A[I.node==j],mus.max=mus.max[I.node==j],
                              pis.A=pis.A[I.node==j],minsplit=minsplit)
          if(is.null(best.H.j)==F && best.H.j$mEy.opt1>output[trunc(j/2),4]+lambda){
            output[j,-1]<-c(best.H.j$X, best.H.j$X.subset, best.H.j$mEy.opt1, NA)
            I.node[I.node==j & H[,best.H.j$X] <= best.H.j$X.subset]<-2*j
            output[2*j,-1]<-c(NA,NA,NA,best.H.j$trt.L)
            I.node[I.node==j & H[,best.H.j$X] > best.H.j$X.subset]<-2*j+1
            output[2*j+1,-1]<-c(NA,NA,NA,best.H.j$trt.R)
          }
        }
      }
      if(sum(is.na(output[(2^(k-1)):(2^k-1),2]))==2^(k-1)) break
    }
  }
  output<-output[!is.na(output[,2]) | !is.na(output[,5]),]
  return(output)
}













