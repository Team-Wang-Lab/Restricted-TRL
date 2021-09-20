#Simulation for Restricted DTR
##a1 = 1 and a2 = 2 is the nonviable arm.
##Naive method (naive T-RL with partial record deletion) is correcting the patients' outcomes in stage 2 as considering 
##only a2 = 1 or a2 = 3 and in stage 1 we correct their pseudo outcomes accordingly.
##Wrong propensity score model used in this simulation.

source("TRL Functions.R")
source('Assignx3.R')
require(randomForest)


results_sim <- function(seed,K=5000){
  #seed = random seed
  #K = sample size    
  set.seed(seed)
  X1 <- rnorm(K,0,1)
  X2 <- rnorm(K,0,1)
  X3 <- rnorm(K,0,1)
  
  
  L <- ifelse((X3>-0.5),1,0)
  A1 <- A.assign.1(X1,X2,L)
  R1 <- R.1(X1,X2,X3,A1)#depend on X3
  A2 <- A.assign.2(X1,X2,R1,A1,L)
  R2 <- R.2(X1,X2,R1,A1,A2,L)
  table(A1,A2)
  
  X <- as.matrix(cbind(X1,X2))
  R <- as.matrix(cbind(R1,R2))
  A <- as.matrix(cbind(A1,A2))
  
  #Estimate Conditional mean and Response at Stage 2
  mu.hat <- Reg.mu(R[,2],A[,2],cbind(X[,2],R1,A[,1],L))[1]$mus.reg
  
  #Fit Tree-Base Treatment Regime Rule for Constrained Data at Stage 2
  #(Does not count people with a1=1 a2=2 when fitting the tree)
  A_res <- (A[,1]==1L&A[,2]!=2L) 
  A_res_all <- (A[,1]==1L)
  X_constrain <- X[A_res,]
  A_constrain <- as.matrix(A[A_res,2],ncol=1)
  R_constrain <- R[A_res,]
  Y_constrain <- R_constrain[,2]
  mu.hat_constrain <- mu.hat[A_res,]
  
  #Propensity scores for stage 2 restricted route
  
  pi.hat_constrain_all <- M.propen(A[A_res_all,2],cbind(exp(X[A_res_all,1]),
                                                        L[A_res_all],R[A_res_all,1]))
  pi.hat_constrain_all <- pi.hat_constrain_all/apply(pi.hat_constrain_all[,c(1,2)],
                                                     MARGIN=1,FUN=sum)
  pi.hat_constrain <- pi.hat_constrain_all[A[A_res_all,2]!=2L,c(1,2)]
  
  #Estimate Optimal DTR
  tree.2.constrain <- DTRtree(Y=Y_constrain,A=A_constrain,H=cbind(X_constrain,R_constrain[,1]),
                              pis.hat=pi.hat_constrain,mus.reg=mu.hat_constrain[,c(1,2)],
                              m.method="AIPW",minsplit=20,lambda.pct=0.05,depth=3)
  
  #Fit Tree-Base Treatment Regime Rule for Unconstrained Data at Stage 2 ( This is for people with a1!=1)
  A_unres <- (A[,1]!=1L)
  X_unconstrain <- X[A_unres,]
  A_unconstrain <- as.matrix(A[A_unres,2],ncol=1)
  R_unconstrain <- R[A_unres,]
  Y_unconstrain <- as.matrix(R_unconstrain[,2],ncol=1)
  mu.hat_unconstrain <- mu.hat[A_unres,]
  pi.hat_unconstrain <- M.propen(A_unconstrain,cbind(exp(X_unconstrain[,1]),
                                                     R_unconstrain[,1]))
  
  tree.2.unconstrain <- DTRtree(Y=Y_unconstrain,A=A_unconstrain,
                                H=cbind(X_unconstrain,A[A_unres,1],R_unconstrain[,1]),
                                pis.hat=pi.hat_unconstrain,m.method="AIPW",
                                mus.reg=mu.hat_unconstrain,minsplit=20,lambda.pct=0.05,depth=3)
  
  #Predict Optimal Treatment for Constrained Data
  predict.2.constrain <- predict.DTR(tree.2.constrain, cbind(X[A_res_all,],R[A_res_all,1]))
  predict.2.unconstrain <- predict.DTR(tree.2.unconstrain, cbind(X_unconstrain,A[A_unres,1],R_unconstrain[,1]))
  
  #Train RandomForest to Estimate R2
  RF.2 <-randomForest(R2~.,data=data.frame(A2,X,A1,R1))
  
  #Rearrange Data
  X.new <- rbind(X[A_res_all,],X[A_unres,])
  R1.new <- as.matrix(c(R[A_res_all,1],R[A_unres,1]), ncol=1)
  A1.new <- as.matrix(c(A[A_res_all,1],A[A_unres,1]), ncol=1)
  A2.new <- as.matrix(c(predict.2.constrain,predict.2.unconstrain), ncol=1)#predicted optimal A2
  A2.old <- as.matrix(c(A[A_res_all,2],A[A_unres,2]), ncol=1)#original A2
  
  R2.new <- as.matrix(c(R[A_res_all,2],R[A_unres,2]), ncol=1)
  L.new <-  as.matrix(c(L[A_res_all],L[A_unres]),ncol=1)
  
  R2.optimal <- predict(RF.2,cbind(A2.new,X.new,A1.new,R1.new))
  R2.old <-  predict(RF.2,cbind(A2.old,X.new,A1.new,R1.new))
  
  
  #Compute Pseudo-Outcome at Stage 1
  R1.pseudo <- R1.new + R2.new + R2.optimal - R2.old 
  
  ###########Naive Method####################################################################################
  res_route<-(!(A[,1]==1L&A[,2]==2L))
  pi.hat_n <- M.propen(A[res_route,2],cbind(exp(X[res_route,2]),X[res_route,1],
                                            R[res_route,1],(L.new[res_route]==1)))#Model to assign A2
  A_n <- as.vector(A[res_route,2])
  R_n <- R[res_route,]
  Y_n <- as.vector(R_n[,2])
  mu.hat_n <- mu.hat[res_route,]
  
  tree.2.naive <- DTRtree(Y=Y_n,A=A_n,H=cbind(X[res_route,],R_n[,1],A[res_route,1]),
                          pis.hat=pi.hat_n,m.method="AIPW",mus.reg=mu.hat_n,
                          minsplit=20,lambda.pct=0.05,depth=3)
  
  #prediction
  predict.2.naive <- predict.DTR(tree.2.naive, 
                                 cbind(X[res_route,],R[res_route,1],A[res_route,1]))
  
  RF.2_n <-randomForest(R2[res_route]~.,data=data.frame(A2=A2[res_route],X=X[res_route,],
                                                        A1=A1[res_route],R1=R1[res_route]))
  
  R1.pseudo.naive<-R1[res_route]+R2[res_route]+
    predict(RF.2_n,newdata=data.frame(A2=predict.2.naive,
                                      X=X[res_route,],A1=A1[res_route],R1=R1[res_route]))-
    predict(RF.2_n,newdata=data.frame(A2=A[res_route,2],
                                      X=X[res_route,],A1=A1[res_route],R1=R1[res_route]))
  
  R2pred_diff = predict(RF.2_n,newdata=data.frame(A2=rep(0,length(A2[!res_route])),
                                                  X=X[!res_route,],
                                                  A1=A1[!res_route],
                                                  R1=R1[!res_route]))-
    predict(RF.2_n,newdata=data.frame(A2=rep(1,length(A2[!res_route])),
                                      X=X[!res_route,],
                                      A1=A1[!res_route],
                                      R1=R1[!res_route]))
  
  R1.pseudo_n_nonviable = R1[!res_route] + 
    predict(RF.2_n,newdata=data.frame(A2=ifelse(R2pred_diff>=0,0,1),
                                      X=X[!res_route,],
                                      A1=A1[!res_route],
                                      R1=R1[!res_route]))
  
  #Rearrange Data
  X.new.naive <- rbind(X[res_route,],X[!res_route,])
  R1.new.naive <- as.matrix(c(R[res_route,1],R[!res_route,1]), ncol=1)
  A1.new.naive <- as.matrix(c(A[res_route,1],A[!res_route,1]), ncol=1)
  A2.new.naive <- as.matrix(c(predict.2.naive,ifelse(R2pred_diff>=0,0,1)), ncol=1)#predicted optimal A2
  A2.old.naive <- as.matrix(c(A[res_route,2],A[!res_route,2]), ncol=1)#original A2
  
  R2.new.naive <- as.matrix(c(R[res_route,2],R[!res_route,2]), ncol=1)
  L.new.naive <-  as.matrix(c(L[res_route],L[!res_route]),ncol=1)
  
  #Compute Pseudo-Outcome at Stage 1
  R1.pseudo_n <- as.matrix(c(R1.pseudo.naive,R1.pseudo_n_nonviable), ncol=1)
  
  ##################################################################################################################
  #Fit Tree-Base Treatment Regime Rule at Stage 1
  X3.new<-c(X3[A_res_all],X3[A_unres])
  pi.hat.1 <- M.propen(A1.new,cbind(X.new[,1],exp(X.new[,2]),(L.new==1)))
  mus.reg.1 <- Reg.mu(R1.pseudo,A1.new,cbind(X.new,X3.new))[1]$mus.reg #True mu regression
  tree.1 <- DTRtree(Y=R1.pseudo,A=A1.new,H=X.new,pis.hat=pi.hat.1,m.method="AIPW",mus.reg=mus.reg.1,minsplit=20,lambda.pct=0.05,depth=3)
  predict.1 <- predict.DTR(tree.1,X.new)
  
  #Stage 1 tree Naive
  pi.hat.1_n <- M.propen(A1,cbind(X.new.naive[,1],exp(X.new.naive[,2]),L.new.naive==1))
  mus.reg.1_n <- Reg.mu(R1.pseudo_n,A1.new.naive,cbind(X.new.naive,L.new.naive))[1]$mus.reg
  tree.1.naive <- DTRtree(Y=R1.pseudo_n,A=A1.new.naive,H=X.new.naive,
                          pis.hat=pi.hat.1_n,m.method="AIPW",
                          mus.reg=mus.reg.1_n,minsplit=20,lambda.pct=0.05,depth=3)
  predict.1.naive <- predict.DTR(tree.1.naive,X.new.naive)
  
  s1.new_m<-numeric(1L)
  s2.new_m<-numeric(1L)
  s1s2.new_m<-numeric(1L)
  s1.naive<-numeric(1L)
  s2.naive<-numeric(1L)
  s1s2.naive<-numeric(1L)
  

  set.seed(seed*50)
  
  K_test <- 1000
  X1.test <- rnorm(K_test,0,1)
  X2.test <- rnorm(K_test,0,1)
  X3.test <- rnorm(K_test,0,1)
  A1.test <- A.assign.1(X1.test,X2.test,(X3.test>-0.5))
  R1.test <- R.1(X1.test,X2.test,X3.test,A1.test)#depend on X3
  A2.test <- A.assign.2(X1.test,X2.test,R1.test,A1.test,(X3.test>-0.5))
  R2.test <- R.2(X1.test,X2.test,R1.test,A1.test,A2.test,(X3.test>-0.5))
  
  A1.predict <- predict.DTR(tree.1,cbind(X1.test,X2.test))
  A1.pred_naive<-predict.DTR(tree.1.naive,cbind(X1.test,X2.test))
  R1.predict <- R.1(X1.test,X2.test,ifelse((X3.test)>-0.5,1,2),A1.predict)
  R1.pred_naive <- R.1(X1.test,X2.test,ifelse((X3.test)>-0.5,1,2),A1.pred_naive)
  
  A2.predict <- A1.predict
  A2.pred_naive<-A1.pred_naive
  for(k in 1:K_test){
    A2.pred_naive[k]<-predict.DTR(tree.2.naive, matrix(c(X1.test[k],X2.test[k],R1.pred_naive[k],A1.pred_naive),nrow=1))
    if(A1.predict[k]==1L) A2.predict[k] <-  predict.DTR(tree.2.constrain, matrix(c(X1.test[k],X2.test[k],R1.predict[k]),nrow=1))
    else A2.predict[k] <- predict.DTR(tree.2.unconstrain, matrix(c(X1.test[k],X2.test[k],A1.predict,R1.predict[k]),nrow=1))
  }
  R2.predict <- R.2(X1=X1.test, X2=X2.test, R1=R1.predict, A1=A1.predict, A2=A2.predict, L=ifelse((X3.test)>-0.5,1,2))
  R2.pred_naive <- R.2(X1=X1.test, X2=X2.test, R1=R1.pred_naive, A1=A1.pred_naive, A2=A2.pred_naive, L=ifelse((X3.test)>-0.5,1,0))
  Y.pred <- R1.predict+R2.predict
  Y.pred_naive <- R1.pred_naive+R2.pred_naive
  Y.observe <- R1.test+R2.test
  Ares <- ifelse(A1.predict==1&A2.predict==2,1,0)
  Ares_naive <- ifelse(A1.pred_naive==1&A2.pred_naive==2,1,0)
  
  A.predict <- cbind(A1.predict,A2.predict)
  A.pred_naive<-cbind(A1.pred_naive,A2.pred_naive)
  A.opt <- Opt(X1.test,X2.test,X3.test)
  R1.true <- R.1(X1.test,X2.test,ifelse((X3.test)>-0.5,1,0),A.opt[,1])
  R2.true <- R.2(X1=X1.test, X2=X2.test, R1=R1.true, A1=A.opt[,1], A2=A.opt[,2], L=ifelse((X3.test)>-0.5,1,0))
  Y.true <- R1.true+R2.true
  
  s1 <- sum(A.opt[,1]==A.predict[,1])
  s2 <- sum(A.opt[,1]==A.predict[,1]&A.opt[,2]==A.predict[,2])
  s21 <- sum(A.opt[,2]==A.predict[,2])
  s3 <- sum(A.opt[,1]==A.pred_naive[,1])
  s4 <- sum(A.opt[,1]==A.pred_naive[,1]&A.opt[,2]==A.pred_naive[,2])
  s22 <- sum(A.opt[,2]==A.pred_naive[,2])
  
  print(s2/K_test)#Opt1+Opt2 are correctly  predicted
  print(s4/K_test)
  
  s1.new_m <-s1/K_test #stage 1 accuracy new method
  s2.new_m <-s21/K_test #stage 2 accuracy new method
  s1s2.new_m <-s2/K_test #stage 1 and 2 accuracy new method
  s1.naive <-s3/K_test #stage 1 accuracy naive method
  s2.naive <-s22/K_test #stage 2 accuracy naive method
  s1s2.naive <-s4/K_test #stage 1 and 2 accuracy naive method
  
  return(list(c(s1.new_m,s2.new_m,s1s2.new_m,s1.naive,s2.naive,s1s2.naive),
              data.frame(Y.observe,Y.pred,Y.pred_naive,Y.true,Ares,Ares_naive)))
}

# Run function iteratively.

p_results3000 <- data.frame()
p_results3000Y <- data.frame()
p_results5000 <- data.frame()
p_results5000Y <- data.frame()

for(i in 1:100){
  print(i)
  temp1 <- results_sim(i,K=3000)
  temp2 <- results_sim(i,K=5000)
  p_results3000<-rbind(p_results3000,temp1[[1]])
  p_results3000Y <- rbind(p_results3000Y,temp1[[2]])
  p_results5000<-rbind(p_results5000,temp2[[1]])
  p_results5000Y <- rbind(p_results5000Y,temp2[[2]])
}#Test about if part of the routes have much lower accuracy

#n = 3000
#opt %
colMeans(p_results3000) 
sapply(p_results3000,sd)*100 #sd

#E(Y|gopt) - E(Y|gobs)
mean(p_results3000Y$Y.pred-p_results3000Y$Y.observe)
mean(p_results3000Y$Y.pred_naive-p_results3000Y$Y.observe)

sd(p_results3000Y$Y.pred-p_results3000Y$Y.observe) 
sd(p_results3000Y$Y.pred_naive-p_results3000Y$Y.observe)

#% of Recommendation with Nonviable Arm
mean(p_results3000Y$Ares_naive)
mean(p_results3000Y$Ares)

#n = 5000

#%opt
colMeans(p_results5000) 
sapply(p_results5000,sd)*100 #sd

#E(Y|gopt) - E(Y|gobs)
mean(p_results5000Y$Y.pred-p_results5000Y$Y.observe)
mean(p_results5000Y$Y.pred_naive-p_results5000Y$Y.observe)

sd(p_results5000Y$Y.pred-p_results5000Y$Y.observe)
sd(p_results5000Y$Y.pred_naive-p_results5000Y$Y.observe)

#% of Recommendation with Nonviable Arm
mean(p_results5000Y$Ares_naive)
mean(p_results5000Y$Ares)