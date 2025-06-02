# Fonction qui calcule 1) les medianes a(T,l,s) avec poids, 2) le minimum s(T,l) de s to f(T,a(T,l,s),l,s), 
#3)le minimum hat{l} de l\in [1,h]\to s(T,l), a(T,l,s(T,l), f(T,a(T,l,s(T,l)),l,s(T,l)) et le minimum hat{l} de f.  
Tmin<-function(data,T_initial,h,tau_1,tau_2,tau_3){
data
##objective function a,l,s
objf<-function(x,l,T){
  return(sum(abs(x[1]*exp(-((1:T)-l)^2/(x[2])^2)-data[1:T]))/T)
}
#mediane avec poids a=a(T,l,s)
med_weighted<-function(l,s,T){
  x<-data[1:T]*exp(((1:T)-l)^2/s^2)
  o<-order(x)
  w<-exp(-((1:T)-l)^2/s^2)
  wo<-w[o[1:T]]
  W<-sum(w)
  y<-c(rep(0,T))
  for (t in 1:T){
    if (sum(wo[1:t])<=(W/2)){y[t]<-t}
  }
  return(x[o[max(y)+1]])
} 
#Objectif medweighted x=(l,s)
objfw<-function(l,s,T){
  return(objf(c(med_weighted(l,s,T),s),l,T))
} 
#################################################################################
hat_l<-function(T){
  out=matrix(0,h,3)
  f <- function(s) return(objfw(l,s,T))
  for (l in 1:h){  
    smin <- optimize(f, c(1,20), tol = 0.0001)
    smin
    out[l,]=c(med_weighted(l,smin$minimum,T),smin$minimum,smin$objective)
  }  
  return(list(out,which.min(out[,3]))) 
}
#############################################################################
###########################################################################
premier_indice <- function(x) {
  N <- length(x)
  for (i in 1:(N-tau_1)) {
    if(all(abs(x[i:(i+tau_1)]-x[i])<=tau_2)) {
      return(i)
    }
  } 
  return(NA)  # on Retourne NA si aucun indice ne satisfait la condition
} 
#############################################################################
deuxieme_indice <- function(x) {
  N <- length(x)
  for (i in (premier_indice(x)+1):(N-tau_1)) {
    if(all(abs(x[i:(N-tau_1)]-x[i])<=tau_3)) {
      return(i)
    }
  } 
  return(NA)  # on Retourne NA si aucun indice ne satisfait la condition
}   
#############################################################################
x=rep(0,h-T_initial+1)
for (T in T_initial:h){x[T-T_initial+1]=hat_l(T)[[2]]}
return(list(x,T_initial-1+premier_indice(x),T_initial-1+deuxieme_indice(x)))
}
####################################################################
####################################################################
data=.....
T_initial=10;h=40;tau_1=3;tau_2=tau_3=1;print(Tmin(data,T_initial,h,tau_1,tau_2,tau_3))
