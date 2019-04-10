#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

Germination.Analysis<-function(germ_dyn, Nv.seed=NULL, cv=1.5, n.seed=20){
  if(is.character(germ_dyn)){
    tot_test_no<-read.csv(germ_dyn, sep=";", h=T)
  }
  else{
    tot_test_no<-germ_dyn
  }

  tot_test<-tot_test_no[order(tot_test_no[,1]),]

  test_type<-tot_test[,1]

  if(is.null(Nv.seed)){
    dyn<-tot_test[,-1]
    vital_seed<-rep(n.seed, NROW(tot_test))
  }else if(!is.null(Nv.seed)) {
    if(length(Nv.seed)==nrow(tot_test)) {
      vital_seed<-rep(Nv.seed, NROW(tot_test))
      dyn<-tot_test[,-c(1)]
    }else if(length(Nv.seed)==1)  {
      vital_seed<-tot_test[,Nv.seed]
      dyn<-tot_test[,-c(1,Nv.seed)]
    } else  stop("error Nv.seed must be a vector of vital seed for all the replicates lenght(Nv.seed)==lenght(germ_dyn[,1]) or the number of column in germ_dyn that contain N. of vital seeds", call. = FALSE)
  }

  Dyn_perc_germ_SP<-matrix(nrow=NROW(dyn), ncol=NCOL(dyn))
  M_SP<-c()
  Z_SP<-c()

  for(i in 1:NROW(dyn)){

    Dyn_perc_germ_SP[i,]<-as.matrix((dyn[i,]*100)/vital_seed[i])

    M_SP[i]<-tail(Dyn_perc_germ_SP[i,], n=1)

    Z_SP[i]<-length(which(Dyn_perc_germ_SP[i,]==0))+1

  }

  # aggregate ordina in ordine alfabetico
  Dyn_perc_germ<-aggregate(Dyn_perc_germ_SP, by=list(test_type), FUN=mean)[,-1]

  M<-aggregate(M_SP, by=list(test_type), FUN=mean)[,-1]
  names(M)<-levels(test_type)

  M_SD<-aggregate(M_SP, by=list(test_type), FUN=sd)[,-1]
  names(M)<-levels(test_type)

  Z<-aggregate(Z_SP, by=list(test_type), FUN=mean)[,-1]
  names(Z)<-levels(test_type)

  Z_SD<-aggregate(Z_SP, by=list(test_type), FUN=sd)[,-1]
  names(Z)<-levels(test_type)

  T1<-c()
  M1<-c()
  T2<-c()
  M2<-c()

  for(i in 1:length(M_SP)){
    if((M_SP[i]/2)==0){
      T1[i]=1
      M1[i]= Dyn_perc_germ_SP[i,T1[i]]
      T2[i]=1
      M2[i]= Dyn_perc_germ_SP[i,T2[i]]
    } else{
      T1[i]= tail(which(Dyn_perc_germ_SP[i,] < (M_SP[i]/2)), n=1)
      M1[i]= Dyn_perc_germ_SP[i,T1[i]]
      T2[i]= head(which(Dyn_perc_germ_SP[i,] > (M_SP[i]/2)), n=1)
      M2[i]= Dyn_perc_germ_SP[i,T2[i]]
    }}

  T50_SP<-c()
  k_SP<-c()
  for(i in 1:length(M_SP)){
    T50_SP[i]= (((M_SP[i]/2)-M1[i])*T2[i]-T1[i])/(M2[i]-M1[i])+T1[i]
    k_SP[i]<- 1/(T50_SP[i]-Z_SP[i])
  }

  T50<-aggregate(T50_SP, by=list(test_type), FUN=mean)[,-1]
  k<-aggregate(k_SP, by=list(test_type), FUN=mean)[,-1]

  T50_SD<-aggregate(T50_SP, by=list(test_type), FUN=sd)[,-1]
  k_SD<-aggregate(k_SP, by=list(test_type), FUN=sd)[,-1]

  Ta<- c(1:dim(dyn)[2])

  if(is.null(cv)){cv<-c(1.5)}

  c<-cv      ##### vedere fattore

  y<-matrix(nrow=length(M), ncol=NCOL(dyn))
  for(i in 1:length(M)){
    for(j in 1:length(Ta)){
      if(Dyn_perc_germ[i,j]<=0){
        y[i,j]=0}
      else{
        y[i,j]= M[i] * (1 - (exp(1) ^ (- (k[i]*(Ta[j]-Z[i]))^c) ))
      }
    }
  }

  y[is.nan(y)]=0
  germ_value<-list(test_type,   Dyn_perc_germ_SP,   M_SP,   Z_SP,   T50_SP,  Dyn_perc_germ,   M,   M_SD,    Z,   Z_SD,   T50,   T50_SD,  k,   y)
  names(germ_value)=c("test_type", "Dyn_perc_germ_SP", "M_SP", "Z_SP", "T50_SP","Dyn_perc_germ", "M", "M_SD",  "Z", "Z_SD", "T50", "T50_SD", "k", "y")
  return(germ_value)
}
