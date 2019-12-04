Germination.niche<-function(germ_dyn, Nv.seed=NULL, n.seed=20, Env_var=NULL){
  if(is.character(germ_dyn)){
    tot_test_no<-read.csv(germ_dyn, sep=";", h=T)
  }
  else{
    tot_test_no<-germ_dyn
  }
  
  if(is.null(Env_var)) stop("error you must specify the column number where the environmental variable is located", call. = FALSE)
  
  
  Env<-germ_dyn[,Env_var]
  if(!is.numeric(Env)) stop("the environmental variable must be numeric", call. = FALSE)
  
  tot_test<-tot_test_no[order(tot_test_no[,1]),]
  
  test_type<-tot_test[,1]
  
  if(is.null(Nv.seed)){
    dyn<-tot_test[,-c(1, Env_var)]
    vital_seed<-rep(n.seed, NROW(tot_test))
    
  }else if(!is.null(Nv.seed)) {
    if(length(Nv.seed)==nrow(tot_test)) {
      vital_seed<-rep(Nv.seed, NROW(tot_test))
      dyn<-tot_test[,-c(1, Env_var)]
    }else if(length(Nv.seed)==1)  {
      vital_seed<-tot_test[,Nv.seed]
      dyn<-tot_test[,-c(1, Nv.seed, Env_var)]
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
  
  Dyn_perc_germ<-aggregate(Dyn_perc_germ_SP, by=list(test_type), FUN=mean,na.rm=TRUE, na.action=NULL)[,-1]
  
  M<-aggregate(M_SP, by=list(test_type), FUN=mean,na.rm=TRUE, na.action=NULL)[,-1]
  names(M)<-levels(test_type)
  
  M_SD<-aggregate(M_SP, by=list(test_type), FUN=sd, na.rm=TRUE)[,-1]
  names(M)<-levels(test_type)
  
  Z<-aggregate(Z_SP, by=list(test_type), FUN=mean,na.rm=TRUE, na.action=NULL)[,-1]
  names(Z)<-levels(test_type)
  
  Z_SD<-aggregate(Z_SP, by=list(test_type), FUN=sd, na.rm=TRUE)[,-1]
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
  
  for(i in 1:length(M_SP)){
    T50_SP[i]=T1[i]+ ((((M_SP[i]/2)-M1[i])*(T2[i]-T1[i]))/(M2[i]-M1[i]))
  }
  
  T50<-aggregate(T50_SP, by=list(test_type), FUN=mean, na.rm=TRUE, na.action=NULL)[,-1]
  
  T50_SD<-aggregate(T50_SP, by=list(test_type), FUN=sd, na.rm=TRUE)[,-1]
  
  M_env<-cbind.data.frame(M_SP,  Env)
  
  Boxplot_M_env<-ggplot(M_env, aes(x=Env, y=M_SP)) +
    geom_boxplot(aes(group=Env, y = M_SP), col="black")+
    geom_smooth(aes(x=Env, y = M_SP), formula=y~x+I(x^2),method = "lm")+
    scale_y_continuous("Final Germination %") +
    scale_x_continuous("Environmental variable")
  
  M_env_lm<-summary(lm(M_SP~Env+I(Env^2), M_env))  
  
ggsave("Boxplot_M_env.tiff", plot = Boxplot_M_env, width = 174, height = 98, units = c("mm"),dpi = 300)
plot(Boxplot_M_env)
  
Z_env<-cbind.data.frame(Z_SP,  Env)
Z_env<-Z_env[-which(M_env==0), ]
Boxplot_Z_env<-ggplot(Z_env, aes(x=Env, y=Z_SP)) +
  geom_boxplot(aes(group=Env, y = Z_SP), col="black")+
  geom_smooth(aes(x=Env, y = Z_SP), formula=y~x+I(x^2),method = "lm")+
  scale_y_continuous("Germination delay - Days") +
  scale_x_continuous("Environmental variable")

Z_env_lm<-summary(lm(Z_SP~Env+I(Env^2), Z_env))

ggsave("Boxplot_Z_env.tiff", plot = Boxplot_Z_env, width = 174, height = 98, units = c("mm"),dpi = 300)
plot(Boxplot_Z_env)

T50_env<-cbind.data.frame( T50_SP,  Env)
T50_env<-T50_env[-which(M_env==0), ]
Boxplot_T50_env<-ggplot(T50_env, aes(x=Env, y= T50_SP)) +
  geom_boxplot(aes(group=Env, y =  T50_SP), col="black")+
  geom_smooth(aes(x=Env, y =  T50_SP), formula=y~x+I(x^2),method = "lm")+
  scale_y_continuous("T50 - Days") +
  scale_x_continuous("Environmental variable")

T50_env_lm<-summary(lm(T50_SP~Env+I(Env^2), T50_env))
ggsave("Boxplot_T50_env.tiff", plot = Boxplot_T50_env, width = 174, height = 98, units = c("mm"),dpi = 300)
plot(Boxplot_T50_env)

  niche_value<-list(M_env, M_env_lm, Boxplot_M_env, Z_env, M_env_lm, Boxplot_Z_env, T50_env, T50_env_lm, Boxplot_T50_env)
  names(niche_value)=c("M_env", "M_env_lm", "Plot_M_env", "Z_env", "M_env_lm" , "Plot_Z_env", "T50_env","T50_env_lm", "Plot_T50_env")
  return(niche_value)
}
