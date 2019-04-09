Weibull<-function(Germ.Analysis.exp, Test.int=levels(Germ.Analysis.exp$test_type)){
  require(plyr)
  require(agricolae)
  require(gdata)
  require(ggplot2)
  require(reshape)
  if(is.null(Test.int))
    Test.int<-levels(Germ.Analysis.exp$test_type)

  lTi<-list()
  vTia<-c()
  for(i in 1:length(Test.int)){
    lTi[[i]]=which(Germ.Analysis.exp$test_type==Test.int[i])
    if(length(lTi[[1]])!=length(lTi[[i]])){
      stop("Number of replicates for each test must be the same")
    }
    vTia[i]=which(levels(Germ.Analysis.exp$test_type)==Test.int[i])
  }
  mTi <- laply(lTi, function(x) x)
  vTi = sort(as.vector(t(mTi)))
  vTia=sort(vTia)

  test<-Germ.Analysis.exp$test_type[vTi]
  test<-drop.levels(test, reorder = FALSE)
  y<-Germ.Analysis.exp$y[vTia,]
  Dyn_perc<-Germ.Analysis.exp$Dyn_perc_germ[vTia,]
  Ta<- c(1:dim(Dyn_perc)[2])

  sel_int<-which(apply(Dyn_perc, 1, function(x)sum(x!=0)) > 0)
  test_int<-levels(test)[sel_int]
  test_int<-drop.levels(test_int, reorder = FALSE)
  Dyn_perc_int<-Dyn_perc[sel_int,]
  y_int<-y[sel_int,]

  wei_int<-cbind.data.frame(t(y_int), Days=Ta)
  m_wei<-melt(wei_int, id="Days")
  m_wei$variable<-rep(test_int, each=length(wei_int[,1]))
  colnames(m_wei)<-c("Days", "Test", "value")
  DyP_int<-cbind.data.frame(t(Dyn_perc_int), Days=Ta)
  m_DyP<-melt(DyP_int, id="Days")
  m_DyP$variable<-rep(test_int, each=length(wei_int[,1]))
  colnames(m_DyP)<-c("Days", "Test", "value")

  plot_wei<-ggplot(data=m_wei, aes(x=Days, y=value, colour=Test)) +
    geom_line() +
    geom_point(data=m_DyP, aes(x=Days, y=value, colour=Test))+
    scale_y_continuous("Germination %") +
    scale_x_continuous("Days")

  ggsave("plot_Weibull_fun.tiff", plot = plot_wei, width = 174, height = 98, units = c("mm"),dpi = 300)
  plot(plot_wei)

  Weibull<-list(y, levels(test), plot_wei)
  names(Weibull)=c("Weibull", "Test", "plot_Weibull")

  return(Weibull)
}
