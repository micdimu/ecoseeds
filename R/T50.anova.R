T50.anova<-function(Germ.Analysis.exp, Test.int=levels(Germ.Analysis.exp$test_type), colour="yes"){
  require(plyr)
  require(agricolae)
  require(gdata)
  require(ggplot2)

  if(is.null(Test.int))
    Test.int<- levels(Germ.Analysis.exp$test_type)

  lTi<-list()
  for(i in 1:length(Test.int)){
    lTi[[i]]=which(Germ.Analysis.exp$test_type==Test.int[i])
    if(length(lTi[[1]])!=length(lTi[[i]])){
      stop("Number of replicates for each test must be the same")
    }
  }
  mTi <- laply(lTi, function(x) x)
  vTi = sort(as.vector(t(mTi)))

  T50<-Germ.Analysis.exp$T50_SP[vTi]
  Fin_germ<-Germ.Analysis.exp$M_SP[vTi]
  test<-Germ.Analysis.exp$test_type[vTi]
  test<-drop.levels(test, reorder = FALSE)

  sel_int<-which(Fin_germ > 0)
  test_int<-test[sel_int]
  test_int<-drop.levels(test_int, reorder = FALSE)
  T50_int<-T50[sel_int]

  A1<-aov(T50_int ~ test_int)
  a<-HSD.test(A1, 'test_int')

  T_T50_ep<-cbind.data.frame(x=sort(test_int), y=T50_int)
  sign <- a$groups[order(row.names(a$groups)),]
  sign_t<-cbind.data.frame(sign=sign$groups, levels(test_int))

  if(colour[1]=="yes"){Legend=T_T50_ep$x
  Boxplot_T50<-ggplot(T_T50_ep, aes(x=x, y=y)) +
    geom_boxplot(aes(group=x, y = y, col=Legend))+
    geom_text(data=sign_t, aes(x=c(1:nlevels(test_int)), label = sign,
                               y=aggregate(T50_int, by=list(test_int), max)[,2]+2),
              size=5)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_y_continuous("T50 Days") +
    scale_x_discrete("Test")

  }  else if(is.vector(colour) & length(colour)>1){Legend=colour
  Legend<-Legend[sel_int]
  Boxplot_T50<-ggplot(T_T50_ep, aes(x=x, y=y)) +
    geom_boxplot(aes(group=x, y = y, col=Legend))+
    geom_text(data=sign_t, aes(x=c(1:nlevels(test_int)), label = sign,
                               y=aggregate(T50_int, by=list(test_int), max)[,2]+2),
              size=5)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_y_continuous("T50 Days") +
    scale_x_discrete("Test")

  } else if(colour[1]=="no"){

  Boxplot_T50<-ggplot(T_T50_ep, aes(x=x, y=y)) +
    geom_boxplot(aes(group=x, y = y))+
    geom_text(data=sign_t, aes(x=c(1:nlevels(test_int)), label = sign,
                               y=aggregate(T50_int, by=list(test_int), max)[,2]+2),
              size=5)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_y_continuous("T50 Days") +
    scale_x_discrete("Test")
  }

  ggsave("Boxplot_T50.anova.tiff", plot = Boxplot_T50, width = 174, height = 98, units = c("mm"),dpi = 300)
  plot(Boxplot_T50)

  T50.V<-list(summary(A1), a, levels(test), T50, Boxplot_T50)
  names(T50.V)=c("Anova", "Tukey", "Test", "T50", "Boxplot_T50")

  if(any(Fin_germ==0)){warning("test with final germination % = 0 are exluded")}
  return(T50.V)
}
