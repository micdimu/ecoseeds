Final.Germ.anova<-function(Germ.Analysis.exp, Test.int=levels(Germ.Analysis.exp$test_type), colour="yes"){
  require(plyr)
  require(agricolae)
  require(gdata)
  require(ggplot2)

  if(is.null(Test.int))
    Test.int<-levels(Germ.Analysis.exp$test_type)

  lTi<-list()
  for(i in 1:length(Test.int)){
    lTi[[i]]=which(Germ.Analysis.exp$test_type==Test.int[i])
    if(length(lTi[[1]])!=length(lTi[[i]])){
      stop("Number of replicates for each test must be the same")
    }
  }
  mTi <- laply(lTi, function(x) x)
  vTi = sort(as.vector(t(mTi)))

  Fin_germ<-Germ.Analysis.exp$M_SP[vTi]
  test<-Germ.Analysis.exp$test_type[vTi]
  test<-drop.levels(test, reorder = FALSE)
  A1<-aov(Fin_germ ~ test)
  a<-HSD.test(A1, 'test')
  T_fg_ep<-cbind.data.frame(x=sort(test), y=Fin_germ)
  sign <- a$groups[order(row.names(a$groups)),]
  sign_t<-cbind.data.frame(sign=sign$groups, levels(test))

  if(colour[1]=="yes"){Legend=T_fg_ep$x

  Boxplot_Fin.Germ<-ggplot(T_fg_ep, aes(x=x, y=y)) +
    geom_boxplot(aes(group=x, y = y, col=Legend))+
    geom_text(data=sign_t, aes(x=c(1:nlevels(test)), label = sign,
                               y=aggregate(Fin_germ, by=list(test), max)[,2]+5),
              size=5)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_y_continuous("Final Germination %") +
    scale_x_discrete("Test")

  }  else if(is.vector(colour) & length(colour)>1){Legend=colour

  Boxplot_Fin.Germ<-ggplot(T_fg_ep, aes(x=x, y=y)) +
    geom_boxplot(aes(group=x, y = y, col=Legend))+
    geom_text(data=sign_t, aes(x=c(1:nlevels(test)), label = sign,
                               y=aggregate(Fin_germ, by=list(test), max)[,2]+5),
              size=5)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_y_continuous("Final Germination %") +
    scale_x_discrete("Test")
  } else if(colour[1]=="no"){

  Boxplot_Fin.Germ<-ggplot(T_fg_ep, aes(x=x, y=y)) +
    geom_boxplot(aes(group=x, y = y))+
    geom_text(data=sign_t, aes(x=c(1:nlevels(test)), label = sign,
                               y=aggregate(Fin_germ, by=list(test), max)[,2]+5),
              size=5)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_y_continuous("Final Germination %") +
    scale_x_discrete("Test")
}
  ggsave("Boxplot_Final.germ.anova.tiff", plot = Boxplot_Fin.Germ, width = 174, height = 98, units = c("mm"),dpi = 300)
  plot(Boxplot_Fin.Germ)


  Fin.G.V<-list(summary(A1), a, levels(test), Fin_germ, Boxplot_Fin.Germ)
  names(Fin.G.V)=c("Anova", "Tukey", "Test", "Fin_germ", "Boxplot_Fin.Germ")

  return(Fin.G.V)
}
