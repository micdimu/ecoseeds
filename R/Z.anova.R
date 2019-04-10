Z.anova<-function(Germ.Analysis.exp, Test.int=levels(Germ.Analysis.exp$test_type), colour="yes"){
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

  Z<-Germ.Analysis.exp$Z_SP[vTi]
  Fin_germ<-Germ.Analysis.exp$M_SP[vTi]
  test<-Germ.Analysis.exp$test_type[vTi]
  test<-drop.levels(test, reorder = FALSE)

  sel_int<-which(Fin_germ > 0)
  test_int<-test[sel_int]
  test_int<-drop.levels(test_int, reorder = FALSE)
  Z_int<-Z[sel_int]

  A1<-aov(Z_int ~ test_int)
  a<-HSD.test(A1, 'test_int')

  T_Z_ep<-cbind.data.frame(x=sort(test_int), y=Z_int)
  sign <- a$groups[order(row.names(a$groups)),]
  sign_t<-cbind.data.frame(sign=sign$groups, levels(test_int))

  if(colour[1]=="yes"){Legend=T_Z_ep$x

  Boxplot_Z<-ggplot(T_Z_ep, aes(x=x, y=y)) +
    geom_boxplot(aes(group=x, y = y, col=Legend))+
    geom_text(data=sign_t, aes(x=c(1:nlevels(test_int)), label = sign,
                               y=aggregate(Z_int, by=list(test_int), max)[,2]+1),
              size=5)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_y_continuous("Z Days") +
    scale_x_discrete("Test")

  }  else if(is.vector(colour) & length(colour)>1){Legend=colour
  Legend<-Legend[sel_int]
  Boxplot_Z<-ggplot(T_Z_ep, aes(x=x, y=y)) +
    geom_boxplot(aes(group=x, y = y, col=Legend))+
    geom_text(data=sign_t, aes(x=c(1:nlevels(test_int)), label = sign,
                               y=aggregate(Z_int, by=list(test_int), max)[,2]+1),
              size=5)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_y_continuous("Z Days") +
    scale_x_discrete("Test")

  } else if(colour[1]=="no"){

  Boxplot_Z<-ggplot(T_Z_ep, aes(x=x, y=y)) +
    geom_boxplot(aes(group=x, y = y))+
    geom_text(data=sign_t, aes(x=c(1:nlevels(test_int)), label = sign,
                               y=aggregate(Z_int, by=list(test_int), max)[,2]+1),
              size=5)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_y_continuous("Z Days") +
    scale_x_discrete("Test")
  }

  ggsave("Boxplot_Z.anova.tiff", plot = Boxplot_Z, width = 174, height = 98, units = c("mm"),dpi = 300)
  plot(Boxplot_Z)


  Z.V<-list(summary(A1), a, levels(test), Z, Boxplot_Z)
  names(Z.V)=c("Anova", "Tukey", "Test", "Z", "Boxplot_Z")

  if(any(Fin_germ==0)){warning("test with final germination % = 0 are exluded")}
  return(Z.V)
}
