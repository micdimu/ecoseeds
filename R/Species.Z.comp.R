Species.Z.comp<-function(Germ.Analysis.exp_sp1, Germ.Analysis.exp_sp2, sp_name=NULL , colour="yes", Test.int=NULL){
  require(plyr)
  require(agricolae)
  require(gdata)
  require(ggplot2)

  sel_M0_sp1<-which(Germ.Analysis.exp_sp1$M_SP > 0)
  sel_M0_sp2<-which(Germ.Analysis.exp_sp2$M_SP > 0)

  test_M0_sp1<-Germ.Analysis.exp_sp1$test_type[sel_M0_sp1]
  test_M0_sp1<-drop.levels(test_M0_sp1, reorder = FALSE)

  test_M0_sp2<-Germ.Analysis.exp_sp2$test_type[sel_M0_sp2]
  test_M0_sp2<-drop.levels(test_M0_sp2, reorder = FALSE)

  Z_M0_sp1<-Germ.Analysis.exp_sp1$Z_SP[sel_M0_sp1]
  Z_M0_sp2<-Germ.Analysis.exp_sp2$Z_SP[sel_M0_sp2]


  int_comm_test<-intersect(test_M0_sp1, test_M0_sp2)

  if(length(int_comm_test)==0) stop('No tests matching between taxa')

  sel_comm_test_sp1<-test_M0_sp1 %in% int_comm_test
  sel_comm_test_sp2<-test_M0_sp2 %in% int_comm_test

  comm_test_sp1<-test_M0_sp1[which(sel_comm_test_sp1==TRUE)]
  comm_test_sp2<-test_M0_sp2[which(sel_comm_test_sp2==TRUE)]

  comm_Z_sp1<-Z_M0_sp1[which(sel_comm_test_sp1==TRUE)]
  comm_Z_sp2<-Z_M0_sp2[which(sel_comm_test_sp2==TRUE)]

  comm_sp1<-cbind.data.frame(Z=comm_Z_sp1,  test=comm_test_sp1)
  comm_sp2<-cbind.data.frame(Z=comm_Z_sp2,  test=comm_test_sp2)

  if(is.null(sp_name)){
    comm_sp1$sp_name<-rep(c("sp_1"), nrow(comm_sp1))
    comm_sp2$sp_name<-rep(c("sp_2"), nrow(comm_sp2))
  } else {
    comm_sp1$sp_name<-rep(sp_name[1], nrow(comm_sp1))
    comm_sp2$sp_name<-rep(sp_name[2], nrow(comm_sp2))
  }

  comm_sp1_sp2<-rbind.data.frame(comm_sp1, comm_sp2)

  comm_sp1_sp2$test<-drop.levels(comm_sp1_sp2$test, reorder = FALSE)



  if(is.null(Test.int))
    Test.int<-levels(comm_sp1_sp2$test)
  lTi<-list()
  for(i in 1:length(Test.int)){
    lTi[[i]]=which(comm_sp1_sp2$test==Test.int[i])
  }

  vTi <- sort(do.call(c, lTi))

  comm_sp1_sp2<-comm_sp1_sp2[vTi,]
  comm_sp1_sp2$test<-drop.levels(comm_sp1_sp2$test, reorder = FALSE)

  comm_sp1_sp2$test_sp<-interaction(comm_sp1_sp2$sp_name,comm_sp1_sp2$test)

  A1<-aov(Z ~ test_sp, data=comm_sp1_sp2)
  a<-HSD.test(A1, 'test_sp')

  sign <- a$groups[order(row.names(a$groups)),]
  sign_t<-cbind.data.frame(sign=sign$groups, levels(comm_sp1_sp2$test_sp))

  if(colour[1]=="yes"){Legend=comm_sp1_sp2$test_sp

  Legend1=comm_sp1_sp2$test_sp

  Boxplot_com_spe<-ggplot(comm_sp1_sp2, aes(x=Legend1, y = Z)) +
    geom_boxplot(aes(group=Legend1, y = Z, col=Legend))+
    geom_text(data=sign_t, aes(x=c(1:nlevels(Legend1)), label = sign,
                               y=aggregate(comm_sp1_sp2$Z, by=list(Legend1), max)[,2]+5),size=5)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_y_continuous("Germination delay") +
    scale_x_discrete("Test")

  }  else if(is.vector(colour) & length(colour)>1){
    Legend1=comm_sp1_sp2$test_sp
    Legend=colour

    Boxplot_com_spe<-ggplot(comm_sp1_sp2, aes(x=Legend1, y = Z)) +
      geom_boxplot(aes(group=Legend1, y = Z, col=Legend))+
      geom_text(data=sign_t, aes(x=c(1:nlevels(Legend1)), label = sign,
                                 y=aggregate(comm_sp1_sp2$Z, by=list(Legend1), max)[,2]+5),size=5)+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      scale_y_continuous("Germination delay") +
      scale_x_discrete("Test")

  } else if(colour[1]=="no"){
    Legend1=comm_sp1_sp2$test_sp

    Boxplot_com_spe<-ggplot(comm_sp1_sp2, aes(x=Legend1, y = Z)) +
      geom_boxplot(aes(group=Legend1, y = Z))+
      geom_text(data=sign_t, aes(x=c(1:nlevels(Legend1)), label = sign,
                                 y=aggregate(comm_sp1_sp2$Z, by=list(Legend1), max)[,2]+5),size=5)+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      scale_y_continuous("Germination delay") +
      scale_x_discrete("Test")
  }
  ggsave("Boxplot_compare_species_Z.tiff", plot = Boxplot_com_spe, width = 174, height = 98, units = c("mm"),dpi = 300)
  plot(Boxplot_com_spe)


  Com.spe.G.V<-list(summary(A1), a, levels(comm_sp1_sp2$test), comm_sp1_sp2$Z, Boxplot_com_spe)
  names(Com.spe.G.V)=c("Anova", "Tukey", "Test", "Germination_delay", "Boxplot_com_spe_Z")

  if(any(Germ.Analysis.exp_sp1$M_SP==0)){warning("test with final germination % = 0 are exluded")}

  return(Com.spe.G.V)
}
