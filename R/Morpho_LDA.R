Morpho_LDA<-function(Morpho_df, N.iter=100, split_size=0.80){
require(plyr)
require(MASS)
require(ggplot2)
require(plotly)
require(devtools)

  if(split_size>1) stop ("split_size must be % from 0 to 1")

  stratified<-function(df, group, size, select = NULL,
           replace = FALSE, bothSets = FALSE) {
    if (is.null(select)) {
      df <- df
    } else {
      if (is.null(names(select))) stop("'select' must be a named list")
      if (!all(names(select) %in% names(df)))
        stop("Please verify your 'select' argument")
      temp <- sapply(names(select),
                     function(x) df[[x]] %in% select[[x]])
      df <- df[rowSums(temp) == length(select), ]
    }
    df.interaction <- interaction(df[group], drop = TRUE)
    df.table <- table(df.interaction)
    df.split <- split(df, df.interaction)
    if (length(size) > 1) {
      if (length(size) != length(df.split))
        stop("Number of groups is ", length(df.split),
             " but number of sizes supplied is ", length(size))
      if (is.null(names(size))) {
        n <- setNames(size, names(df.split))
        message(sQuote("size"), " vector entered as:\n\nsize = structure(c(",
                paste(n, collapse = ", "), "),\n.Names = c(",
                paste(shQuote(names(n)), collapse = ", "), ")) \n\n")
      } else {
        ifelse(all(names(size) %in% names(df.split)),
               n <- size[names(df.split)],
               stop("Named vector supplied with names ",
                    paste(names(size), collapse = ", "),
                    "\n but the names for the group levels are ",
                    paste(names(df.split), collapse = ", ")))
      }
    } else if (size < 1) {
      n <- round(df.table * size, digits = 0)
    } else if (size >= 1) {
      if (all(df.table >= size) || isTRUE(replace)) {
        n <- setNames(rep(size, length.out = length(df.split)),
                      names(df.split))
      } else {
        message(
          "Some groups\n---",
          paste(names(df.table[df.table < size]), collapse = ", "),
          "---\ncontain fewer observations",
          " than desired number of samples.\n",
          "All observations have been returned from those groups.")
        n <- c(sapply(df.table[df.table >= size], function(x) x = size),
               df.table[df.table < size])
      }
    }
    temp <- lapply(
      names(df.split),
      function(x) df.split[[x]][sample(df.table[x],
                                       n[x], replace = replace), ])
    set1 <- do.call("rbind", temp)

    if (isTRUE(bothSets)) {
      set2 <- df[!rownames(df) %in% rownames(set1), ]
      list(SET1 = set1, SET2 = set2)
    } else {
      set1
    }
  }

  confusion <- function(actual, predicted, names = NULL, printit = FALSE,
                        prior = NULL) {
    if (is.null(names))
      names <- levels(actual)
    tab <- table(actual, predicted)
    acctab <- t(apply(tab, 1, function(x) x/sum(x)))
    dimnames(acctab) <- list(Actual = names, "Predicted (cv)" = names)
    if (is.null(prior)) {
      relnum <- table(actual)
      prior <- relnum/sum(relnum)
      acc <- sum(tab[row(tab) == col(tab)])/sum(tab)
    }
    else {
      acc <- sum(prior * diag(acctab))
      names(prior) <- names
    }
    if (printit)
      print(round(c("Overall accuracy" = acc, "Prior frequency" = prior),
                  4))
    if (printit) {
      cat("\nConfusion matrix", "\n")
      print(round(acctab, 4))
    }
    invisible(acctab)
  }

if(is.character(Morpho_df)){morpho<-read.csv(Morpho_df, sep=";", h=T)}
else{morpho<-Morpho_df}
colnames(morpho)[1]<-"group"

mean_train<-c()
mean_test<-c()
conf_train<-list()
conf_test<-list()

for(i in 1:N.iter){

  morpho_train<-stratified(morpho, "group", size=split_size)



  morpho_test<-morpho[-as.numeric(row.names(morpho_train)),]

  morpho.lda <- lda(group ~ ., data =morpho_train)


  pred.mat.train <- predict(morpho.lda, morpho_train)$class
  pred.mat.test <- predict(morpho.lda, morpho_test)$class
  #accuracy on training data
  mean_train[i]<-mean(pred.mat.train == morpho_train$group)
  mean_test[i]<-mean(pred.mat.test == morpho_test$group)
  conf_train[[i]]<-confusion(morpho_train$group ,pred.mat.train )
  conf_test[[i]]<-confusion(morpho_test$group ,pred.mat.test )
}

mtrain<-mean(mean_train)
mtest<-mean(mean_test)

conf_train_mean<-aaply(laply(conf_train, as.matrix), c(2, 3), mean)
conf_test_mean<-aaply(laply(conf_test, as.matrix), c(2, 3), mean)

morpho.lda.values <- predict(morpho.lda)
owd<-data.frame(Legend=morpho.lda.values$class, morpho.lda.values$x)


plot_lda<-ggplot(owd, aes(x=LD1, y=LD2, col=Legend)) + geom_point( size = 1.5, aes(color = Legend))+
  stat_ellipse(type = "t")


if(nlevels(as.factor(morpho$group))>3){
p<-plot_ly(data=owd, x=~LD1, y=~LD2, z=~LD3,  color=~Legend,colors="YlOrRd", type="scatter3d", mode ="markers" )
show(p)}


ggsave("morpho_lda.tiff", plot = plot_lda, width = 174, height = 98, units = c("mm"),dpi = 300)
plot(plot_lda)


lda_group<-list(morpho.lda,  mtrain, mtest, conf_train_mean, conf_test_mean, plot_lda)
names(lda_group)=c("morpho.lda", "accuracy_train", "accuracy_test", "confusion_train", "confusion_test", "plot_lda")

return(lda_group)

}
