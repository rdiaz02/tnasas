####  Copyright (C) 2004, 2005  Ramón Díaz-Uriarte

#### This program is free software; you can redistribute it and/or
#### modify it under the terms of the GNU General Public License
#### as published by the Free Software Foundation; either version 2
#### of the License, or (at your option) any later version.

#### This program is distributed in the hope that it will be useful,
#### but WITHOUT ANY WARRANTY; without even the implied warranty of
#### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#### GNU General Public License for more details.

#### You should have received a copy of the GNU General Public License
#### along with this program; if not, write to the Free Software
#### Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.




## In the future:
##            parallelization? Probably unneeded
##   FUTURE:
##   an R package; include returning function call
##   add bootstrap to evaluate stability of rules (sets of predictors).
##   otros tipos de SVM



##  zz: watch out: use and fix for ever a multtest package.
.First.lib <- function(lib, pkg) {
    library.dynam("Tnasas", pkg, lib) ## beware, must be last
    require(supclust)
    require(e1071) ## for svm
    require(class) ## for knn
    require(multtest) ## for gene selection
    ##dyn.load("dldacpp.so") ## for my dlda
    require(pamr)
    require(randomForest)
    require(Hmisc)
    library.dynam("Tnasas", pkg, lib) ## beware, must be last
    ## move to Namespace and then import only the needed functions
}


dldaC <- function (ls, cll, ts) {
    if(!is.numeric(cll)) stop("cll has to be numeric, starting at 0")
    if(min(cll) != 0) stop("cll has to be numeric, starting at 0")
    if(!is.matrix(ls)) ls <- as.matrix(ls)
    if(!is.matrix(ts)) ts <- as.matrix(ts)
    .C("dlda", as.double(ls), as.integer(cll),
       as.integer(as.vector(table(cll))), as.double(t(ts)),
       as.integer(max(cll) - min(cll) + 1), as.integer(length(cll)),
       as.integer(ncol(ls)), as.integer(nrow(ts)),
       predictions = as.integer(rep(-9, nrow(ts))))$predictions
}
 
###select.k <- function(data, class, max.k = 20) {
##### For knn, choose best number of neighbors by cross-validations
##### But no longer used
###    error.k <- rep(-9, max.k)
###    for (i in 1:max.k) {
###        error.k[i] <- length(which(class != knn.cv(data, class, k = i)))
###    }
###    print(which.min(error.k))
###    return(which.min(error.k))
###}



## F statistic for gene selection
geneSelect.F <- function(data, class) {
    class <- as.numeric(factor(class)) - 1
    tmp   <- mt.maxT(t(data), class, B = 1, test = "f")[, c(1, 2)]
    selected <- tmp[, 1]
    return(selected)
}



geneSelect.Wilcoxon <- function(data, class) {
    ## Wilcoxon for gene selection
    ## Using supclust
    print("Using Wilcoxon for gene selection")
    nclass <- length(levels(class))
    N <- length(class)
    ngenes <- dim(data)[2]
    scores <- matrix(NA, ncol = nclass, nrow = ngenes)
    for(i in 1:nclass) {
        newclass <- rep(0, N)
        newclass[class != levels(class)[i]] <- 1
        sf <- sign.flip(data, newclass)$flipped.matrix
        for(j in 1:ngenes) {
            scores[j, i] <- score(sf[, j], newclass)
        }
    }
    best.scores <- apply(scores, 1, min)
    selected <- order(best.scores)
    return(selected)
}


###geneSelect.Wilcoxon <- function(data, class, size) {
###    ## Wilcoxon for gene selection
###    ## Using multtest; desn't work. I don't know what exactly is returned
###    ## by multtest
###    tdata <- t(data)
###    print("Using Wilcoxon for gene selection")
###    nclass <- length(levels(class))
###    N <- length(class)
###    ngenes <- dim(data)[2]
###    scores <- matrix(NA, ncol = nclass, nrow = ngenes)
###    for(i in 1:nclass) {
###        newclass <- rep(0, N)
###        newclass[class != levels(class)[i]] <- 1
###        tmp <- mt.maxT(tdata, newclass, B = 1, test = "wilcoxon")
###        scores[, i] <- tmp$teststat[order(tmp$index)]
###    }
###    best.scores <- apply(scores, 1, max)
###    selected <- order(best.scores, decreasing = TRUE)[1:size]
###    return(selected)
###}



geneSelect.rF <- function(data, class) {
    ## random Forest for gene selection
    ## Using supclust
    print("Using random forest for gene selection")
    rf1 <- randomForest(data, class, importance = TRUE,
                        ntree = 1000)
    imps <- importance(rf1, type = 1, scale = FALSE)
    selected <- order(imps, decreasing = TRUE)
    return(selected)
}

cv.balance <- function(index, cll, kn) {
    ### CV partitions with equal proportions or each group in
    ### training and testing, as in Dettling et al.
    unlist(tapply(index, cll,
                  function(x) sample(rep(1:kn, length = length(x)),
                                     length(x), replace = FALSE)))
    
}

## to see number of each class in testing:
## cl2 <- factor(c(rep("a", 10), rep("b", 2)))
## table(cl2, cv.balance(1:12, cl2, 10))



### To add predictors, etc, only need to create the appropriate
### internalPredictor.whatever

internalPredictor.svmLinear <- function(x, y, test.data) {
    predict(svm(x, y, kernel = "linear"), test.data)
}

internalPredictor.dldaC <- function(x, y, test.data) {
    yn <- as.integer(unclass(y)) - min(as.integer(unclass(y))) 
    predictions.numeric <- dldaC(x, yn, test.data)
    predictions.factor <- factor(predictions.numeric + min(as.integer(unclass(y))),
                                 levels = seq(along = levels(y)),
                                 labels = levels(y))
    predictions.factor
}

###internalPredictor.knn <- function(x, y, test.data, best.k = 1, ...) {
###    if(is.null(best.k))
###        best.k <- select.k(x, y, max.k = 20)
###    knn(x, test.data, y, best.k)
###}
    

###internalPredictor.knn <- function(x, y, test.data, ...) {
###### so we select best neighbor among odd numbers 1 to 15
##### For knn, choose best number of neighbors by cross-validations
###    error.k <- rep(-9, 8)
###    neighs <- seq(from = 1, to = 15, by = 2)
###    for (i in 1:8) {
###        error.k[i] <- length(which(y != knn.cv(x, y, k = neighs[i])))
###    }
###    best.k <- which.min(error.k)
###    print(best.k)
###    knn(x, test.data, y, best.k)
###}


internalPredictor.knn <- function(x, y, test.data) {
    class::knn1(x, test.data, y)
}


internalPredictor.rF <- function(x, y, test.data, ntree = 1000) {
    randomForest(x, y, xtest = test.data, ntree = ntree)$test$predicted
}



print.cvPred <- function(x) {
    cat("\nError rate:", round(x[[1]], 4), "\n\n")
    cat("Confussion Matrix:\n")
    print(x[[2]])
          cat("\nNumber of predictors that yields minimum error rate:", x[[3]], "\n\n")
          cat("Selected predictor genes: \n")
          print(x$geneNames[x[[4]]])
          cat("\nOOB predictions (and true ---observed--- class): \n")
          print(x[[5]])

          wrong.pred <- which(x[[5]][, 1] != x[[5]][, 2])
          if(length(wrong.pred)) {
                  cat("\n\n Cases with incorrect predictions (errors): \n")
                        print(x[[5]][wrong.pred,])
                }

          ks <- length(x$genesSelected.cv)
          num.selected <- unlist(lapply(x$genesSelected.cv, length))

          cv.names <- paste("CV.run.", 1:ks, sep = "")

          cat("\n\n Stability assessments \n")
          cat(    " ---------------------\n")
      ###    cat("\n Number of genes selected in the cross-validation runs \n")
      ###    print(num.selected)
          cat("\n Genes selected in each of the cross-validation runs \n")
          for(i in 1:ks) {
                    cat(paste("CV run  ", i, " (", num.selected[i], " genes selected):   ", sep = ""), "\n")
                            print(x$geneNames[x$genesSelected.cv[[i]] ])
                            cat("\n---\n")
                  }
          shared.genes <- matrix(NA, nrow = (ks + 1), ncol = (ks + 1))
          tmp.genesSelected <- list()
          tmp.genesSelected[[1]] <- x[[4]]
          tmp.genesSelected <- c(tmp.genesSelected, x$genesSelected.cv)
          for(i in 1:(ks + 1)) {
                    for(j in 1:(ks + 1)) { ## sure, need not be symmetric, but this is fast
                                  shared.genes[i, j] <-
                                                    length(intersect(tmp.genesSelected[[i]],
                                                                                                      tmp.genesSelected[[j]]))
                                }
                  }
          prop.shared <- round(shared.genes/c(x[[3]], num.selected), 3)

          num.selected <- c(x[[3]], num.selected)
          num.selectedS <- paste("(", num.selected, ")", sep = "")
          colnames(shared.genes) <- colnames(prop.shared) <- c("OriginalSample", cv.names)
          rownames(shared.genes) <- rownames(prop.shared) <- paste(c("OriginalSample", cv.names), num.selectedS)

      options(width = 200)
      cat("\n\n Number of shared genes \n")
          print(as.table(shared.genes))

          cat("\n\n Proportion of shared genes (relative to row total) \n")
          print(as.table(prop.shared))

      options(width = 80)
          unlisted.genes.selected <- unlist(x$genesSelected.cv)
          named.unlisted.selected <- x$geneNames[unlisted.genes.selected]


      in.all.data <-
        which(names(table(named.unlisted.selected, dnn = NULL)) %in% x$geneNames[x[[4]]])
      cat("\n\n\n Gene freqs. in cross-validated runs of genes selected in model with all data \n\n")
      print(sort(table(named.unlisted.selected, dnn = NULL)[in.all.data], decreasing = TRUE))
      cat("\n")
      
      
      cat("\n\n Gene frequencies in cross-validated runs \n\n")
      tmp.table <- sort(table(named.unlisted.selected, dnn = NULL),
                        decreasing = TRUE)
      ###    n.tmp.table <- names(tmp.table)
      ###    dim(tmp.table) <- c(dim(tmp.table), 1)
      ###    rownames(tmp.table) <- n.tmp.table
      ###    colnames(tmp.table) <- "Number of times"
      print(tmp.table)
      cat("\n")

    }






linkGene <- function(id) {
    ## Returns a string for use in a web page with a call
    ## to IDClight.
    if ((idtype == "None") | (organism == "None"))
        return(id)
    else
        paste("<a href=\"http://idclight.bioinfo.cnio.es/IDClight.prog",
              "?idtype=", idtype, "&id=", id, "&org=",
              organism,"\" target=\"icl_window\" >",id,"</a>", sep = "")
## target=\"icl_window\"\
}
     


linkGene2 <- function(id) {
    ## Returns a string for use in a web page with a call
    ## to IDClight.
    if ((idtype == "None") | (organism == "None"))
        return(id)
    else 
        paste("http://idclight.bioinfo.cnio.es/IDClight.prog",
              "?idtype=", idtype, "&id=", id, "&org=",
              organism, sep = "")
}

html.data.frame <- function (object, first.col = "Name",
                             file = paste(first.word(deparse(substitute(object))), 
                             "html", sep = "."), append = FALSE, link = NULL, linkCol = 1, 
                             linkType = c("href", "name"), ...) 
{
    linkType <- match.arg(linkType)
    x <- format.df(object, ...)
    adj <- attr(x, "col.just")
    if (any(adj == "r")) 
        for (i in seq(along = adj)[adj == "r"]) x[, i] <- paste("<div align=right>", 
            x[, i], "</div>", sep = "")
    if (length(r <- dimnames(x)[[1]])) 
        x <- cbind(first.col = r, x)
    colnames(x)[1] <- first.col
    cat("<TABLE BORDER>\n", file = file, append = append)
    cat("<tr>", paste("<td>", dimnames(x)[[2]], "</td>", sep = ""), 
        "</tr>\n", sep = "", file = file, append = file != "")
    if (length(link)) 
        x[, linkCol] <- ifelse(link == "", x[, linkCol], paste("<a ", 
            linkType, "=\"", link, "\">", x[, linkCol], "</a>", 
            sep = ""))
    for (i in 1:nrow(x)) cat("<tr>", paste("<td>", x[i, ], "</td>", 
        sep = ""), "</tr>\n", sep = "", file = file, append = file != 
        "")
    cat("</TABLE>\n", file = file, append = file != "")
    structure(list(file = file), class = "html")
}







HTML.cvPred <- function (x, file = "results.txt", append = TRUE) 
{
    sink(file, append)
    cat("<h3>Prediction error (using cross-validation)</h3>")
    cat("<p><b>Error rate</b>", round(x[[1]], 4), "</p>")
    cat("<p>Confussion Matrix:<br><pre>")
    print(x[[2]])
    cat("</pre></p>")
    cat("<br>Number of predictors that yields minimum error rate:", 
        x[[3]], "<br>")
    cat("<h3>Selected predictor genes:</h3>")
    cat("<TABLE frame=\"box\">\n")
    cat("<tr><th width=200>Gene name</th></tr>\n")
    if(length(x[[4]] == 1) && (x[[4]] == "All the genes")) {
      cat("<tr><td> All the genes in the data set</td></tr>\n")
    } else {
      gns <- x$geneNames[x[[4]]]
      for (i in 1:length(gns)) {
        cat("<tr><td>", linkGene(gns[i]), "</td></tr>\n")
      }
    }
    cat("</TABLE>")
    cat("<h3>OOB predictions (and true ---observed--- class): </h3>")
    html.data.frame(x[[5]], file = file, append = TRUE, first.col = "Subject/array")
    wrong.pred <- which(x[[5]][, 1] != x[[5]][, 2])
    if (length(wrong.pred)) {
        cat("<h4>Cases with incorrect predictions (errors): </h4>")
        html.data.frame(x[[5]][wrong.pred, ], file = file, append = TRUE, 
            first.col = "Subject/array")
    }
    ks <- length(x$genesSelected.cv)
    num.selected <- unlist(lapply(x$genesSelected.cv, length))
    cv.names <- paste("CV.run.", 1:ks, sep = "")
    cat("<br><br><h3> Stability assessments </h3>")
    cat("<h4> Genes selected in each of the cross-validation runs </h4>")
    for (i in 1:ks) {
        cat(paste("<h5>CV run  ", i, " (", num.selected[i], " genes selected):   ", 
            sep = ""), "</h5>")
        tmpgs <- x$geneNames[x$genesSelected.cv[[i]]]
        tmpgsl <- linkGene(tmpgs)
        cat("<TABLE frame=\"box\">\n")
        cat("<tr><th width=200>Gene name</th></tr>\n")
        for (i in 1:length(tmpgsl)) {
            cat("<tr><td>", tmpgsl[i], "</td></tr>\n")
        }
        cat("</TABLE>")
    }
    shared.genes <- matrix(NA, nrow = (ks + 1), ncol = (ks + 
        1))
    tmp.genesSelected <- list()
    tmp.genesSelected[[1]] <- x[[4]]
    tmp.genesSelected <- c(tmp.genesSelected, x$genesSelected.cv)
    for (i in 1:(ks + 1)) {
        for (j in 1:(ks + 1)) {
            shared.genes[i, j] <- length(intersect(tmp.genesSelected[[i]], 
                tmp.genesSelected[[j]]))
        }
    }
    prop.shared <- round(shared.genes/c(x[[3]], num.selected), 
        3)
    num.selected <- c(x[[3]], num.selected)
    num.selectedS <- paste("(", num.selected, ")", sep = "")
    colnames(shared.genes) <- colnames(prop.shared) <- c("OriginalSample", 
        cv.names)
    rownames(shared.genes) <- rownames(prop.shared) <- paste(c("OriginalSample", 
        cv.names), num.selectedS)
    options(width = 200)
    cat("<br><br><h4>Number of shared genes</h4>")
    cat("<pre>")
    print(as.table(shared.genes))
    cat("</pre>")
    cat("<h4>Proportion of shared genes (relative to row total)</h4>")
    cat("<pre>")
    print(as.table(prop.shared))
    cat("</pre>")
    options(width = 80)
    unlisted.genes.selected <- unlist(x$genesSelected.cv)
    named.unlisted.selected <- x$geneNames[unlisted.genes.selected]
    cat("<br><br><h4>Gene freqs. in cross-validated runs of genes selected in model with all data</h4>")

    if(length(x[[4]] == 1) && (x[[4]] == "All the genes")) {
      cat("<p> See table below </p>")
    } else {
      in.all.data <- which(names(table(named.unlisted.selected, 
                                       dnn = NULL)) %in% x$geneNames[x[[4]]])
      sgad <- sort(table(named.unlisted.selected, dnn = NULL)[in.all.data], 
                   decreasing = TRUE)
      sgadl <- linkGene(names(sgad))
      cat("<TABLE frame=\"box\">\n")
      cat("<tr><th width=200>Gene name</th><th>Frequency</th></tr>\n")
      for (i in 1:length(sgadl)) {
        cat("<tr><td>", sgadl[i], "</td><td>", sgad[i], "</td></tr>\n")
      }
      cat("</TABLE>")
    } 

    cat("<h4>Gene frequencies in cross-validated runs</h4>")
    tmp.table <- sort(table(named.unlisted.selected, dnn = NULL), 
        decreasing = TRUE)
    tmp.table.l <- linkGene(names(tmp.table))
    cat("<TABLE frame=\"box\">\n")
    cat("<tr><th width=200>Gene name</th><th>Frequency</th></tr>\n")
    for (i in 1:length(tmp.table.l)) {
        cat("<tr><td>", tmp.table.l[i], "</td><td>", tmp.table[i], 
            "</td></tr>\n")
    }
    cat("</TABLE>")
    sink()
}












########################################################
########################################################

###################   PAM



pamr.listgenes.fixed <- function (fit, data, threshold, genenames = FALSE) 
{ ##changed to prevent errors when no genes
    x <- data$x[fit$gene.subset, fit$sample.subset]
    if (genenames) {
        gnames <- data$genenames[fit$gene.subset]
    }
    if (!genenames) {
        gnames <- NULL
    }
    geneid <- data$geneid[fit$gene.subset]
    if (!is.null(fit$y)) {
        nc <- length(fit$y)
    }
    if (is.null(fit$y)) {
        nc <- ncol(fit$proby)
    }
    clabs <- colnames(fit$centroids)
    aa <- pamr.predict(fit, x, threshold = threshold, type = "nonzero")
    cen <- pamr.predict(fit, x, threshold = threshold, type = "cen")
    d <- (cen - fit$centroid.overall)[aa, , drop = FALSE]/fit$sd[aa]
    oo <- order(-apply(abs(d), 1, max))
    d <- round(d, 4)
    g <- gnames[aa]
    g1 <- geneid[aa]
    if (is.null(gnames)) {
        gnhdr <- NULL
    }
    if (!is.null(gnames)) {
        gnhdr <- "name"
    }
##    options(width = 500)
    schdr <- paste(clabs, "score", sep = " ")
    ## add the "drop = FALSE"
    res <- cbind(as.character(g1), g, d)[oo, ,drop = FALSE]
    
    ## somehtings for debugging; not first line

   if(is.null(dim(res))) { 
        the.problem.thing <<- list(g1 = g1, g = g,
                                   d = d, fit = fit, data = data,
                                   threshold = threshold,
                                   genenames = genenames)
        stop("It happened")
    }
    
    if(is.null(res)) { 
        the.problem.thing <<- list(g1 = g1, g = g,
                                   d = d, fit = fit, data = data,
                                   threshold = threshold,
                                   genenames = genenames)
        stop("It happened")
    }
    
    if(length(res) < 1) res <- NULL
    
    return(res)
###     ##    browser()
###     if(dim(res)[1]) { ## I add this
###         dimnames(res) <- list(NULL, c("id", gnhdr, schdr))
###         print(res, quote = FALSE)
###     } else NULL
}




VarSelPAM <- function(x, y, training.rows, testing.rows,
                      select = "size") {
    ## recall we have genes as columns, subjects are rows.
    ## 
    train.pam.data <- list(x = t(x[training.rows,, drop = FALSE]),
##                           y = factor(y[training.rows]),
                           y = y[training.rows],
                           geneid = 1:ncol(x),
                           genenames = NULL)
    if(!is.null(testing.rows)) {
        test.pam.data <- list(x = t(x[testing.rows, , drop = FALSE]),
                          y = NULL,
                          geneid = 1:ncol(x),
                          genenames = NULL)
    }
    trained.pam <- pamr.train(train.pam.data)
    trained.pam.cv <- pamr.cv(trained.pam, data = train.pam.data)

    if(select == "loglik") 
### for automated selection: minimize error and break ties by loglik.
        chosen.threshold <-
            trained.pam.cv$threshold[order(trained.pam.cv$error,
                                           -trained.pam.cv$loglik)[1]]

    if(select == "size") ## minimize number 
        chosen.threshold <-
            trained.pam.cv$threshold[order(trained.pam.cv$error,
                                           -trained.pam.cv$threshold)[1]]
    selected.genes <-
        as.numeric(pamr.listgenes.fixed(trained.pam,
                                  data = train.pam.data,
                                  threshold = chosen.threshold,
                                  genenames = FALSE)[, 1])
    
    if(!is.null(testing.rows)) {
####        predictions.prob <-
####        pamr.predict(fit = trained.pam,
####                     newx = test.pam.data$x,
####                     threshold = chosen.threshold,
####                     type = "posterior")

        predictions.class <-
            pamr.predict(fit = trained.pam,
                         newx = test.pam.data$x,
                         threshold = chosen.threshold,
                         type = "class")
    } else {
        predictions.class <- NULL
    }
    return(list(selected.vars = selected.genes,
                predictions.class = predictions.class,
##                predictions.prob = predictions.prob,
                other =
                list(chosen.treshold = chosen.threshold,
                     order.threshold = order(trained.pam.cv$error,
                                       -trained.pam.cv$loglik)[1],
                     trained.pam.cv = trained.pam.cv)))
}


f.pred.errors.PAM <- function(data, y, 
                              cv.error.k = 5,
                              title.figure = NULL,
                              plot = TRUE,
                              cv.error = TRUE) {

    ## If cv.error = FALSE,
    ## run the predictor "predictor.funct" using
    ## num.genes as the number of genes in the predictor
    ## (num.genes will often be a vector).
    ## Repeat it, for each value of num.genes, num.ks times,
    ## and obtain the cv error rate (with a knumber-fold cross validation).
    ## Plot the values of cv error rates.
    ##
    ## If cv.error = TRUE, in addition return the unbiased, honest
    ## estimate of error rate if we use the rule of choosing the minimim
    ## in the above graphic (and use cv.error.k for cv.error.k-fold measure).
    ##
    ## predictor.funct must return the prediction error of a given
    ## method or rule, with a given number of genes used as predictors.


##    if(!(is.factor(y))) stop("y should be a factor variable")

    gene.names <- colnames(data)
    N <- length(y)
    all.data.run <- VarSelPAM(data, y,
                              training.rows = 1:N,
                              testing.rows = NULL,
                              select = "size")
    the.best.genes <- all.data.run$selected.vars
    
    if(length(all.data.run$selected.vars) == ncol(data)) {
            the.best.genes <- "All the genes"
        }

    if(plot) {
##        par(las = 1)
##        par(mar = c(5, 5, 4, 4) + 0.1)
##        par(cex.lab = 1.5)
##        par(cex.axis = 1.5)
##        par(cex.main = 1.5)
        all.data.errors <- all.data.run$other$trained.pam.cv$error
        ngenes <- all.data.run$other$trained.pam.cv$size
        ##            num.points.plot <- length(ngenes)
###            maxplot <-
###                max(
###                    c(unlist(lapply(object$allBootRuns,
###                                    function(x) max(x$other$trained.pam.cv$error))),
###                      all.data.errors))
###            minplot <-
###                min(
###                    c(unlist(lapply(object$allBootRuns,
###                                    function(x) min(x$other$trained.pam.cv$error))),
###                      all.data.errors))

###            minplot <- minplot * (1 - 0.1)
###            maxplot <- maxplot * (1 + 0.2)
        minplot <- 0; maxplot <- 1
        plot(ngenes, all.data.errors, log = "x",
             type = "l", axes = TRUE, xlab = "Number of genes in predictor",
             ylab = "CV Error rate (but do not quote the smallest one!)",
             ylim = c(minplot, maxplot),
             lty = 1,
             col = "black", lwd = 2,
##             main = "CV Error rate vs. Number of genes in predictor.",
             main = title.figure,
             xlim = c(1, max(ngenes)*1.1))
###             plot(num.points.plot:1,
###                  object$allBootRuns[[1]]$other$trained.pam.cv$error,
###                  type = "l", axes = FALSE, xlab = "Number of genes",
###                  ylab = "CV Error rate", ylim = c(minplot, maxplot), lty = 2,
###                  main = "CV Error rate vs. Number of genes in predictor.")
###            legend(x = 10, y = maxplot,
###                   legend = c("Bootstrap samples", "Original sample"),
###                   lty = c(2, 1), lwd = c(1, 3), col = c("Black", "Red"))
            if(max(ngenes) > 300) axis(1, at = c(1, 2, 3, 5, 8, 15, 20, 25, 35, 50, 75, 150, 200, 300),
                                  labels = c(1, 2, 3, 5, 8, 15, 10, 25, 35, 50, 75, 150, 200, 300))
    }


    
    
    ### If we have asked to return the cross-validated
    ### error of the rule (this takes a while to get done)
   
    rule.cv.error <- NULL
    genesSelected.cv <- list()
    if(cv.error) { ### < if(cv.error)>
        rule.cv.errors <- rep(-99, cv.error.k)
        if(is.null(cv.error.k)) {
            cv.error.k <- length(y)
          }
        N <- length(y)
        OOB.predictions <- y

        index.select <- cv.balance(1:N, y, cv.error.k)
        
        for(sample.number in 1:cv.error.k) {
          print(paste(".... cv.error of rule; cv.sample number = ", sample.number))
###          train.x <- data[index.select != sample.number, , drop = FALSE]
###          train.y <- y[index.select != sample.number]
###          test.x <- data[index.select == sample.number, , drop = FALSE]
###          test.y <- y[index.select == sample.number]

          training.rows <- (1:N)[index.select != sample.number]
          testing.rows <- (1:N)[index.select == sample.number]
          pam.cv.run <- VarSelPAM(data, y,
                                  training.rows = training.rows,
                                  testing.rows = testing.rows,
                                  select = "size")

          if(plot) { ## we add the lines on the fly
              lines(pam.cv.run$other$trained.pam.cv$size,
                    pam.cv.run$other$trained.pam.cv$error,
                    lty = 1, col = "grey", lwd = 1.5)             
          }
            
          
####browser()
          OOB.predictions[testing.rows] <- pam.cv.run$predictions.class
          genesSelected.cv[[sample.number]] <- pam.cv.run$selected.vars
        }
        rule.cv.error <- sum(OOB.predictions != y)/N

        ## check levels are the same in input and output: I think it is OK,
        ## but this is a paranoid check
        ## <paranoid check>
        ## search below for paranoid. check: this is silly
#         oob.tmp <- factor(OOB.predictions)
#         try(paranoid.check <- any(OOB.predictions != oob.tmp))
#         if(class(paranoid.check) == "try-error") {
#           error.msg <- paranoid.check
#           stop(paste("OOOOOPS: try-error in paranoid check", error.msg))
#         } else if(paranoid.check) {
#           stop("OOOOPS: factor error in paranoid check")
#         }
        ## </paranoid check>
        
      } ### </ if(cv.error)>
#    browser()
    ## we also want the confussion matrix

    #### Hey, watch out: this only works if cv.error!!!!! xxxxxx fix it!!!!!!!
    if(cv.error) {
        conf.matrix <- table(y, OOB.predictions)
        overall.error <- (apply(conf.matrix, 1, sum) - diag(conf.matrix))/sum(conf.matrix)
        error.per.class <- (apply(conf.matrix, 1, sum) - diag(conf.matrix)) / apply(conf.matrix, 1, sum)
        conf.matrix <- cbind(conf.matrix, "TotalError" = overall.error,
                             "RelativeErrorPerClass" = error.per.class)
        rownames(conf.matrix) <- paste("Observed", rownames(conf.matrix))
        OOB.predictions <- data.frame(OOB.predictions = OOB.predictions,
                                      TrueClass = y)
        rownames(OOB.predictions) <- rownames(data)
    
        ##<paranoid test 2>
        if(abs(rule.cv.error - sum(overall.error)) > 1e-06)
            stop("OOOOPS: failed second paranoid test")
        ##</paranoid test 2>
    } else {
        OOB.predictions <- NULL
        conf.matrix <- NULL
    }
    
    if(plot) {
        nk <- as.vector(table(y))
        error.max.bet <- 1 - (max(nk)/sum(nk))
        abline(h = error.max.bet, lty = 2, col = "blue", lwd = 1.5)
        axis(4, at = error.max.bet, labels = round(error.max.bet, 3),
             col = "blue", col.axis = "blue", las = 1)        
      if(cv.error) {
        abline(h = rule.cv.error, lty = 2, col = "red", lwd = 1.5)
        axis(4, at = rule.cv.error, labels = round(rule.cv.error, 3),
             col = "red", col.axis = "red", las = 1)
        legend(x = 2, y = 0.95,
               c("Error rate w.o. predictor",
                 "Error rate of predictor building rule",
                 "Original sample",
                 "CV samples"),
               lty = c(2, 3, 1, 1),
               col = c("blue", "red", "black", "grey"),
               cex = 1, lwd = 2)


        
###        legend(x = 2, y = 0.95,
###               c("Error rate w.o. predictor", "Error rate of predictor building rule"),
###                 lty = 2, col = c("blue", "red"), cex = 2, lwd = 2)
      } else {
        legend(x = 2, y = 0.95,
               c("Error rate w.o. predictor", "Original sample"),
               lty = c(2, 1), col = c("blue", "black"),
               cex = 1, lwd = 2)
      }  
    }

    if(plot) {
        ## Paint this line again, just in case.
        lines(ngenes, all.data.errors, col = "black", lwd = 2)
    }

   
    ret.val <- list("Error Rate" = rule.cv.error,
                "ConfussionMatrix" = conf.matrix,
                "Number of predictors that yields minimum error" =
                    length(all.data.run$selected.vars),
                "Selected predictor genes" = the.best.genes,
                "OOB predictions" = OOB.predictions,
                    "genesSelected.cv" = genesSelected.cv,
                    "geneNames" = colnames(data),
                    "subjectNames" = rownames(data))
    class(ret.val) <- "cvPred"
    ret.val
}






########   Other methods


cross.valid.predictor <- function(x, y,
                                  training.rows = NULL,
                                  testing.rows = NULL,
                                  typePredictor = "dldaC",
                                  geneSelection = "F", 
##                                  test.data = NULL, test.class = NULL,
                                  sizeGenes = c(5, 10, 50), knumber = 5) {

    if(!is.null(training.rows)) {
        x.orig <- x
        y.orig <- y
        x <- x[training.rows, , drop = FALSE]
        y <- y[training.rows]
    }
    if(!is.null(testing.rows) & is.null(training.rows))
        stop("If you specify testing rows you MUST specify training rows")
    numgenes <- ncol(x)
    N <- length(y)
    
    fPredictor <-
        eval(parse(text = paste("internalPredictor",
                   typePredictor, sep = "."))) 

    gene.select <-
        eval(parse(text = paste("geneSelect",
                   geneSelection, sep = ".")))

    ## First, selecting the data subsets for cross-validation
    ## data is structured as subjects in rows, genes in columns
    ## (and thus is transposed inside gene.select to be fed to mt.maxT).

###    if(is.null(knumber)) {
###        knumber <- length(y)
###    }
    
    index.select <- cv.balance(1:N, y, knumber) ## cv with balance in
    ## class props. in testing and testing
    
    cv.errors <- matrix(-99, nrow = length(sizeGenes),
                        ncol = knumber)
    
    genesSelected <- list()
    
    
    for(sample.number in 1:knumber) {
        intra.training.rows <- (1:N)[index.select != sample.number]
        intra.testing.rows <- (1:N)[index.select == sample.number]
        gene.subset <- gene.select(x[intra.training.rows, ],
                                   y[intra.training.rows])
        kiter <- 1
        for(k in sizeGenes) {
            intra.test.data <- x[intra.testing.rows, gene.subset[1:k],
                                 drop = FALSE]
            intra.train.data <- x[intra.training.rows, gene.subset[1:k],
                                  drop = FALSE]
            predicted <- fPredictor(intra.train.data,
                                    y[intra.training.rows],
                                    intra.test.data)
            cv.errors[kiter, sample.number] <-
                sum(y[intra.testing.rows] != predicted)
            kiter <- kiter + 1
        }
    }
    if(any(cv.errors < 0)) stop("Something bad happened: some cv.errors < 0")
    cv.error <-  apply(cv.errors, 1, sum)
    cv.error <- cv.error/N
    genes.minimum.error <- sizeGenes[which.min(cv.error)]
    
    if(genes.minimum.error == numgenes) {
        selected.genes <- 1:numgenes
    } else {
        selected.genes <- gene.select(x, y)[1:genes.minimum.error]
    }
    
    if(!is.null(testing.rows)) { ## we are passed a test data and test class
        test.data <- x.orig[testing.rows, selected.genes, drop = FALSE]
        train.data <- x[ ,selected.genes, drop = FALSE]
        predictions.class <- fPredictor(train.data,
                                        y,
                                        test.data)
    } else {
        predictions.class <- NULL
    }
    return(list(selected.vars = selected.genes,
                predictions.class = predictions.class,
                other =
                list(sizeGenes = sizeGenes,
                     cv.error = cv.error)
                ))
}


f.pred.errors <- function(data, y, 
                          sizeGenes,
                          typePredictor = "dldaC",
                          geneSelection = "F",
                          knumber = 5,
                          title.figure = NULL,
                          plot = TRUE,
                          cv.error = TRUE) {

    ## If cv.error = FALSE,
    ## run the predictor "predictor.funct" using
    ## num.genes as the number of genes in the predictor
    ## (num.genes will often be a vector).
    ## Repeat it, for each value of num.genes, num.ks times,
    ## and obtain the cv error rate (with a knumber-fold cross validation).
    ## Plot the values of cv error rates.
    ##
    ## If cv.error = TRUE, in addition return the unbiased, honest
    ## estimate of error rate if we use the rule of choosing the minimim
    ## in the above graphic (and use knumber for knumber-fold measure).
    ##
    ## predictor.funct must return the prediction error of a given
    ## method or rule, with a given number of genes used as predictors.


##    if(!(is.factor(y))) stop("y should be a factor variable")

    gene.names <- colnames(data)
    N <- length(y)
    all.data.run <- cross.valid.predictor(data, y,
                                          training.rows = 1:N,
                                          testing.rows = NULL,
                                          typePredictor = typePredictor,
                                          geneSelection = geneSelection,
                                          sizeGenes = sizeGenes,
                                          knumber = knumber)
    the.best.genes <- all.data.run$selected.vars
    
    if(length(all.data.run$selected.vars) == ncol(data)) {
            the.best.genes <- "All the genes"
        }

    if(plot) {
##        par(las = 1)
##        par(mar = c(5, 5, 4, 4) + 0.1)
##        par(cex.lab = 1.5)
##        par(cex.axis = 1.5)
##        par(cex.main = 1.5)
        
        all.data.errors <- all.data.run$other$cv.error
        ngenes <- all.data.run$other$sizeGenes
        minplot <- 0; maxplot <- 1
        plot(ngenes, all.data.errors, log = "x",
             type = "l", xlab = "Number of genes in predictor",
             ylab = "CV Error rate (but do not quote the smallest one!)",
             ylim = c(minplot, maxplot),
             lty = 1,
             col = "black", lwd = 2,
             main = title.figure,
             ##             main = "CV Error rate vs. Number of genes in predictor.",
             xlim = c(1, max(ngenes)*1.1),
             axes = FALSE)
        box()
        axis(2, las = 1)
##        if(max(ngenes) > 300)
        axis(1, at = ngenes,
             labels = ngenes)
    }
    
    
    
### If we have asked to return the cross-validated
    ### error of the rule (this takes a while to get done)
    
    rule.cv.error <- NULL
    genesSelected.cv <- list()
    if(cv.error) { ### < if(cv.error)>
        rule.cv.errors <- rep(-99, knumber)
        if(is.null(knumber)) {
            knumber <- length(y)
        }
        N <- length(y)
        OOB.predictions <- y
        
        index.select <- cv.balance(1:N, y, knumber)
        
        for(sample.number in 1:knumber) {
            print(paste(".... cv.error of rule; cv.sample number = ", sample.number))
            training.rows <- (1:N)[index.select != sample.number]
            testing.rows <- (1:N)[index.select == sample.number]
            cv.run <- cross.valid.predictor(data, y,
                                            training.rows = training.rows,
                                            testing.rows = testing.rows,
                                            typePredictor = typePredictor,
                                            geneSelection = geneSelection,
                                            sizeGenes = sizeGenes,
                                            knumber = knumber)

          if(plot) { ## we add the lines on the fly
              lines(cv.run$other$sizeGenes,
                    cv.run$other$cv.error,
                   lty = 1, col = "grey", lwd = 1.5)             
          }
            OOB.predictions[testing.rows] <- cv.run$predictions.class
            genesSelected.cv[[sample.number]] <- cv.run$selected.vars
        }
        rule.cv.error <- sum(OOB.predictions != y)/N

        ## check levels are the same in input and output: I think it is OK,
        ## but this is a paranoid check
        ## <paranoid check>
        ## NO!! This is silly: it will fail if OOB.predictions are all identical.
        ## Commented out!!!
#         
#         oob.tmp <- factor(OOB.predictions)
#         try(paranoid.check <- any(OOB.predictions != oob.tmp))
#         if(class(paranoid.check) == "try-error") {
#           error.msg <- paranoid.check
#           stop(paste("OOOOOPS: try-error in paranoid check", error.msg))
#         } else if(paranoid.check) {
#           stop("OOOOPS: factor error in paranoid check")
#         }
        ## </paranoid check>
        
      } ### </ if(cv.error)>
#    browser()
    ## we also want the confussion matrix

    #### Hey, watch out: this only works if cv.error!!!!! xxxxxx fix it!!!!!!!
    if(cv.error) {
        conf.matrix <- table(y, OOB.predictions)
        overall.error <- (apply(conf.matrix, 1, sum) - diag(conf.matrix))/sum(conf.matrix)
        error.per.class <- (apply(conf.matrix, 1, sum) - diag(conf.matrix)) / apply(conf.matrix, 1, sum)
        conf.matrix <- cbind(conf.matrix, "TotalError" = overall.error,
                             "RelativeErrorPerClass" = error.per.class)
        rownames(conf.matrix) <- paste("Observed", rownames(conf.matrix))
        OOB.predictions <- data.frame(OOB.predictions = OOB.predictions,
                                      TrueClass = y)
        rownames(OOB.predictions) <- rownames(data)
    
        ##<paranoid test 2>
        if(abs(rule.cv.error - sum(overall.error)) > 1e-06)
            stop("OOOOPS: failed second paranoid test")
        ##</paranoid test 2>
    } else {
        OOB.predictions <- NULL
        conf.matrix <- NULL
    }
    
    if(plot) {
        nk <- as.vector(table(y))
        error.max.bet <- 1 - (max(nk)/sum(nk))
        abline(h = error.max.bet, lty = 2, col = "blue", lwd = 1.5)
        axis(4, at = error.max.bet, labels = round(error.max.bet, 3),
             col = "blue", col.axis = "blue", las = 1)
        
      if(cv.error) {
        abline(h = rule.cv.error, lty = 2, col = "red", lwd = 1.5)
        axis(4, at = rule.cv.error, labels = round(rule.cv.error, 3),
             col = "red", col.axis = "red", las = 1)
        legend(x = 2, y = 0.95,
               c("Error rate w.o. predictor",
                 "Error rate of predictor building rule",
                 "Original sample",
                 "CV samples"),
               lty = c(2, 3, 1, 1),
               col = c("blue", "red", "black", "grey"),
               cex = 1, lwd = 2)

###        legend(x = 2, y = 0.95,
###               c("Error rate w.o. predictor", "Error rate of predictor building rule"),
###                 lty = 2, col = c("blue", "red"), cex = 2, lwd = 2)
      } else {
        legend(x = 2, y = 0.95,
               c("Error rate w.o. predictor", "Original sample"),
               lty = c(2, 1), col = c("blue", "black"),
               cex = 1, lwd = 2)
      }  
    }

    if(plot) {
        ## Paint this line again, just in case.
        lines(ngenes, all.data.errors, col = "black", lwd = 2)
    }

   
    ret.val <- list("Error Rate" = rule.cv.error,
                    "ConfussionMatrix" = conf.matrix,
                    "Number of predictors that yields minimum error" =
                    length(all.data.run$selected.vars),
                    "Selected predictor genes" = the.best.genes,
                    "OOB predictions" = OOB.predictions,
                    "genesSelected.cv" = genesSelected.cv,
                    "geneNames" = colnames(data),
                    "subjectNames" = rownames(data))
    class(ret.val) <- "cvPred"
    ret.val
}









###### Exporting data for testing and examples:

####colnames(gl.data) <- paste("G", 1:dim(gl.data)[2], sep = "")
####gl.data <- t(gl.data)

####write.table(gl.data, file = "gl.data.txt", row.names = TRUE,
####            col.names = NA, sep = "\t", quote = FALSE)

####colnames(vv.data) <- paste("G", 1:dim(vv.data)[2], sep = "")
####write.table(vv.data, file = "vv.data.txt", row.names = TRUE,
####            col.names = NA, sep = "\t", quote = FALSE)

####colnames(lymphoma.data) <- paste("G", 1:dim(lymphoma.data)[2], sep = "")
####write.table(lymphoma.data, file = "lymphoma.data.txt", row.names = TRUE,
####            col.names = NA, sep = "\t", quote = FALSE)

####colnames(colon.data) <- paste("G", 1:dim(colon.data)[2], sep = "")

####write.table(colon.data, file = "colon.data.txt", row.names = TRUE,
####            col.names = NA, sep = "\t", quote = FALSE)

####gl.cl2 <- factor(gl.class, labels = c("ALL", "AML"))
####write.table(gl.cl2, file = "gl.class.txt", row.names = FALSE,
####            col.names = FALSE, sep = "\t", eol = "\t",
####            quote = FALSE)

####vv.cl2 <- factor(vv.class, labels = c("Metastasis", "DiseaseFree"))
####lymphoma.cl2 <- factor(lymphoma.class,
####                       labels = c("DLBCL", "FL", "CLL"))
####colon.cl2 <- factor(colon.class,
####                    labels = c("Normal", "Tumor"))


####write.table(vv.cl2, file = "vv.class.txt", row.names = FALSE,
####            col.names = FALSE, sep = "\t", eol = "\t",
####            quote = FALSE)
####write.table(colon.cl2, file = "colon.class.txt", row.names = FALSE,
####            col.names = FALSE, sep = "\t", eol = "\t",
####            quote = FALSE)
####write.table(lymphoma.cl2, file = "lymphoma.class.txt", row.names = FALSE,
####            col.names = FALSE, sep = "\t", eol = "\t",
####            quote = FALSE)



######### Examples

####X <- matrix(rnorm(5000 * 50), nrow = 50)
####colnames(X) <- paste("v", 1:ncol(X), sep = "@")
####rownames(X) <- paste("r", 1:nrow(X), sep = ".")
####clases <- c(rep(0, 25), rep(1, 25))
####fac1 <- factor(clases)


####num.genes <- c(5, 10)

####cross.valid.predictor(gl.data[, 1:100], gl.class, sizeGenes = c(5, 10))



####load("all.real.data.RData")


####lcl <- factor(lymphoma.class, labels = c("b", "m", "v"))


####colnames(lymphoma.data) <- paste("v", 1:ncol(lymphoma.data), sep = "@")
####rownames(lymphoma.data) <- paste("r", 1:nrow(lymphoma.data), sep = ".")


####r3 <- f.pred.errors.PAM(lymphoma.data[, 1:100], lcl,
####                        cv.error.k = 2,
####                        cv.error = TRUE,
####                        plot = TRUE, ## all other arguments are irrelevant
####                        typePred = "knn", num.genes,
####                        num.ks = 2, knumber = 2,
####                        title.figure = "DLDA",
####                        best.k = 1)

####r4 <- f.pred.errors.PAM(lymphoma.data, lcl,
####                        cv.error.k = 2,
####                        cv.error = TRUE,
####                        plot = TRUE, ## all other arguments are irrelevant
####                        typePred = "knn", num.genes,
####                        num.ks = 2, knumber = 2,
####                        title.figure = "DLDA",
####                        best.k = 1)
####r5 <- f.pred.errors.PAM(lymphoma.data, lcl,
####                        cv.error.k = 5,
####                        cv.error = TRUE,
####                        plot = TRUE, ## all other arguments are irrelevant
####                        typePred = "knn", num.genes,
####                        num.ks = 2, knumber = 2,
####                        title.figure = "DLDA",
####                        best.k = 1)





####r10b <- f.pred.errors(lymphoma.data, lcl, 
####                      sizeGenes = c(5, 10, 20),
####                      typePredictor = "knn",
####                      geneSelection = "F",
####                      knumber = 5,
####                      title.figure = "DLDA")

####r11b <- f.pred.errors(lymphoma.data, lcl, 
####                      sizeGenes = num.genes,
####                      typePredictor = "knn",
####                      geneSelection = "F",
####                      knumber = 5,
####                      title.figure = "DLDA")

####r11b <- f.pred.errors(lymphoma.data, lcl, 
####                      sizeGenes = num.genes,
####                      typePredictor = "knn",
####                      geneSelection = "F",
####                      knumber = 5,
####                      title.figure = "DLDA")

####rownames(colon.data) <- paste("r", 1:dim(colon.data)[1], sep = "")
####r12b <- f.pred.errors(colon.data, colon.class, 
####                      sizeGenes = c(num.genes, 2000),
####                      typePredictor = "knn",
####                      geneSelection = "F",
####                      knumber = 5,
####                      title.figure = "DLDA")

####r12b <- f.pred.errors(colon.data, colon.class, 
####                      sizeGenes = c(num.genes, 2000),
####                      typePredictor = "knn",
####                      geneSelection = "Wilcoxon",
####                      knumber = 5,
####                      title.figure = "DLDA")

####r13 <- f.pred.errors(colon.data, colon.class, 
####                      sizeGenes = c(num.genes, 2000),
####                      typePredictor = "knn",
####                      geneSelection = "rF",
####                      knumber = 5,
####                      title.figure = "DLDA")

####r14 <- f.pred.errors(colon.data, colon.class, 
####                      sizeGenes = c(num.genes, 2000),
####                      typePredictor = "rF",
####                      geneSelection = "F",
####                      knumber = 5,
####                      title.figure = "DLDA")


