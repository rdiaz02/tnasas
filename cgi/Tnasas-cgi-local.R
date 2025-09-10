## like the cgi, but without it, so can ran locally.
##  create a /tmp directory and ./tmp/arrays

#! /usr/bin/R


source("Tnasas.R")
runif(1)
save.seed <- .Random.seed ## in case there are problems
######   This might be specific to our application

datadir <- "./tmp"

subject.names <- scan(paste(datadir, "arrays", "array_names", sep = "/"),
                      what = "char")
gene.names <- scan(paste(datadir, "arrays", "gene_names", sep = "/"),
                   what = "char")
classes <- factor(scan(paste(datadir, "arrays", "classes", sep = "/"),
                       what = "char"))

datos <- read.table(paste(datadir, "arrays", "predictors", sep = "/"),
                    header = FALSE, sep = "\t", stringsAsFactors = TRUE)
datos <- t(datos)
rownames(datos) <- subject.names
colnames(datos) <- gene.names
num.cols.datos <- ncol(datos)


### for now, this is hard-wired:

num.genes <- c(5, 10, 50, 200, 500)
num.genes <- c(5, 200) ### so that it is faster
if(length(which(num.genes > num.cols.datos)))
    num.genes <- num.genes[-which(num.genes > num.cols.datos)]
cv.error.k <- knumber <- 10
num.ks <- 5
best.k <- 1



######  The actual calls for the CGI


##the.predictor <- switch(formData$model,
##                        svm = "svmLinear",
##                        dlda = "dlda",
##                        knn = "knn",
##                        randomforest = "randomForest"
##                        )
##if(is.null(the.predictor)) stop("the.predictor is NULL")

##title.figure <- switch(formData$model,
##                        svm = "SVM (linear kernel)",
##                        dlda = "DLDA",
##                        knn = "NN",
##                        randomforest = "Random forest"
##                        )


the.predictor <- "dlda"
title.figure <- "figura"



graphDir <- paste(datadir, "/", sep = "")
#url <- "/temp/"
#graphURLroot <- paste(url, formData$pid, "/", sep="")

## Define a better name for the image (avoiding confussions with several runs)
##image.file <- paste("predictor_error_rates_", formData$pid, ".png", sep = "")

## Define the file results (including path cause we're not running
## the script in the temp_process dir)

results.file <- paste(datadir, "/", "results_", ".txt", sep = "")
selected.genes.file <- paste(datadir, "selectedgenes.txt", sep = "/")
nonselected.genes.file <- paste(datadir, "nonselectedgenes.txt", sep = "/")




## utility function
n.print <- function(x) paste(substitute(x), paste(x, collapse = ", "),
                                 sep = " : ")


## minimal error checking
if(min(table(classes)) <  5) {
    postscript()
    plot(c(0,1), c(0,1), type = "n", axes = FALSE, xlab = "", ylab = "")
    box()
    text(0.5, 0.5, "The smallest class must have at least 5 cases", cex = 1.5)
    text(0.5, 0.35, "Yours doesn't", cex = 1.5)
    text(0.5, 0.2, "Please read the help and try again")
    dev.off()
    sink(file = results.file)
    cat("The smallest class must have at least 5 cases\n")
    cat("Yours doesn't\n")
    cat("Please read the help and try again")
    sink()
} else if(length(classes) != nrow(datos)) {
    postscript()
    plot(c(0,1), c(0,1), type = "n", axes = FALSE, xlab = "", ylab = "")
    box()
    text(0.5, 0.5, "The number of cases, as given in the classes")
    text(0.5, 0.35, "of course must equal the number of cases")
    text(0.5, 0.2, "in the data file. Yours doesn't. Please fix your error")
    dev.off()
    sink(file = results.file, append = TRUE)
    cat("The number of cases, as given in the classes")
    cat("of course must equal the number of cases")
    cat("in the data file. Yours doesn't. Please fix your error")
    sink()
} else if((ncol(datos) < 5) | (nrow(datos) == 0)) {
    postscript()
    plot(c(0,1), c(0,1), type = "n", axes = FALSE, xlab = "", ylab = "")
    box()
    text(0.5, 0.5, "You surely should have at least one case and")
    text(0.5, 0.35, "at least five genes.")
    text(0.5, 0.2, "Please, check your data file and try again.")
    dev.off()
    sink(file = results.file, append = TRUE)
    cat("You surely should have at least one case and")
    cat("at least five genes.")
    cat("Please, check your data file and try again.")
    sink()
} else { ## <OK, assuming no other probable mistakes...>

    postscript()
    the.run <- try(f.pred.errors(datos, classes, typePred = the.predictor,
                                 num.genes = num.genes,
                                 num.ks = num.ks,
                                 knumber = knumber,
                                 cv.error.k = cv.error.k,
                                 title.figure = title.figure,
                                 best.k = best.k))
    dev.off()

    write(c("OPTIONS (some might not mean much to you, but they can be",
            "important for tracking errors):",
            n.print(num.genes),
            n.print(cv.error.k),
            n.print(knumber),
            n.print(num.ks),
            n.print(best.k),
            n.print(the.predictor),
            "**********************************************************",
            "**********************************************************",
            " RESULTS: "
            ),
          file = results.file)

    if(class(the.run) == "try-error") {
        sink(file = results.file, append = TRUE)
        cat(" ..... starting traceback \n")
        traceback()
        cat("\n \n The error \n")
        print(the.run)
        cat("\n \n The random number generator seed \n")
        save.seed
        sink()
        save(save.seed, file = paste(datadir, "seed.RData", sep = "/"))
        postscript()
        plot(c(0,1), c(0,1), type = "n", axes = FALSE, xlab = "", ylab = "")
        box()
        text(0.5, 0.5, "An error occurred", cex = 2)
        text(0.5, 0.2, "Please check results file and let us know")
        dev.off()
    } else {##<no errors>
        sink(file = results.file, append = TRUE)
        print(the.run)
        sink()

        write.table(the.run[[4]], sep = "\t", quote = FALSE,
                    row.names = FALSE, col.names = FALSE,
                    file = selected.genes.file)

        write.table(setdiff(colnames(datos), the.run[[4]]), sep = "\t",
                    quote = FALSE,
                    row.names = FALSE, col.names = FALSE,
                    file = nonselected.genes.file)
    }##</no errors>
} ## </OK, assuming no other probable mistakes...>
cat("OK")








#### Algunas pruebecitas:

subject.names <- NULL
gene.names <- NULL

classes <- factor(scan(paste(datadir, "arrays", "golub.class.txt", sep = "/"),
                       what = "char"))
datos <- read.table(paste(datadir, "arrays", "golub.data.txt", sep = "/"),
                    header = FALSE, sep = "\t", stringsAsFactors = TRUE)
datos <- t(datos)
num.cols.datos <- ncol(datos)



postscript()
the.run <- try(f.pred.errors(datos, classes, typePred = "dlda",
                             num.genes = num.genes,
                             num.ks = num.ks,
                             knumber = knumber,
                             cv.error.k = cv.error.k,
                             title.figure = title.figure,
                             best.k = best.k))
dev.off()


postscript()
the.run <- try(f.pred.errors(datos, classes, typePred = "svmLinear",
                             num.genes = num.genes,
                             num.ks = num.ks,
                             knumber = knumber,
                             cv.error.k = cv.error.k,
                             title.figure = title.figure,
                             best.k = best.k))
dev.off()





classes <- factor(scan(paste(datadir, "arrays", "nci.class.txt", sep = "/"),
                       what = "char"))
levels(classes) <- letters[8:1]
datos <- read.table(paste(datadir, "arrays", "nci.data.txt", sep = "/"),
                    header = FALSE, sep = "\t", stringsAsFactors = TRUE)
datos <- t(datos)
num.cols.datos <- ncol(datos)

postscript()
the.run <- try(f.pred.errors(datos, classes, typePred = "dlda",
                             num.genes = num.genes,
                             num.ks = num.ks,
                             knumber = knumber,
                             cv.error.k = cv.error.k,
                             title.figure = title.figure,
                             best.k = best.k))
dev.off()

postscript()
the.run <- try(f.pred.errors(datos, classes, typePred = "svmLinear",
                             num.genes = num.genes,
                             num.ks = num.ks,
                             knumber = knumber,
                             cv.error.k = cv.error.k,
                             title.figure = title.figure,
                             best.k = best.k))
dev.off()



classes <- factor(scan(paste(datadir, "arrays", "vVeer3.class.txt", sep = "/"),
                       what = "char"))
datos <- read.table(paste(datadir, "arrays", "vVeer3.data.txt", sep = "/"),
                    header = FALSE, sep = "\t", stringsAsFactors = TRUE)
datos <- t(datos)
num.cols.datos <- ncol(datos)



postscript()
the.run <- try(f.pred.errors(datos, classes, typePred = "dlda",
                             num.genes = num.genes,
                             num.ks = num.ks,
                             knumber = knumber,
                             cv.error.k = cv.error.k,
                             title.figure = title.figure,
                             best.k = best.k))
dev.off()


postscript()
the.run <- try(f.pred.errors(datos, classes, typePred = "svmLinear",
                             num.genes = num.genes,
                             num.ks = num.ks,
                             knumber = knumber,
                             cv.error.k = cv.error.k,
                             title.figure = title.figure,
                             best.k = best.k))
dev.off()




### for quick tests
num.genes <- c(5, 10) ### so that it is faster
cv.error.k <- knumber <- 5
num.ks <- 2
best.k <- 1

the.run <- try(f.pred.errors(datos, classes, typePred = "dlda",
                             num.genes = num.genes,
                             num.ks = num.ks,
                             knumber = knumber,
                             cv.error.k = cv.error.k,
                             title.figure = title.figure,
                             best.k = best.k))
