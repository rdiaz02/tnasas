####  Copyright (C) 2004-2006  Ram�n D�az-Uriarte and JuanM Vaquerizas

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



rm(list = ls())

.Last <- function() {
    RterminatedOK <- file("RterminatedOK", "w")
    cat("\nNormal termination\n", file = RterminatedOK)
    flush(RterminatedOK)
    close(RterminatedOK)
    ##  try(stopCluster(TheCluster))
  cat("\n\n Normal termination\n")
##  try(system(paste("/http/mpi.log/killLAM.py", lamSESSION, "&")))
##  try(mpi.quit(save = "no"), silent = TRUE)
}

startExecTime <- format(Sys.time())

pid <- Sys.getpid()
write.table(file = "pid.txt", pid,
            row.names = FALSE,
            col.names = FALSE)

##########################################################
######
######   This might be specific to our setup at CNIO
######
##########################################################
system(paste("/http/mpi.log/counterAppsFromR.py Tnasas", getwd()))
## attach pid to name in R.running.procs
hostn <- system("hostname", intern = TRUE)
cat("\n HOSTNAME IS ", hostn, "\n")
new.name1 <- unlist(strsplit(getwd(), "/"))
new.name1 <- paste(new.name1[length(new.name1)], "@", hostn, sep = "")
new.name <- paste("R.", new.name1, "%", pid, sep = "")
new.name1 <- paste("R.", new.name1, sep = "")
system(paste("mv ../../R.running.procs/", new.name1,
             " ../../R.running.procs/", new.name,
             sep = ""))


library(Tnasas)
runif(1)
library(Cairo)
## library(GDD)
                                        #library(CGIwithR)
save.seed <- .Random.seed ## in case there are problems
png.width <- 500
png.height <- 500
##png.res = 144
png.pointsize <- 12
png.family = "Helvetica"
graphDir <- paste(getwd(), "/", sep = "")
##datadir <- paste("../www/tmp", formData$pid, sep = "/") ## zz: delete?


##############################################
##############################################
######                              ##########
######         Error checking       ##########
######                              ##########
##############################################
##############################################


caughtUserError <- function(message) {
    CairoPNG("predictor_error_rates.png", width = png.width,
           height = png.height, 
           ps = png.pointsize)
    plot(x = c(0, 1), y = c(0, 1),
         type = "n", axes = FALSE, xlab = "", ylab = "")
    box()
    text(0.5, 0.7, "There was a PROBLEM with your data.")
    text(0.5, 0.5,
    "Please read carefully the error messages under Results,")
    
    text(0.5, 0.3, "fix the problem, and try again.")
    dev.off()
    sink(file = "results.txt")
    cat(message)
    sink()
    sink(file = "exitStatus")
    cat("Error\n\n")
    cat(message)
    sink()
    quit(save = "no", status = 11, runLast = TRUE)
}

caughtOurError <- function(message) {
    CairoPNG("predictor_error_rates.png", width = png.width,
           height = png.height, 
           ps = png.pointsize)
    plot(x = c(0, 1), y = c(0, 1),
         type = "n", axes = FALSE, xlab = "", ylab = "")
    box()
    text(0.5, 0.7, "There was a PROBLEM with the code.")
    text(0.5, 0.5,
    "Please let us know (send us the URL),")
    
    text(0.5, 0.3, "so that we can fix it.")
    dev.off()
    sink(file = "results.txt")
    cat(message)
    sink()
    sink(file = "exitStatus")
    cat("Error\n\n")
    cat(message)
    sink()
    quit(save = "no", status = 11, runLast = TRUE)
}



### FIXME: move to the Tnasas package later
writeForPaLS <- function(x, file1 = "Selected.genes.txt",
                         file2 = "Selected.and.CV.selected.txt") {

    l1 <- c("#SelectedGenesMainRun", x$geneNames[x[[4]]])
    write(l1, file = file1)

    l2 <- l1
    ks <- length(x$genesSelected.cv)
    for(i in 1:ks) {
        l2 <- c(l2, paste("#CV.run.", i, sep = ""),
                x$geneNames[x$genesSelected.cv[[i]] ])
    }
    write(l2, file = file2)
}




###   CGI must write to tmp directory:
## model, genesel, covariate, arrayNames, Class
## optionally:
##     idtype, organism


### start input

idtype <- try(scan("idtype", what = "", n = 1))
organism <- try(scan("organism", what = "", n = 1))

classModel <- try(scan("model", what = "", n = 1))
geneSel <- try(scan("genesel", what = "", n = 1))

## can make it faster via scan: zz
#### datos <- read.table(paste(datadir, "arrays", "predictors", sep = "/"),
####                     header = FALSE, sep = "\t")
#### subject.names <- scan(paste(datadir, "arrays", "array_names", sep = "/"), sep = "\t",
####                       what = "char", quote = NULL)
#### gene.names <- scan(paste(datadir, "arrays", "gene_names", sep = "/"), sep = "\t",
####                    what = "char", quote = NULL)
#### classes <- factor(scan(paste(datadir, "arrays", "classes", sep = "/"), sep = "\t",
####                        what = "char", quote = NULL))
#### ## to prevent problems with a space at end of classes:
#### if(classes[length(classes)] == "") {
####   classes <- classes[-length(classes)]
####   classes <- factor(classes)
#### }
#### datos <- t(datos)
#### rownames(datos) <- subject.names
#### colnames(datos) <- gene.names
#### num.cols.datos <- ncol(datos)



## For now, I expect covariate file with genes in first column and a separate file
## with array names. And of course, the class file.

tmp <- try(
       	   system("sed 's/\"//g' covariate > cvt; mv cvt covariate")
	   )


num.cols.covariate <- count.fields("covariate", sep = "\t",
                                   quote = "",
                                   comment.char = "#")
if(length(unique(num.cols.covariate)) > 1) {
    emessage <-
    paste("The number of columns in your covariate file\n",
          "is not the same for all rows (genes).\n",
          "We find the following number of columns\n",
          paste(num.cols.covariate, collapse = " "))
    caughtUserError(emessage)
}
tryxdata <- try(
xdata <- read.table("covariate", header = FALSE, sep = "\t",
                    strip.white = TRUE,
                    comment.char = "#",
		    quote = ""))
if(class(tryxdata) == "try-error")
    caughtUserError("The covariate file is not of the appropriate format\n")
geneNames <- xdata[, 1]
xdata <- xdata[, -1]
if(length(unique(geneNames)) < length(geneNames)) {
    dupnames <- which(duplicated(geneNames))
    emessage <- paste("Gene names are not unique.\n",
                      "Please change them so that they are unique.\n",
                      "The duplicated names are in rows", dupnames, "\n")
    caughtUserError(emessage)
}
rownames(xdata) <- geneNames
arrayNames <- scan("arrayNames", sep = "\t", what = "char", quote = "")
if(length(arrayNames) > 0) {
    arrayNames <- arrayNames[-1]
    if(length(unique(arrayNames)) < length(arrayNames)) {
        dupnames <- which(duplicated(arrayNames))
        emessage <- paste("Array names are not unique.\n",
                          "Please change them so that they are unique.\n",
                          "The duplicated names are ", dupnames, "\n")
        caughtUserError(emessage)
    }
    colnames(xdata) <- arrayNames
}
xdata <- t(xdata)
trycl <- try(
             Class <- factor(scan("class", sep = "\t", what = "char", strip.white = TRUE, nlines = 1))
             )
## to prevent problems with a space at end of classes
if(class(trycl) == "try-error")
    caughtUserError("The class file is not of the appropriate format\n")
if(Class[length(Class)] == "") Class <- factor(Class[-length(Class)])
tclass <- table(Class)
if(length(unique(Class)) < 2) {
    caughtUserError(paste("Your data should contain at least two classes\n",
                       "but your data have only one.\n"))
}
if(min(tclass) < 5) {
    caughtUserError(paste("At least one of your classes has less\n",
                       "than 5 cases/subjects/arrays. Although \n",
                       "the programs can deal with this, would you\n",
                       "believe it?\n"))
}
if(length(Class) != dim(xdata)[1]) {
    emessage <- paste("The class file and the covariate file\n",
                      "do not agree on the number of arrays: \n",
                      length(Class), " arrays according to the class file but \n",
                      dim(xdata)[1], " arrays according to the covariate data.\n",
                      "Please fix this problem and try again.\n")
    caughtUserError(emessage)  
}
if(!(is.numeric(xdata))) {
    caughtUserError("Your covariate file contains non-numeric data. \n That is not allowed.\n")
}
if(any(is.na(xdata))) {
    caughtUserError("Your covariate file contains missing values. \n That is not allowed.\n")
}
if(ncol(xdata) < 5) {
    caughtUserError("Your covariate file should contain at least five genes (otherwise, you probably don't need this tool)\n")
}
datos <- xdata
classes <- Class




#####  </end input>



num.cols.datos <- ncol(datos)

### for now, this is hard-wired:
num.genes <- c(2, 5, 10, 20, 35, 50, 75, 120, 200, 500, 1000, 2000)

if(length(which(num.genes > num.cols.datos)))
    num.genes <- num.genes[-which(num.genes > num.cols.datos)]

## to add all genes:
if(max(num.genes) < num.cols.datos)
    num.genes <- c(num.genes, num.cols.datos)


## if at least one class has size 11, we use 10-fold CV
## otherwise, we use the largest fold possible.
## has to be - 2, because I do another round of CV.
cv.error.k <- knumber <- min(10, max(table(classes)) - 2)



######  The actual calls for the CGI



the.predictor <- switch(classModel,
                        svm = "svmLinear",
                        dlda = "dldaC",
                        knn = "knn",
                        randomforest = "rF",
                        PAM = "PAM"
                        )
if(is.null(the.predictor)) stop("the.predictor is NULL")

the.geneselection <- switch(geneSel,
                        Fratio = "F",
                        Wilcoxon = "Wilcoxon",
                        randomforest = "rF"
                        )


the.geneselection.title <- switch(geneSel,
                        Fratio = "F statistic",
                        Wilcoxon = "Wilcoxon statistic",
                        randomforest = "Random Forest"
                        )


title.figure1 <- switch(classModel,
                        svm = "SVM (linear kernel)",
                        dlda = "DLDA",
                        knn = "NN",
                        randomforest = "Random forest",
                        shrunkencentroids = "Shrunken centroids (PAM)"
                        )

if(the.predictor != "PAM") {
    title.figure <- paste(title.figure1, "; gene selection using ",
                           the.geneselection.title, sep = "")
} else title.figure <- "Shrunken centroids (PAM)"

image.file <- "predictor_error_rates.png"
results.file <- "results.txt"
selected.genes.file <- "selectedgenes.txt"
nonselected.genes.file <- "nonselectedgenes.txt"

## utility function
n.print <- function(x) paste(substitute(x), paste(x, collapse = ", "),
                                 sep = " : ")


CairoPNG(image.file, width = 700, height = 500, ps = 10)


par(las = 1)
par(mar = c(5, 5, 4, 4) + 0.1)
par(cex.lab = 1.25)
par(cex.axis = 1)
par(cex.main = 1.25)


    
if(the.predictor != "PAM") {
    the.run <- try(f.pred.errors(datos, classes,
                                 sizeGenes = num.genes,
                                 typePred = the.predictor,
                                 geneSelection = the.geneselection,
                                 knumber = knumber,
                                 title.figure = title.figure,
                                 plot = TRUE,
                                 cv.error = TRUE))
} else {
    the.run <- try(f.pred.errors.PAM(datos, classes,
                                     plot = TRUE,
                                     cv.error = TRUE,
                                     cv.error.k = cv.error.k,
                                     title.figure = title.figure))
    
}

dev.off()


if(class(the.run) == "try-error") {
    save.image()
    caughtOurError(traceback())
}


if(the.predictor == "PAM") num.genes <- "Not available (can differ in each run)"
write(c("<h3>OPTIONS and SETTINGS</h3>",
        paste("<br>Sizes of numbers of genes in predictor examined: ", paste(num.genes, collapse = ", ")),
        paste("<br>Number of folds in cross-validation:             ", knumber),
        paste("<br>Classification algorithm:                        ", title.figure),
##        paste("<br>Gene selection:                                  ", the.geneselection.title),
        "<br><hr /><br>"),
      file = results.file)


HTML.cvPred(the.run, file = results.file, append = TRUE)

writeForPaLS(the.run)


##     These two calls used to be used for sending to FatiGO
write.table(geneNames[the.run[[4]]], sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE,
            file = selected.genes.file)

write.table(setdiff(colnames(datos), geneNames[the.run[[4]]]), sep = "\t",
            quote = FALSE,
            row.names = FALSE, col.names = FALSE,
            file = nonselected.genes.file)
save.image()
cat("\n\n Normal termination\n")



### checkdone is a piece of cake: incorporate figure and results.txt, and compress output.

