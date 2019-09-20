#' Application of Random Forest method
#'
#' Returns the estimated accurasy measures of Sensitivity, Specificity, and Misclassification using Confusion matrix.
#'
#' A random forest model is trained using a training data set and estimate accuracy measures using a test data set. The number of tree is defined by a user and the default value of it is 300.
#'
#' @param train A data frame of training data set
#' @param test A data frame of a testing data
#' @param depvar An outcome variable
#' @param ntree A numeric value of # of trees to build
#' @return A list of RF.fit RF fitted object, RF.err a matrix of the estimated accurasy measures of Sensitivity, Specificity, Misclassification for Training set, Testing set and overall one,
#'   importance_measure as a matrix form of Mean.Decrease.Accuracy and Mean.Decrease.Gini
#' @export

RF.boot.err.Func <- function(train, test, depVar, ntree = 300) {
    variables <- names(train)
    variables <- variables[variables!=depVar]
    f <- as.formula(paste(depVar, paste(variables, collapse = " + "), sep = " ~ "))
    fit <- randomForest(f, data = train, importance = TRUE, ntree = ntree)
    Prediction.test <- predict(fit, test)
    Prediction.train <- predict(fit, train)

    ## For test sets
    AUC.test <- verification::roc.area(as.numeric(test[, which(names(test) %in% depVar)]) - 1, as.numeric(Prediction.test) - 1)$A
    Conf.test.Mat <- confusionMatrix(as.factor(Prediction.test), test[, which(names(test) %in% depVar)])
    Sensitivity.test <- Conf.test.Mat$byClass["Sensitivity"]
    Specificity.test <- Conf.test.Mat$byClass["Specificity"]
    Misclassification.test <- (Conf.test.Mat$table[1, 2] + Conf.test.Mat$table[2, 1])/(sum(Conf.test.Mat$table))

    RF.err.test <- data.frame(Test.Measure = c(AUC.test, Sensitivity.test, Specificity.test, Misclassification.test))
    rownames(RF.err.test) <- c("AUC.test", "Sensitivity.test", "Specificity.test", "Misclassification.test")

    ## For train sets
    AUC.train <- verification::roc.area(as.numeric(train[, which(names(train) %in% depVar)]) - 1, as.numeric(Prediction.train) -
        1)$A
    Conf.train.Mat <- confusionMatrix(as.factor(Prediction.train), train[, which(names(train) %in% depVar)])
    Sensitivity.train <- Conf.train.Mat$byClass["Sensitivity"]
    Specificity.train <- Conf.train.Mat$byClass["Specificity"]
    Misclassification.train <- (Conf.train.Mat$table[1, 2] + Conf.train.Mat$table[2, 1])/(sum(Conf.train.Mat$table))

    RF.err.train <- data.frame(Train.Measure = c(AUC.train, Sensitivity.train, Specificity.train, Misclassification.train))
    rownames(RF.err.train) <- c("AUC.train", "Sensitivity.train", "Specificity.train", "Misclassification.train")

    RF.err <- cbind(RF.err.test, RF.err.train)
    RF.err$overall <- RF.err$Train.Measure * 0.368 + RF.err$Test.Measure * 0.632

    imp.Mat <- merge(importance(fit, type = 1), importance(fit, type = 2), by = "row.names", all = TRUE)  ## Importance measure
    rownames(imp.Mat) <- imp.Mat$Row.names
    imp.Mat$Row.names <- NULL
    names(imp.Mat) <- c("Mean.Decrease.Accuracy", "Mean.Decrease.Gini")
    imp.Mat[, 1:2] <- round(imp.Mat[, 1:2], 4)
    imp.Mat <- imp.Mat[order(imp.Mat$Mean.Decrease.Accuracy, decreasing = TRUE), ]

    return(list(RF.fit = fit, RF.err = RF.err, importance_measure = imp.Mat))
}
