#' Application of SVM method
#'
#' Returns the estimated accurasy measures of Sensitivity, Specificity, and Misclassification using Confusion matrix.
#'
#' A SVM model is trained using a training data set and estimate accuracy measures using a test data set. The number of tree is defined by a user and the default value of it is 300.
#'
#' @inheritParams RF.boot.err.Func A data frame of training data set
#' @param kernel A kernel function to be used.
#' @param cost A cost?
#' @param train ??
#' @param test ??
#' @param depVar ??
#' @return A SVM.err a matrix of the estimated accurasy measures of Sensitivity, Specificity, and Misclassification for Training set, Testing set and overall one.
#' @export

SVM.boot.err.Func <- function(train, test, depVar, kernel = "radial", cost = 5) {
    variables <- names(train)
    variables <- variables[variables!=depVar]
    f <- as.formula(paste(depVar, paste(variables, collapse = " + "), sep = " ~ "))

    fit <- e1071::svm(f, data = train)
    Prediction.test <- predict(fit, test)
    Prediction.train <- predict(fit, train)

    ## For test sets
    AUC.test <- verification::roc.area(as.numeric(test[, which(names(test) %in% depVar)]) - 1, as.numeric(Prediction.test) - 1)$A
    Conf.test.Mat <- caret::confusionMatrix(as.factor(Prediction.test), test[, which(names(test) %in% depVar)])
    Sensitivity.test <- Conf.test.Mat$byClass["Sensitivity"]
    Specificity.test <- Conf.test.Mat$byClass["Specificity"]
    Misclassification.test <- (Conf.test.Mat$table[1, 2] + Conf.test.Mat$table[2, 1])/(sum(Conf.test.Mat$table))

    SVM.err.test <- data.frame(Test.Measure = c(AUC.test, Sensitivity.test, Specificity.test, Misclassification.test))
    rownames(SVM.err.test) <- c("AUC.test", "Sensitivity.test", "Specificity.test", "Misclassification.test")

    ## For train sets
    AUC.train <- verification::roc.area(as.numeric(train[, which(names(train) %in% depVar)]) - 1, as.numeric(Prediction.train) -
        1)$A
    Conf.train.Mat <- caret::confusionMatrix(as.factor(Prediction.train), train[, which(names(train) %in% depVar)])
    Sensitivity.train <- Conf.train.Mat$byClass["Sensitivity"]
    Specificity.train <- Conf.train.Mat$byClass["Specificity"]
    Misclassification.train <- (Conf.train.Mat$table[1, 2] + Conf.train.Mat$table[2, 1])/(sum(Conf.train.Mat$table))

    SVM.err.train <- data.frame(Train.Measure = c(AUC.train, Sensitivity.train, Specificity.train, Misclassification.train))
    rownames(SVM.err.train) <- c("AUC.train", "Sensitivity.train", "Specificity.train", "Misclassification.train")

    SVM.err <- cbind(SVM.err.test, SVM.err.train)
    SVM.err$overall <- SVM.err$Train.Measure * 0.368 + SVM.err$Test.Measure * 0.632

    return(list(SVM.err = SVM.err))
}
