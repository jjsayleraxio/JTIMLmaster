#' Application of Neural Network method
#'
#' Returns the estimated accurasy measures of Sensitivity, Specificity, and Misclassification using Confusion matrix.
#'
#' A Neural Network model is trained using a training data set and estimate accuracy measures using a test data set. The number of tree is defined by a user and the default value of it is 300.
#'
#' @inheritParams RF.boot.err.Func A data frame of training data set
#' @param hidden A numeric value of # of hidden layers
#' @param train ??
#' @param test ??
#' @param depVar ??
#' @param linear.output ??
#' @return A NN.err a matrix of the estimated accurasy measures of Sensitivity, Specificity, and Misclassification for Training set, Testing set and overall one.
#' @export

NN.boot.err.Func <- function(train, test, depVar, hidden = 3, linear.output = FALSE) {

    Outcome <- as.numeric(train[, names(train)%in%depVar]) - 1
    max = apply(train[, names(train) != depVar], 2, max)
    min = apply(train[, names(train) != depVar], 2, min)
    scaled = as.data.frame(scale(train[, names(train) != depVar], center = min, scale = max - min))

    train <- cbind(Outcome, scaled)
    names(train)[1] <- depVar
    n <- names(train)
    n <- n[n!=depVar]
    f <- stats::as.formula(paste(depVar, paste(n, collapse = " + "), sep = " ~ "))
    fit <- neuralnet::neuralnet(f, data = train)

    Outcome <- as.numeric(test[, names(test)%in%depVar]) - 1
    max = apply(test[, names(test) != depVar], 2, max)
    min = apply(test[, names(test) != depVar], 2, min)
    scaled = as.data.frame(scale(test[, names(test) != depVar], center = min, scale = max - min))

    test <- cbind(Outcome, scaled)
    names(test)[1] <- depVar

    Prediction.test <- neuralnet::compute(fit, test[, which(names(test) != depVar)])
    Prediction.train <- neuralnet::compute(fit, train[, which(names(train) != depVar)])

    ## For test sets
    AUC.test <- verification::roc.area(as.numeric(test[, which(names(test) %in% depVar)]), as.numeric(Prediction.test$net.result[,
        1]))$A
    test.out <-
    Conf.test.Mat <- caret::confusionMatrix(factor(ifelse(Prediction.test$net.result[, 1] > 0.5, 1, 0)),
                                     factor(test[, which(names(test) %in% depVar)]))
    Sensitivity.test <- Conf.test.Mat$byClass["Sensitivity"]
    Specificity.test <- Conf.test.Mat$byClass["Specificity"]
    Misclassification.test <- (Conf.test.Mat$table[1, 2] + Conf.test.Mat$table[2, 1])/(sum(Conf.test.Mat$table))

    NN.err.test <- data.frame(Test.Measure = c(AUC.test, Sensitivity.test, Specificity.test, Misclassification.test))
    rownames(NN.err.test) <- c("AUC.test", "Sensitivity.test", "Specificity.test", "Misclassification.test")

    ## For train sets
    AUC.train <- verification::roc.area(as.numeric(train[, which(names(train) %in% depVar)]), as.numeric(Prediction.train$net.result[,
        1]))$A
    Conf.train.Mat <- caret::confusionMatrix(as.factor(ifelse(Prediction.train$net.result[, 1] > 0.5, 1, 0)), as.factor(train[, which(names(train) %in%
        depVar)]))
    Sensitivity.train <- Conf.train.Mat$byClass["Sensitivity"]
    Specificity.train <- Conf.train.Mat$byClass["Specificity"]
    Misclassification.train <- (Conf.train.Mat$table[1, 2] + Conf.train.Mat$table[2, 1])/(sum(Conf.train.Mat$table))

    NN.err.train <- data.frame(Train.Measure = c(AUC.train, Sensitivity.train, Specificity.train, Misclassification.train))
    rownames(NN.err.train) <- c("AUC.train", "Sensitivity.train", "Specificity.train", "Misclassification.train")

    NN.err <- cbind(NN.err.test, NN.err.train)
    NN.err$overall <- NN.err$Train.Measure * 0.368 + NN.err$Test.Measure * 0.632

    return(list(NN.err = NN.err))
}
