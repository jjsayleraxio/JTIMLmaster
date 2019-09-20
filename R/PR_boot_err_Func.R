#' Application of Penalized Regression method
#'
#' Returns the estimated accurasy measures of Sensitivity, Specificity, and Misclassification using Confusion matrix.
#'
#' A Penalized Regression model is trained using a training data set and estimate accuracy measures using a test data set. The number of tree is defined by a user and the default value of it is 300.
#'
#' @param x.train A data frame of training data set for independent variables.
#' @param y.train A data frame of training data set for dependent variable.
#' @param x.test A data frame of testing data set for independent variables.
#' @param y.test A data frame of testing data set for dependent variable.
#' @param alphaLevel User defined Alpha tunning parameter setting.
#' @param family A distribution of the outcome variable.
#' @param type Is this a classification or regression?
#' @return A list of PR.fit PR fitted object. A non0coeff is the estimated non-zerio coefficients. A PR.err.vect a matrix of the estimated accurasy measures of Sensitivity, Specificity, and Misclassification for Training set, Testing set and overall one.
#' @export

PR.boot.err.Func <- function(x.train, y.train, x.test, y.test, alphaLevel, family = "binomial", type = "class") {
    Model <- paste0("Elastic Net(alpha=", alphaLevel, ")")
    Model[which(alphaLevel == 1)] <- "LASSO(alpha=1)"
    x.train <- as.matrix(x.train)
    x.test <- as.matrix(x.test)

    PR.err.vect <- data.frame(Model = Model)
    cvGLMNET <- lapply(1:length(alphaLevel), function(x) {
        glmnet::cv.glmnet(x.train, y.train, type.measure = "mse", alpha = alphaLevel[x], family = family)
    })
    preClassification.test <- lapply(cvGLMNET, function(x) {
        factor(as.vector(predict(x, s = x$lambda.min, newx = x.test, type = type)), levels = c(0, 1))
    })
    preClassification.train <- lapply(cvGLMNET, function(x) {
        factor(as.vector(predict(x, s = x$lambda.min, newx = x.train, type = type)), levels = c(0, 1))
    })

    coefs <- lapply(cvGLMNET, function(x) {
        stats::coef(x, s = "lambda.min")
    })
    non0coeff <- lapply(coefs, function(x) {
        x[which(x[, 1] != 0), ]
    })

    PR.Confusion.test <- lapply(preClassification.test, function(x) {
        caret::confusionMatrix(x, y.test)
    })

    ## For test sets
    PR.err.vect$AUC.test <- unlist(lapply(preClassification.test, function(x) {
        verification::roc.area(as.numeric(y.test) - 1, as.numeric(x) - 1)$A
    }))
    PR.err.vect$sensitivity.test <- unlist(lapply(PR.Confusion.test, function(x) {
        x$byClass["Sensitivity"]
    }))
    PR.err.vect$specificity.test <- unlist(lapply(PR.Confusion.test, function(x) {
        x$byClass["Specificity"]
    }))
    PR.err.vect$misclassification.test <- unlist(lapply(PR.Confusion.test, function(x) {
        (x$table[1, 2] + x$table[2, 1])/(sum(x$table))
    }))

    ## For train sets
    PR.Confusion.train <- lapply(preClassification.train, function(x) {
        caret::confusionMatrix(x, y.train)
    })
    PR.err.vect$AUC.train <- unlist(lapply(preClassification.train, function(x) {
        verification::roc.area(as.numeric(y.train) - 1, as.numeric(x) - 1)$A
    }))
    PR.err.vect$sensitivity.train <- unlist(lapply(PR.Confusion.train, function(x) {
        x$byClass["Sensitivity"]
    }))
    PR.err.vect$specificity.train <- unlist(lapply(PR.Confusion.train, function(x) {
        x$byClass["Specificity"]
    }))
    PR.err.vect$misclassification.train <- unlist(lapply(PR.Confusion.train, function(x) {
        (x$table[1, 2] + x$table[2, 1])/(sum(x$table))
    }))

    PR.err.vect$AUC.overall <- PR.err.vect$AUC.test * 0.632 + PR.err.vect$AUC.train * 0.368
    PR.err.vect$sensitivity.overall <- PR.err.vect$sensitivity.test * 0.632 + PR.err.vect$sensitivity.train * 0.368
    PR.err.vect$specificity.overall <- PR.err.vect$specificity.test * 0.632 + PR.err.vect$specificity.train * 0.368
    PR.err.vect$misclassification.overall <- PR.err.vect$misclassification.test * 0.632 + PR.err.vect$misclassification.train * 0.368

    PR.err.vect$non.zero.num.var <- unlist(lapply(non0coeff, length))
    rownames(PR.err.vect) <- PR.err.vect$Model
    PR.err.vect$Model <- NULL
    PR.err.vect <- as.data.frame(t(PR.err.vect))

    return(list(PR.fit = cvGLMNET, non0coeff = non0coeff, PR.err.vect = PR.err.vect))
}
