#' Application of each ML method to a testing set
#'
#' Returns the estimated accurasy measures of Sensitivity, Specificity, and Misclassification using Confusion matrix.
#'
#' This is a generic function. The candidate predictors are a series of forward-in variables which can be picked by the variables ranked by
#'  importance measures or user-defined ones in `depvar` parameter.
#'  Their coefficients are estimated or the decision tree is built by by the whole training data.
#'
#' @param data A matrix of the whole training data set for coefficients' esimates
#' @param newdata A matrix of a testing data
#' @param indvar A list of candidate predictors
#' @param depvar An outcome variable
#' @param importance_measure A vector of importance measures obtained from Random Forest model
#' @param ML A ML method to apply
#' @param add.vars ??
#' @param rdit ??
#' @return ML.err.vect A data frame including estimated accuracy measures
#' @export

ML_error <- function(data, newdata, importance_measure, add.vars = NULL, rdit = 4, depvar = depVar, indvar = indVar, ML = c("RF", "BOOST", "Logistic")) {
    if (is.null(importance_measure)) {
        importance_measure <- NULL
    } else {
        importance_measure <- lapply(1:length(importance_measure), function(x) {
            importance_measure[1:x]
        })
        ML.err.vect <- data.frame(matrix(rep(0, 1), nrow = length(importance_measure)))
        names(ML.err.vect) <- "Variables"
        ML.Confusion <- vector("list", length(importance_measure))
    }

    if (ML == "Logistic") {
        modelFit <- vector("list", length(importance_measure))
    }
    if (!is.null(add.vars)) {
        importance_measure <- c(importance_measure, add.vars)
        modelFit <- c(importance_measure, vector("list", length(add.vars)))
        ML.err.vect <- data.frame(matrix(rep(0, 1), nrow = (length(add.vars) + length(importance_measure))))
        names(ML.err.vect) <- "Variables"

        ML.Confusion <- vector("list", (length(add.vars) + length(importance_measure)))
    }

    for (i in 1:length(importance_measure)) {
        if (ML == "RF") {
            RF.tempVars <- (importance_measure)[[i]]
            RF.fit <- randomForest(as.factor(PTSD_binary) ~ ., data = data[, which(names(data) %in% c(depVar, RF.tempVars))])
            RF.pred <- predict(RF.fit, newdata = newdata[, -which(names(newdata) %in% depVar)], type = "prob")[, 2]
            ML.err.vect[i, "ML.err"] <- round(roc.area(as.numeric(newdata[, depVar]) - 1, RF.pred)$A, digits = rdit)
            ML.Confusion[[i]] <- confusionMatrix(as.factor(ifelse(RF.pred > 0.5, 1, 0)), newdata[, which(names(newdata) %in% depVar)])
            ML.err.vect[i, 1] <- paste0(RF.tempVars, collapse = ", ")
        }

        if (ML == "BOOST") {
            BOOST.tempVars <- (importance_measure)[1:i]
            BOOST.fit <- gbm((PTSD_binary) ~ ., data = data[, which(names(data) %in% c(depVar, BOOST.tempVars))], distribution = "bernoulli",
                n.trees = 10000, shrinkage = 0.01, interaction.depth = 2)
            BOOST.pred <- predict(BOOST.fit, newdata = newdata[, -which(names(newdata) %in% depVar)], n.trees = 10000, type = "response")
            ML.err.vect[i, "ML.err"] <- round(roc.area(as.numeric(newdata[, depVar]) - 1, Logistic.pred)$A, digits = rdit)
            ML.Confusion[[i]] <- confusionMatrix(as.factor(ifelse(BOOST.pred > 0.5, 1, 0)), newdata[, which(names(newdata) %in% depVar)])
            ML.err.vect[i, 1] <- paste0(BOOST.tempVars, collapse = ", ")
        }

        if (ML == "Logistic") {
            Logistic.tempVars <- (importance_measure)[1:i]
            # Logistic.fit <- glm(as.factor(PTSD_binary)~., data = data[, which(names(data)%in%c(depVar, Logistic.tempVars))] , family =
            # 'binomial' )
            Logistic.fit <- bayesglm(as.factor(PTSD_binary) ~ ., data = data[, which(names(data) %in% c(depVar, indVar, Logistic.tempVars))],
                family = "binomial")
            Logistic.pred <- predict(Logistic.fit, newdata = newdata[, -which(names(newdata) %in% depVar)], type = "response")
            Logistic.pred.binary <- ifelse(Logistic.pred < 0.5, 0, 1)
            ML.err.vect[i, "ML.err"] <- round(roc.area(as.numeric(newdata[, depVar]) - 1, Logistic.pred)$A, digits = rdit)
            ML.Confusion[[i]] <- confusionMatrix(as.factor(ifelse(Logistic.pred > 0.5, 1, 0)), newdata[, which(names(newdata) %in%
                depVar)])
            ML.err.vect[i, 1] <- paste0(Logistic.tempVars, collapse = ", ")

            modelFit[[i]] <- Logistic.fit
        }
    }

    # if(!is.null(add.vars)){ for (i in 1:length(add.vars)){ if(ML=='Logistic'){ Logistic.tempVars <- add.vars[[i]] # Logistic.fit <-
    # glm(as.factor(PTSD_binary)~., data = data[, which(names(data)%in%c(depVar, Logistic.tempVars))] , family = 'binomial' )
    # Logistic.fit <- bayesglm(as.factor(PTSD_binary)~., data = data[, which(names(data)%in%c(depVar, Logistic.tempVars))] , family =
    # 'binomial' ) Logistic.pred <- predict(Logistic.fit, newdata=newdata[, -which(names(newdata)%in%depVar)], type='response')
    # Logistic.pred.binary <- ifelse(Logistic.pred<0.5, 0, 1) ML.err.vect[length(importance_measure)+i, 'ML.err'] <-
    # round(roc.area(as.numeric(newdata[,depVar])-1, Logistic.pred)$A, digits=rdit) ML.Confusion[[length(importance_measure)+i]] <-
    # confusionMatrix(as.factor(ifelse(Logistic.pred>0.5, 1, 0)),newdata[, which(names(newdata)%in%depVar)])
    # ML.err.vect[length(importance_measure)+i,1] <- paste0(Logistic.tempVars, collapse = ', ')
    # modelFit[[length(importance_measure)+i]] <- Logistic.fit } } }
    ML.err.vect$sensitivity <- round(unlist(lapply(ML.Confusion, function(x) {
        x$byClass["Sensitivity"]
    })), digits = rdit)
    ML.err.vect$specificity <- round(unlist(lapply(ML.Confusion, function(x) {
        x$byClass["Specificity"]
    })), digits = rdit)
    ML.err.vect$misclassification <- round(unlist(lapply(ML.Confusion, function(x) {
        (x$table[1, 2] + x$table[2, 1])/(sum(x$table))
    })), digits = rdit)

    if (ML == "Logistic") {
        return(list(ML.err.vect = ML.err.vect, model.fit = modelFit))
    } else {
        return(ML.err.vect)
    }
}
