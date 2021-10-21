source('code/hxl.selectors.r')

oneside_select = function(x,v,n){
    n = min(n,length(x))
    x[sort(order(v)[1:n])]
}

twosides_select = function(x,v,n){
    n = min(n,floor(length(x)/2))
    x[sort(order(v)[c(1:n,(length(x)-n+1):length(x))])] %>% unique
}
glmnet1 = function(X,y,...){
    caret_trainer(
        matrix = X, y = y, 
        trControl = trainControl(
            method = "repeatedcv",number = 10, repeats =3,
            verboseIter=TRUE,returnData=FALSE,classProbs = TRUE,savePredictions=FALSE,sampling='down'), 
        model = 'glmnet', 
        tuneGrid = expand.grid(.alpha=c(0,0.01,0.05,0.1,0.2,0.5,0.8,1),.lambda = c(0,0.01,0.05,0.1,0.2,0.5)), 
        metric = "Kappa",tuneLength =10,...)
}

pc_glmnet=function (X, y, Fitfunc){
    Fitfunc(X, y, preProcess = c("center", "scale", "pca"))
}
glmnet2_oneSE10 = function (X, y,...){
    caret_trainer(matrix = X, y = y, trControl = trainControl(
        method = "repeatedcv", 
        number = 10, repeats = 3, verboseIter = TRUE, returnData = FALSE, 
        classProbs = TRUE, savePredictions = FALSE, sampling = "down",
        summaryFunction = twoClassSummary,
        selectionFunction = 'oneSE10'
    ), 
        model = "glmnet", tuneGrid = expand.grid(.alpha = c(0, 
            0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 1), .lambda = c(0, 
            0.01, 0.05, 0.1, 0.2, 0.5)), tuneLength = 10, metric = "ROC",...)
}

oneSE10 = function(x, metric, maximize){
    oneSE(x, metric, num=10, maximize)
}
coef.pc_glmnet = function(model){
    model %>% with(coef(finalModel,s=bestTune$lambda)) %>% as.matrix %>% .[-1,1]
}

RHS.pc_glmnet = function(model,X){
    trans.X = predict(model$preProcess,X)[,colnames(model$preProcess$rotation)]
    coefs = coef.pc_glmnet(model)
    as.numeric(trans.X %*% coefs)
}
RHS.glmnet = function (model, X) 
{
    if(startsWith(model$coefnames[1],'`') & !startsWith(colnames(X)[1],'`')){
        colnames(X) = sprintf('`%s`',colnames(X))
    }
    X = X[,model$coefnames]
    coefs = coef.pc_glmnet(model)
    if(!is.null(model$preProcess)){
        X = predict(model$preProcess, X)[, names(coefs)]
    }
    as.numeric(X %*% coefs)
}
caret_varImps = function(coefs){
    100*abs(coefs)/max(abs(coefs))
}

varImp.pc_glmnet = function(model){
    caret_varImps(coef2.pc_glmnet(model))
}

coef2.pc_glmnet = function(model){
    (model$preProcess$rotation %*% coef.pc_glmnet(model))[,1]
}