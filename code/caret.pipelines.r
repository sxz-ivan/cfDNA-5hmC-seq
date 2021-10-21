suppressPackageStartupMessages({
  require(caret)
  require(tidyverse)
})

caret_trainer = function(matrix,y,model = 'rf',trControl = trainControl("cv", number = 10),tuneLength = 5,...){
  data = matrix %>% as.data.frame %>% within({y =y })
  {sink("/dev/null");
    model = train(y ~ ., data = data, method = model,trControl = trControl,tuneLength = tuneLength,...);
    sink();}
  model
}

caret_predictor = function(model,X,...){
  predict(model, X,...)
}