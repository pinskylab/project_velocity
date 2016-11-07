#A function for removing info from a model object that is not strictly needed for prediction. This is to reduce size of the saved object and should be done after extracting all relevant summary stats. (Not currently working for GAMs - too many attributes needed for prediction.)

stripGAM = function(cm) {
  cm$y = c()
#  cm$model = c()
  
  cm$residuals = c()
  cm$fitted.values = c()
#  cm$effects = c()
#  cm$qr$qr = c()  
  cm$linear.predictors = c()
  cm$weights = c()
  cm$prior.weights = c()
#  cm$data = c()

  
  cm$family$variance = c()
  cm$family$dev.resids = c()
  cm$family$aic = c()
  cm$family$initialize = c()
#  cm$family$validmu = c()
  cm$family$simulate = c()
  attr(cm$terms,".Environment") = c()
  attr(cm$formula,".Environment") = c()
  
  cm
}
