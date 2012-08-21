## Pretty much at what appears to be a useable state.  Need to manage the class structure/pass through to keep the amended cf object


cfConvert <- function(x, xpre, x2, x2pre, formula){
  xscen <- list(x = NULL, xpre = NULL)
  
  # If only cf data, no formula given (not already a cf object)
  if (!is.null(x) && !is.numeric(x)) stop("counterfactual must be numeric")  # this check needs to be extended across all of them
  

  # Checking & arranging x
  if (is.data.frame(x)) xscen$x <- x
  else if (is.matrix(x)) xscen$x <- as.data.frame(x)
  else if (is.list(x)){
    if (!is.null(x$x) && is.data.frame(x$x)) xscen$x <- x$x
    else if (!is.null(x$x) && is.matrix(x$x)) xscen$x <- as.data.frame(x$x)
    else stop("x must be either a matrix or dataframe of counterfactual scenarios, or must be a list containing a matrix or dataframe of counterfactual scenarios called by $x")
}
    # Arranging xpre
  if (!is.null(xpre) && is.list(x) && !is.null(x$xpre)) stop("xpre is specified twice, in 'x' list object and in 'xpre=' argument in simcf() call. Only one specification allowed.")
  if (is.list(x) && !is.null(x$xpre)) xpre <- x$xpre
  if (is.data.frame(xpre)) xscen$xpre <- xpre
  else if (is.matrix(xpre)) xscen$xpre <- as.data.frame(xpre)
  else if (is.null(xpre)) xscen$xpre <- NULL
  else stop("xpre  must be either a matrix or dataframe of counterfactual scenarios, or there must be a matrix or dataframe included in list 'x' called by x$xpre")
  
  
  #Arranging x2
  if (!is.null(x2) && is.list(x) && !is.null(x$x2)) stop("x2 is specified twice, in 'x' list object and in 'x2=' argument in simcf() call. 
                                           Only one specification allowed.")
  if (is.list(x) && !is.null(x$x2)) x2 <- x$x2
  if (is.data.frame(x2)) xscen$x2 <- x2
  else if (is.matrix(x2)) xscen$x2 <- as.data.frame(x2)
  else if (is.null(x2)) xscen$x2 <- NULL
  else stop("x2  must be either a matrix or dataframe of counterfactual scenarios, or there must be a matrix or dataframe included in list 'x' called by x$x2")
  
  if (!is.null(x2pre) && is.list(x) && !is.null(x$x2pre)) stop("x2pre is specified twice, in 'x' list object and in 'x2pre=' argument in simcf() call. Only one specification allowed.")
  if (is.list(x) && !is.null(x$x2pre)) x2pre <- x$x2pre
  if (is.data.frame(x2pre)) xscen$x2pre <- x2pre
  else if (is.matrix(x2pre)) xscen$x2pre <- as.data.frame(x2pre)
  else if (is.null(x2pre)) xscen$x2pre <- NULL
  else stop("x2pre  must be either a matrix or dataframe of counterfactual scenarios, or there must be a matrix or dataframe included in list 'x' called by x$x2pre")
  
           
  # Checking for congruence between formula, if included, and included cf columns  
  cfFormVarCheck <- function(x,model){
    tms<-terms(model)
    if (attr(tms,"response")==1 && is.null(x[[as.character(attr(tms,"variables")[[2]])]])) x[[as.character(attr(tms,"variables")[[2]])]] <- rep(NA,nrow(x))  
    return(x)
  }
  
  if (!is.null(formula) && !is.null(xscen$x)) xscen$x <- cfFormVarCheck(x=xscen$x, model=formula)
  if (!is.null(formula) && !is.null(xscen$xpre)) xscen$xpre <- cfFormVarCheck(x=xscen$xpre, model=formula)
  if (!is.null(formula) && !is.null(xscen$x2)) xscen$x2 <- cfFormVarCheck(x=xscen$x2, model=formula)
  if (!is.null(formula) && !is.null(xscen$x2pre)) xscen$x2pre <- cfFormVarCheck(x=xscen$x2pre, model=formula)
    
  if (!is.null(formula) && is.null(x)) xscen$x <- model.frame(formula)
                  
  if (!is.null(formula)) xscen$model <- formula
  
  xscen
}
  


validateCF <- function(xscen, qoi, xpre, x2, x2pre,constant,b,...){
  print(qoi)
  if (is.null(xscen$model) && !is.na(constant) && (ncol(xscen$x) != (ncol(b)-1) )) stop("If model formula is not provided, matrix of counterfactuals must include exactly one column for every beta coefficient (excluding the intercept).")
  if (is.null(xscen$model) && is.na(constant) && (ncol(xscen$x) != ncol(b)) ) stop("If model formula is not provided, matrix of counterfactuals must include exactly one column for every beta coefficient (excluding the intercept).")
  
  if (qoi=="fd"){
    if (is.null(xscen$x) || is.null(xscen$xpre)) stop("Both x and xpre scenarios must be provided to calculate first differences")
    else if ((dim(xscen$x)[1] != dim(xscen$spre)[1]) || (dim(xscen$x)[2] != dim(xscen$spre)[2])) stop("x and xpre must have the same dimensions to calculate first differences")
  }
  if (qoi=="rr"){
    if (is.null(xscen$x) || is.null(xscen$xpre)) stop("Both x and xpre scenarios must be provided to calculate relative risk")
    else if ((dim(xscen$x)[1] != dim(xscen$xpre)[1]) || (dim(xscen$x)[2] != dim(xscen$xpre)[2])) stop("x and xpre must have the same dimensions to calculate first differences")
  }
  if (qoi=="dd"){
    if (is.null(xscen$x) || is.null(xscen$xpre)) stop("Both x and xpre scenarios must be provided to calculate difference in differences")
    if (is.null(xscen$x2) || is.null(xscen$x2pre)) stop("Both x2 and x2pre scenarios must be provided (as well as x and xpre) to calculate difference in differences")
    else if ((dim(xscen$x)[1] != dim(xscen$xpre)[1]) || (dim(xscen$x)[2] != dim(xscen$xpre)[2])) stop("x, xpre, x2, and x2pre must all have the same dimensions to calculate difference in differences")
    else if ((dim(xscen$x2)[1] != dim(xscen$x2pre)[1]) || (dim(xscen$x2)[2] != dim(xscen$x2pre)[2])) stop("x, xpre, x2, and x2pre must all have the same dimensions to calculate difference in differences")
    else if ((dim(xscen$x)[1] != dim(xscen$x2)[1]) || (dim(xscen$x)[2] != dim(xscen$x2)[2])) stop("x, xpre, x2, and x2pre must all have the same dimensions to calculate difference in differences")  
  }
  
  if (qoi=="ev"){
    if (is.null(xscen$x)) stop("x scenarios must be provided to calculate expected values.  Include a model formula or matrix/dataframe of counterfactual x-values.")
    }
  
  if (is.null(xscen$model)) warning("Model formula not provided.  Any transformations of independent variables (e.g. logged variables) will not be automatically calculated and should have been done prior to simcf() call.")
  xscen
 }


appendmatrix <- function(x,values,after=ncol(x)) {
  nrowx <- nrow(x)
  ncolx <- ncol(x)
  after <- nrowx*(after-1)
  x <- as.vector(as.matrix(x))
  x <- append(x,values,after)
  matrix(x,nrow=nrowx,ncol=ncolx+1,byrow=FALSE)
}

modelMatrixHelper <- function(cfxs, b, model=NULL, constant, nscen=1, ...) { #cfxs is a single matrix of counterfactuals - x, xpre etc
  if (!is.null(model)){
    mm <- model.matrix(model,cfxs)
  }
  
  else if (!is.matrix(cfxs)) {
    if (is.vector(cfxs)){
      if (is.matrix(b)) {
        mm <- t(cfxs)
        if (!is.na(constant)) {
          mm <- append(mm, 1, constant - 1)
        }
      }
    } else {
      mm <- as.matrix(cfxs)
      if (!is.na(constant)) {
        mm <- appendmatrix(mm, rep(1, nrow(cfxs)), constant)
      }
    }
  } else {
    mm <- cfxs
    if (!is.na(constant)) {
      mm <- appendmatrix(mm, rep(1, nrow(cfxs)), constant)
    }
  }
  mm
}



simcf <- function(x,b, qoi="ev", type="linear", formula=NULL, ci = 0.95, constant = 1, xpre = NULL, x2=NULL, x2pre=NULL, ...){

  if ( !is.na(constant) ) {
    if (length(constant) > 1 || !is.numeric(constant) || constant <= 0)
      stop("constant must be positive scalar (which column of b has model constant?)")
  }
  
  # convert to CF - cfConvert
  if (!is.cf(x)) xscen <- cfConvert(x, xpre, x2, x2pre, formula)
  else xscen<-x
  
  # check input types for each qoi type - function validateCF()
  xscen <- validateCF(xscen,qoi=qoi,constant=constant,b=b)
       
 # make matrices for simulation 
  nscen <- nrow(xscen$x)
  
  x <- modelMatrixHelper(cfxs=xscen$x,b=b,constant=constant, model=xscen$model,nscen=nscen)
  if (!is.null(xscen$xpre)) xpre <-  modelMatrixHelper(cfxs=xscen$xpre,b=b,constant=constant, model=xscen$model, nscen=nscen)
  if (!is.null(xscen$x2)) x2 <-  modelMatrixHelper(cfxs=xscen$x2,b=b,constant=constant, model=xscen$model, nscen=nscen)
  if (!is.null(xscen$x2pre)) x2pre <-  modelMatrixHelper(cfxs=xscen$x2pre,b=b,constant=constant, model=xscen$model, nscen=nscen)
 
 # set up for simulated ys output
  esims <- nrow(as.matrix(b))
  nscen <- nrow(xscen$x)
  nci <- length(ci)
  res <- list(pe = rep(NA, nscen), lower = matrix(NA, nrow = nscen, 
                                                  ncol = nci), upper = matrix(NA, nrow = nscen, ncol = nci))
  
  # Identify type of model simulation calculation

  for (i in 1:nscen) {
    if (qoi=="ev") {
      simy <- simfn(x=(x[i,]),b, type)      
    }
    else if(qoi=="fd"){
      simpre <- simfn(xhyp=(xpre[i,]), b, type)  
      simpost<- simfn(xhyp=(x[i,]), b, type)
      simy <- simpost - simpre      
    } 
    else if (qoi=="rr"){
      simpre <- simfn(xpre[i,],b,type)  
      simpost<- simfn(x[i,],b,type)
      simy <- simpost/simpre
    } 
    else if (qoi=="dd"){
      sim1pre <- simfn(xpre[i,],b,type)  
      sim1post<- simfn(x[i,],b,type)
      sim1y <- sim1post - sim1pre
      sim2pre <- simfn(x2pre[i,],b,type)  
      sim2post<- simfn(x2[i,],b,type)
      sim2y <- sim2post - sim2pre
      simy <- sim2y-sim1y        
    }
    
    res$pe[i] <- mean(simy)
    for (k in 1:nci) {
      cint <- quantile(simy, probs = c((1 - ci[k])/2, (1 - 
        (1 - ci[k])/2)))
      res$lower[i, k] <- cint[1]
      res$upper[i, k] <- cint[2]
    }
  }
  
  res$lower <- drop(res$lower)
  res$upper <- drop(res$upper)
  
  xscen <- append(xscen,res)
  class(xscen) <- c("postsimcf","counterfactual","list")
  xscen
  
}
  