simfn <- function(xhyp,b=b,type=type,...){
  print(paste("length xhyp:", length(xhyp)))
  print(paste("dim b:", dim(b)))
  if ( !is.vector(xhyp) ) stop("x must be a vector of variable values for a single scenario.")
  
  if (type=="linear") {
    print("flag C")
    simmu <- b%*%xhyp
    print(length(simmu))
    print(paste("mean simmu:",mean(simmu)))
    return(simmu)
  }
  
  if (type=="logit") {
    simmu <- b%*%xhyp
    simlt <- 1/(1+ exp(-simmu))
    simlt
  }
}

  
linsim <- function(xhyp,b=b,...){

  simmu <- b%*%x
  simmu
}

lgtsim <- function(x,b=b,...){
  if ( !is.vector(x) ) stop("x must be a vector of variable values for a single scenario.")
  simmu <- x%*%b
  simlt <- 1/(1+ exp(-simmu))
  simlt
}



  