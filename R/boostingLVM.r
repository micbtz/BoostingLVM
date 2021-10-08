

# probabilities of a probit model
probs<-function(param,nodes,D)
{
  alpha<-param[1:D] # discriminations
  beta<-param[D+1] # intercepts
  lp<-matrix(alpha,nrow=1)%*%t(nodes)+beta
  pnorm(lp)
}


# returns unique patterns of responses and relative frequencies
reduce_data<-function(data)
{
  patterns<-apply(data,1,paste, collapse ="_")
  tab<-table(patterns)
  first_patt<-!duplicated(patterns)
  datared<-as.matrix(data[first_patt,])
  patterns_unique<-patterns[first_patt]
  numpatt<-as.matrix(tab[patterns_unique])
  return(list(data=datared,numpatt=numpatt))
}

NormalOgiveLik<-function(par,data,D=1,nodes,weights,fixrotat=FALSE)
{
  patterns<-apply(data,1,paste, collapse ="_")
  tab<-table(patterns)
  first_patt<-!duplicated(patterns)
  datared<-data[first_patt,]
  patterns_unique<-patterns[first_patt]
  numpatt<-as.vector(tab[patterns_unique])
  
  n<-nrow(datared)
  nitems<-ncol(datared)
  param<-matrix(par,nrow=D+1)
  
  if (fixrotat)
  {
    uptri<-upper.tri(matrix(NA,D,D))
    uptri<-cbind(matrix(FALSE,D,nitems-D),uptri)
    uptri<-rbind(uptri,FALSE)
    param[uptri]<-0
  }
  probs<- apply(param,2, FUN=probs, nodes=nodes, D=D)

  logprodj <- matrix(0, n, nrow(nodes))
  for (j in 1:nitems) {
    probsj <- probs[,j]
    probsj <- rbind(1-probsj, probsj)
    logprobsj <- log(probsj)
    xj <- datared[, j]
    na.ind <- is.na(xj)
    logprobsj <- logprobsj[xj+1, ]
    if (any(na.ind))
      logprobsj[na.ind, ] <- 0
    logprodj <- logprodj + logprobsj
  }
  prodj<-exp(logprodj)
  sumq <- (prodj %*% weights)
  mlik <- -sum(log(sumq)*numpatt) # minus log-likelihood
  mlik
}



simdatalvm<-function(n,param,Phi=NULL,type="binary")
{
  nitems<-ncol(param)
  D<-nrow(param)-1
  if (is.null(Phi)) Phi<-diag(1,D)
  Lambda<-t(param[1:D,])
  if (D==1) Lambda<-t(Lambda)
  thresholds<-param[D+1,]
  vary<-Lambda%*%Phi%*%t(Lambda)
  if (any(diag(vary)>1)) stop("variance of y greater than 1")
  diag(vary)<-1
  ys<-rmvnorm(n=n,mean=rep(0,nitems),sigma=vary)
  for (i in 1:nitems) ys[,i]<-ys[,i]-thresholds[i]
  y<-(ys>0)*1
  y<-as.data.frame(y)
  colnames(y)<-paste("I",1:nitems,sep="")
  y
}

fa2<-function(data,D,maxit=5000,eta=0,fixrotat=FALSE)
{
  nitems<-ncol(data)
  bifreq<-compute_bifreq(data)
  fa<-factanal(data,factors = D)
  load<-t(matrix(fa$loadings,ncol=D))
  intercepts<-rep(NA,nitems)
  for (i in 1:nitems) intercepts[i]<- -qnorm(mean(data[,i]))
  param<-rbind(load,intercepts)
  opt<-optim(par=as.vector(param),fn=lik_fa2Rcpp,gr=grad_lik_fa2Rcpp,bifreq=bifreq,nitems=nitems,D=D,control = list(maxit=maxit),method="BFGS",eta=eta,fixrotat=fixrotat)
  opt$par<-matrix(opt$par,D+1,nitems)
  if (fixrotat)
  {
    uptri<-upper.tri(matrix(NA,D,D))
    uptri<-cbind(matrix(FALSE,D,nitems-D),uptri)
    uptri<-rbind(uptri,FALSE)
    opt$par[uptri]<-0
  }
  opt
}


compute_bifreq<-function(data)
{
  nitems<-ncol(data)
  bifreq<-list()
  for(i in 1:(nitems-1))
  {
    bifreq[[i]]<-list()
    for (j in (i+1):nitems)
    {
      tab<-matrix(NA,2,2)
      tab[1,1]<-sum((1-data[,i])*(1-data[,j]))
      tab[1,2]<-sum((1-data[,i])*data[,j])
      tab[2,1]<-sum(data[,i]*(1-data[,j]))
      tab[2,2]<-sum(data[,i]*data[,j])
      bifreq[[i]][[j]]<-tab
    }
  }
  bifreq
}




lvmboost <- function(data=NULL,trace=TRUE,nu=0.05,tau=2,Dmax=ncol(data)-1,kstop=100,
                     obj.lvmboost=NULL,nfreepar=TRUE,eta=0)
{
  if (is.null(data) & is.null(obj.lvmboost)) stop("argument data is NULL")
  if (is.null(obj.lvmboost))
  {
    item_names<-colnames(data)
    nitems<-ncol(data)
    bifreq<-compute_bifreq(data)
    # starting values
    if (kstop<3) message("kstop should be at least 3.")
    D<-1
    start<-rep(NA,nitems)
    for (i in 1:nitems) start[i]<- -qnorm(mean(data[,i]))
    par0<-rep(0,nitems*2)
    par0[seq(2,nitems*2,by=2)]<-start
    
    loglik<-list()
    parall<-list()
    count<-1
    par0mat<-matrix(par0,ncol=nitems)
    colnames(par0mat)<-item_names
    if (trace)
    {
      cat("iteration number: ",count,"\n")
      cat("Parameters: \n")
      print(round(par0mat,2))
    }
    parall[[count]]<-par0mat
    ll<-lik_fa2Rcpp(par=par0,bifreq=bifreq,nitems=nitems,D=D)
    loglik[[count]]<-ll
    
    # =======================================================
    # start with all slopes of first dimension equal to zero
    # and determine which ones move from zero
    # =======================================================
    der<-der_lik_fa2Rcpp(par=par0,bifreq=bifreq,nitems=nitems,D=D,eta=eta)
    gr<-as.vector(der$grad)
    hess<-der$hess
    
    count<-count+1
    
    which_par<-seq(1,nitems*(D+1),by=D+1)
    
    dnc<-DirNegCurvRcpp(which_par=which_par,gr=gr,hess=hess)
    i_par<-dnc$i_par
    j_par<-dnc$j_par
    eigenvector<-dnc$eigenvector
    par1<-par0
    par1[c(i_par,j_par)]<-par0[c(i_par,j_par)]+eigenvector*nu
    par0<-par1
    par0mat<-matrix(par0,ncol=nitems)
    parall[[count]]<-par0mat
    ll<-lik_fa2Rcpp(par=par0,bifreq=bifreq,nitems=nitems,D=D)
    loglik[[count]]<-ll
    

    if (trace)
    {
      cat("\n") ; cat("iteration number: ",count,"\n")
      cat("Parameters: \n")
      print(round(par0mat,2))
    }
    count<-count+1
  }
  else # continuing iterations from a lvmboost object
  {
    item_names<-colnames(obj.lvmboost$parfin)
    nitems<-ncol(obj.lvmboost$parfin)
    bifreq<-obj.lvmboost$bifreq
    data<-obj.lvmboost$data
    nu<-obj.lvmboost$nu
    tau<-obj.lvmboost$tau
    Dmax<-obj.lvmboost$Dmax
    eta<-obj.lvmboost$eta
    loglik<-obj.lvmboost$loglik
    parall<-obj.lvmboost$par
    count<-length(loglik)
    par0mat<-parall[[count]]
    par0<-as.vector(par0mat)
    D<-nrow(par0mat)-1
    kstop<-kstop+count
    count<-count+1
  }
  freepar<-nitems*(nitems-1)/2
  Daug<-TRUE
  noconv<-TRUE
  while(noconv)
  {
    if (D==Dmax) Daug<-FALSE
    # ===========================================================================
    # check adding another dimension
    # ===========================================================================
    if (Daug)
    {
      par0mat<-rbind(0,par0mat)
      par0<-as.vector(par0mat)
    }

    der<-der_lik_fa2Rcpp(par=par0,bifreq=bifreq,nitems=nitems,D=D+Daug,eta=eta)
    gr<-as.vector(der$grad)
    hess<-der$hess
    
    # Newton direction
    nd <- DirNetwRcpp(gr,hess,escl=0)
    nwt_dir<-nd$direction
    val2<-nd$Delta_fz_min # decrese of objective function using the best Newton direction
    i_par_nd<-nd$i_par
    j_par_nd<-nd$j_par
    
    # negative curvature direction
    which_par<-1:length(par0) # all parameters excluding intercepts
    interc<-seq(D+1+Daug,length(par0),by=D+1+Daug) # these are intercepts
    which_par<-which_par[-interc]
    dnc<-DirNegCurvRcpp(which_par=which_par,gr=gr,hess=hess)
    val1<-dnc$Delta_fz_min # decrese of objective function using the best negative curvature direction
    i_par<-dnc$i_par
    j_par<-dnc$j_par
    if (trace) 
      {cat("\n") ; cat("iteration number: ",count,"\n")}

    test<-which.min(c(val1*tau,val2*2))
    if (test==1) # negative curve direction is better
    {
      if (trace) cat("Negative curvature direction chosen \n")
      eigenvector<-dnc$eigenvector
      par1<-par0
      par1[c(i_par,j_par)]<-par0[c(i_par,j_par)]+eigenvector*nu
      Heywood<-TRUE
      frac<-1
      cnt<-1
      while (Heywood)
      {
        frac<-frac*0.1
        par1mat<-matrix(par1,ncol=nitems)
        Lambda<-par1mat[-nrow(par1mat),]
        vy<-t(Lambda)%*%Lambda
        if (any(diag(vy)>1)) 
        {
          par1[c(i_par,j_par)]<-par0[c(i_par,j_par)]+eigenvector*nu*frac
        }
        else Heywood<-FALSE
        cnt<-cnt+1
        if (cnt==4) 
        {
          test<-2
          Heywood<-FALSE
        }
      }
      if (test!=2)
      {
        par0<-par1
        par0mat<-matrix(par0,ncol=nitems)
        # if slopes of new dimension are all zero => delete
        if (all(par0mat[1,]==0)) par0mat<-par0mat[-1,]
        par0<-as.vector(par0mat)
        # if the parameters of the new dimension are linear combination of
        # the parameters of existing dimension => collapse
        if (nrow(par0mat)>2)
        {
          par0mat3<-round(t(par0mat),2)
          par0mat3<-par0mat3[,-ncol(par0mat3)] #remove intercepts
          options(warn=-1)
          correl<-cor(par0mat3)
          options(warn=0)
          correl[is.na(correl)]<-0
          correl<-round(correl,3)
          correl[!lower.tri(correl)]<-999
          if (any(abs(correl)==1)) # if loadings of different dimensions are perfectly correlated, they are are not identifiable => collapse
          {
            lindep<-which(abs(correl)==1,arr.ind = TRUE)
            lindep<-sort(lindep)
            dim1<-lindep[1]
            dim2<-lindep[2]
            y<-par0mat[dim1,]
            x<-par0mat[dim2,]
            alpha<-sum(x*y)/sum(x^2)
            par0mat[dim2,]<-par0mat[dim2,]*sqrt(1+alpha^2)
            par0mat<-par0mat[-dim1,]
          }
          par0<-as.vector(par0mat)
        }
        D<-nrow(par0mat)-1
        llmin1<-lik_fa2Rcpp(par=par0,bifreq=bifreq,nitems=nitems,D=D)
        loglik[[count]]<-llmin1
        parall[[count]]<-par0mat
        
        count<-count+1
      }
    }
    
    if (test==2) # Newton step is better
    {
      if (trace) cat("Newton direction chosen \n")
      cnt1<-1
      escl<-c()
      Heywood1 <- TRUE
      while (Heywood1)
      {
        nwt_dir<-nd$direction
        i_par_nd<-nd$i_par
        j_par_nd<-nd$j_par
        par1<-par0
        par1[c(i_par_nd,j_par_nd)]<-par0[c(i_par_nd,j_par_nd)]+nwt_dir*nu
        Heywood<-TRUE
        frac<-1
        cnt<-1
        while (Heywood)
        {
          frac<-frac*0.1
          par1mat<-matrix(par1,ncol=nitems)
          Lambda<-par1mat[-nrow(par1mat),]
          vy<-t(Lambda)%*%Lambda
          if (any(diag(vy)>1)) 
          {
            par1[c(i_par_nd,j_par_nd)]<-par0[c(i_par_nd,j_par_nd)]+nwt_dir*nu*frac
          }
          else 
          {
            Heywood<-FALSE
            Heywood1<-FALSE
          }
          cnt<-cnt+1
          if (cnt==4) 
          {
            escl<-c(escl,c(i_par_nd,j_par_nd)-1)
            nd <- DirNetwRcpp(gr,hess,escl=escl) #take another Newton direction
            Heywood<-FALSE
          }
          
        }
        cnt1<-cnt1+1
        if (cnt1==4)
        {
          Heywood1<-FALSE
          warnings("Encountered Heywood case")
        }
      }
      par0<-par1
      par0mat<-matrix(par0,ncol=nitems)
      # if slopes of new dimension are all zero => delete
      if (all(par0mat[1,]==0)) par0mat<-par0mat[-1,]
      par0<-as.vector(par0mat)
      if (nrow(par0mat)>2)
      {
        par0mat3<-round(t(par0mat),2)
        par0mat3<-par0mat3[,-ncol(par0mat3)] #remove intercepts
        options(warn=-1)
        correl<-cor(par0mat3)
        options(warn=0)
        correl[is.na(correl)]<-0
        correl<-round(correl,3)
        correl[!lower.tri(correl)]<-999
        if (any(abs(correl)==1)) # if loadings of different dimensions are perfectly correlated, they are are not identifiable => collapse
        {
          lindep<-which(abs(correl)==1,arr.ind = TRUE)
          lindep<-sort(lindep)
          dim1<-lindep[1]
          dim2<-lindep[2]
          y<-par0mat[dim1,]
          x<-par0mat[dim2,]
          alpha<-sum(x*y)/sum(x^2)
          par0mat[dim2,]<-par0mat[dim2,]*sqrt(1+alpha^2)
          par0mat<-par0mat[-dim1,]
        }
        par0<-as.vector(par0mat)
      }
      D<-nrow(par0mat)-1
      parall[[count]]<-par0mat
      ll<-lik_fa2Rcpp(par=par0,bifreq=bifreq,nitems=nitems,D=D)
      
      loglik[[count]]<-ll
      count<-count+1
    }
    if (trace)
    {
      cat("Parameters: \n")
      print(round(par0mat,2))
    }
    if (count==(kstop+1))
    {
      message("max iterations reached \n")
      noconv<-FALSE
    }
    if(mean(abs(gr))<0.00001) {noconv<-FALSE; message("converged at iteration ",count,"\n")}
    if (nfreepar)
    {
      nest<-sum(par0mat[1:D,]!=0)
      if(nest>freepar) 
      {
        noconv<-FALSE
        warning("Algorithm stopped because there are more estimated parameters than free parameters. Consider reducing Dmax.")
      }
    }
  }
  colnames(par0mat)<-item_names
  out<-list(parfin=par0mat,par=parall,loglik=unlist(loglik),D=D,nu=nu,tau=tau,Dmax=Dmax,
            eta=eta,bifreq=bifreq,data=data)
  class(out)<-"lvmboost"
  return(out)
}



cv.lvmboost<-function(obj.lvmboost, K=5, trace=TRUE, seed=1)
{
  data<-obj.lvmboost$data
  n <- nrow(data)
  # the following for excluding subsets with columns all 0 or all 1
  all01 <- TRUE
  set.seed(seed)
  while (all01) {
    all01_tmp <- FALSE
    gr <- split(sample(n, n, replace = FALSE), as.factor(1:K)) # generation of subsets for CROSS VALIDATION
    for (k in 1:K)
    {
      data_k <- data[-gr[[k]],]
      if (any(colMeans(data_k, na.rm = TRUE) == 0 | colMeans(data_k, na.rm = TRUE) == 1)) all01_tmp <- TRUE
    }
    all01 <- all01_tmp
  }
  nitems <- ncol(data)
  kstop <- length(obj.lvmboost$loglik)
  ll<-matrix(NA,kstop,K)
  if (trace) message("folds: ", appendLF = FALSE)
  for (k in 1:K) # k = subset
  {
    if (trace) message(k, appendLF = FALSE)
    training_set <- data[-gr[[k]],]
    validation_set <- data[gr[[k]],]
    bifreq_validation<-compute_bifreq(validation_set)
    suppressMessages(
    res_k<-lvmboost(training_set,kstop=kstop,nu=obj.lvmboost$nu,Dmax=obj.lvmboost$Dmax,
                    tau=obj.lvmboost$tau,eta=obj.lvmboost$eta,
                    trace=FALSE,nfreepar=FALSE)
    )
    
    for (i in 1:length(res_k$par))
    {
      parmat<-res_k$par[[i]]
      D<-nrow(parmat)-1
      # cross-validation error
      ll[i,k]<-lik_fa2Rcpp(par=as.vector(parmat),bifreq=bifreq_validation,nitems=nitems,D=D)
    }
    if (length(res_k$par)<kstop) # if converged before kstop
    {
      ll[(length(res_k$par)+1):kstop,k]<-ll[length(res_k$par),k]
    }
  }
  if (trace) message("\n", appendLF = FALSE)
  cverr<-rowMeans(ll)
  sel<-which.min(cverr)
  par_sel<-obj.lvmboost$par[[sel]]
  colnames(par_sel)<-colnames(data)
  list(cv_err=cverr,sel=sel,par_sel=par_sel)
}




fillzeros<-function(mat,nr)
{
  if(nrow(mat)<nr)
  {
    tmp<-matrix(0,nr,ncol(mat))
    tmp[(nr-nrow(mat)+1):nr,]<-mat
    mat<-tmp
  }
  mat
}


plot.lvmboost<-function(x,which=0,...)
{
  par(mfrow=c(1,2))
  nr<-nrow(x$parfin)
  nc<-ncol(x$parfin)
  par<-sapply(x$par,fillzeros,nr=nr)
  plot(par[1,],type="n",ylim=c(min(par),max(par)),xlab="k",ylab="estimates")
  colors<-rep(c(rep(1,(nr-1)),"grey"),nc)
  if (which==0) which<-1:(nc*nr)
  for (i in which) lines(par[i,],col=colors[i])
  plot(x$loglik,type="l",xlab="k",ylab="negative log-likelihood")
  par(mfrow=c(1,1))
}
