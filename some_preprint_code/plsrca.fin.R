## plsrca v6...

  ## need to allow for components, and MFA weights.
    ## also need R2s and to trim a lot of this down.
plsrca <- function(X, Y, tol=.Machine$double.eps, ncomps=0, weights="mfa"){
  
  if (nrow(X) != nrow(Y)) {
    stop("X and Y must have the same number of rows.")
  }
  ## not needed unless CA weights (which may not even be necessary...)
  R <- t(X) %*% Y
  Or <- R/sum(R)
  wir <- rowSums(Or)
  wjr <- colSums(Or)
  rm(R)
  rm(Or)
  ## I guess if MFA weights, then Es should be derived frm same weights.
  
  Ox <- X/sum(X)
  wix <- rowSums(Ox)
  wjx <- colSums(Ox)
  
  Oy <- Y/sum(Y)
  wiy <- rowSums(Oy)
  wjy <- colSums(Oy)
  

  Ex <- wix %o% wjx
  Ey <- wiy %o% wjy
  
  X.svd <- gsvd(Ox - Ex,1/wix,1/wjx)
  X.trace <- sum(X.svd$d^2)
  X.rank <- length(X.svd$d)
  rm(X.svd)

  Y.svd <- gsvd(Oy - Ey,1/wiy,1/wjy)
  Y.trace <- sum(Y.svd$d^2)
  Y.rank <- length(Y.svd$d)
  rm(Y.svd)
  
  # print(X.trace)
  # print(Y.trace)
  
  if(weights=="mfa"){
    
    global.row.w <- 1/wjx
    global.col.w <- 1/wjy
    
  }else{
    
    global.row.w <- 1/wir
    global.col.w <- 1/wjr
    
  }
  
    ## do I need the sqrt here?
    #### this is actually a mistake... it only works when everything is *strictly* disjunctive.
  # Zx.orig <- Zx.in <- (Ox - Ex) * sqrt(nrow(X))
  # Zy.orig <- Zy.in <- (Oy - Ey) * sqrt(nrow(Y))
  Zx.orig <- Zx.in <- sweep(Ox - Ex, 1, sqrt(wix), "/")
  Zy.orig <- Zy.in <- sweep(Oy - Ey, 1, sqrt(wiy), "/")

  
  #ncomps <- min(X.rank,Y.rank) ## no chance to change this for now.
  #ncomps <- X.rank ## no chance to change this for now.
  if(ncomps > 0){
    ncomps <- ceiling(min(X.rank,ncomps))
  }else{
    ncomps <- X.rank
  }
  
  LX.matrix <- Xt.matrix <- matrix(0,nrow(X),ncomps)
  LY.matrix <- matrix(0,nrow(Y),ncomps)
  
  pred.u.matrix <- fi.matrix <- fi.a.matrix <- p.matrix<- u.matrix <- matrix(0,ncol(X),ncomps)  
  fj.matrix <- fj.a.matrix <- q.matrix<- v.matrix <- matrix(0,ncol(Y),ncomps)
  
  r2.x.cumulative <- r2.y.cumulative <- d.vec <- beta.vec <- vector("numeric",ncomps)
  y.norm.vec <- x.norm.vec <- vector("numeric",ncomps)
  
  X.hats <- X.recs <- X.resids <- array(0,dim=c(nrow(X),ncol(X),ncomps))
  Y.hats <- Y.recs <- Y.resids <- array(0,dim=c(nrow(Y),ncol(Y),ncomps))
  
  #break.out <- F
  for(i in 1:ncomps){
    #print(i)
    
    # gsvd.res <- gsvd(t(Zx.in) %*% Zy.in,global.row.w,global.col.w)
    ## I need a try/catch here so I can safely exit.
    
      ### this doesn't appear to work...
    # tryCatch(
    #   gsvd.res <- gsvd(t(Zx.in) %*% Zy.in,global.row.w,global.col.w,tol = tol), 
    #   error = function(e){warning("gsvd: matrix below tolerance, breaking out of the loop."); }
    #   )
    # 
    # print(break.out)
    # if(break.out){
    #   break
    # }
    
    gsvd.res <- try( gsvd(t(Zx.in) %*% Zy.in,global.row.w,global.col.w,tol = tol), silent = T )
    if(inherits(gsvd.res,"try-error")){
      warning('GSVD::tolerance.svd() stopped. Results will be returned with fewer components than max(ncol(X),ncol(Y)).')
      break
    }else{
      
      u.matrix[,i] <- gsvd.res$u[,1]
      p.matrix[,i] <- gsvd.res$p[,1]
      fi.matrix[,i] <- gsvd.res$fi[,1]
      fi.a.matrix[,i] <- gsvd.res$fi[,1] / gsvd.res$d[1]
      v.matrix[,i] <- gsvd.res$v[,1]
      q.matrix[,i] <- gsvd.res$q[,1]
      fj.matrix[,i] <- gsvd.res$fj[,1]
      fj.a.matrix[,i] <- gsvd.res$fj[,1] / gsvd.res$d[1]
      d.vec[i] <- gsvd.res$d[1]
      
        ## fix this; it should be unnecessary.
      if(length(gsvd.res$d)!=1){
        LX.matrix[,i] <- (Zx.in %*% gsvd.res$fi %*% diag(1/gsvd.res$d))[,1] #gsvd.res$LX[,1]    
        LY.matrix[,i] <- (Zy.in %*% gsvd.res$fj %*% diag(1/gsvd.res$d))[,1] #gsvd.res$LY[,1]
      }else{
        LX.matrix[,i] <- (Zx.in %*% gsvd.res$fi %*% diag(1/as.matrix(gsvd.res$d)))[,1] #gsvd.res$LX[,1]    
        LY.matrix[,i] <- (Zy.in %*% gsvd.res$fj %*% diag(1/as.matrix(gsvd.res$d)))[,1] #gsvd.res$LY[,1]
      }
      
      ## this prediction step is the biggie to get matched.
      Xt.matrix[,i] <- as.matrix(LX.matrix[,i] / sqrt(sum(LX.matrix[,i]^2)))
      beta.vec[i] <- c(t(Xt.matrix[,i]) %*% (Zy.in %*% (fj.matrix[,i] / d.vec[i]) ))      
      pred.u.matrix[,i] <- t(Xt.matrix[,i]) %*% Zx.in
      
      X.rec <- Xt.matrix[,i] %o% pred.u.matrix[,i]
        X.rec[abs(X.rec) < tol] <- 0
        X.recs[,,i] <- X.rec
        X.hats[,,i] <- (X.recs[,,i] + Ex) * sum(X)
      
      Y.rec <- (Xt.matrix[,i] * beta.vec[i]) %o% q.matrix[,i] ###this is much more intuitive (based on above).
        Y.rec[abs(Y.rec) < tol] <- 0
        Y.recs[,,i] <- Y.rec
        Y.hats[,,i] <- (Y.recs[,,i] + Ey) * sum(Y)
  
    
      Zx.in <- Zx.in-X.recs[,,i] #E
        Zx.in[abs(Zx.in) < tol] <- 0
          ### I think this is wrong...; specifically the "/ sqrt(nrow(X))"
        
        ## quick test
        
        # print( paste0("X trace old = ", sum(sweep(sweep(Zx.in / sqrt(nrow(X)),1,sqrt(1/wix),"*"),2,sqrt(1/wjx),"*")^2) ))
        # print( paste0("X trace new2 = ", sum(sweep(sweep(sweep(Zx.in, 1, sqrt(wix), "*"),1,sqrt(1/wix),"*"),2,sqrt(1/wjx),"*")^2) ))
        
        r2.x.cumulative[i] <- (X.trace - sum(sweep(sweep(sweep(Zx.in, 1, sqrt(wix), "*"),1,sqrt(1/wix),"*"),2,sqrt(1/wjx),"*")^2)) / X.trace
        
      Zy.in <- Zy.in-Y.recs[,,i] #F
        Zy.in[abs(Zy.in) < tol] <- 0
          ### I think this is wrong...; specifically the "/ sqrt(nrow(X))"
        
        
        # print(paste0("Y trace old = ", sum(sweep(sweep(Zy.in / sqrt(nrow(Y)),1,sqrt(1/wiy),"*"),2,sqrt(1/wjy),"*")^2)))
        # print(paste0("Y trace new2 = ", sum(sweep(sweep(sweep(Zy.in, 1, sqrt(wiy), "*"),1,sqrt(1/wiy),"*"),2,sqrt(1/wjy),"*")^2)))
        
        r2.y.cumulative[i] <- (Y.trace - sum(sweep(sweep(sweep(Zy.in, 1, sqrt(wiy), "*"),1,sqrt(1/wiy),"*"),2,sqrt(1/wjy),"*")^2)) / Y.trace
      
      # ## could use these with i+1 instead of $Zx and $Zy
        X.resids[,,i] <- Zx.in
        Y.resids[,,i] <- Zy.in
        
          ## and the R2s here...
      # if(sum(abs(Zx.in))<tol){print('Nothing left in X'); break}#else{ r2.x.cumulative[i] <- (X.trace - sum(gsvd(Zx.in  / sqrt(nrow(X)) ,1/wix,1/wjx)$d^2)) / X.trace }
      # if(sum(abs(Zy.in))<tol){print('Nothing left in Y'); break}#else{ r2.y.cumulative[i] <- (Y.trace - sum(gsvd(Zy.in  / sqrt(nrow(Y)),1/wiy,1/wjy)$d^2)) / Y.trace}
        
        
        if(sum(abs(Zx.in))<tol){break}#else{print(paste0("X trace new = ",sum(gsvd(sweep(Zx.in, 1, sqrt(wix), "*"),1/wix,1/wjx)$d^2)))}#else{ r2.x.cumulative[i] <- (X.trace - sum(gsvd(Zx.in  / sqrt(nrow(X)) ,1/wix,1/wjx)$d^2)) / X.trace }
        if(sum(abs(Zy.in))<tol){break}#else{print(paste0("Y trace new = ",sum(gsvd(sweep(Zy.in, 1, sqrt(wiy), "*"),1/wiy,1/wjy)$d^2)))}#else{ r2.y.cumulative[i] <- (Y.trace - sum(gsvd(Zy.in  / sqrt(nrow(Y)),1/wiy,1/wjy)$d^2)) / Y.trace}
      
      test.r <- t(Zx.in) %*% Zy.in
      test.r[abs(test.r) < tol] <- 0
      if(sum(abs(test.r)) < tol){print('Nothing left in R'); break}
      
      print(i)
    }
    
  }	
  
    ## figure this part out... the sqrt(nrow(X)) shouldn't be necessary. but whatever.
  Y.rec <- Xt.matrix %*% diag(beta.vec) %*% t(q.matrix)
  Y.hat <- (Y.rec + Ey) * sum(Y)
    ### ERROR?!?!
  # Y.resid <- ( (Zy.orig - Y.rec) + Ey) * sum(Y)
  # Y.resid <- ( ((Zy.orig - Y.rec)/sqrt(nrow(X))) + Ey) * sum(Y)
    ## I am dum-dum... maybe?
  Y.resid <- ( ((Zy.orig - Y.rec)/sqrt(1/wiy)) + Ey) * sum(Y)
  
  
  rownames(X.recs) <- rownames(X.hats) <- rownames(X.resids) <- rownames(Y.resids) <- rownames(Y.resid) <- rownames(Y.recs) <- rownames(Y.rec) <- rownames(Y.hat) <- rownames(Y.hats) <- rownames(Xt.matrix) <- rownames(LX.matrix) <- rownames(LY.matrix) <- rownames(X)
  colnames(X.recs) <- colnames(X.hats) <- colnames(X.resids) <- rownames(fi.matrix) <- rownames(fi.a.matrix) <- rownames(pred.u.matrix) <- rownames(u.matrix) <- rownames(p.matrix) <- colnames(X)
  colnames(Y.recs) <- colnames(Y.hats) <- colnames(Y.resids) <- rownames(fj.matrix) <- rownames(fj.a.matrix) <- rownames(v.matrix) <- rownames(q.matrix) <- colnames(Y)
  
  
  
  return(list(
    fi.matrix=fi.matrix,
    fj.matrix=fj.matrix,
    fi.a.matrix=fi.a.matrix,
    fj.a.matrix=fj.a.matrix,
    u.matrix=u.matrix, 
    v.matrix=v.matrix,
    p.matrix=p.matrix,
    q.matrix=q.matrix,
    d.vec=d.vec,
    beta.vec=beta.vec,
    pred.u.matrix=pred.u.matrix,
    Xt.matrix=Xt.matrix,
    lx=LX.matrix,
    ly=LY.matrix,
    X.resids=X.resids,
    Y.resids=Y.resids,
    X.recs=X.recs,
    Y.recs=Y.recs,
    X.hats=X.hats,
    Y.hats=Y.hats,
    Ex=Ex,
    Ey=Ey,
    Zx.orig=Zx.orig,
    Zy.orig=Zy.orig,
    wir=wir,
    wjr=wjr,
    wix=wix,
    wiy=wiy,
    wjx=wjx,
    wjy=wjy,
    global.row.w=global.row.w,
    global.col.w=global.col.w,
    Y.hat=Y.hat,
    Y.rec=Y.rec,
    Y.resid=Y.resid,
    r2.x.cumulative=r2.x.cumulative,
    r2.x = diff(c(0,r2.x.cumulative)),
    r2.y.cumulative=r2.y.cumulative,
    r2.y = diff(c(0,r2.y.cumulative))
  )	)
  
}