gpls <- function(ZX,ZY,WX=rep(1,ncol(ZX)),WY=rep(1,ncol(ZY)),pls.type="cor",comps=2, tol=.Machine$double.eps){
  
  if(pls.type=="cor"){
    
    return( gplscor(ZX,ZY,WX,WY,comps=comps, tol=tol) )
    
  }else if(pls.type=="can"){
    
    return( gplscan(ZX,ZY,WX,WY,comps=comps,tol=tol) )
    
  }else if(pls.type=="reg"){
    
    return( gplsreg(ZX,ZY,WX,WY,comps=comps,tol=tol) )
    
  }else{
    
    return( gplscor(ZX,ZY,WX,WY,comps=comps) )
    
  }
  
}