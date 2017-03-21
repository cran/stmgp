
#require(MASS)

stmgplink = function(trainbed,testbed=NULL,gamma=1,taun=NULL,lambda=1,Z=NULL,Zte=NULL,plink="plink --noweb",maf=0.01,hwe=1e-4,geno=0.1,fout="stp",trainfout="train",testfout="test",ll=50,maxal=NULL,alc=NULL){
	return( .STplink(trainbed=trainbed,testbed=testbed,gamma=gamma,lam=taun,ko=lambda,Z=Z,Zte=Zte,plink=plink,maf=maf,hwe=hwe,geno=geno,fout=fout,trainfout=trainfout,testfout=testfout,ll=ll,maxal=maxal,alc=alc,naive=F, omit.tail=F,Sgn=F,probit=F) )
}

.STplink = function(trainbed,testbed=NULL,gamma=1,lam=NULL,ko=1,Z=NULL,Zte=NULL,plink="plink --noweb",maf=0.01,hwe=1e-4,geno=0.1,fout="stp",trainfout="train",testfout="test",ll=50,maxal=NULL,alc=NULL,naive=F, omit.tail=F, Sgn=F, probit=F){
  #trainbed="train"; gamma=1;lam=c(0.5);maxal=NULL;ll=50;Z=NULL;maf=0;r2=0;naive=T; ko=1; omit.tail=F; Sgn=F
  tdir = tempdir()
  fout = paste(tdir,fout,sep="/")
  trainfout = paste(tdir,trainfout,sep="/")
  testfout = paste(tdir,testfout,sep="/")
  nsnp = length(readLines(paste(trainbed,".bim",sep="")))
  FAM = read.table(paste(trainbed,".fam",sep=""))
  y = as.numeric(FAM[,6])
  n = length(y)
  nco = sum(y==1); nca = sum(y==2)
  qua = ifelse(nco+nca==n,"","q")
  if(qua=="") y = y-1
  if(is.null(lam)){
    Tau = n/sqrt(log(n))
  }else{
    Tau = n*lam
  }
  kk = ko
  alo = ifelse(qua=="",3,3)*(n/log(n))/nsnp
  if(is.null(maxal)) maxal = ifelse(qua=="",1,3)*alo
  
  if(is.null(Z)){
    system(paste(plink," --allow-no-sex --bfile ",trainbed," --assoc --pfilter ",maxal," --maf ",maf," --hwe ",hwe," --geno ",geno," --out ",fout,sep=""))
    system(paste(plink," --recodeA --bfile ",trainbed," --extract ",fout,".",qua,"assoc --out ",fout,sep=""))
    system(paste(plink," --bfile ",trainbed," --extract ",fout,".",qua,"assoc --freq --allow-no-sex --out ",fout,sep=""))
    
  }else{
    write.table(cbind(FAM[,1:2],Z),file=paste(fout,".cov",sep=""),row.names=F,col.names=F,quote=F)
    if(qua==""){
      system(paste(plink," --allow-no-sex --bfile ",trainbed," --logistic --covar ",fout,".cov --hide-covar --pfilter ",maxal," --maf ",maf," --hwe ",hwe," --geno ",geno," --out ",fout,sep=""))
      system(paste(plink," --recodeA --bfile ",trainbed," --extract ",fout,".assoc.logistic --out ",fout,sep=""))
      system(paste(plink," --bfile ",trainbed," --extract ",fout,".assoc.logistic --freq --allow-no-sex --out ",fout,sep=""))
    }else{
      system(paste(plink," --allow-no-sex --bfile ",trainbed," --linear --covar ",fout,".cov --hide-covar --pfilter ",maxal," --maf ",maf," --hwe ",hwe," --geno ",geno," --out ",fout,sep=""))
      system(paste(plink," --recodeA --bfile ",trainbed," --extract ",fout,".assoc.linear --out ",fout,sep=""))
      system(paste(plink," --bfile ",trainbed," --extract ",fout,".assoc.linear --freq --allow-no-sex --out ",fout,sep=""))
    }
  }
  
  Qr = read.table(paste(fout,".frq",sep=""),colClasses="character",header=TRUE)[,2:3]
  write.table( Qr, file=paste(fout,".list",sep=""),row.names=F,col.names=F,quote=F )  # For --recode-allele
  D = read.table(paste(fout,".raw",sep=""),header=T)
  X = D[,-(1:6)]
  X = (X==1) + 2*(X==2)
  X = apply(X,2, function(x)ifelse(is.na(x), mean(x, na.rm=TRUE), x))
  ncX = scan(paste(fout,".raw",sep=""),"character",nlines=1)
  sncX = strsplit(ncX[-(1:6)],"_")
  ScX = cbind(Qr[,1],sapply(sncX,tail,1))
  
  ST = .stee(y=y,X=X,Z=Z,Tau=Tau,qua=qua,maxal=maxal,gamma=gamma,ll=ll,naive=naive,kk=kk,alo=alo,grr=1,Sgn=Sgn,alc=alc,probit=probit)
  
  CP = ST$Loss + 2*ST$sig2hato*ST$gdf
  if( which.min(CP)==1 && omit.tail ){
    ilim = which( diff(CP)<0 )[1]
    if(ilim<length(ST$al)) CP[ 1:ilim ] = Inf
  }
  #matplot(-log10(ST$al),CP,type="l")
  lopt = which(CP==min(CP),arr.ind=T)[1,]
  Aopt = which(ST$BA[-(1:length(ST$dz)),lopt[1],lopt[2]]!=0)
  Bopt = ST$BA[-(1:length(ST$dz)),lopt[1],lopt[2]]
  write.table( cbind(ScX[Aopt,],Bopt[Aopt]), file=paste(fout,".score",sep=""),col.names=F,row.names=F,quote=F )


  system(paste(plink," --out ",fout," --bfile ",trainbed," --score ",fout,".score",sep=""))
  system(paste(plink," --out ",trainfout," --bfile ",trainbed," --recodeA --recode-allele ",fout,".list --extract ",fout,".frq",sep=""))
  if(!is.null(testbed)) system(paste(plink," --out ",testfout," --bfile ",testbed," --recodeA --recode-allele ",fout,".list --extract ",fout,".frq",sep=""))

  XG0 = read.table(paste(trainfout,".raw",sep=""),header=TRUE) 
  XG = XG0[,-(1:6)]; XG = (XG==1)+(XG==2)*2
  XG = apply(XG,2, function(x)ifelse(is.na(x), mean(x, na.rm=TRUE), x))

  # Estimated predicted values at optimal al for training data
  if(!is.null(Z)){  # Covariates
    mustp = cbind(Z,XG)%*%ST$BA[-1,lopt[1],lopt[2]] + outer(rep(1,nrow(XG)),ST$BA[1,lopt[1],lopt[2]])
  }else{  # No covariates
    mustp = XG%*%ST$BA[-1,lopt[1],lopt[2]] + outer(rep(1,nrow(XG)),ST$BA[1,lopt[1],lopt[2]])
  }

  if(!is.null(testbed)){
    XG0te = read.table(paste(testfout,".raw",sep=""),header=TRUE)
    XGte = XG0te[,-(1:6)]; XGte = (XGte==1)+(XGte==2)*2
    XGte = apply(XGte,2, function(x)ifelse(is.na(x), mean(x, na.rm=TRUE), x))
    # Estimated predicted values at optimal al for test data
    if(!is.null(Zte)){  # Covariates
      mustpte = cbind(Zte,XGte)%*%ST$BA[-1,lopt[1],lopt[2]] + outer(rep(1,nrow(XGte)),ST$BA[1,lopt[1],lopt[2]])
    }else{  # No covariates
      mustpte = XGte%*%ST$BA[-1,lopt[1],lopt[2]] + outer(rep(1,nrow(XGte)),ST$BA[1,lopt[1],lopt[2]])
    }
  }

  detectSTEE = try(as.character(read.table(paste(fout,".score",sep=""))[,1]),T)  # Non-zero SNPs

  PE = cbind(XG0[,1:6],mustp)
  PEte = NULL
  if(!is.null(testbed)) PEte = cbind(XG0te[,1:6],mustpte)

#  return( list( Muhat=ST$Muhat,gdf=ST$gdf,sig2hat=ST$sig2hat,df=ST$df,al=ST$al,lopt=lopt,BA=ST$BA,Bnaive=ST$Bnaive,RSSnaive=ST$RSSnaive,Munaive=ST$Munaive,Loss=ST$Loss,sig2hato=ST$sig2hato,PE=PE,PEte=PEte,nonzero=detectSTEE ) )

  return( list( Muhat=ST$Muhat,gdf=ST$gdf,sig2hat=ST$sig2hat,df=ST$df,al=ST$al,lopt=lopt,BA=ST$BA,Loss=ST$Loss,sig2hato=ST$sig2hato,CP=CP,PE=PE,PEte=PEte,nonzero=detectSTEE,taun=ST$Tau ) )
}


stmgp = function(y,X,Z=NULL,tau,qb,maxal,gamma=1,ll=50,lambda=1,alc=NULL){
	qua = ifelse(qb=="b","","q")
	ST = .stee(y=y,X=X,Z=Z,Tau=tau,qua=qua,maxal,gamma=gamma,ll=ll,naive=F,kk=lambda,alo=NULL,grr=1,Sgn=F,alc=alc,probit=F)
	CP = ST$Loss + 2*ST$sig2hato*ST$gdf
	lopt = which(CP==min(CP),arr.ind=T)[1,]

	return( list(Muhat=ST$Muhat,gdf=ST$gdf,sig2hat=ST$sig2hat,df=ST$df,al=ST$al,lopt=lopt,BA=ST$BA,Loss=ST$Loss,sig2hato=ST$sig2hato,tau=tau,CP=CP) )
}


.stee = function(y,X,Z,Tau,qua,maxal,gamma=1,ll=50,naive=F,kk=1,alo=NULL,grr=1,Sgn=F,alc=NULL,probit=F){
  n = length(y)
  nT = length(Tau)
  if(is.null(Z)){
    Z = as.matrix( rep(1,n) )
  }else{
    Z = cbind(1,Z)
  }

  groupingXbyCor = function(X,r){
    cX = ( cor(X) >= r )*1
    P1 = cX/colSums(cX)
    return( X%*%P1 )
  }

  X = groupingXbyCor(X,grr)
  alo = ifelse(is.null(alo),maxal,alo)
  nZ = ncol(Z)
  dz = numeric(nZ)
  if(!is.null(alc)) ll = length(alc)
  BA = array(0,c(nZ+ncol(X),ll+1,nT))
  Bnaive = matrix(Inf,nZ+ncol(X),ll+1)
  
  al = c(alo, 10^seq( log10(maxal),log10(5e-8), length=ll) )
  if(!is.null(alc)) al = c(alo,alc)
  
  #print(al)
  
  sig2hat = Loss = matrix(1,ll+1,nT)
  gdf = df = matrix(0,ll+1,nT)
  minpre = Inf
  minal = maxal
  Muhat = array(0,c(n,ll+1,nT))
  Munaive = matrix(0,n,ll+1)
  
  if(qua==""){
    cat("binary phenotype\n")
    muZ = glm(y~Z-1,family="binomial")$fitted.values
    muZd = muZ*(1-muZ)
    GZ = ginv((t(Z)%*%diag(muZd))%*%Z)
    muZdd = muZ*(muZ-1)*(2*muZ-1)
    ZX = (t(Z)%*%diag(muZd))%*%X
    du = X - Z%*%GZ%*%ZX
    dv = Z%*%GZ%*%( (t(Z)%*%diag(muZdd))%*%( du*du ) )
    XX = X;
    
    vv = diag( (t(X)%*%diag(muZd))%*%X - t(ZX)%*%GZ%*%ZX )
    uu = t(X)%*%(y - muZ)
    T2 = as.vector( (uu*uu)/vv ); T2[is.nan(T2)] = 0
    sT2 = sign(uu)
    cr = deT2 = uu
    uuvv = t(matrix( rep( as.vector(uu/vv),n ), length(uu),n ))
    dT2 = 2*du*uuvv - dv*(uuvv*uuvv)
    Th = qchisq(al,df=1,lower.tail=F)
  }else{
    cat("quantitative phenotype\n")
    QZ = diag(1,n,n) - Z%*%ginv(t(Z)%*%Z)%*%t(Z)
    XX = QZ%*%X
    sXX = sqrt(colSums(XX^2))
    XX = XX%*%diag(1/sXX,ncol(X),ncol(X))
    cr = t(XX)%*%y
    sT2 = sign(cr)
    rQy = QZ%*%y
    deT2 = sum(rQy^2) - cr^2
    T2 = (n-nZ-1)*( (cr^2)/deT2 )
    Th = qf(al,1,n-1-nZ,lower.tail=F)
  }
  RSSnaive = rep(Inf,ll+1)
  
  if(qua==""){
    if(probit) qua = "probit"
  }
  
  for(l in 1:(ll+1)){
    Dy = (Th[l]/T2)^((1+gamma)/2)
    A = which(Dy < 1 - 1e-10)
    if(Sgn) A = A[which(sT2[A]>0)]
    nA = length(A)
    
    if(nA>=0){
      dA = nZ + nA
      I = diag(1,dA,dA)
      XA = cbind(Z, X[,A,drop=F] );
      SA = t(XA)%*%XA + kk*I 
      XXA = XX[,A]
      DyA = Dy[A]
      crA = cr[A]
      T2A = T2[A]
      deT2A = deT2[A]
      
      d = c( dz, DyA )
      D = diag(d, dA,dA)
      
      for(it in 1:nT){
        tau = Tau[it]
        bA = as.vector( solve( (I-D)%*%SA + tau*D, (I-D)%*%t(XA)%*%y ) )
        
        if(qua==""){
          rr = 1
          while(rr<20){
            muyA = 1/(1+exp(-as.vector(XA%*%bA)))
            UA = (I-D)%*%(t(XA)%*%(y-muyA) - kk*bA) - tau*D%*%bA
            bA = bA + solve( (I-D)%*%( t(XA)%*%diag(muyA*(1-muyA))%*%XA + kk*I ) + tau*D, UA )
            rr = rr + 1
          }
          BA[c(1:nZ,nZ+A),l,it] = bA
          Muhat[,l,it] = XA%*%bA	
          emu = 1/(1+exp(-as.vector(Muhat[,l,it])))
          if(it==1 && l==1) W0 = diag( emu*(1-emu) )
          mA = -t(X[,A,drop=F])%*%(y-emu) - (tau-kk)*bA[-(1:nZ)]
          if(nA>0){
            sSA = t(matrix( rep( (DyA/T2A)*as.vector(mA),n ), length(mA),n ))
            LA = -((1+gamma)/2)*(dT2[,A]*sSA)
          }else{
            LA = NULL
          }
          LA = cbind( Z*0, LA )
          
          Wa = diag(as.vector(emu*(1-emu)))
          gdf[l,it] = sum(diag( solve( (I-D)%*%( t(XA)%*%Wa%*%XA + kk*I ) + tau*D, (t(LA)%*%W0)%*%XA + (I-D)%*%(t(XA)%*%W0)%*%XA ) ))
          df[l,it] = dA
          Loss[l,it] = -2*(sum(y*Muhat[,l,it] - log(1 + exp(Muhat[,l,it]))))
          
        }else if(qua=="probit"){
          rr = 1
          while(rr<20){
            etaAA = as.vector(XA%*%bA)
            muyA = pnorm(etaAA)
            dmuA = dnorm(etaAA)
            wwA = dmuA/(muyA*(1-muyA))
            UA = (I-D)%*%(t(XA)%*%(wwA*(y-muyA)) - kk*bA) - tau*D%*%bA
            vvsA = as.vector(wwA*( ((1-2*muyA)*wwA + etaAA)*(y-muyA) + dmuA ))
            bA = bA + solve( (I-D)%*%( t(XA)%*%diag(vvsA)%*%XA + kk*I ) + tau*D, UA )
            rr = rr + 1
          }
          BA[c(1:nZ,nZ+A),l,it] = bA
          Muhat[,l,it] = XA%*%bA	
          emu = pnorm( as.vector(Muhat[,l,it]) )
          dmu = dnorm( as.vector(Muhat[,l,it]) )
          wwmu = dmu/( emu*(1-emu) )
          if(it==1 && l==1) W0 = diag( emu*(1-emu) )
          mA = -t(X[,A,drop=F])%*%(wwmu*(y-emu)) - (tau-kk)*bA[-(1:nZ)]
          if(nA>0){
            sSA = t(matrix( rep( (DyA/T2A)*as.vector(mA),n ), length(mA),n ))
            LA = -((1+gamma)/2)*(dT2[,A]*sSA)
          }else{
            LA = NULL
          }
          LA = cbind( Z*0, LA )
          
          Wa = diag(as.vector(wwmu*( ((1-2*emu)*wwmu + Muhat[,l,it])*(y-emu) + dmu )))
          W00 = W0; diag(W00) = diag(W0)*wwmu
          gdf[l,it] = sum(diag( solve( (I-D)%*%( t(XA)%*%Wa%*%XA + kk*I ) + tau*D, (t(LA)%*%W00)%*%XA + (I-D)%*%(t(XA)%*%W00)%*%XA ) ))
          df[l,it] = dA
          Loss[l,it] = -2*sum( y*log(emu) + (1-y)*log(1-emu) )
        }else{
          BA[c(1:nZ,nZ+A),l,it] = bA
          Muhat[,l,it] = XA%*%bA
          if(nA>0){
            mA = SA[-(1:nZ),]%*%bA - (tau-kk)*bA[-(1:nZ)] - t(X[,A,drop=F])%*%y
            ssA = -((1+gamma)/2)*2*(n-1-nZ)*(DyA/T2A)*mA
            cdA = crA/deT2A
            zeA = ((crA/deT2A)^2)*ssA
            CA = diag( as.vector( (1/deT2A + (crA/deT2A)^2)*crA*ssA ), nA,nA )
            LA = XXA%*%CA - rQy%*%t(zeA)
          }else{
            LA = NULL
          }
          LA = cbind( Z*0, LA )
          
          gdf[l,it] = sum(diag( solve( (I-D)%*%SA + tau*D,  t(LA)%*%XA + (I-D)%*%SA ) ))
          df[l,it] = dA
          Loss[l,it] = sum((y - Muhat[,l,it])^2)
          sig2hat[l,it] = sum((y - Muhat[,l,it])^2)/(n-gdf[l,it])
        }
      }
      if(naive==TRUE){
        if(qua==""){
          bnaive = try(glm(y~XA-1,family="binomial")$coef,T)
          if(!is.character(bnaive)){
            Bnaive[,l] = 0
            Bnaive[c(1:nZ,nZ+A),l] = bnaive
            Munaive[,l] = XA%*%bnaive
            RSSnaive[l] = -sum(y*Munaive[,l] - log(1 + exp(Munaive[,l])))
          }
        }else{
          bnaive = try(lm(y~XA-1)$coef,T)
          if(!is.character(bnaive)){
            Bnaive[,l] = 0
            Bnaive[c(1:nZ,nZ+A),l] = bnaive
            Munaive[,l] = XA%*%bnaive
            RSSnaive[l] = sum((y - Munaive[,l])^2)
          }
        }
      }
    }
  }
  
  return(list(Loss=Loss[-1,,drop=F],Muhat=Muhat[,-1,,drop=F],gdf=gdf[-1,,drop=F],sig2hat=sig2hat[-1,,drop=F],df=df[-1,,drop=F],al=al[-1],BA=BA[,-1,,drop=F],Bnaive=Bnaive[,-1,drop=F],RSSnaive=RSSnaive[-1],Munaive=Munaive[,-1,drop=F],dz=dz,sig2hato=sig2hat[1,1],au=Tau))
  
}


## Test
if(0){
  
  D = read.table("stp.raw",header=T)
  X = D[,-(1:6)]
  X = (X==1) + 2*(X==2)
  p = ncol(X)
  n = nrow(X)
  R = 20; ll = 30
  p0 = 20; b0 = log(rep(1.2,p0))
  Op = Oe = GDF = matrix(0,R,ll)
  qua = ""
  probit = F
  sig = 1.2
  
  for(r in 1:R){
    iA0 = sample(1:p,p0)
    Z = as.matrix(cbind(rnorm(n),runif(n)))
    eta = X[,iA0]%*%b0 - 4 + Z%*%c(0.5,0.5)
    if(probit){
      eta = eta - mean(eta)
      mu = pnorm(eta*3-2.5)
      y = rbinom(length(mu),size=1,prob=mu)
      z = Z
      x = X
      N = n
      ica = which(y==1)
      ico = which(y==0)[1:length(ica)]
      y = y[c(ica,ico)]
      mu = mu[c(ica,ico)]
      z = Z[c(ica,ico),]
      x = X[c(ica,ico),]
      N = length(y)
    }else{
      if(qua==""){
        mu = 1/(1+exp(-eta))
        #			mu = exp(eta); mu[which(mu>0.9)] = 0.9
        y = rbinom(n,size=1,prob=mu)
        z = Z
        x = X
        N = n
        #			ica = which(y==1)
        #			ico = which(y==0)[1:length(ica)]
        #			y = y[c(ica,ico)]
        #			mu = mu[c(ica,ico)]
        #			z = Z[c(ica,ico),]
        #			x = X[c(ica,ico),]
        #			N = length(y)
      }else{
        mu = eta
        y = mu + rnorm(n)*sig
        z = Z
        x = X
        N = n
      }
    }
    ST = stee(y,x,z,Tau=N*c(1),qua,maxal=0.1,gamma=1,ll=ll,naive=F,probit=probit)
    MU = outer(as.vector(mu),rep(1,ll))
    Y = outer(as.vector(y),rep(1,ll))
    if(probit){
      Emuhat = pnorm(ST$Muhat[,,1])
      Op[r,] = -2*colSums( MU*log(Emuhat) + (1-MU)*log(1-Emuhat) )
      Oe[r,] = -2*colSums( Y*log(Emuhat) + (1-Y)*log(1-Emuhat) )
      GDF[r,] = ST$gdf
      next
    }
    if(qua==""){
      Op[r,] = -2*colSums( MU*ST$Muhat[,,1] - log(1 + exp(ST$Muhat[,,1])) )
      Oe[r,] = -2*colSums( Y*ST$Muhat[,,1] - log(1 + exp(ST$Muhat[,,1])) )
      GDF[r,] = ST$gdf
    }else{
      Op[r,] = colSums( (MU - ST$Muhat[,,1])^2 ) + N*(sig^2)
      Oe[r,] = colSums( (Y - ST$Muhat[,,1])^2  )
      GDF[r,] = (sig^2)*ST$gdf
    }
    
  }
  
  
  matplot(cbind(colMeans(Op),colMeans(Oe),colMeans(Oe+2*GDF)),type="l")
  
  #matplot(t(Op),type="l")
  #matplot(t(Oe+2*GDF),type="l")
  
  
}


