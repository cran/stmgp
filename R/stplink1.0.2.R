
#require(MASS)

stmgplink = function(trainbed,testbed=NULL,gamma=1,taun=NULL,lambda=1,Z=NULL,Zte=NULL,plink="plink --noweb",maf=0.01,hwe=1e-4,geno=0.1,fout="stp",trainfout="train",testfout="test",ll=50,maxal=NULL,alc=NULL,tdir=NULL,Znames=NULL,trainphenofile=NULL,phenoname=NULL,pSum=NULL){
	return( .STplink(trainbed=trainbed,testbed=testbed,gamma=gamma,lam=taun,ko=lambda,Z=Z,Zte=Zte,plink=plink,maf=maf,hwe=hwe,geno=geno,fout=fout,trainfout=trainfout,testfout=testfout,ll=ll,maxal=maxal,alc=alc,naive=F, omit.tail=F,Sgn=F,tdir=tdir,Znames=Znames,trainphenofile=trainphenofile,phenoname=phenoname,pSum=pSum) )
}


.STplink = function(trainbed,testbed=NULL,gamma=1,lam=NULL,ko=1,Z=NULL,Zte=NULL,plink="plink --noweb",maf=0.01,hwe=1e-4,geno=0.1,fout="stp",trainfout="train",testfout="test",ll=50,maxal=NULL,alc=NULL,naive=F, omit.tail=F, Sgn=F, Nrow=10000, tdir=NULL, Znames=NULL, trainphenofile=NULL, phenoname=NULL, pSum=NULL){
  #trainbed="train"; gamma=1;lam=c(0.5);maxal=NULL;ll=50;Z=NULL;maf=0;r2=0;naive=T; ko=1; omit.tail=F; Sgn=F; tdir=NULL; pSum=NULL


  read.assoc.ch = function(fi,maxal,pSum,nrows=Nrow){
    infile = file(fi, open="r")
    he = read.table(infile, header=FALSE, nrows=1, colClasses="character"); iP = which(he=="P"); iSNP = which(he=="SNP")
    ncols = length(he); temp0 = temp = NULL; temp1 = 1
    while(1){
      temp1 = scan(infile,"character", nlines=nrows)
      if(length(temp1)==0) break
      temp1 = t(matrix(temp1,nrow=ncols))
      rownames(temp1) = temp1[,iSNP]; temp1 = cbind(temp1,temp1[,iP])
      ii = intersect( temp1[,iSNP], rownames(pSum) )
      if(length(ii)>0){
      }
      #temp = rbind(temp,temp1)
      il = which( as.numeric(temp1[,ncol(temp1)]) < maxal )
      if(length(il) > 0) temp0 = rbind(temp0,temp1[il,,drop=F])
    }
    #write.table(temp,file="Aoo.txt")
    colnames(temp0) = c(he,"pFisher")
    close(infile)
    return(temp0)
  }
  #read.assoc.ch(fi="stp.qassoc",maxal=0.001,pSum=s)


  if(is.null(tdir)) tdir = tempdir()
  fout = paste(tdir,fout,sep="/")
  trainfout = paste(tdir,trainfout,sep="/")
  testfout = paste(tdir,testfout,sep="/")
  if(length(trainbed)==1) trainbed = paste0(trainbed,c(".bed",".bim",".fam"))
  nsnp = length(readLines(trainbed[2]))
  FAM = read.table(trainbed[3],na.strings="-9"); rownames(FAM) = as.character(paste(FAM[,1],FAM[,2]))

  if(!is.null(trainphenofile)){
     stopifnot(!is.null(phenoname))
     Phen = read.table(trainphenofile,header=TRUE,na.strings="-9"); rownames(Phen) = as.character(paste(Phen[,1],Phen[,2]))
     cat(paste0("Reading phenotype ",phenoname," from ",trainphenofile,"\n"))
     FAM[,6] = NA; iii = intersect(rownames(FAM),rownames(Phen)); FAM[iii,6] = Phen[iii,phenoname]
  }


  if(!is.null(Z)){
    if(is.character(Z)){  # plink .cov format
      Z = read.table(Z,header=TRUE,na.strings="-9"); rownames(Z) = as.character(paste(Z[,1],Z[,2])); Z = as.matrix(Z[,-(1:2),drop=FALSE]);
      if(!is.null(Znames)){
         cat( paste0( "Used covariates: ",paste(intersect(colnames(Z),Znames),collapse=","),"\n" ) )
         Z = Z[,intersect(colnames(Z),Znames)]
      }
    }
    stopifnot(nrow(Z)==nrow(FAM))
    Z = as.matrix(Z)
    if(is.null(rownames(Z))){
      cat("Assuming rownames of covariate Z are identical to those of .fam file\n")
      rownames(Z) = as.character(paste(FAM[,1],FAM[,2]));
    }
    iZFAM = as.character( intersect(rownames(Z),rownames(FAM)) )  #
    if(length(iZFAM)==0) cat("rownames of covariates (Z) need to be 'FID IID' in .fam file (FID[space]IID)")
    cat(paste("number of individuals in covariate Z overlapping in .fam file is",length(iZFAM),"\n"))
    FAM = FAM[iZFAM,]; Z = Z[iZFAM,,drop=FALSE]
    ia = which(complete.cases(cbind(FAM,Z)))
    FAM = FAM[ia,]; Z = Z[rownames(FAM),,drop=FALSE]
  }else{
    FAM = FAM[which(!is.na(FAM[,6])),]
  }

  write.table(FAM,file=paste0(fout,".keep.fam"),row.names=F,col.names=F,quote=F)

  y = as.numeric(FAM[,6]); names(y) = paste(FAM[,1],FAM[,2])
  n = sum(!is.na(y))
  nco = sum(y==1,na.rm=TRUE); nca = sum(y==2,na.rm=TRUE)
  qua = ifelse(nco+nca==n,"","q")
  if(qua=="") y = y-1
  if(is.null(lam)){
    Tau = n/sqrt(log(n))
  }else{
    Tau = n*lam
  }
  kk = ko

  alo = ifelse(qua=="",3,3)*(n/log(n))/nsnp  # used for sig2hato
  if(is.null(maxal)) maxal = ifelse(qua=="",1,3)*alo

  if(is.null(pSum)){
    if(is.null(Z)){
      system(paste(plink," --allow-no-sex --bed ",trainbed[1]," --bim ",trainbed[2]," --fam ",trainbed[3]," --assoc --pfilter ",maxal," --maf ",maf," --hwe ",hwe," --geno ",geno," --out ",fout," --keep ",fout,".keep.fam",sep=""))
      system(paste(plink," --recodeA --bed ",trainbed[1]," --bim ",trainbed[2]," --fam ",trainbed[3]," --extract ",fout,".",qua,"assoc --out ",fout," --keep ",fout,".keep.fam",sep=""))
      system(paste(plink,"  --bed ",trainbed[1]," --bim ",trainbed[2]," --fam ",trainbed[3]," --extract ",fout,".",qua,"assoc --freq --allow-no-sex --out ",fout," --keep ",fout,".keep.fam",sep=""))

    }else{
      write.table(cbind(FAM[,1:2],Z[as.character(paste(FAM[,1],FAM[,2])),]),file=paste(fout,".cov",sep=""),row.names=F,col.names=F,quote=F)
      if(qua==""){
        system(paste(plink," --allow-no-sex  --bed ",trainbed[1]," --bim ",trainbed[2]," --fam ",trainbed[3]," --logistic --covar ",fout,".cov --hide-covar --pfilter ",maxal," --maf ",maf," --hwe ",hwe," --geno ",geno," --out ",fout," --keep ",fout,".keep.fam",sep=""))
        system(paste(plink," --recodeA --bed ",trainbed[1]," --bim ",trainbed[2]," --fam ",trainbed[3]," --extract ",fout,".assoc.logistic --out ",fout," --keep ",fout,".keep.fam",sep=""))
        system(paste(plink," --bed ",trainbed[1]," --bim ",trainbed[2]," --fam ",trainbed[3]," --extract ",fout,".assoc.logistic --freq --allow-no-sex --out ",fout," --keep ",fout,".keep.fam",sep=""))

      }else{
        system(paste(plink," --allow-no-sex --bed ",trainbed[1]," --bim ",trainbed[2]," --fam ",trainbed[3]," --linear --covar ",fout,".cov --hide-covar --pfilter ",maxal," --maf ",maf," --hwe ",hwe," --geno ",geno," --out ",fout," --keep ",fout,".keep.fam",sep=""))
        system(paste(plink," --recodeA --bed ",trainbed[1]," --bim ",trainbed[2]," --fam ",trainbed[3]," --extract ",fout,".assoc.linear --out ",fout," --keep ",fout,".keep.fam",sep=""))
        system(paste(plink," --bed ",trainbed[1]," --bim ",trainbed[2]," --fam ",trainbed[3]," --extract ",fout,".assoc.linear --freq --allow-no-sex --out ",fout," --keep ",fout,".keep.fam",sep=""))

      }
    }
  }else{

    if( is.null(rownames(pSum)) ) stop( "rownames of pSum must be specified" )

    if(is.null(Z)){
      system(paste(plink," --allow-no-sex --bed ",trainbed[1]," --bim ",trainbed[2]," --fam ",trainbed[3]," --assoc --maf ",maf," --hwe ",hwe," --geno ",geno," --out ",fout," --keep ",fout,".keep.fam",sep=""))
      #Ao = read.table(paste0(fout,".",qua,"assoc"),header=TRUE)  #
      #Ao$pFisher = Ao$P; rownames(Ao) = Ao$SNP
      #ioo = intersect(Ao$SNP,rownames(pSum))
      #Ao[ioo,"pFisher"] = pchisq(-2*log(Ao[ioo,"P"]) - 2*rowSums(log(pSum[ioo,,drop=F]),na.rm=T), df=2*(1+rowSums(!is.na(pSum[ioo,,drop=F]))), lower.tail=F )
      #write.table(Ao,file="Ao.txt")
       write.table(read.assoc.ch(fi=paste0(fout,".",qua,"assoc"),maxal=maxal,pSum=pSum),file=paste0(fout,".",qua,"assoc"),row.names=FALSE,col.names=TRUE,quote=FALSE)  #
      #write.table(Ao[which(-2*log(Ao$P) - 2*rowSums(log(pSum)) > qchisq(maxal,df=2*(1+ncol(pSum)),lower.tail=F)),],file=paste0(fout,".",qua,"assoc"),row.names=F,quote=F)  #
      system(paste(plink," --recodeA --bed ",trainbed[1]," --bim ",trainbed[2]," --fam ",trainbed[3]," --extract ",fout,".",qua,"assoc --out ",fout," --keep ",fout,".keep.fam",sep=""))
      system(paste(plink," --bed ",trainbed[1]," --bim ",trainbed[2]," --fam ",trainbed[3]," --extract ",fout,".",qua,"assoc --freq --allow-no-sex --out ",fout," --keep ",fout,".keep.fam",sep=""))

    }else{
      write.table(cbind(FAM[,1:2],Z[as.character(paste(FAM[,1],FAM[,2])),]),file=paste(fout,".cov",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
      if(qua==""){
        system(paste(plink," --allow-no-sex --bed ",trainbed[1]," --bim ",trainbed[2]," --fam ",trainbed[3]," --logistic --covar ",fout,".cov --hide-covar --maf ",maf," --hwe ",hwe," --geno ",geno," --out ",fout," --keep ",fout,".keep.fam",sep=""))
     write.table(read.assoc.ch(fi=paste0(fout,".assoc.logistic"),maxal=maxal,pSum=pSum),file=paste0(fout,".assoc.logistic"),row.names=FALSE,col.names=TRUE,quote=FALSE)  #
        system(paste(plink," --recodeA --bed ",trainbed[1]," --bim ",trainbed[2]," --fam ",trainbed[3]," --extract ",fout,".assoc.logistic --out ",fout," --keep ",fout,".keep.fam",sep=""))
        system(paste(plink," --bed ",trainbed[1]," --bim ",trainbed[2]," --fam ",trainbed[3]," --extract ",fout,".assoc.logistic --freq --allow-no-sex --out ",fout," --keep ",fout,".keep.fam",sep=""))

      }else{
        system(paste(plink," --allow-no-sex --bed ",trainbed[1]," --bim ",trainbed[2]," --fam ",trainbed[3]," --linear --covar ",fout,".cov --hide-covar --maf ",maf," --hwe ",hwe," --geno ",geno," --out ",fout," --keep ",fout,".keep.fam",sep=""))
        write.table(read.assoc.ch(fi=paste0(fout,".assoc.linear"),maxal=maxal,pSum=pSum),file=paste0(fout,".assoc.linear"),row.names=FALSE,col.names=TRUE,quote=FALSE)  
        system(paste(plink," --recodeA --bed ",trainbed[1]," --bim ",trainbed[2]," --fam ",trainbed[3]," --extract ",fout,".assoc.linear --out ",fout," --keep ",fout,".keep.fam",sep=""))
        system(paste(plink," --bed ",trainbed[1]," --bim ",trainbed[2]," --fam ",trainbed[3]," --extract ",fout,".assoc.linear --freq --allow-no-sex --out ",fout," --keep ",fout,".keep.fam",sep=""))

      }
    }
    
  }

  Qr = read.table(paste(fout,".frq",sep=""),colClasses="character",header=TRUE)[,2:3]
  write.table( Qr, file=paste(fout,".list",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE )  # For --recode-allele
  #D = read.table(paste(fout,".raw",sep=""),header=TRUE)
  ncX = scan(paste(fout,".raw",sep=""),"character",nlines=1)
  D = matrix(scan(paste(fout,".raw",sep=""),"character"),nrow=length(ncX))[,-1]; colnames(D) = paste(D[1,],D[2,])
  X = (t(D[-(1:6),,drop=FALSE])=="1") + 2*(t(D[-(1:6),,drop=FALSE])=="2"); X = X[rownames(FAM),]
  rm(D); gc(); gc();
  X = apply(X,2, function(x)ifelse(is.na(x), mean(x, na.rm=TRUE), x))
  sncX = strsplit(ncX[-(1:6)],"_")
  ScX = cbind(Qr[,1],sapply(sncX,tail,1))
  if(!is.null(pSum)){
    PSum = matrix(NA,ncol(X),ncol(pSum)); rownames(PSum) = Qr$SNP; iQr = intersect(rownames(pSum),rownames(PSum))
    PSum[iQr,] = as.matrix(pSum[iQr,]); pSum = PSum
  }
  ST = .stee(y=y,X=X,Z=Z,Tau=Tau,qua=qua,maxal=maxal,gamma=gamma,ll=ll,naive=naive,kk=kk,alo=alo,grr=1,Sgn=Sgn,alc=alc,pSum=pSum)
  DataTr = list(y=y,X=X,Z=Z)

  CP = ST$Loss + 2*ST$sig2hato*ST$gdf
  if( which.min(CP)==1 && omit.tail ){
    ilim = which( diff(CP)<0 )[1]
    if(ilim<length(ST$al)) CP[ 1:ilim ] = Inf
  }
  #matplot(-log10(ST$al),CP,type="l")
  lopt = which(CP==min(CP),arr.ind=TRUE)[1,]
  Aopt = which(ST$BA[-(1:length(ST$dz)),lopt[1],lopt[2]]!=0)
  Bopt = ST$BA[-(1:length(ST$dz)),lopt[1],lopt[2]]
  allZero = all(Bopt==0)
  names(lopt) = c("opt.al.index","opt.tau.index")


  if(!is.null(testbed) & length(testbed)==1) testbed = paste0(testbed,c(".bed",".bim",".fam"))

  if(!allZero){
    write.table( cbind(ScX[Aopt,],Bopt[Aopt]), file=paste(fout,".score",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE )


    system(paste(plink," --out ",fout," --bed ",trainbed[1]," --bim ",trainbed[2]," --fam ",trainbed[3]," --score ",fout,".score",sep=""))
    system(paste(plink," --out ",trainfout," --bed ",trainbed[1]," --bim ",trainbed[2]," --fam ",trainbed[3]," --recodeA --recode-allele ",fout,".list --extract ",fout,".frq"," --keep ",fout,".keep.fam",sep=""))
    if(!is.null(testbed)){
      system(paste(plink," --out ",testfout," --bed ",testbed[1]," --bim ",testbed[2]," --fam ",testbed[3]," --recodeA --recode-allele ",fout,".list --extract ",fout,".frq",sep=""))
    }

    #XG0 = read.table(paste(trainfout,".raw",sep=""),header=TRUE); XG0[which(XG0[,6]==-9),6] = NA
    #rownames(XG0) = paste(XG0[,1],XG0[,2]); XG0 = XG0[rownames(FAM),]
    ncX0 = scan(paste(trainfout,".raw",sep=""),"character",nlines=1)
    XG0 = t(matrix(scan(paste(trainfout,".raw",sep=""),"character"),nrow=length(ncX0))[,-1]); XG0[which(XG0[,6]=="-9"),6] = NA; rownames(XG0) = paste(XG0[,1],XG0[,2]); XG0 = XG0[rownames(FAM),,drop=FALSE];
    XG = (XG0[,-(1:6),drop=FALSE]=="1") + 2*(XG0[,-(1:6),drop=FALSE]=="2")
    XG = apply(XG,2, function(x)ifelse(is.na(x), mean(x, na.rm=TRUE), x))
  }


  # Estimated predicted values at optimal al for training data
  if(!is.null(Z)){  # Covariates
    if(!allZero){
      mustp = cbind(Z[rownames(FAM),],XG)%*%ST$BA[-1,lopt[1],lopt[2]] + outer(rep(1,nrow(XG)),ST$BA[1,lopt[1],lopt[2]])
    }else{
      mustp = Z%*%ST$BA[1+(1:ncol(Z)),lopt[1],lopt[2]] + outer(rep(1,nrow(Z)),ST$BA[1,lopt[1],lopt[2]])
    }
  }else{  # No covariates
    if(!allZero){
      mustp = XG%*%ST$BA[-1,lopt[1],lopt[2]] + outer(rep(1,nrow(XG)),ST$BA[1,lopt[1],lopt[2]])
    }else{
      mustp = outer(rep(1,nrow(FAM)),ST$BA[1,lopt[1],lopt[2]])
    }
  }

  if(!is.null(testbed)){
    FAMte = read.table(testbed[3]); rownames(FAMte) = paste(FAMte[,1],FAMte[,2])

    if(!allZero){
      #XG0te = read.table(paste(testfout,".raw",sep=""),header=TRUE); XG0te[which(XG0te[,6]==-9),6] = NA; rownames(XG0te) = paste(XG0te[,1],XG0te[,2]); XG0te = XG0te[rownames(FAMte),]
      #XGte = XG0te[,-(1:6)]; XGte = (XGte==1)+(XGte==2)*2;
      #XGte = apply(XGte,2, function(x)ifelse(is.na(x), mean(x, na.rm=TRUE), x))
      ncX0 = scan(paste(testfout,".raw",sep=""),"character",nlines=1)
      XG0te = t(matrix(scan(paste(testfout,".raw",sep=""),"character"),nrow=length(ncX0))[,-1]); XG0te[which(XG0te[,6]=="-9"),6] = NA; rownames(XG0te) = paste(XG0te[,1],XG0te[,2]); XG0te = XG0te[rownames(FAMte),,drop=FALSE]
      XGte = (XG0te[,-(1:6),drop=FALSE]=="1") + 2*(XG0te[,-(1:6),drop=FALSE]=="2")
      XGte = apply(XGte,2, function(x)ifelse(is.na(x), mean(x, na.rm=TRUE), x)); #print(dim(XG0te))
    }

    # Estimated predicted values at optimal al for test data
    if(!is.null(Zte)){  # Covariates
      if(is.character(Zte)){  # plink .cov format
        Zte = read.table(Zte,header=TRUE,na.strings="-9"); rownames(Zte) = as.character(paste(Zte[,1],Zte[,2])); Zte = as.matrix(Zte[,-(1:2)]);
        if(!is.null(Znames)){
           cat( paste0( "Used covariates: ",paste(intersect(colnames(Zte),Znames),collapse=","),"\n" ) )
           Zte = Zte[,intersect(colnames(Zte),Znames)]
        }
      }
      stopifnot(nrow(Zte)==nrow(FAMte))
      if(is.null(rownames(Zte))){
        cat("Assuming rownames of covariate Z are identical to those of .fam file\n")
        rownames(Zte) = as.character(paste(FAMte[,1],FAMte[,2]))
      }

      if(!allZero){
        mustpte = cbind(Zte[rownames(FAMte),],XGte)%*%ST$BA[-1,lopt[1],lopt[2]] + outer(rep(1,nrow(XGte)),ST$BA[1,lopt[1],lopt[2]])
      }else{
        mustpte = Zte[rownames(FAMte),]%*%ST$BA[1+(1:ncol(Zte)),lopt[1],lopt[2]] + outer(rep(1,nrow(Zte)),ST$BA[1,lopt[1],lopt[2]])
      }

    }else{  # No covariates
      if(!allZero){
        mustpte = XGte%*%ST$BA[-1,lopt[1],lopt[2]] + outer(rep(1,nrow(XGte)),ST$BA[1,lopt[1],lopt[2]])
      }else{
        mustpte = outer(rep(1,nrow(FAMte)),ST$BA[1,lopt[1],lopt[2]])
      }
    }
  }

  detectSTEE = try(read.table(paste(fout,".score",sep="")),TRUE)  # Non-zero SNPs

  PE = cbind(FAM,mustp)
  PEte = NULL
  if(!is.null(testbed)) PEte = cbind(FAMte,mustpte)

#  return( list( Muhat=ST$Muhat,gdf=ST$gdf,sig2hat=ST$sig2hat,df=ST$df,al=ST$al,lopt=lopt,BA=ST$BA,Bnaive=ST$Bnaive,RSSnaive=ST$RSSnaive,Munaive=ST$Munaive,Loss=ST$Loss,sig2hato=ST$sig2hato,PE=PE,PEte=PEte,nonzero=detectSTEE ) )

  return( list( Muhat=ST$Muhat,gdf=ST$gdf,sig2hat=ST$sig2hat,df=ST$df,al=ST$al,lopt=lopt,BA=ST$BA,Loss=ST$Loss,sig2hato=ST$sig2hato,CP=CP,PE=PE,PEte=PEte,nonzero=detectSTEE,tau=ST$tau,DataTr=DataTr ) )
}


stmgp = function(y,X,Z=NULL,tau,qb,maxal,gamma=1,ll=50,lambda=1,alc=NULL,pSum=NULL){
	qua = ifelse(qb=="b","","q")
	ST = .stee(y=y,X=X,Z=Z,Tau=tau,qua=qua,maxal,gamma=gamma,ll=ll,naive=FALSE,kk=lambda,alo=NULL,grr=1,Sgn=F,alc=alc,pSum=pSum)
	CP = ST$Loss + 2*ST$sig2hato*ST$gdf
	lopt = which(CP==min(CP),arr.ind=TRUE)[1,]
        names(lopt) = c("opt.al.index","opt.tau.index")
	return( list(Muhat=ST$Muhat,gdf=ST$gdf,sig2hat=ST$sig2hat,df=ST$df,al=ST$al,lopt=lopt,BA=ST$BA,Loss=ST$Loss,sig2hato=ST$sig2hato,tau=ST$tau,CP=CP) )
}


.stee = function(y,X,Z,Tau,qua,maxal,gamma=1,ll=50,naive=FALSE,kk=1,alo=NULL,grr=1,Sgn=FALSE,alc=NULL,pSum=NULL){
  #qb="q";maxal=0.1;gamma=1;ll=50;naive=FALSE;kk=1;alo=NULL;grr=1;Sgn=FALSE;alc=NULL
  n = length(y)
  nT = length(Tau)
  if(is.null(Z)){
    Z = as.matrix( rep(1,n) )
  }else{
    Z = cbind(1,Z)
  }

  if(!is.null(pSum)){
    pSum = as.matrix(pSum)
    stopifnot(nrow(pSum)==ncol(X))
  }

  groupingXbyCor = function(X,r){
    cX = ( abs(cor(X)) >= r )*1
    P1 = cX/colSums(cX)
    return( X%*%P1 )
  }

  if(grr<1) X = groupingXbyCor(X,grr)
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

  dimnames(BA)[[2]] = dimnames(Muhat)[[2]] = dimnames(Loss)[[1]] = dimnames(gdf)[[1]] = dimnames(df)[[1]] = dimnames(sig2hat)[[1]] = al
  dimnames(BA)[[3]] = dimnames(Muhat)[[3]] = dimnames(Loss)[[2]] = dimnames(gdf)[[2]] = dimnames(df)[[2]] = dimnames(sig2hat)[[2]] = Tau

  if(qua==""){
    cat("binary phenotype\n")
    muZ = glm(y~Z-1,family="binomial")$fitted.values
    muZd = muZ*(1-muZ)
    #GZ = ginv((t(Z)%*%diag(muZd))%*%Z)
    GZ = ginv( crossprod(Z,muZd*Z) )
    muZdd = muZ*(muZ-1)*(2*muZ-1)
    #ZX = (t(Z)%*%diag(muZd))%*%X
    ZX = crossprod(Z,muZd*X)
    du = X - Z%*%GZ%*%ZX
    #dv = Z%*%GZ%*%( (t(Z)%*%diag(muZdd))%*%( du*du ) )
    dv = Z%*%GZ%*%( crossprod(Z,muZdd*(du*du)) )
    XX = X;

    #vv = diag( (t(X)%*%diag(muZd))%*%X - t(ZX)%*%GZ%*%ZX )
    vv = diag( crossprod(X,muZd*X) - crossprod(ZX,GZ)%*%ZX )
    uu = crossprod(X,y - muZ)
    T2 = as.vector( (uu*uu)/vv ); T2[is.nan(T2)] = 0
    sT2 = sign(uu)
    cr = deT2 = uu
    uuvv = t(matrix( rep( as.vector(uu/vv),n ), length(uu),n ))
    dT2 = 2*du*uuvv - dv*(uuvv*uuvv)
    if(is.null(pSum)){
      #Th = outer(rep(1,length(T2)),qchisq(al,df=1,lower.tail=FALSE))
      Th = qchisq(al,df=1,lower.tail=FALSE)
    }
    #else{  # Fisher method (abolition)
      #eal = exp( outer( -rowSums(log(pSum),na.rm=TRUE), -0.5*qchisq(al,df=2*rowSums(!is.na(pSum))+2,lower.tail=FALSE) ,"+" ) )
      #eal[eal>1] = 1
      #Th = qchisq(eal,df=1,lower.tail=FALSE)
    #}
  }else{
    cat("quantitative phenotype\n")
#    QZ = diag(1,n,n) - Z%*%ginv(t(Z)%*%Z)%*%t(Z)
#    XX = QZ%*%X
    XX = X - Z%*%ginv(crossprod(Z))%*%(crossprod(Z,X))
    sXX = sqrt(colSums(XX^2))
    #XX = XX%*%diag(1/sXX,ncol(X),ncol(X))
    XX = t(t(XX)/sXX)
    cr = crossprod(XX,y)
    sT2 = sign(cr)
    #rQy = QZ%*%y
    rQy = y - Z%*%ginv(crossprod(Z))%*%(crossprod(Z,y))
    deT2 = sum(rQy^2) - cr^2
    T2 = (n-nZ-1)*( (cr^2)/deT2 )
    if(is.null(pSum)){
      #Th = outer(rep(1,length(T2)),qf(al,1,n-1-nZ,lower.tail=FALSE))
      Th = qf(al,1,n-1-nZ,lower.tail=FALSE)
    }
    #else{  # Fisher method (abolition)
      #eal = exp( outer( -rowSums(log(pSum),na.rm=TRUE), -0.5*qchisq(al,df=2*rowSums(!is.na(pSum))+2,lower.tail=FALSE) ,"+" ) ) #<- mistake
      #eal[eal>1] = 1
      #Th = qf(eal,1,n-1-nZ,lower.tail=FALSE)
    #}
  }
  RSSnaive = rep(Inf,ll+1)

  if(!is.null(pSum)){  # For Fisher method
    slpSum = rowSums(log(pSum),na.rm=TRUE); srpSum = rowSums(!is.na(pSum))  #
  }  #

  for(l in 1:(ll+1)){
    if(is.null(pSum)){
      #Dy = (Th[,l]/T2)^((1+gamma)/2)
      Dy = (Th[l]/T2)^((1+gamma)/2)  #
    }else{  # Fisher method
      eal = exp(-slpSum-0.5*qchisq(al[l],df=2*srpSum+2,lower.tail=FALSE)); eal[eal>1] = 1  #
      if(qua==""){  #
         Dy = (qchisq(eal,df=1,lower.tail=FALSE)/T2)^((1+gamma)/2)  #
      }else{  #
         Dy = (qf(eal,1,n-1-nZ,lower.tail=FALSE)/T2)^((1+gamma)/2)  #
      }  #
    }
    A = which(Dy < 1 - 1e-10)
    if(Sgn) A = A[which(sT2[A]>0)]
    nA = length(A)

    if(nA>=0){
      dA = nZ + nA
#      I = diag(1,dA,dA)
      XA = cbind(Z, X[,A,drop=FALSE] );
      SA = crossprod(XA); diag(SA) = diag(SA) + kk
      XXA = XX[,A]
      DyA = Dy[A]
      crA = cr[A]
      T2A = T2[A]
      deT2A = deT2[A]

      d = c( dz, DyA )
#      D = diag(d, dA,dA)


      for(it in 1:nT){
        tau = Tau[it]
        SAd = (1-d)*SA; diag(SAd) = diag(SAd) + tau*d
        bA = as.vector( solve( SAd, (1-d)*crossprod(XA,y) ) )

        if(qua==""){  # logistic
          rr = 1
          while(rr<20){
            muyA = 1/(1+exp(-as.vector(XA%*%bA)))
            UA = (1-d)*(crossprod(XA,y-muyA) - kk*bA) - tau*d*bA
            #bA = bA + solve( (I-D)%*%( t(XA)%*%diag(muyA*(1-muyA))%*%XA + kk*I ) + tau*D, UA )
            SAd = crossprod(XA, (muyA*(1-muyA))*XA); diag(SAd) = diag(SAd) + kk; SAd = (1-d)*SAd; diag(SAd) = diag(SAd) + tau*d
            bA = bA + solve( SAd, UA )
            rr = rr + 1
          }
          BA[c(1:nZ,nZ+A),l,it] = bA
          Muhat[,l,it] = XA%*%bA
          emu = 1/(1+exp(-as.vector(Muhat[,l,it])))
          #if(it==1 && l==1) W0 = diag( emu*(1-emu) )
          if(it==1 && l==1) W0 = emu*(1-emu)
          mA = -crossprod(X[,A,drop=FALSE],y-emu) - (tau-kk)*bA[-(1:nZ)]
          if(nA>0){
            sSA = t(matrix( rep( (DyA/T2A)*as.vector(mA),n ), length(mA),n ))
            LA = -((1+gamma)/2)*(dT2[,A]*sSA)
          }else{
            LA = NULL
          }
          LA = cbind( Z*0, LA )

          #Wa = diag(as.vector(emu*(1-emu)))
          #gdf[l,it] = sum(diag( solve( (I-D)%*%( t(XA)%*%Wa%*%XA + kk*I ) + tau*D, (t(LA)%*%W0)%*%XA + (I-D)%*%(t(XA)%*%W0)%*%XA ) ))
          Wa = as.vector(emu*(1-emu))
          SAd = crossprod(XA,Wa*XA); diag(SAd) = diag(SAd) + kk; SAd = (1-d)*SAd; diag(SAd) = diag(SAd) + tau*d
          gdf[l,it] = sum(diag( solve( SAd, crossprod(LA,W0*XA) + (1-d)*crossprod(XA,W0*XA) ) ))
          df[l,it] = dA
          Loss[l,it] = -2*(sum(y*Muhat[,l,it] - log(1 + exp(Muhat[,l,it]))))

        }else{  # quantitative
          BA[c(1:nZ,nZ+A),l,it] = bA
          Muhat[,l,it] = XA%*%bA
          if(nA>0){
            mA = SA[-(1:nZ),]%*%bA - (tau-kk)*bA[-(1:nZ)] - crossprod(X[,A,drop=FALSE],y)
            ssA = -((1+gamma)/2)*2*(n-1-nZ)*(DyA/T2A)*mA
            cdA = crA/deT2A
            zeA = ((crA/deT2A)^2)*ssA
            #CA = diag( as.vector( (1/deT2A + (crA/deT2A)^2)*crA*ssA ), nA,nA )
            #LA = XXA%*%CA - rQy%*%t(zeA)
            LA = t(t(XXA)*as.vector( (1/deT2A + (crA/deT2A)^2)*crA*ssA )) - rQy%*%t(zeA)
          }else{
            LA = NULL
          }
          LA = cbind( Z*0, LA )

          SAd = (1-d)*SA; diag(SAd) = diag(SAd) + tau*d
          gdf[l,it] = sum(diag( solve( SAd,  crossprod(LA,XA) + (1-d)*SA ) ))
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

 #print(sig2hat) 
return(list(Loss=Loss[-1,,drop=FALSE],Muhat=Muhat[,-1,,drop=FALSE],gdf=gdf[-1,,drop=FALSE],sig2hat=sig2hat[-1,,drop=FALSE],df=df[-1,,drop=FALSE],al=al[-1],BA=BA[,-1,,drop=FALSE],Bnaive=Bnaive[,-1,drop=FALSE],RSSnaive=RSSnaive[-1],Munaive=Munaive[,-1,drop=FALSE],dz=dz,sig2hato=sig2hat[1,which.min(Tau)],tau=Tau))

}


## Test
if(0){

  D = read.table("stp.raw",header=T)
  X = D[,-(1:6)]
  X = (X==1) + 2*(X==2)
  p = ncol(X)
  n = nrow(X)
  R = 100; ll = 30
  p0 = 10; b0 = log(rep(1.2,p0))
  Op = Oe = GDF = matrix(0,R,ll)
  qua = "q"
  sig = 1.2
  pSum = 1
  crr = numeric(R)

  for(r in 1:R){
    iA0 = sample(1:p,p0)
    if(!is.null(pSum)){
      pSum = runif(p); pSum[iA0] = pchisq(rnorm(length(b0),mean=4*exp(b0),sd=1)^2,df=1,lower.tail=F); pSum = matrix(pSum)
    }
    Z = as.matrix(cbind(rnorm(n),runif(n)))
    eta = X[,iA0]%*%b0 - 4 + Z%*%c(0.1,0.1)/20
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
    ST = .stee(y,x,z,Tau=N*c(1),qua,maxal=0.1,gamma=1,ll=ll,naive=F,pSum=pSum)
    MU = outer(as.vector(mu),rep(1,ll))
    Y = outer(as.vector(y),rep(1,ll))
    if(qua==""){
      Op[r,] = -2*colSums( MU*ST$Muhat[,,1] - log(1 + exp(ST$Muhat[,,1])) )
      Oe[r,] = -2*colSums( Y*ST$Muhat[,,1] - log(1 + exp(ST$Muhat[,,1])) )
      GDF[r,] = ST$gdf
    }else{
      Op[r,] = colSums( (MU - ST$Muhat[,,1])^2 ) + N*(sig^2)
      Oe[r,] = colSums( (Y - ST$Muhat[,,1])^2  )
      GDF[r,] = (sig^2)*ST$gdf
    }

    crr[r] = cor(as.vector(mu),ST$Muhat[,which.min(Oe[r,]+2*GDF[r,]),1])
  }


  mean(crr)

  matplot(cbind(colMeans(Op),colMeans(Oe),colMeans(Oe+2*GDF)),type="l")

  #matplot(t(Op),type="l")
  #matplot(t(Oe+2*GDF),type="l")


}
