##R package for flowPeaks
##Most of them were written in C++, with a couple of R package functions

#remove h0 and h1, use only tol between 0 and 1
#dyn.load("flowPeaks.so")
.ver<-"1.0"
.Kmax<-300
## 2.0 is more robust then version 1.0 but have problems
## in putting neiboring two clusters together

##choices of h of 1.5 or 1.0. 1.5 seems to be too restrict,
##double the volumn, i.e. h=2
#1.0 seems to be too flexible, how to deal with the many comparisons
#of the same data many times, the range test.

replace.levels<-function(x,l.old,l.new)
{
    if(length(l.old)!=length(l.new)){
        stop("The input data for l.old and l.new should have the same size\n")
    }
    res<-.C("Rpack_relevel",as.integer(x),x.new=integer(length(x)),
            as.integer(length(x)),
            as.integer(l.old),as.integer(l.new),as.integer(length(l.old)),
            PACKAGE="flowPeaks")
    invisible(res$x.new)
}

checkxr<-function(x,cl,ch)
{
    if(length(x)!=1 || !is.numeric(x)){
        stop("the setting for ",deparse(substitute(x))," is incorrect, it should be of size 1 numeric vector\n")
    }

    if(x<cl || x>ch){
        stop("the setting for ",deparse(substitute(x))," is incorrect, it should be between ",cl," and ", ch,"\n")
    }
}
adjust.S<-function(S,h,S0,h0,W,K,n)
{
    #S has K columns and p*p rows
    #find the number
    lambda <-K/(W*n+K) ##rep(1,length(W)), to fix
    p<-nrow(S0)
    Snew<-S*h
    S0v<-as.vector(S0)*h0
    for(i in 1:K){
        Snew[,i]<-Snew[,i]*(1-lambda[i])+S0v*lambda[i]
        ##Snew[,i]<-Snew[,i]+S0v*lambda[i]
    }
    invisible(Snew)
}

gen.voronoi<-function(x,y,xybox)
{
    autobox<-1
    argbox<-rep(0,4)
    if(!missing(xybox)){
        argbox<-as.numeric(xybox)
        if(length(argbox)!=4 || sum(is.na(argbox))>0){
            stop("The argument xybox should be of length 4 numerical vector\n");
        }
    }
    n<-length(x)
    if(n<2){
        stop("x should be of vector with length >1\n")
    }
    if(length(y)!=n){
        stop("y should be of the same length as x\n")
    }
    nedgemax<-max(1,3*n-6)
    res<-.C("Rpack_voronoi",as.integer(n),as.double(x),as.double(y),
            as.integer(autobox),as.double(argbox),nedge=integer(1),
            coord=double(nedgemax*4),site=integer(nedgemax*2),
            PACKAGE="flowPeaks")
    nedge<-res$nedge
#    browser()
    cbind(matrix(res$coord[1:(4*nedge)],nedge,4,byrow=TRUE),
          matrix(res$site[1:(2*nedge)],nedge,2,byrow=TRUE))
}

getsummarizecluster<-function(x,clusterid)
{
    x<-data.matrix(x)
    checkx(x)
    n<-nrow(x)
    p<-ncol(x)
    if(nrow(x)!=length(clusterid)){
        stop("The length of clusterid should be the same as the number of rows of x\n")
    }
    K<-max(clusterid) #assuming clusterid is from 1...K
    clusterid<-clusterid-1 #translate to the C++ format
    res<-.C("Rpack_summarize_cluster",as.double(t(x)),as.integer(n),as.integer(p),as.integer(K),
            as.integer(clusterid),nc=integer(K),C=double(K*p),
            S=double(K*p*p),twss=double(1))
    invisible(list(nc=res$nc,C=matrix(res$C,K,p,byrow=TRUE),
                   S=matrix(res$S,p*p,K,byrow=FALSE),twss=res$twss))
}
    
getS0K<-function(x)
{
    p<-ncol(x)
    nc<-sapply(1:p,function(i){nclass.FD(x[,i])})
    K<-ceiling(median(nc))
    K1<-ceiling(nrow(x)^{1/3})
    if(K<K1){
        K<-K1
    }
    if(K>.Kmax){
        K<-.Kmax
    }
    xr<-apply(x,2,max)-apply(x,2,min)
    S0<-diag((xr/K^{1/p})^2,nrow=length(xr))/3 ##assume uniform distribution
    #S0<-diag(apply(x,2,var)) #/(K^{1/p})^2 ##assume with the same marginal variance 
    invisible(list(K=K,S0=S0))
}

unlistS<-function(S)
{
    #change the list into a rectangle matrix
    res<-NULL
    for(i in 1:length(S)){
        res<-cbind(res,S[[i]])
    }
    invisible(res)
}

getpeaks<-function(Cw,Cm,CS,Nb,tol)
{
    K<-nrow(Cm)
    p<-ncol(Cm)
    #be careful to transle to the C conventions for matrix
    res<-.C("Rpack_get_flowpeaks",as.double(Cw),as.double(t(Cm)),as.double(unlistS(CS)),as.integer(K),
            as.integer(p),as.integer(t(Nb)),
            mu=double(K*p),f=double(K), df=double(K*p),
            found=integer(K),cid=integer(K),
            as.double(tol),PACKAGE="flowPeaks")
    invisible(list(mu=matrix(res$mu,K,p, byrow=TRUE),f=res$f,
                   df=matrix(res$df,K,p,byrow=TRUE),
                   found=res$found,cid=res$cid))
}


getpeaks2<-function(Cw,Cm,CS,Nb,tol)
{
    K<-nrow(Cm)
    p<-ncol(Cm)
    #be careful to transle to the C conventions for matrix
    res<-.C("Rpack_get_flowpeaks2",as.double(Cw),as.double(t(Cm)),
            as.double(unlistS(CS)),as.integer(K),
            as.integer(p),as.integer(t(Nb)),
            f=double(K), mu=double(K*p),cid=integer(K),as.double(tol),
            PACKAGE="flowPeaks")
    invisible(list(f=res$f,mu=matrix(res$mu,K,p, byrow=TRUE),
                   cid=res$cid))
}


raster.image<-function(x,y,rawid,nres=400)
{
    raw<-rbind(x,y) #using the C traditions
    n<-length(x)
    res<-.C("Rpack_raster_image",as.double(raw),as.integer(rawid),
            as.integer(n),as.integer(nres),
            grid=double(n*2),grid.id=integer(n),
            ngrid=integer(1),
            PACKAGE="flowPeaks")
    grid<-matrix(res$grid,n,2,byrow=TRUE)
    invisible(cbind(grid[1:res$ngrid,],res$grid.id[1:res$ngrid]))
}
get.kmeans<-function(x,K,stime)
{
    n<-nrow(x)
    p<-ncol(x)
    #future plan
    ##make the transpose in C to reduce the memory footage
    ##so that only one copy of x is stored.
    res<-.C("Rpack_kmeans",
            as.double(t(data.matrix(x))),as.integer(n),as.integer(p),as.integer(K),
            cluster=integer(n),m=double(K*p),nc=integer(K),S=double(K*p*p),
            Nb=integer(K*K),twss=double(1),as.double(stime),
            PACKAGE="flowPeaks") #DUP=FALSE,depreciated
    res<-list(Cw=res$nc/n,Cm=matrix(res$m,K,p,byrow=TRUE),
              CS=matrix(res$S,p*p,K,byrow=FALSE),
              Nb=matrix(res$Nb,K,K,byrow=TRUE),
              cluster=res$cluster+1,twss=res$twss)
    #write(c(K,p),"data.xls",sep="\t",ncol=2)
   # write.table(cbind(res$Cw,res$Cm,t(res$CS)),"data.xls",sep="\t",row.names=FALSE,
   #             col.names=FALSE,append=TRUE)
    invisible(res)
}
    

traditional.kmeans<-function(x,iter.max=20){
    proc0<-proc.time()[3]
    res<-getS0K(x)
    K<-res$K
    S0<-res$S0
    y<-kmeans(x,K,iter.max=iter.max)
    duration<-proc.time()[3]-proc0
    message("Finished kmeans at ",round(duration,digits=3)," sec\n");

    ##obtain information for each groups from the output of the kmeans
    p<-ncol(x)
    Cm<-matrix(0,K,p) #mean
    Cw<-rep(0,K)
    CS<-vector("list",K) #Sigma

    for(i in 1:K){
        Cw[i]<-y$size[i]/nrow(x)
        Cm[i,]<-y$center[i,]
        
        z<-x[y$cluster==i,]
        S<-var(z)
        w.c<-Cw[i]*nrow(x)
        S<-(S*w.c+K*S0)/(w.c+K) 

        CS[[i]]<-S
    }

    ####write the data to an output
    #write(c(K,p),"data.xls",sep="\t",ncol=2)
    res<-NULL
    for(i in 1:length(CS)){
        res<-rbind(res,c(Cw[i],Cm[i,],c(CS[[i]])))
    }
    #write.table(res,"data.xls",sep="\t",row.names=FALSE,
    #            col.names=FALSE,append=TRUE)
}
checkx<-function(x){
    if(!is.matrix(x)){
        stop("The input data for x is not of matrix\n")
    }
    xdim<-dim(x)
    if(xdim[1]<2 || !is.numeric(x) || sum(is.na(x))>0){
        stop("There are problems in reading input data for x: it maybe has only one row, or contains NAs, or contains non-numbers\n")
    }

    xname<-colnames(x)
    if(length(unique(xname))!=ncol(x)){
        stop("The colnames of the data matrix is not unique\n")
    }
}

process.flowPeaks<-function(x,kmeans.info,kmeans.cluster,tol,
                            S0,h0,h,ver)
{
    Cw<-kmeans.info$w
    Cm<-kmeans.info$m
    Nb<-kmeans.info$Nb
    K<-nrow(Cm)
    p<-ncol(Cm)
    CS<-adjust.S(kmeans.info$S,h,S0,h0,Cw,K,nrow(x))*2^{1/p}

    if(ver=="1.0"){
        res<-getpeaks(Cw,Cm,CS,Nb,tol)
    }else{
        res<-getpeaks2(Cw,Cm,CS,Nb,tol)
    }
    
    cid<-res$cid+1
    Cpeaks<-res$mu
    if(ver=="1.0"){
        Cpeaks.fn<-cbind(res$f,res$df)
        Cpeaks.found<-res$found
    }else{
        Cpeaks.fn<-res$f
        for(i in 1:p){
            Cpeaks.fn<-cbind(Cpeaks.fn,NA)
        }
        Cpeaks.found<-rep(NA,nrow(Cpeaks))
    }

    cidw<-tapply(Cw,cid,sum)
    r<-order(cidw,decreasing=TRUE)
    cid.new<-rep(0,length(cid))
    for(i in r){
        cid.new[cid==r[i]]<-i
    }

    y.new<-kmeans.cluster
    y.new<-replace.levels(y.new,1:length(cid.new),cid.new)

    peaks.info<-cbind(cid=cid.new,w=Cw,mu=Cpeaks,fn=Cpeaks.fn,found=Cpeaks.found,kmeans.id=1:K,kmeans.center=Cm)
    id<-order(peaks.info[,1],-peaks.info[,2],peaks.info[,p+3])
    peaks.info<-peaks.info[id,]


    ##summary the Mean and S (nonsmoothed) for each cluster
    peaks<-getsummarizecluster(x,y.new)
    peaks<-list(cid=1:max(cid.new),w=peaks$nc/sum(peaks$nc),
                mu=peaks$C,S=peaks$S)
    fp<-list(peaks.cluster=y.new,
             peaks=peaks,
             kmeans.cluster=kmeans.cluster,
             kmeans=kmeans.info,
             info=peaks.info,
             x=x,
             S0=S0,
             options=list(tol=tol,h=h,h0=h0))
    class(fp)<-"flowPeaks"
    invisible(fp)
}

##furhter manual merging and move the kmeans from one cluster
##to another cluster
##or divide into two clusters
##a two column matrix
## the first is the kmeans id if positive, is the cid if negative
## the second is where it belongs. which can be an exisiting cid: merging the
#kmeans  or the cid to this cid, the old cid should be removed,
#and the kmeans should be removed from the old cid.
#or if the number is negative, it generates a cid

changemember.flowPeaks<-function(fp,z)
{
    if(class(fp)!="flowPeaks"){
        stop("The data input must come from the function flowPeaks\n")
    }

    K<-nrow(fp$info)
    Kp<-max(fp$info[,1])
    cid<-fp$info[,1]
    cidnew<-cid
    for(i in 1:nrow(z)){
        if(z[i,1]<0){
            if(z[i,2]<0){
                stop("the data format is wrong in ",i,"row\n") 
            }
            cidnew[cid==-z[i,1]]<-z[i,2]
        }else{
            cidnew[z[i,1]]<-z[i,2]
        }
    }
    ##reorder the data according to the new weights
    
    
}



adjust.flowPeaks<-function(object,tol,h0,h,...)
{
    if(class(object)!="flowPeaks"){
        stop("The data input for adjust.flowPeaks  must come from the function flowPeaks\n")
    }
    if(missing(tol)){
        tol<-object$options$tol
    }else{
        checktol(tol)
    }
    
    if(missing(h)){
        h<-object$options$h
    }else{
        ##checkxr(h,1.0,10)
    }

    if(missing(h0)){
        h0<-object$options$h0
    }else{
        ##checkxr(h0,1.0,10)
    }
    object<<-process.flowPeaks(object$x,object$kmeans,object$kmeans.cluster,tol,
                               object$S0,h0,h,.ver)
}

checktol<-function(tol)
{
    if(length(tol)!=1 || !is.numeric(tol)){
        stop("the tol setting is incorrect, it should be of size 1 numeric vector\n")
    }
    if(tol<0 || tol>1){
        stop("tol should be between 0 and 1\n")
    }
}

flowPeaks<-function(x,tol=0.1,h0=1,h=1.5)
{
    checktol(tol)
    #checkxr(h,1.5,10)
    #checkxr(h0,1,10)
    ##data checking for x
    ##check if it is a data matrix
    #check if x is coming from data frame
    
    ##if(class(x)==flowFrame){ ###this will introduce the dependency between
                               ###flowCore and flowPeaks
    ##    require(flowCore)
    ##    x<-data.matrix(x@exprs)
    ##    ##remove the Time columns
    ##    x<-x[,toupper(colnames(x))%in% c("TIME","TIMES")]
    ##}else{
        
        x<-data.matrix(x) #x should just be simple matrix

    checkx(x)
    
    message("\nStarting the flow Peaks analysis...\n\n    Task A: compute kmeans...\n");
    proc0<-proc.time()[3]
    n<-nrow(x)
    p<-ncol(x)
    res<-getS0K(x)
    S0<-res$S0
    K<-res$K

    rm("res")
    duration<-proc.time()[3]-proc0
    res<-get.kmeans(x,K,duration)
    
    kmeans.info<-list(kmeans.id=1:K,w=res$Cw,mu=res$Cm,
                      S=res$CS,Nb=res$Nb)
    duration<-proc.time()[3]-proc0
    message("        ...finished summarization at ",round(duration,digits=3),
            " sec\n\n    Task B: find peaks...");
    
    fp<-process.flowPeaks(x,kmeans.info,res$cluster,tol,
                          S0,h0,h,.ver)
    duration<-proc.time()[3]-proc0
    message("finished at ", round(duration,digits=3)," sec\n\n");
    invisible(fp)
}

print.flowPeaks<-function(x,...)
{
    print(summary(x))
}
summary.flowPeaks<-function(object,...)
{
    peaks<-object$peaks
    res<-cbind(peaks$cid,peaks$w,peaks$mu)
    colnames(res)<-c("cluster.id","weight",
                     paste(colnames(object$x),"center",sep="."))
    res
}
plot.flowPeaks<-function(x,idx=c(1,2),drawlab=FALSE,
                         cols=c("red","green3","blue","cyan",
                             "magenta","yellow","gray"),drawvor=TRUE,
                         drawlocalpeaks=FALSE,drawkmeans=TRUE,drawboundary=TRUE,
                         classlab,
                         negcol,negpch,...)
    #white is the background, black is for label and other stuff,
    #cols are recycled only for voronoi boundary
    #needs to take care of the case that p!=px and the order is different
{
    fp<-x
    if(class(fp)!="flowPeaks"){
        stop("The data input for fp must come from the function flowPeaks\n")
    }
    p<-ncol(fp$x)
    if(p==1){
        ##processed differently
        ##idx is ignored
        K<-nrow(fp$info)
        Kc<-max(fp$info[,1])
        cols<-rep(cols,ceiling(Kc/length(cols)))
        hist(fp$x,main="",xlab=colnames(fp$x),nclass=K)
        points(fp$x,rep(0,length(fp$x)),
               col=cols[fp$peaks.cluster],pch=18)
        return(NULL)
    }
    if("black" %in% cols || "white" %in% cols) {
        stop("Please do not use the black and white in the cols specifications\n")}
    if((length(setdiff(idx,1:p))>0) || (length(unique(idx))!=length(idx))){
        stop("The idx assignment is out of arrange or with repeats\n");
    }
    if(ncol(fp$x)==1 || length(idx)<2){
        stop("The idx needs to be at least 2 of length or the data x that is used to generate flowPeaks needs have at least 2 columns\n")
    }
    
    info<-fp$info
    K<-nrow(info)
    Kc<-max(info[,1])
    cols<-rep(cols,ceiling(Kc/length(cols)))
    if(length(cols)<Kc && (p!=2)){
        ##use another setup of the colors
        #jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
        #                                "#7FFF7F", "yellow", "#FF7F00", "red",
        #                                "#7F0000"))
        #cols<-jet.colors(Kc)
        ##change the positions so that they can be more seperated
        #drawlab=TRUE
    }
    if(missing(classlab)){
        classlab<-fp$peaks.cluster
    }
    negid<-classlab<0 |is.na(classlab)
    xlab<-colnames(fp$x)
    for(i in 1:(length(idx)-1)){
        for(j in (i+1):length(idx)){
            xi<-idx[i]
            xj<-idx[j]
            #rimg<-raster.image(fp$x[,xi],fp$x[,xj],fp$peaks.cluster)
            #plot(fp$x[,xi],fp$x[,xj],xlab=xlab[xi],ylab=xlab[xj],pch=".",asp=1)
            plot(range(fp$x[,xi]),range(fp$x[,xj]),xlab=xlab[xi],
                 ylab=xlab[xj],type="n")#,asp=1)
            if(sum(negid)>0){
                if(missing(negcol)){
                    negcol<-"black"
                }
                if(missing(negpch)){
                    negpch<-"."
                }
                points(fp$x[negid,xi],fp$x[negid,xj],pch=negpch,col=negcol)
            }
                
            rimg<-raster.image(fp$x[!negid,xi],fp$x[!negid,xj],classlab[!negid])
            points(rimg[,1],rimg[,2],col=cols[rimg[,3]],pch=18,cex=0.5)
            #pch=18 may cause the file too big, though the most pleasing
            #another choice is *
            ##plot the peaks and cluster
            for(k in 1:K){
                id<-info[k,1]
                if(drawlocalpeaks==TRUE){
                    points(info[k,2+xi],info[k,2+xj],
                           col="black",bg="white",cex=1,pch=24)
                }
                if(drawkmeans==TRUE){
                    points(info[k,2*p+5+xi],info[k,2*p+5+xj],
                           pch=1,cex=0.7,col="black",bg="white")
                }
            }
            #label the cluster by the center
            peaks<-fp$peaks
            if(drawlab==TRUE){
                text(peaks$mu[,xi],peaks$mu[,xj],labels=peaks$cid,cex=1.5,col="black")
            }else{
                points(peaks$mu[,xi],peaks$mu[,xj],col="black",bg="white",
                       pch=10,cex=1.8,lwd=1.2)
            }
            ##plot the voronoi diagram if the dimension is two
            if(p==2){
                xybox<-par()$usr
                ##the round function is required due to a bug in the voroni program.
                vor<-gen.voronoi(round(info[,2*p+5+xi],digits=6),
                                 round(info[,2*p+5+xj],digits=6),xybox)
                #vor<-gen.voronoi(info[,2*p+5+xi],info[,2*p+5+xj],xybox)
                #
                ##first draw the individual edges in light color
                if(drawvor){
                    segments(vor[,1],vor[,2],vor[,3],vor[,4],col="black",lwd=1,lty="longdash")}
                ##find out which edges also belong to different clusters
                if(drawboundary){
                    cid1<-info[vor[,5]+1,1]
                    cid2<-info[vor[,6]+1,1]
                    id<-cid1!=cid2
                    #browser()
                    #write.table(vor,"vor.txt",sep="\t",row.names=FALSE,col.names=FALSE)
                    #write.table(cbind(info[,2*p+5+xi],info[,2*p+5+xj]),"vor_kmeans.txt",row.names=FALSE,col.names=FALSE)
          
                    segments(vor[id,1],vor[id,2],vor[id,3],vor[id,4],
                             col="black",lwd=3)
                }
            }
        }
    }
}
####experimental
assignkmeans<-function(fp,A)
{
    ##good for all assignments
    A<-data.matrix(A)
    p<-as.integer(ncol(A))
    C<-fp$info[,2*p+5+1:p]
    n<-as.integer(nrow(A))
    res<-.C("Rpack_assign_kmeans",as.double(t(A)),n,
            p,as.double(t(C)),as.integer(nrow(C)),
            ret.IC=integer(n),PACKAGE="flowPeaks")
    fp$info[,2*p+5][res$ret.IC+1]
}

assign.flowPeaks<-function(fp,A,tol=0.01,fc=0.8)
{
    ##good for all assignments
    A<-data.matrix(A)
    p<-as.integer(ncol(A))
    Cm<-fp$info[,2*p+5+1:p]
    cid<-fp$info[,1]-1
    Cw<-fp$info[,2]
    CS<-adjust.S(fp$kmeans$S,fp$options$h,
                 fp$S0,fp$options$h0,
                 fp$kmeans$w,
                 nrow(Cm),nrow(fp$x))*2^{1/p}
    CS<-CS[,fp$info[,2*p+5]]
    n<-as.integer(nrow(A))
    res<-.C("Rpack_assign_flowPeaks",as.double(t(A)),n,
            p,as.double(Cw),as.double(t(Cm)),as.double(CS),
            as.integer(nrow(Cm)),
            as.integer(cid),
            as.double(tol),as.double(fc),
            ret.IC=integer(n),PACKAGE="flowPeaks")
    res$ret.IC
}

###potenital exporting terms
#1. assignkmeans, 2. get.kmeans 3 assignflowPeaks
evalCluster<-function(gs,cand,method=c("Rand.index","Fmeasure","Vmeasure"),
                      rm.gs.outliers=TRUE)
{
    measurefunc<-switch(method,Rand.index=randmeasure,
                        Fmeasure=fmeasure,Vmeasure=vmeasure)
    if(is.null(measurefunc)){
        stop("The method should be one of the Rand.index, Fmeasure Vmeasure\n")
    }
    if(length(gs)!=length(cand)){
        stop("the length of gs should be the same as cand\n")
    }
    if(rm.gs.outliers){
        id<- which(gs> -1)
    }else{
        id<-1:length(cand)
    }
    measurefunc(table(gs[id],cand[id]))
}

vmeasure<-function(a,beta=1)
{
    h<-rel.info(a)
    c<-rel.info(t(a))
    (1+beta)*h*c/(beta*h+c)
}
randmeasure<-function(a)
{
    a<-data.matrix(a)
    ac<-apply(a,1,sum)
    ak<-apply(a,2,sum)
    n<-sum(ac)
    x1<-sum(choose(ac,2))
    x2<-sum(choose(ak,2))
    q<-x1*x2/choose(n,2)
    (sum(choose(as.vector(a),2))-q)/((x1+x2)/2.0-q)
}
fmeasure<-function(a)
{
    ac<-apply(a,1,sum)
    ak<-apply(a,2,sum)
    FIJ<-function(i,j){
        r<-a[i,j]/ac[i]
        p<-a[i,j]/ak[j]
        rp<-r+p
        if(rp<1.e-6) return (0)
        2*r*p/rp
    }
    F<-a
    for(i in 1:nrow(a)){
        for(j in 1:ncol(a)){
            F[i,j]<-FIJ(i,j)
        }
    }
    sum(apply(F,1,max)*ac/sum(ac))
}

info<-function(x)
{
    p<-x/sum(x)
    p<-p[p>0]
    -sum(p*log(p))
}
rel.info<-function(a)
    ##a is matrix of CxK
{
    ac<-apply(a,1,sum)
    ak<-apply(a,2,sum)
    HC<-info(ac)
    if(HC==0) return(1)
    HCK<-sum(apply(a,2,info)*ak/sum(ak))
    1-HCK/HC
}
