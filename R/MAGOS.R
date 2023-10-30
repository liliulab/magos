### April 2020.
#Added mag.multiple.run with fold
#mag.exp.var.v3....comment edep=mean(edep.vec)...

### packages:
library(igraph)
library(dplyr)
library(Matrix)


mag.dist.v3<- function(vafs, depths){
  vaf<- mean(vafs[vafs>0.001])
  if(sum(vafs>0.001)==0){
    vaf<- mean(vafs)
  }
  depth<- mean(depths)

  s1<- depth*vaf
  s2<- depth-s1

  #s1/(s1+s2)
  #R<- dbeta(vafs, shape1=s1, shape2=s2, log=T)
  R <-c(dbeta(vafs[vafs!=0.001], shape1=s1, shape2=s2, log=T),
        dbeta(vafs[vafs==0.001], shape1=s1, shape2=s2, log=T)) #log should be F- because if Log=T we could have positive and negative values. We should add log
  #mean(R)
  #plot(c(vafs[vafs!=0.001], vafs[vafs==0.001]),R)
  #rbind(vafs[vafs!=0.001], R[1:58])
  #dbeta(vafs, shape1=s1, shape2=s2, log=T)
  #v<- mad(vafs, center=mean(vaf.test), constant = 1)
  v<-sum((vafs-vaf)^2)/length(vafs)
  #v<- sqrt(v)
  range<- max(vafs)-min(vafs)
  #range<- sqrt(range)
  #range<- range^2
  #cat(v, ' ',range,' ',v*range,'\n')
  #print(range)
  #print(v*range)
  v <- ifelse(v<0.0005, 0.0005, v)
  range<- ifelse(range<0.005,0.005,range)

  if(mean(R)>0){
    R <- R/(range*v)#;mean(R)
  }
  if(mean(R)<0){
    R<- -R*range*v
  }
  dist <- -mean(R)#/exp(5*(weight+0.0001))
  #dist<- ifelse(dist>0, -0.001,dist)

  results<- c(dist=dist, s1=s1, s2=s2)
  #print(paste(results, 'end'))
  return(results)
}



#mARCH 1st
fit.elements.reduce.v3<- function(elements,vaf.data){#, depth.data){
  #  fit.obj<- c()
  # fit.minval<- c()
  #dist.val<-c()
  #dist.shape<- c()
  result <- c();
  for(i in 1:(dim(vaf.data)[2]-1)){
    #print(i)
    variants<- vaf.data[elements, i]
    #depths<- depth.data[elements, i]

    #dist<- mag.dist.v3(variants,depths)
    dist<- max(variants)- min(variants)
    result <- rbind(result, dist)
    # dist.val<- c(dist.val, dist['dist'])
    # dist.shape<- c(dist.shape, dist['s1'])

    #fit1<-optimize(y=variants,exp.depth=exp.depth ,f=optim.obj, interval = c(1,exp.depth-1),tol=0.001)
    #fit.obj<- c(fit.obj, fit1$object)
    #fit.minval<- c(fit.minval, fit1$minimum)
  }

  #return(c(max(fit.obj),fit.minval[which.max(fit.obj)]))
  #return(c(max(dist.val),dist.shape[which.max(dist.val)]))
  #return(result[which.max(result[,'dist']), c('dist', 's1')])
  return(max(result))
}
fit.elements.v3<- function(elements,vaf.data, depth.data){
  #  fit.obj<- c()
  # fit.minval<- c()
  #dist.val<-c()
  #dist.shape<- c()
  result <- c();
  for(i in 1:(dim(vaf.data)[2]-1)){
    #print(i)
    variants<- vaf.data[elements, i]
    depths<- depth.data[elements, i]

    dist<- mag.dist.v3(variants,depths)
    result <- rbind(result, dist)
    # dist.val<- c(dist.val, dist['dist'])
    # dist.shape<- c(dist.shape, dist['s1'])

    #fit1<-optimize(y=variants,exp.depth=exp.depth ,f=optim.obj, interval = c(1,exp.depth-1),tol=0.001)
    #fit.obj<- c(fit.obj, fit1$object)
    #fit.minval<- c(fit.minval, fit1$minimum)
  }

  #return(c(max(fit.obj),fit.minval[which.max(fit.obj)]))
  #return(c(max(dist.val),dist.shape[which.max(dist.val)]))
  return(result[which.max(result[,'dist']), c('dist', 's1')])
}
get.var.v3.1 <-function(x, x.nb, s) {
  x.nb<- x.nb[as.logical(x)]
  a <- mean(x.nb);
  b <- var(x.nb);
  c <- sum(as.logical(x));
  d <- s;
  return(c(a, b, c, d));
}
mag.exp.var.v3<-function(efrq.vec, edep.vec, num, n=500, type='M'){
  expt.var<-c()
  depth.keep<- c()
  vafs<- c()
  #depth.range<- max(edep.vec)-min(edep.vec)
  #if(sd(edep.vec)<10){edep.vec<- mean(edep.vec)}
  a<-quantile(edep.vec,0.1)
  b<-quantile(edep.vec,0.9)
  lowerf.1 <- quantile(efrq.vec,0.15)
  higherf.1 <- quantile(efrq.vec,0.85);


  if(length(edep.vec)>5){
    edep.vec<- edep.vec[edep.vec>=a&
                          edep.vec<=b]
  } else {
    edep.vec<-c(mean(edep.vec), mean(edep.vec))
  }
  if(length(efrq.vec)>5){
    efrq.vec<- efrq.vec[efrq.vec>=lowerf.1&#quantile(efrq.vec,0.05)&
                          efrq.vec<=higherf.1]#quantile(efrq.vec,0.95)]
  } else {
    efrq.vec<-c(mean(efrq.vec), mean(efrq.vec))
  }
  #print('asjdahskdjhasjkd')
  #print(length(edep.vec))
  for(i in 1:n){
    edep<- sample(edep.vec,1, replace=T)
    #edep=mean(edep.vec)
    #efrq<- sample((efrq.vec),1, replace=T)
    efrq<- mean(efrq.vec)
    s1<-efrq*edep
    s2<-(1-efrq)*edep
    x<-rbeta(num,shape1 = s1,shape2 = s2)
    #hist(x, 50)
    #var(x)
    vafs<- c(vafs, s1/(s1+s2))
    expt.var<- c(expt.var, var(x))
    depth.keep<- c(depth.keep, edep)
    i<-i+1
  }
  return(list(exp.vat=(expt.var), depths=depth.keep, vafs=vafs))
}


########################################
########################################          Prep data and other functions

mag.prepdata<- function(arg.data){
  num<- dim(arg.data)[2]/2   # number of samples.     ### in the format of ref, alt, ref, alt
  cat('Number of samples: ', num)
  arg.data<- data.frame(arg.data)
  vaf.data<- c()
  depth.data<- c()
  for(i in 1:num){
    depth<-arg.data[,i*2]+arg.data[,i*2-1]
    vaf<- arg.data[,i*2]/depth

    depth.data<- cbind(depth.data, depth)
    colnames(depth.data)[i]<-paste("depth.", i, sep='')
    vaf.data<- cbind(vaf.data, vaf)
    colnames(vaf.data)[i]<-paste("vaf.", i, sep='')
  }

  ### arg data should have the frequencies in rows and columns for samples. for 3 samples, 100 mutation: 3 columns,100rows.
  vaf.data <- as.matrix(vaf.data)
  vaf.data <- round(vaf.data,3)
  vaf.data <- ifelse(vaf.data<1e-3, 1e-3, vaf.data)
  vaf.data <- ifelse(vaf.data> 1-(1e-3), 1-(1e-3), vaf.data)


  counts.data<- data.frame(arg.data, row.names=NULL)
  vaf.data<- data.frame(vaf.data, row.names=NULL)
  depth.data<- data.frame(depth.data, row.names=NULL)

  ID<- c(1:dim(counts.data)[1])

  counts.data$ID<- ID
  vaf.data$ID<- ID
  depth.data$ID<- ID

  results<- list(counts=counts.data,vafs=vaf.data, depths=depth.data)
  return(results)
}







########################################
########################################          MAG Multiple V3.1


mag.reduce.graph<-function(prep.data,x.el=matrix(), n=-100, threshold=0.03){

  time0<- as.numeric(Sys.time())
  #sampleNum<-1  #added oct 11 2018

  order.var<- apply(prep.data$vafs[,-which(names(prep.data$vafs)=='ID')],1,prod)
  order.ind<- order(order.var)
  sorted.order.var<- order.var[order.ind]

  vaf.data<- prep.data$vafs
  depth.data<- prep.data$depths

  vaf.sorted<- vaf.data[order.ind,]
  depth.sorted<- depth.data[order.ind,]

  IDs<- vaf.sorted[, 'ID']


  # if(dim(x.el)[1]==0){
  #   l<- dim(vaf.data)[1]
  #   x.el<- data.frame(diag(l))
  #   names(x.el)<- IDs
  # }
  #FEB 21 2019
  if(dim(x.el)[1]==1 & dim(x.el)[2]==1){
    l<- dim(vaf.data)[1]
    x.el<- data.frame(diag(l))
    x.el<- sapply(x.el, as.logical)
    colnames(x.el)<- IDs
  }


  mean.frqs<-t(apply( x.el, 1, function(x){mean(sorted.order.var[x])}))
  x.el.sorted<-x.el[order(mean.frqs),]
  x.el.final<- c()


  #x.nb<- data.sorted[,sampleNum]
  #ss.other<- data.sorted[,-c(sampleNum, dim(data.sorted)[2]), drop=F]

  size.threshold<-10

  l<-nrow(x.el)
  full<-floor(l/size.threshold)
  rem<-l%%size.threshold
  row.fold<-0

  if( rem < 5 ){ full<- full-1; rem<- rem+size.threshold}


  # optim.value.bound<-fit.variants.weight(c(0.4,0.46,0.05))[1]


  #optim.value.bound<-fit.elements.v2(c(TRUE,TRUE,TRUE), cbind(c(0.4,0.46, 0.43),c(1,2,3)))[1]


  #depth.bound<- min(apply(depth.sorted[,-dim(depth.sorted)[2]],2,mean))



  # optim.value.bound<- fit.elements.v3(c(TRUE,TRUE), cbind(c(0.4,0.43),c(1,2)),
  #                                    cbind(c(depth.bound,depth.bound),c(1,2)))[1]
  optim.value.bound<-threshold



  ind.full<- c(1:size.threshold)
  combn.el<- t(combn(ind.full,2))

  if(full>0){
    for(i in 1:full){

      fold<- c((1+(i-1)*size.threshold):(i*size.threshold))
      x.el.fold<-x.el.sorted[fold, ]
      #print(i)
      #initial.full<-apply(combn.el, 1, function(x){fit.elements.weight.reduce( x.el.fold[x[1],],
      #                                                                        x.el.fold[x[2],],
      #                                                                       x.nb, ss.other )})

      #### march 1st:
      #size.fold<- ifelse(min(apply(x.el.fold,1,sum))<2, 2,
      #  min(apply(x.el.fold, 1,sum)))
      vafs.fold<- vaf.sorted[apply(x.el.fold, 2, function(x) sum(x)==1),]
      #depths.fold<- depth.sorted[apply(x.el.fold, 2, function(x) sum(x)==1),]


      initial.full<- apply(combn.el, 1, function(x){
        #fit.elements.v2(as.logical(x.el.fold[x[1],]+x.el.fold[x[2],]),data.sorted)})
        #fit.elements.v3(as.logical(x.el.fold[x[1],]+x.el.fold[x[2],]),vaf.sorted, depth.sorted)})
        #fit.elements.v3(x.el.fold[x[1],]|x.el.fold[x[2],],vaf.sorted, depth.sorted)})
        #MARCH 1st
        fit.elements.reduce.v3(x.el.fold[x[1],]|x.el.fold[x[2],],vaf.sorted)})


      tempMatLoop<-matrix(100, nrow=size.threshold, ncol=size.threshold)
      tempMatLoop.2<-lower.tri(diag(100, nrow=size.threshold, ncol=size.threshold))

      #MARCH 1st
      #tempMatLoop[tempMatLoop.2] <- initial.full[1, ]
      tempMatLoop[tempMatLoop.2] <- initial.full
      cat('**********')
      cat(class(tempMatLoop))
      mat.loop<-Matrix::t(tempMatLoop)
      mat.loop.g<- mat.loop
      mat.loop.g[mat.loop <= optim.value.bound]<- 1
      mat.loop.g[mat.loop> optim.value.bound]<- 0
      sum(mat.loop.g==1)
      g1<-graph_from_adjacency_matrix(mat.loop.g, mode = 'undirected')
      complete_points<- max_cliques(g1, min=2)
      if(length(complete_points)>0){
        x.el.graph<-c()
        el.rm.graph<-c()
        for(i in 1:length(complete_points)){
          el.rm<- as.numeric(complete_points[[i]])
          el.rm<- el.rm[!el.rm%in%el.rm.graph]
          if(length(el.rm)>0){
            x.el.graph<- rbind(x.el.graph, as.logical(colSums(x.el.fold[el.rm,, drop=F])))
            el.rm.graph<- c(el.rm.graph, el.rm)
          }}
        x.el.fold<- rbind(x.el.fold[-el.rm.graph, , drop=F],x.el.graph)
      }
      x.el.final<- rbind(x.el.final, x.el.fold)
      {  #while(  min(mat.loop)< optim.value.bound & nrow(mat.loop)>2){#

        # el.rm <- which(mat.loop==min(mat.loop), arr.ind=T)
        #  el.rm <- el.rm[1,]   # april 14

        #el.rm<- unique(as.numeric(el.rm))

        # mat.loop<- mat.loop[-el.rm, -el.rm, drop=F]
        #if(sum(dim(mat.loop))==0){mat.loop=matrix(c(0,0,0,0),nrow = 2)}
        #x.el.fold<- rbind( x.el.fold[-el.rm, , drop=F], colSums(x.el.fold[el.rm,]))
        #colnames(x.el.fold)<-row.names(data.sorted)
        #}
      }
    }}

  rems<- c((l-rem+1):l)
  #print(rems)
  x.el.rem<- x.el.sorted[rems,]
  ind.rem<- c(1:rem)
  combn.rem<- t(combn(ind.rem, 2))
  # initial.rems<-apply(combn.rem, 1, function(x){fit.elements.weight.reduce( x.el.rem[x[1],],
  #                                                                           x.el.rem[x[2],],
  #                                                                           x.nb, ss.other )})
  initial.rems<- apply(combn.rem, 1, function(x){
    #fit.elements.v2(as.logical(x.el.rem[x[1],]+x.el.rem[x[2],]),arg.data)})
    #fit.elements.v3(as.logical(x.el.rem[x[1],] +x.el.rem[x[2],]),vaf.sorted, depth.sorted)})
    #fit.elements.v3(x.el.rem[x[1],]|x.el.rem[x[2],],vaf.sorted, depth.sorted)})
    fit.elements.reduce.v3(x.el.rem[x[1],]|x.el.rem[x[2],],vaf.sorted)})



  ###may 6th
  tempMatLoop<-matrix(100, nrow=rem, ncol=rem)
  tempMatLoop.2<-lower.tri(diag(100, nrow=rem, ncol=rem))
  #tempMatLoop[tempMatLoop.2] <- initial.rems[1, ]
  tempMatLoop[tempMatLoop.2] <- initial.rems
  cat('**********')                                 
  cat(class(tempMatLoop))
  mat.loop<-Matrix::t(tempMatLoop)
  #diag(mat.loop) <- 1e5;
  mat.loop.g<- mat.loop
  mat.loop.g[mat.loop <= optim.value.bound]<- 1
  mat.loop.g[mat.loop> optim.value.bound]<- 0
  sum(mat.loop.g==1)
  g1<-graph_from_adjacency_matrix(mat.loop.g, mode = 'undirected')
  complete_points<- max_cliques(g1, min=2)
  x.el.graph<-c()
  el.rm.graph<-c()

  #print(complete_points)
  if(length(complete_points)>0){
    for(i in 1:length(complete_points)){
      el.rm<- as.numeric(complete_points[[i]])
      el.rm<- el.rm[!el.rm%in%el.rm.graph]
      if(length(el.rm)>0){
        x.el.graph<- rbind(x.el.graph,as.logical(colSums(x.el.rem[el.rm,,drop=F])))
        el.rm.graph<- c(el.rm.graph, el.rm)
      }}
    x.el.rem<- rbind(x.el.rem[-el.rm.graph, , drop=F],x.el.graph)
    colSums(x.el.rem)}
  x.el.final<- rbind(x.el.final, x.el.rem)
  {
    #while(  min(mat.loop)< -10 &nrow(mat.loop)>2){#
    #  el.rm<- which(mat.loop==min(mat.loop), arr.ind=T)
    #  el.rm<- unique(as.numeric(el.rm))

    # mat.loop<- mat.loop[-el.rm, -el.rm, drop=F]
    #  if(sum(dim(mat.loop))==0){mat.loop=matrix(c(0,0,0,0),nrow = 2)}
    #  x.el.rem<- rbind( x.el.rem[-el.rm, , drop=F], colSums(x.el.rem[el.rm,]))
    #  colnames(x.el.fold)<-row.names(data.sorted)
    #}
    #x.el.final<- rbind(x.el.final, x.el.rem)
  }

  results<- list(x.el.red=x.el.final, vaf.data=vaf.data,vaf.sorted=vaf.sorted,
                 depth.data=depth.data,depth.sorted=depth.sorted);
  nn <- nrow(results$x.el.red)
  print("nn")
  print(nn)
  time1<- as.numeric(Sys.time())
  print('TIME:')
  time<- time1-time0
  print(time)
  #jj<-jj+1
  if(nn < size.threshold | nn == n) {
    return(results)
  } else {
    #mag.reduce.2(results, sampleNum=sampleNum,n=nn)
    mag.reduce.graph(prep.data = prep.data,x.el=x.el.final,n=nn, threshold = threshold)
  }



}




### added okay calculation to mag.mult; Also, added the first x.el to x.el.list as the first element, so s should start from 2.
mag.multiple<- function(prep.data,x.el.cut=data.frame()) {


  if(1>dim(prep.data$vafs)[2]-1){print("stop the algorithm and check the sampleNum; run mag.prepdata() first")}

  time0<-as.numeric(Sys.time())
  params<-c();  x.values<-c();  x.steps<-c();  row.loop<-c();x.step<-c()
  x.elements.list <- list()


  ## Sep2019:
  okay.list=list()



  if(ncol(prep.data$vafs) == 2){        ##### change '1' to '2' because the 2nd column is always id
    print('executing single sample clustering....'); flush.console();

  }

  order.var<- apply(prep.data$vafs[,-which(names(prep.data$vafs)=='ID')],1,prod)
  order.ind<- order(order.var)
  sorted.order.var<- order.var[order.ind]

  vaf.data<- prep.data$vafs
  depth.data<- prep.data$depths

  vaf.sorted<- vaf.data[order.ind,]
  depth.sorted<- depth.data[order.ind,]

  IDs<- vaf.sorted[, 'ID']






  #oct 18 2018    I WILL KEEP THE SORTED VERSION.

  n<- dim(vaf.data)[1]


  if( sum(dim(x.el.cut))==0){
    print("Mag.reduce has not been run.")
    x.elements<-data.frame(diag(n))
    #names(x.elements)<-c(1:length(x.nb))   feb 26
    #names(x.elements)<- IDs     APRIL 11 2019
    x.elements<- sapply(x.elements, as.logical)  #APRIL 11 2019
    colnames(x.elements)<- IDs
  }else{
    x.elements<- x.el.cut
  }

  x.elements.list[[1]]=x.elements


  ### Sep 2019
  okay.id= c(1:nrow(x.elements))
  okay.list[[1]]= okay.id



  l<- nrow(x.elements)

  #### new matloop with reduced x.el:
  indx.cut<-c(1:nrow(x.elements))
  combn.cut <- t(combn(indx.cut, 2))

  ###no need to run weights...
  cat('initial fitting ... '); flush.console();
  time1 <- as.numeric(Sys.time())
  #initial.1<-apply(combn.cut, 1, function(x){fit.elements.weight( x.elements[x[1],], x.elements[x[2],],
  #                                                               x.nb, arg.data.other )})
  # x.elements<- sapply(x.elements, as.logical)
  #x.elements <- matrix(as.logical(as.matrix(x.elements)), ncol=ncol(x.elements), byrow=T)

  initial.1<- apply(combn.cut, 1, function(x){
    fit.elements.v3(c(x.elements[x[1],]|x.elements[x[2],]),vaf.sorted, depth.sorted)})

  time2 <- as.numeric(Sys.time())
  cat('took ', time2-time1, ' seconds.\n'); flush.console();

  tempMatLoop<-Matrix(1000000, nrow=l, ncol=l, sparse = T)  #may 6th   DEC 2018: added sparse
  tempMatLoop.2<-lower.tri(diag(1000000, nrow=l, ncol=l))
  tempMatLoop[tempMatLoop.2] <- initial.1[1, ]
  cat('**********')
  cat(class(tempMatLoop))
  mat.loop<-Matrix::t(tempMatLoop)
  #diag(mat.loop) <- 1e5;


  tempPar1Loop<-lower.tri(diag(1000000, nrow=l, ncol=l))
  #tempPar2Loop<-lower.tri(diag(0, nrow=l, ncol=l))
  tempPar1Loop[tempPar1Loop] <- initial.1[2, ]
  #tempPar2Loop[tempPar2Loop] <- initial.1[3, ]

  par.1.loop<-t(tempPar1Loop)
  #par.2.loop<-t(tempPar2Loop)
  #weights.step<-c()
  #var.mean <- c()


  time2<-as.numeric(Sys.time())

  s<-2
  while(nrow(x.elements) > 1){
    cat(s, '...', nrow(x.elements), '...\n'); flush.console();

    x.values[s]<-min(mat.loop)
    el1<-which(mat.loop==min(mat.loop), arr.ind = T)[1,1]
    el2<-which(mat.loop==min(mat.loop), arr.ind = T)[1,2]
    el.rm<-c(el1,el2)
    params<-c(params, par.1.loop[el1,el2])#, par.2.loop[el1,el2]))

    #what happenes at each step

    ### FEB 19 2019
    #x.step <- x.elements[el.rm[1], ] + x.elements[el.rm[2],]   ###this is very clever! I dont need to update the xloop.
    x.step <- x.elements[el.rm[1], ]|x.elements[el.rm[2],]   ###this is very clever! I dont need to update the xloop.



    ###only use x.elements rows for each cluster!
    x.steps <- rbind(x.steps, x.step)
    x.elements <- rbind(x.elements[-el.rm,], x.step)
    x.elements.list[[s]]<-x.elements


    ### Sep2019:
    #update.okay:
    okay.id=c(okay.id[-el.rm],0)
    okay.list[[s]]=okay.id



    #####UPDATE MATRIX HERE
    mat.loop<-mat.loop[-el.rm, -el.rm, drop=F]
    par.1.loop<-par.1.loop[-el.rm, -el.rm, drop=F]


    mat.column.update<-c()
    par1.column.update<-c()
    #par2.column.update<-c()

    x1 <- x.elements[1:(nrow(x.elements)-1), ,drop=F]
    x2 <- x.elements[nrow(x.elements), , drop=F]   #the ones that are put together
    print('here')
    dim(x.elements)
    print(dim(x1))
    print(dim(x2))
    print('there')
    #fit.between <- t(apply(x1, 1, fit.elements.weight, element.2=x2,
    #                      arg.data=x.loop, arg.data.other=arg.data.other))

    #oct 2018
    fit.between<- apply(x1, 1, function(x){
      #fit.elements.v2(as.logical(x+x2), data.sorted)} )
      #fit.elements.v3(as.logical(x+x2), vaf.sorted, depth.sorted)} )
      #FEB 2019 19

      fit.elements.v3(x|x2, vaf.sorted, depth.sorted)} )






    mat.column.update <- fit.between[1,]
    par1.column.update <- fit.between[2,]
    #par2.column.update <- fit.between[, 3]

    #print(mat.loop)
    if(nrow(mat.loop) != 0){
      mat.loop<-cbind(mat.loop, mat.column.update )
      mat.loop<-rbind(mat.loop, rep(1000000,dim(mat.loop)[2]))    #just to keep mat.loop square.
      par.1.loop<-cbind(par.1.loop, par1.column.update)
      par.1.loop<-rbind(par.1.loop, rep(1000000, dim(par.1.loop)[2]))
      # par.2.loop<-cbind(par.2.loop, par2.column.update)
      # par.2.loop<-rbind(par.2.loop, rep(0, dim(par.2.loop)[2]))
      #print("MATLOOP UPDATE SHODE")
      #print(mat.loop)
      #print("END")
    }
    # print(mat.loop)

    s<-s+1
  }
  # colnames(var.mean)<-c("mean", "var","NumberOfPoints", "step", "weights")

  time3 <- as.numeric(Sys.time())
  cat('loop took ', time3 - time2, ' seconds.\n')
  cat('total took ', time3 - time0, ' seconds.\n')

  result <- list(x.el=x.elements.list,vaf.sorted=vaf.sorted,depth.sorted=depth.sorted,
                 prep.data=prep.data, okay.list=okay.list)
  return(result)
}






###Updated SEP2019: okay calculation from MAG
mag.var<- function(mag.out){

  #varmean is geneated based on sorted data
  x.el.list<- mag.out$x.el
  data<- data.frame(mag.out$vaf.sorted)
  #this also assumes that the data was ran on prep.data code so the last column has ID.
  num<- dim(data)[2]-1
  var.mean.list<-list()

  okay.list=mag.out$okay.list



  # SEP 2019: no need used okay.2 instead
  # x.el.red<- x.el.list[[1]]
  #  x.el.red.okay<- apply(x.el.red, 1, function(x) paste(x, collapse = ''))  #april 15
  #
  #
  #
  # SEP19: the x.el.red.IDs should be number of rows, not number of columns.
  # x.el.red.IDs= c(1:nrow(x.el.red))
  #

  for ( i in 1:num){
    var.mean<-c()
    x.data<- data[,i]
    names(data)
    depth.sorted.each<- mag.out$depth.sorted[,i]
    #s<-7
    for( s in 2:length(x.el.list)){
      x.el<- x.el.list[[s]]

      #sep19:
      okay.2=okay.list[[s]]

      temp<- t(apply(x.el, 1, get.var.v3.1, x.nb=x.data, s=s))

      # okay<- apply(x.el, 1, function(x){y=paste(x,collapse = ''); ind=x.el.red.okay%in%y;
      # res=ifelse(sum(ind)==0,0,x.el.red.IDs[ind])
      # return(res)})



      #print(cbind(okay.2, okay))  # SEP 2019: TESTED PERFECTLY! CAN JUST USE OKAY.2

      okay=okay.2

      depths<- apply(x.el,1, function(x) round(mean(depth.sorted.each[as.logical(x)])))



      #temp <- cbind(temp, weights.temp, okay)
      temp <- cbind(temp, depths ,okay)
      var.mean<- rbind(var.mean, temp)
    }

    #colnames(var.mean)<-c("mean", "var","NumberOfPoints", "step","weights","okay")
    colnames(var.mean)<-c("mean", "var","NumberOfPoints", "step","depth","okay")
    var.mean.list[[i]]<- var.mean
    names(var.mean.list)[i]<-paste(names(data)[i],".var.mean",sep="")


  }

  return(list(var.mean.all=var.mean.list, mag.out=mag.out))
}





####add all depths:


### UPDATED PROBABILIES.
cut.off.multiple  <-function(mag.var){
  time0 <- as.numeric(Sys.time())
  #var.mean <- mag.output$var.mean
  mag.var.list<-mag.var$var.mean.all

  ### FEB 13:
  depths.data<- mag.var$mag.out$prep.data$depths
  vafs.data<- mag.var$mag.out$prep.data$vafs

  for( i in 1:length(mag.var.list)){
    temp<- mag.var.list[[i]]
    temp<- data.frame(temp, row.names = NULL)
    temp[is.na(temp)]<-0
    temp$var<- round(temp$var, 6)
    temp$var<- temp$var+1e-7
    mag.var.list[[i]]<-temp
  }
  names(mag.var.list)<- names(mag.var$var.mean.all)
  x.el.list<- mag.var$mag.out$x.el


  flag <- FALSE
  m <- c()
  str <- c()  #clone structure line from temp var.mean matrix
  str.f <- c();  cutstr <- c();cutstr2 <- c()

  i <- max(mag.var.list[[1]]$step)
  # allPoints <- mag.output$x.el[[i]]
  allPoints<- x.el.list[[i]]   #feb28
  ok.points <- allPoints - allPoints  ##### To generate clean slate
  n <- length(allPoints)
  #	print(allPoints)
  id.ok <-c() #clusters with okay variance!
  temp.mv <- c(); temp.points <- c(); temp.step <- c();temp.point.step <- c();

  while(flag == FALSE){
    print(i); flush.console();
    ############
    ## I need to use x.el to get the breaks, because using StepI would cause problem
    ## in identical clusters and mutations.
    ## TO FIX THIS I JUST KEPT MY ORIGINAL METHOD -
    ## I just duplicated the rows in the case of identical mutations
    ## because there should be onlytwo rows at each step
    ## and having only 1 row means identical clusters.
    ############
    v <- i
    stepI.x.el <- x.el.list[[i]]
    stepI.x.el.step <- stepI.x.el[nrow(stepI.x.el), ]  #Only look at the most recent cluster
    stepII.x.el <- x.el.list[[i-1]]
    # stepII.breaks.x.el <- dplyr::setdiff(stepII.x.el, stepI.x.el)
    stepII.breaks.x.el <- dplyr::setdiff(as.data.frame(stepII.x.el), as.data.frame(stepI.x.el))

    temp.points <- rbind(temp.points, c(i, ok.points));
    temp.step <- rbind(temp.step, c(i, stepI.x.el.step));
    temp.point.step <- rbind(temp.point.step, c(i, sum(ok.points & stepI.x.el.step)))

    #print(sum( ok.points & stepI.x.el.step)==0)

    keep<- c(0,0)
    for( j in 1:length(mag.var.list)){
      var.mean<- mag.var.list[[j]]
      stepI <- var.mean[var.mean$step==i,]
      stepII <- var.mean[var.mean$step==i-1, ]  ## previous step

      if(sum(ok.points & stepI.x.el.step) == 0){     ###meaning stepI.x.el is not in ok.points
        #stepII.breaks <- setdiff(stepII[,-c(4,6)], stepI[,-c(4,6)])
        stepII.breaks <- dplyr::setdiff(stepII[,-4], stepI[,-4])
        if(nrow(stepII.breaks)==1){
          # print("why did this happen????")
          stepII.breaks <- rbind(stepII.breaks, stepII.breaks)
        }
        stepII.num <- unique(stepII[, 4])
        stepII.breaks <- cbind(stepII.breaks, step=stepII.num)   #MAY 9th
        # stepII.breaks <- cbind(stepII.breaks, stepII[1:2,c(4,6)])   #just add the step number back to the matrix        #############MARCh 14 add "okay"
        #for( ii in 1:dim(newtemp)[1]){              ###stepII.breaks always has 2 rows
        for(ii in 1:2){
          #sample.sim <- cluster.sim(edep=depth, efrq=stepII.breaks[ii,1], shape=shape, n=10000)
          var.sim <- c()
          size <- ifelse(stepII.breaks$NumberOfPoints[ii]<5, 5, stepII.breaks$NumberOfPoints[ii])
          depth.vec<-depths.data[depths.data$ID%in%names(stepII.breaks.x.el[ii,as.logical(stepII.breaks.x.el[ii,]),drop=F]),j]

          vaf.vec<-vafs.data[vafs.data$ID%in%names(stepII.breaks.x.el[ii,as.logical(stepII.breaks.x.el[ii,]),drop=F]),j]

          #var.sim<- mag.exp.var.v2.1(efrq = stepII.breaks[ii,1], edep.vec =depth.vec , num=size)
          #var.sim<- var.sim$exp.vat

          if(mean(depth.vec)>150){type<-'L'}
          if(mean(depth.vec)<=150){type<-'L'}
          if(mean(depth.vec)<=40){type<-'L'}

          #print(type)

          var.sim<- mag.exp.var.v3(efrq.vec = vaf.vec, edep.vec =depth.vec , num=size, type=type)
          var.sim<- var.sim$exp.vat
          #print('done mag')
          # print(sd(var.sim))


          m <- mean(var.sim) + 1e-6
          vv <- sd(var.sim)
          vmax <- max(var.sim)

          if(length(vaf.vec)>50 & max(vaf.vec)-min(vaf.vec)>0.05){
            ### added the second condition to ignore the clusters with all equal vafs.
            var.check<- var(vaf.vec[-which.max(vaf.vec)[1]])
          }else{
            var.check<- stepII.breaks[ii,2]
          }


          var.threshold<- m+2*vv
          #print(var.threshold)
          if(type=='H'){var.threshold<- m+3*vv}
          #   if(var(vaf.vec)!=0){
          #  test<-varTest(vaf.vec,sigma.squared = m+3*vv)
          # }
          #if(var(vaf.vec)==0){
          #  test<- list(p.value=1)
          #}

          var.threshold<- m+3*vv
          if(var.threshold>0.001){
            var.threshold<- round(var.threshold,3)
          }

          if(stepII.breaks[ii, 'okay'] >0 |
             #floor(var.check*10^4)/10^4 <  ceiling((round(m+3*vv,4)+0.00001)*10^4)/10^4) {
             floor(var.check*10^4)/10^4< m+3*vv  ) {

            keep[ii]<-keep[ii]+1
          }#changed to vmax on oct 14 2018
        }
      }
    }

    for(ii in 1:2){
      if(keep[ii]==length(mag.var.list)){
        cutstr2 <- rbind(cutstr2, stepII.breaks.x.el[ii,])
        #id.ok<-c(id.ok, stepII.breaks$id[ii])
        str <- rbind(str, stepII.breaks[ii,])
        ok.points <- ok.points + stepII.breaks.x.el[ii,]}
    }

    if(sum(ok.points>1)>0){
      print("ERROR")
    }
    if(sum(ok.points) == n){
      flag<-TRUE
    } else {
      i<-(v-1)
    }
    #if(v==i){iii<-TRUE}else{i<-(v-1)}
  }

  time1 <- as.numeric(Sys.time())
  cat("took ", (time1 - time0), " seconds.\n");
  #print(cutstr)
  colors<-colSums(cutstr2*c(1:nrow(cutstr2)))+1
  data<- mag.var$mag.out$vaf.sorted
  final.data<-cbind(data,colors)


  #PROBS:
  t=final.data

  t.1=c()
  for(i in 1:(ncol(t)-2)){
    temp=c()
    for(j in unique(t$colors)){
      #print(j)
      m.1=mean(t[t$colors==j,i])
      v.1=var(t[t$colors==j,i])

      a <- ((1 - m.1) / v.1 - 1 / m.1) * m.1 ^ 2
      b <- a * (1 / m.1 - 1)

      temp=rbind(temp,c(j, a, b))
    }
    colnames(temp)=c('colors',paste0('alpha.',i), paste0('beta.',i))
    t.1=cbind(t.1,temp[,-1])
  }
  t.1=cbind(unique(t$colors), t.1)
  colnames(t.1)[1]='colors'

  #head(t)

  t_means= t%>%select(contains(c('vaf','color')))%>%group_by(colors)%>%summarise_all(mean)
  i=1

  prob.all=c()
  for(i in 1:nrow(t)){
    vf=t[i,1:(ncol(t)-2)]
    #t.1
    prob.s=c()
    #j=1
    for( j in 1:length(vf)){
      #j is the sample number

      probs=c()
      prob=c()
      for(ii in unique(t$colors)){
        v=vf[,j]
        shape1 = t.1[t.1[,1]==ii,1+(2*j-1)]
        shape2=t.1[t.1[,1]==ii, 1+2*j]
        clust_mean= t_means[t_means$colors==ii, j+1]%>%pull()

        prob_flag= ifelse(v<clust_mean, TRUE, FALSE )

        pp=pbeta(v, shape1, shape2 , lower.tail = prob_flag)

        prob= c(prob, round(pp,3))
        #prob=c(prob,dbeta(vf[,j], shape1 = t.1[t.1[,1]==ii,1+(2*j-1)], shape2=t.1[t.1[,1]==ii, 1+2*j]))

      }
      prob.s=rbind(prob.s, prob)
      prob.s=apply(prob.s, 2,min )
    }
    prob.all=rbind(prob.all, prob.s)
  }
  colnames(prob.all)=paste0('colors.', unique(t$colors))
  rownames(prob.all)=NULL
  final.probs=cbind(t, prob.all)






  results<-list(str=str, cutclust=cutstr, cutclust2=cutstr2, colors=colors,probs=final.probs ,final.data=final.data)
  return(results)
}


########################################
########################################          MAG SINGLE V3
mag.single<- function(prep.data){

  if(is.null(prep.data$vafs$ID)){ return(print('Run mag.prepdata on the data'))}
  time0 <- as.numeric(Sys.time())


  vaf.data<- prep.data$vafs
  vaf.data<- round(vaf.data,3)    #change 1 (aug6,2018)
  vaf.data[,-dim(vaf.data)[2]]<- ifelse(vaf.data[,-dim(vaf.data)[2]] < 1e-3, 1e-3, vaf.data[,-dim(vaf.data)[2]])
  vaf.data[,-dim(vaf.data)[2]]<- ifelse(vaf.data[,-dim(vaf.data)[2]] > 1-1e-3, 1-1e-3, vaf.data[,-dim(vaf.data)[2]])

  depth.data<- prep.data$depths


  x.nb <- vaf.data;
  if(!is.null(ncol(vaf.data)) && ncol(vaf.data) > 0) {
    x.nb <- vaf.data$vaf.1
  }
  params<-c();  x.values<-c();  x.steps<-c();  row.loop<-c(); x.step<-c()
  x.elements.list<-list()
  x.elements<-data.frame(diag(length(x.nb)))



  sort.x<-sort(vaf.data$vaf.1)

  #vaf.depth<- merge(prep.data$vafs,prep.data$depths, by="ID")

  #x.loop<-sort.x
  sort.vaf.data <- vaf.data[order(vaf.data$vaf.1),]
  sort.depth.data<- depth.data[order(vaf.data$vaf.1),]


  names(x.elements)<-sort.vaf.data$ID
  #names(x.el.first)<- sort.vaf.data$ID

  ### APRIL 11 2019
  x.elements<- sapply(x.elements, as.logical)
  colnames(x.elements)<- sort.vaf.data$ID


  x.elements.preprc<- c()
  for(i in unique(sort.x)){
    #temp<-colSums(x.elements[sort.x==i,])
    temp<- as.logical(apply(x.elements[sort.x==i,, drop=F],2, sum ))  # APRIL 11 2019
    x.elements.preprc<- rbind(x.elements.preprc, temp)

  }

  #names(x.elements.preprc)<- sort.vaf.data$ID
  colnames(x.elements.preprc)<- sort.vaf.data$ID   # APRIL 11 2019

  ## names checked! correct!



  #### LOOK HERE:       THIS CHANGES TO X.element.preprc
  #x.el.first<- data.frame(diag(length(x.nb)))
  #x.el.first.okay<- apply(x.el.first, 1, function(x) paste(x, collapse = '')) # april 15


  # SEP 2019 commented:
  #x.el.first.okay<- apply(x.elements.preprc, 1, function(x) paste(as.numeric(x), collapse = '')) # april 15
  #####   Sep 2019: I dont even need x.el.first okay. I just need the initial ids. Then use el.rm to update it.
  #####   sep 2019: c

  x.el.first.okay.id<- c(1:dim(x.elements.preprc)[1])
  l<- dim(x.elements.preprc)[1]

  mat.loop<-matrix(1000000, ncol=l, nrow=l)
  par.1.loop<-matrix(100000, ncol=l, nrow=l)

  s.pairs <- cbind(1:(l-1), 2:l)
  s.likls <- t(apply(s.pairs, 1, function(x){
    fit.elements.v3(x.elements.preprc[x[1],]|x.elements.preprc[x[2],], sort.vaf.data, sort.depth.data)    #### APRIL 11 2019
  }))
  #mat.loop<-diag(0, nrow=l, ncol=l)
  for(i in 2:l){
    mat.loop[i-1,i] <- s.likls[i-1, 1]      ### the s.likls is the likelihood between i-1)th element and the i th.
    par.1.loop[i-1,i] <- s.likls[i-1, 2]
    # par.2.loop[i-1,i] <- s.likls[i-1, 3]
  }

  time1 <- as.numeric(Sys.time())
  var.mean<-c()

  ### aug12 add the step 0 to var.mean
  #temp<- t(apply(x.elements.preprc, 1, get.var, x.nb=sort.x, s=0 ))
  temp<- t(apply(x.elements.preprc, 1, get.var.v3.1, x.nb=sort.x, s=0 ))
  #okay<- apply(x.elements.preprc, 1, function(x) sum(paste(x, collapse = '')%in%x.el.first.okay))

  # okay<- apply(x.elements.preprc, 1, function(x) #dec 6 2018
  # {y=paste(as.numeric(x), collapse = '');  #APRIL 11 2019 added as.numeric
  # ind=x.el.first.okay%in%y;
  # res=ifelse(sum(ind)==0, 0, x.el.first.okay.id[ind])
  # return(res)
  # })
  okay=x.el.first.okay.id

  #weights<- rep(0, length(okay))  ##### just to have weight  APRIL 15

  ###add depth to var.mean    dec 6 2018
  #  depths<- apply(x.elements.preprc,1, function(x) round(mean(sort.depth.data[as.logical(x),'depth.1'])))

  depths<- apply(x.elements.preprc,1, function(x) round(mean(sort.depth.data[x,'depth.1'])))   # APRIL 11 2019


  temp<- cbind(temp,depths ,okay)

  var.mean<- rbind(var.mean, temp)
  s<-1
  x.elements<- x.elements.preprc

  m1=as.matrix(x.elements)
  m.preprc=Matrix(m1, sparse = T)   ###just to save the preproces matrix
  #x.elements.list[[1]]<-m2
  ###need to create new "n' for number of mutations. to find neighbours
  n<-length(sort.x)
  while(nrow(x.elements)>= 2){
    # print("INJA")
    # print(s)
    #print(mat.loop)

    #print(min(mat.loop))
    x.values[s]<-min(mat.loop)
    el1<-which(mat.loop==min(mat.loop), arr.ind = T)[1,1]
    el2<-which(mat.loop==min(mat.loop), arr.ind = T)[1,2]
    el.rm<-c(el1,el2)
    # print("Remove ina")
    # print(el.rm)
    params<-rbind(params,par.1.loop[el1,el2])

    #what happenes at each step
    #x.step <- x.elements[el.rm[1],] +x.elements[el.rm[2],]   ###this is very clever! I dont need to update the xloop.
    x.step <- x.elements[el.rm[1],] |x.elements[el.rm[2],]  # APRIL 11 2019



    ###only use x.elements rows for each cluster!
    x.steps<-rbind(x.steps, x.step)
    #print(x.steps)
    x.elements<-rbind(x.elements[-el.rm,], x.step)


    ########## change to sparse matrix saving:

    m1<- as.matrix(x.elements)
    m2<- Matrix(m1, sparse = T)
    x.elements.list[[s]]<- m2

    #### april 14 add "okay"
    temp<- t(apply(x.elements, 1, get.var.v3.1, x.nb=sort.x, s=s ))          ####for oc.v9 I changed x.nb=x.nb to x.nb=
    #####
    #okay<- apply(x.elements, 1, function(x) nrow(merge(t(x), x.el.first)))               ################# APRIL 14 2018   takeeeesss a long time
    ####################### faster version? :
    # okay<- apply(x.elements, 1, function(x) sum(paste(x, collapse = '')%in%x.el.first.okay))
    #print(s)
    # print('1.8')
    # SEP 19: commented the next 7 lines! I just need one line and use el.rm to update the okay.
    #okay<- apply(x.elements, 1, function(x) #dec 6 2018
    #{y=paste(as.numeric(x), collapse = '');   # APRIL 11 2019 added the numeric
    #ind=x.el.first.okay%in%y;
    #res=ifelse(sum(ind)==0, 0, x.el.first.okay.id[ind])
    #res=ifelse(y%in%x.el.first.okay,x.el.first.okay.id[which(x.el.first.okay==y)],0)
    #return(res)
    #})
    okay= c(okay[-el.rm],0)
    #print('1.9')
    # cbind(okay, c(okay.org[-el.rm],0))
    # okay.org=okay
    # weights<- rep(0, length(okay))  ##### just to have weight  APRIL 15


    ###add depth to var.mean    dec 6 2018
    depths<- apply(x.elements,1, function(x) round(mean(sort.depth.data[as.logical(x),'depth.1'])))

    #temp<- cbind(temp,depths,weights ,okay)
    temp<- cbind(temp,depths ,okay)   # removed the weigths APRIL 11 2019

    var.mean<- rbind(var.mean, temp)

    #####UPDATE MATRIX HERE
    mat.loop<-mat.loop[-el.rm, -el.rm]
    par.1.loop<-par.1.loop[-el.rm, -el.rm]
    #par.2.loop<-par.2.loop[-el.rm, -el.rm]

    mat.column.update<-c()
    par1.column.update<-c()
    #par2.column.update<-c()

    max<-max(which(x.step%in%1))
    min<-min(which(x.step%in%1))
    #		print(paste("min", min))
    nei<-c()
    if(min == 1 & max != n){
      neiIND<-which(x.elements[,max+1]==1)
      nei<-x.elements[neiIND,, drop=F]
      #			print("if1")
    }else if(max == n & min != 1){
      neiIND<-which(x.elements[,min-1]==1)
      nei<-x.elements[neiIND,, drop=F]
      #			print("if2")
    }else if( max!=n & min!=1){
      nei1<-which(x.elements[,min-1]==1)
      nei2<-which(x.elements[,max+1]==1)
      neiIND<-c(nei1,nei2)
      nei<-x.elements[neiIND,, drop=F]
      #			print("if3")
    }
    #print(nei)

    if(length(nei) > 0){
      x.last <- x.elements[nrow(x.elements), , drop=F]

      temp <- apply(nei, 1, function(x){
        #fit.elements.v3(as.logical(x+x.last), vaf.data=sort.vaf.data,depth.data=sort.depth.data)
        fit.elements.v3(x|x.last, vaf.data=sort.vaf.data,depth.data=sort.depth.data)   # APRIL 11 2019

      })

      ###try to only look at cluster neighbours
      ###instead of dim i will use sqrt of length. Because mat.loop is always square matrix   WORKS!!! SEP12
      l2<-sqrt(length(mat.loop))
      mat.column.update<-rep(100000,l2)
      mat.column.update[neiIND]<-temp[1,]
      par1.column.update<-rep(100000, l2)
      par1.column.update[neiIND]<-temp[2,]
      #par2.column.update<-rep(0, l2)
      #par2.column.update[neiIND]<-temp[3,]

      mat.loop<-cbind(mat.loop, mat.column.update )
      mat.loop<-rbind(mat.loop, rep(100000,dim(mat.loop)[2]))    #just to keep mat.loop square.
      par.1.loop<-cbind(par.1.loop, par1.column.update)
      par.1.loop<-rbind(par.1.loop, rep(100000, dim(par.1.loop)[2]))
      #par.2.loop<-cbind(par.2.loop, par2.column.update)
      #par.2.loop<-rbind(par.2.loop, rep(0, dim(par.2.loop)[2]))
    }

    s<-s+1
  }

  #var.mean <-cbind(var.mean, 1/(var.mean[,3]))
  colnames(var.mean)<-c("mean", "var","NumberOfPoints", "step","depth", 'okay')

  time2 <- as.numeric(Sys.time())
  cat("while took ", (time2 - time1), " seconds.\n")
  cat("total took ", (time2 - time0), " seconds.\n")

  result <- list(x.el=x.elements.list,x.el.preprocess=m.preprc,
                 params=params, var.mean=var.mean, prep.data=prep.data,
                 vaf.sorted=sort.vaf.data,depth.sorted= sort.depth.data,
                 freq.s= sort.x, time=list(time0,time1, time2))
  return(result)
}


cut.off.single <- function(mag.output){
  time0 <- as.numeric(Sys.time())
  #print(dim(mag.output))
  var.mean <- mag.output$var.mean
  head(var.mean)
  var.mean <- data.frame(var.mean, row.names = NULL)
  #var.mean$id<-round(var.mean$mean*var.mean$NumberOfPoints,5)
  var.mean[is.na(var.mean)] <- 0
  #var.mean.first <- var.mean
  #print(var.mean)
  var.mean$var <- round(var.mean$var, 6)
  var.mean$var <- var.mean$var + 1e-7
  #id.max<-max(var.mean$id)
  flag <- FALSE
  m <- c()
  str <- c()  #clone structure line from temp var.mean matrix
  str.f <- c()
  cutstr <- c()
  cutstr2 <- c()

  i <- max(var.mean$step)
  ##added aug 12:
  if(i==1){
    flag=T
    str<-var.mean
  }
  #print(flag)
  allPoints <- mag.output$x.el[[i]]   ### added +1 with the AUG6th2018 version, because the x.element.list has an extra element.
  allPoints<- as.data.frame(as.matrix(allPoints))
  ok.points <- allPoints - allPoints  ##### To generate clean slate
  n <- length(allPoints)
  #	print(allPoints)
  id.ok <-c() #clusters with okay variance!
  temp.mv <- c();
  temp.points <- c();
  temp.step <- c();
  temp.point.step <- c();



  var.sim.zero<- mag.exp.var.v3(efrq.vec = mag.output$vaf.sorted$vaf.1, edep.vec =mag.output$depth.sorted$depth.1 , num=n)    #### APRIL 11 2019
  var.sim.zero<- var.sim.zero$exp.vat

  m <- mean(var.sim.zero, na.rm=T) + 1e-6
  vv <- sd(var.sim.zero, na.rm=T)
  if(var(mag.output$vaf.sorted$vaf.1) < m+3*vv ){
    flag<-TRUE
    print('Only one clone detected')
    str<- var.mean[var.mean$step==i,]
    colors<- rep(1, n)
  }



  {
    while(flag == FALSE){
      #print(flag)
      # print('umad?')
      # print(flag)
      #print('enter')
      #print(i); flush.console();
      v <- i

      ############
      ## I need to use x.el to get the breaks, because using StepI would cause problem
      ## in identical clusters and mutations.
      ## TO FIX THIS I JUST KEPT MY ORIGINAL METHOD -
      ## I just duplicated the rows in the case of identical mutations
      ## because there should be onlytwo rows at each step
      ## and having only 1 row means identical clusters.
      ############

      stepI <- var.mean[var.mean$step==i,]
      stepI.x.el <- mag.output$x.el[[i]]
      stepI.x.el<- as.data.frame(as.matrix(stepI.x.el))

      stepI.x.el.step <- stepI.x.el[nrow(stepI.x.el), ]  #Only look at the most recent cluster

      stepII <- var.mean[var.mean$step==i-1, ]  ## previous step
      #stepII$var<-stepII$var,)
      if(i>1){
        stepII.x.el <- mag.output$x.el[[i-1]]
        stepII.x.el<- as.data.frame(as.matrix(stepII.x.el))
      }
      if(i==1){
        stepII.x.el <- mag.output$x.el.preprocess
        stepII.x.el<- as.data.frame(as.matrix(stepII.x.el))

      }
      #print("II")
      #print(stepII.x.el)

      #if( !(id.sum%in%id.ok)){    ###this means if the division on this step is not happening on the "ok" clusters. so the okay clusters should still exist in the temp.
      #if(sum(abs(id.ok-id.sum)<0.1)==0){
      #newtemp<-temp[!(temp$id%in%id.ok),]
      #if(!id.I%in%id.ok){
      temp.points <- rbind(temp.points, c(i, ok.points));
      temp.step <- rbind(temp.step, c(i, stepI.x.el.step));
      temp.point.step <- rbind(temp.point.step, c(i, sum(ok.points & stepI.x.el.step)))


      #print(ok.points)
      #print(sum(ok.points & stepI.x.el.step))
      if(sum(ok.points & stepI.x.el.step) == 0){     ###meaning stepI.x.el is not in ok.points

        #stepII.breaks<- stepII[!stepII$id%in%stepI,]
        #print("stepI ")
        #print(stepI )
        #print("stepII ")
        #print(stepII )
        #stepII.breaks <- setdiff(stepII[,-c(4,7)], stepI[,-c(4,7)])
        stepII.breaks <- dplyr::setdiff(stepII[,-4], stepI[,-4])   #dec 6 2018
        if(nrow(stepII.breaks)==1){
          #print("why did this happen????")
          stepII.breaks <- rbind(stepII.breaks, stepII.breaks)
        }
        stepII.num <- stepII[, 4]
        stepII.breaks <- cbind(stepII.breaks, stepII[1:2,4])   #just add the step number back to the matrix
        stepII.breaks.x.el <- dplyr::setdiff(stepII.x.el, stepI.x.el)

        #for( ii in 1:dim(newtemp)[1]){              ###stepII.breaks always has 2 rows
        for(ii in 1:2){
          #sample.sim <- cluster.sim(edep=depth, efrq=stepII.breaks[ii,1], shape=shape, n=10000)
          var.sim <- c()
          size <- ifelse(stepII.breaks$NumberOfPoints[ii]<5, 5, stepII.breaks$NumberOfPoints[ii])
          #print(size)
          #for(j in 1:500){
          # sample.sim.2 <- sample.sim[sample(1:10000, size=size),]   #temp[1,3]
          #  var.sim <- c(var.sim, var(sample.sim.2[, 4]))
          #}
          #var.sim<- mag.exp.var.v2(efrq=stepII.breaks[ii,1], edep=stepII.breaks[ii,'depth'], num=size)

          #print(stepII.breaks.x.el)

          vafs<- mag.output$vaf.sorted$vaf.1[as.matrix(stepII.breaks.x.el[ii,])]                 #### APRIL 11 2019
          depths<- mag.output$depth.sorted$depth.1[as.matrix(stepII.breaks.x.el[ii,])]

          var.sim<- mag.exp.var.v3(efrq.vec = (vafs) , edep.vec = (depths) , num=size)    #### APRIL 11 2019
          var.sim<- var.sim$exp.vat
          #print(var.sim)

          ###MAY 2019:
          if(length(vafs)>20 & max(vafs)-min(vafs)>0.05){
            #print('UMAD TOOOOOOOOO')
            ### added the second condition to ignore the clusters with all equal vafs.
            var.check<- var(vafs[-which.max(vafs)[1]])
          }else{
            var.check<- stepII.breaks[ii,2]}




          m <- mean(var.sim, na.rm=T) + 1e-6
          vv <- sd(var.sim, na.rm=T)
          vmax <- max(var.sim, na.rm=T)
          #temp.mv <- rbind(temp.mv, c(i, ii, m, vv, vmax, stepII.breaks[ii,2]))

          #if(stepII.breaks[ii, 'NumberOfPoints'] < 11 | stepII.breaks[ii,2]*(1+stepII.breaks[ii, "weights"]) < m+vv) {        #okay cluster
          #if(stepII.breaks[ii, 'okay'] ==1 | (stepII.breaks[ii,2])*(1+stepII.breaks[ii, "weights"]) < m+vv) {        #okay cluster   APRIL 15
          # print(stepII.breaks)
          #print(vv)
          if(stepII.breaks[ii, 'okay'] >0 |
             # stepII.breaks[ii,2] < m+3*vv ) {   #dec 6 2018
             var.check<m+3*vv){
            str.f <- c(str.f, stepII.breaks[ii,1])
            #print("here")
            #print(str)
            str <- rbind(str, stepII.breaks[ii,])
            #print(stepII.breaks.x.el)
            cutstr2 <- rbind(cutstr2, stepII.breaks.x.el[ii,])
            #print(str)

            #id.ok<-c(id.ok, stepII.breaks$id[ii])
            ok.points <- ok.points + stepII.breaks.x.el[ii,]
            #print(ok.points)
          }
        }
      }

      ###stop condition
      #if(abs(sum(str$id)-id.max)<0.000001){flag<-TRUE}else{i<-(v-1)}
      #print(i)      ###new condition: sum of id.okay== id.max

      #if(abs(sum(str$id)-id.max)<0.000001){flag<-TRUE}else{i<-(v-1)}

      ####new condition    oct 31
      #print("LAST PRINT")
      if(sum(ok.points>1)>0){
        print("ERROR")
      }
      if(sum(ok.points) == n){
        flag<-TRUE
      } else {
        i<-(v-1)
      }
      #if(v==i){iii<-TRUE}else{i<-(v-1)}
    }

    for( i in 1:nrow(str)){
      #print(i)
      if(str$step[i]!=0){
        xel <- mag.output$x.el[[str$step[i]]]}
      if(str$step[i]==0){
        xel <- mag.output$x.el.preprocess}

      #print(xel)
      for(j in 1:nrow(xel)){
        m1 <- mean(mag.output$freq.s[xel[j,]==1])    ##changed to sorted. after oc.9#### OCT 30: THis should be changed for OC11 since there is no sorting in that
        if(abs(m1-str[i,1]) < 0.001){
          cutstr<-rbind(cutstr, xel[j,])  ## what's the difference between cutstr and cutstr2????????
        }
      }
    }


    time1 <- as.numeric(Sys.time())
    cat("took ", (time1 - time0), " seconds.\n");
    #print(cutstr)
    colors<-colSums(cutstr*c(1:nrow(cutstr)))+1
  }

  final.data<-cbind(mag.output$vaf.sorted,colors)
  #print('hi')
  #  print(final.data)

  #PROBS:
  t=final.data


  t.1=c()
  for(i in 1:(ncol(t)-2)){
    temp=c()
    for(j in unique(t$colors)){
      #print(j)
      m.1=mean(t[t$colors==j,i])
      v.1=var(t[t$colors==j,i])

      a <- ((1 - m.1) / v.1 - 1 / m.1) * m.1 ^ 2
      b <- a * (1 / m.1 - 1)

      temp=rbind(temp,c(j, a, b))
    }
    colnames(temp)=c('colors',paste0('alpha.',i), paste0('beta.',i))
    t.1=cbind(t.1,temp[,-1, drop=F])
  }
  t.1=cbind(unique(t$colors), t.1)
  colnames(t.1)[1]='colors'

  #head(t)

  t_means= t%>%select(contains(c('vaf','color')))%>%group_by(colors)%>%summarise_all(mean)
  i=1

  prob.all=c()
  for(i in 1:nrow(t)){
    vf=t[i,1:(ncol(t)-2), drop=F]
    #t.1
    prob.s=c()
    #j=1
    for( j in 1:length(vf)){
      #j is the sample number

      probs=c()
      prob=c()
      for(ii in unique(t$colors)){
        v=vf[,j]
        shape1 = t.1[t.1[,1]==ii,1+(2*j-1)]
        shape2=t.1[t.1[,1]==ii, 1+2*j]
        clust_mean= t_means[t_means$colors==ii, j+1]%>%pull()

        prob_flag= ifelse(v<clust_mean, TRUE, FALSE )

        pp=pbeta(v, shape1, shape2 , lower.tail = prob_flag)

        prob= c(prob, round(pp,3))
        #prob=c(prob,dbeta(vf[,j], shape1 = t.1[t.1[,1]==ii,1+(2*j-1)], shape2=t.1[t.1[,1]==ii, 1+2*j]))

      }
      prob.s=rbind(prob.s, prob)
      prob.s=apply(prob.s, 2,min )
    }
    prob.all=rbind(prob.all, prob.s)
  }
  colnames(prob.all)=paste0('colors.', unique(t$colors))
  rownames(prob.all)=NULL
  final.probs=cbind(t, prob.all)





  results<-list(str=str, cutclust=cutstr, cutclust2=cutstr2, colors=colors, final.data=final.data, probs=final.probs)
  return(results)
}






####################################
####################################        MAG RUN

mag.single.run<- function(input.data, fold=F){
  data.prep<- mag.prepdata(input.data)
  purity<- 1
  mag<- mag.single(data.prep)
  cut<- cut.off.single(mag)
  temp<- merge(cut$final.data, data.prep$depths, by='ID')
  #plot(temp$vaf.1, temp$depth.1, col=temp$colors, pch=19)
  sum1<-temp%>%group_by(colors)%>%summarize(meanVAF=mean(vaf.1))
  purity<- 2*max(sum1$meanVAF)
  if(purity>1 & fold==F){purity='There are clusters with frequency higher than 0.5- consider folding.'}
  if(fold==T & sum(sum1$meanVAF>0.5)>0){
    fold.colors<- sum1$colors[sum1$meanVAF>0.5]
    data.folded<- temp
    data.folded$vaf.1[data.folded$colors%in%fold.colors]<- 1-data.folded$vaf.1[data.folded$colors%in%fold.colors]
    data.prep.fold<- data.prep
    data.prep.fold$vafs$vaf.1<- data.folded$vaf.1
    mag<- mag.single(data.prep.fold)
    cut<- cut.off.single(mag)
    temp<- merge(cut$final.data, data.prep$depths, by='ID')
    #plot(temp$vaf.1, temp$depth.1, col=temp$colors, pch=19)
    sum2<- temp%>%group_by(colors)%>%summarize(meanVAF=mean(vaf.1))
    purity<- 2*max(sum2$meanVAF)
  }
  temp<- temp[,c(2,4,3,1)]
  res<- list(mag=mag, cut=cut, results=temp, fold=fold, purity=purity)
  return(res)
}

# APRIL 2020: added fold.
mag.multiple.run<- function(input.data, reduce=T, threshold=0.03, fold= F){
  prep= mag.prepdata(input.data)
  purity=c()

  if(reduce){
    red=mag.reduce.graph(prep, threshold = threshold)
    mag=mag.multiple(prep, x.el.cut = red$x.el.red)
    mv=mag.var(mag)
    cut=cut.off.multiple(mv)
  }
  if(!reduce){
    #red=mag.reduce.graph(prep)
    mag=mag.multiple(prep)#, x.el.cut = red$x.el.red)
    mv=mag.var(mag)
    cut=cut.off.multiple(mv)
  }

  temp=cut$final.data
  #plot(temp$vaf.1, temp$vaf.2, col=temp$colors, xlim=c(0,1), ylim=c(0,1))
  sum1=temp%>%group_by(colors)%>%summarise_all(mean)
  purity=apply(sum1, 2, max)*2
  purity=purity[c(-1, -length(purity))]

  fs=which(sum1[,-c(1, ncol(sum1))]>0.5, arr.ind = T)
  print('what')
  print(purity)
  if(fold==F & any(purity>1)){purity='There are clusters with frequency higher than 0.5- consider folding.'
  print(purity)}

  if(fold==T & nrow(fs)>0){
    fold.prep= prep
    for( i in 1:nrow(fs)){
      bad_col= as.numeric(sum1[fs[i,1],1 ])
      bad_sam= fs[i,2]

      bad_ID=temp$ID[temp$colors==bad_col]
      head(temp)

      fold.prep$vafs[fold.prep$vafs$ID%in%bad_ID,bad_sam]=1-fold.prep$vafs[fold.prep$vafs$ID%in%bad_ID,bad_sam]
    }

    #plot(fold.prep$vafs$vaf.1, fold.prep$vafs$vaf.2, xlim=c(0,1), ylim=c(0,1))


    if(reduce){
      red=mag.reduce.graph(fold.prep, threshold = threshold)
      mag=mag.multiple(fold.prep, x.el.cut = red$x.el.red)
      mv=mag.var(mag)
      cut=cut.off.multiple(mv)
    }
    if(!reduce){
      #red=mag.reduce.graph(prep)
      mag=mag.multiple(fold.prep)#, x.el.cut = red$x.el.red)
      mv=mag.var(mag)
      cut=cut.off.multiple(mv)
    }
  }
  temp=cut$final.data
  #plot(temp$vaf.1, temp$vaf.2, col=temp$colors, xlim=c(0,1), ylim=c(0,1))
  sum1=temp%>%group_by(colors)%>%summarise_all(mean)
  purity=apply(sum1, 2, max)*2
  purity=purity[c(-1, -length(purity))]

  res=list(results=temp, purity=purity, reduce=reduce, mag=mag, cut=cut)
  return(res)

}





####################################
####################################        CNV ASSIGNMENT




## determine the best cluster of each SNV
## input:   cl.cnv is the output from mag.cn.block()
##
mag.cn.call.tcn=function(cut, cn.input, cl.cnv, upper=10){
  cn.output=c()
  temp=c()
  temp2=data.frame(cut$final.data%>%group_by(colors)%>%summarise_all(funs(mean)))

  for(m in 1:nrow(cn.input[[1]])){  ## each snv in CNV region
    res=c()
    for(c in unique(temp2$colors)){  ## each cluster containing a snv
      for(s in 1:length(cn.input)){  ## each sample
        cn.sample=cn.input[[s]][m, ]
        cn.tmp=cn.sample$cn
        cn.vaf=cn.sample$vaf
        cn.blk=cn.sample$block
        blk.cl = cl.cnv[which(cl.cnv$block == cn.blk), 'cl.cnv']
        vaf.tcn=temp2[temp2$colors==blk.cl, s+1]
        for(k in 1:max(min(upper, cn.tmp), 1)){
          vaf.sc=temp2[temp2$colors==c, s+1]
          exp.vaf=2*k*vaf.sc/(2*vaf.tcn*cn.tmp+(1-2*vaf.tcn)*2)
          res=rbind(res,c(abs(cn.vaf-exp.vaf),s,c,blk.cl,k))
        }
      }
    }
    colnames(res)=c('error', 'sample', 'cl.snv', 'cl.cnv', 'k')
    agg = aggregate(error ~ cl.snv + cl.cnv + k, data=res, 'sum')
    best = agg[which.min(agg$error), ]
    cn.output = rbind(cn.output, best)
  }
  if(nrow(cut$final.data)<nrow(cn.input[[1]])){
    print('Warning: number of CNV events is greater than SNVs. CNV assignment may not be reliable.')}
  return(data.frame(cn.output))
}

## determine the best cluster a CNV belongs to
## it is the majority vote of all SNVs in a CNV block
## input: cut is from mag
## input: cn.input is a list of data frames.
## each data frame has 3 columns: vaf, cn that is tcn, and block that is the CNV block ID.
## input: upper is the max cnv it searches
## output: a data frame with 3 columns: cluster a cnv belongs to, block ID of a CNV, number of SNVs that support the CNV cluster assignment
mag.cn.block=function(cut, cn.input, upper=10){
  cn.output=c()
  temp=c()
  temp2=data.frame(cut$final.data%>%group_by(colors)%>%summarise_all(funs(mean)))

  for(m in 1:nrow(cn.input[[1]])){  ## each snv in CNV region
    res=c()
    for(c in unique(temp2$colors)){  ## each cluster containing a snv
      for(r in unique(temp2$colors)){  ## each cluster containing a cnv region
        for(s in 1:length(cn.input)){  ## each sample
          cn.sample=cn.input[[s]][m, ]
          cn.tmp=cn.sample$cn
          cn.blk=cn.sample$block
          cn.vaf=cn.sample$vaf
          for(k in 1:max(min(upper, cn.tmp), 1)){
            vaf.sc=temp2[temp2$colors==c, s+1]
            vaf.tcn=temp2[temp2$colors==r, s+1]
            exp.vaf=2*k*vaf.sc/(2*vaf.tcn*cn.tmp+(1-2*vaf.tcn)*2)
            res=rbind(res,c(abs(cn.vaf-exp.vaf),s,c,r,k,cn.blk))
          }
        }
      }
    }
    colnames(res)=c('error', 'sample', 'cl.snv', 'cl.cnv', 'k', 'block')
    agg = aggregate(error ~ cl.snv + cl.cnv + k + block, data=res, 'sum')
    best = agg[which.min(agg$error), ]
    cn.output = rbind(cn.output, best)
  }
  block.best.cl = do.call(rbind, lapply(split(cn.output, cn.output$block), function(x) {
    a = aggregate(cl.snv ~ cl.cnv + block, data=x, 'length')
    a[which.max(a$cl.snv), ]
  }))

  if(nrow(cut$final.data)<nrow(cn.input[[1]])){
    print('Warning: number of CNV events is greater than SNVs. CNV assignment may not be reliable.')}
  return(block.best.cl)
}





