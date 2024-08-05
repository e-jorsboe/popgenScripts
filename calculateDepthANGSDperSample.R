## For more on depth see:
## evernote: Admixture proportions from NGS data

args<-commandArgs(TRUE)

if(length(args)==0){

    print("Arguments have to be supplied: ")
    print("depthSample=(single depthSample file)")
    q()
}

getArgs<-function(x,l){
  unlist(strsplit(grep(paste("^",x,"=",sep=""),args,val=T),"="))[2]
}

depthSample<-getArgs("depthSample",args)

print("depthSample:")
print(depthSample)

print("IF ANGSD IS NOT RUNG WITH A .sites FILE IT IS DEPTH OF COVERAGE - IF .sites FILE IT IS DEPTH (OF PRE-DEFINED COVERAGE)!")
      

depth<-read.table(depthSample,as.is=T,h=F)

depthPerSample<-apply(depth,1,function(x) (sum(0:100*x)/sum(x)))

print("Depth per sample:")
print(depthPerSample)

print("Mean depth for all samples:")
print(mean(depthPerSample))
