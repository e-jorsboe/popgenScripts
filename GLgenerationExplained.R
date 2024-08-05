getLikes<-function(x,d=2,e=0.01,norm=FALSE){                                                                                                                                                                        
  n<-length(x)                                                                                                                                                                                                      
  ## simulates depth - poission distributed                                                                                                                                                                         
  dep<-rpois(n,d)                                                                                                                                                                                                   
  ## x denotes number of 's B (biallelic A and B), then this is number of A with genotype and errorrate                                                                                                             
  nA<-rbinom(n,dep,c(e,0.5,1-e)[x+1])                                                                                                                                                                               
  ## then gets the probability or likelihood                                                                                                                                                                        
  res<-rbind(dbinom(nA,dep,e),dbinom(nA,dep,0.5),dbinom(nA,dep,1-e))                                                                                                                                                
  if(norm){                                                                                                                                                                                                         
    res<-t(t(res)/colSums(res))                                                                                                                                                                                     
  }                                                                                                                                                                                                                 
  return(res)                                                                                                                                                                                                       
} 

getLikes2<-function(x,d=2,e=0.01,norm=FALSE,haploid=F){                                                                                                                                                                        
  n<-length(x)                                                                                                                                                                                                      
  ## simulates depth - poission distributed
  dep<-rpois(n,d)                                                                                                             
  if(haploid){
    ## x denotes number of B's (biallelic A and B), then this is number of A with genotype and errorrate
    nA<-rbinom(n,dep,c(e,1-e)[x+1])
    ## then gets the probability or likelihood
    res<-rbind(dbinom(nA,dep,e),dbinom(nA,dep,1-e))
  } else{
    ## x denotes number of 's B (biallelic A and B), then this is number of A with genotype and errorrate
    ## because when geno == 0 meaning BB, then we only get A when there is an error
    ## this generate reads so either A or not A, depends on underlying genotype (assumed true genotype)
    nA<-rbinom(n,dep,c(e,0.5,1-e)[x+1])
    ## then gets the probability (same as bionomial density) 
    ## of observing this many As, given the genotype, or prob of observing this many As
    ## given binomial with depth flips and true genotype and error rate
    ## order BB, BA, AA
    res<-rbind(dbinom(nA,dep,e),dbinom(nA,dep,0.5),dbinom(nA,dep,1-e))
  }
  if(norm){
    res<-t(t(res)/colSums(res))
  }
  return(res)
} 
