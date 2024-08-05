
args<-commandArgs(trailingOnly=T)

if(length(args)==0){
    print("Converts a gzipped .gen file (.gen.gz) to a recessive BIMBAM .mgf file, where minor allele is counted as recessive!")
    print("Output will be prefixREC.mgf")
    print("1st argument is .gen.gz file")
    q()
}

fileName<-args[1]

outConRec<-file(sub("\\.gen.gz","\\REC.mgf",fileName),"w")
outConAdd<-file(sub("\\.gen.gz","\\ADD.mgf",fileName),"w")
outConAnn<-file(sub("\\.gen.gz","\\.ann",fileName),"w")

con <- gzcon(file(fileName, "rb"))

n<-0

####################################################################################################
####################################################################################################
## converts .gen file to .mgf file ADD and REC with dosage in direction of minor allele
## the .gen file is expected to have 1:54676_T_C in first column
## EXAMPLE
## .gen: ... A G p(AA) p(AG) pGG) 
## if ( p(AG) + 2*p(GG) >= 0.5 ) --> A minor, G major --> mgf: ...,A,G,(p(AG) + 2*p(AA)) - what assumed by GEMMA
## else --> G minor, A major --> mgf: ...,G,A,(p(AG) + 2*p(GG)) - FLIPPING
## in .mgf file minor allele should be first (2nd col in file) and is the one that is counted in the dosage
####################################################################################################
####################################################################################################


##con <- gzfile(fileName, "r")
while ( TRUE ) {
    line <- readLines(con, n = 1)

    if ( length(line) == 0 ) {
        break
    }
     
    line1<-unlist(strsplit(line," "))

    ## 5 first cols are chr, id, pos A0, A1
    gprobs<-as.numeric(line1[6:length(line1)])
    
    if(length(gprobs) %% 3 != 0){
        print("Something is wrong with .gen file gprobs not in pairs of three (WT,HE,HO)")
        q()
    }

    ## dose is in direction of A1 (5th col)
    dose<-0
    ## how many indis with data (might have NA as value)
    indis<-0
    for(i in 1:(length(gprobs)/3)){

        if(is.na(gprobs[(i*3-2)])){
            next
        }
        dose<-dose+gprobs[(i*3-1)]+gprobs[(i*3)]*2
        indis<-indis+1
        
    }

    ## get frequency by dividing dosage by number of chr
    freq<-dose/(2*indis)
    
    ## if major allele is second allele then flip!
    if(freq>=0.5){        
        ## alleles have to be swapped around       
        add<-as.numeric(line1[seq(from=6,to=length(line1),by=3)])*2+as.numeric(line1[seq(from=7,to=length(line1),by=3)])
        add[is.na(add)]<-NA
        
        rec<-as.numeric(line1[seq(from=6,to=length(line1),by=3)])
        rec[is.na(rec)]<-NA

        write(paste(c(line1[c(2,4,5)],add),collapse=","), outConAdd,append=TRUE)        
        write(paste(c(line1[c(2,4,5)],rec),collapse=","), outConRec,append=TRUE)
        chr<-unlist(strsplit(line1[1],":"))[1]
        pos<-unlist(strsplit(unlist(strsplit(line1[1],":"))[2],"_"))[1]
        id<-line1[2]
        write(paste0(id,", ",pos,", ",chr), outConAnn,append=TRUE)

        
    } else{
        add<-as.numeric(line1[seq(from=8,to=length(line1),by=3)])*2+as.numeric(line1[seq(from=7,to=length(line1),by=3)])
        add[is.na(add)]<-NA
                      
        rec<-as.numeric(line1[seq(from=8,to=length(line1),by=3)])
        rec[is.na(rec)]<-NA
                
        write(paste(c(line1[c(2,5,4)],add),collapse=","), outConAdd,append=TRUE)        
        write(paste(c(line1[c(2,5,4)],rec),collapse=","), outConRec,append=TRUE)

        chr<-unlist(strsplit(line1[1],":"))[1]
        pos<-unlist(strsplit(unlist(strsplit(line1[1],":"))[2],"_"))[1]
        id<-line1[2]
        write(paste0(id,", ",pos,", ",chr), outConAnn,append=TRUE)
        
    }

}    

close(con)
close(outConRec)
close(outConAdd)
