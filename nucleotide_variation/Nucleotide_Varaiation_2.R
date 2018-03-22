########################WORKSPACE SETUP####################################

setwd("~/Documents/UAB/Module3/Molecular Population Genetics")
source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings")                 
require(Biostrings)
#library(plyr)
bseqs<-readDNAStringSet("AM.fasta")
split.bseqs<-strsplit(as.character(x=bseqs),split="")
df.split<-as.data.frame(split.bseqs)

########################FILE_WRITING######################################
############Specify file names before running code #######################

thetaFile="theta_CO_50000.wig"
cat("fixedStep  chrom=2L  start=800001  step=50000",file=thetaFile,append=FALSE,sep="\n")
SFile="S_CO_50000.wig"
cat("fixedStep  chrom=2L  start=800001  step=50000",file=SFile,append=FALSE,sep="\n")
sfsFile="SFS_CO_50000.txt"
cat("Site Frequency Spectrum", file=sfsFile,append=FALSE,sep="\n")
piFile="Pi_50000.wig"
cat("fixedStep  chrom=2L  start=800001  step=50000",file=piFile,append=FALSE,sep="\n")

############################REQUIRED INPUTS################################

#The window, total number of nucleotides, total number of sequences and frames 
#are all required inputs for the following code to run.
s<-50000                                
total.s<-length(bseqs[[1]])                      
n<-length(bseqs)                                
u<-total.s/s                             
start<-seq(1,total.s,s)
end<-seq(s,total.s,s)
f<-cbind(start,end)

################################FUNCTION##################################

#theta.count() function takes start and end position per window to compute 
#minor allele frequencies computed across all n sequences; NAs excluded

theta.count<-function(start,end){               
        SNP.minor.na<-c()
        SNP.minor<-c()
        for(i in start:end){
                SNP.freq<-round(nucleotideFrequencyAt(bseqs,at=i)/sum(nucleotideFrequencyAt(bseqs,at=i)),digits=2)
                if (1 %in% SNP.freq == FALSE) {
                        SNP<-min(SNP.freq[SNP.freq>0])
                        SNP.minor.na<-append(SNP.minor.na,SNP)
                }
        }
        SNP.minor<-SNP.minor.na[!is.na(SNP.minor.na)]
        return (SNP.minor)
}

#theta.func takes minor allele frequency count (SNP.minor) to compute theta. Output is a wig file.
#theta.func also computes S number of segregating sites per nucleotide from minor allele frequency count.
#S output is a wig file.

theta.func<-function(count, v){                 
        seg.sites<-length(count)
        ss<-seg.sites/s                         
        wat.num<-rep(1,n-1)
        for (i in 1:(n-1)){
                wat.num[i]<-1/(i)
        }
        wat.num<-sum(wat.num)
        theta<-ss/wat.num                                       
        cat(theta,file=thetaFile,append=TRUE,sep="\n")          
        cat(seg.sites,file=SFile,append=TRUE,sep="\n")          
}

#SFS parses minor allele frequencies in the appropriate bin, computing SNP site frequency distribution
#output is a txt file

SFS<-function(vector,v) {
        a<-sum(vector<=0.1)                                     
        b<-sum(vector<=0.2 & vector>0.1)
        c<-sum(vector<=0.3 & vector>0.2)
        d<-sum(vector<=0.4 & vector>0.3)
        e<-sum(vector<=0.5 & vector>0.4)
        newrow<-c(a,b,c,d,e)
        cat(newrow,file=sfsFile,append=TRUE,sep="\n")
}

#get.cnt performs efficient pattern matching to compare sequences pairwise and identify when bases do 
#not match, indicating a SNP.

get.cnt <- function(seq1,seq2){                                  
        base <- "ATCG"
        tmp.cnt <- 0
        for(i in seq(1:length(seq1))){
                if(grepl(seq1[i],base) == TRUE & grepl(seq1[i],base)){
                        if(seq1[i] != seq2[i]){
                                tmp.cnt = tmp.cnt + 1
                        }
                }
        }
        return (tmp.cnt)
}


#getdifs takes in sequences pairwise, (seq1 and seq2) sequence names extracted from 
#colnames of df.split, and a windowsize. Sequence identifier is retrived for every window.
get.difs <- function(seq1,seq2,s){
        results <- c(rep(0,u))
        counter <- 1
        for(i in start){
                results[[counter]] <- get.cnt(as.vector(df.split[seq1][i:((i+s)-1),]),as.vector(df.split[seq2][i:((i+s)-1),]))
                counter <- counter + 1
        }
        return (results)
}


#pairs computes all combns (combinations in pairs) of sequences. There are 435 when n=30. To save memory, only the column name 
#is sent to get.difs

pairs <- combn(names(df.split),2)
final <- c(rep(0,u))

# the for loop adds the counts from each pairing into the running total stored in final

for(i in seq(1:ncol(pairs))){
        final <- mapply("+",(get.difs(pairs[,i][1],pairs[,i][2],s)),final)
}

#Pi is computed from get.count and get.difs output stored in final. Pair.seq 
#is the number of combinations resulting from pairwise comparisons. Pi computes nucleotide diversity.
#theta.count theta.func and SFS are called within the specified window

for(i in seq(1:length(final))){
        pair.seq<-(n*(n-1))/2 
        pi.avg<-final[i]/pair.seq 
        pi <- pi.avg/s 
        cat(pi,file=piFile,append=TRUE,sep="\n")  
}

for (v in 1:(length(f)/2)){
        start.t=as.numeric(f[v,1])
        end=as.numeric(f[v,2])
        temp<-theta.count(start.t,end)
        theta.func(temp,v)
        SFS(temp,v)        
}











