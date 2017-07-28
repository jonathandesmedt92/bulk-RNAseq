####################################################
# RNA sequencing analysis protocol              ####
# R implementation for Linux bash commands      ####
# Version 14042017                              ####
####################################################


######################################
# Code chunk 0: Set parameters     ###             
######################################

source("https://bioconductor.org/biocLite.R")


# Single-end or paired-end?
paired=TRUE
# Adapter sequences?
adapter = NULL
adapter_rev = NULL
#How many low-quality bases should be trimmed at the 5' end?
trim.5 = 0
#How many low-quality bases should be trimmed at the 3' end?
trim.3 = 10
#What is the minimum read length required?
minlength = 35
#What are the number of CPU cores to be used for parallel computing?  (a numeric or "max")
cores = "max"
#Human genome index directory
hgi = "/mnt/nfs/data/database/human/star2.5.1index/"
#twopassmode?
twopassMode = "Basic"
#alignIntronMax?
alignIntronMax = 1000000
#alignIntronMin?
alignIntronMin = 21



######################################
# Code chunk 1: Quality control    ###
# Program used: FastQC             ###
######################################

files<-list.files()
files<-files[grep(pattern = "*fastq.gz",files)]

for (i in 1:length(files)){
  system(paste("fastqc ",files[i],sep = ""))
  print(paste("Quality control completed on fastq file ",i," of ",length(files),sep = ""))
}


#################################################################
# Code chunk 2: Adapter removal and low-quality base removal  ###
# Program used: Cutadapt                                      ###
#################################################################

removeext<-function(x){
  tmp<-NULL
  for(i in 1:length(x)){
    tmp[i]<-strsplit(x[i],"\\.")
  }
  tmp<-data.frame(tmp, stringsAsFactors = F)
  tmp<-as.character(tmp[1,])
  return(tmp)
}

cutadapt <-
  function( fastqfiles, adapter = "ATCG" , adapter_rev = NULL, quality.cutoff=c(0,0) , minlength = 10 , cores="max"){
    library(parallel)
    if(cores=="max"){cores=detectCores()-1}
    if(length(fastqfiles) < cores){cores=length(fastqfiles)}
    outnames<-paste0(basename(removeext(fastqfiles)),"_clipped.fastq")
    if(quality.cutoff[1]==0&quality.cutoff[2]==0){
      quality<-""
    }
    if(quality.cutoff[1]==0&quality.cutoff[2]!=0){
      quality<-paste(" -q ",quality.cutoff[2],sep="")
    }
    if(quality.cutoff[1]!=0&quality.cutoff[2]==0){
      quality<-paste(" -q ",quality.cutoff[1],",0",sep="")
    }  
    if(quality.cutoff[1]!=0&quality.cutoff[2]!=0){
      quality<-paste(" -q ",quality.cutoff[1],",",quality.cutoff[2],sep="")
    }  
    if(is.null(adapter)&is.null(adapter_rev)){
      adapters<-""
    }  
    if(is.null(adapter)&!is.null(adapter_rev)){
      adapters<-paste(" -A ", adapter_rev, sep = "")
    }
    if(!is.null(adapter)&is.null(adapter_rev)){
      adapters<-paste(" -a ", adapter, sep = "")
    }
    if(!is.null(adapter)&!is.null(adapter_rev)){
      adapters<-paste(" -a ", adapter," -A ", adapter_rev, sep = "")
    }
    a<-mclapply(1:length(fastqfiles) , function(x) system(paste("cutadapt -O 1"," -m ",minlength,quality ,adapters," ",fastqfiles[x],">",outnames[x], sep = "")), mc.cores=cores , mc.preschedule=F)
    return(outnames)
  }

cutadapt(fastqfiles = files, adapter = adapter, adapter_rev = adapter_rev, quality.cutoff = c(trim.5, trim.3), minlength = minlength, cores = cores)



###############################################
# Code chunk 3: Alignment                   ###
# Program used: STAR                        ###
###############################################

#system("nohup nice -n 19 STAR --runMode genomeGenerate --runThreadN 16 --genomeDir ./ --genomeFastaFiles hg.GRCh38.headerFormatted.fa --sjdbGTFfile gencode.v24.annotation.gtf --sjdbOverhang 100 &")
#is already done but need to check how to do it again

clipped.files<-list.files()
clipped.files<-clipped.files[grep(pattern = "_clipped.fastq",clipped.files)]

for(i in 1:length(clipped.files)){
  system(paste("gzip ",clipped.files[i],sep = ""))
}

clipped.files<-list.files()
clipped.files<-clipped.files[grep(pattern = "_clipped.fastq.gz",clipped.files)]

#Next loop takes a while; when it's finished check the directory as process continues in background
#approx 3 minutes
for (i in 1:length(clipped.files)){
  system(paste("nohup nice -n 19 STAR --genomeDir ",hgi,"  --runThreadN 30 --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outFilterType BySJout --quantMode GeneCounts --twopassMode ",twopassMode," --outFilterMultimapNmax 20 --alignIntronMax ",alignIntronMax," --alignMatesGapMax 1000000 --alignIntronMin ",alignIntronMin," --chimSegmentMin 0 --outFilterMismatchNoverLmax 0.1 --outFilterMismatchNmax 8 --alignSJoverhangMin 8 --alignSJDBoverhangMin 3 --outReadsUnmapped None --readFilesCommand zcat --readFilesIn ",clipped.files[i]," --outFileNamePrefix ",removeext(clipped.files[i]),"_ > star.running.log &",sep = ""))
  print(paste("Fastq file ",i," of ",length(clipped.files)," aligned to genome.",sep = ""))
  Sys.sleep(600)
}

###############################################
# Code chunk 4: Bam-to-Sam conversion       ###
# Program used: Rsamtools                   ###
###############################################

library(Rsamtools)

aligned.files<-list.files()
aligned.files<-aligned.files[grep(pattern = "Aligned.sortedByCoord",aligned.files)]

for(i in 1:length(aligned.files)){
  system(paste("samtools view -h -o ",removeext(aligned.files[i]),".sam ",aligned.files[i], sep = ""))
  print(paste("Bam-file ",i," of ", length(aligned.files)," converted to sam-file.",sep = ""))
}


#######################################################
# Code chunk 5: Generating count matrix             ###
# Program used: featureCounts (Rsubread package)    ###
#######################################################

#not advised local install of packages in case of no writing permission (lib path needs to be specified)
#install.packages("/mnt/nfs/data/jonathan/S-J-inc/ToSend/2_Quality_trimming/1_RawFastqFiles/Rsubread_1.24.2.tar.gz", repos = NULL, type = "source", lib = "/mnt/nfs/data/jonathan/S-J-inc/ToSend/2_Quality_trimming/1_RawFastqFiles/")
library(Rsubread, lib.loc = "/mnt/nfs/data/jonathan/S-J-inc/ToSend/2_Quality_trimming/1_RawFastqFiles/")

sam.files<-list.files()
sam.files<-sam.files[grep(pattern = ".sam",sam.files)]

counts<-featureCounts(files = sam.files,
                      annot.inbuilt = "hg38",
                      strandSpecific = 2)
#Extract counts
count_mat<-data.frame(counts[[1]])

#Summmate lane splits
count_mat2<-NULL

for(i in 1:12){
  l1<-i*4-3
  l2<-i*4-2
  l3<-i*4-1
  l4<-i*4
  new<-as.numeric(count_mat[,l1])+as.numeric(count_mat[,l2])+as.numeric(count_mat[,l3])+as.numeric(count_mat[,l4])
  count_mat2<-cbind(count_mat2,new)
  rownames(count_mat2)<-rownames(count_mat)
  colnames(count_mat2)[i]<-colnames(count_mat)[l1]
}

count_mat2<-data.frame(count_mat2, stringsAsFactors = F)

write.table(count_mat2, "RawReadsCount.txt", sep="\t")

##########################################################################
# Code chunk 6: Differential gene expression and quality assessment    ###
# Program used: EDASeq, EdgeR                                          ###
##########################################################################

# Remove absent genes (less than 1 counts-per-million)

exp.tresholds<-unname(colSums(count_mat2)/1000000)
absent.genes<-count_mat2

for(i in 1:ncol(absent.genes)){
  absent.genes[,i]<-as.numeric(absent.genes[,i]<exp.tresholds[i])
}
absent.genes[,"absent"]<-NULL
for(i in 1:nrow(absent.genes)){
  absent.genes[i,"absent"]<-as.logical(prod(absent.genes[i,1:12]))
  print(i/nrow(absent.genes))
}

present.counts<-count_mat2[!absent.genes[,"absent"],]

# Within-sample normalisation (GC-content)

#install.packages("/mnt/nfs/data/jonathan/S-J-inc/ToSend/2_Quality_trimming/1_RawFastqFiles/EDASeq_2.8.0.tar.gz", repos = NULL, type = "source", lib = "/mnt/nfs/data/jonathan/S-J-inc/ToSend/2_Quality_trimming/1_RawFastqFiles/")


# Between-sample normalisation (library size and RNA composition)


# EdgeR

# Convert ENTREZ to SYMBOL

library(org.Hs.eg.db)

mappedIDs<-data.frame(cbind(rownames(present.counts),mapIds(org.Hs.eg.db, keys=rownames(present.counts), column="SYMBOL", keytype = "ENTREZID")))

present.counts[,"GeneSymbol"]<-mappedIDs[,2]
present.counts<-present.counts[!is.na(present.counts$GeneSymbol),]



##########################################################################
# Code chunk 6: Differential gene expression and quality assessment    ###
# Program used: EDASeq, EdgeR                                          ###
##########################################################################


