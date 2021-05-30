library(groHMM)
library(rtracklayer) #Import bed files
library(BRGenomics) #To merge GRanges

##################################Useful###################################################
rm(arabisopsis)
##################################Functions########################################

plot_hist<-function(ranges,title){
  hist(width(ranges),breaks=100,main=title,xlim=c(0,150000))
}


extract_characteristics<-function(ranges,genomeLength){
  print(paste("Transcript count:",length(width(ranges))))
  show(quantile(width(ranges)))
  print(paste("Transcript average length:",mean(width(ranges))))
  print(paste("Transcript median length:",median(width(ranges))))
  print(paste("Transcript density:",sum(width(ranges)/genomeLength)))
}

preprocess_drosophila<-function(gRange){
  processed<-gRange
  mcols(processed)['score']<- as.integer(mcols(processed)$name)
  mcols(processed)['name'] <- "*"
  return (processed)
}

write_wiggle<-function(gRange,path){
  writeWiggle(reads=gRange, file=paste(path,"_Plus.wig",sep=""), 
              fileType="wig", strand="+", reverse=FALSE)
  writeWiggle(reads=gRange, file=paste(path,"_Minus.wig",sep=""), 
              fileType="wig", strand="-", reverse=TRUE)
}

write_bigWig<-function(gRange,path){
  seqinfo(gRange)<-NULL
  writeWiggle(reads=gRange, file=paste(path,"_Plus.bw",sep=""), 
              fileType="BigWig", strand="+", reverse=FALSE)
  writeWiggle(reads=gRange, file=paste(path,"_Minus.bw",sep=""), 
              fileType="BigWig", strand="-", reverse=TRUE)
}

add_strand<-function(gRange){
  temp<-gRange
  strand(temp[score(temp)>0])<-'+'
  strand(temp[score(temp)<0])<-'-'
  score(temp)<-abs(score(temp))
  return(temp)
}

export_hmm<-function(hmm,path){
  transcripts<-hmm$transcripts
  export.bed(transcripts,paste(path,".bed",sep=""))#Save to file
  export.bed(transcripts[strand(transcripts)=="+"], paste(path,"_Plus.bed",sep=""))#Save to file
  export.bed(transcripts[strand(transcripts)=="-"], paste(path,"_Minus.bed",sep=""))#Save to file
}
##################################MCF_7 Processing###########################################
MCF7Path = "../data/MCF-7/"
MCF7LiftedPath = paste(MCF7Path,"lifted/",sep="")
#MUST concatenate data with duplicate ranges if duplicate scores

MCF7VehicleRep1 <- import(paste(MCF7LiftedPath,"GRO-seq_Vehicle_rep1.bed",sep=""), format="bed", genome = 'hg19')
MCF7VehicleRep2 <- import(paste(MCF7LiftedPath,"GRO-seq_Vehicle_rep2.bed",sep=""), format="bed", genome = 'hg19')
MCF7E2_10mRep1 <- import(paste(MCF7LiftedPath,"GRO-seq_E2_10m_rep1.bed",sep=""), format="bed", genome = 'hg19')
MCF7E2_10mRep2 <- import(paste(MCF7LiftedPath,"GRO-seq_E2_10m_rep2.bed",sep=""), format="bed", genome = 'hg19')
MCF7E2_40mRep1 <- import(paste(MCF7LiftedPath,"GRO-seq_E2_40m_rep1.bed",sep=""), format="bed", genome = 'hg19')
MCF7E2_40mRep2 <- import(paste(MCF7LiftedPath,"GRO-seq_E2_40m_rep2.bed",sep=""), format="bed", genome = 'hg19')
MCF7E2_160mRep1 <- import(paste(MCF7LiftedPath,"GRO-seq_E2_160m_rep1.bed",sep=""), format="bed", genome = 'hg19')
MCF7E2_160mRep2 <- import(paste(MCF7LiftedPath,"GRO-seq_E2_160m_rep2.bed",sep=""), format="bed", genome = 'hg19')

MCF7 <- mergeGRangesData(MCF7VehicleRep1,MCF7VehicleRep2,
                          MCF7E2_10mRep1,MCF7E2_10mRep2,
                          MCF7E2_40mRep1,MCF7E2_40mRep2,
                          MCF7E2_160mRep1,MCF7E2_160mRep2,
                          field = "score",exact_overlaps = FALSE)

MCF7 <- keepStandardChromosomes(MCF7,pruning.mode="coarse")
genome(MCF7)<-'hg19'
export.bed(MCF7, paste(MCF7Path,"MCF-7.bed",sep=""))#It is then processed using panda

############################MCF_7 groHMM#########################################
#After processing using panda
MCF7Duplicated <- import(paste(MCF7Path,"MCF-7_duplicated.bed",sep=""), format="bed", genome = 'hg19')
MCF7Duplicated<-keepStandardChromosomes(MCF7Duplicated,pruning.mode="coarse")

mcols(MCF7Duplicated)$name<-NULL
mcols(MCF7Duplicated)$score<-NULL
show(MCF7Duplicated)

#Double the number of reads according to the paper
MCF7Duplicated2n<-sort(c(MCF7Duplicated,MCF7Duplicated))

seqlevels(MCF7Duplicated2n,pruning.mode="coarse")<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8"
                                  ,"chr9","chr10","chr11","chr12","chr13","chr14","chr15",
                                  "chr16","chr17","chr18","chr19","chr20","chr21","chrX")

show(MCF7Duplicated2n)

hmmMCF7<-detectTranscripts(MCF7Duplicated2n, LtProbB=-350, UTS=30, threshold=0.1)
hmmMCF7Transcripts<-hmmMCF7$transcripts
export.bed(hmmMCF7Transcripts, "data/MCF-7/groHMM_MCF7_2n.bed")#Save to file
export.bed(hmmMCF7Transcripts[strand(hmmMCF7Transcripts)=="+"], paste(MCF7Path,"groHMM_MCF7_2nPlus.bed",sep=""))#Save to file
export.bed(hmmMCF7Transcripts[strand(hmmMCF7Transcripts)=="-"], paste(MCF7Path,"groHMM_MCF7_2nMinus.bed",sep=""))#Save to file

####################### Wiggle file  visualization##############
write_wiggle(MCF7Duplicated2n,paste(MCF7Path,"MCF_7_2n",sep=""))

########################MCF7 Analysis###########################################
HUMAN_GENOME_LENGTH<-2859970000

show(hmmMCF7Transcripts)
plot_hist(ranges(hmmMCF7Transcripts),"Human groHMM predicted transcription region sizes")
extract_characteristics(hmmMCF7Transcripts,HUMAN_GENOME_LENGTH)

##################################DROSOPHILA Processing########################################
DROSOPHILIA_GENOME_LENGTH<-137688000

DROSOPHILA_PATH = "../data/Drosophila/"
DROSOPHILA_PREPROCESSED_PATH = paste(DROSOPHILA_PATH,"preprocessed/",sep="")

drosophilia_2h_Minus <- import(paste(DROSOPHILA_PREPROCESSED_PATH,"GSE41611_GROseq_wt_emb_2-2.5_both_Minus.bedgraph",sep=""), format="bed")
drosophilia_2h_Plus <- import(paste(DROSOPHILA_PREPROCESSED_PATH,"GSE41611_GROseq_wt_emb_2-2.5_both_Plus.bedgraph",sep=""), format="bed")
drosophilia_3h_Minus <- import(paste(DROSOPHILA_PREPROCESSED_PATH,"GSE41611_GROseq_wt_emb_3-3.5_both_Minus.bedgraph",sep=""), format="bed")
drosophilia_3h_Plus <- import(paste(DROSOPHILA_PREPROCESSED_PATH,"GSE41611_GROseq_wt_emb_3-3.5_both_Plus.bedgraph",sep=""), format="bed")

strand(drosophilia_2h_Minus)<-"-"
strand(drosophilia_3h_Minus)<-"-"
strand(drosophilia_2h_Plus)<-"+"
strand(drosophilia_3h_Plus)<-"+"

drosophilia_2h_Minus<-preprocess_drosophila(drosophilia_2h_Minus)
drosophilia_2h_Plus<-preprocess_drosophila(drosophilia_2h_Plus)
drosophilia_3h_Minus<-preprocess_drosophila(drosophilia_3h_Minus)
drosophilia_3h_Plus<-preprocess_drosophila(drosophilia_3h_Plus)

drosophila <- mergeGRangesData(drosophilia_2h_Minus,drosophilia_2h_Plus,
                               drosophilia_3h_Plus,drosophilia_3h_Minus,
                              field = "score",exact_overlaps = FALSE)

score(drosophila)<-abs(score(drosophila))

export.bed(drosophila, paste(DROSOPHILA_PATH,"drosophila.bed",sep=""))
#Must now process it to duplicate it using panda and python

############################Drosophila groHMM#########################################

drosophila_duplicated<- import(paste(DROSOPHILA_PATH,"drosophila_duplicated.bed",sep=""), format="bed")

show(drosophila_duplicated)
mcols(drosophila_duplicated)$name<-NULL
mcols(drosophila_duplicated)$score<-NULL

hmmDrosophila <- detectTranscripts(drosophila_duplicated, LtProbB=-50, UTS=50, threshold=0.09)
hmmDrosophilaTranscript<-hmmDrosophila$transcripts

export.bed(hmmDrosophilaTranscript, "groHMM_Drosophilia.bed")
# print(hmmDrosophilia$emisParams)
# print(hmmDrosophilia$transParams)

export.bed(hmmDrosophilaTranscript, paste(DROSOPHILA_PATH,"groHMM_drosophila.bed",sep=""))#Save to file
export.bed(hmmDrosophilaTranscript[strand(hmmDrosophilaTranscript)=="+"], paste(DROSOPHILA_PATH,"groHMM_drosophila_Plus.bed",sep=""))#Save to file
export.bed(hmmDrosophilaTranscript[strand(hmmDrosophilaTranscript)=="-"], paste(DROSOPHILA_PATH,"groHMM_drosophila_Minus.bed",sep=""))#Save to file

hmmDrosophilaRanges<- ranges(hmmDrosophilaTranscript)
extract_characteristics(hmmDrosophilaRanges,DROSOPHILIA_GENOME_LENGTH)
plot_hist(hmmDrosophilaRanges,"Drosophilia groHMM predicted transcription region sizes")


#######################Wiggle file visualization################################
write_wiggle(drosophila_duplicated,paste(DROSOPHILA_PATH,"drosophila_duplicated",sep=""))
##################################Arabidopsis processing########################################
ARABIDOPSIS_GENOME_LENGTH<-119669000

ARABIDOPSIS_PATH<-"../data/Arabidopsis/"

arabidopsis_N_Rep1<-import(paste(ARABIDOPSIS_PATH,"raw/GSM3682879_N_GRO_rep1.bedgraph",sep=""))
arabidopsis_N_Rep2<-import(paste(ARABIDOPSIS_PATH,"raw/GSM3682880_N_GRO_rep2.bedgraph",sep=""))
arabidopsis_T_Rep1<-import(paste(ARABIDOPSIS_PATH,"raw/GSM3682881_T_GRO_rep1.bedgraph",sep=""))
arabidopsis_T_Rep2<-import(paste(ARABIDOPSIS_PATH,"raw/GSM3682882_T_GRO_rep2.bedgraph",sep=""))


score(arabidopsis_N_Rep1)<-round(score(arabidopsis_N_Rep1)/161.05)
score(arabidopsis_N_Rep2)<-round(score(arabidopsis_N_Rep2)/156.111)
score(arabidopsis_T_Rep1)<-round(score(arabidopsis_T_Rep1)/151.98)
score(arabidopsis_T_Rep2)<-round(score(arabidopsis_T_Rep2)/152.95)

arabidopsis_N_Rep1<-add_strand(arabidopsis_N_Rep1)
arabidopsis_N_Rep2<-add_strand(arabidopsis_N_Rep2)
arabidopsis_T_Rep1<-add_strand(arabidopsis_T_Rep1)
arabidopsis_T_Rep2<-add_strand(arabidopsis_T_Rep2)

arabidopsis<-mergeGRangesData(arabidopsis_N_Rep1,arabidopsis_N_Rep2,
                              arabidopsis_T_Rep1,arabidopsis_T_Rep2,
                              field = "score",exact_overlaps = FALSE)

export.bed(arabidopsis, paste(ARABIDOPSIS_PATH,"arabidopsis.bed",sep=""))
#Must now process it to duplicate it using panda and python

############################Arabidopsis groHMM#########################################

arabidopsis_duplicated<- import(paste(ARABIDOPSIS_PATH,"arabidopsis_duplicated.bed",sep=""), format="bed")
mcols(arabidopsis_duplicated)$name<-NULL
mcols(arabidopsis_duplicated)$score<-NULL

show(arabidopsis_duplicated)

hmmArabidopsis_drosophila_param <- detectTranscripts(arabidopsis_duplicated, LtProbB=-50, UTS=50,threshold=0.09)
hmmArabidopsis_human_param <- detectTranscripts(arabidopsis_duplicated, LtProbB=-350, UTS=30,threshold=0.1)

export_hmm(hmmArabidopsis_drosophila_param,paste(ARABIDOPSIS_PATH,"groHMM_arabidopsis_50_50",sep=""))
export_hmm(hmmArabidopsis_human_param,paste(ARABIDOPSIS_PATH,"groHMM_arabidopsis_350_30",sep=""))

hmmArabidopsisTranscript<-hmmArabidopsis$transcripts
export.bed(hmmArabidopsisTranscript, "groHMM_Arabidopsis.bed")

#hmmArabidopsisTranscript<-import("groHMM_Arabidopsis.bed")
hmmArabidopsisRanges<- ranges(hmmArabidopsisTranscript)

plot_hist(hmmArabidopsisRanges,"Arabidopsis groHMM predicted transcription region sizes")
extract_characteristics(hmmArabidopsisRanges,ARABIDOPSIS_GENOME_LENGTH)

#######################Wiggle file visualization################################
write_wiggle(arabidopsis_duplicated,paste(ARABIDOPSIS_PATH,"arabidopsis_duplicated",sep=""))
test<-arabidopsis_duplicated
show(test)
show(seqinfo(test))
isCircular(test)<-rep(FALSE,times=5)
genome(test)<-rep("TAIR10.1",times=5)
seqlengths(test)<-rep(0,times=5)
show(test)
wigToBigWig(paste(ARABIDOPSIS_PATH,"arabidopsis_duplicated_Plus",sep=""))

