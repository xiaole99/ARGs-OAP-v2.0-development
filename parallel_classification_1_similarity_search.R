getwd()
################# give cut off ################################################
cutoff <- 80
cutoff <- 90
cutoff <- 70

################# load files  #################################################################
coredatabase<- read.delim(file="str_subtype_type_for_classification.txt",sep="\t",header=T,stringsAsFactors = F)
blast6out<- as.data.frame(read.delim("blast6out_nr.txt",sep="\t",header = TRUE,stringsAsFactors = F))
classify_keyword<- as.data.frame(read.delim("90-100_type_subtype.txt",sep="\t",header = F,stringsAsFactors = F))

################# just run #############################################################################
### similarity search
blast6out<- blast6out[,c(1,2,3)]
colnames(blast6out)<- c("core_name","ncbi_name","identity")
blast6out <-blast6out[order(blast6out$ncbi_name,-blast6out$identity),]
blast6out <- subset(blast6out,identity <100 )
# best hit for each ncbi_id
dereplication <- blast6out[!duplicated(blast6out$ncbi_name),]
# filter by cut off
dereplication_filter <- subset(dereplication,identity >= cutoff)
classify_genelevel <- merge(dereplication_filter,coredatabase,by.x="core_name",by.y="sequenceid",all=F)
classify_genelevel <- classify_genelevel[,c(2,1,3,4)]
colnames(classify_genelevel)[4]<-"genelevel"
classify_genelevel$mark <- 1:nrow(classify_genelevel)
rm(dereplication_filter)
### keywords search
classify_keyword$keywordssearch <- paste(classify_keyword$V2, classify_keyword$V3,sep="__")
classify_keyword <- classify_keyword[,c(1,4)]
### merge
classify_both <- merge(classify_genelevel,classify_keyword,by.x="ncbi_name",by.y="V1",all.y=T)
tmp <- classify_both[is.na(classify_both$genelevel),]
tmptmp<-merge(tmp,coredatabase,by.x="ncbi_name",by.y="sequenceid",all=F)
tmp <- tmptmp[,c(1,2,3,7,5,6)]
rm(tmptmp)
colnames(tmp)[4]<-"genelevel"
tmp$mark<- paste("c",1:nrow(tmp),sep="")
tmp2 <- classify_both[!(is.na(classify_both$genelevel)),]
classify_both <- rbind(tmp,tmp2)
rm(tmp)
rm(tmp2)
classify_both <-classify_both[order(classify_both$ncbi_name,-classify_both$identity),]
# check if keywordsearch column has NA
classify_both[is.na(classify_both$keywordssearch),]
### find match
matchbesthit <-classify_both[which(classify_both$keywordssearch== classify_both$genelevel),]
matchbesthit <-matchbesthit[,c(1,4)]
write.table(matchbesthit,file="matchbesthit.txt",quote=F,sep="\t",col.names=F,row.names=F)

