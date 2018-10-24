rm(list=ls(all=TRUE))

###http://davetang.org/muse/2011/07/27/edger-vs-deseq-using-pnas_expression-txt/
###http://cgrlucb.wikispaces.com/edgeRWikiVersion
###source("http://bioconductor.org/biocLite.R")
###biocLite("edgeR")
#specify input directory  
#illumina rest shrimp
# setwd("/home/jan/Data/projects/mirna_expression/miRNA_marek/CLL/new_results/res_detailed_v0.23/shrimp/illumina_rest")
#illumina rest miranalyzer
#setwd("/home/jan/Projects/mirna_expression/miRNA_marek/CLL/new_results/res_detailed_v0.23/miranalyzer/illumina_rest")
#illumina main shrimp
# setwd("/home/jan/Data/projects/mirna_expression/miRNA_marek/CLL/new_results/res_detailed_v0.23/shrimp/illumina")
#ilumina main mirnalyzer
setwd("/home/jan/Projects/mirna_expression/miRNA_marek/CLL/new_results/res_detailed_v0.23/miranalyzer/illumina")
#specify output directory
datasetDir<-"/home/jan/Desktop/mirna_expression_illumina_miranalyzer_main"
#specify name of a table withou counts in the input directory
#illumina rest shrimp
# dataset<-"RCtableLimSub_rest_mature.txt"
#illumina main shrimp
# dataset<-"RCtableLimSub_illumina_mature.txt"
#illumina rest mirnalyzer
#dataset<-"RCtableLimSub_rest_mature_unique.txt"
#illumina main miranalyzer
dataset<-"RCtableLimSub_mature_unique.txt"
#load library and calculate statistics
############TO DO COMMENTS#################################################################################################################
library(edgeR)
RCtableEdgeR<-read.table(dataset,sep=";",dec=".",header=T)

dir.create(datasetDir,showWarnings=T)
setwd(datasetDir)

#rownames(RCtableEdgeR)<-RCtableEdgeR[,1]
#RCtableEdgeR<-RCtableEdgeR[,-1]

sink("library_sizes.txt")
colSums(RCtableEdgeR)
sink()

### illumina main
RCtableEdgeR<-RCtableEdgeR[c("lane5IV1","lane5IV5","lane5IV7","lane5IV9","lane5IV11","lane5IV2","lane5IV6","lane5IV8","lane5IV10","lane5IV12")]
### illumina rest
#RCtableEdgeR<-RCtableEdgeR[c("lane5IV3","lane5IV4","lane5IV17","lane5IV18","lane5IV19","lane5IV20","lane5IV21","lane5IV22","lane5IV23","lane5IV24")]

group<-c(rep("Bf",5),rep("Af",5))
d<-DGEList(counts=RCtableEdgeR,group=group)
#d<-DGEList(counts=as.matrix(RCtableEdgeR),group=group)
dim(d)

# filter uninformative genes
#m <- 1e6 * t(t(d$counts) / d$samples$lib.size)
#ridx <- rowSums(m > 1) >= 2
#d <- d[ridx,]

sink("genes_RClower_than10.txt")
nrow(RCtableEdgeR[rowSums(RCtableEdgeR)<10,])
RCtableEdgeR[rowSums(RCtableEdgeR)<10,]
sink()

sink("genes_RClower_than20.txt")
nrow(RCtableEdgeR[rowSums(RCtableEdgeR)<20,])
RCtableEdgeR[rowSums(RCtableEdgeR)<20,]
sink()

sink("samples_control.txt")
d$samples
sink()

d<-calcNormFactors(d)
d

png(file="MDS.png",height=600,width=600)
plotMDS(d,main="MDSPlot",xlim=c(-1,1),labels=c("Bf1","Bf3","Bf4","Bf5","Bf6","Af1","Af3","Af4","Af5","Af6"),xlab="LeadinglogFCdim1",ylab="LeadinglogFCdim2")
dev.off()

#treatment<-factor(group)
treatment<-factor(c(rep("Bf",5),rep("Af",5)),levels=c("Bf","Af"))
patient<-factor(c(1,3,4,5,6,1,3,4,5,6))
design <- model.matrix(~patient+treatment)
rownames(design) <- colnames(d)

sink("design_controll")
data.frame(Sample=colnames(d),patient,treatment)
sink()


contrasts(patient)<-contr.sum(5)
contrasts(treatment)<-contr.sum(2)

d <- estimateGLMCommonDisp(d,design)

d <- estimateGLMTrendedDisp(d,design)

d <- estimateGLMTagwiseDisp(d,design)

fit_com <- glmFit(d, design, dispersion=d$common.dispersion)
lrt_com <- glmLRT(fit_com)
topTags(lrt_com,adjust.method="BH",sort.by="PValue")


fit_trd <- glmFit(d, design, dispersion=d$trended.dispersion)
lrt_trd <- glmLRT(fit_trd)
topTags(lrt_trd,adjust.method="BH",sort.by="PValue")

fit_tgw <- glmFit(d, design, dispersion=d$tagwise.dispersion)
lrt_tgw <- glmLRT(fit_tgw)
topTags(lrt_tgw,adjust.method="BH",sort.by="PValue")

o_com <- order(lrt_com$table$PValue)

sink("counts_per_million_byPValue_common.txt")
cpm(d)[o_com,]
#cpm(d)[o_com[1:10],]
sink()

summary(de_com <- decideTestsDGE(lrt_com))

detags_com <- rownames(d)[as.logical(de_com)]

#png(file="logFC_logCPM_common.png",height=600,width=800)
#plotSmear(lrt_com, de.tags=detags_com)
#abline(h=c(-1, 1), col="blue")
#dev.off()

FDR_com <- p.adjust(lrt_com$table$PValue, method="BH")
sum(FDR_com < 0.1)

com=summary(decideTestsDGE(lrt_com,adjust.method="BH",p.value=0.05))
com_total=com[1]+com[3]

sink("topTags_com.txt")
topTags(lrt_com,n=com_total,adjust.method="BH")
sink()


o_trd <- order(lrt_trd$table$PValue)

sink("counts_per_million_byPValue_trended.txt")
#cpm(d)[o_trd[1:10],]
cpm(d)[o_trd,]
sink()

summary(de_trd <- decideTestsDGE(lrt_trd))

detags_trd <- rownames(d)[as.logical(de_trd)]

#png(file="logFC_logCPM_trended.png",height=600,width=800)
#plotSmear(lrt_trd, de.tags=detags_trd)
#abline(h=c(-1, 1), col="blue")
#dev.off()

FDR_trd <- p.adjust(lrt_trd$table$PValue, method="BH")
sum(FDR_trd < 0.1)

trd=summary(decideTestsDGE(lrt_trd,adjust.method="BH",p.value=0.05))
trd_total=trd[1]+trd[3]


sink("topTags_trd")
topTags(lrt_trd,n=trd_total,adjust.method="BH")
sink()


o_tgw <- order(lrt_tgw$table$PValue)

sink("counts_per_million_byPValue_tagwise.txt")
#cpm(d)[o_tgw[1:10],]
cpm(d)[o_tgw,]
sink()

summary(de_tgw <- decideTestsDGE(lrt_tgw))

detags_tgw <- rownames(d)[as.logical(de_tgw)]

#png(file="logFC_logCPM_tagwise.png",height=600,width=800)
#plotSmear(lrt_tgw, de.tags=detags_tgw)
#abline(h=c(-1, 1), col="blue")
#dev.off()


FDR_tgw <- p.adjust(lrt_tgw$table$PValue, method="BH")
sum(FDR_tgw < 0.1)

tgw=summary(decideTestsDGE(lrt_tgw,adjust.method="BH",p.value=0.05))
tgw_total=tgw[1]+tgw[3]

sink("topTags_tgw")
topTags(lrt_tgw,n=tgw_total,adjust.method="BH")
sink()

#isDE <- as.logical(de)
#DEnames <- rownames(d)[isDE]

png(file="Library_size_mill.png",height=600,width=800)
barplot(d$samples$lib.size*1e-6, names=1:10, ylab="Library size (millions)")
dev.off()

resultsTbl.com<-topTags(lrt_com,n=nrow(lrt_com$table))$table
resultsTbl.trd<-topTags(lrt_trd,n=nrow(lrt_trd$table))$table
resultsTbl.tgw<-topTags(lrt_tgw,n=nrow(lrt_tgw$table))$table

wh.rows.com <- match( rownames( resultsTbl.com ) , rownames( d$counts ) )
wh.rows.trd <- match( rownames( resultsTbl.trd ) , rownames( d$counts ) )
wh.rows.tgw <- match( rownames( resultsTbl.tgw ) , rownames( d$counts ) )

combResults.com<-cbind(resultsTbl.com,"Com.Disp"=d$common.dispersion[wh.rows.com],"UpDown.Com"=decideTestsDGE(lrt_com,adjust.method="BH",p.value=0.05)[wh.rows.com],d$counts[wh.rows.com,])
combResults.trd<-cbind(resultsTbl.trd,"trd.Disp"=d$trended.dispersion[wh.rows.trd],"UpDown.trd"=decideTestsDGE(lrt_trd,adjust.method="BH",p.value=0.05)[wh.rows.trd],d$counts[wh.rows.trd,])
combResults.tgw<-cbind(resultsTbl.tgw,"tgw.Disp"=d$tagwise.dispersion[wh.rows.tgw],"UpDown.tgw"=decideTestsDGE(lrt_tgw,adjust.method="BH",p.value=0.05)[wh.rows.tgw],d$counts[wh.rows.tgw,])


write.table(combResults.com,file="commonResultsDetailed.csv",sep=",",row.names=TRUE)
write.table(combResults.trd,file="trendedResultsDetailed.csv",sep=",",row.names=TRUE)
write.table(combResults.tgw,file="tagwiseResultsDetailed.csv",sep=",",row.names=TRUE)


png(file="mean_variation_tgw.png",height=600,width=800)
mv<-plotMeanVar(d,show.raw.vars=TRUE,show.tagwise.vars=TRUE,NBline=TRUE,main="Mean Variation Tagwise")
dev.off()


de.genes.com<-rownames(resultsTbl.com)[resultsTbl.com$FDR<=0.1]
de.genes.trd<-rownames(resultsTbl.trd)[resultsTbl.trd$FDR<=0.1]
de.genes.tgw<-rownames(resultsTbl.tgw)[resultsTbl.tgw$FDR<=0.1]


png(file="MA_top_tags.png",height=1200,width=1000)
par(mfrow=c(3,1))

plotSmear(d,de.tags=de.genes.com,main="Common",pair=c("Bf","Af"),cex=.35,xlab="Log Concentration",ylab="LogFold-Change")
abline(h=c(-1,1),col="dodgerblue")

plotSmear(d,de.tags=de.genes.trd,main="Tagwise",pair=c("Bf","Af"),cex=.35,xlab="Log Concentration",ylab="LogFold-Change")
abline(h=c(-1,1),col="dodgerblue")

plotSmear(d,de.tags=de.genes.tgw,main="Poisson",pair=c("Bf","Af"),cex=.35,xlab="Log Concentration",ylab="LogFold-Change")
abline(h=c(-1,1),col="dodgerblue")

dev.off()

png(file="tgw_disp_vs_log2_conc_common.png",height=600,width=800)
plotBCV(d,cex=0.4)
abline(h=d$common.dispersion,col="dodgerblue",lwd=3)
dev.off()

sink("logFC_1_-1_common_trended_tagwise.txt")
topNUMcom<-topTags(lrt_com,n=com_total,adjust.method="BH")
sum(topNUMcom$table$logFC>1)
sum(topNUMcom$table$logFC<-1)

topNUMtrd<-topTags(lrt_trd,n=trd_total,adjust.method="BH")
sum(topNUMtrd$table$logFC>1)
sum(topNUMtrd$table$logFC<-1)

topNUMtgw<-topTags(lrt_tgw,n=tgw_total,adjust.method="BH")
sum(topNUMtgw$table$logFC>1)
sum(topNUMtgw$table$logFC<-1)
sink()

detagsNUMcom<-rownames(topNUMcom$table)
png(file="MA_plot_com_disp.png",height=600,width=800)
plotSmear(lrt_com,de.tags=detagsNUMcom,main="FC plot using common dispersion")
abline(h=c(-1,1),col="dodgerblue")
dev.off()

detagsNUMtrd<-rownames(topNUMtrd$table)
png(file="MA_plot_trd_disp.png",height=600,width=800)
plotSmear(lrt_trd,de.tags=detagsNUMtrd,main="FC plot using trended dispersion")
abline(h=c(-1,1),col="dodgerblue")
dev.off()

detagsNUMtgw<-rownames(topNUMtgw$table)
png(file="MA_plot_tgw_disp.png",height=600,width=800)
plotSmear(lrt_tgw,de.tags=detagsNUMtgw,main="FC plot using tagwise dispersion")
abline(h=c(-1,1),col="dodgerblue")
dev.off()

#save.image(file="marek_edgeR_paired_img.R")
#savehistory(file="marek_edgeR_paired_his.R")

##################NEW################################
##################NEW################################
##################NEW################################

#normalized counts with tagwise dispersion estimate
counts.per.m <- cpm(d,normalized.lib.sizes=T,log=F)
counts.per.m.log <- cpm(d,normalized.lib.sizes=T,log=T)#for logarithmic transformation

#specify selected genes
genes_select1<-c("hsa-miR-34a-5p","hsa-miR-155-3p","hsa-miR-29b-1-5p","hsa-miR-4521","hsa-miR-1248","hsa-miR-1246")
genes_select<-NULL
for (name in 1:length(genes_select1)){
if ((sum((rownames(RCtableEdgeR)==genes_select1[name])==T))==1){
  genes_select<-rbind(genes_select,genes_select1[name])
}
}

genes_select<-as.vector(genes_select)

counts.per.m <- as.data.frame(counts.per.m)
#create sub group of only selected genes
# counts.per.m.sub<-counts.per.m[genes_select,]

counts.per.m <- data.matrix(counts.per.m)
# counts.per.m.sub <- data.matrix(counts.per.m.sub)
#select only those rows and data which are complete
#counts.per.m.sub<-counts.per.m.sub[complete.cases(counts.per.m.sub),]
#see quantiles from counts.per.m
quantile(rowSums(counts.per.m))
quantile(rowSums(counts.per.m[genes_select,]))
################################################
#load libraries
library(gplots)#heatmap
library(RColorBrewer)#color pallete

#heatmap.2(counts.per.m.sub,na.rm=T)
#####HEAT MAP START#######
#####WORKS BUT HAVE TO APPLIED TO NORMALIZED VALUES####
###install.packages("gplots")

# png(file="test.png",height=4000,width=6000)
# heatmap.2(as.matrix(counts.per.m), col=redgreen(75), scale="row", key=T, keysize=1,
#           density.info="none", trace="none",cexCol=1.5,cexRow=1,dendrogram = "column",margins=c(10,10))#, labRow=NA
# dev.off()

### do not remove "scale="row"" in this calculation - the result is non-sense picture
#
#heamap for all mirnas - log2cpm (default in edgeR calculation)
png(file="heatmap_logcpm_all_mirnas.png",height=3000,width=4500)
heatmap.2(counts.per.m.log,col=brewer.pal(11,"RdBu"),cexCol=2,cexRow=0.8, scale="row", trace="none", dendrogram="column",keysize=0.8,margins=c(12,12), na.rm = TRUE)
dev.off()
#heamap for selected mirnas - log2cpm (default in edgeR calculation)
png(file="heatmap_logcpm_selected_mirnas.png",height=1000,width=1500)
heatmap.2(counts.per.m.log[genes_select,],col=brewer.pal(11,"RdBu"),cexCol=1.5,cexRow=1.5, scale="row", trace="none", dendrogram="column",keysize=1,margins=c(12,12), na.rm = TRUE)
dev.off()
#heamap for selected mirnas - cpm
png(file="heatmap_cpm_all_mirnas.png",height=3000,width=4500)
heatmap.2(counts.per.m,col=brewer.pal(11,"RdBu"),cexCol=2,cexRow=0.8, scale="row", trace="none", dendrogram="column",keysize=0.8,margins=c(12,12), na.rm = TRUE)
dev.off()
#heamap for all mirnas - cpm
png(file="heatmap_cpm_selected_mirnas.png",height=1000,width=1500)
heatmap.2(counts.per.m[genes_select,],col=brewer.pal(11,"RdBu"),cexCol=1.5,cexRow=1.5, scale="row", trace="none", dendrogram="column",keysize=1,margins=c(12,12), na.rm = TRUE)
dev.off()
#####HEAT MAP END#######
#################################################
#####PERCENTAGE START#######
#load libraries
library(reshape2)
library(plyr)

#extract normalized counts
tmp_cpm<-NULL
tmp_cpm<-cpm(d,normalized.lib.sizes=T,log=F)
#calculate frequencies for all samples
#create empty objects
counts.normalized<-NULL;i<-NULL
#create empty dataframe
counts.normalized<-as.data.frame(matrix(, nrow = nrow(tmp_cpm), ncol = ncol(tmp_cpm)))#number of rows=number of rows in the input ad number of columns=number of column in input
for (i in 1:ncol(tmp_cpm)){#calculate frequencies of all samples
  counts.normalized[,i]<-(tmp_cpm[,i]/sum(tmp_cpm[,i]))*100#results are in percentages
}
#rename columns and rows based on the input object
colnames(counts.normalized)<-colnames(tmp_cpm)
rownames(counts.normalized)<-rownames(tmp_cpm)
################################################################
#extract information about top 20 expressed miRNA in all samples and sort them from the most expressed based on normalized counts
#create empty object to do it
i<-NULL;

for (i in 1:ncol(tmp_cpm))#calculate it for all samples
{
  a<-NULL;tmp2<-NULL;newdata<-NULL;selection.data<-NULL;rest<-NULL;sample<-NULL;melted.counts<-NULL
  tmp2<-tmp_cpm[,i]#copy current sample
  tmp2<-as.data.frame(tmp2)
  tmp2[,2]<-rownames(tmp_cpm)#copy rownames to second column
  colnames(tmp2)<-c("norm_expr","mirna")#rename colnames
  newdata<-tmp2[with(tmp2, order(-norm_expr)), ]#sort them by "norm_expr"=frequencies
  rownames(newdata)<-newdata[,2]#copy rownames from second column
  sample<-colnames(tmp_cpm)[i]#save current sample name
  #select top 20 most expressed mirnas
  selection.data<-newdata[1:20,]
  selection.data$mirna<-NULL#remove mirna column
  selection.data$counts<-d$counts[rownames(selection.data),i]#add raw counts
  #write the table  
  write.table(selection.data,file=paste(sample,"_norm_expr_top20.csv",sep=""),sep=";",col.names=NA)
}

################################################################
#load libraries
library(ggplot2)
#extract sample one by one and select top frequent miRNAs
#the same thing as previously but based on frequency
i<-NULL;
selection.data.total<-as.data.frame(matrix(nrow=40,ncol=0))
selection.data.total.tmp<-as.data.frame(matrix(nrow=21,ncol=1))
selection.data.total.names<-as.data.frame(matrix(ncol=1,nrow=0))
selection.data.total.names.tmp<-matrix(ncol=1,nrow=20)

for (i in 1:ncol(counts.normalized))#calculate it for all 
{
a<-NULL;tmp2<-NULL;newdata<-NULL;selection.data<-NULL;rest<-NULL;sample<-NULL;melted.counts<-NULL
tmp2<-counts.normalized[,i]
tmp2<-as.data.frame(tmp2)
tmp2[,2]<-rownames(tmp_cpm)
colnames(tmp2)<-c("counts","mirna")
newdata<-tmp2[with(tmp2, order(-counts)), ]
rownames(newdata)<-newdata[,2]
#select only chosen mirnas
selection.data3<-NULL
selection.data3<-newdata[genes_select,]
colnames(selection.data3)<-c("frequency","mirna")
counts_tmp12<-NULL
counts_tmp12<-as.data.frame(tmp_cpm[genes_select,][,i])
colnames(counts_tmp12)<-"norm_count"
counts_tmp12$mirna<-rownames(counts_tmp12)
selection.data3<-merge(selection.data3,counts_tmp12)

counts_tmp13<-NULL
counts_tmp13<-as.data.frame(d$counts[genes_select,][,i])
colnames(counts_tmp13)<-"raw_count"
counts_tmp13$mirna<-rownames(counts_tmp13)
selection.data3<-merge(selection.data3,counts_tmp13)

sample<-colnames(tmp_cpm)[i]

write.table(selection.data3[order(as.numeric(selection.data3$frequency),decreasing=T),],file=paste(sample,"_selected_mirnas.csv",sep=""),sep=";",
            dec=".",col.names=NA)

#select top 20 most expressed mirnas
selection.data<-newdata[1:20,]
# selection.data<-newdata[newdata[,1]>1,]#select only miRNAs with more that 1% of total "expression"

rest<-100-sum(selection.data[,1])

selection.data[nrow(selection.data)+1,1]<-rest
selection.data[nrow(selection.data),2]<-"rest"
rownames(selection.data)[nrow(selection.data)]<-"rest"

sample<-colnames(counts.normalized)[i]

selection.data<-selection.data[with(selection.data, sort(mirna,decreasing=T)), ]

selection.data2<-selection.data
colnames(selection.data2)<-c("frequency","mirna")
names_tmp<-NULL
names_tmp<-rownames(selection.data2)
names_tmp<-names_tmp[names_tmp != "rest"]
counts_tmp10<-NULL
counts_tmp10<-as.data.frame(tmp_cpm[names_tmp,][,i])
colnames(counts_tmp10)<-"norm_count"
counts_tmp10$mirna<-rownames(counts_tmp10)
selection.data2<-merge(selection.data2,counts_tmp10)

counts_tmp11<-NULL
counts_tmp11<-as.data.frame(d$counts[names_tmp,][,i])
colnames(counts_tmp11)<-"raw_count"
counts_tmp11$mirna<-rownames(counts_tmp11)
selection.data2<-merge(selection.data2,counts_tmp11)

selection.data2[21,]<-c("rest",selection.data["rest",1],(sum(tmp_cpm[,i])-sum(as.numeric(selection.data2$norm_count[1:20]),na.rm=T)),
                               (sum(d$count[,i])-sum(as.numeric(selection.data2$raw_count[1:20]),na.rm=T)))



write.table(selection.data2[order(as.numeric(selection.data2$frequency),decreasing=T),],file=paste(sample,"_top20_selected.csv",sep=""),sep=";",
            dec=".",col.names=NA)

colnames(selection.data)<-c("counts","mirna")

melted.counts<-melt(selection.data)
melted.counts = rename(melted.counts, c(mirna="mirna", variable=sample))

Palette1 <- rainbow(nrow(selection.data))

a <- ggplot(melted.counts, aes(x = sample, y = value, fill = mirna)) 
a <- a + geom_bar(stat = "identity",  position = "stack", colour="black")
a <- a + scale_fill_manual(values=Palette1)

png(file= paste(sample,"logfreq.png",sep=""),height=800,width=600)
print(a)
dev.off()
}

##################TEST
average_cpm<-NULL
average_cpm<-rowMeans(tmp_cpm,na.rm=T)
average_cpm<-as.data.frame(average_cpm)

average.normalized<-NULL;
average.normalized<-as.data.frame(matrix(, nrow = nrow(average_cpm), ncol = 2))
rownames(average.normalized)<-rownames(average_cpm)
average.normalized[,1]<-(average_cpm[,1]/sum(average_cpm[,1]))*100
average.normalized[,2]<-rownames(average.normalized)

colnames(average.normalized)<-c("counts","mirna")
average.normalized<-average.normalized[with(average.normalized, order(-counts)), ]

selection.data.average<-average.normalized[1:20,]

rest2<-100-sum(average.normalized[,1])

selection.data.average[nrow(selection.data.average)+1,1]<-rest2
selection.data.average[nrow(selection.data.average),2]<-"rest"
rownames(selection.data.average)[nrow(selection.data.average)]<-"rest"

selection.data.average<-selection.data.average[with(selection.data.average, sort(mirna,decreasing=T)), ]

melted.counts.average<-melt(selection.data.average)
sample<-"averaged"
melted.counts.average = rename(melted.counts.average, c(mirna="mirna", variable=sample))

Palette1 <- rainbow(nrow(selection.data.average))


a <- ggplot(melted.counts.average, aes(x = sample, y = value, fill = mirna)) 
a <- a + geom_bar(stat = "identity",  position = "stack", colour="black")
a <- a + scale_fill_manual(values=Palette1)

png(file= paste("averaged_freq.png",sep=""),height=800,width=600)
print(a)
dev.off()
##################TEST



#####PERCENTAGE END#######
##########################
######EXTRACT LOGFC FOR SELECTED GENES##################
logFC_selection<-combResults.tgw[genes_select,1]
logFC_selection<-as.data.frame(logFC_selection)
rownames(logFC_selection)<-genes_select
colnames(logFC_selection)<-"logFC"

barplot_tmp<-logFC_selection[,1]

all_tmp<-NULL
all_tmp<-as.data.frame(combResults.tgw[,1])
colnames(all_tmp)<-"logFC"
rownames(all_tmp)<-rownames(combResults.tgw)

write.table(all_tmp,file="logFC_all.csv",sep=";",dec=".", col.names=NA)

write.table(logFC_selection,file="logFC_selected_mirnas.csv",sep=";",dec=".", col.names=NA)

png(file= paste("logFC.png",sep=""),height=600,width=800)
barplot(barplot_tmp, names.arg=genes_select, main="logFC between groups", ylim=c(-2,2))
abline(h=0)
box()
dev.off()

logFC_table<-NULL;

###replace zero with medians
###better is to handle values without any replacement
# for (i in 1:nrow(counts.per.m.sub))
# {
# counts.per.m.sub[i,(counts.per.m.sub[i,]==0)]<-median(counts.per.m.sub[i,])
# }

###CHECK WHETHER IT IS NOT BETTER TO CALCULATE log2(counts.per.m.sub...)/log2(counts.per.m.sub...)
####Illumina main
logFC_table<-log2(counts.per.m[,"lane5IV2"]/counts.per.m[,"lane5IV1"])
logFC_table<-as.data.frame(logFC_table)
colnames(logFC_table)<-"lane5IV2/lane5IV1"
logFC_table$"lane5IV6/lane5IV5"<-log2(counts.per.m[,"lane5IV6"]/counts.per.m[,"lane5IV5"])
logFC_table$"lane5IV8/lane5IV7"<-log2(counts.per.m[,"lane5IV8"]/counts.per.m[,"lane5IV7"])
logFC_table$"lane5IV10/lane5IV9"<-log2(counts.per.m[,"lane5IV10"]/counts.per.m[,"lane5IV9"])
logFC_table$"lane5IV12/lane5IV11"<-log2(counts.per.m[,"lane5IV12"]/counts.per.m[,"lane5IV11"])

###Illumina rest
# logFC_table<-log2(counts.per.m[,"lane5IV4"]/counts.per.m[,"lane5IV3"])
# logFC_table<-as.data.frame(logFC_table)
# colnames(logFC_table)<-"lane5IV4/lane5IV3"
# logFC_table$"lane5IV18/lane5IV17"<-log2(counts.per.m[,"lane5IV18"]/counts.per.m[,"lane5IV17"])
# logFC_table$"lane5IV20/lane5IV19"<-log2(counts.per.m[,"lane5IV20"]/counts.per.m[,"lane5IV19"])
# logFC_table$"lane5IV22/lane5IV21"<-log2(counts.per.m[,"lane5IV22"]/counts.per.m[,"lane5IV21"])
# logFC_table$"lane5IV24/lane5IV23"<-log2(counts.per.m[,"lane5IV24"]/counts.per.m[,"lane5IV23"])





is.na(logFC_table) <- sapply(logFC_table, is.infinite) #would replace infinite values with NAs
                                                         #this can happen if we have zero expression at some genes and we do not replace them or handle them
                                                        #testing showed this is the closest results (=not replace zeros)
barplot_tmp2<-as.data.frame(apply(logFC_table, 1, mean, na.rm = TRUE))
logFC_table$mean<-barplot_tmp2[,1]
write.table(logFC_table,file="manual_calculated_logFC.csv",sep=";",dec=".",col.names=NA)
write.table(logFC_table[genes_select,],file="manual_calculated_logFC_selected_mirnas.csv",sep=";",dec=".",col.names=NA)


png(file= paste("logFC_manual.png",sep=""),height=600,width=800)
barplot(barplot_tmp2[genes_select,], names.arg=genes_select, main="logFC between groups", ylim=c(-2,2))
abline(h=0)
box()
dev.off()

logFC_table_select<-NULL
logFC_table_select<-logFC_table[genes_select,]

stand_dev<-NULL
stand_dev<-sd(logFC_table_select[1,],na.rm=T)
stand_dev<-as.data.frame(stand_dev)

for (i in 2:nrow(logFC_table_select))
{
stand_dev[i,]<-sd(logFC_table_select[i,], na.rm=T)
}
rownames(stand_dev)<-genes_select
############################

logFC_table_select<-as.matrix(logFC_table_select)

rowMeans(logFC_table_select,na.rm=T)

means<-rowMeans(logFC_table_select,na.rm=T)

png(file= paste("logFC_manual_withSD.png",sep=""),height=600,width=800)
mp <- barplot(means, axes=FALSE, axisnames=FALSE, ylim=c(-3,3),main="logFC between groups")

axis(1, labels=genes_select, at = mp)

axis(2, at=seq(-3 , 3, by=0.5))

stDevs <- as.matrix(stand_dev)

segments(mp, means - stDevs, mp, means + stDevs, lwd=2)

segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)

segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)
abline(h=0)
box()
dev.off()

######################################################################
tmp_cpm2<-tmp_cpm

### illumina main
tmp_cpm2<-tmp_cpm2[,c("lane5IV1","lane5IV2","lane5IV5","lane5IV6","lane5IV7","lane5IV8","lane5IV9","lane5IV10","lane5IV11","lane5IV12")]
# ### illumina rest
# tmp_cpm2<-tmp_cpm2[,c("lane5IV3","lane5IV4","lane5IV17","lane5IV18","lane5IV19","lane5IV20","lane5IV21","lane5IV22","lane5IV23","lane5IV24")]


for (i in 1:length(genes_select))
{
png(file= paste(genes_select[i],"expression_before_after.png",sep=""),height=600,width=800)
barplot(tmp_cpm2[genes_select[i],],main=genes_select[i],axes=F)
axis(2,at=seq(0,max(tmp_cpm2[genes_select[i],]),(max(tmp_cpm2[genes_select[i],]))/5),
     labels=round(seq(0,max(tmp_cpm2[genes_select[i],]),(max(tmp_cpm2[genes_select[i],]))/5),0))
#box()
dev.off()
}
######################################################################
###CALCULATION OF LOGFC BY EDGER
###log2(mean(cpm(d,normalized.lib.sizes=T)["hsa-miR-34a-5p",]))
###is the same as
###tmp10<-NULL
###tmp10<-d$AveLogCPM
###tmp10<-as.data.frame(tmp10)
###rownames(tmp10)<-rownames(d$counts)
###tmp10["hsa-miR-34a-5p",]
# tmp_cpm10<-NULL
# tmp_cpm10<-cpm(d,normalized.lib.sizes=T)
# 
# for (i in genes_select)
# {
#   print(i)
#   print(log2(mean(tmp_cpm10[i,])))
# }
# 
# order(tmp_cpm10)

mirna_norm_count<-NULL;norm_count_ave<-NULL;total_table<-NULL
mirna_norm_count<-as.data.frame(tmp_cpm)
norm_count_ave<-as.data.frame(rowMeans(mirna_norm_count))
mirna_norm_count$average<-norm_count_ave[,1]
write.table(mirna_norm_count,"all_mirnas_norm_counts.csv",sep=";",dec=".",col.names=NA)
write.table(mirna_norm_count[genes_select,],"selected_mirnas_norm_counts.csv",sep=";",dec=".",col.names=NA)

colnames(mirna_norm_count)<-paste(colnames(mirna_norm_count),"norm_count",sep="_")

total_table<-merge(combResults.tgw,mirna_norm_count,by="row.names",all.x=T)
colnames(total_table)[9:(ncol(combResults.tgw)+1)]<-paste(colnames(total_table)[9:(ncol(combResults.tgw)+1)],"raw_counts",sep="_")
rownames(total_table)<-total_table[,1]
total_table<-subset(total_table, select = -c(Row.names) )

mirna_norm_count<-NULL;norm_count_ave<-NULL;tmp_cpm2<-NULL
mirna_norm_count<-as.data.frame(log2(tmp_cpm))
is.na(mirna_norm_count) <- sapply(mirna_norm_count, is.infinite)
mirna_norm_count[is.na(mirna_norm_count)==T]<-0
norm_count_ave<-as.data.frame(rowMeans(mirna_norm_count))
mirna_norm_count$average<-norm_count_ave[,1]
write.table(mirna_norm_count,"all_mirnas_log2_norm_counts.csv",sep=";",dec=".",col.names=NA)
write.table(mirna_norm_count[genes_select,],"selected_mirnas_log2_norm_counts.csv",sep=";",dec=".",col.names=NA)

colnames(mirna_norm_count)<-paste(colnames(mirna_norm_count),"log2_norm_count",sep="_")

total_table<-merge(total_table,mirna_norm_count,by="row.names",all.x=T)
rownames(total_table)<-total_table[,1]
total_table<-subset(total_table, select = -c(Row.names) )


mirna_norm_count<-NULL;norm_count_ave<-NULL;tmp_cpm2<-NULL
mirna_norm_count<-as.data.frame(counts.normalized)
# is.na(mirna_norm_count) <- sapply(mirna_norm_count, is.infinite)
# mirna_norm_count[is.na(mirna_norm_count)==T]<-0
norm_count_ave<-as.data.frame(rowMeans(mirna_norm_count))
mirna_norm_count$average<-norm_count_ave[,1]
write.table(mirna_norm_count,"all_mirnas_frequencies.csv",sep=";",dec=".",col.names=NA)
write.table(mirna_norm_count[genes_select,],"selected_mirnas_frequencies.csv",sep=";",dec=".",col.names=NA)

colnames(mirna_norm_count)<-paste(colnames(mirna_norm_count),"frequencies",sep="_")

total_table<-merge(total_table,mirna_norm_count,by="row.names",all.x=T)
rownames(total_table)<-total_table[,1]
total_table<-subset(total_table, select = -c(Row.names) )

write.table(d$counts,"all_mirnas_counts.csv",sep=";",dec=".",col.names=NA)
write.table(d$counts[genes_select,],"selected_mirnas_counts.csv",sep=";",dec=".",col.names=NA)

mirna_norm_count<-NULL;norm_count_ave<-NULL;tmp_cpm2<-NULL
mirna_norm_count<-as.data.frame(logFC_table)
colnames(mirna_norm_count)<-paste(colnames(mirna_norm_count),"manual_logFC",sep="_")
total_table<-merge(total_table,mirna_norm_count,by="row.names",all.x=T)
rownames(total_table)<-total_table[,1]
total_table<-subset(total_table, select = -c(Row.names) )

write.table(total_table,file="all_results.csv",sep=";",dec=".",col.names=NA)
write.table(total_table[genes_select,],file="all_results_selected_mirnas.csv",sep=";",dec=".",col.names=NA)


# #manual calculation of logFC should be done as following
# # log2FC = log2(expression in mutant backround) - log2(expression in wildtype)
# # or
# # log10FC = log10(expression in mutant backround) - log10(expression in wildtype)
# manual_logFC<-NULL
# manual_logFC<-as.data.frame(matrix(ncol=ncol(d$counts),nrow=nrow(total_table)))
# ### illumina main
# colnames(manual_logFC)<-c("lane5IV2/lane5IV1","lane5IV6/lane5IV5","lane5IV8/lane5IV7","lane5IV10/lane5IV9","lane5IV12/lane5IV11")
# manual_logFC$"lane5IV2/lane5IV1"<-total_table$lane5IV2_log2_norm_count/total_table$lane5IV1_log2_norm_count

# tmp_cpm2<-NULL
# tmp_cpm2<-tmp_cpm
# tmp_cpm2<-log2(tmp_cpm2)
# is.na(tmp_cpm2) <- sapply(tmp_cpm2, is.infinite)
# tmp_cpm2[is.na(tmp_cpm2)]<-0.000000000000000001
# heatmap.2(tmp_cpm2,scale="row")


tmp_avelog<-NULL
tmp_avelog<-as.data.frame(2^d$AveLogCPM)
rownames(tmp_avelog)<-rownames(d$counts)
colnames(tmp_avelog)<-"average_cpm"
tmp_avelog<-tmp_avelog[order(-tmp_avelog$average_cpm), , drop = FALSE]
write.table(tmp_avelog,file="average_cpm_all.csv",sep=";",col.names=NA,dec=".")

tmp_avelog<-as.data.frame(tmp_avelog)
tmp_avelog<-as.data.frame(tmp_avelog[genes_select,])
rownames(tmp_avelog)<-genes_select
colnames(tmp_avelog)<-"average_cpm"
write.table(tmp_avelog,file="average_cpm_selected.csv",sep=";",col.names=NA,dec=".")


# ####TESTING OF MERGING
# selection_data_total<-NULL
# merge(x=selection_data_total,y=selection.data,by="mirna",all=T)
# #if it doesn't work
# library(sqldf)
# sqldf("SELECT * FROM selection_data_total LEFT JOIN selection.data USING mirna")

order_select<-NULL
order_select<-as.data.frame(total_table[genes_select,(8+ncol(RCtableEdgeR)):(8+(2*(ncol(RCtableEdgeR))))])
rownames(order_select)<-genes_select
order_select$mirna<-rownames(order_select)
order_select<-arrange(order_select,desc(average_norm_count))
rownames(order_select)<-order_select$mirna
order_select$mirna<-NULL

tmp_colnames<-NULL
tmp_colnames<-colnames(order_select)
colnames(order_select)<-gsub("_norm_count", "", colnames(order_select)) #remove "_norm_count" from the name

# png(file= "heatmap_selected_mirnas_norm_counts.png",height=600,width=800)
# heatmap.2(as.matrix(order_select[,1:(ncol(order_select))-1]),scale="row",margins=c(12,12),dendrogram="column",Rowv=F)
# dev.off()

png(file= "heatmap_selected_mirnas_norm_counts.png",height=1200,width=1600)
heatmap.2(as.matrix(order_select[,1:(ncol(order_select))-1]),margins=c(12,12),dendrogram="none",scale="row",Colv=F,Rowv=F,trace="none")
dev.off()

colnames(order_select)<-tmp_colnames
###########
order_select<-NULL
order_select<-as.data.frame(total_table[genes_select,(8+(2*(ncol(RCtableEdgeR))+1)):(8+(3*(ncol(RCtableEdgeR))+1))])
rownames(order_select)<-genes_select
order_select$mirna<-rownames(order_select)
order_select<-arrange(order_select,desc(average_log2_norm_count))
rownames(order_select)<-order_select$mirna
order_select$mirna<-NULL

tmp_colnames<-NULL
tmp_colnames<-colnames(order_select)
colnames(order_select)<-gsub("_log2_norm_count", "", colnames(order_select)) #remove "_norm_count" from the name

png(file= "heatmap_selected_mirnas_log_norm_counts.png",height=600,width=800)
heatmap.2(as.matrix(order_select[,1:(ncol(order_select))-1]),scale="row",margins=c(14,12),dendrogram="column",Rowv=F)
dev.off()

colnames(order_select)<-tmp_colnames
#############################
order_select<-NULL
order_select<-as.data.frame(total_table[,(8+ncol(RCtableEdgeR)):(8+(2*(ncol(RCtableEdgeR))))])
rownames(order_select)<-rownames(total_table)
order_select$mirna<-rownames(order_select)
order_select<-arrange(order_select,desc(average_norm_count))
rownames(order_select)<-order_select$mirna
order_select$mirna<-NULL

tmp_colnames<-NULL
tmp_colnames<-colnames(order_select)
colnames(order_select)<-gsub("_norm_count", "", colnames(order_select)) #remove "_norm_count" from the name

png(file= "heatmap_all_mirnas_norm_counts.png",height=600,width=800)
heatmap.2(as.matrix(order_select[,1:(ncol(order_select))-1]),scale="row",margins=c(12,12),dendrogram="column",Rowv=F)
dev.off()

colnames(order_select)<-tmp_colnames
#############################
order_select<-NULL
order_select<-as.data.frame(total_table[,(8+(2*(ncol(RCtableEdgeR))+1)):(8+(3*(ncol(RCtableEdgeR))+1))])
rownames(order_select)<-rownames(total_table)
order_select$mirna<-rownames(order_select)
order_select<-arrange(order_select,desc(average_log2_norm_count))
rownames(order_select)<-order_select$mirna
order_select$mirna<-NULL

tmp_colnames<-NULL
tmp_colnames<-colnames(order_select)
colnames(order_select)<-gsub("_log2_norm_count", "", colnames(order_select)) #remove "_norm_count" from the name

png(file= "heatmap_all_mirnas_log_norm_counts.png",height=600,width=800)
heatmap.2(as.matrix(order_select[,1:(ncol(order_select))-1]),scale="row",margins=c(14,12),dendrogram="column",Rowv=F)
dev.off()

colnames(order_select)<-tmp_colnames
# ####TEST
selection.data.total.names.tmp<-as.data.frame(selection.data2$mirna)
selection.data.total.names<-rbind(selection.data.total.names,selection.data.total.names.tmp)
selection.data.total.names<-unique(selection.data.total.names)

rownames(selection.data.total.tmp)<-selection.data.total.names[,1]
selection.data.total.tmp[,1]<-selection.data2$norm_count
colnames(selection.data.total.tmp)<-sample
###TEST

total_table_tmp10<-NULL
total_table_tmp10<-total_table[,(8+(2*(ncol(RCtableEdgeR))+1)):(8+(3*(ncol(RCtableEdgeR)))+1)]
total_table_tmp10$mirnas<-rownames(total_table_tmp10)

library(plyr)
total_table_tmp10<-arrange(total_table_tmp10,desc(average_log2_norm_count))
rownames(total_table_tmp10)<-total_table_tmp10$mirnas
colnames(total_table_tmp10)
total_table_tmp10$mirnas<-NULL

tmp_colnames<-NULL
tmp_colnames<-colnames(total_table_tmp10)
colnames(total_table_tmp10)<-gsub("_log2_norm_count", "", colnames(total_table_tmp10)) #remove "_log2_norm_count" from the name

png(file= "heatmap_top20_log_norm_counts_scaled_column_bw.png",height=1200,width=1600)
heatmap.2(as.matrix(total_table_tmp10[1:20,]),margins=c(14,12),dendrogram="column",scale="column",Rowv=F, col=colorpanel(9, "grey10", "white"),tracecol="black",trace="none")
dev.off()

png(file= "heatmap_top20_log_norm_counts_scaled_none_bw.png",height=1200,width=1600)
heatmap.2(as.matrix(total_table_tmp10[1:20,]),margins=c(14,12),dendrogram="column",scale="none",Rowv=F, col=colorpanel(9, "grey10", "white"),tracecol="black",trace="none")
dev.off()

png(file= "heatmap_top20_log_norm_counts_scaled_row_bw.png",height=1200,width=1600)
heatmap.2(as.matrix(total_table_tmp10[1:20,]),margins=c(14,12),dendrogram="column",scale="row",Rowv=F, col=colorpanel(9, "grey10", "white"),tracecol="black",trace="none")
dev.off()
##############################3
png(file= "heatmap_top20_log_norm_counts_scaled_column.png",height=1200,width=1600)
heatmap.2(as.matrix(total_table_tmp10[1:20,]),margins=c(14,12),dendrogram="column",scale="column",Rowv=F,trace="none")
dev.off()

png(file= "heatmap_top20_log_norm_counts_scaled_none.png",height=1200,width=1600)
heatmap.2(as.matrix(total_table_tmp10[1:20,]),margins=c(14,12),dendrogram="column",scale="none",Rowv=F,trace="none")
dev.off()

png(file= "heatmap_top20_log_norm_counts_scaled_row.png",height=1200,width=1600)
heatmap.2(as.matrix(total_table_tmp10[1:20,]),margins=c(14,12),dendrogram="column",scale="row",Rowv=F,trace="none")
dev.off()

colnames(total_table_tmp10)<-tmp_colnames
##############################
##############################reordered
total_table_tmp10_reor<-NULL
#illumina main #order made by clustering order
total_table_tmp10_reor<-total_table_tmp10[c("lane5IV1_log2_norm_count","lane5IV9_log2_norm_count","lane5IV7_log2_norm_count","lane5IV5_log2_norm_count","lane5IV11_log2_norm_count",
                                            "lane5IV2_log2_norm_count","lane5IV10_log2_norm_count","lane5IV6_log2_norm_count","lane5IV8_log2_norm_count","lane5IV12_log2_norm_count","average_log2_norm_count")]
#original order illumina main
## total_table_tmp10_reor<-total_table_tmp10[c("lane5IV1_log2_norm_count","lane5IV5_log2_norm_count","lane5IV7_log2_norm_count","lane5IV9_log2_norm_count","lane5IV11_log2_norm_count",
##                                             "lane5IV2_log2_norm_count","lane5IV6_log2_norm_count","lane5IV8_log2_norm_count","lane5IV10_log2_norm_count","lane5IV12_log2_norm_count","average_log2_norm_count")]
#######################################################
#illumina rest #order made by clustering order
# total_table_tmp10_reor<-total_table_tmp10[c("lane5IV17_log2_norm_count","lane5IV23_log2_norm_count","lane5IV21_log2_norm_count","lane5IV19_log2_norm_count","lane5IV3_log2_norm_count",
#                                             "lane5IV18_log2_norm_count","lane5IV24_log2_norm_count","lane5IV22_log2_norm_count","lane5IV20_log2_norm_count","lane5IV4_log2_norm_count","average_log2_norm_count")]
#original order illumina rest
## total_table_tmp10_reor<-total_table_tmp10[c("lane5IV3_log2_norm_count","lane5IV17_log2_norm_count","lane5IV19_log2_norm_count","lane5IV21_log2_norm_count","lane5IV23_log2_norm_count",
##                                             "lane5IV4_log2_norm_count","lane5IV18_log2_norm_count","lane5IV20_log2_norm_count","lane5IV22_log2_norm_count","","average_log2_norm_count")]

colnames(total_table_tmp10_reor)<-gsub("_log2_norm_count", "", colnames(total_table_tmp10_reor)) #remove "_log2_norm_count" from the name

png(file= "heatmap_top20_log_norm_counts_scaled_column_bw_reor.png",height=1200,width=1600)
heatmap.2(as.matrix(total_table_tmp10_reor[1:20,]),margins=c(14,12),dendrogram="none",scale="column",Rowv=F,Colv=F, col=colorpanel(9, "grey10", "white"),tracecol="black",trace="none", add.expr = abline(v = c(5.5,10.5), lwd = 3))
dev.off()

png(file= "heatmap_top20_log_norm_counts_scaled_none_bw_reor.png",height=1200,width=1600)
heatmap.2(as.matrix(total_table_tmp10_reor[1:20,]),margins=c(14,12),dendrogram="none",scale="none",Rowv=F,Colv=F, col=colorpanel(9, "grey10", "white"),tracecol="black",trace="none", add.expr = abline(v = c(5.5,10.5), lwd = 3))
dev.off()

png(file= "heatmap_top20_log_norm_counts_scaled_row_bw_reor.png",height=1200,width=1600)
heatmap.2(as.matrix(total_table_tmp10_reor[1:20,]),margins=c(14,12),dendrogram="none",scale="row",Rowv=F,Colv=F, col=colorpanel(9, "grey10", "white"),tracecol="black",trace="none", add.expr = abline(v = c(5.5,10.5), lwd = 3))
dev.off()
##############################3
png(file= "heatmap_top20_log_norm_counts_scaled_column_reor.png",height=1200,width=1600)
heatmap.2(as.matrix(total_table_tmp10_reor[1:20,]),margins=c(14,12),dendrogram="none",scale="column",Rowv=F,Colv=F,trace="none", add.expr = abline(v = c(5.5,10.5), lwd = 3))
dev.off()

png(file= "heatmap_top20_log_norm_counts_scaled_none_reor.png",height=1200,width=1600)
heatmap.2(as.matrix(total_table_tmp10_reor[1:20,]),margins=c(14,12),dendrogram="none",scale="none",Rowv=F,Colv=F,trace="none", add.expr = abline(v = c(5.5,10.5), lwd = 3))
dev.off()

png(file= "heatmap_top20_log_norm_counts_scaled_row_reor.png",height=1200,width=1600)
heatmap.2(as.matrix(total_table_tmp10_reor[1:20,]),margins=c(14,12),dendrogram="none",scale="row",Rowv=F,Colv=F,trace="none", add.expr = abline(v = c(5.5,10.5), lwd = 3))
dev.off()





