library(ggplot2)


#Reading in and fixing tables

data <- read.table("all_counts_fixed.txt", sep="\t")
key <- read.csv("KeySequencing_BC011.csv", header=TRUE)


head(key)
key$X <- as.factor(key$X)
key$Mouse <- as.factor(key$Mouse)
key$Route <- as.factor(key$Route)
key$Genotype <- as.factor(key$Genotype)
key$Antiviral <- as.factor(key$Antiviral)
key$Day <- as.factor(key$Day)
key$Experiment <- as.factor(key$Experiment)
str(key)

#Modified the above to introduce other factors not included in my original

data2 <- data[order(data$V2, data$V1),]
colnames(key)[1] <- c("Sample")
colnames(data2) <- c("Barcode_sequence","Sample","Reads")
data2$Barcode_sequence <- as.factor(data2$Barcode_sequence)
total_reads <- read.table("total_reads_fixed.txt")
colnames(total_reads) <- c("Total_Reads", "Sample")
total_reads$Total_Reads <- total_reads$Total_Reads/4

#This will get you to the first, generic plot showing total abundances of each barcode, as compared to the total reads per sample
data3 <- merge(key,data2,by.x="Sample",by.y="Sample")
data3 <- merge(data3,total_reads,by.x="Sample",by.y="Sample")
data3$Reads[is.na(data3$Reads)] <- 0
data4 <- matrix(data3$Reads, nrow=length(unique(data3$Barcode_sequence)))
rownames(data4) <- unique(data3$Barcode_sequence)
colnames(data4) <- unique(data3$Sample)
data5 <- matrix(data3$Total_Reads, nrow=length(unique(data3$Barcode_sequence)))
rownames(data5) <- unique(data3$Barcode_sequence)
colnames(data5) <- unique(data3$Sample)
library(viridis)
palette <- viridis(24)
par(mfrow=c(2,2))
barplot(data4, col=palette, main="Total abundance per barcode", ylab="Number of reads")
barplot(data5[1,], col=palette, main="Total read number", ylab="Number of reads")
percent_barcoded <- seq(1:ncol(data4))
for (i in (1:ncol(data4)))
{
  percent_barcoded[i] <- (sum(data4[,i])/data5[1,i])
}
boxplot(percent_barcoded, main="% barcoded", ylab="% of reads containing barcodes")
hist(data3$Total_Reads[seq(1, length(data3$Total_Reads), by=length(unique(data3$Barcode_sequence)))], xlab="Total number of reads", main="Total reads")


library(viridis)
palette <- viridis(24)


#Calculating relative abundances of barcodes and graphing

data3$Relative_abund <- data3$Reads/data3$Total_Reads
data6 <- matrix(data3$Relative_abund, nrow=length(unique(data3$Barcode_sequence)))
colnames(data6) <- unique(data3$Sample)
rownames(data6) <- unique(data3$Barcode_sequence)
par(bg="white", mfrow=c(1,1))
barplot(data6,col=palette, ylab="Relative abundance of barcode in reads", main="Relative abundance per barcode")

##Everything above here can be run at once and produce useful initial sanity-checking sort of graphs

#Now let's remove the unwanted rows

data3 <- data3[data3$Barcode_sequence!="GACGGATGGTAC",]
data3 <- data3[data3$Barcode_sequence!="GACGGAAGTTGTTGGTAC",]
data3 <- data3[data3$Barcode_sequence!="GACGGACTCGATTGGTAC",]
data3 <- data3[data3$Barcode_sequence!="GACGGAGTCAATTGGTAC",]
data3 <- data3[data3$Barcode_sequence!="GACGGAACCTTTTGGTAC",]
data3 <- data3[data3$Barcode_sequence!="GACGGAAACCGTTGGTAC",]

#Now Count the #reads per sample
library(dplyr)
data3A = data3 %>% group_by (SampleID) %>% summarize(Total_Reads_Mod = sum(Reads))

#merge Total_mod._reads in data3
data3 <- merge(data3,data3A,by.x ="SampleID",by.y = "SampleID")

#Calculating Mod. Rel. abundance
data3$Mod.Relative_abund <- data3$Reads/data3$Total_Reads_Mod


#Set read threshold >= 10 and calculate richness
read_cutoff = 10
data3$richness <- 0
data3$richness[data3$Reads >= read_cutoff] = 1

#Count #barcodes per sampleID

data3b = data3 %>% group_by (SampleID) %>% summarize(Total_richness = sum(richness))


#merge Total_richness in data3
data3 <- merge(data3,data3b,by.x ="SampleID",by.y = "SampleID")
data3ba <- data3[unique(data3$SampleID),]
library(reshape2)
data3ba <- dcast(data3[,c('Mouse','Genotype','Tissue','Barcode_ID','Total_richness', 'PercentBarcoded')])

#Calculate %barcoded reads
data3$PercentBarcoded <- data3$Total_Reads_Mod/data3$Total_Reads
data3$PercentBarcoded <- data3$PercentBarcoded*100

#Now let's subset things and pull out one experiment

data3$Day <- factor(data3$Day, levels=c('D-5', 'D-14', 'D-21'))

tuft_data <- data3[data3$Experiment=="Tuft",]
tuft_data_no <- tuft_data[tuft_data$IL.4=="None",]
p <- ggplot(tuft_data_no) + aes(fill=Barcode_sequence, y=Relative_abund, x=Mouse) + geom_bar(position="Stack", stat="identity")
p + facet_grid(cols=vars(Genotype), scale="free")

#And the other experiment

dissem_data <- data3[data3$Experiment=="Dissemination",]
g <- ggplot(dissem_data) + aes(fill=Barcode_sequence, y=Relative_abund, x=Mouse) + geom_bar(position="Stack", stat="identity")
g + facet_grid(rows=vars(Route), cols=vars(Genotype), scale="free")

#I haven't run anything below this for Somya's data, so I have no idea if it works as-is.

shannon <- as.data.frame(matrix(ncol=5,nrow=(length(unique(data3$Mouse)) * length(unique(data3$Tissue)))))
colnames(shannon) <- c("Mouse", "Tissue", "Genotype", "Virus", "Shannon")
for (q in (1:length(unique(data3$Tissue))))
{
  day <- data3[data3$Tissue==as.character(unique(data3$Tissue)[q]),]
  for (i in (1:length(unique(day$Mouse))))
  {
    shannon$Mouse[(((q-1)*11)+i)] <- as.vector(unique(day$Mouse))[i]
    shannon$Genotype[(((q-1)*11)+i)] <- as.character(unique(day$Genotype[day$Mouse==unique(day$Mouse)[i]]))
    shannon$Virus[(((q-1)*11)+i)] <- as.character(unique(day$Virus[day$Mouse==unique(day$Mouse)[i]]))
    shannon$Shannon[(((q-1)*11)+i)] <- diversity(day$Reads[day$Mouse==unique(day$Mouse)[i]])
    shannon$Tissue[(((q-1)*11)+i)] <- as.character(unique(day$Tissue))
  }
}


shannon$Mouse <- as.factor(shannon$Mouse)
shannon$Genotype <- as.factor(shannon$Genotype)
shannon$Tissue <- as.factor(shannon$Tissue)
shannon$Virus <- as.factor(shannon$Virus)

shannon2 <- shannon[shannon$Virus=="Pool",]

par(mfrow=c(1,1))
boxplot(shannon2$Shannon ~ shannon2$Genotype * shannon2$Tissue)
