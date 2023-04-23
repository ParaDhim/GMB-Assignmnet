# load the libraries
library(RUVSeq)
library(zebrafishRNASeq)

# load the data
data(zfGenes)

# glimpse of the data
head(zfGenes)
tail(zfGenes)

# filter out the genes with little to no expression
filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=2)
head(filter)

filtered <- zfGenes[filter,]

# seperate data into genes and spikes:
genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]
spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]

head(genes)

# convert to S4 object:
x <- as.factor(rep(c("Ctl", "Trt"), each=3))
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x, row.names=colnames(filtered)))
set


# Check for variability in the replicates:
library(RColorBrewer)

colors <- brewer.pal(3, "Set2")

plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)


# Upper quantile normalization to reduce batch effects:
set <- betweenLaneNormalization(set, which="upper")

plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

# Spike ins to detect unwanted variations:
set1 <- RUVg(set, spikes, k=1)
pData(set1)

plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set1, col=colors[x], cex=1.2)

# Differential Gene Expression Analysis ####
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts(set1),
                              colData = pData(set1),
                              design = ~ W_1 + x)

dds <- DESeq(dds)
res <- results(dds)
res

# convert to dataframe
output <- as.data.frame(res)
#write.csv(output,"D:\\iiit delhi\\sem 3\\Assignment\\GMB\\ASSMT2\\assmt.csv",row.names = TRUE)


# volcano plots:
final <- data.frame(row.names(output))
final[,c(2,3)] <- output[,c(2,5)]


names(final)[1] <- "Gene"

# Make a basic volcano plot
with(final, plot(log2FoldChange, -log10(pvalue), pch = 20, main = "ParasDhiman_2021482"))
ff <- subset(final, pvalue < .05 & log2FoldChange > 1)
ff2 <- ff[1:50,]
write.csv(ff2,"D:\\iiit delhi\\sem 3\\Assignment\\GMB\\ASSMT2\\assmt.csv",row.names = TRUE)

ff1 <- subset(final, pvalue < .05 & log2FoldChange < -1)
ff3 <- ff1[1:50,]
write.csv(ff3,"D:\\iiit delhi\\sem 3\\Assignment\\GMB\\ASSMT2\\assmt1.csv",row.names = TRUE)

with(subset(final, pvalue < .05 & abs(log2FoldChange) > 1), points(log2FoldChange, -log10(pvalue), pch = 20, col = "blue"))



