setwd("~/NexAmp/framecorr/")

shared <- read.csv("DvR.QCD.forStackedBar.shared.csv", header = TRUE, check.names = FALSE) # mothur shared file that I modified slightly, will be available in repo
rownames(shared) <- shared$Group
shared <- shared[,-c(1,2,3)]

#####################################

# This section is for the NMDS figure #

nmds <- read.table("DvR.QCD.sorclass.0.03.lt.nmds.axes", header = TRUE) # this was a file output by mothur
metadata <- read.table("DvR.design", header = FALSE)   # my design file from mothur
colnames(metadata) <- c('group', 'set')
nmds <- merge(metadata, nmds)

DNA <- nmds[which(nmds$set == "DNA"),]
RNA <- nmds[which(nmds$set == "RNA"),]



plot(RNA$axis1, RNA$axis2, pch = c(7, 8, 9, 10, 11, 12), xlim = c(-.5,.5), ylim = c(-.5,.5),
     cex =2, lwd = 3, col = 'red',
     xlab = 'nmds axis 1', ylab = 'nmds axis 2', main = "Community membership \nNMDS, Sorenson Disimilarity, DNA v RNA", sub = "AMOVA p = 0.10, F = 1.96")
points(DNA$axis1, DNA$axis2, pch = c(7, 8, 9, 10, 11, 12), cex = 2, lwd = 3, col = 'blue')
segments(x0 = RNA$axis1, y0 = RNA$axis2 , x1 = DNA$axis1 , y1 = DNA$axis2, lty = 2)
rpoly <- chull(RNA$axis1, RNA$axis2)
dpoly <- chull(DNA$axis1, DNA$axis2)
polygon(RNA$axis1[rpoly], RNA$axis2[rpoly], col = 'tomato', density = 9)
polygon(DNA$axis1[dpoly], DNA$axis2[dpoly], col = 'skyblue2', density = 9)
legend("topright", c('RNA', 'DNA'), fill = c('red', 'blue'))
text(x = 0.4, y= -0.4, labels = "Stress = 0.107 \nRsq = 0.900")


# stacked barplot stuff #
## there are multiple different versions of this chart below, the last one is the one that was submitted for publication
# 25 most abundant OTUs for bar graph
library(RColorBrewer)
display.brewer.all()
reds <- brewer.pal(9,"Reds")
reds <- reds[4:7]
purps <- rev(brewer.pal(9,"Purples"))
purps <- purps[4:7]
orang <-brewer.pal(9,"Oranges")
orang <- orang[4:7]
greens <- brewer.pal(9,"Greens")
greens<- rev(greens[4:7])
blues <- brewer.pal(9,"Blues")
blues <- blues[4:7]
last <- brewer.pal(9,"YlOrRd")
last <- rev(last[4:7])
lastest <- brewer.pal(9, "PuRd")
lastest <- lastest[4:7]

coltab <- rbind(reds, purps, orang, greens, blues, last, lastest)
colvec <- c(coltab[,4], coltab[,3], coltab[,2], coltab[,1])
colvec <- c(colvec, 'blue', 'darkgrey')
colvec[2] <- "limegreen"
colvec[4] <- 'gold'
colvec[9] <- 'darkmagenta'
colvec[6] <- 'cyan'
colvec[27] <- 'red'

shared <- read.csv("DvR.QCD.forStackedBar.shared.csv", header = TRUE, check.names = FALSE)
rownames(shared) <- shared$Group
shared <- shared[,-c(1,2,3)]


corectorder <- colnames(shared)

abundOTUs <- sort(colSums(shared), decreasing = TRUE)[1:25]
names(abundOTUs)
barshared <- shared[,names(abundOTUs)]


otherOTUs <- !(colnames(shared) %in% names(abundOTUs))
otherOTUs <- rowSums(shared[,otherOTUs])
barshared$others <- otherOTUs

nameees <- sort(colnames(barshared))
length(nameees)
nameees <- nameees[c(2:26,1)]
barshared <- barshared[,nameees]
barshared <- barshared[c(1,7,2,8,3,9,4,10,5,11,6,12),]
rownames(barshared) <- c("18.DNA", "18.RNA", "62.DNA","62.RNA","7.DNA","7.RNA","71.DNA","71.RNA","85.DNA","85.RNA","87.DNA","87.RNA")
barshared <- barshared/1696
barshared <- barshared*100

par(mar=c(5.1,4.1,4.1,17.5), xpd = TRUE, font = 1)
barplot(t(barshared), col = colvec, las=2, ylab = "Percent Relative Abundance", space = c(04,0,0.4,0,0.4,0,0.4,0,0.4,0,0.4,0))
legend("topright", fill= rev(colvec), legend = rev(rownames(t(barshared))), cex = .7, inset = c(-.55,0))

length(colvec)
colvec <- colvec[1:26]
colvec[26] <- 'darkgrey'
colvec[8] <- 'red'
colvec[25] <- 'green'
colvec[22] <- 'brown'
colvec[15] <- 'blue'
colvec[17] <- 'darkgreen'
colvec[24] <- 'darkgoldenrod'
colvec[16] <- 'magenta'




par(mar=c(5.1,4.1,4.1,17.5), xpd = TRUE, font = 1)
barplot(t(barshared), col = colvec, las=2, ylab = "Percent Relative Abundance", space = c(04,0,0.4,0,0.4,0,0.4,0,0.4,0,0.4,0))
legend("topright", fill= rev(colvec), legend = rev(rownames(t(barshared))), cex = .7, inset = c(-.55,0))

################################################################################################################
# this is the code for the final version of this figure
# if you are running this in Rstudio, you'll need to enlarge the viewing window so the legend is in the right spot,
# either that or you'll need to tweak some settings in barplot() or par()

colvec <- c('firebrick1','chartreuse','blue','yellow','cyan','green','aquamarine3','darkgreen',
               'darkmagenta','darkorange','cornflowerblue','tan2','darkslateblue','darkslategray4','deepskyblue','magenta',
               'midnightblue','red','orangered','plum3','seagreen','saddlebrown','purple','khaki','darkgoldenrod')

colvec[26] <- 'darkgrey'

par(mar=c(5.1,4.1,4.1,19.5), xpd = TRUE, font = 3)
barplot(t(barshared), col = colvec, las=2, ylab = "Percent Relative Abundance", space = c(0.3,0,0.3,0,0.3,0,0.3,0,0.3,0,0.3,0), cex.axis = 1.4, cex.names = 1.4, cex.lab = 1.4)
legend("topright", fill= rev(colvec), legend = rev(rownames(t(barshared))), cex = .97, inset = c(-.67,0), y.intersp = 1.45, x.intersp = .3)



######################## number of metagenome sequence ####################
# these numbers are for the results section #

colnames(shared)

splits <- strsplit(colnames(shared), split = '-')

percents <- c()
for (x in splits){
  y <- x[length(x)]
  percents <- c(percents, y)
}

percents <- gsub('%', '', percents)
percents <- gsub('[*]', '', percents)
percents <- as.numeric(percents)
sum(percents < 90)
sum(percents < 80)
sum(percents < 70)
max(percents)
min(percents)
############## 50 most abundant OTUs DNA vs RNA ###################

shared <- barshared
DNA <- shared[1:6,]  
RNA <- shared[7:12,]

DNA[DNA > 0] <- 1
RNA[RNA > 0] <- 1


DNAsums <- colSums(DNA)
RNAsums <- colSums(RNA)

notDNA <- DNAsums[DNAsums == 0]
notRNA <- RNAsums[RNAsums == 0]

everyDNA <- DNAsums[DNAsums == 6]
everyRNA <- RNAsums[RNAsums == 6]

onlyDNA <- DNAsums[DNAsums > 0 & RNAsums == 0]
onlyRNA <- RNAsums[RNAsums > 0 & DNAsums == 0]

everyone <- DNAsums[DNAsums == 6 & RNAsums == 6]



