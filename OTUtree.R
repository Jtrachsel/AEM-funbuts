setwd("~/NexAmp/framecorr/figglin/coproadd/")


library(ggtree)
library(ape)


############## ############### ################### ################ #############


metadata <- read.csv("reftreeMETA.csv", header = TRUE, stringsAsFactors = FALSE)
shared <- read.csv("DvR.QCD.OTUsOrdered.shared.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(shared) <- shared[,2]

shared <- shared[,4:95]
abundance <- colSums(shared)

abundance <- data.frame(abundance)
rownames(metadata) <- metadata[,1]
rownames(abundance)
rownames(metadata)
metadata <- merge(metadata, abundance, by = 0, all.x = TRUE)
rownames(metadata)<- metadata[,1]
metadata <- metadata[,2:8]

#####################

short <- read.raxml("RAxML_bipartitionsBranchLabels.shorts")
s <- ggtree(short)

tiff
s <- s %<+% metadata

s <- s + geom_hilight(node = 126, fill = 'green', alpha = .4, extendto = 5.1) + geom_hilight(node = 226, fill = 'yellow', alpha = .5, extendto= 5.1) +
  geom_hilight(node = 118, fill = 'red', alpha = .4, extendto = 5.1) + scale_color_continuous(name="Total Reads", limits=c(1,500),oob=scales::squish, low='white', high = 'black') +
  theme(legend.position=c(.175,.5), legend.text = element_text(size = 14), legend.title = element_text(size = 18, family = 'bold'),
        legend.background = element_rect(fill = 'gray90', size = .5), legend.key.size = unit(1.5,'cm'))+ geom_tippoint(aes(color=abundance), size=5)

x1 <- s$data$x[s$data$functCONF == 'TRUE']
y1 <- s$data$y[s$data$functCONF == 'TRUE']
t1 <- s$data$label[s$data$functCONF == 'TRUE']
confirmed <- as.data.frame(cbind(x1,y1,t1))

x2 <- s$data$x[s$data$functOTHER == 'TRUE']
y2 <- s$data$y[s$data$functOTHER == 'TRUE']
t2 <- s$data$label[s$data$functOTHER == 'TRUE']
other <- as.data.frame(cbind(x2,y2,t2))

x3 <- s$data$x[s$data$functUNKNOWN == "TRUE"]
y3 <- s$data$y[s$data$functUNKNOWN == 'TRUE']
t3 <- s$data$label[s$data$functUNKNOWN =='TRUE']
unknown <- as.data.frame(cbind(x3,y3,t3))

xb <- s$data$x[s$data$bootstrap >50]
yb <- s$data$y[s$data$bootstrap >50]
lb <- s$data$bootstrap[s$data$bootstrap >50]
boots <- as.data.frame(cbind(xb,yb,lb))


s+ geom_point(data=confirmed, x=x1, y=y1, color="darkgreen", shape=20, size=6) + geom_point(data=other, x=x2,y=y2, color="red", shape=20, size=6) +
 geom_text(data=confirmed, x=x1+.05, y=y1, label=t1, color="darkgreen", size = 5, fontface='bold.italic', hjust = 'left') + geom_text(data=other, x=x2+.05,y=y2, label=t2, color="red", size = 5, fontface='bold.italic', hjust = 'left')+
  geom_text(data=unknown, x=x3+.15, y=y3, label=t3, color='black', size=3.5)+
  geom_text(data=boots, x=xb, y=yb, label=lb, color='black', size=3, hjust=.25) + ggplot2::xlim(0,5.1)


#################################################################################################3
# this section calculates how many sequences are in the main clade with all functional reference seqs
s + geom_text2(aes(subset=!isTip, label=node))
# node125 contains all good seqs#
# node 126 all good, node 226 somegood node 118 no good

library('caper')

good.tips <- clade.members(125, short@phylo)

goods <- sum(s$data$abundance[s$data$node %in% good.tips], na.rm = TRUE)

bads <- sum(s$data$abundance[!(s$data$node %in% good.tips)], na.rm = TRUE)

tot <- goods + bads

goods/tot * 100

bads/tot *100





######################### supplemental reference tree  ###########################################

raxml <- read.raxml("~/NexAmp/framecorr/figglin/coproadd/RAxML_bipartitionsBranchLabels.reftreeO2")
r <- ggtree(raxml)


r <- r %<+% metadata

x1 <- r$data$x[r$data$functCONF == 'TRUE']
y1 <- r$data$y[r$data$functCONF == 'TRUE']
t1 <- r$data$label[r$data$functCONF == 'TRUE']
confirmed <- as.data.frame(cbind(x1,y1,t1))


x2 <- r$data$x[r$data$functOTHER == 'TRUE']
y2 <- r$data$y[r$data$functOTHER == 'TRUE']
t2 <- r$data$label[r$data$functOTHER == 'TRUE']
other <- as.data.frame(cbind(x2,y2,t2))



xb <- r$data$x[r$data$bootstrap >50]
yb <- r$data$y[r$data$bootstrap >50]
lb <- r$data$bootstrap[r$data$bootstrap >50]
boots <- as.data.frame(cbind(xb,yb,lb))

par(mar = c(1,1,2,10))

r + geom_point(data=confirmed, x=x1, y=y1, color="darkgreen", shape=20, size=6) + geom_point(data=other, x=x2,y=y2, color="red", shape=20, size=6) +
  geom_text(data=confirmed, x=x1+.05, y=y1, label=t1, color="darkgreen", size = 5, fontface='bold.italic', hjust = 'left') + geom_text(data=other, x=x2+.05,y=y2, label=t2, color="red", size = 5, fontface='bold.italic', hjust = 'left')+
  geom_text(data=boots, x=xb, y=yb, label=lb, color='black', size=3, hjust=.25) + ggplot2::xlim(0,6)

