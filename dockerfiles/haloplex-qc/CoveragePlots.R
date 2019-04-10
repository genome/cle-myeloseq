require(scales)

sample.name <- commandArgs(T)[1]

pdf(height=8.5,width=11,file=paste0(sample.name,".coverage_qc.pdf"))
layout(matrix(c(1,1,2:5),nrow=3,byrow = T),heights = c(.2,1,1))

op <- par(mar=c(2,2,2,2))
plot.new()
text(0,.5,labels = sample.name,cex=2,font=2,xpd=T,pos=4)
par(op)

cov <- read.table(paste0(sample.name,".coverage.txt"),stringsAsFactors = F)
cov$V5 <- gsub("_.+","",cov$V5)
genes <- levels(factor(cov$V5))
gene.chunks <- split(genes, ceiling(seq_along(genes)/10))
invisible(lapply(gene.chunks,
       function(g){
         x <- subset(cov,cov$V5 %in% g)
         stats <- cbind(aggregate(x$V4 ~ x$V5,FUN = length),
                        aggregate(x$V4 ~ x$V5,FUN = mean),
                        aggregate(x$V4 ~ x$V5,FUN = median),
                        aggregate(x$V4 ~ x$V5,FUN = function(X){ length(which(X<50)) }))[,c(1,2,4,6,8)]
         rownames(stats) <- stats[,1]
         apply(stats[,1:5],1,function(X){ cat(c("COVERAGE",X,"\n"),sep="\t"); })
         op <- par(mar=c(4.5,4,4.5,2))
         boxplot(log2(x$V4+1) ~ x$V5,las=2,axes=F,col="light gray")
         axis(1,at=1:nrow(stats),labels=rownames(stats),las=2)
         yAx <- c(0,1,20,50,100,500,seq(1000,max(x$V4,1000),by=1000))
         axis(2,at=log2(yAx+1),labels = yAx)
         mtext("Coverage depth (by position)",side=2,line=2.5)
         box()
         abline(h=log2(c(20,50)+1),lty=2,lwd=3,col="red")
         mtext(c("Mean",round(stats[,3],0)),at=c(-.25,1:nrow(stats)),line=2.75,cex=.8)
         mtext(c("Median",round(stats[,4],0)),at=c(-.25,1:nrow(stats)),line=1.5,cex=.8)
         mtext(c("<50x",round(stats[,5],0)),at=c(-.25,1:nrow(stats)),line=0.25,cex=.8)
         invisible(par(op))
  }))
invisible(dev.off())
invisible(require(scales))
pdf(height=11,width=8.5,file=paste0(sample.name,".gc_length_qc.pdf"))
layout(matrix(1:3,nrow=3,byrow = T),heights = c(.2,1,1))

op <- par(mar=c(2,2,2,2))
plot.new()
text(0,.5,labels = sample.name,cex=2,font=2,xpd=T,pos=4)
par(op)

x <- read.table(paste0(sample.name,".amplicon_counts.txt"),stringsAsFactors = F)
op <- par(mar=c(4,4,2,1),cex=1.25)
plot(x$V3-x$V2,log2(x$V5),pch=16,cex=.4,col=alpha("black",.5),xlab="Amplicon length",ylab="Read counts (log2)")
length.loess <- loess.smooth(x$V3-x$V2,log2(x$V5))
lines(length.loess,col="blue",lwd=3)
abline(h=log2(50),col="red",lwd=3,lty=3)
abline(h=log2(20),col="red",lwd=3,lty=3)
cat(paste0(c("LENGTHFITX\t",paste0(length.loess$x,collapse=","))))
cat("\n")
cat(paste0(c("LENGTHFITY\t",paste0(length.loess$y,collapse=","))))
cat("\n")

plot(x$V6*100,log2(x$V5),pch=16,cex=.4,col=alpha("black",.5),xlab="Amplicon GC%",ylab="Read counts (log2)")
gc.loess <- loess.smooth(x$V6*100,log2(x$V5))
lines(gc.loess,col="blue",lwd=3)
abline(h=log2(50),col="red",lwd=3,lty=3)
abline(h=log2(20),col="red",lwd=3,lty=3)

cat(paste0(c("GCFITX\t",paste0(gc.loess$x,collapse=","))))
cat("\n")
cat(paste0(c("GCFITY\t",paste0(gc.loess$y,collapse=","))))
cat("\n")
par(op)
invisible(dev.off())
