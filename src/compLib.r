#! /usr/bin/Rscript

args = commandArgs(trailingOnly=T)
input = args[1]
output = args[2]

x = read.table(input, header = T, row.names = 1)

lib_sum = sort(apply(x, 2, sum), decreasing=T)

pdf(paste(output, 'figure1.pdf', sep=''),width=5,heigh=5)
par(mar=c(7,4,4,2))
bp = barplot(lib_sum, xaxt = 'n', las = 2, ylab = 'Number of circRNA')
text(bp, lib_sum/2, lib_sum)
text(bp+0.5, rep(-max(lib_sum)/20,length(bp)), names(lib_sum), srt = 45, xpd = T,pos=2)
dev.off()

circ_occurs = apply(x, 1, sum)
circ_occurs_summary = table(circ_occurs)

pdf(paste(output, 'figure2.pdf', sep=''),width=5,heigh=5)
par(mar=c(5,4,4,2))
bp = barplot(circ_occurs_summary,xlab='Number of library',ylab='Number of circRNA', las = 1, main = paste('Total number of circRNA:', nrow(x))) 
text(bp, circ_occurs_summary/2, circ_occurs_summary)
dev.off()

#pdf(paste(output, 'figure3.pdf', sep=''),width=5,heigh=7)
#heatmap(as.matrix(x), scale="none", labRow=F, col=c('white','red'))
#dev.off()


x_filt = x[circ_occurs>=3,]
x_filt_sum = apply(x_filt, 2, sum)[names(lib_sum)]
prop = rbind(x_filt_sum/lib_sum, (lib_sum-x_filt_sum)/lib_sum)
prop = prop[, order(prop[1, ])]


#pdf(paste(output, 'figure4.pdf', sep=''),width=5,heigh=5)
#par(mar=c(7,4,4,2))
#bp = barplot(prop, xaxt = 'n', las = 2, ylab = 'Proportion', main = 'circRNA in >= 3 libraries')
#text(bp+0.5, rep(-max(prop)/20,length(bp)), colnames(prop), srt = 45, xpd = T,pos=2)
#text(bp, prop[1, ]/2, x_filt_sum[colnames(prop)], col = 'white')
#text(bp, prop[1, ]+prop[2, ]/2, (lib_sum-x_filt_sum)[colnames(prop)])
#dev.off()

