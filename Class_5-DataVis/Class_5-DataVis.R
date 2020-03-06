
setwd("~/Desktop/All_Things_School/UCSD/Courses/BGGN213/BGGN213/bimm143_05_rstats")

# chart of baby weight to age
babies <- read.table("weight_chart.txt", header=T, sep="\t")
plot(x=babies$Age, y=babies$Weight, type="o", pch=15)

# counts of genomic features
f_c <- read.table("feature_counts.txt", header=T, sep="\t")

old.par <- par()$mar
par(mar=c(3.1, 11.1, 4.1, 2.1))
barplot(f_c$Count, names.arg=f_c$Feature, horiz=T, las=1)
dotchart(f_c$Count, labels=f_c$Feature)

# male and female counts with rainbow colors in plot

mf_counts <- read.table("male_female_counts.txt", header=T, sep="\t")

par(mar=c(7.1, 4.1, 4.1, 2.1))
barplot(mf_counts$Count, names.arg=mf_counts$Sample, col=rainbow(nrow(mf_counts)), las=2)
par(mar=old.par)

# gene expression profile
u_d_expr <- read.table("up_down_expression.txt", header=T, sep="\t")
head(u_d_expr)
nrow(u_d_expr)

table(u_d_expr$State)
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(u_d_expr$Condition1, u_d_expr$Condition2, col=u_d_expr$State)
palette(c("blue", "grey", "red"))

# methylation data
meth <- read.table("expression_methylation.txt", header=T, sep="\t")
head(meth)
nrow(meth)

plot(meth$gene.meth, meth$expression)
dcols <- densCols(meth$gene.meth, meth$expression)
plot(meth$gene.meth, meth$expression, col=dcols)

inds <- meth$expression > 0
dcols <- densCols(meth$gene.meth[inds], meth$expression[inds])
plot(meth$gene.meth[inds], meth$expression[inds], col=dcols, pch=20)

c_r <- colorRampPalette(c("blue", "green", "red", "yellow"))
dcols <- densCols(meth$gene.meth[inds], meth$expression[inds], 
                  colramp = colorRampPalette(c("blue", "green", "red", "yellow"))
                  )
plot(meth$gene.meth[inds], meth$expression[inds], col=dcols, pch=20)
smoothScatter(meth$gene.meth[inds], meth$expression[inds])
