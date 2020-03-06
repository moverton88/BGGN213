title: 'Class10: Machine Learning Project'
author: "Michael Overton"
date: "2/7/2020"
output: github_document

## Using PCA and heirarchical clustering to investigate data from 
```{r, echo=F}
fna.data <- "Lecture10_data/WisconsinCancer.csv"

wisc.df <- read.csv(fna.data)
# View(wisc.df)
class(wisc.df)
dim(wisc.df)
sum(wisc.df$diagnosis == "M")

wisc.data <- as.matrix(wisc.df[ ,3:32])
rownames(wisc.data) <- wisc.df$id
head(wisc.data)

diagnosis <- as.vector(wisc.df$diagnosis)

```

Q5. Why do you think we are using the indices 3:32 here?
  The specific sample IDs are converted into the rownames, and we want this to be unsupervised, so we do not
  want to bias our analysis with the final diagnoses, we want to be able to predict malignancy from the data
  
```{r, echo=F}
length(grep("_mean", colnames(wisc.data)))

colMeans(wisc.data)
apply(wisc.data,2,sd)

wisc.pr <- prcomp(wisc.data, scale. = T)

summary(wisc.pr)


```


```{r, echo=F}
biplot(wisc.pr)

```
Q10. What stands out to you about this plot? Is it easy or difficult to understand? Why?
The plot is extremely cluttered, as ID names are used as points rather than some marker. As well, the various factors influencing the data point locations are numerous and the names are cluttered.

```{r}

plot(wisc.pr$x, col = (diagnosis=="M")+1, 
     xlab = "PC1", ylab = "PC2")

```


Q11. Generate a PC1 vs PC2 plot as described above colored by diagnosis. What do the points in this plot represent? What are the red points? Where did this coloring come from and why were these colors chosen by R?
  Each point is a particular biopsy result, the red points are malignant diagnoses. By using the conversion of logical to numeric, R looks up the associated colors in a palette to the values 1 and 0 (or actually, 2 and 1, because we add 1 to the logical) which are 1=black and 2=red.
  

```{r}
plot(wisc.pr$x[,c(1,3)], col = (diagnosis=="M")+1, 
     xlab = "PC1", ylab = "PC3")
```

Q12. Generate a similar plot for principal components 1 and 3. What do you notice about these plots?
  PC3 is less informative than PC2. There is less partitioning of the datapoints in PC3
  
```{r}
pr.var <- wisc.pr$sdev^2
pve <- pr.var / sum(pr.var)


plot(pve, xlab = "Principal Component", 
     ylab = "Proportion of Variance Explained", 
     ylim = c(0, 1), type = "o")

barplot(pve, ylab = "Precent of Variance Explained",
     names.arg=paste0("PC",1:length(pve)), las=2, axes = FALSE)
axis(2, at=pve, labels=round(pve,2)*100 )
```
  
```{r}
library(factoextra)

fviz_eig(wisc.pr, addlabels = TRUE)
```


```{r}
wisc.pr$rotation[,1][["radius_mean"]]
wisc.pr$rotation[,1][["smoothness_se"]]
max(abs(wisc.pr$rotation[,1]))
wisc.pr$rotation[,1]
```


```{r}
# Scale the wisc.data data: data.scaled
library(scales)
data.scaled <- apply(wisc.data, MARGIN = 2, FUN=rescale)
head(data.scaled)


data.dist <- dist(data.scaled)

wisc.hclust <- hclust(data.dist, method="complete")

plot(wisc.hclust)
```

```{r}

wisc.hclust.k2 <- cutree(wisc.hclust, k=2)
table(wisc.hclust.k2, diagnosis)

wisc.hclust.k3 <- cutree(wisc.hclust, k=3)
table(wisc.hclust.k3, diagnosis)

wisc.hclust.k4 <- cutree(wisc.hclust, k=4)
table(wisc.hclust.k4, diagnosis)

wisc.hclust.k5 <- cutree(wisc.hclust, k=5)
table(wisc.hclust.k5, diagnosis)

wisc.hclust.k10 <- cutree(wisc.hclust, k=10)
table(wisc.hclust.k10, diagnosis)
```

# clustering with PCA

```{r}


wisc.hc.pr <- hclust(dist(wisc.pr$x[,1:2]), method="ward.D2")
#plot(wisc.hc.pr) 
grps <- cutree(wisc.hc.pr, k=2)
table(grps, diagnosis)
plot(wisc.pr$x[ ,1:2], col=grps, pch=20+(diagnosis=="M"))
```






