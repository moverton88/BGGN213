Lecture\_9-Unsup\_Learn
================
Michael Overton
02.05.2020

K-means clustering
------------------

We will try to cluster data using the kmeans() function

``` r
tmp_data <- c(rnorm(30,-3), rnorm(30,0,0.4), rnorm(30,3))
x <- cbind(x=tmp_data, y=rev(tmp_data))
plot(x)
```

![](Lecture_9-Unsup_Learning_files/figure-markdown_github/unnamed-chunk-1-1.png)

``` r
km <- kmeans(x, 2, 20)
km_3 <- kmeans(x, 3, 50)
```

Use the kmeans() function setting k to 2 and nstart=20 Inspect/print the results Q. What is in the output object?

``` r
attributes(km)
```

    ## $names
    ## [1] "cluster"      "centers"      "totss"        "withinss"    
    ## [5] "tot.withinss" "betweenss"    "size"         "iter"        
    ## [9] "ifault"      
    ## 
    ## $class
    ## [1] "kmeans"

Q. How many points are in each cluster? 30 and 30 Q. What ‘component’ of your result object details - cluster size?

``` r
km$size
```

    ## [1] 30 60

      - cluster assignment/membership?

``` r
km$cluster
```

    ##  [1] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
    ## [36] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1
    ## [71] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

      - cluster center?    

``` r
km$centers
```

    ##           x         y
    ## 1  3.259821 -3.309319
    ## 2 -1.684162  1.600408

Plot x colored by the kmeans cluster assignment and add cluster centers as blue points

``` r
km_df <- data.frame(x=x[ ,1], y=x[ ,2], k=km$cluster)

plot(x=km_df$x, y=km_df$y, col=km_df$k)
```

![](Lecture_9-Unsup_Learning_files/figure-markdown_github/unnamed-chunk-6-1.png)

![](Lecture_9-Unsup_Learning_files/figure-markdown_github/unnamed-chunk-7-1.png)

Hierarchical clustering with distance matrix, hclust, and dendrogram
--------------------------------------------------------------------

![](Lecture_9-Unsup_Learning_files/figure-markdown_github/unnamed-chunk-8-1.png)

    ## 
    ##  1  2  3 
    ## 29 32 29

    ## [1] 32 29 29

![](Lecture_9-Unsup_Learning_files/figure-markdown_github/unnamed-chunk-9-1.png)![](Lecture_9-Unsup_Learning_files/figure-markdown_github/unnamed-chunk-9-2.png)

Q. Use the dist(), hclust(), plot() and cutree() functions to return 2 and 3 clusters

![](Lecture_9-Unsup_Learning_files/figure-markdown_github/unnamed-chunk-10-1.png)

    ##   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 2
    ##  [36] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
    ##  [71] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1
    ## [106] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 1 2 2
    ## [141] 2 2 2 2 2 2 2 2 2 1

    ##   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 2
    ##  [36] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
    ##  [71] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 1
    ## [106] 3 3 3 3 3 3 3 3 3 3 3 3 2 3 3 2 3 2 3 3 3 3 2 3 3 2 2 1 2 2 3 2 1 3 3
    ## [141] 3 3 2 3 3 3 3 3 3 1

Q. How does this compare to your known 'col' groups?

![](Lecture_9-Unsup_Learning_files/figure-markdown_github/unnamed-chunk-11-1.png)![](Lecture_9-Unsup_Learning_files/figure-markdown_github/unnamed-chunk-11-2.png)![](Lecture_9-Unsup_Learning_files/figure-markdown_github/unnamed-chunk-11-3.png)

    ## 
    ##  1  2  3 
    ## 52 62 36

    ## 
    ##  1  2  3 
    ## 50 50 50

PCA analysis
------------

Q1. How many rows and columns are in your new data frame named x? What R functions could you use to answer this questions?

    ## [1] 17  5

    ##                      X England Wales Scotland N.Ireland
    ## 1               Cheese     105   103      103        66
    ## 2        Carcass_meat      245   227      242       267
    ## 3          Other_meat      685   803      750       586
    ## 4                 Fish     147   160      122        93
    ## 5       Fats_and_oils      193   235      184       209
    ## 6               Sugars     156   175      147       139
    ## 7      Fresh_potatoes      720   874      566      1033
    ## 8           Fresh_Veg      253   265      171       143
    ## 9           Other_Veg      488   570      418       355
    ## 10 Processed_potatoes      198   203      220       187
    ## 11      Processed_Veg      360   365      337       334
    ## 12        Fresh_fruit     1102  1137      957       674
    ## 13            Cereals     1472  1582     1462      1494
    ## 14           Beverages      57    73       53        47
    ## 15        Soft_drinks     1374  1256     1572      1506
    ## 16   Alcoholic_drinks      375   475      458       135
    ## 17      Confectionery       54    64       62        41

    ##                  X         England           Wales           Scotland     
    ##  Alcoholic_drinks : 1   Min.   :  54.0   Min.   :  64.0   Min.   :  53.0  
    ##  Beverages        : 1   1st Qu.: 156.0   1st Qu.: 175.0   1st Qu.: 147.0  
    ##  Carcass_meat     : 1   Median : 253.0   Median : 265.0   Median : 242.0  
    ##  Cereals          : 1   Mean   : 469.6   Mean   : 503.9   Mean   : 460.2  
    ##  Cheese           : 1   3rd Qu.: 685.0   3rd Qu.: 803.0   3rd Qu.: 566.0  
    ##  Confectionery    : 1   Max.   :1472.0   Max.   :1582.0   Max.   :1572.0  
    ##  (Other)          :11                                                     
    ##    N.Ireland     
    ##  Min.   :  41.0  
    ##  1st Qu.: 135.0  
    ##  Median : 209.0  
    ##  Mean   : 429.9  
    ##  3rd Qu.: 586.0  
    ##  Max.   :1506.0  
    ## 

    ##                X England Wales Scotland N.Ireland
    ## 1         Cheese     105   103      103        66
    ## 2  Carcass_meat      245   227      242       267
    ## 3    Other_meat      685   803      750       586
    ## 4           Fish     147   160      122        93
    ## 5 Fats_and_oils      193   235      184       209
    ## 6         Sugars     156   175      147       139

    ##                England Wales Scotland N.Ireland
    ## Cheese             105   103      103        66
    ## Carcass_meat       245   227      242       267
    ## Other_meat         685   803      750       586
    ## Fish               147   160      122        93
    ## Fats_and_oils      193   235      184       209
    ## Sugars             156   175      147       139

    ## [1] 17  4

    ## [1] 17  4

Q2. Which approach to solving the ‘row-names problem’ mentioned above do you prefer and why? Is one approach more robust than another under certain circumstances? Ascribing row names inside the read.csv function is shorter and more robust, since the x\[ ,-1\] argument could be used multiple times and delete useful data from x.

![](Lecture_9-Unsup_Learning_files/figure-markdown_github/unnamed-chunk-14-1.png)

Q3: Changing what optional argument in the above barplot() function results in the following plot? change "beside=" argument to FALSE

![](Lecture_9-Unsup_Learning_files/figure-markdown_github/unnamed-chunk-15-1.png)

Q5: Generating all pairwise plots may help somewhat. Can you make sense of the following code and resulting figure? What does it mean if a given point lies on the diagonal for a given plot?

This function produces scatterplots of the values of each food for each pair of countries. If a value lies on the diagonal, it means that the value is equal for both countries.

![](Lecture_9-Unsup_Learning_files/figure-markdown_github/unnamed-chunk-16-1.png)

Q6. What is the main differences between N. Ireland and the other countries of the UK in terms of this data-set? Pairwise comparisons of England, Wales, and Scotland give datapoints that lie very close to the diagonal. The comparisons with N. Ireland show a few datapoints that differ significantly, with higher consumption of fresh potatos and lower consumption of fresh fruit. However, a holistic comparison is difficult.

    ## Importance of components:
    ##                             PC1      PC2      PC3       PC4
    ## Standard deviation     324.1502 212.7478 73.87622 4.189e-14
    ## Proportion of Variance   0.6744   0.2905  0.03503 0.000e+00
    ## Cumulative Proportion    0.6744   0.9650  1.00000 1.000e+00

Q7. Complete the code below to generate a plot of PC1 vs PC2. The second line adds text labels over the data points.

![](Lecture_9-Unsup_Learning_files/figure-markdown_github/unnamed-chunk-18-1.png)

Q8. Customize your plot so that the colors of the country names match the colors in our UK and Ireland map and table at start of this document.

![](Lecture_9-Unsup_Learning_files/figure-markdown_github/unnamed-chunk-19-1.png)

![](Lecture_9-Unsup_Learning_files/figure-markdown_github/unnamed-chunk-20-1.png)
