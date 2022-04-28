# Poolseq replicate variance
## Evan Durland
## 4/28/2022

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
```

## A visualization of variance in estimates of minor allele frequency between independent culture replicates

One of the central assumptions of analyses using pooled DNA sequencing is that the allele frequencies are representative of the broader population. This assumption is vulnerable to several sources of error from sampling bias to PCR duplications to sequencing error. 

Ideally, technical replicates should be taken from any single population in order to ascertain the extent of the technical variation with this approach. In practice this may not be feasible due to constraints on individual samples or cost of analyses (library prep and sequencing).

In our recent study (https://doi.org/10.1098/rspb.2020.3223) we used pooled DNA samples from five replicate cultures of larval oysters to evaluate changes in allele frequencies across larval development. Lacking technical replicates begs the question whether estimates of minor allele frequencies (MAF) in each culture are accurate to the 'true' value of that population. 

In this scenario, we analyzed variance in MAF across all five biological replicates at several time points. This approach folds the technical variation of the pooled DNA sequencing approach in with the biological variance between replicate cultures. Any variance between the culture replicates at a given sampling (time) point, then, encompasses natural biological differences as well as variance arising from unintentional technical error. 

This script explores the nature of this variation in our dataset to more fully describe the significance of the dynamic trajectories of allele frequency we observed.

```{r}
#load the data:
df <- read.delim("will_avg_sepAF.txt", sep="\t")
head(df)
```

this dataframe contains the following information:
   $loc: the ID of the SNP locus
   $af_eggs --> $afspat: mean allele frequency for the 6 sampling points (eggs, D-larvae, day 6, day 10, eyed, spat)
   $flop2: categorical type of change ("none"=gradual, "linear"=uni-directional and "flip"= bi-directional)
   $clust: kmeans cluster assignment.
   $WA1.D2af --> WA5.D22Baf: alle frequencies of indiviual cultures (WA1-WA5) at each sampling point(D2-D22).

for plotting, we have to select the allele frequency (af) from each replicate (e.g. WA1, WA2)

```{r}
af_cols <- grep("af$",colnames(df))

#keep the 'loc' columns for plotting and pivot the rest to a melted dataframe:
m_af <- df[,c(1,af_cols)]%>%
  pivot_longer(,cols=c(2:23),names_to = "sample",values_to="af")%>%
  as.data.frame()
head(m_af)

#designate ages:
m_af$age<-apply(m_af,1,function(x){
  smpl <- x[2]
  day <- strsplit(smpl,split="\\.")%>%
    unlist()%>%
    .[2]
  return(day)
})
head(m_af)
m_af$age <- factor(m_af$age,levels=c("D2af","D6af","D10af","D16af","D22Baf"))
```
The entire dataset is **473** loci with **10 406** datapoints in total.
Much too large to plot in full, individually.
We can evaluate the nature of the variation in small random samples:
First create a subset of the loci:

```{r}
# a random selection:
rand <- sample(unique(df$loc),16, replace=FALSE)

# alternatively, we can select loci with varying levels of change.  
# In the original analysis we used a k-means clustering algorighm to group
# loci by their trajectory of change in allele frequency.  
# These clusters (n=5) break down as follows: 
xtabs(~clust,data=df)
# cluster 'name' (1-5) is randomly allocated in the algorithm 
# but we know from the original analysis that the size of each cluster 
# is inversely correlated to it's 'dynamic' pattern of change (e.g. more loci in more boring patterns) 

# if we want to select the most extreme examples:
xtrm <- df %>%
  filter(clust==4)%>%
  .$loc %>%
  sample(.,15,replace=FALSE)

# or their boring counterparts:

boring <- df %>%
  filter(clust==1)%>%
  .$loc %>%
  sample(.,16, replace=FALSE)
```
From here we can plot what the variance looks like:

```{r}
#plot:
m_af%>%
  filter(loc%in% rand)%>%
  ggplot(.,aes(age,af))+
  geom_point()+
  stat_smooth(aes(group=loc))+
  facet_wrap(~loc)

# In the above code chunk, changing the filter on the dataset (rand, xtrm or boring)
# will plot the associated random sample
```
![](https://github.com/E-Durland/poolseq_variance/blob/main/rep_var.png)
