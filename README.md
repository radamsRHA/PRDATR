
---
# PRDATR: PRobabilistic Distances under models of Adaptive Trait evolution in R
**NOTE See the file https://github.com/radamsRHA/PRDATR/blob/master/PRDATR_Manual.pdf for detailed instructions**

## Installing R package PRDATR from github
The R package PRDATR is freely available to download and distribute from github <https://github.com/radamsRHA/PRDATR/>. To install and load PRDATR, you must first install the R package `devtools`, Additionally, make sure the most updated version of R is installed 

```
install.packages("devtools")
```
Now using devtools we can install `PRDATR` from github:

```
library(devtools)
install_github("radamsRHA/PRDATR")
library(PRDATR) # Load package 
```
`PRDATR` also requires the following dependencies to be installed:

* install.packages('geiger')  
* install.packages('ape')  
* install.packages('gaussDiff')  


To begin using `PRDATR` try using the examples associated with each function. 

## Example: compute probabilistic distances between evolutionary models  

We can use `Function_ComputeDistances` to compute the probabilistic distance between a BM and an OU model for an example dataset. First, let's load the R package `PRDATR` and its dependancies:

```
################
# Load depends #
################
library(ape)
library(geiger)
library(gaussDiff)
library(PRDATR)
```

Now, let's specify the phylogenetic tree we want to use as a framework (here we are using the tree from Felsenstein 1985)

```
############################################################################################
# Specifity example tree (Fig. 8, Felsenstein 1985) used for demonstrating model distances #
############################################################################################
string.Figure01_Felsenstein1985_Tree <- "(((Species_7:1.635983031,Species_8:0.8079384444):1.801052391,(Species_6:0.4510394335,Species_5:1.208543249):0.4146682462):0.4476024657,((Species_4:0.9434278477,Species_3:0.2480806489):2.642993143,(Species_2:2.588686978,Species_1:0.908678702):0.3223280154):1.399150111):1;"
handle.Figure01_Felsenstein1985_Tree <- read.tree(text = string.Figure01_Felsenstein1985_Tree)
```

Next, we specify two lists, each representing the two models:

```
#################################
# Construct a list for Model 01 #
#################################
vector.Model_01_Theta <- c(1) # using a rate of 1
names(vector.Model_01_Theta) <- c("Sig2")

list.Model01_BM <- list(handle.Phylogeny = handle.Figure01_Felsenstein1985_Tree,
                        string.Model = "BM",
                        vector.Z = rep(0, length(handle.Figure01_Felsenstein1985_Tree$tip.label)), # mean (ancestral root value) for all tips in the tree
                        vector.Theta = vector.Model_01_Theta)
                        


#################################
# Construct a list for Model 02 #
#################################
vector.Model_02_Theta <- c(1, 1) # rate of 1, alpha of 1
names(vector.Model_02_Theta) <- c("Sig2", "alpha")

list.Model02_OU <- list(handle.Phylogeny = handle.Figure01_Felsenstein1985_Tree,
                        string.Model = "OU",
                        vector.Z = rep(0, length(handle.Figure01_Felsenstein1985_Tree$tip.label)), # mean (ancestral root value) for all tips in the tree
                        vector.Theta = vector.Model_02_Theta)

```

After the two models have been specified in their respective lists (`Model01_BM` and `list.Model02_OU`), we can compute the probabilistic distances between the two:

```
Function_ComputeDistances(list.Model_01 = list.Model01_BM, list.Model_02 = list.Model02_OU)

       dH        dKL 
 0.9991598 15.5200103 
```

As you can see, the two distances (dH and dKL) are 0.9991598 and 15.5200103, respectively. In practice, the specific values of the model parameters can be specified as shown above, to represent maximum likelihood or Bayesian inferences, for example. 
