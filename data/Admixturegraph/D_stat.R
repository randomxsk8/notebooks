#---------------------------
'''
ADMIXTUREGRAPH

In this notebook there is a try to fitting the ABBA-BABA Results with the Treemix graph of the 3L-free region

'''

install.packages("admixturegraph")
setwd('/home/randomx/admixturegraph/')
library("admixturegraph")

### BUILDING THE GRAPH
leaves <- c("chri", "KE", "FRgam","GAgam","GQgam","UGgam", "BFgam", "GNgam", "GHgam", "CMgam", "AOcol","GNcol", "BFcol", "CIcol", "GHcol", "GM", "GW")
inner_nodes <- c("R", "r", "u", "k","g", "gg", "c", "cc")
edges <- parent_edges(c(edge("chri", "r"),
                        edge("FRgam", "k"),
                        edge("KE", "k"),
                        edge("GAgam", "g"),
                        edge("GQgam", "g"),
                        edge("UGgam", "gg"),
                        edge("BFgam", "gg"),
                        edge("GNgam", "gg"),
                        edge("GHgam", "gg"),
                        edge("CMgam", "gg"),
                        edge("AOcol", "c"),
                        edge("GNcol", "cc"),
                        edge("BFcol", "cc"),
                        edge("CIcol", "cc"),
                        edge("GHcol", "cc"),
                        edge("r", "R"),
                        edge("u", "R"),
                        edge("k","u"),
                        edge("g", "u"),
                        edge("gg", "g"),
                        edge("c","u"),
                        edge("cc", "c")
))
leaf <- c("chri", "KE", "FRgam","GAgam","GQgam","UGgam", "BFgam", "GNgam", "GHgam", "CMgam", "AOcol", "GNcol", "BFcol", "CIcol", "GHcol")
an_graph <- agraph(leaves, inner_nodes, edges)

### PLOTTING THE GRAPH
plot(an_graph, leaf_order = leaf) 

## LOADING ABBA-BABA Results file
d_stat <- read.csv('ABBA-BABA_Kenya_chrom_3.txt', sep="\t")
plot(f4stats(d_stat))

## FITTING THE TREE WITH THE Results
an_fit <- fit_graph(d_stat, an_graph)
summary(an_fit)
plot(an_fit)

### MCM model?
mcmc1 <- make_mcmc_model(an_graph, d_stat)
initial1 <- rep(0.5, length(mcmc1$parameter_names))
chain1 <- run_metropolis_hasting(mcmc1, initial1, iterations = 10000, verbose = FALSE)
head(chain1)
plot(chain1[, "edge_u_k"])
thinned_1 <- thinning(burn_in(chain1, 4000), 100)
plot(thinned_1[, "edge_u_k"])
hist(thinned_1[, "edge_u_k"])
model_likelihood(thinned_1[, "likelihood"])
model_likelihood_n(thinned_1[, "likelihood"], 100)
  
#-----------------------------------------------------