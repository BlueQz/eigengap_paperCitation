
###############################################################################
#
#  This is the R code to reprocdue the results in the comment by Song Wang 
# and Karl Rohe titled  "Don't Mind the (Eigen) Gap", which is an invited comment
# on the paper titled "Coauthorship and Citation Networks for Statisticians", 
# to appear in the Annals of Applied Statistics.
#
################################################################################




rm(list=ls())   # clean up the space

## load the required packages

library(igraph)  # graph visualization
library(tm) # text mining
library(Matrix)  ## sparse matrix computation
library(rARPACK)  ## fast eigenvector computation
library(irlba)   ## fast svd computation
library(xtable) #generate table in latex code


##read data into R

p2p = as.matrix(read.table("../data/paperCitAdj.txt", header=F)) ## paper cite paper
p2p = Matrix(p2p)  # transforms to sparse matrix
absk = read.table("../data/paperList_Abstracts_Keyword.txt", sep = ",",
                  colClass = c("character","integer", "character", "character"), 
                  header = T) # paper - abstract - keywords


## processing on the graph:  graph with core >= 4
g = graph.adjacency(adjmatrix = p2p,mode = "undire")
V(g)$name = 1:length(V(g))
core = coreness(g,mode = "all")
coreID4 <- which( core >= 4 )
g2 = induced_subgraph(g, coreID4)
A = get.adjacency(g2)


## bag of words matrix text analysis
text <- absk[coreID4,3]  # use specific set
vc <- VCorpus( VectorSource(text) ) # just change the vector of strings to corpus
ctrl <- list(removePunctuation = list(preserve_intra_word_dashes = TRUE)
             , stopwords = TRUE
             ,removeNumbers = TRUE
             #, stemming = TRUE                    # remove prefix or postfix
             #, bounds = list(global= c(15,Inf))   # remove low-frequency/high frequency words
             #, wordLengths = c(4, 20) # remove short words like 'a' 
             #, weighting = function(x) weightSMART(x, spec = "nnn")
)
tdm <- TermDocumentMatrix(vc, control =  ctrl)  ## term - ducoment matrix
terms <- tdm$dimnames$Terms
print ( sprintf( "after initial cleaning: %s words remains in %s docuemnts",
     dim(tdm)[1], dim(tdm)[2], '\n') )  


B = spMatrix(i = tdm$i, j = tdm$j, x = tdm$v, nrow = tdm$nrow, ncol  = tdm$ncol)         # frequency count
rownames(B)  = tdm$dimnames$Terms

# remove 's' 
# for word ends with 's', whether the word without 's' is in terms. 
# like designs is in, check the posiition of design, all the locations of design alike are returned
# some are NA, means like "boss" exists, but "bos" not.
idx <- match( gsub(pattern = '(.*)s$', replacement = '\\1', x= terms[grep('s$',terms)]), terms)
idx1 <- match(paste0(terms[idx[!is.na(idx)]],'s'),terms)    # location of plural terms
idx2 <- match(terms[idx[!is.na(idx)]], terms)   #location of single terms with out s
B[idx1,] <- B[idx1,]+B[idx2,]
terms <- terms[-idx1];  B<- B[terms,]; #update terms, tdm

# remvoe 'ed'
idx <- match( gsub(pattern = '(.*)ed$', replacement = '\\1', x= terms[grep('ed$',terms)]), terms)
idx1 <- match(paste0(terms[idx[!is.na(idx)]],'ed'),terms)
idx2 <- match(terms[idx[!is.na(idx)]], terms)
B[idx1,] <- B[idx1,]+B[idx2,]
terms <- terms[-idx1];  B<- B[terms,]; #update terms, tdm
print (sprintf( "after combining 's','ed' cleaning: %s words remains in %s docuemnts",
         dim(B)[1], dim(B)[2], '\n') )


## remove words that appears in less than 10 document
B = (B>0)+0  # converte counts in term document matrix to  {0,1}
rownames(B)  = terms
kw = B
kw = kw[rowSums(kw)>=10,]   
print(dim(kw))
kw = t(kw)  # document term matrix



### Partioning the graph by Regularized spectral clustering with 11 clusters

nr = nrow(A); nc = ncol(A)
rs=rowSums(A); cs = colSums(A)
E = sum(A)
tau = sqrt(E/nr)  

Drow = Diagonal(n = nr, x = 1/sqrt(rs + tau))
Dcol = Diagonal(n = nc, x = 1/sqrt(cs + tau))
tmp = Drow%*%A
L = tmp %*% Dcol
e = irlba(L,nu = 50)

## Fig 1
#par(mfcol = c(1,1), mar = c(5,4,4,2)+0.1)
#pdf('../latex/imsart-ims/screeplot.pdf',width = 5,height = 5)
par(mfcol = c(1,1),mar=c(3,3,1,1))
plot(1:25,e$d[1:25], cex=2, pch=15, type="b", lty=2, 
     lwd=3, col="red", bg="red", xlab="", ylab="")
#dev.off()

e1 = eigen(L,symmetric = TRUE)
#pdf('../latex/imsart-ims/histogram.pdf',width = 5,height = 5)
par(mfcol = c(1,1),mar=c(3,3,1,1))
hist(e1$values,breaks = 50, main ="",xlab = "", ylab="", probability = T)
lines(density(e1$values, bw =0.05),col="red")
#dev.off()


nclust <- 11 
#nclust <- 20
X <- e$u[,1:nclust]
nX= t(apply(X, 1, function(x,nr) return(x/sqrt(sum(x^2)+length(x)/(100*nr))), nr = nr))
set.seed(10)
km_11 <- kmeans(nX,nclust,iter.max = 1000, nstart = 1000)
clust <- km_11$cluster
sort(table(clust))  

#table 1, summerize the clustering results
hv <- Matrix(0,nrow(nX), nclust)
for(i in 1:nclust) {hv[clust==i,i] <- 1}
ord <- order(km_11$size,decreasing = T)
hv <- hv[,ord]
cat("volumn within and cross clusters:","\n")
vol_clust <- t(hv)%*%A%*%hv
print(vol_clust)
print(cbind(diag(vol_clust),rowSums(vol_clust)-diag(vol_clust)))


## moving to interpret clusters where ncluster = 11,

bestWords = matrix("", nrow = nclust, ncol = 20)
#pdf('../latex/imsart-ims/bestWords.pdf', width = 3.5, height =3.5)
for(i in 1:nclust){
    clust_i = which(hv[,i]==1)
    # variance statblizing 
    vals =  asin(sqrt(colMeans(kw[clust_i,])))- asin(sqrt(colMeans(kw[-clust_i,])))
    bestWords[i,] = (names(vals)[order(-vals)[1:20]])
}
bestWords[,1:10]
#dev.off()


#Table 2
topwords <- apply(bestWords[,1:10],MARGIN = 1, function(x) paste0(x,collapse = ", "))
results <- list()
results$id = 1:nclust
name11 <- c("Vari Selection","Mutiple Testing", "Semi/Non Para", "Functional Data", 
                 "Cov Matrix",'Sliced Inverse Regr', "Spatial", "Classification",
                 "Bayesian","Learn Theory",'Den Estimation')
if(nclust == 11)  results$name <- name11
results$size <-  km_11$size[ord]
results$In <- as.integer(diag(vol_clust))
results$Out <- as.integer(rowSums(vol_clust)-diag(vol_clust))
results$topwords <- topwords[ord]                      
results <- data.frame(results,row.names = NULL)  

#generate latex file for tables
if(nclust ==11){
    print(xtable(    cbind(results[1:5,1:5],results[7:11,1:5])    ))
    print(xtable(    results[,c(2,6)]    ))
}


## Fig 2: visualization of clusters with 11 clusters
if(nclust == 11){
    colours <- c("red","yellow","green","blue","pink",'brown',
                 "cyan","grey","darkcyan", 'deepskyblue','black')
    vColor <- character(vcount(g2))      
    for( i in 1:nclust){
        vColor[which(clust == i)] <- colours[i]
    }
    
    #pdf('../latex/imsart-ims/graph_visualization_11.pdf', width = 14, height = 10 )
    set.seed(10)
    plot(g2, vertex.color = vColor, 
         vertex.label = "",
         #vertex.label.cex=3, vertex.label.color='red', 
         #vertex.label.dist = 0,
         vertex.size=sqrt(degree(g2)),  asp=.7, 
         layout = layout.auto,    # same as fruchterman.reingold
         # equivalent to <100, layout.kamada.kawai,
         # layout.fruchterman.reingold, 
         # layout.drl
         edge.color='grey',edge.width = 0.25,  margin=c(0,0,0,0.5)
    )
    ord <- order(table(clust),decreasing = T)
    legend(x= c(1,1), y = c(1,1),legend = name11,  col = colours[ord], pch = 16, cex = 1.5)
    #dev.off()
}


### The results when asssuming there are 20 clusters

nclust <- 20
X <- e$u[,1:nclust]
nX= t(apply(X, 1, function(x,nr) return(x/sqrt(sum(x^2)+length(x)/(100*nr))), nr = nr))
set.seed(10)
km_20 <- kmeans(nX,nclust,iter.max = 1000, nstart = 100)
#membership matrix hv
hv <- Matrix(0,nrow(nX), nclust)
for(i in 1:nclust) { hv[km_20$cluster ==i,i] <- 1 }
ord <- order(km_20$size,decreasing = T)
hv <- hv[,ord]
cat("volumn within and cross clusters:","\n")

bestWords = matrix("", nrow = nclust, ncol = 20)
#pdf('../latex/imsart-ims/bestWords.pdf', width = 3.5, height =3.5)
for(i in 1:nclust){
    clust_i = which(hv[,i]==1)
    #variance stablizing 
    vals =  asin(sqrt(colMeans(kw[clust_i,])))- asin(sqrt(colMeans(kw[-clust_i,])))
    bestWords[i,] = (names(vals)[order(-vals)[1:20]])
}
#bestWords[,1:10]

topwords <- apply(bestWords[,1:10],MARGIN = 1, function(x) paste0(x,collapse = ", "))
results <- list()
results$id = 1:nclust
name20 <- c("Multiple Testing","Lasso I","FDA","Cov Esitmation",'Dim Reduction',
                 'Lasso II',"Longitudinal", "Forecast(in other fields)","Bayesian nonpara","Non-parametric I",
                 "Bayes(others)",'Spatial(bayesian)','Quantile regreesion','Learning Theory I', 'Learning Theory II',
                 'Classification','Non-parametric II','Spatial II(frequentist)','Designs','Semiparametric')

#table 3,4
results$name <- name20

vol_clust <- t(hv)%*%A%*%hv
print(vol_clust)
print(cbind(diag(vol_clust),rowSums(vol_clust)-diag(vol_clust)))

results$size <-  km_20$size[ord]
results$In <- as.integer(diag(vol_clust))
results$Out <- as.integer(rowSums(vol_clust)-diag(vol_clust))
results$topwords <- topwords[ord]                      
results <- data.frame(results,row.names = NULL)                     
head(results)
# latex file
print(xtable(    cbind(results[1:10,1:5],results[11:20,1:5]) ))
print(xtable(   results[,c(2,6)]  ))









#################################################################################
##   Extra tables and visualization for nclust = 11 or 20
#################################################################################

# visualization for the core in the whole plot

#pdf("../latex/imsart-ims/graph_visualization_2.pdf", width =12,  height = 6)
Vcolor <- character(nrow(p2p))
Vcolor[] <-"white"
Vcolor[core>3] <- vColor
plot(g, vertex.color = Vcolor, vertex.label= "",
     vertex.label.cex=1.2,
     vertex.size=log(degree(g)+1), vertex.label.color='red', asp=.5,
     layout = layout.drl, #layout.kamada.kawai , #layout.fruchterman.reingold, #
     edge.color='green', edge.width = 0.1, margin=c(0,0,0,0), 
     main = "3 core in the whole plot"
)
#dev.off()


# plot of the sorted adjcency matrix based on clustering results

ord <- order(km_11$size, decreasing =  T)
ord_node <- NULL
for(i in ord){
    ord_node <- c(ord_node, which ( km_11$cluster==i))
}
idx <- which(as.matrix(A)==1,arr.ind = T)
#pdf('./A_plot.pdf',width = 4, height =4)
plot(idx[,1],idx[,2],col= "blue", pch =".", xlab="", ylab="")
#dev.off()
A.sort <- as.matrix(A)[ord_node,ord_node]
idx.sort <- which(A.sort ==1, arr.ind =  1)
#pdf('./sbm_plot.pdf',width = 4, height =4)
plot(idx.sort[,1],idx.sort[,2],col= "blue", pch =".",xlab="", ylab="")
#dev.off()



# cluster interpretaion of 11, 20, for 635 papers in the core >=4.
paper_core4 <- absk[coreID4,]
dict <- data.frame(list(name = name11, number = order(km_11$size,decreasing = T)))
clust_name <- as.character( dict$name[order(dict$number)] )
paper_core4$clust11 <- km_11$cluster
for(i in 1:11){
    paper_core4$clust11[which(km_11$cluster==i)] <- clust_name[i]
} 

paper_core4$clust20 <- km_20$cluster
dict <- data.frame(list(name = name20, number = order(km_20$size,decreasing = T)))
clust_name <- as.character( dict$name[order(dict$number)] )
for (i in 1:20){
    paper_core4$clust20[which(km_20$cluster==i)] <- clust_name[i]
}
# table(paper_core4[c(6,5)])
#first 20 samples
print(paper_core4[1:20,c(4,5,6)])

