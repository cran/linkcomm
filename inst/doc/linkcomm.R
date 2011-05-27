### R code from vignette source 'linkcomm.Rnw'

###################################################
### code chunk number 1: linkcomm.Rnw:56-58 (eval = FALSE)
###################################################
## source("load_linkcomm.R")
## source("load_OCG.R")


###################################################
### code chunk number 2: linkcomm.Rnw:65-66 (eval = FALSE)
###################################################
## install.packages("linkcomm")


###################################################
### code chunk number 3: linkcomm.Rnw:71-72
###################################################
library(linkcomm)


###################################################
### code chunk number 4: linkcomm.Rnw:77-78 (eval = FALSE)
###################################################
## demo(topic = "linkcomm", package = "linkcomm")


###################################################
### code chunk number 5: linkcomm.Rnw:92-93
###################################################
head(lesmiserables)


###################################################
### code chunk number 6: linkcomm.Rnw:109-111 (eval = FALSE)
###################################################
## yeast_pp <- read.table("pp_rnapol.txt", header = FALSE)
## head(yeast_pp)


###################################################
### code chunk number 7: linkcomm.Rnw:114-117
###################################################
load(file = "networks.rda")
zz<-data.frame(V1=pp_rnapol[,1],V2=pp_rnapol[,2])
head(zz)


###################################################
### code chunk number 8: linkcomm.Rnw:122-123 (eval = FALSE)
###################################################
## lc <- getLinkCommunities(yeast_pp, hcmethod = "single")


###################################################
### code chunk number 9: linkcomm.Rnw:126-127
###################################################
cat("   Removing loops...\n   Removing edge duplicates...\n   Calculating edge similarities for 449 edges... 100.00%\n   Hierarchical clustering of edges...\n   Calculating link densities... 100.00%\n   Partition density maximum height =  0.33333\n   Finishing up...4/4... 100%\n   Plotting...\n   Colouring dendrogram... 100%\n")


###################################################
### code chunk number 10: linkcomm.Rnw:141-142
###################################################
lc <- getLinkCommunities(pp_rnapol, hcmethod = "single")


###################################################
### code chunk number 11: linkcomm.Rnw:145-146
###################################################
print(lc)


###################################################
### code chunk number 12: linkcomm.Rnw:175-177 (eval = FALSE)
###################################################
## plot(lc, type = "graph", layout = layout.fruchterman.reingold)
## plot(lc, type = "graph", layout = "spencer.circle")


###################################################
### code chunk number 13: linkcomm.Rnw:182-183 (eval = FALSE)
###################################################
## plot(lc, type = "graph", layout = "spencer.circle", shownodesin = 3)


###################################################
### code chunk number 14: linkcomm.Rnw:188-189 (eval = FALSE)
###################################################
## plot(lc, type = "graph", shownodesin = 2, node.pies = TRUE)


###################################################
### code chunk number 15: linkcomm.Rnw:195-196 (eval = FALSE)
###################################################
## plot(lc, type = "members")


###################################################
### code chunk number 16: linkcomm.Rnw:201-202 (eval = FALSE)
###################################################
## plot(lc, type = "summary")


###################################################
### code chunk number 17: linkcomm.Rnw:207-208 (eval = FALSE)
###################################################
## plot(lc, type = "dend")


###################################################
### code chunk number 18: linkcomm.Rnw:216-217
###################################################
plot(lc, type = "members", nodes = head(names(lc$numclusters),10))


###################################################
### code chunk number 19: linkcomm.Rnw:229-230
###################################################
nc <- getAllNestedComm(lc)


###################################################
### code chunk number 20: linkcomm.Rnw:233-234 (eval = FALSE)
###################################################
## getAllNestedComm(lc)


###################################################
### code chunk number 21: linkcomm.Rnw:237-238
###################################################
nc


###################################################
### code chunk number 22: linkcomm.Rnw:243-244 (eval = FALSE)
###################################################
## getNestedHierarchies(lc, clusid = 9)


###################################################
### code chunk number 23: linkcomm.Rnw:249-250 (eval = FALSE)
###################################################
## plot(lc, type = "graph", clusterids = c(9,11))


###################################################
### code chunk number 24: linkcomm.Rnw:257-258
###################################################
plot(lc, type = "graph", clusterids = c(9,11), frame = TRUE)


###################################################
### code chunk number 25: linkcomm.Rnw:268-269
###################################################
cr <- getClusterRelatedness(lc, hcmethod = "ward")


###################################################
### code chunk number 26: linkcomm.Rnw:272-273 (eval = FALSE)
###################################################
## cr <- getClusterRelatedness(lc, hcmethod = "ward")


###################################################
### code chunk number 27: linkcomm.Rnw:292-293 (eval = FALSE)
###################################################
## cutDendrogramAt(cr, cutat = 1.2)


###################################################
### code chunk number 28: linkcomm.Rnw:296-297
###################################################
tt <- cutDendrogramAt(cr, cutat = 1.2)


###################################################
### code chunk number 29: linkcomm.Rnw:300-301
###################################################
print(tt)


###################################################
### code chunk number 30: linkcomm.Rnw:320-321
###################################################
cc <- getCommunityCentrality(lc)


###################################################
### code chunk number 31: linkcomm.Rnw:324-325 (eval = FALSE)
###################################################
## cc <- getCommunityCentrality(lc)


###################################################
### code chunk number 32: linkcomm.Rnw:328-329
###################################################
head(sort(cc, decreasing = TRUE))


###################################################
### code chunk number 33: linkcomm.Rnw:335-336
###################################################
head(lc$numclusters)


###################################################
### code chunk number 34: linkcomm.Rnw:351-353 (eval = FALSE)
###################################################
## cm <- getCommunityConnectedness(lc, conn = "modularity")
## plot(lc, type = "commsumm", summary = "modularity")


###################################################
### code chunk number 35: linkcomm.Rnw:360-361
###################################################
plot(lc, type = "commsumm", summary = "modularity", verbose = FALSE)


###################################################
### code chunk number 36: linkcomm.Rnw:372-373
###################################################
lc2 <- newLinkCommsAt(lc, cutat = 0.4)


###################################################
### code chunk number 37: linkcomm.Rnw:376-377 (eval = FALSE)
###################################################
## lc2 <- newLinkCommsAt(lc, cutat = 0.4)


###################################################
### code chunk number 38: linkcomm.Rnw:382-383
###################################################
print(lc2)


###################################################
### code chunk number 39: linkcomm.Rnw:391-392
###################################################
getNodesIn(lc, clusterids = c(4,5))


###################################################
### code chunk number 40: linkcomm.Rnw:402-403 (eval = FALSE)
###################################################
## lc <- getLinkCommunities(yeast_pp, directed = TRUE, dirweight = 0.8)


###################################################
### code chunk number 41: linkcomm.Rnw:408-409
###################################################
head(weighted)


###################################################
### code chunk number 42: linkcomm.Rnw:427-428 (eval = FALSE)
###################################################
## lc <- getLinkCommunities(yeast_pp, edglim = 10)


###################################################
### code chunk number 43: linkcomm.Rnw:438-439 (eval = FALSE)
###################################################
## linkcomm2cytoscape(lc, interaction = "pp", ea = "linkcomms.ea")


###################################################
### code chunk number 44: linkcomm.Rnw:457-458 (eval = FALSE)
###################################################
## lm <- getLinkCommunities(lesmiserables, plot = FALSE)


###################################################
### code chunk number 45: linkcomm.Rnw:463-466 (eval = FALSE)
###################################################
## nf <- graph.feature(lm, type = "nodes", indices = which(V(lm$igraph)$name == "Valjean"), 
##                     features = 30, default = 5)
## plot(lm, type = "graph", vsize = nf, vshape = "circle", shownodesin = 4)


###################################################
### code chunk number 46: linkcomm.Rnw:471-474 (eval = FALSE)
###################################################
## nf <- graph.feature(lm, type = "nodes", indices = getNodesIn(lm, clusterids = 1, 
##                     type = "indices"), features = 30, default = 5)
## plot(lm, type = "graph", vsize = nf, vshape = "circle", vlabel = FALSE)


###################################################
### code chunk number 47: linkcomm.Rnw:479-482 (eval = FALSE)
###################################################
## ef <- graph.feature(lm, type = "edges", indices = getEdgesIn(lm, clusterids = 14), 
##                     features = 5, default = 1)
## plot(lm, type="graph", ewidth = ef)


###################################################
### code chunk number 48: linkcomm.Rnw:487-492 (eval = FALSE)
###################################################
## ef <- graph.feature(lm, type = "edges", indices = getEdgesIn(lm, nodes = "Myriel"), 
##                     features = 5, default = 1)
## nf <- graph.feature(lm, type = "nodes", indices = which(V(lm$igraph)$name == "Myriel"), 
##                     features = 30, default = 5)
## plot(lm, type = "graph", vsize = nf, ewidth = ef, vshape = "circle", vlabel = FALSE)


###################################################
### code chunk number 49: linkcomm.Rnw:508-509
###################################################
oc <- getOCG.clusters(lesmiserables)


###################################################
### code chunk number 50: linkcomm.Rnw:512-513
###################################################
cat("Calculating Initial class System....Done\nNb. of classes 43\nNb. of edges not within the classes 19\nNumber of initial classes 43\nRunning....\nRemaining classes: None\nReading OCG data...\nExtracting cluster sizes... 100%\n")


###################################################
### code chunk number 51: linkcomm.Rnw:518-519
###################################################
print(oc)


###################################################
### code chunk number 52: linkcomm.Rnw:524-525 (eval = FALSE)
###################################################
## plot(oc, type = "graph", shownodesin = 7, scale.vertices = 0.1)


