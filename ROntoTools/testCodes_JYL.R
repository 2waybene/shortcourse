##=================================================
##  testCodes_JYL.R
### R code from vignette source 'rontotools.Rnw'

###################################################
### code chunk number 1: rontotools.Rnw:55-58
###################################################
require(graph)
require(ROntoTools)
kpg <- keggPathwayGraphs("hsa", verbose = FALSE)


###################################################
### code chunk number 2: rontotools.Rnw:62-63 (eval = FALSE)
###################################################
## kpg <- keggPathwayGraphs("hsa", updateCache = TRUE, verbose = TRUE)

length(keggPathwayNames(organism = "hsa", updateCache = FALSE, verbose = TRUE))
#[1] 337

###################################################
### code chunk number 3: rontotools.Rnw:69-70
###################################################
head(names(kpg))
length(kpg)
#[1] 221
###################################################
### code chunk number 4: rontotools.Rnw:74-77
###################################################
kpg[["path:hsa04110"]]
head(nodes(kpg[["path:hsa04110"]]))
head(edges(kpg[["path:hsa04110"]]))


###################################################
### code chunk number 5: rontotools.Rnw:81-82
###################################################
head(edgeData(kpg[["path:hsa04110"]], attr = "subtype"))


###################################################
### code chunk number 6: rontotools.Rnw:86-90
###################################################
kpg <- setEdgeWeights(kpg, edgeTypeAttr = "subtype", 
                      edgeWeightByType = list(activation = 1, inhibition = -1, 
                                              expression = 1, repression = -1), 
                      defaultWeight = 0)


###################################################
### code chunk number 7: rontotools.Rnw:94-95
###################################################
head(edgeData(kpg[["path:hsa04110"]], attr = "weight"))


###################################################
### code chunk number 8: rontotools.Rnw:99-101
###################################################
kpn <- keggPathwayNames("hsa")
head(kpn)

length(kpn)
#[1] 337
kpn[which(names(kpn) %in% "path:hsa05216")]
###################################################
### code chunk number 9: rontotools.Rnw:109-111
###################################################
load(system.file("extdata/E-GEOD-21942.topTable.RData", package = "ROntoTools"))
head(top)
str(top)

###################################################
### code chunk number 10: rontotools.Rnw:127-136
###################################################
fc <- top$logFC[top$adj.P.Val <= .01]
names(fc) <- top$entrez[top$adj.P.Val <= .01]

pv <- top$P.Value[top$adj.P.Val <= .01]
names(pv) <- top$entrez[top$adj.P.Val <= .01]

head(fc)

head(pv)


###################################################
### code chunk number 11: rontotools.Rnw:140-145
###################################################
fcAll <- top$logFC
names(fcAll) <- top$entrez

pvAll <- top$P.Value
names(pvAll) <- top$entrez


###################################################
### code chunk number 12: rontotools.Rnw:149-151
###################################################
ref <- top$entrez
head(ref)


###################################################
### code chunk number 13: rontotools.Rnw:163-165
###################################################
kpg <- setNodeWeights(kpg, weights = alphaMLG(pv), defaultWeight = 1)
head(nodeWeights(kpg[["path:hsa04110"]]))


###################################################
### code chunk number 14: rontotools.Rnw:177-178
###################################################
peRes <- pe(x = fc, graphs = kpg, ref = ref,  nboot = 200, verbose = FALSE)


###################################################
### code chunk number 15: rontotools.Rnw:182-186
###################################################
head(Summary(peRes))

head(Summary(peRes, pathNames = kpn, totalAcc = FALSE, totalPert = FALSE,
             pAcc = FALSE, pORA = FALSE, comb.pv = NULL, order.by = "pPert"))


###################################################
### code chunk number 16: peRes_twoway1
###################################################
plot(peRes)


###################################################
### code chunk number 17: peRes_twoway2
###################################################
plot(peRes, c("pAcc", "pORA"), comb.pv.func = compute.normalInv, threshold = .01)


###################################################
### code chunk number 18: fig1
###################################################
plot(peRes)


###################################################
### code chunk number 19: fig2
###################################################
plot(peRes, c("pAcc", "pORA"), comb.pv.func = compute.normalInv, threshold = .01)


###################################################
### code chunk number 20: pePathway_twoway_Acc
###################################################
plot(peRes@pathways[["path:hsa05216"]], type = "two.way")


###################################################
### code chunk number 21: pePathway_boot_Acc
###################################################
plot(peRes@pathways[["path:hsa05216"]], type = "boot")


###################################################
### code chunk number 22: fig3
###################################################
plot(peRes@pathways[["path:hsa05216"]], type = "two.way")


###################################################
### code chunk number 23: fig4
###################################################
plot(peRes@pathways[["path:hsa05216"]], type = "boot")


###################################################
### code chunk number 24: pePathway_graph_Pert
###################################################
p <- peRes@pathways[["path:hsa05216"]]
g <- layoutGraph(p@map, layoutType = "dot")
graphRenderInfo(g) <- list(fixedsize = FALSE)
edgeRenderInfo(g) <- peEdgeRenderInfo(p)
nodeRenderInfo(g) <- peNodeRenderInfo(p)
renderGraph(g)


###################################################
### code chunk number 25: pePathway_graph_Pert2
###################################################
p <- peRes@pathways[["path:hsa04660"]]
g <- layoutGraph(p@map, layoutType = "dot")
graphRenderInfo(g) <- list(fixedsize = FALSE)
edgeRenderInfo(g) <- peEdgeRenderInfo(p)
nodeRenderInfo(g) <- peNodeRenderInfo(p)
renderGraph(g)


###################################################
### code chunk number 26: fig5
###################################################
p <- peRes@pathways[["path:hsa05216"]]
g <- layoutGraph(p@map, layoutType = "dot")
graphRenderInfo(g) <- list(fixedsize = FALSE)
edgeRenderInfo(g) <- peEdgeRenderInfo(p)
nodeRenderInfo(g) <- peNodeRenderInfo(p)
renderGraph(g)


###################################################
### code chunk number 27: fig6
###################################################
p <- peRes@pathways[["path:hsa04660"]]
g <- layoutGraph(p@map, layoutType = "dot")
graphRenderInfo(g) <- list(fixedsize = FALSE)
edgeRenderInfo(g) <- peEdgeRenderInfo(p)
nodeRenderInfo(g) <- peNodeRenderInfo(p)
renderGraph(g)


###################################################
### code chunk number 28: rontotools.Rnw:300-303
###################################################
require(graph)
require(ROntoTools)
kpg <- keggPathwayGraphs("hsa", verbose = FALSE)


###################################################
### code chunk number 29: rontotools.Rnw:307-308 (eval = FALSE)
###################################################
## kpg <- keggPathwayGraphs("hsa", updateCache = TRUE, verbose = TRUE)


###################################################
### code chunk number 30: rontotools.Rnw:313-315
###################################################
kpn <- keggPathwayNames("hsa")
head(kpn)
length(kpn)

###################################################
### code chunk number 31: rontotools.Rnw:321-323
###################################################
load(system.file("extdata/E-GEOD-21942.topTable.RData", package = "ROntoTools"))
head(top)


###################################################
### code chunk number 32: rontotools.Rnw:339-348
###################################################
fc <- top$logFC[top$adj.P.Val <= .01]
names(fc) <- top$entrez[top$adj.P.Val <= .01]

pv <- top$P.Value[top$adj.P.Val <= .01]
names(pv) <- top$entrez[top$adj.P.Val <= .01]

head(fc)

head(pv)


###################################################
### code chunk number 33: rontotools.Rnw:352-357
###################################################
fcAll <- top$logFC
names(fcAll) <- top$entrez

pvAll <- top$P.Value
names(pvAll) <- top$entrez


###################################################
### code chunk number 34: rontotools.Rnw:361-363
###################################################
ref <- top$entrez
head(ref)


###################################################
### code chunk number 35: rontotools.Rnw:375-376
###################################################
pDisRes <- pDis(x = fc, graphs = kpg, ref = ref,  nboot = 200, verbose = FALSE)


###################################################
### code chunk number 36: rontotools.Rnw:380-384
###################################################
head(Summary(pDisRes))

head(Summary(pDisRes, pathNames = kpn, totalpDis = FALSE,
             pORA = FALSE, comb.pv = NULL, order.by = "ppDis"))

