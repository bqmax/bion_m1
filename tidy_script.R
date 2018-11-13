### required packages
require(limma)
require(tidyverse)
require(AnnotationDbi)
require(org.Mm.eg.db)
require(Vennerable)
require(FactoMineR)
require(factoextra)
require(dendextend)
require(gage)
require(gageData)
require(Biobase)
require(biomaRt)
require(GOexpress)
require(GO.db)
require(superheat)
require(GOplot)

### colors for visualizations
col <- c('blue', 'green4', 'orange')
colvs <- c('#FEFFB3', '#8DD3C8', '#BEBBDA')

### load data
targets <- readTargets('targets.txt')
targets$Condition <- factor(targets$Condition, 
                            levels = unique(targets$Condition))
x <- read.maimages(targets, path = './scans', 
                   source = 'agilent', green.only = TRUE)

### intensities before preprocessing
dim(x) ## 44397
boxplot(log2(x$E), range=0, ylab="log2 intensity")
plotDensities(x, legend = F)

### preprocessing
y <- backgroundCorrect(x, method="normexp") ## background correction
plotDensities(y, legend = F)
y <- normalizeBetweenArrays(y, method="quantile") ## normalisation
plotDensities(y, legend = F)

### filtering and averaging
sum(y$genes$ControlType != 0) ## control probes
y <- y[y$genes$ControlType == 0, ]
dim(y) ## 43020

##### Ensembl annotation
y$genes$ensembl <- ifelse(grepl('ENSMUS', y$genes$GeneName), y$genes$GeneName, 
                          mapIds(org.Mm.eg.db, keys = y$genes$GeneName, 
                                 column = "ENSEMBL",
                                 keytype = "SYMBOL", 
                                 multiVals = 'first'))
sum(!is.na(y$genes$ensembl)) ## 32106

##### Entrez annotation
y$genes$entrez <- mapIds(org.Mm.eg.db, keys = y$genes$GeneName, 
                         column="ENTREZID",
                         keytype="SYMBOL", 
                         multiVals = 'first')
sum(!is.na(y$genes$entrez)) # 33368

##### Filtering items having ensembl ID
y <- y[!is.na(y$genes$ensembl), ]

##### Averaging by ensembl ID
y.ave <- avereps(y, ID=y$genes$ensembl)
y.ave$targets$Condition <- as.factor(y.ave$targets$Condition)
dim(y.ave) ## 20644

### Linear model fitting
##### Making design
f <- targets$Condition
design <- model.matrix(~ 0 + f)
colnames(design) <- levels(f)

##### Fitting
fit <- lmFit(y.ave, design)
summary(fit)
contrast.matrix <- makeContrasts(
    "Space-Control", "Recovery-Control", "Recovery-Space", 
    levels = design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

hist(fit$coefficients, main = 'Coefficients', xlab = '', col = 'wheat')
hist(p.adjust(fit$p.value, method = 'BH'), 
     main = 'FDR', xlab = '', col = 'wheat')
hist(fit$p.value, main = 'P-values', xlab = '', col = 'wheat')

output <- decideTests(fit, adjust.method = "BH", p.value = 0.05, lfc = 2)
summary(decideTests(fit, adjust.method = "BH", p.value = 0.05, lfc = 2))

vennDiagram(output, include=c("up", "down"), 
            circle.col= c('red', 'blue', 'green'))

### differential expression
##### space vs control
sc <- topTable(fit, adjust = "BH", coef = "Space-Control", 
               genelist = y.ave$genes, number = Inf)

s_sc <- subset(sc, abs(logFC) >= 2 & adj.P.Val < 0.05)

plot(sc$logFC, -log10(sc$adj.P.Val), 
     xlab = 'log2(FC)', ylab = '-log10(p-value)', main = 'S vs C')
points(s_sc$logFC, -log10(s_sc$adj.P.Val), col = 'red')

##### recovery vs control
rc <- topTable(fit, adjust = "BH", coef = "Recovery-Control", 
               genelist = y.ave$genes, number = Inf)
s_rc <- subset(rc, abs(logFC) >= 2 & adj.P.Val < 0.05)
s_rc$Regulation <- ifelse(s_rc$logFC > 0, 'UP', 'DOWN')
names(s_rc)[11:12] <- stringr::str_to_title(names(s_rc)[11:12])

plot(rc$logFC, -log10(rc$P.Value), 
     xlab = 'log2(FC)', ylab = '-log10(p-value)', main = 'Recovery vs Control')
points(s_rc$logFC, -log10(s_rc$P.Value), col = 'red')

##### recovery-space
rs <- topTable(fit, adjust = "BH", coef = "Recovery-Space", 
               genelist = y.ave$genes, number = Inf)
s_rs <- subset(rs, abs(logFC) >= 2 & adj.P.Val < 0.05)
s_rs$Regulation <- ifelse(s_rs$logFC > 0, 'UP', 'DOWN')
names(s_rs)[11:12] <- stringr::str_to_title(names(s_rs)[11:12])

plot(rs$logFC, -log10(rs$P.Value), 
     xlab = 'log2(FC)', ylab = '-log10(p-value)', main = 'Recovery vs Space')
points(s_rs$logFC, -log10(s_rs$P.Value), col = 'red')

discuss <- rs[rs$GeneName %in% c('Syp', 'Dlg4', 'Chat', 'Nefl', 'Nefm', 'Nefh'),
              c('GeneName', 'logFC')]
discuss$logFC <- 2^abs(discuss$logFC)

##### up & down
VennRaw <- Venn(list(`S vs C` = rownames(s_sc),
                     `R vs S` = rownames(s_rs),
                     `R vs C` = rownames(s_rc)))
plot(VennRaw, doWeights = F)


##### up
VennRaw <- Venn(list(`S vs C` = rownames(s_sc),
                     `R vs S` = rownames(s_rs[s_rs$Regulation == 'UP', ]),
                     `R vs C` = rownames(s_rc[s_rc$Regulation == 'UP', ])))
plot(VennRaw, doWeights = F)

##### down
VennRaw <- Venn(list(`S vs C` = rownames(s_sc),
                     `R vs S` = rownames(s_rs[s_rs$Regulation == 'DOWN', ]),
                     `R vs C` = rownames(s_rc[s_rc$Regulation == 'DOWN', ])))
plot(VennRaw, doWeights = F)

##### list of differentially expressed genes
deg <- union(s_rc$GeneName, s_rs$GeneName)

##### g:profiler overrepresentation analysis 
bp <- gProfileR::gprofiler(deg, organism = "mmusculus", 
                           src_filter = 'GO:BP', correction_method = 'fdr')
bp <- bp[order(bp$p.value),]
mf <- gProfileR::gprofiler(deg, organism = "mmusculus", 
                           src_filter = 'GO:MF', correction_method = 'fdr')
mf <- mf[order(mf$p.value),]
cc <- gProfileR::gprofiler(deg, organism = "mmusculus", 
                           src_filter = 'GO:CC', correction_method = 'fdr')
cc <- cc[order(cc$p.value),]
hp <- gProfileR::gprofiler(deg, organism = "mmusculus", 
                           src_filter = 'HP', correction_method = 'fdr')
hp <- hp[order(hp$p.value),]
hp$term.name <- tolower(hp$term.name)

##### mapping therms from previous steps to word clouds
wordcloud::wordcloud(hp$term.name[1:15], 1/-log10(hp$p.value[1:15]),
                     scale=c(0.2, 1), random.order=F, rot.per=0,
                     colors="black", fixed.asp = T)

### PCA
t <- y.ave$E
row.names(t) <- y.ave$genes$GeneName
res.pca <- PCA(t(t), ncp = Inf, graph = FALSE)

##### biplot
fviz_pca_biplot(res.pca, geom = 'point', axes = 1:2,
                habillage = y.ave$targets$Condition,
                palette = col, pointsize = 4,
                pointshape = 19, alpha = 0.6, mean.point = F, invisible = 'var',
                title = '') +
    theme(axis.title = element_text(size = 20, colour = 'black', face = 'bold'),
          axis.text = element_text(size = 18, colour = 'black'),
          legend.title = element_blank(),
          legend.text = element_text(size = 18, colour = 'black'),
          legend.position = 'bottom')
ggsave('t.svg')

##### scree plot
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 75),
         barcolor = 'black', main = '') +
    theme(axis.text = element_text(size = 14, color = "black"), 
          legend.text = element_text(size = 14), 
          axis.title = element_text(size = 16), 
          legend.title = element_blank())
ggsave('scree.svg', scale = 2)

##### Variables - PCs
fviz_pca_var(res.pca, select.var = list(contrib = 5), repel = T)

##### list of top 20 genes - contributors to PC1
gen <- rownames(res.pca$var$contrib[order(res.pca$var$contrib[, 1], 
                                          decreasing = T), ][1:20, ])

##### g:profiler overrepresentation analysis
go <- gProfileR::gprofiler(gen, organism = "mmusculus")

##### mapping therms from previous steps to word cloud
wordcloud::wordcloud(go$term.name, 1/-log10(go$p.value),
          scale=c(0.3, 1), random.order=F, rot.per=0,
          colors="black", fixed.asp = T)

##### subsets of these genes from each comparison
d1 <- subset(sc, GeneName %in% gen)
d1$Condition <- 1
d2 <- subset(rc, GeneName %in% gen)
d2$Condition <- 2
d3 <- subset(rs, GeneName %in% gen)
d3$Condition <- 3
d <- rbind(d1, d2, d3)

##### barplots for subsets
d$Condition <- factor(d$Condition, 
                      labels = c(' S vs C  ', ' R vs C  ', ' R vs S'))
d$FDR <- ifelse(d$adj.P.Val < 0.05, 1, 0)
d$FDR <- factor(d$FDR, labels = c('  p > 0.05  ', '  p < 0.05'))
ggplot(d, aes(fill = Condition, y = logFC, x = GeneName, colour = FDR))+
    geom_bar(stat = 'identity', position = position_dodge(), size = 0.7)+
    scale_color_manual(values = c('red2', 'green4'))+
    scale_fill_manual(values = colvs)+
    theme_bw()+
    scale_y_continuous(name = expression(paste(log[2],"(FC)",sep="")), 
                       breaks = seq(-13, 10, 2))+
    theme(axis.ticks.y = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_text(size = 16, colour = 'black'),
          axis.text = element_text(size = 14, colour = 'black'),
          legend.title = element_blank(),
          legend.text = element_text(size = 14, colour = 'black'),
          legend.position = 'bottom')+
    xlab('Gene name')+
    coord_flip()

##### heatmap for subsets

### hierarchical clustering
dend <- t(t) %>%  scale %>% 
    dist %>% hclust %>% as.dendrogram
par(mar = c(8,5,1,0))
dend %>% set("leaves_pch", 19) %>% 
    set("leaves_cex", 1.5) %>%
    set("labels_cex", 1.5) %>% 
    set("leaves_col", c('blue', 'green4', 'blue', 'green4', "orange", "orange")) %>%
    plot(cex.axis = 1.3, ylab = 'Distance', cex.lab = 1.5, font.lab = 2)

### KEGG enrichment
data("kegg.sets.mm")
data("sigmet.idx.mm")
kegg.sets.mm <- kegg.sets.mm[sigmet.idx.mm]

##### SPACE - CONTROL
sc <- topTable(fit, adjust="BH", coef="Space-Control", 
               genelist=y.ave$genes, number=Inf)
lfc <- sc$logFC
names(lfc) <- sc$entrez
keggres <- gage(lfc, gsets=kegg.sets.mm, same.dir=F)
sc_kegg <- as.data.frame(keggres$greater[, 1:5])
sc_kegg <- sc_kegg[sc_kegg$p.val < 0.05 & !(is.na(sc_kegg$p.val)), ]
# write.csv(sc_kegg, 'keggSC.csv')


##### RECOVERY - CONTROL
rc <- topTable(fit, adjust="BH", coef="Recovery-Control", 
               genelist=y.ave$genes, number=Inf)
lfc <- rc$logFC
names(lfc) <- rc$entrez
keggres <- gage(lfc, gsets=kegg.sets.mm, same.dir=F)
rc_kegg <- as.data.frame(keggres$greater[, 1:5])
rc_kegg <- rc_kegg[rc_kegg$q.val < 0.05 & !(is.na(rc_kegg$q.val)), ]
# write.csv(rc_kegg, 'keggRC.csv')

##### RECOVERY - SPACE
rs <- topTable(fit, adjust="BH", coef="Recovery-Space", 
               genelist=y.ave$genes, number=Inf)
lfc <- rs$logFC
names(lfc) <-rs$entrez
keggres <- gage(lfc, gsets=kegg.sets.mm, same.dir=F)
rs_kegg <- as.data.frame(keggres$greater[, 1:5])
rs_kegg <- rs_kegg[rs_kegg$q.val < 0.05 & !(is.na(rs_kegg$q.val)), ]
# write.csv(rs_kegg, 'keggRS.csv')
# nrow(rc_kegg)
# nrow(rs_kegg)
# setdiff(rownames(rs_kegg), rownames(rc_kegg))


##### Visualization
VennRaw <- Venn(list(`S vs C` = rownames(sc_kegg),
                     `R vs S` = rownames(rs_kegg),
                     `R vs C` = rownames(rc_kegg)))
plot(VennRaw, doWeights = F )

### GO enrichment
##### Preparing Expression Set
minimalSet <- y.ave$E
pData <- y.ave$targets
pData$Condition <- as.factor(pData$Condition)
phenoData <- new("AnnotatedDataFrame",
                 data=pData[,3, drop = F])
eset <- ExpressionSet(assayData=minimalSet, phenoData=phenoData)

##### Preparing local DBs
tt <- listAttributes(ensembl75)

ensembl75 <- useMart(
    host='feb2014.archive.ensembl.org',
    biomart='ENSEMBL_MART_ENSEMBL', dataset='mmusculus_gene_ensembl')

allgenes.Ensembl <- getBM(
    attributes=c('ensembl_gene_id', 'external_gene_id', 'description'),
    mart=ensembl75)
colnames(allgenes.Ensembl)[1] = 'gene_id'

allGO.Ensembl <- getBM(
    attributes=c('go_id', 'name_1006', 'namespace_1003'),
    mart=ensembl75)

GOgenes.Ensembl <- getBM(
    attributes=c('ensembl_gene_id', 'go_id'),
    mart=ensembl75)
colnames(GOgenes.Ensembl)[1] <- 'gene_id'

GOgenes.Ensembl <- GOgenes.Ensembl[GOgenes.Ensembl$go_id != '', ]
GOgenes.Ensembl <- GOgenes.Ensembl[GOgenes.Ensembl$gene_id != '', ]

##### Results of random forest
set.seed(4543)
resgo <- GO_analyse(
    eSet=eset, f="Condition",
    GO_genes=GOgenes.Ensembl,
    all_GO=allGO.Ensembl,
    all_genes=allgenes.Ensembl)

##### Results of permutations
pVal <- pValue_GO(result=resgo, N=1000) ## N=1000!!!

##### filtering
bp <- subset_scores(result = pVal,
                    namespace = "biological_process",
                    total = 10,
                    p.val=0.01)

##### mapping matrix
m <- table(bp$mapping$gene_id, bp$mapping$go_id)
attributes(m)$class <- 'matrix'
m <- m[rownames(bp$genes)[1:1000], ]
superheat(m, legend = F)

a <- hclust(dist(t(m), method = 'binary'))
b <- hclust(dist(m, method = 'binary'))

m <- m[b$order, a$order]
superheat(m, legend = F)

t <- m[867:894, 269:274]

superheat(t, legend = F,
          bottom.label.text.size = 8,
          bottom.label.text.angle = 90,
          left.label.text.size = 8,
          left.label.text.angle = 0,
          bottom.label.size = 0.2,
          left.label.size = 0.6,
          padding = 0.2)

superheat(m, legend = F,
          # yt = goscores * 1000,
          # yt.plot.type = 'scatter',
          # yt.point.size = 0.3,
          # yr = genscores * 1000,
          # yr.plot.type = 'scatter',
          # yr.point.size = 0.3,
          clustering.method = 'hierarchical',
          dist.method = 'binary',
          pretty.order.rows = T,
          pretty.order.cols = T,
          yt.axis = T,
          yr.axis = T,
          left.label.size = 0.05,
          bottom.label.size = 0.05,
          force.left.label = T,
          force.bottom.label = T
)


gind <- subset(bp$genes, rownames(bp$genes) %in% rownames(t))
gind <- gind[rownames(t), ]
t <- cbind(t, gind$Score)
rownames(t) <- mapIds(org.Mm.eg.db, keys = rownames(t), column = "SYMBOL",
                      keytype = "ENSEMBL", multiVals = 'first')
colnames(t)[ncol(t)] <- 'logFC'
colnames(t)
length(rownames(t))

##### visualization
GOChord(t, gene.size = 9, nlfc = 0)+
    theme(legend.text = element_text(size = 26),
          legend.title = element_text(size = 28, face = 'bold'))

# GO:0051387
heatmap_GO(go_id = 'GO:0007015', result = pVal, eSet = eset)

### DEG based
deg2 <- union(s_rc$Ensembl, s_rs$Ensembl)
deg3 <- intersect(deg2, rownames(bp$genes)[1:1000])
m <- table(bp$mapping$gene_id, bp$mapping$go_id)
attributes(m)$class <- 'matrix'

m <- m[rownames(m) %in% deg3, ]
superheat(m, legend = F)

a <- hclust(dist(t(m), method = 'binary'))
b <- hclust(dist(m, method = 'binary'))

m <- m[b$order, a$order]
superheat(m, legend = F)
superheat(m[1442:1479, 281:286], legend = F)

gind <- subset(bp$genes, rownames(bp$genes) %in% rownames(m[1442:1479, 281:286]))
gind <- gind[rownames(m[1442:1479, 281:286]), ]
t <- cbind(m[1442:1479, 281:286], gind$Score)
rownames(t) <- mapIds(org.Mm.eg.db, keys = rownames(t), column = "SYMBOL",
                      keytype = "ENSEMBL", multiVals = 'first')
colnames(t)[ncol(t)] <- 'logFC'
colnames(t)
length(rownames(t))

##### visualization
GOChord(t, nlfc = 1, gene.size = 9, lfc.min = 0, lfc.max = 0.005, 
        gene.order = "logFC")+
    theme(legend.text = element_text(size = 26),
          legend.title = element_text(size = 28, face = 'bold'))

