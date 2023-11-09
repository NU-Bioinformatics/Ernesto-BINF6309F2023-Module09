BiocManager::install("bioformat")
BiocManager::install("metagenomeSeq")
library(metagenomeSeq)
library(biomformat)
dataDirectory <- system.file("extdata", package = "metagenomeSeq")
lung = loadMeta(file.path(dataDirectory, "CHK_NAME.otus.count.csv"))
dim(lung$counts)
taxa = read.delim(file.path(dataDirectory, "CHK_otus.taxonomy.csv"),
                  stringsAsFactors = FALSE)
clin = loadPhenoData(file.path(dataDirectory, "CHK_clinical.csv"),
                     tran = TRUE)
ord = match(colnames(lung$counts), rownames(clin))
clin = clin[ord, ]
head(clin[1:6, ])
phenotypeData = AnnotatedDataFrame(clin)
phenotypeData
OTUdata = AnnotatedDataFrame(taxa)
OTUdata
obj = newMRexperiment(lung$counts,phenoData=phenotypeData,featureData=OTUdata)
obj
data(lungData)
lungData
data(mouseData)
mouseData
condition = mouseData$diet
mouseData = wrenchNorm(mouseData, condition = condition)
data(lungData)
p = cumNormStatFast(lungData)
lungData = cumNorm(lungData, p = p)
mat = MRcounts(lungData, norm = TRUE, log = TRUE)[1:5, 1:5]
exportMat(mat, file = file.path(dataDirectory, "tmp.tsv"))
exportStats(lungData[, 1:5], file = file.path(dataDirectory,
                                              "tmp.tsv"))
head(read.csv(file = file.path(dataDirectory, "tmp.tsv"), sep = "\t"))
data(lungData)
lungData = lungData[, -which(is.na(pData(lungData)$SmokingStatus))]
lungData = filterData(lungData, present = 30, depth = 1)
lungData <- cumNorm(lungData, p = 0.5)
pd <- pData(lungData)
mod <- model.matrix(1 + SmokingStatus, data = pd)
lungres1 = fitFeatureModel(lungData, mod)
head(MRcoefs(lungres1))

data(lungData)
controls = grep("Extraction.Control", pData(lungData)$SampleType)
lungTrim = lungData[, -controls]
rareFeatures = which(rowSums(MRcounts(lungTrim) > 0) < 10)
lungTrim = lungTrim[-rareFeatures, ]
lungp = cumNormStat(lungTrim, pFlag = TRUE, main = "Trimmed lung data")

lungTrim = cumNorm(lungTrim, p = lungp)
smokingStatus = pData(lungTrim)$SmokingStatus
bodySite = pData(lungTrim)$SampleType
normFactor = normFactors(lungTrim)
normFactor = log2(normFactor/median(normFactor) + 1)
mod = model.matrix(smokingStatus + bodySite + normFactor)
settings = zigControl(maxit = 10, verbose = TRUE)
fit = fitZig(obj = lungTrim, mod = mod, useCSSoffset = FALSE,
             control = settings)

settings = zigControl(maxit = 1, verbose = FALSE)
mod = model.matrix(bodySite)
colnames(mod) = levels(bodySite)
# fitting the ZIG model
res = fitZig(obj = lungTrim, mod = mod, control = settings)
# The output of fitZig contains a list of various useful
# items. hint: names(res). Probably the most useful is the
# limma 'MLArrayLM' object called fit.
zigFit = slot(res, "fit")
finalMod = slot(res, "fit")$design
contrast.matrix = makeContrasts(BAL.A - BAL.B, OW - PSB, levels = finalMod)
fit2 = contrasts.fit(zigFit, contrast.matrix)
fit2 = eBayes(fit2)
topTable(fit2)
classes = pData(mouseData)$diet
res = fitDO(mouseData[1:100, ], cl = classes, norm = FALSE, log = FALSE)
head(res)
cors = correlationTest(mouseData[55:60, ], norm = FALSE, log = FALSE)
head(cors)
obj = aggTax(mouseData, lvl = "phylum", out = "matrix")
head(obj[1:5, 1:5])
cors = correlationTest(mouseData[55:60, ], norm = FALSE, log = FALSE)
head(cors)

trials = pData(mouseData)$diet
heatmapColColors = brewer.pal(12, "Set3")[as.integer(factor(trials))]
heatmapCols = colorRampPalette(brewer.pal(9, "RdBu"))(50)
plotMRheatmap(obj = mouseData, n = 200, cexRow = 0.4, cexCol = 0.4,
              trace = "none", col = heatmapCols, ColSideColors = heatmapColColors)
