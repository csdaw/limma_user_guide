library(limma)
load("~/Downloads/Apoa1.RData")
RG2 <- backgroundCorrect(RG, method = "normexp")
MA <- normalizeWithinArrays(RG2)
MA2 <- normalizeBetweenArrays(MA, method = "Aq")
targets <- data.frame(
  Cy3 = I(rep("Pool", 16)),
  Cy5 = I(rep(c("WT", "KO"), each = 8))
)
targets
targets.sc <- targetsA2C(targets)
targets.sc
targets.sc$Target <- factor(targets.sc$Target, levels = c("Pool", "WT", "KO"))
targets.sc
design <- model.matrix(~Target, data = targets.sc)
#debugonce(intraspotCorrelation)
corfit <- intraspotCorrelation(MA, design)
fit <- lmscFit(MA, design, correlation = corfit$consensus)
cont.matrix <- cbind(KOvsWT = c(0, -1, 1))
fit2 <-contrasts.fit(fit, cont.matrix)
fit3 <- eBayes(fit2)
topTable(fit3, adjust.method = "fdr")
