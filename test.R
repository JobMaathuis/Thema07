library(DESeq2)

annotation <- data.frame(Patient = factor(rep(1:5, times = 3), 
                                          labels = c("m29", "m42", "m43", "s10", "s11")),
                         Condition = factor(rep(c(rep("mild", times = 3),
                                                  rep("severe", times = 2)), times = 3),
                                            levels = c('severe', 'mild')),
                         Repeat = factor(rep(rep(1:3, times = 1), each = 5),
                                         labels = c("R1", "R2", "R3")))

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = annotation,
                              design = ~ Condition)
dds <- DESeq(dds, betaPrior = FALSE)
dds.collapse <- collapseReplicates(dds, groupby = annotation$Patient)
resultsNames(dds.collapse)

# res <- results(dds.collapse, alpha = 0.1, contrast = c("Condition", "mild", "severe"))
res <- lfcShrink(dds.collapse, type="apeglm", coef = "Condition_mild_vs_severe")

## Note: input data is the corrected DESeq2 output using the 'lfcShrink' function (see chapter 4)
deseq.volcano(res = res, datasetName = 'Corona')

