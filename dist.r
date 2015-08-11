library("ctc")

df <- read.table("distance.txt", sep="\t", col.names = c('a', 'b', 'sim'))
items <- with(df, unique(c(as.character(a), as.character(b))))
distm <- with(df, structure(df$sim,
                            Size = length(items),
                            Labels = items,
                            Diag = FALSE,
                            Upper = FALSE,
                            method = "user",
                            class = "dist"))
hc <- hclust(distm, method = "ward.D")
newick <- hc2Newick(hc)
write(newick, file = "dendo.newick")
