library("ctc")

## input should be Taxa1\tTaxa2\tDistance

f <- file("stdin")
open(f, blocking = T)
df <- read.table(f, sep="\t", col.names = c('a', 'b', 'sim'))
close(f)

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

write(newick, stdout())
