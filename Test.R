fasta = Biostrings::readDNAStringSet("../DATASETS/JASPAR/MA0003.4.fasta.masked.dust")
fasta = as.character(fasta[!grepl("N", fasta)])
k = 14
n = length(fasta)
t = unname(nchar(fasta[1]))
m = unname(t - k + 1)

alpha = matrix (
  data = c(
  4335,	2763,	2842,	1088,	280,	 116,	  278,	  530,	 13567,	186,	  154,	   1221,	4220,	5204,
  4092,	4620,	2634,	3615,	15321, 15351,	1173,	  13536, 739,	  217,	  234,	   9445,	4783,	3915,
  4072,	3374,	3791,	9961,	187,	 212,	  750,	  1406,	 1229,	15432,	15422,	 4098,	3051,	3841,
  3469,	5211,	6701,	1304,	180,	 289,	  13767,	496,	 433,	  133,	  158,	   1204,	3914,	3008
  ), nrow = 4, ncol = k, byrow = T
)


alpha = alpha / colSums(alpha)
alpha

alpha = t(extraDistr::rdirichlet(k, c(1,1,1,1)))
beta = createMarkovChain(fasta, 1)
w = rep(1/m, m)
w[round(t/2)] = 0.1
new_alpha = oops(fasta = fasta, alpha = alpha, beta = beta, 1e-3, 1e3, w)
UTILS::plot_data(new_alpha)
