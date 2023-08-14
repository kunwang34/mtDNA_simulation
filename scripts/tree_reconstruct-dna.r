library(ape)
library(phangorn)
library(getopt)

spec <- matrix(c('filename', 'f', 2, 'character', 'asfd'), byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
fn <- opt$filename


seq_dna = read.dna(paste('../results/dna_mut0.5_', fn, '.phy',sep=''), format='sequential')

dist_mat_dna = dist.dna(seq_dna, model='raw')

tree_dna_nj = nj(dist_mat_dna)

seq2_dna <- as.phyDat(seq_dna)

tree_dna_mp <- optim.parsimony(tree_dna_nj, seq2_dna)
tree_dna_mp <- nnls.phylo(tree_dna_mp, dist_mat_dna)

write.tree(tree_dna_nj, file=paste('../results/dna_mut0.5_', fn, '_nj.nwk',sep=''))
write.tree(tree_dna_mp, file=paste('../results/dna_mut0.5_', fn, '_mp.nwk',sep=''))
