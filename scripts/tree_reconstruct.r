library(ape)
library(phangorn)
library(getopt)

spec <- matrix(c('filename', 'f', 2, 'character', 'asfd'), byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
fn <- opt$filename

seq_mt = read.dna(paste('../results/mt_constant_mut_', fn, '.phy',sep=''), format='sequential')
seq_dna = read.dna(paste('../results/dna_constant_mut_', fn, '.phy',sep=''), format='sequential')
dist_mat_mt = dist.dna(seq_mt, model='raw')
dist_mat_dna = dist.dna(seq_dna, model='raw')
tree_mt_nj = nj(dist_mat_mt)
tree_dna_nj = nj(dist_mat_dna)

seq2_mt <- as.phyDat(seq_mt)
seq2_dna <- as.phyDat(seq_dna)

tree_mt_mp <- optim.parsimony(tree_mt_nj, seq2_mt)
tree_mt_mp <- nnls.phylo(tree_mt_mp, dist_mat_mt)

tree_dna_mp <- optim.parsimony(tree_dna_nj, seq2_dna)
tree_dna_mp <- nnls.phylo(tree_dna_mp, dist_mat_dna)

write.tree(tree_mt_nj, file=paste('../results/mt_constant_mut_', fn, '_nj.nwk',sep=''))
write.tree(tree_mt_mp, file=paste('../results/mt_constant_mut_', fn, '_mp.nwk',sep=''))

write.tree(tree_dna_nj, file=paste('../results/dna_constant_mut_', fn, '_nj.nwk',sep=''))
write.tree(tree_dna_mp, file=paste('../results/dna_constant_mut_', fn, '_mp.nwk',sep=''))