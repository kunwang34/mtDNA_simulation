library(ape)
library(phangorn)
library(getopt)

spec <- matrix(c('filename', 'f', 2, 'character', 'asfd',
                 'path', 'p', 2, 'character', 'asdffdsd'
                ), byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
fn <- opt$filename
path <- opt$path

files <- list.files(path)
files <- files[grep(paste(fn, '.phy$', sep=''), files)]

for (file in files){
    seq = read.dna(paste(path, '/', file,sep=''), format='sequential')
    dist_mat = dist.dna(seq, model='raw')
    tree_nj = nj(dist_mat)
    seq2 <- as.phyDat(seq)
    tree_mp <- optim.parsimony(tree_nj, seq2)
    write.tree(tree_nj, file=paste(path, '/', file, '_nj.nwk',sep=''))
    write.tree(tree_mp, file=paste(path, '/', file, '_mp.nwk',sep=''))
}
