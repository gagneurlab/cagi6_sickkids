file=$1

# restrict only to chromosomes 1-22,Y,X
regions=`echo "chr"{1..22} "chr"{Y,X} {1..22} Y X | tr ' ' ','`

bcftools view -O v --trim-alt-alleles -r $regions $file
