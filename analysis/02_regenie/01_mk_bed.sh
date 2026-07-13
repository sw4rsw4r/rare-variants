DATA='../../data/'

awk '{print $1,$2,$3,$4,$5,$6}' $DATA/LCT_region_geno.raw > $DATA/LCT_region_geno.fam
