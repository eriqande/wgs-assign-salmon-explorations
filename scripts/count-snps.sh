# pass this a gzipped beagle file.  It assumes
# that the first three columns are: chrom_pos,
# ref alt.  So individuals start at position 4.
gunzip -c $1 | awk -v out1=$2 -v out2=$3 '

function abs(x) {return ((x < 0.0) ? -x : x)}

NR > 1 {
  nindiv = (NF - 3) / 3
  nHaveData=0;
  indiv=0
  for(i=4;i<=NF;i+=3) {
    indiv++
    a=$(i);
    b=$(i+1);
    c=$(i+2);
    # if they are not all 0.333 then there were reads at the site
    if(!(abs(a - 0.3333) < 0.01 && abs(b - 0.3333) < 0.01 &&  abs(c - 0.3333) < 0.01)) {
      nHaveData++
      indiv_snp_counts[indiv]++
    }
  }
  snps[nHaveData]++
}

END {
  # now, we need to print these things out
  # first do the distribution of the number of individuals for each SNP
  for(i=0;i<=nindiv;i++) printf("%d\t%d\n", i, snps[i]) > out1;

  # and now print out the number of SNPs for each of the individuals
  for(i=1;i<=nindiv;i++) printf("%d\t%d\n", i, indiv_snp_counts[i]) > out2;
}
'


