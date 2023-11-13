# pass this a gzipped beagle file.  It assumes
# that the first three columns are: chrom_pos,
# ref alt.  So individuals start at position 4.
gunzip -c $1 | awk -v out1=$2 '


NR > 1 {
  nindiv = (NF - 3) / 3
  indiv=0
  for(i=4;i<=NF;i+=3) {
    indiv++
    a=$(i);
    b=$(i+1);
    c=$(i+2);
    # if they are not all "0.333333" then we will say there were reads at the site
    if(!(a == "0.333333" && b == "0.333333" && c == "0.333333")) {
      indiv_snp_counts[indiv]++
    }
  }
}

END {
  # now, we need to print these things out
  # first do the distribution of the number of individuals for each SNP
  for(i=1;i<=nindiv;i++) printf("%d\t%d\n", i, indiv_snp_counts[i]) > out1;
}
'


