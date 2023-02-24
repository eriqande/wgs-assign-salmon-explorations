



rule thin_bam:
	input:
		bam="BAMs/full-depth/{samp}.bam",
		bai="BAMs/full-depth/{samp}.bam.bai",
		cov="BAMs/coverages.tsv"
	output:
		bam="BAMs/{cov}X/{samp}.bam",
		#bai="BAMs/{cov}X/{samp}.bam.bai",
	conda:
		"envs/samtools.yaml"
	shell:
		"FC=$(awk '/{samp}/ {{print $2}}') {input.cov} > {output.bam} "