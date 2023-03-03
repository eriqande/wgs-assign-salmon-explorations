
SAMPS=[
"DPCh_plate1_A03_S3",
"DPCh_plate1_A04_S4",
"DPCh_plate1_B03_S15",
"DPCh_plate1_C01_S25",
"DPCh_plate1_D03_S39",
"DPCh_plate1_D04_S40",
"DPCh_plate1_E01_S49",
"DPCh_plate1_F01_S61",
"DPCh_plate1_G03_S75",
"DPCh_plate1_G04_S76",
"DPCh_plate1_G08_S80",
"DPCh_plate1_H03_S87",
"DPCh_plate1_H04_S88",
"DPCh_plate1_H06_S90",
"DPCh_plate2_A01_S97",
"DPCh_plate2_A02_S98",
"DPCh_plate2_A04_S100",
"DPCh_plate2_B01_S105",
"DPCh_plate2_B02_S106",
"DPCh_plate2_C01_S113",
"DPCh_plate2_C02_S114",
"DPCh_plate2_D03_S123",
"DPCh_plate2_E03_S131",
"DPCh_plate2_H03_S155"
]



rule all:
	input:
		expand("results/BAMs/{cov}X/rep_{rep}/{s}.bam", cov=[1.0, 0.1, 0.05, 0.01, 0.005, 0.001], rep = [1,2,3,4,5], s=SAMPS)


rule thin_bam:
	input:
		bam="BAMs/full-depth/{samp}.rmdup.bam",
		bai="BAMs/full-depth/{samp}.rmdup.bam.bai",
		dps="BAMs/coverages.tsv"
	output:
		bam="results/BAMs/{cov}X/rep_{rep}/{samp}.bam",
		bai="results/BAMs/{cov}X/rep_{rep}/{samp}.bam.bai",
	conda:
		"envs/samtools.yaml"
	shell:
		" OPT=$(awk '/{wildcards.samp}/ {{ fract = {wildcards.cov} / $2; if(fract < 1) print fract; else print \"NOSAMPLE\"; }}' {input.dps});  "
		" if [ $OPT = \"NOSAMPLE\" ]; then "
		"     ln  {input.bam} {output.bam}; "
		"     ln  {input.bai} {output.bai}; " 
		" else "
		"     samtools view --subsample $OPT --subsample-seed {wildcards.rep}  -b {input.bam} > {output.bam}; "
		"     samtools index {output.bam}; "
		" fi "


# in the following, "mprun" picks out the different vcfs like filt_snps05 and filt_snps05_miss30
rule get_sites:
	input:
		vcf="mega-post-bcf-exploratory-snakeflows/results/bcf_cal_chinook/{mprun}/all/thin_0_0/main.bcf"
	output:
		sites="results/sites/{mprun}.txt"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' {input.vcf} | sed 's/_/-/g;' > {output.sites}"

rule make_bamlist:
	input: 
		bam=expand("results/BAMs/{{cov}}X/rep_{{rep}}/{s}.bam", s=SAMPS)
	output:
		bamlist="results/bamlists/{cov}X/rep_{rep}/bamlist.txt"
	shell:
		"for i in {input.bam}; do echo $i; done > {output.bamlist} "

# this rule creates genotype likelihoods using ANGSD from the sites that we want
rule angd_likes:
	input:
		sites="results/sites/{mprun}.txt"
		
		bamlist="results/BAMs/{{cov}}X/rep_{rep}"
	output:
		sites="results/angsd_beagle/{cov}X/rep_{rep}/{samp}.bam"
	conda:
		"envs/bcftools.yaml"
	shell:
		"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' {input.vcf} > {output.sites}"

#		"envs/angsd.yaml"