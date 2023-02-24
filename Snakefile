
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
		expand("BAMs/{cov}X/{s}.bam", cov=[1, 0.1, .05, .01, .005, .001], s=SAMPS)
rule thin_bam:
	input:
		bam="BAMs/full-depth/{samp}.rmdup.bam",
		bai="BAMs/full-depth/{samp}.rmdup.bam.bai",
		dps="BAMs/coverages.tsv"
	output:
		bam="BAMs/{cov}X/{samp}.bam",
		#bai="BAMs/{cov}X/{samp}.bam.bai",
	conda:
		"envs/samtools.yaml"
	shell:
		" OPT=$(awk '/{wildcards.samp}/ {{ fract = {wildcards.cov} / $2; if(fract < 1) print \" --subsample \", fract; else print \"\"; }}' {input.dps});  "
		" samtools view $OPT {input.bam} > {output.bam}; "
		" samtools index {output.bam} "