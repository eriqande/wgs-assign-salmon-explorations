# This is a Snakefile to run WGSassign multiple times
# on data from BAM files downsampled to various degrees.

#### INPUT VARIABLES ####


# Names of the samples for downsampling. For this to run
# there must be a bam named BAMs/full_depth/{sample_name}.rmdup.bam
# where {sample_name} is the name of each of the samples in SAMPS.
# Those BAM files should be indexed, too.
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

# A list of average read depths you would like to 
# downsample the SAMPS bams to.  FD does them at
# full depth.
COVIES=["FD", 1.0, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001]

# list of replicate numbers.  [1,2,3,4,5] will do 5 reps, numbered
# 1 through 5.
REPLIST=[1,2,3,4,5]


# Specifier for the initial condition of the reference BEAGLE file
# that you want to use. The beagle file must be in a directory
# that is named for this variable. 
MPRUN="filt_snps05_miss30"


#### INPUT FILES ####

# here are the files that have to be in place for this to work:

# The big beagle.gz file that holds the genotype likelihoods
# of the reference samples.  {mprun} let's you set different
# directories for different data sets
# --> reference/{mprun}/reference-beagle-gl.gz

# the TSV file that tells us what pop each sample in the
# reference sample file is from
# --> reference/{mprun}/ref_pops.tsv

# the BAM files named according to SAMPS + rmdup.bam and
# also their indexes

# A TSV file giving information about the mixture samples
# with columns: vcf_name NMFS_DNA_ID ref_pop group original_depth
# --> BAMs/sample_info.txt
# that looks like this:
# vcf_name            NMFS_DNA_ID ref_pop group   original_depth
# DPCh_plate1_G03_S75 T145116     Klamath mixture         1.33  
# DPCh_plate1_A03_S3  T145109     Klamath mixture         0.195 
# DPCh_plate1_H03_S87 T145117     Klamath mixture         0.149 
# ...




#### Snakefile Rules   #####

localrules: get_sites, index_sites, make_bamlist

rule all:
	input:
		"results/collated_mixture_likes.txt",
		"results/chinook-logls-and-depths-fig.pdf",
		"results/snps_in_indivs.tsv",
		"results/indivs_with_reads.tsv"
		#expand("results/BAMs/{cov}X/rep_{rep}/{s}.bam", cov=COVIES, rep = REPLIST, s=SAMPS),
		#expand("results/angsd_beagle/{mprun}/{cov}X/rep_{rep}/ref.beagle.gz", 
		#	mprun=["filt_snps05_miss30"], cov=COVIES, rep=REPLIST)


rule thin_bam:
	input:
		bam="BAMs/full-depth/{samp}.rmdup.bam",
		bai="BAMs/full-depth/{samp}.rmdup.bam.bai",
		dps="BAMs/sample_info.tsv"
	output:
		bam="results/BAMs/{cov}X/rep_{rep}/{samp}.bam",
		bai="results/BAMs/{cov}X/rep_{rep}/{samp}.bam.bai",
	conda:
		"envs/samtools.yaml"
	shell:
		" OPT=$(awk '/{wildcards.samp}/ {{ wc = \"{wildcards.cov}\"; if(wc == \"FD\") {{print \"NOSAMPLE\"; exit}} fract = wc / $NF; if(fract < 1) print fract; else print \"NOSAMPLE\"; }}' {input.dps});  "
		" if [ $OPT = \"NOSAMPLE\" ]; then "
		"     ln -sr {input.bam} {output.bam}; "
		"     ln -sr {input.bai} {output.bai}; " 
		" else "
		"     samtools view --subsample $OPT --subsample-seed {wildcards.rep}  -b {input.bam} > {output.bam}; "
		"     samtools index {output.bam}; "
		" fi "


# in the following, "mprun" picks out the different vcfs like filt_snps05 and filt_snps05_miss30
# rule get_sites:
# 	input:
# 		vcf="mega-post-bcf-exploratory-snakeflows/results/bcf_cal_chinook/{mprun}/all/thin_0_0/main.bcf"
# 	output:
# 		sites="results/sites/{mprun}.txt",
# 		chroms="results/sites/{mprun}.chroms"
# 	conda:
# 		"envs/bcftools.yaml"
# 	shell:
# 		"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' {input.vcf} > {output.sites}; "
# 		" cut -f1 {output.sites} | sort | uniq > {output.chroms}"

# this gets a list of the sites from the reference beagle file
# that ANGSD will do genotype likelihood calculations at from the
# BAMs
rule get_sites:
	input:
		beag="reference/{mprun}/reference-beagle-gl.gz"
	output:
		sites="results/sites/{mprun}.txt",
		chroms="results/sites/{mprun}.chroms"
	shell:
		"zcat {input.beag} | awk 'NR>1 {{print $1, $2, $3}}' | sed 's/_/ /g; s/-/_/g;' | "
        " awk '{{printf(\"%s\\t%s\\t%s\\t%s\\n\", $1, $2, $3, $4);}}' > {output.sites}; "
		" cut -f1 {output.sites} | sort | uniq > {output.chroms} "

rule index_sites:
	input:
		sites="results/sites/{mprun}.txt"
	output:
		idx="results/sites/{mprun}.txt.idx",
		bn="results/sites/{mprun}.txt.bin"
	conda:
		"envs/angsd.yaml"
	shell:
		" angsd sites index {input.sites} "


rule make_bamlist:
	input: 
		bam=expand("results/BAMs/{{cov}}X/rep_{{rep}}/{s}.bam", s=SAMPS)
	output:
		bamlist="results/bamlists/{cov}X/rep_{rep}/bamlist.txt"
	shell:
		"for i in {input.bam}; do echo $i; done > {output.bamlist} "

# this rule creates genotype likelihoods using ANGSD from the sites that we want
rule angsd_likes:
	input:
		sites="results/sites/{mprun}.txt",
		idx="results/sites/{mprun}.txt.idx",
		bn="results/sites/{mprun}.txt.bin",
		chroms="results/sites/{mprun}.chroms",
		bamlist="results/bamlists/{cov}X/rep_{rep}/bamlist.txt"
	output:
		beagle="results/angsd_beagle/{mprun}/{cov}X/rep_{rep}/out.beagle.gz"
	conda:
		"envs/angsd.yaml"
	resources:
		time = "3-00:00:00",
		mem_mb = "38400"
	threads: 8
	shell:
		" OUTPRE=$(echo {output.beagle} | sed 's/\.beagle\.gz//g;'); "
		" angsd -out $OUTPRE -GL 1 -rf {input.chroms}          "
		"    -nThreads {threads} -doGlf 2 -doMajorMinor 3 "
		"    -sites {input.sites} -bam {input.bamlist} "

# F&^%ck!  Angsd drops sites that don't have any data.
# So, now we need to make the "reference" beagle file concordant
# with that.  Part of that is also going to be changing the chromosome
# names to not have underscores
rule concordify_beagle_files:
	input:
		beagle="results/angsd_beagle/{mprun}/{cov}X/rep_{rep}/out.beagle.gz",
		big_ref="reference/{mprun}/reference-beagle-gl.gz"  # should update this to depend on mprun, but is OK for the miss30's
	output:
		ref_beagle="results/angsd_beagle/{mprun}/{cov}X/rep_{rep}/ref.beagle.gz",
		mix_beagle="results/angsd_beagle/{mprun}/{cov}X/rep_{rep}/mix.beagle.gz"
	shell:
		" " # first, deal with the underscores in the chromosome names   
		" zcat {input.beagle} | sed 's/^NC_/NC-/g; s/^NW_/NW-/g;' | gzip -c > {output.mix_beagle}; "
		" " # then pick out from the reference beagle sites only those found in the downsampled ones
		" (zcat {output.mix_beagle} | awk 'NR>1 {{print $1}}' ; zcat {input.big_ref}) | "
		" awk 'BEGIN {{OFS=\"\\t\"}} NF==1 {{g[$1]++; next}} /^marker/ || ($1 in g) {{print}}' | gzip -c > {output.ref_beagle} "
	


# this rule is here to address a referee comment.  They were
# curious about how many SNPs were involved in the different
# levels of downsampling.  I am going to go one further than that---
# I want to record both the total number of SNPs, as well as, for each
# individual, the number of SNPs that actually have data---i.e. the number
# that don't have flat likelihoods.  Maybe a better summary would be the number
# of SNPs that have only 1, 2, 3, individuals observed at them.  That would
# be an easier summary.  WE can grab that with an awk script.
rule count_snps_in_mixture_indivs:
	input: 
		mix_beagle="results/angsd_beagle/{mprun}/{cov}X/rep_{rep}/mix.beagle.gz"
	output:
		snps_in_indivs="results/angsd_beagle/{mprun}/{cov}X/rep_{rep}/snps_in_indivs.tsv",
		indivs_with_reads="results/angsd_beagle/{mprun}/{cov}X/rep_{rep}/indivs_with_reads.tsv",
	shell:
		"./scripts/count-snps.sh {input} {output.snps_in_indivs} {output.indivs_with_reads} " 


rule get_wgs_assign_installed:
	output:
		"results/wgs_assign_test.txt"
	conda:
		"envs/wgsassign.yaml"
	shell:
		"echo Yay > {output}"


rule test_WGSassign_in_snakemake:
    output:
    	"results/wgs-assign-spew.txt"
    conda:
    	"envs/wgsassign.yaml"
    shell:
    	"WGSassign > {output}"


rule get_reference_af:
	input:
		ref_beagle="results/angsd_beagle/{mprun}/{cov}X/rep_{rep}/ref.beagle.gz",
		IDs="reference/{mprun}/ref_pops.tsv"
	output:
		ref_af="results/angsd_beagle/{mprun}/{cov}X/rep_{rep}/reference.pop_af.npy"
	conda:
		"envs/wgsassign.yaml"
	threads: 4
	log:
		"results/logs/get_reference_af/{mprun}/{cov}X/rep_{rep}/log.txt"
	shell:
		" WGSassign --beagle {input.ref_beagle} --pop_af_IDs {input.IDs} "
		" --get_reference_af --out $(dirname {input.ref_beagle})/reference --threads {threads} > {log} 2>&1 "


rule infer_mixture_fish:
	input:
		ref_npy="results/angsd_beagle/{mprun}/{cov}X/rep_{rep}/reference.pop_af.npy",
		mix_beagle="results/angsd_beagle/{mprun}/{cov}X/rep_{rep}/mix.beagle.gz",
	output:
		mixfish="results/angsd_beagle/{mprun}/{cov}X/rep_{rep}/mixfish.pop_like.txt",
	conda:
		"envs/wgsassign.yaml"
	threads: 4
	log:
		"results/logs/infer_mixture_fish/{mprun}/{cov}X/rep_{rep}/log.txt"
	shell:
		" WGSassign --beagle {input.mix_beagle}  --pop_af_file {input.ref_npy} "
		"  --get_pop_like --out $(dirname {output.mixfish})/mixfish --threads {threads} > {log} 2>&1 "

rule collate_mixture_likes:
	input:
		expand("results/angsd_beagle/{mprun}/{cov}X/rep_{rep}/mixfish.pop_like.txt",
			mprun=["filt_snps05_miss30"], cov=COVIES, rep=REPLIST)
	output:
		"results/collated_mixture_likes.txt"
	shell:
		" for i in {input}; do awk -v file=$i '{{print file, ++n, $0}}' $i; done > {output} "


rule collate_snp_counts:
	input:
		snps_in_indivs=expand("results/angsd_beagle/{mprun}/{cov}X/rep_{rep}/snps_in_indivs.tsv",
			mprun=["filt_snps05_miss30"], cov=COVIES, rep=REPLIST),
		indivs_with_reads=expand("results/angsd_beagle/{mprun}/{cov}X/rep_{rep}/indivs_with_reads.tsv",
			mprun=["filt_snps05_miss30"], cov=COVIES, rep=REPLIST)
	output:
		snps_in_indivs="results/snps_in_indivs.tsv",
		indivs_with_reads="results/indivs_with_reads.tsv"
	shell:
		" for i in {input.snps_in_indivs}; do awk -v file=$i 'BEGIN {{OFS=\"\t\"; print \"file\tindiv_idx\tnum_snps\" }} {{print file, $0}}' $i; done > {output.snps_in_indivs}; "
		" for i in {input.indivs_with_reads}; do awk -v file=$i 'BEGIN {{OFS=\"\t\"; print \"file\tnum_indivs_with_data\tnum_snps\" }} {{print file, $0}}' $i; done > {output.indivs_with_reads}; "



# this is a stupid little rule to get the order of the BAMs in the
# output.
rule get_sample_order:
	params: 
		names = expand("{s}", s=SAMPS)
	output:
		"results/sample_order.txt"
	shell:
		"for i in {params.names}; do echo $i; done > {output} "


rule make_plot:
	input:
		sample_info="BAMs/sample_info.tsv",
		sample_order="results/sample_order.txt",
		sim_output="results/collated_mixture_likes.txt"
	output:
		outfig="results/chinook-logls-and-depths-fig.pdf"
	envmodules: "R/4.0.3"
	script:
		"scripts/make-plots.R"