'''
maps fastq files to a reference genome, converts them to bam format, sorts them,
and indexes the sorted bam file.
'''
configfile: 'modular_config.yaml'
output_folder=config['output_folder']+'/map_fastqs_paired_outputs'

rule all:
	input:
		sorted_bam=expand(output_folder+'/bam_main_files/{sample}_sorted.bam', sample=config['samples']),
		sorted_bam_index=expand(output_folder+'/bam_main_files/{sample}_sorted.bam.bai', sample=config['samples'])

rule align_paired_sam:
	'''
	the second fastq file is passed in as a parameter because if you have single
	end sequencing data, then the second fastq file is not a file that 'needs'
	to exist in order for the program to run. 2nd file empty quotes becomes
	extra whitespace.
	'''
	input:
		indexed_genome=config['indexed_genome_path'],
		fastq_file1=config['fastq_folder']+'/{sample}'+config['mate1_suffix'],
		fastq_file2=config['fastq_folder']+'/{sample}'+config['mate2_suffix']
	output:
		sample_sam=output_folder+'/sam_main_files/{sample}.sam'
	shell:
		'''
		bwa mem -o {output.sample_sam} {input.indexed_genome} {input.fastq_file1} {input.fastq_file2}
		'''

rule make_bam:
	input:
		sample_sam=output_folder+'/sam_main_files/{sample}.sam'
	output:
		sample_bam=temp(output_folder+'/bam_temp_files/{sample}.bam')
	shell:
		'samtools view -b -o {output.sample_bam} {input.sample_sam}'

rule sort_bam:
	input:
		sample_bam=output_folder+'/bam_temp_files/{sample}.bam'
	output:
		sorted_bam=output_folder+'/bam_main_files/{sample}_sorted.bam'
	shell:
		'samtools sort -o {output.sorted_bam} {input.sample_bam}'

rule index_bam:
	input:
		sorted_bam=output_folder+'/bam_main_files/{sample}_sorted.bam'
	output:
		bam_index=output_folder+'/bam_main_files/{sample}_sorted.bam.bai'
	shell:
		'samtools index {input.sorted_bam}'
