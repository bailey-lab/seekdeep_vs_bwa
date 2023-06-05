'''
counts the number of reads that overlap each targeted region for each sample and
summarizes them in a table
'''
configfile: 'modular_config.yaml'

output_folder=config['output_folder']+'/count_bwa_reads_outputs'

rule all:
	input:
		regional_sam=expand(output_folder+'/primer_overlapping_reads/{primer}_{sample}.sam', primer=config['primers'], sample=config['samples']),
		summarized=output_folder+'/primer_by_sample/primer_by_sample_raw_counts.tsv',
#		summary='count_bwa_reads_outputs/'+config['count_bwa_reads_output_folder']+'/sample_by_primer/primer_by_sample_valid_counts.tsv',
		summary_alt=output_folder+'/primer_by_sample/primer_by_sample_valid_counts.tsv'

rule extract_regions:
	input:
		bam_file=config['output_folder']+'/map_fastqs_paired_outputs/bam_main_files/{sample}_sorted.bam',
		bam_index = config['output_folder']+'/map_fastqs_paired_outputs/bam_main_files/{sample}_sorted.bam.bai',
		coords_file=config['amplicon_location_size']
	params:
		primer='{primer}'
	output:
		regional_sam=output_folder+'/primer_overlapping_reads/{primer}_{sample}.sam',
	script:
		'scripts/extract_regions.py'

rule summarize_all_reads:
	input:
		all_regions=expand(output_folder+'/primer_overlapping_reads/{primer}_{sample}.sam', primer=config['primers'], sample=config['samples'])
	output:
		summary=output_folder+'/primer_by_sample/primer_by_sample_raw_counts.tsv'
	script:
		'scripts/summarize_all_reads.py'

rule summarize_valid_reads:
	input:
		all_regions=expand(output_folder+'/primer_overlapping_reads/{primer}_{sample}.sam', primer=config['primers'], sample=config['samples']),
		coords_file=config['amplicon_location_size']
	params:
		score_threshold=config['score_threshold']
	output:
		summary=output_folder+'/primer_by_sample/primer_by_sample_valid_counts.tsv'
	script:
		'scripts/summarize_valid_reads.py'
