'''
counts the number of reads that overlap each targeted region for each sample and
summarizes them in a table
'''
configfile: 'modular_config.yaml'

output_folder=config['output_folder']+'/obs_seekdeep_vs_exp_bwa'

rule all:
	input:
		edited_extraction_profile=output_folder+'/allExtractionProfile_edited.tab.txt',
		fraction_of_expected=output_folder+'/primer_by_sample_fraction_of_expected.tsv'

rule compare_extraction_profile_to_valid_counts:
	input:
		extraction_profile = config['extraction_profile'],
		valid_counts_file = config['output_folder']+'/count_bwa_reads_outputs/primer_by_sample/primer_by_sample_valid_counts.tsv'
	output:
		extraction_profile_edited = output_folder+'/allExtractionProfile_edited.tab.txt',
		valid_counts_file_edited = output_folder+'/primer_by_sample_fraction_of_expected.tsv'
	script:
		'scripts/extraction_summary_vs_mapped.py'
