configfile: "modular_config.yaml"

module map_fastqs_paired:
	#prefix: 'map_fastqs_paired_inputs'
	snakefile: "snakemake_pipelines/map_fastqs_paired/map_fastqs_paired.smk"
	#config: "map_fastqs_paired_inputs/map_fastqs_paired.yaml"
	config: config
use rule * from map_fastqs_paired as map_fastqs_paired_*

module count_bwa_reads:
	snakefile: "snakemake_pipelines/count_bwa_reads/count_bwa_reads.smk"
	config: config
use rule * from count_bwa_reads as count_bwa_reads_*

module analyze_extraction_profile:
	snakefile: "snakemake_pipelines/seekdeep_vs_bwa/seekdeep_vs_bwa.smk"
	config: config
use rule * from analyze_extraction_profile as analyze_extraction_profile_*

rule all:
	input:
		rules.map_fastqs_paired_all.input,
		rules.count_bwa_reads_all.input,
		rules.analyze_extraction_profile_all.input,
		snakefile=config['output_folder']+'/snakemake_parameters/modular_seekdeep_vs_bwa.smk'
	default_target: True

rule copy_files: #needs to be edited and integrated
	'''
	copies snakemake and yaml files to the output folder so users and
	collaborators can see how the data was produced.
	'''
	input:
		snakefile='modular_seekdeep_vs_bwa.smk',
		configfile='modular_config.yaml',
		scripts='snakemake_pipelines'
	output:
		snakefile=config['output_folder']+'/snakemake_parameters/modular_seekdeep_vs_bwa.smk',
		configfile=config['output_folder']+'/snakemake_parameters/modular_config.yaml',
		scripts=directory(config['output_folder']+'/snakemake_pipelines')
	shell:
		'''
		cp {input.snakefile} {output.snakefile}
		cp {input.configfile} {output.configfile}
		cp -r {input.scripts} {output.scripts}
		'''
