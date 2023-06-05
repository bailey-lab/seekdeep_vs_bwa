import subprocess
bam_file=snakemake.input['bam_file']
coords_file=snakemake.input['coords_file']
primer=snakemake.params['primer']
regional_sam=snakemake.output['regional_sam']

coord_dict={}
for line in open(coords_file):
	primer_name, chrom, start, end, size=line.strip().replace('"', '').split(',')
	if len(chrom)==1:
		chrom='0'+chrom
	coord_dict[primer_name]=f'Pf3D7_{chrom}_v3:{start}-{end}'

region=coord_dict[primer]
subprocess.call(f'samtools view {bam_file} {region} >{regional_sam}', shell=True)

#print(f'samtools view {bam_file} {region} >{regional_sam}')
