all_regions=snakemake.input['all_regions']
summary=snakemake.output['summary']
score_threshold=snakemake.params['score_threshold']
coords_file=snakemake.input['coords_file']
#valid_sam=open(snakemake.output['valid_sam'], 'w')

def make_coord_dict(coords_file):
	coord_dict={}
	for line_number, line in enumerate(open(coords_file)):
		if line_number>0:
			primer_name, chrom, start, end, size=line.strip().replace('"', '').split(',')
			if len(chrom)==1:
				chrom='0'+chrom
			chrom=f'Pf3D7_{chrom}_v3'
			coord_dict[primer_name]=[chrom, int(start), int(end)]
	return coord_dict

def calculate_score(CIGAR):
	'''
	takes an input CIGAR string and calculates a very rough alignment score for
	filtering purposes - insertions, deletions, and matches all count as 'good'
	'''
	digits=''
	score,size,clip=0,0,0
	for char in CIGAR:
		if char.isdigit():
			digits+=char
		elif char=='M' or char=='D':
			size+=int(digits)
			score+=int(digits)
			digits=''
		elif char=='I':
			score+=int(digits)
			digits=''
		elif char=='S' or char=='H':
			clip+=int(digits)
		else:
			digits=''
	return score, size, clip

coord_dict=make_coord_dict(coords_file)
summary_dict={}
all_reads=set([])
for region in all_regions:
	print('region is', region)
	primer, sample=region.split('/')[-1].split('_')
	target_chrom, target_start, target_end=coord_dict[primer]
	summary_dict.setdefault(sample, {})
	summary_dict[sample].setdefault(primer, set([]))
	read_dict={}
	region_file=open(region)
	for line in region_file:
		line=line.strip().split()
		read, flag, chrom, pos, CIGAR=line[0], int(line[1]), line[2], int(line[3]), line[5]
		score, size, clip=calculate_score(CIGAR)
		read_start, read_end=pos, pos+abs(size)
#		if read=='FS10001583:22:BSB09425-1416:1:1102:13980:1970':
#			print(line)
#			print('target was', target_chrom, target_start, target_end)
#			print('pos was', pos, 'size was', size, 'end pos was', pos+abs(size))
#		if chrom==target_chrom and pos>=target_start and (pos+abs(size))<=target_end:
		if score>score_threshold and flag<256 and chrom==target_chrom and (abs(read_start-target_start)<10 or abs(read_end-target_end)<10):
			read_dict[read]=read_dict.setdefault(read, 0)+1
#		if score>score_threshold and clip<10 and flag<256 and chrom==target_chrom and (abs(read_start-target_start)<10 or abs(read_end-target_end)<10):
#			pass
#		elif score>score_threshold:
#			print('\n', line)
#			print('target was', target_chrom, target_start, target_end)
#			print('pos was', pos, 'size was', size, 'end pos was', pos+abs(size))
	for read in read_dict:
		if read_dict[read]==2:
			summary_dict[sample][primer].add(read)
			if read in all_reads:
				print(read, flag, pos, CIGAR, size, score)
			all_reads.add(read)

#for line in open(region):
#	split_line=line.strip().split()
#	if split_line[0] in all_reads:
#		valid_sam.write(line)


for sample in summary_dict:
	print(sample)
	print(summary_dict[sample].keys())

all_samples=list(summary_dict.keys())
all_primers=list(summary_dict[all_samples[0]].keys())

summary=open(summary, 'w')
summary.write('sample'+'\t'+'\t'.join(all_primers)+'\n')

for sample in all_samples:
	summary_line=[sample]
	for primer in all_primers:
		summary_line.append(str(len(summary_dict[sample][primer])))
	summary.write('\t'.join(summary_line)+'\n')
