all_regions=snakemake.input['all_regions']
summary=snakemake.output['summary']
score_threshold=snakemake.params['score_threshold']

def calculate_score(CIGAR):
	'''
	takes an input CIGAR string and calculates a very rough alignment score for
	filtering purposes - insertions, deletions, and matches all count as 'good'
	'''
	digits=''
	score=0
	for char in CIGAR:
		if char.isdigit():
			digits+=char
		elif char in 'DIM':
			score+=int(digits)
			digits=''
		else:
			digits=''
	return score

summary_dict={}
for region in all_regions:
	primer, sample=region.split('/')[-1].split('_')
	read_dict={}
	for line in open(region):
		line=line.strip().split()
		read, flag, pos, CIGAR, size=line[0], line[1], line[3], line[5], line[8]
		score=calculate_score(CIGAR)
		if score>score_threshold:
			read_dict[read]=read_dict.setdefault(read, 0)+1
	good_read_count=0
	for read in read_dict:
		if read_dict[read]==2:
			good_read_count+=1
	summary_dict.setdefault(sample, {})
	summary_dict[sample].setdefault(primer, good_read_count)

all_samples=list(summary_dict.keys())
all_primers=list(summary_dict[all_samples[0]].keys())

summary=open(summary, 'w')
summary.write('sample'+'\t'+'\t'.join(all_primers)+'\n')

for sample in all_samples:
	summary_line=[sample]
	for primer in all_primers:
		summary_line.append(str(summary_dict[sample][primer]))
	summary.write('\t'.join(summary_line)+'\n')
