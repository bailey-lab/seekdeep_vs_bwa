all_regions=snakemake.input['all_regions']
summary=snakemake.output['summary']

summary_dict={}
for region in all_regions:
	primer, sample=region.split('/')[-1].split('_')
	reads=set([])
	for line in open(region):
		reads.add(line.split()[0])
	summary_dict.setdefault(sample, {})
	summary_dict[sample].setdefault(primer, len(reads))

all_samples=list(summary_dict.keys())
all_primers=list(summary_dict[all_samples[0]].keys())

summary=open(summary, 'w')
summary.write('sample'+'\t'+'\t'.join(all_primers)+'\n')

for sample in all_samples:
	summary_line=[sample]
	for primer in all_primers:
		summary_line.append(str(summary_dict[sample][primer]))
	summary.write('\t'.join(summary_line)+'\n')
