extraction_profile = snakemake.input['extraction_profile']
valid_counts_file = snakemake.input['valid_counts_file']
output_file_1 = snakemake.output['extraction_profile_edited']
output_file_2 = snakemake.output['valid_counts_file_edited']

def get_extracted_counts(extraction_profile):
	extracted_dict = {}
	for line in open(extraction_profile):
		if 'inputName' not in line:
			line = line.strip().split()
			sample = line[0].split('_')[0]
			primer = line[1].split('MID')[0]
			totalMatching = line[2]
			good = int(line[3].split('(')[0])
			total = int(line[2])
			if sample not in extracted_dict:
				extracted_dict[sample]={primer:[good]}
			if sample in extracted_dict:
				extracted_dict[sample][primer]=[good]
	return extracted_dict

def get_expected_counts(valid_counts_file):
	expected_list = []
	for line in open(valid_counts_file):
		if 'sample' in line:
			line = line.strip().split()
			primer_list = line[1:]
		if 'sample' not in line:
			line = line.strip().split()
			sample = line[0].split('.')[0]
			expected_counts = line[1:]
			for number in range(len(primer_list)):
				expected_list.append([sample,primer_list[number],int(expected_counts[number])])
	return expected_list,primer_list

def create_fraction_dict(expected_list, extracted_dict):
	fraction_dict = extracted_dict
	for expected in expected_list:
		if expected[1] in extracted_dict[expected[0]]:
			fraction_dict[expected[0]][expected[1]].append(expected[2])
		if expected[1] not in extracted_dict[expected[0]]:
			fraction_dict[expected[0]][expected[1]] = [0,expected[2]]
	return fraction_dict

def write_primer_by_sample_output(fraction_dict,output_file,primer_list):
	open(output_file,'w').write('\t'.join(['sample_name']+primer_list)+'\n')
	for sample in fraction_dict:
		output_line = [sample]
		for primer in primer_list:
			numerator,denominator = fraction_dict[sample][primer]
			if denominator>0:
				output_line.append(f'{numerator}/{denominator} ({round(100* numerator/denominator)} %)')
			else:
				output_line.append(f'{numerator}/{denominator} (0.00) %')

		open(output_file,'a').write('\t'.join(output_line)+'\n')

def create_edited_extraction_profile(output_file, extraction_profile, fraction_dict):
	for line in open(extraction_profile):
		line = line.strip()
		if 'inputName' in line:
			open(output_file,'w').write(f'{line}\tpercent of expected\n')
		if 'inputName' not in line:
				row = line.strip().split()
				sample = row[0].split('_')[0]
				primer = row[1].split('MID')[0]
				totalMatching = row[2]
				good = int(row[3].split('(')[0])
				total = int(row[2])
				numerator,denominator = fraction_dict[sample][primer]
				if denominator>0:
					output_item = f'{line}\t{numerator}/{denominator} ({round(100*numerator/denominator)} %)'
				else:
					output_item = f'{line}\t{numerator}/{denominator} (0.00) %)'
				open(output_file,'a').write(output_item+'\n')
extracted_dict = get_extracted_counts(extraction_profile)
expected_list,primer_list = get_expected_counts(valid_counts_file)[0:]
fraction_dict = create_fraction_dict(expected_list, extracted_dict)
primer_list.sort()
write_primer_by_sample_output(fraction_dict,output_file_2,primer_list)
create_edited_extraction_profile(output_file_1,extraction_profile,fraction_dict)

