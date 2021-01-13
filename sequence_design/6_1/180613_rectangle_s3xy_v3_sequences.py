import sys
import pickle
import copy
import inspect
import random
import os

def pv(name):
    record=inspect.getouterframes(inspect.currentframe())[1]
    frame=record[0]
    val=eval(name,frame.f_globals,frame.f_locals)
    print('{0}: {1}'.format(name, val))


def self_dG_FN(seq):
	energy = float(os.popen('hybrid-ss-min -q --tmin=45 --tmax=45 --sodium=0.003 --magnesium=0.015 --NA=DNA ' + seq).read())
	if energy > 0:
		energy = 0
	return energy

def duplex_dG_FN(seq):
	energy = float(os.popen('hybrid-ss-min -q --tmin=45 --tmax=45 --sodium=0.003 --magnesium=0.015 --NA=DNA ' + seq + 'LLL' + comp_seq_FN(seq)).read())
	if energy > 0:
		energy = 0
	return energy


def total_score_FN(word_ra, strand_pointer_ra):
	score = 0
	for sub_ra in strand_pointer_ra:
		strand_seq = strand_seq_FN(word_ra, sub_ra)
		score += self_dG_FN(strand_seq)
	return score


def comp_seq_FN(raw_sequence):
    uppercase = {'a':'A', 'A':'A', 'c':'C', 'C':'C', 'g':'G', 'G':'G', 't':'T', 'T':'T'}
    complement = {'a':'T', 'A':'T', 'c':'G', 'C':'G', 'g':'C', 'G':'C', 't':'A', 'T':'A'}
    antisense_seq = ''
    for letter in raw_sequence:
	if letter in uppercase:
            antisense_seq = complement[letter] + antisense_seq
    return antisense_seq

def stap_color_string_FN(stap_color_int):
	return color_dc[stap_color_int]


def break_strand_FN(vstrands, parity, bp):
	[helix_num, helix_base_num] = bp
	if bp != null_bp:
		[next_helix_num, next_helix_base_num] = vstrands[vsn[helix_num]][parity][helix_base_num][2:]
		if [next_helix_num, next_helix_base_num] != null_bp:
			vstrands[vsn[helix_num]][parity][helix_base_num][2:] = null_bp
			vstrands[vsn[next_helix_num]][parity][next_helix_base_num][:2] = null_bp
	else:
		print "No strand present at", parity, bp, "for breakage."
	return

def force_xover_FN(vstrands, parity, first_bp, second_bp):
	[first_helix_num, first_helix_base_num] = first_bp
	[second_helix_num, second_helix_base_num] = second_bp
	vstrands[vsn[first_helix_num]][parity][first_helix_base_num][2:] = second_bp
	vstrands[vsn[second_helix_num]][parity][second_helix_base_num][:2] = first_bp
	return


def strand_seq_FN(word_dc, pointer_ra):
	seq = ''
	if pointer_ra[0][1] == pointer_ra[1][1]:
		strand_type = 'y_strand'
	else:
		strand_type = 'x_strand'
	for [y, x] in pointer_ra:
		word = word_dc[y][x]
		if strand_type == 'y_strand':
			seq += comp_seq_FN(word)
		else:
			seq += word
	return seq


def report_FN(swap_attempt_num, path_dc, word_dc):
	self_dG_score = 0
	stack_dG_score = 0
	for group in ['g1-y', 'g2-x', 'g3-y']:
		for i in range(12):
			seq  = strand_seq_FN(word_dc, path_dc[group][i])
			self_dG_score += self_dG_FN(seq)
			if group == 'g2-x':
				stack_dG_score -= duplex_dG_FN(seq)
			if group == 'g3-y':
				stack_dG_score += duplex_dG_FN(seq)
	total_dG_score = self_dG_score + stack_dG_score
	print swap_attempt_num, 'self_dG:', self_dG_score, ', stack_dG:', stack_dG_score, ', total_dG_score:', total_dG_score
	return


#################################################################
#adjustable parameters
square_lattice = False #Set to False if honeycomb lattice is used
on_lattice_xover_scaf_loop_length = 2 #ssDNA loops at scaffold crossovers placed on the square lattice points?
off_lattice_xover_scaf_loop_length = 0 #ssDNA loops at scaffold crossovers placed off the square lattice points?
no_loop_exception_ra = [[17, 24], [16, 247]] #no scaffold loops allowed after the indicated base pointer positions
caDNAno_json_filename = '170204_axs_queen_v0.json'
maxiscaf_seq_filename = 'p8064.txt'
prefix = '180613_rectangle_s3xy_v2'	#Annotation prefix
#################################################################





#Read files
##Need to read in designed sequences for every length of scaffold strand in the caDNAno json file,
##otherwise an error will occur
input_file = file(caDNAno_json_filename, 'r')
all_file = eval(input_file.read())
input_file.close()
vstrands = all_file['vstrands']
num_vstrands = len(vstrands)
vsn = {}
for vstrand_num in range(num_vstrands):
	vsn[vstrands[vstrand_num]['num']] = vstrand_num
vsn[-1] = -1


null_bp = [-1, -1]


#Import maxiscaf seq
second_hairpin_seq = 'GGGTGATGGTTCACGTAGTGGGCCATCGCCC'
input_file = file(maxiscaf_seq_filename, 'r')
maxiscaf_seq = input_file.read()
input_file.close()
maxiscaf_seq = comp_seq_FN(comp_seq_FN(maxiscaf_seq))
end_hairpins_position = maxiscaf_seq.find(second_hairpin_seq) + len(second_hairpin_seq)
maxiscaf_seq = maxiscaf_seq[end_hairpins_position:] + maxiscaf_seq[:end_hairpins_position]	


pickled_designed_seq_dc_filename = '130126_1127_designed_seqs_05-69.txt'
input_file = file(pickled_designed_seq_dc_filename, 'r')
seq_dc = pickle.load(input_file)
input_file.close()
seq_dc = seq_dc.copy()



#Initialization of variables
null_bp = [-1, -1]
num_vstrands = len(vstrands)
num_helix_bases = len(vstrands[0]['scaf'])
stap_color_dc  = {}
stap_color_dc[13369344] = 'red'
stap_color_dc[16204552] = 'red orange'
stap_color_dc[16225054] = 'light orange'
stap_color_dc[11184640] = 'olive'
stap_color_dc[5749504]  = 'light green'
stap_color_dc[29184]    = 'dark green'
stap_color_dc[243362]   = 'cyan'
stap_color_dc[1507550]  = 'blue'
stap_color_dc[7536862]  = 'purple'
stap_color_dc[12060012] = 'magenta'
stap_color_dc[3355443]  = 'dark gray'
stap_color_dc[8947848]  = 'light gray'

#Generate arrays for stap paths and scaf paths
scaf_path_ra = []
stap_path_ra = []
for vstrand_num in range(num_vstrands):
	for helix_base_num in range(num_helix_bases):
		for [parity, path_ra] in [['scaf', scaf_path_ra], ['stap', stap_path_ra]]:
			[prev_helix_num, prev_helix_base_num, next_helix_num, next_helix_base_num] = vstrands[vstrand_num][parity][helix_base_num]
			if ([prev_helix_num, prev_helix_base_num] == null_bp) and ([next_helix_num, next_helix_base_num] != null_bp):
				sub_ra = []
				curr_helix_num, curr_helix_base_num = vstrands[vstrand_num]['num'], helix_base_num
				end_of_strand = False
				while not end_of_strand:
					sub_ra.append([curr_helix_num, curr_helix_base_num])
					[curr_helix_num, curr_helix_base_num] = vstrands[vsn[curr_helix_num]][parity][curr_helix_base_num][2:]
					end_of_strand = [curr_helix_num, curr_helix_base_num] == null_bp
				path_ra.append(sub_ra)

print "scaf_path_ra length", len(scaf_path_ra)
print "stap_path_ra length", len(stap_path_ra)


#Determine length of maxiscaf and miniscafs
seq_counter_dc = {}
maxiscaf_length = 0
for sub_ra in scaf_path_ra[:1]:
	for [helix_num, helix_base_num] in sub_ra:
		maxiscaf_length += 1 + vstrands[vsn[helix_num]]['skip'][helix_base_num] + vstrands[vsn[helix_num]]['loop'][helix_base_num]
		if vstrands[vsn[helix_num]]['scaf'][helix_base_num][2] not in [-1, helix_num]:
			maxiscaf_length += 2
print "maxiscaf_length", maxiscaf_length
seq_counter_dc[maxiscaf_length] = 0
seq_dc[maxiscaf_length] = [maxiscaf_seq]
seq_counter_dc[48] = 0
seq_counter_dc[56] = 0

#Assign scaf base sequences
sub_ra = ['.' for i in range(num_helix_bases)]
scaf_base_seq_dc = {}
for vstrand in vstrands:
	helix_num = vstrand['num']
	scaf_base_seq_dc[helix_num] = sub_ra[:]
for sub_ra in scaf_path_ra:
	seq_length = 0
	for [helix_num, helix_base_num] in sub_ra:
		seq_length += 1 + vstrands[vsn[helix_num]]['skip'][helix_base_num] + vstrands[vsn[helix_num]]['loop'][helix_base_num]
		if vstrands[vsn[helix_num]]['scaf'][helix_base_num][2] not in [-1, helix_num]:
			seq_length += 2
	seq = seq_dc[seq_length][seq_counter_dc[seq_length]]
	seq_counter_dc[seq_length] += 1
	seq_pointer = 0
	for [helix_num, helix_base_num] in sub_ra:
		base_length = 1 + vstrands[vsn[helix_num]]['skip'][helix_base_num] + vstrands[vsn[helix_num]]['loop'][helix_base_num]
		scaf_base_seq_dc[helix_num][helix_base_num] = seq[seq_pointer:seq_pointer + base_length]
		seq_pointer += base_length
		if vstrands[vsn[helix_num]]['scaf'][helix_base_num][2] not in [-1, helix_num]:
			seq_pointer += 2


#Print scaf_base_sequences
for vstrand in vstrands:
	helix_num = vstrand['num']
	print helix_num, ':', ''.join(scaf_base_seq_dc[helix_num])



meta_stap_color_dc = {}
for vstrand in vstrands:
	helix_num = vstrand['num']
	meta_stap_color_dc[helix_num] = [-1 for helix_base_num in range(num_helix_bases)]
	for [helix_base_num, stap_color_int] in vstrand['stap_colors']:
		meta_stap_color_dc[helix_num][helix_base_num] = stap_color_dc[stap_color_int]


#Generate staple strand output
stap_output_ra = []
for sub_ra in stap_path_ra:
	seq = ''
	for [helix_num, helix_base_num] in sub_ra:
		seq += comp_seq_FN(scaf_base_seq_dc[helix_num][helix_base_num])
	start_pointer = sub_ra[0]
	end_pointer = sub_ra[-1]
	[first_helix_num, first_helix_base_num] = start_pointer
	stap_color = meta_stap_color_dc[first_helix_num][first_helix_base_num]
	stap_output_ra.append([stap_color, len(seq), start_pointer, end_pointer, seq])


#Sort staple strands according to caDNAno color
sorted_stap_output_ra = sorted(stap_output_ra, key = lambda stap_output:stap_output[0])

print
print


#Print sorted staple strand sequences with annotations
num_additions = 0
for sub_ra in sorted_stap_output_ra:
	seq = sub_ra[-1]
	if (sub_ra[3][1] == 127) and (sub_ra[3][0] in [13, 15, 17, 19, 31, 33, 35, 37, 53, 55, 57, 59, 73, 75, 77]):
		seq += 'A'
		num_additions += 1
	note  = prefix + ' stap strand, ' + str(sub_ra[0]) + ', ' + str(len(seq)) + 'mer, '
	note += str(sub_ra[2]) + ' start, ' + str(sub_ra[3]) + ' end'
	print seq + '\t\t' + note
print
print

#Add terminal A to 56mer miniscaf strands
#Print miniscaf strands and corresponding stap strands
scaf_output_ra = []
for sub_ra in scaf_path_ra[1:]:
	seq = ''
	for [helix_num, helix_base_num] in sub_ra:
		seq += scaf_base_seq_dc[helix_num][helix_base_num]
		if len(seq) == 56:
			seq += 'A'
	start_pointer = sub_ra[0]
	end_pointer = sub_ra[-1]
	[first_helix_num, first_helix_base_num] = start_pointer
	scaf_output_ra.append([len(seq), start_pointer, end_pointer, seq])

for sub_ra in scaf_output_ra:
	seq = sub_ra[-1]
	note  = prefix + ' miniscaf strand, ' + str(sub_ra[0]) + 'mer, '
	note += str(sub_ra[1]) + ' start, ' + str(sub_ra[2]) + ' end'
	print seq + '\t\t' + note
	

print
print
print









#v1_66_g1f_x=0-7_d0


#read in word_dc
input_file = file('180613_1508_pickled_rectangle_s3xy_v3_word_dc.txt', 'r')
word_dc = pickle.load(input_file)
input_file.close()

#Initialize path_dc and map_dc
path_dc = {}

path_dc['g1-y'] = []
for x in range(12):
	sub_ra = []
	for y in range(0 + x, -12 + x, -1):
		if x%2 == 0:
			sub_ra.append([y, x])
		else:
			sub_ra.insert(0, [y, x])
	path_dc['g1-y'].append(sub_ra)

path_dc['g2-x'] = []
for y in range(1, 13):
	sub_ra = []
	for x in range(y, y + 12):
		if y%2 == 0:
			sub_ra.append([y, x])
		else:
			sub_ra.insert(0, [y, x])
	path_dc['g2-x'].append(sub_ra)

path_dc['g3-y'] = []
for x in range(12, 24):
	sub_ra = []
	for y in range(0 + x, -12 + x, -1):
		if x%2 == 0:
			sub_ra.append([y, x])
		else:
			sub_ra.insert(0, [y, x])
	path_dc['g3-y'].append(sub_ra)


#score by self-energy, heterogeneous-energy, x-versus-y-stacking energy
report_FN(0, path_dc, word_dc)

strand_dc = {}
for i in range(12):
	seq = strand_seq_FN(word_dc, path_dc['g1-y'][i])
	print seq + '\t\t' + '180613_rectangle_s3xy_v3, g1-y, x = ', i
for i in range(12):
	seq = strand_seq_FN(word_dc, path_dc['g2-x'][i])
	seq = 'T'*8 + seq + 'T'*8
	print seq + '\t\t' + '180613_rectangle_s3xy_v3, g2-x, y = ', i + 1
for i in range(12):
	seq = strand_seq_FN(word_dc, path_dc['g3-y'][i])
	print seq + '\t\t' + '180613_rectangle_s3xy_v3, g3-y, x = ', i + 12















