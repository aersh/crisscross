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
	energy = float(os.popen('hybrid-ss-min -q --tmin=20 --tmax=20 --sodium=0.003 --magnesium=0.030 --NA=DNA ' + seq).read())
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


def strand_seq_FN(word_ra, pointer_ra):
	seq = ''
	if pointer_ra[0][1] == pointer_ra[1][1]:
		strand_type = 'y_strand'
	else:
		strand_type = 'x_strand'
	for [z, y, x] in pointer_ra:
		word = word_ra[z][y][x]
		if strand_type == 'y_strand':
			seq += comp_seq_FN(word)
		else:
			seq += word
	word_length = len(word_ra[0][0][0])
	T_brush = 'T'*word_length
	seq = T_brush + seq[word_length:-word_length] + T_brush
	return seq



def y_strand_pointer_FN(block_ra):
	strand_pointer_ra = []
	for x in range(1, 9):
		sub_ra = []
		for z in block_ra:
			for y in range(10):
				if x%2 == 1:
					sub_ra.append([z, y, x])
				else:
					sub_ra.insert(0, [z, y, x])
		strand_pointer_ra.append(sub_ra)
	return strand_pointer_ra


def x_strand_pointer_FN(block_ra):
	strand_pointer_ra = []
	for y in range(1, 9):
		sub_ra = []
		for z in block_ra:
			for x in range(10):
				if y%2 == 1:
					sub_ra.append([z, y, x])
				else:
					sub_ra.insert(0, [z, y, x])
		strand_pointer_ra.append(sub_ra)
	return strand_pointer_ra






#################################################################
#adjustable parameters
square_lattice = False #Set to False if honeycomb lattice is used
on_lattice_xover_scaf_loop_length = 2 #ssDNA loops at scaffold crossovers placed on the square lattice points?
off_lattice_xover_scaf_loop_length = 0 #ssDNA loops at scaffold crossovers placed off the square lattice points?
no_loop_exception_ra = [[17, 24], [16, 247]] #no scaffold loops allowed after the indicated base pointer positions
caDNAno_json_filename = '160915_pxs_v5.1_queen.json'
maxiscaf_seq_filename = 'p8064.txt'
prefix = '160915_pxs_v5.1_queen'	#Annotation prefix
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


#Create queen_base_seq_dc
queen_base_seq_dc = {}
for helix_num in range(8):
	queen_base_seq_dc[helix_num] = scaf_base_seq_dc[helix_num][224:352]


#Generate strands
for word_length in [5, 6]:
	word_cache_ra = copy.deepcopy(seq_dc[word_length])
	queen_offset = word_length
	print '\n\n\nword length', word_length, '\tqueen offset', queen_offset
	#Initialize word lattice
	sub_word_ra = [['.'*word_length]*10 for y in range(10)]
	word_ra = [copy.deepcopy(sub_word_ra) for z in range(5)]

	#Fill in z=0 (queen) words
	for y in range(1, 9):
		for x in range(1, 9):
			word = ''
			for base_num in range(word_length):
				base = queen_base_seq_dc[y - 1][queen_offset + (x - 1)*word_length + base_num]
				if y%2 == 1:
					word += base
				else:
					word = base + word
			word_ra[0][y][x] = word
	counter = 0
	for x in range(1, 9):
		word_ra[0][0][x] = word_cache_ra.pop(0)
		word_ra[0][9][x] = word_cache_ra.pop(0)
	
	#Fill other layer words
	for z in range(1, 5):
		for y in range(10):
			for x in range(10):
				if (y in range(1, 9)) or (x in range(1, 9)):
					word_ra[z][y][x] = word_cache_ra.pop(0)
	print len(word_cache_ra)
#	for z in range(5):
#		print '\nlayer ', z
#		for y in range(10):
#			print word_ra[z][y]

	#Define all strands in terms of the words
	strand_pointer_ra = []
	strand_pointer_ra += y_strand_pointer_FN([2, 1, 0]) #g1
	strand_pointer_ra += x_strand_pointer_FN([1, 3, 4]) #g2a
	strand_pointer_ra += x_strand_pointer_FN([4, 3, 2]) #g2b
	strand_pointer_ra += y_strand_pointer_FN([2, 1, 3]) #g3a
	strand_pointer_ra += y_strand_pointer_FN([4, 1, 2]) #g3b
	num_strands = len(strand_pointer_ra)
	print "num strands", num_strands

	subsub_ra = [[] for x in range(10)]
	sub_ra = [copy.deepcopy(subsub_ra) for y in range(10)]
	word_pointer_ra = [copy.deepcopy(sub_ra) for z in range(5)]
	for strand_num in range(num_strands):
		for [z, y, x] in strand_pointer_ra[strand_num]:
			word_pointer_ra[z][y][x].append(strand_num)

#	for z in range(5):
#		for y in range(10):
#			print z, y, word_pointer_ra[z][y]







#	Swaps on layers 1-4
	swap_attempt_num = 0
	num_swaps = 0
	num_swap_attempts = 3000
	old_word_ra = word_ra
	old_word_cache_ra = word_cache_ra[:]
	print "before swaps score", total_score_FN(old_word_ra, strand_pointer_ra)
	while swap_attempt_num < num_swap_attempts:
		new_word_ra = copy.deepcopy(old_word_ra)
		new_word_cache_ra = copy.deepcopy(old_word_cache_ra)
		z = random.choice(range(1, 5))
		y = 0
		x = 0
		while (x in [0, 9]) and (y in [0, 9]):
			x = random.choice(range(10))
			y = random.choice(range(10))
		new_word_cache_ra.append(new_word_ra[z][y][x])
		new_word_ra[z][y][x] = new_word_cache_ra.pop(0)
		score_dc = {'old':0, 'new':0}
		for [version, word_ra] in [['old', old_word_ra], ['new', new_word_ra]]:
			for strand_num in word_pointer_ra[z][y][x]:
				strand_seq = strand_seq_FN(word_ra, strand_pointer_ra[strand_num])
				score_dc[version] += self_dG_FN(strand_seq)
		if score_dc['new'] > score_dc['old']:
			num_swaps += 1
			old_word_ra = new_word_ra
			old_word_cache_ra = new_word_cache_ra
		if swap_attempt_num%100 == 0:
			print 'swap_attempt_num', swap_attempt_num, '\tnum_swaps', num_swaps
		swap_attempt_num += 1
	print "after swaps score", total_score_FN(old_word_ra, strand_pointer_ra)

	strand_dc = {'g1':[], 'g2a':[], 'g2b':[], 'g3a':[], 'g3b':[]}
	splint_dc = {'g1':[], 'g2a':[], 'g2b':[], 'g3a':[], 'g3b':[]}

	#Generate strands
	for strand_num in range(8):
		strand_dc['g1']  += [strand_seq_FN(word_ra, strand_pointer_ra[ 0 + strand_num])]
		strand_dc['g2a'] += [strand_seq_FN(word_ra, strand_pointer_ra[ 8 + strand_num])]
		strand_dc['g2b'] += [strand_seq_FN(word_ra, strand_pointer_ra[16 + strand_num])]
		strand_dc['g3a'] += [strand_seq_FN(word_ra, strand_pointer_ra[24 + strand_num])]
		strand_dc['g3b'] += [strand_seq_FN(word_ra, strand_pointer_ra[32 + strand_num])]
	
	
	sub_length = len(strand_dc['g1'][0])/3
	for subset in ['g1', 'g2a', 'g2b', 'g3a', 'g3b']:
		splint_dc[subset] = []
		for strand_num in range(8):
			strand_seq = strand_dc[subset][strand_num]
			for i in [1, 2]:
				splint_seq = comp_seq_FN(strand_seq[sub_length*i - 10:sub_length*i + 20])
				splint_dc[subset].append(splint_seq)
	
	for subset in ['g1', 'g2a', 'g2b', 'g3a', 'g3b']:
		print '\n' + subset
		for i in range(8):
			strand_seq = strand_dc[subset][i]
			for j in range(3):
				print strand_seq[sub_length*j:sub_length*(j + 1)]
		for i in range(16):
			print splint_dc[subset][i]
	
	
	
	

























