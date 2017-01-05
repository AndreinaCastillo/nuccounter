####Import library
import fileinput
#import os, sys
import glob, os, sys
import collections
import itertools
from os.path import basename
from operator import itemgetter
from itertools import islice
from itertools import groupby
import numpy as np
# from plotly import __version__
# from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
# from plotly.graph_objs import Bar, Scatter, Figure, Layout
# init_notebook_mode(connected=True)



####Global variable 
window_size = 2000



##File input/path

fasta_line_source = fileinput.FileInput(os.path.expanduser('chromosome_16911_5.faa'))
file_name = (os.path.splitext('chromosome_16911_5.faa')[0])
file_name = file_name + '_out.txt'
# print(file_name)



###Read each line on the input file and splits sequence headers and sequenced (es a string) into a dictionary

def lines_to_contigs(lines):
    contig_dictionary = None
    current_contig_name = None    
    current_contig_sequence = []
    lines_with_eof_token = itertools.chain(lines, ['>EOF'])
    for line in lines_with_eof_token:
        if line.startswith('>'):
            #print('line: {}'.format(line))
            if current_contig_name:
                #print(line)
                #print('current_contig_name: {}'.format(current_contig_name))
                contig_sequence_string = ''.join(current_contig_sequence)
                contig_sequence_string = contig_sequence_string.upper( )
                contig_dictionary = {'name' : current_contig_name, 'sequence' : contig_sequence_string}
                yield contig_dictionary
                current_contig_sequence = []
            current_contig_name = line[1:].strip()
        else:  
            current_contig_sequence.append(line.strip())
            
generator_lines_to_contigs = lines_to_contigs(fasta_line_source)
# for i in generator_lines_to_contigs:
#     print(i)


###Groups A and T into AT, and C and G into GC. Other nucleotides (including ambiguous characters are added to their own category) 

def group_nucleotides(nucleotide_seqs):
    nucleotide_group_map = {'A': 'AT','T': 'AT','G': 'GC','C': 'GC','N': 'N','X': 'X', 
                            'M':'M','R':'R','Y':'Y','S':'S','W':'W','K':'K','B':'B','V':'V','D':'D','H':'H'}
    for nucleotide_seq in nucleotide_seqs:
        nucleotide_groups = (nucleotide_group_map[base] for base in nucleotide_seq['sequence'])
        yield {'name': nucleotide_seq['name'], 'sequence': nucleotide_groups}

nucleotides_seq = group_nucleotides(lines_to_contigs(fasta_line_source))
# for i in nucleotides_seq:
#     print(list(islice(i['sequence'], 10)))


### Creates a sliding window (size of the window is a global variable). 
### Calculates %AT, %GC, %N and %X inside the window. It also calculates %of ambiguous nucleotides if found.
### After first window, the first nucleotide inside the window is droped and the next in the sequence is added. 
### counter updates to discount the removed nucleotide and count the new nucleotide
### Percentages are recalculated.
### Each window is yielded to a dictionary.

def sliding_percentages(nucleos, window_size, unique_nucleotides):
    for nucleo in nucleos:
        sliding_window = collections.deque(itertools.islice(nucleo['sequence'], window_size), maxlen=window_size)
        counter = collections.Counter({nucleotides: 0 for nucleotides in unique_nucleotides}) 
        counter.update(sliding_window)
        window_percentage = {nucleotides: (counter / window_size)*100 for nucleotides, counter in counter.items()}
        w_start=0
        for base in nucleo['sequence']:
            w_start = w_start+1
            itertools.islice(nucleo['sequence'], window_size, None)
            trailing_nucleotide = sliding_window.popleft()
            counter.subtract([trailing_nucleotide])
            w_end = w_start+window_size
            sliding_window.append(base)
            counter.update([base])
            window_percentage = {nucleotides: (counter / window_size)*100 for nucleotides, counter in counter.items()}
            yield {'name': nucleo['name'], 'sequence': window_percentage, 'start': w_start, 'end': w_end}
        
nucleotides4 = sliding_percentages(group_nucleotides(lines_to_contigs(fasta_line_source)),window_size,('AT', 'GC', 'N', 'X'))
# for i in nucleotides4:
#     print(i)



### Turns all values calculated on the previous function and found on the dictionary into list
### Output list as a tab delimited format into a .txt file

def contig_list_creator(nucleos):    
    #sorted_iterator_list = sorted(nucleos, key=itemgetter('name'))

    calc_out_file = open(file_name, 'w')
    name = ''
    for nucleo in nucleos:
        if nucleo['name'] != name:
            calc_out_file.write('#Chromosome/contig name, window_start_position, window_end_position, %AT, %GC, %N, %X\n')
            name = nucleo['name']
        calc_out_file.write(name + '\t')
        calc_out_file.write(str(nucleo['start']) + '\t')
        calc_out_file.write(str(nucleo['end']) + '\t')
        calc_out_file.write(str(nucleo['sequence']['AT']) + '\t')
        calc_out_file.write(str(nucleo['sequence']['GC']) + '\t')
        calc_out_file.write(str(nucleo['sequence']['N']) + '\t')
        calc_out_file.write(str(nucleo['sequence']['X']) + '\n')
        
    calc_out_file.close()
                
    return

contig_list_creator(nucleotides4)
# nucleotides5 = contig_list_creator(sliding_percentages(group_nucleotides(sequence_split(lines_to_contigs(fasta_line_source))),window_size,('AT', 'GC', 'N', 'X')))
# for i in nucleotides5:
#     print(i)
