#!/usr/bin/env python

import argparse
import gzip
import os
import pipes
import string

from Bio.SeqIO.QualityIO import FastqGeneralIterator


def make_strict_match(read_index_length, expected_indexes):

    sample_index_lengths = map(len, expected_indexes)
    sample_index_lengths = set(sample_index_lengths)

    if (len(sample_index_lengths) == 1):

        sample_index_length = sample_index_lengths.pop()

        if (read_index_length == sample_index_length):

            def strict_match_same_length(read_index, sample_indexes):

                if (read_index in sample_indexes):
                    return read_index
                else:
                    return None

            return strict_match_same_length
   
        elif (read_index_length > sample_index_length):

            def strict_match_different_length(read_index, sample_indexes):
       
                read_index_trimmed = read_index[:sample_index_length]
 
                if (read_index_trimmed in sample_indexes):
                    return read_index_trimmed
                else:
                    return None
 
            return strict_match_different_length

        else:
            sys.exit('sample_index_length greater than read_index_length')
    else:

        def strict_match_variable_length(read_index, sample_indexes):

            matches = [ ]

            for sample_index in sample_indexes:

                is_a_match = 1

                for pos in (range(len(sample_index))):
                    if (read_index[pos] != sample_index[pos]):
                        is_a_match = 0
                        break

                if (is_a_match):
                    matches.append(sample_index)
                
            if (len(matches) != 1):
                return None
            else:
                return matches.pop()
        
        return strict_match_variable_length

def make_fuzzy_match(mismatches_allowed, n_penalty):

    def fuzzy_match(read_index, sample_indexes):

        matches = [ ]

        for sample_index in sample_indexes:

            mismatches = 0

            for pos in (range(len(sample_index))): 
                if (read_index[pos] != sample_index[pos]):
                    if (read_index[pos] == 'N'):
                        mismatches += n_penalty
                    else:
                        mismatches += 1

            if (mismatches <= mismatches_allowed):
                matches.append(sample_index) 
    
        if (len(matches) != 1):
            return None
        else:
            return matches.pop()          

    return fuzzy_match

def lookup_index_cycles(index_fn):

    iterator = FastqGeneralIterator(gzip.open(args.index_read_file))
        
    name, seq, qual = iterator.next()

    return len(seq)


parser = argparse.ArgumentParser(description='FASTQ demultiplexer')

parser.add_argument('--read1-file',
                    '-read1-file',
                    '-r1',
                    required=True,
                    help='read1 file',
                   )

parser.add_argument('--read2-file',
                    '-read2-file',
                    '-r2',
                    help='read2 file',
                   )

parser.add_argument('--read3-file',
                    '-read3-file',
                    '-r3',
                    help='read3 file',
                   )

parser.add_argument('--index-read-file',
                    '-index-read-file',
                    '-i',
                    required=True,
                    help='index read file',
                   )

parser.add_argument('--tag-file',
                    '-tag-file',
                    '-t',
                    required=True,
                    help='file with expected indexes',
                   )

parser.add_argument('--mismatches',
                    '-mismatches',
                    '-m',
                    type=int,
                    default=0,
                    help='mismatches allowed')

parser.add_argument('--n-penalty',
                    '-n-penalty',
                    '-n',
                    type=int,
                    default=1,
                    help='mismatch penalty for Ns in index')

args = parser.parse_args()

indexes  = set()

tag_fh              = open(args.tag_file, 'r')
sample_index_length = None

for line in tag_fh:
    index               = line.rstrip(os.linesep)
    sample_index_length = len(index)
    indexes.add(index)
 
read1_index_out_fh = { }
read2_index_out_fh = { }
read3_index_out_fh = { }

pipe_template = pipes.Template()
pipe_template.append('gzip -c', '--')

for index in indexes.union(set(['unknown'])): 

    read1_index_out_fn        = string.join([index, '1', 'fq', 'gz'], '.')
    read1_index_out_fh[index] = pipe_template.open(read1_index_out_fn, 'w')  

    if (args.read2_file is not None):
        read2_index_out_fn        = string.join([index, '2', 'fq', 'gz'], '.')
        read2_index_out_fh[index] = pipe_template.open(read2_index_out_fn, 'w')
    if (args.read3_file is not None):
        read3_index_out_fn        = string.join([index, '3', 'fq', 'gz'], '.')
        read3_index_out_fh[index] = pipe_template.open(read3_index_out_fn, 'w')

read_index_length = lookup_index_cycles(args.index_read_file)

iteratori = FastqGeneralIterator(gzip.open(args.index_read_file))

iterator2 = None
iterator3 = None

if (args.read2_file != None):
    iterator2 = FastqGeneralIterator(gzip.open(args.read2_file))
if (args.read3_file != None):
    iterator3 = FastqGeneralIterator(gzip.open(args.read3_file))

func = None

if (args.mismatches > 0):
    func = make_fuzzy_match(args.mismatches, args.n_penalty)
else:
    func = make_strict_match(read_index_length, indexes)

for rname1, seq1, qual1 in FastqGeneralIterator(gzip.open(args.read1_file)):

    rnamei, seqi, quali = iteratori.next()

    out_index = func(seqi, indexes)
 
    if (out_index is None):
        out_index = 'unknown'
    
    if (iterator2 is not None):

        rname2, seq2, qual2 = iterator2.next()

    if (iterator3 is not None):

        rname3, seq3, qual3 = iterator3.next()

        read3_index_out_fh[out_index].write('@' + rname3 + ":" + seq3 + "\n")
        read3_index_out_fh[out_index].write(seq3 + "\n")
        read3_index_out_fh[out_index].write('+' + "\n")
        read3_index_out_fh[out_index].write(qual3 + "\n")

    nm1, stuff1=rname1.split();
    nm2, stuff2=rname2.split();

    read1_index_out_fh[out_index].write('@' + nm1 + ":" + seq3 + "\n")
    read1_index_out_fh[out_index].write(seq1 + "\n")
    read1_index_out_fh[out_index].write('+' + "\n")
    read1_index_out_fh[out_index].write(qual1 + "\n")

    read2_index_out_fh[out_index].write('@' + nm2 + ":" + seq3 + "\n")
    read2_index_out_fh[out_index].write(seq2 + "\n")
    read2_index_out_fh[out_index].write('+' + "\n")
    read2_index_out_fh[out_index].write(qual2 + "\n")



for index in read1_index_out_fh:
    read1_index_out_fh[index].close()   

for index in read2_index_out_fh:
    read2_index_out_fh[index].close() 

