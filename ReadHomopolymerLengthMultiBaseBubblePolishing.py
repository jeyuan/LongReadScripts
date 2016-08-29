# -*- coding: utf-8 -*-
"""
Created on Sun Aug 28 17:48:02 2016

@author: jeffrey_yuan
"""

import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file',action='store')
    parser.add_argument('out_file',action='store')
    parser.add_argument('-v','--verbose',action='store_true',dest='v',default=False)
    parser.add_argument('-m','--minlen',action='store',dest='min_hopo_len',type=int,default=2)
    parser.add_argument('-p','--minprop',action='store',dest='min_prop',type=float,default=0.1)
    parser.add_argument('--minunit',action='store',dest='min_unit',type=int,default=2)
    parser.add_argument('--maxunit',action='store',dest='max_unit',type=int,default=3)
    
    args = parser.parse_args()
    unit_lens = range(args.min_unit,args.max_unit+1)
    
    all_heads,all_reads,all_cands = getBubblesFromAllMisha(args.input_file)
    with open(args.out_file,'w') as of:        
        for i,(head,reads,cand) in enumerate(zip(all_heads,all_reads,all_cands)):
            old_cand = cand
            for num_nuc in unit_lens:
                new_cand = read_homopolymer_length_multi_base_polishing(i,old_cand,reads,num_nuc,args.min_hopo_len,args.min_prop,args.v)
                old_cand = new_cand
            of.write('%s\n%s\n' % (head,new_cand))
    if args.v:
        print 'Results written to %s' % args.out_file

def test_polisher():
    '''Tests a set of example bubbles with errors to show corrections'''
    input_file = 'bls_polished_0802_4_error_bubbles_dinuc_examples.txt'
    out_file = 'test_error_examples.txt'
    v = True
    min_hopo_len = 2
    min_prop = 0.1
    min_unit = 2
    max_unit = 3
    unit_lens = range(min_unit,max_unit+1)
    
    out_lst = readErrorBubbles(input_file)
    
    with open(out_file,'w') as of:
        for r_i,bubble in enumerate(out_lst):
            header,errors,known,reads = bubble
            err_i,bub_i = header
            orig,hopo,gen = known        
            cand = hopo
            head = '>circular_0 %d' % bub_i
            i = bub_i
        
            old_cand = cand
            for num_nuc in unit_lens:
                new_cand = read_homopolymer_length_multi_base_polishing(i,old_cand,reads,num_nuc,min_hopo_len,min_prop,v)
                old_cand = new_cand
            of.write('%s\n%s\n' % (head,new_cand))
    if v:
        print 'Results written to %s' % out_file
    

def read_homopolymer_length_multi_base_polishing(bub_i,cand,reads,num_nuc,min_hopo_len,min_prop,v):
    '''Given a candidate sequence and a set of read segments corresponding
    to that sequence, polish all of the num_nuc homopolymer sequences
    in cand at least min_hopo_len long based on the longest length of that
    homopolymer that occurs in at least min_prop of reads. Output the new
    polished cand. v is a verbose option.
    Note: If a hompolymer of the same bases occurs multiple times in cand,
    only the LONGEST homopolymer will be polished.'''
    
    hopo_inds = find_homopolymer_indices_multiple_nucs(cand,min_hopo_len,num_nuc)
    if not hopo_inds:
        if v:
            print 'Bubble %d: No %d-unit homopolymers found.' % (bub_i,num_nuc)
        return cand
    else:
        prev_end = 0
        out_cand = []
        
        filtered_hopo_inds = filter_overlap_homopolymers(filter_dup_homopolymers(hopo_inds))
        for i,(start,end,unit) in enumerate(filtered_hopo_inds):
            len_counts = {}
            total_num = len(reads)
            
            for r in reads:
                longest_in_read = longest_hopo_length(r,num_nuc,unit)
                if longest_in_read not in len_counts:
                    len_counts[longest_in_read] = 1
                else:
                    len_counts[longest_in_read] += 1
            
            longest_hopo = ''
            if v:
                longest_len = 0
            for length in sorted(len_counts,reverse=True):
                if len_counts[length] >= total_num*min_prop:
                    longest_hopo = unit*length
                    if v:
                        longest_len = length
                    break
                
            out_cand.append(cand[prev_end:start])
            out_cand.append(longest_hopo)
            prev_end = end
            
            if v:
                cand_hopo_len = (end-start)/num_nuc
                if cand_hopo_len != longest_len:
                    print 'Bubble %d: %d-unit homopolymer %d changed from %d%s to %d%s' % (bub_i,num_nuc,i,cand_hopo_len,unit,longest_len,unit)
            
        out_cand.append(cand[prev_end:])
        return ''.join(out_cand)
    
def test_filter_overlap():
    seqs = ['','AAAAAAATAAAAA','ACACACACACACACAC','TACACACTGACACACACGACAC','TACACACGACACACAT','ACACGCGCGTGTGTG']
    min_hopo_len = 2
    num_nucs = [1,1,2,2,2,2]
    for (num_nuc,seq) in zip(num_nucs,seqs):
        hopo_inds = find_homopolymer_indices_multiple_nucs(seq,min_hopo_len,num_nuc)
        print seq,num_nuc,hopo_inds
        filtered_hopo_inds = filter_overlap_homopolymers(hopo_inds)
        print filtered_hopo_inds

def filter_overlap_homopolymers(hopo_inds):
    '''Given a list of hopo_inds in the format, [(start,end,unit),...],
    remove the shorter or later occurring of all homopolymers that are
    overlapping with each other.'''
    output_inds = []
    prev_hopo = ()
    prev_len = 0
    for start,end,unit in sorted(hopo_inds):
        if not prev_hopo:
            prev_hopo = (start,end,unit)
            prev_len = end-start
        else:
            if start >= prev_hopo[1]:
                output_inds.append(prev_hopo)
                prev_hopo = (start,end,unit)
                prev_len = end-start
            else:
                if prev_len >= end-start:
                    pass
                else:
                    prev_hopo = (start,end,unit)
                    prev_len = end-start
    if prev_hopo:
        output_inds.append(prev_hopo)
    return output_inds

def test_filter_dup():
    seqs = ['','AAAAAAATAAAAA','ACACACACACACACAC','TACACACTGACACACACGACAC','TACACACGACACACAT']
    min_hopo_len = 2
    num_nucs = [1,1,2,2,2]
    for (num_nuc,seq) in zip(num_nucs,seqs):
        hopo_inds = find_homopolymer_indices_multiple_nucs(seq,min_hopo_len,num_nuc)
        print seq,num_nuc,hopo_inds
        filtered_hopo_inds = filter_dup_homopolymers(hopo_inds)
        print filtered_hopo_inds

def filter_dup_homopolymers(hopo_inds):
    '''Given a list of hopo_inds in the format, [(start,end,unit),...],
    filter all but the longest homopolymers that have identical strings
    for unit and output filtered_hopo_inds.'''
    max_unit_len = {}
    max_unit_ind = {}
    for i,(start,end,unit) in enumerate(hopo_inds):
        if unit not in max_unit_len:
            max_unit_len[unit] = end-start
            max_unit_ind[unit] = i
        else:
            if end-start > max_unit_len[unit]:
                max_unit_len[unit] = end-start
                max_unit_ind[unit] = i
    
    return [(start,end,unit) for i,(start,end,unit) in enumerate(hopo_inds) if i == max_unit_ind[unit]]

def find_homopolymer_indices_multiple_nucs(seq,min_hopo_len,num_nuc):
    '''Given a DNA sequence seq, find the indices of all homopolymers that consist
    of >= min_hopo_len repeats of num_nuc units.
    Output: [(start,end,unit),...] for all homopolymers in seq.
    N.B. if num_nuc > 1, does not allow for hompolymers of the same single base.'''
    hopo_indices = []
    curr_units = ['' for b in range(num_nuc)]
    curr_lengths = [1 for b in range(num_nuc)]
    for i in range(len(seq)):
        j = i%num_nuc
        if seq[i:i+num_nuc] == curr_units[j]:
            curr_lengths[j] += 1
        else:
            if curr_units[j]:
                if curr_lengths[j] >= min_hopo_len and (num_nuc == 1 or curr_units[j] != len(curr_units[j])*curr_units[j][0]):
                    hopo_indices.append((i-(curr_lengths[j]*num_nuc),i,curr_units[j]))
                curr_lengths[j] = 1
            curr_units[j] = seq[i:i+num_nuc]
    j = len(seq)%num_nuc
    if curr_lengths[j] >= min_hopo_len and (num_nuc == 1 or curr_units[j] != len(curr_units[j])*curr_units[j][0]):
        hopo_indices.append((len(seq)-(curr_lengths[j]*num_nuc),len(seq),curr_units[j]))
    return hopo_indices

def longest_hopo_length(seq,num_nuc,unit):
    '''Given a DNA sequence seq and a string unit consisting of num_nuc nucleotides,
    output the length of the longest homopolymer in seq consisting of repeats of unit.
    For example, if seq = ACGCGCGT, num_nuc = 2 and unit = CG, output 3.'''
    longest = 0
    curr = 0
    i = 0
    while i < len(seq)-num_nuc:
        if seq[i:i+num_nuc] == unit:
            curr += 1
            i += num_nuc
        else:
            if curr > longest:
                longest = curr
            curr = 0
            i += 1
    if curr > longest:
        longest = curr
    return longest

def getBubblesFromAllMisha(all_file):
    '''Reads (potentially unsorted) bubbles in the form:
    >circular 0 30
    "consensus_sequence_0"
    >0
    "read 0"
    >1
    "read 1"
    ...
    >29
    "read 29"
    >circular 1 25
    "consensus_sequence_1"
    ...
    Returns the headers, bubbles, and consensus sequences in sorted order'''
    headers,segments = readAllLongFASTA(all_file)
    curr_ai = 0
    all_heads = []
    all_cands = []
    all_segs = []
    while curr_ai < len(headers):
        curr_head = headers[curr_ai]
        if curr_head[:7] == '>linear' or curr_head[:9] == '>circular':
            b_i,num_seg = map(int,curr_head.split(' ')[1:])
            cand = segments[curr_ai]
            read_seqs = segments[curr_ai+1:curr_ai+1+num_seg]
            all_heads.append(headers[curr_ai])
            all_cands.append(cand)
            all_segs.append(read_seqs)
            curr_ai += num_seg+1
        else:
            raise IOError('getBubblesFromAllMisha failed - curr_ai=%d is not a bubble head' % curr_ai)
    
    inds = [int(x.split(' ')[1]) for x in all_heads]
    ordered = sorted(zip(inds,all_heads,all_cands,all_segs))
    out_inds,out_heads,out_cands,out_segs = zip(*ordered)
    
    return out_heads,out_segs,out_cands

def readAllLongFASTA(file_name):
    '''Given a FASTA file of reads where the header starts with a '>' and the
    corresponding read can span multiple succeeding lines, return a tuple of
    headers,seqs.
    Note: reads all lines, even if they are blank'''
    seq_file = open(file_name,'r')
    headers = []
    seqs = []
    curr_seq = ''
    for i,line in enumerate(seq_file):
        line = line.strip()
        if line and line[0] == '>':
            headers.append(line)
            if i > 0:
                seqs.append(curr_seq.upper())
            curr_seq = ''
        else:
            curr_seq += line
    seqs.append(curr_seq.upper())
    
    seq_file.close()
    return headers,seqs

def readErrorBubbles(err_bubble_file):
    #(err_i,bub_i),[(offset,con_i,type,base,hopo,contig,genome),...],(cand,lik,hopo,gen),[seq1,seq2,...]
    out_lst = []
    first = True
    header = ()
    errors = []
    known = []
    bub_seqs = []
    with open(err_bubble_file,'r') as ebf:
        line = ebf.readline().strip()
        while line:
            if line[0] == '>':
                if not first:
                    out_lst.append((header,errors,known,bub_seqs))
                else:
                    first = False
                parts = line.split('|')
                header = (int(parts[1]),int(parts[3]))
                errors = []
                known = []
                bub_seqs = []
            else:
                parts = line.split('\t')
                if len(parts) > 2:
                    if parts[0] != 'Offset':
                        errors.append((int(parts[0]),int(parts[1]),parts[2],parts[3],parts[4],parts[5],parts[6]))
                elif len(parts) == 2:
                    if parts[0] in ['Cand','Lik','Hopo','Gen']:
                        known.append(parts[1])
                    else:
                        bub_seqs.append(parts[1])
            line = ebf.readline().strip()
        out_lst.append((header,errors,known,bub_seqs))
    return out_lst

if __name__ == '__main__':
    main()