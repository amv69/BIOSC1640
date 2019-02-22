"""
General utilities for processing biological data

Author: Miler Lee
Last modified: 18 Feb 2012
"""

import random
import gzip
import bz2

# Config settings
#import config
#DATA_ROOT = config.DATA_ROOT
#GENOMES = config.GENOMES

#GENOMIC_SEQ_FILES = {'danRer7': GENOMES + 'danRer7.fa', 'danrer7': GENOMES + 'danRer7.fa', 'xentro3': GENOMES + 'xenTro3.fa', 'xenTro3': GENOMES + 'xenTro3.fa', 'orylat2': GENOMES + 'oryLat2.fa', 'oryLat2': GENOMES + 'oryLat2.fa'}
GENOMIC_SEQ_FILES = {'danRer7': '/data/users/miler/data/genomes/full/danRer7.fa'}

##OLD
GENOMIC_SEQ_ROOTS = {'rn4':'/data15/public.repo/genomes/rn4/', \
                     'mm9': '/data15/public.repo/genomes/mm9/'}
GENOMIC_RMSK_SEQ_ROOTS = {'rn4': '/data15/public.repo/genomes/rn4_rmsk/'}

#######################
# Nucleotide processing
#######################

rna_bases = ['A', 'C', 'G', 'U']

def random_rna_base():
    """
    Returns an upper case rna base uniformly.  Faster than random_base()
    """

    return random.choice(rna_bases)


def random_base(upper=False, rna=False, exclude=''):
    """
    Returns a uniform pseudorandom DNA or RNA nucleotide.  exclude is a
    user-specified string of nucleotide characters, none of which will
    be returned as a random base.
    """

    import random
    if rna:
        fourth_base = 'u'
    else:
        fourth_base = 't'

    bases = ['a', 'c', 'g', fourth_base]

    try:
        if not exclude:
            bases.remove(exclude.lower())
    except:
        pass

    new_base = random.choice(bases)

    if upper:
        return new_base.upper()
    else:
        return new_base


def random_base_p(probs = [0.25, 0.25, 0.25, 0.25]):
    """
    Returns a random rna base according to the specified distribution,
    encoded as a list of probabilities for A, C, and G, and U in that order.
    Default is uniform probability for all four bases, i.e. identical
    to random_base()

    Probabilities are normalized to 1 and missing probs are assumed to be 0

    Due to precision issues, it is possible but rare for base selection to
    fail, in which case 'n' is returned.
    """

    padded_probs = probs + [0] * (4-len(probs))
    cumulative = sum(padded_probs)

    prob_bases = list(zip(padded_probs, ['a', 'c', 'g', 'u']))

    prob = random.random() * cumulative
    cum_sum = 0
    for base_prob, base in prob_bases:
        if not base_prob:
            continue
        cum_sum += base_prob
        if prob <= cum_sum:
            return base

    return 'n'
    

###string.maketrans('acgturyACGTURY', 'tgcaayrTGCAAYR')
DNA_TRANS = '\x00\x01\x02\x03\x04\x05\x06\x07\x08\t\n\x0b\x0c\r\x0e\x0f\x10\x11\x12\x13\x14\x15\x16\x17\x18\x19\x1a\x1b\x1c\x1d\x1e\x1f !"#$%&\'()*+,-./0123456789:;<=>?@TBGDEFCHIJKLMNOPQYSAAVWXRZ[\\]^_`tbgdefchijklmnopqysaavwxrz{|}~\x7f\x80\x81\x82\x83\x84\x85\x86\x87\x88\x89\x8a\x8b\x8c\x8d\x8e\x8f\x90\x91\x92\x93\x94\x95\x96\x97\x98\x99\x9a\x9b\x9c\x9d\x9e\x9f\xa0\xa1\xa2\xa3\xa4\xa5\xa6\xa7\xa8\xa9\xaa\xab\xac\xad\xae\xaf\xb0\xb1\xb2\xb3\xb4\xb5\xb6\xb7\xb8\xb9\xba\xbb\xbc\xbd\xbe\xbf\xc0\xc1\xc2\xc3\xc4\xc5\xc6\xc7\xc8\xc9\xca\xcb\xcc\xcd\xce\xcf\xd0\xd1\xd2\xd3\xd4\xd5\xd6\xd7\xd8\xd9\xda\xdb\xdc\xdd\xde\xdf\xe0\xe1\xe2\xe3\xe4\xe5\xe6\xe7\xe8\xe9\xea\xeb\xec\xed\xee\xef\xf0\xf1\xf2\xf3\xf4\xf5\xf6\xf7\xf8\xf9\xfa\xfb\xfc\xfd\xfe\xff'
RNA_TRANS = '\x00\x01\x02\x03\x04\x05\x06\x07\x08\t\n\x0b\x0c\r\x0e\x0f\x10\x11\x12\x13\x14\x15\x16\x17\x18\x19\x1a\x1b\x1c\x1d\x1e\x1f !"#$%&\'()*+,-./0123456789:;<=>?@UBGDEFCHIJKLMNOPQYSAAVWXRZ[\\]^_`ubgdefchijklmnopqysaavwxrz{|}~\x7f\x80\x81\x82\x83\x84\x85\x86\x87\x88\x89\x8a\x8b\x8c\x8d\x8e\x8f\x90\x91\x92\x93\x94\x95\x96\x97\x98\x99\x9a\x9b\x9c\x9d\x9e\x9f\xa0\xa1\xa2\xa3\xa4\xa5\xa6\xa7\xa8\xa9\xaa\xab\xac\xad\xae\xaf\xb0\xb1\xb2\xb3\xb4\xb5\xb6\xb7\xb8\xb9\xba\xbb\xbc\xbd\xbe\xbf\xc0\xc1\xc2\xc3\xc4\xc5\xc6\xc7\xc8\xc9\xca\xcb\xcc\xcd\xce\xcf\xd0\xd1\xd2\xd3\xd4\xd5\xd6\xd7\xd8\xd9\xda\xdb\xdc\xdd\xde\xdf\xe0\xe1\xe2\xe3\xe4\xe5\xe6\xe7\xe8\xe9\xea\xeb\xec\xed\xee\xef\xf0\xf1\xf2\xf3\xf4\xf5\xf6\xf7\xf8\xf9\xfa\xfb\xfc\xfd\xfe\xff'


def rc_dna_fast(sequence):
    """
    Fast reverse complement assuming DNA nucleotides
    """

    result = sequence.translate(DNA_TRANS)
    return result[::-1]

def rc_rna_fast(sequence):
    """
    Fast reverse complement assuming RNA nucleotides
    """

    result = sequence.translate(RNA_TRANS)
    return result[::-1]


def rc(sequence):
    """
    Alias for reversecomp
    """

    return reversecomp(sequence)

def reversecomp(sequence):
    """
    Returns DNA reverse complement of input string.

    Operates by reversing the string and
    calling watsoncrick on every element
    """

    result = list(map(watsoncrick, sequence))
    result.reverse()
    return "".join(result)

def complement(sequence):
    """
    Similar to reversecomp() but does not reverse.
    Resulting sequence is thus read 3'-5'
    """

    result = list(map(watsoncrick, sequence))
    return ''.join(result)


def watsoncrick(base, rna=False):
    """
    Returns the watson crick complementary base to the input.

    A<->T, C<->G, U->A, Y<->R, all others->N
    case of the original base is maintained
    """

    lc = base.islower()
    base = base.upper()

    if base == 'A':
        if rna:
            rcbase = 'U'
        else:
            rcbase = 'T'
    elif base == 'T':
        rcbase = 'A'
    elif base == 'U':
        rcbase = 'A'
    elif base == 'C':
        rcbase = 'G'
    elif base == 'G':
        rcbase = 'C'
    elif base == 'Y':
        rcbase = 'R'
    elif base == 'R':
        rcbase = 'Y'
    elif base == 'S':
        rcbase = 'W'
    elif base == 'W':
        rcbase = 'S'
    elif base == 'K':
        rcbase = 'M'
    elif base == 'M':
        rcbase = 'K'
    elif base == 'B':
        rcbase = 'V'
    elif base == 'D':
        rcbase = 'H'
    elif base == 'H':
        rcbase = 'D'
    elif base == 'V':
        rcbase = 'B'
    else:
        rcbase = 'N'
        
    if lc:
        rcbase = rcbase.lower()

    return rcbase


def rna(seq):
    """
    Converts all T's into U's
    """

    seq = seq.replace('T', 'U')
    seq = seq.replace('t', 'u')    
    return seq

def dna(seq):
    """
    Converts all U's into T's
    """

    seq = seq.replace('U', 'T')
    seq = seq.replace('u', 't')    
    return seq

def disambig(seq, dna = False):
    """
    Given the input sequence, disambiguates all IUPAC wildcard
    characters with a randomly consistent nucleotide
    """
    if dna:
        t = 'T'
    else:
        t = 'U'

    equiv = {'R': ['A', 'G'], 'Y': ['C', t], 'M': ['A', 'C'], 'K': ['G', t],\
             'S': ['C', 'G'], 'W': ['A', t], 'B': ['C', 'G', t], \
             'D': ['A', 'G', t], 'H': ['A', 'C', t], 'V': ['A', 'C', 'G'],\
             'N': ['A', 'C', 'G', t]}

    new_seq = []
    for nt in seq:
        nt_u = nt.upper()
        if nt_u in equiv:
            new_nt = random.choice(equiv[nt_u])
            if nt == nt_u:
                new_seq.append(new_nt)
            else:
                new_seq.append(new_nt.lower())
        else:
            new_seq.append(nt)
    return ''.join(new_seq)


genetic_code = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S', 'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I', 'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TGA': '*', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C', 'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'}


def translate(dna_sequence, frame = 0):
    """
    Standard genetic code nt to amino acid. Stop codon becomes *.
    Frame = 1 means skip the first nt, 2 means skip the first 2.
    Extra nts at the end of the sequence are ignored. Ambiguous
    codons (e.g., containing N characters) are translated to ?

    Not a particularly efficient implementation
    """

    dna_sequence = dna_sequence.upper()
    protein = []
    
    for i in range(frame, len(dna_sequence)-2, 3):
        codon = dna_sequence[i:(i+3)]
        try:
            aa = genetic_code[codon]
        except:
            aa = '?'
        protein.append(aa)
    return ''.join(protein)


def reading_frame(dna_sequence):
    """
    Returns the reading frame of the sequence, assuming that it
    encodes a maximal ORF. Frame 0 = starts at the first
    nt, 1 = skip the first nt, 2 = skip the first 2. -1 means
    the sequence has no maximal ORF. No checking for ATGs
    is performed.

    If there is more than one valid ORF, the lowest frame
    is favored.
    """
    
    stop_codons = {'TAG': 1, 'TGA': 1, 'TAA': 1}
    
    dna_sequence = dna_sequence.upper()
    
    for frame in range(3):
        for k in range(frame, len(dna_sequence)-3, 3): #skip the last 3 nts always
            codon = dna_sequence[k:(k+3)]
            if codon in stop_codons:
                break
        else:
            return frame
    return -1
    

#######################
# FASTA and FASTQ file processing
#######################

def read_fasta(filename):
    """
    Returns the contents of a fasta file in a list of (id, sequence)
    tuples. Empty list returned if there are no fasta sequences in the file
    """

    a = fasta_reader(filename)
    seqs = []
    while a.has_next():
        seqs.append(next(a))
    return seqs

def count_fasta(filename):
    """
    Returns the number of sequences in the fasta file
    """
    count = 0
    a = fasta_reader(filename)
    while a.has_next():
        next(a)
        count += 1
    return count
    
def length_fasta(filename, outfile='', print_lengths = False, write_ids = True):
    """
    Prints the nt length of each entry in the fasta file,
    optionally to an output file; returns the lengths in an array
    """

    length_array = []

    a = fasta_reader(filename)
    wr = None
    if outfile != '':
        wr = open(outfile, 'w')
    while a.has_next():
        len, id = a.next_length(True)
        length_array.append(len)
        if wr != None:
            if write_ids:
                wr.write('%s\t%d\n' % (id, len))
            else:
                wr.write(str(len) + '\n')
        if print_lengths:
            print(id, len)        
    a.close()
    if wr != None:
        wr.close()
    return length_array


def length_fasta_ids(filename, as_hash = False, id_first_word_only = False):
    """
    Same as length_fasta() except returns lengths tupled with ids

    If as_hash = True, the return is a hashtable keyed by id. No
    checking for duplicate ids.

    If id_first_word_only, the ID is truncated to be the first word
    """
    
    a = fasta_reader(filename)
    
    if as_hash:
        length_array = {}
        while a.has_next():
            len, id = a.next_length(True)
            if id_first_word_only:
                id = id.split()[0]
            
            length_array[id] = len
    else:
        length_array = []

        while a.has_next():
            len, id = a.next_length(True)
            if id_first_word_only:
                id = id.split()[0]
            length_array.append((id, len))
            
    a.close()
    return length_array



def fasta_ids(filename, first_word_only = True):
    """
    Processes a fasta file and returns a list of the identifiers
    in order.  If first_word_only is True, then only the first
    word of the identifier is returned (everything before the
    first white space character)
    """

    entries = read_fasta(filename)
    if first_word_only:
        return [x[0].split()[0] for x in entries]
    else:
        return [x[0] for x in entries]


def shuffle_fasta(infile, outfile):
    """
    Writes a new fasta file that is a random permutation of the
    original.  Potentially very slow and memory intensive.
    """

    seqs = read_fasta(infile)
    wr = fasta_writer(outfile)
    while len(seqs) > 0:
        index = random.randint(0, len(seqs)-1)
        entry = seqs.pop(index)
        wr.write(entry)
    wr.close()


def random_subset_fasta(infile, outfile, n):
    """
    Writes a new fasta file that is a random subset
    of size min(n, size(infile)) sequences from infile
    """ 

    seqs = read_fasta(infile)
    wr = fasta_writer(outfile)
    if n > len(seqs):
        n = len(seqs)
    indices = random.sample(range(len(seqs)), n)
    for index in indices:
        wr.write(seqs[index])
    wr.close()


def map_fasta(infile, outfile, mapfn):
    """
    Given a fasta file, outputs it applying the mapfn
    to each nucleotide in each sequence
    """

    a = fasta_reader(infile)
    wr = fasta_writer(outfile)
    a.set_filter(mapfn)
    while a.has_next():
        wr.write(next(a))
    a.close()
    wr.close()


def remove_n_fasta(infile, outfile):
    """
    Rewrites the input fasta file concatenating all the non-N bases
    """
    map_fasta(infile, outfile, remove_n_filter_fn)


def disambig_fasta(infile, outfile, dna = False):
    """
    Replaces all IUPAC wildcard nucleotides with consistent canonical
    nucleotides
    """

    a = fasta_reader(infile)
    wr = fasta_writer(outfile)
    while a.has_next():
        id, seq = next(a)
        wr.write(id, disambig(seq))
    a.close()
    wr.close()


def filter_fasta(infile, outfile, idfile):
    """
    The id file contains identifiers.  Entries in the fasta infile are output
    to the outfile if the id is contained in the idfile (ie, if the id
    string is found anywhere in the identifier...be careful of prefixes)
    """

    file = open(idfile, 'r')
    ids = {}
    try:
        for line in file:
            ids[line.rstrip()] = 1
    except:
        pass
    file.close()

    a = fasta_reader(infile)
    wr = fasta_writer(outfile)
    while a.has_next():
        next = next(a)
        for key in list(ids.keys()):
            if next[0].find(key) > -1:
                ids[key] = 0
                wr.write(next)
                break
    a.close()
    wr.close()
    for key in list(ids.keys()):
        if ids[key] == 1:
            print(key)


def extract_fasta_low_mem(infile, outfile, idfile):
    """
    Entire ID must match
    """

    file = open(idfile, 'r')
    ids = {}
    try:
        for line in file:
            id = line.rstrip()
            if id == '':
                break
            if id in ids:
                print("repeat: " + id)
            ids[id] = 0
    except:
        pass
    file.close()

    wr = fasta_writer(outfile)
    fr = fasta_reader(infile)

    while fr.has_next():
        header, seq = next(fr)
        if header in ids:
            wr.write(header, seq)
            ids[header] = 1
    fr.close()
    wr.close()


def extract_fasta(infile, outfile, idfile, prefix=True):
    """
    For each id, looks for the corresponding identifier in the in fasta
    file, and outputs that entry into the outfile, otherwise prints
    not found.  prefix = True if the id (plus a white space) only needs to
    match the beginning of the identifier.  Slow.
    """

    file = open(idfile, 'r')
    ids = []
    try:
        for line in file:
            id = line.rstrip()
            if id == '':
                break
            ids.append([id, '', ''])
    except:
        pass
    file.close()

    wr = fasta_writer(outfile)
    fr = fasta_reader(infile)
    while fr.has_next():
        header, seq = next(fr)


        for tuple in ids:
            id = tuple[0]
            if prefix == True:
                index = header.find(id)
                if index > -1:
                    if index + len(id) >= len(header) or \
                           header[index+len(id)].isspace():
                        if tuple[1] != '':
                            print('Collision: ' + id)
                        tuple[1] = header
                        tuple[2] = seq
            else:
                if header.find(id) > -1:
                    if tuple[1] != '':
                        print('Collision: ' + id)
                    tuple[1] = header
                    tuple[2] = seq
    fr.close()

    for id, header, seq in ids:
        if header != '':
            wr.write(header, seq)
        else:
            print(id)
    wr.close()



def filter_fasta_by_id(infile, outfile, keyword, inverse=False):
    """
    Writes fasta entries from infile to outfile if keyword is
    present in the identifier.

    If inverse is True, then entries are written to outfile if
    keyword IS NOT present in the identifier
    """

    a = fasta_reader(infile)
    wr = fasta_writer(outfile)
    while a.has_next():
        next = next(a)
        if next[0].find(keyword) > -1:
            if not inverse:
                wr.write(next)
        else:
            if inverse:
                wr.write(next)
    a.close()
    wr.close()


def clean_up_rfam(infile, outfile, identifier=''):
    """
    Specific method that takes a fasta seed alignment file from rfam
    and outputs it minus gaps and prepending the rna class name in
    front of the identifier
    """

    seqs = read_fasta(infile)

    fw = fasta_writer(outfile)
    for seq in seqs:
        fw.write(identifier + ';' + seq[0], seq[1].replace('.', ''))
    fw.close()


def fasta_repeats(infile, outfile):
    """
    Removes repeated entries in the fasta file.  The first
    occurrence of an id is kept.  No sequence checking is performed
    """

    seqs = read_fasta(infile)
    added = {}
    fw = fasta_writer(outfile)
    for id, seq in seqs:
        if id not in added:
            added[id] = 1
            fw.write(id, seq)
    fw.close()


def fasta_remains(infile1, infile2, outfile):
    """
    Outputs all entries in infile1 not in infile2 based on id equality only
    """

    seqs1 = read_fasta(infile1)
    seqs2 = read_fasta(infile2)

    seq2_hash = dict(seqs2)

    fw = fasta_writer(outfile)
    for id, seq in seqs1:
        if id not in seq2_hash:
            fw.write(id, seq)
    fw.close()


def fasta_intersection(file1, file2):
    """
    Returns a list of ids contained in both fasta files.  No
    sequence checking is performed.
    """

    ids1 = fasta_ids(file1, False)
    ids2 = fasta_ids(file2, False)

    id1_hash = {}
    repeats = []
    for id in ids1:
        id1_hash[id] = 1
    for id in ids2:
        if id in id1_hash:
            repeats.append(id)
    return repeats


def fasta_union(infile_list, outfile):
    """
    Writes a single fasta file joining the
    input files together. Any read present in
    more that none input file is listed only once in the
    outfile. IDs are assumed to be consistent
    across the files
    """
    seqs = {}
    for infile in infile_list:
        f = fasta_reader(infile)
        while f.has_next():
            id, seq = next(f)
            seqs[id] = seq
        f.close()

    fw = open(outfile)
    
    for id, seq in list(seqs.items()):
        fw.write(id, seq)
    fw.close()
    

def split_fasta(fastafile, seqs_per_file = 1000000):
    """
    Divides a fasta file into many files with at most seqs_per_file sequences
    per file.  Filenames are of the form <ORIG_FASTAFILE_NAME>.px where x is
    the xth partition (starting at 0).
    """

    fnames = []
    outfilename = fastafile + '.p%d'
    curr_fw_id = -1
    fw = None
    
    fr = fasta_reader(fastafile)
    counter = 0
    while fr.has_next():
        if counter % seqs_per_file == 0:
            if fw:
                fw.close()
            curr_fw_id += 1
            fw = fasta_writer(outfilename % curr_fw_id)
            fnames.append(outfilename % curr_fw_id)

        fw.write(next(fr))
        counter += 1
    fr.close()
    return fnames


def balance_fasta(fastafile, npartitions = 2, single_line = False):
    """
    Reads in the fasta file, approximately evenly splits
    it into npartitions parts using a greedy approach.
    Each partition is output to a fasta file with suffix
    '.px' where x is the xth partition

    if single_line, ensures each fasta entry has no line breaks
    """

    seqs = read_fasta(fastafile)
    nt_count = 0
    for id, seq in seqs:
        nt_count += len(seq)

    partition_size = nt_count / npartitions
    current_count = 0
    partition_number = 0
    for id, seq in seqs:
        if current_count == 0:
            if single_line:
                width = 0
            else:
                width = 60
            fw = fasta_writer(fastafile + '.p' + str(partition_number), width)

        current_count += len(seq)
        fw.write(id, seq)

        if current_count > partition_size:
            current_count = 0
            partition_number += 1
            fw.close()


def trim_fasta(fasta_in, fasta_out, seq_trim_amt):
    """
    Given a fasta file, trims seq_trim_amt number of nts
    from each sequence and outputs to the fasta_out file,
    with no changes to the IDs.  Sequences shorter than
    or equal to seq_trim_amt are converted to empty strings.

    To specify final length of the sequence rather than
    number of nts to trim, use truncate_fasta()
    """
    trim_index = -1 * seq_trim_amt

    fr = fasta_reader(fasta_in)
    fw = fasta_writer(fasta_out)

    while fr.has_next():
        id, seq = next(fr)
        seq = seq[:trim_index]
        fw.write(id, seq)
    fr.close()
    fw.close()


def truncate_fasta(fasta_in, fasta_out, max_seq_len = 50, \
                   discard_shorter = False):
    """
    Given a fasta file, trims each sequence to a maximum
    length of max_seq_len nucleotides.  Sequences shorter
    than max_seq_len are unaffected unless discard_shorter is True, then
    they are discarded.  Output to the fasta_out file with no changes to
    the IDs.

    To specify a fixed number of nts to trim from each sequence
    regardless of length, use trim_fasta()
    """

    fr = fasta_reader(fasta_in)
    fw = fasta_writer(fasta_out)
    while fr.has_next():
        id, seq = next(fr)
        seq = seq[:max_seq_len]
        if len(seq) == max_seq_len or not discard_shorter:
            fw.write(id, seq)
    fr.close()
    fw.close()

    
def join_fasta(fasta_file_list, fasta_out):
    """
    Simply combines the fasta files together into a new fasta file,
    but renames each sequence identifier to be the filename of the
    original file (if multiple seqs in a file, _2, _3, etc are appended)
    """
    
    of = fasta_writer(fasta_out)
    
    for f in fasta_file_list:
        seqs = read_fasta(f)
        new_id = f.split('/')[-1].split('.')[0]
        
        for i, (id, seq) in enumerate(seqs):
            if i:            
                of.write(new_id + '_%d' % (i+1), seq)
            else:
                of.write(new_id, seq)                
    of.close()
                
class fasta_reader:
    """
    Lightweight class for incrementally reading fasta files.

    Supports reading directly from properly named
    gzipped (.gz or .z) or bzip2ed (.bz2) files.
    """

    file = None
    nextheader=''
    filter_fn = None
    
    def __init__(self, filename):
        try:
            if filename.endswith('.gz') or filename.endswith('.z'):
                self.file = gzip.open(filename, 'rb')
            elif filename.endswith('.bz2'):
                self.file = bz2.BZ2File(filename, 'rb')
            else:
                self.file = open(filename, 'r')

            # fast forward to the first entry
            while 1:
                line = self.file.readline()
                if line == '':
                    self.close()
                    return
                elif line[0] == '>':
                    self.nextheader = line[1:].rstrip()
                    return
        except IOError:
            print('No such file', filename)
            raise
            
    def set_filter(self, filter_fn):
        self.filter_fn = filter_fn

    def filter(self, sequence_list):
        """
        Applies the filter_fn on each character of the sequence
        """

        new_list = []
        for sequence in sequence_list:
            new_list.append(''.join(map(self.filter_fn, sequence)))
        return new_list
        
    def has_next(self):
        """
        Returns true if there are still fasta entries
        """
        
        return len(self.nextheader) > 0

    def __next__(self):
        """
        Returns an (id, sequence) tuple, or () if file is finished
        """

        #if global nextheader is empty, return empty
        #otherwise, the header is the nextheader
        try:
            identifier = self.nextheader
            total = []
            while 1:
                line = self.file.readline()
                if line == '' or line[0] == '>':  #EOF, end of entry
                    break
                total.append(line.rstrip())

            if self.filter_fn != None:
                total = self.filter(total)
            sequence = ''.join(total)

            if len(line) > 0:
                self.nextheader = line[1:].rstrip()
            else:
                self.nextheader = ''
                self.close()

            return (identifier, sequence)

        except:
            self.nextheader=''
            self.close()
            return ()

    def next_length(self, return_id = False):
        """
        Returns the length of the next sequence.
        Should not be combined with next, as seqs get discarded
        """

        a = next(self)
        if len(a) == 2:
            if return_id:
                return len(a[1]), a[0]
            else:
                return len(a[1])

    def close(self):
        """
        Close the fasta file
        """
        
        self.file.close()


def write_fasta(filename, id_or_list, seq='', width=60, gzip_compress = False):
    """
    Writes a fasta file with the sequence(s)
    version 1: write_fasta(myfilename, 'seq1_id', 'AAAAA')
    version 2: write_fasta(myfilename, [('seq1_id', 'AAAAA'),
    ('seq2_id', BBBBB)])
    """

    a = fasta_writer(filename, width=width, gzip_compress = gzip_compress)
    a.write(id_or_list, seq)
    a.close()
    

class fasta_writer:
    """
    Rudimentary fasta file writer

    Supports writing out to a gzipped file. If the passed in filename
    does not end with .gz or .z, .gz is appended.
    """

    file = None
    width = 0
    
    def __init__(self, filename, width=60, gzip_compress = False):
        self.width = width
        try:
            if gzip_compress:
                if not filename.endswith('.gz') and not filename.endswith('.z'):
                    filename += '.gz'
                self.file = gzip.open(filename, 'wb')
            else:            
                self.file = open(filename, 'w')
        except IOError:
            print('Can\'t open file.')

    def write(self, id, seq=''):
        """
        Supports an id and a sequence, an (id, seq) tuple, or
        a list of sequence tuples
        """

        if type(id) == type([]):
            list(map(self.writeone, id))
        else:
            self.writeone(id, seq)

    def writeone(self, id, seq=''):
        """
        Internal method.
        """

        if type(id) == type((0,0)):  ###hack!
            seq = id[1]
            id = id[0]

        line_width = self.width
        if self.width == 0:
            line_width = len(seq)
        self.file.write(">" + id + "\n")
        i = 0

        while i < len(seq):
            self.file.write(seq[i:i+line_width] + "\n")
            i+=line_width
            
    def close(self):
        """
        Closes the fasta file.
        """

        self.file.close()


class fasta_struct_writer:
    """
    Writes a combined sequence/structure fasta file

    >id
    sequence
    structure

    """

    file = ''
    width = 0
    
    def __init__(self, filename, width=60):
        self.width = width
        try:
            self.file = open(filename, 'w')
        except IOError:
            print('Can\'t open file.')

    def write(self, id, seq='', struct=''):
        """
        Supports an id and a sequence, an (id, seq, struct) tuple, or
        a list of such tuples
        """

        if type(id) == type([]):
            list(map(self.writeone, id))
        else:
            self.writeone(id, seq, struct)

    def writeone(self, id, seq='', struct=''):
        """
        Internal method
        """

        if type(id) == type((0,0,0)):  ###hack!
            struct = id[2]
            seq = id[1]
            id = id[0]

        line_width = self.width
        if self.width == 0:
            line_width = len(seq)
        self.file.write(">" + id + "\n")
        i = 0
        while i < len(seq):
            self.file.write(seq[i:i+line_width] + "\n")
            i+=line_width

        i = 0
        while i < len(struct):
            self.file.write(struct[i:i+line_width] + "\n")
            i+=line_width


    def close(self):
        """
        Closes the fasta file.
        """
        
        self.file.close()


class struct_fasta_reader:
    """
    Incrementally reads a fasta seq/struct file as output by Vienna.
    """

    file = None
    nextheader=''
    mfe = True
    
    def __init__(self, filename, mfe=True):
        self.mfe = mfe
        try:
            self.file = open(filename, 'r')
            # fast forward to the first entry
            while 1:
                line = self.file.readline()
                if line == '':
                    self.close()
                    return
                elif line[0] == '>':
                    self.nextheader = line[1:].rstrip()
                    return
        except IOError:
            print('No such file', filename)

    def has_next(self):
        """
        True if there are still entries in the fasta file
        """
        
        return len(self.nextheader) > 0

    def __next__(self):
        """Returns an (id, sequence, structure, mfe) tuple,
        or () if file is finished
        """

        #if global nextheader is empty, return empty
        #otherwise, the header is the nextheader
        try:
            identifier = self.nextheader
            seq_total = []
            struct_total = []

            while True:
                line = self.file.readline()
                if line == '' or line[0] in '>(.)':
                    break
                seq_total.append(line.rstrip())

            while True:
                if line == '' or line[0] == '>':
                    break
                struct_total.append(line.rstrip())
                line = self.file.readline()

            sequence = ''.join(seq_total)
            structure = ''.join(struct_total)

            if len(line) > 0:
                self.nextheader = line[1:].rstrip()
            else:
                self.nextheader = ''
                self.close()

            if self.mfe:
                mfe_space = structure.rfind('(')
                mfe = structure[mfe_space+1:-1]
                structure = structure[0:mfe_space]

                return (identifier, sequence, structure, mfe)
            else:
                return (identifier, sequence, structure)

        except:
            self.nextheader=''
            self.close()
            return ()

    def close(self):
        """
        Closes the fasta file.
        """

        self.file.close()


class fastq_reader:
    """
    Lightweight class for incrementally reading fastq files as
    output by Illumina sequencing runs.  Each record is on four
    lines; line wrapping is not supported (i.e., sequences must not
    contain line breaks)

    Supports reading directly from properly named
    gzipped (.gz or .z) or bzip2ed (.bz2) files.
    """
    
    file = None
    nextheader = ''

    def __init__(self, filename):

        try:
            if filename.endswith('.gz') or filename.endswith('.z'):
                self.file = gzip.open(filename, 'rb')
            elif filename.endswith('.bz2'):
                self.file = bz2.BZ2File(filename, 'rb')
            else:
                self.file = open(filename, 'r')
            line = self.file.readline()
            if line == '' or line[0] != '@':
                self.close()
                return
            self.nextheader = line[1:].rstrip()

        except IOError:
            print('No such file', filename)
            raise
            
    def has_next(self):
        """
        True if there are still entries in the fastq file.
        """
        
        return len(self.nextheader) > 0
        

    def __next__(self):
        """
        Returns a 4-tuple consisting of each of the four lines
        minus the @ character for the first line and the + character
        for the third line.  Returns empty tuple if no next
        """

        identifier = self.nextheader
        try:
            seq = self.file.readline().rstrip()
            third = self.file.readline().rstrip()
            quals = self.file.readline().rstrip()

            line = self.file.readline().rstrip()
            if len(line) > 0:
                self.nextheader = line[1:].rstrip()
            else:
                self.nextheader = ''
                self.close()

            return (identifier, seq, third[1:], quals)
        except:
            self.nextheader = ''
            self.close()
            return ()

    def close(self):
        """
        Closes the fastq file.
        """

        self.file.close()
            

class fastq_writer:
    """
    Writes a standard fastq file, four lines per entry, no line wrapping.
    First line prepends a @, third line prepends a +

    Supports writing out to a gzipped file. If the passed in filename
    does not end with .gz or .z, .gz is appended.
    """

    file = None
    
    def __init__(self, filename, gzip_compress = False):
        try:
            if gzip_compress:
                if not filename.endswith('.gz') and not filename.endswith('.z'):
                    filename += '.gz'
                self.file = gzip.open(filename, 'wb')
            else:
                self.file = open(filename, 'w')
        except IOError:
            print('Can\'t open file.')

    def write(self, id, seq='', desc='', quals = ''):
        """
        Supports individual fields, an (id, seq, desc, qual) tuple, or
        a list of such tuples
        """

        if type(id) == type([]):
            list(map(self.writeone, id))
        else:
            self.writeone(id, seq, desc, quals)

    def writeone(self, id, seq='', desc='', quals = ''):
        """
        Internal method
        """
        
        if type(id) == type((0,0,0,0)):  ###hack!
            quals = id[3]
            desc = id[2]
            seq = id[1]
            id = id[0]
        self.file.write('@' + id + '\n')
        self.file.write(seq + '\n')
        self.file.write('+' + desc + '\n')
        self.file.write(quals + '\n')


    def close(self):
        """
        Closes the fastq file.
        """
        
        self.file.close()

        
def read_fastq(filename, abridged = False, as_hash = False):
    """
    Returns the contents of a fastq file in a list of \
    (id, sequence, description, quals) tuples by default or
    (id, sequence) if abridged is True. Empty list returned if there are
    no fastq sequences in the file
    """

    if as_hash:
        return read_fastq_as_hash(filename)
    
    a = fastq_reader(filename)
    seqs = []
    while a.has_next():
        next = next(a)
        if abridged:
            seqs.append((next[0], next[1]))
        else:
            seqs.append(next)
    a.close()
    return seqs


def read_fastq_as_hash(filename):
    a = fastq_reader(filename)
    seqs = {}
    while a.has_next():
        next_seq = next(a)
        seqs[next_seq[0]] = next_seq[1]
    a.close()
    return seqs


def ilmn_fastq_info(ilmn_fastq):
    """
    Queries the first entry in the fastq file and reports
    information, assuming illumina semantics:
    flowcell, lane, pair number, read length
    
    @DCM97JN1:188:C11R8ACXX:5:1101:1208:2164

    This will generally give unreliable results
    if the ID is not a well formed Illumina ID
    """

    a = fastq_reader(ilmn_fastq)
    id, seq, comment, quals = next(a)
    a.close()

    identifier, meta = id.split()
    id_fields = identifier.split(':')
    flowcell = ':'.join(id_fields[:-4])
    lane = int(id_fields[-4])

    if meta:
        pair = int(meta[0])
    else:
        pair = int(id_fields[-1][0])

    readlen = len(seq)

    return flowcell, lane, pair, readlen
    

def count_fastq(filename):
    """
    Returns the number of sequences in the fastq file
    """
    count = 0
    a = fastq_reader(filename)
    while a.has_next():
        next(a)
        count += 1
    return count


def trim_fastq(fastq_in, fastq_out, seq_trim_amt):
    """
    Given a fastq file, trims seq_trim_amt number of nts
    from each sequence and quality string and outputs to the
    fastq_out file, with no changes to the IDs or descriptions.
    Sequences shorter than or equal to seq_trim_amt are converted
    to empty strings.

    To specify final length of the sequence rather than
    number of nts to trim, use truncate_fastq()
    """

    trim_index = -1 * seq_trim_amt

    fr = fastq_reader(fastq_in)
    fw = fastq_writer(fastq_out)

    while fr.has_next():
        id, seq, desc, quals = next(fr)
        seq = seq[:trim_index]
        quals = quals[:trim_index]
        fw.write(id, seq, desc, quals)
    fr.close()
    fw.close()


def trim_fastq_qual(fastq_in, fastq_out, qfail_char = 'B', min_seq_length = 1):
    """
    Given a fastq file, trims each sequence according to the
    quality scores.  For Illumina sequencing, the 'B' character
    is a special quality score indicating unreliable base.  All nts are
    trimmed from the 3' end until the quality string does not end
    with a 'B'.

    The resulting fastq file will contain sequences of different length.
    To impose a minimum length for trimmed sequence, set min_seq_length
    to a value > 1.
    """

    fr = fastq_reader(fastq_in)
    fw = fastq_writer(fastq_out)

    while fr.has_next():
        id, seq, desc, quals = next(fr)
        quals = quals.rstrip(qfail_char)
        lq = len(quals)
        if lq >= min_seq_length:
            seq = seq[:lq]
            fw.write(id, seq, desc, quals)
    fr.close()
    fw.close()


def truncate_fastq(fastq_in, fastq_out, max_seq_len = 50, discard_shorter = False):
    """
    Given a fastq file, trims each sequence and quality string
    to a maximum length of max_seq_len nucleotides.  Sequences shorter
    than max_seq_len are unaffected, unless discard_shorter is True,
    then they are discarded.  Output to the fastq_out
    file with no changes to the IDs or descriptions.

    To specify a fixed number of nts to trim from each sequence
    regardless of length, use trim_fastq()
    """

    fr = fastq_reader(fastq_in)
    fw = fastq_writer(fastq_out)
    while fr.has_next():
        id, seq, desc, quals = next(fr)
        seq = seq[:max_seq_len]
        quals = quals[:max_seq_len]
        if len(seq) == max_seq_len or not discard_shorter:
            fw.write(id, seq, desc, quals)
    fr.close()
    fw.close()


    
def fastq2fasta(fastq_filename, fasta_filename = ''):
    """
    Converts a fastq file to a fasta file as specified.
    Quals and comments are discarded
    if not fasta_filename, then the fastq_filename is used
    replacing the trailing .fq with .fa
    """

    if not fasta_filename:
        if fastq_filename[-3:] == '.fq':
            fasta_filename = fastq_filename[:-3] + '.fa'
        else:
            fasta_filename = fastq_filename + '.fa'            

    fw = fasta_writer(fasta_filename, width=200)
    fq = fastq_reader(fastq_filename)
    while fq.has_next():
        id, seq, comment, quals = next(fq)
        fw.write(id, seq)
    fq.close()
    fw.close()


def fastq2split_fasta(fastq_filename, fasta_filename = '', seqs_per_file = 4000000):
    """
    Converts a fastq file to a series of fasta files each with
    no more than seqs_per_file sequences.

    Quals and comments are discarded.

    if not fasta_filename, then the fastq_filename is used
    replacing the trailing .fq with .fa.px where x is the xth partition
    (starting at 0).
    """

    if not fasta_filename:
        if fastq_filename[-3:] == '.fq':
            fasta_filename = fastq_filename[:-3] + '.fa'
        else:
            fasta_filename = fastq_filename + '.fa'            
    fasta_filename = fasta_filename + '.p%d'

    fnames = []
    curr_fw_id = -1
    fw = None
    
    fr = fastq_reader(fastq_filename)
    counter = 0
    while fr.has_next():
        if counter % seqs_per_file == 0:
            if fw:
                fw.close()
            curr_fw_id += 1
            fw = fasta_writer(fasta_filename % curr_fw_id)
            fnames.append(fasta_filename % curr_fw_id)

        id, seq, comment, quals = next(fr)
        fw.write(id, seq)
        counter += 1
    fr.close()
    try:
        fw.close()
    except:
        pass
    return fnames



def fasta2fastq(fasta_filename, fastq_filename = '', default_qual = 'b'):
    """
    Converts a fasta file to a fastq file with dummy quality scores
    """
    if not fastq_filename:
        if fasta_filename[-3:] == '.fa':
            fastq_filename = fasta_filename[:-3] + '.fq'
        else:
            fastq_filename = fasta_filename + '.fq'

    fq = fastq_writer(fastq_filename)
    fr = fasta_reader(fasta_filename)
    while fr.has_next():
        id, seq = next(fr)
        fq.write(id, seq, '', ''.join([default_qual] * len(seq)))
    fq.close()
    fr.close()

    return fastq_filename


def fastq_ids(fastq_filename, first_word_only = True, strip_slash = True):
    """
    Processes a fastq file and returns a hash of identifiers
    contained within. If first_word_only is True, then only the first
    word of the identifier is returned (everything before the
    first white space character). If strip_slash is True, then all
    characters from a trailing / onward are removed.
    """

    fq_ids = {}
    fq = fastq_reader(fastq_filename)
    while fq.has_next():
        id, seq, comment, qual = next(fq)
        if first_word_only:
            id = id.split()[0]
        if strip_slash:
            a = id.rfind('/')
            if a > -1:
                id = id[:a]
        fq_ids[id] = 1
    fq.close()
    
    return fq_ids


def fastq_strip_comment(fastq_filename, outfile):
    """
    Removes the comment in the input fastq file
    """
    fq = fastq_reader(fastq_filename)
    fw = fastq_writer(outfile)

    while(fq.has_next()):
        id, seq, comment, qual = next(fq)
        fw.write(id, seq, '', qual)
    fq.close()
    fw.close()
    

def fasta_id_compress(fasta_filename, outfile_name = ''):
    """
    Clears the ID lines for the fasta file and replaces
    with a numbered scheme 'S1', 'S2', ... to facilitate
    smaller file sizes
    """

    if not outfile_name:
        outfile_name = fasta_filename[:-3] + '_compress.fa'
    fw = fasta_writer(outfile_name)

    counter = 0    
    f = fasta_reader(fasta_filename)    
    while f.has_next():
        id, seq = next(f)
        fw.write('S%d' % counter, seq)
        counter += 1
    f.close()
    fw.close()
    

#####################
# Sequence processing
#####################

def remove_n_filter_fn(ch):
    """
    Function that detects 'N' or 'n' characters.  If the character
    passed in is 'N', the empty string is returned; otherwise the
    character is returned unmodified.
    """

    if ch == 'n' or ch == 'N':
        return ''
    else:
        return ch


def gc_content(seq):
    """
    Returns the fraction of G or C nts in the sequence.
    """

    count = 0
    for base in seq:
        if base == 'g' or base == 'G' or base == 'C' or base == 'c':
            count += 1
    return 1.0 * count / len(seq)


def create_windows(fasta_file, rc=True, window_size=60, overlap=30, \
                   preserve_shorter = True, outfile='out1.fa'):

    """
    Given a fasta file, outputs a new fasta file with the
    sequence separated into windows; if rc, then the reverse
    complement is also windowed and output
    """

    a = fasta_reader(fasta_file)
    fasta_out = fasta_writer(outfile)
    directions = ['+', '-']

    while a.has_next():
        (id, seq) = next(a);
        #iterate through the sequence
        if len(seq) < window_size:
            if preserve_shorter:
                fasta_out.write(id, seq)
            else:
                continue
        elif len(seq) == window_size:
            new_id = id + ' (' + str(i) + ' - ' + str(end_index-1) + ') ' \
                     + directions[0]
            fasta_out.write(new_id, seq)
            new_id = id + ' (' + str(i) + ' - ' + str(end_index-1) + ') ' \
                     + directions[1]
            fasta_out.write(new_id, reversecomp(seq))
        
        else:
            for strand in range(2):
                i = 0
                done = False

                while i < len(seq) and not done:
                    end_index = i+window_size
                    if end_index > len(seq):
                        i = len(seq) - window_size
                        end_index = len(seq)
                        done = True
                    subseq = seq[i:end_index]
                    new_id = id + ' (' + str(i) + ' - ' + str(end_index-1) \
                             + ') ' + directions[strand]
                    #write fasta
                    fasta_out.write(new_id, subseq)
                    i = i + overlap
                if rc == False:
                    break
                if strand == 0:
                    seq = reversecomp(seq)

    fasta_out.close()
    a.close()


def random_segments(fasta_file, n, length, print_seqs=True):
    """
    Generates uniformly n random length-segments from the
    sequence contained in the fasta file.  Not recommended for large
    sequence files.  Sequences are returned as a list.
    """

    segments = []

    seqs = read_fasta(fasta_file)
    seq = seqs[0][1]
    total_length = len(seq)
    for i in range(n):
        index = random.randint(0, total_length-length)
        segments.append(seq[index:index+length])

    if print_seqs:
        i=0
        for segment in segments:
            print('>segment'+str(i))
            print(segment)
            i=i+1
    else:
        return segments


def random_segment(seq, k, left_bound = 0, right_bound = -1):
    """
    Given a nucleotide sequence, returns a random k-length
    interval between [left_bound, right_bound).

    right_bound = -1 indicates the right boundary is the end of the sequence.
    All other negative numbers are right offsets

    Empty string is returned if right_bound - left_bound < k
    """
    
    if right_bound < 0:
        right_bound = len(seq) + 1 + right_bound
    if right_bound - left_bound < k:
        return ''

    index = random.randint(left_bound, right_bound - k)
    return seq[index:index+k]


def indices_of_motif(seq, motif):
    """
    Wrapper around string find function; returns list of all indices
    where the motif is found, or the empty list if the motif is not found.
    Case sensitive.
    """

    locations = []
    index = 0
    while True:
        index = seq.find(motif, index)
        if index == -1:
            break
        else:
            locations.append(index)
            index += 1
    return locations



################
# Genomic DNA retrieval
################

class genomic_sequence_engine:
    """
    Random access for the genomic sequence in one chromosome at a time.
    """

    db_roots = GENOMIC_SEQ_ROOTS
    rmsk_roots = GENOMIC_RMSK_SEQ_ROOTS

    master_seq = ''

    def __init__(self, chr, species, rmsk = False, preserve_case = False):
        """
        Initialization to specify chromosome number (e.g., chr1) and species
        (currently only 'rn4' and 'mm9' are supported).  If rmsk, then
        the sequence returned is repeatmasked with 'N' characters.  Otherwise
        repeatmasked nucleotides are lower case and normal nts are upper
        case unless preserve_case = False (default), in which case all
        nts are lower case.
        """

        if chr.find('chr') == -1:
            chr = 'chr' + chr
        if rmsk:
            filename = chr + '.fa.masked'
            if species in self.rmsk_roots:
                id, seq = read_fasta(self.rmsk_roots[species] + '/' + filename)[0]
                if not preserve_case:
                    self.master_seq = seq.lower()
                else:
                    self.master_seq = seq
        else:
            filename = chr + '.fa'
            if species in self.db_roots:
                id, seq = read_fasta(self.db_roots[species] + '/' + filename)[0]
                if not preserve_case:
                    self.master_seq = seq.lower()
                else:
                    self.master_seq = seq
            else:
                return None

    def get_sequence(self, start, end, reverse_complement=False):
        """
        For the given start, end range (inclusive) in 1-based coordinates,
        returns the sequence at that position for the chromosome that was
        specified during initialization.  No error checking for coordinate
        ranges beyond the length of the chromosome.

        If reverse_complement, the reverse complement sequence is returned.
        """

        if reverse_complement:
            return rc(self.master_seq[start-1:end])
        else:
            return self.master_seq[start-1:end]


    def unique_hit(self, target_seq):
        """
        Simple sequence search (exact match) in the current chromosome.
        0 = target_seq not present, 1 = present exactly once, 2 = present
        >1 locations.
        """
        
        return unique_sequence(target_seq, self.master_seq)


class GenomicSequenceEngineFull:
    """
    Random access for genomic sequence across all chromosomes simultaneously.
    Designed for small genomes that are reported in a single file
    """

    fa_files = GENOMIC_SEQ_FILES
    seqs = {}

    def __init__(self, species, preserve_case = False, repeat_mask = False, fa_file = None):
        if fa_file:
            chr_seqs = read_fasta(fa_file)
        else:
            chr_seqs = read_fasta(self.fa_files[species])
        if repeat_mask:
            for id, seq in chr_seqs:
                for nt in ['a', 'c', 'g', 't']:
                    seq = seq.replace(nt, 'N')
                self.seqs[id] = seq
        
        elif not preserve_case:
            for id, seq in chr_seqs:
                self.seqs[id] = seq.lower()
        else:
            for id, seq in chr_seqs:
                self.seqs[id] = seq

                
    def get_sequence(self, chr, start, end, reverse_complement=False):
        """
        For the given start, end range (inclusive) in 1-based coordinates,
        returns the sequence at that position for the chromosome
        No error checking for coordinate
        ranges beyond the length of the chromosome.

        Chromosome must be specified in UCSC style with the 'chr' prefix.
        
        If reverse_complement, the reverse complement sequence is returned.
        """

        if reverse_complement:
            return rc(self.seqs[chr][start-1:end])
        else:
            return self.seqs[chr][start-1:end]

    def get_sequence_from_segments(self, chr, segment_list, reverse_complement=False, padding = 0, three_padding = 0):
        """
        Rather than specifying a single start and end, a list of (start,end)
        tuples is provided as an ordered segment_list. Sequences for each
        segment are concatenated together. If reverse_complement,
        the final concatenated product is reverse complemented.

        4Jun12 -- fixed bug: segment_list is no longer mutated, copy
        is made first.
        """

        segment_list = segment_list[:]
        
        if padding:
            left_coord = segment_list[0][0]
            if left_coord > 1:
                segment_list.insert(0, (max(1, left_coord - padding), left_coord - 1))            
            right_coord = segment_list[-1][1]
            segment_list.append((right_coord + 1, right_coord + padding))

        if three_padding:
            if reverse_complement:
                left_coord = segment_list[0][0]
                if left_coord > 1:
                    segment_list.insert(0, (max(1, left_coord - three_padding), left_coord - 1))
            else:
                right_coord = segment_list[-1][1]
                segment_list.append((right_coord + 1, right_coord + three_padding))
                            
        seq = []
        for start, end in segment_list:
            seq.append(self.seqs[chr][start-1:end])

        seq = ''.join(seq)
        if reverse_complement:
            return rc(seq)
        else:
            return seq
            

def genomic_sequence(chr, start, end, species):
    """
    Returns the sequence specified by the coordinate range
    passed in.  Start and end are with respect to 1-indexing
    and are inclusive.  Wrapper for a single query to the
    genomic_sequence_engine class; for multiple queries, use
    the engine explicitly.

    Currently supported species are rn4 and mm9
    """

    #hack
    if species == 'danRer7':
        gse = GenomicSequenceEngineFull(species)
        return gse.get_sequence(chr, start, end)
    else:  
        gse = genomic_sequence_engine(chr, species)
        return gse.get_sequence(start,end)
    

def unique_sequence(target_seq, source_seq):
    """
    Simple sequence search (exact match) in the current chromosome.
    0 = target_seq not present, 1 = present exactly once, 2 = present
    >1 locations.
    """

    left = source_seq.find(target_seq)
    if left == -1:
        return 0
    right = source_seq.rfind(target_seq)
    if left == right:
        return 1
    else:
        return 2



################
# Sequence simulation
################

def random_dna_sequence(length, joined=True):
    """
    Given a specified length, returns a uniformly random DNA sequence.
    If joined == False, the sequence is returned as a list, otherwise
    a string is returned.
    """

    seq = []
    for i in range(length):
        seq.append(random_base(True))
    if joined:
        return ''.join(seq)
    else:
        return seq


def random_rna_sequence(length, uniform=True, hyperuniform=False,\
                        base_proportion = []):
    """
    Given a specified length, returns by default a uniformly random
    RNA sequence as a string.  If uniform = False, then the base proportions
    will deviate from 0.25.  Keeping the default will cause the base
    proportions to be selected from a normal distribution with mean = 1
    std = 0.25, then normalized to sum to 1.  Setting hyperuniform = True
    causes the proportions to be selected from a uniform distribution
    and normalized to 1.  Alternatively, specify the base_proportion
    as a list of proportions for A, C, and G, and T (which will be normalized
    if they don't sum to 1).
    """    

    seq = []
    if uniform:
        for i in range(length):
            seq.append(random.choice(rna_bases))
    else:
        if hyperuniform:
            p1 = random.random()
            p2 = random.random()
            p3 = random.random()
            p4 = random.random()
        else:
            if base_proportion == []:
                p1 = max(0.0001, random.gauss(1, 0.25))
                p2 = max(0.0001, random.gauss(1, 0.25))            
                p3 = max(0.0001, random.gauss(1, 0.25))
                p4 = max(0.0001, random.gauss(1, 0.25))
                total = p1+p2+p3+p4
                p1 = p1/total
                p2 = p2/total
                p3 = p3/total                
                p4 = p4/total
            else:
                p1 = base_proportion[0]
                p2 = base_proportion[1]
                p3 = base_proportion[2]
                p4 = base_proportion[3]

        for i in range(length):
            seq.append(random_base_p([p1,p2,p3,p4]))

    return ''.join(seq)


def mutate_sequence(seq, mu=0.001, upper = False, rna=True):
    """
    Mutates the input sequence with mutation rate mu.  Returns
    the mutated sequence
    """

    new_seq = []
    for base in seq:
        if random.random() < mu:
            new_seq.append(random_base(upper=upper, rna=rna))
        else:
            new_seq.append(base)
    return ''.join(new_seq)


class random_ma_creator:
    """
    Creates a rudimentary random multiple alignment and outputs
    to a specified file.  Note that the output is not actually
    aligned (i.e., no gap characters are inserted), but rather
    consists of a fasta file with alignable sequences.
    """

    indel_rate = .1
    deletion_rate = .125     # P(deletion|indel)
    insertion_rate = .125
    gap_extend_rate = .3
    substitution_rate = .75

    def __init__(self, length = None, num_seqs = None, seed='', outfile=''):
        """
        Initializes the class.  Keep all the input parameters as default
        to use this class for multiple MA files (by calling generate()
        several times).  Alternatively, to create one MA, specify the
        parameters, which will be passed to a single call to generate()
        """

        if length != None:
            self.generate(length, num_seqs, seed, outfile)

    
    def generate(self, length, num_seqs, seed='', outfile=''):
        """
        Generate the ma of starting length with num_seqs sequences
        and write to outfile.  Seed is an optional starting sequence
        to guide the MA sequences; otherwise a random dna sequence
        is used.  All entries in the file will be random mutations
        of the original sequence.
        """

        ma = []
        if seed == '':
            seed = random_dna_sequence(length, joined=False)
        else:
            seed = list(seed)
        ma.append(''.join(seed))
        for i in range(num_seqs-1):
            ma.append(''.join(map(self.mutate, seed)))
        
        if outfile == '':
            return ma
        else:
            fw = fasta_writer(outfile)
            for i in range(len(ma)):
                fw.write('seq_' + str(i), ma[i])
            fw.close()


    def mutate(self, base):
        """
        Mutates the base with some probability.  Called during generate()
        """

        if random.random() < self.indel_rate:
            indel_choice = random.random()
            if indel_choice < self.deletion_rate:
                return ''
            elif indel_choice < self.deletion_rate + self.insertion_rate:
                return self.insert() + base
            else:
                return random_base(upper=True, exclude=base)
        else:
            return base

            
    def insert(self):
        """
        Nt insertion scheme with separate probabilities for insertion
        site and length of insertion. Called during generate()
        """

        seq = []
        extend_chance = 0
        while extend_chance < self.gap_extend_rate:
            seq.append(random_base(True))            
            extend_chance = random.random()
        return ''.join(seq)



################
# Genomic coordinate processing
################

def gb_coord(chr, start, stop, prepend_chr = True, intelligent_parse = False):
    """
    Returns a string of the form 'chrXX:start-stop'

    If intelligent_parse, the chromosome name is queried, and chr is
    automatically prepended when appropriate
    """

    coord = ['chr', chr, ':', str(start), '-', str(stop)]
    
    if intelligent_parse:
        if chr[0].isdigit() or chr in ['X', 'Y', 'M']:
            return ''.join(coord)
        else:
            return ''.join(coord[1:])
    else:
        if prepend_chr and not chr.startswith('chr'):
            return ''.join(coord)
        else:
            return ''.join(coord[1:])

#    if not prepend_chr or chr.find('chr') == 0:
#        return chr + ':' + str(start) + '-' + str(stop)
#    else:
#        return 'chr' + chr + ':' + str(start) + '-' + str(stop)
        

def parse_gb_coord_comma(gb_string, dir=False):
    """
    Accomodates commas in the numbers.  Separate method due
    to efficiency
    """
    return parse_gb_coord(gb_string.replace(',', ''), dir)


def parse_gb_coord_fast(gb_string):
    """
    Implements a parse_gb_coord() with dir = False
    and include_chr = True
    """

    chr, indices = gb_string.split(':')
    start, end = indices.split('-')
    return chr, int(start), int(end)


def parse_gb_coord(gb_string, dir = False, include_chr = False):
    """
    Parses a genome browser coordinate and returns the chromosome
    (without 'chr' prefix by default), start, and stop indices.

    If dir=True, then the last character is assumed to be a single-symbol
    direction (+/-/F/R/etc) and is returned as well

    If include_chr = True then chromosome names are returned as
    'chrX' rather than just 'X'

    If the coordinate doesn't parse a ValueError is thrown
    """

    chr, indices = gb_string.split(':')
    if not include_chr:
        chr = chr[3:]

    start, end = indices.split('-')

    if dir:
        return chr, int(start), int(end[:-1]), end[-1]
    else:
        return chr, int(start), int(end)


def parse_gb_coord_list(gb_string, sep=';', include_chr = False):
    """
    Parses lists of genome browser coordinates, each separated
    by the sep character -- e.g., 'chr1:12345-12445;chr1:45567-45600'.

    This is the most flexible of the gb_coord functions.
    """

    fields = gb_string.split(sep)
    coords = [parse_gb_coord(x, include_chr = include_chr) for x in fields]
    return coords


def parse_gb_coord_list_fast(gb_string, sep=';'):
    """
    Faster implementation of parse_gb_coord_list using calls to
    parse_gb_coord_fast for each coord in the list.
    """

    fields = gb_string.split(sep)
    coords = list(map(parse_gb_coord_fast, fields))
    return coords


def coord_span(gb_string, zero_indexing = False):
    """
    Calculates the number of nts spanned by the input coord string,
    which is processed using parse_gb_coord_list_fast.
    """

    if zero_indexing:
        offset = 0
    else:
        offset = 1
        
    coords = parse_gb_coord_list_fast(gb_string)
    span = sum([x[2]-x[1]+offset for x in coords])
    return span


def coord_boundaries(gb_string):
    """
    Returns the start and end coordinate of the coord list,
    which is parsed using parse_gb_coord_list_fast
    """

    coords = parse_gb_coord_list_fast(gb_string)
    start = coords[0][1]
    end = coords[-1][2]
    return start, end


def spanning_coord(gb_string):
    """
    Returns a gb coord that spans the entire coord list.
    No check in made that the same chromosome applies
    to all members
    """

    positions = []
    coords = parse_gb_coord_list(gb_string, include_chr = True)
    for chr, start, end in coords:
        positions.append(start)
        positions.append(end)
    return gb_coord(chr, min(positions), max(positions), prepend_chr = False)

    

################
# Sequence alignment
################

def edit_distance(s1, s2, match=1, gap=-1, mismatch=-1, seq_len_equal=False):
    """
    Given two sequences, runs NW and returns the global score.
    Uses a backpointer-less NW implementation, thus the trace is not
    returned
    """

    len_s1 = len(s1)
    len_s2 = len(s2)

    #initialize the dp matrix with zeros
    dp_matrix = []
    row = []
    for j in range(len_s2+1):
        row.append(0)
    for i in range(len_s1+1):
        dp_matrix.append(row[:])

    #initialize the top row/top column with appropriate gap penalty
    for i in range(len_s2+1):
        dp_matrix[0][i] = gap * i
    for i in range(len_s1+1):
        dp_matrix[i][0] = gap * i

    #go for it, row by row
    for i in range(len_s1+1):
        if i == 0:
            continue
        row = dp_matrix[i]
        for j in range(len_s2+1):
            if j == 0:
                continue
            current = 0
            if s1[i-1] == s2[j-1]:
                current = match
            else:
                current = mismatch
            row[j] = max(dp_matrix[i-1][j]+gap, dp_matrix[i][j-1]+gap, dp_matrix[i-1][j-1]+current)

    return dp_matrix[len_s1][len_s2]


def edit_distance_optim(s1, s2, seq_len, match=1, gap=-1, mismatch=-1):
    """
    Fast edit distance calculation, assuming equal length sequences.
    """

    matrix_size = seq_len + 1

    col_prev = [gap*i for i in range(matrix_size)]
    col_curr = [0] * matrix_size
    
    for i in range(1, matrix_size):  #cols
        col_curr[0] = gap * i
        for j in range(1, matrix_size):
            if s1[i-1] == s2[j-1]:
                current = match + col_prev[j-1]
            else:
                current = mismatch + col_prev[j-1]
            gap1 = col_prev[j]+gap
            gap2 = col_curr[j-1]+gap
            if current > gap1:
                if current > gap2:
                    col_curr[j] = current
                else:
                    col_curr[j] = gap2
            else:
                if gap1 > gap2:
                    col_curr[j] = gap1
                else:
                    col_curr[j] = gap2

        temp_col = col_prev
        col_prev = col_curr
        col_curr = temp_col

    return col_prev[seq_len]


def hamming_distance(s1, s2, strict = False):
    """
    Compare s1 and s2 and return the no. of characters in which they differ.
    If strict, assumes len(s1) = len(s2) and returns None otherwise.
    Otherwise, truncates the longer
    """

    if strict and len(s1) != len(s2):
        return None
    r = 0
    for i in range(min(len(s1), len(s2))):
        if s1[i] != s2[i]: r+=1
    return r


def hamming_distance_indices(s1, s2):
    """
    Compare s1 and s2 and return a list of indices where the
    sequences differ.  Assumes lengths are the same.  If not, None
    is returned.
    """

    if len(s1) != len(s2):
        return None

    diff_list = []

    for i in range(len(s1)):
        if s1[i] != s2[i]: diff_list.append(i)
    return diff_list


def smith_water(s1, s2, match=1, gap=-1, mismatch=-1):
    """
    Returns the maximum score of a local alignment between the two
    sequences.  Alignment is not returned.
    """

    max_score = 0

    #initialize the dp matrix with zeros
    dp_matrix = []
    row = []
    for j in range(len(s2)+1):
        row.append(0)
    for i in range(len(s1)+1):
        dp_matrix.append(row[:])

    #go for it, row by row
    for i in range(1, len(s1)+1):
        row = dp_matrix[i]
        for j in range(1, len(s2)+1):
            current = 0
            if s1[i-1] == s2[j-1]:
                current = match
            else:
                current = mismatch
            row[j] = max(dp_matrix[i-1][j]+gap, dp_matrix[i][j-1]+gap, dp_matrix[i-1][j-1]+current, 0)
            if row[j] > max_score:
                max_score = row[j]

    return max_score


def sw_memoize(s1, s2, min_pathlength=2, count_diagonal=False):
    """
    Specialized function that sums the Smith-Waterman matrix scores
    over the entire alignment, as a proxy for degree of global repetitive
    similarity.
    
    Records the length of the path taken at a particular cell, and then
    sums only cells at the end of paths greater than min_pathlength.  Set
    min_pathlength = 0 to sum the entire SW matrix.

    count_diagonal = False ignores scores occurring on the diagonal; this
    is useful if s1 = s2 and you don't want to count trivial similarity.
    """

    gap = -6
    mismatch = -2
    match = 1

    match_p = 2
    match_d = 1

    #initialize the dp matrix with zeros
    dp_matrix = []
    row = []
    for j in range(len(s2)+1):
        row.append((0, 0))
    for i in range(len(s1)+1):
        dp_matrix.append(row[:])


    #go for it, row by row
    for i in range(len(s1)+1):
        if i == 0:
            continue
        row = dp_matrix[i]
        for j in range(len(s2)+1):
            if j == 0:
                continue
            current = 0
            if s1[i-1] == s2[j-1]:
                if s1[i-1] == '.':
                    current = match_d
                else:
                    current = match_p
            else:
                current = mismatch

            gap_i = dp_matrix[i-1][j][0]+gap
            gap_j = dp_matrix[i][j-1][0]+gap
            no_gap = dp_matrix[i-1][j-1][0]+current

            if gap_i > gap_j and gap_i > no_gap and gap_i > 0:
                row[j] = (gap_i, dp_matrix[i-1][j][1] + 1)
#                dp_matrix[i-1][j] = (dp_matrix[i-1][j][0], True)
            elif gap_j > no_gap and gap_j > 0:
                row[j] = (gap_j, dp_matrix[i][j-1][1] + 1)
#                dp_matrix[i][j-1] = (dp_matrix[i][j-1][0], True)
            elif no_gap > 0:
                row[j] = (no_gap, dp_matrix[i-1][j-1][1] + 1)
#                dp_matrix[i-1][j-1] = (dp_matrix[i-1][j-1][0], True)
            else:
                row[j] = (0, 0)

#            row[j] = max(dp_matrix[i-1][j]+gap, dp_matrix[i][j-1]+gap, dp_matrix[i-1][j-1]+current, 0)

#    for row in dp_matrix:
#        print row

    sum = 0
    for i in range(len(dp_matrix)):
        for j in range(len(dp_matrix[i])):
            if i == j:
                if not count_diagonal:
                    continue
            if dp_matrix[i][j][1] >= min_pathlength:
                sum += dp_matrix[i][j][0]                

#    for row in dp_matrix:
#        for elem in row:
#            sum += elem[0]
    return sum


def nw_endfree(s1, s2, match=5, gap_open=-10, gap_extend=-0.5, mismatch=-2):
    """
    Needleman wunsch global alignment with free end gaps.  Alignment is
    returned.
    """

    nrows = len(s1) + 1
    ncols = len(s2) + 1

    # score, matrix, row, col
    match_matrix = [[(0,-1,-1,-1)]*ncols for i in range(nrows)]
    gap_matrix = [[(0,-1,-1,-1)]*ncols for i in range(nrows)]

    score_matrix = [match_matrix, gap_matrix]

    for j in range(1, ncols):
        for i in range(1, nrows):
            if s1[i-1] == s2[j-1]:
                match_value = match
            else:
                match_value = mismatch

            #UPDATE MATCH MATRIX
            if match_matrix[i-1][j-1][0] > gap_matrix[i-1][j-1][0]:
                match_matrix[i][j] = (match_matrix[i-1][j-1][0] + match_value, 0, i-1, j-1)
            else:
                match_matrix[i][j] = (gap_matrix[i-1][j-1][0] + match_value, 1, i-1, j-1)


            #UPDATE GAP MATRIX
            max_score = (-32768, -1, -1, -1)

            if match_matrix[i][j-1][0] + gap_open > max_score[0]:
                max_score = (match_matrix[i][j-1][0] + gap_open, 0, i, j-1)
            if match_matrix[i-1][j][0] + gap_open > max_score[0]:
                max_score = (match_matrix[i-1][j][0] + gap_open, 0, i-1, j)
            if gap_matrix[i][j-1][0] + gap_extend > max_score[0]:
                max_score = (gap_matrix[i][j-1][0] + gap_extend, 1, i, j-1)
            if gap_matrix[i-1][j][0] + gap_extend > max_score[0]:
                max_score = (gap_matrix[i-1][j][0] + gap_extend, 1, i-1, j)

            gap_matrix[i][j] = max_score


    # MAX COORDINATE IN LAST ROW/LAST COL
    max_score = -32768
    max_coord = (0, -1, -1, -1)

    for i in range(nrows):
        if match_matrix[i][-1][0] > max_score:
            max_score = match_matrix[i][-1][0]
            max_coord = (0, 0, i, ncols-1)
    for j in range(ncols):
        if match_matrix[-1][j][0] > max_score:
            max_score = match_matrix[-1][j][0]
            max_coord = (0, 0, nrows-1, j)


    # TRACEBACK FROM MAX COORDINATE
    curr_coord = max_coord
    prev_coord = None

    s1_align = []
    s2_align = []

    while curr_coord != (0, -1, -1, -1):
        if prev_coord == None:
            s1_align.append(s1[curr_coord[2] - 1])
            s2_align.append(s2[curr_coord[3] - 1])
        elif curr_coord[2] == 0 or curr_coord[3] == 0:
            break
        elif prev_coord[2] == curr_coord[2]:
            s1_align.append('-')
            s2_align.append(s2[curr_coord[3] - 1])
        elif prev_coord[3] == curr_coord[3]:
            s1_align.append(s1[curr_coord[2] - 1])
            s2_align.append('-')
        else:
            s1_align.append(s1[curr_coord[2] - 1])            
            s2_align.append(s2[curr_coord[3] - 1])                            

        prev_coord = curr_coord
        curr_coord = score_matrix[prev_coord[1]][prev_coord[2]][prev_coord[3]]

    s1_align.reverse()
    s2_align.reverse()

#    print prev_coord[2], max_coord[2], ''.join(s1_align)
#    print prev_coord[3], max_coord[3], ''.join(s2_align)

    return prev_coord[2], max_coord[2], ''.join(s1_align), prev_coord[3], max_coord[3], ''.join(s2_align)



def nw_gapfree_5prime_anchored(s1, s2, match=5, mismatch=-2, s1_anchor = False, s2_anchor = True):
    """
    Needleman wunsch global alignment.
    Free end gaps on the 3' end, but one of the 5' ends must align,
    by default the second sequence.  If both s1_anchor and s2_anchor
    are True, then the 5' ends of both sequences must align.

    Interior gaps are not allowed.

    Return is a 5-tuple: (start1, end1, start2, end2, score) where
    start and end coords refer to aligned portions of s1 and s2
    """

    nrows = len(s1) + 1
    ncols = len(s2) + 1

    max_row = -1
    max_col = -1
    max_score = 0

    # score, matrix, row, col
    score_matrix = [[0]*ncols for i in range(nrows)]

    for j in range(1, ncols):
        for i in range(1, nrows):
            if s1[i-1] == s2[j-1]:
                match_value = match
            else:
                match_value = mismatch

            #UPDATE SCORE MATRIX
            score = score_matrix[i-1][j-1] + match_value
            score_matrix[i][j] = score
            if score > max_score:
                triangle = i - j
                if (triangle == 0 and s1_anchor and s2_anchor) or \
                   (triangle > 0 and s2_anchor and not s1_anchor) or \
                   (triangle < 0 and s1_anchor and not s2_anchor):

                    max_score = score
                    max_row = i
                    max_col = j

    # TRACE COORD BACK TO CORRESPONDING FIRST ROW/FIRST COL
    if max_row <= max_col:
        start_row = 0
        start_col = max_col - max_row
    else:
        start_row = max_row - max_col
        start_col = 0
    return start_row, max_row - 1, start_col, max_col - 1, max_score

def strand_reverse(strand):
    """
    + maps to -, - maps to +, all others map to '.'
    """

    opposite = {'+':'-', '-': '+', '.': '.'}
    
    try:
        return opposite[strand]
    except:
        return '.'
