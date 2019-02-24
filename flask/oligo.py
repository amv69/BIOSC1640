import math
import bioutil

iupac_wc = {'A': ('A',),
            'C': ('C',),
            'G': ('G',),
            'T': ('T',),
    'R': ('A', 'G'),
    'Y': ('C', 'T'),
    'S': ('C', 'G'),
    'W': ('A', 'T'),
    'K': ('G', 'T'),
    'M': ('A', 'C'),
    'B': ('C', 'G', 'T'),
    'D': ('A', 'G', 'T'),
    'H': ('A', 'C', 'T'),
    'V': ('A', 'C', 'G'),
    'N': ('A', 'C', 'G', 'T')}

iupac_wc_rev = {('A', 'G'): 'R',
               ('C', 'T'): 'Y',
               ('C', 'G'): 'S',
               ('A', 'T'): 'W',
               ('G', 'T'): 'K',
               ('A', 'C'): 'M',
               ('C', 'G', 'T'): 'B',
               ('A', 'G', 'T'): 'D',
               ('A', 'C', 'T'): 'H',
               ('A', 'C', 'G'): 'V',
               ('A', 'C', 'G', 'T'): 'N'}

#Source: Sugimoto 1995 (RNA/DNA)
# listed with respect to the RNA target
RNADNAdeltaH = {'AA': 7.8,
          'AC': 5.9,
          'AG': 9.1,
          'AT': 8.3,
          'CA': 9.0,
          'CC': 9.3,
          'CG': 16.3,
          'CT': 7.0,
          'GA': 5.5,
          'GC': 8.0,
          'GG': 12.8,
          'GT': 7.8,
          'TA': 7.8,
          'TC': 8.6,
          'TG': 10.4,
          'TT': 11.5
          }
           
RNADNAdeltaS = {'AA': 0.0219,
          'AC': 0.0123,
          'AG': 0.0235,
          'AT': 0.0239,
          'CA': 0.0261,
          'CC': 0.0232,
          'CG': 0.0471,
          'CT': 0.0197,
          'GA': 0.0135,
          'GC': 0.0171,
          'GG': 0.0319,
          'GT': 0.0216,
          'TA': 0.0232,
          'TC': 0.0229,
          'TG': 0.0284,
          'TT': 0.0364
          }
          

#Source: Sugimoto 1996 (DNA/DNA)
DNAdeltaH = {'AA': 8.0,
          'TT': 8.0,
          'AT': 5.6,
          'TA': 6.6,
          'CA': 8.2,
          'TG': 8.2,
          'GT': 9.4,
          'AC': 9.4,
          'CT': 6.6,
          'AG': 6.6,
          'GA': 8.8,
          'TC': 8.8,
          'CG': 11.8,
          'GC': 10.5,
          'GG': 10.9,
          'CC': 10.9
           }

#kcal/K*mol
DNAdeltaS = {'AA': 0.0219,
          'TT': 0.0219,
          'AT': 0.0152,
          'TA': 0.0184,
          'CA': 0.0210,
          'TG': 0.0210,
          'GT': 0.0255,
          'AC': 0.0255,
          'CT': 0.0164,
          'AG': 0.0164,
          'GA': 0.0235,
          'TC': 0.0235,
          'CG': 0.0290,
          'GC': 0.0264,
          'GG': 0.0284,
          'CC': 0.0284
          }

#Source: SantaLucia 1998
RNAdeltaH = {'AA': 6.82,
          'AC': 11.4,
          'AG': 10.48,
          'AT': 9.38,
          'CA': 10.44,
          'CC': 13.39,
          'CG': 10.64,
          'CT': 10.48,
          'GA': 12.44,
          'GC': 14.88,
          'GG': 13.39,
          'GT': 11.4,
          'TA': 7.69,
          'TC': 12.44,
          'TG': 10.44,
          'TT': 6.82,
          }
           
RNAdeltaS = {'AA': 0.019,
          'AC': 0.0295,
          'AG': 0.0271,
          'AT': 0.0267,
          'CA': 0.0269,
          'CC': 0.0327,
          'CG': 0.0267,
          'CT': 0.0271,
          'GA': 0.0325,
          'GC': 0.0369,
          'GG': 0.0327,
          'GT': 0.0295,
          'TA': 0.0205,
          'TC': 0.0325,
          'TG': 0.0269,
          'TT': 0.0190
          }


RNA_DNA_HELIX_INIT = -3.1 #deltaG
DNA_HELIX_INIT = -3.4
RNA_HELIX_INIT = -4.09
R = 0.001987 #kcal/K*mol #gas constant


def expand_wildcards(seq):
    """
    Makes all combination of sequences with wildcards replaced
    by unambiguous bases. Returns a list
    """
    bases = list(seq.upper())
    seqs = [bases]
    for i, base in enumerate(bases):
        expansion = iupac_wc.get(base, ('X',))
#        if not expansion:
#            continue
        new_seqs = []
        for temp_seq in seqs:
            for x in expansion:
                new_seq = temp_seq[:]
                new_seq[i] = x
                new_seqs.append(new_seq)
        seqs = new_seqs
    seqs = [ ''.join(x) for x in seqs ]
        
    return seqs


def wildcard_stats(seq):
    """
    Returns number of wildcard positions and number of possible
    expansions
    """

    wc_count = 0
    expansions = 1
    for base in seq:
        if base not in ['A', 'C', 'G', 'T', 'X']:
            wc_count += 1
            expansions *= len(iupac_wc.get(base, ()))
    return seq.count('X'), wc_count, expansions
    

def calculate_tm(dh, ds, helix_init, oligo_concen, na_concen):
    return round( (dh + helix_init)/ (ds + R * math.log((1.0/oligo_concen), math.e)) - 273.15 + 16.6*math.log(na_concen, 10) , 1)


def do_tm(seq, C = 0.00000005, Na_concen = 0.2, molec = 'dnarna'):
    """
    Formula given by oligocalc
    http://biotools.nubic.northwestern.edu/OligoCalc.html
    (they use a K to C conversion of -272.9 for some reason,
    so the value will not exactly match)
    (Also, they use the DNA helix initiation energy for RNA)
    
    Na_concen in mol/L
    C in mol/L concen of oligo
    """
    if molec == 'dnarna':
        dH = RNADNAdeltaH
        dS = RNADNAdeltaS
        helix_init = RNA_DNA_HELIX_INIT
    elif molec == 'rna':
        dH = RNAdeltaH
        dS = RNAdeltaS
        helix_init = RNA_HELIX_INIT
    else:
        dH = DNAdeltaH
        dS = DNAdeltaS
        helix_init = DNA_HELIX_INIT
        
    sum_dh = 0
    sum_ds = 0

    seq = seq.upper()
    for i in range(len(seq)-1):
        dnt = seq[i:i+2]
        sum_dh += dH[dnt]
        sum_ds += dS[dnt]

#    print(sum_dh, sum_ds)
    return calculate_tm(sum_dh, sum_ds, helix_init, C, Na_concen)


def iterative_tms(seq, start_at = 0, min_len = 40, max_len = 50, C = 0.00000005, Na_concen = 0.2, molec = 'dnarna', min_tm = 70, max_tm = 80, strict = True):
    """
    Calculates Tm for all prefixes, retaining only the
    ones with acceptable length.
    """

    if molec == 'dnarna':
        dH = RNADNAdeltaH
        dS = RNADNAdeltaS
        helix_init = RNA_DNA_HELIX_INIT
    elif molec == 'rna':
        dH = RNAdeltaH
        dS = RNAdeltaS
        helix_init = RNA_HELIX_INIT
    else:
        dH = DNAdeltaH
        dS = DNAdeltaS
        helix_init = DNA_HELIX_INIT

    sub_seq = seq[start_at:start_at+max_len]
        
    sum_dh = 0
    sum_ds = 0
    tms = []
    acceptable_lens = []
    
    if len(sub_seq) < 8:
        return tms, acceptable_lens
    
    for i in range(1, len(sub_seq)):
        dnt = sub_seq[i-1:i+1]
        sum_dh += dH[dnt]
        sum_ds += dS[dnt]
        tm = calculate_tm(sum_dh, sum_ds, helix_init, C, Na_concen)
        if i+1 >= min_len:
            if min_tm <= tm <= max_tm:
                acceptable_lens.append(i+1)
                tms.append((start_at, start_at + i+1, tm))
            elif not strict:
                tms.append((start_at, start_at + i+1, tm))
    return tms, acceptable_lens


def iterative_tms_wc(seq, start_at = 0, min_len = 40, max_len = 50, C = 0.00000005, Na_concen = 0.2, molec = 'dnarna', min_tm = 70, max_tm = 80, strict = True, max_num_wildcards = 2):
    """
    Calculates Tm for all prefixes, retaining only the
    ones with acceptable length. Support for wildcard bases,
    which induces a Tm range for each of the possible bases
    Returns a list of Tms: (start_pos, end_pos, Tm_low, Tm_high)
    """

    if molec == 'dnarna':
        dH = RNADNAdeltaH
        dS = RNADNAdeltaS
        helix_init = RNA_DNA_HELIX_INIT
    elif molec == 'rna':
        dH = RNAdeltaH
        dS = RNAdeltaS
        helix_init = RNA_HELIX_INIT
    else:
        dH = DNAdeltaH
        dS = DNAdeltaS
        helix_init = DNA_HELIX_INIT

    sub_seq = seq[start_at:start_at+max_len]
        
    tms = []
    oligo_lens = []
    
    if len(sub_seq) < 8:
        return [], []

    n_wildcards = 0
    prev_params = [(0, 0, '')]
    for i in range(1, len(sub_seq)):
        dnt = sub_seq[i-1:i+1]
        if 'X' in dnt:
            break
        dnts = expand_wildcards(dnt)
        n_wildcards += len(dnts) - 1
        tm_range = []
        curr_params = []
        for dnt in dnts:
            for sum_dh, sum_ds, base in prev_params:
                if not base or base == dnt[0]:
                    sum_dh += dH[dnt]
                    sum_ds += dS[dnt]
                    curr_params.append((sum_dh, sum_ds, dnt[1]))
                    tm_range.append(calculate_tm(sum_dh, sum_ds, helix_init, C, Na_concen))
        prev_params = curr_params
            
        if i+1 >= min_len and math.ceil(n_wildcards/2) <= max_num_wildcards:
            tm_low = min(tm_range)
            tm_high = max(tm_range)
            if (tm_low >= min_tm and tm_high <= max_tm) or not strict:
                tms.append((start_at, start_at + i+1, tm_low, tm_high))
                oligo_lens.append(i+1)
    return tms, oligo_lens


def generate_tms(seq, min_tm = 70, max_tm = 80, min_len = 40, max_len = 50, C = 0.00000005, Na_concen = 0.2, molec = 'dnarna', strict = True):
    """
    Calculates all possible tms iteratively starting at each position.
    Also calculates the median length of oligo within the tm range

    If strict, only oligos with Tm within range are returned    
    """

    tms = []
    oligo_lens = []

    for i in range(len(seq)):
        tm_list, acceptable_lens = iterative_tms_wc(seq, i, min_len, max_len, C, Na_concen, molec, min_tm, max_tm, strict = strict)
        tms.append(tm_list)
        oligo_lens += acceptable_lens

    oligo_lens.sort()

    try:
        return tms, oligo_lens[len(oligo_lens)//2]
    except:
        return tms, min_len

    
def tile_oligos(seq, min_tm = 70, max_tm = 80, min_len = 40, max_len = 50, C = 0.00000005, Na_concen = 0.2, molec = 'dnarna'):
    """
    TODO: store intermediate solutions, propose overlap
    """
    
    tms, median_oligo_len = generate_tms(seq, min_tm, max_tm, min_len, max_len, C, Na_concen, molec)
        
    curr_tile_paths = [ [x] for x in tms[0] ]

    while True:
        new_tp = []
        for tp in curr_tile_paths:
            if tp[-1][1] == len(seq): #Success
                return tp
            for cand in tms[tp[-1][1]]:
                new_tp.append(tp + [cand])
        curr_tile_paths = new_tp
        if not curr_tile_paths:
            break

    return []


def tm_dev(tm, min_tm, max_tm):
    if min_tm <= tm <= max_tm:
        return 0
    return min(abs(tm-min_tm), abs(tm-max_tm))

def tm_range_dev(tm_low, tm_hi, min_tm, max_tm):
    return max(tm_dev(tm_low, min_tm, max_tm), tm_dev(tm_hi, min_tm, max_tm))

    
def balanced_oligo_score_fn(tm, left_bound, right_bound, min_len, max_len, min_tm, max_tm, max_untiled_len):
    """
    Favors short oligos with balanced upstream/downstream gap sizes.
    Secondarily favors oligos with Tms within the tm range
    """
    
    if tm[1]-tm[0] < min_len or tm[1]-tm[0] > max_len:
        return 0
    
    left_untiled_len = tm[0] - left_bound
    right_untiled_len = right_bound - tm[1]
    if left_untiled_len > max_untiled_len or right_untiled_len > max_untiled_len or left_untiled_len < 0 or right_untiled_len < 0:
        return 0

    balance_score = left_untiled_len * right_untiled_len

    tm_dev_score = tm_range_dev(tm[2], tm[3], min_tm, max_tm)
    if tm_dev_score == 0:
        if balance_score == max_untiled_len * max_untiled_len:
            return 1000000        
        return balance_score
    else:
        return -1 * (tm_dev_score + 0.1) #prevent 0

    
def unlimited_untiled_oligo_score_fn(tm, left_bound, right_bound, min_len, max_len, min_tm, max_tm, max_untiled_len):
    """
    Does not impose max_untiled_len restriction
    """
    return balanced_oligo_score_fn(tm, left_bound, right_bound, min_len, max_len, min_tm, max_tm, max_untiled_len = float("inf"))
    

def find_oligo(tms, left_bound, right_bound, min_len, max_len, max_untiled_len, min_tm, max_tm, score_fn = balanced_oligo_score_fn, return_all = False):
    """
    Finds an oligo preferentially in the center of the bounds with
    optimal Tm. On failure, returns the best Tm satisfying length
    requirements. Returns () if there is no solution
    """

    #Determine the minimum acceptable size
    target_min = right_bound - left_bound - 2*max_untiled_len

    if target_min < min_len or target_min > max_len:
        target_min = min_len
        
#    if target_min < min_len or target_min > max_len:
#        print ("FIND: UNACCEPTABLE SIZE")
#        return () #Unacceptable size given constraints

    #Choose a starting point in the middle
    start = (left_bound + right_bound)//2 - min_len//2

    #Search positions radiating out from the middle. If an
    #optimal solution is found, return immediately, otherwise
    #return the best (shortest, most centered), which is optimized
    #by the max of left_untiled_len * right_untiled_len
    candidates = []
    
    sign = -1
    counter = 1
    
    while left_bound <= start <= right_bound:
        for tm in tms[start]:
            score = score_fn(tm, left_bound, right_bound, target_min, max_len, min_tm, max_tm, max_untiled_len)
            if score == 1000000 and not return_all:
                return tm
            elif score:
                candidates.append((score, tm))
                
        start += sign*counter
        counter += 1
        sign *= -1
        
    #Return the max
    if not candidates:
#        print('CANDIDATE NOT FOUND FOR ', left_bound, right_bound)
        return ()
    if return_all:
        candidates.sort(key = lambda x: x[0], reverse=True)
        return [ x[1] for x in candidates ]
    else:
        return max(candidates, key = lambda x: x[0])[1]


####
# Overall strategy: partition the sequence into blocks. Go left to right
# finding acceptable oligos that satisfy length requirements, but
# not necessarily tm.
# On refinement, if there is a suboptimal tm, attempt to re-find
# oligos in the subsection of +1/-1 the suboptimal, starting with the
# suboptimal.


def find_3oligos(tms, oligo_block, min_tm, max_tm, min_len, max_len, max_untiled_len):
    """
    Find all optimal oligos in an extended window in the middle, and
    attempt to find compatible left and right oligos

    At least two oligos are required

    oligo_block: [-2, -1, 0, 1, 2]
    """

    left_bound = oligo_block[0][1] + min_len

    if len(oligo_block) > 4:
        right_bound = oligo_block[4][0] - min_len
    else:
        right_bound = oligo_block[3][0]
        
    oligo_sets = []
    oligo_cache = {}

    middle_oligos = find_oligo(tms, left_bound, right_bound, min_len, max_len, max_untiled_len, min_tm, max_tm, score_fn = unlimited_untiled_oligo_score_fn, return_all = True)

    for oligo in middle_oligos:
        tm_diff = tm_range_dev(oligo[2], oligo[3], min_tm, max_tm)
#        print(oligo, tm_diff)
        
        #Find left oligo consistent with middle oligo
        ll = oligo_block[0][1]
        lr = oligo[0]

        if (ll,lr) not in oligo_cache:
            oligo_cache[(ll,lr)] = find_oligo(tms, ll, lr, min_len, max_len, max_untiled_len, min_tm, max_tm, return_all = False)

        left_oligo = oligo_cache[(ll,lr)]
        if not left_oligo:
            if 0 <= oligo[0] - oligo_block[1][1] <= max_untiled_len:
                tm_diff_l = 1000 #big
            else:
                tm_diff_l = 10000
#            continue
        else:
            tm_diff_l = tm_range_dev(left_oligo[2], left_oligo[3], min_tm, max_tm)
#        print(left_oligo, oligo)
        #Find right
        if len(oligo_block) > 4:
            rl = oligo[1]
            rr = oligo_block[4][0]
            
            if (rl,rr) not in oligo_cache:
                oligo_cache[(rl,rr)] = find_oligo(tms, rl, rr, min_len, max_len, max_untiled_len, min_tm, max_tm, return_all = False)

            right_oligo = oligo_cache[(rl,rr)]
            if not right_oligo:
                if 0 <= oligo_block[3][0] - oligo[1] <= max_untiled_len:
                    tm_diff_r = 1000 #big
                else:
                    tm_diff_r = 10000
#                continue
            else:
                tm_diff_r = tm_range_dev(right_oligo[2], right_oligo[3], min_tm, max_tm)
        else:
            #Check the middle oligo with respect to end coord
            if right_bound - oligo[1] > max_untiled_len: #illegal bound
                continue
            
            tm_diff_r = 0
            right_oligo = ()
            
        #Evaluate oligos
        if tm_diff == 0 and tm_diff_l == 0 and tm_diff_r == 0:
            return (left_oligo, oligo, right_oligo)
        else:
            score = 1000 * (tm_diff + 0.1) * (tm_diff_l + 0.1) * (tm_diff_r + 0.1)
            oligo_sets.append((score, (left_oligo, oligo, right_oligo)))
            
    #Return the best
    if oligo_sets:
        return min(oligo_sets, key = lambda x: x[0])[1]
    else:
        return ()


def oligo_set_score(oligos, min_tm, max_tm, min_len, max_len, max_untiled_len):
    length_violations = 0
    tm_deviations = 1

    for i, oligo in enumerate(oligos[1:-1]):
        if oligo[1] - oligo[0] < min_len or oligo[1] - oligo[0] > max_len:
            length_violations += 1
        if oligos[i][1] - oligo[0] > max_untiled_len or oligos[i+1][0] - oligo[1] > max_untiled_len:
            length_violations += 1

        tm_deviations *= (tm_range_dev(oligo[2], oligo[3], min_tm, max_tm) + 0.1)
        
    return length_violations, round(tm_deviations ** (1.0/(len(oligos)-2)), 3)

        
def find_oligo_set(seq, tms, target_oligo_len, min_tm = 70, max_tm = 80, min_len = 40, max_len = 50, max_untiled_len = 25):

    #Estimate the number of oligos
    min_n_oligos = math.ceil((len(seq) - max_untiled_len) / (max_len + max_untiled_len))
    max_n_oligos = len(seq) // min_len

    candidate_oligo_sets = []
    
    #Iterate from low to high number of oligos, until a good solution is found
    for n_oligos in range(min_n_oligos, max_n_oligos + 1): 
        max_untiled_allowed = (n_oligos+1) * max_untiled_len
        expected_untiled = min(len(seq) - (n_oligos * target_oligo_len), max_untiled_allowed)
        target_untiled_len = expected_untiled // (n_oligos+1)
        remaining_nts = expected_untiled % (n_oligos+1)

        #Estimate oligo coords
        oligos = [(0,0,-1,-1)]
        for i in range(n_oligos):
            if i == 0:
                start = target_untiled_len
            else:
                start = oligos[-1][1] + target_untiled_len
            if remaining_nts:
                start += 1
                remaining_nts -= 1
            oligos.append((start, start+target_oligo_len, -1, -1))
        oligos.append((len(seq),len(seq),-1,-1))

        #Oligos are generated in blocks of 3: a candidate oligo in the middle
        #is found, then the left and right oligos are generated to
        #satsify constraints.

        if n_oligos > 1:
            for i in range(2, max(3, len(oligos)-2)):
                oligo_block = find_3oligos(tms, oligos[i-2:i+3], min_tm, max_tm, min_len, max_len, max_untiled_len)        
                for j, oligo in enumerate(oligo_block):
                    if oligo:
                        oligos[i-1+j] = oligo

        #Validate the solution -- if not, try again with one more oligo
        len_score, tm_score = oligo_set_score(oligos, min_tm, max_tm, min_len, max_len, max_untiled_len)
        if len_score == 0:
           #Final refinement
            if True:
                for i in range(1, len(oligos)-1):
                    left_bound = oligos[i-1][1]
                    right_bound = oligos[i+1][0]
                    new_oligo = find_oligo(tms, left_bound, right_bound, min_len, max_len, max_untiled_len, min_tm, max_tm, return_all = False)
                    if new_oligo:
                        oligos[i] = new_oligo
                        
            #Check again
            len_score, tm_score = oligo_set_score(oligos, min_tm, max_tm, min_len, max_len, max_untiled_len)
            if tm_score == 0.1:
                return oligos
            else:
                candidate_oligo_sets.append((tm_score, oligos))
        else:
            print("Length violation", oligos)
                
    if candidate_oligo_sets:
        return min(candidate_oligo_sets, key = lambda x: x[0])[1]
    else:
        print("No acceptable solution found")


def tile_oligos_with_gaps(seq, min_tm = 70, max_tm = 80, min_len = 40, max_len = 50, C = 0.00000005, Na_concen = 0.2, molec = 'dnarna', max_untiled_len = 25):
    """
    Returns an oligo list tuple: (start, end, Tm)
    """
    if max_untiled_len > len(seq):
        max_untiled_len = len(seq)-1
    
    #Check if a solution exists based on input parameters
    max_n_oligos = len(seq) // min_len
    remaining_nts = len(seq) % min_len
    possible_slop = max_untiled_len + max_n_oligos * (max_len - min_len + max_untiled_len)
    if possible_slop < remaining_nts:
        print("No possible solutions")
        return ()

    #Generate Tms for candidate oligos
    tms, median_oligo_len = generate_tms(seq, min_tm, max_tm, min_len, max_len, C, Na_concen, molec, strict = False)
        
    #Attempt to find an optimal oligo set
    oligos = find_oligo_set(seq, tms, median_oligo_len, min_tm, max_tm, min_len, max_len, max_untiled_len)

    return oligos


def pretty_print_oligos(seq, oligos):
    """
    Also count number of wildcard positions, expansions, oligo seq
    """
    print('#Target sequence: %d nts' % (len(seq)))
    print('\t'.join(['Start', 'End', 'Length', 'Tm_low', 'Tm_high', 'X_pos', 'Ambig_pos', 'Num_targets', 'Target_seq', 'Antisense_oligo']))
    for oligo in oligos[1:-1]:
        subseq = seq[oligo[0]:oligo[1]]
        x_count, wc_count, expand_count = wildcard_stats(subseq)
        print ('\t'.join([str(oligo[0]+1), str(oligo[1]), str(oligo[1]-oligo[0]), str(oligo[2]), str(oligo[3]), str(x_count), str(wc_count), str(expand_count), subseq, bioutil.rc(subseq)]))

    

#print(tm('AAAAACCCCCGGGGGTTTTT', dna_rna = False))
#print(tm('AAAAACCCCCGGGGGTTTTT', molec = 'dnarna'))
#print(tm('AAAAACCCCCGGGGGT', molec = 'dnarna'))
#print(Tm('TGGCTTAATCTTTGAGACAAGCATATGCTACCTGGCAGGATCAACCAGGT'))
#print(Tm_v2('TGGCTTAATCTTTGAGACAAGCATATGCTACCTGGCAGGATCAACCAGGT'))
#print(Tm_v2('CGCGTACGCGTACGCG', dna_rna = False))
#print(Tm_v2('TTGTAATCCATT', .00000005, .2))
#print(expand_wildcards('AGTAAATTATGC'))

#print(tm('AAAAACCCCC'))
#print(tm('GGGGGTTTT'))


#tt = generate_tms('AAAAACCCCCGGGGGTTTTT', min_tm = 10, max_tm = 40, max_len = 15)
#for t in tt:
#    print(t)

#print(tile_oligos('TGGCTTAATCTTTGAGACAAGCATATGCTACCTGGCAGGATCAACCAGGT', min_tm = 10, max_tm = 30, max_len = 15, C = 0.00000005, Na_concen = 0.2, molec = 'dnarna'))

#tile_oligos_with_gaps('TGGCTTAATCTTTGAGACAAGCATATGCTACCTGGCAGGATCAACCAGGT', min_tm = 10, max_tm = 30, max_len = 15, C = 0.00000005, Na_concen = 0.2, molec = 'dnarna', max_untiled_len = 7)

#print(tile_oligos_with_gaps('GACTCTTAGCGGTGGATCACTCGGCTCGTGCGTCGATGAAGAACGCAGCTAGCTGCGAGAATTAATGTGAATTGCAGGACACATTGATCATCGACACTTCGAACGCACTTGCGGCCCCGGGTTCCTCCCGGGGCTACGCCTGTCTGAGCGTCGGTTG', min_tm=70, max_tm=80,max_untiled_len = 25)) #mouse 5.8S
#seq = 'GACTCTTAGCGGTGGATCACTCGGCTCGTGCGTCGATGAAGAACGCAGCTAGCTGCGAGAATTAATGTGAATTGCAGGACACATTGATCATCGACACTTCGAACGCACTTGCGGCCCCGGGTTCCTCCCGGGGCTACGCCTGTCTGAGCGTCGGTTG'
seq = 'GACTCTTAGCRGYGGATXACTCGGCTCGTGCGTCGATGAAGAACGCAGCTAGCTGCGAGAATTAATGTGAATTGCAGGACACATTGATCATCGACACTTCGAACGCACTTGCGGCCCCGGGTTCCTCCCGGGGCTACGCCTGTCTGAGCGTCGGTTG'
pretty_print_oligos(seq,tile_oligos_with_gaps(seq, min_len = 40, max_len = 50, min_tm=70, max_tm=80,max_untiled_len = 25)) #mouse 5.8S


# 1. from a multiple sequence alignment, generate a consensus sequence
# 2. Find oligos that work with the consensus
# 3. For each individual sequence, mask the found oligos and attempt
# to generate additional oligos with the pieces


frag_28s = 'atggatggcgctggagcgtcgggcccatacccggccgtcg'
frag_28s = 'ccgggttaaggcgcccgatgccgacgctcatcagacccca'

#print(do_tm(frag_28s, C = 0.00000005, Na_concen = 0.2, molec = 'dnarna'))

#seqs = bioutil.read_fasta('zfm_28s.fa')
#pretty_print_oligos(seqs[0][1],tile_oligos_with_gaps(seqs[0][1], min_len = 40, max_len = 50, min_tm=70, max_tm=80,max_untiled_len = 25))


