#!/usr/local/cs/bin/python3
import random as rand
import sys
import bisect
import time
import tracemalloc # for memory usage profiling

############################################################################
######## Instruction #######################################################
############################################################################

'''
This script can be run from shell.

First make sure the first line of this script is invoking Python from the
correct directory in your environment.

Then give execution permission using commands like:

    `chmod 744 ./pj4.py`

Then the script can be run with the following syntax.
    `./pj4.py [NUM_READS] [LEN_READ] [LEN_GENOME] [WHICH_ALGORITHM]`
    where:
        NUM_READS       - number of reads
        LEN_READ        - length of each read
        LEN_GENOME      - length of the genome
        WHICH_ALGORITHM - which algorithm to perform
            > '0' for PatternMatchingWithSuffixArray
            > '1' for BetterBWMatching
            > '2' for checking if the two algorithms's results agree
    e.g.
        `./pj4.py 100000 100 500000 0` means apply PatternMatchingWithSuffixArray
        to a genome of length 500,000 with 100,000 reads all of which are of
        length 100.
'''


############################################################################
######## General helper functions for both #################################
######## PatternMatchingWithSuffixArray and BetterBWMatching ###############
############################################################################

# Input: genome length
# Output: a randomly generated genome of length `glen` with 4 possible
#           nucleotides "ACTG"
def generate_genome(glen):
    genome = ''
    nucleotides = 'ACTG'
    for i in range(glen):
        genome += nucleotides[rand.randrange(4)]
    return genome

# Input: 1. rnum - number of reads to generate 
#        2. rlen - length of each read
#        3. genome - genome as a string 
#        4. glen - length of the genome
# Output: a list of reads represented using strings
def get_reads(rnum, rlen, genome, glen): 
    if rlen > glen:
        sys.stderr.write("Error: Read length must be smaller than genome length.\n")
        exit(1)
    last_index = glen - rlen # the largest possible index of the beginning of any read
    reads = []
 
    # generate one read each cycle
    for i in range(rnum):
        rand_index = rand.randrange(last_index+1)
        reads.append(genome[rand_index : (rand_index+rlen)])
    
    return reads
    

############################################################################
######### Functions for PatternMatchingWithSuffixArray #####################
############################################################################

# Input:    1. a string representing the whole genome
#           2. length of the genome
# Output: suffix array as a list
def get_suffixarray(genome, glen):
    class Suffix:
        def __init__(self, suffix, index):
            self.suffix = suffix
            self.index = index
        def __gt__(self, other):
            return self.suffix > other.suffix

    # generate suffixes and keep them in alphabetical order
    suffix_table = []
    for i in range(glen):
        bisect.insort(suffix_table, Suffix(genome[i:], i))
    # get the indices from the suffixes
    sufarray = []
    for suf in suffix_table:
        sufarray.append(suf.index)

    return sufarray

# Output: a tuple containing indices to the first and last match
#           between `read` and suffixes of `genome` in `sufarray` 
def pattern_matching_suffix_array(genome, read, sufarray):
    glen = len(genome)
    rlen = len(read)
    min_index = 0
    max_index = glen - 1
    mid_index = None
    first = None
    last = None
    genome_fract = ''

    # find first match
    while min_index <= max_index:
        mid_index = int((min_index + max_index) / 2)
        genome_fract =  genome[sufarray[mid_index] : min(glen, sufarray[mid_index]+rlen)]
        if read > genome_fract:
            min_index = mid_index + 1
            mid_index += 1
        
        elif mid_index == 0:
            # `read` smaller than first suffix in `sufarray`
            break
        else:
            max_index = mid_index - 1

    genome_fract =  genome[sufarray[mid_index] : min(glen, sufarray[mid_index]+rlen)]
    if read == genome_fract:
        # first match found
        first = min_index
    else:
        # no match found
        sys.stderr.write("Error: read not found in genome\n")
        exit(1)
        return None

    # find last match
    min_index = first
    max_index = glen - 1
    while min_index <= max_index:
        mid_index = int((min_index + max_index) / 2)
        if read == genome[sufarray[mid_index] : sufarray[mid_index]+rlen]:
            min_index = mid_index + 1
        else:
            max_index = mid_index - 1
    last = max_index

    return (first, last)

############################################################################
######### Functions for BetterBWMatching ###################################
############################################################################

# Name: get_lc_and_fo - get last column and first occurence
# Output: a tuple containing the last column as a list and the
#           first occurence as a dictionary
def get_lc_and_fo(genome, glen):
    class Table:
        def __init__(self, text, index):
            self.before = text[index:] + text[:index]
            self.last = text[index-1]
        def __gt__(self, other):
            return self.before > other.before
    
    # get all cyclic rotations of `genome`
    tables = []
    for i in range(glen):
        bisect.insort(tables, Table(genome, i))
    
    # get last column
    last_column = []
    for table in tables:
        last_column.append(table.last)
    
    # get first occurence
    first_occur = {}
    cur_symbol = None
    for i, table in enumerate(tables):
        if table.before[0] != cur_symbol:
            cur_symbol = table.before[0]
            first_occur[cur_symbol] = i

    return (last_column, first_occur)

# Output: the count matrix for BetterBWMatching
def get_count(last_column):
    # initialize the count matrix
    count = {}
    unique_symbols = set(last_column)
    for symbol in unique_symbols:
        count[symbol] = [0]
    
    # fill in the count matrix
    for i, target in enumerate(last_column):
        for symbol in unique_symbols:
            count[symbol].append(count[symbol][i])
        count[target][i + 1] =  1 + count[target][i + 1]

    return count

def better_bw_matching(first_occurence, last_column, read, count):
    top = 0
    bottom = len(last_column) - 1
    symbol = None
    while top <= bottom:
        if read != '':
            symbol = read[-1]
            read = read[0:-1]
            if symbol in last_column[top:bottom+1]:
                top = first_occurence[symbol] + count[symbol][top]
                bottom = first_occurence[symbol] + count[symbol][bottom+1] - 1
            else:
                return 0
        else:
            return bottom - top + 1

    sys.stderr.write("Error: reached unreachable region in better_bw_matching()\n")
    exit(1)

############################################################################
######### Main #############################################################
############################################################################

def main():
    # get arguments from command line
    if len(sys.argv) != 5:
        sys.stderr.write\
        (
"Error: Need four artuments:\n\
    1. number of reads\n\
    2. length of each reads\n\
    3. length of the genome\n\
    4. which algorithm to perform\n\
        - '0' for PatternMatchingWithSuffixArray\n\
        - '1' for BetterBWMatching\n\
        - '2' for checking if the two algorithms's results agree\n"
        )
        exit(1)
    rnum = int(sys.argv[1])
    rlen = int(sys.argv[2])
    glen = int(sys.argv[3])
    genome = generate_genome(glen)
    reads = get_reads(rnum, rlen, genome, glen)


    # PatternMatchingWithSuffixArray
    if sys.argv[4] == '0':
        tracemalloc.start()
        sufarray = get_suffixarray(genome, glen)
        start_time = time.time()
        num_match = 0
        for read in reads:
            res = pattern_matching_suffix_array(genome, read, sufarray)
            num_match += res[1] - res[0] + 1
        current, peak = tracemalloc.get_traced_memory()
        print("PatternMatchingWithSuffixArray:")
        print("- runtime:", str(time.time() - start_time), "seconds")
        print("- results:", str(num_match), "matches found")
        print(f"- current memory usage: {current / 10**6}MB; peak: {peak / 10**6}MB")
        # tracemalloc.stop()

    # BetterBWMatching
    elif sys.argv[4] == '1':
        tracemalloc.start()
        last_column, first_occurence = get_lc_and_fo(genome, len(genome))
        count = get_count(last_column)
        start_time = time.time()
        num_match = 0
        for read in reads:
            res = better_bw_matching(first_occurence, last_column, read, count)
            num_match += res
        current, peak = tracemalloc.get_traced_memory()
        print("BetterBWMatching:")
        print("- runtime:", str(time.time() - start_time), "seconds")
        print("- results:", str(num_match), "matches found")
        print(f"- current memory usage: {current / 10**6}MB; peak: {peak / 10**6}MB")
        tracemalloc.stop()
    
    # check if the two algorithms agree
    elif sys.argv[4] == '2':
        genome = generate_genome(glen)
        reads = get_reads(rnum, rlen, genome, glen)

        # PatternMatchingWithSuffixArray
        sufarray = get_suffixarray(genome, glen)
        num_match_array = 0
        for read in reads:
            res = pattern_matching_suffix_array(genome, read, sufarray)
            num_match_array += res[1] - res[0] + 1

        # BetterBWMatching
        last_column, first_occurence = get_lc_and_fo(genome, len(genome))
        count = get_count(last_column)
        maxi = 1
        num_match_bwm = 0
        for read in reads:
            res = better_bw_matching(first_occurence, last_column, read, count)
            maxi = max(maxi, res)
            num_match_bwm += res
        if num_match_bwm == num_match_array:
            print("Success. Results from the two algorithms matches")

    # error handling
    else:
        sys.stderr.write(
"""Error: invalid input for the fourth/last argument.
    - '0' for PatternMatchingWithSuffixArray
    - '1' for BetterBWMatching
    - '2' for checking if the two algorithms's results agree
"""     )
        exit(1)


if __name__ == "__main__":
    main()