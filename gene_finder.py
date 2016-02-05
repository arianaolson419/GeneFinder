# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Ariana Olson

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide

        the program should not take other letters

    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('Z')
    False
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    elif nucleotide == 'T':
        return 'A'
    else:
        print nucleotide + ' is not a valid nucleotide' #error message if the nucleotide is not a valid letter
        return False
    


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
        This seems to be sufficient because the get_complement function checked for non-nucleotides

    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    index = 1
    complement = ''
    reverse_complement = ''
    #creates sting of all complement nucleotides
    for nucleotide in dna:
        complement = complement + get_complement(nucleotide)
    #reverses order of the complements
    while index <= len(complement):
        reverse_complement = reverse_complement + complement[-index]
        index += 1
    return reverse_complement
    
get_reverse_complement('ATGCCGTTT')

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        need doctest to check that whole string is returned if the sequence begins with a stop codon
        need doctest to check if full string is returned if no in frame stop codon

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("TGAACGCGG")
    'TGAACGCGG'
    >>> rest_of_ORF('CGAAAAT')
    'CGAAAAT'
    """
    k = 0
    snippet = '' 
    #checks for stop codon at beggining. returns full string if True
    if dna[0:3] == 'TAG' or dna[0:3]=='TAA' or dna[0:3]=='TGA': 
        return dna
    #searches for stop codon, returns everything before stop, longest 
    #string with length a multiple of 3 if stop DNE
    else:
        while k + 2 < len(dna):
            #stops loop if stop codon present
            if dna[k:k+3] == 'TAG' or dna[k:k+3] == 'TGA' or dna[k:k+3] == 'TAA':
                break
            #adds next three letters and advances check to next chunk of three letters
            else:
                snippet = snippet + dna[k:k+3]
                k = k +3
        #prints the extra letters at the end if no stop codon found
        if float(k +1/3) != (k+1)/3.0 and dna[k:k+3] != 'TAG' and dna[k:k+3]!= 'TAA' and dna[k:k+3] != 'TGA':
            snippet = snippet + dna[k:]
        return snippet

    

def find_all_ORFs_oneframe(dna, start_at_index):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        needed a test for empty list if ORF could not start on a multiple of 3
        needed a test to confirm ORFs not nested
    >>> find_all_ORFs_oneframe('TATGCCCC', 0)
    []
    >>>find_all_ORFs_oneframe('ATGCGAATGTAGCAT', 0)
    [ATGCGAATG]
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC", 0)
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    list_orf = []
    i = 0 + start_at_index
    while i < len(dna):
        #checks that the initial position minus the start offset is a multiple of three
        if (i - start_at_index)/3.0 == (i - start_at_index)/3:
            #starts the ORF at the first ATG on a multiple of 3
            if dna[i:i+3] == 'ATG':
                #'ATG' + the sequence until the next stop codon in the strand is found
                orf = dna[i:i+3] + rest_of_ORF(dna[i+3:])
                #advances position, should be a multiple of three from the start
                i = i + len(orf)
                #adds the ORF to the list
                list_orf.append(orf)
            else:
                #advances three positions forward to check for 'ATG' in the next chunk
                i = i + 3
    return list_orf
    


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        The doctest in the find_all_ORFs_oneframe function checked for non-nested frames 
        which should be sufficient. No further doctests needed for this function

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    j = 0
    all_frames = []
    while j < 3:
        one_frame = find_all_ORFs_oneframe(dna, j)
        for i in one_frame:
            all_frames.append(rest_of_ORF(i))
        j = j + 1
    return all_frames
    
find_all_ORFs('ATGCATATGTAGGATGC')

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        again, non-nesting has been checked for so no more doctrsts are needed

    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    reverse_complement = get_reverse_complement(dna)
    forwards = find_all_ORFs(dna)
    backwards = find_all_ORFs(reverse_complement)
    return forwards + backwards

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string

        needs a test to make sure the function can pick the largest element 
        out of lists containing more than two entries
        needs to return an empty string for an empty list

    >>> longest_ORF('ATGCGAATGTAGCATCAAAATGCAT')
    'ATGCATTTTGATGCTACATTCGCAT'
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    >>> longest_ORF('AAAAAAAAA')
    ''
    """
    ORF_list = find_all_ORFs_both_strands(dna)
    #each iteration compares the 'first' element of the list with the next
    #the varible first is reassigned to the larger entry of the two
    #moves on to the next entry to compare with 'first'
    #once all element have been cycled through, the largest element in the list is assigned to 'first'
    
    if len(ORF_list) > 0:
        i = 0
        first = ORF_list[0]
        for index in ORF_list:
            if len(ORF_list) > i + 1:
                second = ORF_list[i+1]
                if len(first) < len(second):
                    first = second  
                    i = i + 1
                else:
                    i = i + 1
            #stops comparison is only one list element
            else:
                break
        return first
    else:
        return ''
  
import random

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles returns: the maximum length longest ORF 
        It's hard to do a doctest for a random outcome, but I think that the more shuffles, 
        the higher the probability will be that the longest ORF will be as long as the dna strand
        I would also like my function to handle strings of zero lengths
        >>> longest_ORF_noncoding('', 1)
        0
        >>> longest_ORF_noncoding('ATGCGAATGTAGCATCAAAATGCAT',200)
        25
        """
    strands = []
    longest_orfs = []
    lengths = []
    i = 0
    #generates a list of all orfs of all trials
    while i < num_trials:
        #shuffles original sequence
        rand_strand = shuffle_string(dna)
        #finds all the ORFs of that string, creates a list
        strands.append(rand_strand)
        #adds the longest orf to a list of the longest orfs in each trial
        longest_orfs.append(longest_ORF(rand_strand))
        i += 1
    #returns the length of the longest string
    return len(max(longest_orfs, key=len))

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment
        I don't think anything needs to be tested besides returning a string of 
        amino acids and being able to handle dna strands that aren't a multiple
        of 3 long

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    aa_sequence = ''
    i = 0
    #shortens dna until its length is a multiple of 3
    while len(dna)%3 != 0:
        dna = dna[:len(dna)-1]
    #looks up the codon and adds an amino acid to the sequence
    while i<len(dna):
        codon = dna[i:i+3]
        amino_acid = aa_table[codon]
        aa_sequence= aa_sequence + amino_acid
        i = i+3

    return aa_sequence



def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.

        I can't do much to verify that my amino acid squence is right, but my gene_finder 
        should return a list

        It is heavily reliant on other previously tested functions, so no other doctest 
        are necessary

        >>> type(gene_finder(dna))
        <type 'list'>
    """
    genes = []
    amino_acid_sequence = []
    threshold = longest_ORF_noncoding(dna, 1500)
    ORFs = find_all_ORFs_both_strands(dna)
    for ORF in ORFs:
        if len(ORF) > threshold:
            genes.append(ORF)
    for gene in genes:
        amino_acid_sequence.append(coding_strand_to_AA(gene))
    return amino_acid_sequence

from load import load_seq
dna=load_seq("./data/X73525.fa")

print gene_finder(dna)
print'done'
if __name__ == "__main__":
    import doctest
    doctest.run_docstring_examples(gene_finder, globals())
