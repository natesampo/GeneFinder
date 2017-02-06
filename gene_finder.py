# -*- coding: utf-8 -*-
"""
GENE FINDER FTW

@author: Nate Sampo

"""

import random
import math
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
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    if nucleotide == ('A'):
        return 'T'
    if nucleotide == ('C'):
        return 'G'
    if nucleotide == ('T'):
        return 'A'
    if nucleotide == ('G'):
        return 'C'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    reverse_dna_complement = ''
    for i in range(0,len(dna)):
        reverse_dna_complement += get_complement(dna[(len(dna)-1)-i])
    return reverse_dna_complement


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    if dna.count('TAG') > 0 and dna.index('TAG') % 3 == 0:
        return dna[:dna.index('TAG')]
    elif dna.count('TAA') > 0 and dna.index('TAA') % 3 == 0:
        return dna[:dna.index('TAA')]
    elif dna.count('TGA') > 0 and dna.index('TGA') % 3 == 0:
        return dna[:dna.index('TGA')]
    else:
        return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    all_ORFs = []
    counter = -1
    for i in range(0,len(dna)):
        if dna[i:i+3] == 'ATG' and i % 3 == 0:
            counter += 1
            all_ORFs.append(str(rest_of_ORF(dna[i:]))[0:])
            dna = dna.replace(all_ORFs[counter],' ' * len(all_ORFs[counter]))
    return all_ORFs



def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    all_ORFs = []
    for i in range(0,3):
        if find_all_ORFs_oneframe(dna[i:]) != '':
            this_frame = find_all_ORFs_oneframe(dna[i:])
            for t in range(0,len(this_frame)):
                all_ORFs.append(this_frame[t])
    return all_ORFs


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    return find_all_ORFs(dna) + find_all_ORFs(get_reverse_complement(dna))


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    big_ORF = ''
    all_ORFs = find_all_ORFs_both_strands(dna)
    for i in range(0, len(all_ORFs)):
        if len(all_ORFs[i]) > len(big_ORF):
            big_ORF = all_ORFs[i]
    return big_ORF


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    random_ORF_length = 0
    for i in range(0,num_trials):
        temp = len(longest_ORF(shuffle_string(dna)))
        if random_ORF_length < temp:
            random_ORF_length = temp
    return random_ORF_length


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    amino_acid = ''
    for i in range(0,math.floor(len(dna)/3)):
        if sum(dna[i*3:(i*3)+3].count(x) for x in ('A','T','G','C')) == 3:
            amino_acid += aa_table[dna[i*3:(i*3)+3]]
    return amino_acid


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    the_answer = []
    threshold = longest_ORF_noncoding(dna, 1500)/2
    all_of_the_ORFs = find_all_ORFs_both_strands(dna)
    for i in range(0,len(all_of_the_ORFs)):
        if len(all_of_the_ORFs[i]) >= threshold:
            the_answer.append(coding_strand_to_AA(all_of_the_ORFs[i]))
    return the_answer


if __name__ == "__main__":
    from load import load_seq
    dna = load_seq("./data/X73525.fa")
    gene = gene_finder(dna)
    for g in gene:
        print(g)
        print('-'*50)
