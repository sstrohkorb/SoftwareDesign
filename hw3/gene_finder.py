# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 11:24:42 2014

@author: Sarah Strohkorb
"""

# you may find it useful to import these variables (although you are not required to use them)
from amino_acids import aa, codons

def collapse(L):
    """ Converts a list of strings to a string by concatenating all elements of the list """
    output = ""
    for s in L:
        output = output + s
    return output

def two_lists_contain_same_elements(list1, list2):
    if len(list1) != len(list2):
        return False
    else: 
        for list_item in list1:
            if list_item in list2:
                continue
            else:
                return False
        return True


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).
        
        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment
    """
    amino_acid_output = ""
    for i in xrange(0,len(dna),3):
        my_codon = dna[i:i+3]
        for j in range(len(codons)):
            if my_codon in codons[j]:
                amino_acid_output += aa[j]
                break
    return amino_acid_output


def coding_strand_to_AA_unit_tests():
    """ Unit tests for the coding_strand_to_AA function """

    # DNA input strands 
    dna_input1 = "ACTGCCCC"
    dna_input2 = "AGCTGAGGGTGTTTTGGA"
    dna_input3 = "CAGGCTTGCGGCTTCTTAA"

    # Expected output amino acid strands
    e_output1 = "TA"
    e_output2 = "S|GCFG"
    e_output3 = "QACGFL"

    # Actual output amino acid strands 
    a_output1 = coding_strand_to_AA(dna_input1)
    a_output2 = coding_strand_to_AA(dna_input2)
    a_output3 = coding_strand_to_AA(dna_input3)

    test1_result = (e_output1 == a_output1)
    test2_result = (e_output2 == a_output2)
    test3_result = (e_output3 == a_output3)

    if not test1_result:
        print "Test 1 FAILED: " + str(a_output1) + " != " + str(e_output1)
    if not test2_result:
        print "Test 2 FAILED: " + str(a_output2) + " != " + str(e_output2)
    if not test3_result:
        print "Test 3 FAILED: " + str(a_output3) + " != " + str(e_output3)

    if test1_result and test2_result and test3_result:
        return True
    else:
        return False
    

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
    
        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    """
    reverse_complement_output = ""
    for i in range(len(dna)):
        base = dna[len(dna)-(i+1)]
        complementary_base = ""
        if base == 'T':
            complementary_base = 'A'
        elif base == 'A':
            complementary_base = 'T'
        elif base == 'G':
            complementary_base = 'C'
        elif base == 'C':
            complementary_base = 'G'
        reverse_complement_output += complementary_base
    return reverse_complement_output
    
def get_reverse_complement_unit_tests():
    """ Unit tests for the get_complement function """

    # DNA input strands 
    dna_input1 = "ACTGCCCC"
    dna_input2 = "AGCTGAGGGTGTTTTGGA"
    dna_input3 = "CAGGCTTGCGGCTTCTTAA"

    # Expected output DNA strands
    e_output1 = "GGGGCAGT"
    e_output2 = "TCCAAAACACCCTCAGCT"
    e_output3 = "TTAAGAAGCCGCAAGCCTG"

    # Actual output DNA strands
    a_output1 = get_reverse_complement(dna_input1)
    a_output2 = get_reverse_complement(dna_input2)
    a_output3 = get_reverse_complement(dna_input3)

    test1_result = (e_output1 == a_output1)
    test2_result = (e_output2 == a_output2)
    test3_result = (e_output3 == a_output3)

    if not test1_result:
        print "Test 1 FAILED: " + str(a_output1) + " != " + str(e_output1)
    if not test2_result:
        print "Test 2 FAILED: " + str(a_output2) + " != " + str(e_output2)
    if not test3_result:
        print "Test 3 FAILED: " + str(a_output3) + " != " + str(e_output3)

    if test1_result and test2_result and test3_result:
        return True
    else:
        return False


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start codon and returns
        the sequence up to but not including the first in frame stop codon.  If there
        is no in frame stop codon, returns the whole string.
        
        dna: a DNA sequence
        returns: the open reading frame represented as a string
    """
    ORF_output = ""
    for i in xrange(0,len(dna),3):
        my_codon = dna[i:i+3]
        if my_codon == 'TAA' or my_codon == 'TAG' or my_codon == 'TGA':
            return ORF_output
        else:
            ORF_output += my_codon
    return ORF_output

def rest_of_ORF_unit_tests():
    """ Unit tests for the rest_of_ORF function """

    # DNA input strands 
    dna_input1 = "ATGTTTAAGCCGTAAATGAAACCGGGC"
    dna_input2 = "ATGTAG"
    dna_input3 = "ATGCCGATAGCCTGACCGATAAAATTG"
    dna_input4 = "ATGC"

    # Expected output DNA strands
    e_output1 = "ATGTTTAAGCCG"
    e_output2 = "ATG"
    e_output3 = "ATGCCGATAGCC"
    e_output4 = "ATGC"

    # Actual output DNA strands
    a_output1 = rest_of_ORF(dna_input1)
    a_output2 = rest_of_ORF(dna_input2)
    a_output3 = rest_of_ORF(dna_input3)
    a_output4 = rest_of_ORF(dna_input4)

    test1_result = (e_output1 == a_output1)
    test2_result = (e_output2 == a_output2)
    test3_result = (e_output3 == a_output3)

    if not test1_result:
        print "Test 1 FAILED: " + str(a_output1) + " != " + str(e_output1)
    if not test2_result:
        print "Test 2 FAILED: " + str(a_output2) + " != " + str(e_output2)
    if not test3_result:
        print "Test 3 FAILED: " + str(a_output3) + " != " + str(e_output3)

    if test1_result and test2_result and test3_result:
        return True
    else:
        return False

        
def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence and returns
        them as a list.  This function should only find ORFs that are in the default
        frame of the sequence (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    """
    all_ORFs_output = []
    dna_length = len(dna)
    current_dna_location = 0
    while(current_dna_location < dna_length):
        my_codon = dna[current_dna_location : current_dna_location + 3]
        if my_codon == 'ATG':
            discovered_ORF = rest_of_ORF(dna[current_dna_location : dna_length])
            all_ORFs_output.append(discovered_ORF)
            current_dna_location += len(discovered_ORF)
        else:
            current_dna_location += 3
    return all_ORFs_output
     
def find_all_ORFs_oneframe_unit_tests():
    """ Unit tests for the find_all_ORFs_oneframe function """

    # DNA input strands 
    dna_input1 = "ATGCATGAATGTAGATAGATGTGCCC"
    dna_input2 = "ATGCCCGGGTATCCGGAAATAG"
    dna_input3 = "CCGTTTATGCCGTAGTTAGACATGCCCGAGTAAGCGATGTTTATAGGGC"

    # Expected output list of ORFs
    e_output1 = ['ATGCATGAATGTAGA', 'ATGTGCCC']
    e_output2 = ['ATGCCCGGGTATCCGGAAATAG']
    e_output3 = ['ATGCCG', 'ATGCCCGAG', 'ATGTTTATAGGGC']

    # Actual output list of ORFs
    a_output1 = find_all_ORFs_oneframe(dna_input1)
    a_output2 = find_all_ORFs_oneframe(dna_input2)
    a_output3 = find_all_ORFs_oneframe(dna_input3)

    test1_result = two_lists_contain_same_elements(e_output1, a_output1)
    test2_result = two_lists_contain_same_elements(e_output2, a_output2)
    test3_result = two_lists_contain_same_elements(e_output3, a_output3)

    if not test1_result:
        print "Test 1 FAILED: " + str(a_output1) + " != " + str(e_output1)
    if not test2_result:
        print "Test 2 FAILED: " + str(a_output2) + " != " + str(e_output2)
    if not test3_result:
        print "Test 3 FAILED: " + str(a_output3) + " != " + str(e_output3)

    if test1_result and test2_result and test3_result:
        return True
    else:
        return False

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in all 3
        possible frames and returns them as a list.  By non-nested we mean that if an
        ORF occurs entirely within another ORF and they are both in the same frame,
        it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    """
    dna_frame_0 = dna
    dna_frame_1 = dna[1:]
    dna_frame_2 = dna[2:]

    ORFs_frame_0 = find_all_ORFs_oneframe(dna_frame_0)
    ORFs_frame_1 = find_all_ORFs_oneframe(dna_frame_1)
    ORFs_frame_2 = find_all_ORFs_oneframe(dna_frame_2)

    return ORFs_frame_0 + ORFs_frame_1 + ORFs_frame_2


def find_all_ORFs_unit_tests():
    """ Unit tests for the find_all_ORFs function """

    # DNA input strands 
    # tests offset by 2 and 0 ORFs, multiple ORFs
    dna_input1 = "ATGCATGAATGTAGATAGATGTGCCC"
    # tests offset by two and only one ORF
    dna_input2 = "GCATGCCCGGGTATCCGGAAATAG"
    # nested reading frame, that won't be included
    dna_input3 = "GATGTTTATGGGCTAGTAG"

    # Expected output list of ORFs
    e_output1 = ['ATGCATGAATGTAGA', 'ATGTGCCC', 'ATGAATGTAGATAGATGTGCCC', 'ATG']
    e_output2 = ['ATGCCCGGGTATCCGGAAATAG']
    e_output3 = ['ATGTTTATGGGC']

    # Actual output list of ORFs
    a_output1 = find_all_ORFs(dna_input1)
    a_output2 = find_all_ORFs(dna_input2)
    a_output3 = find_all_ORFs(dna_input3)

    test1_result = two_lists_contain_same_elements(e_output1, a_output1)
    test2_result = two_lists_contain_same_elements(e_output2, a_output2)
    test3_result = two_lists_contain_same_elements(e_output3, a_output3)

    if not test1_result:
        print "Test 1 FAILED: " + str(a_output1) + " != " + str(e_output1)
    if not test2_result:
        print "Test 2 FAILED: " + str(a_output2) + " != " + str(e_output2)
    if not test3_result:
        print "Test 3 FAILED: " + str(a_output3) + " != " + str(e_output3)

    if test1_result and test2_result and test3_result:
        return True
    else:
        return False


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    """
    reversed_dna = get_reverse_complement(dna)
    print reversed_dna

def find_all_ORFs_both_strands_unit_tests():
    """ Unit tests for the find_all_ORFs_both_strands function """

    # DNA input strands 
    # reverse of 1 is 
    dna_input1 = "ATGCATGAATGTAGATAGATGTGCCC"
    dna_input2 = "GCATGCCCGGGTATCCGGAAATAG"
    dna_input3 = "GATGTTTATGGGCTAGTAG"

    # Expected output list of ORFs
    e_output1 = ['ATGCATGAATGTAGA', 'ATGTGCCC', 'ATGAATGTAGATAGATGTGCCC', 'ATG']
    e_output2 = ['ATGCCCGGGTATCCGGAAATAG']
    e_output3 = ['ATGTTTATGGGC']

    # Actual output list of ORFs
    a_output1 = find_all_ORFs_both_strands(dna_input1)
    a_output2 = find_all_ORFs_both_strands(dna_input2)
    a_output3 = find_all_ORFs_both_strands(dna_input3)

    test1_result = two_lists_contain_same_elements(e_output1, a_output1)
    test2_result = two_lists_contain_same_elements(e_output2, a_output2)
    test3_result = two_lists_contain_same_elements(e_output3, a_output3)

    if not test1_result:
        print "Test 1 FAILED: " + str(a_output1) + " != " + str(e_output1)
    if not test2_result:
        print "Test 2 FAILED: " + str(a_output2) + " != " + str(e_output2)
    if not test3_result:
        print "Test 3 FAILED: " + str(a_output3) + " != " + str(e_output3)

    if test1_result and test2_result and test3_result:
        return True
    else:
        return False

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string"""

    # YOUR IMPLEMENTATION HERE

def longest_ORF_unit_tests():
    """ Unit tests for the longest_ORF function """

    # YOUR IMPLEMENTATION HERE

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence
        
        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """

    # YOUR IMPLEMENTATION HERE

def gene_finder(dna, threshold):
    """ Returns the amino acid sequences coded by all genes that have an ORF
        larger than the specified threshold.
        
        dna: a DNA sequence
        threshold: the minimum length of the ORF for it to be considered a valid
                   gene.
        returns: a list of all amino acid sequences whose ORFs meet the minimum
                 length specified.
    """

    # YOUR IMPLEMENTATION HERE

if __name__ == "__main__":
    print str(coding_strand_to_AA_unit_tests()) + " - DNA to AA strand conversion"
    print str(get_reverse_complement_unit_tests()) + " - reverse complement of DNA strand"
    print str(rest_of_ORF_unit_tests()) + " - get ORF that starts with a start codon"
    print str(find_all_ORFs_oneframe_unit_tests()) + ' - find all ORFs one frame'
    print str(find_all_ORFs_unit_tests()) + ' - find all ORFs all frames'



