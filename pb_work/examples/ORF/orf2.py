"""=================================================================================================
Given a DNA Sequence, what are all the open reading frames (ORF) longer than L
    lengths in nucleotides
    lengths in residues
    amino acid sequence

Michael Gribskov     15 January 2021
================================================================================================="""
codon2aa = {"AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N",
            "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
            "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S",
            "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I",

            "CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H",
            "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
            "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
            "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",

            "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D",
            "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
            "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
            "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",

            "TAA": "*", "TAC": "Y", "TAG": "*", "TAT": "T",
            "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
            "TGA": "*", "TGC": "C", "TGG": "W", "TGT": "C",
            "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F"}
# --------------------------------------------------------------------------------------------------
# Main program
# --------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    # Read nucleotide file
    # for now just a dummy sequence string
    seq = 'CATATCATCGATCTACGTACGTATCGAGCTATCGGACTCGATCGATCGATCGTCGATCGTAGCATGCTAGCTAGCTAG'

    # check codon table

    # Translate each reading frame and save ORFs

    # ncodon = 0
    #     # for codon in codon2aa:
    #     #     ncodon += 1
    #     #     print(codon, codon2aa[codon])
    #     #
    #     # print('number of codons:', ncodon)
    start = 0
    while start < len(seq):
        codon = seq[start:start+3]
        print(start, codon, codon2aa[codon])

        start+= 3

    # Report results

    exit(0)
