""" IMPORT LIBRARIES """
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62
import re
import time

"""" DEFINE FUNCTIONS """

def trimNucleotide(seq):
    """
    Trim Nucleotide Sequence.
    Remove the overhanging bases so the sequence length is a multiple of three
    :param seq:
    :return: seq
    """
    leftover = len(seq) % 3
    return seq[:-leftover]

# ==================================================================================================
# main
# ==================================================================================================
# DEFINE VARIABLES
t = []  # Define Timer List.
Start, Stop = [], []  # Define Start and Stop Codon Location Variables.

# OPEN FILES
# Open Isoform List.
print("\nLoading Isoform List...")
Start_time = time.time()
# with open("IsoformsTopTen.txt", "r") as Isoforms:	# Top Ten Test Code.
with open("isoforms.txt", "r") as Isoforms:
    IsoformList = Isoforms.read().splitlines()
print("\tLoaded", len(IsoformList), "Isoforms in %1.3f seconds.\n" % (time.time() - Start_time))

"""# Index UniRef FASTA (Only need to do once, then load .idx file with SeqIO.index_db).
print("Indexing UniRef90 database...")
Start_time = time.time()
#ReferenceDB = SeqIO.index_db("uniref90TopTen.idx", "uniref90TopTen.fasta", "fasta")	# Top Ten Test Code
ReferenceDB = SeqIO.index_db("uniref90.idx", "uniref90.fasta", "fasta")
print("\tIndexed UniRef90 in %1.3f seconds.\n" % (time.time() - Start_time))"""

# Open Indexed UniRef90 FASTA.
print("Loading UniRef90 Database...")
Start_time = time.time()
# Reference = SeqIO.index_db("uniref90TopTen.idx")	# Top Ten Test Code.
Reference = SeqIO.index_db("uniref90.idx")
print("\tLoaded UniRef90 in %1.3f seconds.\n" % (time.time() - Start_time))

# Index Sequences From FASTA FILE.
print("Indexing Acomys Transcripts...")
Start_time = time.time()
# seq_record = SeqIO.index("acomysTopTen.transcripts.fa", "fasta")	# Top Ten Test Code.
seq_record = SeqIO.index("acomys_0.transcripts.fa", "fasta")
print("\tLoaded", len(seq_record),
      "Acomys Transcripts in %1.3f seconds." % (time.time() - Start_time))

# MAIN LOOP
for iso in range(1137, len(IsoformList)):  # For all sequences in Isoforms list,
    # Define Open Reading Frame Variables
    n_ORF = [0, 0, 0, 0, 0, 0]  # Define/Empty list of nucleotide open reading frames.
    p_ORF = [0, 0, 0, 0, 0, 0]  # Define/Empty list of peptide open reading frames.
    a_ORF = [0, 0, 0, 0, 0,
             0]  # Define/Empty list of longest AA sequence in each open reading frame.
    Alignment = [0, 0, 0, 0, 0, 0]  # Define/Empty list of pairwise2 alignment scores for each ORF.
    lens = [0, 0, 0, 0, 0, 0]  # Define/Empty list of length of longest AA sequence.

    # Load current Isoform.
    currentID = IsoformList[iso]  # Current Acomys Transcript ID.
    # Print which Transcript is being analyzed and the progress.
    print("\nCurrently Analyzing", currentID, "(%d/%d)..." % (iso + 1, len(IsoformList)))
    Start_time = time.time()

    # Get matched Reference ID from BLASTx results.
    # print("Searching BLASTx Results for Reference ID for %s..." % currentID)
    # Start_time = time.time()
    # BLASTxResult = open("dmndTopTen.f6.blastx")	# Top Ten Test Code.
    BLASTxResult = open("aocken_022318.uniref90.dmnd.f6.blastx")  # Open BLASTx file.
    currentResult = [line for line in BLASTxResult if
                     currentID in line]  # Search for current transcript.
    BLASTxResult.close()  # Close BLASTx file.

    UniRef_beg = str(currentResult).find("UniRef90_")  # Search string for "UniRef90_".
    UniRef_end = str(currentResult).find("\\", UniRef_beg)  # Search string for next "\".
    UniRef = str(currentResult)[UniRef_beg:UniRef_end]  # Define Reference ID.
    # print("\tLoaded", UniRef, "reference sequence in %1.3f seconds.\n" % (time.time() - Start_time))

    # Load Reference Protein Sequence.
    # Start_time = time.time()
    # print("Loading", UniRef, " Protein Sequence...")
    currentRef_seq = Reference[UniRef].seq  # Load Reference Sequence.
    currentRef_description = Reference[UniRef].description  # Load Reference Description.
    # print("\tLoaded in %1.3f seconds.\n" % (time.time() - Start_time))

    # Load Acomys Nucleotide Sequence.
    # Start_time = time.time()
    # print("Loading", currentID, "Nucleotide Sequence...")
    currentSeq = seq_record[currentID].seq  # Load current Transcript Sequence.
    currentLen = len(currentSeq)  # Length of current Transcript Sequence.
    revComplement = currentSeq.reverse_complement()  # Reverse complement of Current Sequence.
    # print("\tLoaded in %1.3f seconds.\n" % (time.time() - Start_time))

    # Generate Open Reading Frames
    # print("Generating ORFs for: %s..." % currentID)
    # Start_time = time.time()
    for n in range(0, 6):  # For index between zero and six,
        if (n < 3):  # if index is less than three,
            n_ORF[n] = trimNucleotide(
                currentSeq[n:currentLen])  # generate forward Open Reading Frames.
        else:  # Otherwise,
            n_ORF[n] = trimNucleotide(
                revComplement[(n - 3):currentLen])  # generate reverse Open Reading Frames.
    # print("\tGenerated ORFs for", currentID, "transcript in %1.3f seconds.\n" % (time.time() - Start_time))

    # Translate Open Reading Frames, Locate Codons, and List AA Sequence in ORFs
    for o in range(0, 6):  # For each Open Reading Frame,
        # translate the Open Reading Frame,
        p_ORF[o] = n_ORF[o].translate(stop_symbol="-")  # Translate the Open Reading Frame.

        # locate the Start and Stop Codons,
        currentORF = str(p_ORF[o])  # Current Open Reading Frame sequence.
        Start = [m.start() for m in re.finditer("M", currentORF)]  # Locate Start Codons.
        Stop = [s.start() for s in re.finditer("-", currentORF)]  # Locate Stop Codons.

        # and list the non-redundant AA sequences for the current Open Reading Frame.
        a = 0  # Define/Empty Index for number of non-redundant AA sequences.
        p_AA = []  # Define/Empty List of AA sequences.
        if (len(Start) == 0):  # If no Start Codons were identified,
            # print("ORF", o+1, "has no coding AA Sequences.")		# return no AA sequences.
            a = a
        elif (len(Stop) == 0):  # Otherwise, if no Stop Codons were identified,
            p_AA.append(p_ORF[o][Start[0]:len(
                p_ORF[o])])  # return the single non-redundant AA sequence identified,
            a_ORF[o] = p_AA[a]  # and define it as the AA sequence for current ORF.
        # print("ORF", o+1, "has 1 non-redundant AA Sequence:")
        # print(p_AA[a], "- Longest")
        else:  # Otherwise,
            Stop.append(len(p_ORF[o]))  # consider the entire length to be coding,
            for i in range(0, len(Start)):  # and for all Start Codons
                for j in range(0, len(Stop)):  # and for all Stop Codons,
                    if (Start[i] < Stop[j]):  # if Start comes before Stop,
                        p_AA.append(p_ORF[o][Start[i]:Stop[j]])  # list the AA sequence,
                        # print("ID'd:", p_AA[a])	# print the AA sequence identified,
                        a = a + 1  # and increase the list index,

            # Remove redundant AA sequences.
            # Remove AA sequences with Stop Codon in the middle.
            a = 0
            while (a < (len(p_AA))):  # For all AA sequences identified,
                if (str(p_AA[a]).find("-") != -1):  # if it contains a stop Codon,
                    # print("Removed: ", p_AA[a])
                    del p_AA[a]  # remove it from the list,
                    a = a - 1  # and recycle the index.
                else:  # Otherwise,
                    a = a + 1  # continue through the list.

            # Remove AA sequences with Start Codon
            a = 1  # Starting at the second sequence identified,
            while (a < (len(p_AA))):  # for all AA sequences identified,
                if (p_AA[a - 1].count("M") > 1):  # if the previous sequence contained an "M",
                    while (str(p_AA[a - 1]).find(
                            str(p_AA[a])) != -1):  # while the current sequence is redundant,
                        # print("Removed: ", p_AA[a])
                        del p_AA[a]  # remove it from the list,
                        a = a - 1  # and recycle the index.
                        if (len(p_AA) == 1):  # If its the last sequence,
                            break  # end the loop.
                a = a + 1

            # Print the number of non-redundant AA sequence identified.
            # print("ORF", o+1, "has", a, "non-redundant AA Sequence(s):")
            for a in range(0, len(p_AA)):  # For all non-redundant AA sequences identified,
                # Determine the longest non-redundant AA sequence in the ORF.
                if (len(p_AA[a]) == len(max(p_AA, key=len))):  # if it is the longest AA sequences,
                    # print("\t", p_AA[a], "- Longest")
                    a_ORF[o] = p_AA[a]  # define it as the AA sequence for current ORF.
            # else:
            # print(p_AA[a])

        # Print longest AA sequence and Align each Open Reading Frame with the Reference Sequence
        if (a_ORF[o] == 0):  # If ORF has no coding sequence,
            # print("ORF", o+1, "has no coding AA Sequences.")
            a = a  # do nothing.
        else:  # Otherwise,
            # Align the longest AA sequence with the Reference Sequence.
            Alignment[o] = pairwise2.align.localds(a_ORF[o], currentRef_seq, blosum62, -10, -0.5,
                                                   score_only=True)
        # print("ORF", o+1, ":", a_ORF[o], "\n\tAlignment Score:", Alignment[o])

        # Return Longest AA sequence for the Transcript and Alignment Score
        if (o == 5):  # If on the last ORF,
            for l in range(0, 6):  # for all ORF sequences,
                lens[l] = len(str(a_ORF[l]))  # return the length.
            """# Longest: 
            # I am not sure if this is necessary, my tests always confirm higher score is more accruate than longest.
            longest = a_ORF[lens.index(max(lens))]
            l_score = Alignment[a_ORF.index(longest)]
            print("\tLongest:\tORF %d: Length: %d, Alignment Score: %1.1f" % (a_ORF.index(longest)+1, len(longest), l_score))
            if (len(longest) > 75):
                print("\t\t\t%s...%s" % (longest[0:70], longest[-4:-1]))
            else:
                print("\t\t\t%s" % longest)
            #print("\t\t\t%s" % longest)"""
            # Highest Score:
            highest = max(Alignment)  # Determine the best score of all ORF sequences.
            h_seq = a_ORF[
                Alignment.index(highest)]  # Return the sequence of the best scoring sequence.
            print("\tHighest Score:\tORF %d: Length: %d, Alignment Score: %1.1f " % (
            Alignment.index(highest) + 1, len(h_seq), highest))
            if (len(h_seq) > 75):
                print("\t\t\t%s...%s" % (h_seq[0:70], h_seq[-3:-0]))
            else:
                print("\t\t\t%s" % h_seq)
            # print("\t\t\t%s" % h_seq)
            t.append(time.time() - Start_time)  # Update the timer.
            print("\tAnalyzed", currentID, "in %1.3f seconds." % (time.time() - Start_time))

            # Write Result to File
            num = str(iso + 1)  # Derive a unique identifier.
            while (len(num) != 6):  # While the length of the unique identifier is less than six,
                num = "0" + num  # append a "0" to the beginning of the identifier.
            # output = open("acomysTopTen.peptides.fasta", "a+")		# Top Ten Test Code.
            output = open("acomys_0.peptides.fasta", "a+")  # Open output file,
            # include FASTA header with >"initials"|"unique identifier"|"transcript ID"|"BLASTx reference ID & description"
            output.write(">aro|acomys_D0_%s|%s|%s\n" % (num, currentID, currentRef_description))
            output.write(str(h_seq) + "\n")  # write AA sequence to output file,
            output.close()  # and close the output file.
tSUM = sum(t)  # Calculate total time.
tAVG = tSUM / len(t)  # Calculate average time per transcript.
print("\nAnalyzed", len(IsoformList),
      "transcripts in %1.3f seconds (%1.3f seconds per transcript)\n" % (tSUM, tAVG))