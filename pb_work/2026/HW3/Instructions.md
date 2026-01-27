# Homework 3

The focus of homework 3 is to practice using both for and while loops to accomplish the same thing, 
and to practice some basic string manipulation.
This assignment uses the same data as HW2, scaffolds_short.fa, it is provided in the HW3 folder 
for your convenience, but it is identical.

Make sure to include headers and follow PEP8 coding conventions.

## Read the sequences
 * Read the sequence file only once, and store the sequences
 * Use the same file opening and line-by-line reading method as in HW2
 * Print the number of sequences and the total number of bases read: e.g.,
_35 sequences with 16934 bases read from scaffolds_short.fa_
## For each sequence in scaffolds_short.fa
 * Loop over all sequences using a while loop (outer loop)
 * For each sequence. Use an inner loop to break the sequence each at each 'A' in the sequence
 * Each segment should either end in 'A' or at the end of the sequence
 * Output for each sequence (see hw3_example.output)
   * For each sequence, report the sequence of each segment in fasta format
   * Add a sequentially incrementing suffix (beginning at _0) to the sequence name indicating 
   which segment it is from
   * Drop the documentation part of the original title line </br>
   * Example of sequence output: this is FastA format with no comment. See the hw3_example_output for details.</br>
_&gt;test0_1</br>
A_

    * After the complete set of sequence segments, print a line indicating how many total segments 
were found, for instance:</br> _5793 subsequences with 16934 bases read from 35 sequences_
    *  If your code is correct, the number of bases in the original sequences should be the same as in the segments
    * Check the spacing in your output to make sure it is exactly the same as the example. PyCharm can
   compare two files and see if there are differences
 * Using a for loop for the inner loop, repeat the above. Do not re-read the sequences.