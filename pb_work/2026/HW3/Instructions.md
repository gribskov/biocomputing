# Homework 3

The focus of homework 3 is to practice using both for and while loops to accomplish the same thing, 
and to practice some basic string manipulation.
This assignment uses the same data as HW2, scaffolds_short.fa, it is provided in the HW3 folder 
for your convenience, but it is identical.

Make sure to include headers and follow PEP8 coding conventions.

## For each sequence in scaffolds_short.fa

 * Read the sequence file only once, and store the sequences
 * Use the same file opening and line-by-line reading method as in HW2
 * Using a for loop, break each sequence into the segments that begin with 'A'.
   * Each segment should contain the initial A so the minimum length will be 1 base
   * Output (see hw3_example.output)
     * Print the number of sequences read: e.g.,35 sequences read from scaffolds_short.fa
     * For each sequence, report the sequence of each segment in fasta format
       * add a sequentially incrementing suffix (beginning at _1) to the sequence name indicating which segment it is
       from, for instance </br>
_&gt;test0_1</br>
A_
        * Drop the documentation part of the original title line </br>
        * See the hw3_example_output for details.
      * After the complete set of sequence segments, print a line indicating how many total segments 
were found, for instance:</br> 5793 subsequences read from 35 sequences
      * Check the spacing in your output to make 
sure it is the same as the example
 * Using a while loop, do the same as above