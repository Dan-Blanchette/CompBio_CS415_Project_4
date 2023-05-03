# CompBio_CS415_Project_4
Project 4 for my computational biology and sequence analysis class


## Project 4:

(You may do this project in small groups.)

In Project 3 you were given series of aligned sequences to analyze. In this project we will start, but not complete, the process of building a profile HMM from a new set of sequences. 

The data file project4data-1.txt Download project4data-1.txt  contains 8 sequences generated using a GA similar to the one that generated the sequences from Project 3. (Note: these sequences are variable length and they may wrap oddly in the Canvas interface, so it's safer to download the file to work with). The first step for most progressive alignment algorithms used to build a multiple sequence alignment is to calculate the diagonal matrix of all pairwise alignments (pg. 145-147 of the text). We can build these pairwise alignments using our global alignment algorithm from the beginning of the course (possibly used in Project 2).

We also need a substitution matrix to generate the global alignments. Test two options: the one you generated in Project 3 and the Blosum50 matrix. You will probably need to experiment with the gap penalty - if it's too large gaps will never be inserted, if it's too small gaps will be preferred where misaligning similar symbols/residues/amino acids is more appropriate. Ideally, you will find a value where the gaps are not all placed at one end or the other of the alignments. 

Write-up: This may be a slightly shorter project write-up than the previous one. 

- Abstract: a summary of what you did and the results.
- Methods: a short description of what you did.
- Results: the results section should include two tables showing the pairwise scores, one using your substitution matrix, one using Blosum50. Also, include samples of the global alignments - i.e. examples of two aligned sequences. You may want to include the best (highest scoring) and worst (lowest scoring) sequences.
- Conclusions: a summary of your results. Include a discussion of the following questions:
  - How 'good' do you think the results were, i.e. do the alignments seem plausible?
  - Where the result of from the two substitution matrices the same? If not, what was the difference and which seemed better?
  - Based on the alignment scores do you think the sequences could be clustered into two or more subsets?
