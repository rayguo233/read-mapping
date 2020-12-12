# Description of files
- **report.docx** - the report of this project
- **pj4.py** - the implementation of PatternMatchingWithSuffixArray and BetterBWMatching; instruction to run can be found bellow
- **plot code.pdf** - the code used to plot figures used in the report (the code is written in Jupyter Notebook environment)

# Instruction to Run pj4.py

The script "pj4.py" can be run from shell.  

First edit the first line of this script (`#!/usr/local/cs/bin/python3`) to make sure that it is invoking Python from the correct directory in your environment. 

Then give execution permission using commands like:  
`chmod 744 ./pj4.py`

Then the script can be run with the following syntax.  
`./pj4.py [NUM_READS] [LEN_READ] [LEN_GENOME] [WHICH_ALGORITHM]`
where:
```
NUM_READS - number of reads
LEN_READ - length of each read
LEN_GENOME - length of the genome
WHICH_ALGORITHM - which algorithm to perform:
	'0' for PatternMatchingWithSuffixArray
	'1' for BetterBWMatching
	'2' for checking if the two algorithms's results agree
```
For example, `./pj4.py 100000 100 500000 0` means apply PatternMatchingWithSuffixArray to a genome of length 500,000 with 100,000 reads all of which are of length 100.
