# jellyphy
This program is unfinished. This program consists of two scripts for generating a distance matrix for phylogenetic tree building programs (e.g. PHYLIP) from large whole genome datasets using Jellyfish software (https://github.com/gmarcais/Jellyfish) to count k-mer (DNA words) frequencies and a custom python script used to compare these k-mer profiles to get a phylogenetic distance matrix. 
There is a shell script (jellyphy.sh) which takes contig files as input and then calls jellyfish software to count k-mer frequencies within each genome and prints this in column format to a newly created folder. Upon completion of this step a second script written in Python is called (jellyphy.py) which takes all newly created kmer count files and generates a phylogenetic distance matrix using the Jensen-Shannon divergence distance. The resulting output file can then be input into most commonly used tree building programs like PHYLIP. 

To run the program:

bash jellyphy.sh /inputDir -o outputfile

Python file can be ran seperately if wanted with the following command:

python jellyphy.py --inputDir /inputDir > outputfile
