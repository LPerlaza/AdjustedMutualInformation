# AdjustedMutualInformation

AdjMI.c is a program that calculates the AMI for each pair of SNP sites in an alignment. Takes two SNP sites from an alignment per time and calculates the Adjusted Mutual Information. It is optimised for sides that are biallelic (only two variations in each site). It has into account Ns, and uses only the intersection of the two sites that have no Ns. Only will have into account SNPs that have a minor allele frequency of 5%, thus SNP sites with really low variability will not be included. This program is also optimised to be run in parallel and use shared memory to read the input file. This two features will help you to improve the speed of your calculations and reduce the computational resources used for calculating AMI.

![AdjMI] (https://github.com/LPerlaza/AdjustedMutualInformation/blob/master/AdjMI_updated.pdf)

1. Download and install the MI dependency from here: https://github.com/Craigacp/MIToolbox/

2. C code has to be compiled. Please compile using this command:

```gcc -O3 AdjMI.c -lMIToolbox -lm -o AdjMI.exe```

3. You have to run the program as:

```./AdjMI.exe <SNPs_InputFile> <Intervals_File> <Threshold>```
 

The *SNPs_InputFile*: The input file has to be formatted like the attached example (Example_inputFile.txt). The first line should have a header for the chromosome (CHROM) and Position (POS). and then a number for each accession separate by commas. The second line will have the corresponding chromosome for the first SNP and the corresponding position, and then the nucleotide corresponding to each accession in the same order that is in the first line. This time NO commas. and so on. I know this format is not conventional but facilitates the reading of big files. You should be able to transform VCF files into this format using R (if your tables are not too big) or linux commands to transpose. Let me know if you need some help with this.

The *intervals_File*:  The intervals-file is intended to guide AdjMI.exe to spread equitably the number of comparisons along the different parallel jobs (cores). The function behind Intervals-file divides efficiently the number of computations. AdjMI.exe make pairwise comparisons visiting each line of the SNPs-file, thus considering that at the beginning of the file, the number of comparisons to perform is larger than towards the end of the file. Then, Intervals-file is a file with the number of the lines each parallel job needs to perform to spend an equitable amount of time per each processors, helping to make the calculations as fast as possible. You can use Intervals_fileGenerator.r for this or use your own way to divide the load equitably. 

The *Threshold* is the minimum number of accession of the intersection you want to consider to have unambiguous variant calls per pairwise comparison. Imagen you have 1000 genomes in your alignment, when doing the intersection of the nonambiguous variants what is the minimun number the vector should be?.

4. If you manage to do all previous steps you will have the output files which are named as:  inputfile_interval.start_interval.end.txt . Due to the large number of comparisons that are needed and splitting the processing into a number of cores you are going to have several output files that you can concatenate as you want, or process separately, as you prefer. The format is simple:

EMI: Expected Mutual Information

MI: Mutual Information

AMI: Adjusted Mutual Information

NMI: Normalised Mutual Information

EntropyA: Entropy of SNP_A

EntropyB: Entropy of SNP_B

CHR_A: Chromosome location of SNP_A

POS_A: Position in the chromosome of SNP_A

FreqA1: Frequency of the allele_1 in SNP_A

FreqA2: Frequency of the allele_2 in SNP_A

CHR_B: Chromosome location of SNP_B

POS_B: Position in the chromosome of SNP_B

FreqB1: Frequency of the allele_1 in SNP_B

FreqB2: Frequency of the allele_2 in SNP_B

Length: Length of the vectors intersection (both have NO Ns)

AccCode: The order of the accessions numbers that were included in the AMI calculation


