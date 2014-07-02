DEPENDENCIES
---------------------------------------------------------------------
classMarker.py - class definition for autosomal and X chromosome markers

gsl.py - gsl wrapper

USAGE
---------------------------------------------------------------------
python extFBAT.py -fm=<feature matrix path> -phenotype=<phenotypeFile.txt> -gender=<genderFile.txt> [other optional arguments]


INPUT 
-----------------------------------------------------
1. Required: feature matrix (command line option -fm)
Example: sampleFeatureMatrix.txt
Row 1: trio ids
Row 2: trio member type (father = 1, mother = 2, offspring = 3)
Row 3 onwards is genotype data with marker id as first column. 
marker id examples : chr13:328658 or chrX:83769 or chr23:83769 (chrY and chrM are not analyzed)
Genotype coding : allele 1 homozygous = 0, heterozygous = 1, allele2 homozygous = 2, half calls or no calls = NA

2. Required: phenotype file (command line option -pheno)
Example: samplePhenotype.txt
2 column file with trio ids in first column and the offspring's affectation status in second column (1 = control, 2 = case, NA = unknown)
Trios with unknown phenotype are not analyzed.

3. Required: gender file (command line option -gender)
Example: sampleGenderFile.txt
2 column file with trio ids in first column and the offspring's gender in second column (1 = control, 2 = case, NA = unknown). Used for chrX markers.
Trios with unknown gender are not analyzed.

4. Optional: pedigree list file (command line option -pedigree)
Example: samplePedList.txt
list of trio ids to be analyzed. If omitted, all the trios in the feature matrix are analyzed.

5. Optional: test type (command line option -test)
Options: tdt or fbat or omit to run both tests


6. Optional: test version (command line option -version)
Options: ext or std (extended or standard) or omit to run both


7. Optional: genetic model list(command line option -models)
Options: subset of [a,d,r] or omit to run only additive model. 'a' stands for additive model, 'd' for dominant model, 'r' for recessive model


8. Optional: trait offset (command line option -offset)
Any number in [0,1]. If omitted, default value of 0.5 is used. This option is used only for FBAT to account for unaffected families.

9. Optional: output file name (command line option -out)
If omitted, default output file is tdt.out.gz in the current working directory.


OUTPUT
------------------------------------------------------------------------
1. Standard output: The software writes log messages to standard output which can be redirected to a file for future reference.
2. Results file: 'tdt.out.gz' or a user specified file with -out command line option
Columns:

MarkerID: <chr#:position>

MAF: minor allele frequency computed from the data set

n[0/0,0/1,1/1,./.]: count of each genotype aggregated over all samples at a marker

nCompleteInformative[Cases/Controls]_[MaleNB/FemaleNB]: counts of informative trio types with complete genotypes within case trios and control trios. Trios with male and female offspring are counted separately. These counts are used for computing the standard TDT and FBAT statistics.

nCompleteNonInformative[Cases/Controls]: counts of non-informative trio types with complete genotypes within case trios and control trios. These trios do not contribute to the statistic, but are reported for debugging purposes. 

nIncompleteInformative[Cases/Controls]_[MaleNB/FemaleNB]: counts of informative trio types with incomplete genotypes within case trios and control trios. Trios with male and female offspring are counted separately. These counts are used for computing the robustTDT and extFBAT statistics.

nIncompleteNonInformative[Cases/Controls]: counts of non-informative trio types with incomplete genotypes within case trios and control trios. These trios do not contribute to the statistic, but are reported for debugging purposes.

nMIE: count of trios with a Mendelian Inheritance error.

Depending on the test type, version, and genetic models chosen by the user, additional columns are printed to the output file to report scores and p-values. 

TRIO TYPES 
---------------------------------------------------------------------------
Complete Informative Trio Types for autosomal chromosomes: [set(Parent1,Parent2),offspring]
[set(['0','1']),'0']
[set(['0','1']),'1']
[set(['1','1']),'0']
[set(['1','1']),'1']
[set(['1','1']),'2']
[set(['1','2']),'1']
[set(['1','2']),'2']

Complete Informative Trio Types for ChrX for Male offspring: [Father,Mother,Offspring]
['0','1','0']
['0','1','1']
['1','1','0']
['1','1','1']


Complete Informative Trio Types for ChrX for Female offspring: [Father,Mother,Offspring]
['0','1','0']
['0','1','1']
['1','1','1']
['1','1','2']

Incomplete informative trio types  [set(Parent1,Parent2),offspring]
[set(['NA','NA']),'NA']
[set(['NA','NA']),'0']
[set(['NA','NA']),'1']
[set(['NA','NA']),'2']
[set(['NA','0']),'NA']
[set(['NA','0']),'0']
[set(['NA','0']),'1']
[set(['NA','1']),'NA']
[set(['NA','1']),'0']
[set(['NA','1']),'1']
[set(['NA','1']),'2']
[set(['NA','2']),'NA']
[set(['NA','2']),'1']
[set(['NA','2']),'2']
[set(['0','1']),'NA']
[set(['1','1']),'NA']
[set(['1','2']),'NA']

Incomplete informative Trio Type for ChrX for Male offspring: [Father,Mother,Offspring] 
['NA','NA','NA']
['NA','NA','0']
['NA','NA','1']
['NA','1','NA']
['NA','1','0']
['NA','1','1']
['0','NA','NA']
['0','NA','0']
['0','NA','1']
['1','NA','NA']
['1','NA','0']
['1','NA','1']
['0','1','NA']
['1','1','NA']

Incomplete informative Trio Type for ChrX for Female offspring: [Father,Mother,Offspring] 
['NA','NA','NA']
['NA','NA','0']
['NA','NA','1']
['NA','NA','2']
['NA','1','NA']
['NA','1','0']
['NA','1','1']
['NA','1','2']
['0','NA','NA']
['0','NA','0']
['0','NA','1']
['1','NA','NA']
['1','NA','1']
['1','NA','2']
['0','1','NA']
['1','1','NA']


Non-informative mating types for autosomal chromosomes: set(Parent1, Parent2)
set(['0','0'])
set(['0','2'])
set(['2','2'])

Non-informative mother genotypes for chrX
['0','2']  #regardless of father's genotype, if mother is homozygous, the mating type is non-informative for chrX


 
