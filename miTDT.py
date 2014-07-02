import sys
import string
import gzip
import numpy as np
from classMarker import AutosomalMarker,ChrXMarker

#USAGE:  python miTDT.py 
#-fm=<featurematrix.txt> <required>
#-phenotype=<phenotype.txt> <required>
#-pedigrees=<pedID.txt> OR DEFAULT to all trios in the fm 
#-offset=<a number between 0 and 1> DEFAULT is 0.5 to give equal and opposite weightage to cases and controls (used only for FBAT)
#-test=tdt OR fbat OR omit to run both
#-version=ext OR std OR omit to run both (NOTE: standard score is computed and reported even with extended version)
#-models=a(dditive) OR d(ominant) OR r(ecessive) OR DEFAULT to additive 
#-out=<output file path> or DEFAULT to tdt.out.gz
#-gender=<NB gender file path> <required> (column 1 is pedID of the NB, column 2 is '1' for male, '2' for female)

#INPUT FORMATS
#feature matrix - first row contains trio ids, second row indicates member type: 1 (father), 2(mother), 3(offspring). 
#Cell values set to 0 (ref homozygous), 1 (heterozygous), 2 (non ref homozygous), or NA (missing). First column is row ids. For marker rows, it is <chr:position>
#phenotype file - 2 column file with pedigree ids in first column, affection status in second column (1 = control, 2 = case, NA = unknown)
#pedigree file - pedigree ids to analyze


#ASSUMPTIONS:
#For trios only, multiple offsprings and larger pedigrees are not handled
#NOTE: TDT statistic is always positive, so ref/alt annotation doesn't matter. FBAT computes over or under transmission of the alternate allele(coded as 1)
#NOTE: phenoFile may contain NA too. those pedigrees will not be used in the analysis
#####################################
##GLOBAL VARIABLES and LOOK-UP TABLES
headerColumns = []
pedIDs = []
pedMemberIndices = {}
pedMemberType = {}
pedPhenoDict = {}
pedNBGender = {}
idColumns = []
memberTypeColumns = []

FM_FILENAME=""                    #required
PHENO_FILENAME=""                 #required
PED_FILENAME=""                   #optional
OFFSET = 0.5			  #optional	
TEST=""                           #optional
VERSION=""                        #optional
MODELS=[]                         #optional
OUTPUT_FILENAME = "tdt.out.gz"    #optional
GENDER_FILENAME = ""

fmFile = None
phenoFile = None
pedIDFile = None
outputFile = None
genderFile = None

############################################
#####FUNCTION DEFINITIONS
def readInputArguments():
	global FM_FILENAME
	global PHENO_FILENAME
	global PED_FILENAME
	global OFFSET
	global TEST
	global VERSION
	global MODELS
	global OUTPUT_FILENAME
	global GENDER_FILENAME

        #PARSE COMMAND LINE ARGS. (3 of the above arguments are required)
        assert (len(sys.argv) >=4 ), 'Insufficient number of arguments.'

        while len(sys.argv) > 1:
                # Command line arguments containing '=' are implicitly options.
                thisArg = sys.argv.pop(1)
                if thisArg.find("=")==-1:
                        print 'Unrecognised argument: '+thisArg
                        sys.exit(1)
                else:
                        name,value = thisArg.split("=")  #split name value pairs
			if __debug__: # ...just to see what's going on.
                                print( "{},{}".format( name, value ) )
                        
			name = name.lower().strip("- ")  #strip hyphens and white spaces
                        if name == "fm":
                                FM_FILENAME = value.strip(" ")
                        elif name == "phenotype":
                                PHENO_FILENAME = value.strip(" ")
                        elif name == "pedigrees":
                                PED_FILENAME = value.strip(" ")
			elif name == "offset":
				OFFSET = float(value.strip(" "))
                        elif name == "test":
                                TEST = value.lower().strip(" ")
                        elif name == "version":
                                VERSION = value.lower().strip(" ")
                        elif name == "model":
                                MODELS = value.lower().strip(" ").split(",")
                        elif name == "out":
                                OUTPUT_FILENAME = value.strip(" ")
			elif name == "gender":
				GENDER_FILENAME = value.strip(" ")
                        else:
                                print "unrecognized option:", name
                                sys.exit(1)

        assert (FM_FILENAME <> ""), 'Feature matrix path was not provided'
        assert (PHENO_FILENAME <> ""), 'Phenotype file path was not provided'
	assert (GENDER_FILENAME <> ""), 'Gender file path was not provided'

def createFileObjects():
	global FM_FILENAME
        global PHENO_FILENAME
        global PED_FILENAME
	global VERSION
        global OUTPUT_FILENAME
	global GENDER_FILENAME

	global fmFile
	global phenoFile
	global pedIDFile
	global outputFile
	global genderFile
	
	#INPUT: feature matrix - first row contains pedigree ids, second row indicates member type: 1 (father), 2(mother), 3(newborn). (not assuming that pedIDs are part of sample IDs)
	#Cell values set to 0 (ref homozygous), 1 (heterozygous), 2 (non ref homozygous), or NA (missing)
	if FM_FILENAME.endswith("gz"):
	        fmFile = gzip.open(FM_FILENAME,"r")
	else:
        	fmFile = open(FM_FILENAME,"r")

	#INPUT: phenotype file - 2 column file with pedigree ids in first column, affection status in second column (1 = control, 2 = case)
	phenoFile = open(PHENO_FILENAME,"r")

	#INPUT: pedigree ids to analyze
	if PED_FILENAME <> "":
       		pedIDFile = open(PED_FILENAME,"r")

	#INPUT: GENDER FILE
	if GENDER_FILENAME <> "":
                genderFile = open(GENDER_FILENAME,"r")

	#OUTPUT: output file
	if not OUTPUT_FILENAME.endswith(".gz"):
	        OUTPUT_FILENAME = OUTPUT_FILENAME+".gz"

	outputFile = gzip.open(OUTPUT_FILENAME,"w")

	

def getPedIDs():
        global pedIDFile
	global pedIDs
	global idColumns

	#READ pedigree ids to analyze. If no pedIDFile has been provided, pedIDs vector contains all the pedigree ids from the header line of the feature matrix
        if pedIDFile:
                pedIDs = pedIDFile.read().split()
        else:
                pedIDs = list(set(idColumns))
	#NOTE:pedIDs does not retain the order of predigrees in the input feature matrix
        #DEBUG
#       print pedIDs,len(pedIDs)


def createPedPhenoDict():
	global phenoFile
	global pedIDs
	global pedPhenoDict

        #get these pedigrees from phenotype file into a dictionary
        for line in phenoFile:
                columns = line.strip().split()
                if columns[0] in pedIDs: 
                        if columns[1].isdigit():
				pedPhenoDict[columns[0]] = int(columns[1])
			else:
				print "Missing phenotype for ",columns[0],". This pedigree will not be analysed."
				pedIDs.remove(columns[0])

def getPedNBGender():
	global genderFile
        global pedIDs
	global pedNBGender
	
        #get these pedigrees from phenotype file into a dictionary
        for line in genderFile:
                columns = line.strip().split()
                if columns[0] in pedIDs:
                        if columns[1].isdigit():
                                pedNBGender[columns[0]] = int(columns[1])
                        else:
                                print "Missing gender for ",columns[0],". This pedigree will not be analysed for chrX."
                                

def getPedMemberIndicesAndType():
	global pedIDs
	global idColumns
	global pedMemberIndices
	global pedMemberType
        for thisPed in pedIDs:
		#DEBUG
#		print thisPed
#		print idColumns
                memberIndices = [x for x in range(len(idColumns)) if idColumns[x]==thisPed]
		pedMemberIndices[thisPed]=memberIndices
                pedMemberType[thisPed]=[memberTypeColumns[x] for x in memberIndices]

#	        #DEBUG
#	        print pedMemberIndices
#	        print (idColumns[x] for x in memberIndices)
#	        print pedMemberType
#	        raw_input('continue')

#	print len(pedMemberIndices),len(pedMemberType)
		
###################################################
######### PROCESSING STARTS HERE ##################

#parse input arguments
readInputArguments()

###open input files for reading and output files for writing
createFileObjects()

#PRINT HEADER line according to options selected by user
outputColumns = ['MarkerID','MAF','n[0/0,0/1,1/1,./.]','nCompleteInformativeCases_MaleNB','nCompleteInformativeCases_FemaleNB','nCompleteInformativeControls_MaleNB','nCompleteInformativeControls_FemaleNB','nCompleteNonInformativeCases','nCompleteNonInformativeControls','nIncompleteInformativeCases_MaleNB','nIncompleteInformativeCases_FemaleNB','nIncompleteInformativeControls_MaleNB','nIncompleteInformativeControls_FemaleNB','nIncompleteNonInformativeCases','nIncompleteNonInformativeControls','nMIE']
if TEST == "" or TEST == "tdt": #both or TDT
	if not MODELS or "a" in MODELS:
		outputColumns.extend(['ChiSq_TDT_Additive','P-value_TDT_Additive'])
		if VERSION == "ext" or VERSION == "":
			outputColumns.extend(['min_ChiSq_rTDT_Additive','min_P-value_rTDT_Additive','max_ChiSq_rTDT_Additive','max_P-value_rTDT_Additive'])
	if "d" in MODELS:    
		outputColumns.extend(['ChiSq_TDT_Dominant','P-value_TDT_Dominant'])
		if VERSION == "ext" or VERSION == "":
                        outputColumns.extend(['min_ChiSq_rTDT_Dominant','min_P-value_rTDT_Dominant','max_ChiSq_rTDT_Dominant','max_P-value_rTDT_Dominant'])	
	if "r" in MODELS:
		outputColumns.extend(['ChiSq_TDT_Recessive','P-value_TDT_Recessive'])
		if VERSION == "ext" or VERSION == "":	
			outputColumns.extend(['min_ChiSq_rTDT_Recessive','min_P-value_rTDT_Recessive','max_ChiSq_rTDT_Recessive','max_P-value_rTDT_Recessive'])	
	if MODELS and "a" not in MODELS and "d" not in MODELS and "r" not in MODELS: 
		print "Unrecognised model: "+MODELS
       		sys.exit(1)


if TEST == "" or TEST == "fbat": #both or FBAT
	if not MODELS or "a" in MODELS:
                outputColumns.extend(['Z_FBAT_Additive','P-value_FBAT_Additive'])
		if VERSION == "ext" or VERSION == "":
                        outputColumns.extend(['min_Z_extFBAT_Additive','min_P-value_extFBAT_Additive','max_Z_extFBAT_Additive','max_P-value_extFBAT_Additive'])
        if "d" in MODELS:
                outputColumns.extend(['Z_FBAT_Dominant','P-value_FBAT_Dominant'])
		if VERSION == "ext" or VERSION == "":
                        outputColumns.extend(['min_Z_extFBAT_Dominant','min_P-value_extFBAT_Dominant','max_Z_extFBAT_Dominant','max_P-value_extFBAT_Dominant'])
        if "r" in MODELS:
                outputColumns.extend(['Z_FBAT_Recessive','P-value_FBAT_Recessive'])
		if VERSION == "ext" or VERSION == "":
                        outputColumns.extend(['min_Z_extFBAT_Recessive','min_P-value_extFBAT_Recessive','max_Z_extFBAT_Recessive','max_P-value_extFBAT_Recessive'])
        if not MODELS and "a" not in MODELS and "d" not in MODELS and "r" not in MODELS:
                print "Unrecognised model: "+MODELS
                sys.exit(1)	

	
outputFile.write('\t'.join(outputColumns)+'\n')

	
#READ header line containing pedIDs and 2nd row containing member type of each sample
idColumns = fmFile.readline().strip().split('\t')[1:]  #1st row of the feature matrix - pedigree ids (skip 1st column)
memberTypeColumns = fmFile.readline().strip().split('\t')[1:]  #second row of the feature matrix (skip 1st column)

#assert that all membertype assignments are numeric - 1/2/3 for F/M/NB
assert(all(v.isdigit() for v in memberTypeColumns[1:len(memberTypeColumns)]))   

#get predigree ids 
getPedIDs()
#generate pedigree-phenotype dictionary
createPedPhenoDict()
#pedigree-NBgender dictionary
getPedNBGender()
#get pedigree member indices and member types from 1st two lines of the feature matrix
getPedMemberIndicesAndType()

## ALL GLOBAL CHECKS and ASSERTS HERE
## TDT must be provided with only case trios, FBAT must have both case and controls
if TEST == "tdt" and 1 in pedPhenoDict.values():
        print 'Pedigree file must contain only case pedigrees for TDT.'
        sys.exit(1)
if TEST == "fbat" and 1 not in pedPhenoDict.values():
        print 'Pedigree file must contain control pedigrees for FBAT.'
        sys.exit(1)



#read feature matrix one line at a time. First two lines have been read above for pedigree ids and member type
for line in fmFile:
	#create new marker object
	#chrM and chrY are not tested
	if line.startswith('chrM') or line.startswith('chr25') or line.startswith('chrY') or line.startswith('chr24'):
		continue	
	elif line.startswith('chrX') or line.startswith('chr23'): 
		thisMarker = ChrXMarker()
	else:  #TODO: more thorough check for validity of data format, like chr numbers??
		thisMarker = AutosomalMarker()
	
	#get vcf columns
 	vcfValues = line.strip().split('\t')

	#assert that all genotype values are numeric #TODO: verify that this assert works
	assert(all(v.isdigit() or v=="NA" for v in vcfValues[1:]))   #1st column is variant id, 2nd onwards are sample genotypes

	
	#set this marker object's sample values 
	thisMarker.markerID = vcfValues[0]
	thisMarker.getSampleGenotypes(vcfValues[1:],pedMemberIndices)
	
	#Note: thought of checking if chrX marker has heterozygous males, but it's not possible since the feature matrix comes in with encoded genotypes. 
	#So all you can check is whether autosomal chrs have any genotypes other than 0/1/2/NA and chrX has any genotypes other that 0/1/NA for males, 0/1/2/NA for females
	#for chrX
	#TODO: verify that hasValidGenotypes() works
	if isinstance(thisMarker,ChrXMarker) and not thisMarker.hasValidGenotypes(pedNBGender,pedMemberType):
		print 'Invalid genotype found at ',thisMarker.markerID,'. This marker will not be tested.'
		continue
	#for autosomal chromosomes
	elif not isinstance(thisMarker,ChrXMarker) and not thisMarker.hasValidGenotypes():
		print 'Invalid genotype found at ',thisMarker.markerID,'. This marker will not be tested.'
                continue
		
	#COMPUTE ALLELE FREQUENCY
	#compute regardless of whether 'mi' has been selected or not, because allele frequency will be reported in the output
	if isinstance(thisMarker,ChrXMarker):
		thisMarker.computeMAF(pedNBGender,pedMemberType)
	else:
		thisMarker.computeMAF()
	
	try:
		#assert that MAF is always positive
		assert(thisMarker.maf >= 0)
	except(AssertionError):
		print 'Marker ',thisMarker.markerID,', MAF=',thisMarker.maf
		exit(1)	

        #DEBUG
#	print 'Ref and Alt Frequencies'
#       print vcfValues[0],thisMarker.maf
#       raw_input('continue')
	
	#get variant distribution
	if isinstance(thisMarker,ChrXMarker):	
		thisMarker.getVariantDistribution(pedNBGender,pedMemberType)
	else:
		thisMarker.getVariantDistribution()

	#count complete and incomplete case and control trio types and populate corresponding vectors
	thisMarker.populateTrioTypeCountVectors(pedMemberType,pedPhenoDict,pedNBGender)
	
	#concatenate output string and print to output file
        outputColumns = [thisMarker.markerID,str(thisMarker.maf),str(thisMarker.nVariantType),str(thisMarker.nCompleteInformativeCaseTrio_MaleNB),str(thisMarker.nCompleteInformativeCaseTrio_FemaleNB),str(thisMarker.nCompleteInformativeControlTrio_MaleNB),str(thisMarker.nCompleteInformativeControlTrio_FemaleNB),str(thisMarker.nCompleteNonInformativeCaseTrio),str(thisMarker.nCompleteNonInformativeControlTrio),str(thisMarker.nIncompleteInformativeCaseTrio_MaleNB),str(thisMarker.nIncompleteInformativeCaseTrio_FemaleNB),str(thisMarker.nIncompleteInformativeControlTrio_MaleNB),str(thisMarker.nIncompleteInformativeControlTrio_FemaleNB),str(thisMarker.nIncompleteNonInformativeCaseTrio),str(thisMarker.nIncompleteNonInformativeControlTrio),str(thisMarker.nMIE)]			
	# run appropriate tests based on options selected by the user, add appropriate output columns	
	#**************TDT***************************************************************************
	if TEST == "" or TEST == "tdt":
		#ADDITIVE std. TDT-------------------------------------------------------------------
		if not MODELS or "a" in MODELS:
			thisMarker.stdTDT("a")
			outputColumns.extend([str(thisMarker.chiSq_StdTDT),str(thisMarker.pValue_StdTDT)])
			#ADDITIVE mi-TDT--------------------------------------
			if VERSION == "ext" or VERSION == "":
				thisMarker.extendedTDT("a")
				outputColumns.extend([str(thisMarker.minChiSq_rTDT),str(thisMarker.minPValue_rTDT),str(thisMarker.maxChiSq_rTDT),str(thisMarker.maxPValue_rTDT)])

		#DOMINANT std. TDT----------------------------------------------------------------------------------
		if "d" in MODELS:   
			thisMarker.stdTDT("d")
			outputColumns.extend([str(thisMarker.chiSq_StdTDT),str(thisMarker.pValue_StdTDT)])
			#DOMINANT mi-TDT----------------------------------------------------------------------------
			if VERSION == "ext" or VERSION == "":
				thisMarker.extendedTDT("d")
				outputColumns.extend([str(thisMarker.minChiSq_rTDT),str(thisMarker.minPValue_rTDT),str(thisMarker.maxChiSq_rTDT),str(thisMarker.maxPValue_rTDT)])
		
		#RECESSIVE std. TDT----------------------------------------------------------------------------------
                if "r" in MODELS:   
                        thisMarker.stdTDT("r")
			outputColumns.extend([str(thisMarker.chiSq_StdTDT),str(thisMarker.pValue_StdTDT)])
                        #RECESSIVE mi-TDT----------------------------------------------------------------------------
			if VERSION == "ext" or VERSION == "":	
                                thisMarker.extendedTDT("r")
                                outputColumns.extend([str(thisMarker.minChiSq_rTDT),str(thisMarker.minPValue_rTDT),str(thisMarker.maxChiSq_rTDT),str(thisMarker.maxPValue_rTDT)])
	
	if TEST == "" or TEST == "fbat":
		#ADDITIVE std. FBAT-------------------------------------------------------------------
                if not MODELS or "a" in MODELS:
			thisMarker.stdFBAT("a",OFFSET)
                        outputColumns.extend([str(thisMarker.Z_stdFBAT),str(thisMarker.pValue_stdFBAT)])
                        #ADDITIVE ext-FBAT------------------------------
			if VERSION == "ext" or VERSION == "":
				thisMarker.extendedFBAT("a",OFFSET)
				outputColumns.extend([str(thisMarker.minZ_extFBAT),str(thisMarker.minPValue_extFBAT),str(thisMarker.maxZ_extFBAT),str(thisMarker.maxPValue_extFBAT)])
		
		#DOMINANT std. FBAT-------------------------------------------------------------------
                if "d" in MODELS:
                        thisMarker.stdFBAT("d",OFFSET)
                        outputColumns.extend([str(thisMarker.Z_stdFBAT),str(thisMarker.pValue_stdFBAT)])
                        #DOMINANT ext-FBAT------------------------------
			if VERSION == "ext" or VERSION == "":
                                thisMarker.extendedFBAT("d",OFFSET)
                                outputColumns.extend([str(thisMarker.minZ_extFBAT),str(thisMarker.minPValue_extFBAT),str(thisMarker.maxZ_extFBAT),str(thisMarker.maxPValue_extFBAT)])	

		#RECESSIVE std. FBAT-------------------------------------------------------------------
                if "r" in MODELS:
                        thisMarker.stdFBAT("r",OFFSET)
                        outputColumns.extend([str(thisMarker.Z_stdFBAT),str(thisMarker.pValue_stdFBAT)])
                        #RECESSIVE ext-FBAT-----------------------------
			if VERSION == "ext" or VERSION == "":
                                thisMarker.extendedFBAT("r",OFFSET)
                                outputColumns.extend([str(thisMarker.minZ_extFBAT),str(thisMarker.minPValue_extFBAT),str(thisMarker.maxZ_extFBAT),str(thisMarker.maxPValue_extFBAT)])


	#print to outputfile
	outputFile.write('\t'.join(outputColumns)+'\n')

#TODO: close all files
fmFile.close()
phenoFile.close()
if pedIDFile:
	pedIDFile.close()
outputFile.close()

