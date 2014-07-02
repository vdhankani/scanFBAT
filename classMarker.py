import sys
import string
import math
from scipy import stats

##LOOK-UP TABLES
#Informative mating types
#Mating type 2: set(1,2)
#Mating type 4: set(1,1)
#Mating type 5: set(0,1)
#trio type notation: [X_50, X_51, X_40, X_41, X_42, X_21, X_22]
informativeTrioType = [[set(['0','1']),'0'],[set(['0','1']),'1'],[set(['1','1']),'0'],[set(['1','1']),'1'],[set(['1','1']),'2'],[set(['1','2']),'1'],[set(['1','2']),'2']]
#separate for chrX and for male and female NB. Note: in chrX, mating type is not a set but an ordered vector [F,M,NB] because a mating type is informative only if M is heterozygous. 
#trio type notation: [X_50, X_51, X_40, X_41]
informativeTrioType_ChrX_MaleNB = [['0','1','0'],['0','1','1'],['1','1','0'],['1','1','1']]
#trio type notation: [X_50, X_51, X_41, X_42]
informativeTrioType_ChrX_FemaleNB = [['0','1','0'],['0','1','1'],['1','1','1'],['1','1','2']]

#Non-informative mating types
#for autosomal chromosomes
nonInformativeMatingType = [set(['0','0']),set(['0','2']),set(['2','2'])]
#for ChrX
nonInformativeMotherType_ChrX = ['0','2']  #regardless of F's genotype, if M is homozygous, the mating type is non-informative

#stdTDT look-up tables
informativeTrioBinc_Additive = [1,0,2,1,0,1,0]
informativeTrioCinc_Additive = [0,1,0,1,2,0,1]
informativeTrioBinc_Dominant = [1,0,1,0,0,0,0]
informativeTrioCinc_Dominant = [0,1,0,1,1,0,0]
informativeTrioBinc_Recessive = [0,0,1,1,0,1,0]
informativeTrioCinc_Recessive = [0,0,0,0,1,0,1]

informativeTrioBinc_AllModels_ChrX_MaleNB = [1,0,1,0]
informativeTrioCinc_AllModels_ChrX_MaleNB = [0,1,0,1]

informativeTrioBinc_Additive_ChrX_FemaleNB = [1,0,1,0]
informativeTrioCinc_Additive_ChrX_FemaleNB = [0,1,0,1]
informativeTrioBinc_Dominant_ChrX_FemaleNB = [1,0,0,0]
informativeTrioCinc_Dominant_ChrX_FemaleNB = [0,1,0,0]
informativeTrioBinc_Recessive_ChrX_FemaleNB = [0,0,1,0] 
informativeTrioCinc_Recessive_ChrX_FemaleNB = [0,0,0,1]
 
#stdFBAT look-up tables
informativeTrioU_Additive = [-0.5,0.5,-1,0,1,-0.5,0.5]
informativeTrioU_Dominant = [-0.5,0.5,-2/float(3),1/float(3),1/float(3),0,0]
informativeTrioU_Recessive = [0,0,-1/float(3),2/float(3),2/float(3),-0.5,0.5]
informativeTrioVarU_Additive = [0.25,0.25,0.5,0.5,0.5,0.25,0.25]
informativeTrioVarU_Dominant = [0.25,0.25,2/float(9),2/float(9),2/float(9),0,0]
informativeTrioVarU_Recessive = [0,0,2/float(9),2/float(9),2/float(9),0.25,0.25]

informativeTrioU_AllModels_ChrX_MaleNB = [-0.5,0.5,-0.5,0.5]
informativeTrioVarU_AllModels_ChrX_MaleNB = [0.25,0.25,0.25,0.25]

informativeTrioU_Additive_ChrX_FemaleNB = [-0.5,0.5,-0.5,0.5]
informativeTrioVarU_Additive_ChrX_FemaleNB = [0.25,0.25,0.25,0.25]
informativeTrioU_Dominant_ChrX_FemaleNB = [-0.5,0.5,0,0]
informativeTrioVarU_Dominant_ChrX_FemaleNB = [0.25,0.25,0,0]
informativeTrioU_Recessive_ChrX_FemaleNB = [0,0,-0.5,0.5]
informativeTrioVarU_Recessive_ChrX_FemaleNB = [0,0,0.25,0.25]


#Incomplete trio look-up table  [set('F','M'),'NB']
incompleteTrioType = [[set(['NA','NA']),'NA'],[set(['NA','NA']),'0'],[set(['NA','NA']),'1'],[set(['NA','NA']),'2'],[set(['NA','0']),'NA'],[set(['NA','0']),'0'],[set(['NA','0']),'1'],[set(['NA','1']),'NA'],[set(['NA','1']),'0'],[set(['NA','1']),'1'],[set(['NA','1']),'2'],[set(['NA','2']),'NA'],[set(['NA','2']),'1'],[set(['NA','2']),'2'],[set(['0','1']),'NA'],[set(['1','1']),'NA'],[set(['1','2']),'NA']]

incompleteTrioType_ChrX_MaleNB = [['NA','NA','NA'],['NA','NA','0'],['NA','NA','1'],['NA','1','NA'],['NA','1','0'],['NA','1','1'],['0','NA','NA'],['0','NA','0'],['0','NA','1'],['1','NA','NA'],['1','NA','0'],['1','NA','1'],['0','1','NA'],['1','1','NA']]

incompleteTrioType_ChrX_FemaleNB = [['NA','NA','NA'],['NA','NA','0'],['NA','NA','1'],['NA','NA','2'],['NA','1','NA'],['NA','1','0'],['NA','1','1'],['NA','1','2'],['0','NA','NA'],['0','NA','0'],['0','NA','1'],['1','NA','NA'],['1','NA','1'],['1','NA','2'],['0','1','NA'],['1','1','NA']]

#rTDT look-up tables
incompleteTrioBinc_Additive = [2,2,1,0,1,1,0,2,2,1,0,1,1,0,1,2,1]
incompleteTrioCinc_Additive = [2,0,1,2,1,0,1,2,0,1,2,1,0,1,1,2,1]
incompleteTrioBinc_Dominant = [1,1,0,0,1,1,0,1,1,0,0,0,0,0,1,1,0]
incompleteTrioCinc_Dominant = [1,0,1,1,1,0,1,1,0,1,1,0,0,0,1,1,0]
incompleteTrioBinc_Recessive = [1,1,1,0,0,0,0,1,1,1,0,1,1,0,0,1,1]
incompleteTrioCinc_Recessive = [1,0,0,1,0,0,0,1,0,0,1,1,0,1,0,1,1]

incompleteTrioBinc_AllModels_ChrX_MaleNB = [1,1,0,1,1,0,1,1,0,1,1,0,1,1]
incompleteTrioCinc_AllModels_ChrX_MaleNB = [1,0,1,1,0,1,1,0,1,1,0,1,1,1]

incompleteTrioBinc_Additive_ChrX_FemaleNB = [1,1,1,0,1,1,1,0,1,1,0,1,1,0,1,1]
incompleteTrioCinc_Additive_ChrX_FemaleNB = [1,0,1,1,1,0,1,1,1,0,1,1,0,1,1,1]
incompleteTrioBinc_Dominant_ChrX_FemaleNB = [1,1,0,0,1,1,0,0,1,1,0,0,0,0,1,0]
incompleteTrioCinc_Dominant_ChrX_FemaleNB = [1,0,1,0,1,0,1,0,1,0,1,0,0,0,1,0]
incompleteTrioBinc_Recessive_ChrX_FemaleNB = [1,0,1,0,1,0,1,0,0,0,0,1,1,0,0,1]
incompleteTrioCinc_Recessive_ChrX_FemaleNB = [1,0,0,1,1,0,0,1,0,0,0,1,0,1,0,1]

#rFBAT look-up tables
incompleteCaseTrioUmin_Additive = [-1,-1,-0.5,0,-0.5,-0.5,0,-1,-1,-0.5,0.5,-0.5,-0.5,0,-0.5,-1,-0.5]
incompleteCaseTrioVarUmin_Additive = [0.5,0.5,0.25,0,0.25,0.25,0,0.5,0.5,0.25,0.25,0.25,0.25,0,0.25,0.5,0.25]
incompleteControlTrioUmin_Additive = [1,0,0.5,1,0.5,0,0.5,1,-0.5,0.5,1,0.5,0,0.5,0.5,1,0.5]
incompleteControlTrioVarUmin_Additive = [0.5,0,0.25,0.5,0.25,0,0.25,0.5,0.25,0.25,0.5,0.25,0,0.25,0.25,0.5,0.25]
incompleteCaseTrioUmax_Additive = [1,0,0.5,1,0.5,0,0.5,1,-0.5,0.5,1,0.5,0,0.5,0.5,1,0.5]
incompleteCaseTrioVarUmax_Additive = [0.5,0,0.25,0.5,0.25,0,0.25,0.5,0.25,0.25,0.5,0.25,0,0.25,0.25,0.5,0.25]
incompleteControlTrioUmax_Additive = [-1,-1,-0.5,0,-0.5,-0.5,0,-1,-1,-0.5,0.5,-0.5,-0.5,0,-0.5,-1,-0.5]
incompleteControlTrioVarUmax_Additive = [0.5,0.5,0.25,0,0.25,0.25,0,0.5,0.5,0.25,0.25,0.25,0.25,0,0.25,0.5,0.25]

incompleteCaseTrioUmin_Dominant = [-2/float(3),-2/float(3),0,0,-0.5,-0.5,0,-2/float(3),-2/float(3),0,0,0,0,0,-0.5,-2/float(3),0]
incompleteCaseTrioVarUmin_Dominant = [2/float(9),2/float(9),0,0,0.25,0.25,0,2/float(9),2/float(9),0,0,0,0,0,0.25,2/float(9),0]
incompleteControlTrioUmin_Dominant = [0.5,0,0.5,1/float(3),0.5,0,0.5,0.5,-0.5,0.5,1/float(3),0,0,0,0.5,1/float(3),0]
incompleteControlTrioVarUmin_Dominant = [0.25,0,0.25,2/float(9),0.25,0,0.25,0.25,0.25,0.25,2/float(9),0,0,0,0.25,2/float(9),0]
incompleteCaseTrioUmax_Dominant = [0.5,0,0.5,1/float(3),0.5,0,0.5,0.5,-0.5,0.5,1/float(3),0,0,0,0.5,1/float(3),0]
incompleteCaseTrioVarUmax_Dominant = [0.25,0,0.25,2/float(9),0.25,0,0.25,0.25,0.25,0.25,2/float(9),0,0,0,0.25,2/float(9),0]
incompleteControlTrioUmax_Dominant = [-2/float(3),-2/float(3),0,0,-0.5,-0.5,0,-2/float(3),-2/float(3),0,0,0,0,0,-0.5,-2/float(3),0]
incompleteControlTrioVarUmax_Dominant = [2/float(9),2/float(9),0,0,0.25,0.25,0,2/float(9),2/float(9),0,0,0,0,0,0.25,2/float(9),0]

incompleteCaseTrioUmin_Recessive = [-0.5,-1/float(3),-0.5,0,0,0,0,-0.5,-1/float(3),-0.5,0.5,-0.5,-0.5,0,0,-1/float(3),-0.5]
incompleteCaseTrioVarUmin_Recessive = [0.25,2/float(9),0.25,0,0,0,0,0.25,2/float(9),0.25,0.25,0.25,0.25,0,0,2/float(9),0.25]
incompleteControlTrioUmin_Recessive = [2/float(3),0,2/float(3),2/float(3),0,0,0,2/float(3),0,2/float(3),2/float(3),0.5,0,0.5,0,2/float(3),0.5]
incompleteControlTrioVarUmin_Recessive = [2/float(9),0,2/float(9),2/float(9),0,0,0,2/float(9),0,2/float(9),2/float(9),0.25,0,0.25,0,2/float(9),0.25]
incompleteCaseTrioUmax_Recessive = [2/float(3),0,2/float(3),2/float(3),0,0,0,2/float(3),0,2/float(3),2/float(3),0.5,0,0.5,0,2/float(3),0.5]
incompleteCaseTrioVarUmax_Recessive = [2/float(9),0,2/float(9),2/float(9),0,0,0,2/float(9),0,2/float(9),2/float(9),0.25,0,0.25,0,2/float(9),0.25]
incompleteControlTrioUmax_Recessive = [-0.5,-1/float(3),-0.5,0,0,0,0,-0.5,-1/float(3),-0.5,0.5,-0.5,-0.5,0,0,-1/float(3),-0.5]
incompleteControlTrioVarUmax_Recessive = [0.25,2/float(9),0.25,0,0,0,0,0.25,2/float(9),0.25,0.25,0.25,0.25,0,0,2/float(9),0.25]

incompleteCaseTrioUmin_AllModels_ChrX_MaleNB = [-0.5,-0.5,0,-0.5,-0.5,0.5,-0.5,-0.5,0,-0.5,-0.5,0,-0.5,-0.5]
incompleteCaseTrioVarUmin_AllModels_ChrX_MaleNB = [0.25,0.25,0,0.25,0.25,0.25,0.25,0.25,0,0.25,0.25,0,0.25,0.25]
incompleteControlTrioUmin_AllModels_ChrX_MaleNB = [0.5,0,0.5,0.5,-0.5,0.5,0.5,0,0.5,0.5,0,0.5,0.5,0.5]
incompleteControlTrioVarUmin_AllModels_ChrX_MaleNB = [0.25,0,0.25,0.25,0.25,0.25,0.25,0,0.25,0.25,0,0.25,0.25,0.25]
incompleteCaseTrioUmax_AllModels_ChrX_MaleNB = [0.5,0,0.5,0.5,-0.5,0.5,0.5,0,0.5,0.5,0,0.5,0.5,0.5]
incompleteCaseTrioVarUmax_AllModels_ChrX_MaleNB = [0.25,0,0.25,0.25,0.25,0.25,0.25,0,0.25,0.25,0,0.25,0.25,0.25]
incompleteControlTrioUmax_AllModels_ChrX_MaleNB = [-0.5,-0.5,0,-0.5,-0.5,0.5,-0.5,-0.5,0,-0.5,-0.5,0,-0.5,-0.5]
incompleteControlTrioVarUmax_AllModels_ChrX_MaleNB = [0.25,0.25,0,0.25,0.25,0.25,0.25,0.25,0,0.25,0.25,0,0.25,0.25]

incompleteCaseTrioUmin_Additive_ChrX_FemaleNB = [-0.5,-0.5,-0.5,0,-0.5,-0.5,-0.5,0.5,-0.5,-0.5,0,-0.5,-0.5,0,-0.5,-0.5]
incompleteCaseTrioVarUmin_Additive_ChrX_FemaleNB = [0.25,0.25,0.25,0,0.25,0.25,0.25,0.25,0.25,0.25,0,0.25,0.25,0,0.25,0.25]
incompleteControlTrioUmin_Additive_ChrX_FemaleNB = [0.5,0,0.5,0.5,0.5,-0.5,0.5,0.5,0.5,0,0.5,0.5,0,0.5,0.5,0.5]
incompleteControlTrioVarUmin_Additive_ChrX_FemaleNB = [0.25,0,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.25,0.25,0,0.25,0.25,0.25]
incompleteCaseTrioUmax_Additive_ChrX_FemaleNB = [0.5,0,0.5,0.5,0.5,-0.5,0.5,0.5,0.5,0,0.5,0.5,0,0.5,0.5,0.5]
incompleteCaseTrioVarUmax_Additive_ChrX_FemaleNB = [0.25,0,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0,0.25,0.25,0,0.25,0.25,0.25]
incompleteControlTrioUmax_Additive_ChrX_FemaleNB = [-0.5,-0.5,-0.5,0,-0.5,-0.5,-0.5,0.5,-0.5,-0.5,0,-0.5,-0.5,0,-0.5,-0.5]
incompleteControlTrioVarUmax_Additive_ChrX_FemaleNB = [0.25,0.25,0.25,0,0.25,0.25,0.25,0.25,0.25,0.25,0,0.25,0.25,0,0.25,0.25]

incompleteCaseTrioUmin_Dominant_ChrX_FemaleNB = [-0.5,-0.5,0,0,-0.5,-0.5,0,0,-0.5,-0.5,0,0,0,0,-0.5,0]
incompleteCaseTrioVarUmin_Dominant_ChrX_FemaleNB = [0.25,0.25,0,0,0.25,0.25,0,0,0.25,0.25,0,0,0,0,0.25,0]
incompleteControlTrioUmin_Dominant_ChrX_FemaleNB = [0.5,0,0.5,0,0.5,-0.5,0.5,0,0.5,0,0.5,0,0,0,0.5,0]
incompleteControlTrioVarUmin_Dominant_ChrX_FemaleNB = [0.25,0,0.25,0,0.25,0.25,0.25,0,0.25,0,0.25,0,0,0,0.25,0]
incompleteCaseTrioUmax_Dominant_ChrX_FemaleNB = [0.5,0,0.5,0,0.5,-0.5,0.5,0,0.5,0,0.5,0,0,0,0.5,0]
incompleteCaseTrioVarUmax_Dominant_ChrX_FemaleNB = [0.25,0,0.25,0,0.25,0.25,0.25,0,0.25,0,0.25,0,0,0,0.25,0]
incompleteControlTrioUmax_Dominant_ChrX_FemaleNB = [-0.5,-0.5,0,0,-0.5,-0.5,0,0,-0.5,-0.5,0,0,0,0,-0.5,0]
incompleteControlTrioVarUmax_Dominant_ChrX_FemaleNB = [0.25,0.25,0,0,0.25,0.25,0,0,0.25,0.25,0,0,0,0,0.25,0]

incompleteCaseTrioUmin_Recessive_ChrX_FemaleNB = [-0.5,0,-0.5,0,-0.5,0,-0.5,0.5,0,0,0,-0.5,-0.5,0,0,-0.5]
incompleteCaseTrioVarUmin_Recessive_ChrX_FemaleNB = [0.25,0,0.25,0,0.25,0,0.25,0.25,0,0,0,0.25,0.25,0,0,0.25]
incompleteControlTrioUmin_Recessive_ChrX_FemaleNB = [0.5,0,0,0.5,0.5,0,0,0.5,0,0,0,0.5,0,0.5,0,0.5]
incompleteControlTrioVarUmin_Recessive_ChrX_FemaleNB = [0.25,0,0,0.25,0.25,0,0,0.25,0,0,0,0.25,0,0.25,0,0.25]
incompleteCaseTrioUmax_Recessive_ChrX_FemaleNB = [0.5,0,0,0.5,0.5,0,0,0.5,0,0,0,0.5,0,0.5,0,0.5]
incompleteCaseTrioVarUmax_Recessive_ChrX_FemaleNB = [0.25,0,0,0.25,0.25,0,0,0.25,0,0,0,0.25,0,0.25,0,0.25]
incompleteControlTrioUmax_Recessive_ChrX_FemaleNB = [-0.5,0,-0.5,0,-0.5,0,-0.5,0.5,0,0,0,-0.5,-0.5,0,0,-0.5]
incompleteControlTrioVarUmax_Recessive_ChrX_FemaleNB = [0.25,0,0.25,0,0.25,0,0.25,0.25,0,0,0,0.25,0.25,0,0,0.25]


##----------------------------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------------------------
class AutosomalMarker:
	
	###METHODS
	def __init__(self):
		#MEMBER VARIABLES
	        self.markerID = ""
	        self.sampleGenotypes = {}
	        self.maf = 0.0
	        self.nVariantType = [0]*4   #[refHomozygous, heterozygous, nonRefHomozygous, missing]
	        self.nCompleteInformativeCaseTrio_MaleNB = [0]*7
		self.nCompleteInformativeCaseTrio_FemaleNB = [0]*7
		self.nCompleteInformativeControlTrio_MaleNB = [0]*7
	        self.nCompleteInformativeControlTrio_FemaleNB = [0]*7
	        self.nCompleteNonInformativeCaseTrio = [0]*3
		self.nCompleteNonInformativeControlTrio = [0]*3
	        self.nIncompleteInformativeCaseTrio_MaleNB = [0]*17
		self.nIncompleteInformativeCaseTrio_FemaleNB = [0]*17
		self.nIncompleteInformativeControlTrio_MaleNB = [0]*17
		self.nIncompleteInformativeControlTrio_FemaleNB = [0]*17
		self.nIncompleteNonInformativeCaseTrio = [0]*3
		self.nIncompleteNonInformativeControlTrio = [0]*3
		self.nMIE = 0
		#I considered a class for TDT itself, but then TDT is always associated with a marker; the class will have no use on its own. So in favor of keeping things simple...
	        #These values will be different for different genetic models(a,d,r), but I don't need to save them all (i'm sending them to print directly), so overwriting these variable is fine.
	        self.bComplete = 0   #no. of ref allele transmissions within complete cases
	        self.cComplete = 0   #no. of alt allele transmissions within complete cases
	        self.bMaxIncrement = 0
	        self.cMaxIncrement = 0
	        self.bMin = 0
	        self.cMin = 0
	        self.bMax = 0
	        self.cMax = 0
	        #stdTDT
	        self.chiSq_StdTDT = 0.0
	        self.pValue_StdTDT = 0.0
	        #extendedTDT
	        self.minChiSq_rTDT = 0.0
	        self.minPValue_rTDT = 0.0
	        self.maxChiSq_rTDT = 0.0
	        self.maxPValue_rTDT = 0.0
        	#miTDT
	        self.minChiSq_miTDT = 0.0
	        self.minPValue_miTDT = 0.0
	        self.maxChiSq_miTDT = 0.0
	        self.maxPValue_miTDT = 0.0
		#stdFBAT
		self.Z_stdFBAT = 0.0
		self.pValue_stdFBAT = 0.0
		self.caseU = 0.0
		self.caseVarU = 0.0
		self.controlU = 0.0
		self.controlVarU = 0.0
		#rFBAT
		self.minZ_extFBAT = 0.0
		self.minPValue_extFBAT = 0.0
		self.maxZ_extFBAT = 0.0
                self.maxPValue_extFBAT = 0.0		
	
	###Std TDT
	def stdTDT(self,MODEL):
        	#stdTDT statistic (chiSq) is always positive. It doesn't matter whether you test ref vs. alt or major vs. minor
	        if MODEL.lower() == "a":
        	        #Note: nInformativeTrio[3] counts in both b and c. adding or removing it does not affect the numerator b-c, but does affect the denominator b+c
               		#so we keep it, because this trio transmits 1 REF allele and 1 ALT allele
        	        self.cComplete = sum([a*b for a,b in zip(informativeTrioCinc_Additive,self.nCompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioCinc_Additive,self.nCompleteInformativeCaseTrio_FemaleNB)]) #no. of alt alleles transmitted
	                self.bComplete = sum([a*b for a,b in zip(informativeTrioBinc_Additive,self.nCompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioBinc_Additive,self.nCompleteInformativeCaseTrio_FemaleNB)]) #no. of ref alleles transmitted
	        elif MODEL.lower() == "d":
        	        self.cComplete = sum([a*b for a,b in zip(informativeTrioCinc_Dominant,self.nCompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioCinc_Dominant,self.nCompleteInformativeCaseTrio_FemaleNB)]) #no. of alt alleles transmitted
	                self.bComplete = sum([a*b for a,b in zip(informativeTrioBinc_Dominant,self.nCompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioBinc_Dominant,self.nCompleteInformativeCaseTrio_FemaleNB)]) #no. of ref alleles transmitted
	        elif MODEL.lower() == "r":
        	        self.cComplete = sum([a*b for a,b in zip(informativeTrioCinc_Recessive,self.nCompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioCinc_Recessive,self.nCompleteInformativeCaseTrio_FemaleNB)]) #no. of alt alleles transmitted
	                self.bComplete = sum([a*b for a,b in zip(informativeTrioBinc_Recessive,self.nCompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioBinc_Recessive,self.nCompleteInformativeCaseTrio_FemaleNB)]) #no. of ref alleles transmitted

        	try:
	                self.chiSq_StdTDT = (self.bComplete-self.cComplete)**2/float(self.bComplete+self.cComplete)
        	        self.pValue_StdTDT = 1 - stats.chi2.cdf(self.chiSq_StdTDT, 1)   #multiply by 2 for two-sided test??? No, since TDT statistic is always positive, we want to compute the probability of observing a score higher than this
	        except ZeroDivisionError:
        	        self.chiSq_StdTDT = 'NA'
	                self.pValue_StdTDT = 'NA'

	
	### rTDT
	def extendedTDT(self,MODEL):

                if MODEL.lower() == "a":
                        #compute maximum increments to b and c
                        self.bMaxIncrement = sum([a*b for a,b in zip(incompleteTrioBinc_Additive,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteTrioBinc_Additive,self.nIncompleteInformativeCaseTrio_FemaleNB)])
                        self.cMaxIncrement = sum([a*b for a,b in zip(incompleteTrioCinc_Additive,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteTrioCinc_Additive,self.nIncompleteInformativeCaseTrio_FemaleNB)])

                        self.bMin = self.bComplete + self.nIncompleteInformativeCaseTrio_MaleNB[8] + self.nIncompleteInformativeCaseTrio_FemaleNB[8]
                        self.cMin = self.cComplete + self.nIncompleteInformativeCaseTrio_MaleNB[10] + self.nIncompleteInformativeCaseTrio_FemaleNB[10]
                elif MODEL.lower() == "d":
                        #compute maximum increments to b and c
                        self.bMaxIncrement = sum([a*b for a,b in zip(incompleteTrioBinc_Dominant,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteTrioBinc_Dominant,self.nIncompleteInformativeCaseTrio_FemaleNB)])
                        self.cMaxIncrement = sum([a*b for a,b in zip(incompleteTrioCinc_Dominant,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteTrioCinc_Dominant,self.nIncompleteInformativeCaseTrio_FemaleNB)])

                        self.bMin = self.bComplete + self.nIncompleteInformativeCaseTrio_MaleNB[8] + self.nIncompleteInformativeCaseTrio_FemaleNB[8]
                        self.cMin = self.cComplete #none of the incomplete trio types exclusively and always increase c without increasing b
                elif MODEL.lower() == "r":
                        #compute maximum increments to b and c
                        self.bMaxIncrement = sum([a*b for a,b in zip(incompleteTrioBinc_Recessive,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteTrioBinc_Recessive,self.nIncompleteInformativeCaseTrio_FemaleNB)])
                        self.cMaxIncrement = sum([a*b for a,b in zip(incompleteTrioCinc_Recessive,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteTrioCinc_Recessive,self.nIncompleteInformativeCaseTrio_FemaleNB)])

                        self.bMin = self.bComplete #none of the incomplete trio types exclusively and always increase b without increasing c
                        self.cMin = self.cComplete + self.nIncompleteInformativeCaseTrio_MaleNB[10] + self.nIncompleteInformativeCaseTrio_FemaleNB[10]
		
		self.bMax = self.bComplete + self.bMaxIncrement
	        self.cMax = self.cComplete + self.cMaxIncrement	
		[self.minChiSq_rTDT,self.minPValue_rTDT,self.maxChiSq_rTDT,self.maxPValue_rTDT] = self.getMinMaxStatistic() 


	##Assign min and max test statistic and p-value based on bMin, cMin, bMax and cMax
	def getMinMaxStatistic(self):
		try:
                        extTDT1 = (self.bMin-self.cMax)**2/float(self.bMin+self.cMax)
                        pValue1 = 1 - stats.chi2.cdf(extTDT1, 1)
                except ZeroDivisionError:
                        extTDT1 = 'NA'
                        pValue1 = 'NA'
                try:
                        extTDT2 = (self.bMax-self.cMin)**2/float(self.bMax+self.cMin)
                        pValue2 = 1 - stats.chi2.cdf(extTDT2, 1)
                except ZeroDivisionError:
                        extTDT2 = 'NA'
                        pValue2 = 'NA'
                if self.bMin >= self.cMax:
			return [extTDT1,pValue1,extTDT2,pValue2]
                elif self.bMax <= self.cMin:
			return [extTDT2,pValue2,extTDT1,pValue1]
                else:  #overlapping ranges of b and c
                        return [0,1,max(extTDT1,extTDT2),min(pValue1,pValue2)]

	###Get sample genotypes, only for pedigrees that the user wants to analyze
	def getSampleGenotypes(self,vcfValues,pedMemberIndices):
		for ped in pedMemberIndices.keys():
			memberIndices = pedMemberIndices[ped]
			self.sampleGenotypes[ped]=[vcfValues[x] for x in memberIndices]

	###Check validity of genotypes
        def hasValidGenotypes(self):
                for thisPed in self.sampleGenotypes.keys():
                        thisPedGenotypes = self.sampleGenotypes[thisPed]
			try:
				assert(all(member in ['0','1','2','NA'] for member in thisPedGenotypes))
			except(AssertionError):
				return 0
		return 1

	###Compute Minor Allele Frequency
	def computeMAF(self):
		refCount = 0
		altCount = 0
	        #COMPUTE ALLELE FREQUENCY from feature matrix (Note: missing values might bias this computation)
		for thisPed in self.sampleGenotypes.keys():  #assuming all diploid calls for autosomal markers
			thisPedGenotypes = self.sampleGenotypes[thisPed]
			refCount += sum((2-int(member)) for member in thisPedGenotypes if member.isdigit())
	       		altCount += sum(int(member) for member in thisPedGenotypes if member.isdigit()) 

        	if refCount < altCount:
                	self.maf = refCount/float(refCount+altCount)
	        else:
        		self.maf = altCount/float(refCount+altCount)

        ### Get variant distribution
	#NOTE: same method for chrX. So, '0' for males in chrX is counted as ref.homo., '1' is counted as het.
        def getVariantDistribution(self):
                #get variant distribution
		for ped in self.sampleGenotypes.keys():
			pedGenotypes = self.sampleGenotypes[ped]
			for genotype in pedGenotypes:
				if genotype.isdigit():
					self.nVariantType[int(genotype)] += 1
				else:  #NA genotype
					self.nVariantType[3] +=1 

	###Populate Trio Type vectors
	def populateTrioTypeCountVectors(self,pedMemberType,pedPhenoDict,pedNBGender):
	        #loop through each pedigree from the list pedIDs
	        for thisPed in self.sampleGenotypes.keys():
			thisPedGenotypes = self.sampleGenotypes[thisPed]
	                thisMemberType = pedMemberType[thisPed]
			thisNBGender = pedNBGender[thisPed]
			thisNBPheno = pedPhenoDict[thisPed]
	                # get numeric genotype code
	                genoF = thisPedGenotypes[thisMemberType.index('1')]  #Father
	                genoM = thisPedGenotypes[thisMemberType.index('2')]  #Mother
	                genoNB = thisPedGenotypes[thisMemberType.index('3')]  #Newborn

#        	        print 'reading fm',self.markerID
#	                print pedIDs[loopIndex],thisPedIndices,thisMemberType
#	                print genoF,genoM,genoNB

        	        #populate trio type count vectors
	                #if complete trio info available
	                if genoF <> "NA" and genoM <> "NA" and genoNB <> "NA":
        	                try:
                	                #increment counts of case and control complete trios
	                                if thisNBPheno == 2:
						if thisNBGender == 2:
							self.nCompleteInformativeCaseTrio_FemaleNB[informativeTrioType.index([set([genoF,genoM]),genoNB])] += 1
						elif thisNBGender == 1:
							self.nCompleteInformativeCaseTrio_MaleNB[informativeTrioType.index([set([genoF,genoM]),genoNB])] += 1
	                                elif thisNBPheno == 1:
						if thisNBGender == 2:
							self.nCompleteInformativeControlTrio_FemaleNB[informativeTrioType.index([set([genoF,genoM]),genoNB])] += 1
						elif thisNBGender == 1:
							self.nCompleteInformativeControlTrio_MaleNB[informativeTrioType.index([set([genoF,genoM]),genoNB])] += 1

                	        except (KeyError,ValueError):    # trio not informative
	                                try:
						if thisNBPheno == 2:
							self.nCompleteNonInformativeCaseTrio[nonInformativeMatingType.index(set([genoF,genoM]))] += 1
						elif thisNBPheno == 1:
                                                        self.nCompleteNonInformativeControlTrio[nonInformativeMatingType.index(set([genoF,genoM]))] += 1
	                                except (KeyError,ValueError):
        	                                print 'Unmatched trio type at ',self.markerID,': ',thisPed,' [',genoF,genoM,genoNB,']. Counting as MIE' 
                	                        self.nMIE += 1
	                else:     #atleast one member genotype missing; incomplete trio
        	                try:
                	                if thisNBPheno == 2:
						if thisNBGender == 2:
	                        	                self.nIncompleteInformativeCaseTrio_FemaleNB[incompleteTrioType.index([set([genoF,genoM]),genoNB])] += 1
						elif thisNBGender == 1:
                                                	self.nIncompleteInformativeCaseTrio_MaleNB[incompleteTrioType.index([set([genoF,genoM]),genoNB])] += 1
	                                elif thisNBPheno == 1:
						if thisNBGender == 2:
        	                                	self.nIncompleteInformativeControlTrio_FemaleNB[incompleteTrioType.index([set([genoF,genoM]),genoNB])] += 1
						elif thisNBGender == 1:
                                                	self.nIncompleteInformativeControlTrio_MaleNB[incompleteTrioType.index([set([genoF,genoM]),genoNB])] += 1
                	        except (KeyError,ValueError):
					try:
						if thisNBPheno == 2:
							self.nIncompleteNonInformativeCaseTrio[nonInformativeMatingType.index(set([genoF,genoM]))] += 1
						elif thisNBPheno == 1:
							self.nIncompleteNonInformativeControlTrio[nonInformativeMatingType.index(set([genoF,genoM]))] += 1
					except (KeyError,ValueError):
	                        	        print 'Unmatched trio type at ',self.markerID,': ',thisPed,' [',genoF,genoM,genoNB,']. Counting as MIE'
		                                self.nMIE += 1

	def stdFBAT(self,MODEL,OFFSET):
		#Default offset is 0.5 
		#FBAT statistic has equal magnitude but opposite signs for the two alleles. 
		if MODEL.lower() == "a":
	                #compute U
	                self.caseU = (1-OFFSET) * (sum([a*b for a,b in zip(informativeTrioU_Additive,self.nCompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioU_Additive,self.nCompleteInformativeCaseTrio_FemaleNB)]))
	                self.controlU = (0-OFFSET) * (sum([a*b for a,b in zip(informativeTrioU_Additive,self.nCompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioU_Additive,self.nCompleteInformativeControlTrio_FemaleNB)]))
	                #compute Variance(U); FBAT tech report says Var(U) = sum(Var(Si)) assuming that all Si's are independent
	                self.caseVarU = ((1-OFFSET)**2) * (sum([a*b for a,b in zip(informativeTrioVarU_Additive,self.nCompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioVarU_Additive,self.nCompleteInformativeCaseTrio_FemaleNB)]))
		        self.controlVarU = ((0-OFFSET)**2) * (sum([a*b for a,b in zip(informativeTrioVarU_Additive,self.nCompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioVarU_Additive,self.nCompleteInformativeControlTrio_FemaleNB)]))
		elif MODEL.lower() == "d":
			#compute U
                        self.caseU = (1-OFFSET) * (sum([a*b for a,b in zip(informativeTrioU_Dominant,self.nCompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioU_Dominant,self.nCompleteInformativeCaseTrio_FemaleNB)]))
                        self.controlU = (0-OFFSET) * (sum([a*b for a,b in zip(informativeTrioU_Dominant,self.nCompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioU_Dominant,self.nCompleteInformativeControlTrio_FemaleNB)]))
                        #compute Variance(U); FBAT tech report says Var(U) = sum(Var(Si)) assuming that all Si's are independent
                        self.caseVarU = ((1-OFFSET)**2) * (sum([a*b for a,b in zip(informativeTrioVarU_Dominant,self.nCompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioVarU_Dominant,self.nCompleteInformativeCaseTrio_FemaleNB)]))
                        self.controlVarU = ((0-OFFSET)**2) * (sum([a*b for a,b in zip(informativeTrioVarU_Dominant,self.nCompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioVarU_Dominant,self.nCompleteInformativeControlTrio_FemaleNB)]))
		elif MODEL.lower() == "r":
			#compute U
                        self.caseU = (1-OFFSET) * (sum([a*b for a,b in zip(informativeTrioU_Recessive,self.nCompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioU_Recessive,self.nCompleteInformativeCaseTrio_FemaleNB)]))
                        self.controlU = (0-OFFSET) * (sum([a*b for a,b in zip(informativeTrioU_Recessive,self.nCompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioU_Recessive,self.nCompleteInformativeControlTrio_FemaleNB)]))
                        #compute Variance(U); FBAT tech report says Var(U) = sum(Var(Si)) assuming that all Si's are independent
                        self.caseVarU = ((1-OFFSET)**2) * (sum([a*b for a,b in zip(informativeTrioVarU_Recessive,self.nCompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioVarU_Recessive,self.nCompleteInformativeCaseTrio_FemaleNB)]))
                        self.controlVarU = ((0-OFFSET)**2) * (sum([a*b for a,b in zip(informativeTrioVarU_Recessive,self.nCompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioVarU_Recessive,self.nCompleteInformativeControlTrio_FemaleNB)]))

                U = self.caseU + self.controlU
                varU = self.caseVarU + self.controlVarU
                try:
                	self.Z_stdFBAT = U / float(math.sqrt(varU))
                        self.pValue_stdFBAT = (1 - stats.norm.cdf(abs(self.Z_stdFBAT)))*2.0
                except ZeroDivisionError:
                	self.Z_stdFBAT = 'NA'
                        self.pValue_stdFBAT = 'NA'

	#NOTE: this is just rFBAT. allele frequencies have not been incorporated yet.
	def extendedFBAT(self,MODEL,OFFSET):
		if MODEL.lower() == "a":
			#min statistic
	                caseUmin = self.caseU + (1-OFFSET) * (sum([a*b for a,b in zip(incompleteCaseTrioUmin_Additive,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteCaseTrioUmin_Additive,self.nIncompleteInformativeCaseTrio_FemaleNB)]))
                        caseVarUmin = self.caseVarU + ((1-OFFSET)**2) * (sum([a*b for a,b in zip(incompleteCaseTrioVarUmin_Additive,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteCaseTrioVarUmin_Additive,self.nIncompleteInformativeCaseTrio_FemaleNB)]))
                        controlUmin = self.controlU + (0-OFFSET) * (sum([a*b for a,b in zip(incompleteControlTrioUmin_Additive,self.nIncompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteControlTrioUmin_Additive,self.nIncompleteInformativeControlTrio_FemaleNB)]))
                        controlVarUmin = self.controlVarU + ((0-OFFSET)**2) * (sum([a*b for a,b in zip(incompleteControlTrioVarUmin_Additive,self.nIncompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteControlTrioVarUmin_Additive,self.nIncompleteInformativeControlTrio_FemaleNB)]))

                        #max statistic
                        caseUmax =  self.caseU + (1-OFFSET) * (sum([a*b for a,b in zip(incompleteCaseTrioUmax_Additive,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteCaseTrioUmax_Additive,self.nIncompleteInformativeCaseTrio_FemaleNB)]))
                        caseVarUmax = self.caseVarU + ((1-OFFSET)**2) * (sum([a*b for a,b in zip(incompleteCaseTrioVarUmax_Additive,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteCaseTrioVarUmax_Additive,self.nIncompleteInformativeCaseTrio_FemaleNB)]))
                        controlUmax =  self.controlU + (0-OFFSET) * (sum([a*b for a,b in zip(incompleteControlTrioUmax_Additive,self.nIncompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteControlTrioUmax_Additive,self.nIncompleteInformativeControlTrio_FemaleNB)]))
                        controlVarUmax = self.controlVarU + ((0-OFFSET)**2) * (sum([a*b for a,b in zip(incompleteControlTrioVarUmax_Additive,self.nIncompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteControlTrioVarUmax_Additive,self.nIncompleteInformativeControlTrio_FemaleNB)]))

		elif MODEL.lower() == "d":
			#min statistic
                        caseUmin = self.caseU + (1-OFFSET) * (sum([a*b for a,b in zip(incompleteCaseTrioUmin_Dominant,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteCaseTrioUmin_Dominant,self.nIncompleteInformativeCaseTrio_FemaleNB)]))
                        caseVarUmin = self.caseVarU + ((1-OFFSET)**2) * (sum([a*b for a,b in zip(incompleteCaseTrioVarUmin_Dominant,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteCaseTrioVarUmin_Dominant,self.nIncompleteInformativeCaseTrio_FemaleNB)]))
                        controlUmin = self.controlU + (0-OFFSET) * (sum([a*b for a,b in zip(incompleteControlTrioUmin_Dominant,self.nIncompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteControlTrioUmin_Dominant,self.nIncompleteInformativeControlTrio_FemaleNB)]))
                        controlVarUmin = self.controlVarU + ((0-OFFSET)**2) * (sum([a*b for a,b in zip(incompleteControlTrioVarUmin_Dominant,self.nIncompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteControlTrioVarUmin_Dominant,self.nIncompleteInformativeControlTrio_FemaleNB)]))

                        #max statistic
                        caseUmax =  self.caseU + (1-OFFSET) * (sum([a*b for a,b in zip(incompleteCaseTrioUmax_Dominant,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteCaseTrioUmax_Dominant,self.nIncompleteInformativeCaseTrio_FemaleNB)]))
                        caseVarUmax = self.caseVarU + ((1-OFFSET)**2) * (sum([a*b for a,b in zip(incompleteCaseTrioVarUmax_Dominant,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteCaseTrioVarUmax_Dominant,self.nIncompleteInformativeCaseTrio_FemaleNB)]))
                        controlUmax =  self.controlU + (0-OFFSET) * (sum([a*b for a,b in zip(incompleteControlTrioUmax_Dominant,self.nIncompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteControlTrioUmax_Dominant,self.nIncompleteInformativeControlTrio_FemaleNB)]))
                        controlVarUmax = self.controlVarU + ((0-OFFSET)**2) * (sum([a*b for a,b in zip(incompleteControlTrioVarUmax_Dominant,self.nIncompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteControlTrioVarUmax_Dominant,self.nIncompleteInformativeControlTrio_FemaleNB)]))

		elif MODEL.lower() == "r":
			#min statistic
                        caseUmin = self.caseU + (1-OFFSET) * (sum([a*b for a,b in zip(incompleteCaseTrioUmin_Recessive,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteCaseTrioUmin_Recessive,self.nIncompleteInformativeCaseTrio_FemaleNB)]))
                        caseVarUmin = self.caseVarU + ((1-OFFSET)**2) * (sum([a*b for a,b in zip(incompleteCaseTrioVarUmin_Recessive,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteCaseTrioVarUmin_Recessive,self.nIncompleteInformativeCaseTrio_FemaleNB)]))
                        controlUmin = self.controlU + (0-OFFSET) * (sum([a*b for a,b in zip(incompleteControlTrioUmin_Recessive,self.nIncompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteControlTrioUmin_Recessive,self.nIncompleteInformativeControlTrio_FemaleNB)]))
                        controlVarUmin = self.controlVarU + ((0-OFFSET)**2) * (sum([a*b for a,b in zip(incompleteControlTrioVarUmin_Recessive,self.nIncompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteControlTrioVarUmin_Recessive,self.nIncompleteInformativeControlTrio_FemaleNB)]))

                        #max statistic
                        caseUmax =  self.caseU + (1-OFFSET) * (sum([a*b for a,b in zip(incompleteCaseTrioUmax_Recessive,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteCaseTrioUmax_Recessive,self.nIncompleteInformativeCaseTrio_FemaleNB)]))
                        caseVarUmax = self.caseVarU + ((1-OFFSET)**2) * (sum([a*b for a,b in zip(incompleteCaseTrioVarUmax_Recessive,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteCaseTrioVarUmax_Recessive,self.nIncompleteInformativeCaseTrio_FemaleNB)]))
                        controlUmax =  self.controlU + (0-OFFSET) * (sum([a*b for a,b in zip(incompleteControlTrioUmax_Recessive,self.nIncompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteControlTrioUmax_Recessive,self.nIncompleteInformativeControlTrio_FemaleNB)]))
                        controlVarUmax = self.controlVarU + ((0-OFFSET)**2) * (sum([a*b for a,b in zip(incompleteControlTrioVarUmax_Recessive,self.nIncompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteControlTrioVarUmax_Recessive,self.nIncompleteInformativeControlTrio_FemaleNB)]))
			

		Umin = caseUmin + controlUmin
                varUmin = caseVarUmin + controlVarUmin
                try:
                	self.minZ_extFBAT = Umin/float(math.sqrt(varUmin))
                        self.minPValue_extFBAT = (1 - stats.norm.cdf(abs(self.minZ_extFBAT)))*2.0
                except:
                	self.minZ_extFBAT = 'NA'
                        self.minPValue_extFBAT = 'NA'

		Umax = caseUmax + controlUmax
                varUmax = caseVarUmax + controlVarUmax
                try:
                	self.maxZ_extFBAT = Umax/float(math.sqrt(varUmax))
                        self.maxPValue_extFBAT = (1 - stats.norm.cdf(abs(self.maxZ_extFBAT)))*2.0
                except:
                	self.maxZ_extFBAT = 'NA'
                        self.maxPValue_extFBAT = 'NA'
	

##----------------------------------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------------------------------
class ChrXMarker(AutosomalMarker):
	
	 ###METHODS
        def __init__(self):
		AutosomalMarker.__init__(self)
                #MEMBER VARIABLES 
                self.nCompleteInformativeCaseTrio_MaleNB = [0]*4
		self.nCompleteInformativeCaseTrio_FemaleNB = [0]*4
		self.nCompleteInformativeControlTrio_MaleNB = [0]*4
		self.nCompleteInformativeControlTrio_FemaleNB = [0]*4
                self.nCompleteNonInformativeCaseTrio = [0]*2
		self.nCompleteNonInformativeControlTrio = [0]*2
                self.nIncompleteInformativeCaseTrio_MaleNB = [0]*14
                self.nIncompleteInformativeControlTrio_MaleNB = [0]*14
		self.nIncompleteInformativeCaseTrio_FemaleNB = [0]*16
                self.nIncompleteInformativeControlTrio_FemaleNB = [0]*16
		self.nIncompleteNonInformativeCaseTrio = [0]*2
		self.nIncompleteNonInformativeControlTrio = [0]*2
		
	###Check validity of genotypes
	def hasValidGenotypes(self,pedNBGender,pedMemberType):
		for thisPed in self.sampleGenotypes.keys():
                        thisPedGenotypes = self.sampleGenotypes[thisPed]
			#check parents first
			if thisPedGenotypes[pedMemberType[thisPed].index('1')] not in ['0','1','NA'] or thisPedGenotypes[pedMemberType[thisPed].index('2')] not in ['0','1','2','NA']:
				return 0
			#check NB conditional on gender
			elif (pedNBGender[thisPed]==1 and thisPedGenotypes[pedMemberType[thisPed].index('3')] not in ['0','1','NA']) or (pedNBGender[thisPed]==2 and thisPedGenotypes[pedMemberType[thisPed].index('3')] not in ['0','1','2','NA']):
				return 0
		return 1

	###Compute Minor Allele Frequency
        def computeMAF(self,pedNBGender,pedMemberType):
		altCount = 0
		refCount = 0
                #COMPUTE ALLELE FREQUENCY from feature matrix (Note: missing values might bias this computation)
		for thisPed in self.sampleGenotypes.keys():
			thisPedGenotypes = self.sampleGenotypes[thisPed]
			genoF = thisPedGenotypes[pedMemberType[thisPed].index('1')]
			genoM = thisPedGenotypes[pedMemberType[thisPed].index('2')]
			genoNB = thisPedGenotypes[pedMemberType[thisPed].index('3')]

			if genoF.isdigit():
				refCount += 1-int(genoF)
				altCount += int(genoF)
			if genoM.isdigit():
				refCount += 2-int(genoM)
				altCount += int(genoM)
			if genoNB.isdigit():
				if pedNBGender[thisPed] == 1:
					refCount += 1-int(genoNB)
	                                altCount += int(genoNB)	 
				elif pedNBGender[thisPed] == 2:
					refCount += 2-int(genoNB)
                                        altCount += int(genoNB)
		if refCount < altCount:
                        self.maf = refCount/float(refCount+altCount)
                else:
                        self.maf = altCount/float(refCount+altCount)
	

	def getVariantDistribution(self,pedNBGender,pedMemberType):
        	#get variant distribution
                for thisPed in self.sampleGenotypes.keys():
                	thisPedGenotypes = self.sampleGenotypes[thisPed]
                        for index in range(len(thisPedGenotypes)):
				genotype = thisPedGenotypes[index]
				if genotype.isdigit():
					if pedMemberType[thisPed][index]=='1' or (pedMemberType[thisPed][index]=='3' and pedNBGender[thisPed]==1):
						if genotype == '0':
							self.nVariantType[0] += 1
						elif genotype == '1':
							self.nVariantType[2] += 1
					elif pedMemberType[thisPed][index]=='2' or (pedMemberType[thisPed][index]=='3' and pedNBGender[thisPed]==2):
						self.nVariantType[int(genotype)] += 1	
                                else:  #NA genotype
                                        self.nVariantType[3] +=1
			
	
	###Populate Trio Type vectors
        def populateTrioTypeCountVectors(self,pedMemberType,pedPhenoDict,pedNBGender):
                #loop through each pedigree from the list pedIDs
                for thisPed in self.sampleGenotypes.keys():
			thisPedGenotypes = self.sampleGenotypes[thisPed]
                        thisMemberType = pedMemberType[thisPed]
			thisNBGender = pedNBGender[thisPed]
			thisNBPheno = pedPhenoDict[thisPed]
                        # get numeric genotype code
                        genoF = thisPedGenotypes[thisMemberType.index('1')]  #Father
                        genoM = thisPedGenotypes[thisMemberType.index('2')]  #Mother
                        genoNB = thisPedGenotypes[thisMemberType.index('3')]  #Newborn

#                       print 'reading fm',self.markerID
#                       print pedIDs[loopIndex],thisPedIndices,thisMemberType
#                       print genoF,genoM,genoNB
                        #populate trio type count vectors
                        #if complete trio info available
                        if genoF <> "NA" and genoM <> "NA" and genoNB <> "NA":
                                try:
                                        #increment counts of case and control complete trios
                                        if thisNBPheno == 2:
						if thisNBGender == 2:  #female
	                                                self.nCompleteInformativeCaseTrio_FemaleNB[informativeTrioType_ChrX_FemaleNB.index([genoF,genoM,genoNB])] += 1
						elif thisNBGender == 1: #male
							self.nCompleteInformativeCaseTrio_MaleNB[informativeTrioType_ChrX_MaleNB.index([genoF,genoM,genoNB])] += 1
                                        elif thisNBPheno == 1:
						if thisNBGender == 2:  #female
                                                        self.nCompleteInformativeControlTrio_FemaleNB[informativeTrioType_ChrX_FemaleNB.index([genoF,genoM,genoNB])] += 1
                                                elif thisNBGender == 1: #male
                                                        self.nCompleteInformativeControlTrio_MaleNB[informativeTrioType_ChrX_MaleNB.index([genoF,genoM,genoNB])] += 1

                                except (KeyError,ValueError):    # trio not informative
                                        try:
						if thisNBPheno == 2:
	                                                self.nCompleteNonInformativeCaseTrio[nonInformativeMotherType_ChrX.index(genoM)] += 1
						elif thisNBPheno == 1:
                                                        self.nCompleteNonInformativeControlTrio[nonInformativeMotherType_ChrX.index(genoM)] += 1	
                                        except (KeyError,ValueError):
                                                print 'Unmatched trio type at ',self.markerID,': ',thisPed,' [',genoF,genoM,genoNB,']. Counting as MIE'
                                                self.nMIE += 1
                        else:     #atleast one member genotype missing; incomplete trio
                                try:
                                        if pedPhenoDict[thisPed] == 2:
						if thisNBGender == 2:
							self.nIncompleteInformativeCaseTrio_FemaleNB[incompleteTrioType_ChrX_FemaleNB.index([genoF,genoM,genoNB])]+=1
						elif thisNBGender == 1:
							self.nIncompleteInformativeCaseTrio_MaleNB[incompleteTrioType_ChrX_MaleNB.index([genoF,genoM,genoNB])]+=1
                                        elif pedPhenoDict[thisPed] == 1:
						if thisNBGender == 2:
                                                        self.nIncompleteInformativeControlTrio_FemaleNB[incompleteTrioType_ChrX_FemaleNB.index([genoF,genoM,genoNB])]+=1
                                                elif thisNBGender == 1:
                                                        self.nIncompleteInformativeControlTrio_MaleNB[incompleteTrioType_ChrX_MaleNB.index([genoF,genoM,genoNB])]+=1
                                except (KeyError,ValueError):
                                        try:
						if thisNBPheno == 2:
	                                                self.nIncompleteNonInformativeCaseTrio[nonInformativeMotherType_ChrX.index(genoM)] += 1
						elif thisNBPheno == 1:
                                                        self.nIncompleteNonInformativeControlTrio[nonInformativeMotherType_ChrX.index(genoM)] += 1
                                        except (KeyError,ValueError):
                                                print 'Unmatched trio type at ',self.markerID,': ',thisPed,' [',genoF,genoM,genoNB,']. Counting as MIE'
                                                self.nMIE += 1	
		
	
	###Std TDT
        def stdTDT(self,MODEL):
                #stdTDT statistic (chiSq) is always positive. It doesn't matter whether you test ref vs. alt or major vs. minor
                if MODEL.lower() == "a":
                        self.cComplete = sum([a*b for a,b in zip(informativeTrioCinc_AllModels_ChrX_MaleNB,self.nCompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioCinc_Additive_ChrX_FemaleNB,self.nCompleteInformativeCaseTrio_FemaleNB)]) #no. of alt alleles transmitted
                        self.bComplete = sum([a*b for a,b in zip(informativeTrioBinc_AllModels_ChrX_MaleNB,self.nCompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioBinc_Additive_ChrX_FemaleNB,self.nCompleteInformativeCaseTrio_FemaleNB)]) #no. of ref alleles transmitted
                elif MODEL.lower() == "d":
                        self.cComplete = sum([a*b for a,b in zip(informativeTrioCinc_AllModels_ChrX_MaleNB,self.nCompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioCinc_Dominant_ChrX_FemaleNB,self.nCompleteInformativeCaseTrio_FemaleNB)]) #no. of alt alleles transmitted
                        self.bComplete = sum([a*b for a,b in zip(informativeTrioBinc_AllModels_ChrX_MaleNB,self.nCompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioBinc_Dominant_ChrX_FemaleNB,self.nCompleteInformativeCaseTrio_FemaleNB)]) #no. of ref alleles transmitted
                elif MODEL.lower() == "r":
			self.cComplete = sum([a*b for a,b in zip(informativeTrioCinc_AllModels_ChrX_MaleNB,self.nCompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioCinc_Recessive_ChrX_FemaleNB,self.nCompleteInformativeCaseTrio_FemaleNB)]) #no. of alt alleles transmitted
                        self.bComplete = sum([a*b for a,b in zip(informativeTrioBinc_AllModels_ChrX_MaleNB,self.nCompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioBinc_Recessive_ChrX_FemaleNB,self.nCompleteInformativeCaseTrio_FemaleNB)]) #no. of ref alleles transmitted

                try:
                        self.chiSq_StdTDT = (self.bComplete-self.cComplete)**2/float(self.bComplete+self.cComplete)
                        self.pValue_StdTDT = 1 - stats.chi2.cdf(self.chiSq_StdTDT, 1)   #TODO: multiply by 2 for two-sided test??? (TODO for all tdt tests below)
                except ZeroDivisionError:
                        self.chiSq_StdTDT = 'NA'
                        self.pValue_StdTDT = 'NA'


        def extendedTDT(self,MODEL):

                if MODEL.lower() == "a":
                        #compute maximum increments to b and c
                        self.bMaxIncrement = sum([a*b for a,b in zip(incompleteTrioBinc_AllModels_ChrX_MaleNB,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteTrioBinc_Additive_ChrX_FemaleNB,self.nIncompleteInformativeCaseTrio_FemaleNB)])
                        self.cMaxIncrement = sum([a*b for a,b in zip(incompleteTrioCinc_AllModels_ChrX_MaleNB,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteTrioCinc_Additive_ChrX_FemaleNB,self.nIncompleteInformativeCaseTrio_FemaleNB)])

                        self.bMin = self.bComplete + self.nIncompleteInformativeCaseTrio_MaleNB[4] + self.nIncompleteInformativeCaseTrio_FemaleNB[5]
                        self.cMin = self.cComplete + self.nIncompleteInformativeCaseTrio_MaleNB[5] + self.nIncompleteInformativeCaseTrio_FemaleNB[7]
                elif MODEL.lower() == "d":
                        #compute maximum increments to b and c
			self.bMaxIncrement = sum([a*b for a,b in zip(incompleteTrioBinc_AllModels_ChrX_MaleNB,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteTrioBinc_Dominant_ChrX_FemaleNB,self.nIncompleteInformativeCaseTrio_FemaleNB)])
                        self.cMaxIncrement = sum([a*b for a,b in zip(incompleteTrioCinc_AllModels_ChrX_MaleNB,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteTrioCinc_Dominant_ChrX_FemaleNB,self.nIncompleteInformativeCaseTrio_FemaleNB)])

                        self.bMin = self.bComplete + self.nIncompleteInformativeCaseTrio_MaleNB[4] + self.nIncompleteInformativeCaseTrio_FemaleNB[5]
                        self.cMin = self.cComplete + self.nIncompleteInformativeCaseTrio_MaleNB[5]#none of the incomplete trio types with female NB exclusively and always increase c without increasing b
                elif MODEL.lower() == "r":
                        #compute maximum increments to b and c
			self.bMaxIncrement = sum([a*b for a,b in zip(incompleteTrioBinc_AllModels_ChrX_MaleNB,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteTrioBinc_Recessive_ChrX_FemaleNB,self.nIncompleteInformativeCaseTrio_FemaleNB)])
                        self.cMaxIncrement = sum([a*b for a,b in zip(incompleteTrioCinc_AllModels_ChrX_MaleNB,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteTrioCinc_Recessive_ChrX_FemaleNB,self.nIncompleteInformativeCaseTrio_FemaleNB)])

                        self.bMin = self.bComplete + self.nIncompleteInformativeCaseTrio_MaleNB[4]#none of the incomplete trio types with female NB  exclusively and always increase b without increasing c 
                        self.cMin = self.cComplete + self.nIncompleteInformativeCaseTrio_MaleNB[5]+ self.nIncompleteInformativeCaseTrio_FemaleNB[7]

                self.bMax = self.bComplete + self.bMaxIncrement
                self.cMax = self.cComplete + self.cMaxIncrement
                [self.minChiSq_rTDT,self.minPValue_rTDT,self.maxChiSq_rTDT,self.maxPValue_rTDT] = self.getMinMaxStatistic()


        def stdFBAT(self,MODEL,OFFSET):
                #TODO: the comment below about directionality seems incorrect for Z score. you might have to take care of major vs minor allele
                #Since the statistic is always positive (does not have directionality), comparing ref and alt counts instead of minor and major counts does not affect the statistic.
                if MODEL.lower() == "a":
                        #compute U
                        self.caseU = (1-OFFSET) * (sum([a*b for a,b in zip(informativeTrioU_AllModels_ChrX_MaleNB,self.nCompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioU_Additive_ChrX_FemaleNB,self.nCompleteInformativeCaseTrio_FemaleNB)]))
                        self.controlU = (0-OFFSET) * (sum([a*b for a,b in zip(informativeTrioU_AllModels_ChrX_MaleNB,self.nCompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioU_Additive_ChrX_FemaleNB,self.nCompleteInformativeControlTrio_FemaleNB)]))
                        #compute Variance(U); FBAT tech report says Var(U) = sum(Var(Si)) assuming that all Si's are independent
                        self.caseVarU = ((1-OFFSET)**2) * (sum([a*b for a,b in zip(informativeTrioVarU_AllModels_ChrX_MaleNB,self.nCompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioVarU_Additive_ChrX_FemaleNB,self.nCompleteInformativeCaseTrio_FemaleNB)]))
                        self.controlVarU = ((0-OFFSET)**2) * (sum([a*b for a,b in zip(informativeTrioVarU_AllModels_ChrX_MaleNB,self.nCompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioVarU_Additive_ChrX_FemaleNB,self.nCompleteInformativeControlTrio_FemaleNB)]))
                elif MODEL.lower() == "d":
                        #compute U
			self.caseU = (1-OFFSET) * (sum([a*b for a,b in zip(informativeTrioU_AllModels_ChrX_MaleNB,self.nCompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioU_Dominant_ChrX_FemaleNB,self.nCompleteInformativeCaseTrio_FemaleNB)]))
                        self.controlU = (0-OFFSET) * (sum([a*b for a,b in zip(informativeTrioU_AllModels_ChrX_MaleNB,self.nCompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioU_Dominant_ChrX_FemaleNB,self.nCompleteInformativeControlTrio_FemaleNB)]))
                        #compute Variance(U); FBAT tech report says Var(U) = sum(Var(Si)) assuming that all Si's are independent
                        self.caseVarU = ((1-OFFSET)**2) * (sum([a*b for a,b in zip(informativeTrioVarU_AllModels_ChrX_MaleNB,self.nCompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioVarU_Dominant_ChrX_FemaleNB,self.nCompleteInformativeCaseTrio_FemaleNB)]))
                        self.controlVarU = ((0-OFFSET)**2) * (sum([a*b for a,b in zip(informativeTrioVarU_AllModels_ChrX_MaleNB,self.nCompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioVarU_Dominant_ChrX_FemaleNB,self.nCompleteInformativeControlTrio_FemaleNB)]))
                elif MODEL.lower() == "r":
			#compute U
                        self.caseU = (1-OFFSET) * (sum([a*b for a,b in zip(informativeTrioU_AllModels_ChrX_MaleNB,self.nCompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioU_Recessive_ChrX_FemaleNB,self.nCompleteInformativeCaseTrio_FemaleNB)]))
                        self.controlU = (0-OFFSET) * (sum([a*b for a,b in zip(informativeTrioU_AllModels_ChrX_MaleNB,self.nCompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioU_Recessive_ChrX_FemaleNB,self.nCompleteInformativeControlTrio_FemaleNB)]))
                        #compute Variance(U); FBAT tech report says Var(U) = sum(Var(Si)) assuming that all Si's are independent
                        self.caseVarU = ((1-OFFSET)**2) * (sum([a*b for a,b in zip(informativeTrioVarU_AllModels_ChrX_MaleNB,self.nCompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioVarU_Recessive_ChrX_FemaleNB,self.nCompleteInformativeCaseTrio_FemaleNB)]))
                        self.controlVarU = ((0-OFFSET)**2) * (sum([a*b for a,b in zip(informativeTrioVarU_AllModels_ChrX_MaleNB,self.nCompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(informativeTrioVarU_Recessive_ChrX_FemaleNB,self.nCompleteInformativeControlTrio_FemaleNB)]))

                U = self.caseU + self.controlU
                varU = self.caseVarU + self.controlVarU
                try:
                        self.Z_stdFBAT = U / float(math.sqrt(varU))
                        self.pValue_stdFBAT = (1 - stats.norm.cdf(abs(self.Z_stdFBAT)))*2.0
                except ZeroDivisionError:
                        self.Z_stdFBAT = 'NA'
                        self.pValue_stdFBAT = 'NA'


	
	def extendedFBAT(self,MODEL,OFFSET):
                if MODEL.lower() == "a":
                        #min statistic
                        caseUmin = self.caseU + (1-OFFSET) * (sum([a*b for a,b in zip(incompleteCaseTrioUmin_AllModels_ChrX_MaleNB,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteCaseTrioUmin_Additive_ChrX_FemaleNB,self.nIncompleteInformativeCaseTrio_FemaleNB)]))
                        caseVarUmin = self.caseVarU + ((1-OFFSET)**2) * (sum([a*b for a,b in zip(incompleteCaseTrioVarUmin_AllModels_ChrX_MaleNB,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteCaseTrioVarUmin_Additive_ChrX_FemaleNB,self.nIncompleteInformativeCaseTrio_FemaleNB)]))
                        controlUmin = self.controlU + (0-OFFSET) * (sum([a*b for a,b in zip(incompleteControlTrioUmin_AllModels_ChrX_MaleNB,self.nIncompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteControlTrioUmin_Additive_ChrX_FemaleNB,self.nIncompleteInformativeControlTrio_FemaleNB)]))
                        controlVarUmin = self.controlVarU + ((0-OFFSET)**2) * (sum([a*b for a,b in zip(incompleteControlTrioVarUmin_AllModels_ChrX_MaleNB,self.nIncompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteControlTrioVarUmin_Additive_ChrX_FemaleNB,self.nIncompleteInformativeControlTrio_FemaleNB)]))

                        #max statistic
			caseUmax = self.caseU + (1-OFFSET) * (sum([a*b for a,b in zip(incompleteCaseTrioUmax_AllModels_ChrX_MaleNB,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteCaseTrioUmax_Additive_ChrX_FemaleNB,self.nIncompleteInformativeCaseTrio_FemaleNB)]))
                        caseVarUmax = self.caseVarU + ((1-OFFSET)**2) * (sum([a*b for a,b in zip(incompleteCaseTrioVarUmax_AllModels_ChrX_MaleNB,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteCaseTrioVarUmax_Additive_ChrX_FemaleNB,self.nIncompleteInformativeCaseTrio_FemaleNB)]))
                        controlUmax = self.controlU + (0-OFFSET) * (sum([a*b for a,b in zip(incompleteControlTrioUmax_AllModels_ChrX_MaleNB,self.nIncompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteControlTrioUmax_Additive_ChrX_FemaleNB,self.nIncompleteInformativeControlTrio_FemaleNB)]))
                        controlVarUmax = self.controlVarU + ((0-OFFSET)**2) * (sum([a*b for a,b in zip(incompleteControlTrioVarUmax_AllModels_ChrX_MaleNB,self.nIncompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteControlTrioVarUmax_Additive_ChrX_FemaleNB,self.nIncompleteInformativeControlTrio_FemaleNB)]))

                elif MODEL.lower() == "d":
			#min statistic
                        caseUmin = self.caseU + (1-OFFSET) * (sum([a*b for a,b in zip(incompleteCaseTrioUmin_AllModels_ChrX_MaleNB,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteCaseTrioUmin_Dominant_ChrX_FemaleNB,self.nIncompleteInformativeCaseTrio_FemaleNB)]))
                        caseVarUmin = self.caseVarU + ((1-OFFSET)**2) * (sum([a*b for a,b in zip(incompleteCaseTrioVarUmin_AllModels_ChrX_MaleNB,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteCaseTrioVarUmin_Dominant_ChrX_FemaleNB,self.nIncompleteInformativeCaseTrio_FemaleNB)]))
                        controlUmin = self.controlU + (0-OFFSET) * (sum([a*b for a,b in zip(incompleteControlTrioUmin_AllModels_ChrX_MaleNB,self.nIncompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteControlTrioUmin_Dominant_ChrX_FemaleNB,self.nIncompleteInformativeControlTrio_FemaleNB)]))
                        controlVarUmin = self.controlVarU + ((0-OFFSET)**2) * (sum([a*b for a,b in zip(incompleteControlTrioVarUmin_AllModels_ChrX_MaleNB,self.nIncompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteControlTrioVarUmin_Dominant_ChrX_FemaleNB,self.nIncompleteInformativeControlTrio_FemaleNB)]))

                        #max statistic
                        caseUmax = self.caseU + (1-OFFSET) * (sum([a*b for a,b in zip(incompleteCaseTrioUmax_AllModels_ChrX_MaleNB,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteCaseTrioUmax_Dominant_ChrX_FemaleNB,self.nIncompleteInformativeCaseTrio_FemaleNB)]))
                        caseVarUmax = self.caseVarU + ((1-OFFSET)**2) * (sum([a*b for a,b in zip(incompleteCaseTrioVarUmax_AllModels_ChrX_MaleNB,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteCaseTrioVarUmax_Dominant_ChrX_FemaleNB,self.nIncompleteInformativeCaseTrio_FemaleNB)]))
                        controlUmax = self.controlU + (0-OFFSET) * (sum([a*b for a,b in zip(incompleteControlTrioUmax_AllModels_ChrX_MaleNB,self.nIncompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteControlTrioUmax_Dominant_ChrX_FemaleNB,self.nIncompleteInformativeControlTrio_FemaleNB)]))
                        controlVarUmax = self.controlVarU + ((0-OFFSET)**2) * (sum([a*b for a,b in zip(incompleteControlTrioVarUmax_AllModels_ChrX_MaleNB,self.nIncompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteControlTrioVarUmax_Dominant_ChrX_FemaleNB,self.nIncompleteInformativeControlTrio_FemaleNB)]))

		elif MODEL.lower() == "r":
			#min statistic
                        caseUmin = self.caseU + (1-OFFSET) * (sum([a*b for a,b in zip(incompleteCaseTrioUmin_AllModels_ChrX_MaleNB,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteCaseTrioUmin_Recessive_ChrX_FemaleNB,self.nIncompleteInformativeCaseTrio_FemaleNB)]))
                        caseVarUmin = self.caseVarU + ((1-OFFSET)**2) * (sum([a*b for a,b in zip(incompleteCaseTrioVarUmin_AllModels_ChrX_MaleNB,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteCaseTrioVarUmin_Recessive_ChrX_FemaleNB,self.nIncompleteInformativeCaseTrio_FemaleNB)]))
                        controlUmin = self.controlU + (0-OFFSET) * (sum([a*b for a,b in zip(incompleteControlTrioUmin_AllModels_ChrX_MaleNB,self.nIncompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteControlTrioUmin_Recessive_ChrX_FemaleNB,self.nIncompleteInformativeControlTrio_FemaleNB)]))
                        controlVarUmin = self.controlVarU + ((0-OFFSET)**2) * (sum([a*b for a,b in zip(incompleteControlTrioVarUmin_AllModels_ChrX_MaleNB,self.nIncompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteControlTrioVarUmin_Recessive_ChrX_FemaleNB,self.nIncompleteInformativeControlTrio_FemaleNB)]))

                        #max statistic
                        caseUmax = self.caseU + (1-OFFSET) * (sum([a*b for a,b in zip(incompleteCaseTrioUmax_AllModels_ChrX_MaleNB,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteCaseTrioUmax_Recessive_ChrX_FemaleNB,self.nIncompleteInformativeCaseTrio_FemaleNB)]))
                        caseVarUmax = self.caseVarU + ((1-OFFSET)**2) * (sum([a*b for a,b in zip(incompleteCaseTrioVarUmax_AllModels_ChrX_MaleNB,self.nIncompleteInformativeCaseTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteCaseTrioVarUmax_Recessive_ChrX_FemaleNB,self.nIncompleteInformativeCaseTrio_FemaleNB)]))
                        controlUmax = self.controlU + (0-OFFSET) * (sum([a*b for a,b in zip(incompleteControlTrioUmax_AllModels_ChrX_MaleNB,self.nIncompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteControlTrioUmax_Recessive_ChrX_FemaleNB,self.nIncompleteInformativeControlTrio_FemaleNB)]))
                        controlVarUmax = self.controlVarU + ((0-OFFSET)**2) * (sum([a*b for a,b in zip(incompleteControlTrioVarUmax_AllModels_ChrX_MaleNB,self.nIncompleteInformativeControlTrio_MaleNB)]) + sum([a*b for a,b in zip(incompleteControlTrioVarUmax_Recessive_ChrX_FemaleNB,self.nIncompleteInformativeControlTrio_FemaleNB)]))
			

                Umin = caseUmin + controlUmin
                varUmin = caseVarUmin + controlVarUmin
                try:
                        self.minZ_extFBAT = Umin/float(math.sqrt(varUmin))
                        self.minPValue_extFBAT = (1 - stats.norm.cdf(abs(self.minZ_extFBAT)))*2.0
                except:
                        self.minZ_extFBAT = 'NA'
                        self.minPValue_extFBAT = 'NA'

                Umax = caseUmax + controlUmax
                varUmax = caseVarUmax + controlVarUmax
                try:
                        self.maxZ_extFBAT = Umax/float(math.sqrt(varUmax))
                        self.maxPValue_extFBAT = (1 - stats.norm.cdf(abs(self.maxZ_extFBAT)))*2.0
                except:
                        self.maxZ_extFBAT = 'NA'
                        self.maxPValue_extFBAT = 'NA'
