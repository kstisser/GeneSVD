import numpy as np
from scipy.linalg import svd
from enum import Enum

class GeneType(Enum):
    EYECOLOR = 1
    HAIRCOLOR = 2
    METABOLISM = 3
    DIGESTION = 4
    TEST = 5
    PERO = 6
    HEMOGLOBIN = 7
    FOX = 8
    OLFACTORY = 9

class GenePurpose(Enum):
    CONTRIBUTER = 1
    CONTROLLER = 2
    UNKNOWN = 3

class Gene:
    def __init__(self, name, description, sequence, geneType, purpose):
        self.name = name
        self.description = description
        self.sequence = sequence
        self.geneType = GeneType(geneType)
        self.purpose = GenePurpose(purpose)
        
        #convert to one hot encoding
        self.numericSequence = self.convertSequenceToNumeric(sequence)
        #print(self.numericSequence)

        #convert to 1234 vector representation
        self.vectorSequence = self.convertSequenceToNumberVector(sequence)
        
        #compute SVD, store eigenNucleotides
        self.eigenNucleotide, self.s, vh = np.linalg.svd(self.numericSequence, full_matrices=True)   
        print("Adding gene: ", self.name, " with sequence length: ", str(len(self.sequence)), " and eigen nucleotide length: ", str(len(self.eigenNucleotide)))
        if(name == 'AAAA' or name == 'ACTG'):
            print("Name: " + str(name) + " eigennucleotide: ")
            print(self.eigenNucleotide)
            print("v")
            print(vh) 
        print("S:")
        print(self.s)
        self.singularRatio = self.s[0]/self.s[-1]
        self.ARatio = float(self.sequence.count('A'))/float(len(self.sequence))
        self.TRatio = float(self.sequence.count('T'))/float(len(self.sequence))
        self.CRatio = float(self.sequence.count('C'))/float(len(self.sequence))
        self.GRatio = float(self.sequence.count('G'))/float(len(self.sequence))
        self.geneRatios = [self.ARatio, self.TRatio, self.CRatio, self.GRatio]
        mn = min(self.ARatio,min(self.TRatio,min(self.CRatio,self.GRatio)))
        mx = max(self.ARatio,max(self.TRatio,max(self.CRatio,self.GRatio)))
        self.nucDiff = mn/mx
        print("Singular ratio: " + str(self.singularRatio))
        print("Ratios: A(" + str(self.ARatio) + ") T(" + str(self.TRatio) + ") C(" + str(self.CRatio) + ") G(" + str(self.GRatio) + ")")
        
    def convertSequenceToNumeric(self, sequence):
        numericSequence = np.zeros((len(sequence),4), dtype=float)
        for i in range(len(sequence)):
            nucleotide = sequence[i]
            if nucleotide == 'A':
                numericSequence[i,0] = 1.0
            elif nucleotide == 'T':
                numericSequence[i,1] = 1.0
            elif nucleotide == 'C':
                numericSequence[i,2] = 1.0
            elif nucleotide == 'G':
                numericSequence[i,3] = 1.0
            else:
                print("Error! Nucleotide does not match known nucleotides!", str(nucleotide))
        return numericSequence

    def convertSequenceToNumberVector(self, sequence):
        numericSequence = np.zeros((len(sequence),1), dtype=float)
        for i in range(len(sequence)):
            nucleotide = sequence[i]
            if nucleotide == 'A':
                numericSequence[i] = 64.0
            elif nucleotide == 'T':
                numericSequence[i] = 128.0
            elif nucleotide == 'C':
                numericSequence[i] = 192.0
            elif nucleotide == 'G':
                numericSequence[i] = 255.0
            else:
                print("Error! Nucleotide does not match known nucleotides!", str(nucleotide))
        return numericSequence
