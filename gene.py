import numpy as np
from scipy.linalg import svd
from enum import Enum

class GeneType(Enum):
    EYECOLOR = 1
    HAIRCOLOR = 2
    METABOLISM = 3
    DIGESTION = 4

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
        
        #compute SVD, store eigenNucleotides
        self.eigenNucleotide, s, vh = np.linalg.svd(self.numericSequence, full_matrices=True)

        print("Adding gene: ", self.name, " with sequence length: ", str(len(self.sequence)), " and eigen nucleotide length: ", str(len(self.eigenNucleotide)))

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
