import gene
import numpy as np

class GeneHolder():
    def __init__(self):
        self.genes = []
        self.maxNucleotideLength = 0
        self.allOneHot = np.empty((0,4))
        self.addAllGenes()
        self.compileNumberedNucleotides()

    def getStrongestEigenGene(self):
        eigenGene, s, vh = np.linalg.svd(self.allOneHot, full_matrices=True)
        return eigenGene[:,0]

    def addGene(self, name, description, purpose, sequence, geneType):
        self.genes.append(gene.Gene(name, description, sequence, geneType, purpose))
        self.allOneHot = np.vstack((self.allOneHot, self.genes[-1].numericSequence))
        if(len(sequence) > self.maxNucleotideLength):
            self.maxNucleotideLength = len(sequence)

    def compileNumberedNucleotides(self, printout):
        self.allNucleotides = np.full((self.maxNucleotideLength, len(self.genes)),0)
        count = 0
        for g in self.genes:
            print(g.vectorSequence.shape)
            self.allNucleotides[0:len(g.vectorSequence),count] = (g.vectorSequence).reshape(len(g.vectorSequence))
            count = count + 1
        [self.numU, self.numS, self.numV] = np.linalg.svd(self.allNucleotides, full_matrices=False)
        print(printout)
        print(self.numS)    
