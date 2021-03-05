import gene
import numpy as np

class TestGenes():
    def __init__(self):
        self.genes = []
        self.addAllTestGenes()

    def getStrongestEigenGene(self):
        eigenGene, s, vh = np.linalg.svd(self.allOneHot, full_matrices=True)
        return eigenGene[:,0]

    def addGene(self, name, description, purpose, sequence):
        self.genes.append(gene.Gene(name, description, sequence, gene.GeneType.TEST, purpose))

    def addAllTestGenes(self):
        self.addGene("AAAA", "", gene.GenePurpose.UNKNOWN, "AAAA")
        self.addGene("CCCC", "", gene.GenePurpose.UNKNOWN, "CCCC")
        self.addGene("TTTT", "", gene.GenePurpose.UNKNOWN, "TTTT")
        self.addGene("GGGG", "", gene.GenePurpose.UNKNOWN, "GGGG")
        self.addGene("ATAT", "", gene.GenePurpose.UNKNOWN, "ATAT")
        self.addGene("CGCG", "", gene.GenePurpose.UNKNOWN, "CGCG")
        self.addGene("AACC", "", gene.GenePurpose.UNKNOWN, "AACC")
        self.addGene("TTGG", "", gene.GenePurpose.UNKNOWN, "TTGG")
        self.addGene("ACTG", "", gene.GenePurpose.UNKNOWN, "ACTG")
        self.addGene("GTCA", "", gene.GenePurpose.UNKNOWN, "GTCA")
