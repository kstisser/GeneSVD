import gene
import numpy as np

class TestGenes():
    def __init__(self):
        self.genes = []
        st = "ACTG"
        length = len(st)
        data = [""] * (length)
        #self.addAllTestGenes(st, data, (length - 1), 0)
        self.simpleAddAll()

    def getStrongestEigenGene(self):
        eigenGene, s, vh = np.linalg.svd(self.allOneHot, full_matrices=True)
        return eigenGene[:,0]

    def addGene(self, name, description, purpose):
        self.genes.append(gene.Gene(name, description, name, gene.GeneType.TEST, purpose))

    def addAllTestGenes(self, st, data, last, index):
        length= len(st)
        for i in xrange(length):
            data[index] = st[i]
            if index == last:
                sequence = ''.join(data)
                self.addGene(sequence, "", gene.GenePurpose.UNKNOWN)
            else:
                self.addAllTestGenes(st, data, last, index+1)

    def simpleAddAll(self):
        self.addGene("AAAA", "", gene.GenePurpose.UNKNOWN)
        self.addGene("CCCC", "", gene.GenePurpose.UNKNOWN)
        self.addGene("TTTT", "", gene.GenePurpose.UNKNOWN)
        self.addGene("GGGG", "", gene.GenePurpose.UNKNOWN)
        self.addGene("CCCA", "", gene.GenePurpose.UNKNOWN)
        self.addGene("TTTA", "", gene.GenePurpose.UNKNOWN)
        self.addGene("GGGA", "", gene.GenePurpose.UNKNOWN)
        self.addGene("AAAT", "", gene.GenePurpose.UNKNOWN)
        self.addGene("CCCT", "", gene.GenePurpose.UNKNOWN)
        self.addGene("GGGT", "", gene.GenePurpose.UNKNOWN)
        self.addGene("AAAC", "", gene.GenePurpose.UNKNOWN)
        self.addGene("TTTC", "", gene.GenePurpose.UNKNOWN)
        self.addGene("GGGC", "", gene.GenePurpose.UNKNOWN)
        self.addGene("AAAG", "", gene.GenePurpose.UNKNOWN)
        self.addGene("CCCG", "", gene.GenePurpose.UNKNOWN)
        self.addGene("TTTG", "", gene.GenePurpose.UNKNOWN)
        self.addGene("ACCC", "", gene.GenePurpose.UNKNOWN)
        self.addGene("ATTT", "", gene.GenePurpose.UNKNOWN)
        self.addGene("AGGG", "", gene.GenePurpose.UNKNOWN)
        self.addGene("TGGG", "", gene.GenePurpose.UNKNOWN)        
        self.addGene("ATAT", "", gene.GenePurpose.UNKNOWN)
        self.addGene("CGCG", "", gene.GenePurpose.UNKNOWN)
        self.addGene("AACC", "", gene.GenePurpose.UNKNOWN)
        self.addGene("TTGG", "", gene.GenePurpose.UNKNOWN)
        self.addGene("ACTG", "", gene.GenePurpose.UNKNOWN)
        self.addGene("GTCA", "", gene.GenePurpose.UNKNOWN)
 

