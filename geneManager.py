import gene
import digestionGenes
import eyeGenes
import hairGenes
import metabolismGenes
import geneAgglomerator as ga
import dataComparator as dc


class GeneManager:
    def __init__(self):
        print("Starting the Gene Manager, compiling all data")
        self.compileAllData()

    def compileAllData(self):
        geneAgg = ga.GeneAgglomerator()
        dataComparator = dc.DataComparator(geneAgg)
        dataComparator.dotTestGenes()
        dataComparator.dotAllGenePairs()
        dataComparator.plotDottedData()
        dataComparator.dotEigenGenes()
        
        print("Read in ", str(len(geneAgg.metabolismGenes.genes)), " metabolism genes")
        print("Read in ", str(len(geneAgg.digestionGenes.genes)), " digestion genes")
        print("Read in ", str(len(geneAgg.eyeGenes.genes)), " eye genes")
        print("Read in ", str(len(geneAgg.hairGenes.genes)), " hair genes")
                


if __name__ == "__main__":
    gm = GeneManager()
