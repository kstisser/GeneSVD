import gene
import digestionGenes
import eyeGenes
import hairGenes
import metabolismGenes
import foxTranscriptionGenes
import olfactoryReceptorGenes
import hemoglobinGenes
import peroxiredoxinGenes
import geneAgglomerator as ga
import dataComparator as dc
import statsFinder as sf


class GeneManager:
    def __init__(self):
        print("Starting the Gene Manager, compiling all data")
        self.compileAllData()

    def compileAllData(self):
        geneAgg = ga.GeneAgglomerator()
        stats = sf.StatsFinder(geneAgg)
        stats.plotSingularValueRatiosByNumNucleotides()
        stats.plotLetterRatiosByCategory()
        stats.showSingularValuesOfAll()
        dataComparator = dc.DataComparator(geneAgg)
        dataComparator.dotTestGenes()
        dataComparator.dotAllGenePairs()
        dataComparator.plotDottedData()
        dataComparator.dotAgainstOneself()
        dataComparator.dotEigenGenes()
        dataComparator.dotNumericRepresentations()               


if __name__ == "__main__":
    gm = GeneManager()
