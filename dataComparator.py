import numpy as np
import geneAgglomerator
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

class DataComparator:
    def __init__(self, agglomerator):
        self.agglomerator = agglomerator

    def dotNumericRepresentations(self):
        hemoglobinEG = self.agglomerator.hemoglobinGenes.numU[:,0]
        peroEG = self.agglomerator.peroGenes.numU[:,0]
        olfactoryEG = self.agglomerator.olfactoryGenes.numU[:,0]
        foxEG = self.agglomerator.foxGenes.numU[:,0]
        self.dotGenes(hemoglobinEG, peroEG, olfactoryEG, foxEG, 'Numeric Eigen gene dot products')        

    def dotTestGenes(self):
        testDF = self.getDotProductDataFrame(self.agglomerator.testGenes.genes, self.agglomerator.testGenes.genes)
        ax = sns.heatmap(testDF, vmin=0.0, vmax=1.0)
        plt.title('Test dot products')
        plt.show()

    def dotAgainstOneself(self):
        hemoglobinDF = self.getDotProductDataFrame(self.agglomerator.hemoglobinGenes.genes, self.agglomerator.hemoglobinGenes.genes)
        peroDF = self.getDotProductDataFrame(self.agglomerator.peroGenes.genes, self.agglomerator.peroGenes.genes)
        foxDF = self.getDotProductDataFrame(self.agglomerator.foxGenes.genes, self.agglomerator.foxGenes.genes)
        olfactoryDF = self.getDotProductDataFrame(self.agglomerator.olfactoryGenes.genes, self.agglomerator.olfactoryGenes.genes)

        fig, ax = plt.subplots(2,2)
        fig.suptitle('Eigen Nucleotides Against Itself')
        g1 = sns.heatmap(hemoglobinDF, ax = ax[0,0], vmin=0.0, vmax=1.0)
        ax[0,0].set_title('Hemoglobin')
        g2 = sns.heatmap(peroDF, ax = ax[0,1], vmin=0.0, vmax=1.0)
        ax[0,1].set_title('Pero')
        g3 = sns.heatmap(foxDF, ax = ax[1,0], vmin=0.0, vmax=1.0)
        ax[1,0].set_title('Fox')
        g4 = sns.heatmap(olfactoryDF, ax = ax[1,1], vmin=0.0, vmax=1.0)
        ax[1,1].set_title('Olfactory')
        plt.show()

    def dotEigenGenes(self):
        hemoglobinEG = self.agglomerator.hemoglobinGenes.getStrongestEigenGene()
        peroEG = self.agglomerator.peroGenes.getStrongestEigenGene()
        olfactoryEG = self.agglomerator.olfactoryGenes.getStrongestEigenGene()
        foxEG = self.agglomerator.foxGenes.getStrongestEigenGene()
        self.dotGenes(hemoglobinEG, peroEG, olfactoryEG, foxEG, 'One Hot Eigen gene dot products')

    def dotGenes(self, hemoglobinEG, peroEG, olfactoryEG, foxEG, title):
        #dot each combination and plot
        labels = ["Hemoglobin", "Pero", "Fox", "Olfactory"]
        eigenDots = np.zeros((4,4))
        eigenDots[0,0] = np.dot(hemoglobinEG, hemoglobinEG)
        minLength = min(len(hemoglobinEG), len(peroEG))
        eigenDots[0,1] = np.dot(hemoglobinEG[:minLength], peroEG[:minLength])
        minLength = min(len(hemoglobinEG), len(foxEG))
        eigenDots[0,2] = np.dot(hemoglobinEG[:minLength], foxEG[:minLength])
        minLength = min(len(hemoglobinEG), len(olfactoryEG))
        eigenDots[0,3] = np.dot(hemoglobinEG[:minLength], olfactoryEG[:minLength])
        eigenDots[1,0] = eigenDots[0,1]
        eigenDots[1,1] = np.dot(peroEG, peroEG)
        minLength = min(len(peroEG), len(foxEG))
        eigenDots[1,2] = np.dot(peroEG[:minLength], foxEG[:minLength])
        minLength = min(len(olfactoryEG), len(peroEG))
        eigenDots[1,3] = np.dot(peroEG[:minLength], olfactoryEG[:minLength])
        eigenDots[2,0] = eigenDots[0,2]
        eigenDots[2,1] = eigenDots[1,2]
        eigenDots[2,2] = np.dot(foxEG, foxEG)
        minLength = min(len(foxEG), len(olfactoryEG))
        eigenDots[2,3] = np.dot(foxEG[:minLength], olfactoryEG[:minLength])
        eigenDots[3,0] = eigenDots[0,3]
        eigenDots[3,1] = eigenDots[1,3]
        eigenDots[3,2] = eigenDots[2,3]
        eigenDots[3,3] = np.dot(olfactoryEG, olfactoryEG)

        eigenDF = pd.DataFrame(eigenDots, columns = labels, index = labels)
        ax = sns.heatmap(eigenDF, vmin=0.0, vmax=1.0)
        plt.title(title)
        plt.show()

    def dotAllGenePairs(self):
        #GENERATE DATAFRAMES FOR ALL DOT PRODUCTS BETWEEN BIGGEST EIGENNUCLEOTIDE
        #generate dataframe for hemoglobin and pero
        self.hemoglobinPeroNucEigDF = self.getDotProductDataFrame( self.agglomerator.hemoglobinGenes.genes, self.agglomerator.peroGenes.genes)
        
        #generate dataframe for hemoglobin and olfactory
        self.hemoglobinOlfactoryEigNucDF = self.getDotProductDataFrame( self.agglomerator.hemoglobinGenes.genes, self.agglomerator.olfactoryGenes.genes)
        
        #generate dataframe for hemoglobin and fox
        self.hemoglobinFoxEigNucDF = self.getDotProductDataFrame( self.agglomerator.hemoglobinGenes.genes, self.agglomerator.foxGenes.genes)

        #generate dataframe for pero and olfactory
        self.peroOlfactoryEigNucDF = self.getDotProductDataFrame( self.agglomerator.peroGenes.genes, self.agglomerator.olfactoryGenes.genes)

        #generate dataframe for pero and fox
        self.peroFoxEigNucDF = self.getDotProductDataFrame( self.agglomerator.peroGenes.genes, self.agglomerator.foxGenes.genes)

        #generate dataframe for olfactory and fox
        self.olfactoryFoxEigNucDF = self.getDotProductDataFrame( self.agglomerator.olfactoryGenes.genes, self.agglomerator.foxGenes.genes)

    def getDotProductDataFrame(self, geneListA, geneListB):
        data = np.zeros((len(geneListA), len(geneListB)))
        columnNames = []
        rowNames = []
        rowCount = 0
        for aGene in geneListA:
            rowNames.append(aGene.name)
            columnCount = 0
            for bGene in geneListB:
                if rowCount == 0:
                    columnNames.append(bGene.name)
                minLength = min(len(aGene.eigenNucleotide[:,0]), len(bGene.eigenNucleotide[:,0]))
                dp = np.dot(aGene.eigenNucleotide[:minLength,0], bGene.eigenNucleotide[:minLength,0])
                data[rowCount, columnCount] = dp
                columnCount = columnCount + 1
            rowCount = rowCount + 1
        df = pd.DataFrame(data, columns=columnNames, index=rowNames)   
        return df

    def plotDottedData(self):       
        #make 2 x 3 plots showing eigen nucleotide relationships between genes of different types
        fig, ax = plt.subplots(2,3)
        fig.suptitle('Eigen Nucleotides')
        g1 = sns.heatmap(self.hemoglobinPeroNucEigDF, ax = ax[0,0], vmin=0.0, vmax=1.0)
        ax[0,0].set_title('Hemoglobin vs Pero')
        g2 = sns.heatmap(self.hemoglobinOlfactoryEigNucDF, ax = ax[0,1], vmin=0.0, vmax=1.0)
        ax[0,1].set_title('Hemoglobin vs Olfactory')
        g3 = sns.heatmap(self.hemoglobinFoxEigNucDF, ax = ax[0,2], vmin=0.0, vmax=1.0)
        ax[0,2].set_title('Hemoglobin vs Fox')
        g4 = sns.heatmap(self.peroOlfactoryEigNucDF, ax = ax[1,0], vmin=0.0, vmax=1.0)
        ax[1,0].set_title('Pero vs Olfactory')
        g5 = sns.heatmap(self.peroFoxEigNucDF, ax = ax[1,1], vmin=0.0, vmax=1.0)
        ax[1,1].set_title('Pero vs Fox')
        g6 = sns.heatmap(self.olfactoryFoxEigNucDF, ax = ax[1,2], vmin=0.0, vmax=1.0)
        ax[1,2].set_title('Olfactory vs Fox')
        plt.show()
