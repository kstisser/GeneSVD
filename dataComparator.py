import numpy as np
import geneAgglomerator
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

class DataComparator:
    def __init__(self, agglomerator):
        self.agglomerator = agglomerator

    def dotEigenGenes(self):
        eyeEG = self.agglomerator.eyeGenes.getStrongestEigenGene()
        hairEG = self.agglomerator.hairGenes.getStrongestEigenGene()
        digestionEG = self.agglomerator.digestionGenes.getStrongestEigenGene()
        metabolismEG = self.agglomerator.metabolismGenes.getStrongestEigenGene()

        #dot each combination and plot
        labels = ["Eye", "Hair", "Digestion", "Metabolism"]
        eigenDots = np.zeros((4,4))
        eigenDots[0,0] = np.dot(eyeEG, eyeEG)
        minLength = min(len(eyeEG), len(hairEG))
        eigenDots[0,1] = np.dot(eyeEG[:minLength], hairEG[:minLength])
        minLength = min(len(eyeEG), len(digestionEG))
        eigenDots[0,2] = np.dot(eyeEG[:minLength], digestionEG[:minLength])
        minLength = min(len(eyeEG), len(metabolismEG))
        eigenDots[0,3] = np.dot(eyeEG[:minLength], metabolismEG[:minLength])
        eigenDots[1,0] = eigenDots[0,1]
        eigenDots[1,1] = np.dot(hairEG, hairEG)
        minLength = min(len(hairEG), len(digestionEG))
        eigenDots[1,2] = np.dot(hairEG[:minLength], digestionEG[:minLength])
        minLength = min(len(eyeEG), len(hairEG))
        eigenDots[1,3] = np.dot(hairEG[:minLength], metabolismEG[:minLength])
        eigenDots[2,0] = eigenDots[0,2]
        eigenDots[2,1] = eigenDots[1,2]
        eigenDots[2,2] = np.dot(digestionEG, digestionEG)
        minLength = min(len(digestionEG), len(metabolismEG))
        eigenDots[2,3] = np.dot(digestionEG[:minLength], metabolismEG[:minLength])
        eigenDots[3,0] = eigenDots[0,3]
        eigenDots[3,1] = eigenDots[1,3]
        eigenDots[3,2] = eigenDots[2,3]
        eigenDots[3,3] = np.dot(metabolismEG, metabolismEG)

        eigenDF = pd.DataFrame(eigenDots, columns = labels, index = labels)
        ax = sns.heatmap(eigenDF, vmin=0.0, vmax=1.0)
        plt.title('Eigen gene dot products')
        plt.show()

    def dotAllGenePairs(self):
        #GENERATE DATAFRAMES FOR ALL DOT PRODUCTS BETWEEN BIGGEST EIGENNUCLEOTIDE
        #generate dataframe for eyes and hair
        self.eyeHairNucEigDF = self.getDotProductDataFrame( self.agglomerator.eyeGenes.genes, self.agglomerator.hairGenes.genes)
        
        #generate dataframe for eyes and metabolism
        self.eyeMetabolismEigNucDF = self.getDotProductDataFrame( self.agglomerator.eyeGenes.genes, self.agglomerator.metabolismGenes.genes)
        
        #generate dataframe for eyes and digestion
        self.eyeDigestionEigNucDF = self.getDotProductDataFrame( self.agglomerator.eyeGenes.genes, self.agglomerator.digestionGenes.genes)

        #generate dataframe for hair and metabolism
        self.hairMetabolismEigNucDF = self.getDotProductDataFrame( self.agglomerator.hairGenes.genes, self.agglomerator.metabolismGenes.genes)

        #generate dataframe for hair and digestion
        self.hairDigestionEigNucDF = self.getDotProductDataFrame( self.agglomerator.hairGenes.genes, self.agglomerator.digestionGenes.genes)

        #generate dataframe for metabolism and digestion
        self.metabolismDigestionEigNucDF = self.getDotProductDataFrame( self.agglomerator.metabolismGenes.genes, self.agglomerator.digestionGenes.genes)

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
                #print((aGene.eigenNucleotide[0]))
                minLength = min(len(aGene.eigenNucleotide[:,0]), len(bGene.eigenNucleotide[:,0]))
                dp = np.dot(aGene.eigenNucleotide[:minLength,0], bGene.eigenNucleotide[:minLength,0])
                data[rowCount, columnCount] = dp
                columnCount = columnCount + 1
            rowCount = rowCount + 1
        df = pd.DataFrame(data, columns=columnNames, index=rowNames)   
        return df

    def plotDottedData(self):
        fig, ax = plt.subplots(2,3)
        fig.suptitle('Eigen Nucleotides')
        g1 = sns.heatmap(self.eyeHairNucEigDF, ax = ax[0,0])
        ax[0,0].set_title('Eye vs Hair')
        g2 = sns.heatmap(self.eyeMetabolismEigNucDF, ax = ax[0,1])
        ax[0,1].set_title('Eye vs Metabolism')
        g3 = sns.heatmap(self.eyeDigestionEigNucDF, ax = ax[0,2])
        ax[0,2].set_title('Eye vs Digestion')
        g4 = sns.heatmap(self.hairMetabolismEigNucDF, ax = ax[1,0])
        ax[1,0].set_title('Hair vs Metabolism')
        g5 = sns.heatmap(self.hairDigestionEigNucDF, ax = ax[1,1])
        ax[1,1].set_title('Hair vs Digestion')
        g6 = sns.heatmap(self.metabolismDigestionEigNucDF, ax = ax[1,2])
        ax[1,2].set_title('Metabolism vs Digestion')
        plt.show()
        
        #make 2 x 3 plots showing eigen nucleotide relationships between genes of different types
        fig, ax = plt.subplots(2,3)
        fig.suptitle('Eigen Nucleotides')
        g1 = sns.heatmap(self.eyeHairNucEigDF, ax = ax[0,0], vmin=0.0, vmax=1.0)
        ax[0,0].set_title('Eye vs Hair')
        g2 = sns.heatmap(self.eyeMetabolismEigNucDF, ax = ax[0,1], vmin=0.0, vmax=1.0)
        ax[0,1].set_title('Eye vs Metabolism')
        g3 = sns.heatmap(self.eyeDigestionEigNucDF, ax = ax[0,2], vmin=0.0, vmax=1.0)
        ax[0,2].set_title('Eye vs Digestion')
        g4 = sns.heatmap(self.hairMetabolismEigNucDF, ax = ax[1,0], vmin=0.0, vmax=1.0)
        ax[1,0].set_title('Hair vs Metabolism')
        g5 = sns.heatmap(self.hairDigestionEigNucDF, ax = ax[1,1], vmin=0.0, vmax=1.0)
        ax[1,1].set_title('Hair vs Digestion')
        g6 = sns.heatmap(self.metabolismDigestionEigNucDF, ax = ax[1,2], vmin=0.0, vmax=1.0)
        ax[1,2].set_title('Metabolism vs Digestion')
        plt.show()
