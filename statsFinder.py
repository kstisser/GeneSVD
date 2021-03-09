import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

class StatsFinder:
    def __init__(self, agg):
        self.data = agg

    def showSingularValuesOfAll(self):
        hemoRows, hemoCols = (self.data.hemoglobinGenes.allNucleotides).shape
        peroRows, peroCols = (self.data.peroGenes.allNucleotides).shape
        olfactoryRows, olfactoryCols = (self.data.olfactoryGenes.allNucleotides).shape
        foxRows, foxCols = (self.data.foxGenes.allNucleotides).shape
        maxLength = max(hemoRows, max(peroRows, max(olfactoryRows, foxRows)))
        combinedNucleotides = np.full((maxLength, (hemoCols + peroCols + olfactoryCols + foxCols)), 0)
        combinedNucleotides[0:hemoRows, 0:hemoCols] = self.data.hemoglobinGenes.allNucleotides
        combinedNucleotides[0:peroRows, hemoCols:(hemoCols + peroCols)] = self.data.peroGenes.allNucleotides
        combinedNucleotides[0:olfactoryRows, (hemoCols + peroCols):(hemoCols + peroCols + olfactoryCols)] = self.data.olfactoryGenes.allNucleotides
        combinedNucleotides[0:foxRows, (hemoCols + peroCols + olfactoryCols):(hemoCols + peroCols + olfactoryCols + foxCols)] = self.data.foxGenes.allNucleotides

        [u, s, v] = np.linalg.svd(combinedNucleotides, full_matrices=True)
        
        '''print("*********************************************************")
        print("Singular values of all data:")
        print(s)
        print("*********************************************************")'''
        s = (s)/255.0

        print("*********************************************************")
        print("Normalized Singular values of all data:")
        print(s)
        print("*********************************************************")
        for g in self.data.hemoglobinGenes.genes:
            plt.plot(g.s, 'r')

        for g in self.data.peroGenes.genes:
            plt.plot(g.s, 'g')

        for g in self.data.foxGenes.genes:
            plt.plot(g.s, 'b')

        for g in self.data.olfactoryGenes.genes:
            plt.plot(g.s, 'c')
        
        plt.plot(s, 'ko')
        blue_patch = mpatches.Patch(color='blue', label='Fox')
        red_patch = mpatches.Patch(color='red', label='Hemoglobin')
        green_patch = mpatches.Patch(color='green', label='Pero')
        cyan_patch = mpatches.Patch(color='cyan', label='Olfactory')
        black_patch = mpatches.Patch(color='black', label='All')
        plt.legend(handles=[blue_patch, red_patch, green_patch, cyan_patch, black_patch])        
        plt.title('Singular Values')
        plt.show()
        

    def plotSingularValueRatiosByNumNucleotides(self):
        self.addSVRatioPlots('bo', self.data.foxGenes.genes)
        self.addSVRatioPlots('ro', self.data.hemoglobinGenes.genes)
        self.addSVRatioPlots('go', self.data.peroGenes.genes)
        self.addSVRatioPlots('co', self.data.olfactoryGenes.genes)
        plt.title("Singular Value Ratio to # Nucleotides")
        plt.xlabel("Singular Value Ratio (min/max)")
        plt.ylabel("Number of nucleotides")
        blue_patch = mpatches.Patch(color='blue', label='Fox')
        red_patch = mpatches.Patch(color='red', label='Hemoglobin')
        green_patch = mpatches.Patch(color='green', label='Pero')
        cyan_patch = mpatches.Patch(color='cyan', label='Olfactory')
        plt.legend(handles=[blue_patch, red_patch, green_patch, cyan_patch])
        plt.show()        

    def plotLetterRatiosByCategory(self):
        fig, ax = plt.subplots(1,2)
        #AT
        ax[0].set_title("A vs T ratio")

        #CG
        ax[1].set_title("C vs G ratio")

        self.addLetterPlots(ax, 'bo', self.data.foxGenes.genes)
        self.addLetterPlots(ax, 'ro', self.data.hemoglobinGenes.genes)
        self.addLetterPlots(ax, 'go', self.data.peroGenes.genes)
        self.addLetterPlots(ax, 'co', self.data.olfactoryGenes.genes)
        blue_patch = mpatches.Patch(color='blue', label='Fox')
        red_patch = mpatches.Patch(color='red', label='Hemoglobin')
        green_patch = mpatches.Patch(color='green', label='Pero')
        cyan_patch = mpatches.Patch(color='cyan', label='Olfactory')
        plt.legend(handles=[blue_patch, red_patch, green_patch, cyan_patch])        
        plt.show()

    def addLetterPlots(self, ax, color, genes):
        for g in genes:
            ax[0].plot(g.ARatio, g.TRatio, color)
            ax[1].plot(g.CRatio, g.GRatio, color)
        
    def addSVRatioPlots(self, color, genes):
        for g in genes:
            plt.plot(g.singularRatio, len(g.sequence), color)
