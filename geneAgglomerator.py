import numpy as np
import pandas as pd
import testGenes
import eyeGenes
import hairGenes
import metabolismGenes
import digestionGenes
import foxTranscriptionGenes
import olfactoryReceptorGenes
import hemoglobinGenes
import peroxiredoxinGenes
import gene
import enum

#This class collects all of the data in one location
class GeneAgglomerator:
    def __init__(self):
        #cols = ['Name', 'GeneType', 'GenePurpose', 'Description', 'Sequence', 'NumericSequence', 'EigenNucleotides']
        self.testGenes = testGenes.TestGenes()
        self.hemoglobinGenes = hemoglobinGenes.HemoglobinGenes()
        self.peroGenes = peroxiredoxinGenes.PeroGenes()
        self.foxGenes = foxTranscriptionGenes.FoxTranscriptionGenes()
        self.olfactoryGenes = olfactoryReceptorGenes.OlfactoryGenes()
        self.allData = np.concatenate((self.hemoglobinGenes.genes, self.peroGenes.genes, self.foxGenes.genes, self.olfactoryGenes.genes))


        contributerCount = 0
        controllerCount = 0
        for g in self.allData:
            #count contributer vs controller
            if (g.purpose) == (gene.GenePurpose.CONTRIBUTER):
                contributerCount = contributerCount + 1
            elif (g.purpose) == (gene.GenePurpose.CONTROLLER):
                controllerCount = controllerCount + 1

        print("#################SUMMARY###########################")
        print("We have ", str(contributerCount), " known contributer genes")
        print("We have ", str(controllerCount), " known controller genes")
        print("###################################################")
        print("We have ", len(self.allData), " total genes to analyze")
        print("###################################################")
        print("We have ", len(self.testGenes.genes), " total test genes")
        print("We have ", len(self.foxGenes.genes), " total Fox genes")
        print("We have ", len(self.hemoglobinGenes.genes), " total Hemoglobin genes")
        print("We have ", len(self.peroGenes.genes), " total Pero genes")
        print("We have ", len(self.olfactoryGenes.genes), " total Olfactory genes")
        print("###################################################")



