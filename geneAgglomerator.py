import numpy as np
import pandas as pd
import eyeGenes
import hairGenes
import metabolismGenes
import digestionGenes
import gene
import enum

class GeneAgglomerator:
    def __init__(self):
        #cols = ['Name', 'GeneType', 'GenePurpose', 'Description', 'Sequence', 'NumericSequence', 'EigenNucleotides']
        self.eyeGenes = eyeGenes.EyeGenes()
        self.hairGenes = hairGenes.HairGenes()
        self.digestionGenes = digestionGenes.DigestionGenes()
        self.metabolismGenes = metabolismGenes.MetabolismGenes()
        allData = np.concatenate((self.eyeGenes.genes, self.hairGenes.genes, self.digestionGenes.genes, self.metabolismGenes.genes))


        contributerCount = 0
        controllerCount = 0
        for g in allData:
            #count contributer vs controller
            if (g.purpose) == (gene.GenePurpose.CONTRIBUTER):
                contributerCount = contributerCount + 1
            elif (g.purpose) == (gene.GenePurpose.CONTROLLER):
                controllerCount = controllerCount + 1

        print("#################SUMMARY###########################")
        print("We have ", str(contributerCount), " known contributer genes")
        print("We have ", str(controllerCount), " known controller genes")
        print("###################################################")
        print("We have ", len(allData), " total genes to analyze")
        print("###################################################")



