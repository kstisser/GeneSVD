import geneHolder
import gene

#Genes affecting smell
class OlfactoryGenes(geneHolder.GeneHolder, object):
    def __init__(self):
        super(OlfactoryGenes, self).__init__()      

    def addGene(self, name, description, purpose, sequence):
        super(OlfactoryGenes, self).addGene(name, description, purpose, sequence, gene.GeneType.OLFACTORY)

    def compileNumberedNucleotides(self):
        super(OlfactoryGenes, self).compileNumberedNucleotides("Olfactory singular values:")                                     

    def addAllGenes(self):
        self.addGene("OR10C1", "Smell gene- 2-Ethylfenchol", gene.GenePurpose.CONTRIBUTER, "ATGAGTGCAAACACCTCCATGGTGACTGAGTTTCTTCTTCTCGGCTTCTCCCACCTGGCCGACCTCCAGG\
GCTTGCTCTTCTCTGTCTTTCTCACTATCTACCTGCTGACCGTGGCAGGCAATTTCCTCATTGTGGTGCT\
GGTCTCCACTGATGCTGCCCTCCAGTCCCCTATGTACTTCTTCCTGCGCACCCTCTCGGCCTTGGAGATT\
GGCTATACGTCTGTCACGGTCCCCCTGCTACTTCACCACCTCCTTACTGGCCGGCGCCACATCTCTCGCT\
CTGGATGTGCTCTCCAGATGTTCTTCTTCCTCTTCTTTGGCGCCACGGAGTGCTGCCTCCTGGCAGCCAT\
GGCCTATGACCGCTATGCAGCCATCTGTGAACCCCTCCGCTACCCACTGCTGCTGAGCCACCGGGTGTGT\
CTACAGCTAGCTGGGTCGGCGTGGGCCTGTGGGGTGCTGGTGGGGCTGGGCCACACCCCTTTCATCTTCT\
CTTTGCCCTTCTGCGGCCCCAATACCATCCCGCAGTTCTTCTGTGAGATCCAGCCTGTCCTGCAGCTGGT\
ATGTGGAGACACCTCGCTTAATGAACTGCAGATTATCCTGGCAACAGCCCTCCTCATCCTCTGCCCCTTT\
GGCCTCATCCTGGGCTCCTACGGGCGTATCCTCGTTACCATCTTCCGGATCCCATCTGTTGCGGGCCGCC\
GCAAGGCCTTCTCCACCTGCTCCTCCCACCTGATCATGGTCTCCCTCTTCTATGGCACCGCACTCTTTAT\
CTATATTCGCCCTAAGGCCAGCTACGATCCGGCCACTGACCCTCTGGTGTCCCTCTTCTATGCTGTGGTC\
ACCCCCATCCTCAACCCCATCATCTACAGCCTGCGGAACACAGAGGTCAAAGCTGCCCTAAAGAGAACCA\
TCCAGAAAACGGTGCCTATGGAGATTTGA")
        self.addGene("OR11A1", "Smell gene- 2-Ethylfenchol", gene.GenePurpose.CONTRIBUTER, "ATGGAAATTGTCTCCACAGGAAACGAAACTATTACTGAATTTGTCCTCCTTGGCTTCTATGACATCCCTG\
AACTGCATTTCTTGTTTTTTATTGTATTCACTGCTGTCTATGTCTTCATCATCATAGGGAATATGCTGAT\
TATTGTAGCAGTGGTTAGCTCCCAGAGGCTCCACAAACCCATGTATATTTTCTTGGCGAATCTGTCCTTC\
CTGGATATTCTCTACACCTCCGCAGTGATGCCAAAAATGCTGGAGGGCTTCCTGCAAGAAGCAACTATCT\
CTGTGGCTGGTTGCTTGCTCCAGTTCTTTATCTTCGGCTCTCTAGCCACAGCTGAATGCTTACTGCTGGC\
TGTCATGGCATATGACCGCTACCTGGCAATTTGCTACCCACTCCACTACCCACTCCTGATGGGGCCCAGA\
CGGTACATGGGGCTGGTGGTCACAACCTGGCTCTCTGGATTTGTGGTAGATGGACTGGTTGTGGCCCTGG\
TGGCCCAGCTGAGGTTCTGTGGCCCCAACCACATTGACCAGTTTTACTGTGACTTTATGCTTTTCGTGGG\
CCTGGCTTGCTCGGATCCCAGAGTGGCTCAGGTGACAACTCTCATTCTGTCTGTGTTCTGCCTCACTATT\
CCTTTTGGACTGATTCTGACATCTTATGCCAGAATTGTGGTGGCAGTGCTGAGAGTTCCTGCTGGGGCAA\
GCAGGAGAAGGGCTTTCTCCACATGCTCCTCCCACCTAGCTGTAGTGACCACATTCTATGGAACGCTCAT\
GATCTTTTATGTTGCACCCTCTGCTGTCCATTCCCAGCTCCTCTCCAAGGTCTTCTCCCTGCTCTACACT\
GTGGTCACCCCTCTCTTCAATCCTGTGATCTATACCATGAGGAACAAGGAGGTGCATCAGGCACTTCGGA\
AGATTCTCTGTATCAAACAAACTGAAACACTTGATTGA")
        self.addGene("OR10R2", "Smell gene- Diacetyl", gene.GenePurpose.CONTRIBUTER, "ATGCCCCAAATTCTTATATTCACATACCTGAATATGTTTTACTTCTTTCCCCCTTTGCAGATCTTGGCAG\
AAAACCTCACCATGGTCACCGAATTCCTGTTGCTGGGTTTTTCCAGCCTTGGTGAAATTCAGCTGGCCCT\
CTTTGTAGTTTTTCTTTTTCTGTATCTAGTCATTCTTAGTGGCAATGTCACCATTATCAGTGTCATCCAC\
CTGGATAAAAGCCTCCACACACCAATGTACTTCTTCCTTGGCATTCTCTCAACATCTGAGACCTTCTACA\
CCTTTGTCATTCTACCCAAGATGCTCATCAATCTACTTTCTGTGGCCAGGACAATCTCCTTCAACTGTTG\
TGCTCTTCAAATGTTCTTCTTCCTTGGTTTTGCCATTACCAACTGCCTGCTATTGGGTGTGATGGGTTAT\
GATCGCTATGCTGCCATTTGTCACCCTCTGCATTACCCCACTCTTATGAGCTGGCAGGTGTGTGGAAAAC\
TGGCAGCTGCCTGTGCAATTGGTGGCTTCTTGGCCTCTCTTACAGTAGTAAATTTAGTTTTCAGCCTCCC\
TTTTTGTAGCGCCAACAAAGTCAATCATTACTTCTGTGACATCTCAGCAGTCATTCTTCTGGCTTGTACC\
AACACAGATGTTAACGAATTTGTGATATTCATTTGTGGAGTTCTTGTACTTGTGGTTCCCTTTCTGTTTA\
TCTGTGTTTCTTATCTCTGCATTCTGAGGACTATCCTGAAGATTCCCTCAGCTGAGGGCAGACGGAAAGC\
GTTTTCCACCTGCGCCTCTCACCTCAGTGTTGTTATTGTTCATTATGGCTGTGCTTCCTTCATCTACCTG\
AGGCCTACAGCAAACTATGTGTCCAACAAAGACAGGCTGGTGACGGTGACATACACGATTGTCACTCCAT\
TACTAAACCCCATGGTTTATAGCCTCAGAAACAAGGATGTCCAACTTGCTATCAGAAAAGTGTTGGGCAA\
GAAAGGTTCTCTAAAACTATATAATTGA")
        self.addGene("OR10Z1", "Smell gene- Diacetyl", gene.GenePurpose.CONTRIBUTER, "ATGGGGCAGACCAACGTAACCTCCTGGAGGGATTTTGTCTTCCTGGGCTTCTCCAGTTCTGGGGAGTTGC\
AGCTCCTTCTCTTTGCCTTGTTCCTCTCTCTGTATCTAGTCACTCTGACCAGCAATGTCTTCATTATCAT\
AGCCATCAGGCTGGATAGCCATCTGCACACCCCCATGTACCTCTTCCTTTCCTTCCTATCCTTCTCTGAG\
ACCTGCTACACTTTGGGCATCATCCCTAGAATGCTCTCTGGCCTGGCTGGGGGGGACCAGGCTATCTCCT\
ATGTGGGCTGTGCTGCCCAGATGTTCTTTTCTGCCTCATGGGCCTGTACTAACTGCTTCCTTCTGGCTGC\
CATGGGCTTTGACAGATATGTGGCCATCTGTGCTCCACTCCACTATGCCAGCCACATGAATCCTACCCTC\
TGTGCCCAGCTGGTCATTACTTCCTTCCTGACTGGATACCTCTTTGGACTGGGAATGACACTAGTTATTT\
TCCACCTCTCATTCTGCAGCTCCCATGAAATCCAGCACTTTTTTTGTGACACGCCACCTGTGCTGAGCCT\
AGCCTGTGGAGATACAGGCCCGAGTGAGCTGAGGATCTTTATCCTCAGTCTTTTGGTCCTCTTGGTCTCC\
TTCTTCTTCATCACCATCTCCTACGCCTACATCTTGGCAGCAATACTGAGGATCCCCTCTGCTGAGGGGC\
AGAAGAAGGCCTTCTCCACTTGTGCCTCGCACCTTACAGTGGTCATTATTCATTATGGCTGTGCTTCCTT\
CGTGTACCTGAGGCCCAAAGCCAGCTACTCTCTTGAGAGAGATCAGCTTATTGCCATGACCTATACTGTA\
GTGACCCCCCTCCTTAATCCCATTGTTTATAGTCTAAGGAATAGGGCTATACAGACAGCTCTGAGGAATG\
CTTTCAGAGGGAGATTGCTGGGTAAAGGATGA")
        self.addGene("OR6Y1", "Smell gene- Diacetyl", gene.GenePurpose.CONTRIBUTER, "ATGACCACCATAATTCTGGAAGTAGATAATCATACAGTGACAACACGTTTCATTCTTCTGGGGTTTCCAA\
CACGACCAGCCTTCCAGCTTCTCTTTTTCTCCATTTTCCTGGCAACCTATCTGCTGACACTGCTGGAGAA\
TCTTCTTATCATCTTAGCTATCCACAGTGATGGGCAGCTGCATAAGCCCATGTACTTCTTCTTGAGCCAC\
CTCTCCTTCCTGGAGATGTGGTATGTCACAGTCATCAGCCCCAAGATGCTTGTTGACTTCCTCAGTCATG\
ACAAGAGTATTTCCTTCAATGGCTGCATGACTCAACTTTACTTTTTTGTGACCTTTGTCTGCACTGAGTA\
CATCCTTCTTGCTATCATGGCCTTTGACCGCTATGTAGCCATTTGTAATCCACTACGCTACCCAGTCATC\
ATGACCAACCAGCTCTGTGGCACACTGGCTGGAGGATGCTGGTTCTGTGGACTCATGACTGCCATGATTA\
AGATGGTTTTTATAGCACAACTTCACTACTGTGGCATGCCTCAGATCAATCACTACTTTTGTGATATCTC\
TCCACTCCTTAACGTCTCCTGTGAGGATGCCTCACAGGCTGAGATGGTGGACTTCTTCTTGGCCCTCATG\
GTCATTGCTATTCCTCTTTGTGTTGTGGTGGCATCCTACGCTGCTATCCTTGCCACCATCCTCAGGATCC\
CTTCTGCTCAGGGCCGCCAAAAGGCATTCTCCACCTGTGCCTCCCACCTGACCGTCGTAATTCTCTTCTA\
TTCCATGACACTTTTCACCTATGCCCGTCCCAAACTCATGTATGCCTACAATTCCAACAAAGTGGTATCT\
GTTCTCTACACTGTCATTGTTCCACTCCTCAACCCCATCATTTACTGTCTGAGGAACCATGAAGTAAAGG\
CAGCCCTCAGAAAGACCATACATTGCAGAGGAAGTGGGCCCCAGGGAAATGGGGCTTTCAGTAGTTAA")


