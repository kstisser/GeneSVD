import geneHolder
import gene

class FoxTranscriptionGenes(geneHolder.GeneHolder, object):
    def __init__(self):
        super(FoxTranscriptionGenes, self).__init__()      

    def addGene(self, name, description, purpose, sequence):
        super(FoxTranscriptionGenes, self).addGene(name, description, purpose, sequence, gene.GeneType.FOX)

    def compileNumberedNucleotides(self):
        super(FoxTranscriptionGenes, self).compileNumberedNucleotides("Fox singular values:")                                     

    def addAllGenes(self):
        self.addGene("FOXA1", "Transcription factor that plays an important role in regulating \
the expression of genes involved in cell growth, proliferation, differentiation, \
and longevity", gene.GenePurpose.CONTROLLER, "ATGTTAGGAACTGTGAAGATGGAAGGGCATGAAACCAGCGACTGGAACAGCTACTACGCAGACACGCAGG\
AGGCCTACTCCTCCGTCCCGGTCAGCAACATGAACTCAGGCCTGGGCTCCATGAACTCCATGAACACCTA\
CATGACCATGAACACCATGACTACGAGCGGCAACATGACCCCGGCGTCCTTCAACATGTCCTATGCCAAC\
CCGGGCCTAGGGGCCGGCCTGAGTCCCGGCGCAGTAGCCGGCATGCCGGGGGGCTCGGCGGGCGCCATGA\
ACAGCATGACTGCGGCCGGCGTGACGGCCATGGGTACGGCGCTGAGCCCGAGCGGCATGGGCGCCATGGG\
TGCGCAGCAGGCGGCCTCCATGAATGGCCTGGGCCCCTACGCGGCCGCCATGAACCCGTGCATGAGCCCC\
ATGGCGTACGCGCCGTCCAACCTGGGCCGCAGCCGCGCGGGCGGCGGCGGCGACGCCAAGACGTTCAAGC\
GCAGCTACCCGCACGCCAAGCCGCCCTACTCGTACATCTCGCTCATCACCATGGCCATCCAGCAGGCGCC\
CAGCAAGATGCTCACGCTGAGCGAGATCTACCAGTGGATCATGGACCTCTTCCCCTATTACCGGCAGAAC\
CAGCAGCGCTGGCAGAACTCCATCCGCCACTCGCTGTCCTTCAATGACTGCTTCGTCAAGGTGGCACGCT\
CCCCGGACAAGCCGGGCAAGGGCTCCTACTGGACGCTGCACCCGGACTCCGGCAACATGTTCGAGAACGG\
CTGCTACTTGCGCCGCCAGAAGCGCTTCAAGTGCGAGAAGCAGCCGGGGGCCGGCGGCGGGGGCGGGAGC\
GGAAGCGGGGGCAGCGGCGCCAAGGGCGGCCCTGAGAGCCGCAAGGACCCCTCTGGCGCCTCTAACCCCA\
GCGCCGACTCGCCCCTCCATCGGGGTGTGCACGGGAAGACCGGCCAGCTAGAGGGCGCGCCGGCCCCCGG\
GCCCGCCGCCAGCCCCCAGACTCTGGACCACAGTGGGGCGACGGCGACAGGGGGCGCCTCGGAGTTGAAG\
ACTCCAGCCTCCTCAACTGCGCCCCCCATAAGCTCCGGGCCCGGGGCGCTGGCCTCTGTGCCCGCCTCTC\
ACCCGGCACACGGCTTGGCACCCCACGAGTCCCAGCTGCACCTGAAAGGGGACCCCCACTACTCCTTCAA\
CCACCCGTTCTCCATCAACAACCTCATGTCCTCCTCGGAGCAGCAGCATAAGCTGGACTTCAAGGCATAC\
GAACAGGCACTGCAATACTCGCCTTACGGCTCTACGTTGCCCGCCAGCCTGCCTCTAGGCAGCGCCTCGG\
TGACCACCAGGAGCCCCATCGAGCCCTCAGCCCTGGAGCCGGCGTACTACCAAGGTGTGTATTCCAGACC\
CGTCCTAAACACTTCCTAG")
        self.addGene("FOXA2", "Transcription factor that plays an important role in regulating \
the expression of genes involved in cell growth, proliferation, differentiation, \
and longevity", gene.GenePurpose.CONTROLLER, "ATGCTGGGAGCGGTGAAGATGGAAGGGCACGAGCCGTCCGACTGGAGCAGCTACTATGCAGAGCCCGAGG\
GCTACTCCTCCGTGAGCAACATGAACGCCGGCCTGGGGATGAACGGCATGAACACGTACATGAGCATGTC\
GGCGGCCGCCATGGGCAGCGGCTCGGGCAACATGAGCGCGGGCTCCATGAACATGTCGTCGTACGTGGGC\
GCTGGCATGAGCCCGTCCCTGGCGGGGATGTCCCCCGGCGCGGGCGCCATGGCGGGCATGGGCGGCTCGG\
CCGGGGCGGCCGGCGTGGCGGGCATGGGGCCGCACTTGAGTCCCAGCCTGAGCCCGCTCGGGGGGCAGGC\
GGCCGGGGCCATGGGCGGCCTGGCCCCCTACGCCAACATGAACTCCATGAGCCCCATGTACGGGCAGGCG\
GGCCTGAGCCGCGCCCGCGACCCCAAGACCTACAGGCGCAGCTACACGCACGCAAAGCCGCCCTACTCGT\
ACATCTCGCTCATCACCATGGCCATCCAGCAGAGCCCCAACAAGATGCTGACGCTGAGCGAGATCTACCA\
GTGGATCATGGACCTCTTCCCCTTCTACCGGCAGAACCAGCAGCGCTGGCAGAACTCCATCCGCCACTCG\
CTCTCCTTCAACGACTGTTTCCTGAAGGTGCCCCGCTCGCCCGACAAGCCCGGCAAGGGCTCCTTCTGGA\
CCCTGCACCCTGACTCGGGCAACATGTTCGAGAACGGCTGCTACCTGCGCCGCCAGAAGCGCTTCAAGTG\
CGAGAAGCAGCTGGCGCTGAAGGAGGCCGCAGGCGCCGCCGGCAGCGGCAAGAAGGCGGCCGCCGGAGCC\
CAGGCCTCACAGGCTCAACTCGGGGAGGCCGCCGGGCCGGCCTCCGAGACTCCGGCGGGCACCGAGTCGC\
CTCACTCGAGCGCCTCCCCGTGCCAGGAGCACAAGCGAGGGGGCCTGGGAGAGCTGAAGGGGACGCCGGC\
TGCGGCGCTGAGCCCCCCAGAGCCGGCGCCCTCTCCCGGGCAGCAGCAGCAGGCCGCGGCCCACCTGCTG\
GGCCCGCCCCACCACCCGGGCCTGCCGCCTGAGGCCCACCTGAAGCCGGAACACCACTACGCCTTCAACC\
ACCCGTTCTCCATCAACAACCTCATGTCCTCGGAGCAGCAGCACCACCACAGCCACCACCACCACCAACC\
CCACAAAATGGACCTCAAGGCCTACGAACAGGTGATGCACTACCCCGGCTACGGTTCCCCCATGCCTGGC\
AGCTTGGCCATGGGCCCGGTCACGAACAAAACGGGCCTGGACGCCTCGCCCCTGGCCGCAGATACCTCCT\
ACTACCAGGGGGTGTACTCCCGGCCCATTATGAACTCCTCTTAA\
")
        self.addGene("FOXA3", "Transcription factor that plays an important role in regulating\
the expression of genes involved in cell growth, proliferation, differentiation, \
and longevity", gene.GenePurpose.CONTROLLER, "ATGCTGGGCTCAGTGAAGATGGAGGCCCATGACCTGGCCGAGTGGAGCTACTACCCGGAGGCGGGCGAGG\
TCTACTCGCCGGTGACCCCAGTGCCCACCATGGCCCCCCTCAACTCCTACATGACCCTGAATCCTCTAAG\
CTCTCCCTATCCCCCTGGGGGGCTCCCTGCCTCCCCACTGCCCTCAGGACCCCTGGCACCCCCAGCACCT\
GCAGCCCCCCTGGGGCCCACTTTCCCAGGCCTGGGTGTCAGCGGTGGCAGCAGCAGCTCCGGGTACGGGG\
CCCCGGGTCCTGGGCTGGTGCACGGGAAGGAGATGCCGAAGGGGTATCGGCGGCCCCTGGCACACGCCAA\
GCCACCGTATTCCTATATCTCACTCATCACCATGGCCATCCAGCAGGCGCCGGGCAAGATGCTGACCTTG\
AGTGAAATCTACCAGTGGATCATGGACCTCTTCCCTTACTACCGGGAGAATCAGCAGCGCTGGCAGAACT\
CCATTCGCCACTCGCTGTCTTTCAACGACTGCTTCGTCAAGGTGGCGCGTTCCCCAGACAAGCCTGGCAA\
GGGCTCCTACTGGGCCCTACACCCCAGCTCAGGGAACATGTTTGAGAATGGCTGCTACCTGCGCCGCCAG\
AAACGCTTCAAGCTGGAGGAGAAGGTGAAAAAAGGGGGCAGCGGGGCTGCCACCACCACCAGGAACGGGA\
CAGGGTCTGCTGCCTCGACCACCACCCCCGCGGCCACAGTCACCTCCCCGCCCCAGCCCCCGCCTCCAGC\
CCCTGAGCCTGAGGCCCAGGGCGGGGAAGATGTGGGGGCTCTGGACTGTGGCTCACCCGCTTCCTCCACA\
CCCTATTTCACTGGCCTGGAGCTCCCAGGGGAGCTGAAGCTGGACGCGCCCTACAACTTCAACCACCCTT\
TCTCCATCAACAACCTAATGTCAGAACAGACACCAGCACCTCCCAAACTGGACGTGGGGTTTGGGGGCTA\
CGGGGCTGAAGGTGGGGAGCCTGGAGTCTACTACCAGGGCCTCTATTCCCGCTCTTTGCTTAATGCATCC\
TAG")
        self.addGene("FOXB1", "Transcription factor that plays an important role in regulating\
the expression of genes involved in cell growth, proliferation, differentiation, \
and longevity", gene.GenePurpose.CONTROLLER, "ATGCCTCGGCCCGGCCGCAACACGTACAGCGACCAGAAGCCGCCCTACTCGTACATCTCGCTGACCGCTA\
TGGCCATCCAGAGCTCTCCCGAGAAGATGCTGCCGCTGAGCGAGATCTACAAGTTCATCATGGACCGCTT\
CCCCTACTACAGGGAGAACACGCAGCGCTGGCAGAACAGTCTGCGCCACAACCTCTCCTTCAACGACTGC\
TTCATCAAGATCCCGCGGCGGCCGGACCAGCCAGGCAAGGGCAGCTTCTGGGCGCTGCACCCAAGCTGCG\
GGGACATGTTCGAGAACGGCAGCTTCCTGCGGCGCCGCAAGCGCTTCAAGGTGCTTAAGTCCGACCACCT\
GGCGCCCAGCAAGCCAGCCGACGCGGCGCAGTACCTGCAGCAGCAGGCCAAGCTGCGGCTCAGCGCGCTG\
GCGGCCTCGGGCACGCACCTGCCACAGATGCCCGCCGCCGCCTACAACTTGGGCGGCGTGGCGCAGCCCT\
CGGGCTTCAAGCACCCCTTCGCCATCGAGAACATCATCGCGCGGGAATACAAGATGCCTGGGGGGCTGGC\
CTTCTCCGCCATGCAGCCGGTGCCCGCTGCCTACCCGCTCCCCAACCAGTTGACTACCATGGGCAGCTCG\
CTGGGCACCGGCTGGCCACACGTGTATGGCTCCGCCGGCATGATCGACTCGGCCACCCCCATCTCCATGG\
CGAGTGGCGACTACAGCGCCTACGGCGTGCCGTTGAAGCCGCTGTGCCACGCGGCGGGCCAAACGCTGCC\
CGCCATCCCCGTGCCCATTAAGCCCACGCCGGCCGCCGTGCCCGCGCTGCCTGCGCTGCCAGCGCCCATC\
CCCACCTTGCTCTCGAACTCGCCGCCCTCGCTCAGCCCCACGTCCTCGCAAACAGCCACCAGCCAAAGCA\
GCCCCGCCACCCCCAGCGAAACGCTCACCAGCCCGGCCTCCGCCTTGCACTCGGTGGCGGTGCACTGA")
        self.addGene("FOXB2", "Transcription factor that plays an important role in regulating\
the expression of genes involved in cell growth, proliferation, differentiation, \
and longevity", gene.GenePurpose.CONTROLLER, "ATGCCGCGGCCGGGGAAGAGCTCGTACAGCGACCAAAAACCGCCCTACTCTTACATCTCGCTGACCGCCA\
TGGCAATCCAGCACTCGGCCGAGAAGATGCTGCCGCTGAGCGACATCTACAAGTTCATCATGGAGCGCTT\
CCCCTACTACCGCGAGCACACACAGCGCTGGCAGAACAGCCTGCGCCACAACCTCTCCTTCAACGACTGC\
TTCATCAAGATTCCGCGGAGGCCCGACCAGCCTGGCAAGGGTAGCTTCTGGGCGCTGCACCCCGACTGCG\
GGGACATGTTCGAGAACGGCAGCTTCCTGCGGCGTCGCAAGCGCTTCAAGGTGCTGCGCGCCGACCATAC\
TCACTTGCACGCGGGAAGCACCAAGAGCGCGCCGGGCGCCGGTCCGGGAGGGCACCTTCACCCCCATCAC\
CACCACCACCCCCACCACCACCATCATCACCACGCTGCCGCACACCACCACCATCACCACCACCCACCCC\
AGCCGCCGCCGCCGCCGCCCCCGCCGCCGCCGCACATGGTACACTATTTCCATCAGCAACCGCCTACTGC\
TCCGCAGCCGCCTCCGCACCTCCCGTCACAGCCCCCGCAGCAACCGCCCCAGCAGTCGCAGCCTCAGCAG\
CCGTCTCACCCCGGCAAGATGCAGGAGGCGGCGGCCGTGGCGGCGGCGGCGGCGGCGGCCGCGGCAGCCG\
CGGTGGGCAGCGTGGGACGCCTGTCTCAGTTCCCACCCTACGGGCTGGGCTCGGCCGCCGCCGCTGCCGC\
CGCGGCCGCGGCGTCCACGTCAGGCTTCAAGCACCCCTTTGCCATTGAGAACATTATTGGCCGGGACTAC\
AAGGGCGTGCTGCAGGCTGGAGGGCTGCCCTTGGCGTCCGTCATGCACCACCTGGGCTACCCCGTGCCCG\
GCCAGCTTGGCAACGTCGTCAGCTCCGTGTGGCCGCACGTTGGCGTCATGGATTCGGTGGCCGCCGCCGC\
GGCCGCCGCAGCCGCAGCCGGAGTCCCTGTAGGCCCGGAGTATGGGGCCTTCGGGGTCCCGGTCAAGTCC\
CTGTGCCACTCGGCAAGCCAGAGCCTGCCTGCCATGCCGGTGCCCATCAAGCCCACGCCTGCGCTGCCGC\
CCGTGTCCGCGCTGCAGCCGGGGCTCACTGTCCCCGCGGCTTCGCAGCAGCCTCCGGCGCCATCCACCGT\
GTGCTCCGCGGCCGCGGCCTCGCCCGTTGCCTCTCTGCTGGAGCCCACAGCCCCTACCTCGGCCGAAAGC\
AAGGGCGGCTCCTTGCACTCGGTGCTAGTGCACTCCTAG")
        self.addGene("FOXC1", "Transcription factor that plays an important role in regulating\
the expression of genes involved in cell growth, proliferation, differentiation, \
and longevity", gene.GenePurpose.CONTROLLER, "ATGCAGGCGCGCTACTCCGTGTCCAGCCCCAACTCCCTGGGAGTGGTGCCCTACCTCGGCGGCGAGCAGA\
GCTACTACCGCGCGGCGGCCGCGGCGGCCGGGGGCGGCTACACCGCCATGCCGGCCCCCATGAGCGTGTA\
CTCGCACCCTGCGCACGCCGAGCAGTACCCGGGCGGCATGGCCCGCGCCTACGGGCCCTACACGCCGCAG\
CCGCAGCCCAAGGACATGGTGAAGCCGCCCTATAGCTACATCGCGCTCATCACCATGGCCATCCAGAACG\
CCCCGGACAAGAAGATCACCCTGAACGGCATCTACCAGTTCATCATGGACCGCTTCCCCTTCTACCGGGA\
CAACAAGCAGGGCTGGCAGAACAGCATCCGCCACAACCTCTCGCTCAACGAGTGCTTCGTCAAGGTGCCG\
CGCGACGACAAGAAGCCGGGCAAGGGCAGCTACTGGACGCTGGACCCGGACTCCTACAACATGTTCGAGA\
ACGGCAGCTTCCTGCGGCGGCGGCGGCGCTTCAAGAAGAAGGACGCGGTGAAGGACAAGGAGGAGAAGGA\
CAGGCTGCACCTCAAGGAGCCGCCCCCGCCCGGCCGCCAGCCCCCGCCCGCGCCGCCGGAGCAGGCCGAC\
GGCAACGCGCCCGGTCCGCAGCCGCCGCCCGTGCGCATCCAGGACATCAAGACCGAGAACGGTACGTGCC\
CCTCGCCGCCCCAGCCCCTGTCCCCGGCCGCCGCCCTGGGCAGCGGCAGCGCCGCCGCGGTGCCCAAGAT\
CGAGAGCCCCGACAGCAGCAGCAGCAGCCTGTCCAGCGGGAGCAGCCCCCCGGGCAGCCTGCCGTCGGCG\
CGGCCGCTCAGCCTGGACGGTGCGGATTCCGCGCCGCCGCCGCCCGCGCCCTCCGCCCCGCCGCCGCACC\
ATAGCCAGGGCTTCAGCGTGGACAACATCATGACGTCGCTGCGGGGGTCGCCGCAGAGCGCGGCCGCGGA\
GCTCAGCTCCGGCCTTCTGGCCTCGGCGGCCGCGTCCTCGCGCGCGGGGATCGCACCCCCGCTGGCGCTC\
GGCGCCTACTCGCCCGGCCAGAGCTCCCTCTACAGCTCCCCCTGCAGCCAGACCTCCAGCGCGGGCAGCT\
CGGGCGGCGGCGGCGGCGGCGCGGGGGCCGCGGGGGGCGCGGGCGGCGCCGGGACCTACCACTGCAACCT\
GCAAGCCATGAGCCTGTACGCGGCCGGCGAGCGCGGGGGCCACTTGCAGGGCGCGCCCGGGGGCGCGGGC\
GGCTCGGCCGTGGACGACCCCCTGCCCGACTACTCTCTGCCTCCGGTCACCAGCAGCAGCTCGTCGTCCC\
TGAGTCACGGCGGCGGCGGCGGCGGCGGCGGGGGAGGCCAGGAGGCCGGCCACCACCCTGCGGCCCACCA\
AGGCCGCCTCACCTCGTGGTACCTGAACCAGGCGGGCGGAGACCTGGGCCACTTGGCGAGCGCGGCGGCG\
GCGGCGGCGGCCGCAGGCTACCCGGGCCAGCAGCAGAACTTCCACTCGGTGCGGGAGATGTTCGAGTCAC\
AGAGGATCGGCTTGAACAACTCTCCAGTGAACGGGAATAGTAGCTGTCAAATGGCCTTCCCTTCCAGCCA\
GTCTCTGTACCGCACGTCCGGAGCTTTCGTCTACGACTGTAGCAAGTTTTGA")
        self.addGene("FOXD4", "Transcription factor that plays an important role in regulating\
the expression of genes involved in cell growth, proliferation, differentiation, \
and longevity", gene.GenePurpose.CONTROLLER, "ATGAACTTGCCAAGAGCTGAGCGCCTTCGCTCCACACCGCAGCGCAGCCTCCGGGACTCCGATGGGGAAG\
ACGGTAAAATCGATGTCCTGGGAGAGGAGGAAGATGAAGACGAGGAGGAGGCGGCGAGCCAGCAGTTCCT\
AGAGCAGTCGCTCCAGCCGGGGCTGCAGGTGGCCCGGTGGGGCGGGGTTGCGCTTCCCCGAGAGCACATC\
GAGGGCGGCGGCGGCCCGAGCGACCCCTCAGAGTTTGGCACCGAGTTCAGGGCACCGCCAAGGTCTGCGG\
CGGCCTCTGAAGATGCCCGGCAGCCGGCAAAGCCCCCCTCCTCGTACATCGCGCTCATCACCATGGCCAT\
CCTGCAAAGCCCGCACAAGCGCCTCACGCTCAGCGGCATCTGCGCCTTCATTAGTGACCGCTTCCCCTAC\
TACCGCCGCAAGTTCCCCGCCTGGCAGAACAGCATCCGCCACAACCTCTCGCTGAACGACTGCTTCGTCA\
AGATCCCCCGCGAGCCGGGCCGCCCAGGCAAGGGCAACTACTGGAGCCTGGACCCCGCCTCCCAGGACAT\
GTTCGACAATGGCAGCTTTCTCCGGCGTAGGAAGCGTTTCCAGCGCCACCAACCGACCCCGGGAGCCCAC\
CTGCCCCACCCCTTCCCTCTACCTGCTGCACACGCCGCCCTGCACAACCCCCGCCCAGGCCCTCTGCTTG\
GGGCCCCTGCCCCGCCGCAGCCAGTCCCGGGGGCCTACCCCAACACCGGCCCCGGGAGACGCCCTTACGC\
TCTGCTGCACCCGCATCCTCCTCGCTACCTACTGCTCTCGGCCCCCGCCTATGCCGGGGCACCGAAGAAA\
GCAGAAGGCGCGGACCTGGCGACCCCGGCACCCTTCCCGTGCTGCAGCCCTCACTTGGTCCTCAGCCTTG\
GGAGGAGGGCAAGGGTCTGGCGTCGCCACCGGGAGGCGGATGCATCTCTTTCAGCATTGAGAGTATCATG\
CAAGGGGTCAGGGGAGCGGGTACAGGGGCTGCGCAGAGTTTGTCCCCGACCGCGTGGAGCTACTGCCCCC\
TGCTCCAGCGACCGTCAAGCCTGTCGGACAATTTTGCAGCAACAGCAGCGGCATCAGGAGGAGGACTGCG\
CCAACGGCTGCGCTCCCACCAAGGGCGCGGTGCTGGGCGGGCACCTGTCGGCCGCGTCGGCGCTGCTGCG\
GTATCAGGCGGTGGCAGAGGGCTCTGGGCTGACATCGCTGGCCGCCCCTTTGGGCGGAGAGGGGACCTCA\
CCAGTTTTTTTAGTATCGCCCACGCCCAGTTCCCTGGCCGAGTCCGCAGGGCCCTCCTAG")
        self.addGene("FOXD2", "Transcription factor that plays an important role in regulating\
the expression of genes involved in cell growth, proliferation, differentiation, \
and longevity", gene.GenePurpose.CONTROLLER, "ATGACCCTGGGCAGCTGCTGCTGCGAGATCATGTCCTCCGAGAGCTCCCCGGCCGCGCTGTCCGAGGCCG\
ACGCAGACATAGACGTGGTGGGCGGCGGCAGCGGCGGGGGGGAGCTCCCAGCTCGCTCCGGGCCCCGCGC\
CCCCCGGGACGTGCTCCCCCACGGCCACGAGCCTCCCGCGGAGGAAGCCGAGGCAGACTTAGCCGAGGAC\
GAGGAGGAGTCTGGTGGCTGCTCGGACGGCGAGCCCCGCGCTCTGGCGTCCCGGGGGGCGGCGGCCGCAG\
CGGGGAGCCCGGGGCCAGGCGCCGCGGCGGCCCGCGGCGCAGCGGGGCCCGGGCCGGGACCGCCGTCGGG\
GGGCGCGGCGACGCGGAGCCCGCTGGTGAAGCCGCCCTACTCGTACATCGCGCTCATCACCATGGCCATC\
CTGCAGAGCCCCAAGAAGCGGCTGACGTTGAGCGAGATCTGCGAGTTCATCAGCGGCCGCTTCCCCTACT\
ACCGGGAGAAGTTCCCCGCCTGGCAGAACAGCATCCGCCACAACCTCTCTCTCAACGACTGCTTCGTCAA\
GATCCCCCGCGAGCCGGGCAACCCGGGCAAGGGCAACTACTGGACGCTGGACCCGGAGTCGGCCGACATG\
TTCGACAACGGCAGCTTCCTGCGGCGTCGCAAGCGCTTCAAGCGGCAGCCCCTGCCGCCGCCGCACCCAC\
ACCCGCACCCTCACCCGGAGCTGCTGCTGCGTGGCGGGGCCGCGGCGGCGGGGGATCCCGGCGCTTTCCT\
GCCCGGCTTCGCTGCCTACGGCGCCTACGGCTACGGCTACGGGCTGGCTCTCCCGGCCTACGGCGCACCC\
CCGCCGGGGCCGGCCCCGCATCCGCACCCGCACCCGCACGCCTTCGCTTTCGCCGCGGCAGCCGCCGCCG\
CTCCTTGCCAGCTGTCGGTACCCCCAGGCCGCGCCGCCGCGCCTCCACCCGGACCTCCGACGGCCTCGGT\
GTTCGCAGGCGCGGGATCGGCCCCAGCTCCTGCGCCTGCCTCAGGCTCGGGCCCGGGCCCGGGCCCCGCA\
GGCCTGCCCGCCTTCCTGGGCGCGGAGCTGGGCTGCGCCAAAGCCTTCTACGCGGCGTCCCTGAGTCCTC\
CCGCAGCCGGCACCGCGGCGGGTCTGCCCACCGCACTTCTGCGCCAGGGCCTCAAGACGGACGCGGGCGG\
TGGTGCAGGCGGCGGGGGCGCCGGGGCAGGGCAGAGGCCTTCCTTCTCTATAGACCACATCATGGGCCAC\
GGTGGCGGCGGGGCAGCACCCCCGGGCGCCGGCGAGGGCTCTCCGGGACCGCCATTCGCGGCAGCCGCGG\
GTCCTGGGGGCCAAGCCCAGGTCTTGGCCATGCTGACTGCTCCGGCCCTGGCTCCCGTTGCTGGCCACAT\
TCGCCTCTCGCATCCCGGGGACGCGCTGCTGTCCTCAGGGTCCCGGTTTGCCAGCAAAGTCGCCGGCCTT\
AGTGGCTGCCACTTCTGA")
