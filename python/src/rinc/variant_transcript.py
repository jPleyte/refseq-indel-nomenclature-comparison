class VariantTranscript:
    def __init__(self, chromosome, position, reference, alt, cdna_transcript=None):
        self.chromosome = chromosome
        self.position = position
        self.reference = reference
        self.alt = alt

        self.exon = None
        self.gene = None
        self.c_dot = None
        self.p_dot1 = None
        self.p_dot2 = None

        self.cdna_transcript = cdna_transcript
        self.protein_transcript = None
        
        self.strand = None
        self.notes = []

