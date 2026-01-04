class VariantTranscript:
    def __init__(self, chromosome, position, ref, alt, cdna_transcript=None):
        self.chromosome = chromosome
        self.position = position
        self.ref = ref
        self.alt = alt

        self.exon = None
        self.gene = None
        self.c_dot = None
        self.p_dot = None

        self.cdna_transcript = None
        self.protein_transcript = None

