class VariantTranscript:
    def __init__(self, chromosome, position, reference, alt, cdna_transcript=None):
        self.chromosome = chromosome
        self.position = position
        self.reference = reference
        self.alt = alt

        self.g_dot = None
        
        self.exon = None
        self.gene = None
        self.c_dot = None
        self.p_dot1 = None
        self.p_dot3 = None
        
        self.cdna_transcript = cdna_transcript
        self.protein_transcript = None

        self.strand = None
        self.notes = []
        
        self.genomic_region_type = None
        self.protein_variant_type = None
        

    def __str__(self):
        return f"{self.chromosome}-{self.position}-{self.reference}-{self.alt} {self.cdna_transcript}:{self.c_dot}"