from dataclasses import dataclass, field

@dataclass(slots=True)
class VariantTranscript:
    chromosome: str
    position: int
    reference: str
    alt: str
    cdna_transcript: str
    
    g_dot: str = None
    exon: int = None
    gene: str = None
    c_dot: str = None
    p_dot1: str = None
    p_dot3: str = None
    protein_transcript: str = None
    strand: int = None
    
    # exon, intron, utr, or intergenic
    genomic_region_type: str = None
    
    # Synonymous, missense, start loss, stop gain, etc
    protein_variant_type: str = None
    
    notes: list[str] = field(default_factory=list)
    additional_fields: dict = field(default_factory=dict)

    def __str__(self):
        return f"{self.chromosome}-{self.position}-{self.reference}-{self.alt} {self.cdna_transcript}:{self.c_dot}"
