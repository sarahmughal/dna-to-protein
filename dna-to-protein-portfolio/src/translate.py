from typing import Dict

CODON_TABLE: Dict[str, str] = {
    # U
    'UUU':'F','UUC':'F','UUA':'L','UUG':'L',
    'UCU':'S','UCC':'S','UCA':'S','UCG':'S',
    'UAU':'Y','UAC':'Y','UAA':'*','UAG':'*',
    'UGU':'C','UGC':'C','UGA':'*','UGG':'W',
    # C
    'CUU':'L','CUC':'L','CUA':'L','CUG':'L',
    'CCU':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAU':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'CGU':'R','CGC':'R','CGA':'R','CGG':'R',
    # A
    'AUU':'I','AUC':'I','AUA':'I','AUG':'M',
    'ACU':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAU':'N','AAC':'N','AAA':'K','AAG':'K',
    'AGU':'S','AGC':'S','AGA':'R','AGG':'R',
    # G
    'GUU':'V','GUC':'V','GUA':'V','GUG':'V',
    'GCU':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAU':'D','GAC':'D','GAA':'E','GAG':'E',
    'GGU':'G','GGC':'G','GGA':'G','GGG':'G',
}

def dna_to_rna(dna: str) -> str:
    return dna.upper().replace('T','U')

def translate_rna(rna: str, frame: int = 0, stop_behavior: str = "truncate") -> str:
    """Translate RNA to amino acids.
    stop_behavior: 'truncate' (stop at first *), 'keep' (include *), or 'ignore' (skip *).
    """
    rna = rna.upper()
    if frame not in (0,1,2):
        raise ValueError("frame must be 0,1,2")
    aa = []
    for i in range(frame, len(rna)-2, 3):
        codon = rna[i:i+3]
        amino = CODON_TABLE.get(codon, 'X')  # unknown -> X
        if amino == '*':
            if stop_behavior == 'truncate':
                break
            elif stop_behavior == 'ignore':
                continue
        aa.append(amino)
    return ''.join(aa)

def translate_dna(dna: str, frame: int = 0, stop_behavior: str = 'truncate') -> str:
    return translate_rna(dna_to_rna(dna), frame=frame, stop_behavior=stop_behavior)
