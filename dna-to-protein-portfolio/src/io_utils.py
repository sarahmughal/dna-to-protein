from typing import Dict, List, Tuple, Iterable

def read_fasta(path: str) -> Dict[str, str]:
    """Read a FASTA file into a dict {header: sequence}. Supports multiline sequences."""
    records: Dict[str, List[str]] = {}
    current = None
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].strip()
                if current in records:
                    raise ValueError(f"Duplicate FASTA header: {current}")
                records[current] = []
            else:
                if current is None:
                    raise ValueError("FASTA missing header before sequence lines")
                records[current].append(line.upper())
    return {h: "".join(seq) for h, seq in records.items()}

def write_fasta(path: str, records: Dict[str, str], width: int = 80) -> None:
    """Write records {header: sequence} to FASTA with line wrapping."""
    with open(path, "w", encoding="utf-8") as f:
        for header, seq in records.items():
            f.write(f">{header}\n")
            for i in range(0, len(seq), width):
                f.write(seq[i:i+width] + "\n")

def iter_fasta(path: str) -> Iterable[Tuple[str, str]]:
    """Stream FASTA records as (header, sequence)."""
    with open(path, "r", encoding="utf-8") as f:
        header = None
        chunks = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.upper())
        if header is not None:
            yield header, "".join(chunks)
