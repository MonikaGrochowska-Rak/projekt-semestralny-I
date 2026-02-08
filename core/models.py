from dataclasses import dataclass

@dataclass
class SequenceRecord:
    id: str
    description: str
    sequence: str