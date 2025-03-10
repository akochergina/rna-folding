import warnings
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from modules.rna_structure import *

class RNAMolecule:
    """
    Represents an RNA sequence and its corresponding structure.
    
    Attributes:
        seq (str): aligned RNA sequence (e.g., "AUGC").
        structure (str, optional): RNA structure in dot-bracket notation (e.g., "..((..))..").
        
    Methods:
    """
    @staticmethod
    def preprocess_rna_seq(seq: str) -> str:
        """
        Preprocess RNA sequence by uppercasing and converting T to U

        Args:
            seq: RNA sequence

        Returns:
            preprocessed RNA sequence    
        """
        return seq.upper().replace('T', 'U')

    def __init__(self, seq, structure=None):
        """
        Initializes an RNAMolecule object.
        
        Args:
            seq (str): RNA sequence.
            structure (str, optional): RNA structure in dot-bracket notation. Defaults to None.
        """
        self.seq = self.preprocess_rna_seq(seq)
        if structure is not None:
            self.structure = structure
        else:
            self.structure = None

    def __str__(self):
        """
        Returns a string representation of the RNA molecule.
        
        Returns:
            str: Sequence and structure in the format "Sequence: <seq>, Structure: <structure>".
        """
        return f"Sequence: {self.seq}\nStructure: {self.structure or 'None'}"
    
    @staticmethod
    def edges_to_dot_bracket(edges, seq_len):
        """
        Convert edges to dot bracket notation

        Args:
            edges: List of edges in the alignment graph
            seq_len: Length of the sequence
            
        Returns:
            Dot bracket notation of the alignment graph

        Example:
            edges_to_dot_bracket([(2, 8), (3, 7), (9, 13)], 13) -> '.((...))(...)'
        """
        dot_bracket = ['.'] * seq_len
        for i, j in edges:
            assert i!=j, f"Edges must be of the form (i, i) where i != j, i={i}, j={j}"
            assert i <= seq_len and j <= seq_len, f"Edges must be within the sequence length, i={i}, j={j}, seq_len={seq_len}"
            assert i < j, f"Edges must be of the form (i, j) where i < j, i={i}, j={j}"
            dot_bracket[i-1] = '('
            dot_bracket[j-1] = ')'
        return ''.join(dot_bracket)

    def parse_rna_structure(self) -> list:
        """
        Parse RNA structure from dot bracket notation

        Args:
            dot_bracket: Dot bracket notation of the RNA structure
            
        Returns:
            List of base pairs

        Example:
            parse_rna_structure('.((...))(...)') -> [(2, 8), (3, 7), (9, 13)]
        """
        dot_bracket = self.structure
        stack = []
        base_pairs = []
        for i, c in enumerate(dot_bracket):
            if c == '(':
                stack.append(i+1)
            elif c == ')':
                assert len(stack) > 0, "Unmatched brackets )"
                base_pairs.append((stack.pop(), i+1))
        base_pairs = sorted(base_pairs, key=lambda x: x[0])
        assert len(stack) == 0, "Unmatched brackets ("
        return base_pairs
    
    @staticmethod
    def is_canonical_base_pair(x, y):
        """
        Check if two bases form a canonical base pair: (A,U),(C,G),(G,C),(G,U),(U,A), or (U,G)

        Args:
            x: base 1
            y: base 2
            
        Returns:
            True if x, y form a canonical base pair, False otherwise
        """
        assert x in ['A', 'C', 'G', 'U'], f"Invalid base {x}"
        assert y in ['A', 'C', 'G', 'U'], f"Invalid base {y}"
        return (x == 'A' and y == 'U') or (x == 'U' and y == 'A') or (x == 'G' and y == 'C') or (x == 'C' and y == 'G') or (x == 'G' and y == 'U') or (x == 'U' and y == 'G')
    
    def is_canonical(self) -> bool:
        """
        Check if the RNA structure is canonical

        Args:
            dot_bracket: Dot bracket notation of the RNA structure
            seq: RNA sequence
            
        Returns:
            True if the RNA structure is canonical, False otherwise
        """
        seq = self.seq
        dot_bracket = self.structure
        if dot_bracket is None:
            warnings.warn("Structure is not provided")
            return False
        base_pairs = self.parse_rna_structure()
        for i, j in base_pairs:
            assert i < j, f"Base pair ({i}, {j}) is invalid, i must be less than j"
            assert i > 0 and j > 0, f"Base pair ({i}, {j}) is invalid, positions must be positive"
            assert i <= len(seq) and j <= len(seq), f"Base pair ({i}, {j}) is invalid, sequence length is {len(seq)}"
            if not is_canonical_base_pair(seq[i-1], seq[j-1]):
                return False
                    
        return True
    
    def rna_structure_projection(self):
        """ 
        Perform RNA structure projection. The projected structure does not contain any base pair that involved gaps in the alignment string.

        Args:
            aligned_sequence: RNA sequence
            dot_bracket_structure: RNA structure in dot-bracket notation
        
        Returns:
            (sequence, projected_structure): Tuple of RNA sequence and projected structure
        
        Example:
            rna_structure_projection('GCC-CUUAG-U-GAAUCCAGC', '((.((...))(((...)))))') == ('GCCCUUAGUGAAUCCAGC', '((.(...)((...).)))')
        """
        assert len(self.seq) == len(self.structure), "Sequences must have the same length"
        edges = parse_rna_structure(self.structure)
        sequence = []
        projected_structure = []
        to_convert = [] # if first bracket is deleted, we need to delete the second bracket as well, so to convert second one in a dot
        for i, (s, b) in enumerate(zip(self.seq, self.structure)):
            if i in to_convert:
                if b == ')':
                    b = '.'
                to_convert.remove(i)
            if s == '-':
                if b == '(':
                    # we need to find in edges if there is a pair with i
                    found = False
                    for edge in edges:
                        if edge[0] == i+1:
                            found = True
                            break
                    if found:
                        to_convert.append(edge[1]-1)
            else:
                sequence.append(s)
                projected_structure.append(b)
        return ''.join(sequence), ''.join(projected_structure)