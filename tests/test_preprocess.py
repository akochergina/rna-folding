import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from modules.preprocess_rna import *

def test_preprocess_rna_seq():
    """
    Test function for preprocess_rna_seq

    Returns:
        String All tests passed if all tests are successful, otherwise raises an exception
    """
    assert preprocess_rna_seq('ATGC') == 'AUGC', "Error in test case 1: T to U"
    assert preprocess_rna_seq('atgc') == 'AUGC', "Error in test case 2: lower case and T to U"
    assert preprocess_rna_seq('ATGCATGC') == 'AUGCAUGC', "Error in test case 3: T to U"
    assert preprocess_rna_seq('ATGCTTGC') == 'AUGCUUGC', "Error in test case 4: T to U"
    assert preprocess_rna_seq('ATGCTTGCATGC') == 'AUGCUUGCAUGC', "Error in test case 5: T to U"
    assert preprocess_rna_seq('ATGCTTGCatgc') == 'AUGCUUGCAUGC', "Error in test case 6: lower case and T to U"
    return "Preprocess RNA sequence tests passed"

def test_sequence_identity():
    """
    Test function for sequence_identity

    Returns:
        String All tests passed if all tests are successful, otherwise raises an exception
    """
    assert sequence_identity("ACGT", "ACGT", False) == "100.00%", f"No preprocessing. Identical sequences test failed, percentage should be 100 and not {sequence_identity('ACGT', 'ACGT', False)}"
    assert sequence_identity("ACGT", "AC-A", False) == "50.00%", f"No preprocessing. Half identical sequences test failed, percentage should be 50 and not {sequence_identity('ACGT', 'ACGT', False)}"
    assert sequence_identity("ACGT", "----", False) == "0.00%", f"No preprocessing. No identical sequences test failed, percentage should be 0 and not {sequence_identity('ACGT', 'ACGT', False)}"
    assert sequence_identity("ACGT", "ACGU", True) == "100.00%", f"Preprocessing T to U. Identical sequences test failed, percentage should be 100 and not {sequence_identity('ACGT', 'ACGU', True)}"
    assert sequence_identity("ACGT", "AC-U", True) == "75.00%", f"Preprocessing T to U. 3/4 identical sequences test failed, percentage should be 75 and not {sequence_identity('ACGT', 'AC-U', True)}"
    assert sequence_identity("ACGT", "----", True) == "0.00%", f"Preprocessing T to U. No identical sequences test failed, percentage should be 0 and not {sequence_identity('ACGT', '----', True)}"
    assert sequence_identity("acgu", "ACGU", True) == "100.00%", f"Preprocessing lower case. Identical sequences test failed, percentage should be 100 and not {sequence_identity('acgu', 'ACGU', True)}"

    return "Sequence identity percentage tests passed"

def all_preprocess_tests():
    print(test_preprocess_rna_seq())
    print(test_sequence_identity())
    print("Example of using the preprocess_rna module:")
    print(f"Preprocessed RNA sequence aTGc: {preprocess_rna_seq('ATGC')}")
    print(f"Sequence identity of ATGC and AU-C: {sequence_identity('ATGC', 'AU-C')}")
