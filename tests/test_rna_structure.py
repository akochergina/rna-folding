import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from modules.rna_structure import *

def test_alignment_edges():
    """
    Test function for alignment_edges

    Returns:
        String All tests passed if all tests are successful, otherwise raises an exception
    """
    assert alignment_edges("ACGT", "ACGT") == [(1,1), (2, 2), (3, 3), (4, 4)], f"Identical sequences test failed, edges should be [(1,1), (2, 2), (3, 3), (4, 4)] and not {alignment_edges('ACGT', 'ACGT')}"
    assert alignment_edges("ACGT", "AC-T") == [(1,1), (2, 2), (4, 3)], f"Half identical sequences test failed, edges should be [(1,1), (2, 2), (4, 3)] and not {alignment_edges('ACGT', 'AC-T')}"
    assert alignment_edges("ACGT", "----") == [], f"No identical sequences test failed, edges should be [] and not {alignment_edges('ACGT', '----')}"
    return "Alignment edges tests passed"

def test_edges_to_dot_bracket():
    """
    Test function for edges_to_dot_bracket

    Returns:
        String All tests passed if all tests are successful, otherwise raises an exception
    """
    assert edges_to_dot_bracket([(2, 8), (3, 7), (9, 13)], 13) == '.((...))(...)', f"Test case 1 failed, expected .((...))(...) but got {edges_to_dot_bracket([(2, 8), (3, 7), (9, 13)], 13)}"
    assert edges_to_dot_bracket([],6) == '......', f"Test case 2 failed, expected ...... but got {edges_to_dot_bracket([],6)}"
    assert edges_to_dot_bracket([(1, 6)],6) == '(....)', f"Test case 3 failed, expected (....) but got {edges_to_dot_bracket([(1, 6)],6)}"
    try:
        edges_to_dot_bracket([(1, 1)],1)
    except AssertionError:
        pass
    else:
        raise AssertionError(f"Test case 4 failed, expected an assertion error but got {edges_to_dot_bracket([(1, 1)],1)}")
    try:
        edges_to_dot_bracket([(1, 2)],1)
    except AssertionError:
        pass
    else:
        raise AssertionError(f"Test case 5 failed, expected an assertion error but got {edges_to_dot_bracket([(1, 2)],1)}")
    try:
        edges_to_dot_bracket([(2, 1)],2)
    except AssertionError:
        pass
    else:
        raise AssertionError(f"Test case 6 failed, expected an assertion error but got {edges_to_dot_bracket([(2, 1)],2)}")

    return "Edges to dot bracket tests passed"

def test_parse_rna_structure():
    """
    Test function for parse_rna_structure

    Returns:
        String All tests passed if all tests are successful, otherwise raises an exception
    """
    assert parse_rna_structure('.((...))(...)') == [(2, 8), (3, 7), (9, 13)], f"Test case 1 failed, expected [(2, 8), (3, 7), (9, 13)] but got {parse_rna_structure('.((...))(...)')}"
    assert parse_rna_structure('......') == [], f"Test case 2 failed, expected [] but got {parse_rna_structure('......')}"
    assert parse_rna_structure('(....)') == [(1, 6)], f"Test case 3 failed, expected [(1, 6)] but got {parse_rna_structure('(....)')}"
    try:
        parse_rna_structure('..(..')
    except AssertionError:
        pass
    else:
        raise AssertionError(f"Test case 4 failed, expected an assertion error but got {parse_rna_structure('..(..')}")
    try:
        parse_rna_structure('..)..')
    except AssertionError:
        pass
    else:
        raise AssertionError(f"Test case 5 failed, expected an assertion error but got {parse_rna_structure('..)..')}")
    
    return "Parse RNA structure tests passed"

def test_is_canonical_base_pair():
    """
    Test function for is_canonical_base_pair

    Returns:
        String All tests passed if all tests are successful, otherwise raises an exception
    """
    assert is_canonical_base_pair('A', 'U') == True, f"Test case 1 failed, expected True but got {is_canonical_base_pair('A', 'U')}"
    assert is_canonical_base_pair('U', 'A') == True, f"Test case 2 failed, expected True but got {is_canonical_base_pair('U', 'A')}"
    assert is_canonical_base_pair('G', 'C') == True, f"Test case 3 failed, expected True but got {is_canonical_base_pair('G', 'C')}"
    assert is_canonical_base_pair('C', 'G') == True, f"Test case 4 failed, expected True but got {is_canonical_base_pair('C', 'G')}"
    assert is_canonical_base_pair('G', 'U') == True, f"Test case 5 failed, expected True but got {is_canonical_base_pair('G', 'U')}"
    assert is_canonical_base_pair('U', 'G') == True, f"Test case 6 failed, expected True but got {is_canonical_base_pair('U', 'G')}"
    assert is_canonical_base_pair('A', 'C') == False, f"Test case 7 failed, expected False but got {is_canonical_base_pair('A', 'C')}"
    assert is_canonical_base_pair('C', 'A') == False, f"Test case 8 failed, expected False but got {is_canonical_base_pair('C', 'A')}"
    assert is_canonical_base_pair('G', 'A') == False, f"Test case 9 failed, expected False but got {is_canonical_base_pair('G', 'A')}"
    assert is_canonical_base_pair('A', 'G') == False, f"Test case 10 failed, expected False but got {is_canonical_base_pair('A', 'G')}"
    assert is_canonical_base_pair('U', 'C') == False, f"Test case 11 failed, expected False but got {is_canonical_base_pair('U', 'C')}"
    assert is_canonical_base_pair('C', 'U') == False, f"Test case 12 failed, expected False but got {is_canonical_base_pair('C', 'U')}"
    assert is_canonical_base_pair('A', 'A') == False, f"Test case 13 failed, expected False but got {is_canonical_base_pair('A', 'A')}"
    assert is_canonical_base_pair('C', 'C') == False, f"Test case 14 failed, expected False but got {is_canonical_base_pair('C', 'C')}"
    assert is_canonical_base_pair('G', 'G') == False, f"Test case 15 failed, expected False but got {is_canonical_base_pair('G', 'G')}"
    assert is_canonical_base_pair('U', 'U') == False, f"Test case 16 failed, expected False but got {is_canonical_base_pair('U', 'U')}"
    try:
        is_canonical_base_pair('P', 'T')
    except AssertionError:
        pass
    else:
        raise AssertionError(f"Test case 17 failed, expected an assertion error but got {is_canonical_base_pair('A', 'T')}")
    return "Is canonical base pair tests passed"

def test_is_canonical():
    """
    Test function for is_canonical

    Returns:
        String All tests passed if all tests are successful, otherwise raises an exception
    """
    assert is_canonical('......', 'ACGUAC') == True, f"Test case 1 failed, expected True but got {is_canonical('......', 'ACGUAC')}"
    assert is_canonical('(....)', 'ACGUAC') == False, f"Test case 2 failed, expected False but got {is_canonical('(....)', 'ACGUAC')}"
    assert is_canonical('(....)', 'ACGUAU') == True, f"Test case 3 failed, expected True but got {is_canonical('(....)', 'ACGUAC')}"
    assert is_canonical('((....))', 'ACGGGUAC') == False, f"Test case 4 failed, expected False but got {is_canonical('((....))', 'ACGGGUAC')}"
    assert is_canonical('((....))', 'ACGGGUGU') == True, f"Test case 5 failed, expected True but got {is_canonical('((....))', 'ACGGGUGU')}"
    assert is_canonical('.((...))(...)', 'ACGUACCGACGGU') == True, f"Test case 6 failed, expected True but got {is_canonical('.((...))(...)', 'ACGUACGGACGGU')}"
    assert is_canonical('.((...))(...)', 'ACGUACGUACGGC') == False, f"Test case 7 failed, expected False but got {is_canonical('.((...))(...)', 'ACGUACGUACGGC')}"
    try:
        is_canonical('.((...))(...)', 'ACGUACG')
    except AssertionError:
        pass
    else:
        raise AssertionError(f"Test case 8 failed, expected an assertion error but got {is_canonical('.((...))(...)', 'ACGUACGUACGGCU')}")
   
    return "RNA structure is canonical tests passed"

