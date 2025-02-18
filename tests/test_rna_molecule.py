import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from modules.rna_molecule import *



def test_class_rna():
    """
    Test function for RNAMolecule class

    Returns:
        String All tests passed if all tests are successful, otherwise raises an exception
    """
    rna_test = RNAMolecule("ACGUAC")
    assert rna_test.seq == "ACGUAC", f"Test case 1 failed, expected ACGUAC but got {rna_test.seq}"
    assert rna_test.structure == None, f"Test case 2 failed, expected None but got {rna_test.structure}"

    rna_test = RNAMolecule("ACGUAC", "......")
    assert rna_test.seq == "ACGUAC", f"Test case 3 failed, expected ACGUAC but got {rna_test.seq}"
    assert rna_test.structure == "......", f"Test case 4 failed, expected ...... but got {rna_test.structure}"

    rna_test = RNAMolecule("ACGUAC", "(....)")
    assert rna_test.seq == "ACGUAC", f"Test case 5 failed, expected ACGUAC but got {rna_test.seq}"
    assert rna_test.structure == "(....)", f"Test case 6 failed, expected (....) but got {rna_test.structure}"
    assert str(rna_test) == "Sequence: ACGUAC\nStructure: (....)", f"Test case 7 failed, expected Sequence: ACGUAC\nStructure: (....) but got {str(rna_test)}"
    print("Printing test:\n", rna_test, sep="")

    rna_test = RNAMolecule("acguac", "(....)")
    assert rna_test.seq == "ACGUAC", f"Test case 8 failed, expected ACGUAC but got {rna_test.seq}"
    assert rna_test.structure == "(....)", f"Test case 9 failed, expected (....) but got {rna_test.structure}"
    assert rna_test.is_canonical() == False, f"Test case 10 is canonical failed, expected False but got {rna_test.is_canonical()}"
    assert rna_test.parse_rna_structure() == [(1, 6)], f"Test case 11 failed, expected [(1, 6)] but got {rna_test.parse_rna_structure()}"

    rna_test = RNAMolecule('AC--G-T', "((..).)")
    assert rna_test.rna_structure_projection() == ('ACGU', '(())'), f"Expected ('ACGT', '(())') but got {rna_structure_projection('AC--G-T', '((..).)')}"

    rna_test = RNAMolecule('--CA---T', '(.)((.))') 
    assert rna_test.rna_structure_projection() == ('CAU', '.()'), f"Expected ('CAT', '.()') but got {rna_structure_projection('--CA---T', '(.)((.))')}"

    rna_test = RNAMolecule('GCC-CUUAG-U-GAAUCCAGC', '((.((...))(((...)))))')
    assert rna_test.rna_structure_projection() == ('GCCCUUAGUGAAUCCAGC', '((.(...)((...).)))'), f"Expected ('GCCCUUAGUGAAUCCAGC', '((.(...)((...).)))') but got {rna_structure_projection('GCC-CUUAG-U-GAAUCCAGC', '((.((...))(((...)))))')}"
    return "All tests passed"