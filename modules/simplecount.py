complementary_bases=["AT","CG","AU","UG"]

def simple_count_cost (base1, base2):
    """ computes wether two RNA bases (a,u,c,g or t) are complementary
    IN : two car base 1 and base 2, upper case
    Returns : int, 1 if complementary, 0 otherwise
    """
    if base1+base2 in complementary_bases or base2+base1 in complementary_bases:
        return 1
    return 0
