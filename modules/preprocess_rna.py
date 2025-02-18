def preprocess_rna_seq(seq: str) -> str:
    """
    In publicly available RNA sequence data, one will often find small letter symbols and even T instead of U. 
    Preprocess RNA sequence by uppercasing and converting T to U.

    Args:
        seq: RNA sequence

    Returns:
        preprocessed RNA sequence    
    """
    return seq.upper().replace('T', 'U')

def sequence_identity(xp:str, yp:str, need_preprocess = True) -> str:
    """
    Sequence identity of two aligned sequences

    Args:
      xp: alignment string of sequence x
      yp: alignment string of sequence y
      need_preprocess: whether to preprocess the sequences
    
    Returns:
      Percentage of identical alignment columns
    """
    assert len(xp) == len(yp), "Sequences must have the same length"
    if need_preprocess:
        xp = preprocess_rna_seq(xp)
        yp = preprocess_rna_seq(yp)

    xp_symbols = [x for x in xp]
    yp_symbols = [y for y in yp]
    return f"{sum([1 for x, y in zip(xp_symbols, yp_symbols) if x == y]) / len(xp) * 100:.2f}%"

