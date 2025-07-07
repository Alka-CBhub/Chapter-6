"""
Functios to convert the "implicit equations" produced by SINDy-PI into explicit ODEs for derivatives:

-----------
Workflow:
-----------
1.  feature_names     <--  model.get_feature_names()
2.  coeffs            <--   model.coefficients()
3.  left_coeff        <--  identity matrix  (makes each row start with its own feature)
4.  eq_list, symbols  <--   generate_symbolic_equations_both_sides(...)
"""

import re
import numpy as np
import sympy as sp


# -------------------------------------------------------------------------
# 1.  Token-level parsing helpers
# -------------------------------------------------------------------------
def parse_feature(feature: str):
    """
    Split a raw feature-string such as 'x0x0x1_dot'
    into a tuple of 'tokens'  ('x0','x0','x1_dot') that preserves
    any derivative suffix (_dot or _t).
    """
    # locate every “…_dot”  or  “…_t”            e.g. x1_dot or x1_t
    matches = list(re.finditer(r"([a-zA-Z]\d*)(_dot|_t)", feature))

    if not matches:                                # plain monomial like 'x0x1'
        return tuple(re.findall(r"[a-zA-Z]\d*", feature))

    tokens, last = [], 0
    for m in matches:
        # prefix  = text before this _dot token
        prefix = feature[last:m.start()]
        if prefix:
            tokens.extend(re.findall(r"[a-zA-Z]\d*", prefix))

        base, deriv = m.groups()                   # e.g. ('x1', '_dot')
        tokens.append(base + deriv)                # keep derivative suffix
        last = m.end()

    # any trailing part after the final _dot / _t
    if last < len(feature):
        tokens.extend(re.findall(r"[a-zA-Z]\d*", feature[last:]))

    return tuple(tokens)

# -------------------------------------------------------------------------
# Above works for variables like a to z, A to Z, each alphabet with numeric 
#.....indices as state variables can take any variable name.

# If the string contains no “_dot”/“_t”, it simply extracts every alphanumeric
#.....  variable block.

################################################################################################################################


# -------------------------------------------------------------------------
# 2.  Build a vocabulary of unique tokens and a map  feature -> token list
# -------------------------------------------------------------------------
def extract_distinct_features(feature_names):
    """
    Parameters
    ----------
    feature_names : list[str]
        Raw strings from model.get_feature_names().

    Returns
    -------
    unique_tokens : list[str]
        Sorted list of every distinct variable / derivative token.
    feature_map   : dict[str, tuple[str]]
        For each feature name, its ordered tuple of tokens.
    """
    unique_tokens, feature_map = set(), {}
    for feat in feature_names:
        toks = parse_feature(feat)
        feature_map[feat] = toks
        unique_tokens.update(toks)
    return sorted(unique_tokens), feature_map


################################################################################################################################
# -------------------------------------------------------------------------
# 3.  Convert a token list -> SymPy product, automatically handling powers
# -------------------------------------------------------------------------
def consolidate_product(tokens, token_symbols):
    """
    ('x0','x0','x1_dot')  →  x0**2 * x1_dot
    """
    counts = {}
    for t in tokens:
        counts[t] = counts.get(t, 0) + 1

    return sp.Mul(*(
        token_symbols[t] ** p for t, p in counts.items()
    ), evaluate=False)                            # keep symbolic power

################################################################################################################################
# -------------------------------------------------------------------------
# 4.  Main routine: produce SymPy Eq objects for every row of the implicit SINDy model
# -------------------------------------------------------------------------
def generate_symbolic_equations_both_sides(feature_names,
                                           left_coeff,
                                           right_coeff):
    """
    Given the *coefficient matrices* from SINDy-PI, build:

        Eq(expression, 0)

     left_coeff  : e.g. identity matrix  -> LHS is that row’s own feature.
     right_coeff : np.ndarray  (same shape) of sparse coefficients.

    Returns
    -------
    eq_list      : list[sympy.Eq],   length = n_features
    token_symbols: dict[token -> SymPy Symbol]
                   so that we can later access token_symbols['x_dot'] etc.
    """
    distinct, fmap = extract_distinct_features(feature_names)
    token_symbols  = {tok: sp.Symbol(tok) for tok in distinct}
    n = len(feature_names)

    L = sp.Matrix(left_coeff)
    R = sp.Matrix(right_coeff)

    eq_list = []
    for i in range(n):
        lhs = sum(L[i, j] * consolidate_product(fmap[feature_names[j]], token_symbols)
                  for j in range(n))
        rhs = sum(R[i, j] * consolidate_product(fmap[feature_names[j]], token_symbols)
                  for j in range(n))
        eq_list.append(sp.Eq(lhs - rhs, 0))              # implicit row i

    return eq_list, token_symbols


# -------------------------------------------------------------------------
# 5.  Produce an identity of proper size
# -------------------------------------------------------------------------
def generate_identity_matrix(features):
    """
    Convenience:  np.eye(len(features)) so that
    each implicit equation starts with its “own” feature on the LHS.
    """
    return np.eye(len(features), dtype=int)

# -------------------------------------------------------------------------
# 6.  Produce reformatted features
# -------------------------------------------------------------------------
def get_reformatted_feature_names(feature_names):
    """
    Return reformatted feature names as simplified symbolic strings
    """
    distinct_tokens, feature_map = extract_distinct_features(feature_names)
    token_symbols = {token: sp.Symbol(token) for token in distinct_tokens}

    reformatted = []
    for feat in feature_names:
        tokens = feature_map[feat]
        sym_expr = consolidate_product(tokens, token_symbols)
        reformatted.append(str(sym_expr))
    return reformatted



































