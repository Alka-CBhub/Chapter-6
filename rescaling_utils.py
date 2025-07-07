# Import necessary libraries
import sympy as sp
import pandas as pd
from sympy.printing.latex import latex
from typing import Optional, List, Dict, Union, Tuple

###################################################################
# Utility 1: Rational Simplification
###################################################################
def rational_simplify(expression, expand_denominator=False):
    combined = sp.cancel(sp.together(expression))  # Combine and cancel
    numerator, denominator = combined.as_numer_denom()
    if expand_denominator:
        denominator = sp.expand(denominator)
    return numerator, denominator


def _expanded_fraction(expr):
    num, den = sp.fraction(expr)
    return sp.expand(num) / sp.expand(den)

###################################################################
# Utility 2: Drop Small Terms
###################################################################
def drop_small_terms(expr: sp.Expr, tol: float = 1e-6) -> sp.Expr:
    num, den = expr.as_numer_denom()

    # Expand first
    num = sp.expand(num)
    den = sp.expand(den)

    def filter_terms(terms):
        new_terms = []
        for t in terms:
            coeff, _ = t.as_coeff_Mul()
            if coeff.is_number and not coeff.is_zero and abs(coeff.evalf()) > tol:
                new_terms.append(t)
        return new_terms

    num_terms = filter_terms(num.as_ordered_terms())
    den_terms = filter_terms(den.as_ordered_terms())

    # Optional fallback to preserve model
    if not num_terms and num.as_ordered_terms():
        num_terms = [max(num.as_ordered_terms(), key=lambda t: abs(t.as_coeff_Mul()[0].evalf()))]
    if not den_terms and den.as_ordered_terms():
        den_terms = [max(den.as_ordered_terms(), key=lambda t: abs(t.as_coeff_Mul()[0].evalf()))]

    num_clean = sp.Add(*num_terms)
    den_clean = sp.Add(*den_terms)

    if den_clean.is_zero:
        return sp.zoo  # invalid model

    return _expanded_fraction(num_clean / den_clean)



###################################################################
# Utility 3: Rescaling
###################################################################
def rescale_expression(expr: sp.Expr,
                       *,
                       target: Optional[sp.Expr] = None,
                       gens: Optional[Union[List[sp.Symbol], tuple]] = None,
                       verbose: bool = False) -> sp.Expr:
    expr = sp.together(expr)
    N, D = expr.as_numer_denom()

    if D.is_number:
        if verbose:
            print("[Skip] Denominator is constant.")
        return expr

    gens = tuple(gens) if gens else tuple(sorted(D.free_symbols, key=lambda s: s.name))
    scaling_coeff = None
    used_monomial = None

    try:
        poly = sp.Poly(D, *gens)

        # === Try to use target if provided ===
        if target is not None:
            coeff = poly.coeff_monomial(sp.expand(target))
            if (
                coeff.is_number and
                not coeff.has(sp.nan) and
                not coeff.is_zero and
                abs(coeff.evalf()) > 1e-12
            ):
                scaling_coeff = coeff
                used_monomial = target
                if verbose:
                    print(f"[Target Found] Using monomial {target} with coefficient {scaling_coeff}")
            else:
                if verbose:
                    print(f"[Target Invalid] Target {target} missing or zero. Falling back to leading term.")

        # === Fallback to leading monomial ===
        if scaling_coeff is None:
            monoms = poly.monoms()
            if not monoms:
                if verbose:
                    print("[Fail] No monomials in denominator.")
                return expr

            total_degrees = [sum(m) for m in monoms]
            max_deg = max(total_degrees)
            lead = sorted([m for m, d in zip(monoms, total_degrees) if d == max_deg])[0]
            leading_term = sp.prod(v**e for v, e in zip(poly.gens, lead))
            coeff = poly.coeff_monomial(leading_term)
            if (
                coeff.is_number and
                not coeff.has(sp.nan) and
                not coeff.is_zero and
                abs(coeff.evalf()) > 1e-12
            ):
                scaling_coeff = coeff
                used_monomial = leading_term
                if verbose:
                    print(f"[Leading Term] Using monomial {leading_term} with coefficient {scaling_coeff}")
            else:
                if verbose:
                    print("[Abort] Leading term coefficient is invalid.")
                return expr

    except Exception as e:
        if verbose:
            print(f"[Error] Poly failed: {e}")
        return expr

    # Final rescaling
    if verbose:
        print(f"[Success] Rescaling by coefficient {scaling_coeff}")
    rescaled_expr = (N / scaling_coeff) / (D / scaling_coeff)
    return _expanded_fraction(rescaled_expr)



###################################################################
# Utility 4: Print Rational Model
###################################################################
def print_rational_eq(i: int, var: sp.Symbol, expr: sp.Expr) -> None:
    num, den = sp.fraction(expr)
    num, den = sp.expand(num), sp.expand(den)
    block = [
        "=" * 80,
        f"Model {i}:",
        " " * (len(str(var)) + 4) + str(num),
        f"{var} = " + "-" * max(len(str(num)), len(str(den))),
        " " * (len(str(var)) + 4) + str(den),
        f"[Terms: numerator={len(num.as_ordered_terms())}, denominator={len(den.as_ordered_terms())}]",
        "=" * 80
    ]
    print("\n".join(block))

###################################################################
# Utility 5: Solve, Simplify, and Export Models
###################################################################
def print_and_store_models(eq_list: List[sp.Eq],
                           xdot_symbol: sp.Symbol,
                           *,
                           target: Optional[sp.Expr] = None,
                           gens: Optional[Union[List[sp.Symbol], tuple]] = None,
                           sig_digits: int = 4,
                           tol: float = 1e-6,
                           latex_path: Optional[str] = None,
                           csv_path: Optional[str] = None
                          ) -> Tuple[Dict[int, sp.Expr], Dict[int, sp.Expr], Dict[int, sp.Expr], Dict[int, sp.Expr]]:

    raw_models = {}
    rescaled_models = {}
    cleaned_models = {}
    final_models = {}

    failed_models = 0
    successful_models = 0
    rows_for_csv = []
    latex_file = open(latex_path, "w") if latex_path else None

    for i, eq in enumerate(eq_list):
        sols = sp.solve(eq, xdot_symbol)
        if not sols:
            print(f"Model {i}: No solution.")
            failed_models += 1
            continue

        raw = sols[0]
        raw_expanded = sp.together(sp.expand(raw))
        raw_models[i] = raw_expanded # store raw

        num, _ = sp.fraction(raw_expanded)
        if getattr(num, "is_zero", False):
            print(f"Model {i}: Solver failed (zero numerator).")
            failed_models += 1
            continue

        # Step 1: Rescale
        rescaled = rescale_expression(raw_expanded, target=target, gens=gens)
        rescaled_models[i] = rescaled   # store rescaled
    

        # Step 2: Drop small terms BEFORE rounding
        cleaned = drop_small_terms(rescaled, tol=tol)
        num_clean, den_clean = sp.fraction(cleaned)
        if getattr(num_clean, "is_zero", False):
            print(f"Model {i}: All numerator terms dropped (zero).")
            failed_models += 1
            continue

        # Step 3: Round remaining terms
        final_expr = _expanded_fraction(cleaned.evalf(sig_digits))

        # Store all surviving stages
        cleaned_models[i]  = cleaned
        final_models[i]    = final_expr
        successful_models += 1

        # Print and export
        print_rational_eq(i, xdot_symbol, final_expr)

        if latex_file:
            latex_expr = sp.latex(sp.Eq(xdot_symbol, final_expr))
            latex_file.write(f"% Model {i}\n\\[\n{latex_expr}\n\\]\n\n")

        has_negative = any(
            term.as_coeff_Mul()[0] < 0
            for term in den_clean.as_ordered_terms()
            if term.as_coeff_Mul()[0].is_number
        )

        rows_for_csv.append({
            "Model": i,
            "Equation": f"{xdot_symbol} = {final_expr}",
            "Numerator Terms": len(num_clean.as_ordered_terms()),
            "Denominator Terms": len(den_clean.as_ordered_terms()),
            "Any Negative in Denominator": has_negative
        })

    if latex_file:
        latex_file.close()

    if csv_path:
        df = pd.DataFrame(rows_for_csv)
        df.to_csv(csv_path, index=False)

    print("\n======================")
    print(f" Total models processed: {len(eq_list)}")
    print(f" Models obtained      : {successful_models}")
    print(f" Models failed        : {failed_models}")
    print("======================\n")

    return raw_models, rescaled_models, cleaned_models, final_models
