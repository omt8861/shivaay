# Stoichiometry Solver (Pure Python, Fully Working)
# Author: ChatGPT (GPT-5)

import re
from fractions import Fraction

# --- Parse formula into element counts ---
def parse_formula(formula):
    tokens = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
    composition = {}
    for el, num in tokens:
        composition[el] = composition.get(el, 0) + (int(num) if num else 1)
    return composition

# --- Build element matrix ---
def build_matrix(reactants, products):
    elements = sorted({el for cmpd in reactants + products for el in parse_formula(cmpd)})
    matrix = []
    for el in elements:
        row = []
        for r in reactants:
            row.append(parse_formula(r).get(el, 0))
        for p in products:
            row.append(-parse_formula(p).get(el, 0))
        matrix.append(row)
    return matrix

# --- Gaussian elimination (fraction-safe) ---
def gaussian_elimination(matrix):
    matrix = [list(map(Fraction, row)) for row in matrix]
    rows, cols = len(matrix), len(matrix[0])
    for i in range(min(rows, cols)):
        # find pivot
        if matrix[i][i] == 0:
            for j in range(i+1, rows):
                if matrix[j][i] != 0:
                    matrix[i], matrix[j] = matrix[j], matrix[i]
                    break
        pivot = matrix[i][i]
        if pivot == 0:
            continue
        # normalize row
        matrix[i] = [x / pivot for x in matrix[i]]
        # eliminate others
        for j in range(rows):
            if j != i:
                factor = matrix[j][i]
                matrix[j] = [a - factor*b for a, b in zip(matrix[j], matrix[i])]
    return matrix

# --- Extract integer coefficients ---
def extract_coeffs(matrix):
    rows, cols = len(matrix), len(matrix[0])
    coeffs = [0]*cols
    coeffs[-1] = 1  # free variable = 1
    for i in range(rows-1, -1, -1):
        s = sum(matrix[i][j]*coeffs[j] for j in range(i+1, cols))
        coeffs[i] = -s / matrix[i][i] if matrix[i][i] != 0 else 0
    # Convert to integers
    lcm = 1
    for c in coeffs:
        if isinstance(c, Fraction):
            lcm = lcm * c.denominator // gcd(lcm, c.denominator)
    coeffs = [int(c*lcm) for c in coeffs]
    # normalize sign
    if all(c <= 0 for c in coeffs):
        coeffs = [-c for c in coeffs]
    return coeffs

def gcd(a, b):
    while b:
        a, b = b, a % b
    return abs(a)

# --- Balance equation ---
def balance_equation(reactants, products):
    matrix = build_matrix(reactants, products)
    reduced = gaussian_elimination(matrix)
    coeffs = extract_coeffs(reduced)
    return coeffs

# --- Stoichiometry calculator ---
def stoichiometry():
    print("\n--- Stoichiometry Solver ---")
    reactants = input("Enter reactants (comma-separated): ").replace(" ", "").split(",")
    products = input("Enter products (comma-separated): ").replace(" ", "").split(",")

    coeffs = balance_equation(reactants, products)
    lhs = " + ".join(f"{coeffs[i]} {reactants[i]}" for i in range(len(reactants)))
    rhs = " + ".join(f"{coeffs[len(reactants)+j]} {products[j]}" for j in range(len(products)))
    print(f"\n✅ Balanced Equation:\n{lhs} -> {rhs}")

    all_cmpds = reactants + products
    given = input("\nEnter known compound: ").strip()
    if given not in all_cmpds:
        print("❌ Compound not found in equation.")
        return
    target = input("Enter target compound: ").strip()
    if target not in all_cmpds:
        print("❌ Compound not found in equation.")
        return
    moles_given = float(input(f"Enter moles of {given}: "))

    given_idx = all_cmpds.index(given)
    target_idx = all_cmpds.index(target)
    ratio = coeffs[target_idx] / coeffs[given_idx]
    moles_target = moles_given * ratio
    print(f"\n➡️ {moles_given} mol of {given} → {moles_target:.3f} mol of {target}")

# --- Main ---
if __name__ == "__main__":
    stoichiometry()
