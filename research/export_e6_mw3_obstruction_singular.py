#!/usr/bin/env python3
"""Emit certified E6/MW3 obstruction ideals in Singular syntax."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

import hensel_lift_e6_mw3_gf29_flexible_slice as flexible
import solve_e6_mw3_gf29_second_order as second


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "certificates" / "e6_mw3_gf29_independence.json"


def singular_expression(expression):
    return str(expression).replace("**", "^")


def emit(
    input_path,
    p3_index,
    add_field_equations,
    ring_order,
    variable_order,
    print_basis,
    specialization,
):
    data = json.loads(input_path.read_text())
    record = second.obstruction_model(data, p3_index, system=flexible)
    variables, polynomials = second.sympy_polynomials(record)
    substitutions = {}
    if specialization:
        by_name = {str(variable): variable for variable in variables}
        for item in specialization.split(","):
            name, value = item.split("=", 1)
            substitutions[by_name[name]] = int(value) % 29
        polynomials = [sp.expand(poly.subs(substitutions)) for poly in polynomials]
        variables = tuple(variable for variable in variables if variable not in substitutions)
    names = variable_order or ",".join(map(str, variables))
    equations = [singular_expression(poly) for poly in polynomials]
    if add_field_equations:
        equations.extend(f"{variable}^29-{variable}" for variable in variables)
    print(f"ring r=29,({names}),{ring_order};")
    print("option(redSB);")
    print("ideal I=" + ",\n".join(equations) + ";")
    print("int start_timer=timer;")
    print("ideal G=std(I);")
    print('print("elapsed_ms="); print(timer-start_timer);')
    print('print("basis_size="); print(size(G));')
    print('print("dimension="); print(dim(G));')
    print('print("independent_set="); print(indepSet(G));')
    if add_field_equations:
        print('print("vector_space_dimension="); print(vdim(G));')
    print('print("remainder_of_one="+string(reduce(1,G)));')
    if print_basis:
        print('print("basis="); print(G);')
    print("quit;")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--p3", type=int, required=True)
    parser.add_argument("--field-equations", action="store_true")
    parser.add_argument("--ring-order", default="dp")
    parser.add_argument("--variable-order")
    parser.add_argument("--print-basis", action="store_true")
    parser.add_argument("--specialize")
    args = parser.parse_args()
    emit(
        args.input,
        args.p3,
        args.field_equations,
        args.ring_order,
        args.variable_order,
        args.print_basis,
        args.specialize,
    )


if __name__ == "__main__":
    main()
