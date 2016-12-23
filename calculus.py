#!/usr/bin/python
from sympy import diff, Eq, init_printing, pprint, symbols, simplify
init_printing()
C, Ca, Cb, f, F, rhoS  = symbols('C Ca Cb f F rhoS')

chemenergy = rhoS * (C - Ca)**2 * (Cb - C)**2

pprint(Eq(F,simplify(chemenergy)))
pprint(Eq(f,simplify(diff(chemenergy, C))))
