"""qmsd -- discovery of qudit triorthogonal codes for magic state distillation.

Reimplements and extends the search of arXiv:2510.10852 (Saha & Prakash, 2025):
given a prime qudit dimension p, build triorthogonal codes by puncturing
maximal-degree Reed-Muller codes over F_p and score them by the yield parameter
gamma = log(n/k)/log(d) and by single-round distillation cost C.

Authoritative math reference: arXiv:2510.10852
Software design:               ../IMPLEMENTATION_BLUEPRINT.md

This package keeps the submodules import-light on purpose (no heavy re-exports
at package import time) so partial builds stay importable.
"""

__version__ = "0.1.0"
