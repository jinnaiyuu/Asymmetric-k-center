# Asymmetric-k-center

Implementation of an O(log* k) approximation algorithm (Archer 2001) for asymmetric k-center problem.

(Archer 2001) Archer, Aaron. "Two O (log* k)-approximation algorithms for the asymmetric k-center problem." International Conference on Integer Programming and Combinatorial Optimization. Springer, Berlin, Heidelberg, 2001.


# Usage

AKC_APPROX(D, k) returns the approximate solution for the asymmetric k-center problem using Archer 2001 algorithm.

AKC_OPT(D, k) returns the exact optimal assignment of the asymmetric k-center problem using a CSP solver, interfaced by Google's ortools. To run this you need to install ortools.

```
	pip install numpy scipy ortools
```

You do not need ortools to run AKC_APPROX.


# Author

Yuu Jinnai <ddyuudd@gmail.com>
