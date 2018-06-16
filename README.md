# Asymmetric-k-center

Implementation of an O(log* k) approximation algorithm (Archer 2001) for asymmetric k-center problem.

The input to the problem is a directed graph G and a distance function d: V x V -> R.
The distance function must satisfy triangle inequality.
The goal is to find a set C of k vertices, or centers, such that the maximal distance of a vertex from its center is minimized.
That is, we wish to minimize the radius R:

R = max {s in V} min {t in C} d(s, t)

The problem is a generalization of the k-center problem but the distance function can be asymmetric.
It is known to be an NP-hard problem (Lewis 1983) even for a (symmetric) k-center problem.
The problem is hard (NP-hard) to approximate up to a factor of log* n - \theta(1), where n is the number of nodes (Chuzhoy 2004).
Thus, the algorithm by Archer (2001) is the best known algorithm in terms of time complexity.

The algorithm finds a solution with radius R within r*(3 log*k + O(1)), where r* is the optimum.

# Usage

AKC_APPROX(D, k) returns the approximate solution for the asymmetric k-center problem using Archer's algorithm.

AKC_OPT(D, k) returns the exact optimal assignment of the asymmetric k-center problem using a CSP solver, interfaced by Google's ortools. To run this you need to install ortools.

```
	pip install numpy scipy ortools
```

You do not need ortools to run AKC_APPROX.

# Literature

(Archer 2001) Archer, Aaron. "Two O (log* k)-approximation algorithms for the asymmetric k-center problem." International Conference on Integer Programming and Combinatorial Optimization. Springer, Berlin, Heidelberg, 2001.

(Lewis 1983) Lewis, Harry R. "Computers and Intractability. A Guide to the Theory of NP-Completeness." (1983): 498-500.

(Chuzhoy 2004) Chuzhoy, Julia, et al. "Asymmetric k-center is log* n-hard to approximate." Journal of the ACM (JACM) 52.4 (2005): 538-551.

(Panigrahy 1998) Panigrahy, Rina, and Sundar Vishwanathan. "An O(log* n) Approximation Algorithm for the Asymmetric p-Center Problem." Journal of Algorithms 27.2 (1998): 259-268.

# Author

Yuu Jinnai <ddyuudd@gmail.com>
