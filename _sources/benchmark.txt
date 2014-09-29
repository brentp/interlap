Benchmarking
============

Here we will benchmark `InterLap`

**TLDR** `InterLap` is:

- 40 times faster than naive method with no memory overhead
- 20 times slower than the bx-python interval tree written in cython


Create 50K random intervals of random lengths between 2KB and 20KB across
100MB of chromosome:

.. ipython::

    In [1]: import random, sys

    In [1]: from interlap import InterLap

    In [2]: random.seed(42)

    In [3]: sites = random.sample(range(22, 100000000, 12), 50000)

    In [4]: sites = [(s, s + random.randint(2000, 20000)) for s in sites]

    In [4]: %timeit InterLap(sites)

    In [4]: inter = InterLap(sites)

    In [5]: len(inter)

Now we have the searchable object `inter` with 50K intervals. Since `InterLap` is 
a sorted list, time to create is simply time to create a sorted list.

We can time to see how fast we can test every site in the tree

.. ipython::

    In [6]: %timeit assert all(r in inter for r in sites)

And make sure we dont see intervals at positions < 0:

.. ipython::

    In [7]: %timeit assert not any((-r[0], -r[1]) in inter for r in sites)


Now do 100K (10K * 10) queries.

.. ipython::

    In [8]: def times(inter, test_intervals, n_times):
       ...:     for i in xrange(n_times):
       ...:         a = sum(t in inter for t in test_intervals)
       ...:
    

    In [9]: test_sites = random.sample(range(21, 100000000, 12), 10000)

    In [10]: test_sites = [(t, t + random.randint(2000, 200000)) for t in test_sites]

    In [11]: %timeit times(inter, test_sites, 10)

Compared to Naive Method
------------------------

We can compare to the naive method (only 1 iteration since it will be slow):

.. ipython::

   In [1]: def naivetimes(query_intervals, test_intervals):
      ...:     for t in test_intervals:
      ...:         isin = any(q[1] > t[0] and q[0] < t[1] for q in query_intervals)
      ...:
   
   In [2]: %timeit naivetimes(sites, test_sites)

So `InterLap` is over 40 times (remember we did only 1 rep here instead of 10)
as fast as the naive method with *no memory overhead*.

Compared to bx-python interval-tree
-----------------------------------

We can compare this to bx-python which uses a tree structure written in cython.

.. ipython::

    In [1]: from bx.intervals import IntervalTree

    In [2]: %timeit IntervalTree(sites)

    In [3]: tree = IntervalTree(sites)

    In [4]: def bxtimes(tree, test_intervals, n_times):
       ...:     for i in xrange(n_times):
       ...:         a = sum(len(tree.find(ts[0], ts[1])) > 0 for ts in test_intervals)
       ...:

    In [5]: %timeit bxtimes(tree, test_sites, 10)

As expected, bx-python is much faster, but `InterLap` does perform quite well
at ~50K queries per second.


`InterLap` will come within 2X of bx-python if most of the query intervals are
found in the database.

