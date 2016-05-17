InterLap: simple, fast interval overlap testing
-----------------------------------------------
[ ![Codeship Status for brentp/interlap](https://www.codeship.io/projects/b4029ac0-2a1a-0132-a72c-3a1c9f74327f/status)](https://www.codeship.io/projects/38165)

InterLap is >20 times faster than doing a naive search (see: https://brentp.github.io/interlap/benchmark.html)
with **no memory overhead** because it operates on a sorted list. It is pure python and has no
dependencies.

It uses binary search and a knowledge of the longest observed interval to quickly query datasets
with 100's of thousands of intervals.

See the documentation at [https://brentp.github.io/interlap/](https://brentp.github.io/interlap)

Usage
-----

Interlap takes tuples or lists where the first 2 elements are start, end and the remaining
elements can be anything.


```Python
>>> from interlap import InterLap
>>> inter = InterLap()

#Here create 10K random intervals:

>>> import random
>>> random.seed(42)
>>> sites = random.sample(range(22, 100000000, 12), 50000)
>>> ranges = [(i, i + random.randint(2000, 20000)) for i in sites]

>>> inter.update(ranges)
>>> inter.add((20, 22, {'info': 'hi'}))

>>> [20, 21] in inter
True

>>> next(inter.find((21, 21)))
(20, 22, {'info': 'hi'})

>>> inter.add((2, 3, {'info': 'hello'}))

*NOTE*: below shows how edge-cases are handled.
>>> list(inter.find((2, 2)))
[(2, 3, {'info': 'hello'})]
>>> list(inter.find((3, 3)))
[(2, 3, {'info': 'hello'})]

Test every item in the InterLap. These 50K queries take < 0.5 seconds:

>>> for s, e in ranges:
...     assert (s, e) in inter

InterLap objects are iterable:

>>> for seo in inter:
...     assert (seo[0], seo[1]) in inter

```

Installation
------------

```Shell

pip install interlap

```

Example
-------

In general, we will want one interlap per chromosome for genomic data.
The snippet below shows a simple way to do that in the process of creating
and then querying some intervals.

```Python

from interlap import InterLap
from collections import defaultdict
import sys

# use defaultdict to key by chromosome.
inter = defaultdict(InterLap)

for toks in (x.rstrip().split("\t") for x in open(sys.argv[1])):
    start, end = map(int, toks[1:3])
    inter[toks[0]].add((start, end, toks))

# now find overlaps in another file:

for toks in (x.rstrip().split("\t") for x in open(sys.argv[2])):
    start, end = map(int, toks[1:3])

    print list(inter[toks[0]].find((start, end)))

```

Why
---

I am aware of bx-python's interval tree (I wrote the cython version)
but for some projects it is nice to have a simple dependency (or no
dependency since this can be included as a single file or 30 lines
of code) when we just need something a bit better than naive overlap
testing.

In my testing, the method implemented here, using a sorted list and keeping
track of the longest observed interval is the fastest *pure python* method
*as long as the longest observed interval is does not cover a substantial 
fraction of intervals in the set*.


IntervalSet Operations
----------------------

As of version 0.2.0 Interlap also includes an `Interval` class that behaves
like what is normally called an interval set.

```python

# note how it merges overlapping sub-intervals.
>>> Interval([(1, 95), (95, 100)]).add(Interval([(90, 100)]))
Interval([(1, 100)])

# it also has a fairly specialize 'split' function to split an interval-set
# by another set of intervals:
>>> Interval([(1, 50), (60, 80)]).split([(45, 65), (70, 74), (76, 78)])
[Interval([(1, 45)]), Interval([(65, 70)]), Interval([(74, 76)]), Interval([(78, 80)])]

>>> Interval([(45, 65), (70, 74), (76, 78)]).split([(1, 50), (60, 80)])
[Interval([(50, 60)])]

```

See the doctests under the Interval class for more details
