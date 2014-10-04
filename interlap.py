"""
`InterLap` does fast interval overlap testing with a simple python data
structure.

It works well on the types of querying done in genomic datasets where we have
10's of thousands of intervals and we check for overlap millions of times. It
is very simple and has no dependencies.

It takes tuples or lists where the first 2 elements are start, end and the
remaining elements can be anything.

>>> from interlap import InterLap
>>> inter = InterLap()

Here create 10K random intervals:

>>> import random
>>> sites = random.sample(range(22, 100000000, 12), 50000)
>>> ranges = [(i, i + random.randint(2000, 20000)) for i in sites]

Add them to the interval tree (this takes < 0.5 seconds):

>>> inter.update(ranges)

We can also add one at a time:

>>> inter.add((20, 22, {'info': 'hi'}))

Now do overlap testing:

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
"""

__all__ = ['InterLap']

__version__ = '0.1.1'

try:
    int_types = (int, long)
except NameError:
    int_types = (int,)

######################################
def binsearch_left_start(intervals, x, lo, hi):
    while lo < hi:
        mid = (lo + hi)//2
        f = intervals[mid]
        if f[0] < x: lo = mid + 1
        else: hi = mid
    return lo

# like python's bisect_right find the _highest_ index where the value x 
# could be inserted to maintain order in the list intervals
def binsearch_right_end(intervals, x, lo, hi):
    while lo < hi:
        mid = (lo + hi)//2
        f = intervals[mid]
        if x < f[0]: hi = mid
        else: lo = mid + 1
    return lo

class InterLap(object):

    """Create an Interlap object. (See module docstring)."""

    def __init__(self, ranges=()):
        """The ranges are tuples of (start, end, *whatever)."""
        self._iset = sorted(ranges)
        self._maxlen = max(r[1] - r[0] for r in (ranges or [[0, 0]]))

    def add(self, ranges):
        r"""Add a single (or many) [start, end, \*] item to the tree."""
        if len(ranges) and isinstance(ranges[0], int_types):
            ranges = [ranges]
        iset = self._iset
        self._maxlen = max(self._maxlen, max(r[1] - r[0] + 1 for r in ranges))

        if len(ranges) > 30 or len(iset) < len(ranges):
            iset.extend(ranges)
            iset.sort()
        else:
            for o in ranges:
                iset.insert(binsearch_left_start(iset, o[0], 0, len(iset)), o)

    update = add

    def __len__(self):
        """Return number of intervals."""
        return len(self._iset)

    def find(self, other):
        """Return an interable of elements that overlap other in the tree."""
        iset = self._iset
        l = binsearch_left_start(iset, other[0] - self._maxlen, 0, len(iset))
        r = binsearch_right_end(iset, other[1], 0, len(iset))
        iopts = self._iset[l:r]
        iiter = (s for s in iopts if s[0] <= other[1] and s[1] >= other[0])
        for o in iiter: yield o

    def __contains__(self, other):
        """Indicate whether `other` overlaps any elements in the tree."""
        iset = self._iset
        l = binsearch_left_start(iset, other[0] - self._maxlen, 0, len(iset))
        # since often the found interval will overlap, we short cut that
        # case of speed.
        max_search = 8
        if l == len(iset): return False
        for left in iset[l:l + max_search]:
            if left[1] >= other[0] and left[0] <= other[1]: return True
            if left[0] > other[1]: return False

        r = binsearch_right_end(iset, other[1], 0, len(iset))
        return any(s[0] <= other[1] and s[1] >= other[0]
                   for s in iset[l + max_search:r])

if __name__ == "__main__":
    import time
    t0 = time.time()
    import doctest
    print(doctest.testmod(verbose=0,
          optionflags=doctest.REPORT_ONLY_FIRST_FAILURE))
    print(time.time() - t0)

