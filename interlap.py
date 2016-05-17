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

>>> for i, se in enumerate(inter):
...     if i > 10: break
...     assert se[0] < se[1]
...

>>> list(inter.closest((2, 2)))
[(2, 3, {'info': 'hello'})]

>>> list(inter.closest((2, 21)))
[(2, 3, {'info': 'hello'}), (20, 22, {'info': 'hi'})]

>>> list(inter.closest((2, 21)))
[(2, 3, {'info': 'hello'}), (20, 22, {'info': 'hi'})]

>>> list(inter.closest((11, 12)))
[(2, 3, {'info': 'hello'}), (20, 22, {'info': 'hi'})]

>>> for i in range(10):
...     inter.add((20, 21))
>>> len(list(inter.closest((10, 13))))
12

>>> for i in range(10):
...     inter.add((18, 21))
>>> len(list(inter.closest((10, 13))))
10

>>> for i in range(10):
...     inter.add((5, 5))
>>> len(list(inter.closest((10, 13))))
20

>>> for i in range(10):
...     inter.add((3, 5))
>>> len(list(inter.closest((10, 13))))
30

>>> inter.add((9, 9))
>>> list(inter.closest((10, 13)))
[(9, 9)]

"""
from itertools import groupby
from operator import itemgetter

__all__ = ['InterLap']

__version__ = '0.2.2'

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
        iopts = iset[l:r]
        iiter = (s for s in iopts if s[0] <= other[1] and s[1] >= other[0])
        for o in iiter: yield o

    def closest(self, other):
        iset = self._iset
        l = binsearch_left_start(iset, other[0] - self._maxlen, 0, len(iset))
        r = binsearch_right_end(iset, other[1], 0, len(iset))
        l, r = max(0, l - 1), min(len(iset), r + 2)

        while r < len(iset) and iset[r - 1][0] == iset[r][0]:
            r += 1

        while l > 1 and iset[l][1] == iset[l + 1][1]:
            l -= 1
        iopts = iset[l:r]
        ovls = [s for s in iopts if s[0] <= other[1] and s[1] >= other[0]]
        if ovls:
            for o in ovls: yield o
        else:
            iopts = sorted([(min(abs(i[0] - other[1]), abs(other[0] - i[1])), i) for i in iopts])
            for dist, g in groupby(iopts, itemgetter(0)):
                # only yield the closest intervals
                for d, ival in g:
                    yield ival
                break

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

    def __iter__(self):
        return iter(self._iset)


def overlaps(s1, e1, s2, e2):
    """
    >>> overlaps(2, 4, 3, 5)
    True
    >>> overlaps(2, 4, 4, 5)
    False
    >>> overlaps(2, 200, 3, 5)
    True
    >>> overlaps(3, 5, 2, 200)
    True
    >>> overlaps (3, 5, 5, 6)
    False
    >>> overlaps (5, 6, 3, 5)
    False
    """
    return not (e1 <= s2 or s1 >= e2)

def reduce(args):
    """
    >>> reduce([(2, 4), (4, 9)])
    [(2, 4), (4, 9)]

    >>> reduce([(2, 6), (4, 10)])
    [(2, 10)]
    """
    if len(args) < 2: return args
    args.sort()
    ret = [args[0]]
    for next_i, (s, e) in enumerate(args, start=1):
        if next_i == len(args):
            ret[-1] = ret[-1][0], max(ret[-1][1], e)
            break

        ns, ne = args[next_i]
        if e > ns or ret[-1][1] > ns:
            ret[-1] = ret[-1][0], max(e, ne, ret[-1][1])
        else:
            ret.append((ns, ne))
    return ret


class Interval(object):
    """
    >>> i = Interval([(2, 10), (8, 20), (30, 40)])
    >>> i
    Interval([(2, 20), (30, 40)])

    >>> i.add([(20, 22)])
    >>> i
    Interval([(2, 20), (20, 22), (30, 40)])

    >>> i.add([(10, 31)])
    >>> i
    Interval([(2, 40)])

    >>> i.add([(55, 65), (65, 75), (75, 85), (85, 95), (95, 100)])
    >>> i
    Interval([(2, 40), (55, 65), (65, 75), (75, 85), (85, 95), (95, 100)])

    >>> i.add([(1, 95)])
    >>> i
    Interval([(1, 95), (95, 100)])

    >>> i.add(Interval([(90, 100)]))
    >>> i
    Interval([(1, 100)])

    >>> Interval([(1, 100)]).split([(55, 65), (75, 85)])
    [Interval([(1, 55)]), Interval([(65, 75)]), Interval([(85, 100)])]

    >>> Interval([(1, 50), (60, 80)]).split([(45, 65), (75, 85)])
    [Interval([(1, 45)]), Interval([(65, 75)])]

    >>> Interval([(1, 50), (60, 80)]).split([(45, 65), (70, 74), (76, 78)])
    [Interval([(1, 45)]), Interval([(65, 70)]), Interval([(74, 76)]), Interval([(78, 80)])]

    >>> Interval([(45, 65), (70, 74), (76, 78)]).split([(1, 50), (60, 80)])
    [Interval([(50, 60)])]

    >>> Interval([(45, 65), (70, 95)]).split([(66, 67)])
    [Interval([(45, 65)]), Interval([(70, 95)])]

    >>> Interval()
    Interval([])

    """

    __slots__ = ('_vals', '_fixed')

    def __init__(self, args=None):
        self._vals = []
        if args is None: return
        assert isinstance(args, list)
        if len(args) > 0:
            assert isinstance(args[0], tuple)
            assert isinstance(args[0][0], (int, long))
            self._vals = reduce(args)

    def _as_tuples(self, args):
        vals = []
        if isinstance(args, Interval):
            vals.extend(args._vals)
        else:
            for a in args:
                if isinstance(a, Interval):
                    vals.extend(a._vals)
                else:
                    vals.append(a)
        return vals

    def add(self, args):
        self._vals = reduce(self._vals + self._as_tuples(args))

    def __repr__(self):
        return "%s(%r)" % (self.__class__.__name__, self._vals)

    def split(self, others):
        others = sorted(self._as_tuples(others))
        old = self._vals

        ret = []
        last = []
        n_greater = 0
        for s in self._vals:
            os = [o for o in others if overlaps(s[0], s[1], o[0], o[1])]
            greater = sum(1 for o in others if s[0] >= o[1])
            inew = False
            if greater != n_greater:
                inew = True
                n_greater = greater
            if os:
                if last:
                    ret.append(Interval(last))
                    last = []
                # split or truncate the current s interval
                # truncate right-end of interval
                start = s[0]
                for i, o in enumerate(os):
                    if s[0] < o[0]:
                        last.append((start, min(s[1], o[0])))
                        #print("start, s, o:", start, s, o, file=sys.stderr)
                        ret.append(Interval(last))
                        last = []
                    if s[1] > o[1]:
                        if last:
                            ret.append(Interval(last))
                            last = []
                        last.append((max(s[0], o[1]), s[1]))
                        if i < len(os) - 1:
                            if os[i + 1][0] < last[-1][1]:
                                last[-1] = last[-1][0], os[i + 1][0]

                    start = o[1]
            else:
                last.append(s)
            # one of the splitters fell between the most recently added interval
            # and what preceded it.
            if inew:
                if len(last) > 1:
                    a, last = last[:-1], last[-1:]
                    ret.append(Interval(a))

        if last:
            ret.append(Interval(last))
        return ret


if __name__ == "__main__":
    import time
    t0 = time.time()
    import doctest
    print(doctest.testmod(verbose=0,
          optionflags=doctest.REPORT_ONLY_FIRST_FAILURE))
    print(time.time() - t0)

