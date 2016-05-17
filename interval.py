from __future__ import print_function

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


    """

    __slots__ = ('_vals', '_fixed')

    def __init__(self, args):
        assert isinstance(args, list)
        assert isinstance(args[0], tuple)
        assert isinstance(args[0][0], (int, long))
        self._vals = []
        self._vals = reduce(args)
        self._fixed = True

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

if __name__ == "__main__":
    import doctest
    print(doctest.testmod())

