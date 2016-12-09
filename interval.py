from __future__ import print_function
import sys

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

    >>> Interval([(11, 18), (22, 28), (32, 39), (42, 48)]).split([(12, 18), (22, 26),(43, 44),(43, 48)])
    [Interval([(11, 12)]), Interval([(26, 28)]), Interval([(32, 39)]), Interval([(42, 43)])]

    >>> Interval([(11, 18), (22, 28), (32, 39), (42, 48)]).split([(12, 18), (22, 26),(43, 44),(44, 48)])
    [Interval([(11, 12)]), Interval([(26, 28)]), Interval([(32, 39)]), Interval([(42, 43)])]
    
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

    """

    __slots__ = ('_vals', '_fixed')

    def __init__(self, args):
        assert isinstance(args, list)
        self._vals = []
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
                        if min(s[1], o[0])>start:
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
                            if os[i + 1][0] < last[-1][1] and last[-1][0] < os[i + 1][0]:
                                last[-1] = last[-1][0], os[i + 1][0]
                            elif last[-1][0] >= os[i + 1][0]:
                                last.pop()

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
