# Time:  O(n^2)
# Space: O(n)

class Solution:
    def gcd(self, a, b):
        return gcd(b % a, a) if a != 0 else b
			
    def maxPoints(self, points: List[List[int]]) -> int:
        if not points:
            return 0

        pointsInLine = {}
		
        for i in range(len(points)):
            for j in range(i+1, len(points)):
                x1, y1 = points[i]
                x2, y2 = points[j]
                # line: ax + by + c = 0
                if x1 == x2:
                    # equation: x = x1
                    a, b, c = 1, 0, -x1
                else:
                    # not using usual linear equation y = kx + b
                    # because of the precision problem of float
                    # we need to make sure all numbers are integer
                    # so use ax + by + c = 0 instead
                    # value of a, b, c can be deduced from the following
                    # simultaneous equations:
                    # 1) ax1 + by1 + c = 0
                    #    ax2 + by2 + c = 0
                    # using basic algebra, we get:
                    # 2) a/b = -(y2-y1)/(x2-x1), set a to y2-y1, b to x1-x2
                    # substitute a and b into the previous simultaneous equations
                    # we get: 3) c = x2y1-x1y2
                    # To make sure the same (a, b, c) are calculated for all
                    # same lines. We need to set some rules. Because values of
                    # a and b are only bounded by a ratio (eq 2), and a and b
                    # are picked arbitrarily.
                    # first: make sure a is positive
                    # 		reason: -ax - by - c = 0 represents the same line
                    # 				as ax + by + c = 0
                    # second: reduce fraction of a, b and c
                    # 		reason: 2ax + 2by + 2c = 0 represents the same line
                    # 				as ax + by + c = 0				
                    a, b, c = y1 - y2, x2 - x1, x1 * y2 - x2 * y1
                    # get gcd of a, b and c
                    rate = self.gcd(self.gcd(abs(a), abs(b)), abs(c))
                    #rate = reduce(self.gcd, (abs(a), abs(b), abs(c)))
                    # reduce fraction
                    a, b, c = a / rate, b / rate, c / rate
            line = (a, b, c)
            pointsInLine.setdefault(line, set())
            pointsInLine[line].add(i)
            pointsInLine[line].add(j)
        # edge case: return 0 for 0 point, 1 for 1 point
        return max(map(len, pointsInLine.values())) if pointsInLine else len(points)
