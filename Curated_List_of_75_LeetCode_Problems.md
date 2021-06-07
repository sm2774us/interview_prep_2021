# LeetCode Curated List of 75 Problems.

## Math

## Array

## Bit Manipulation

## Dynamic Programming

## Graph

## Interval

## Linked List

## Matrix

## String

## Tree

## Heap

## Recursion

## Sliding Window

## Greedy and Backtracking


## Math
| #     | Title	               | url                                                 | Time     | Space | Difficulty | Tag	      | Note                   |
| ----- | -------------------- | --------------------------------------------------- | -------- | ----- | ---------- | ---------- | ---------------------- |
| 0149  | Max Points on a Line | https://leetcode.com/problems/max-points-on-a-line/ | O(n^2)	| O(n)  | Hard       |            |                        |
# The Math Behind the Solution
```
Q1> What is the equation of a line, in the form:
ax + by + c = 0, 
with gradient -2 through the point (4, -6) ?

A1> 
First, we should know that slope of linear equation is :
m = (y1 − y2)/(x1 − x2)
and we can form the equation by this formula.

In this case, we have gradient (slope) = −2 and the 
point (4, -6).
We can just simply substitute the things we know into 
the above equation.

So, the equation will be:
−2 = (y − (−6)) / (x − 4)
−2(x − 4) = y + 6
−2x + 8 = y + 6

And we can change it in the form ax + by + c = 0, which is
−2x −y + 2 = 0

Q2> Finding a linear equation ax + by + c = 0 given 2 points ?

A2>
If you have 2 points, for example:

P( 2, 1 ); Q( 5, 7 )

You can find the linear equation of the line that passes through those points in the form:

Ax + By + C = 0

in one step by simply using the formula:

(y1 – y2)x + (x2 – x1)y + (x1y2 – x2y1) = 0

OR, we can write this as (easier to read!)
(py – qy)x + (qx – px)y + (px*qy – qx*py) = 0

Let’s try it:

Take Point1=( 2, 1 )
Take Point2=( 5, 7 )

Find the LINEAR EQUATION of the line that passes through the points (2,1) and (5,7). Your answer must be in the form of Ax + By + C = 0.

Using the equation:
(y1 – y2)x + (x2 – x1)y + (x1*y2 – x2*y1) = 0

We’ll just plug numbers in:
(1 – 7)x + (5 – 2)y + ( (2 x 7) – (5 x 1) ) = 0
-6x + 3y + (14 – 5) = 0
-6x + 3y + 9 = 0

Factoring a -3 out:
-3( 2x – y – 3 ) = 0

Dividing both sides by -3:
2x – y – 3 = 0

And that’s the answer.
```
```
Solution Explanation
-----------------------
1) Remember from Euclidean Geometry that the general equation of a line is represented as: 
ax + by + c = 0
or
ax + by = c

2) if x1 == x2:
equation: x = x1
       o: a, b, c = 1, 0, -x1

3) Not using usual linear equation y = kx + b
because of the precision problem of float.
We need to make sure all numbers are integer
so use ax + by + c = 0 instead
values of a, b, c can be deduced from the following simultaneous equations:
3.1) ax1 + by1 + c = 0
3.2) ax2 + by2 + c = 0
Using basic algebra, we get:
3.3) a/b = (y1-y2)/(x2-x1), set a to y1-y2, b to x2-x1
Substitute a and b into the previous simultaneous equations
we get:
3.4) c = x1y2-x2y1

4)
To make sure the same (a, b, c) are calculated for all
same lines. We need to set some rules. Because values of
a and b are only bounded by a ratio (eq 2), and a and b
are picked arbitrarily.

  first: make sure a is positive
 reason: -ax - by - c = 0 represents the same line
     as: ax + by + c = 0

 second: reduce fraction of a, b and c
 reason: 2ax + 2by + 2c = 0 represents the same line
     as: ax + by + c = 0				

a, b, c = y1 - y2, x2 - x1, x1 * y2 - x2 * y1

5) get gcd of a, b and c

6) reduce fraction

7) handle edge case: return 0 for 0 point, 1 for 1 point

8) get max on the length of values (x and y coordinates of each point that constitutes a line) in the result dictionary -> pointsInLine 
   and that is your answer
```
```
#Complexity Analysis
# ---------------------
TC: O(N^2) ... Test all pairs: O(N^2)
SC: O(N)   ... Storing the result in a dictionary <key=<int, int>, value=count>.
               Considering all points are distinct we will have n entries in dictionary.
```
```python
# Solution 1>
# ---------------------
# Prefer this over Solution-2 in an interview setting.
from typing import List

class Solution:
    def gcd(self, a: int, b: int) -> int:
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
```
```python
# Solution 2>
# ---------------------
# Easier solution using Python built-ins:
# 1) itertools.combinations
# 2) collections.defaultdict
# 3) collections.Counter
#		- most_common
import collections
from itertools import combinations
from math import gcd
from typing import List

class Solution:
    def maxPoints(self, points: List[List[int]]) -> int:
        # Represent line uniquely as (a,b,c) where ax+by+c=0 and
        # a>0 or (a==0 and b>0) or (a==0 and b==0).
        if not points:
            return 0

        points = [(p[0], p[1]) for p in points]

        counter, points, lines = collections.Counter(points), set(points), collections.defaultdict(set)

        for p1, p2 in combinations(points, 2):
            (x1, y1), (x2, y2) = p1, p2
            a, b, c = y1 - y2, x2 - x1, x1 * y2 - x2 * y1
            gcd_ = gcd(gcd(abs(a), abs(b)), abs(c))
            lines[(a // gcd_, b // gcd_, c // gcd_)] |= {p1, p2}

        counter_most_common = counter.most_common(1)[0][1]
        if not lines:
            return counter_most_common
        return max(counter_most_common, max(
            sum(counter[p] for p in ps)
            for ps in lines.values()
        ))
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

## String
| #     | Title	                                         | url                                                                           | Time   | Space | Difficulty | Tag	      | Note                   |
| ----- | ---------------------------------------------- | ----------------------------------------------------------------------------- | ------ | ----- | ---------- | ------------ | ---------------------- |
| 0149  | Longest Substring Without Repeating Characters | https://leetcode.com/problems/longest-substring-without-repeating-characters/ | O(n)	  | O(1)  | Medium     |              |                        |
| 0424  | Longest Repeating Character Replacement        | https://leetcode.com/problems/longest-repeating-character-replacement/        | O(n)	  | O(1)  | Medium     |              | Sliding Window         |
| 0076  | Minimum Window Substring                       | https://leetcode.com/problems/minimum-window-substring/                       | O(n)   | O(k)  | Hard       |              |                        |
| 0242  | Valid Anagram                                  | https://leetcode.com/problems/valid-anagram/                                  | O(n)   | O(1)  | Easy       |              |                        |
| 0049  | Group Anagrams                                 | https://leetcode.com/problems/group-anagrams/                                 | O(n)   | O(1)  | Easy       |              |                        |

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>
