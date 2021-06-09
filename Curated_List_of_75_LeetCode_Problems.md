# LeetCode Curated List of 75 Problems.

## Algorithms

* [Math](#math)
* [Array](#array)
* [Bit Manipulation](#bit-manipulation)
* [Dynamic Programming](#dynamic-programming)
* [Graph](#graph)
* [Interval](#interval)
* [Linked List](#linked-list)
* [Matrix](#matrix)
* [String](#string)
* [Tree](#tree)
* [Heap](#heap)
* [Recursion](#recursion)
* [Sliding Window](#sliding-window)
* [Greedy and Backtracking](#greedy-and-backtracking)

## Math
| #     | Title	               | url                                                 | Time                  | Space        | Difficulty | Tag	                               | Note                   |
| ----- | -------------------- | --------------------------------------------------- | --------------------- | ------------ | ---------- | ----------------------------------- | ---------------------- |
| 0012  | [Integer to Roman](#lc-12integer-to-roman) | https://leetcode.com/problems/integer-to-roman/     | _O(n)_                | _O(1)_       | Medium     |                                     |                        |
| 0013  | [Roman to Integer](#lc-13roman-to-integer) | https://leetcode.com/problems/roman-to-integer/     | _O(n)_                | _O(1)_       | Easy       |                                     |                        |
| 0149  | [Max Points on a Line](#lc-149max-points-on-a-line) | https://leetcode.com/problems/max-points-on-a-line/ | _O(n^2)_	             | _O(n)_       | Hard       |                                     | Linear Equation `ax + by + c = 0` |
| 0204  | [Count Primes](#lc-204count-primes) | https://leetcode.com/problems/count-primes/         | _O( N*Log(Log(N)) )_  | _O(N)_       | Easy       |                                     | Sieve of Eratosthenes  |
| 0372  | [Super Pow](#lc-372super-pow) | https://leetcode.com/problems/super-pow/            | _O(n)_                | _O(1)_       | Medium     |                                     |                        |
| 0509  | [Fibonacci Number](#lc-509fibonacci-number) | https://leetcode.com/problems/fibonacci-number/     | _O(logn)_             | _O(1)_       | Easy       | variant of [Climbing Stairs](https://leetcode.com/problems/climbing-stairs/) | Matrix Exponentiation, Binet's Formula |
| 1390  | [Four Divisors](#solution-1--sieve-of-eratosthenes-) | https://leetcode.com/problems/four-divisors/        | _O(N+K*Log(Log(K)))_, where, N = max(nums), M=len(nums), K is len(primes) | _O(N+M+K^2)_, where, N = max(nums), M=len(nums), K is len(primes) | Medium     |                                     | Sieve of Eratosthenes  |
| 1390  | [Four Divisors](#solution-2--recursion-) | https://leetcode.com/problems/four-divisors/        | _O(N * sqrt(M))_, where, N = length of nums and M = nums[i]  | _O(1)_       | Medium     |                                     | Recursion              |
#### [LC-12:Integer to Roman](https://leetcode.com/problems/integer-to-roman/)
##### Solution Explanation
```
Idea:

Just like Roman to Integer, this problem is most easily solved using a lookup table (dictionary) for the conversion between digit and numeral. 
In this case, we can easily deal with the values in descending order and insert the appropriate numeral (or numerals) 
as many times as we can while reducing the our target number (N) by the same amount.

Once N runs out, we can return ans.
```
##### Complexity Analysis
```
TC: O(13) = O(1), iterate the length of dictionary keys
SC: O(13) = O(1), one hash map (dictionary)
```
```python
class Solution(object):
    def __init__(self):
        self.value_map = {
			1000: 'M',
			900: 'CM',
			500: 'D',
			400: 'CD',
			100: 'C',
			90: 'XC',
			50: 'L',
			40: 'XL',
			10: 'X',
			9: 'IX',
            5: 'V',
			4: 'IV',
			1: 'I'
        }

    def intToRoman(self, num):
        """
        :type num: int
        :rtype: str
        """

        #d = {1000: 'M', 900: 'CM', 500: 'D', 400: 'CD', 100: 'C', 90: 'XC', 50: 'L', 40: 'XL', 10: 'X', 9: 'IX',
        #     5: 'V', 4: 'IV', 1: 'I'}
		##d = {1: 'I', 4: 'IV', 5: 'V', 9: 'IX', 10: 'X', 40: 'XL', 50: 'L', 90: 'XC', 100: 'C', 400: 'CD',
		##     500: 'D', 900: 'CM', 1000: 'M'}
        roman = ""

		for v in sorted(self.value_map.keys(), reverse=True):
        #for v in sorted(d.keys(), reverse=True):
		##for v in d.keys():
            roman += (num // v) * self.value_map[v]
            num -= (num // v) * v

        return roman
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-13:Roman to Integer](https://leetcode.com/problems/roman-to-integer/)
##### Solution Explanation
```
Approach:

Use dictionary for fast and easy lookup of numeral to integer value.
Go through each numeral in the input string
If numeral is smaller than the next numeral in the input we have a value like IV so subtract the current numeral from the next numeral.
Else add the value of the numeral to our result.
```
##### Complexity Analysis
```
TC: O(N)
SC: O(1)
```
```python
class Solution:
    
    def __init__(self):
        self.value_map = {
            'I': 1,
            'V': 5,
            'X': 10,
            'L': 50,
            'C': 100,
            'D': 500,
            'M': 1000
        }
    
    def romanToInt(self, s: str) -> int:
        result = 0
        index = 0
        length = len(s)
        
        while index < length:
            current_num = self.value_map[s[index]]
            
            # Check next value to see if it is larger. If
            # it is that means that the we are dealing with
            # something like 4 'IV' and see need to subtract
            # the current numeral from next
            if index+1 < length:
                next_num = self.value_map[s[index+1]]
                if next_num > current_num:
                    current_num = next_num - current_num
                    # Skip ahead an additional index since we combined two numerals
                    index += 1
            
            result += current_num
            index += 1
            
        return result
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-149:Max Points on a Line](https://leetcode.com/problems/max-points-on-a-line/)
##### The Math Behind the Solution
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
##### Solution Explanation
```
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
##### Complexity Analysis
```
TC: O(N^2) ... Test all pairs: O(N^2)
SC: O(N)   ... Storing the result in a dictionary <key=<int, int>, value=count>.
               Considering all points are distinct we will have n entries in dictionary.
```
##### Solution-1 ( Prefer this solution in an interview setting ).
```python
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
##### Solution-2 ( Uses Python built-ins `math.gcd`, `itertools.combinations` and `collections.Counter.most_common` and data structures `collections.defaultdict` and `collections.Counter` ).
```python
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

#### [LC-204:Count Primes](https://leetcode.com/problems/count-primes/)
##### The Math Behind the Solution
```
Reference: https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes

In mathematics, the sieve of Eratosthenes is an ancient algorithm for finding all prime numbers up to any given limit.

A prime number is a natural number that has exactly two distinct natural number divisors: the number 1 and itself.

To find all the prime numbers less than or equal to a given integer n by Eratosthenes' method:

  1. Create a list of consecutive integers from 2 through n: (2, 3, 4, ..., n).
  2. Initially, let p equal 2, the smallest prime number.
  3. Enumerate the multiples of p by counting in increments of p from 2p to n, and mark them in the list (these will be 2p, 3p, 4p, ...; the p itself should not be marked).
  4. Find the smallest number in the list greater than p that is not marked. If there was no such number, stop. Otherwise, let p now equal this new number (which is the next prime), and repeat from step 3.
  5. When the algorithm terminates, the numbers remaining not marked in the list are all the primes below n.
  6. The main idea here is that every value given to p will be prime, because if it were composite it would be marked as a multiple of some other, smaller prime. Note that some of the numbers may be marked more than once (e.g., 15 will be marked both for 3 and 5).

As a refinement, it is sufficient to mark the numbers in step 3 starting from p2, as all the smaller multiples of p will have already been marked at that point. This means that the algorithm is allowed to terminate in step 4 when p2 is greater than n.[1]

Another refinement is to initially list odd numbers only, (3, 5, ..., n), and count in increments of 2p from p2 in step 3, thus marking only odd multiples of p. This actually appears in the original algorithm.[1] This can be generalized with wheel factorization, forming the initial list only from numbers coprime with the first few primes and not just from odds (i.e., numbers coprime with 2), and counting in the correspondingly adjusted increments so that only such multiples of p are generated that are coprime with those small primes, in the first place.[6]
-------------------------
Example
To find all the prime numbers less than or equal to 30, proceed as follows.

First, generate a list of integers from 2 to 30:
```
> 2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
```
The first number in the list is 2; cross out every 2nd number in the list after 2 by counting up from 2 in increments of 2 (these will be all the multiples of 2 in the list):
```
> 2  3  ~~4~~  5  ~~6~~  7  ~~8~~  9  ~~10~~ 11 ~~12~~ 13 ~~14~~ 15 ~~16~~ 17 ~~18~~ 19 ~~20~~ 21 ~~22~~ 23 ~~24~~ 25 ~~26~~ 27 ~~28~~ 29 ~~30~~
```
The next number in the list after 2 is 3; cross out every 3rd number in the list after 3 by counting up from 3 in increments of 3 (these will be all the multiples of 3 in the list):
```
> 2  3  ~~4~~  5  ~~6~~  7  ~~8  9  10~~ 11 ~~12~~ 13 ~~14 15 16~~ 17 ~~18~~ 19 ~~20 21 22~~ 23 ~~24~~ 25 ~~26 27 28~~ 29 ~~30~~
```
The next number not yet crossed out in the list after 3 is 5; cross out every 5th number in the list after 5 by counting up from 5 in increments of 5 (i.e. all the multiples of 5):
```
> 2  3  ~~4~~  5  ~~6~~  7  ~~8  9  10~~ 11 ~~12~~ 13 ~~14 15 16~~ 17 ~~18~~ 19 ~~20 21 22~~ 23 ~~24 25 26 27 28~~ 29 ~~30~~
```
The next number not yet crossed out in the list after 5 is 7; the next step would be to cross out every 7th number in the list after 7, but they are all already crossed out at this point, as these numbers (14, 21, 28) are also multiples of smaller primes because 7 × 7 is greater than 30. The numbers not crossed out at this point in the list are all the prime numbers below 30:
```
> 2  3     5     7           11    13          17    19          23                29
```
-------------------------
Pseudocode
-------------------------
The sieve of Eratosthenes can be expressed in pseudocode, as follows:
```
> algorithm Sieve of Eratosthenes is
>     input: an integer n > 1.
>     output: all prime numbers from 2 through n.
> 
>     let A be an array of Boolean values, indexed by integers 2 to n,
>     initially all set to true.
>     
>     for i = 2, 3, 4, ..., not exceeding √n do
>         if A[i] is true
>             for j = i2, i2+i, i2+2i, i2+3i, ..., not exceeding n do
>                 A[j] := false
> 
>     return all i such that A[i] is true.
> 
```
This algorithm produces all primes not greater than n. It includes a common optimization, which is to start enumerating t
he multiples of each prime i from i^2. The time complexity of this algorithm is O( N* Log(Log(N)) ),
provided the array update is an O(1) operation, as is usually the case.
```

##### Solution Explanation
```
Algorithm:
-------------------------------
  1. Create a list of consecutive integers from 2 through n: (2, 3, 4, ..., n).
  2. Initially, let p equal 2, the smallest prime number.
  3. Enumerate the multiples of p by counting in increments of p from 2p to n, and mark them in the list (these will be 2p, 3p, 4p, ...; the p itself should not be marked).
  4. Find the smallest number in the list greater than p that is not marked. If there was no such number, stop. Otherwise, let p now equal this new number (which is the next prime), and repeat from step 3.
  5. When the algorithm terminates, the numbers remaining not marked in the list are all the primes below n.
  6. The main idea here is that every value given to p will be prime, because if it were composite it would be marked as a multiple of some other, smaller prime. Note that some of the numbers may be marked more than once (e.g., 15 will be marked both for 3 and 5).

Make your code faster:
-------------------------------
  * The code line:
    
	lst[m * m: n: m] = [0] *((n-m*m-1)//m + 1) 
	
	is key to reduce the run time.
	You could write a loop like this (but it would be very expensive):
	
    for i in range(m * m, n,  m):
        lst[i] = 0

  * After marking all the even indices in the first iteration, I do not check even numbers again, 
    and will only check odd numbers in the remaining iterations.

  * I created a list with numeral elements, instead of boolean elements.
  
  * Do not use function sqrt, because it is expensive [do not use: m < sqrt(n)].
    Instead, use: m * m < n.
```
##### Complexity Analysis
```
TC: O( N* Log(Log(N)) )

The time complexity is O(n/2 + n/3 + n/5 + n/7 + n/11 + ....) which is equivalent to O( N* Log(Log(N)) ).

SC: O(N)
```

```python
class Solution:
    def countPrimes(self, n: int) -> int:
        if n < 3: return 0     ###// No prime number less than 2
        lst = [1] * n          ###// create a list for marking numbers less than n
        lst[0] = lst[1] = 0    ###// 0 and 1 are not prime numbers
        m = 2
        while m * m < n:       ###// we only check a number (m) if its square is less than n
            if lst[m] == 1:    ###// if m is already marked by 0, no need to check its multiples.
			
			    ###// If m is marked by 1, we mark all its multiples from m * m to n by 0. 
			    ###// 1 + (n - m * m - 1) // m is equal to the number of multiples of m from m * m to n
                lst[m * m: n: m] = [0] *(1 + (n - m * m - 1) // m)
				
			###// If it is the first iteration (e.g. m = 2), add 1 to m (e.g. m = m + 1; 
			### // which means m will be 3 in the next iteration), 
            ###// otherwise: (m = m + 2); This way we avoid checking even numbers again.	
            m += 1 if m == 2 else 2
        return sum(lst)
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-372:Super Pow](https://leetcode.com/problems/super-pow/)
##### Solution-1 ( Euler's theorem ).
##### The Math Behind the Solution
> ![equation](https://latex.codecogs.com/png.image?\dpi{150}%20a^{\Phi(n)}\equiv1)  `mod n`
> where, ![equation](https://latex.codecogs.com/png.image?\dpi{150}%20\Phi(n)), is Euler's totient function.
>
```
Let exp denote the exponent extracted from input b

Goal = (a ^ exp) mod 1337
= a ^ ( exp mod φ(1337) ) mod 1337
= a ^ ( exp mod 1140) mod 1337
use the formula derived above to reduce computation loading and to have higher performance.

Remark:

φ( 1337 )
= φ( 7 x 191 ) where 7 and 191 are prime factor of 1337's factor decomposition
= φ( 7 ) x φ( 191 )
= ( 7 - 1 ) x ( 191 - 1 )
= 6 x 190
=1140

[Euler's totient function φ( n )](https://en.wikipedia.org/wiki/Euler's_totient_function) counts the positive integers up to a given integer n that are relatively prime to n
```
##### References
```
= Euler's Theorem/Fermat's Little Theorem
  - Euler's theorem, [[https://en.wikipedia.org/wiki/Euler's_theorem](https://en.wikipedia.org/wiki/Euler's_theorem](https://en.wikipedia.org/wiki/Euler's_theorem](https://en.wikipedia.org/wiki/Euler's_theorem)\)
  - Fermat's little theorem, [[https://en.wikipedia.org/wiki/Fermat's_little_theorem](https://en.wikipedia.org/wiki/Fermat's_little_theorem](https://en.wikipedia.org/wiki/Fermat's_little_theorem](https://en.wikipedia.org/wiki/Fermat's_little_theorem)\)
  - Euler's totient function, [[https://en.wikipedia.org/wiki/Euler's_totient_function](https://en.wikipedia.org/wiki/Euler's_totient_function](https://en.wikipedia.org/wiki/Euler's_totient_function](https://en.wikipedia.org/wiki/Euler's_totient_function)\)
  - Examples, http://mathonline.wikidot.com/examples-using-euler-s-theorem
  - Three cases for reducing the power
    + a is multiple of 1337, result is 0
    + a is coprime of 1337, then a ^ b % 1337 = a ^ (b % phi(1337)) % 1337
    + a has divisor of 7 or 191,
  - phi(1337) = phi(7 * 191) = 6 * 190 = 1140
  - a ^ b mod 1337 = a ^ x mod 7, where x = b mod 1140
```
##### Solution Explanation
```
Hint & Reference:

1. First step is to extract exponent from input parameter b, either by string operation, or reduce with lambda.

2. Second step is to compute power with speed up by Euler theorem.
```
##### Complexity Analysis
```
Time complexity O(n)
Space complexity O(1)
```
```python
#Implementation_#1
#Extract exponent by string operation
class Solution:
    def superPow(self, a: int, b: List[int]) -> int:
        
        
        if a % 1337 == 0:
            return 0
        
        else:
            # Euler theorem:
            #
            # ( a ^ phi(n) ) mod n = 1
            #  
            # => ( a ^ b ) mod n = ( a ^ ( b mod phi(n) ) ) mod n
            
            exponent = int( ''.join( map( str, b) ) )
			
			return pow( base = a, exp = exponent % 1140, mod = 1337 )
			

#Implementation_#2
#Extract exponent by reduce with lambda
from functools import reduce

class Solution:
    def superPow(self, a: int, b: List[int]) -> int:
        
        
        if a % 1337 == 0:
            return 0
        
        else:
            # Euler theorem:
            #
            # ( a ^ phi(n) ) mod n = 1
            #  
            # => ( a ^ b ) mod n = ( a ^ ( b mod phi(n) ) ) mod n
            
            decimal = lambda x,y : 10*x+y
            exponent = reduce( decimal, b )
			
            return pow( base = a, exp = exponent % 1140, mod = 1337 )
```
##### Solution-2 ( Modular Exponentiation ).
##### Solution Explanation
```
i)   It's useful in the field of public-key cryptography, https://en.wikipedia.org/wiki/Modular_exponentiation
ii)  c mod m = (a * b) mod m = [(a mod m) * (b mod m)] mod m
iii) Starting from the front, we do the powMod the base with the corresponding digit
iv)  Then, we need to remember to powMod the previous result with 10

pow_mod is quite similar to ordinary pow except computing remainder instead of product. 
Since b could be larger than 2^31 so superPow is added.
Otherwise pow_mod(a, c) is enough with O(lgN), 
where N is the value of c.
Inside superPow, I only call pow_mod(a, 10) with O(1). So the overall time complexity is O(N), 
where N is the length of b.
```
##### Complexity Analysis
```
TC: O(N)
SC: O(1)
```
```python
class Solution:
    def superPow(self, a, b):
		k = 1337
        if not b:
            return 1
        
        a = a % k
        ans = pow(a, b[len(b) - 1]) % k
        for i in range(len(b)-2, -1, -1):
            a = self.pow_mod(a, 10)
            ans = ans * pow(a, b[i]) % k
        return ans
    
    def pow_mod(self, a, b):
		k = 1337
        ans = 1
        while b > 0:
            if b & 1 == 1:
                ans = ans * a % k
            a = a * a % k
            b = b >> 1
        return ans
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-509:Fibonacci Number](https://leetcode.com/problems/fibonacci-number/)
##### Solution-1 ( Using Binet's Formula or the Golden Ratio  ).

##### The Math Behind the Solution
* Calculating terms of the Fibonacci sequence can be tedious when using the recursive formula, especially when finding terms with a large n. 
* Luckily, a mathematician named Leonhard Euler discovered a formula for calculating any Fibonacci number.
* This formula was lost for about 100 years and was rediscovered by another mathematician named Jacques Binet.
* The original formula, known as Binet’s formula, is shown below :
* **Binet’s Formula**: The nth Fibonacci number is given by the following formula:
![equation](https://latex.codecogs.com/png.image?\dpi{150}%20f_{n}=\frac{\left[\left(\frac{1+\sqrt{5}}{2}\right)^{n}-\left(\frac{1-\sqrt{5}}{2}\right)^{n}\right]}{\sqrt{5}})

* Another interesting fact arises when looking at the ratios of consecutive Fibonacci numbers.
* It appears that these ratios are approaching a number.
* The number that these ratios are getting closer to is a special number called the Golden Ratio which is denoted by  (the Greek letter phi). You have seen this number in Binet’s formula.
* The Golden Ratio:
![equation](https://latex.codecogs.com/png.image?\dpi{150}%20\frac{1+\sqrt{5}}{2})
* The Golden Ratio has the decimal approximation of ϕ =1.6180339887.

##### Complexity Analysis
```
TC: O(1)
SC: O(1)
```

```python
# Reference: https://math.libretexts.org/Bookshelves/Applied_Mathematics/Book%3A_College_Mathematics_for_Everyday_Life_(Inigo_et_al)/10%3A_Geometric_Symmetry_and_the_Golden_Ratio/10.04%3A_Fibonacci_Numbers_and_the_Golden_Ratio#:~:text=The%20number%20that%20these%20ratios,this%20number%20in%20Binet%27s%20formula.&text=The%20Golden%20Ratio%20has%20the,for%20a%20variety%20of%20reasons.
class Solution:
    def fib(self, N: int) -> int:
        phi = round((1 + 5 ** 0.5) / 2, 6)
        return round((phi ** N - (-phi) ** (-N)) / (5 ** 0.5))
```

##### Solution-2 ( Using Tail Recursion  ).

##### Complexity Analysis
```
TC: O(N)
SC: O(1)
```

```python
class Solution:
    def fib(self, n: int) -> int:
        def fib_helper(n: int, a: int =0, b: int =1) -> int:
            # tail recursion
            # ie. input = 5
            # 4, 1, 1
            # 3, 1, 2
            # 2, 2, 3
            # 1, 3, 5
            if n == 0:
                return a
            elif n == 1:
                return b
            else:
                return fib_helper(n-1, b, a+b)
        
        return fib_helper(n)
```

##### Solution-3 ( Using Matrix Exponentiation - DP ).

##### Solution Explanation
```
Usually if any problem is solvable in constant space using iterative dp,
then we can apply matrix exponentiation and convert an O(n) problem to O(logn)

We will try to create a matrix out of the recursive relations we know:

Relation 1:      fib(n)   =   1*fib(n-1) +  1*fib(n-2)
Relation 2:      fib(n-1) =   1*fib(n-1) +  0*fib(n-2)

Matrix M:  [1,1]
           [1,0]
				  
|fib (n) | =   |1 1| |fib(n-1)|
|fib(n-1)|     |1 0| |fib(n-2)|

or 

|fib (n) |  =  |1 1| |1 1||fib(n-2)|
|fib(n-1)|     |1 0| |1 0||fib(n-3)|

or

|fib (n) |  =  |1 1| |1 1| |1 1||fib(n-3)|
|fib(n-1)|     |1 0| |1 0| |1 0||fib(n-4)|

or

|fib (n) |  =  |1 1| ^ (n-2) |fib(2)|
|fib(n-1)|     |1 0|         |fib(1)|
---------------------------------------------------------------------
Intuition
---------------------------------------------------------------------
The core idea behind this is to evaluate the equation F[n] = F[n-1] + F[n-2] 
and represent this as somehow a power of a matrix. 
We know that a^n can be calculated in O(log n) time using binary exponentiation.
If we can somehow represent the above recurrence relation as a power of matrix and the base values (F[1], F[0]),
then we can get the F[n] in O(log n) time.
---------------------------------------------------------------------
Algorithm
---------------------------------------------------------------------
So consider this :- [F[n], F[n-1]] = [[1,1], [1,0]] * [F[n-1], F[n-2]].
So what I have done is just make the matrix on the right side a 2 X 1 matrix 
so that it can be represented in the terms of this matrix [[1,1],[1,0]]. 
This is done because if you see the left side, then you can put 
[F[n-1], F[n-2]] = [[1,1],[1,0]] * [F[n-2],F[n-3]]. 

So you if you keep doing this then the above recurrence relation boils down to
[F[n], F[n-1]] = [[1,1], [1,0]] ^ (n-1) * [F[1], F[0]].

So what we have to do is calculate the (n-1)th power of the matrix [[1,1],[1,0]]. 
So we can calculate that just like we calculate the a^n in O(log n) time using binary exponentiation.
Just that a would be the matrix [[1,1],[1,0]] instead of an integer.
Following is the code for this :-
```

##### Complexity Analysis
```
TC: O(log(N))
SC: O(1)
```

```python
class Solution:
    def fib(self, n: int) -> int:
        
        def get_mat_mult(mat, other_mat):
            res = [[0 for _ in range(len(mat[0]))] for _ in range(len(mat))]
            for i in range(len(mat)):
                for j in range(len(mat[i])):
                    for k in range(len(other_mat[i])):
                        res[i][j] += mat[i][k] * other_mat[k][j]
            return res
        
        if n == 0 or n == 1:
            return n
        
        final_mat = [[1,0],[0,1]]
        start_mat = [[1,1], [1,0]]
        n -= 1
        while(n):
            if (n & 1):
                final_mat = get_mat_mult(start_mat, final_mat)
            start_mat = get_mat_mult(start_mat, start_mat) 
            n >>= 1
        return final_mat[0][0]
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-1390:Four Divisors](https://leetcode.com/problems/four-divisors/)

##### Solution-1 ( Sieve of Eratosthenes ).
##### Solution Explanation
```
Based on the property of prime number to count the divisors:

n = p1 ^ (a1) * p2 ^ (a2) *... * pn ^ (an) ; (with p1,p2,..,pn is the prime)

= > the total divisors of number S = (a1+1) * (a2 + 1) * .. *(an +1)

***So if one number have 4 divisors we have 2 cases:

case 1: S = (a1+1) = (3 + 1) = 4 -> N = p1^3

case 2: S = (a1+1) * (a2 + 1) = 2 * 2 = 4 -> N = p1 ^ 1 * p2 ^ 1 ; (with p1,p2 is prime and p2!=p1 )
```
##### Complexity Analysis
```
Time  : O(N+K*log(log(K)))  where N = max(nums), K is len(primes) 
Space : O(N+M+K^2)          where N = max(nums), M=len(nums), K is len(primes) 
```
```python
class Solution:
	def sieve(self, n):
		mark = [True] * (n + 1)
		mark[0] = False
		mark[1] = False
		for i in range(2, int(n **(0.5)) + 2):
			for j in range(i*i, n + 1, i):
				mark[j]  = False
		return mark

    def sumFourDivisors(self, nums: List[int]) -> int:		
		mark = sieve(max(nums))
		
		def f(num, mark):
			s = 1 + num
			for i in range(2, int(num**(0.5))+1):
				if(i ** 3 == num):                     # case 1 :  n = p1 ^ 3 
					return s + i + num//i
				if(num % i == 0):
					s += i
					num = num//i
					if(mark[num] == True and num!=i):  # case 2 : n = p1 ^ 1 * p2 ^ 1 -> the "p2" must be prime and different with the "p1" 
						return s + num
					return 0
			return 0
```

##### Solution-2 ( Recursion ).
##### Solution Explanation
```
Need to check up to floor(sqrt(num)) = s (inclusive) only because 
for any divisor d < s, there is another divisor num/d > s.

E.g. 81 has divisors 1, 3, 9, 27, 81. sqrt(81) = 9 has divisors 1, 3, 9.
81/1 = 81, 81/3 = 27, 81/9 = 9. So if we only check for divisors up to 9 
and account for 81/divisor, we reduce time complexity by sqrt(num).
```
##### Complexity Analysis
```
Time  : O(N * sqrt(M)) where N = length of nums and M = nums[i] 
Space : O(1) since the length of set will not be more than 4
```
```python
from functools import cache
from math import isqrt

class Solution:
    def sumFourDivisors(self, nums: List[int]) -> int:
        #@lru_cache(None)
		#Note that @cache was introduced in Python 3.9, and it's a shorthand for @lru_cache(maxsize=None)
		@cache
        def f(x):
            s = set()
			#for i in range(1, floor(sqrt(x)) + 1):
            for i in range(1, isqrt(x)+1):
                if not x % i:
                    s.add(i)
                    s.add(x//i)
                    if len(s) > 4: break
            return sum(s) if len(s) == 4 else 0
        return sum(map(f, nums))
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

## Array
| #     | Title	                                         | url                                                                           | Time   | Space  | Difficulty | Tag	       | Note                   |
| ----- | ---------------------------------------------- | ----------------------------------------------------------------------------- | ------ | ------ | ---------- | ------------ | ---------------------- |
| 0001  | [Two Sum](#lc-1two-sum) | https://leetcode.com/problems/two-sum/                                        | _O(n)_ | _O(n)_ | Easy       |              |                        |
| 0121  | [Best Time to Buy and Sell Stock](#lc-121best-time-to-buy-and-sell-stock) | https://leetcode.com/problems/best-time-to-buy-and-sell-stock/                | _O(n)_ | _O(1)_ | Easy       |              |                        |
| 0217  | [Contains Duplicate](#lc-217contains-duplicate) | https://leetcode.com/problems/contains-duplicate/                             | _O(n)_ | _O(n)_ | Easy       |              |                        |
| 0238  | [Product of Array Except Self](#lc-238product-of-array-except-self) | https://leetcode.com/problems/product-of-array-except-self/                   | _O(n)_ | _O(1)_ | Medium     | LintCode     |                        |
| 0053  | [Maximum Subarray](#lc-53maximum-subarray) | https://leetcode.com/problems/maximum-subarray/                               | _O(n)_ | _O(1)_ | Medium     |              | `Kadane's Algorithm`   |
| 0152  | [Maximum Product Subarray](#lc-152maximum-product-subarray) | https://leetcode.com/problems/maximum-product-subarray/                       | _O(n)_ | _O(1)_ | Medium     |              |                        |
| 0153  | [Find Minimum in Rotated Sorted Array](#lc-153find-minimum-in-rotated-sorted-array) | https://leetcode.com/problems/find-minimum-in-rotated-sorted-array/           | _O(logn)_ | _O(1)_ | Medium  |              |                        |
| 0033  | [Search in Rotated Sorted Array](#lc-33search-in-rotated-sorted-array) | https://leetcode.com/problems/search-in-rotated-sorted-array/                 | _O(logn)_ | _O(1)_ | Medium  | CTCI         |                        |
| 0015  | [3 Sum](#lc-153-sum) | https://leetcode.com/problems/3sum/)                                          | _O(n^2)_  | _O(1)_ | Medium  |              | Two Pointers           |
| 0011  | [Container With Most Water](#lc-11container-with-most-water) | https://leetcode.com/problems/container-with-most-water/                      | _O(n)_ | _O(1)_ | Medium     |              |                        |

#### [LC-1:Two Sum](https://leetcode.com/problems/two-sum/)
##### Solution Explanation
```
Approach: One-pass Hash Table

While we iterate and insert elements into the table : 
 - we also look back to check if current element's complement already exists in the table. 
 - If it exists, we have found a solution and return immediately.
```
##### Complexity Analysis:
```
Time  : O(N)
========================
We traverse the list containing N elements only once. Each look up in the table costs only O(1) time.

Space : O(N)
========================
The extra space required depends on the number of items stored in the hash table, which stores at most N elements.
```
```python
def twoSum(nums: List[int], target: int) -> List[int]:
    dic = {}
    for i,num in enumerate(nums):
        if target-num in dic:
            return [dic[target-num], i]
        dic[num]=i

if __name__ == "__main__":
    #Input: nums = [2,7,11,15], target = 9
    #Output: [0,1]
    #Explanation: Because nums[0] + nums[1] == 9, we return [0, 1].
    nums = [2,7,11,15]
    target = 9
    print(twoSum(nums))
```
```kotlin
fun twoSum(nums: IntArray, target: Int): IntArray {
    val complements = HashMap<Int, Int>()
    nums.forEachIndexed { index, num ->
        with(target - num) { with(complements[this]) { this?.let { return intArrayOf(this, index) } } }
        complements[num] = index
    }
    throw IllegalArgumentException("No two sum solution")
}

fun main(args: Array<String>) {
    //Input: nums = [2,7,11,15], target = 9
    //Output: [0,1]
    //Explanation: Because nums[0] + nums[1] == 9, we return [0, 1].
    val nums = intArrayOf(2, 7, 11, 15)
    val target = 9
    //contentToString() deprecated from kotlin version 1.4 onwards
    //println(twoSum(nums, target).contentToString())
    println(twoSum(nums, target).joinToString(","))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-121:Best Time to Buy and Sell Stock](https://leetcode.com/problems/best-time-to-buy-and-sell-stock/)
##### Solution Explanation
```
Approach: Dynamic Programming + State Machine
=================================================================================================================================================================

- It is impossible to have stock to sell on first day, so -infinity is set as initial value for cur_hold and cur_not_hold is set to 0
- Iterate through the list of prices
  + Either keep in hold, or just buy today with stock price
  + Either keep in not holding, or just sell today with stock price
- Max profit must be in not-hold state

                                  Sell
                                  + stock price to balance
                      +------------------------------------------------+
                     /                                                  \
   ___    _________|/_                                                   \_________    ___
  /  _\| /         \                                                     /         \ |/_  \
 |      + Not Hold  +                                                   +    Hold   +      | Keep in hold
  \___/  \_________/                                                     \_________/  \___/
                    \                                                    /
 Keep in not holding \                                                  /
                      +------------------------------------------------+
                                  Buy
                                  - stock price from balance
```
##### Complexity Analysis:
```
Time  : O(N)
========================
We traverse the list containing n elements only once.

Space : O(1)
========================
No need for a dp array. We have replaced dp array with a variable called cur_not_hold.
```
```python
from typing import List

def maxProfit(prices: List[int]) -> int:
    if not price: return 0

    # It is impossible to have stock to sell on first day, so -infinity is set as initial value
    cur_hold, cur_not_hold = -float('inf'), 0
            
    for stock_price in prices:        
        prev_hold, prev_not_hold = cur_hold, cur_not_hold
            
        # either keep in hold, or just buy today with stock price
        cur_hold = max(prev_hold, -stock_price)
            
        # either keep in not holding, or just sell today with stock price
        cur_not_hold = max(prev_not_hold, prev_hold + stock_price)
        
    # max profit must be in not-hold state
    return cur_not_hold if prices else 0

if __name__ == "__main__":
    #Input: prices = [7,1,5,3,6,4]
    #Output: 5
    #Explanation: Buy on day 2 (price = 1) and sell on day 5 (price = 6), profit = 6-1 = 5.
    prices = [7,1,5,3,6,4]
    print(maxProfit(prices))
```
```kotlin
fun maxProfit(prices: IntArray): Int {
    //if (prices?.isEmpty() ?: true) return 0
    if (prices.isEmpty()) return 0
    // It is impossible to have stock to sell on first day, so -infinity is set as initial value
    var curHold = Int.MIN_VALUE
    var curNotHold = 0

    for (stockPrice in prices) {
        val prevHold = curHold
        val prevNotHold = curNotHold
            
        // either keep in hold, or just buy today with stock price
        curHold = maxOf(prevHold, stockPrice.unaryMinus())
            
        // either keep in not holding, or just sell today with stock price
        curNotHold = maxOf(prevNotHold, prevHold + stockPrice)
    }
    // max profit must be in not-hold state
    return curNotHold
}

fun main(args: Array<String>) {
    //Input: prices = [7,1,5,3,6,4]
    //Output: 5
    //Explanation: Buy on day 2 (price = 1) and sell on day 5 (price = 6), profit = 6-1 = 5.
    var prices = intArrayOf(7, 1, 5, 3, 6, 4)
    println(maxProfit(prices))
    //Input: prices = [7,6,4,3,1]
    //Output: 0
    //Explanation: In this case, no transactions are done and the max profit = 0.
    prices = intArrayOf(7, 6, 4, 3, 1)
    println(maxProfit(prices))
    //Input: prices = []
    //Output: 0
    prices = intArrayOf()
    println(maxProfit(prices))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-217:Contains Duplicate](https://leetcode.com/problems/contains-duplicate/)
##### Solution Explanation ( 2 possible solutions w/ their trade-offs explained ):
```
- For a sufficiently large value of N, choose Approach 2 ( Hash Table based approach ).
- When N is not sufficiently large, choose Approach 1 ( Sorting based approach ).

Two Possible Solutions (with their trade-offs explained)
=================================================================================================================================================================
Approach 1: Sorting
=================================================================================================================================================================
Intuition
-------------------------
If there are any duplicate integers, they will be consecutive after sorting.

Algorithm
-------------------------
- This approach employs sorting algorithm.
- Since comparison sorting algorithm like heapsort is known to provide O(N*log(N)) worst-case performance, sorting is often a good preprocessing step. 
- After sorting, we can sweep the sorted array to find if there are any two consecutive duplicate elements.

=================================================================================================================================================================
Complexity Analysis:
=================================================================================================================================================================

Time complexity : O(N*log(N))
========================
Sorting is O(N*log(N)) and the sweeping is O(N).
The entire algorithm is dominated by the sorting step, which is O(N*log(N)).

Space complexity : O(1)
========================
Space depends on the sorting implementation which, usually,
costs O(1) auxiliary space if heapsort is used.

Note
-------------------------
+ The implementation here modifies the original array by sorting it. 
+ In general, it is not a good practice to modify the input unless it is clear to the caller that the input will be modified. 
+ One may make a copy of nums and operate on the copy instead.

=================================================================================================================================================================
Approach 2: Hash Table
=================================================================================================================================================================
Intuition
-------------------------
Utilize a dynamic data structure that supports fast search and insert operations.

Algorithm
-------------------------
From Approach #1 we know that search operations is O(N) in an unsorted array and we did so repeatedly.
Utilizing a data structure with faster search time will speed up the entire algorithm.

There are many data structures commonly used as dynamic sets such as Binary Search Tree and Hash Table.
The operations we need to support here are search() and insert().
For a self-balancing Binary Search Tree (TreeSet or TreeMap in Java), search() and insert() are both O(log(N)) time.
For a Hash Table (HashSet or HashMap in Java), search() and insert() are both O(1) on average.
Therefore, by using hash table, we can achieve linear time complexity for finding the duplicate in an unsorted array.

=================================================================================================================================================================
Complexity Analysis:
=================================================================================================================================================================

Time complexity : O(N)
========================
We do search() and insert() for N times and each operation takes constant time.

Space complexity : O(N)
========================
The space used by a hash table is linear with the number of elements in it.

Note
-------------------------
+ For certain test cases with not very large N, the runtime of this method can be slower than Approach #2.
+ The reason is hash table has some overhead in maintaining its property.
+ One should keep in mind that real world performance can be different from what the Big-O notation says.
+ The Big-O notation only tells us that for sufficiently large input, one will be faster than the other.
+ Therefore, when N is not sufficiently large, an O(N) algorithm can be slower than an O(N*log(N)) algorithm.
```
```python
from typing import List

# Approach 1: Sorting [ TC: O(N*log(N)) ; SC: O(1) ]
def containsDuplicateApproachOne(nums: List[int]) -> bool:
    if len(nums) == 0 or len(nums) == 1: return False
    nums.sort()
    for i in range(len(nums) - 1):
        if nums[i] == nums[i + 1]:
            return True
    return False

# Approach 2: Hash Table [ TC: O(N) ; SC: O(N) ]
def containsDuplicateApproachTwo(nums: List[int]) -> bool:
    if len(nums) == 0 or len(nums) == 1: return False
    s = set()
    for num in nums:
        if num in s: return True
        s.add(num)
    return False

if __name__ == "__main__":
    #Input: nums = [1,2,3,1]
    #Output: true
    nums = [1,2,3,1]
    print(containsDuplicateApproachOne(nums))
    nums = [1,2,3,1]
    print(containsDuplicateApproachTwo(nums))
```
```kotlin
fun containsDuplicateApproachOne(nums: IntArray): Boolean {
    //if ( (nums?.isEmpty() ?: true) || (nums.size == 1) ) return false
    if ( (nums.isEmpty()) || (nums.size == 1) ) return false
    nums.sort()
    for (i in 0 until nums.size) {
        if (nums[i] == nums[i + 1]) {
            return true
        }
    }
    return false
}

fun containsDuplicateApproachTwo(nums: IntArray): Boolean {
    //if ( (nums?.isEmpty() ?: true) || (nums.size == 1) ) return false
    if ( (nums.isEmpty()) || (nums.size == 1) ) return false
    val numsSet = hashSetOf<Int>()
    //nums.forEach {value -> 
    //    if (numsSet.contains(value)) {
    //        return true
    //    }
    //    numsSet.add(value)
    //}
    for (num in nums) {
        if (numsSet.contains(num)) {
            return true
        }
        numsSet.add(num)
    }
    return false
}

fun main(args: Array<String>) {
    //Input: nums = [1,2,3,1]
    //Output: true
    var nums = intArrayOf(1, 2, 3, 1)
    println(containsDuplicateApproachOne(nums))
    nums = intArrayOf(1, 2, 3, 1)
    println(containsDuplicateApproachTwo(nums))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-238:Product of Array Except Self](https://leetcode.com/problems/product-of-array-except-self/)
##### Solution Explanation
```
Approach: Rolling Twice Algorithm
=================================================================================================================================================================

Thought Process:
---------------------------
- Look at the example array -> [1, 2, 3, 4]
- Somehow realize that:
  + output[i] = product([i+1:]) * product([:i])
  + which means output[i] = "product of all elements left of i" * "product of all elements right of i"
- Realize that we can find find "left products" and "right products" by iterating back and forth across the array.
  + Realize that we can store the left products in output (counts as constant space as per problem constraints) and compute the running right product and store it in a variable right (allows for constant space) as we iterate from the back of the input array to the front.
- Come up with Pseudocode:
  + Iterate through nums while computing "left products" and save them in the output array
  + Then, while going from right to left across the arrays:
    => store the running "right product" in a variable right
    => compute output[i] as we iterate
```
##### Complexity Analysis:
```
Time complexity : O(N) [ Technically O(2N) ]
========================
We traverse the list containing N elements twice. Each look up in the list costs only O(1) time.

Space complexity : O(1) [ As per problem, the output array does not count as extra space for space complexity analysis. ]
========================
Constant space since we only create a single output array to store the results.
```
```python
from typing import List

def productExceptSelf(nums: List[int]) -> List[int]:
    if not nums: return []
    size = len(nums)
        
    # create a list for output
    product_excluding_self = [ 1 for i in range(size) ]
        
    # Step_#1
    # record product of terms on the left hand side        
    for i in range( 1, size ):
        product_excluding_self[i] = product_excluding_self[i-1] * nums[i-1]
     
    # Step_#2
    # Update array elements as the product of ( product of left hand side ) * ( produt of right hand side )
    product_of_right_hand_side = 1
    for j in reversed( range( size) ):
        product_excluding_self[j] *= product_of_right_hand_side
        product_of_right_hand_side *= nums[j]
        
    return product_excluding_self

if __name__ == "__main__":
    #Input: nums = [1,2,3,4]
    #Output: [24,12,8,6]
    nums = [1,2,3,4]
    print(productExceptSelf(nums))
    #Input: nums = [-1,1,0,-3,3]
    #Output: [0,0,9,0,0]
    nums = [-1,1,0,-3,3]
    print(productExceptSelf(nums))
```
```kotlin
fun productExceptSelf(nums: IntArray): IntArray {
    // create a list for output
    val results = IntArray(nums.size)
    // edge case
    if (nums.isEmpty()) return results

    // Step_#1
    // record product of terms on the left hand side        
    results[0] = 1
    for (i in 1 until nums.size) {
        results[i] = results[i - 1] * nums[i - 1]
    }

    // Step_#2
    // Update array elements as the product of ( product of left hand side ) * ( produt of right hand side )
    var right = 1
    for (i in nums.size - 1 downTo 0) {
        results[i] *= right
        right *= nums[i]
    }
    return results
}

fun main(args: Array<String>) {
    //Input: nums = [1,2,3,4]
    //Output: [24,12,8,6]
    var nums = intArrayOf(1,2,3,4)
    //contentToString() deprecated from kotlin version 1.4 onwards
    //println(productExceptSelf(nums).contentToString())
    println(productExceptSelf(nums).joinToString(","))

    //Input: nums = [-1,1,0,-3,3]
    //Output: [0,0,9,0,0]
    nums = intArrayOf(-1,1,0,-3,3)
    //contentToString() deprecated from kotlin version 1.4 onwards
    //println(productExceptSelf(nums).contentToString())
    println(productExceptSelf(nums).joinToString(","))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-53:Maximum Subarray](https://leetcode.com/problems/maximum-subarray/)
##### Solution Explanation ( 2 possible solutions w/ their trade-offs explained ):
```
Approach 1: Kadane's Algorithm
-------------------------------
The largest subarray is either:
  - the current element
  - sum of previous elements

If the current element does not increase the value a local maximum subarray is found.

If local > global replace otherwise keep going

Problem is called Kadane's algorithm.

Reference: https://www.youtube.com/watch?v=86CQq3pKSUw
           https://medium.com/@rsinghal757/kadanes-algorithm-dynamic-programming-how-and-why-does-it-work-3fd8849ed73d


=================================================================================================================================================================
Complexity Analysis:
=================================================================================================================================================================

Time complexity : O(N)
========================
The time complexity of Kadane’s algorithm is O(N) because there is only one for loop which scans the entire array exactly once.

Space complexity : O(1)
========================
Kadane’s algorithm uses a constant space. So, the space complexity is O(1).


Approach 2: Divide and Conquer
-------------------------------
Explanation :

In the approach, we follow a divide and conquer approach similar to merge sort.

We use a helper functionhelper for this, wherein we pass in a starting index and the ending index to look at. The idea is to use this helper function recurseively.

Within the helper function, for a given start and end index, we find the mid of the array and split the array into two parts. Part 1 being the start ... mid and part 2 being mid+1 ... end. For each of the parts, we return 4 pieces of information.

1. The best possible answer within the subArray ==> ans
2. The maxSubarraySum starting at the beginning of the subArray ==> maxFromBegging
3. The maxSubarraySum ending at the end of the subArray ==> maxFromEnd
4. The total sum of the subArray ==> totalSum

With these four pieces of information for the two split parts, it is possible to combine them to generate a similar four pieces of information for the aggregated array. The trick to get an O(n) solution is to combine the information in a constant time.

The details of combing the information to curry up the information is as follows. (For each of 1-4 above, I will be prefexing left_ or right_ to denote they came from left/right subArray)

1. For part (1), we need to take the maximum of left_ans (best answer in left subarray), right_ans (nest answer in right subarray) and the crossover. The crossover is simply left_maxFromEnd + right_maxFromBeginning. The max these 3 components gives us the ans. For the highest level recursion, this is all we need.

2. For maxFromBeginning, we take the maximum among left_maxFromBeginning and left_totalSum + right_maxFromBeginning. This logically makes sense i.e either we want to take the result of left part or take the entire left part and merge it with the result from the right part.

3. For maxFromEnd, it is similar to above and we take the maximum among right_maxFromEnd and right_totalSum + left_maxFromEnd

4. For totalSum, we add the left_totalSum and right_totalSum.

Time Complexity : T(n) = 2*T(n/2) + 1 ==> O(n)

See this for proof https://youtu.be/OynWkEj0S-s?t=273

=================================================================================================================================================================
Complexity Analysis:
=================================================================================================================================================================

Time complexity : O(N)
========================
Time Complexity : T(n) = 2*T(n/2) + 1 ==> O(N) .. See this for proof https://youtu.be/OynWkEj0S-s?t=273.

Space complexity : O(log(N))
========================
Space Complexity : O(log(N)) since we are recursing and the call stack/number of recursive calls is of the order of log(N).
```
##### Solution-1 ( Kadane's Algorithm ):
```python
# If we are only interested in returning the sum of max sub-array
def maxSubArray(nums: List[int]) -> int:
    # largest subarray found in entire problem
    maxGlobal = nums[0]
    # current maximum subarray that is reset
    maxCurrent = nums[0]
        
    for i in range(1, len(nums)):
        maxCurrent  = max(nums[i], maxCurrent + nums[i])
        maxGlobal = max(maxCurrent, maxGlobal)
        
    return maxGlobal

if __name__ == "__main__":
    #Input: nums = [-2,1,-3,4,-1,2,1,-5,4]
    #Output: 6
    #Explanation: [4,-1,2,1] has the largest sum = 6.
    nums = [-2,1,-3,4,-1,2,1,-5,4]
    print(maxSubArray(nums))


# Variant: If we are interested in returning the actual max sub-array

# We can easily solve this problem in linear time using Kadane's algorithm.
# The idea is to maintain a maximum (positive-sum) subarray "ending" at each index of the given array.
# 	- The subarray is either empty (in which case its sum is zero) or 
#	- The subarray consists of one more element than the maximum subarray ending at the previous index.
 
#Variant (Also print the list)
# - Modify Kadane's algorithm which outputs only the sum of contiguous subarray with the largest sum but
#   does not print the subarray itself.
# - Keep track of the maximum subarray's starting and ending indices.

# Function to print contiguous sublist with the largest sum
# in a given set of integers
def maxSubArray(A: List[int]) -> List[int]:
 
    # stores maximum sum sublist found so far
    maxSoFar = 0
 
    # stores the maximum sum of sublist ending at the current position
    maxEndingHere = 0
 
    # stores endpoints of maximum sum sublist found so far
    start = end = 0
 
    # stores starting index of a positive-sum sequence
    beg = 0
 
    # traverse the given list
    for i in range(len(A)):
 
        # update the maximum sum of sublist "ending" at index `i`
        maxEndingHere = maxEndingHere + A[i]
 
        # if the maximum sum is negative, set it to 0
        if maxEndingHere < 0:
            maxEndingHere = 0        # empty sublist
            beg = i + 1
 
        # update result if the current sublist sum is found to be greater
        if maxSoFar < maxEndingHere:
            maxSoFar = maxEndingHere
            start = beg
            end = i
 
    #print(f"The sum of contiguous sublist with the largest sum is: {maxSoFar}")
    #print(f"The contiguous sublist with the largest sum is: {A[start: end + 1]}")
    return A[start: end + 1]
	
if __name__ == '__main__':
    A = [-2, 1, -3, 4, -1, 2, 1, -5, 4]
    print(maxSubArray(A))
```
```kotlin
fun maxSubArray(nums: IntArray): Int {
    var maxGlobal = nums[0] // largest subarray found in entire problem
    var maxCurrent = nums[0] // current maximum subarray that is reset

    for (i in 1 until nums.size) {
        maxCurrent = maxOf(nums[i], maxCurrent + nums[i])
        maxGlobal = maxOf(maxCurrent, maxGlobal)
    }

    return maxGlobal
}

fun main(args: Array<String>) {
    //Input: nums = [-2,1,-3,4,-1,2,1,-5,4]
    //Output: 6
    //Explanation: [4,-1,2,1] has the largest sum = 6.
    val nums = intArrayOf(-2,1,-3,4,-1,2,1,-5,4)
    println(maxSubArray(nums))
}
```
##### Solution-2 ( Divide and Conquer Algorithm ):
```python
def maxSubArray(nums: List[int]) -> int:
    def helper(nums, start, end):
        if start == end:
            return nums[start], nums[start], nums[start], nums[start]
        else:
            mid = start + (end - start)//2
                
        left_ans , left_maxFromBeginning , left_maxFromEnd , left_totalSum  = helper(nums, start, mid)
        right_ans, right_maxFromBeginning, right_maxFromEnd, right_totalSum = helper(nums, mid+1, end)
                
        ans = max(left_ans, right_ans, left_maxFromEnd + right_maxFromBeginning)
        maxFromBeginning = max(left_maxFromBeginning, left_totalSum + right_maxFromBeginning)
        maxFromEnd = max(right_maxFromEnd, right_totalSum + left_maxFromEnd)
        totalSum = left_totalSum + right_totalSum
                
        return (ans, maxFromBeginning, maxFromEnd, totalSum)

    ans, _, _, _ = helper(nums, 0, len(nums)-1)
    return ans

if __name__ == "__main__":
    #Input: nums = [-2,1,-3,4,-1,2,1,-5,4]
    #Output: 6
    #Explanation: [4,-1,2,1] has the largest sum = 6.
    nums = [-2,1,-3,4,-1,2,1,-5,4]
    print(maxSubArray(nums))
```
```kotlin
//The max subarray sum of an array of length 1 is equal to the unique element in that array.

//For arrays of size 2 or greater, the pivot index is used to divide the array into two halves.

//Case 1: The max subarray is divided between the two halves. In this case, the max subarray sum can be found by considering the sum of the largest sum found by iterating backward over the left half and the largest sum found by iterating forward over the right half.

//Case 2: The max subarray is not divided between the two halves. In this case, the max subarray sum of the overall array is equal to the maximum of the max subarray sums for the two array halves. This is the divide and conquer step.

fun maxLeftSumFromPivot(pivot: Int, startInclusive: Int): Int {
    var maxLeftSum = nums.get(pivot)
    var leftSum = nums.get(pivot)
    for (index in pivot - 1 downTo startInclusive) {
        leftSum += nums.get(index)
        maxLeftSum = Math.max(maxLeftSum, leftSum)
    }
    return maxLeftSum
}
    
fun maxRightSumFromPivot(pivot: Int, endExclusive: Int): Int {
    var maxRightSum = nums.get(pivot)
    var rightSum = nums.get(pivot)
    for (index in pivot + 1 until endExclusive) {
        rightSum += nums.get(index)
        maxRightSum = maxOf(maxRightSum, rightSum)
    }
    return maxRightSum
}

fun maxSubArray(startInclusive: Int, endExclusive: Int): Int {
    val diameter = endExclusive - startInclusive
    if (diameter == 1) {
        return nums.get(startInclusive)
    }
    val pivot = startInclusive + diameter / 2
    val leftMaxSubArray = maxSubArray(startInclusive, pivot)
    val rightMaxSubArray = maxSubArray(pivot, endExclusive)
    val maxRecursiveSubArray = Math.max(leftMaxSubArray, rightMaxSubArray)
    val maxLeftSumFromPivot = maxLeftSumFromPivot(pivot, startInclusive)
    val maxRightSumFromPivot = maxRightSumFromPivot(pivot, endExclusive)
    val maxSumFromPivot = maxLeftSumFromPivot + maxRightSumFromPivot - nums.get(pivot)
    return maxOf(maxRecursiveSubArray, maxSumFromPivot)
}

fun maxSubArray(nums: IntArray): Int {
    this.nums = nums
    return maxSubArray(0, nums.size)
}

fun main(args: Array<String>) {
    //Input: nums = [-2,1,-3,4,-1,2,1,-5,4]
    //Output: 6
    //Explanation: [4,-1,2,1] has the largest sum = 6.
    val nums = intArrayOf(-2,1,-3,4,-1,2,1,-5,4)
    println(maxSubArray(nums))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-152:Maximum Product Subarray](https://leetcode.com/problems/maximum-product-subarray/)
##### Solution Explanation
```
Kadane's Algorithm
=================================================================================================================================================================
In this solution we keep track of both the minimum and the maximum subarray encountered so far. 
This is done in a Kadane like fashion, where the updates for both the current minimal streak and maximal streak depend on the other, 
while the maximum encountered so far only depends on the current maximum, and the updated current maximal streak.
```
##### Complexity Analysis:
```
Time  : O(N)
========================
The time complexity of Kadane’s algorithm is O(N) because there is only one for loop which scans the entire array exactly once.

Space : O(1)
========================
Kadane’s algorithm uses a constant space. So, the space complexity is O(1).
```
```python
from typing import List

def maxProduct(nums: List[int]) -> int:
    ## RC ##
    ## APPROACH : KADANES ALGORITHM ##

    ## TIME COMPLEXITY : O(N) ##
    ## SPACE COMPLEXITY : O(1) ##

    # 1. Edge Case : Negative * Negative = Positive
    # 2. So we need to keep track of minimum values also, as they can yield maximum values.
    global_max = prev_max = prev_min = nums[0]
    for i in range(1, len(nums)):
        curr_min = min(prev_max*nums[i], prev_min*nums[i], nums[i])
        curr_max = max(prev_max*nums[i], prev_min*nums[i], nums[i])
        global_max = max(global_max, curr_max)
        prev_max = curr_max
        prev_min = curr_min
    return global_max

if __name__ == "__main__":
    #Input: nums = [2,3,-2,4]
    #Output: 6
    #Explanation: [2,3] has the largest product 6.    nums = [-2,1,-3,4,-1,2,1,-5,4]
    nums = [2,3,-2,4]
    print(maxProduct(nums))
```
```kotlin
fun maxProduct(nums: IntArray): Int {
    var maxGlobal = nums[0]
    var prevMax = nums[0]
    var prevMin = nums[0]
    for (i in 1 until nums.size) {
        val currMin = minOf(prevMax*nums[i], prevMin*nums[i], nums[i])
        val currMax = maxOf(prevMax*nums[i], prevMin*nums[i], nums[i])

        maxGlobal = maxOf(maxGlobal, currMax)
        prevMax = currMax
        prevMin = currMin
    }

    return maxGlobal
}

fun main(args: Array<String>) {
    //Input: nums = [2,3,-2,4]
    //Output: 6
    //Explanation: [2,3] has the largest product 6.    nums = [-2,1,-3,4,-1,2,1,-5,4]
    val nums = intArrayOf(2,3,-2,4)
    println(maxProduct(nums))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-153:Find Minimum in Rotated Sorted Array](https://leetcode.com/problems/find-minimum-in-rotated-sorted-array/)
##### Solution Explanation:
```
Binary Search Algorithm
=================================================================================================================================================================

Algorithm

1. Find the mid element of the array.

2. If mid element > first element of array this means that we need to look for the inflection point on the right of mid.

3. If mid element < first element of array this that we need to look for the inflection point on the left of mid.

          6 > 4
     +-------------+
     |             |
    \|/            |
 +---*--+------+---*--+------+------+------+
 |   4  |   5  |   6  |   7  |   2  |   3  |
 +------+------+------+------+------+------+
   Left           Mid                 Right
                     ---------------------->

In the above example mid element 6 is greater than first element 4. Hence we continue our search for the inflection point to the right of mid.

4 . We stop our search when we find the inflection point, when either of the two conditions is satisfied:

nums[mid] > nums[mid + 1] Hence, mid+1 is the smallest.

nums[mid - 1] > nums[mid] Hence, mid is the smallest.

                          +------+
                          |      |
                         \|/     |
 +------+------+------+---*--+---*--+------+
 |   4  |   5  |   6  |   7  |   2  |   3  |
 +------+------+------+------+------+------+
                        Left    Mid   Right

In the above example. With the marked left and right pointers. 
The mid element is 2. The element just before 2 is 7 and 7>2 i.e. nums[mid - 1] > nums[mid]. 
Thus we have found the point of inflection and 2 is the smallest element.

Detailed Algorithm
-----------------------

1) set left and right bounds
2) left and right both converge to the minimum index; DO NOT use left <= right because that would loop forever
  2.1) find the middle value between the left and right bounds (their average);
       can equivalently do: mid = left + (right - left) // 2,
       if we are concerned left + right would cause overflow (which would occur
       if we are searching a massive array using a language like Java or C that has
       fixed size integer types)
  2.2) the main idea for our checks is to converge the left and right bounds on the start
       of the pivot, and never disqualify the index for a possible minimum value.
  2.3) in normal binary search, we have a target to match exactly,
       and would have a specific branch for if nums[mid] == target.
       we do not have a specific target here, so we just have simple if/else.
  2.4) if nums[mid] > nums[right]
    2.4.1) we KNOW the pivot must be to the right of the middle:
           if nums[mid] > nums[right], we KNOW that the
           pivot/minimum value must have occurred somewhere to the right
           of mid, which is why the values wrapped around and became smaller.
    2.4.2) example:  [3,4,5,6,7,8,9,1,2]
           in the first iteration, when we start with mid index = 4, right index = 9.
           if nums[mid] > nums[right], we know that at some point to the right of mid,
           the pivot must have occurred, which is why the values wrapped around
           so that nums[right] is less then nums[mid]
    2.4.3) we know that the number at mid is greater than at least
           one number to the right, so we can use mid + 1 and
           never consider mid again; we know there is at least
           one value smaller than it on the right
  2.5) if nums[mid] <= nums[right]
    2.5.1) here, nums[mid] <= nums[right]:
           we KNOW the pivot must be at mid or to the left of mid:
           if nums[mid] <= nums[right], we KNOW that the pivot was not encountered
           to the right of middle, because that means the values would wrap around
           and become smaller (which is caught in the above if statement).
           this leaves the possible pivot point to be at index <= mid.

    2.5.2) example: [8,9,1,2,3,4,5,6,7]
           in the first iteration, when we start with mid index = 4, right index = 9.
           if nums[mid] <= nums[right], we know the numbers continued increasing to
           the right of mid, so they never reached the pivot and wrapped around.
           therefore, we know the pivot must be at index <= mid.

    2.5.3) we know that nums[mid] <= nums[right].
           therefore, we know it is possible for the mid index to store a smaller
           value than at least one other index in the list (at right), so we do
           not discard it by doing right = mid - 1. it still might have the minimum value.

3) at this point, left and right converge to a single index (for minimum value) since
   our if/else forces the bounds of left/right to shrink each iteration:

4) when left bound increases, it does not disqualify a value
   that could be smaller than something else (we know nums[mid] > nums[right],
   so nums[right] wins and we ignore mid and everything to the left of mid).

5) when right bound decreases, it also does not disqualify a
   value that could be smaller than something else (we know nums[mid] <= nums[right],
   so nums[mid] wins and we keep it for now).

6) so we shrink the left/right bounds to one value,
   without ever disqualifying a possible minimum.
```
##### Complexity Analysis:
```
Time  : O(log(N))
========================
Same as Binary Search O(log(N))

Space : O(1)
========================
```
```python
from typing import List

def findMin(nums: List[int]) -> int:
    """
    :type nums: List[int]
    :rtype: int
    """
    # set left and right bounds
    left, right = 0, len(nums)-1

    # left and right both converge to the minimum index;
    # DO NOT use left <= right because that would loop forever
    while left < right:
        # find the middle value between the left and right bounds (their average);
        # can equivalently do: mid = left + (right - left) // 2,
        # if we are concerned left + right would cause overflow (which would occur
        # if we are searching a massive array using a language like Java or C that has
        # fixed size integer types)
        #mid = (left + right) // 2
        mid = left + (right - left) // 2
            
        # the main idea for our checks is to converge the left and right bounds on the left
        # of the pivot, and never disqualify the index for a possible minimum value.

        # in normal binary search, we have a target to match exactly,
        # and would have a specific branch for if nums[mid] == target.
        # we do not have a specific target here, so we just have simple if/else.
            
        if nums[mid] > nums[right]:
            # we KNOW the pivot must be to the right of the middle:
            # if nums[mid] > nums[right], we KNOW that the
            # pivot/minimum value must have occurred somewhere to the right
            # of mid, which is why the values wrapped around and became smaller.

            # example:  [3,4,5,6,7,8,9,1,2] 
            # in the first iteration, when we left with mid index = 4, right index = 9.
            # if nums[mid] > nums[right], we know that at some point to the right of mid,
            # the pivot must have occurred, which is why the values wrapped around
            # so that nums[right] is less then nums[mid]

            # we know that the number at mid is greater than at least
            # one number to the right, so we can use mid + 1 and
            # never consider mid again; we know there is at least
            # one value smaller than it on the right
            left = mid + 1
        else:
            # here, nums[mid] <= nums[right]:
            # we KNOW the pivot must be at mid or to the left of mid:
            # if nums[mid] <= nums[right], we KNOW that the pivot was not encountered
            # to the right of middle, because that means the values would wrap around
            # and become smaller (which is caught in the above if statement).
            # this leaves the possible pivot point to be at index <= mid.
                
            # example: [8,9,1,2,3,4,5,6,7]
            # in the first iteration, when we left with mid index = 4, right index = 9.
            # if nums[mid] <= nums[right], we know the numbers continued increasing to
            # the right of mid, so they never reached the pivot and wrapped around.
            # therefore, we know the pivot must be at index <= mid.

            # we know that nums[mid] <= nums[right].
            # therefore, we know it is possible for the mid index to store a smaller
            # value than at least one other index in the list (at right), so we do
            # not discard it by doing right = mid - 1. it still might have the minimum value.
            right = mid

    # at this point, left and right converge to a single index (for minimum value) since
    # our if/else forces the bounds of left/right to shrink each iteration:

    # when left bound increases, it does not disqualify a value
    # that could be smaller than something else (we know nums[mid] > nums[right],
    # so nums[right] wins and we ignore mid and everything to the left of mid).

    # when right bound decreases, it also does not disqualify a
    # value that could be smaller than something else (we know nums[mid] <= nums[right],
    # so nums[mid] wins and we keep it for now).

    # so we shrink the left/right bounds to one value,
    # without ever disqualifying a possible minimum
    return nums[left]

if __name__ == "__main__":
    #Input: nums = [3,4,5,1,2]
    #Output: 1
    #Explanation: The original array was [1,2,3,4,5] rotated 3 times.
    nums = [3,4,5,1,2]
    print(findMin(nums))

### Uncommented concise solution
from typing import List

def findMin(nums: List[int]) -> int:
    left, right = 0, len(nums)-1
    while left < right:
        mid = left + (right - left) // 2            
        if nums[mid] > nums[right]:
            left = mid + 1
        else:
            right = mid

    return nums[left]

if __name__ == "__main__":
    #Input: nums = [3,4,5,1,2]
    #Output: 1
    #Explanation: The original array was [1,2,3,4,5] rotated 3 times.
    nums = [3,4,5,1,2]
    print(findMin(nums))
```
```kotlin
fun findMin(nums: IntArray): Int {
    var left = 0
    var right = nums.size - 1
    while (left < right) {
        //var mid = (right + left) / 2
        var mid = left + (right - left) / 2
        if (nums[mid] >= nums[left] && nums[mid] > nums[right]) {
            left = mid + 1
        } else {
            right = mid
        }
    }
    return nums[left]
}

fun main(args: Array<String>) {
    //Input: nums = [3,4,5,1,2]
    //Output: 1
    //Explanation: The original array was [1,2,3,4,5] rotated 3 times.
    val nums = intArrayOf(3,4,5,1,2)
    println(findMin(nums))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-33:Search in Rotated Sorted Array](https://leetcode.com/problems/search-in-rotated-sorted-array/)
##### Solution Explanation:
```
Binary Search Algorithm
=================================================================================================================================================================
Idea:
--------------------------
We have an ascending array, which is rotated at some pivot.
Let's call the rotation the inflection point. (IP)
One characteristic the inflection point holds is: arr[IP] > arr[IP + 1] and arr[IP] > arr[IP - 1]
So if we had an array like: [7, 8, 9, 0, 1, 2, 3, 4] the inflection point, IP would be the number 9.

One thing we can see is that values until the IP are ascending. And values from IP + 1 until end are also ascending (binary search, wink, wink).
Also the values from [0, IP] are always bigger than [IP + 1, n].

Intuition:
--------------------------
We can perform a Binary Search.
If A[mid] is bigger than A[left] we know the inflection point will be to the right of us, meaning values from a[left]...a[mid] are ascending.

So if target is between that range we just cut our search space to the left.
Otherwise go right.

The other condition is that A[mid] is not bigger than A[left] meaning a[mid]...a[right] is ascending.
In the same manner we can check if target is in that range and cut the search space correspondingly.
```
##### Complexity Analysis:
```
Time Complexity : O(log(N))
========================
Same as Binary Search O(log(N))

Space Complexity : O(1)
========================
```
```python
from typing import List

def search(nums: List[int], target: int) -> int:
    n = len(nums)
    left, right = 0, n - 1
    if n == 0: return -1
        
    while left <= right:
        mid = left + (right - left) // 2
        if nums[mid] == target: return mid
            
        # inflection point to the right. Left is strictly increasing
        if nums[mid] >= nums[left]:
            if nums[left] <= target < nums[mid]:
                right = mid - 1
            else:
                left = mid + 1
                    
        # inflection point to the left of me. Right is strictly increasing
        else:
            if nums[mid] < target <= nums[right]:
                left = mid + 1
            else:
                right = mid - 1
            
        return -1

if __name__ == "__main__":
    #Input: nums = [4,5,6,7,0,1,2], target = 0
    #Output: 4
    nums = [4,5,6,7,0,1,2]
    target = 0
    print(search(nums, target))
```
```kotlin
fun search(nums: IntArray, target: Int): Int {
    var left = 0
    var right = nums.size - 1
    while (left <= right) {
        val mid = left + (right - left) / 2
        when {
            nums[mid] == target -> return mid
            nums[left] <= nums[mid] -> if (target in nums[left] .. nums[mid]) right = mid - 1 else left = mid + 1
            nums[mid] <= nums[right] -> if (target in nums[mid] .. nums[right]) left = mid + 1 else right = mid - 1
        }
    }
    return -1
}

fun main(args: Array<String>) {
    //Input: nums = [4,5,6,7,0,1,2], target = 0
    //Output: 4
    val nums = intArrayOf(4,5,6,7,0,1,2)
    val target = 0
    println(search(nums, target))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-15:3Sum](https://leetcode.com/problems/3sum/)
##### Solution Explanation:
```
Two Pointers
=================================================================================================================================================================
1) Let n be the size of the list.
2) First sort the list.
3) Iterate through the list from index => i = 0 to i < n-2.
   3.1) Since the list is sorted, if nums[i] > 0, then all 
        nums[j] with j > i are positive as well, and we cannot
        have three positive numbers sum up to 0. Return immediately.
   3.2) The nums[i] == nums[i-1] condition helps us avoid duplicates.
        E.g., given [-1, -1, 0, 0, 1], when i = 0, we see [-1, 0, 1]
        works. Now at i = 1, since nums[1] == -1 == nums[0], we avoid
        this iteration and thus avoid duplicates. The i > 0 condition
        is to avoid negative index, i.e., when i = 0, nums[i-1] = nums[-1]
        and you don't want to skip this iteration when nums[0] == nums[-1]
   3.3) Set left pointer to i+1 and right pointer to n-1
   3.4) While left pointer is less than right pointer ( Now we have a Classic two pointer solution )
        3.4.1) calculate sum => sum = nums[i] + nums[left] + nums[right]
        3.4.2) sum too small ( sum < 0 ), move left ptr
        3.4.3) sum too large ( sum > 0 ), move right ptr
        3.4.4) sum == 0
               3.4.4.1) we need to skip elements that are identical to our
                        current solution, otherwise we would have duplicated triples
```
##### Complexity Analysis:
```
Time  : O(N^2) [ Quadratic ]
========================
2SUM time complexity is O(N).
Every higher SUM (k-sum) will have complexity => O(N^k-1) ... So for 3SUM with k=3 we get O(N^(3-1)) = O(N^2).

reason for k-1: last 2 depths are solved by 2 sum sorted

Details
-------
python built-in sort method is using timSort. Therefore bestcase time complexity is O(N), otherwise O(NlogN).
The rest of the algorithm takes O(N^2).

Space : O(N)
========================
If we do not consider the result list, the space complexity is bounded by O(N).

If we consider the result list, then space complexity becomes => O(N^2).
```
```python
from typing import List

def threeSum(nums: List[int]) -> List[List[int]]:
    if len(nums) < 3: return []
    res = []        # Triples
    n = len(nums)   # Length of the list
    nums.sort()     # We need to sort the list first!
        
    for i in range(n-2):            
        # Since the list is sorted, if nums[i] > 0, then all 
        # nums[j] with j > i are positive as well, and we cannot
        # have three positive numbers sum up to 0. Return immediately.
        if nums[i] > 0:
            break
                
        # The nums[i] == nums[i-1] condition helps us avoid duplicates.
        # E.g., given [-1, -1, 0, 0, 1], when i = 0, we see [-1, 0, 1]
        # works. Now at i = 1, since nums[1] == -1 == nums[0], we avoid
        # this iteration and thus avoid duplicates. The i > 0 condition
        # is to avoid negative index, i.e., when i = 0, nums[i-1] = nums[-1]
        # and you don't want to skip this iteration when nums[0] == nums[-1]
        if i > 0 and nums[i] == nums[i-1]:
            continue
                
        # Classic two pointer solution
        left, right = i + 1, n - 1
        while left < right:
            sumOfNums = nums[i] + nums[left] + nums[right]
            if sumOfNums < 0: # sum too small, move left ptr
                left += 1
            elif sumOfNums > 0: # sum too large, move right ptr
                right -= 1
            else:
                res.append([nums[i], nums[left], nums[right]])
                    
                # we need to skip elements that are identical to our
                # current solution, otherwise we would have duplicated triples
                while left < right and nums[left] == nums[left+1]:
                    left += 1
                while left < right and nums[right] == nums[right-1]:
                    right -= 1
                left += 1
                right -= 1
    return res

if __name__ == "__main__":
    #Input: nums = [-1,0,1,2,-1,-4]
    #Output: [[-1,-1,2],[-1,0,1]]
    nums = [-1,0,1,2,-1,-4]
    print(threeSum(nums))

# Concise w/o comments solution
from typing import List

def threeSum(nums: List[int]) -> List[List[int]]:
    if len(nums) < 3: return []
    res = []        # Triples
    n = len(nums)   # Length of the list
    nums.sort()     # We need to sort the list first!
        
    for i in range(n-2):            
        if nums[i] > 0:
            break
                
        if i > 0 and nums[i] == nums[i-1]:
            continue
                
        left, right = i + 1, n - 1
        while left < right:
            sumOfNums = nums[i] + nums[left] + nums[right]
            if sumOfNums < 0: # sum too small, move left ptr
                left += 1
            elif sumOfNums > 0: # sum too large, move right ptr
                right -= 1
            else:
                res.append([nums[i], nums[left], nums[right]])
                    
                while left < right and nums[left] == nums[left+1]:
                    left += 1
                while left < right and nums[right] == nums[right-1]:
                    right -= 1
                left += 1
                right -= 1
    return res

if __name__ == "__main__":
    #Input: nums = [-1,0,1,2,-1,-4]
    #Output: [[-1,-1,2],[-1,0,1]]
    nums = [-1,0,1,2,-1,-4]
    print(threeSum(nums))
```
```kotlin
fun threeSum(nums: IntArray): List<List<Int>> {
    val ans:MutableList<List<Int>> = mutableListOf()
    if (nums.size < 3) return ans
    nums.sort()
    for (i in 0 until nums.size - 2) {
        if (nums[i] > 0) break
        if ((i > 0) && (nums[i] == nums[i-1])) {
            continue
        }
        var left = i+1
        var right = nums.size -1
        while (left < right) {
            val sumOfNums = nums[i] + nums[left] + nums[right]
            if (sumOfNums < 0) {
                left++
            } else if (sumOfNums > 0) {
                right--
            } else {
                ans.add(listOf(nums[i],nums[left],nums[right]))
                while ((left < right) && (nums[left] == nums[left+1])) left++
                while ((left < right) && (nums[right] == nums[right-1])) right--
                left++
                right--
            }
        }
    }
    return ans.toList()        
}

fun main(args: Array<String>) {
    //Input: nums = [-1,0,1,2,-1,-4]
    //Output: [[-1,-1,2],[-1,0,1]]
    val nums = intArrayOf(-1,0,1,2,-1,-4)
    println(threeSum(nums))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-11:Container With Most Water](https://leetcode.com/problems/container-with-most-water/)
##### Formal Proof than an O(N) solution exists:
```
Formal Proof than an O(N) solution exists:
=================================================================================================================================================================
Problem Description:
Condition: We have two pointers at i & j, suppose h[i] <= h[j].
Goal to Prove: If there is a better answer within the sub-range of [i, j], then the range [i + 1, j] must contain that optimal sub-range. (This doesn't mean range [i, j - 1] can't contain it, we just want to prove range [i + 1, j] will contain it).

Proof:
Since we assume there is a better answer in the sub-range of [i, j], then this optimal range must be contained by either [i + 1, j] or [i, j - 1], or both.

Let's assume [i + 1, j] doesn't contain the optimal range, but [i, j - 1] contains it. Then this means two things:

the optimal range is not in [i + 1, j - 1], otherwise, [i + 1, j] will contain it.
The optimal range contains the block [i, i + 1] (since this is the part which exists in [i, j - 1] but not in [i+1, j]).
However, notice that, len(i, j - 1) < len(i, j), and in the range [i, j], the area is constrained by the height of h[i] (even in the case h[i] == h[j]). Thus, in the range [i, j - 1], even all pillar are no shorter than h[j], the maximum area is smaller than the area formed by i & j, which contradicts our assumption there is a better answer in the sub-range of [i, j]. This contradiction suggests [i + 1, j] contains the optimal sub-range, if such sub-range exists.

According to above theorem, we can design the algorithm, whenever h[i] < h[j], we advance the pointer i.
```
##### Solution Explanation:
```
=================================================================================================================================================================
Solution Explanation:
Two Pointers
=================================================================================================================================================================
- O(N) solution which is explained in editorial. Why the solution works needs some thought.
- Use two pointers start and end initialized at 0 and N-1
- Now compute the area implied by these pointers as (end-start) * min (height[start], height[end])
- if height[start] < height[end], start = start + 1 else end = end -1
- Why? Imagine height[start] < height[end]. Then is there any need to compare height[end-1] with height[start]? 
  There is no way we can now get a larger area using height[start] as one of the end points. We should therefore move start.
```
##### Complexity Analysis:
```
Time  : O(N)
========================
We traverse the list containing N elements only once. Each look up in the table costs only O(1) time.

Space : O(1)
========================
The constant space required by the output variable max_water.
```
```python
from typing import List

def maxArea(height: List[int]) -> int:
    """
    :type height: List[int]
    :rtype: int
    """
    start, end = 0, len(height)-1
    max_water = -1
    while start < end:
        max_water = max(max_water, (end-start)*min(height[start], height[end]))
        if height[start] < height[end]:
            start = start+1
        else:
            end = end-1
    return max_water

if __name__ == "__main__":
    #Input: height = [1,8,6,2,5,4,8,3,7]
    #Output: 49
    #Explanation: The above vertical lines are represented by array [1,8,6,2,5,4,8,3,7]. 
    #             In this case, the max area of water (blue section) the container 
    #             can contain is 49.
    height = [1,8,6,2,5,4,8,3,7]
    print(maxArea(height))
```
```kotlin
fun maxArea(height: IntArray): Int {    
    val ans:MutableList<List<Int>> = mutableListOf()
    var start = 0
    var end = height.size - 1
    var maxWater = -1
    while (start < end) {
        maxWater = maxOf(maxWater, (end-start)*minOf(height[start], height[end]))
        if (height[start] < height[end]) {
            start++
        } else {
            end--
        }
    }
    return maxWater
}

fun main(args: Array<String>) {
    //Input: height = [1,8,6,2,5,4,8,3,7]
    //Output: 49
    //Explanation: The above vertical lines are represented by array [1,8,6,2,5,4,8,3,7]. 
    //             In this case, the max area of water (blue section) the container 
    //             can contain is 49.
    val height = intArrayOf(1,8,6,2,5,4,8,3,7)
    println(maxArea(height))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

## Bit Manipulation
| #     | Title	                                              | url                                                                           | Time     | Space    | Difficulty | Tag	        | Note                   |
| ----- | --------------------------------------------------- | ----------------------------------------------------------------------------- | -------- | -------- | ---------- | ------------ | ---------------------- |
| 0371  | [Sum of Two Integers](#lc-371sum-of-two-integers)   | https://leetcode.com/problems/sum-of-two-integers/                            | _O(1)_   | _O(1)_   | Easy       | LintCode     |                        |
| 0191  | [Number of 1 Bits](#lc-191number-of-1-bits)         | https://leetcode.com/problems/number-of-1-bits/                               | _O(1)_   | _O(1)_   | Easy       |              |                        |
| 0338  | [Counting Bits](#lc-338counting-bits)               | https://leetcode.com/problems/counting-bits/                                  | _O(n)_   | _O(n)_   | Medium     |              |                        |
| 0268  | [Missing Number](#lc-268missing-number)             | https://leetcode.com/problems/missing-number/                                 | _O(n)_   | _O(1)_   | Medium     | LintCode     |                        |
| 0190  | [Reverse Bits](#lc-190reverse-bits)                 | https://leetcode.com/problems/reverse-bits/                                   | _O(1)_   | _O(1)_   | Easy       |              |                        |

#### [LC-371:Sum of Two Integers](https://leetcode.com/problems/sum-of-two-integers/)
##### Solution Explanation:
```
=====================================
Java/Kotlin
=====================================
For this problem, the main crux is that, we are dividing the task of adding 2 numbers into two parts -

Let a = 13 and b = 10. Then we want to add them, a+b
In binary, it would look as follows -
a      =  1 1 0 1
b      =  1 0 1 0
--------------------
a+b = (1) 0 1 1 1
Now we can break the addition into two parts, one is simple addition without taking care of carry, and other is taking the carry.
With that strategy in mind, we have the following -
simpleAddition(a, b):
a      =  1 1 0 1
b      =  1 0 1 0
--------------------
a+b    =  0 1 1 1  

carry(a, b):
carry = 1 0 0 0 0
*shift left by one, because we add the carry next left step :P
Now if we can add the carry to our simpleAddition result, we can get our final answer. So add them using the simpleAddition method, and take care of the new carry again, unless carry becomes zero.
This simpleAddition is performed by XOR ^operator, and the carry is performed by AND &operator.
So our final answer would look as -
(a^b) =    0 1 1 1 
+carry = 1 0 0 0 0 
-----------------------------
ans =    1 0 1 1 1
-----------------------------
which is 23

References:
===============================================================
https://en.wikipedia.org/wiki/Adder_%28electronics%29#Half_adder

The half adder adds two single binary digits A and B. It has two outputs, sum (S) and carry (C).
The carry signal represents an overflow into the next digit of a multi-digit addition. 
The value of the sum is 2C + S. 

=====================================
Python
=====================================
Python doesn't respect this int boundary that other strongly typed languages like Java and C++ have defined.
So we need to use a mask.

=====================================
Let's recall the rule for taking two's complements: Flip all the bits, then plus one.

So, to take the two's complement of -20 in the 32 bits sense. We flip all the 32 bits of 0xFFFFFFEC and add 1 to it.
Note that here we cannot use the bit operation ~ because it will flip infinite many bits, not only the last 32.
Instead, we xor it with the mask 0xFFFFFFFF. Recall that xor with 1 has the same effect as flipping.
This only flips the last 32 bits, all the 0's to the far left remains intact.
Then we add 1 to it to finish the two's complement and produce a valid 20

(0xFFFFFFEC^mask)+1 == 0x14 == 20
Next, we take the two's complement of 20 in the Python fashion. Now we can direcly use the default bit operation

~20+1 == -20
Write these two steps in one line

~((0xFFFFFFEC^mask)+1)+1 == -20 == 0x...FFFFFFFFFFFFFFEC
Wait a minute, do you spot anything weird? We are not supposed to use + in the first place, right? Why are there two +1's?
Does it mean this method won't work? Hold up and let me give the final magic of today:

for any number x, we have

~(x+1)+1 = ~x
(Here the whole (0xFFFFFFEC^mask) is considered as x).

In other words, the two +1's miracly cancel each other! so we can simly write

~(0xFFFFFFEC^mask) == -20
To sum up, since Python allows arbitary length for integers, we first use a mask 0xFFFFFFFF to restrict the lengths.
But then we lose information for negative numbers, so we use the magical formula ~(a^mask) to convert the result to Python-interpretable form.
```
##### Complexity Analysis:
```
Time  : O(1)
Space : O(1)
```
```python
def getSum(a: int, b: int) -> int:
    """
    :type a: int
    :type b: int
    :rtype: int
    """
    mask = 0xffffffff
    while b:
        sum = (a^b) & mask
        carry = ((a&b)<<1) & mask
        a = sum
        b = carry

    if (a>>31) & 1: # If a is negative in 32 bits sense
        return ~(a^mask)
    return a

if __name__ == "__main__":
    #Input: a = 1, b = 2
    #Output: 3
    a = 1
    b = 2
    print(getSum(a, b))
```
```kotlin
fun getSum(a: Int, b: Int): Int {
    var a = a
    var b = b
    while (b != 0) {
        var carry = a.and(b) // carry contains common set bits
        a = a.xor(b) // sum of bits where at least 1 common bit is not set
        carry = carry.shl(1) // carry needs to be added 1 place left side
        b = carry
    }
    return a
}

fun main(args: Array<String>) {
    //Input: a = 1, b = 2
    //Output: 3
    val a = 1
    val b = 2
    println(getSum(a, b))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-191:Number of 1 Bits](https://leetcode.com/problems/number-of-1-bits/)
##### Solution Explanation:
```
=================================================================================================================================================================
Approach-1 ( Using masking )
=================================================================================================================================================================
Solution Explanation:
Using masking
=================================================================================================================================================================
The best solution for this problem is to use "divide and conquer" to count bits:

- First, count adjacent two bits, the results are stored separatedly into two bit spaces;
- Second is to count the results from each of the previous two bits and store results to four bit spaces;
- Repeat those steps and the final result will be sumed.
- Check the following diagram from Hack's Delight book to understand the procedure:


x = (x & 0x55555555) + ((x >> 1) & 0x55555555);
x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
x = (x & 0x0F0F0F0F) + ((x >> 4) & 0x0F0F0F0F);
x = (x & 0x00FF00FF) + ((x >> 8) & 0x00FF00FF);
x = (x & 0x0000FFFF) + ((x >> 16) & 0x0000FFFF);

The first line uses (x >> 1) & 0x55555555 rather than the perhaps more natural (x & 0xAAAAAAAA) >> 1,
because the code shown avoids generating two large constants in a register. This would cost an
instruction if the machine lacks the and not instruction. A similar remark applies to the other lines.
Clearly, the last and is unnecessary, and other and’s can be omitted when there is no danger that a
field’s sum will carry over into the adjacent field. Furthermore, there is a way to code the first line
that uses one fewer instruction. This leads to the simplification shown in Figure 5–2, which executes
in 21 instructions and is branch-free.

Resource: https://doc.lagout.org/security/Hackers%20Delight.pdf (Chapter 5)


=================================================================================================================================================================
Approach-2 ( Using Brian Kernighan Algorithm )
=================================================================================================================================================================
Solution Explanation:
Kernighan way
=================================================================================================================================================================
If we can somehow get rid of iterating through all the 32 bits and only iterate as many times as there are 1's, wouldn't that be better?
Below is a solution that does this. It's based on Kernighan's number of set bits counting algorithm.
=================================================================================================================================================================
To solve this problem efficiently one must be familiar with Brian Kernighans Bit Manipulation Algorithm 
which is used to count the number of set bits in a binary representation of a number k. (a bit is considered set if it has the value of 1)

Resource: https://www.techiedelight.com/brian-kernighans-algorithm-count-set-bits-integer/
=================================================================================================================================================================

Using Brian Kernighan Algorithm, we will not check/compare or loop through all the 32 bits present
but only count the set bits which is way better than checking all the 32 bits

Suppose we have a number 00000000000000000000000000010110 (32 bits).

Now using this algorithm we will skip the 0's bit and directly jump to set bit(1's bit) 
and we don't have to go through each bit to count set bits i.e. the loop will be executed 
only for 3 times in the mentioned example and not for 32 times.


Assume we are working for 8 bits for better understanding, but the same logic apply for 32 bits
So, we will take a number having 3 set bits.
n = 00010110
n - 1 = 00010101 (by substracting 1 from the number, all the bits gets flipped/toggled after the MSB(most significant right bit) including the MSB itself
After applying &(bitwise AND) operator on n and n - 1 i.e. (n & n - 1), the righmost set bit will be turned off/toggled/flipped

Let's understand step by step:
===============================
* 1st Iteration
     00010110 --> (22(n) in decimal)
  &  00010101 --> (21(n - 1) in decimal i.e. flipping all the bits of n(22) after MSB including the MSB)
  -----------
     00010100 --> (20(n & n - 1) in decimal i.e after applying bitwise AND(&), the MSB will be turned off)

After applying bitwise AND(&) ,assign this number to n i.e. n = n & n - 1
n = 00010100(20 in decimal)
and increase the count
bitCount++ (Initial bitCount = 0. By incrementing it, the bitCount = 1)
-------------------------------------------------------------------------------------------------------------------------------
* 2nd Iteration
     00010100 --> (20(n) in decimal)
  &  00010011 --> (19(n - 1) in decimal i.e. flipping all the bits of n(20) after MSB including the MSB)
  -----------
     00010000 --> (16(n & n - 1) in decimal i.e after applying bitwise AND(&), the MSB will be turned off)

After applying bitwise AND(&) ,assign this number to n i.e. n = n & n - 1
n = 00010000(16 in decimal)
and increase the count
bitCount++ (previous bitCount = 1. By incrementing it, the bitCount = 2)
-------------------------------------------------------------------------------------------------------------------------------
* 3rd Iteration
     00010000 --> (16(n) in decimal)
  &  00001111 --> (15(n - 1) in decimal i.e. flipping all the bits of n(16) after MSB including the MSB)
  -----------
     00000000 --> (0(n & n - 1) in decimal i.e after applying bitwise AND(&), the MSB will be turned off)

After applying bitwise AND(&) ,assign this number to n i.e. n = n & n - 1
n = 00000000 (0 in decimal)
and increase the count
bitCount++ (previous bitCount = 2. By incrementing it, the bitCount = 3)
-------------------------------------------------------------------------------------------------------------------------------

Now, since the n = 0, there will be no furthur iteration as the condition becomes false, so it will come 
out of the loop and return bitCount which is 3 which is desired output.
```
##### Complexity Analysis:
```
For both approaches:

Time  : O(1)
Space : O(1)
```
```python
# Approach-1 ( Using masking )
def hammingWeight(n: int) -> int:
    mask_sum_2bit = 0x55555555
    mask_sum_4bit = 0x33333333
    mask_sum_8bit = 0x0F0F0F0F
    mask_sum_16bit = 0x00FF00FF
    mask_sum_32bit = 0x0000FFFF

    n = (n & mask_sum_2bit) + ((n >> 1) & mask_sum_2bit)
    n = (n & mask_sum_4bit) + ((n >> 2) & mask_sum_4bit)
    n = (n & mask_sum_8bit) + ((n >> 4) & mask_sum_8bit)
    n = (n & mask_sum_16bit) + ((n >> 8) & mask_sum_16bit)
    n = (n & mask_sum_32bit) + ((n >> 16) & mask_sum_32bit)

    return n

if __name__ == "__main__":
    #Input: n = 00000000000000000000000000001011
    #Output: 3
    #Explanation: The input binary string 00000000000000000000000000001011 has a total of three '1' bits.
    n = 0b00000000000000000000000000001011
    print(hammingWeight(n))

# Approach-2 ( Using Brian Kernighan Algorithm )
def hammingWeight(n: int) -> int:
    count = 0
    while n:
        count += 1
        n = n & (n - 1)    
    return count


if __name__ == "__main__":
    #Input: n = 00000000000000000000000000001011
    #Output: 3
    #Explanation: The input binary string 00000000000000000000000000001011 has a total of three '1' bits.
    n = 0b00000000000000000000000000001011
    print(hammingWeight(n))
```
```kotlin
// Approach-1 ( Using masking )
fun hammingWeight(n:Int):Int {
    var num = n
    val mask_sum_2bit: Int = 0x55555555.toInt()
    val mask_sum_4bit: Int = 0x33333333.toInt()
    val mask_sum_8bit: Int = 0x0F0F0F0F.toInt()
    val mask_sum_16bit: Int = 0x00FF00FF.toInt()
    val mask_sum_32bit: Int = 0x0000FFFF.toInt()

    
    num = ((num and 0xAAAAAAAA.toInt()) ushr 1) + (num and mask_sum_2bit)
    num = ((num and 0xCCCCCCCC.toInt()) ushr 2) + (num and mask_sum_4bit)
    num = ((num and 0xF0F0F0F0.toInt()) ushr 4) + (num and mask_sum_8bit)
    num = ((num and 0xFF00FF00.toInt()) ushr 8) + (num and mask_sum_16bit)
    num = ((num and 0xFFFF0000.toInt()) ushr 16) + (num and mask_sum_32bit)
    return num
}

fun main(args: Array<String>) {
    //Input: n = 00000000000000000000000000001011
    //Output: 3
    //Explanation: The input binary string 00000000000000000000000000001011 has a total of three '1' bits.
    val n = 0b00000000000000000000000000001011
    println(hammingWeight(n))
}

// Approach-2 ( Using Brian Kernighan Algorithm )
fun hammingWeight(n:Int):Int {
    var num = n
    var count = 0
    while (num != 0) {
        num = num and (num - 1)
        count++      
    }
   
   return count
}

fun main(args: Array<String>) {
    //Input: n = 00000000000000000000000000001011
    //Output: 3
    //Explanation: The input binary string 00000000000000000000000000001011 has a total of three '1' bits.
    val n = 0b00000000000000000000000000001011
    println(hammingWeight(n))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-338:Counting Bits](https://leetcode.com/problems/counting-bits/)
##### Prerequisite:
```
=================================================================================================================================================================
To solve this problem efficiently one must be familiar with Brian Kernighans Bit Manipulation Algorithm 
which is used to count the number of set bits in a binary representation of a number k. (a bit is considered set if it has the value of 1)

Resource: https://www.techiedelight.com/brian-kernighans-algorithm-count-set-bits-integer/
```
##### Solution Explanation:
```
Intuition:
=================================================================================================================================================================

This problem can be solved using dynamic programming combined with a bit manipulation technique.

Overview
------------------
Recall the goal is to count the number of set bits in the binary representation of all integers 0 to n and record the set bits count of each number.
This can be done in O(N) time and using O(N) space.

To solve this problem efficiently one must be familiar with Brian Kernighans Bit Manipulation Algorithm which is used to count the number 
of set bits in a binary representation of a number k. (a bit is considered set if it has the value of 1)

Intially it seems just to be aware of Kernighans algorithm is enough to solve this problem, but this is not the case. 
It would be naive to calculate the set bit count of each number from 1 to n and record the result for each number.
This is because Kernighans Algorithm has a runtime of of O(S) where S is the number of set bits in a number. 
This is because in each iteration a set bits in the original number is changed from 1 to 0 until there are no more set bits. 
This must be done for all N numbers.
The runtime of this approach is O(N*S) time. We can do better.


Optimization
We can use dynamic programming to eliminate uncessary work.

The uncessary work is calculating the set bit count for each number from 0 to n.

if we use a cache and leverage the heart of kernighans algorithm, we can avoid explictly calculating the set bit count for each number reducing the runtime to O(N)

The heart of kernighans algorithm uses a bit manipulation techinuque to turn off (set 1 to 0) the least significant bit (rightmost) in a number. an AND operation is perfomerd be between the binary representations of k and k-1. this results in a number m who's set bits count is one less than k. The algorithm performs this operation on each iteration, eahc time updating k to the value of m. in python this operation is k = k & (k - 1)

The key observation for optimization is that, if we count the set bits of each number in sorted order and we cache the results. for any given number k which no set bit count has been recorded, we can always lookup a number m who has a set bits count that is one less than k. if we know the set bit count of m, we can get the set bits count of k by adding one to the set bits count of m. lookup is possible because we are going in sorted order and it is the case m <= k.

if we use kernighans bit manipulation technique, m can always be found in constant time.
if we start at 0 we can build our way up to n. obtaining all counts along the way.

if n = 4
0 -> 0 set bits by default
thus dp[0] = 0

Now how do we count the set bits of the number 1?
Using kernighans technique we can find a number that has one less set bit than the binary representation of 1. that number is zero.

1 -> 1 & (1 - 1) = 0

You can think of 1 ask and 0 asm from above explanation

Thus if we add 1 + 0 we get the set bit count for the binary representation of one.
count set bits in binary representation of 1 dp[1] = 1 + dp[0] = 1 + 0 = 1

count set bits in binary representation of  2 
dp = [0, 1, 0, 0, 0]
2 -> 2 & (2 - 1) = 0
base 2: 0010  & 0001 = 0000
dp[2] = 1 + dp[0] =  0 + 1 = 1 

count set bits in binary representation of  3 
dp = [0, 1, 1, 0, 0]
3 -> 3 & (3 - 1) = 3 & 2 = 2
base 2: 0011  & 0010 = 0010 
dp[3] = 1 + dp[2] =  1 + 1 = 2

count set bits in binary representation of  4 
dp = [0, 1, 1, 2, 0]
4 -> 4 & (4 - 1) = 4 & 3 = 2
base 2: 0100  & 0011 = 0100
dp[4] = 1 + dp[4] =  1 + 0 = 1

result 
dp = [0, 1, 1, 2, 1]

the pattern from these examples is 
dp[i]  = 1 + dp[i & (i - 1)]
```
##### Complexity Analysis:
```
Time  : O(N)
Space : O(N)
```
```python
from typing import List

def countBits(n: int) -> List[int]:
    if num < 0: return []
    dp = [0]*(num+1)
    for i in range(1, num+1):
        # Use kernighans algorithm for bit manipulation techinuque to turn off (set 1 to 0) the least significant bit (rightmost) in a number.
        # i & i - 1,  yeilds a number m , where m <= i, and has one less set bit than i.
        dp[i] = dp[i & (i-1)] + 1
    return dp

if __name__ == "__main__":
    #Input: n = 2
    #Output: [0,1,1]
    #Explanation:
    #0 --> 0
    #1 --> 1
    #2 --> 10
    n = 2
    print(countBits(n))
```
```kotlin
fun countBits(n: Int): IntArray {
    if (n < 0) return intArrayOf()
    val dp = IntArray(n+1)
    dp[0] = 0
    for (i in 1 until n+1) {
        dp[i] = dp[i and (i-1)]+1;
    }
    return dp
}

fun main(args: Array<String>) {
    //Input: n = 2
    //Output: [0,1,1]
    //Explanation:
    //0 --> 0
    //1 --> 1
    //2 --> 10
    val n = 2
    println(countBits(n).joinToString(","))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-268:Missing Number](https://leetcode.com/problems/missing-number/)
##### Solution Explanation:
```
Bitwise XOR Operation
=================================================================================================================================================================

- The basic idea is to use XOR operation.
- We all know that a^b^b =a, which means two xor operations with the same number will eliminate the number and reveal the original number.
- In this solution, I apply XOR operation to both the index and value of the array. 
- In a complete array with no missing numbers, the index and value should be perfectly corresponding( nums[index] = index), 
  so in a missing array, what left finally is the missing number.
```
##### Complexity Analysis:
```
Time  : O(N)
Space : O(1)
```
```python
from typing import List

def missingNumber(nums: List[int]) -> int:
    missing_number = len(nums)
    for i in range(len(nums)):
        missing_number ^= nums[i] ^ i
    
    return missing_number

if __name__ == "__main__":
    #Input: nums = [3,0,1]
    #Output: 2
    #Explanation: n = 3 since there are 3 numbers, so all numbers are in the range [0,3].
    #2 is the missing number in the range since it does not appear in nums.
    nums = [3,0,1]
    print(missingNumber(nums))
```
```kotlin
fun missingNumber(nums: IntArray): Int {
    var result = nums.size
    for (i in 0 until nums.size) {
        //result = result xor nums[i]
        //result = result xor i
        result = result.xor(nums[i]).xor(i)
    }
    return result    
}

fun main(args: Array<String>) {
    //Input: nums = [3,0,1]
    //Output: 2
    //Explanation: n = 3 since there are 3 numbers, so all numbers are in the range [0,3].
    //2 is the missing number in the range since it does not appear in nums.
    var nums = intArrayOf(3,0,1)
    println(missingNumber(nums))
    //Input: nums = [0,1]
    //Output: 2
    //Explanation: n = 2 since there are 2 numbers, so all numbers are in the range [0,2].
    //2 is the missing number in the range since it does not appear in nums.
    nums = intArrayOf(0,1)
    println(missingNumber(nums))
    //Input: nums = [9,6,4,2,3,5,7,0,1]
    //Output: 8
    //Explanation: n = 9 since there are 9 numbers, so all numbers are in the range [0,9].
    //8 is the missing number in the range since it does not appear in nums.
    nums = intArrayOf(9,6,4,2,3,5,7,0,1)
    println(missingNumber(nums))
    //Input: nums = [0]
    //Output: 1
    //Explanation: n = 1 since there is 1 number, so all numbers are in the range [0,1].
    //1 is the missing number in the range since it does not appear in nums.
    nums = intArrayOf(0)
    println(missingNumber(nums))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-190:Reverse Bits](https://leetcode.com/problems/reverse-bits/)
##### Solution Explanation:
```
Bitwise XOR Operation
=================================================================================================================================================================
Approach 1: Bit Manipulation and Bit-wise XOR operation
=================================================================================================================================================================

- In each loop, use logical AND operation n & 1 to get the least significant bit and add it to the ans. 
- To reverse the bit, shift n and ans in opposite directions.
- What needs attention is that, after we add the last bit of n to ans at the 32nd loop,
  the following left shift of ans is no longer needed.
- The most significant bit of n has already at the right most position after previous 31 loops.
=================================================================================================================================================================

=================================================================================================================================================================
Approach 2: Using Masking
=================================================================================================================================================================
For 8 bit binary number A B C D E F G H, the process is like: A B C D E F G H --> E F G H A B C D --> G H E F C D A B --> H G F E D C B A
```
##### Complexity Analysis:
```
For both approaches
=================================================================================================================================================================
Time  : O(log(N))
Space : O(1)
```
```python
# Approach 1: Bitwise XOR operation
def reverseBits(n: int) -> int:
    result = 0
    for i in range(32):
        result <<= 1
        result |= n & 1
        n >>= 1
    return result

if __name__ == "__main__":
    #Input       : n = 00000010100101000001111010011100
    #Output      : 964176192 (00111001011110000010100101000000)
    #Explanation : The input binary string 00000010100101000001111010011100 represents the unsigned integer 43261596,
    #so return 964176192 which its binary representation is 00111001011110000010100101000000.
    nums = 0b00000010100101000001111010011100
    print(reverseBits(nums))

# Approach 2: Using Masking
def reverseBits(n: int) -> int:
    n = (n >> 16) | (n << 16)
    n = ((n & 0xff00ff00) >> 8) | ((n & 0x00ff00ff) << 8)
    n = ((n & 0xf0f0f0f0) >> 4) | ((n & 0x0f0f0f0f) << 4)
    n = ((n & 0xcccccccc) >> 2) | ((n & 0x33333333) << 2)
    n = ((n & 0xaaaaaaaa) >> 1) | ((n & 0x55555555) << 1)
    return n

if __name__ == "__main__":
    #Input       : n = 00000010100101000001111010011100
    #Output      : 964176192 (00111001011110000010100101000000)
    #Explanation : The input binary string 00000010100101000001111010011100 represents the unsigned integer 43261596,
    #so return 964176192 which its binary representation is 00111001011110000010100101000000.
    nums = 0b00000010100101000001111010011100
    print(reverseBits(nums))
```
```kotlin
// Approach 1: Bitwise XOR operation
fun reverseBits(n:Int):Int {
    var num = n
    var result = 0

    for (i in 0 until 32) {
        result = result.shl(1)
        result += (num and 1)
        num = num ushr 1
    }

    return result
}

fun main(args: Array<String>) {
    //Input       : n = 00000010100101000001111010011100
    //Output      : 964176192 (00111001011110000010100101000000)
    //Explanation : The input binary string 00000010100101000001111010011100 represents the unsigned integer 43261596,
    //so return 964176192 which its binary representation is 00111001011110000010100101000000.
    val nums = 0b00000010100101000001111010011100
    println(reverseBits(nums))
}

// Approach 2: Using Masking
fun reverseBits(n:Int):Int {
    var num = n
    num = (num ushr 16) or (num shl 16);
    num = ((num and 0xFF00FF00.toInt()) ushr 8) or ((num and 0x00FF00FF.toInt()) shl 8)
    num = ((num and 0xF0F0F0F0.toInt()) ushr 4) or ((num and 0x0F0F0F0F.toInt()) shl 4)
    num = ((num and 0xCCCCCCCC.toInt()) ushr 2) or ((num and 0x33333333.toInt()) shl 2)
    num = ((num and 0xAAAAAAAA.toInt()) ushr 1) or ((num and 0x55555555.toInt()) shl 1)
    return num
}

fun main(args: Array<String>) {
    //Input       : n = 00000010100101000001111010011100
    //Output      : 964176192 (00111001011110000010100101000000)
    //Explanation : The input binary string 00000010100101000001111010011100 represents the unsigned integer 43261596,
    //so return 964176192 which its binary representation is 00111001011110000010100101000000.
    val nums = 0b00000010100101000001111010011100
    println(reverseBits(nums))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

## Dynamic Programming
| #     | Title	                                              | url                                                                           | Time       | Space   | Difficulty | Tag	         | Note                       |
| ----- | --------------------------------------------------- | ----------------------------------------------------------------------------- | ---------- | ------- | ---------- | ------------ | -------------------------- |
| 0070  | [Climbing Stairs](#lc-70climbing-stairs)            | https://leetcode.com/problems/climbing-stairs/                                | _O(logn)_  | _O(1)_  | Easy       |              | DP and Fibonacci Sequence  |
| 0322  | [Coin Change](#lc-322coin-change)                   | https://leetcode.com/problems/coin-change/                                    | _O(n * k)_ | _O(k)_  | Medium     |              |                            |
| 0300  | [Longest Increasing Subsequence (LIS)](#lc-300longest-increasing-subsequence) | https://leetcode.com/problems/longest-increasing-subsequence/                 | _O(nlogn)_ | _O(n)_  | Medium     | CTCI, LintCode | Patience Sorting using Binary Search, DP |
| 1143  | [Longest Common Subsequence (LCS)](#lc-1143longest-common-subsequence) | https://leetcode.com/problems/longest-common-subsequence/                     | _O(m * n)_ | _O(min(m, n))_ | Medium |           |                            |
| 0139  | [Word Break](#lc-139word-break)                     | https://leetcode.com/problems/word-break/                                     | _O(n * l^2)_ | _O(n)_ | Medium    |              |                            |
| 0377  | [Combination Sum IV](#lc-377combination-sum-iv)     | https://leetcode.com/problems/combination-sum-iv/                             | _O(nlogn + n * t)_ | _O(t)_ | Medium |           |                            |
| 0198  | [House Robber](#lc-198house-robber)                 | https://leetcode.com/problems/house-robber/                                   | _O(n)_     | _O(1)_  | Easy       |              |                            |
| 0213  | [House Robber II](#lc-213house-robber-ii)           | https://leetcode.com/problems/house-robber-ii/                                | _O(n)_     | _O(1)_  | Medium     |              |                            |
| 0091  | [Decode Ways](#lc-91decode-ways)                    | https://leetcode.com/problems/decode-ways/                                    | _O(n)_     | _O(1)_  | Medium     |              |                            |
| 0062  | [Unique Paths](#lc-62unique-paths)                  | https://leetcode.com/problems/unique-paths/                                   | _O(m * n)_ | _O(m + n)_ | Medium  |              |                            |
| 0055  | [Jump Game](#lc-55jump-game)                        | https://leetcode.com/problems/jump-game/                                      | _O(n)_     | _O(1)_  | Medium     |              |                            |

#### [LC-70:Climbing Stairs](https://leetcode.com/problems/climbing-stairs/)
##### Solution Explanation:
```
Intuition
---------
Solution to this problem makes a Fibonacci sequence. We can understand it better if we start from the end. 
To reach to Step N, you can either reach to step N-1 and take 1 step from there or take 2 step from N - 2.
Therefore it can be summarized as:
F(N) = F(N-1) + F(N-2)

Once you have recognized the pattern, it is very easy to write the code:
---------

Solution Approach:
DP and Fibonacci Sequence
=================================================================================================================================================================
To reach a specific stair x, we can either climb 1 stair from x-1, or 2 stairs from x-2. 
Therefore, suppose dp[i] records the number of ways to reach stair i, dp[i] = dp[i-1]+dp[i-2]. 
And it's a Fibonacci Array.
The base case is to reach the first stair, we only have one way to do it so dp[1] = 1.

Besides, since only dp elements we used is most recent two elements, we can use two pointers to save using of dp array. 
So space complexity is O(1).
```
##### Complexity Analysis:
```
Time  : O(N)
Space : O(1)
```
```python
def climbStairs(n: int) -> int:
    prev, curr = 0, 1
    for _ in range(n):
        curr = prev + curr
        prev = curr
    return curr

if __name__ == "__main__":
    #Input: n = 2
    #Output: 2
    #Explanation: There are two ways to climb to the top.
    #1. 1 step + 1 step
    #2. 2 steps
    n = 2
    print(climbStairs(n))
```
```kotlin
fun climbStairs(n: Int): Int {
    var prev = 0
    var curr = 1
    for (i in 0 until n) {
        curr = prev + curr
        prev = curr
    }
    return curr        
}

fun main(args: Array<String>) {
    //Input: n = 2
    //Output: 2
    //Explanation: There are two ways to climb to the top.
    //1. 1 step + 1 step
    //2. 2 steps
    val n = 2
    println(climbStairs(n))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-322:Coin Change](https://leetcode.com/problems/coin-change/)
##### Solution Explanation:
```
Solution Approach:
DP
=================================================================================================================================================================

Assume dp[i] is the fewest number of coins making up amount i, then for every coin in coins, dp[i] = min(dp[i - coin] + 1).

The time complexity is O(amount * coins.length) and the space complexity is O(amount).
```
##### Complexity Analysis:
```
TIME COMPLEXITY : O(amount * coins.length)
SPACE COMPLEXITY : O(amount)
```
```python
from typing import List

def coinChange(coins: List[int], amount: int) -> int:
    dp = [float('Inf')]*(amount+1)
    dp[0] = 0
    for i in range(1, amount+1):
        for coin in coins:
            if i - coin >= 0:
                dp[i] = min(dp[i], dp[i-coin] + 1)
    return dp[amount] if dp[amount] != float('Inf') else -1

if __name__ == "__main__":
    #Input: coins = [1,2,5], amount = 11
    #Output: 3
    #Explanation: 11 = 5 + 5 + 1
    coins = [1,2,5]
    amount = 11
    print(coinChange(coins, amount))
```
```kotlin
fun coinChange(coins: IntArray, amount: Int): Int {
    val dp = IntArray(amount+1)
    for (i in 1 until amount+1) {
        dp[i] = Int.MAX_VALUE
        for (coin in coins) {
            if (i - coin >= 0 && dp[i - coin] != Int.MAX_VALUE) {
                dp[i] = minOf(dp[i], dp[i - coin] + 1)
            }
        }
    }
    return if (dp[amount] == Int.MAX_VALUE) -1 else dp[amount]
}

fun main(args: Array<String>) {
    //Input: coins = [1,2,5], amount = 11
    //Output: 3
    //Explanation: 11 = 5 + 5 + 1
    val coins = intArrayOf(1,2,5)
    val amount = 11
    println(coinChange(coins, amount))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-300:Longest Increasing Subsequence](https://leetcode.com/problems/longest-increasing-subsequence/)
##### Solution Explanation:
```
=================================================================================================================================================================
Approach 1: Patience Sorting using Binary Search
=================================================================================================================================================================
This algorithm is actually Patience sorting. 
It might be easier for you to understand how it works if you think about it as piles of cards instead of tails.
The number of piles is the length of the longest subsequence.
For more info see Princeton lecture.

1) Initially, there are no piles. The first card dealt forms a new pile consisting of the single card.
2) Each subsequent card is placed on the leftmost existing pile whose top card has a value greater 
   than or equal to the new card's value, or to the right of all of the existing piles, thus forming a new pile.
3) When there are no more cards remaining to deal, the game ends.

Detailed Algorithm:
-----------------------
piles is an array storing the smallest tail of all increasing subsequences with length i+1 in piles[i].
For example, say we have nums = [4,5,6,3], then all the available increasing subsequences are:

len = 1   :      [4], [5], [6], [3]   => piles[0] = 3
len = 2   :      [4, 5], [5, 6]       => piles[1] = 5
len = 3   :      [4, 5, 6]            => piles[2] = 6
We can easily prove that piles is a increasing array. Therefore it is possible to do a binary search in piles array to find the one needs update.

Each time we only do one of the two:

(1) if num is larger than all piles, append it, increase the size by 1
(2) if piles[i-1] < num <= piles[i], update piles[i]
(3) Doing so will maintain the piles invariant. The the final answer is just the size.


-----------------------
References:
https://en.wikipedia.org/wiki/Patience_sorting
Priceton Lecture on LIS: https://www.cs.princeton.edu/courses/archive/spring13/cos423/lectures/LongestIncreasingSubsequence.pdf

=================================================================================================================================================================
Approach 2: DP
=================================================================================================================================================================
1) Check the base case, if nums has size less than or equal to 1, then return length of nums
2) Create a 'dp' array of size nums.length to track the longest sequence length
3) Fill each position with value 1 in the array
4) Mark one pointer at i. For each i, start from j=0.
   4.1.1) If, nums[j] < nums[i], it means next number contributes to increasing sequence. 
      4.1.1.1) But increase the value only if it results in a larger value of the sequence than dp[i].
               It is possible that dp[i] already has larger value from some previous j'th iteration.
5) Find the maximum length from the array that we just generated.

-----------------------
References:
https://www.youtube.com/watch?v=CE2b_-XfVDk
```
##### Interview Tips:
```
NOTES: In an Interview Situation - Choose Approach 1 ( Patience Sorting using Binary Search ) over Approach 2 ( DP ),
       since it is an O(N*log(N)) TC [ DP is worse off at O(N^2) ].
```
##### Complexity Analysis:
```
=================================================================================================================================================================
Approach 1: Patience Sorting using Binary Search
=================================================================================================================================================================
TIME COMPLEXITY : O(N*log(N))
SPACE COMPLEXITY : O(N)
=================================================================================================================================================================
Approach 2: DP
=================================================================================================================================================================
TIME COMPLEXITY : O(N^2)
SPACE COMPLEXITY : O(N)
```
```python
#=================================================================================================================================================================
#Approach 1 ( Patience Sorting using Binary Search ) ... answer for the Follow-Up question
#TC: O(N*log(N))
#SC: O(N)
#=================================================================================================================================================================
from typing import List

def binary_search(nums: List[int], target: int) -> int:
    lo, hi = 0, len(nums)
    while lo < hi:
        mid = (lo + hi) // 2
        if nums[mid] < target:
            lo = mid + 1
        else:
            hi = mid
    return lo

def lengthOfLIS(nums: List[int]) -> int:
    piles = []
    for num in nums:
        if not piles or num > piles[-1]:
            piles.append(num)
        else:
            pos = binary_search(piles, num)
            piles[pos] = num
    return len(piles)

def lengthOfLIS(nums: List[int]) -> int:
    if len(nums) == 0:
        return 0
    piles = [0] * len(nums)
    longest = 0
    for num in nums:
        i, j = 0, longest
        while i != j:
            m = (i + j) // 2
            if piles[m] < x:
                i = m + 1
            else:
                j = m
        piles[i] = x
        longest = max(i + 1, longest)
    return longest

if __name__ == "__main__":
    #Input: nums = [10,9,2,5,3,7,101,18]
    #Output: 4
    #Explanation: The longest increasing subsequence is [2,3,7,101], therefore the length is 4.
    nums = [10,9,2,5,3,7,101,18]
    print(lengthOfLIS(nums))
	
#=================================================================================================================================================================
#Approach 2 ( DP )
#TC: O(N^2)
#SC: O(N)
#=================================================================================================================================================================
from typing import List

def lengthOfLIS(nums: List[int]) -> int:
    if len(nums) == 0:
        return 0
    nums_length = len(nums)
    dp = [1] * nums_length
    longest = 1
    for i in range(nums_length):
        for j in range(i):
            if nums[i] < nums[j]:
                dp[i] = max(dp[i], dp[j]+1)
        longest = max(longest, dp[i])
    return longest

if __name__ == "__main__":
    #Input: nums = [10,9,2,5,3,7,101,18]
    #Output: 4
    #Explanation: The longest increasing subsequence is [2,3,7,101], therefore the length is 4.
    nums = [10,9,2,5,3,7,101,18]
    print(lengthOfLIS(nums))
```
```kotlin
//=================================================================================================================================================================
//Approach 1 ( Patience Sorting using Binary Search ) ... answer for the Follow-Up question
//TC: O(N*log(N))
//SC: O(N)
//=================================================================================================================================================================
fun binarySearch(nums: IntArray, target: Int): Int {
    var lo = 0
    var hi = nums.size
    while (lo < hi) {
        val mid = (lo + hi) / 2
        if (nums[mid] < target) {
            lo = mid + 1
        } else {
            hi = mid
        }
    }
    return lo
}

fun lengthOfLIS(nums: IntArray): Int {
    // sanity check
    if (nums.isEmpty()) return 0

    var piles: MutableList<Int> = mutableListOf()
    for (num in nums) {
        if ( num > piles?.lastOrNull() ?: -1 ) {
            piles.add(num)
        } else {
            val pos = binarySearch(piles.toIntArray(), num)
            piles[pos] = num
        }
    }
    return piles.size
}

fun lengthOfLIS(nums: IntArray): Int {
    // sanity check
    if (nums.isEmpty()) return 0
    val piles = IntArray(nums.size)
    var longest = 0
    for (x in nums) {
        var i = 0
        var j = longest
        while (i != j) {
            var m = (i + j) / 2
            if (piles[m] < x) {
                i = m + 1
            } else {
                j = m
            }
        }
        piles[i] = x
        if (i == longest) ++longest
    }
    return longest
}

fun main(args: Array<String>) {
    //Input: nums = [10,9,2,5,3,7,101,18]
    //Output: 4
    //Explanation: The longest increasing subsequence is [2,3,7,101], therefore the length is 4.
    val nums = intArrayOf(10,9,2,5,3,7,101,18)
    println(lengthOfLIS(nums))
}

//=================================================================================================================================================================
//Approach 2 ( DP )
//TC: O(N^2)
//SC: O(N)
//=================================================================================================================================================================
fun lengthOfLIS(nums: IntArray): Int {
    // sanity check
    if (nums.isEmpty()) return 0

    val nums_length = nums.size

    val dp = IntArray(nums_length) { 1 }
    var longest = 1

    for (i in 0 until nums_length) {
        for (j in 0 until i) {
            if (nums[i] > nums[j]) {
                dp[i] = maxOf(dp[i], dp[j] + 1)
            }
        }
        longest = maxOf(longest, dp[i])
    }
    return longest
}

fun main(args: Array<String>) {
    //Input: nums = [10,9,2,5,3,7,101,18]
    //Output: 4
    //Explanation: The longest increasing subsequence is [2,3,7,101], therefore the length is 4.
    val nums = intArrayOf(10,9,2,5,3,7,101,18)
    println(lengthOfLIS(nums))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-1143:Longest Common Subsequence](https://leetcode.com/problems/longest-common-subsequence/)
##### Solution Explanation:
```
Solution Approach:
=================================================================================================================================================================
DP with Memoization and 1D array for "Space Optimization"
-----------------------------------------------------------

Find LCS;
Let X be “XMJYAUZ” and Y be “MZJAWXU”. The longest common subsequence between X and Y is “MJAU”. 
The following table shows the lengths of the longest common subsequences between prefixes of X and Y.
The ith row and jth column shows the length of the LCS between X_{1..i} and Y_{1..j}.

+-------+---+---+---+---+---+---+---+---+
|       | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|       +---+---+---+---+---+---+---+---+
|       | 0 | M | Z | J | A | W | X | U |
+---+---+---+---+---+---+---+---+---+---+
| 0 | 0 |*0*| 0 | 0 | 0 | 0 | 0 | 0 | 0 |
+---+---+---+---+---+---+---+---+---+---+
| 1 | X | 0 | 0 | 0 | 0 | 0 | 0 | 1 | 1 |
+---+---+---+---+---+---+---+---+---+---+
| 2 | M | 0 |*1*| 1 | 1 | 1 | 1 | 1 | 1 |
+---+---+---+---+---+---+---+---+---+---+
| 3 | J | 0 | 1 | 1 |*2*| 2 | 2 | 2 | 2 |
+---+---+---+---+---+---+---+---+---+---+
| 4 | Y | 0 | 1 | 1 | 2 | 2 | 2 | 2 | 2 |
+---+---+---+---+---+---+---+---+---+---+
| 5 | A | 0 | 1 | 1 | 2 |*3*| 3 | 3 | 3 |
+---+---+---+---+---+---+---+---+---+---+
| 6 | U | 0 | 1 | 1 | 2 | 3 | 3 | 3 |*4*|
+---+---+---+---+---+---+---+---+---+---+
| 7 | Z | 0 | 1 | 2 | 2 | 3 | 3 | 3 | 4 |
+---+---+---+---+---+---+---+---+---+---+

References:
-----------------
https://en.m.wikipedia.org/wiki/Longest_common_subsequence_problem
https://www.ics.uci.edu/~eppstein/161/960229.html
```
##### Complexity Analysis:
```
TIME COMPLEXITY : O(M*N)
SPACE COMPLEXITY : O(MIN(M,N))

where:
------
M = length of string text1
N = length of string text2
```
```python
#Q & A:
#---------------
#Q1: What is the difference between [[0] * m * n] and [[0] * m for _ in range(n)]? 
#    Why does the former update all the rows of that column when I try to update one particular cell ?
#A1: [[0] * m * n] creates n references to the exactly same list objet: [0] * m; 
#    In contrast: [[0] * m for _ in range(n)] creates n different list objects that have same value of [0] * m.
#
def longestCommonSubsequence(text1: str, text2: str) -> int:
    m, n = map(len, (text1, text2))
    if m < n:
        return self.longestCommonSubsequence(text2, text1)
    dp = [0] * (n + 1)
    for c in text1:
        prevRow, prevRowPrevCol = 0, 0
        for j, d in enumerate(text2):
            prevRow, prevRowPrevCol = dp[j + 1], prevRow
            dp[j + 1] = prevRowPrevCol + 1 if c == d else max(dp[j], prevRow)
    return dp[-1]

if __name__ == "__main__":
    #Input: text1 = "abcde", text2 = "ace" 
    #Output: 3  
    #Explanation: The longest common subsequence is "ace" and its length is 3.
    text1 = "abcde"
    text2 = "ace"
    print(longestCommonSubsequence(text1,text2))
```
```kotlin
fun longestCommonSubsequence(text1: String, text2: String): Int {
    val m = text1.length
    val n = text2.length
    if (m < n) {
        return longestCommonSubsequence(text2, text1);
    }
    val dp = IntArray(n+1)
    for (c in text1) {
        var prevRow = 0
        var prevRowPrevCol = 0
        for ((j, d) in text2.withIndex()) {
            prevRow = dp[j + 1]
            prevRowPrevCol = prevRow
            dp[j + 1] = if (c.equals(d)) prevRowPrevCol + 1 else maxOf(dp[j], prevRow)
        }
    }
    return dp.last()
}

fun main(args: Array<String>) {
    //Input: text1 = "abcde", text2 = "ace" 
    //Output: 3  
    //Explanation: The longest common subsequence is "ace" and its length is 3.
    val text1 = "abcde"
    val text2 = "ace"
    println(longestCommonSubsequence(text1,text2))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-139:Word Break](https://leetcode.com/problems/word-break/)
##### Solution Explanation:
```
Solution Approach:
DP
=================================================================================================================================================================
(1) dp[hi] = s[0:hi] is segmentable ( or breakable ).
(2) Considering all possible substrings of s.
(3) If s[0:lo] is segmentable and s[lo:hi] is segmentable, then s[0:hi] is segmentable. Equivalently, if dp[lo] is True and s[lo:hi] is in the wordDict, then dp[hi] is True.
(4) Our goal is to determine if dp[hi] is segmentable, and once we do, we don't need to consider anything else. This is because we want to construct dp.
(5) dp[len(s)] tells us if s[0:len(s)] (or equivalently, s) is segmentable.
```
##### Complexity Analysis:
```
TIME COMPLEXITY  : O(N*L) [ O(N*L*N) or O(L*N^2) is substr is considered ]. 
SPACE COMPLEXITY : O(N)

where:
------
L = size of wordDict
N = length of string s
```
```python
#(1) dp[i] = s[0:i] is breakable
#(2) Considering all possible substrings of s.
#(3) If s[0:j] is breakable and s[j:i] is breakable, then s[0:i] is breakable.
#    Equivalently, if dp[j] is True and s[j:i] is in the wordDict, then dp[i] is True.
#(4) Our goal is to determine if dp[i] is breakable, and once we do, we don't need to consider anything else. This is because we want to construct dp.
#(5) dp[len(s)] tells us if s[0:len(s)] (or equivalently, s) is breakable.
def wordBreak(s: str, wordDict: List[str]) -> bool:
    if not s: return False
    dp = [False for i in range(len(s) + 1)] #(1)
    dp[0] = True
    
    for hi in range(len(s) + 1): #(2)
        for lo in range(hi):
            if dp[lo] and s[lo:hi] in wordDict: #(3)
                dp[hi] = True
                break #(4)
        
    return dp[len(s)] #(5)

if __name__ == "__main__":
    #Input: s = "leetcode", wordDict = ["leet","code"]
    #Output: true
    #Explanation: Return true because "leetcode" can be segmented as "leet code".
    s = "leetcode"
    wordDict = ["leet","code"]
    print(wordBreak(s, wordDict))
```
```kotlin
//(1) dp[i] = s[0:i] is breakable
//(2) Considering all possible substrings of s.
//(3) If s[0:j] is breakable and s[j:i] is breakable, then s[0:i] is breakable.
//     Equivalently, if dp[j] is True and s[j:i] is in the wordDict, then dp[i] is True.
//(4) Our goal is to determine if dp[i] is breakable, and once we do, we don't need to consider anything else. This is because we want to construct dp.
//(5) dp[len(s)] tells us if s[0:len(s)] (or equivalently, s) is breakable.
fun wordBreak(s: String, wordDict: List<String>): Boolean {
    // sanity check
    if(s.isEmpty()) return false

    val len = s.length
    val wordSet = HashSet(wordDict)

    val dp = BooleanArray(len + 1) //(1)
    dp[0] = true

    for (hi in 1..len) { //(2)
        for (lo in 0..hi) {
            if (dp[lo] && wordSet.contains(s.substring(lo, hi))) { //(3)
                dp[hi] = true
                break //(4)
            }
        }
    }
    
    return dp[len] //(5)
}

fun main(args: Array<String>) {
    //Input: s = "leetcode", wordDict = ["leet","code"]
    //Output: true
    //Explanation: Return true because "leetcode" can be segmented as "leet code".
    val s = "leetcode"
    val wordDict: List<String> = arrayListOf("leet","code")
    print(wordBreak(s, wordDict))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-377:Combination Sum IV](https://leetcode.com/problems/combination-sum-iv/)
##### Solution Explanation:
```
Solution Approach:
DP
=================================================================================================================================================================

Idea:
With this problem, we can easily imagine breaking up the solution into smaller pieces that we can use as stepping stones 
towards the overall answer. For example, if we're searching for a way to get from 0 to our target number (T), 
and if 0 < x < y < T, then we can see that finding out how many ways we can get 
from y to T will help us figure out how many ways we can get from x to T, 
all the way down to 0 to T. 

This is a classic example of a top-down (memoization) dyanamic programming (DP) solution.

Of course, the reverse is also true, and we could instead choose to use a bottom-up (tabulation) DP solution with the same result.

Top-Down DP Approach: 
---------------------

Our DP array (dp) will contain cells (dp[i]) where i will represent the remaining space left before T 
and dp[i] will represent the number of ways the solution (dp[T]) can be reached from i.
At each value of i as we build out dp we'll iterate through the different nums in our number array (N) 
and consider the cell that can be reached with each num (dp[i-num]). 
The value of dp[i] will therefore be the sum of the results of each of those possible moves.

We'll need to seed dp[0] with a value of 1 to represent the value of the completed combination, 
then once the iteration is complete, we can return dp[T] as our final answer.

Bottom-Up DP Approach: 
---------------------

Our DP array (dp) will contain cells (dp[i]) where i will represent the current count as we head towards T 
and dp[i] will represent the number of ways we can reach i from the starting point (dp[0]).
This means that dp[T] will represent our final solution.

At each value of i as we build out dp we'll iterate through the different nums in our number array (N)
and update the value of the cell that can be reached with each num (dp[i+num]) by adding the result 
of the current cell (dp[i]). If the current cell has no value, then we can continue without 
needing to iterate through N.

We'll need to seed dp[0] with a value of 1 to represent the value of the common starting point, 
then once the iteration is complete, we can return dp[T] as our final answer.
```
##### Complexity Analysis:
```
In both the top-down and bottom-up DP solutions, the time complexity is O(N * T) and the space complexity is O(T).


TIME COMPLEXITY  : O(N*T)
SPACE COMPLEXITY : O(T)
```
```python
# Top-Down DP
def combinationSum4_TopDownDP(nums: List[int], target: int) -> int:
    dp = [0] * (target + 1)
    dp[0] = 1
    for i in range(1, target+1):
        for num in nums:
            if num <= i: dp[i] += dp[i-num]
    return dp[target]

# Bottom-Up DP
def combinationSum4_BottomUpDP(nums: List[int], target: int) -> int:
    dp = [0] * (target + 1)
    dp[0] = 1
    for i in range(target):
        if not dp[i]: continue
        for num in nums:
            if num + i <= target: dp[i+num] += dp[i]
    return dp[target]

if __name__ == "__main__":
    #Input: nums = [1,2,3], target = 4
    #Output: 7
    #Explanation:
    #The possible combination ways are:
    #(1, 1, 1, 1)
    #(1, 1, 2)
    #(1, 2, 1)
    #(1, 3)
    #(2, 1, 1)
    #(2, 2)
    #(3, 1)
    #Note that different sequences are counted as different combinations.
    nums = [1,2,3]
    target = 4
    print(combinationSum4_TopDownDP(nums, target))
    print(combinationSum4_BottomUpDP(nums, target))
```
```kotlin
// Top-Down DP
fun combinationSum4_TopDownDP(nums: IntArray, target: Int): Int {
    val dp = IntArray(target + 1) { if (it == 0) 1 else 0 }
    for (i in 1..target) {
        for (num in nums) {
            if (num <= i) {
                dp[i] += dp.getOrNull(i - num) ?: 0
            }
        }
    }
    //return dp.last()
    return dp[target]
}

// Bottom-Up DP
fun combinationSum4_BottomUpDP(nums: IntArray, target: Int): Int {
    val dp = IntArray(target + 1) { if (it == 0) 1 else 0 }
    for (i in 0..target-1) {
        if (dp[i] == 0) continue
        for (num in nums) {
            if ((num + i) <= target) {
                dp[i+num] += dp.getOrNull(i) ?: 0
            }
        }
    }
    //return dp.last()
    return dp[target]
}


fun main(args: Array<String>) {
    //Input: nums = [1,2,3], target = 4
    //Output: 7
    //Explanation:
    //The possible combination ways are:
    //(1, 1, 1, 1)
    //(1, 1, 2)
    //(1, 2, 1)
    //(1, 3)
    //(2, 1, 1)
    //(2, 2)
    //(3, 1)
    //Note that different sequences are counted as different combinations.
    val nums = intArrayOf(1,2,3)
    val target = 4
    println(combinationSum4_TopDownDP(nums, target))
    println(combinationSum4_BottomUpDP(nums, target))
    //Input: nums = [9], target = 3
    //Output: 0
    val nums1 = intArrayOf(9)
    val target1 = 3
    println(combinationSum4_TopDownDP(nums1, target1))
    println(combinationSum4_BottomUpDP(nums1, target1))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-198:House Robber](https://leetcode.com/problems/house-robber/)
##### Solution Explanation:
```
Solution Approach:
-----------------------------------
Bottom-Up DP 2 ways:
-----------------------------------
1) Iterative + memo (bottom-up)
2) Iterative + N variables (bottom-up)

NOTE: For interview situation use method (2) => rob_iteratively_with_variables()
=================================================================================================================================================================

This particular problem and most of others can be approached using the following sequence:

  1. Find recursive relation
  2. Recursive (top-down)
  3. Recursive + memo (top-down)
  4. Iterative + memo (bottom-up)
  5. Iterative + N variables (bottom-up)

Step 1> Figure out recursive relation.
------------------------------------------------
A robber has 2 options: a) rob current house i; b) don't rob current house.
If an option "a" is selected it means she can't rob previous i-1 house but can safely proceed to the one before previous i-2 and gets all cumulative loot that follows.
If an option "b" is selected the robber gets all the possible loot from robbery of i-1 and all the following buildings.
So it boils down to calculating what is more profitable:

robbery of current house + loot from houses before the previous
loot from the previous house robbery and any loot captured before that


rob(i) = Math.max( rob(i - 2) + currentHouseValue, rob(i - 1) )

Step 2. Recursive (top-down)
------------------------------------------------

public int rob(int[] nums) {
    return rob(nums, nums.length - 1);
}
private int rob(int[] nums, int i) {
    if (i < 0) {
        return 0;
    }
    return Math.max(rob(nums, i - 2) + nums[i], rob(nums, i - 1));
}

Step 3. Recursive + memo (top-down).
------------------------------------------------

int[] memo;
public int rob(int[] nums) {
    memo = new int[nums.length + 1];
    Arrays.fill(memo, -1);
    return rob(nums, nums.length - 1);
}

private int rob(int[] nums, int i) {
    if (i < 0) {
        return 0;
    }
    if (memo[i] >= 0) {
        return memo[i];
    }
    int result = Math.max(rob(nums, i - 2) + nums[i], rob(nums, i - 1));
    memo[i] = result;
    return result;
}

Much better, this should run in O(n) time. Space complexity is O(n) as well, because of the recursion stack, let's try to get rid of it.


Step 4. Iterative + memo (bottom-up)
------------------------------------------------
public int rob(int[] nums) {
    if (nums.length == 0) return 0;
    int[] memo = new int[nums.length + 1];
    memo[0] = 0;
    memo[1] = nums[0];
    for (int i = 1; i < nums.length; i++) {
        int val = nums[i];
        memo[i+1] = Math.max(memo[i], memo[i-1] + val);
    }
    return memo[nums.length];
}


Step 5. Iterative + 2 variables (bottom-up)
------------------------------------------------
We can notice that in the previous step we use only memo[i] and memo[i-1], so going just 2 steps back. We can hold them in 2 variables instead.
This optimization is met in Fibonacci sequence creation and some other problems.

/* the order is: prev2, prev1, num  */
public int rob(int[] nums) {
    if (nums.length == 0) return 0;
    int prev1 = 0;
    int prev2 = 0;
    for (int num : nums) {
        int tmp = prev1;
        prev1 = Math.max(prev2 + num, prev1);
        prev2 = tmp;
    }
    return prev1;
}
```
##### Complexity Analysis:
```
For both the solutions mentioned in Steps 4 (rob_iteratively_using_memo) and 5 (rob_iteratively_with_variables)

TIME COMPLEXITY  : O(N)
SPACE COMPLEXITY : O(1)
```
```python
#Bottom-Up DP 2 ways:
#-----------------------------------
#1) Iterative + memo (bottom-up)
#2) Iterative + N variables (bottom-up)
#
#NOTE: For interview situation use method (2) => rob_iteratively_with_variables()
"""
Note to self:
rerun "pylint house_robber.py; python3 house_robber.py"

I use pyright to do static type checking in VSCode
"""

# pylint: disable = too-few-public-methods, no-self-use
from typing import List # so you can do List[int]

# https://leetcode.com/problems/house-robber/discuss/156523/From-good-to-great.-How-to-approach-most-of-DP-problems.
# my take-away:
# 1) recursive
# 2) recursive memo
# 3) iterative memo
# 4) iterative "pointer" variables

class Solution:
    """solution for 'House Robber' on leetcode"""

    def __init__(self):
        self.memo = {}

    def rob(self, houses: List[int]) -> int:
        """
        Gets the max loot from non-adjacent houses.
        Stores data in the given list of houses.
        """
        # handle trivial cases:
        if not houses:
            return 0
        if len(houses) == 1:
            return houses[0]
        if len(houses) == 2:
            return max(houses[0], houses[1])
        # return self.rob_recursively(houses, 0) # works
        # return self.rob_iteratively(houses) # better
        return self.rob_iteratively_with_variables(houses) # even better

    def rob_recursively(self, houses: List[int], i: int) -> int:
        """
        TO FIND A RECURSIVE SOLUTION: THINK OF THE BASE CASES!!!
        Either (1) loot this house (and then must skip next one),
        or (2) don't loot this house (and loot the next one).
        If you don't loot this house and not the next one either,
        then you're being silly and are better off with case (1) above anyways.
        (BTW: going left to right recursively until base case i >= len.)
        """
        if i >= len(houses):
            return 0
        if i in self.memo:
            return self.memo[i]
        loot_this_house_and_next_house = houses[i] + self.rob_recursively(houses, i + 2)
        loot_next_house = self.rob_recursively(houses, i + 1)
        max_loot = max(loot_this_house_and_next_house, loot_next_house)
        self.memo[i] = max_loot
        return max_loot

    def rob_iteratively_using_memo(self, houses: List[int]) -> int:
        """
        Iterative solution: loot up to the NEXT house =
        either (1) loot from this house + house 2 ago,
        or (2) loot previous house.
        (BTW: going left to right with O(n).)
        """
        self.memo[0] = houses[0]
        # loop starting at the 2nd house:
        for i in range(1, len(houses)):
            if i - 2 < 0: # as if only getting 2nd house (otherwise invalid index)
                loot_this_house_and_2_ago = houses[i]
            else: # otherwise i - 2 is a valid index
                loot_this_house_and_2_ago = houses[i] + self.memo[i - 2]
            loot_previous_house = self.memo[i - 1]
            self.memo[i] = max(loot_this_house_and_2_ago, loot_previous_house)
        return self.memo[len(houses) - 1]

    def rob_iteratively_with_variables(self, houses: List[int]) -> int:
        """
        (Iterative solution that uses "pointer" variables instead of memo.)
        Improve on the iterative memo solution by noticing:
            self.memo[i - 2] = houses[i - 2]
            self.memo[i - 1] = houses[i - 1]
        so:
            loot_this_house_and_2_ago = this house + two ago
            loot_previous_house = one ago
        (BTW: going left to right with O(n).)
        """
        prev = 0
        curr = 0
        for house in houses:
            # curr: current house = either previous house, or this house + two ago
            # prev: just moves one to the next position
            loot_this_house_and_2_ago = house + prev
            loot_previous_house = curr
            curr = max(loot_previous_house, loot_this_house_and_2_ago)
            prev = loot_previous_house
        return curr # = current house = either previous house, or this house + two ago

if __name__ == "__main__":
    def check_answer(houses, correct):
        """helper function"""
        answer = Solution().rob(houses)
        assert answer == correct, f'{answer} should be {correct}'
        print(answer, 'ok' if answer == correct else 'error')
    check_answer(houses=[], correct=0) # empty
    check_answer(houses=None, correct=0) # invalid
    check_answer(houses=[1], correct=1) # simple
    check_answer(houses=[111], correct=111) # simple
    check_answer(houses=[1, 2], correct=2)
    check_answer(houses=[1, 2, 3], correct=4)
    check_answer(houses=[1, 2, 3, 1], correct=4)
    check_answer(houses=[2, 7, 9, 3, 1], correct=12)
    check_answer(houses=[9, 1, 1, 9], correct=18)
    check_answer(houses=[0], correct=0)
    check_answer(houses=[1, 0, 0, 0], correct=1)
    check_answer(houses=[0, 1, 0, 0, 0], correct=1)
    check_answer(houses=[0, 1, 0, 0, 0], correct=1)
    check_answer(houses=[1, 0, 1, 0, 0, 1], correct=3)
    check_answer(houses=[0, 1, 0, 1, 0, 0, 1], correct=3)
    check_answer(houses=[155, 44, 52, 58, 250, 225, 109, 118, 211, \
        73, 137, 96, 137, 89, 174, 66, 134, 26, 25, 205, 239, 85, 146, \
        73, 55, 6, 122, 196, 128, 50, 61, 230, 94, 208, 46, 243, 105, \
        81, 157, 89, 205, 78, 249, 203, 238, 239, 217, 212, 241, 242, \
        157, 79, 133, 66, 36, 165], correct=4517) # requires fast algorithm
    check_answer(houses=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], correct=30)
    check_answer(houses=[10, 9, 8, 7, 6, 5, 4, 3, 2, 1], correct=30)
```
```kotlin
fun rob_iteratively_using_memo(houses: IntArray): Int {
    //
    //Iterative solution: loot up to the NEXT house =
    //either (1) loot from this house + house 2 ago,
    //or     (2) loot previous house.
    //(BTW: going left to right with O(n).)
    //
    val len = houses.size
    val memo = IntArray(houses.size)
    memo[0] = houses[0]
    var loot_this_house_and_2_ago = -1
    var loot_previous_house = -1 
    // loop starting at the 2nd house:
    for (i in 1 until len) { 
        if ((i - 2) < 0) { // as if only getting 2nd house (otherwise invalid index)
           loot_this_house_and_2_ago = houses[i]
        } else { // otherwise i - 2 is a valid index
            loot_this_house_and_2_ago = houses[i] + memo[i - 2]
            loot_previous_house = memo[i - 1]
            memo[i] = maxOf(loot_this_house_and_2_ago, loot_previous_house)
        }
    }
    return memo[len - 1]
}

fun rob_iteratively_with_variables(houses: IntArray): Int {
    //
    //(Iterative solution that uses "pointer" variables instead of memo.)
    //Improve on the iterative memo solution by noticing:
    //    self.memo[i - 2] = houses[i - 2]
    //    self.memo[i - 1] = houses[i - 1]
    //so:
    //    loot_this_house_and_2_ago = this house + two ago
    //    loot_previous_house = one ago
    //(BTW: going left to right with O(n).)
    //
    var prev = 0
    var curr = 0
    for (house in houses) {
        // curr: current house = either previous house, or this house + two ago
        // prev: just moves one to the next position
        val loot_this_house_and_2_ago = house + prev
        val loot_previous_house = curr
        curr = maxOf(loot_previous_house, loot_this_house_and_2_ago)
        prev = loot_previous_house
    }
    return curr // = current house = either previous house, or this house + two ago
}

fun main(args: Array<String>) {
    //Input: houses = [1,2,3,1]
    //Output: 4
    //Explanation: Rob house 1 (money = 1) and then rob house 3 (money = 3).
    //Total amount you can rob = 1 + 3 = 4.
    val houses = intArrayOf(1,2,3,1)
    println(rob_iteratively_using_memo(houses))
    println(rob_iteratively_with_variables(houses))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-213:House Robber II](https://leetcode.com/problems/house-robber-ii/)
##### Solution Explanation:
```
Solution Approach:
-----------------------------------
Bottom-Up DP => Iterative + N variables (bottom-up)

Variant of [ House Robber | LeetCode Problem 198 | https://leetcode.com/problems/house-robber/ ] ... can be solved by just calling the solution for LC-198 twice.
=================================================================================================================================================================
This problem can be seen as follow-up question for problem 198. House Robber. 

If thief choose to rob first house, then thief cannot rob last house. Simiarly, if choose to rob last house, then cannot rob first.
So if we are given houses [1,3,4,5,6,7], then we are taking max from:
1.[1,3,4,5,6] OR
2.[3,4,5,6,7]

whichever's max value is larger.
So we just do 2 pass of dp and take the larger one. And the problem is reduced to House Robber I [ House Robber | LeetCode Problem 198 | https://leetcode.com/problems/house-robber/ ].
```
##### Complexity Analysis:
```
Time Complexity: time complexity is O(n), because we use dp problem with complexity O(n) twice. 
Space complexity is O(1), because in python lists passed by reference and space complexity of House Robber problem is O(1).

TIME COMPLEXITY  : O(N)
SPACE COMPLEXITY : O(1)
```
```python
from typing import List

class Solution:

    def __init__(self):
        self.memo = {}

    def rob_dp(self, houses: List[int], start: int, end: int) -> int:
        prev = 0
        curr = 0
        for i in range(start,end+1):
            loot_this_house_and_2_ago = houses[i] + prev
            loot_previous_house = curr
            curr = max(loot_previous_house, loot_this_house_and_2_ago)
            prev = loot_previous_house
        return curr

    def rob(self, houses: List[int]) -> int:
        return max(self.rob_dp(houses, 0, len(houses)-2), self.rob_dp(houses, 1, len(houses)-1))

if __name__ == "__main__":
    #Input: nums = [2,3,2]
    #Output: 3
    #Explanation: You cannot rob house 1 (money = 2) and then rob house 3 (money = 2), because they are adjacent houses.
    houses = [2,3,2]
    solution = Solution()
    print(solution.rob(houses))
```
```kotlin
fun rob_dp(houses: IntArray, start: Int, end: Int): Int {
    var prev = 0
    var curr = 0
    for (i in start..end) {
        val loot_this_house_and_2_ago = houses[i] + prev
        val loot_previous_house = curr
        curr = maxOf(loot_previous_house, loot_this_house_and_2_ago)
        prev = loot_previous_house
    }
    return curr
}

fun rob(houses: IntArray): Int {
    return maxOf(rob_dp(houses, 0, houses.size-2), rob_dp(houses, 1, houses.size-1))
}

fun main(args: Array<String>) {
    //Input: nums = [2,3,2]
    //Output: 3
    //Explanation: You cannot rob house 1 (money = 2) and then rob house 3 (money = 2), because they are adjacent houses.
    val houses = intArrayOf(2,3,2)
    println(rob(houses))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-91:Decode Ways](https://leetcode.com/problems/decode-ways/)
##### Solution Explanation:
```
DP Approach 1>
------------------
Use a dp array of size n + 1 to save subproblem solutions.

dp[0] means an empty string will have one way to decode,
dp[1] means the way to decode a string of size 1.

Check one digit and two digit combination and save the results along the way.

In the end, dp[n] will be the end result.

For example:
s = "231"
index 0: extra base offset. dp[0] = 1
index 1: # of ways to parse "2" => dp[1] = 1
index 2: # of ways to parse "23" => "2" and "23", dp[2] = 2
index 3: # of ways to parse "231" => "2 3 1" and "23 1" => dp[3] = 2

DP Approach 2> ( SPACE OPTIMIZATION - Constant Space )
------------------
We can use two variables to store the previous results.
Since we only use dp[i-1] and dp[i-2] to compute dp[i]. 
Why not just use two variable prev1, prev2 instead?
This can reduce the space to O(1)
```
##### Complexity Analysis:
```
For DP Approach 1>

TIME COMPLEXITY  : O(N)
SPACE COMPLEXITY : O(N)

For DP Approach 2>

TIME COMPLEXITY  : O(N)
SPACE COMPLEXITY : O(1)
```
```python
#DP Approach 1>
#------------------
def numDecodings(s: str) -> int:
    if not s or s[0] == '0':
        return 0

    dp = [0 for x in range(len(s) + 1)]
    dp[0] = 1
    dp[1] = 1 if 0 < int(s[0]) <= 9 else 0

    for i in range(2, len(s) + 1):
        first = int(s[i-1:i])
        second = int(s[i-2:i])
        if 1 <= first <= 9:
            dp[i] += dp[i - 1]
        if 10 <= second <= 26:
            dp[i] += dp[i - 2]
    return dp[len(s)]

if __name__ == "__main__":
    #Input: s = "12"
    #Output: 2
    #Explanation: "12" could be decoded as "AB" (1 2) or "L" (12).
    s = "12"
    print(numDecodings(s))

#DP Approach 2> ( SPACE OPTIMIZATION - Constant Space )
#------------------
def numDecodings(s: str) -> int:
    if not s or s[0] == '0':
        return 0

    # pre represents dp[i-1]
    pre = 1
    # ppre represents dp[i-2]
    ppre = 0
    for i in range(1, len(s) + 1):
        temp = pre
        if s[i - 1] == "0":
            pre = 0
        if i > 1 and 10 <= int(s[i-2:i]) <= 26:
            pre += ppre
        ppre = temp
    return pre

if __name__ == "__main__":
    #Input: s = "12"
    #Output: 2
    #Explanation: "12" could be decoded as "AB" (1 2) or "L" (12).
    s = "12"
    print(numDecodings(s))
```
```kotlin
//DP Approach 1>
//------------------
fun numDecodings(s: String): Int {
    if (s.isEmpty() || s[0] == '0') return 0

    val n = s.length
    val dp = IntArray(n + 1)
        
    /*
     * dp[0] is set to 1 only to get the result for dp[2].
     * For example, you have a string "12" , "12" could be decoded as "AB" (1 2) or "L" (12).
     * Now if you select "12" , then dp[2] += dp[0]. If dp[0] is 0, you wont count '12' as a way to decode. Hence dp[0]
     * needs to be 1.
     */
    dp[0] = 1 // To handle the case like "12"
    dp[1] = if ((s[0] > '0') and (s[0] <= '9')) 1 else 0
        
    for (i in 2..n) {
        val first = s.substring(i - 1, i).toInt()
        val second = s.substring(i - 2, i).toInt()
            
        if (first in 1..9) {
            dp[i] += dp[i - 1]
        }
        if (second in 10..26) {
            dp[i] += dp[i -2]
        }
    }
    return dp[n]        
}

fun main(args: Array<String>) {
    //Input: s = "12"
    //Output: 2
    //Explanation: "12" could be decoded as "AB" (1 2) or "L" (12).
    val s = "12"
    println(numDecodings(s))
}

//DP Approach 2> ( SPACE OPTIMIZATION - Constant Space )
//------------------
fun numDecodings(s: String): Int {
    if (s.isEmpty() || s[0] == '0') return 0
    val n = s.length
    // pre represents dp[i-1]
    var pre = 1
    // ppre represents dp[i-2]
    var ppre = 0
    for (i in 1..n) {
        val temp = pre
        if (s[i - 1] == '0') {
            pre = 0
        }
        if ( i > 1 ) {
            if (s.substring(i - 2, i).toInt() in 10..26) pre += ppre
        }
        ppre = temp        
    }

    return pre
}

fun main(args: Array<String>) {
    //Input: s = "12"
    //Output: 2
    //Explanation: "12" could be decoded as "AB" (1 2) or "L" (12).
    val s = "12"
    println(numDecodings(s))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-62:Unique Paths](https://leetcode.com/problems/unique-paths/)
##### Solution Explanation:
```
Solution Approach:
=================================================================================================================================================================
DP without Recursion + Space Optimization

- path[i,j] = Number of paths from [0,0] to [i,j].
- path[0,j] = 1 and path[i,0] = 1
- path[i,j] = path[i,j-1] + path[i-1,j]
- return path[m-1, n-1]
- We can start from row 1 and column 1 after initializing the path matrix to 1.

Time and Space complexity: O(MN)

- Space Optimization: Instead of 2D matrix, a single array can do the job and reduce space complexity to O(N)
```
##### Complexity Analysis:
```
TIME COMPLEXITY  : O(N)
SPACE COMPLEXITY : O(1)
```
```python
def uniquePaths(m: int, n: int) -> int:
    if m == 0 or n == 0:
        return 0
    dp = [1]*n
    for i in range(1,m):
        for j in range(1,n):
            dp[j] = dp[j-1] + dp[j]
    return dp[-1]

if __name__ == "__main__":
    #Input: m = 3, n = 2
    #Output: 3
    #Explanation:
    #From the top-left corner, there are a total of 3 ways to reach the bottom-right corner:
    #1. Right -> Down -> Down
    #2. Down -> Down -> Right
    #3. Down -> Right -> Down
    m = 3
    n = 2
    print(uniquePaths(m, n))
```
```kotlin
fun uniquePaths(m: Int, n: Int): Int {
    if ((m == 0) or (n == 0)) {
        return 0
    }

    var dp = IntArray(n){1}

    for (i in 1 until m) {
        for (j in 1 until n){
            dp[j] += dp[j-1]
        }
    }

    return dp[n-1]        
}

fun main(args: Array<String>) {
    //Input: m = 3, n = 2
    //Output: 3
    //Explanation:
    //From the top-left corner, there are a total of 3 ways to reach the bottom-right corner:
    //1. Right -> Down -> Down
    //2. Down -> Down -> Right
    //3. Down -> Right -> Down
    val m = 3
    val n = 2
    println(uniquePaths(m, n))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

#### [LC-55:Jump Game](https://leetcode.com/problems/jump-game/)
##### Solution Explanation:
```
Solution Approach:
=================================================================================================================================================================
Greedy Algorithm

Greedy -- Reach means last num's maximum reach position.
If current index > reach, then means can't reach current position from last number.

1. We start travering the array from start
2. While traversing, we keep a track on maximum reachable index and update it accordingly. 
3. If we cannot reach the maxium reachable index we get out of loop ( If current index > reach, then means can't reach current position from last number ).
```
##### Complexity Analysis:
```
TIME COMPLEXITY  : O(N)
SPACE COMPLEXITY : O(1)
```
```python
from typing import List

def canJump(nums: List[int]) -> bool:
    reachable_ind = 0
    for ind, val in enumerate(nums):
        if ind > reachable_ind:
            return False
        reachable_ind = max(reachable_ind, ind + val) 
            
    return True

if __name__ == "__main__":
    #Input: nums = [2,3,1,1,4]
    #Output: true
    #Explanation: Jump 1 step from index 0 to 1, then 3 steps to the last index.
    nums = [2,3,1,1,4]
    print(canJump(nums))
```
```kotlin
fun canJump(nums: IntArray): Boolean {
    var reachable_ind = 0
    for ((index, value) in nums.withIndex()) {
        if (index > reachable_ind) {
            return false
        }
        reachable_ind = maxOf(reachable_ind, index + value)
    }            
    return true
}

fun main(args: Array<String>) {
    //Input: nums = [2,3,1,1,4]
    //Output: true
    //Explanation: Jump 1 step from index 0 to 1, then 3 steps to the last index.
    val nums = intArrayOf(2,3,1,1,4)
    println(canJump(nums))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

## String
| #     | Title	                                         | url                                                                           | Time   | Space   | Difficulty | Tag	        | Note                   |
| ----- | ---------------------------------------------- | ----------------------------------------------------------------------------- | ------ | ------- | ---------- | ------------ | ---------------------- |
| 0003  | Longest Substring Without Repeating Characters | https://leetcode.com/problems/longest-substring-without-repeating-characters/ | _O(n)_ | _O(1)_  | Medium     |              |                        |
| 0424  | Longest Repeating Character Replacement        | https://leetcode.com/problems/longest-repeating-character-replacement/        | _O(n)_ | _O(1)_  | Medium     |              | Sliding Window         |
| 0076  | Minimum Window Substring                       | https://leetcode.com/problems/minimum-window-substring/                       | _O(n)_ | _O(k)_  | Hard       |              |                        |
| 0242  | Valid Anagram                                  | https://leetcode.com/problems/valid-anagram/                                  | _O(n)_ | _O(1)_  | Easy       |              |                        |
| 0049  | Group Anagrams                                 | https://leetcode.com/problems/group-anagrams/                                 | _O(n)_ | _O(1)_  | Easy       |              |                        |
| 0678  | Valid Parenthesis String                       | https://leetcode.com/problems/valid-parenthesis-string/                       | _O(n)_ | _O(1)_  | Medium     |              |                        |
| 0125  | Valid Palindrome                               | https://leetcode.com/problems/valid-palindrome/                               | _O(n)_ | _O(1)_  | Easy       |              |                        |
| 0005  | Longest Palindromic Substring                  | https://leetcode.com/problems/longest-palindromic-substring/                  | _O(n)_ | _O(n)_  | Medium     |              | `Manacher's Algorithm` |
| 0647  | Palindromic Substrings                         | https://leetcode.com/problems/palindromic-substrings/                         | _O(n)_ | _O(n)_  | Medium     |              | `Manacher's Algorithm` |
| 0271  | Encode and Decode Strings                      | https://leetcode.com/problems/encode-and-decode-strings/                      | _O(n)_ | _O(1)_  | Medium     | 🔒           |                        |

####  [LC-3:Longest Substring Without Repeating Characters](https://leetcode.com/problems/longest-substring-without-repeating-characters/)
##### Solution Explanation:
```
Approach: Sliding Window Algorithm.
=================================================================================================================================================================

From the input, we can gather the following information -

1. Given data structure is a string which is a linear data structure.
2. The output must be a substring - a part of the given string.
3. Naive solution is to check for each combination of characters in the string

Are you thinking what I am thinking 🤔 ? 
Yes, this is a classic example of a problem that can be solved using the legendary technique - Sliding Window Algorithm.

Following are the steps that we will follow -

1. Have two pointers which will define the starting index start and ending index end of the current window. Both will be 0 at the beginning.
2. Declare a Set that will store all the unique characters that we have encountered.
3. Declare a variable maxLength that will keep track of the length of the longest valid substring.
4. Scan the string from left to right one character at a time.
5. If the character has not encountered before i.e., not present in the Set the we will add it and increment the end index. The maxLength will be the maximum of Set.size() and existing maxLength.
6. If the character has encounter before, i.e., present in the Set, we will increment the start and we will remove the character at start index of the string.
7. Steps #5 and #6 are moving the window.
8. After the loop terminates, return maxLength.
```
##### Complexity Analysis:
```
a) Time  : O(N)
============================
We are scanning the string from left to right only once, hence the time complexity will be O(n).

b) Space : O(1)
============================
We are using Set as a data structure to facilitate our computation, therefore, the space complexity should also be O(n), right? Wrong!

WHY?

The problem clearly states that the string contains only English letters, digits, symbols and spaces and are covered in 256 code points.
Therefore, a string will be made up of a combination of these characters.

Since a Set can contain only unique elements, at any point the size of Set cannot be more than 256.

What does this mean? This means that the size of set is a function independent of the size of the input.
It is a constant.
Therefore, the space complexity will be O(1) (let me know in comments if you think otherwise).
```
```python
def lengthOfLongestSubstring(s: str) -> int:
    # Base condition
    if s == "":
        return 0
    # Starting index of window
    start = 0
    # Ending index of window
    end = 0
    # Maximum length of substring without repeating characters
    maxLength = 0
    # Set to store unique characters
    unique_characters = set()
    # Loop for each character in the string
    while end < len(s):
        if s[end] not in unique_characters:
            unique_characters.add(s[end])
            end += 1
            maxLength = max(maxLength, len(unique_characters))
        else:
            unique_characters.remove(s[start])
            start += 1
    return maxLength

if __name__ == "__main__":
    #Input: s = "abcabcbb"
    #Output: 3
    #Explanation: The answer is "abc", with the length of 3.
    s = "abcabcbb"
    print(lengthOfLongestSubstring(s))
```
```kotlin
fun lengthOfLongestSubstring(s: String): Int {
    // Base condition
    if (s == "") {
        return 0
    }
    // Starting window index
    var start = 0
    // Ending window index
    var end = 0
    // Maximum length of substring
    var maxLength = 0
    // This set will store the unique characters
    val uniqueCharacters: MutableSet<Char> = HashSet()
    // Loop for each character in the string
    while (end < s.length) {
        if (uniqueCharacters.add(s[end])) {
            end++
            maxLength = maxLength.coerceAtLeast(uniqueCharacters.size)
        } else {
            uniqueCharacters.remove(s[start])
            start++
        }
    }
    return maxLength
}

fun main(args: Array) {
    //Input: s = "abcabcbb"
    //Output: 3
    //Explanation: The answer is "abc", with the length of 3.
    val s = "abcabcbb"
    println(lengthOfLongestSubstring(s))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

####  [LC-424:Longest Repeating Character Replacement](https://leetcode.com/problems/longest-substring-without-repeating-characters/)
##### Solution Explanation:
```
Approach: Sliding Window Algorithm.
=================================================================================================================================================================
- The idea here is to find a window that satisfies the condition -
- count of most repeatable character + no. of allowed replacements <= length of the window
- Since the no. of allowed replacements is fixed, then the window size is directly proportional to the count of the most repeating character.
- Initially the window keeps growing from the end, until all the allowed replacements are added up in the window until it reaches the max size.
- The moment the condition is not satisfied (i.e., count of most repeatable character + no. of allowed replacements > size of the window), 
  then we need to slide the window (not shrink) 
  to the right and decrement the frequency of the character that is moved out of the window.
- If the next character coming in is the most repeating character, then the window grows or else it simply slides again.
```
##### Complexity Analysis:
```
Time  : O(N)
Space : O(N)
```
```python
import collections
def characterReplacement(s, k):        
    # Base condition
    if s == "":
        return 0
    longest_window = 0
    window_counts = collections.defaultdict(int)
    start = 0
    for end in range(len(s)):
        window_counts[s[end]] += 1
        while ( (end - start + 1) - max(window_counts.values()) ) > k:
            window_counts[s[start]] -= 1 
            start += 1
        longest_window = max(longest_window, end - start + 1)
    return longest_window

if __name__ == "__main__":
    #Input: s = "ABAB", k = 2
    #Output: 4
    s = "ABAB"
    k = 2
    print(characterReplacement(s,k))
```
```kotlin
fun characterReplacement(s: String, k: Int): Int {
    // Base condition
    if (s == "") {
        return 0
    }
    var mostFreqCharCount = 0; var start = 0; var max=0
    val map = mutableMapOf<Char, Int>()
        
    for (end in 0 until s.length){
        map.put(s[end], map.getOrDefault(s[end], 0) + 1)
        mostFreqCharCount = Math.max(map.get(s[end])!!, mostFreqCharCount)
        if ( ( (end - start + 1) - mostFreqCharCount ) > k ) {
            map.put(s[start], map.get(s[start])!! - 1)
            start++                
        }
        max = Math.max(max, end - start + 1)
    }
    return max
}

fun main(args: Array) {
    //Input: s = "ABAB", k = 2
    //Output: 4
    val s = "ABAB"
    val k = 2
    println(characterReplacement(s,k))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

####  [LC-76:Minimum Window Substring](https://leetcode.com/problems/minimum-window-substring/)
##### Solution Explanation:
```
Approach: Sliding Window Algorithm.
=================================================================================================================================================================
- The idea is we use a variable-length sliding window which is gradually applied across the string.
- We use two pointers: start and end to mark the sliding window.
- We start by fixing the start pointer and moving the end pointer to the right.
- The way we determine the current window is a valid one is by checking if all the target letters have been found in the current window.
- If we are in a valid sliding window, we first make note of the sliding window of the most minimum length we have seen so far.
- Next we try to contract the sliding window by moving the start pointer.
- If the sliding window continues to be valid, we note the new minimum sliding window. 
- If it becomes invalid (all letters of the target have been bypassed), 
  we break out of the inner loop and go back to moving the end pointer to the right.
```
##### Complexity Analysis:
```
Time  : O(N)
Space : O(N)
```
```python
import collections
def minWindow(s: str, t: str) -> str:
    # Base condition
    if (s == "" || t == "" || len(s) < len(t))  return ""
    char_freq_in_target = collections.Counter(target)
    start = 0
    end = 0
    shortest = ""
    target_len = len(target)
        
    for end in range(len(s)):
        # If we see a target letter, decrease the total target letter count
        if char_freq_in_target[s[end]] > 0:
            target_len -= 1

        # Decrease the letter count for the current letter
        # If the letter is not a target letter, the count just becomes -ve
        char_freq_in_target[s[end]] -= 1

        # If all letters in the target are found:
        while target_len == 0:
            window_len = end - start + 1
            if not shortest or window_len < len(shortest):
                # Note the new minimum window
                shortest = s[start : end + 1]

            # Increase the letter count of the current letter
            char_freq_in_target[s[start]] += 1

            # If all target letters have been seen and now, a target letter is seen with count > 0
            # Increase the target length to be found. This will break out of the loop
            if char_freq_in_target[s[start]] > 0:
                target_len += 1
                    
            start+=1
                
    return shortest

if __name__ == "__main__":
    #Input: s = "ADOBECODEBANC", t = "ABC"
    #Output: "BANC"
    s = "ADOBECODEBANC"
    t = "ABC"
    print(minWindow(s,k))
```
```kotlin
fun minWindow(s: String, t: String): String {
    // Base condition
    if (s.isEmpty() || t.isEmpty() || s.length < t.length)  return ""
    //val charFreqInTarget = IntArray(128){ 0 }
    //for(ch in t){
    //    ++charFreqInTarget[ch.toInt()]
    //}
    val charFreqInTarget = t.groupingBy { it }.eachCount().toMutableMap()        
    var start = 0
    var end = 0
    var shortest = ""
    var lengthOfTarget = t.length

    for (end in 0..s.length - 1) {
        //if (charFreqInTarget[s[end].toInt()]-- > 0) --lengthOfTarget
        if (charFreqInTarget.contains(s[end].toInt())) --lengthOfTarget
            
        while (lengthOfTarget == 0){
            if (shortest.isEmpty() || end - start + 1 < shortest.length){
                shortest = s.substring(start, end + 1)
            }
                
            //if (++charFreqInTarget[s[start].toInt()] > 0) ++lengthOfTarget
            if (charFreqInTarget.contains(s[start].toInt())) ++lengthOfTarget
            ++start
        }
    }
        
    return shortest
}

fun main(args: Array) {
    //Input: s = "anagram", t = "nagaram"
    //Output: true
    val s = "ADOBECODEBANC"
    val t = "ABC"
    println(minWindow(s,k))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

####  [LC-242:Valid Anagram](https://leetcode.com/problems/valid-anagram/)
##### Solution Explanation:
```
Approach: HashMap
=================================================================================================================================================================
Algorithm
----------
- Simple question. Build a frequency map for s. Now check t against this frequency map.
- Read the editorial about the followup about unicode characters.
- Unicode has 4 bytes per character. So 2^32 or 4 billion characters
- Using an array so big is not good. Use a hash-table.
```
##### Complexity Analysis:
```
Time  : O(N)
Space : O(N)
```
```python
def isAnagram(s: str, t: str) -> bool:
    if (len(s) != len(t)) return False
    """
    Use dict to check whether the count of every element is the same or not.
    """
    hashMap = {}
        
    for i in s:
        hashMap[i] = hashMap.get(i, 0) + 1
            
    for j in t:
        if j not in hashMap:
            return False
        else:
            #hashMap[j] -= 1
            hashMap[j] = hashMap.get(j, 0) - 1
        
    #for v in hashMap.values():
    #    if v != 0:
    #        return False
    #return True
    return False not in [hashMap[char] == 0 for char in hashMap]

if __name__ == "__main__":
    #Input: s = "anagram", t = "nagaram"
    #Output: true
    s = "anagram"
    t = "nagaram"
    print(isAnagram(s, t))
```
```kotlin
fun isAnagram(s: String, t: String): Boolean {
    if (s.length != t.length) return false
    val hashMap = HashMap<String, Int>()
    for (i in s)
        hashMap[i.toString()] = (hashMap[i.toString()] ?: 0) + 1
    for (j in t) {
        if (hashMap[j.toString()] == null)
            return false
        //hashMap[j.toString()] = hashMap[j.toString()]!! - 1
        hashMap[j.toString()] = (hashMap[i.toString()] ?: 0) - 1
    }
    return hashMap.values.all { it == 0 }
}

fun main(args: Array<String>) {
    //Input: s = "anagram", t = "nagaram"
    //Output: true
    val s = "anagram"
    val t = "nagaram"
    println(isAnagram(s, t))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

####  [LC-49:Group Anagrams](https://leetcode.com/problems/group-anagrams/)
##### Solution Explanation:
```
Approach: HashMap
=================================================================================================================================================================
Algorithm
----------
Explanation
We loop through each input string and determine the frequency of each letter in it, considering all 26 English lowercase letters.
We transform the frequency list into a tuple which will then be used as a key to access the list of anagrams.
The given input string will then be appended to it.

For the python solution we use a dict where:
key: tuple of the frequency of 26 letters, value: [string].

For the kotlin soluton we use a Pair where:
encodeToPair encode strings to a pair of array of letters' frequencies and the string itself.
```
##### Complexity Analysis:
```
Time  : O(N)
Space : O(N)

Definitions
n: The total number of letters.

Runtime Complexity
O(n) for examining each input letter once. The dictionary operations are expected to be in O(1).

Space Complexity
O(n) for storing each input string.
```
```python
import collections
from typing import List

def groupAnagrams(strs: List[str]) -> List[List[str]]:
    if not strs: 
        return []    
    if len(strs) < 2:
        return [strs]
    #key: tuple of the frequency of 26 letters, value: [string]
    d = collections.defaultdict(list)
    # only the frequency of each letter matters
    for s in strs:
        arr = [0] * 26
        for c in s:
            arr[ord(c) - ord('a')] += 1
        d[tuple(arr)].append(s)
    # turn the values of dict into a list
    return list(d.values())


if __name__ == "__main__":
    #Input: strs = ["eat","tea","tan","ate","nat","bat"]
    #Output: [["bat"],["nat","tan"],["ate","eat","tea"]]
    #Any order is acceptable, so the below is also correct:
    #[[eat, tea, ate], [tan, nat], [bat]]
    strs = ["eat","tea","tan","ate","nat","bat"]
    print(groupAnagrams(strs))
```
```kotlin
private val ALPHABET_LENGTH = 26
private val ASCII_OF_LOWERCASE_A = 97
    
private fun encodeToPair(str: String): Pair<String, String> {
    var theList = IntArray(ALPHABET_LENGTH) { 0 }
    for (char in str) ++theList[char.toInt() - ASCII_OF_LOWERCASE_A]
    return Pair(theList.joinToString(), str)
}
    
fun groupAnagrams(strs: Array<String>): List<List<String>> {
    if (strs.isEmpty()) {
        return emptyList()
    }
    if (strs.size < 2) {
        return emptyList()
    }
    return strs
        .map { encodeToPair(it) }
        .groupBy { it.first }
        .toList()
        .map{ it.second.map{ it.second } };
}

fun main(args: Array<String>) {
    //Input: strs = ["eat","tea","tan","ate","nat","bat"]
    //Output: [["bat"],["nat","tan"],["ate","eat","tea"]]
    //Any order is acceptable, so the below is also correct:
    //[[eat, tea, ate], [tan, nat], [bat]]
    val strs = arrayOf("eat","tea","tan","ate","nat","bat")
    println(groupAnagrams(strs))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

####  [LC-678:Valid Parenthesis String](https://leetcode.com/problems/valid-parenthesis-string/)
##### Solution Explanation:
```
Approach: Stack
=================================================================================================================================================================
Algorithm
----------
- We can use a stack to keep track of the order of brackets seen in s and know whether the current bracket matches the last one seen.
- The stack will only contain openBrackets and we will use the openToCloseBracket mapping to check that the latest open bracket has a matching close bracket.
- I.e. The top of the stack (stack[-1]) should be an open bracket that matches the current close bracket, if the current bracket is a close bracket.
```
##### Complexity Analysis:
```
Time  : O(N)
Space : O(N)

O(n) time, 
O(n) space for when each of the characters are "(", "{", "[". This would result in a stack of size n
```
```python
def isValid(s: str) -> bool:
    stack = []
    openBrackets = {'(', '[', '{'}
    openToCloseBracket = {
        '(': ')',
        '[': ']',
        '{': '}'
    }
        
    for bracket in s:
        # Only add open brackets to the stack.
        if bracket in openBrackets:
            stack.append(bracket)
        # The stack is non-empty, and the last open bracket seen matches the close bracket.
        elif stack and openToCloseBracket[stack[-1]] == bracket:
            # Remove the matching open bracket from the stack.
            stack.pop()
        # The rules for valid bracket matching have been violated.
        else:
            return False
                
    # All brackets must be paired up, so the stack must be empty by the end of the string.
    return len(stack) == 0

if __name__ == "__main__":
    #Input: s = "()[]{}"
    #Output: True
    s = "()[]{}"
    print(isValid(s))
```
```kotlin
import java.util.ArrayDeque

fun isValid(s: String): Boolean {
    val stack = ArrayDeque<Char>()

    s.forEach {

        when (it) {
            '(', '[', '{' -> stack.push(it)
            else -> {

                val end: Char = when (it) {
                    '}' -> '{'
                    ']' -> '['
                    ')' -> '('
                    else -> throw RuntimeException("Unknown char $it")
                }

                if (end != stack.poll()) return false
            }
        }
    }

    return stack.isEmpty()
}

fun main(args: Array<String>) {
    //Input: s = "()[]{}"
    //Output: True
    val s = "()[]{}"
    println(isValid(s))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

####  [LC-125:Valid Palindrome](https://leetcode.com/problems/valid-palindrome/)
##### Solution Explanation:
```
Approach: Two-Pointers
=================================================================================================================================================================
Algorithm
----------
- Normalize the string and convert to lowercase.
- Use 2 pointers start and end to compare characters.
- Skip non alphanumeric characters. Return False if characters do not match.
```
##### Complexity Analysis:
```
Time  : O(N)
Space : O(1)
=============================
Time Complexity  : O(N)
Since your checking every letter in the string. Because we move beg only to the right and end only to the left, until they meet.
Space Complexity : O(1)
We just use a couple of additional variables.
```
```python
def isPalindrome(s: str) -> bool:
    start, end = 0, len(s) - 1
    while start < end:
        if not s[start].isalnum():
            start += 1
        elif not s[end].isalnum():
            end -= 1
        elif s[start].lower() != s[end].lower():
            return False
        else:
            start += 1
            end -= 1
    return True

if __name__ == "__main__":
    #Input: s = "A man, a plan, a canal: Panama"
    #Output: true
    #Explanation: "amanaplanacanalpanama" is a palindrome.    s = "()[]{}"
    s = "A man, a plan, a canal: Panama"
    print(isPalindrome(s))
```
```kotlin
fun isPalindrome(s: String): Boolean {
    var start = 0
    var end = s.length - 1
    while (true) {
        if (start >= end) return true
        if (!isAlnum(s[start])) {
            start++
        } else if (!isAlnum(s[end])) {
            end--
        } else if (!equal(s[start], s[end])) {
            return false
        } else {
            start++
            end--
        }
    }
	return true
}

fun equal(char1: Char, char2: Char): Boolean {
    return char1.toLowerCase() == char2.toLowerCase()
}

fun isAlnum(char: Char): Boolean {
    return char in '0'..'9' || char in 'a'..'z' || char in 'A'..'Z'
}

fun main(args: Array<String>) {
    //Input: s = "A man, a plan, a canal: Panama"
    //Output: true
    //Explanation: "amanaplanacanalpanama" is a palindrome.
    val s = "A man, a plan, a canal: Panama"
    println(isPalindrome(s))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

####  [LC-5:Longest Palindromic Substring](https://leetcode.com/problems/longest-palindromic-substring/)
##### Solution Explanation:
```
Approach: Manacher's Algorithm
=================================================================================================================================================================
There are many ways to solve this problem. Most common way is to treat each character of the string as the center and expand left and right.
Keep track of their lengths and return the string with maximum length.

So, what’s the problem 🤔? The problem is the time complexity - it will be O(n2). Not so good, right?

Let’s see what’s hurting us. We are expanding left and right treating each character as the center. 
What if we only expand only at the necessary characters instead of expanding at each character?

Can we do that 🤔? Yes, we can using the Manacher’s Algorithm. This algorithm intelligently uses characteristics 
of a palindrome to solve the problem in linear O(n) time -

1. The left side of a palindrome is a mirror image of its right side.
2. Odd length palindrome will be centered on a letter and even length palindrome will be centered in between 
the two letters (thus there will be total 2n + 1 centers).

Manacher’s Algorithm deals with the problem of finding the new center.
Below are the steps -

1. Initialize the lengths array to the number of possible centers.
2. Set the current center to the first center.
3. Loop while the current center is valid:
 (a) Expand to the left and right simultaneously until we find the largest palindrome around this center.
 (b) Fill in the appropriate entry in the longest palindrome lengths array.
 (c) Iterate through the longest palindrome lengths array backwards and 
     fill in the corresponding values to the right of the entry for the current center 
	 until an unknown value (as described above) is encountered.
 (d) set the new center to the index of this unknown value.
4. Return the longest substring.
```
##### Complexity Analysis:
```
Time Complexity  : O(N)
========================
Note that at each step of the algorithm we’re either incrementing our current position in the input string or filling in an entry 
in the lengths array. Since the lengths array has size linear in the size of the input array, 
the algorithm has worst-case linear O(N) running time.

Space Complexity : O(N)
========================
Since we are using the palindrome array to store the length of palindromes centered at each character,
the space complexity will also be O(N).
```
```python
def get_updated_string(s):
    sb = ''
    for i in range(0, len(s)):
        sb += '#' + s[i]
    sb += '#'
    return sb

# Manacher's Algorithm
def longestPalindrome(s: str) -> str:
    # Update the string to put hash "#" at the beginning, end and in between each character
    updated_string = get_updated_string(s)
    # Length of the array that will store the window of palindromic substring
    length = 2 * len(s) + 1
    # List to store the length of each palindrome centered at each element
    p = [0] * length
    # Current center of the longest palindromic string
    c = 0
    # Right boundary of the longest palindromic string
    r = 0
    # Maximum length of the substring
    maxLength = 0
    # Position index
    position = -1
    for i in range(0, length):
        # Mirror of the current index
        mirror = 2 * c - i
        # Check if the mirror is outside the left boundary of current longest palindrome
        if i < r:
            p[i] = min(r - i, p[mirror])
        # Indices of the characters to be compared
        a = i + (1 + p[i])
        b = i - (1 + p[i])
        # Expand the window
        while a < length and b >= 0 and updated_string[a] == updated_string[b]:
            p[i] += 1
            a += 1
            b -= 1
        # If the expanded palindrome is expanding beyond the right boundary of
        # the current longest palindrome, then update c and r
        if i + p[i] > r:
            c = i
            r = i + p[i]
        if maxLength < p[i]:
            maxLength = p[i]
            position = i
    offset = p[position]
    result = ''
    for i in range(position - offset + 1, position + offset):
        if updated_string[i] != '#':
            result += updated_string[i]
    return result

if __name__ == "__main__":
    #Input: s = "babad"
    #Output: "bab"
    #Note: "aba" is also a valid answer.
    s = "babad"
    print(longestPalindrome(s))
```
```kotlin
fun getUpdatedString(s: String): String {
    val sb = StringBuilder()
    for (element in s) {
        sb.append("#").append(element)
    }
    sb.append("#")
    return sb.toString()
}

// Manacher's Algorithm
fun longestPalindrome(s: String): String {
    // Update the string to put hash "#" at the beginning, end and in between each character
    val updatedString = getUpdatedString(s)
    // Length of the array that will store the window of palindromic substring
    val length = 2 * s.length + 1
    // Array to store the length of each palindrome centered at each element
    val p = IntArray(length)
    // Current center of the longest palindromic string
    var c = 0
    // Right boundary of the longest palindromic string
    var r = 0
    // Maximum length of the substring
    var maxLength = 0
    // Position index
    var position = -1
    for (i in 0 until length) {
        // Mirror of the current index
        val mirror = 2 * c - i
        // Check if the mirror is outside the left boundary of current longest palindrome
        if (i < r) {
            p[i] = (r - i).coerceAtMost(p[mirror])
        }
        // Indices of the characters to be compared
        var a = i + (1 + p[i])
        var b = i - (1 + p[i])
        // Expand the window
        while (a < length && b >= 0 && updatedString[a] == updatedString[b]) {
            p[i]++
            a++
            b--
        }
        // If the expanded palindrome is expanding beyond the right boundary of
        // the current longest palindrome, then update c and r
        if (i + p[i] > r) {
            c = i
            r = i + p[i]
        }
        if (maxLength < p[i]) {
            maxLength = p[i]
            position = i
        }
    }
    val offset = p[position]
    val result = StringBuilder()
    for (i in position - offset + 1 until position + offset) {
        if (updatedString[i] != '#') {
            result.append(updatedString[i])
        }
    }
    return result.toString()
}

fun main(args: Array<String>) {
    //Input: s = "babad"
    //Output: "bab"
    //Note: "aba" is also a valid answer.    
    val s = "babad"
    println(longestPalindrome(s))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

####  [LC-647:Palindromic Substrings](https://leetcode.com/problems/palindromic-substrings/)
##### Solution Explanation:
```
Approach: Manacher's Algorithm
=================================================================================================================================================================
ALGO:

1. Define formatted_string by adding # between characters, and @ at the beginning.
2. An array(palindrome_count) is used to mark the radius of the largest odd-length palindromic substring centered at index
3. Traverse throughout the length of the formatted string, and check for 3 basic conditions:
   i.   if the mirror is still in a valid position, update palindrome count for that index as minimum of current radius and expected mirror radius
   ii.  if the mirror is valid for all characters from the centre to the respective boundaries, keep updating the palindromic count
   iii. Shift the mirror if palindrome is found
```
##### Complexity Analysis:
```
Time Complexity  : O(N)
========================
Note that at each step of the algorithm we’re either incrementing our current position in the input string or filling in an entry 
in the lengths array. Since the lengths array has size linear in the size of the input array, 
the algorithm has worst-case linear O(N) running time.

Space Complexity : O(N)
========================
Since we are using the palindrome array to store the length of palindromes centered at each character,
the space complexity will also be O(N).
```
```python
def countSubstrings(s: str) -> int:
    # Pre-processed for Manacher's Algorithm
    formatted_string = '@#' + '#'.join(s) + '#$'
    palindrome_count  = [0] * len(formatted_string)
    maxRight          = 0 # The most-right position ever touched by sub-strings
    center            = 0 # The center for the sub-string touching the maxRight
    for i in range(1, len(formatted_string) - 1):
        if i < maxRight:
            palindrome_count[i] = min(maxRight - i, palindrome_count[2 * center - i])
        while formatted_string[i + palindrome_count[i] + 1] == formatted_string[i - palindrome_count[i] - 1]:
            palindrome_count[i] += 1
        if i + palindrome_count[i] > maxRight:
            center = i
            maxRight = i + palindrome_count[i]
    return sum((v+1)//2 for v in palindrome_count)

if __name__ == "__main__":
    #Input: s = "abc"
    #Output: 3
    #Explanation: Three palindromic strings: "a", "b", "c".
    s = "abc"
    print(countSubstrings(s))
```
```kotlin
fun countSubstrings(input: String): Int {
    val formatted_string = CharArray(2 * input.length() + 3)
    formatted_string[0] = '@'
    formatted_string[1] = '#'
    formatted_string[formatted_string.size - 1] = '$'
    var t = 2
    for (c in input.toCharArray()) {
        formatted_string[t++] = c
        formatted_string[t++] = '#'
    }
    val palindrome_count = IntArray(formatted_string.size)
    // The center for the sub-string touching the maxRight
    var center = 0
    // The most-right position ever touched by sub-strings
    var maxRight = 0
    for (index in 1 until palindrome_count.size - 1) {
        if (index < maxRight) {
            palindrome_count[index] = Math.min(maxRight - index, palindrome_count[2 * center - index]) //min of (current mirror radius, expected mirror radius)
        }
        //mirror
        while (formatted_string[index + palindrome_count[index] + 1] == formatted_string[index - palindrome_count[index] - 1]) {
            palindrome_count[index]++
        }
        if (index + palindrome_count[index] > maxRight) { //shift the mirror if palindrome found
            center = index
            maxRight = index + palindrome_count[index]
        }
    }
    var ans = 0
    for (v in palindrome_count) ans += (v + 1) / 2
    return ans
}

fun main(args: Array<String>) {
    //Input: s = "abc"
    //Output: 3
    //Explanation: Three palindromic strings: "a", "b", "c".
    val s = "abc"
    println(countSubstrings(s))
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

####  [LC-271:Encode and Decode Strings](https://leetcode.com/problems/encode-and-decode-strings/)
##### Solution Explanation:
```
References: http://leetcode.libaoj.in/encode-and-decode-strings.html
            https://developer.mozilla.org/en-US/docs/Web/HTTP/Headers/Transfer-Encoding
			https://en.wikipedia.org/wiki/Chunked_transfer_encoding
=================================================================================================================================================================
Approach: Chunked Transfer Encoding ( Similar to encoding used in HTTP v1.1 )
=================================================================================================================================================================
Pay attention to this approach because last year Google likes to ask that sort of low-level optimisation.
Serialize and deserialize BST problem is a similar example.

This approach is based on the encoding used in HTTP v1.1 . It doesn't depend on the set of input characters, 
and hence is more versatile and effective than Approach 1 [ Approach 1: Non-ASCII Delimiter ].

Data stream is divided into chunks. Each chunk is preceded by its size in bytes.

Encoding Algorithm
====================================================

        +---+---+---+ +---+---+ 
Input	| t | b | c | | d | e |
        +---+---+---+ +---+---+ 

Encode          Each chunk is preceded 
                by its 4-bytes size

+---+---+---+---+---+---+---+---+---+---+---+---+---+ 
| 0 | 0 | 0 | 3 | t | b | c | 0 | 0 | 0 | 2 | d | e |
+---+---+---+---+---+---+---+---+---+---+---+---+---+ 
       \         /                 \          /
        \       /                   \        /
    size of next chunk               \      /
                                 size of next chunk

- Iterate over the array of chunks, i.e. strings.
  + For each chunk compute its length, and convert that length into 4-bytes string.
  + Append to encoded string :
    [] 4-bytes string with information about chunk size in bytes.
    [] Chunk itself.

- Return encoded string.

Decoding Algorithm
====================================================

        +---+---+---+ +---+---+ 
Input	| t | b | c | | d | e |
        +---+---+---+ +---+---+ 

----------------------------------------------------------------

Encode          Each chunk is preceded 
                by its 4-bytes size

+---+---+---+---+---+---+---+---+---+---+---+---+---+ 
| 0 | 0 | 0 | 3 | t | b | c | 0 | 0 | 0 | 2 | d | e |
+---+---+---+---+---+---+---+---+---+---+---+---+---+ 
       \         /                 \          /
        \       /                   \        /
    size of next chunk               \      /
                                 size of next chunk

----------------------------------------------------------------

Decode             1. Read next chunk length
                   2. Read chunk itself and add it to output


        +---+---+---+ +---+---+ 
	| t | b | c | | d | e |
        +---+---+---+ +---+---+ 

- Iterate over the encoded string with a pointer i initiated as 0. While i < n :
  + Read 4 bytes s[i: i + 4] . It's chunk size in bytes. Convert this 4-bytes string to integer length .
  + Move the pointer by 4 bytes i += 4 .
  + Append to the decoded array string s[i: i + length] .
  + Move the pointer by length bytes i += length .
- Return decoded array of strings.
```
##### Complexity Analysis:
```
Time Complexity  : O(N)
========================
O(N) both for encode and decode, where N is a number of strings in the input array.

Space Complexity : O(1)
========================
O(1) for encode to keep the output, since the output is one string.O(N) for decode keep the output, since the output is an array of strings.
```
```python
class Codec:
    def len_to_str(self, x):
        """
        Encodes length of string to bytes string
        """
        x = len(x)
        bytes = [chr(x >> (i * 8) & 0xff) for i in range(4)]
        bytes.reverse()
        bytes_str = ''.join(bytes)
        return bytes_str
    
    def encode(self, strs):
        """Encodes a list of strings to a single string.
        :type strs: List[str]
        :rtype: str
        """
        # encode here is a workaround to fix BE CodecDriver error
        return ''.join(self.len_to_str(x) + x.encode('utf-8') for x in strs)
        
    def str_to_int(self, bytes_str):
        """
        Decodes bytes string to integer.
        """
        result = 0
        for ch in bytes_str:
            result = result * 256 + ord(ch)
        return result
    
    def decode(self, s):
        """Decodes a single string to a list of strings.
        :type s: str
        :rtype: List[str]
        """
        i, n = 0, len(s)
        output = []
        while i < n:
            length = self.str_to_int(s[i: i + 4])
            i += 4
            output.append(s[i: i + length])
            i += length
        return output

if __name__ == "__main__":
    input = "Hello World"
    print(f'Original word: {input}')
    strs = [s.strip() for s in input.split(' ')]
    codec = Codec()
    encodedInput = codec.encode(strs)
    decodedInput = codec.decode(encodedInput)
    result = ' '.join(decodedInput)
    print(f'Decoded word: {result}')
```
```kotlin
class Codec {
    // Encodes string length to bytes string
    fun intToString(s: String): String {
        val x: Int = s.length
        val bytes = CharArray(4)
        for (i in 3 downTo -1 + 1) {
            bytes[3 - i] = (x shr i * 8 and 0xff).toChar()
        }
        return String(bytes)
    }

    // Encodes a list of strings to a single string.
    fun encode(strs: List<String>): String {
        val sb = StringBuilder()
        for (s in strs) {
            sb.append(intToString(s))
            sb.append(s)
        }
        return sb.toString()
    }

    // Decodes bytes string to integer
    fun stringToInt(bytesStr: String): Int {
        var result = 0
        for (b in bytesStr.toCharArray()) result = (result shl 8) + b.code
        return result
    }

    // Decodes a single string to a list of strings.
    fun decode(s: String): List<String> {
        var i = 0
        val n = s.length
        val output = mutableListOf<String>()
        while (i < n) {
            val length = stringToInt(s.substring(i, i + 4))
            i += 4
            output.add(s.substring(i, i + length))
            i += length
        }
        return output
    }

}

fun main(args: Array<String>) {
    val input = "Hello World"
    println("Original word: [$input]")
    var strs: List<String> = input.split(",").map { it.trim() }
    val codec = Codec()
    val encodedInput = codec.encode(strs)
    val decodedInput = codec.decode(encodedInput)
    val result = decodedInput.joinToString(" ")
    println("Decoded word: [$result]")
}
```

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>

