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

## String
| #     | Title	                                         | url                                                                           | Time   | Space   | Difficulty | Tag	        | Note                   |
| ----- | ---------------------------------------------- | ----------------------------------------------------------------------------- | ------ | ------- | ---------- | ------------ | ---------------------- |
| 0149  | Longest Substring Without Repeating Characters | https://leetcode.com/problems/longest-substring-without-repeating-characters/ | _O(n)_ | _O(1)_  | Medium     |              |                        |
| 0424  | Longest Repeating Character Replacement        | https://leetcode.com/problems/longest-repeating-character-replacement/        | _O(n)_ | _O(1)_  | Medium     |              | Sliding Window         |
| 0076  | Minimum Window Substring                       | https://leetcode.com/problems/minimum-window-substring/                       | _O(n)_ | _O(k)_  | Hard       |              |                        |
| 0242  | Valid Anagram                                  | https://leetcode.com/problems/valid-anagram/                                  | _O(n)_ | _O(1)_  | Easy       |              |                        |
| 0049  | Group Anagrams                                 | https://leetcode.com/problems/group-anagrams/                                 | _O(n)_ | _O(1)_  | Easy       |              |                        |
| 0678  | Valid Parenthesis String                       | https://leetcode.com/problems/valid-parenthesis-string/                       | _O(n)_ | _O(1)_  | Medium     |              |                        |
| 0125  | Valid Palindrome                               | https://leetcode.com/problems/valid-palindrome/                               | _O(n)_ | _O(1)_  | Easy       |              |                        |
| 0005  | Longest Palindromic Substring                  | https://leetcode.com/problems/longest-palindromic-substring/                  | _O(n)_ | _O(n)_  | Medium     |              | `Manacher's Algorithm` |
| 0647  | Palindromic Substrings                         | https://leetcode.com/problems/palindromic-substrings/                         | _O(n)_ | _O(n)_  | Medium     |              | `Manacher's Algorithm` |
| 0271  | Encode and Decode Strings                      | https://leetcode.com/problems/encode-and-decode-strings/                      | _O(n)_ | _O(1)_  | Medium     | 🔒           |                        |

<br/>
<div align="right">
    <b><a href="#algorithms">⬆️ Back to Top</a></b>
</div>
<br/>
