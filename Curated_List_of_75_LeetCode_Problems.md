# LeetCode Curated List of 75 Problems.

## Algorithms

* [Math](#math)

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
| #     | Title	               | url                                                 | Time                  | Space        | Difficulty | Tag	                               | Note                   |
| ----- | -------------------- | --------------------------------------------------- | --------------------- | ------------ | ---------- | ----------------------------------- | ---------------------- |
| 0012  | Integer to Roman     | https://leetcode.com/problems/integer-to-roman/     | _O(n)_                | _O(1)_       | Medium     |                                     |                        |
| 0013  | Roman to Integer     | https://leetcode.com/problems/roman-to-integer/     | _O(n)_                | _O(1)_       | Easy       |                                     |                        |
| 0029  | Divide Two Integers  | https://leetcode.com/problems/divide-two-integers/  | _O(1)_                | _O(1)_       | Medium     |                                     |                        |
| 0050  | Pow(x, n)            | https://leetcode.com/problems/powx-n/               | _O(1)_                | _O(1)_       | Medium     |                                     |                        |
| 0149  | Max Points on a Line | https://leetcode.com/problems/max-points-on-a-line/ | _O(n^2)_	             | _O(n)_       | Hard       |                                     | Linear Equation `ax + by + c = 0` |
| 0204  | Count Primes         | https://leetcode.com/problems/count-primes/         | _O( N*Log(Log(N)) )_  | _O(N)_       | Easy       |                                     | Sieve of Eratosthenes  |
| 0372  | Super Pow            | https://leetcode.com/problems/super-pow/            | _O(n)_                | _O(1)_       | Medium     |                                     |                        |
| 0509  | Fibonacci Number     | https://leetcode.com/problems/fibonacci-number/     | _O(logn)_             | _O(1)_       | Easy       | variant of [Climbing Stairs](https://leetcode.com/problems/climbing-stairs/) | Matrix Exponentiation, Binet's Formula |
| 1390  | Four Divisors        | https://leetcode.com/problems/four-divisors/        | _O(N+K*Log(Log(K)))_, where, N = max(nums), M=len(nums), K is len(primes) | _O(N+M+K^2)_, where, N = max(nums), M=len(nums), K is len(primes) | Medium     |                                     | Sieve of Eratosthenes  |
| 1390  | Four Divisors        | https://leetcode.com/problems/four-divisors/        | _O(N * sqrt(M))_, where, N = length of nums and M = nums[i]  | _O(1)_       | Medium     |                                     | Recursion              |
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
