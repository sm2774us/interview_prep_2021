# Solution-1 ( Using Binet's Formula or the Golden Ratio  ).
# Time: O(1)
# Space: O(1)
# Reference: https://math.libretexts.org/Bookshelves/Applied_Mathematics/Book%3A_College_Mathematics_for_Everyday_Life_(Inigo_et_al)/10%3A_Geometric_Symmetry_and_the_Golden_Ratio/10.04%3A_Fibonacci_Numbers_and_the_Golden_Ratio#:~:text=The%20number%20that%20these%20ratios,this%20number%20in%20Binet%27s%20formula.&text=The%20Golden%20Ratio%20has%20the,for%20a%20variety%20of%20reasons.
class Solution1:
    def fib(self, N: int) -> int:
        phi = round((1 + 5 ** 0.5) / 2, 6)
        return round((phi ** N - (-phi) ** (-N)) / (5 ** 0.5))

# Solution-2 ( Using Tail Recursion  ).
# Time: O(N)
# Space: O(1)
class Solution2:
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

##### Solution-3 ( Using Matrix Exponentiation - DP ).

##### Solution Explanation
# Time:  O(log(N))
# Space: O(1)
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
