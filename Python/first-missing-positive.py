from typing import List

# Time:  O(n)
# Space: O(1)

#
# Notice, that first missing positive will always be in range 1,2,...n,n+1,
# where n is length of nums. Let us rearrange numbers, putting each number in its place:
# number i on place i-1 (because indexes start with 0): let us iterate over our numbers
# and change two numbers if one of them not on its place: we break if number not in
# range 1,...,n or if we are trying to put number on the place,
# which is already occupied with this place (because we have infinite loop in this case).
#
# When we iterate all numbers we find for number which is not on its place,
# using i == nums[i] - 1. It can happen that all numbers between 1 and n are here,
# so I add [0] to the end. Finally, I found index of False in this array: it will be our number.
#
# Complexity: even though we have while loop inside for loop, complexity will be O(n):
#             on each step we put at least one number to its proper place.
#             Additional space complexity is O(1), however we modify our nums.
#
# This solution is very similar to problem 442. Find All Duplicates in an Array,
# see my detailed solution here:
#
# https://leetcode.com/problems/find-all-duplicates-in-an-array/discuss/775738/Python-2-solutions-with-O(n)-timeO(1)-space-explained
#
#
class Solution(object):
    # @param A, a list of integers
    # @return an integer
    def firstMissingPositive(self, nums: List[int]) -> int:
        n = len(nums)
        for i in range(n):
            while nums[i] - 1 in range(n) and nums[i] != nums[nums[i] - 1]:
                nums[nums[i] - 1], nums[i] = nums[i], nums[nums[i] - 1]

        return next((i + 1 for i, num in enumerate(nums) if num != i + 1), n + 1)
