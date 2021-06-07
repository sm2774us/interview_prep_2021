from typing import List
# Time:  O(n)
# Space: O(1)

#
# Solution 1: rearrange numbers
#
# If we are not allowed to use additional memory, the only memory we can use is our original list nums.
# Let us try to put all numbers on they own places, by this I mean, we try to up number k to place k-1.
# Let us iterate through numbers, and change places for pairs of them, if it is needed.
#
# Consider an example [4, 3, 2, 7, 8, 2, 3, 1]:
#
# 1. We look at first number 4 and change it with number, which is on the place with number 4-1,
#    so we have [7, 3, 2, 4, 8, 2, 3, 1], i = 0.
# 2. We still look at first number and now we put it to place 7-1: [3, 3, 2, 4, 8, 2, 7, 1], i = 0.
# 3. We still look at first number and now we put it to place 3-1: [2, 3, 3, 4, 8, 2, 7, 1], i = 0.
# 4. We still look at first number and now we put it to place 2-1: [3, 2, 3, 4, 8, 2, 7, 1], i = 0.
#    Now, if we look at first place, we see number 3, which we need to put on place, where we already have number 3.
#    So, we stop with i = 0.
# 5. Continue with i = 1, 2, 3 and see that number is already on its place, so we do nothin.
# 6. For i = 4, we do: [3, 2, 3, 4, 1, 2, 7, 8].
# 7. Still, i = 4, we do: [1, 2, 3, 4, 3, 2, 7, 8], we stop here, because nums[4] = 3 but number 3 is already on place 3-1.
# 8. We continue with i = 5, 6, 7 and do nothing
#
# Finally, what we need to do with obtained array: find all places with wrong numbers and return these places,
# here it is numbers 3 and 2: [1, 2, 3, 4, 3, 2, 7, 8].
#
# Complexity:
#
# time complexity O(n), because with each swap on numbers, we put at least one of them on its place.
# Additional space complexity is O(1).
#
class Solution(object):
    def findDuplicates(self, nums: List[int]) -> List[int]:
        """
        :type nums: List[int]
        :rtype: List[int]
        """
        for i in range(len(nums)):
            while i != nums[i] - 1 and nums[i] != nums[nums[i] - 1]:
                nums[nums[i] - 1], nums[i] = nums[i], nums[nums[i] - 1]

        return [nums[it] for it in range(len(nums)) if it != nums[it] - 1]


# Time:  O(n)
# Space: O(1)

#
# Solution 2: hash numbers
#
# Actually, there is more clean, but much more tricky solution, where we again change our nums,
# but in different way.
# Let us again go through the same example and see how it works:
#
# 1. First, we look at number 4, look at place number 3, see, that number there is positive, it means,
#    we did not met number 4 yet, so we just change sign and have: nums = [4, 3, 2, -7, 8, 2, 3, 1], ans = [].
# 2. Similarly, we look at next number 3 and change sign of number in cell 3-1:
#    nums = [4, 3, -2, -7, 8, 2, 3, 1], ans = [].
# 3. Continue with next number, which is -2 now. Number is negative, but what we are interested in is nums[2-1],
#    which is 3, so we change its sign and have nums = [4, -3, -2, -7, 8, 2, 3, 1], ans = [].
# 4. Next number is -7, so we have nums = [4, -3, -2, -7, 8, 2, -3, 1], ans = [].
# 5. Next number is 8, so we have nums = [4, -3, -2, -7, 8, 2, -3, -1], ans = [].
# 6. Next number is 2, we look at place with number 2-1 = 1 and see, that we have there negative number.
#    It means, that we already seen number 2 before, so it is duplicate,
#    and we add it to our answer: nums = [4, -3, -2, -7, 8, 2, -3, -1], ans = [2].
# 7. Next number is -3, so we look at cell with number 3-1, where we see negative number, so we see number 3 before,
#    we update: nums = [4, -3, -2, -7, 8, 2, -3, -1], ans = [2, 3].
# 8. Finally, we see number -1, so we change sign of number with index 1-1:
#    ans = [-4, -3, -2, -7, 8, 2, -3, -1], ans = [2, 3].
#
# Complexity: again, time complexity is O(n), because we iterate once over our nums, space complexity is O(1).
#
class Solution2(object):
    def findDuplicates(self, nums: List[int]) -> List[int]:
        """
        :type nums: List[int]
        :rtype: List[int]
        """
        ans = []
        for num in nums:
            if nums[abs(num) - 1] < 0:
                ans.append(abs(num))
            else:
                nums[abs(num) - 1] *= -1
        return ans
