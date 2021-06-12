import unittest

# Definition for singly-linked list.
class ListNode:
    def __init__(self, val=0, next=None):
        self.val = val
        self.next = next

    @classmethod
    def initList(self, nums):
        if not nums:
            return None
        head = None
        current = None

        for n in nums:
            if not head:
                head = ListNode(n)
                current = head
            else:
                node = ListNode(n)
                current.next = node
                current = node
        return head

    @classmethod
    def linkedListToList(self, head):
        if not head:
            return []

        pointer = head
        sll_list = []
        while pointer:
            sll_list.append(pointer.val)
            pointer = pointer.next
        return sll_list


class Solution:
    def removeNthFromEndUsingOnePointerTwoPass(self, head: ListNode, n: int) -> ListNode:
	    ptr, length = head, 0
	    while ptr:
		    ptr, length = ptr.next, length + 1
	    if length == n : return head.next
	    ptr = head
	    for i in range(1, length - n):
		    ptr = ptr.next
	    ptr.next = ptr.next.next
	    return head

    def removeNthFromEndUsingTwoPointersOnePass(self, head: ListNode, n: int) -> ListNode:
	    fast = slow = head
	    for i in range(n):
		    fast = fast.next
	    if not fast: return head.next
	    while fast.next:
		    fast, slow = fast.next, slow.next
	    slow.next = slow.next.next
	    return head

class Test(unittest.TestCase):
    def setUp(self) -> None:
        pass

    def tearDown(self) -> None:
        pass

    def test_removeNthFromEnd(self) -> None:
        sol = Solution()
        for head, n, solution in (
            [
                [1, 2, 3, 4, 5], 2,
                [1, 2, 3, 5],
            ],
            [
                [1], 1,
                []
            ],
            [
                [1, 2], 1,
                [1]
            ]
        ):
            self.assertEqual(
                ListNode.linkedListToList(sol.removeNthFromEndUsingOnePointerTwoPass(ListNode.initList(head), n)),
                solution
            )
            self.assertEqual(
                ListNode.linkedListToList(sol.removeNthFromEndUsingTwoPointersOnePass(ListNode.initList(head), n)),
                solution
            )

if __name__ == "__main__":
    unittest.main()