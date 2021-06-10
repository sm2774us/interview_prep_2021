import unittest

# Definition for a Node.
class Node:
    def __init__(self, val = 0, neighbors = None):
        self.val = val
        self.neighbors = neighbors if neighbors is not None else []

    def flatten_into_str(self):
        return '{} neighbors: {}'.format(
            self.val,
            "".join([x.__repr__() for x in self.neighbors])
        )
			
    def __repr__(self):
        return '{} neighbors: {}'.format(
            self.val,
            [x.val for x in self.neighbors]
        )

class Solution:
    def cloneGraph(self, node: 'Node') -> 'Node':
        if not node: return
        d = {node: Node(node.val)}
        stack = [node]
        while stack:
            curNode = stack.pop()
            for nei in curNode.neighbors:
                if nei not in d:
                    d[nei] = Node(nei.val)
                    stack.append(nei)
                d[curNode].neighbors.append(d[nei])
        return d[node] # return the value of the original node which is a copy of that original node

class Test(unittest.TestCase):
    def setUp(self):
        pass
	
    def tearDown(self):
        pass
		
    def test_cloneGraph(self):
        solution = Solution()
        root = Node(1, [Node(2, [Node(1),Node(3, [Node(2), Node(4)])]), Node(4, [Node(1),Node(3, [Node(2), Node(4)])])])
        self.assertEqual(root.flatten_into_str(), solution.cloneGraph(root).flatten_into_str())

if __name__ == "__main__":
    unittest.main()
