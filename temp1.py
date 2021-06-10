from collections import deque
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
        # map original nodes to their clones
        d = {node : Node(node.val)}
        q = deque([node])        
        while q:
            for i in range(len(q)):
                currNode = q.popleft()
                for nei in currNode.neighbors:
                    if nei not in d:
                        # store copy of the neighboring node
                        d[nei] = Node(nei.val)
                        q.append(nei)
                    # connect the node copy at hand to its neighboring nodes (also copies) -------- [1]
                    d[currNode].neighbors.append(d[nei])
        # return copy of the starting node ------- [2]
        return d[node]

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