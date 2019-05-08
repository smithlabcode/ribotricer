# Part of ribotricer software
#
# Copyright (C) 2019 Wenzheng Li, Saket Choudhary and Andrew D Smith
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.


class ITNode:
    """Class for interval tree node"""

    def __init__(self, interval, left=None, right=None):
        self.interval = interval
        # upper is the max end of the subtree rooted at the current node
        self.upper = interval.end
        self.left = left
        self.right = right


class IntervalTree:
    """Class for interval tree
       All the intervals used in this project is 1-based and closed
    """

    def __init__(self):
        self.root = None

    def insert(self, interval):
        self.root = self._add(self.root, interval)

    def _add(self, root, interval):
        if root is None:
            root = ITNode(interval)
        else:
            low = root.interval.start
            if interval.start < low:
                root.left = self._add(root.left, interval)
            else:
                root.right = self._add(root.right, interval)
            if root.upper < interval.end:
                root.upper = interval.end
        return root

    def search(self, start, end):
        results = []
        self._find(self.root, start, end, results)
        return results

    def _find(self, root, start, end, results):
        if root is None:
            return
        if not (root.interval.start > end or root.interval.end < start):
            results.append(root.interval)
        if root.left is not None and root.left.upper >= start:
            self._find(root.left, start, end, results)
        if root.right is not None and root.interval.start <= end:
            self._find(root.right, start, end, results)
