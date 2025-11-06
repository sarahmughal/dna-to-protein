from typing import List, Tuple
from dataclasses import dataclass
import math

@dataclass
class Node:
    name: str = None
    left: 'Node' = None
    right: 'Node' = None
    height: float = 0.0  # UPGMA height (ultrametric)

def upgma(names: List[str], D: List[List[float]]) -> Node:
    n = len(names)
    clusters = {i: Node(name=names[i], height=0.0) for i in range(n)}
    dist = { (i,j): D[i][j] for i in range(n) for j in range(n) if i<j }
    sizes = { i: 1 for i in range(n) }
    next_id = n
    while len(clusters) > 1:
        # find closest pair
        (i,j), mind = min(dist.items(), key=lambda kv: kv[1])
        # merge
        hi = clusters[i]; hj = clusters[j]
        new_h = mind/2
        new_node = Node(name=f"C{next_id}", left=hi, right=hj, height=new_h)
        # update structures
        clusters[next_id] = new_node
        size_i, size_j = sizes[i], sizes[j]
        # recompute distances to new cluster
        for k in list(clusters.keys()):
            if k in (i,j,next_id): 
                continue
            dik = dist[tuple(sorted((i,k)))]
            djk = dist[tuple(sorted((j,k)))]
            d_new = (size_i*dik + size_j*djk) / (size_i + size_j)
            dist[(min(k,next_id), max(k,next_id))] = d_new
        # remove old distances
        keys_to_del = [k for k in dist.keys() if i in k or j in k]
        for k in keys_to_del:
            del dist[k]
        # finalize
        del clusters[i]; del clusters[j]
        sizes[next_id] = size_i + size_j
        next_id += 1
    # return the sole cluster
    return next(iter(clusters.values()))

def _branch_length(parent_h, child: Node) -> float:
    return parent_h - child.height

def to_newick(node: Node) -> str:
    if node.left is None and node.right is None:
        return node.name
    l = to_newick(node.left)
    r = to_newick(node.right)
    bl_l = _branch_length(node.height, node.left)
    bl_r = _branch_length(node.height, node.right)
    return f"({l}:{bl_l:.4f},{r}:{bl_r:.4f})" + (node.name or '')

