import sys
import numpy as np


INF = sys.maxsize - 1


# DP algorithm adapted from Langmead's notebooks.
def trace(D: np.ndarray, p: str, t: str) -> int:
    ''' Backtrace edit-distance matrix D for strings x and y '''
    i, j = len(p), len(t)
    
    while i > 0:
        diag, vert, horz = INF, INF, INF
        if i > 0 and j > 0:
            delt = 0 if p[i-1] == t[j-1] else 1
            diag = D[i-1, j-1] + delt
        if i > 0:
            vert = D[i-1, j] + 1
        if j > 0:
            horz = D[i, j-1] + 1

        if diag <= vert and diag <= horz:  # diagonal was best
            i -= 1
            j -= 1
        elif vert <= horz: # vertical was best; this is an insertion in x w/r/t y
            i -= 1
        else: # horizontal was best
            j -= 1
   
    return j # offset of the first (leftmost) character of t involved in the alignment


def kEditDp(p: str, t: str, max_shift: int = INF, k: int = INF) -> tuple[int, int, int] | None:
    ''' Find the alignment of p to a substring of t with the fewest edits, up to k,
        and shifted by at most max_shift characters.
        Return the edit distance and the coordinates of the substring. '''
    len_p, len_t = len(p), len(t)
    D = np.zeros((len_p + 1, len_t + 1), dtype=int) + INF
    D[0, :] = np.zeros(len_t + 1)
    D[:, 0] = range(len_p + 1)

    for i in range(1, len_p + 1):
        start_j = max(1, i - max_shift)
        end_j = min(len_t + 1, i + max_shift + 1)
        current_min = INF

        for j in range(start_j, end_j):
            delta = 1 if p[i - 1] != t[j - 1] else 0
            D[i, j] = min(D[i - 1, j - 1] + delta, D[i - 1, j] + 1, D[i, j - 1] + 1)
            current_min = min(current_min, D[i, j])

        if current_min > k:
            return None
    
    # Find minimum edit distance in last row:
    start_j = max(1, len_p - max_shift)
    end_j = min(len_t + 1, len_p + max_shift + 1)

    best_j = start_j + int(np.argmin(D[len_p, start_j:end_j]))
    best_score = D[len_p, best_j]
    
    # Backtrace; note: stops as soon as it gets to first row
    off = trace(D, p, t[:best_j])
    
    return best_score, off, best_j
