from __future__ import annotations
from array import array
from collections import Counter


ALPHABET = "$ACGT"


# Karkkainen-Sanders algorithm
# from https://code.google.com/archive/p/pysuffix/

def radixpass(a, b, r, n, k) :
  c = array("i", [0]*(k+1))
  for i in range(n) :
    c[r[a[i]]]+=1

  somme = 0
  for i in range(k+1):
    freq, c[i] = c[i], somme
    somme += freq

  for i in range(n) :
    b[c[r[a[i]]]] = a[i]
    c[r[a[i]]] += 1

def kark_sort(s, SA, n, K):
  n0  = (n+2) // 3
  n1  = (n+1) // 3
  n2  = n // 3
  n02 = n0 + n2
      
  SA12 = array('i', [0]*(n02+3))
  SA0  = array('i', [0]*n0)

  s12 = [i for i in range(n+(n0-n1)) if i%3] 
  s12.extend([0]*3)
  s12 = array('i', s12)

  radixpass(s12, SA12, s[2:], n02, K)
  radixpass(SA12, s12, s[1:], n02, K)
  radixpass(s12, SA12, s, n02, K)

  name = 0
  c0, c1, c2 = -1, -1, -1
  for i in range(n02) :
    if s[SA12[i]] != c0 or s[SA12[i]+1] != c1 or s[SA12[i]+2] != c2 :
      name += 1
      c0 = s[SA12[i]]
      c1 = s[SA12[i]+1]
      c2 = s[SA12[i]+2]
    if SA12[i] % 3 == 1 :
      s12[SA12[i]//3] = name
    else :
      s12[SA12[i]//3 + n0] = name

  if name < n02 :
    kark_sort(s12, SA12, n02, name+1)
    for i in range(n02) :
      s12[SA12[i]] = i+1
  else :
    for i in range(n02) :
      SA12[s12[i]-1] = i

  s0 = array('i',[SA12[i]*3 for i in range(n02) if SA12[i]<n0])
  radixpass(s0, SA0, s, n0, K)
  
  p = j = k = 0
  t = n0 - n1
  while k < n :
    i = SA12[t]*3+1 if SA12[t]<n0 else (SA12[t] - n0)*3 + 2
    j = SA0[p] if p < n0 else 0

    if SA12[t] < n0 :
      test = (s12[SA12[t]+n0] <= s12[j//3]) if(s[i]==s[j]) else (s[i] < s[j])
    elif(s[i]==s[j]) :
      test = s12[SA12[t]-n0+1] <= s12[j//3 + n0] if(s[i+1]==s[j+1]) else s[i+1] < s[j+1]
    else :
      test = s[i] < s[j]

    if(test) :
      SA[k] = i
      t += 1
      if t == n02 :
        k += 1
        while p < n0 :
          SA[k] = SA0[p]
          p += 1
          k += 1
        
    else : 
      SA[k] = j
      p += 1
      if p == n0 :
        k += 1
        while t < n02 :
          SA[k] = (SA12[t] * 3) + 1 if SA12[t] < n0 else ((SA12[t] - n0) * 3) + 2
          t += 1
          k += 1
    k += 1

def suffixArray(s) -> array:
  k = len(ALPHABET)
  n = len(s)
  t = {c: i for i, c in enumerate(ALPHABET)}
  SA = array('i', [0]*(n+3))
  kark_sort(array('i', [t[c] for c in s]+[0]*3), SA, n, k)
  return SA[:n]

# ----------------------------------------------------------------------
# FM adapted from Langmead's notebooks.

def bwtFromSa(t: str, sa: array) -> str:
    return ''.join(t[i-1] if i else '$' for i in sa)

class FmCheckpoints(object):
    ''' Manages rank checkpoints and handles rank queries, which are
        O(1) time, with the checkpoints taking O(m) space, where m is
        length of text. '''
    
    def __init__(self, bw: str, cpIval: int = 4):
        ''' Scan BWT, creating periodic checkpoints as we go '''
        self.cps = {c: [] for c in ALPHABET}  # checkpoints
        self.cpIval = cpIval                  # spacing between checkpoints
        tally = {c: 0 for c in ALPHABET}      # tally so far

        for i, c in enumerate(bw):
            tally[c] += 1 # up to *and including*
            if i % cpIval == 0:
                for c in tally.keys():
                    self.cps[c].append(tally[c])
    
    def rank(self, bw: str, c: str, row: int) -> int:
        ''' Return # c's there are in bw up to and including row '''
        if row < 0:
            return 0
        i, nocc = row, 0
        # Always walk to left (up) when calculating rank
        while (i % self.cpIval) != 0:
            if bw[i] == c:
                nocc += 1
            i -= 1
        return self.cps[c][i // self.cpIval] + nocc


class FmIndex():
    ''' O(m) size FM Index, where checkpoints and suffix array samples are
        spaced O(1) elements apart.  Queries like count() and range() are
        O(n) where n is the length of the query.  Finding all k
        occurrences of a length-n query string takes O(n + k) time.
        
        Note: The spacings in the suffix array sample and checkpoints can
        be chosen differently to achieve different bounds. '''
        
    def __init__(self, t: str, cpIval: int = 4):
        self.ssa = suffixArray(t)
        self.bwt = bwtFromSa(t, self.ssa)
        self.slen = len(self.bwt)
        
        self.cps = FmCheckpoints(self.bwt, cpIval)
        tots = Counter(self.bwt)
        # Calculate concise representation of first column
        self.first = {}
        totc = 0
        for c, count in sorted(tots.items()):
            self.first[c] = totc
            totc += count
    
    def count(self, c: str) -> int:
        ''' Return number of occurrences of characters < c '''
        return self.first[c]
    
    def range(self, p: str) -> tuple[int, int]:
        ''' Return range of BWM rows having p as a prefix '''
        l, r = 0, self.slen - 1 # closed (inclusive) interval
        for i in range(len(p)-1, -1, -1): # from right to left
            l = self.cps.rank(self.bwt, p[i], l-1) + self.count(p[i])
            r = self.cps.rank(self.bwt, p[i], r)   + self.count(p[i]) - 1
            if r < l:
                break
        return l, r+1
    
    def resolve(self, row: int) -> int:
        ''' Given BWM row, return its offset w/r/t T '''
        return self.ssa[row]
    
    def hasSubstring(self, p: str) -> bool:
        ''' Return true if and only if p is substring of indexed text '''
        l, r = self.range(p)
        return r > l
    
    def hasSuffix(self, p: str) -> bool:
        ''' Return true if and only if p is suffix of indexed text '''
        l, r = self.range(p)
        off = self.resolve(l)
        return r > l and off + len(p) == self.slen - 1
    
    def occurrences(self, p: str) -> list[int]:
        ''' Return offsets for all occurrences of p, in no particular order '''
        l, r = self.range(p)
        return [self.resolve(x) for x in range(l, r)]
