import argparse
import time

from Bio import SeqIO
from typing import Any

from bwt import FmIndex
from kEditDp import kEditDp

K = 110             
MAX_DIFF = 20        
PART_K = 34            
STEP = 5          


def main():
    # Reading data
    reference_file, reads_file, output_file = read_args()
    text, reads = read_data(reference_file, reads_file)

    t0 = time.time()
    fm = bulid_fm_index(text)
    
    t1 = time.time()
    lines = align_all(reads, text, fm)

    t2 = time.time()
    with open(output_file, "w") as f:
        f.writelines(lines)

    print(f"Output written to {output_file}.")
    print(f"Total time: {t2-t0:.2f}s.")
    print(f"Preprocessing time: {t1-t0:.2f}s.")
    print(f"Alignment time: {t2-t1:.2f}s.")


def read_args() -> tuple[str, str, str]:
    parser = argparse.ArgumentParser(description="Map reads to reference using FM-index.")
    parser.add_argument("reference", type=str, help="Path to reference file.")
    parser.add_argument("reads", type=str, help="Path to reads file.")
    parser.add_argument("output", type=str, help="Path to output file.")
    args = parser.parse_args()
    return args.reference, args.reads, args.output


def read_data(reference: str, reads: str) -> tuple[str, Any]:
    seq_rec_list = [seq_record for seq_record in SeqIO.parse(reference, "fasta")]
    text = seq_rec_list[0].seq + "$"
    del seq_rec_list
    return text, SeqIO.parse(reads, "fasta")


def bulid_fm_index(text: str) -> FmIndex:
    return FmIndex(text)


def align_all(reads: Any, text: str, fm_index: FmIndex) -> list[str]:
    lines = []
    i = 0
    for read in reads:
        r = str(read.seq)
        hits = align(r, text, fm_index)
        if hits is not None:
            lines.append("{}\t{}\t{}\n".format(read.id, hits[1], hits[2]))
        i += 1
    return lines

    
def align(pattern: str, text: str, fm_index: FmIndex) -> tuple[int, int, int] | None:
    intervals, part_to_pat = process_seeds(pattern, fm_index)
    
    if not intervals:
        return None

    mems = set()
    text_scores = {}
    for part, inter in intervals.items():
        for it in inter:
            for p in part_to_pat[part]:
                score, mem, start, end = extend_seed(pattern, text, p, it)
                if mem not in mems:
                    text_scores[(start, end)] = text_scores.get((start, end), 0) + score
                    mems.add(mem)
    
    text_scores = text_scores.items()

    scores = []
    for (start, end), _ in text_scores:
        res = local_align(pattern, text, start, end)
        if res is not None:
            scores.append(res)

    if not scores:
        return None

    return min(scores)


def process_seeds(pattern: str, fm_index: FmIndex) -> tuple[Any, Any]:
    parts, part_to_pat = eff_part_overlap(pattern, STEP)
    intervals = find_offsets(fm_index, parts)
    return intervals, part_to_pat 


def eff_part_overlap(pattern: str, step: int = 1
                    ) -> tuple[list[str], dict[str, list[tuple[int, int]]]]:
    m = len(pattern)
    q = m // (PART_K+1)
    part_to_pat = {}
    parts = [pattern[i:i+q] for i in range(0, m-q+1, step)]

    for i in range(0, m-q+1, step):
        if pattern[i:i+q] in part_to_pat:
            part_to_pat[pattern[i:i+q]].append((i,i+q))
        else:
            part_to_pat[pattern[i:i+q]] = [(i,i+q)]
            
    return parts, part_to_pat


def find_offsets(fm_index: FmIndex, kmers: list[str]) -> dict[str, list[tuple[int, int]]]:
    occurences = {part:fm_index.occurrences(part) for part in kmers if fm_index.hasSubstring(part)}
    intervals = {}
    for part, occ in occurences.items():
        m = len(part)
        for ind in occ:
            if part in intervals:
                intervals[part].append((ind, ind + m))
            else:
                intervals[part] = [(ind, ind + m)]
                
    return intervals


def extend_seed(
        pattern: str, text: str, 
        interval_p: tuple[int, int], 
        interval_t: tuple[int, int]
        ) -> tuple[int, tuple[int, int, int, int], int, int]:
    p1, p2, t1, t2 = extend_to_mem(pattern, text, interval_p, interval_t)
    m, n = len(pattern), len(text)
    start = max(t1 - p1 - MAX_DIFF//2, 0)
    end = min(t2 + m - p2 + MAX_DIFF//2, n)
    
    return p2 - p1, (p1, p2, t1, t2), start, end


def local_align(pattern: str, text: str, start: int, end: int) -> tuple[int, int, int] | None:
    res = kEditDp(pattern, text[start:end], MAX_DIFF, K)
    if res is None:
        return None

    score, soff, eoff = res
    return score, start + soff, start + eoff


def extend_to_mem(
        pattern: str, text: str, 
        interval_p: tuple[int, int], 
        interval_t: tuple[int, int]
        ) -> tuple[int, int, int, int]:
    m, n = len(pattern), len(text)
    p1, p2 = interval_p
    t1, t2 = interval_t

    while p1 > 0 and t1 > 0 and text[t1-1] == pattern[p1-1]:
        t1 -=1
        p1 -=1

    while p2 < m-1 and t2 < n-1 and text[t2+1] == pattern[p2+1]:
        t2 +=1
        p2 +=1
    
    return p1, p2 + 1, t1, t2 + 1

    
if __name__ == '__main__':
    main()
