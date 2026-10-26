def brute_force(text, p):
    l = len(p)
    result = []
    for i in range(len(text)-l+1):
        if text[i:i+l] == p:
            result.append(i)
    return result


def suffix_array(text):
    SA = []
    for i in range(len(text)):
        SA.append([i, text[i:]])
    SA.sort(key=lambda x: x[1])
    SA = [row[0] for row in SA]
    return SA


def find_block(text, sa, p):
    lo, hi = 0, len(sa)
    # binary search for low
    while hi > lo:
        mid = int(lo + (hi-1-lo)/2)
        idx = sa[mid]
        if text[idx:idx+len(p)] < p:
            lo = mid + 1
        else:
            hi = mid
    low = lo
    # binary search for high
    hi = len(sa)
    while hi > lo:
        mid = int(lo + (hi - 1 - lo) / 2)
        idx = sa[mid]
        if text[idx:idx + len(p)] > p:
            hi = mid
        else:
            lo = mid + 1
    high = hi
    return low, high


def bwt(text, sa):
    BWT_str = ''
    for idx in sa:
        BWT_str += text[idx-1]
    return BWT_str


def first_occurrence(letters):
    sorted_list = [l for l in letters]
    sorted_list.sort()
    sorted = ''
    for l in sorted_list:
        sorted += l
    keys = set(letters)
    first_occ = {}
    for key in keys:
        idx = sorted.index(key)
        first_occ[key] = idx
    return first_occ


def inverse_bwt(bwt):
    first_occ = first_occurrence(bwt)
    k = []
    count = {key: 0 for key in first_occ.keys()}
    for letter in bwt: # preparing the k
        count[letter] += 1
        k.append(count[letter])
    result = ''
    row = 0
    step = 1 # counting steps
    while step <= len(bwt):
        x = bwt[row]
        result += x
        row = first_occ[x] + k[row] - 1
        step += 1
    return result


class FMIndex:
    def __init__(self, text, s):
        sa = suffix_array(text)
        self.bwt = bwt(text, sa)
        self.first_occurrence = first_occurrence(self.bwt)
        self.counts = {key: [0] for key in ['A', 'C', 'G', 'T']}
        for letter in self.bwt:
            for key in self.counts:
                self.counts[key].append(self.counts[key][-1])
                if key == letter:
                    self.counts[key][-1] += 1
        self.sampled = {}
        for row in range(len(sa)):
            if (row % s) == 0:
                self.sampled[row] = sa[row]
        self.s = s



if __name__ == "__main__":
    print(brute_force('AAAA', 'AA'))
    text = 'AAGGCGGCTTCA$'
    SA = suffix_array(text)
    print(SA)
    block_of_GGC = find_block(text, SA, 'GGC')
    print(block_of_GGC)
    BWT1 = bwt(text, SA)
    occurrence = first_occurrence(BWT1)
    print(occurrence)
    print(inverse_bwt(BWT1))

    # --- backward search ---
    print('--- backward search ---')
    toy = FMIndex(text, 4)
    print(toy.bwt)
    print(toy.first_occurrence)
    print(toy.counts)
    print(toy.sampled)