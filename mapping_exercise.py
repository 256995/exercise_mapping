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
    return SA


def find_block(text, sa, p):
    lo, hi = 0, len(sa)
    # binary search for low
    while hi > lo:
        mid = int(lo + (hi-1-lo)/2)
        if sa[mid][1][:len(p)] < p:
            lo = mid + 1
        else:
            hi = mid
    low = lo
    # binary search for high
    hi = len(sa)
    while hi > lo:
        mid = int(lo + (hi - 1 - lo) / 2)
        if sa[mid][1][:len(p)] > p:
            hi = mid
        else:
            lo = mid + 1
    high = hi
    return low, high


if __name__ == "__main__":
    print(brute_force('AAAA', 'AA'))
    SA = suffix_array('AAGGCGGCTTCA$')
    print([row[0] for row in SA])
    block_of_GGC = find_block('AAGGCGGCTTCA$', SA, 'GGC')
    print(block_of_GGC)