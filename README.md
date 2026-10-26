# Read Mapping

## From a genome to an index

In Exercise 05 you built a genome from reads without a reference. Today the genome
*is* the reference: you build the index from the lecture yourselves and place reads on
it. Each task reuses the previous one.

Structure of the session:
1. **Suffix array** (theory) → **Task 1**
2. **Burrows-Wheeler transform** (theory) → **Task 2**
3. **FM-index** (theory) → **Task 3**
4. **Mapping reads** (theory) → **Task 4** (code only; finish at home if needed)

The paper parts of Tasks 1–3 (hand traces on a small toy text, and the questions that
go with them) are on the printed handout. This README has the code.

**Where your answers go.** Write the answers to the paper tasks and to the handout's
*Think about it* questions on the handout. Write your own code from scratch in a single
file, `mapping_exercise.py`, and answer the *Think about it* questions of this README in
a comment at the end of that file. Fork this repository to your own GitHub account,
clone your fork, and commit and push your work to it as you go (e.g. one commit per
finished task).

**How you check yourself.**
Every function has self-checks (properties that must hold for any input), and from
Task 1 on you have a slow but obviously correct function, `brute_force`, to compare with.
Run the self-checks on **many random texts**, not only on the examples.

Reading: Jones & Pevzner, *An Introduction to Bioinformatics Algorithms*, Chapter 9
(suffix trees and pattern matching); Compeau & Pevzner, *Bioinformatics Algorithms: An
Active Learning Approach*, 2nd ed., Vol. II, Chapter 9 ("How Do We Locate
Disease-Causing Mutations?"): suffix arrays, the Burrows-Wheeler transform, backward
search, the sampled suffix array.

### Words used below

- Every text ends with **`$`**, a character that occurs nowhere else. In Python
  `"$" < "A"`, so `$` sorts before every DNA letter.
- **Suffix** `i` of `text` = `text[i:]`. Positions are **0-based**, as in Python.
- **Row**: the position of a suffix in the alphabetically sorted list of all suffixes
  (row 0 is always the suffix `"$"`). In code a row is simply an index: `sa[row]`,
  `bwt[row]` and `counts[X][row]` are the columns of that sorted table.
- **Half-open interval** `[a, b)`: the rows a, a + 1, …, b − 1, exactly like
  `range(a, b)` or `sa[a:b]`. It is empty when a = b.
- **Block** of a pattern: the rows whose suffixes start with that pattern.
- **F** and **L**: the first and the last column of the sorted rotation matrix (the
  lecture's First and Last column). L is the BWT; F is its letters sorted. In code, L
  is the string `bwt` (`self.bwt` in Task 3).
- **Occurrence number** k of a letter in L: how many equal letters stand above it in L,
  + 1 (the lecture writes A1, A2, …). It is not the row.
- **FirstOccurrence(X)** and **Count<sub>X</sub>(i)** of the lecture and the handout are
  `first_occurrence[X]` and `counts[X][i]` in your code.
- **LF** (Last-to-First): from a row to the row of the suffix that starts one position
  earlier in the text.
- **Random texts** for self-checks: write a helper such as
  `"".join(random.choice("ACGT") for _ in range(n)) + "$"` (`import random`). Use short
  texts (n from 1 to about 40), and also an alphabet of only `"AC"`: it gives many
  repeats, which is where index code goes wrong.

---

## Task 1: Suffix array and binary search

First do the paper part of Task 1 on the handout.

#### `brute_force(text, p)`
- **Input:** `text`, `p`, strings.
- **Output:** `list[int]`, every start position of `p` in `text`, in increasing order,
  overlapping occurrences included.
- **Example:** `brute_force("AAAA", "AA")` → `[0, 1, 2]`.
- **Hint:** write it first. It is the reference you compare everything else with today.

#### `suffix_array(text)`
- **Input:** `text`, a string ending with `$`.
- **Output:** `list[int]`, the start positions of the suffixes of `text` in alphabetical
  order of the suffixes.
- **Hint:** one call of `sorted` with a `key`. This copies every suffix, so it needs
  memory proportional to n², fine for a few thousand letters, not for a human genome.
- **Self-checks:**
  1. `sorted(sa) == list(range(len(text)))`: every position exactly once.
  2. Neighbouring suffixes increase: `text[sa[i]:] < text[sa[i + 1]:]` for every i.

#### `find_block(text, sa, p)`
- **Input:** `text`; `sa`, its suffix array; `p`, a non-empty pattern.
- **Output:** `tuple[int, int]`, the half-open interval `(a, b)` of rows whose suffixes
  start with `p`. If `p` does not occur, `a == b` (an empty block at the row where `p`
  would be).
- **Hint:** two binary searches over `lo, hi = 0, len(sa)`. Both compare only the first
  `len(p)` letters of the suffix with `p`: the start `a` is the first row where they are
  ≥ `p`, the end `b` the first row where they are > `p`. Do not scan the rows one by
  one: the point is about log₂ n steps.
- **Self-check:** `sorted(sa[a:b]) == brute_force(text, p)` for hundreds of random texts
  and random patterns of length 1–5, including patterns that do not occur and patterns
  longer than some suffixes.

---

## Task 2: The Burrows-Wheeler transform and its inverse

First do the paper part of Task 2 on the handout.

#### `bwt(text, sa)`
- **Input:** `text`, a string ending with `$`; `sa`, its suffix array.
- **Output:** `str`, the BWT: the letter **before** each suffix, in row order,
  `text[sa[row] - 1]`.
- **Example:** `bwt("CAT$", [3, 1, 0, 2])` → `"TC$A"`.
- **Hint:** for `sa[row] == 0`, Python's `text[-1]` is already the `$` you need. No
  rotation matrix in code.
- **Self-checks:** the same letters as `text` (`sorted(bwt(text, sa)) == sorted(text)`);
  exactly one `$`.

#### `first_occurrence(letters)`
- **Input:** `letters`, a string (here, always a BWT).
- **Output:** `dict[str, int]`, for each letter X of `letters`, the row where X starts in
  F = the number of letters in `letters` that are smaller than X.
- **Example:** `first_occurrence("TC$A")` → `{"$": 0, "A": 1, "C": 2, "T": 3}`.
- **Self-check:** the value of `$` is 0, and the values increase in alphabetical order
  of the letters.

#### `inverse_bwt(bwt)`
- **Input:** `bwt`, the BWT (L) of some text ending with `$`.
- **Output:** `str`, that text.
- **Hint:** First-Last property: the k-th X in L is the k-th X in F. First compute, in
  one pass, the occurrence number k of every letter of `bwt`. Then start in row 0 (the
  row of `"$"`) and repeat: X = `bwt[row]` with occurrence number k; collect X; the next
  row is `first_occurrence(bwt)[X] + k - 1`. You collect the text backwards. Stop after
  `len(bwt) - 1` letters, reverse them, and append the `$`.
- **Self-check:** `inverse_bwt(bwt(t, suffix_array(t))) == t` for hundreds of random
  texts.

*Think about it:* count the **runs** (maximal blocks of equal neighbouring letters; e.g.
`AACCCA` has 3) in the reference text (`reference.fasta` + `"$"`), in its BWT, and in
the BWT of the reference written **twice** in a row (+ `"$"`). What happens to each
number when the text is doubled, and why?

---

## Task 3: The `FMIndex` class

First do the paper part of Task 3 on the handout.

The index is built once and then used for every read: that is why it is a class (just as
`MotifFinder` in Exercise 03 kept its precomputed windows).

#### `FMIndex(text, s)`
- **Input:** `text`, a string ending with `$`; `s`, a positive int (the sampling step).
- **Stores:**
  - `self.bwt`: the BWT of `text`;
  - `self.first_occurrence`: `first_occurrence(self.bwt)`;
  - `self.counts`: `dict[str, list[int]]`, for every letter X a list of length
    `len(text) + 1` with `counts[X][i]` = the number of X in `self.bwt[:i]`;
  - `self.sampled`: `dict[int, int]`, part of the suffix array: row → `sa[row]`, kept
    only where `sa[row]` is a multiple of `s`;
  - `self.s`: `s`.
- **Hint:** compute the full suffix array inside `__init__` (Task 1), use it for the
  BWT and `sampled`, and do not keep it. Build all `counts` lists in one pass over the
  BWT.
- **Example:** `FMIndex("CAT$", 2)` stores `bwt = "TC$A"`, `counts["T"] = [0, 1, 1, 1, 1]`,
  `sampled = {2: 0, 3: 2}`.

#### `lf(row)`
- **Input:** `row`, an int with 0 ≤ row < `len(self.bwt)`.
- **Output:** `int`, the row of the suffix that starts one position **earlier** in the
  text: `first_occurrence[X] + counts[X][row]` with X = `self.bwt[row]`.
- **Example:** for `FMIndex("CAT$", 2)`, `lf(1)` → `2` (row 1 is `AT$`, row 2 is
  `CAT$`).

#### `count(p)`
- **Input:** `p`, a non-empty string.
- **Output:** `tuple[int, int]`, the interval `(top, bottom)` of rows that start with
  `p`, found by **backward search**: start with `(0, len(self.bwt))`, then for each
  letter X of `p` from the last one: `top = first_occurrence[X] + counts[X][top]`, and
  the same for `bottom`.
- **Hint:** stop as soon as the interval is empty (any `(t, t)` is a correct answer
  then). A letter that does not occur in the text gives an empty interval, not a
  `KeyError`.
- **Example:** for `FMIndex("CAT$", 2)`, `count("AT")` → `(1, 2)`.

#### `locate(p)`
- **Input:** `p`, a non-empty string.
- **Output:** `list[int]`, the text positions of all occurrences of `p` (one per row of
  `count(p)`, in row order).
- **Hint:** from each row, apply `lf` until you reach a row in `self.sampled`; the
  position is the stored value + the number of steps you took.
- **Example:** for `FMIndex("CAT$", 2)`, `locate("AT")` → `[1]` (row 1 → row 2, stored
  0, one step).
- **Self-checks** (on hundreds of random texts and patterns):
  1. `bottom - top == len(locate(p)) == len(brute_force(text, p))`.
  2. `sorted(locate(p)) == brute_force(text, p)`, for s = 1, 4 and 32.
  3. If `p` occurs, `count(p)` gives the same interval as `find_block` from Task 1; if
     it does not, both intervals are empty.
  4. Every stored value is a multiple of `s`, and there is one for each such position:
     `sorted(self.sampled.values()) == list(range(0, len(text), s))`.
  5. *(optional)* Count the LF steps inside `locate`: there are always fewer than `s`.

*Think about it:* how many integers do your `counts` hold for the reference (one list
per letter, `$` included: 5 lists), compared with the number of letters in its BWT?
How would storing `counts` only at every K-th row (checkpoints, as in the lecture)
change the memory and the work per `lf`?

*(Optional)* Implement the checkpoints: store `counts[X][i]` only for i divisible by K
and count the rest directly in the BWT. All self-checks must still pass.

---

## Task 4: Mapping reads

`reference.fasta` holds one genome (the one from Exercise 05). `reads.fasta` holds
315 reads of 100 bp. Reads come from **either strand**, have a few sequencing errors
(substitutions only, no insertions or deletions), and some come from somewhere else
entirely. For every read, find **all** places on either strand where it fits with at
most d = 2 mismatches, or report it as unmapped.

**Pigeonhole (lecture):** cut the read into d + 1 pieces (seeds). A read with at most
d mismatches has at least one seed without a mismatch, so that seed occurs exactly in
the reference.

#### `reverse_complement(seq)`
- Copy yours from Exercise 02.
- **Example and self-check:** `reverse_complement("AACG")` → `"CGTT"` (reversed *and*
  complemented); applying it twice returns the input.

#### `hamming(a, b)`
- **Input:** two strings of the same length.
- **Output:** `int`, the number of positions where they differ.
- **Example:** `hamming("ACGT", "AGGA")` → `2`.
- **Hint:** check that the lengths are equal (`assert`): `zip` silently stops at the
  shorter string, which hides bugs in Task 4.

#### `seeds(read, n)`
- **Input:** `read`, a string; `n`, the number of seeds.
- **Output:** `list[tuple[int, str]]`, `n` consecutive, nearly equal pieces that cover
  the whole read, each with its offset in the read: piece i starts at
  `i * len(read) // n`.
- **Example:** `seeds("ACGTACGTAC", 3)` → `[(0, "ACG"), (3, "TAC"), (6, "GTAC")]`.

#### `map_read(read, fm, ref, d=2, n_seeds=None)`
- **Input:** `read`; `fm`, an `FMIndex` of `ref + "$"`; `ref`, the reference string;
  `d`, the maximum number of mismatches; `n_seeds`, the number of seeds (`None` means
  d + 1).
- **Output:** `tuple[list[tuple[str, int, int]], int]`: the sorted list of hits
  `(strand, position, mismatches)`, strand `"+"` for the read and `"-"` for its reverse
  complement, position 0-based in `ref` (for `"-"`: the leftmost position where the
  reverse complement matches); and the number of **candidates**: the distinct starts
  inside the reference that you verified with `hamming`, summed over both strands.
- **Hint:** for each strand separately: for every seed `(o, seed)`, every position p in
  `fm.locate(seed)` gives a candidate start p − o. Collect the candidates in a `set`
  (several seeds often vote for the same start), drop starts that would make the read
  stick out of the reference, and verify each one with `hamming` on the whole read.

#### `brute_map(read, ref, d=2)`
- **Input:** `read`, `ref` and `d` as in `map_read` (no index).
- **Output:** the sorted list of hits, as in `map_read`, found the slow way: try every
  start on both strands.

### Map all reads

1. Read both files with Biopython (`from Bio import SeqIO`; `SeqIO.read` for the one
   reference, `SeqIO.parse` for the reads; `str(record.seq)` for the sequence).
2. Build `FMIndex(ref + "$", 32)` once.
3. Map every read with `map_read`. Print how many reads are mapped, how many are
   unmapped, and how many have more than one hit, and print the hits of the first five
   reads (name, strand, position, mismatches).
4. **Self-checks:**
   1. For every hit, 0 ≤ position ≤ `len(ref) - 100`, and recomputing the mismatches of
      the read (or its reverse complement) against `ref[position:position + 100]` gives
      the same number, at most d.
   2. For 30 reads chosen with `random.sample`, `map_read` gives exactly the same hits
      as `brute_map` (brute force takes a fraction of a second per read).
   3. No unmapped read has a hit in `brute_map` (check them all).
   4. The ends of the reference: `ref[:100]` maps to `("+", 0, 0)`, and the reverse
      complement of `ref[-100:]` maps to `("-", len(ref) - 100, 0)`.

*Think about it:* run the mapping again with `n_seeds` = 10, 12 and 16. How long are
the seeds, how many candidates per read (total candidates / number of reads) does each
give, and do the hits change? Why is a seed of 6 letters a bad seed for this genome,
and what would happen with a human genome?

*Think about it:* this genome contains two different 12-bp sequences that each occur
twice. Does that make any read map to two places? Why or why not? How long would a
repeat have to be to cause it?
