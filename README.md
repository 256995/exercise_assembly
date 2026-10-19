# De Novo Genome Assembly

## From graphs to genomes

In the lecture you saw two ways to turn reads into a graph. In the **overlap graph**,
the assembler looks for a path through every *node* (Hamiltonian). In the **de Bruijn
graph**, it looks for a path through every *edge* (Eulerian). In Exercise 04 you decided
whether an Eulerian path *exists*. Today you find one, read a genome off it, and see
where that goes wrong.

Structure of the session:
1. **The overlap graph** → **Task 1** (paper only)
2. **The de Bruijn graph** → **Task 2**
3. **Hierholzer's algorithm** (theory) → **Task 3**
4. **Repeats and contigs** (theory) → **Task 4**
5. **A mystery genome** → **Task 5** (code only)

The paper parts of Tasks 1–4 (drawing graphs, hand traces, and the questions that go
with them) are on the printed handout. The handout also gives the three toy examples
A, B and C. This README has the code. Parts marked *(optional)* are for students who
finish early.

**Where your answers go.** Write the answers to the paper tasks and to the handout's
*Think about it* questions on the handout. Write your own code from scratch in a single
file, `assembly_exercise.py`. Fork this
repository to your own GitHub account, clone your fork, and commit and push your work
to it as you go (e.g. one commit per finished task).

Assumptions for the whole exercise (the "ideal data" of the lecture): reads have **no
sequencing errors**, they all come from the **same strand**, and every k-mer of the
genome is covered by at least one read.

Reading: Compeau & Pevzner, *Bioinformatics Algorithms: An Active Learning Approach*,
2nd ed., Vol. I, Chapter 3 ("How Do We Assemble Genomes?"): the de Bruijn graph, Euler's
theorem, the algorithm for an Eulerian cycle used here, and contigs (maximal
non-branching paths).

### Words used below

- **node** = vertex (the lecture and Exercise 04 use both words).
- **prefix** and **suffix** of a k-mer: its first k − 1 and its last k − 1 letters (the
  lecture's *left* and *right* (k−1)-mer). Example: the prefix of `GATT` is `GAT`, its
  suffix is `ATT`.
- **balanced** node: in-degree = out-degree.
- **walk**: a list of nodes in which every two neighbours are joined by an edge. An
  **Eulerian path** is a walk that uses every edge exactly once (the lecture's
  "Eulerian walk"); if it ends where it started, it is an **Eulerian cycle** (circuit).
- **`Counter`** (`from collections import Counter`) counts how often each item occurs:
  `Counter(["AC", "CG", "AC"])` is `{"AC": 2, "CG": 1}`. Two Counters are equal (`==`)
  exactly when they contain the same items the same number of times, in any order.
  The self-checks below use it to compare collections with repeats.

---

## Task 1: The overlap graph

Paper only, on the handout (toy A).

---

## Task 2: The de Bruijn graph

First do the paper part of Task 2 on the handout (the de Bruijn graph of toy A).

### Code

A de Bruijn graph is a **directed** graph, stored as an adjacency list, exactly as in
Exercise 04: a `dict` mapping each node to the list of nodes its edges point to.

**Copy your own Exercise 04 functions** into `assembly_exercise.py`: `in_out_degree`,
`dfs`, `connected_components`, `is_weakly_connected` and `eulerian_type`. You will run
them on de Bruijn graphs (always with `directed=True`). If one of them was not finished
or did not pass its self-checks, fix it first: everything below builds on them.

**No recursion.** Python stops a recursion after about 1 000 nested calls, and the
graphs in Task 5 have thousands of nodes. If your `dfs` from Exercise 04 calls itself,
rewrite it with a loop and an explicit stack (a list) now. The same holds for every
function below.

#### `kmer_composition(text, k)`
- **Input:** `text`, a string; `k`, an int with 1 ≤ k ≤ `len(text)`.
- **Output:** `list[str]`, every substring of length `k` of `text`, in the order they
  occur, repeats included.
- **Example:** `kmer_composition("GATTACA", 4)` → `["GATT", "ATTA", "TTAC", "TACA"]`.
- **Self-check:** the list has exactly `len(text) - k + 1` entries.

#### `debruijn_graph(kmers)`
- **Input:** `kmers`, a list of strings, all of the same length k ≥ 2.
- **Output:** `dict[str, list[str]]`, the de Bruijn graph. Every k-mer is one edge from
  its prefix to its suffix. Every prefix and every suffix is a key of the dict, even if
  its list ends up empty.
- **Hint:** the same k-mer twice means two parallel edges: two entries in the list, as
  for the multigraphs of Exercise 04. Do not remove duplicates here.
- **Self-checks:**
  1. The total length of all neighbour lists equals `len(kmers)`: one edge per k-mer.
  2. The keys are exactly the prefixes and suffixes:
     `set(graph) == {km[:-1] for km in kmers} | {km[1:] for km in kmers}`.

#### `spell_path(nodes)`
- **Input:** `nodes`, a non-empty list of strings of length k − 1, where each one
  overlaps the next by k − 2 letters (a walk in a de Bruijn graph).
- **Output:** `str`, the string this walk spells: the first node, followed by the last
  letter of every later node.
- **Example:** `spell_path(["GAT", "ATT", "TTA"])` → `"GATTA"`.
- **Self-check:** for any string `text` and any `m` with 1 ≤ m ≤ `len(text)`,
  `spell_path(kmer_composition(text, m)) == text`. Try it on a few strings of your own.

Now build the de Bruijn graph of toy A's 3-mers (from the handout) in code and compare
it with your drawing: the same nodes, the same edges, the same in- and out-degrees, and
the same `eulerian_type` verdict as by hand. If anything disagrees, recheck both sides.

**The toys have no repeated k-mer, so they never test parallel edges.** Run all
self-checks of this task once more on the 3-mers of a string you make up that contains
a repeated 3-mer (e.g. a short tandem repeat such as `CATCATCAT`).

*Think about it:* take **any** string, build the de Bruijn graph of its k-mer
composition and run `eulerian_type` on it. Can the answer ever be `"none"`? When exactly
is it `"circuit"` rather than `"path"`? Make up strings to test your answer.

---

## Task 3: Hierholzer's algorithm

`eulerian_type` tells you whether an Eulerian path exists. Now you find one. The
algorithm from today's slides, with fixed rules so that a hand trace and your code
always agree step by step:

1. From the start node, walk along **unused** edges. Whenever there is a choice, take
   the unused edge whose end node comes **first alphabetically**. Mark each edge as used
   when you walk it. Stop when the current node has no unused edge left.
2. If every edge is used, you are done.
3. Otherwise go along the **current walk** (walk 1 with every earlier walk already
   spliced in) from its first node, and stop at the **first position** whose node still
   has an unused edge. Do not pick by alphabet, and do not look only at walk 1. From
   that node, walk by the same rule as in step 1 until you are stuck.
4. **Splice:** the new walk starts and ends at that node. Replace the node at **that
   position** (only that one occurrence) by the whole new walk. Go back to step 2.

First do the hand trace of toy B on the handout, and **show it to the instructor before
you start coding this task.**

The version of this algorithm in most code online (e.g. Rosalind solutions) or from a chatbot uses
a stack and builds the path backwards. It also returns a valid Eulerian path, but often
a **different** one, so it fails self-check 5 below. Implement the rules above.

#### `hierholzer(graph, start, verbose=False)`
- **Input:** `graph`, a de Bruijn graph for which `eulerian_type` is not `"none"`;
  `start`, the node to start from (see `eulerian_path` for which node that must be);
  `verbose`, a bool.
- **Output:** `list[str]`, the nodes of the final walk, in order, using every edge
  exactly once. If `verbose` is true, also print every walk (walk 1, then each new walk
  with the node it restarted at) and the walk after each splice.
- **Ties:** exactly the rules above.
- **Hint:** do not change `graph` itself. Keep a separate `dict` of unused edges, for
  example a sorted copy of every neighbour list, and remove an edge from it when you
  use it. Use loops, not recursion.

#### `eulerian_path(graph, verbose=False)`
- **Input:** `graph`, a de Bruijn graph for which `eulerian_type` is not `"none"`;
  `verbose`, a bool, passed on to `hierholzer`.
- **Output:** `list[str]`, the walk returned by `hierholzer(graph, start, verbose)`. The `start`
  is the node with out − in = +1 if there is one (`eulerian_type` is `"path"`; any other
  start would miss edges); otherwise every node is balanced and `start` is the
  alphabetically first node, `min(graph)`.
- **Self-checks** (for any graph, not only the toys):
  1. The walk has exactly (number of edges) + 1 nodes.
  2. Every edge is used exactly once: the `Counter` of consecutive pairs
     `zip(walk, walk[1:])` equals the `Counter` of all edges `(u, v)` of `graph`.
  3. **The key invariant of the whole exercise:** the string spelled by the walk has
     the same k-mer composition as the k-mers you built the graph from (compare the two
     `Counter`s). It must hold for every Eulerian path, whether or not it is the true
     genome.
  4. `graph` is unchanged after the call (compare it with a copy taken before).
  5. On toy B, `verbose=True` prints exactly the walks, restart nodes and splices of
     your hand trace.
  6. The result does not depend on the order of the input: `eulerian_path` of the
     graph built from `kmers` and of the graph built from `kmers[::-1]` (reversed list)
     are equal. (If not, you forgot to sort somewhere.)

*Think about it:* in step 1, starting from the node with out − in = +1, why can the
first walk only get stuck at the node with in − out = +1? And why does every later walk
(step 3) get stuck exactly at the node where it started?

---

## Task 4: Repeats and contigs

First do the paper part of Task 4 on the handout (toy C). Then run `eulerian_path` on
toy C with `verbose=True`: it must return the genome you predicted on the handout
(*Think about it*, "which genome would Hierholzer's algorithm return?"), and each restart
must be the first node on the current walk that still has an unused edge.

Hierholzer always returns *one* Eulerian path. If there are several, it cannot tell
which one is the real genome. An assembler therefore reports only the pieces it can be
sure of: **contigs**, the maximal non-branching paths.

A node is **1-in-1-out** if its in-degree and its out-degree are both 1. A path is
**non-branching** if every node *inside* it (not its first or last node) is 1-in-1-out.
It is **maximal** if it cannot be extended at either end without becoming branching.

#### `contigs(graph)`
- **Input:** `graph`, a de Bruijn graph.
- **Output:** `list[str]`, the spelled strings of all maximal non-branching paths,
  sorted.
- **Hint:** start at every node that is **not** 1-in-1-out. For each of its outgoing
  edges, follow edges as long as the node you reach is 1-in-1-out; the path ends at the
  first node that is not (include that node). Parallel edges give the same contig
  twice: keep both. Afterwards, some edges may still be left: they form cycles in which
  every node is 1-in-1-out (nothing enters or leaves them). Each such cycle is one more
  contig. Spell it as a closed walk that starts and ends at its alphabetically first
  node (the first node repeated at the end, as `spell_path` expects), so the string
  repeats its first k − 1 letters at its end.
- **Self-checks** (k = length of a node + 1):
  1. Every edge of `graph` lies in exactly one contig. For each contig, build the de
     Bruijn graph of its k-mers and take the `Counter` of its edges. Add these Counters
     together. The sum must equal the `Counter` of the edges of `graph`.
  2. Make up a string in which no (k − 1)-mer occurs twice. The de Bruijn graph of its
     k-mers must give exactly one contig: the string itself.
  3. Make up a string whose first k − 1 letters equal its last k − 1 letters and in
     which no other (k − 1)-mer repeats (e.g. `ACGTTGCAAC` with k = 3). It must give
     exactly one contig, of the same length (a rotation of your string).
  4. On toy C, your contigs are the ones you marked on the handout.
  5. The result does not depend on the order of the input k-mers (as in Task 3).

*Think about it:* why must a contig stop at a node with two outgoing edges, even if one
of the two continuations is the true genome?

*(Optional)* Design your own toy genome (k = 3, 12–16 letters) for which the rules of
Task 3 return a genome **different** from the one you made up. Check it with your code
and explain on the handout why it happens.

---

## Task 5: A mystery genome

`reads.fasta` contains error-free reads of length 100 from one small genome, all from
the same strand, with high coverage. Assemble it and find out what it is.

**Before you run anything,** write on the handout your prediction: will a larger k give
fewer or more contigs, and why?

1. Read the file with Biopython (`from Bio import SeqIO`, `SeqIO.parse(..., "fasta")`;
   use `str(record.seq)` for the sequence).
2. For each k in **9, 12, 13, 14, 21**:
   - collect the **distinct** k-mers of all reads in a `set` (you will see why in
     *Think about it* below), and pass `sorted(...)` of it to `debruijn_graph`;
   - record |V| (number of nodes), |E| (number of edges) and `eulerian_type`;
   - compute `contigs`: how many, and how long is the longest;
   - if `eulerian_type` is not `"none"`, run `eulerian_path` and spell the walk.
3. **Self-checks** for every k, on the strings exactly as your functions return them:
   the set of k-mers of every spelled string equals your set of distinct k-mers, and
   the spelled string has |E| + k − 1 letters; the edges of all contigs together are
   exactly the edges of the graph (Task 4, self-check 1).
4. **Trim cycles.** A walk or a contig that ends at the node where it started repeats
   its first k − 1 letters at its end. Cut those letters off (only after step 3) before
   you report a length or compare strings.
5. Fill in your own table on the handout (one row per k):

   | k | \|V\| | \|E\| | `eulerian_type` | contigs | longest contig (trimmed) | spelled string (length after trimming, or —) |
   |---|---|---|---|---|---|---|

6. Compare the trimmed spelled strings you got for different k. A circular sequence
   has no fixed starting point, so two strings `s` and `t` describe the same circle
   exactly when `len(s) == len(t)` and `s in t + t`. Do all of your spelled strings
   describe the same circle?
7. **BLAST.** At [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), choose
   *Nucleotide BLAST* (blastn) and keep the default database (*core_nt*). Submit:
   (a) the trimmed single contig of the smallest k at which `contigs` returned exactly
   one contig; (b) the trimmed spelled string of every smaller k at which you got one;
   (c) *(optional)* the contigs of that smaller k. What organism is it, and how long is
   its genome? For every query, look at *Query cover* and at the *Graphic Summary*:
   does **one** alignment cover the whole query, or does it split into several pieces?
   Note: the database stores a circular genome as a linear string that starts at an
   arbitrary point. A correct assembly may therefore split **once**, where the database
   string starts. More breaks than that mean a misassembly.

*Think about it:* the reads are error-free and cover every k-mer, yet for some k
`eulerian_type` says `"none"`. Which assumption behind your answer to the Task 2
*Think about it* breaks once you keep only distinct k-mers?

*Think about it:* if two of your spelled strings in step 6 describe different circles,
at most one can be the genome, yet both pass self-check 3 of Task 3. Which of your
results would you trust at that k, the spelled string or the contigs? What does BLAST
say? What does that tell you about why assemblers output contigs?

*Think about it:* why is `eulerian_type` `"circuit"` from some k on? Where does this
genome "start"?

*Think about it (optional):* why must you remove duplicate k-mers here? Build the graph
once from **all** k-mers of all reads, duplicates included, and run `eulerian_type`.
What does it say, and why?

*Think about it:* look up what the organism BLAST found is used for in sequencing labs.
