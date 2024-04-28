# GWFA-TB

GWFA (Graph WaveFront Alignment) is an algorithm designed to align a sequence with a sequence graph. It is an adaptation of the [WFA algorithm][wfa] algorithm for graphs and this repository presents an improved implementation of the [original GWFA *proof-of-concept*][gwfa]. This implementation supports alignment backtracking and computes the edit distance between a sequence and graph, reporting the complete traceback in the standard GAF format.

It's important to note that the algorithm assumes the start of the sequence to be aligned with the start of the first segment in the graph and requires the query sequence to be fully aligned. This makes it similar to the SHW mode of [edlib][edlib]. However, it does not support semi-global alignment and is therefore not intended for read mapping.

GWFA is optimized for graphs consisting of long segments and like WFA, it is fast when the edit distance is small.

## Usage

```sh
git clone https://github.com/sgarlata/gwfa-tb.git
cd gwfa-tb && make
./gwfa-tb test/C4-90.gfa.gz test/C4-NA19240.1.fa.gz
./gwfa-tb test/C4-NA19240.1.fa.gz test/C4-NA19240.2.fa.gz
```

Output GAF files are named in the following format:
`gwfa-tb_2023-12-19_17:22:00.gaf`.

[gwfa]: https://github.com/lh3/gwfa
[wfa]: https://github.com/smarco/WFA
[edlib]: https://github.com/Martinsos/edlib
