# model_align

## General
Scan sequences for alignment(s) to a specified model
of a DNA sequence.

**Sequences must be defined as character vectors containing only the character
set `ACTUGactug`. Other characters will not cause a crash, but will result in
nonsensical alignments.**

The function takes as an input one or more sequences and a model that the
sequences should be aligned to.

The model contains a set of motifs that need to be matched in a
given order and where we have little constraint as to the distances between
adjacent motifs. This can be achieved by setting the sequence insertion
penalty to 0 at the end of each motif. It would be better to define
an equation or function that is based upon the known distribution of
motif distances, but in order to implement this I would need a better idea of
the relevant distributions. Having a free gap insertion between motifs means that the
functions cannot currently be used to scan long sequences.

The model is specified as a matrix where columns indicate the positions of
the model. Each column contains the scores associated with:

1. A
2. C
3. T
4. G
5. Insertion (nucleotides in the sequence that don't match to the model)
6. Deletion  (nucleotides present in the model but not the sequence)

For example, the model for the consensus:

A,C,A|G,N<sub>x</sub>,T,A,T,A

(where A|G indicates that an A or G are allowed, and N<sub>x</sub> indicates a
sequence of any length separating the two motifs) could be defined by the
following table:

|       |    1 |    2 |    3 |    4 |    5 |    6 |    7 |
|     -:|    -:|    -:|    -:|    -:|    -:|    -:|    -:|
| **A** | **4**| -4   | **2**| -4   | **4**| -4   | **4**|
| **C** | -4   | **4**| -4   | -4   | -4   | -4   | -4   |
| **T** | -4   | -4   | -4   | **4**| -4   | **4**| -4   |
| **G** | -4   | -4   | **2**| -4   | -4   | -4   | -4   |
| _I_   |-10   |-10   | 0    |-10   |-10   |-10   |-10   |
| _D_   |-10   |-10   |-10   |-10   |-10   |-10   |-10   |


Where *I* and *D* represent the cost of inserting or deleting a residue from
the scanned sequence.

The sequence can then be scanned against a model by a simple dynamic
programming algorithm as in Needlemen-Wunsch or Smith-Waterman. As I
understand it, this is more or less an implementation of a Hidden Markov Model
(HMM), though we do not specify that the scores need be related to emission
likelihoods.

Given the similarity to an HMM, it would probably be better to simply use
existing HMM software (eg. hmmer), or at least to properly implement an HMM
scanning algortithm. The reasons for me to write these functions are
primarily:

1. As a logical exercise
2. To provide a framework with R where the statistical consequences of
   parameter choice can be explored


## Compilation

To compile the code, enter the `src` directory and call:

```sh
R CMD SHLIB model_aligner_r.c model_aligner.c
```

This will create a file called `model_aligner_r.so`. This file needs to be
included in your R session using `dyn.load()`.

## Usage

The shared object `model_aligner_r.so` defines a function `align_model` that
can be called through the `.Call()` interface:

```R
dyn.load( "<path to model_aligner>/src/model_aligner_r.so" )

## example model:
model <- cbind(c(4,-4,-4,-4,-10, -10), c(-4, 4, -4, -4, -10, -10),
               c(2,-4,-4,-2, 0,  -10), c(-4, -4, 4, -4, -10, -10),
               c(4,-4,-4,-4,-10, -10), c(-4, -4, 4, -4, -10, -10),
               c(4, -4, -4, -4, -10, -10))

## make sure the model is represented by integer values
model <- matrix( as.integer(model), nrow=nrow(model), ncol=ncol(model))
seq <- "ACATTAGATAGGACTATA"

## use the .Call() interface to call the code
alignment <- .Call("align_model", seq, model)
```

`align_model` returns a list of of lists with one element for each sequence
defined by the sequences passed to the function. Each element is a list
containing the following components:

1. The maximum score for the last column of the score table; this is
   equivalent to the global alignment score for the sequence.
2. The scoring table used to trace the alignment.
3. A character vector of length two that gives the alignment. The locations
   where the model was aligned are given by `M` in the second of the two
   character elements.
4. A NULL value; may be removed in future versions.

Clearly we would not normally require the full score matrix, but it is
returned here to enable debugging and to illustrate the alignment
process. Future versions will make the return of this matrix optional.


