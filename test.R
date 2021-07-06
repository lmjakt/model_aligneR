dyn.load( "src/model_aligner_r.so" )

## try to call "align_model" with no arguments. Should give an error:
.Call("align_model");

## good. Lets call with the wrong type of arguments:
.Call("align_model", 3, 4)

## good.
.Call("align_model", "ACATTAGATAGGACTATA", cbind( c(4,-4,-4,-4,-10,-10), c(-4, 4, -4, -4, -10, -10),
                                                 c(4,-4,-4,-4,-10, 0), c(-4, -4, 4, -4, -10, -10),
                                                 c(4,-4,-4,-4,-10, 0), c(-4, -4, 4, -4, -10, -10),
                                                 c(4, -4, -4, -4, -10, -10) ) )
## Note that the model consists of 6 rows where the rows indicate the following:
## 1-4:  scores for A, C, T|U, G
## 5:    The cost of an insertion into the sequence (i.e. opening a gap between motifs)
## 6:    The cost of a deletion from the sequence (i.e. a nucleotide missing in the sequence
##       but present in the model)
align.model <- function(seq, model){
    model=matrix(as.integer(model), nrow=nrow(model), ncol=ncol(model))
    .Call("align_model", seq, model)
}

## that does something without warnings, but seems to return something not so reasonable.
model <- cbind( c(4,-4,-4,-4,-10,-10), c(-4, 4, -4, -4, -10, -10),
               c(4,-4,-4,-4,0, -10), c(-4, -4, 4, -4, -10, -10),
               c(4,-4,-4,-4,-10, -10), c(-4, -4, 4, -4, -10, -10),
               c(4, -4, -4, -4, -10, -10))

model <- matrix( as.integer(model), nrow=nrow(model), ncol=ncol(model))
image(t(model[1:4,]))
seq <- "ACATTAGATAGGACTATA"

align.model(seq, model)

