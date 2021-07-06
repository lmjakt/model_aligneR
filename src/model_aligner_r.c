#include <R.h>
#include <Rinternals.h>
#include "model_aligner.h"

// R wrapper functions for model_aligner.

// seq_r    : a character vector of length > 0
// model_r  : a matrix of n columns and 6 rows. The rows indicate
//            scores associated with: A, C, T/U, G and
//            cost of inserting a vertical gap (position present in sequence, not model)
//                                horizontal gap (position present in model, missing from sequence)
// At the moment there is no sane way of defining the cost of the distance between two adjacent motifs
// defined by the model, but for short sequences it should be OK to allow free vertical gaps at the end
// of motifs. Ideally we would want it to be possible to define an R function that takes the distance
// and calculates a context specific gap.
SEXP align_model_r(SEXP seq_r, SEXP model_r){
  if(TYPEOF(seq_r) != STRSXP || length(seq_r) < 1)
    error("seq should be a character vector of length > 0");
  if(TYPEOF(model_r) != INTSXP)
    error("the model should be defined as a matrix of integers");
  SEXP model_dims_r = PROTECT( getAttrib( model_r, R_DimSymbol ));
  if(length(model_dims_r) != 2)
    error("The model should be a matrix (i.e. have two dimensions)");
  int *model_dims = INTEGER(model_dims_r);
  int model_row_n = model_dims[0];
  int model_column_n = model_dims[1];
  if(model_row_n != 6)
    error("The model matrix should have 6 rows.");
  if(model_column_n < 1)
    error("The model has a length of 0.");
  int *model = INTEGER(model_r);
  // Then we can go ahead and make the required data structures.
  // We want to return a list of alignment data, with one element for each
  // sequence.
  // Each element of this list is also a list containning:
  // 1. The score of the best alignment
  // 2. The score matrix
  // 3. A character vector containing the gapped alignments
  // [4. And maybe one day a two dimensional matrix containing the row and column
  //    positions of the traceback. ]
  SEXP ret_list = PROTECT(allocVector(VECSXP, length(seq_r)));
  for(int i=0; i < length(seq_r); ++i){
    SET_VECTOR_ELT( ret_list, i, allocVector(VECSXP, 4) );
    SEXP ret_i = VECTOR_ELT( ret_list, i );
    SET_VECTOR_ELT( ret_i, 0, allocVector(INTSXP, 1));
    SEXP seq_i = STRING_ELT(seq_r, i);
    int seq_len = length(seq_i);
    const unsigned char *seq = (const unsigned char*)CHAR(seq_i);
    SET_VECTOR_ELT( ret_i, 1, allocMatrix( INTSXP, (1 + seq_len), (1 + model_column_n) ) );
    int *score_table = INTEGER( VECTOR_ELT(ret_i, 1) );
    // we don't return the pointer table, so use malloc to allocate and then free
    unsigned char *ptr_table = malloc( (1 + seq_len) * (1 + model_column_n) * sizeof(unsigned char) );
    score_model( seq, seq_len, model, model_column_n, model_row_n, score_table, ptr_table );
    int gapped_l = seq_len + model_column_n;
    char *gapped_seq = malloc( 1 + gapped_l * sizeof( char ));
    char *gapped_model = malloc( 1 + gapped_l * sizeof( char ));
    int gapped_offset = trace_path( seq, seq_len, ptr_table, score_table, 1 + model_column_n, 1 + seq_len,
				    gapped_seq, gapped_model, gapped_l,
				    INTEGER(VECTOR_ELT(ret_i, 0)));
    SET_VECTOR_ELT( ret_i, 2, allocVector(STRSXP, 2));
    SEXP gapped_r = VECTOR_ELT(ret_i, 2 );
    SET_STRING_ELT(gapped_r, 0, mkChar( gapped_seq + gapped_offset ));
    SET_STRING_ELT(gapped_r, 1, mkChar( gapped_model + gapped_offset ));
    free( gapped_seq );
    free( gapped_model );
    free( ptr_table );
  }
  UNPROTECT(2);
  return( ret_list );
}


static const R_CallMethodDef callMethods[] = {
  {"align_model", (DL_FUNC)&align_model_r, 2},
  {NULL, NULL, 0}
};

void R_init_model_aligner_r(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
