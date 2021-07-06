#ifndef MODEL_ALIGNER_H
#define MODEL_ALIGNER_H

// seq: should contain a sequence of only ACTGactg
// Nucleotides are mapped to offsets by:
// (c & 6) >> 1
// This maps, A => 0, C => 1, T => 2, U => 2, G => 3
//      other characters will cause nonsensical results.
// seq_length; the length of the sequence
// model: a matrix of width seq_l with 6 rows. The 6 rows give
//        the position weight matrix for:
//            A, C, T, G, insertion, deletion
// model_height: the number of rows in model. Should be 6
//        with respect to the sequence (i.e. insertion means
//        characters present in the sequence but not the model).
// score_table: a preallocated score table with (seq_length + 1) rows
//              and (1 + model_length) columns.
// ptr_table: a preallocated table containing the reverse pointers giving
//            the prior cell from which the score was calculated.
//            0 => vertical, 1 => horizontal, 2 => diagonal (i.e. match)
//
// Note that we use a column major order as in R for all tables.
void score_model( const unsigned char *seq, int seq_length,
		  int *model, int model_length, int model_height,
		  int *score_table, unsigned char *ptr_table );

// gapped_seq and gapped model should both be of length:
// t_width + t_height + 1 and be 0 terminated.
// The function will fill in the gaps of the traced
// sequences; for the sequence it will give the sequence characters
// or gap characters (for deletions)
// for the gapped model, it will return simply, M or - indicating
// either an insertion (-) or a match to the model.
// The length of the gapped sequencs is returned.
int trace_path(const unsigned char *seq, int seq_length,
	       unsigned char *ptr_table, int *score_table,
	       int t_width, int t_height,
	       char *gapped_seq, char *gapped_model, int gapped_l,
	       int *max_score);
  

#endif
