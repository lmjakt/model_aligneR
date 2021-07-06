#include <string.h>
#include "model_aligner.h"

int which_max(int *scores, int l){
  int max = scores[0];
  int max_i = 0;
  for(int i=1; i < l; ++i){
    if(max < scores[i]){
      max = scores[i];
      max_i = i;
    }
  }
  return(max_i);
}

void score_model(  const unsigned char *seq, int seq_length,
		   int *model, int model_length, int model_height,
		   int *score_table, unsigned char *ptr_table )
{
  // Set the score of the first column to all 0s since we do not
  // care where the model starts to align to the sequence
  if(model_height != 6)  // should provide a sutiable error message here
    return;
  //  int t_width = 1 + model_length;
  int t_height = 1 + seq_length;
  memset( score_table, 0, sizeof(int) * t_height );
  memset( ptr_table, 0, sizeof(unsigned char) * t_height);
  // Set the initial row of scores to the cumulative sum of the the
  // deletion gap scores (indicating a loss of the model components
  // from the sequence.
  for(int i=0; i < model_length; ++i){
    int t_offset = (i + 1) * t_height;
    int model_offset = i * model_height;
    score_table[ t_offset ] = score_table[ t_offset - t_height ] + model[model_offset + 5];
    ptr_table[ t_offset ] = 1;
    for(int j=0; j < seq_length; ++j){
      ++t_offset;
      ptr_table[ t_offset ] = 0;
      int scores[3] = {0, 0, 0}; // insertion, deletion, match
      scores[0] = score_table[ t_offset - 1 ] + model[ model_offset + 4 ];
      scores[1] = score_table[ t_offset - t_height ] + model[ model_offset + 5];
      scores[2] = score_table[ t_offset - (t_height + 1) ] + model[ model_offset + ((seq[j] & 6) >> 1) ];
      int max_k = which_max(scores, 3);
      score_table[ t_offset ] = scores[ max_k ];
      ptr_table[ t_offset ] = max_k;
    }
  }
  // and that is basically it.
}

// memory for gapped_seq and gapped_model must have been allocated
// properlly.
// max_score is a pointer to an integer whose value will be set by the function.
int trace_path(const unsigned char *seq, int seq_length,
		unsigned char *ptr_table, int *score_table,
		int t_width, int t_height,
	       char *gapped_seq, char *gapped_model, int gapped_l,
	       int *max_score){
  int gapped_o = gapped_l;
  gapped_seq[gapped_o] = 0;
  gapped_model[gapped_o] = 0;
  // First find the maximum score in the last column of the table;
  int t_offset = (t_width -1) * t_height;
  *max_score = score_table[t_offset];
  int max_offset = 0;
  for(int t_offset = (t_width -1) * t_height; t_offset < (t_width * t_height); ++t_offset){
    if( score_table[ t_offset ] > *max_score ){
      max_offset = t_offset;
      *max_score = score_table[ t_offset ];
    }
  }
  // then follow the arrows backwards:
  t_offset = max_offset;
  int row = t_offset % t_height;
  int column = t_offset / t_height;
  while(row && column){
    --gapped_o;
    switch( ptr_table[ t_offset ] ){
    case 0:
      gapped_seq[gapped_o] = seq[row-1];
      gapped_model[gapped_o] = '-';
      --t_offset;
      break;
    case 1:
      gapped_seq[gapped_o] = '-';
      gapped_model[gapped_o] = 'M';
      t_offset -= t_height;
      break;
    case 2:
      gapped_seq[gapped_o] = seq[row-1];
      gapped_model[gapped_o] = 'M';
      t_offset -= (t_height + 1);
      break;
      //    default:
      // do nothing.
    }
    row = t_offset % t_height;
    column = t_offset / t_height;
  }
  return(gapped_o);
}
