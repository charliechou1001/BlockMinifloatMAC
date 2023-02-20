#ifndef TYPEDEF_
#define TYPEDEF_

//W = 1+EL+ML = 1+ER+MR
#define W 8
//bit-length of data, 8 or 4 bit
#define WI 5
//extra bit for accumulator to avoid overflow

//left input
#define EL 2
#define ML 5

//right input
#define ER 4
#define MR 3

//output mode
#define ResE 3
#define ResM 4

#define SinM 16
#define SinN 16
#define SinK 16

const unsigned Kadd = (1+((1<<EL)-1+ML)+((1<<ER)-1+MR)+WI);

typedef ap_int<8> SEXP_T;


#endif
