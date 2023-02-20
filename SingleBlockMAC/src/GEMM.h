// this is the GEMM for arithmetic test
#include "MAC.h"

#define Ebit 4


void SingleTileGEMM(ap_uint<W> A[SinM][SinK], ap_uint<W> B[SinK][SinN],SEXP_T betaA, SEXP_T betaB, SEXP_T &betaC, ap_uint<W> C[SinM][SinN]){
	ap_int<Kadd> accum[SinM][SinN];
	ap_int<Kadd-WI> mult;
	for(int i=0;i<SinM;i++){
		for(int j=0;j<SinN;j++){
			for(int k=0;k<SinK;k++){
				if(k==0)
					accum[i][j] = 0;
				BMMul(A[i][k], B[k][j], mult);
				accum[i][j] += mult;
			}
		}
	}

	ap_uint<W> Z_min=(1<<W)-1;
	for(int i=0;i<SinM;i++){
		for(int j=0;j<SinN;j++){
			CheckZmin<W>(accum[i][j],Z_min);
		}
	}

	betaC = betaA + betaB;
	ap_int<W> Zshift=0;
	BiasAdjust<W,W>(betaC, Z_min, Zshift);
	for(int i=0;i<SinM;i++){
		for(int j=0;j<SinN;j++){
			Normalization<W,W,4>(accum[i][j],Zshift,C[i][j]);
		}
	}
}

