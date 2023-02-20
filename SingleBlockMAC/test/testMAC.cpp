#include "ap_int.h"
#include "pydata.h"
#include "GEMM.h"
#include "singleBMToFP.h"

void singleGEMM_test(){
	SEXP_T betaC;
	ap_uint<W> C[SinM][SinN];

	SingleTileGEMM(SA00,SB00,betaSA00,betaSB00,betaC,C);
	SingleMatrixCompare(C, SC_0_0, betaC, betaSC_0_0);

}


int main(){
	singleGEMM_test();
}
