#include <cmath>
#include "ap_int.h"
#include "typedef.h"

float BMToFP(ap_uint<W> a, SEXP_T shared_exp){
	float two_exp;
	float frac;
	float sexp = std::pow(2,int(shared_exp));
	float res;

	ap_uint<ResM> man = a(ResM-1,0);
	ap_int<2> sgn = a[W-1]?-1:1;

	ap_uint<1> denorm = a(W-2,W-ResE-1)==0? ap_uint<1>(1):ap_uint<1>(0);
	ap_int<ResE+1> exp = a(W-2,W-ResE-1) - (1<<(ResE-1))+1 + denorm;
	two_exp = std::pow(2,exp.to_int());
	frac = man.to_float()*std::pow(2,-ResM) + (~denorm);
	res = sgn*frac*two_exp*sexp;

	return res;
}

void SingleMatrixCompare(ap_uint<W> A[SinM][SinN], ap_uint<W> B[SinM][SinN], SEXP_T  betaA, SEXP_T betaB){
	float sum =0;
	for(int i=0;i<SinM;i++){
		for(int j=0;j<SinN;j++){
			float a = BMToFP(A[i][j],betaA);
			float b = BMToFP(B[i][j],betaB);
			if(a!=b){
				printf("i:%d, j:%d, A: %s, B: %s, a:%f, b:%f\n",i,j,A[i][j].to_string(2).c_str(), B[i][j].to_string(2).c_str(), a,b);
				sum += (a-b)*(a-b);
			}
		}
	}
	sum = sum/(SinM*SinN);
	printf("average MSE is: %f\n",sum);
}

