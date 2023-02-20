#include "ap_int.h"
#include "typedef.h"
#include "newCLZ.h"

void BMMul(ap_uint<W> a, ap_uint<W> b, ap_int<Kadd-WI> &mult){

	const unsigned EA = (1<<EL)-1;
	const unsigned EB = (1<<ER)-1;

	ap_uint<1> sgnA = a[W-1];
	ap_uint<1> sgnB = b[W-1];

	ap_uint<EL> expA = a(W-2,W-2-EL+1);
	ap_uint<1> denormA = (expA == 0)? ap_uint<1>(1) : ap_uint<1>(0);

	ap_uint<ER> expB = b(W-2,W-2-ER+1);
	ap_uint<1> denormB = (expB == 0)? ap_uint<1>(1) : ap_uint<1>(0);

	ap_uint<ML+1> fracA = ( ~denormA, a(ML-1,0));
	ap_int<1+EA+ML> A = ((1-2*sgnA)*ap_int<1+EA+ML>(fracA)) << (expA -1 + denormA);//(expA -1 + denormA) is the number for shift, -1 for e=0 & e=1 case
	ap_uint<MR+1> fracB = ( ~denormB, b(MR-1,0));
	ap_int<1+EB+MR> B = ((1-2*sgnB)*ap_int<1+EB+MR>(fracB)) << (expB -1 + denormB);

	mult = ap_int<Kadd-WI>(A*B);
}


template<int Zm=W, int Zsh=W>
void BiasAdjust(SEXP_T &beta, ap_uint<Zm> Z_min, ap_int<Zsh> &Zshift){

	ap_uint<ML> manA = ML;
	ap_uint<MR> manB = MR;

	ap_uint<EL> biasA = ap_uint<EL>((1<<(EL-1))-1);
	ap_uint<ER> biasB = ap_uint<ER>((1<<(ER-1))-1);
	ap_uint<ResE> ResBias = ap_uint<ResE>((1<<(ResE-1))-1);

	ap_uint<8> CntOverflow = (Kadd-1) - (biasA-1+manA + biasB-1+manB) - ( (1<<ResE)-1 -ResBias + 1);
	ap_int<Zm> emax_shift = CntOverflow - Z_min;

	Zshift = emax_shift;
	beta =  beta + Zshift;
}

template<int Zm=W>
void CheckZmin(ap_int<Kadd> mac_out, ap_uint<Zm> &Z_min){
	ap_uint<1> signR = mac_out[Kadd-1];
	ap_uint<Kadd-1> sgn_rst = signR? (~ap_uint<Kadd-1>(mac_out(Kadd-2,0))+1) : mac_out(Kadd-2,0);

	ap_uint<Zm> count = CLZ64((sgn_rst,ap_uint<64-Kadd+1>(1)));

	if(Z_min > count){
		Z_min = count;
	}
}


template<int Zsh=W, int Zm=W, int mTR=4>
void Normalization(ap_int<Kadd> mac_out, ap_int<Zsh> Zshift, ap_uint<W> &R){

	ap_uint<ML> manA = ML;
	ap_uint<MR> manB = MR;

	ap_uint<EL> biasA = ap_uint<EL>((1<<(EL-1))-1);
	ap_uint<ER> biasB = ap_uint<ER>((1<<(ER-1))-1);
	ap_uint<ResE> ResBias = ap_uint<ResE>((1<<(ResE-1))-1);

	ap_uint<8> expZeroPoint = (Kadd -1) - (biasA -1 + manA + biasB -1 + manB);
	ap_uint<8> CntDenorm = (Kadd-1) - (biasA-1+manA + biasB-1+manB) +(ResBias-1);

	ap_uint<1> signR = mac_out[Kadd-1];
	ap_uint<Kadd-1> sgn_rst = signR? (~ap_uint<Kadd-1>(mac_out(Kadd-2,0))+1) : mac_out(Kadd-2,0);
	sgn_rst = sgn_rst >> Zshift;

	ap_uint<Zm> count = CLZ64((sgn_rst,ap_uint<64-Kadd+1>(1)));

	ap_uint<ResE> exp;
	if(count >= CntDenorm)//denorm and underflow
		exp = 0;
	else
		exp = expZeroPoint - (count + 1) + ResBias;//1 for implicit bit

	ap_uint<Zm> count_adjust;
	ap_uint<ResM> man;
	ap_uint<1> rand;
	if(count < CntDenorm)
		count_adjust = count;
	else
		count_adjust = ap_uint<8>(CntDenorm-1);
	sgn_rst = sgn_rst << count_adjust;

	man = sgn_rst(Kadd-3,Kadd-2-ResM);
	ap_uint<1> guard = sgn_rst[Kadd-2-ResM];
	ap_uint<1> round = sgn_rst[Kadd-3-ResM];
	ap_uint<1> sticky = sgn_rst(Kadd-4-ResM,0)==0? ap_uint<1>(0) : ap_uint<1>(1);
	rand = round & sticky | guard&round&(~sticky);
	if(man < (1<<ResM)-1)
		man = man + rand;
	else if(rand==1 && exp < (1<<ResE)-1){
		exp ++;
		man = 0;
	}

	R[W-1] = signR;
	R(W-2,W-1-ResE) = exp;
	R(ResM-1,0) = man;

}
