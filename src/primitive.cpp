/*
 * primitive.cpp
 * 
 * Copyright (C) 2013 Louis-Francois Handfield
 * e-mail: lfhandfield@gmail.com
 *
 * This program is free software; upon notification by email to the licensor
 * of the licencee identity and nature of use, the licencee can redistribute
 * this program and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 
 * of the License, or (at the licencee option) any later version. As such,
 * no further notifications are required.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "primitive.h"

namespace LFHPrimitive{
	LFH_USE_MEMLEAK_HELP(LFHCONCAT2(myHashmap<unsigned int, unsigned int> memleak_helper;))
	
	// 8x12 font for 96 characters
	unsigned int def_font[] = {
	0x00000000,0x00000000,0x00000000,// space
	0x00180000,0x18181818,0x00181818,// !
	0x00000000,0x00000000,0x003C3C3C,
	0x12120000,0x2424247F,0x004848FE,
	0x3C180000,0x063C6066,0x00183C66,
	0xF3600000,0x30180C66,0x0006CF66,
	0x33FE0000,0x0C0EDB7B,0x001C0606,
	0x00000000,0x00000000,0x00181818,
	0x0C181830,0x0C0C0C0C,0x00301818,
	0x3018180C,0x30303030,0x000C1818,
	0x00000000,0x3CFF3C66,0x00000066,
	0x18180000,0x1818FF18,0x00000018,
	0x00180C00,0x00000000,0x00000000,
	0x00000000,0x0000FF00,0x00000000,
	0x00180000,0x00000000,0x00000000,
	0x03000000,0x30180C06,0x0000C060,
	
	0x663C0000,0x66666666,0x003C6666,
	0x187E0000,0x18181818,0x00181E18,
	0x067E0000,0x6030180C,0x003C6660,
	0x663C0000,0x60386060,0x003C6660,
	0x30780000,0x3C3C367E,0x00303838,
	0x663C0000,0x063E6060,0x007E0606,
	0x663C0000,0x063E6666,0x00380C06,
	0x0C0C0000,0x30301818,0x007E6660,
	0x663C0000,0x663C6666,0x003C6666,
	0x301C0000,0x667C6060,0x003C6666,
	0x00180000,0x18000000,0x00000000,
	0x00180C00,0x18000000,0x00000000,
	0x18300000,0x180C060C,0x00000030,
	0x00000000,0x7E007E00,0x00000000,
	0x180C0000,0x18306030,0x0000000C,
	0x00180000,0x60301818,0x003C6660,
	
	0x063C0000,0xA5A5BDDB,0x007CC6BB,
	0x66E70000,0x3C24667E,0x00101818,
	0xC67F0000,0xC67EC6C6,0x007FC6C6,
	0xC67C0000,0x03030303,0x007CC603,
	0x663F0000,0xC6C6C6C6,0x003F66C6,
	0xC6FF0000,0x363E3606,0x00FFC606,
	0x061F0000,0x363E3606,0x00FFC606,
	0xC67C0000,0x0303F3C3,0x007CC603,
	0x66E70000,0x667E6666,0x00E76666,
	0x187E0000,0x18181818,0x007E1818,
	0x331E0000,0x30303033,0x007C3030,
	0x66EF0000,0x361E1E36,0x00EF6636,
	0xCCFF0000,0x0C0C0C0C,0x003F0C0C,
	0x66E70000,0x7E5A5A42,0x00E7667E,
	0x666F0000,0x6E7E7676,0x00F7666E,
	0x663C0000,0xC3C3C3C3,0x003C66C3,
	
	0x061F0000,0xC67E0606,0x007FC6C6,
	0x663CFC00,0xC3C3C3C3,0x003C66C3,
	0x66EF0000,0xC67E3636,0x007FC6C6,
	0xC37E0000,0x037EC0C0,0x007EC303,
	0x183C0000,0x18181818,0x00FFDB18,
	0x663C0000,0x66666666,0x00FF6666,
	0x18180000,0x663C3C3C,0x00E76666,
	0x24240000,0x425A7E7E,0x00E76666,
	0x66E70000,0x3C183C3C,0x00E7663C,
	0x183C0000,0x3C3C1818,0x00E76666,
	0xC6FF0000,0x30180C0C,0x00FF6330,
	0x0C0C0C3C,0x0C0C0C0C,0x003C0C0C,
	0xC0000000,0x0C183060,0x00000306,
	0x3030303C,0x30303030,0x003C3030,
	0x00000000,0x00000000,0x183C6600,
	0x000000FF,0x00000000,0x00000000,
	
	0x00000000,0x00000000,0x000C1800,
	0x63FE0000,0x3E607E63,0x00000000, // a
	0xC67F0000,0x7EC6C6C6,0x00070606,
	0xC37E0000,0x7EC30303,0x00000000,
	0x63FE0000,0x7E636363,0x00706060,
	0xC37E0000,0x7EC3FF03,0x00000000, // e
	0x0C3E0000,0x3E0C0C0C,0x00380C0C,
	0x7E60603E,0xFE636363,0x00000000,
	0x66EF0000,0x3E6E6666,0x00070606,
	0x187E0000,0x1E181818,0x00180000,
	0x3030301E,0x3E303030,0x00300000,
	0x36E70000,0x76361E1E,0x00070606,
	0x187E0000,0x18181818,0x001C1818,
	0x4ADB0000,0x3F6A6A6A,0x00000000,
	0x66E70000,0x3F6E6666,0x00000000,
	0xC37E0000,0x7EC3C3C3,0x00000000,
	
	0xC67E060F,0x7FC6C6C6,0x00000000,
	0x637E60F0,0xFE636363,0x00000000,
	0x0C3F0000,0x7FDC0C0C,0x00000000,
	0xC37E0000,0x7EC31E70,0x00000000,
	0x6C380000,0x3E0C0C0C,0x00000C0C,
	0x76FC0000,0x77666666,0x00000000,
	0x3C180000,0xE766663C,0x00000000,
	0x3C240000,0xC3665A7E,0x00000000,
	0x66C30000,0xC3663C3C,0x00000000,
	0x3C18180E,0xC366663C,0x00000000,
	0x667F0000,0x7F33180C,0x00000000,
	0x18181870,0x18180E18,0x00701818,
	0x18181818,0x18181818,0x18181818,
	0x1818180E,0x18187018,0x000E1818,
	0x00000000,0x73000000,0x0000CEDB,
	0x00000000,0x00000000,0x00000000,
	
	0xC37E0000,0x7EC3FF03,0x60FC6000,
	0xC37E0000,0x7EC3FF03,0x00FC0000};
	
	
	myHashmap< unsigned int , pair< void*, unsigned int> > central_alias_bank;
	myHashmap< void* , unsigned int > central_alias_bank_back;
	
	char* cloneString(const char* const what){
		unsigned int i = strlen(what) +1;
		char* _out = new char[i];
		memcpy(_out,what,sizeof(char)*i);
		return(_out);
	}
	
LFH_GOLD double lngamma(double xx){ double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091, -1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
    int j; y=x=xx; tmp=x+5.5f;tmp-=(x+0.5)*log(tmp);ser= 1.000000000190015;
    for(j=0;j<=5;j++)ser += cof[j]/++y;
    return(-tmp+log(2.5066282746310005*ser/x));
}

LFH_GOLD double logdensity_chisquarre(double x, double free){ double fout = log(0.25f) - x + (free - 2.0f) *log(x * 0.5f); fout *= 0.5f; fout -= lngamma(free * 0.5); return fout;}

#ifdef GNU_SCIENTIFIC_LIBRARY
#include "GSLfunc.hpp"	
	
LFH_GOLD double Pvalue_chisquarre_Ptail(double x, double free){return incgamma_frac(free* 0.5f, x * 0.5f);}
LFH_GOLD double Pvalue_chisquarre_Ntail(double x, double free){return 1.0f - incgamma_frac(free* 0.5f, x* 0.5f);}
LFH_GOLD double Pvalue_Gamma_Ptail(double x, double k, double theta){return incgamma_frac(k, x / theta);}
LFH_GOLD double Pvalue_Gamma_Ntail(double x, double k, double theta){return 1.0f - incgamma_frac(k, x / theta);}
LFH_GOLD double Pvalue_GammaRate_Ptail(double x, double alpha, double beta){return incgamma_frac(alpha, x* beta);}
LFH_GOLD double Pvalue_GammaRate_Ntail(double x, double alpha, double beta){return 1.0f - incgamma_frac(alpha, x* beta);}
LFH_GOLD double Pvalue_SumLogPvalues_Ptail(double sum, int nbpvals){ return incgamma_frac((double) nbpvals , -sum); }
	
LFH_GOLD double Pvalue_Beta_Ntail(double x, double a, double b){gsl_sf_result_struct res;gsl_sf_beta_inc_e(b,a,1.0-x,&res);	return res.val;}
LFH_GOLD double Pvalue_Beta_Ptail(double x, double a, double b){gsl_sf_result_struct res;gsl_sf_beta_inc_e(a,b,x,&res);	return res.val;}
LFH_GOLD double Pvalue_Fdistrib_Ntail(double x, double a, double b){gsl_sf_result_struct res;gsl_sf_beta_inc_e(0.5f*a,0.5f*b,a*x/(a*x+b),&res); return res.val;}
LFH_GOLD double Pvalue_Fdistrib_Ptail(double x, double a, double b){gsl_sf_result_struct res;gsl_sf_beta_inc_e(0.5f*b,0.5f*a,b/(a*x+b),&res);	return res.val;}

// log are all natural logarithms
	
LFH_GOLD double LogPvalue_chisquarre_Ptail(double x, double free){return incgamma_log_frac(free* 0.5f, x * 0.5f);}
LFH_GOLD double LogPvalue_chisquarre_Ntail(double x, double free){return incgamma_log_1minus_frac(free* 0.5f, x* 0.5f);}
LFH_GOLD double LogPvalue_Gamma_Ptail(double x, double k, double theta){return incgamma_log_frac(k, x / theta);}
LFH_GOLD double LogPvalue_Gamma_Ntail(double x, double k, double theta){return incgamma_log_1minus_frac(k, x / theta);}
LFH_GOLD double LogPvalue_GammaRate_Ptail(double x, double alpha, double beta){return incgamma_log_frac(alpha, x* beta);}
LFH_GOLD double LogPvalue_GammaRate_Ntail(double x, double alpha, double beta){return incgamma_log_1minus_frac(alpha, x* beta);}
LFH_GOLD double LogPvalue_SumLogPvalues_Ptail(double sum, int nbpvals){ return incgamma_log_frac((double) nbpvals ,-sum); }
	
LFH_GOLD double LogPvalue_Beta_Ntail(double x, double a, double b){gsl_sf_result_struct res;gsl_sf_log_beta_inc_e(b,a,1.0-x,&res);	return res.val;}
LFH_GOLD double LogPvalue_Beta_Ptail(double x, double a, double b){gsl_sf_result_struct res;gsl_sf_log_beta_inc_e(a,b,x,&res);	return res.val;}
LFH_GOLD double LogPvalue_Fdistrib_Ntail(double x, double a, double b){gsl_sf_result_struct res;gsl_sf_log_beta_inc_e(0.5f*a,0.5f*b,a*x/(a*x+b),&res); return res.val;}
LFH_GOLD double LogPvalue_Fdistrib_Ptail(double x, double a, double b){gsl_sf_result_struct res;gsl_sf_log_beta_inc_e(0.5f*b,0.5f*a,b/(a*x+b),&res); return res.val;}
	
double Gcummul_HotelingTTest(double stat, unsigned int nb_dims, unsigned int nbsample_A, unsigned int nbsample_b){
	// stat = [(NaNb) * (Na+Nb-D-1) / (D*(Na+Nb))] *  ( ((sum_allA X)/Na) - ((sum_allB X)/Nb))^T ( sum_allA&B XX^T )^-1 (((sum_allA X)/Na) - ((sum_allB X)/Nb))	
}
#endif
	
	
void GPParamSearch(const Vector<KeyElem<double, double> >& vec, double &mean, double &sigma_s, double &sigma_n, double &scale){
	CurvatureSearchScope cs; cs.init(0.0000001f, 3);
	double param[3]; param[0] = 0; param[1] =0; param[2] = 0;
	double deriv[3]; double tmp, like;
	Trianglix<double> covar; covar.setSize(vec.getSize());
	unsigned int i,j,k;
	
    // after normalization, noise <= 1 noise = 1 / 1+ exp(param)
	// MarjLikeli = 
	// DMarjLikeli/Da = tr(( K^-1yyTK^-1 - K^-1)  * DK/Da)
	Tuple<double> yyy; yyy.setSize(vec.getSize());
	tmp =0.0f;like =0.0f;scale =0.0f;sigma_n=0.0f;
	for(i=0;i<vec.getSize();i++) {yyy[i] = vec[i].d; tmp += vec[i].d; like += vec[i].d * vec[i].d; sigma_n +=vec[i].k;  scale += vec[i].k*vec[i].k;}
	tmp /= vec.getSize(); sigma_n /= vec.getSize();
	like /= vec.getSize(); scale /= vec.getSize();
	like -= tmp*tmp; sigma_n = like;
	like = pow(like,-0.5f);
	yyy -= tmp; yyy *= like; // centering!
    mean = tmp;
    param[2] = -3.0f-log( 2.0f*(scale - sigma_n));
	do{
        scale = 1.0625f / (1.0f + exp(param[1]));
		for(k=0,i=0;i<vec.getSize();i++){
			for(j=0;j<i;j++) {tmp = vec[i].k - vec[j].k; covar.data[k++] = exp(param[0] - tmp*tmp * exp(param[2]));}
			covar.data[k++] = exp(param[0]) + scale;
		}
		like = (covar.Xformed_inner_product_of_inverse(yyy) + covar.log_determinant() + log(2.0f* M_PI)) * -1.0f;
		Trianglix<double> dainverse = covar.mkinverse();
		Trianglix<double> factrig = dainverse.Xformed_outer_product(yyy) - dainverse;
		covar.cell(0,0) -= scale;
		deriv[1] = factrig.cell(0,0);
        
		for(i=1;i<vec.getSize();i++) {
			deriv[1] += factrig.cell(i,i);
			covar.cell(i,i) -= scale;
		} deriv[1] *= -1.0625f * exp(param[1]) / (1.0f + 2.0f * exp(param[1]) +  exp(2.0f* param[1]));
		deriv[0] = factrig.trace_of_product(covar) * exp(param[0]);
		for(k=0,i=0;i<vec.getSize();i++){
			for(j=0;j<i;j++) {tmp = vec[i].k - vec[j].k; covar.data[k++] *= -tmp*tmp;}
			covar.data[k++] = 0.0f;
		}
		deriv[2] = factrig.trace_of_product(covar) * exp(param[2]);
        printf("f[%e,%e,%e] = %e   deriv = {%e,%e,%e}\n", param[0], param[1], param[2], like,deriv[0],deriv[1],deriv[2]);
	}while(cs.updateAscent(like,param,deriv) > 0.0000001f);
	sigma_s = exp(param[0]) * sigma_n;
	sigma_n *= 1.0625f / (1.0f + exp(param[1]));
	scale = exp(param[2]);
}
	
// the inputs are SQUARED distances
double distanceToEllipse(double d1, double d2, double w, double f){
	double tmpb;
    //double distb;

	tmpb = d1 - d2;
	w *= 0.5f;
	w *= d1 + d2;
	f *= 0.5f;

	//distb = -2.0f * tmpb * tmpb;

	tmpb = d1 - d2;
	d1 += d2;
	w *= 0.5f;
	d1 = 4.0f * w * (d1+d2);
	d2 = -2.0f * tmpb * tmpb;


	f *= 0.5f;

/*	LFHPrimitive::PolyThing<double> poly;
	poly.setOrder(4);

	poly.coef[0] = w*w*(-8.0f*w*w+ 15*w*f-8.0f*f*f+d1) + (w-f)*(w-f)*d2;
	poly.coef[1] = 2.0f*(w*w*(-16.0f*w*w+ 23*w*f-8.0f*f*f+d1) + (w-f)*d2);
	poly.coef[2] = -48.0f*w*w+ 47*w*f-8.0f*f*f+ d1 + d2;
	poly.coef[3] = 16.0f * (f -2*w);
	poly.coef[4] = -8.0f;

	w =poly.newtonRoot(0.0f, f - w, LFHPrimitive::ExCo<double>::max());
	*/
	return(w);
}


double CubicRealRoot(double* coef, bool wantmiddle){
	double tmp = 3 *coef[3];
	double q = (coef[1] - (coef[2]*coef[2] / tmp) ) / tmp;
	double r = (-1.5f*coef[0] - ( ((coef[2]*coef[2]*coef[2] / tmp) - 1.5f*coef[2]*coef[1] )/tmp) )/tmp;
	tmp = q*q*q + r*r;
	//	printf("r:%f\td:%f\n",r,q);
	double s,t,ac,as;
	if ( tmp >= 0){ // only 1 real
		//		printf("preroot: s=%f\tt=%f\n",r + sqrt(q),r - sqrt(q));
		q = sqrt(tmp);
		tmp = r + q; s = pow(fabs(tmp),1/3.0f) * ((tmp <0) ? -1.0f:1.0f);
		tmp = r - q; t = pow(fabs(tmp),1/3.0f) * ((tmp <0) ? -1.0f:1.0f);
		//		printf("s=%f\tt=%f\n",s,t);
		return(s + t - coef[2] / (3*coef[3]));

		}else{ // all 3 are real, return the smallest one (closest to 0)
			tmp = r*r - tmp;
			t = acos(r/sqrt(tmp))/3;
			q = sqrt(-q); // q >= 0 is guarrantied!


			ac = cos(t);
			as = sin(t);
			if (wantmiddle){

				if (fabs(ac)*3.0f > fabs(as)) return((-ac +as* ((ac < 0.0f) ? sqrt(3.0f): -sqrt(3.0f) ))*q- coef[2] / (3*coef[3]));
				else return(ac*q*2 - coef[2] / (3*coef[3]));


			}else{

			s = ac*q*2 - coef[2] / (3*coef[3]);
			tmp = (-ac +as*sqrt(3.0f) )*q - coef[2] / (3*coef[3]);
			if (fabs(tmp) < fabs(s)) s =tmp;
			tmp = (-ac -as*sqrt(3.0f) )*q - coef[2] / (3*coef[3]);
			if (fabs(tmp) < fabs(s)) s =tmp;
			return(s);
			}
		}

}

// point at -3,-1,1,3 return closest to center, valid root
double CubicInterpolationRoot(double* pt , bool is_monotonic){

	double coef[4];
	double q,r;
	q = (pt[3] - pt[0])/ (3 * 16.0f);
	r = (pt[2] - pt[1])/16.0f;
	coef[1] = -q + 9.0f *r;
	coef[3] = q - r;
	q = (pt[3] + pt[0])/ 16.0f;
	r = (pt[2] + pt[1]) /(16.0f);
	coef[0] = -q + 9.0f *r;
	coef[2] = q - r;



	double tmp = 3 *coef[3];
	q = (coef[1] - (coef[2]*coef[2] / tmp) ) / tmp;
	r = (-1.5f*coef[0] - ( ((coef[2]*coef[2]*coef[2] / tmp) - 1.5f*coef[2]*coef[1] )/tmp) )/tmp;
	tmp = q*q*q + r*r;
	//	printf("r:%f\td:%f\n",r,q);
	double s,t,ac,as;
	if ( tmp >= 0){ // only 1 real
		//		printf("preroot: s=%f\tt=%f\n",r + sqrt(q),r - sqrt(q));
		q = sqrt(tmp);
		tmp = r + q; s = pow(fabs(tmp),1/3.0f) * ((tmp <0) ? -1.0f:1.0f);
		tmp = r - q; t = pow(fabs(tmp),1/3.0f) * ((tmp <0) ? -1.0f:1.0f);
		//		printf("s=%f\tt=%f\n",s,t);
		return(s + t - coef[2] / (3*coef[3]));

	}else{ // all 3 are real, return the smallest one (closest to 0)
		tmp = r*r - tmp;
		t = acos(r/sqrt(tmp))/3;
		q = sqrt(-q); // q >= 0 is guarrantied!


		ac = cos(t);
		as = sin(t);

//		if (is_monotonic){
			// ignores illegal roots!
//		}else{
			s = ac*q*2 - coef[2] / (3*coef[3]);
			tmp = (-ac +as*sqrt(3.0f) )*q - coef[2] / (3*coef[3]);
			if (fabs(tmp) < fabs(s)) s =tmp;
			tmp = (-ac -as*sqrt(3.0f) )*q - coef[2] / (3*coef[3]);
			if (fabs(tmp) < fabs(s)) s =tmp;
		return(s);
//		}



	}
}

double sinc(double xx){

	if (fabs(xx) < pow(0.5f,23.0f)) return(1.0f);
	else return(sin(xx) / xx);

	}

double sinc_d(double xx){
	// (sin x - x cos x) / x^2

		if (fabs(xx) < pow(0.5f,23.0f)) return((-2.0f  / 3.0f) * xx);
	else return(((sin(xx) / xx) - cos(xx))/ xx);

	}

double sampleGaussian(){return(sqrt(-2 * log((1 + rand()) / ((float)RAND_MAX))) * cos(M_PI *2 * rand() / ((float)RAND_MAX)));}
	
SetComparison::SetComparison(){}
SetComparison::SetComparison(int v) : comp(v){}

WarH<LFH_WARNINGS> static_warning_handdle;
	
void ArgumentParser::readList(char* const text, vector<unsigned int> & _out){

	char* cur = text;
	bool range = false;
	unsigned int tmpint;
	unsigned int tmpint_ite;
	while(true){
		switch(*cur){
			case '\0': return;
			case ' ': break;
			case '-': range = true;
			default:
				if (((*cur) >= '0')&&((*cur) <= '9')){
					tmpint = atoi(cur);
					while(((*cur) >= '0')&&((*cur) <= '9')) cur++;
					cur--;
					if (range){
						tmpint_ite =  _out[_out.size()-1];
						if (tmpint >= tmpint_ite) for(;tmpint_ite <=tmpint ;tmpint_ite++) _out.push_back(tmpint_ite);
						else for(;tmpint_ite >=tmpint ;tmpint_ite--) _out.push_back(tmpint_ite);
					}else _out.push_back(tmpint);
			}	}
			cur++;
	}	}
	void ArgumentParser::readList(char* const text, Vector<unsigned int> & _out){
		
		char* cur = text;
		bool range = false;
		unsigned int tmpint;
		unsigned int tmpint_ite;
		while(true){
			switch(*cur){
				case '\0': return;
				case ' ': break;
				case '-': range = true;
				default:
				if (((*cur) >= '0')&&((*cur) <= '9')){
					tmpint = atoi(cur);
					while(((*cur) >= '0')&&((*cur) <= '9')) cur++;
					cur--;
					if (range){
						tmpint_ite =  _out[_out.getSize()-1];
						if (tmpint >= tmpint_ite) for(;tmpint_ite <=tmpint ;tmpint_ite++) _out.push_back(tmpint_ite);
						else for(;tmpint_ite >=tmpint ;tmpint_ite--) _out.push_back(tmpint_ite);
					}else _out.push_back(tmpint);
			}	}
			cur++;
	}	}

int ArgumentParser::operator()(int nbargs, char * const * args){
	int min, max,i,j;
	nbargs--;
	char * const * cargs = args +1;
	int min_default;
	
	char emptyfun[] = "";
	max =0;nbaddtoken(emptyfun, min_default, max);
	
	while((nbargs > min_default)&&(cargs[0][0] == '-')){
		
		if (cargs[0][1] == '\0') {nbargs--; cargs++; break;}
		else if (cargs[0][1] == 'h') {
			if (cargs[0][2] == 'L')	printf("* Copyright (C) 2013 Louis-Francois Handfield\n* e-mail: lfhandfield@gmail.com\n* \n* This program is free software; uppon the notification to the licensor\n* of the licencee identity and nature of use, the licencee can redistribute this\n* program and/or modify it under the terms of the GNU General Public\n* License as published by the Free Software Foundation; either version 2\n* of the License, or (at the licencee option) any later version. As such,\n* no further notification of any kind are required.\n*\n* This program is distributed in the hope that it will be useful, but\n* WITHOUT ANY WARRANTY; without even the implied warranty of\n* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU\n* General Public License for more details.\n*\n* You should have received a copy of the GNU General Public License\n* along with this program; if not, write to the Free Software\n* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.\n");
			else {
				help();
				printf("\n-hL for source code licence for this software\n");
				printf("Compiled on the %s, at %s\n",__DATE__, __TIME__);
			}
			exit(0);
		}
		max =0;nbaddtoken((*cargs) + 1, min, max);
		
		if ((min == max)&&(max != 0)){
			store(cargs,nbargs);
			return defstore(cargs,nbargs);
		}else{
			if (max < min) max =min;
			for(i=1;i<nbargs;i++) {
				if (i >= nbargs - min_default) break;
				if (cargs[i][0] == '-') { //number? ignore
					for(j=1; cargs[i][j] == ' ';j++);
					if (((cargs[i][j] < '0')||(cargs[i][j] > '9'))&&(cargs[i][j] != '.')) break;
				}
			}
			if ( i-1 > max){
				if (i >= nbargs - min_default) i = max+1;
				else{
				if (max > min) fprintf(stderr,"flag %s needs %i to %i arguemnts, found %i\n",(*cargs) + 1 , min, max,i-1);
				else fprintf(stderr,"flag %s needs %i arguments, found %i\n", (*cargs) + 1 , min,i-1);
				fprintf(stderr,"for more details, try %s -h\n", *args);
				exit(1);
				}
			}
			if (i-1 < min) {
				if (max > min) fprintf(stderr,"flag %s needs %i to %i arguemnts, found %i\n",(*cargs) + 1 , min, max,i-1);
				else fprintf(stderr,"flag %s needs %i arguments, found %i\n", (*cargs) + 1 , min,i-1);
				fprintf(stderr,"for more details, try %s -h\n", *args);
				exit(1);
			}
			store(cargs,i-1);
			cargs += i;
			nbargs -= i;
	}	}
	
	for(i=0;i<nbargs;i++){
		if (cargs[i][0] == '-'){
			if (cargs[i][1] == '\0') {nbargs--; cargs++; break;}
			else if (cargs[i][1] == 'h') {
				if (cargs[i][2] == 'L')	printf("* Copyright (C) 2013 Louis-Francois Handfield\n* e-mail: lfhandfield@gmail.com\n* \n* This program is free software; uppon the notification to the licensor\n* of the licencee identity and nature of use, the licencee can redistribute this\n* program and/or modify it under the terms of the GNU General Public\n* License as published by the Free Software Foundation; either version 2\n* of the License, or (at the licencee option) any later version. As such,\n* no further notification of any kind are required.\n*\n* This program is distributed in the hope that it will be useful, but\n* WITHOUT ANY WARRANTY; without even the implied warranty of\n* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU\n* General Public License for more details.\n*\n* You should have received a copy of the GNU General Public License\n* along with this program; if not, write to the Free Software\n* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.\n");
				else {
					help();
					printf("\n-hL for source code licence for this software\n");
					printf("Compiled on the %s, at %s\n",__DATE__, __TIME__);
				}
				return 0;
			}
			max =0;	nbaddtoken((*cargs) + 1, min, max);
			if ((min == max)&&(max != 0)){
				store(cargs,nbargs);
				return defstore(cargs,nbargs);
	}	}	}
	// no flag = default flag \0, (last flag too)
	max =0; nbaddtoken(emptyfun, min, max);
	if (max < min) max =min;
	if ((nbargs < min)||(nbargs > max)){
		if (max > min) fprintf(stderr,"process needs %i to %i arguemnts, found %i\n", min, max,nbargs);
		else fprintf(stderr,"process needs %i arguments, found %i\n", min,nbargs);
		fprintf(stderr,"for more details, try %s -h\n", *args);
		exit(1);
	}
	return defstore(cargs,nbargs);
}
	

	unsigned int state;
	unsigned int lenght;
	void ProgressBarPrint::start(const char* text){
		unsigned int i;
		for(i=0;i<lenght;i++) putc(' ', stdout);
		printf(": %s", text);
		for(i=lenght + strlen(text)+2;i != 0;i--) putc('\b', stdout);fflush(stdout);
		state =0;lasttime = clock();
	}
	void ProgressBarPrint::update(double fraction){
		if ((fraction * lenght) > state) {
		while ((fraction * lenght) > state) {state++;printf("=");}
			if ((clock() - lasttime) & 0xFFFFFC00 ) {fflush(stdout); lasttime =clock();} 
		}
	}
	void ProgressBarPrint::finish(){
		while (lenght > state) {state++;printf("=");}
		printf("\n");fflush(stdout);
	}

	TableReader::TableReader(const char* path, const Vector<const char *> &_collumns){
		Vector< KeyElem<unsigned int, unsigned int> > loc_ibo;
		unsigned int i,j,k,l;
		myHashmap<string, unsigned int> all_needed_colname;
		char buffer[65536];
		char sep;
		unsigned int extra_read =0; unsigned int found=0;
		for(j=0;j<_collumns.size();j++) {
			if ((_collumns[j][0] >= '0')&&(_collumns[j][0] <= '9')){
				for(i= strlen(_collumns[j])-1;i!=0;i--) if ((_collumns[j][i] < '0')||(_collumns[j][i] > '9')) break;
			} else i = 1;
			if (i ==0) loc_ibo.push_back(KeyElem<unsigned int, unsigned int>(atoi(_collumns[j]), j));
			else{
				k=0;
				for(i=0;_collumns[j][i] != '\0';i++){
					if ((_collumns[j][i] == '-') || (_collumns[j][i] == '+')|| (_collumns[j][i] == '/')|| (_collumns[j][i] == '*')|| (_collumns[j][i] == '(')|| (_collumns[j][i] == ')')){
						memcpy(buffer, _collumns[j] + k, i-k);
						buffer[i-k] = '\0';
						l = all_needed_colname.find(string(buffer));
						if (l == 0xFFFFFFFF) {all_needed_colname[string(buffer)] = 0xFFFFFFFF; extra_read++;}
						
						k=i+1;
					}
				}
				if (k == 0){
					l = all_needed_colname.find(string(_collumns[j]));
					if (l == 0xFFFFFFFF) all_needed_colname[string(_collumns[j])] = j;					
					else if (all_needed_colname.deref(l) == 0xFFFFFFFF) {all_needed_colname.deref(l) = j; extra_read--;}
					else found++; // needs copying
				}else{
					l = all_needed_colname.find(string(_collumns[j] + k));
					if (l == 0xFFFFFFFF) {all_needed_colname[string(_collumns[j] + k)] = 0xFFFFFFFF;extra_read++;}
					found++;
				}
			}
		}
		

		
		loc_ibo.sort();
		f = fopen(path, "r+"); if (f ==NULL) {fprintf(stderr, "Could not open file %s!\n", path); exit(1);}
		nbcols = _collumns.size(); nbread =nbcols + extra_read - found;
	//	printf("special collumns: %i, extra reads: %i nb collumns: %i\n",found,extra_read,nbcols );
		read_info = new unsigned int[nbread*2];
		col_name = new string[nbcols + extra_read];
		// fill col_name
		unsigned int cur_num_index=0;
		
		for(i=0;i<_collumns.size();i++){
			if ((cur_num_index >= loc_ibo.size())||(loc_ibo[cur_num_index].k != k)) col_name[i] = string(_collumns[i]);
			else col_name[i] = string("+"); // imposible to match
		}
		j = _collumns.size();
		for(i=0;i<all_needed_colname.heap.size();i++){
			if (all_needed_colname.heap[i].first.d == 0xFFFFFFFF){
				col_name[j] = all_needed_colname.heap[i].first.k;
				all_needed_colname.heap[i].first.d = j++;
			}
		}

		found=0;
		
		read_info[0] =0;k=0;cur_num_index=0;j=0;
		while(2 == fscanf(f,"%[^\t\n\r]%c",buffer, &sep)){
			if ((cur_num_index < loc_ibo.size())&&(loc_ibo[cur_num_index].k == k)) {i = loc_ibo[cur_num_index++].d;col_name[found] = string(buffer);}
			else for(i=0;i<nbcols + extra_read;i++) if (strcmp(buffer, col_name[i].c_str()) == 0) break;
			
			if (i == nbcols + extra_read) read_info[(found << 1)]++;
			else{
				read_info[(found << 1) | 1] = i;
				found++;
				
				if (found == nbread) break;
				read_info[(found << 1)] =0;
			}
			k++;
			
			if ((sep == '\n')||(sep == '\r')) break;
		}
		if (found != nbread){
			printf("missing collumns, found %i out of %i!\n", found , nbread );
			exit(1);
		}
		if (j < loc_ibo.size()) {fprintf(stderr,"table file %s does not a %ith collumn!\n", loc_ibo[j].k);exit(1);}

		while((sep != '\n')&&(sep != '\r'))	fscanf(f,"%[^\t\n\r]%c",buffer, &sep); // flush the resto of the row
		
		current_row = new string[nbcols + extra_read];

		// operattions!
		
		for(j=0;j<_collumns.size();j++) {
			if ((_collumns[j][0] >= '0')&&(_collumns[j][0] <= '9')){
				for(i= strlen(_collumns[j])-1;i!=0;i--) if ((_collumns[j][i] < '0')||(_collumns[j][i] > '9')) break;
			} else i = 1;
			if (i !=0){
				k=0;sep = '\0';
				for(i=0;_collumns[j][i] != '\0';i++){
					if ((_collumns[j][i] == '-') || (_collumns[j][i] == '+')|| (_collumns[j][i] == '/')|| (_collumns[j][i] == '*')|| (_collumns[j][i] == '(')|| (_collumns[j][i] == ')')){
						memcpy(buffer, _collumns[j] + k, i-k);
						buffer[i-k] = '\0';
						oper.push_back(all_needed_colname[string(buffer)]);
						switch(sep){
							case '+': oper.push_back(0x80000001);break;
							case '-': oper.push_back(0x80000002);break;
							case '*': oper.push_back(0x80000003);break;
							case '/': oper.push_back(0x80000004);break;
						}
						sep = _collumns[j][i];
						k = i+1;
					}
				}
				if (k != 0){
					oper.push_back(all_needed_colname[string(_collumns[j] + k)]);
					switch(sep){
						case '+': oper.push_back(0x80000001);break;
						case '-': oper.push_back(0x80000002);break;
						case '*': oper.push_back(0x80000003);break;
						case '/': oper.push_back(0x80000004);break;
					}
					oper.push_back(0x40000000 | j);
				}else{
					if (all_needed_colname[string(_collumns[j])] != j){ // multiple uses... copy!
						oper.push_back(all_needed_colname[string(_collumns[j])]);
						oper.push_back(0x40000000 | j);
					}
				}
			}
		}
		
	//	printf("operations!:\n");for(i=0;i<oper.size();i++) printf("%X\n", oper[i]);fflush(stdout);
	}
	TableReader::~TableReader(){fclose(f);
		delete[](read_info);
		delete[](current_row);
		delete[](col_name);
	}
	
	bool TableReader::nextRow(){
		double dbuffer[16];
		unsigned int i,j;
		char buffer[65536];
		char sep;
		for(j=0;j<nbread;j++){
			for(i=0;i<=read_info[(j << 1)];i++) {
				if (2 != fscanf(f,"%[^\t\n\r]%c",buffer, &sep)){
					if (1 == fscanf(f,"%c", &sep)){
						buffer[0] = '\0';
						if ((sep == '\n')||(sep == '\r')) break;
					}else return false;
				}
			}
			if (i<read_info[(j << 1)]) return false;
			current_row[read_info[(j << 1) | 1]] = string(buffer);
			if ((i == read_info[(j << 1)])&&(j != nbcols-1)) return false;
		}
		while((sep != '\n')&&(sep != '\r'))	if (2 != fscanf(f,"%[^\t\n\r]%c",buffer, &sep)) fscanf(f,"%c", &sep);


		i=0;
		for(j=0;j<oper.size();j++){
			if (oper[j] &0x80000000){
				switch(oper[j] &0x7FFFFFFF){
					case 1: i--;dbuffer[i-1] += dbuffer[i]; break;
					case 2: i--;dbuffer[i-1] -= dbuffer[i]; break;
					case 3: i--;dbuffer[i-1] *= dbuffer[i]; break;
					case 4: i--;dbuffer[i-1] /= dbuffer[i]; break;
				}
			}else if (oper[j] &0x40000000){ // write
				i--;//printf("w\n");fflush(stdout);
				sprintf(buffer, "%.32e", dbuffer[i]); 
				//printf("%s!\n",buffer);
				current_row[(oper[j] &0x3FFFFFFF)] = string(buffer);
			}else{//printf("r\n");fflush(stdout);
				dbuffer[i++] = atof(current_row[oper[j]].c_str()); //printf("%e read!\n",dbuffer[i-1]);
			}
		}
	//	printf("end\n");fflush(stdout);
		return true;
	}

	TiffFile::TiffFile(const char* path, bool writeonly) : curfp_pos(4), endfile_pos(0){
		if (path != NULL){
			f = (writeonly) ? NULL : fopen(path,"rb+");
	char buffer[256];

	if (f == NULL){ // it's a new file, init header!
		f = fopen(path,"wb+");
		if (f ==NULL) {fprintf(stderr,"could not open/create %s!\n", path); return;}
		buffer[0] = 'I';
		buffer[1] = 'I';
		*(short*)(buffer + 2) = 42;
		*(short*)(buffer + 2) = 42;
		*(int*)(buffer + 4) = 0;
		fwrite(buffer,sizeof(char),8,f);

	}else{ // filed exists!
		fread(buffer,sizeof(char),4,f);
		if (buffer[0] == 'M'){
			if (buffer[1] != 'M') {fprintf(stderr,"tried to open %s, with is not a tiff file!\n", path); exit(1);}
			inv =true;
			if (*((short*)(buffer + 2)) != 10752) {fprintf(stderr,"tried to open %s, with is not a tiff file!\n", path); exit(1);}
			} else if (buffer[0]  == 'I'){
			if (buffer[1] != 'I') {fprintf(stderr,"tried to open %s, with is not a tiff file!\n", path); exit(1);}
			inv=false;
		if (*((short*)(buffer + 2)) != 42) {fprintf(stderr,"tried to open %s, with is not a tiff file!\n", path); exit(1);}
		}else fprintf(stderr,"tried to open %s, with is not a tiff file!\n", path);

	}
	} else f = NULL;
}

TiffFile::~TiffFile(){
	if (f != NULL) fclose(f);
}

bool TiffFile::gotoNext(){
	unsigned int i;
	unsigned short s;
	fseek(f, curfp_pos, SEEK_SET);
	fread(&i,sizeof(unsigned int),1,f);
	if (inv) ExOp::bytereverse(i);
	if (i == 0) {fseek(f, -4 * sizeof(char), SEEK_CUR); return(false);}
	fseek(f, i, SEEK_SET);
	fread(&s,sizeof(unsigned short),1,f);
	if (inv) ExOp::bytereverse(s);
//	printf("move to %i, found %i flags\n", i,s);
	curflaglist.setSize(s*12);
	fread(&(curflaglist[0]),sizeof(unsigned char),s*12,f);
	curfp_pos = i + (2+s*12)*sizeof(unsigned char);
	return(true);
}

int TiffFile::flagType(int flagindex){
		unsigned short type = (*(unsigned short*)&(curflaglist[flagindex*12+2]));
		if (inv) ExOp::bytereverse(type);
		int out;
		switch(type){
			case 1:
			case 2:
			case 6:
			case 7:
				out = (*(unsigned int*)&(curflaglist[flagindex*12+4]));
				if (inv) ExOp::bytereverse(out);
				if (out <= 4) out =  type | 0x00000100;
				else out = type;
			break;
			case 3:
			case 8:
				out = (*(unsigned int*)&(curflaglist[flagindex*12+4]));
				if (inv) ExOp::bytereverse(out);
				if (out <= 2) out =  type | 0x00000100;
				else out = type;
			break;
			case 4: // int
			case 9:
			case 11:
				out = (*(unsigned int*)&(curflaglist[flagindex*12+4]));
				if (inv) ExOp::bytereverse(out);
				if (out == 1) out =  type | 0x00000100;
				else out = type;
			break;
		}
	return(out);
	}

template< >
int TiffFile::getValue<int>(int flagindex){
	int type = flagType(flagindex);
	switch(type){
		case 0x00000004:
			fseek(f, *(unsigned int*)&(curflaglist[flagindex*12+ 8]), SEEK_SET);
		break;
		case 0x00000101:
		case 0x00000102: {unsigned char c = *(unsigned char*)&(curflaglist[flagindex*12+ 8]); if (inv) ExOp::bytereverse(c); return(c);}
		case 0x00000103: {unsigned short c = *(unsigned short*)&(curflaglist[flagindex*12+ 8]); if (inv) ExOp::bytereverse(c); return(c);}
		case 0x00000104: {unsigned int c = *(unsigned int*)&(curflaglist[flagindex*12+ 8]); if (inv) ExOp::bytereverse(c); return(c);}
		case 0x00000106: {char c = *(char*)&(curflaglist[flagindex*12+8]); if (inv) ExOp::bytereverse(c); return(c);}
		case 0x00000108: {short c = *(short*)&(curflaglist[flagindex*12+8]); if (inv) ExOp::bytereverse(c); return(c);}
		case 0x00000109: {int c = *(int*)&(curflaglist[flagindex*12+ 8]); if (inv) ExOp::bytereverse(c); return(c);}
	}
	return(0);
}

template< >
unsigned int TiffFile::getValue<unsigned int>(int flagindex){
		int type = flagType(flagindex);
		unsigned int i_out;
		switch(type){
			case 0x00000003: {//:
				unsigned short i_buf;
				fseek(f, *(unsigned int*)&(curflaglist[flagindex*12+8]), SEEK_SET);
				fread(&i_buf,sizeof(unsigned short), 1, f);
				i_out =i_buf;
				}//:
				return(i_out);
			case 0x00000004:
				fseek(f, *(unsigned int*)&(curflaglist[flagindex*12+8]), SEEK_SET);
				fread(&i_out,sizeof(unsigned int), 1, f);
				return(i_out);
			case 0x00000101:
		case 0x00000102: {unsigned char c = *(unsigned char*)&(curflaglist[flagindex*12+8]); if (inv) ExOp::bytereverse(c); return(c);}
		case 0x00000103: {unsigned short c = *(unsigned short*)&(curflaglist[flagindex*12+8]); if (inv) ExOp::bytereverse(c); return(c);}
		case 0x00000104: {unsigned int c = *(unsigned int*)&(curflaglist[flagindex*12+8]); if (inv) ExOp::bytereverse(c); return(c);}
		case 0x00000106: {char c = *(char*)&(curflaglist[flagindex*12+8]); if (inv) ExOp::bytereverse(c); return(c);}
		case 0x00000108: {short c = *(short*)&(curflaglist[flagindex*12+8]); if (inv) ExOp::bytereverse(c); return(c);}
		case 0x00000109: {int c = *(int*)&(curflaglist[flagindex*12+8]); if (inv) ExOp::bytereverse(c); return(c);}
		}
		return(0);
	}

	template< >
	vector<unsigned int> TiffFile::getValue< vector< unsigned int> >(int flagindex){
		int type = flagType(flagindex);
		vector< unsigned int> f_out;
		unsigned int tmptmp;
		unsigned short tmptmps;
		unsigned int i = *(unsigned int*)&(curflaglist[flagindex*12+8]);
		unsigned int arraysize = (*(unsigned int*)&(curflaglist[flagindex*12+4]));
		if (inv) {ExOp::bytereverse(i);ExOp::bytereverse(arraysize);}
		switch(type){
			case 0x00000004:
				fseek(f, i, SEEK_SET);
				for(i=0;i< arraysize ;i++) {
					fread(&tmptmp,sizeof(unsigned int), 1, f);
					if (inv) ExOp::bytereverse(tmptmp);
					f_out.push_back(tmptmp);
				}
				break;
			case 0x00000003:
				fseek(f, i, SEEK_SET);
				for(i=0;i< arraysize ;i++) {
					fread(&tmptmps,sizeof(unsigned short), 1, f);
					if (inv) ExOp::bytereverse(tmptmps);
					f_out.push_back(tmptmps);
				}
				break;
			case 0x00000104:
				f_out.push_back(i);
				break;
			case 0x00000103:
				tmptmps = *(unsigned short*)&(curflaglist[flagindex*12+8]);
				if (inv) ExOp::bytereverse(tmptmps);
				f_out.push_back((unsigned int) tmptmps);
				if ( arraysize  > 1) {
					tmptmps = *(unsigned short*)&(curflaglist[flagindex*12+10]);
					if (inv) ExOp::bytereverse(tmptmps);
					f_out.push_back((unsigned int) tmptmps);
				}
				break;
			case 0x00000102:
				f_out.push_back(*(unsigned char*)&(curflaglist[flagindex*12+8]));
				if ( arraysize > 1) f_out.push_back(*(unsigned char*)&(curflaglist[flagindex*12+9]));
				if ( arraysize  > 2) f_out.push_back(*(unsigned char*)&(curflaglist[flagindex*12+10]));
				if ( arraysize  > 3) f_out.push_back(*(unsigned char*)&(curflaglist[flagindex*12+11]));
				break;
		}
		return(f_out);
	}

	void TiffFile::addopt_Description(Vector< char* > & opt, char const * const descript) const{
		unsigned int j = opt.getSize();
		unsigned int i = strlen(descript);
		opt.push_back(new char[i+5]);
		memcpy(opt[j], "DESC", sizeof(char)*4);
		memcpy(opt[j]+4, descript, sizeof(char)*(i+1));
	}

	void TiffFile::addopt_Xscale(Vector< char* > & opt, double min_x, double max_x) const{
		unsigned int j = opt.getSize();
		opt.push_back(new char[4 + sizeof(double)*2]);
		memcpy(opt[j], "RNGX", sizeof(char)*4);
		memcpy(opt[j]+4, &min_x, sizeof(char)*4);
		memcpy(opt[j]+4+sizeof(double),  &max_x, sizeof(char)*4);
	}
	void TiffFile::addopt_Yscale(Vector< char* > & opt, double min_y, double max_y) const{
		unsigned int j = opt.getSize();
		opt.push_back(new char[4 + sizeof(double)*2]);
		memcpy(opt[j], "RNGY", sizeof(char)*4);
		memcpy(opt[j]+4, &min_y, sizeof(char)*4);
		memcpy(opt[j]+4+sizeof(double),  &max_y, sizeof(char)*4);
	}

    void TiffFile::startFrameWrite(){
        current = new WriteScope();
    }

    void TiffFile::endFrameWrite(){
        delete[](current);

    }
    void TiffFile::savePLYasTIFF(const char * path){
        FILE* f = fopen(path, "r+");
        char buffer[65536];
        unsigned int nbvertex;
        unsigned int nbface;
        // parse header
        while(1 == fscanf(f," %s",buffer)){
            if (strcmp(buffer, "element") == 0){
                fscanf(f," %s",buffer);
                if (strcmp(buffer, "vertex") == 0) fscanf(f," %i", &nbvertex);
                else if (strcmp(buffer, "face") == 0)fscanf(f," %i", &nbface);
            } else if (strcmp(buffer, "end_header") == 0) break;
            }
        
        Tuple<unsigned int, 2> vertex_coor;
        vertex_coor[0] = 8;
        vertex_coor[1] = nbvertex;
        DataGrid<float, 2> vertexmap(vertex_coor);
        unsigned int i,j;
        float* buf;
        float maxcoor =0.0f;
        for(i=0;i<nbvertex;i++) {
            buf = vertexmap.data + i * 8;
            fscanf(f," %f %f %f %f %f %f %f %f", buf+1, buf+2, buf, buf+4, buf+5, buf+3, buf+6, buf+7);
            for(j=0;j<3;j++) if (maxcoor < fabs(buf[j])) maxcoor = fabs(buf[j]);
        }

        for(i=0;i<nbvertex;i++) {
            buf = vertexmap.data + i * 8;
            for(j=0;j<3;j++) buf[j] /= maxcoor;
            /*
            if (fabs(buf[4]) >= 0.9f){
                buf[6] = (buf[0] + 1.0f) * ((buf[4]  > 0.0f) ? 0.5f : -0.5f);
                buf[7] = (buf[2] + 1.0f) * 0.5f;
            }else if (fabs(buf[5]) >= 0.999f){
                buf[6] = (buf[0] + 1.0f) * ((buf[5]  > 0.0f) ? 0.5f : -0.5f);
                buf[7] = (buf[1] + 1.0f) *0.5f;
            }else if (fabs(buf[3]) >= 0.999f){
                buf[6] = (buf[2] + 1.0f) * ((buf[3]  > 0.0f) ? 0.5f : -0.5f);
                buf[7] = (buf[1] + 1.0f) *0.5f;
            }else if (fabs(buf[4]) < sqrt(0.095f)){
                buf[6] = ((float)rand()) / RAND_MAX;
                buf[7] = (buf[1] + 1.0f) *0.5f;
            }*/
        }
        this->put(vertexmap, (float) -1.0f, (float) 1.0f);

        Vector<unsigned int> triangles;
        Vector<unsigned int> quads;
        for(i=0;i<nbface;i++) {
            fscanf(f," %i", &(vertex_coor[0]));
            switch(vertex_coor[0]){
                case 3: for(vertex_coor[1]=0;vertex_coor[1]<3;vertex_coor[1]++) {fscanf(f," %i", &(vertex_coor[0])); triangles.push_back(vertex_coor[0]); } break;
                case 4: for(vertex_coor[1]=0;vertex_coor[1]<4;vertex_coor[1]++) {fscanf(f," %i", &(vertex_coor[0])); quads.push_back(vertex_coor[0]); } break;
                default: fscanf(f,"%*[^\n]"); // flush line
            }
        }
        DataGrid<unsigned int, 2> index_map;
        if (triangles.getSize() != 0){
            vertex_coor[0] =3;
            vertex_coor[1] =triangles.getSize() / 3;
            index_map.setSizes(vertex_coor);
            memcpy(index_map.data, triangles.darray, triangles.getSize() * sizeof(unsigned int) );
            this->put(index_map, (unsigned int) 0, (unsigned int) nbvertex-1);
        }
        if (quads.getSize() != 0){
            vertex_coor[0]= 4;
            vertex_coor[1]= quads.getSize() / 4;
            index_map.setSizes(vertex_coor);
            memcpy(index_map.data, quads.darray, quads.getSize() * sizeof(unsigned int) );
            this->put(index_map, (unsigned int) 0, (unsigned int) nbvertex-1);
        }
        fclose(f);
        printf("Conversion done!, %i vertice, %i triangles %i quads\n", nbvertex, triangles.getSize()/3, quads.getSize()/4);
    }



	void Bluestein::setSize(unsigned int n_size){
		tsize= n_size;
		pre_pow2 =0;
		if (n_size == 0) {static_warning_handdle << LFH_WARNING_UNEXPECTED_INPUT; return;}
		while((tsize & 1) == 0) {pre_pow2++;tsize >>= 1;}
		if (tsize == 1){ // pure power of 2!
			post_pow2 =0;
			mult = 0;
		}else{
			mult = tsize;
			unsigned int t = (tsize << 1) + 1;
			for(post_pow2=0;t != 0;post_pow2++) t >>=1;


			post_pow2 += pre_pow2;
			mult <<= pre_pow2;
			pre_pow2 =0;

			bluewindow.setSize(1 << post_pow2);

			unsigned int i;
			
			for(i=0; i < mult;i++) {
				double ang = i*i*M_PI/ mult;
				bluewindow[i] =  complex(cos(ang),sin(ang));
			}
			for(;((i + mult-1)  >> post_pow2) == 0;i++) bluewindow[i] =  ExCo< complex >::zero();

			for(;(i >> post_pow2) == 0;i++) {
				double ang = ((1 << post_pow2)-i)*((1 << post_pow2)-i)*M_PI/ mult;
				bluewindow[i] =  complex(cos(ang),sin(ang));
			}
			pow2_FFT_routine(bluewindow.darray, post_pow2);
		}
		tsize= n_size;
	}

	unsigned int Bluestein::getBufferSize()const{
		if (mult == 0) return tsize;
		else return tsize - mult + (1 << post_pow2);
	}

	void Bluestein::show(FILE* f, int level)const{
		fprintf(f, "Bluestein! %i,%i,%i,%i;\n", pre_pow2, mult, post_pow2, getBufferSize());
	}

	FoldedGaussianDistribution::FoldedGaussianDistribution() {}
	FoldedGaussianDistribution::FoldedGaussianDistribution(double _mean, double _std): mean(_mean),i_std(sqrt(0.5f)/_std){}

	void FoldedGaussianDistribution::operator()(double &_out_prob, double &obs) const{

		double tmp = (mean - obs) * i_std;
		double tmp2 = (mean + obs) * i_std;
		_out_prob = i_std * (exp(-tmp * tmp) + exp(-tmp2*tmp2)) * (1.0f / sqrt(M_PI));

	}


	double FoldedGaussianDistribution::LL(const double &obs) const{
		double tmp = (mean - obs) * i_std;
		double tmp2 = (mean + obs) * i_std;
		return log(i_std / sqrt(M_PI)) * log(exp(-tmp * tmp) + exp(-tmp2*tmp2));
		// X = log(COSH X + SINH X)

	}



	void FoldedGaussianDistribution::EMinit(){
	}

	void FoldedGaussianDistribution::EMAlphainit(double){
	}

	void FoldedGaussianDistribution::EMregist( const double &instance, const double prob){
	}

	double FoldedGaussianDistribution::EMfinit(){
		return(0.0f);
	}

	FoldedGaussianDistribution FoldedGaussianDistribution::EMnext_exec(double alpha) const{
		FoldedGaussianDistribution ha;
		return ha;
	}

	void FoldedGaussianDistribution::EMclear(){
	}

	void FoldedGaussianDistribution::show(FILE* o) const{
	}



	void GaussianDistributionV::setSizes(unsigned int s){
		if (s == nbsizes) return;
		if (nbsizes != 0){
			delete[](mean);
		}
		mean = new double[s*2]; // mean and var eigevals
		unsigned int coor[2];
		coor[0] =s;
		coor[1] =s;
		eigenvec.setSizes(coor);
		if (stat_Scope != NULL) {
			delete[](stat_Scope);
			stat_Scope = NULL;
		}
		nbsizes = s;
	}


	void GaussianDistributionV::operator()(double & prob, double * & dd) const{

	}

	double GaussianDistributionV::LL( double * const & dd) const{
		return(0.0f);
	}
	void GaussianDistributionV::EMinit(){
		if (stat_Scope == NULL) stat_Scope = new double[2 +  ((nbsizes*(nbsizes+3))/2)];
		for(unsigned int i=0;i<2 + ((nbsizes*(nbsizes+3))/2) ;i++) stat_Scope[i] =0.0f;
	}
	void GaussianDistributionV::EMAlphainit(double){
		if (stat_Scope == NULL) stat_Scope = new double[2 + ((nbsizes*(nbsizes+3))/2)];
		for(unsigned int i=0;i<2 + ((nbsizes*(nbsizes+3))/2) ;i++) stat_Scope[i] =0.0f;
	}
	void GaussianDistributionV::EMregist( double * const  &instance, const double prob){
		unsigned int i,j;

		for(i=0;i<nbsizes;i++) if (!ExOp::isValid(instance[i])) break;
		if (i<nbsizes) return;

		stat_Scope[0] += prob;
		stat_Scope[1] += prob*prob;

		for(i=0;i<nbsizes;i++) stat_Scope[2+i] += instance[i] * prob;
		unsigned int k = nbsizes+2;
		for(i=0;i<nbsizes;i++){
			for(j=i;j<nbsizes;j++,k++) stat_Scope[k] += instance[i] * instance[j] * prob;
		}
	}


	double GaussianDistributionV::EMfinit(){
		unsigned int i;
		unsigned int coor[2];
		unsigned int acoor[2];
	//	for(i=0;i< 2+((nbsizes*(nbsizes+3))/2) ;i++) printf("%f\t", stat_Scope[i]);
		double num = (stat_Scope[1] - stat_Scope[0]*stat_Scope[0]) /2;

		double tmp = 1.0f / stat_Scope[0];
		for(i=0;i<nbsizes;i++) {
			mean[i] = stat_Scope[i+2] * tmp;
		}

	//	DataGrid<double, 2> mat;
		coor[0] = nbsizes;
		coor[1] = nbsizes;
		eigenvec.setSizes(coor);

		i = 2 + nbsizes;
		for(coor[0]=0;coor[0]<nbsizes;coor[0]++){
		//	ExOp::show(coor);
			tmp = stat_Scope[0] * stat_Scope[i] - stat_Scope[2+coor[0]] * stat_Scope[2+coor[0]];
			if (tmp <= 0) tmp = 0.0f;
			coor[1]=coor[0];
			eigenvec(coor) = tmp;
			acoor[1] =coor[0];
			for(coor[1]++,i++;coor[1]<nbsizes;coor[1]++, i++){
		//		ExOp::show(coor);
				tmp = stat_Scope[0] * stat_Scope[i] - stat_Scope[2+coor[1]] * stat_Scope[2+coor[0]];
				eigenvec(coor) = tmp;
				acoor[0] =coor[1];
		//		ExOp::show(acoor);
				eigenvec(acoor) = tmp;
			}
		}
//		(mat / -2.0f * num) is the covariance Matrix!!!

		eigenvec *= (-0.5f /  num);

	//	Vector<double> eig;
	//	eigenvec.makeDiagonalizer_ofinverse(mat, eig);

	//	ExOp::show(mat);
	//	for(i=0;i<nbsizes;i++) {
	//		mean[i+nbsizes] = eig[i] * num;
	//	}

		return (0.0f); // todo
	}

	void GaussianDistributionV::EMclear(){
		if (stat_Scope != NULL) {
			delete[](stat_Scope);
			stat_Scope = NULL;
		}
	}

	void GaussianDistributionV::show(FILE* o) const{
		fprintf(o,"mean: ");
		unsigned int i;
		for(i=0;i<nbsizes;i++) fprintf(o,"%f%c", mean[i], i +1 == nbsizes ? '\n' : '\t');
		fprintf(o,"EigenVar: ");
		for(i=0;i<nbsizes;i++) fprintf(o,"%f%c", mean[i+nbsizes], i +1 == nbsizes ? '\n' : '\t');


		ExOp::show(eigenvec,o,0);

	}


	void GaussianDistributionV::save(FILE* f) const{
		//		fwrite(&normfactor,sizeof(double),1,f);
		//		fwrite(&mean,sizeof(Tuple<double, nbchannels>),1,f);
		//		fwrite(&ihvar,sizeof(Matrix<double, nbchannels, nbchannels>),1,f);
	}

	void GaussianDistributionV::load(FILE* f, unsigned int size){
		//		fread(&normfactor,sizeof(double),1,f);
		//		fread(&mean,sizeof(Tuple<double, nbchannels>),1,f);
		//		fread(&ihvar,sizeof(Matrix<double, nbchannels, nbchannels>),1,f);
	}

	GradientSearchScope::GradientSearchScope(const GradientSearchScope& in) : size(in.size),lastvalue(in.lastvalue), expdiff(in.expdiff), log_alpha(in.log_alpha){
		if (in.lastderiv) {lastderiv = new double[size]; memcpy(lastderiv,in.lastderiv,sizeof(double)*size);}
		else lastderiv = NULL;
	}

	GradientSearchScope& GradientSearchScope::operator=(const GradientSearchScope& in){
		size = in.size;
		if (in.lastderiv) {lastderiv = new double[size]; memcpy(lastderiv,in.lastderiv,sizeof(double)*size);}
		else lastderiv = NULL;
		lastvalue = in.lastvalue;
		expdiff = in.expdiff;
		log_alpha = in.log_alpha;
		return(*this);
	}

	double GradientSearchScope::operator()()const{return(exp(log_alpha));}
	void GradientSearchScope::punish(double log_mag){log_alpha -= log_mag; expdiff *= exp(-log_mag);}

	void GradientSearchScope::init(double initalpha, int in_size){
		size = in_size;
		if (lastderiv != NULL) delete[](lastderiv);
		lastderiv = NULL;
		log_alpha = log(initalpha);
	}

	void GradientSearchScope::registAscent(const double &value, double const * const deriv, double mod_fact){
		int i;
		if (lastderiv == NULL){ // first time, dont update alpha
			lastderiv = new double[size];
			lastvalue = value;
			expdiff  = 0;
			memcpy(lastderiv,deriv,sizeof(double)*size);
			for(i=0;i<size;i++) expdiff  += deriv[i]*deriv[i];
			expdiff = lastvalue + exp(log_alpha) * expdiff;
		}else{
			// update alpha!
			double sum, suma ,sumb, tmp, fact;
			tmp = (value - expdiff); sum = (tmp * tmp);
			tmp = (lastvalue - expdiff); suma = (tmp * tmp); // magnitude of prediction
			tmp = (value - lastvalue); sumb = tmp * tmp;
			fact = exp(- 0.7f * sum * pow(suma*sumb,-0.5f));
			sum =0;
			suma =0;
			sumb =0;
			lastvalue = value;
			expdiff = 0;
			for(i=0;i<size;i++) {tmp = lastderiv[i] - deriv[i]; sum += tmp*tmp; suma += deriv[i]*deriv[i]; sumb += lastderiv[i] * lastderiv[i]; expdiff += deriv[i]*deriv[i];}
			log_alpha += fact + exp(- 0.7f * sum * pow(suma*sumb,-0.5f)) - 1.0f ; // multiply by exp(1) or exp(-1) in extreme cases
			expdiff = lastvalue + exp(log_alpha) *expdiff;
		}
	}
	void GradientSearchScope::registDescent(const double &value, double const * const deriv, double mod_fact){
		int i;
		if (lastderiv == NULL){ // first time, dont update alpha
			lastderiv = new double[size];
			expdiff = 0;
			for(i=0;i<size;i++) expdiff += deriv[i]*deriv[i];
		}else{
			double sump, modif;
			
			sump =0;
			expdiff = 0;
			for(i=0;i<size;i++) {expdiff += deriv[i] * deriv[i]; sump += deriv[i] * lastderiv[i];}

	//		printf("F[x] = %e; dF[x]/dx = %e; dF[x]/dy = %e; F[x+s]-F[x] = %e; GR = %e; old deriv = %e; new deriv = %e \n", lastvalue, lastderiv[0], lastderiv[1], value - lastvalue, -lexpdiff * exp(log_alpha) * mod_fact, lexpdiff , sump);
	//		printf("curvature! = %e, stepsize %e ,log alpha incr : %e\n", 1.0f / fabs(lexpdiff - sump), exp(log_alpha) * (sqrt(lexpdiff) * mod_fact), 0.5f * log(lexpdiff) + log(mod_fact) - log(fabs(lexpdiff - sump)) );


			//		modif = -0.5f * log( fabs( 1.0f - (sump / expdiff) )); // geometric mean update
					modif = -0.5f * log( 0.1f + fabs( 1.0f - (sump / expdiff) )); // geometric mean update

			if (ExCo<double>::isValid(modif)){
				log_alpha +=  modif;
			}
	//		printf("curvature! = %e   %e\n",log_alpha, modif);

			// F[-S/2] = lastvalue, F[S/2] = value F'[-S] = lexpdiff  F[S] = expdiff

			// F[-S/2] = lastvalue, F[S/2] = value F'[-S] = lexpdiff  F[S] = expdiff

			//  0.5f * (lastvalue + value) = A + C * S^2 // useless!
			//  0.5f * (-lastvalue + value) = B*S + D * S^3
			//  0.5f * (lexpdiff + expdiff) = B * 3D*S^2
			//   (-lexpdiff + expdiff) = C*S //!

			// 1/2C ~=  2.0f / (-lexpdiff + expdiff)




		}
		lastvalue = value;
		memcpy(lastderiv,deriv,sizeof(double)*size);
	}

	double GradientSearchScope::updateAscent(const double &value, double * guess,  double const * const deriv){// returns a safe bound on the parameter changes;
		int i;
		double sums,  sump, modif;
		if (lastderiv == NULL){ // first time, dont update alpha
			lastderiv = new double[size];
			expdiff = 0.0f;
			for(i=0;i<size;i++) {expdiff += deriv[i] * deriv[i];}
		}else if (value <= lastvalue){
               // reject point!, use cubic fit to estimate minimum between previous point an this point F[0] = lastvalue F[S] = value F'[0]
                sump =0;
                for(i=0;i<size;i++) {sump += deriv[i] * lastderiv[i];}
         //       printf("super scope! f0 = %f f1 = %f d0 = %f d1 = %f\n", lastvalue, value, expdiff, sump);
                // distance between points:   alpha *expdiff

                double tmp[3];
                sums = exp(log_alpha);
                tmp[0] = -lastvalue -  sums* expdiff + value; //objective
                tmp[1] = sump - expdiff; //objective // sqrt(expdiff) premulti
          //      printf("%f,%f\n", tmp[0],tmp[1]);

                tmp[2] = -(3.0f * tmp[0] - tmp[1]* sums);
                tmp[0] = 3.0f * (2.0f * tmp[0]- tmp[1]* sums); // sqrt(expdiff) premulti
           //     printf("x2 =%f  x3 = %f in absolute direction\n", tmp[2], tmp[3]);
                
          //      printf("deriv discriminant, be positive :/ = %f\n", tmp[4]);
                tmp[1] = -tmp[2] / tmp[0];
                if (fabs(tmp[1]) > 1000000000.0f){ // pivot is too large!
                    tmp[0] = sums * expdiff / (tmp[2] * 2.0f);
                    
                }else{
                tmp[0] = tmp[1] + (sqrt(tmp[2] * tmp[2] + tmp[0] * expdiff * sums) / tmp[0]);
           //     printf("altroot: %f\n", tmp[5] - (sqrt(tmp[4]) / tmp[3]));
                }
                printf("in [0,1] range, desired root = %f\n", tmp[0]); // pick same sign as x3
               
          //      if ((tmp[0] < 0.0f)||(tmp[0]> 1.0f)) printf("NO way!!!!\n");

                //printf("recover guess: ");
                sums *= (1.0f - tmp[0]); 
                for(i=0;i<size;i++) guess[i] -= sums * lastderiv[i];
                // alpha = alpha * tmp8 -> logalpha += log(tmp[8])
                log_alpha += log(tmp[0]);
                return(ExCo<double>::maximum());
        }else{
			sums = expdiff;
			sump =0;
			expdiff = 0;
			for(i=0;i<size;i++) {expdiff += deriv[i] * deriv[i]; sump += deriv[i] * lastderiv[i];}
			//		modif = -0.5f * log( fabs( 1.0f - (sump / sums) )); // geometric mean update
			modif = -0.5f * log( 0.1f + fabs( 1.0f - (sump / sums) )); // geometric mean update
			if (ExCo<double>::isValid(modif)){
				log_alpha +=  modif;
			}
		}
		lastvalue = value;
		memcpy(lastderiv,deriv,sizeof(double)*size);
		sums = exp(log_alpha);
		modif = 0.0f;
		for(i=0;i<size;i++) {sump = guess[i] ; guess[i] += sums * deriv[i]; modif += (fabs(guess[i]) > 1.2e-200) ? fabs((guess[i] - sump)/ guess[i]) : 0.0f ;}
		return(modif /size);
	}

	double GradientSearchScope::updateDescent(const double &value, double * guess , double const * const deriv){ // returns a safe bound on the parameter changes;
		int i;
		double sums,  sump, modif;
		if (lastderiv == NULL){ // first time, dont update alpha
			lastderiv = new double[size];
			expdiff = 0.0f;
			for(i=0;i<size;i++) {expdiff += deriv[i] * deriv[i];}
		}else if (value >= lastvalue) {
               // reject point!, use cubic fit to estimate minimum between previous point an this point F[0] = lastvalue F[S] = value F'[0]
                sump =0;
                for(i=0;i<size;i++) {sump += deriv[i] * lastderiv[i];}
         //       printf("super scope! f0 = %f f1 = %f d0 = %f d1 = %f\n", lastvalue, value, expdiff, sump);
                // distance between points:   alpha *expdiff

                double tmp[3];
                sums = exp(log_alpha);
                tmp[0] = lastvalue -  sums* expdiff - value; //objective
                tmp[1] = sump - expdiff; //objective // sqrt(expdiff) premulti
                
                tmp[2] = -(3.0f * tmp[0] - tmp[1]* sums);
                tmp[0] = 3.0f * (2.0f * tmp[0]- tmp[1]* sums); // sqrt(expdiff) premulti
           //     printf("x2 =%f  x3 = %f in absolute direction\n", tmp[2], tmp[3]);
                
          //      printf("deriv discriminant, be positive :/ = %f\n", tmp[4]);
                tmp[1] = -tmp[2] / tmp[0];
                if (fabs(tmp[1]) > 1000000000.0f){ // pivot is too large!
                    tmp[0] = sums * expdiff / (tmp[2] * 2.0f);
                }else{
                tmp[0] = tmp[1] + (sqrt(tmp[2] * tmp[2] + tmp[0] * expdiff * sums) / tmp[0]);
           //     printf("altroot: %f\n", tmp[5] - (sqrt(tmp[4]) / tmp[3]));
                }
               // printf("in [0,1] range, desired root = %f\n", tmp[0]); // pick same sign as x3
                
                if ((tmp[0] < 0.0f)||(tmp[0]> 1.0f)) printf("NO way!!!!\n");

                //printf("recover guess: ");
                sums *= (1.0f - tmp[0]); 
                for(i=0;i<size;i++) guess[i] += sums * lastderiv[i];
                // alpha = alpha * tmp8 -> logalpha += log(tmp[8])
                log_alpha += log(tmp[0]);
                return(ExCo<double>::maximum());
        }else{
			sums = expdiff;
			sump =0;
			expdiff = 0;
			for(i=0;i<size;i++) {expdiff += deriv[i] * deriv[i]; sump += deriv[i] * lastderiv[i];}
			//		modif = -0.5f * log( fabs( 1.0f - (sump / sums) )); // geometric mean update
	//		printf("%e expected, %e obtained\n", );
			modif = -0.5f * exp(log_alpha) - ((value - lastvalue) / (sump + sums)); // cubic error
			modif = -0.5f * log( fabs( modif) + fabs( 1.0f - (sump / sums) )); // geometric mean update
			if (ExCo<double>::isValid(modif)){
				log_alpha +=  modif;
			}

		}
		lastvalue = value;
		memcpy(lastderiv,deriv,sizeof(double)*size);
		sums = exp(log_alpha);
		modif = 0.0f;
		for(i=0;i<size;i++) {sump = guess[i] ; guess[i] -= sums * deriv[i]; modif += (fabs(guess[i]) > 1.2e-200) ? fabs((guess[i] - sump)/ guess[i]) : 0.0f ;}
		return(modif /size);
	}



	CurvatureSearchScope::CurvatureSearchScope(const CurvatureSearchScope& in) : size(in.size),lastvalue(in.lastvalue), expdiff(in.expdiff){
		if (in.lastderiv) {lastderiv = new double[size]; memcpy(lastderiv,in.lastderiv,sizeof(double)*size);
            curvature = new double[size & 1 ? ((size >> 1)+1)*(size) : (size >> 1)*(size+1) ]; memcpy(lastderiv,in.lastderiv,sizeof(double)*size);
        }
		else lastderiv = NULL;
	}

	CurvatureSearchScope& CurvatureSearchScope::operator=(const CurvatureSearchScope& in){
        size = in.size;
		if (in.lastderiv) {lastderiv = new double[size]; memcpy(lastderiv,in.lastderiv,sizeof(double)*size);}
		else lastderiv = NULL;
		lastvalue = in.lastvalue;
		expdiff = in.expdiff;
		log_alpha = in.log_alpha;
shrinkcount = in.shrinkcount;
		return(*this);
	}


	void CurvatureSearchScope::init(double initalpha, unsigned int in_size){
		if ((curvature)&&(size != in_size)) {delete[](curvature); curvature = NULL;}
        size = in_size;
		if (lastderiv != NULL) delete[](lastderiv);
		lastderiv = NULL;
		if (curvature == NULL) curvature = new double[size & 1 ? ((size >> 1)+1)*(size) : (size >> 1)*(size+1) ];
        unsigned int i,j,k;
        for(i=0,k=0;i<size;i++) {
            curvature[k] = initalpha;
            for(j=i+1,k++;j<size;j++,k++) curvature[k] = 0.0f;
        } 
shrinkcount =0;
	}

	double CurvatureSearchScope::updateAscent(const double &value, double * guess,  double const * const deriv){// returns a safe bound on the parameter changes;
		unsigned int i,j,k;
		double sums,  sump, modif;

		if (lastderiv == NULL) lastderiv = new double[size*3];
		else{

			// check Wolfe Conditions!
			
			// value - lastvalue ≤ 0.0001(guess - oldguess) \cdot oldderiv 
			// (guess - oldguess) \cdot deriv ≥ 0.9 (guess - oldguess) \cdot oldderiv
            sums =0; sump =0;
            for(i=0;i<size;i++) {sums -= lastderiv[i] * (guess[i] - lastderiv[i+size]);
                                     sump -= deriv[i] * (guess[i] - lastderiv[i+size]);}
         //   printf("wolfe 1: %e <= %e  %c\n",  lastvalue - value, 0.0001f * sums, (lastvalue-value <= 0.0001f * sums) ? 'Y' : 'N' );
       //     printf("wolfe 2: %e >= %e  %c\n", sump, 0.9f * sums , (sump>= 0.9f * sums) ? 'Y' : 'N' );
        if (value - lastvalue  < 0.0001f * sums){  // step too long
            // scaling inverse hessian by 0.25
			// sums is the denominator!
            
		//	printf("Shrink!\n");
            shrinkcount++;
            
            
            sums = -0.75f / fabs(sums);
            sump=0;
            for(i=0;i<size;i++) guess[i] -= lastderiv[i+size];
			
            for(i=0,k=0;i<size;i++) {
                curvature[k] += sums * guess[i] * guess[i]; 
                for(j=i+1,k++;j<size;j++,k++){
                    curvature[k] += sums * guess[i] * guess[j];
                }
            }
			//        if ((shrinkcount & 7) == 7) {for(i=0;i<size;i++) lastderiv[i] = -lastderiv[i];}
			/*
			 for(i=0,k=0;i<size;i++) {
			 curvature[k] *= 0.25f; 
			 for(j=i+1,k++;j<size;j++,k++){
			 curvature[k] *= 0.25f; 
			 }
			 }*/
            memcpy(guess,lastderiv+size,sizeof(double)*size);
			//   printf("curvature:");
			
            modif = 0.0f;
			
            for(i=0,k=0;i<size;i++) {
                guess[i] += curvature[k] * lastderiv[i];//printf("%e\t",curvature[k] );
                for(j=i+1,k++;j<size;j++,k++){
                    guess[j] += curvature[k] * lastderiv[i];//printf("%e\t",curvature[k] );
                    guess[i] += curvature[k] * lastderiv[j];
					
                }
                modif += (fabs(guess[i]) > 1.2e-100) ? fabs((guess[i] - lastderiv[i+size])/ guess[i]) : 0.0f ;
            }//printf("\n");
			
			
            return( modif );
        }else if (sump < 0.9f * sums){ // step too short
            // scaling inverse hessian by 4
            // sums is the denominator!
			//     printf("Expand!\n");
            sums = 3.0f / fabs(sums);
            shrinkcount =0;
            for(i=0;i<size;i++) lastderiv[i+size] -= guess[i];
			
            for(i=0,k=0;i<size;i++) {
                curvature[k] += sums * (lastderiv[i+size] * lastderiv[i+size]);// printf("%f\t", curvature[k]);
                for(j=i+1,k++;j<size;j++,k++){
                    curvature[k] += sums * (lastderiv[j+size] * lastderiv[i+size]); //printf("%f\t", curvature[k]);
                }
            }
			
			
			memcpy(lastderiv + size ,guess,sizeof(double)*size); 
			//    printf("curvature:");
			for(i=0,k=0;i<size;i++) {
				guess[i] += curvature[k] * deriv[i];//printf("%e\t",curvature[k] );
				for(j=i+1,k++;j<size;j++,k++){
					guess[j] += curvature[k] * deriv[i];
					guess[i] += curvature[k] * deriv[j];//printf("%e\t",curvature[k] );
				}
				
			}// printf("\n");
			lastvalue = value;
			memcpy(lastderiv,deriv,sizeof(double)*size); 
			
			return(ExCo<double>::maximum());
        }else { shrinkcount =0;
            // if condition1 fails, the step is too long: it is impossible to update the hessian
			// if condition2 fails, the step is too short: it is impossible to update the hessian
     //   printf("AlRight!\n");

        for(i=0;i<size;i++) lastderiv[i+size] -= guess[i];
        i=0;k=0;

        lastderiv[i+size*2] = curvature[k] * deriv[i];

        for(j=i+1,k++;j<size;j++,k++){
            lastderiv[j+size*2] = curvature[k] * deriv[i];
            lastderiv[i+size*2] += curvature[k] * deriv[j];
   //         printf("%f\t", curvature[k]);
        }
        for(i++;i<size;i++) {
            lastderiv[i+size*2] += curvature[k] * deriv[i];
            for(j=i+1,k++;j<size;j++,k++){
                lastderiv[j+size*2] += curvature[k] * deriv[i];
                lastderiv[i+size*2] += curvature[k] * deriv[j];
            }
        }

        sums -= sump;
        sump = lastderiv[size*2] * (deriv[0] - lastderiv[0]);
        for(i=1;i<size;i++) sump += lastderiv[i+size*2] * (deriv[i] - lastderiv[i]);
        sums = -1.0f / sums;
        sump *= sums;
   //     printf("K2 = %f\tK1 = %f\n", sump,sums);
        for(i=0,k=0;i<size;i++) {
            curvature[k] += sums * (lastderiv[i+size] * (sump * lastderiv[i+size] - 2.0f * lastderiv[i+size*2]));
            for(j=i+1,k++;j<size;j++,k++){
                curvature[k] += sums *(sump * lastderiv[i+size]* lastderiv[j+size] - (lastderiv[i+size]* lastderiv[j+size*2] + lastderiv[i+size*2]* lastderiv[j+size])); 
            }
        }

		} 
        }
        // has updated curvature by now
        modif = (fabs(value) > 1.2e-100) ? fabs((value - lastvalue) / (value)) : 0.0f;
        memcpy(lastderiv + size ,guess,sizeof(double)*size); 
    //   printf("curvature:");
        for(i=0,k=0;i<size;i++) {
            guess[i] += curvature[k] * deriv[i];//printf("%e\t",curvature[k] );
            for(j=i+1,k++;j<size;j++,k++){
                guess[j] += curvature[k] * deriv[i];
                guess[i] += curvature[k] * deriv[j];//printf("%e\t",curvature[k] );
            }
            modif += (fabs(guess[i]) > 1.2e-100) ? fabs((guess[i] - lastderiv[i+size])/ guess[i]) : 0.0f ;
        } //printf("\n");
        lastvalue = value;
        memcpy(lastderiv,deriv,sizeof(double)*size); 
    //    printf("%e\n",modif /size);
        
		return(modif /size);


	}

	double CurvatureSearchScope::updateDescent(const double &value, double * guess , double const * const deriv){ // returns a safe bound on the parameter changes;
		unsigned int i,j,k;
		double sums,  sump, modif;
		if (lastderiv == NULL) lastderiv = new double[size*3];
		else{

			// check Wolfe Conditions!
			
			// value - lastvalue ≤ 0.0001(guess - oldguess) \cdot oldderiv 
			// (guess - oldguess) \cdot deriv ≥ 0.9 (guess - oldguess) \cdot oldderiv
            sums =0; sump =0;
            for(i=0;i<size;i++) {sums += lastderiv[i] * (guess[i] - lastderiv[i+size]);
                                     sump += deriv[i] * (guess[i] - lastderiv[i+size]);}
        //    printf("wolfe 1: %e <= %e  %c\n",  value - lastvalue, 0.0001f * sums, (value - lastvalue <= 0.0001f * sums) ? 'Y' : 'N' );
       //     printf("wolfe 2: %e >= %e  %c\n", sump, 0.9f * sums , (sump>= 0.9f * sums) ? 'Y' : 'N' );
       //      printf("wolfe 1: %e <= %e  %c\n", value - lastvalue, 0.0001f * sums, (value - lastvalue <= 0.0001f * sums) ? 'Y' : 'N' );
       //     printf("wolfe 2: %e >= %e  %c\n", sump, 0.9f * sums , (sump>= 0.9f * sums) ? 'Y' : 'N' );

        if (value - lastvalue > 0.0001f * sums){  // step too long
            // scaling inverse hessian by 0.25
			// sums is the denominator!
            
			//   printf("desc Shrink!\n");
            sums = -0.75f / fabs(sums);
            sump=0;
            for(i=0;i<size;i++) guess[i] -= lastderiv[i+size];
			
            for(i=0,k=0;i<size;i++) {
                curvature[k] += sums * guess[i] * guess[i];// printf("%f\t", curvature[k]);
                for(j=i+1,k++;j<size;j++,k++){
                    curvature[k] += sums * guess[i] * guess[j]; //printf("%f\t", curvature[k]);
                }
            }
			
            memcpy(guess,lastderiv+size,sizeof(double)*size);
			//     printf("curvature:");
            for(i=0,k=0;i<size;i++) {
                guess[i] -= curvature[k] * lastderiv[i];//printf("%e\t",curvature[k] );
                for(j=i+1,k++;j<size;j++,k++){
                    guess[j] -= curvature[k] * lastderiv[i];//printf("%e\t",curvature[k] );
                    guess[i] -= curvature[k] * lastderiv[j];
                }
            }//printf("\n");
			
			
            return(ExCo<double>::maximum());
        }else if (sump < 0.9f * sums){ // step too short
            // scaling inverse hessian by 4
            // sums is the denominator!
			//     printf("desc Expand!\n");
            sums = 3.0f / fabs(sums);
            
            for(i=0;i<size;i++) lastderiv[i+size] -= guess[i];
			
            for(i=0,k=0;i<size;i++) {
                curvature[k] += sums * (lastderiv[i+size] * lastderiv[i+size]);// printf("%f\t", curvature[k]);
                for(j=i+1,k++;j<size;j++,k++){
                    curvature[k] += sums * (lastderiv[j+size] * lastderiv[i+size]); //printf("%f\t", curvature[k]);
                }
            }
			
			
			memcpy(lastderiv + size ,guess,sizeof(double)*size); 
			//    printf("curvature:");
			for(i=0,k=0;i<size;i++) {
				guess[i] -= curvature[k] * deriv[i];//printf("%e\t",curvature[k] );
				for(j=i+1,k++;j<size;j++,k++){
					guess[j] -= curvature[k] * deriv[i];
					guess[i] -= curvature[k] * deriv[j];//printf("%e\t",curvature[k] );
				}
				
			} //printf("\n");
			lastvalue = value;
			memcpy(lastderiv,deriv,sizeof(double)*size); 
			
			return(ExCo<double>::maximum());
        }else {
            // if condition1 fails, the step is too long: it is impossible to update the hessian
			// if condition2 fails, the step is too short: it is impossible to update the hessian

//printf("desc Right!\n");
        for(i=0;i<size;i++) lastderiv[i+size] -= guess[i];
        i=0;k=0;

        lastderiv[i+size*2] = curvature[k] * deriv[i];

        for(j=i+1,k++;j<size;j++,k++){
            lastderiv[j+size*2] = curvature[k] * deriv[i];
            lastderiv[i+size*2] += curvature[k] * deriv[j];
   //         printf("%f\t", curvature[k]);
        }
        for(i++;i<size;i++) {
            lastderiv[i+size*2] += curvature[k] * deriv[i];
            for(j=i+1,k++;j<size;j++,k++){
                lastderiv[j+size*2] += curvature[k] * deriv[i];
                lastderiv[i+size*2] += curvature[k] * deriv[j];
            }
        }

        sums -= sump;
        sump = lastderiv[size*2] * (deriv[0] - lastderiv[0]);
        for(i=1;i<size;i++) sump += lastderiv[i+size*2] * (deriv[i] - lastderiv[i]);
        sums = 1.0f / sums;
        sump *= sums;
   //     printf("K2 = %f\tK1 = %f\n", sump,sums);
        for(i=0,k=0;i<size;i++) {
            curvature[k] += sums * (lastderiv[i+size] * (sump * lastderiv[i+size] - 2.0f * lastderiv[i+size*2]));
            for(j=i+1,k++;j<size;j++,k++){
                curvature[k] += sums *(sump * lastderiv[i+size]* lastderiv[j+size] - (lastderiv[i+size]* lastderiv[j+size*2] + lastderiv[i+size*2]* lastderiv[j+size])); 
            }
        }

		} 
        }
        modif = (fabs(value) > 1.2e-100) ? fabs((value - lastvalue) / (value)) : 0.0f;
        memcpy(lastderiv + size ,guess,sizeof(double)*size); 
     //    printf("curvature:");
        for(i=0,k=0;i<size;i++) {
            guess[i] -= curvature[k] * deriv[i];//printf("%e\t",curvature[k] );
            for(j=i+1,k++;j<size;j++,k++){
                guess[j] -= curvature[k] * deriv[i];
                guess[i] -= curvature[k] * deriv[j];//printf("%e\t",curvature[k] );
            }
            modif += (fabs(guess[i]) > 1.2e-100) ? fabs((guess[i] - lastderiv[i+size])/ guess[i]) : 0.0f ;
        } //printf("\n");
        lastvalue = value;
        memcpy(lastderiv,deriv,sizeof(double)*size); 
    //    printf("%e\n",modif /size);
        
		return(modif /size);
	}
	
} // end of namespace
