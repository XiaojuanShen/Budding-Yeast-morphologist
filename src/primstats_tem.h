/*
 * primstats_tem.hpp
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


	template<class C, unsigned int nbstates>
	void Classifier<C,nbstates>::setUnknownProbability(double val){
		if ((val == 0.0f)^(unknown_scope == NULL)){
			if (val == 0) {delete[](unknown_scope); unknown_scope = NULL;}
			else{
				unknown_scope = new double[8];
				unknown_scope[0] = 0.0f;
				unknown_scope[1] = val;
				unknown_scope[2] = 0.0f; // step_size
			}
		}
	}

	template<class C, unsigned int nbstates>
	double Classifier<C,nbstates>::UnknownProbability(C &instance, const Tuple<double, nbstates>& prob){
		int i;
		double sum =0.0f;
		double tmp;
		for(i=0;i<nbstates;i++){(*classes[i])(tmp,instance);
			if (tmp != 0.0f) sum += prob[i] * log(tmp);
		}
		return(unknown_scope[0] / (exp(sum) + unknown_scope[0]));
	}


	template<class C, unsigned int nbstates>
	void Classifier<C,nbstates>::operator()(Tuple<double, nbstates>& _out, C& _in )const{
		Tuple<double, nbstates> tmp;
		int i;
		for(i=0;i<nbstates;i++){
			(*classes[i])(tmp[i],_in);
			/*	if (!ExCo<double>::isValid(tmp[i])) {

			 ((GaussianDistribution<2> *)classes[0])->show(stdout);
			 ((GaussianDistribution<2> *)classes[1])->show(stdout);

			 printf("%f\t%f\n", _in[0], _in[1]);
			 printf("%f\n", tmp[i]);
			 printf("AN EARLY ERROR!\n");exit(1);
			 }*/
		}
		// tmp[nbstates] = unknown_scope == NULL ? 0.0f : unknown_scope[0];
		_out = tmp.normalize();


		for(i=0;i<nbstates;i++){


		}

	}

	template<class C, unsigned int nbstates>
	void Classifier<C,nbstates>::EMinit(){
		unsigned int i; for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMinit();
		if (unknown_scope != NULL){
			memset(unknown_scope+3,'\0',sizeof(double)*5);
		}
	}
	template<class C, unsigned int nbstates>
	void Classifier<C,nbstates>::EMAlphainit(double alp){
		unsigned int i; for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMAlphainit(alp);
		if (unknown_scope != NULL){
			memset(unknown_scope+3,'\0',sizeof(double)*5);
		}
	}
	template<class C, unsigned int nbstates>
	void Classifier<C,nbstates>::EMregist(C &instance, const Tuple<double, nbstates> prob){
		int i;
		double pix[8];
		if (unknown_scope == NULL){
			for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMregist(instance, prob[i]);
		}else{
			// the probability of unknown is only defined by the anterior
			double tmp,sum,tmp2;
			sum =0.0f;
			for(i=0;i<nbstates;i++){(*classes[i])(tmp,instance);
				if (tmp > 0.0f) sum += prob[i] * log(tmp);
			}

			sum = exp(sum);

			tmp = sum / (sum + unknown_scope[0]);





			/*
			 unknown_scope[4] += 1.0f / (1.0f + exp(3.0f * unknown_scope[2]) * sum);
			 unknown_scope[5] += 1.0f / (1.0f + exp(unknown_scope[2]) * sum);
			 unknown_scope[6] += 1.0f / (1.0f + exp(-unknown_scope[2]) * sum);
			 unknown_scope[7] += 1.0f / (1.0f + exp(-3.0f * unknown_scope[2]) * sum);
			 */
			if (unknown_scope[2] > 0.125f){
				for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMregist(instance, prob[i] * tmp);
				// far approach
				tmp2 = sum / unknown_scope[0];
				sum = unknown_scope[0] / sum;
				if ((!ExCo<double>::isValid(sum))||(!ExCo<double>::isValid(tmp2))){
					unknown_scope[3] += 1.0f;
					unknown_scope[4] += 2.0f;
					unknown_scope[6] += 2.0f;
				}else {

					tmp = 2.0f + (sum + tmp2) / cosh(3.0f*unknown_scope[2]);
					tmp2 = 2.0f + (sum + tmp2) / cosh(unknown_scope[2]);

					if ((ExCo<double>::isValid(tmp))&&(ExCo<double>::isValid(tmp2))){

						if (ExCo<double>::isValid(2.0f * sinh(unknown_scope[2]) / tmp2)){

							pix[0] = (1.0f + (sum / cosh(3.0f*unknown_scope[2]))) / tmp;
							pix[1] = tanh(3.0f*unknown_scope[2]) / tmp;
							pix[2] = (1.0f + (sum/cosh(unknown_scope[2]))) / tmp2;
							pix[3] = tanh(unknown_scope[2]) / tmp2;

							if ((ExCo<double>::isValid(pix[0]))&&(ExCo<double>::isValid(pix[1]))&&(ExCo<double>::isValid(pix[2]))&&(ExCo<double>::isValid(pix[3]))){
								unknown_scope[4] += pix[0];
								unknown_scope[5] += pix[1];
								unknown_scope[6] += pix[2];
								unknown_scope[7] += pix[3];
								unknown_scope[3] += 1.0f;
							}
						}
					}}
			}else if (unknown_scope[2] == 0.0f){
				// passive guess
				for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMregist(instance, prob[i]);
				if (ExCo<double>::isValid(log(sum))){
					unknown_scope[3] += 1.0f;
					unknown_scope[4] += log(sum);

				}
			}else {
				//close approach
				for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMregist(instance, prob[i] * tmp);
				sum /= unknown_scope[0];

				tmp = 1.0f + sum;

				pix[0] = 1.0f / tmp;
				double tmp2 = 2.0f + sum + (1.0f / sum);
				pix[1] = 1.0f / tmp2;
				pix[2] = 1.0f / (tmp2 * tmp);
				pix[3] = 1.0f / (tmp2 * tmp * tmp);

				if ((ExCo<double>::isValid(pix[0]))&&(ExCo<double>::isValid(pix[1]))&&(ExCo<double>::isValid(pix[2]))&&(ExCo<double>::isValid(pix[3]))){

					unknown_scope[3] += 1.0f;
					unknown_scope[4] += pix[0];

					unknown_scope[5] += pix[1];
					unknown_scope[6] += pix[2];
					unknown_scope[7] += pix[3];
				}
			}
		}
	}
	template<class C, unsigned int nbstates>
	void Classifier<C,nbstates>::EMclear(){
		int i; for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->clear();
	}

	template<class C, unsigned int nbstates>
	double Classifier<C,nbstates>::EMfinit(){
		int i;
		double LL_out = 0.0f;
		if (unknown_scope != NULL){
			double coef[4];

			if (unknown_scope[2] > 0.125f){
				// far approach
				//			printf("Far Step:\n");
				//		printf("%e\t%e\t%e\t%e\n", (unknown_scope[4]-unknown_scope[5])/(unknown_scope[3]), (unknown_scope[6]-unknown_scope[7])/(unknown_scope[3]),(unknown_scope[6]+unknown_scope[7])/(unknown_scope[3]),(unknown_scope[4]+unknown_scope[5])/(unknown_scope[3]));

				coef[0] =  (-unknown_scope[4] + 9.0f * unknown_scope[6]) - 8.0f * unknown_scope[3] * unknown_scope[1];
				coef[1] =  ((unknown_scope[5]/-3.0f) + 9.0f * unknown_scope[7]) ;
				coef[2] =  (unknown_scope[4] - unknown_scope[6]);
				coef[3] =  ((unknown_scope[5]/3.0f) - unknown_scope[7]);
				//		printf("%e\t%e\t%e\t%e\n", coef[0], coef[1], coef[2], coef[3]);
				double shift = CubicRealRoot(coef, false);
				//printf("Far Step: F[%e] = %e!\n",unknown_scope[0], unknown_scope[6]/(unknown_scope[3]));
				//		printf("%e\t%e\n", shift, unknown_scope[2]);



				shift *=antioss;
				if (!ExCo<double>::isValid(shift)){
					unknown_scope[0] *= exp((unknown_scope[6] > unknown_scope[3]* unknown_scope[1] ? -3.0f : 3.0f)* unknown_scope[2]);

				}else if (fabs(shift) > 3.0f){ // outside!
					unknown_scope[0] *= exp((shift < 0 ? -3.0f : 3.0f)* unknown_scope[2]);
					unknown_scope[2] *= 4.0f;

				}else{
					unknown_scope[0] *= exp(shift * unknown_scope[2]);
					unknown_scope[2] *= fabs(unknown_scope[1] -(unknown_scope[5]+ unknown_scope[6])/(unknown_scope[3]));
				}
				//		printf("end: %e\t    %e\n", unknown_scope[0], unknown_scope[2]);
				relerr = log(unknown_scope[6]) - log(unknown_scope[3] * unknown_scope[1]);

				if (relerr < 1.0f){ // dont update if too much data is in unknown class
					for(i=0;i<nbstates;i++) LL_out += ((ProbabilityDistribution<C>*)classes[i])->EMfinit();
					LL_out /= (1.0f - unknown_scope[1] * exp(relerr));
				} else {unknown_scope[0] *= exp(-3);antioss *=0.9f; }
			}else if (unknown_scope[2] == 0.0f){
				unknown_scope[0] = exp(unknown_scope[4] / unknown_scope[3] );
				unknown_scope[2] = 20.0f;
				relerr = 20.0f;
				antioss =1.0f;
				for(i=0;i<nbstates;i++) LL_out +=((ProbabilityDistribution<C>*)classes[i])->EMfinit();
			}else{
				// close approach
				double tmp;
				//	printf("Close Step:\n");
				relerr = log(unknown_scope[4]) - log(unknown_scope[1] * unknown_scope[3]);
				if ((ExCo<double>::isValid(unknown_scope[6]))&&(ExCo<double>::isValid(unknown_scope[7]))&&(unknown_scope[6] != 0.0f)&&(unknown_scope[7] != 0.0f)&&(fabs(1.0f / unknown_scope[6]) != 0.0f) &&(fabs(1.0f / unknown_scope[7]) != 0.0f) ){



					coef[0] = unknown_scope[4] - unknown_scope[1] * unknown_scope[3];
					coef[1] = unknown_scope[5];
					coef[2] = unknown_scope[6];
					coef[3] = unknown_scope[7];
					//		poly[2] = 2 * statbuf[3] + statbuf[2];
					//		poly[3] = 6 * statbuf[4] + 6 * statbuf[3] + statbuf[2];
					//		printf("%e\t%e\t%e\t%e\n",coef[0], coef[1], coef[2] , coef[3]  );

					tmp = -CubicRealRoot(coef,false);

					//		printf("%e poly eval\n", coef[0] -tmp *(coef[1] -tmp * (coef[2] -tmp * coef[3])));

					if ((!ExCo<double>::isValid(tmp))||(fabs(tmp) > 2.0f)){

						tmp = (unknown_scope[4] - unknown_scope[1] * unknown_scope[3]) / unknown_scope[5];
						//				printf("%e\t%e\n",(unknown_scope[4] - unknown_scope[1] * unknown_scope[3]) ,  unknown_scope[5]);
						if (!ExCo<double>::isValid(tmp)) tmp = 2.0f;

					}

				}else if (unknown_scope[3] > 0.0f){

					tmp = (unknown_scope[4] - unknown_scope[1] * unknown_scope[3]) / unknown_scope[5];
					//		printf("wrongwrong F[%e] = %e ! %e\n",tmp ,unknown_scope[3],unknown_scope[4]);
					//		printf("%e\t%e\n",(unknown_scope[4] - unknown_scope[1] * unknown_scope[3]) ,  unknown_scope[5]);
					if (!ExCo<double>::isValid(tmp)) tmp = 2.0f;


				} else {tmp = 2.0f; unknown_scope[0] *= exp(-3);antioss *=0.9f; relerr = 4.0f;}
				//printf("Near Step: F[%e] = %e !\n",unknown_scope[0], unknown_scope[4] / unknown_scope[3]);



				tmp *=antioss;
				if (fabs(tmp) >= 2.0f){
					unknown_scope[0] *=  (tmp < 0.0f) ? exp(2.0f) : exp(-2.0f);
					unknown_scope[2] *= 1.25f;
					if (unknown_scope[2] > 20.0f) unknown_scope[2] = 0.0f;
				}else{
					unknown_scope[2] *= fabs(tmp) /2.0f;
					unknown_scope[0] *= exp(-tmp);
				}
				//		printf("step = %e, new guess %e\n", tmp,unknown_scope[0]);

				if ((relerr < 1.0f)&&(ExCo<double>::isValid(relerr))){ // dont update if too much data is in unknown class
					for(i=0;i<nbstates;i++) LL_out += ((ProbabilityDistribution<C>*)classes[i])->EMfinit();
					LL_out /= (1.0f - unknown_scope[1] * exp(relerr));
				}else {unknown_scope[0] *= exp(-3);antioss *=0.9f; } // want to converge from below
			}

		}else for(i=0;i<nbstates;i++) LL_out += ((ProbabilityDistribution<C>*)classes[i])->EMfinit();
		return(LL_out);
	}

	template<class C, unsigned int nbstates>
	Classifier<C,nbstates> Classifier<C,nbstates>::EMfinit(double alpha) const{
		int i;
		Classifier<C,nbstates> f_out;
		if (unknown_scope != NULL){
			double coef[4];

			if (unknown_scope[2] > 0.125f){
				// far approach
				//			printf("Far Step:\n");
				//		printf("%e\t%e\t%e\t%e\n", (unknown_scope[4]-unknown_scope[5])/(unknown_scope[3]), (unknown_scope[6]-unknown_scope[7])/(unknown_scope[3]),(unknown_scope[6]+unknown_scope[7])/(unknown_scope[3]),(unknown_scope[4]+unknown_scope[5])/(unknown_scope[3]));

				coef[0] =  (-unknown_scope[4] + 9.0f * unknown_scope[6]) - 8.0f * unknown_scope[3] * unknown_scope[1];
				coef[1] =  ((unknown_scope[5]/-3.0f) + 9.0f * unknown_scope[7]) ;
				coef[2] =  (unknown_scope[4] - unknown_scope[6]);
				coef[3] =  ((unknown_scope[5]/3.0f) - unknown_scope[7]);
				//		printf("%e\t%e\t%e\t%e\n", coef[0], coef[1], coef[2], coef[3]);
				double shift = CubicRealRoot(coef, false);
				//printf("Far  Step: F[%e] = %e !\n",unknown_scope[0], unknown_scope[6]/(unknown_scope[3]));
				//		printf("%e\t%e\n", shift, unknown_scope[2]);



				shift *=antioss;
				if (!ExCo<double>::isValid(shift)){
					unknown_scope[0] *= exp((unknown_scope[6] > unknown_scope[3]* unknown_scope[1] ? -3.0f : 3.0f)* unknown_scope[2]);

				}else if (fabs(shift) > 3.0f){ // outside!
					unknown_scope[0] *= exp((shift < 0 ? -3.0f : 3.0f)* unknown_scope[2]);
					unknown_scope[2] *= 4.0f;

				}else{
					unknown_scope[0] *= exp(shift * unknown_scope[2]);
					unknown_scope[2] *= fabs(unknown_scope[1] -(unknown_scope[5]+ unknown_scope[6])/(unknown_scope[3]));
				}
				//		printf("end: %e\t    %e\n", unknown_scope[0], unknown_scope[2]);
				relerr = log(unknown_scope[6]) - log(unknown_scope[3] * unknown_scope[1]);

				if (relerr < 1.0f){ // dont update if too much data is in unknown class
					for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMfinit();

				} else {unknown_scope[0] *= exp(-3);antioss *=0.9f; }
			}else if (unknown_scope[2] == 0.0f){
				unknown_scope[0] = exp(unknown_scope[4] / unknown_scope[3] );
				unknown_scope[2] = 20.0f;
				relerr = 20.0f;
				antioss =1.0f;
				for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMfinit();
			}else{
				// close approach
				double tmp;
				//	printf("Close Step:\n");
				relerr = log(unknown_scope[4]) - log(unknown_scope[1] * unknown_scope[3]);
				if ((ExCo<double>::isValid(unknown_scope[6]))&&(ExCo<double>::isValid(unknown_scope[7]))&&(unknown_scope[6] != 0.0f)&&(unknown_scope[7] != 0.0f)&&(fabs(1.0f / unknown_scope[6]) != 0.0f) &&(fabs(1.0f / unknown_scope[7]) != 0.0f) ){



					coef[0] = unknown_scope[4] - unknown_scope[1] * unknown_scope[3];
					coef[1] = unknown_scope[5];
					coef[2] = unknown_scope[6];
					coef[3] = unknown_scope[7];
					//		poly[2] = 2 * statbuf[3] + statbuf[2];
					//		poly[3] = 6 * statbuf[4] + 6 * statbuf[3] + statbuf[2];
					//		printf("%e\t%e\t%e\t%e\n",coef[0], coef[1], coef[2] , coef[3]  );

					tmp = -CubicRealRoot(coef,false);

					//		printf("%e poly eval\n", coef[0] -tmp *(coef[1] -tmp * (coef[2] -tmp * coef[3])));

					if ((!ExCo<double>::isValid(tmp))||(fabs(tmp) > 2.0f)){

						tmp = (unknown_scope[4] - unknown_scope[1] * unknown_scope[3]) / unknown_scope[5];
						//				printf("%e\t%e\n",(unknown_scope[4] - unknown_scope[1] * unknown_scope[3]) ,  unknown_scope[5]);
						if (!ExCo<double>::isValid(tmp)) tmp = 2.0f;

					}

				}else if (unknown_scope[3] > 0.0f){

					tmp = (unknown_scope[4] - unknown_scope[1] * unknown_scope[3]) / unknown_scope[5];
					//		printf("wrongwrong F[%e] = %e ! %e\n",tmp ,unknown_scope[3],unknown_scope[4]);
					//		printf("%e\t%e\n",(unknown_scope[4] - unknown_scope[1] * unknown_scope[3]) ,  unknown_scope[5]);
					if (!ExCo<double>::isValid(tmp)) tmp = 2.0f;


				} else {tmp = 2.0f; unknown_scope[0] *= exp(-3);antioss *=0.9f; relerr = 4.0f;}
				//printf("Near Step: F[%e] = %e !\n",unknown_scope[0], unknown_scope[4] / unknown_scope[3]);



				tmp *=antioss;
				if (fabs(tmp) >= 2.0f){
					unknown_scope[0] *=  (tmp < 0.0f) ? exp(2.0f) : exp(-2.0f);
					unknown_scope[2] *= 1.25f;
					if (unknown_scope[2] > 20.0f) unknown_scope[2] = 0.0f;
				}else{
					unknown_scope[2] *= fabs(tmp) /2.0f;
					unknown_scope[0] *= exp(-tmp);
				}
				//		printf("step = %e, new guess %e\n", tmp,unknown_scope[0]);

				if ((relerr < 1.0f)&&(ExCo<double>::isValid(relerr))){ // dont update if too much data is in unknown class
					for(i=0;i<nbstates;i++) ((ProbabilityDistribution<C>*)classes[i])->EMfinit();

				}else {unknown_scope[0] *= exp(-3);antioss *=0.9f; } // want to converge from below
			}

		}else for(i=0;i<nbstates;i++) f_out.push_back(ExOp::new_class(classes[i]->EMfinit(alpha)));
		return f_out;
	}


template<class C, unsigned int nbstates> void  Classifier<C,nbstates>::performEM(Vector< C > &obs, int nbstep){
	// EM on independent observations
	Tuple<double, nbstates> likeli;
	int i , j, x;
	double LLhood;
	double LLhood_old;
	double alpha = 1.0f;

	Vector< ProbabilityDistribution<C>* > swaps_classes;

	for(i=0;i<nbstep;i++){

		this->EMinit();
		for(x = 0;x< obs.size();x++) {
			for(j=0;j< classes; j++) classes[j](likeli[j], obs[i]);
			likeli.normalize();
			this->EMregist(obs[i] , likeli);
		}

		LLhood = this->EMfinit();

		if (LLhood_old > LLhood) {i--;
			printf("Failed: Log-likelihood = %f\n", LLhood);
			alpha *= 0.1f;
			for(j=0;j< classes; j++) delete(classes[j]);
		}else{
			printf("Step %i: Log-likelihood = %f\n", LLhood);
			LLhood_old = LLhood;
			alpha = 1.0f;
			if (swaps_classes.size() == nbstates){
			for(j=0;j< classes; j++) {delete(swaps_classes[j]); swaps_classes[j] =classes [j];}
			}else for(j=0;j< classes; j++) swaps_classes.push_back(classes [j]);
		}
		for(j=0;j< classes; j++) classes[j] = swaps_classes[j]->EMnext(alpha);
	}

}

#undef LFHTEMP
#define LFHTEMP template<unsigned int DIMS>


LFHTEMP EuclidianEllipse<DIMS>::EuclidianEllipse():fixvars(0), updateGreedy(false) {	alpha_search.init(1.0f, DIMS*2);}
LFHTEMP EuclidianEllipse<DIMS>::EuclidianEllipse(const Tuple<double, 2> &in_center): center(in_center), fixvars(0), updateGreedy(false){ alpha_search.init(1.0f, DIMS*2); ExOp::toZero(eccent);}
// LFHTEMP EuclidianEllipse<DIMS>::EuclidianEllipse(const EuclidianEllipse<DIMS>&o) : center(o.center),eccent(o.eccent),updateGreedy(o.updateGreedy), fixvars(o.fixvars),width(o.width),nihvar(o.nihvar),fact(o.fact),EMscope(o.EMscope),alpha_search(o.alpha_search){}



LFHTEMP void EuclidianEllipse<DIMS>::operator()(double & f_out , pair< Tuple<double, DIMS> , double> &instance) const{
	double error = ExOp::norm( (instance.first) - center + eccent) + ExOp::norm( (instance.first) - center - eccent) + instance.second - width;


	// error is the expected value for
	//	double lamba = pow(-0.5, 0.5 * nihvar); // std
	//		lamba = 1.0f / ((error < lamba) ? lamba : error);


	//	f_out = inside_factor + (1 - inside_factor) * erf( nihvar * error ); // probability for the
	//	f_out *= inside_factor + lamba *exp(lamba * error);

	f_out = fact *exp( error*nihvar*error);

}

LFHTEMP double EuclidianEllipse<DIMS>::LL(const pair< Tuple<double, DIMS> , double> &instance) const{
	double error = ExOp::norm( (instance.first) - center + eccent) + ExOp::norm( (instance.first) - center - eccent) + instance.second - width;
	return (log(fact) + error*nihvar*error);
}

LFHTEMP	void EuclidianEllipse<DIMS>::EMinit(){ // uses local gradient!
	ExOp::toZero(EMscope);
}
LFHTEMP	void EuclidianEllipse<DIMS>::EMAlphainit(double alpha){
	//		if (!EMscope) {
	//			EMscope = new pair<INPUT, Tuple<double,4> >();
	//			ExOp::zero(EMscope->first);
	//			ExOp::zero(EMscope->second);
	//		}else{
	//			EMscope->first *= alpha;
	//			EMscope->second *= alpha;
	//		}
}

LFHTEMP	void EuclidianEllipse<DIMS>::EMpreregist(unsigned int step, const pair< Tuple<double, DIMS> ,double> &instance, const double prob){
	// find center of mass!
	EMscope[0] += prob; // sw
	int i;
	for(i=0;i<DIMS;i++){
		EMscope[3 +i ] -= prob * instance.first[i];
	}

}

LFHTEMP void EuclidianEllipse<DIMS>::EMprefinit(unsigned int step){
	int i;
	for(i=0;i<DIMS;i++){
		center[i] = EMscope[3 +i] / EMscope[0];
		EMscope[3 +i] =0.0f;
	}
	EMscope[0] = 0.0f;
}


LFHTEMP	void EuclidianEllipse<DIMS>::EMregist(const pair< Tuple<double, DIMS> ,double> &instance, const double prob){


	Tuple<double, DIMS> tmppt = (instance.first) - center + eccent;
	double norma = tmppt.norm();
	tmppt = (instance.first) - center - eccent;
	double normb = tmppt.norm();

	double errorval = norma + normb + instance.second  ; // deviation!

	EMscope[2] += prob * errorval * errorval; // swx2
	errorval *= prob;
	EMscope[1] += errorval; // swx
	EMscope[0] += prob; // sw
	// sufficien stat for radius update, conditionn

	// likelifunc += sum prob * errorval /  )

	int i;
	for(i=0;i<DIMS;i++){
		EMscope[3 +i ] -= errorval * (((instance.first[i] - center[i] + eccent[i]) / norma) + ((instance.first[i] - center[i] - eccent[i]) / normb));
		EMscope[3 +i + DIMS] += errorval * (((instance.first[i] - center[i] + eccent[i]) / norma) - ((instance.first[i] - center[i] - eccent[i]) / normb));
		EMscope[3 +i + DIMS*2] -= prob * (((instance.first[i] - center[i] + eccent[i]) / norma) + ((instance.first[i] - center[i] - eccent[i]) / normb));
		EMscope[3 +i + DIMS*3] += prob * (((instance.first[i] - center[i] + eccent[i]) / norma) - ((instance.first[i] - center[i] - eccent[i]) / normb));
	}

}

LFHTEMP double EuclidianEllipse<DIMS>::EMfinit(){

	if (EMscope[0] == 0.0f) return(0.0f); // did not get any points, dont update

	if (updateGreedy){
		return 0.0f;



	}else{
	// updates radius and nihval only!
	width = EMscope[1] / EMscope[0];
	if (!this->has_fix_stddev()){
	nihvar = -0.5f * EMscope[0] / (EMscope[2] - width *  EMscope[1]);
	fact = sqrt(nihvar / -M_PI);
	}
	double derivative[DIMS*2];

	// cached next calculation
	int i;
	for(i=0;i<DIMS;i++){
		EMscope[3 + i] =        derivative[i]        = nihvar * (EMscope[3 + i]       - width *  EMscope[3 + i +DIMS*2]) / EMscope[0] ;
		EMscope[3 + i +DIMS]  = derivative[i + DIMS] = nihvar * (EMscope[3 + i +DIMS] - width *  EMscope[3 + i +DIMS*3]) / EMscope[0] ;
	}

	alpha_search.registAscent( (log(fact) - 0.5f) ,derivative);
	EMscope[2] = alpha_search(); // stores alpha
	//printf("wnabe alpha! %f\n", EMscope[2]);

	switch (fixvars & 3){
		case 0: return EMscope[0] * (log(fact) - 0.5f);
		case 2: return EMscope[0] * log(fact) - 0.5f * nihvar * (EMscope[2] - EMscope[1] * width); // log likelihood!
		default: return EMscope[0] * log(fact) - 0.5f * nihvar * (EMscope[2] - 2*EMscope[1] * width + EMscope[0] * width * width); // log likelihood!
		}
	}
}


LFHTEMP	double EuclidianEllipse<DIMS>::EMloglikelihood() const{
	switch (fixvars & 3){
		case 0: return EMscope[0] * (log(fact) - 0.5f);
		case 2: return EMscope[0] * log(fact) - 0.5f * nihvar * (EMscope[2] - EMscope[1] * width); // log likelihood!
		default: return EMscope[0] * log(fact) - 0.5f * nihvar * (EMscope[2] - 2*EMscope[1] * width + EMscope[0] * width * width); // log likelihood!
	}
}




LFHTEMP	EuclidianEllipse<DIMS> EuclidianEllipse<DIMS>::EMnext_exec(double alpha) const{
	// updates the focal points! alpha updates!

	EuclidianEllipse<DIMS> f_out;
	f_out.width = width;
	f_out.nihvar = nihvar;
	f_out.fact = fact;
	f_out.alpha_search = alpha_search;
	int i;
	for(i=0;i<DIMS;i++){
		f_out.center[i] = center[i] + alpha * (EMscope[2] * EMscope[3 + i] + ((0.01f * width * rand()) / RAND_MAX));
		f_out.eccent[i] = eccent[i] + alpha * (EMscope[2] * EMscope[3 + i + DIMS] + ((0.01f * width * rand()) / RAND_MAX)) ;
	}

	return f_out;
}

LFHTEMP	void EuclidianEllipse<DIMS>::EMclear(){
	//		if (EMscope) {delete(EMscope); EMscope = NULL;}
}





LFHTEMP	void EuclidianEllipse<DIMS>::show(FILE* o) const{
	fprintf(o, "DistanceError Distribution: rad=%f +- %f\tcenter :[", width, pow(-2.0f * nihvar, -0.5f ));
	ExOp::show(center,o,2);
	fprintf(o, "]+-[");
	ExOp::show(eccent,o,2);
	fprintf(o, "]\n");
}





template<class C, unsigned int DIMS> DataGrid<double,2> MetaPixels::performEM(unsigned int nbstep, const DataGrid<C, DIMS> &obs, ClassifierV<C> &model, const DataGrid<unsigned int, DIMS> &meta, unsigned int nbmetapix, unsigned int first_val){
	Vector< C > *pixlists = new Vector< C >[nbmetapix];
	unsigned int i,j;

	typename DataGrid<C, DIMS>::KeyIterator ite = obs.getKeyIterator();

	if (ite.first()) do{
		j = meta(ite());
		if ((j >= first_val)&&(j < first_val + nbmetapix)) pixlists[j -  first_val].push_back(obs(ite()));
	} while ( ite.next());


	unsigned int nbstates = model.nbstates();


	double *likeli = new double[nbstates];
	double* prior_freq = new double[nbstates];
	double* prior_freq_learn = new double[nbstates];

	for(j=0;j< nbstates; j++) prior_freq[j] = 1.0f / nbstates;



	int r;
	unsigned int x,y;

	unsigned int totnbpix =0;for(y = 0;y< nbmetapix;y++) totnbpix += pixlists[y].size();

	double LLhood;
	double LLhood_old;

	double alpha = 1.0f;
	double sum,saf;


	bool rand_init = false;
	for(i=0;i<nbstep;i++){

		model.EMinit();
		memset(likeli,'\0',sizeof(double)*nbstates);
		memset(prior_freq_learn,'\0',sizeof(double)*nbstates);
		if ((i==0)&&(rand_init)){
			for(y = 0;y< nbmetapix;y++) {
				r = rand() % nbstates;
				for(j=0;j< nbstates; j++) likeli[j] = (r == j) ? 1.0f : 0.0f;
			for(x = 0;x< pixlists[y].size();x++) model.EMregist(pixlists[y][x] , likeli);
			}
		}else{
			for(y = 0;y< nbmetapix;y++) {
				memset(likeli,'\0',sizeof(double)*nbstates);
				for(x = 0;x< pixlists[y].size();x++) for(j=0;j< nbstates; j++) likeli[j] += model[j]->LL(pixlists[y][x]);
				saf = likeli[0];
				for(j=1;j< nbstates; j++) if (saf < likeli[j]) saf = likeli[j];
				sum =0;
				for(j=0;j< nbstates; j++) {likeli[j] = prior_freq[j] * exp(likeli[j] - saf); sum += likeli[j];}
				for(j=0;j< nbstates; j++) likeli[j] /= sum;
				for(x = 0;x< pixlists[y].size();x++) model.EMregist(pixlists[y][x] , likeli);
				for(j=0;j< nbstates; j++) prior_freq_learn[j] += likeli[j] * pixlists[y].size();
			}
		}
		LLhood = model.EMfinit();
		for(j=0;j< nbstates; j++) prior_freq[j] = (ExCo<double>::epsilon() + prior_freq_learn[j]) / ( (ExCo<double>::epsilon() *nbstates) + totnbpix);

		if (((LLhood_old > LLhood)||(!ExCo<double>::isValid(LLhood)))&&(i != 0)) {
			printf("Failed: Log-likelihood = %e\n", LLhood);
			alpha *= 0.1f;
			model.EMaccept(false);
		}else{
			printf("Step %i: Log-likelihood = %e\n", i , LLhood);
			LLhood_old = LLhood;
			alpha = 1.0f;
			model.EMaccept();

		}
		model.EMnext(alpha);
	}
	model.EMclear();

	DataGrid<double,2> fout;
	unsigned int coor[2];
	coor[0] = nbstates;
	coor[1] = nbmetapix;
	fout.setSizes(coor);


	for(coor[1] = 0;coor[1]< nbmetapix;coor[1]++) {
		memset(likeli,'\0',sizeof(double)*nbstates);
		for(x = 0;x< pixlists[coor[1] ].size();x++) for(j=0;j< nbstates; j++) likeli[j] += model[j]->LL(pixlists[coor[1]][x]);
		saf = likeli[0];
		for(j=1;j< nbstates; j++) if (saf < likeli[j]) saf = likeli[j];
		sum =0;
		for(j=0;j< nbstates; j++) {likeli[j] = prior_freq[j] * exp(likeli[j] - saf); sum += likeli[j];}
		for(coor[0]=0;coor[0]< nbstates; coor[0]++) fout(coor) = likeli[ coor[0] ] / sum;
	}



	delete[](likeli);
	delete[](prior_freq);
	delete[](prior_freq_learn);





	delete[](pixlists);
	return fout;
}

#undef LFHTEMP
#define LFHTEMP template<class C, unsigned int SIZE>


LFHTEMP typename ExCo<C>::LMUL_TYPE GaussianObj<C,SIZE>::computeProj(const Tuple<C,SIZE> &) const{

}

LFHTEMP void GaussianObj<C,SIZE>::operator()(double &fout , Tuple<C,SIZE> &fin) const{
    typename ExCo<C>::LMUL_TYPE mag = computeProj(fin);
    fout = factor * exp(-ExOp::norm(mag));
}

LFHTEMP double GaussianObj<C,SIZE>::LL(const C& fin) const{
    typename ExCo<C>::LMUL_TYPE mag = computeProj(fin);
    return log(factor) -ExOp::norm(mag);
}

LFHTEMP void GaussianObj<C,SIZE>::EMinit(){
 //   delete[](dascope); dascope = new GOscope();
}
LFHTEMP void GaussianObj<C,SIZE>::EMAlphainit(double a){
 //   delete[](dascope); dascope = new GOscope();
}

LFHTEMP double GaussianObj<C,SIZE>::EMfinit(){
	mean = scope.getMean();
	ihvar = scope.cov.inverse() * -0.5f;
    return(0.0f);
}

LFHTEMP void GaussianObj<C,SIZE>::EMregist(const Tuple<C,SIZE> &instance, const double prob){
	scope += GaussElem<C,SIZE>(instance, prob);
}

LFHTEMP GaussianObj<C,SIZE>* GaussianObj<C,SIZE>::EMnext(double alpha) const{
    return new GaussianObj<C,SIZE>(*this);
}

LFHTEMP void GaussianObj<C,SIZE>::EMclear(){
//delete[](dascope);dascope = NULL;
}

LFHTEMP void GaussianObj<C,SIZE>::show(FILE *f, unsigned int level)const{
 //   fprintf(f,"Gaussian Object: mean = "); ExOp::show(mean, f, level+1);
 //   fprintf(f,"\tProject = ");ExOp::show(proj, f, level+1);
 //   fprintf(f,"\n");
}


#undef LFHTEMP
#define LFHTEMP template<unsigned int NBPARAM, unsigned int NBINPUT>

	
LFHTEMP	void GaussianProcessProbSimple<NBPARAM,NBINPUT>::EMinit(){

	noise.EMinit();
	EM_obs.clear();
    

/*	
	if (resample_to_fixgrid){

		Tuple<double, NBINPUT> daz; ExOp::toZero(daz);
		for(unsigned int i=0;i<registered.size();i++){registered[i].first.d = daz; registered[i].second = 0.01f;}
	}else{
		registered.clear();
	}
*/

}
LFHTEMP	void GaussianProcessProbSimple<NBPARAM,NBINPUT>::EMAlphainit(double a){
	noise.EMAlphainit(a);
    noise.EMinit();
	EM_obs.clear();
}


LFHTEMP	void GaussianProcessProbSimple<NBPARAM,NBINPUT>::EMregist(const KeyElem< Tuple<double, NBPARAM> , Tuple<double, NBINPUT> > &data,  const double prob){
    if (resample_to_fixgrid){
    }else EM_obs.push_back(KeyElem< double, KeyElem< Tuple<double, NBPARAM> , Tuple<double, NBINPUT> > >(prob, data));
	
	noise.EMregist(data.d, prob);
	
}


LFHTEMP	void GaussianProcessProbSimple<NBPARAM,NBINPUT>::EMclear(){
	EM_obs.clear();
}


LFHTEMP double GaussianProcessProbSimple<NBPARAM,NBINPUT>::setDerivative_routine(double* i_deriv, bool showmat){
    Trianglix<double, 0u> damatrix;
    Trianglix<double, 0u> damatrix_nonoise;
   
	Trianglix<double, 0u> daderivative;
	Trianglix<double, 0u> dainverse;
    Trianglix<double, 0u> dainner;

	inner_param.setSize(EM_obs.size());
	
	Tuple<double,0u> datuple;
Matrix<double> test;


    double ttmptmp,tmp_dev;
    
 //   DataGrid<double, 2> da_vectors;
   // DataGrid<double, 1> inner_prods;
 //   unsigned int coor[2];
 //   coor[0] = 1;
    unsigned int s = EM_obs.size();
 //   coor[1] = s;
 //   da_vectors.setSizes(coor);
    unsigned int i,j,k;
 //   printf("%e\t%e\t%e\n",scope[0],scope[1],scope[2]);
				damatrix_nonoise.setSize(s);damatrix.setSize(s);datuple.setSize(s);
				daderivative.setSize(s);dainverse.setSize(s);dainner.setSize(s);
					for(i=0,k=0;i<s;i++){
               //         coor[1] = i;
               //         for(coor[0] = 0 ; coor[0] < 1 ; coor[0]++) da_vectors(coor) = EM_obs[i].d.d[coor[0]];
						datuple[i] = EM_obs[i].d.d[0];
						inner_param[i].k = EM_obs[i].d.k;
				//		printf("%e\t%e\n",EM_obs[i].d.k[0],EM_obs[i].d.d[0]);
						for(j=0;j<i;j++,k++){
							double d = ExOp::pnorm(EM_obs[i].d.k[0] -EM_obs[j].d.k[0]); 
							//damatrix.data[k] = exp(scope[1] - exp(scope[0]) * d*d);
							//daderivative.data[k] = damatrix.data[k] * (-exp(scope[0]) * d*d);
							damatrix_nonoise.data[k] =  exp(scope[1] - exp(scope[0]) * d);
							daderivative.data[k] = damatrix_nonoise.data[k] * (-d);
						}
						damatrix_nonoise.data[k] = exp(scope[1]);
						daderivative.data[k] = 0.0f;
						k++;
					}
					damatrix = damatrix_nonoise;
             //       for(i=0,j=0;j<s;j++) {damatrix.data[i] += exp(scope[2]+ scope[1]);i += j+2;}
                    for(i=0,j=0;j<s;j++) {damatrix.data[i] += exp(scope[2]);i += j+2;}
	if (showmat)	{Matrix<double>(damatrix).show(); printf("%e\t%e\n",scope[0],exp(scope[0]));}
	
		dainverse = damatrix.mkinverse();
	if (!ExOp::isValid(dainverse.data[0])) dainverse = damatrix.inverse_MK2();
	
	
	
	//Matrix<double>(damatrix_nonoise).show();
	
		
	//	(Matrix<double>(damatrix) * Matrix<double>(dainverse)).show();
                //    if (fabs(test.data[6] - 1.0f) > 0.001f) {damatrix.show(); printf("OMG!\n"); }
                    //test.show();

				//	datuple.show();
                 //   dainverse.show();
				//	(dainverse * datuple).show();
//	ExOp::show(datuple);
					dainner = dainverse.Xformed_outer_product(datuple);
                 //   printf("davec: ");
                 //   ExOp::show(dainverse * datuple );

                //    printf("damat: ");dainner.show();
					dainner -= dainverse;
					
				/*	dainner *= exp(-scope[1]);
					dainverse *= exp(-scope[1]);
					damatrix *= exp(scope[1]);
					damatrix_nonoise *= exp(scope[1]);
					daderivative *= exp(scope[1]);*/
				//	ExOp::show( dainverse.Xformed_inner_product(datuple) );
				//	ExOp::show(dainner.trace()); 
				//	ExOp::show(dainner.trace_of_product(damatrix));
					i_deriv[2] = 0.5f * exp(scope[2]) * dainner.trace(); // +scope[1] in exp
                    i_deriv[1] = 0.5f * dainner.trace_of_product(damatrix_nonoise) ; // + i_deriv[2]
					i_deriv[0] = 0.5f * exp(scope[0]) * dainner.trace_of_product(daderivative);
                //    printf("%e\t%e\t%e    .. %e   %e \n", deriv[0], deriv[1], deriv[2],dainverse.Xformed_inner_product(datuple),damatrix.log_determinant());
					
              //   ttmptmp = -0.5f * (dainverse.Xformed_inner_product(datuple) + damatrix.log_determinant());    
               
               	  // ttmptmp = dainverse.Xformed_inner_product(datuple);
               	  ttmptmp = damatrix.Xformed_inner_product_of_inverse(datuple);
	
	
               	  if (ttmptmp < 0.0f) ExOp::toMin(ttmptmp);
               	  else ttmptmp = -0.5f * (damatrix.log_determinant() + ttmptmp); 
	
	       //         printf("damat: %e\n", damatrix.log_determinant());
              // 	  ttmptmp = damatrix.log_determinant();
              //  if (ttmptmp > 500) printf("%e\t%e\n",dainverse.Xformed_inner_product(datuple), damatrix.log_determinant() ); //	(Matrix<double>(damatrix) * Matrix<double>(dainverse)).show();
	for(i=0;i<NBPARAM*2+1;i++){
		tmp_dev = (scope[i] -  param_prior_mean[i]) / param_prior_std[i];
		ttmptmp -= 0.5f * tmp_dev * tmp_dev; i_deriv[i] -= tmp_dev;
	}
	
	Tuple<double,0u> datuple2 = dainverse * datuple;
	
	for(i=0;i<s;i++) inner_param[i].d[0] = datuple2[i];
	
		return ttmptmp;
    }

LFHTEMP	double GaussianProcessProbSimple<NBPARAM,NBINPUT>::EMfinit(){
	
	
//	if (EMstep_count == 0){
		
//		llikelyhood = noise.EMfinit();
//		scope[1] = log(noise.getVar(0));
//		scope[2] = log(noise.getVar(0));
//	}
//		noise.EMfinit();
	// WHERE THE FUN HAPPENS !!!


        double eps = 0.000001f;
        
        llikelyhood = setDerivative_routine(deriv);

//	KeyElem< Tuple<double, NBPARAM> , Tuple<double, NBINPUT> > * inner_param;
	

  //     double fake[3];double tmp;
  //      for(unsigned int i =0;i<3;i++){ tmp = scope[i]; scope[i] += eps; printf("%i: Rel Err = %e  (%e)\n",i, 1.0f  - ((setDerivative_routine(fake)- llikelyhood) / (deriv[i] *eps) ), deriv[i] ); scope[i] = tmp;}


       


/*
					
					dainverse = damatrix.inverse();
					
				//	(Matrix<double>(damatrix) * Matrix<double>(dainverse) ).show();
					datuple.show();
					(dainverse * datuple).show();
					dainner = dainverse.Xformed_outer_product(datuple);
					dainner -= dainverse;
					
					dainner.show();	printf("%f\n", dainner.trace());
					
					deriv[0] = -0.5f * dainner.trace();
					deriv[1] = -0.5f * (dainner.trace_of_product(damatrix) - dainner.trace());
					deriv[2] = -0.5f * dainner.trace_of_product(daderivative);
					value = -0.5f * (dainverse.Xformed_inner_product(datuple) + damatrix.log_determinant());

					printf("deriv = [%e,%e,%e]\n", deriv[0]*0.000001f,deriv[1]*0.000001f,deriv[2]*0.000001f);
					printf("F[%e,%e,%e] = %e\n", guess[0],guess[1],guess[2], value);
*/

    return llikelyhood;
}

LFHTEMP	void GaussianProcessProbSimple<NBPARAM,NBINPUT>::EMnext(double alpha, GaussianProcessProbSimple<NBPARAM, NBINPUT> &fout) const{
     fout.alpha_search = alpha_search;
     memcpy(fout.scope, scope, sizeof(scope) );
     fout.alpha_search.updateAscent(llikelyhood, fout.scope, deriv);
}
LFHTEMP	void GaussianProcessProbSimple<NBPARAM,NBINPUT>::EMstep(double alpha){

  //   printf("step %e\n", alpha_search.updateAscent(llikelyhood, scope, deriv));
	
	printf("bef: %e\t%e\t%e\n", scope[0], scope[1], scope[2]);
    printf("%e \t", alpha_search.updateAscent(llikelyhood, scope, deriv));
	printf("aft: %e\t%e\t%e\n", scope[0], scope[1], scope[2]);
}
	
	LFHTEMP	Tuple<double,NBINPUT> GaussianProcessProbSimple<NBPARAM,NBINPUT>::getSignalMean(const Tuple<double, NBPARAM>& whe) const{ Tuple<double,NBINPUT> fout;
		// needs to compute distance to all key points
		unsigned int i,j,k,l;
		Tuple<double, NBPARAM> dif;
		for(i=0 ; i < inner_param.getSize();i++){
			dif = whe - inner_param[i].k;
			// multiply my upper triangular and computes self-innerproduct:
			if (i == 0) fout = inner_param[0].d * exp(scope[1] - exp(scope[0]) * ExOp::pnorm(dif) );
			else fout +=  inner_param[i].d * exp(scope[1] - exp(scope[0]) * ExOp::pnorm(dif) );
		}
		return(fout);
	}
	

LFHTEMP	void GaussianProcessProbSimple<NBPARAM,NBINPUT>::setGrid(const Vector< Tuple<double, NBPARAM> > &ingr){

	if (ingr.size() == 0) {static_warning_handdle << LFH_WARNING_UNEXPECTED_INPUT;return;}
	resample_to_fixgrid = true;
	nb_param = ingr.size();
	inner_param.setSize(nb_param);
	KeyElem< Tuple<double,NBPARAM>,  Tuple<double,NBINPUT> > a;
	ExOp::toZero(a);

	for(int i=0;i<ingr.size();i++) {a.k = ingr[i]; inner_param[i] = a;}

	}

LFHTEMP	void GaussianProcessProbSimple<NBPARAM,NBINPUT>::show(FILE* f, int level)const{
    
}
