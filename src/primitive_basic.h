/*
 * primitive_basic.h
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

// contains KeyElem, Tuple, TMatrix, matriangle..., Complex, Quaternion, WeightElem

	template<class KEY, class DATA,bool isMax> ExtremumScope<KEY,DATA,isMax>& ExtremumScope<KEY,DATA,isMax>::toZero(){ExOp::toMin(best); return(*this);}
template<class KEY, class DATA,bool isMax> void ExtremumScope<KEY,DATA,isMax>::init(const KEY& key, const DATA& data){if (ExOp::isValid(data)) {best  =data; best_key = key;}else ExOp::toMin(best);}
template<class KEY, class DATA,bool isMax> void ExtremumScope<KEY,DATA,isMax>::regist(const KEY& key, const DATA& data){if ((ExOp::isValid(data))&&(best < data)) {best  =data; best_key = key;}}
template<class KEY, class DATA> ExtremumScope<KEY,DATA,false>& ExtremumScope<KEY,DATA,false>::toZero(){ExOp::toMax(best);return(*this);}
template<class KEY, class DATA> void ExtremumScope<KEY,DATA,false>::init(const KEY& key, const DATA& data){if (ExOp::isValid(data)) {best  =data; best_key = key;}else ExOp::toMax(best);}
template<class KEY, class DATA> void ExtremumScope<KEY,DATA,false>::regist(const KEY& key, const DATA& data){if ((ExOp::isValid(data))&&(best > data)) {best  =data; best_key = key;}}

	// class KeyElem
#undef LFHTEMP
#define LFHTEMP template <class C, class B>
	
LFHTEMP KeyElem<C,B>::KeyElem(){}
LFHTEMP KeyElem<C,B>::KeyElem(const C& _k, const B& _n): k(_k),d(_n){}
LFHTEMP bool KeyElem<C,B>::isValid() const{return ( ExOp::isValid(k) && ExOp::isValid(d) );}
	
#undef LFHTEMP
#define LFHTEMP template <class C,unsigned int TSIZE,Tuple_flag Cflag>
LFHTEMP double Tuple<C,TSIZE,Cflag>::weight(){return(data[TSIZE]);}


LFHTEMP void Tuple<C,TSIZE,Cflag>::zero(){for(unsigned int i=0;i<TSIZE;i++) ExOp::toZero(data[i]);}
LFHTEMP void Tuple<C,TSIZE,Cflag>::random(){for(unsigned int i=0;i<TSIZE;i++) ExOp::toRand(data[i]);}

LFHTEMP template<class D> Tuple<C,TSIZE,Cflag>::operator Tuple<D,TSIZE> () const{
		Tuple<D,TSIZE> out;
		for(unsigned int i=0;i<TSIZE;i++) out.data[i] = (D)data[i];
		return(out);
		}

LFHTEMP double Tuple<C,TSIZE,Cflag>::pnorm() const{
		double sum = ExOp::pnorm(data[0]);
		for(unsigned int i =1;i<TSIZE;i++) sum += ExOp::pnorm(data[i]);
		return(sum);
	}

LFHTEMP double Tuple<C,TSIZE,Cflag>::norm() const{
		double sum = ExOp::pnorm(data[0]);
		for(unsigned int i =1;i<TSIZE;i++) sum += ExOp::pnorm(data[i]);
		return(sqrt(sum));
	}

LFHTEMP void Tuple<C,TSIZE,Cflag>::save( FILE *f) const {
	if (ExCo<C>::IsPOD) fwrite(data,sizeof(C),TSIZE,f);
	else for(unsigned int i=0;i<TSIZE;i++) ExOp::save(data[i],f);

}
LFHTEMP void Tuple<C,TSIZE,Cflag>::load( FILE *f, unsigned int ch_TSIZE) {
	if (ExCo<C>::IsPOD) fread(data,sizeof(C),TSIZE,f);
	else for(unsigned int i=0;i<TSIZE;i++) ExOp::load(data[i],f, ch_TSIZE/TSIZE);
}

LFHTEMP bool Tuple<C,TSIZE,Cflag>::isValid() const{
unsigned int i;
for(i=0;i<TSIZE;i++) if (!ExOp::isValid(data[i])) break;
return(i == TSIZE);
}


LFHTEMP void Tuple<C,TSIZE,Cflag>::show( FILE *f_out, int l) const{
	switch(l){
		case 1:
		case 0:
			for(unsigned int i=0;i<TSIZE;i++) {
				if (i != 0) fprintf(f_out,"\t");
				ExOp::show(data[i],f_out,2);
			}
			if (l==0) fprintf(f_out,"\n");
		break;
		case 2:
			fprintf(f_out,"[");
			for(unsigned int i=0;i<TSIZE;i++) {
				if (i != 0) fprintf(f_out,";");
				ExOp::show(data[i],f_out,3);
			}
			fprintf(f_out,"]");
			break;
		default:
			fprintf(f_out,"(");
			for(unsigned int i=0;i<TSIZE;i++) {
				if (i != 0) fprintf(f_out,",");
				ExOp::show(data[i],f_out,4);
			}
			fprintf(f_out,")");
			break;
	}
}

LFHTEMP string Tuple<C,TSIZE,Cflag>::type_tostring() const{char buffer[256]; sprintf(buffer,"%u",TSIZE); return string("Tuple<") + ExOp::type_tostring(data[0]) +string(",")+ ExOp::type_tostring(buffer) + string(">");}

LFHTEMP GaussScope< Tuple<C,TSIZE> , typename MT_IFTYPE<ExCo<C>::IS_COMMUTATIVE::ans, Tuple<C, TEMPLATE_TRIANGLE_NUMBER<TSIZE, 2>::ans> , TMatrix<C,TSIZE,TSIZE> >::TYPE >  Tuple<C,TSIZE,Cflag>::mkgaussstat(double &w) const{
    typename MT_IFTYPE<ExCo<C>::IS_COMMUTATIVE::ans, Tuple<C, TEMPLATE_TRIANGLE_NUMBER<TSIZE, 2>::ans> , TMatrix<C,TSIZE,TSIZE> >::TYPE outer;
    unsigned int i,j,k;
    if (ExCo<C>::IS_COMMUTATIVE::ans){
        for(i=0,k=0; i< TSIZE;i++) for(j=i;j<TSIZE;j++,k++) outer[k] = data[i] * data[j];
    }else{
        for(i=0; i< TSIZE;i++) for(j=0;j<TSIZE;j++) outer[i + j * TSIZE] = data[i] * data[j];
    }
    return GaussScope< Tuple<C,TSIZE> , typename MT_IFTYPE<ExCo<C>::IS_COMMUTATIVE::ans, Tuple<C, TEMPLATE_TRIANGLE_NUMBER<TSIZE, 2>::ans> , TMatrix<C,TSIZE,TSIZE> >::TYPE >(*this, outer ,w);
    }

LFHTEMP static double overlap(const GaussScope< Tuple<C,TSIZE> , typename MT_IFTYPE<ExCo<C>::IS_COMMUTATIVE::ans, Tuple<C, TEMPLATE_TRIANGLE_NUMBER<TSIZE, 2>::ans> , TMatrix<C,TSIZE,TSIZE> >::TYPE >&a, const GaussScope< Tuple<C,TSIZE> , typename MT_IFTYPE<ExCo<C>::IS_COMMUTATIVE::ans, Tuple<C, TEMPLATE_TRIANGLE_NUMBER<TSIZE, 2>::ans> , TMatrix<C,TSIZE,TSIZE> >::TYPE >&b){

    return 0;
}


//LFHTEMP Tuple<C,TSIZE,Cflag>::Tuple(C const & other){	int i;	for(i = TSIZE-1;i>=0;i--) data[i] = other;}

//LFHTEMP
//Tuple<C,TSIZE,Cflag>::Tuple(Tuple<C,TSIZE,Cflag> const & clonefrom){memcpy(data,clonefrom.data, TSIZEof(C)*(TSIZE));}

LFHTEMP C& Tuple<C,TSIZE,Cflag>::operator[](int const pos) {return(data[pos]);}
LFHTEMP const C& Tuple<C,TSIZE,Cflag>::operator[](int const pos) const {return(data[pos]);}

LFHTEMP template<class O,unsigned int oTSIZE, Tuple_flag Oflag>
char Tuple<C,TSIZE,Cflag>::compare(const Tuple<O,oTSIZE, Oflag>& other) const{
	unsigned int i=0;
	unsigned int m = (TSIZE < oTSIZE) ? TSIZE : oTSIZE;
	for(;i<m;i++) if ((*this)[i] != other[i]) break;
	if (i == m) return( (TSIZE == oTSIZE) ? 0 : ((TSIZE < oTSIZE)  ? 2 : 1) );
	else return( ((*this)[i] < other[i]) ? 2 : 1 );
}

// Super Operators!
LFHTEMP
template<class I> void Tuple<C,TSIZE,Cflag>::operator() (Oper1<I> const & op){ // not a match
for(unsigned int i = 0;i < TSIZE;i++) data[i](op);
}
LFHTEMP 
void Tuple<C,TSIZE,Cflag>::operator() (Oper1< C> const & op){ // match
for(unsigned int i = 0;i < TSIZE;i++)  op(data[i], data[i]);
}
LFHTEMP
template<class A_1, class A_2, class C_2, unsigned int TSIZE_2> void Tuple<C,TSIZE,Cflag>::operator() (Oper2<A_1,A_2> const & op, Tuple<C_2, TSIZE_2,Cflag> const & a_2 ){ // not a match
unsigned int i = (TSIZE > TSIZE_2) ? TSIZE_2 : TSIZE;
for(i--;i!=ExCo<unsigned int>::maximum();i--) data[i](op,a_2.data[i]);
}

LFHTEMP
template<class C_2, unsigned int TSIZE_2> void Tuple<C,TSIZE,Cflag>::operator() (Oper2<C,C_2> const & op, Tuple<C_2, TSIZE_2,Cflag> const & a_2){ // match
unsigned int i = (TSIZE > TSIZE_2) ? TSIZE_2 : TSIZE;
for(i--;i!=ExCo<unsigned int>::maximum();i--) op(data[i], a_2.data[i]);
}
LFHTEMP
template<class A_1, class A_2, class C_2, unsigned int TSIZE_2> void Tuple<C,TSIZE,Cflag>::operator() (Oper2<A_1,A_2> const & op, Tuple<C_2, TSIZE_2,Cflag> & a_2 ){ // not a match
unsigned int i = (TSIZE > TSIZE_2) ? TSIZE_2 : TSIZE;
for(i--;i!=ExCo<unsigned int>::maximum();i--) data[i](op,a_2.data[i]);
}
LFHTEMP
template<class C_2, unsigned int TSIZE_2> void Tuple<C,TSIZE,Cflag>::operator() (Oper2<C,C_2> const & op, Tuple<C_2, TSIZE_2,Cflag> & a_2){ // match
unsigned int i = (TSIZE > TSIZE_2) ? TSIZE_2 : TSIZE;
for(i--;i!=ExCo<unsigned int>::maximum();i--) op(data[i], a_2.data[i]);
}

	LFHTEMP
	template<class A_1, class A_2, class A_3, class C_2, unsigned int TSIZE_2, class C_3, unsigned int TSIZE_3> void Tuple<C,TSIZE,Cflag>::operator() (Oper3<A_1,A_2,A_3> const & op, Tuple<C_2, TSIZE_2,Cflag> const & a_2, Tuple<C_3, TSIZE_3,Cflag> const & a_3 ){ // not a match
		int i = (TSIZE > TSIZE_2) ? (TSIZE_2 > TSIZE_3 ? TSIZE_3-1 : TSIZE_2 -1) : (TSIZE > TSIZE_3 ? TSIZE_3-1 : TSIZE - 1);
		for(;i>=0;i--) data[i](op,a_2.data[i],a_3.data[i]);
	}
	LFHTEMP
	template<class A_1, class A_2, class A_3, class C_2, unsigned int TSIZE_2, class C_3, unsigned int TSIZE_3> void Tuple<C,TSIZE,Cflag>::operator() (Oper3<A_1,A_2,A_3> const & op, Tuple<C_2, TSIZE_2,Cflag> & a_2, Tuple<C_3, TSIZE_3,Cflag> const & a_3 ){ // not a match
		int i = (TSIZE > TSIZE_2) ? (TSIZE_2 > TSIZE_3 ? TSIZE_3-1 : TSIZE_2 -1) : (TSIZE > TSIZE_3 ? TSIZE_3-1 : TSIZE - 1);
		for(;i>=0;i--) data[i](op,a_2.data[i],a_3.data[i]);
	}
	LFHTEMP
	template<class A_1, class A_2, class A_3, class C_2, unsigned int TSIZE_2, class C_3, unsigned int TSIZE_3> void Tuple<C,TSIZE,Cflag>::operator() (Oper3<A_1,A_2,A_3> const & op, Tuple<C_2, TSIZE_2,Cflag> & a_2, Tuple<C_3, TSIZE_3,Cflag> & a_3 ){ // not a match
		int i = (TSIZE > TSIZE_2) ? (TSIZE_2 > TSIZE_3 ? TSIZE_3-1 : TSIZE_2 -1) : (TSIZE > TSIZE_3 ? TSIZE_3-1 : TSIZE - 1);
		for(;i>=0;i--) data[i](op,a_2.data[i],a_3.data[i]);
	}
	LFHTEMP
	template<class C_2, unsigned int TSIZE_2, class C_3, unsigned int TSIZE_3> void Tuple<C,TSIZE,Cflag>::operator() (Oper3<C,C_2,C_3> const & op, Tuple<C_2, TSIZE_2,Cflag> const & a_2, Tuple<C_3, TSIZE_3,Cflag> const & a_3){ // match
		int i = (TSIZE > TSIZE_2) ? (TSIZE_2 > TSIZE_3 ? TSIZE_3-1 : TSIZE_2 -1) : (TSIZE > TSIZE_3 ? TSIZE_3-1 : TSIZE - 1);
		for(;i>=0;i--) op(data[i], a_2.data[i], a_3.data[i]);
	}
	LFHTEMP
	template<class C_2, unsigned int TSIZE_2, class C_3, unsigned int TSIZE_3> void Tuple<C,TSIZE,Cflag>::operator() (Oper3<C,C_2,C_3> const & op, Tuple<C_2, TSIZE_2,Cflag> & a_2, Tuple<C_3, TSIZE_3,Cflag> const & a_3){ // match
		int i = (TSIZE > TSIZE_2) ? (TSIZE_2 > TSIZE_3 ? TSIZE_3-1 : TSIZE_2 -1) : (TSIZE > TSIZE_3 ? TSIZE_3-1 : TSIZE - 1);
		for(;i>=0;i--) op(data[i], a_2.data[i], a_3.data[i]);
	} 
	
	LFHTEMP
	template<class C_2, unsigned int TSIZE_2, class C_3, unsigned int TSIZE_3> void Tuple<C,TSIZE,Cflag>::operator() (Oper3<C,C_2,C_3> const & op, Tuple<C_2, TSIZE_2,Cflag> & a_2, Tuple<C_3, TSIZE_3,Cflag> & a_3){ // match
		int i = (TSIZE > TSIZE_2) ? (TSIZE_2 > TSIZE_3 ? TSIZE_3-1 : TSIZE_2 -1) : (TSIZE > TSIZE_3 ? TSIZE_3-1 : TSIZE - 1);
		for(;i>=0;i--) op(data[i], a_2.data[i], a_3.data[i]);
	}
	

	
LFHTEMP
	void Tuple<C,TSIZE,Cflag>::fourierTransform_routine(){
		// assumes TSIZE is a power of two, and C can be multiplied by a complex number
		int step;
		int cur;
		int icur;
		complex multip[TSIZE];
		multip[0] = complex(1.0f,0.0f);
		for(step=1;step<TSIZE;step<<=1){
			for(icur=1;icur<(step<<1);icur++){
				double tmpdouble =  ((double)icur * M_PI) / step;
				multip[icur] = complex(cos(tmpdouble),sin(tmpdouble));
			//	ExOp::show(multip[icur]);
			}
			for(cur =0;cur<TSIZE;cur+= (step<< 1)){
				for(icur=0;icur<step;icur++){
	//					printf("%f\t%f\t",((complex)data[cur | step | icur])[0],((complex)data[cur | step | icur])[1]);
						C tmp = data[cur | step | icur];
					//	ExOp::show(tmp);
					//	ExOp::show(multip[step | icur]);
						tmp *= multip[step | icur];
					//	printf("gives\n");
					//	ExOp::show(tmp);
			//			printf("%f\t%f\t",tmp[0],tmp[1]);
						tmp += data[cur | icur];
				//		printf("%f\t%f\t",tmp[0],tmp[1]);
						data[cur | step | icur] *= multip[icur];
						data[cur | icur] += data[cur | step | icur];
						data[cur | step | icur] = tmp;
//						printf("%f\t%f\n",((complex)data[cur | step | icur])[0],((complex)data[cur | step | icur])[1]);
				}
			}
		}
	}

	LFHTEMP
	void Tuple<C,TSIZE,Cflag>::invfourierTransform_routine(){
		// assumes TSIZE is a power of two, and C can be multiplied by a complex number
		int step;
		int cur;
		int icur;
		complex multip[TSIZE];

		multip[0] = complex(1.0f,0.0f);
		for(step=1;step<TSIZE;step<<=1){
			for(icur=1;icur<(step<<1);icur++){
				double tmpdouble =  -((double)icur * M_PI) / step;
				multip[icur] = complex(cos(tmpdouble),sin(tmpdouble));
			}
			for(cur =0;cur<TSIZE;cur+= (step<< 1)){
				for(icur=0;icur<step;icur++){
	//					printf("%f\t%f\t",((complex)data[cur | step | icur])[0],((complex)data[cur | step | icur])[1]);
						C tmp = data[cur | step | icur];
						tmp *= multip[step | icur];
			//			printf("%f\t%f\t",tmp[0],tmp[1]);
			//			data[cur | icur] *= half;
						tmp += data[cur | icur];
				//		printf("%f\t%f\t",tmp[0],tmp[1]);
						data[cur | step | icur] *= multip[icur];
						data[cur | icur] += data[cur | step | icur];
						data[cur | step | icur] = tmp;
//						printf("%f\t%f\n",((complex)data[cur | step | icur])[0],((complex)data[cur | step | icur])[1]);
				}
			}
		}
		complex half = complex(1/((double)TSIZE),0.0f);
		for(cur =0;cur<TSIZE;cur++) data[cur] *= half;

	}

	LFHTEMP
	void Tuple<C,TSIZE,Cflag>::fourierTransform_routine2(){
		// assumes TSIZE is a power of two, and C can be multiplied by a complex number
		int step;
		int icur;
        step = TSIZE >> 1;
            for(icur=0;icur<step;icur++){
                    C tmp = data[icur] - data[step | icur];
                    data[ icur] += data[step | icur];
                    data[step | icur] = tmp;
            }
		if (TSIZE <= 2) return;
		complex damultip;
		double tmpdouble;
		for(step=step >> 1;step != 0 ;step=step >> 1){
            // run
			for(icur =0;icur<TSIZE;icur+= step){
                tmpdouble = (M_PI * icur) / TSIZE;
			    damultip = complex(cos(tmpdouble),sin(tmpdouble));

				for(;(icur & step) == 0 ; icur++){
                        data[ step | icur] *= damultip;
                        C tmp = data[icur];
                        tmp -= data[ step | icur];
						data[icur] += data[step | icur];
						data[step | icur] = tmp;
				}
			}
		}
	}

	LFHTEMP
	void Tuple<C,TSIZE,Cflag>::invfourierTransform_routine2(){
		// assumes TSIZE is a power of two, and C can be multiplied by a complex number
				int step;
		int icur;

		if (TSIZE > 2) {
		complex damultip;
		double tmpdouble;
		for(step= 1 ; (step << 1) < TSIZE ;step=step << 1){
            // run
			for(icur =0;icur<TSIZE;icur+= step){
                tmpdouble = (-M_PI * icur) / TSIZE;
			    damultip = complex(cos(tmpdouble),sin(tmpdouble));
				for(;(icur & step) == 0 ; icur++){
                        C tmp = data[icur];
                        tmp -= data[ step | icur];
						data[icur] += data[step | icur];
						data[step | icur] = tmp;
						data[ step | icur] *= damultip;
				}
			}
		}
		}
		step = TSIZE >> 1;
            for(icur=0;icur<step;icur++){
                    C tmp = data[icur] - data[step | icur];
                    data[ icur] += data[step | icur];
                    data[step | icur] = tmp;
            }
		complex half = complex(1/((double)TSIZE),0.0f);
		for(icur =0;icur<TSIZE;icur++) data[icur] *= half;
	}

	LFHTEMP
	void Tuple<C,TSIZE,Cflag>::invfourierTransform_routine(complex* array){
		// assumes TSIZE is a power of two, and C can be multiplied by a complex number
		int step;
		int cur;
		int icur;
		complex multip[TSIZE];

		multip[0] = complex(1.0f,0.0f);
		for(step=1;step<TSIZE;step<<=1){
			for(icur=1;icur<(step<<1);icur++){
				double tmpdouble =  -((double)icur * M_PI) / step;
				multip[icur] = complex(cos(tmpdouble),sin(tmpdouble));
			}
			for(cur =0;cur<TSIZE;cur+= (step<< 1)){
				for(icur=0;icur<step;icur++){
						C tmp = array[cur | step | icur];
						tmp *= multip[step | icur];
						tmp += array[cur | icur];
						array[cur | step | icur] *= multip[icur];
						array[cur | icur] += array[cur | step | icur];
						array[cur | step | icur] = tmp;
				}
			}
		}
		complex half = complex(1/((double)TSIZE),0.0f);
		for(cur =0;cur<TSIZE;cur++) array[cur] *= half;

	}



LFHTEMP
	void Tuple<C,TSIZE,Cflag>::fourierTransform_routine(complex* array){
		// assumes TSIZE is a power of two, and C can be multiplied by a complex number
		int step;
		int cur;
		int icur;
		complex multip[TSIZE];
		multip[0] = complex(1.0f,0.0f);
		for(step=1;step<TSIZE;step<<=1){
			for(icur=1;icur<(step<<1);icur++){
				double tmpdouble =  ((double)icur * M_PI) / step;
				multip[icur] = complex(cos(tmpdouble),sin(tmpdouble));
			}
			for(cur =0;cur<TSIZE;cur+= (step<< 1)){
				for(icur=0;icur<step;icur++){
						C tmp = array[cur | step | icur];
						tmp *= multip[step | icur];
						tmp += array[cur | icur];
						array[cur | step | icur] *= multip[icur];
						array[cur | icur] += array[cur | step | icur];
						array[cur | step | icur] = tmp;
				}
			}
		}
	}


LFHTEMP
Tuple<complex, TSIZE,Cflag> Tuple<C,TSIZE,Cflag>::bluesteinWindow(int subTSIZE){ // the uutput should be cached...
	Tuple<C, TSIZE,Cflag> out;
	int i;

	int j,s;
	j=0;
			for(i=0;i<subTSIZE;i++) {
				double ang = i*i*M_PI/ subTSIZE;
				out[j] =  complex(cos(ang),sin(ang));
	//			printf("b:%i: %f\t%f\n",i,out[j][0],out[j][1]);
				for(s = TSIZE >> 1;(j & s);s >>=1) j ^= s;
				j |= s;
			}
			for(;i<TSIZE-subTSIZE+1;i++) {
				out[j] =  ExCo<C>::zero();
	//			printf("b:%i: %f\t%f\n",i,out[j][0],out[j][1]);
				for(s = TSIZE >> 1;(j & s);s >>=1) j ^= s;
				j |= s;
			}
			for(;i<TSIZE;i++) {
				double ang = (TSIZE-i)*(TSIZE-i)*M_PI/ subTSIZE;
				out[j] =  complex(cos(ang),sin(ang));
	//			printf("b:%i: %f\t%f\n",i,out[j][0],out[j][1]);
				for(s = TSIZE >> 1;(j & s);s >>=1) j ^= s;
				j |= s;
			}

		out.fourierTransform_routine();

	return(out);
	}

LFHTEMP
void Tuple<C,TSIZE,Cflag>::bluesteinWindow(complex *& fout, int subTSIZE){ // the uutput should be cached...
	int i;
    fout = new complex[TSIZE];
	int j,s;
	j=0;
			for(i=0;i<subTSIZE;i++) {
				double ang = i*i*M_PI/ subTSIZE;
				fout[j] =  complex(cos(ang),sin(ang));
	//			printf("b:%i: %f\t%f\n",i,out[j][0],out[j][1]);
				for(s = TSIZE >> 1;(j & s);s >>=1) j ^= s;
				j |= s;
			}
			for(;i<TSIZE-subTSIZE+1;i++) {
				fout[j] =  ExCo<C>::zero();
	//			printf("b:%i: %f\t%f\n",i,out[j][0],out[j][1]);
				for(s = TSIZE >> 1;(j & s);s >>=1) j ^= s;
				j |= s;
			}
			for(;i<TSIZE;i++) {
				double ang = (TSIZE-i)*(TSIZE-i)*M_PI/ subTSIZE;
				fout[j] =  complex(cos(ang),sin(ang));
	//			printf("b:%i: %f\t%f\n",i,out[j][0],out[j][1]);
				for(s = TSIZE >> 1;(j & s);s >>=1) j ^= s;
				j |= s;
			}

		Tuple<C,TSIZE,Cflag>::fourierTransform_routine(fout);

	}

LFHTEMP
template<unsigned int superTSIZE> Tuple<C, TSIZE,Cflag> Tuple<C,TSIZE,Cflag>::fourierTransform(const Tuple<complex, superTSIZE,Cflag>& bluewindow) const{

		// implement Bluestein algorithm!

		// B windows : 2 {1.0+1.0i, 1.0-1.0i }
		// B windows : 4 {0.0+0.0i, 2.0-0.0i,(1 + i) * sqrt(2)/2, (1 + i) * sqrt(2)/2}


		Tuple<C, superTSIZE,Cflag> tmp;
		int i;
		int j=0;
		int s;

		for(i=0;i<TSIZE;i++){
				double ang = -i * i * M_PI / TSIZE;
				complex factor = complex(cos(ang),sin(ang));
				tmp[j] = data[i];
				tmp[j] *= factor;
				for(s = superTSIZE >> 1;(j & s);s >>=1) j ^= s;
				j |= s;
			}
		for(;i<superTSIZE;i++){
				tmp[j] = ExCo<C>::zero();
				for(s = superTSIZE >> 1;(j & s);s >>=1) j ^= s;
				j |= s;
			}

//				for(i=0;i<superTSIZE;i++) printf("1:%i: %f\t%f\n",i,tmp[i][0],tmp[i][1]);
		tmp.fourierTransform_routine();
//				for(i=0;i<superTSIZE;i++) printf("1:%i: %f\t%f\n",i,tmp[i][0],tmp[i][1]);

		Tuple<C, superTSIZE,Cflag> tmp2;
//		Tuple<C, superTSIZE,Cflag> tmpb = bluesteinWindow<superTSIZE>();

//		for(i=0;i<superTSIZE;i++) printf("b:%i: %f\t%f\n",i,tmpb[i][0],tmpb[i][1]);

		for(i=0;i<superTSIZE;i++){
				tmp2[j]	= tmp[i] * bluewindow[i];
				for(s = superTSIZE >> 1;(j & s);s >>=1) j ^= s;
				j |= s;
			}

//		for(i=0;i<superTSIZE;i++) printf("2:%i: %f\t%f\n",i,tmp2[i][0],tmp2[i][1]);
		tmp2.invfourierTransform_routine();
//		for(i=0;i<superTSIZE;i++) printf("2:%i: %f\t%f\n",i,tmp2[i][0],tmp2[i][1]);

		Tuple<C, TSIZE,Cflag> out;
		for(i=0;i<TSIZE;i++) {
			 double ang = -i * i * M_PI / TSIZE;
			complex factor = complex(cos(ang),sin(ang));
			 out[i] = tmp2[i];
			 out[i] *= factor;

			 }
		return(out);
	}

LFHTEMP
template<unsigned int superTSIZE> Tuple<C, TSIZE,Cflag> Tuple<C,TSIZE,Cflag>::invfourierTransform(const Tuple<complex, superTSIZE,Cflag>& bluewindow) const{

		// implement Bluestein algorithm!

		// B windows : 2 {1.0+1.0i, 1.0-1.0i }
		// B windows : 4 {0.0+0.0i, 2.0-0.0i,(1 + i) * sqrt(2)/2, (1 + i) * sqrt(2)/2}


		Tuple<C, superTSIZE,Cflag> tmp;
		int i;
		int j=0;
		int s;

		for(i=0;i<TSIZE;i++){
				double ang = i * i * M_PI / TSIZE;
				complex factor = complex(cos(ang),sin(ang));
				tmp[j] = data[i];
				tmp[j] *= factor;
				for(s = superTSIZE >> 1;(j & s);s >>=1) j ^= s;
				j |= s;
			}
		for(;i<superTSIZE;i++){
				tmp[j] = ExCo<C>::zero();
				for(s = superTSIZE >> 1;(j & s);s >>=1) j ^= s;
				j |= s;
			}

//				for(i=0;i<superTSIZE;i++) printf("1:%i: %f\t%f\n",i,tmp[i][0],tmp[i][1]);
		tmp.fourierTransform_routine();
//				for(i=0;i<superTSIZE;i++) printf("1:%i: %f\t%f\n",i,tmp[i][0],tmp[i][1]);

		Tuple<C, superTSIZE,Cflag> tmp2;
//		Tuple<C, superTSIZE,Cflag> tmpb = bluesteinWindow<superTSIZE>();

//		for(i=0;i<superTSIZE;i++) printf("b:%i: %f\t%f\n",i,tmpb[i][0],tmpb[i][1]);

		for(i=0;i<superTSIZE;i++){
				tmp2[j]	= tmp[i];
				complex factor =bluewindow[i];
				factor[1] = -factor[1];
				tmp2[j]	*= factor;
				for(s = superTSIZE >> 1;(j & s);s >>=1) j ^= s;
				j |= s;
			}

//		for(i=0;i<superTSIZE;i++) printf("2:%i: %f\t%f\n",i,tmp2[i][0],tmp2[i][1]);
		tmp2.invfourierTransform_routine();
//		for(i=0;i<superTSIZE;i++) printf("2:%i: %f\t%f\n",i,tmp2[i][0],tmp2[i][1]);

		Tuple<C, TSIZE,Cflag> out;
		for(i=0;i<TSIZE;i++) {
			 double ang = i * i * M_PI / TSIZE;
			complex factor = complex(cos(ang)/TSIZE,sin(ang)/TSIZE);
			 out[i] = tmp2[i];
			 out[i] *= factor;

			 }
		return(out);
	}

	LFHTEMP
	C Tuple<C,TSIZE,Cflag>::max() const{
		C _out = data[0];
		int i;
		for(i=1;i<TSIZE;i++) if (data[i] > _out) _out = data[i];
		return(_out);
	}
	LFHTEMP
	C Tuple<C,TSIZE,Cflag>::min() const{
		C _out = data[0];
		int i;
		for(i=1;i<TSIZE;i++) if (data[i] < _out) _out = data[i];
		return(_out);
	}

	LFHTEMP
	Tuple<C, TSIZE, Cflag> Tuple<C,TSIZE,Cflag>::normalize() const{
		C sum = data[0];
		Tuple<C, TSIZE, Cflag> _out;
		int i;
		for(i=1;i<TSIZE;i++) sum += data[i];
//		C sumi = ExCo<C>::invert(sum);
		for(i=0;i<TSIZE;i++) if (!ExCo<C>::isValid(_out[i] = data[i] / sum)) break;
		if (i == TSIZE) return(_out);
		// sum or sumi is NAN, normalize by the max beforehand

		C cax = max();
		C cix = min();
		sum = (cax > -cix) ? cax : -cix;
		if (sum  == ExCo<C>::zero()) {
			// prior guess
			_out[i] =  ExCo<C>::one() * (1.0f / TSIZE);
			for(i=1;i<TSIZE;i++) _out[i] = _out[0];
			return(_out);// all entries are zeroes
		}
		for(i=0;i<TSIZE;i++) _out[i] = data[i] / sum;
		sum = _out[0];
		for(i=1;i<TSIZE;i++) sum += _out[1];
		for(i=0;i<TSIZE;i++) _out[i] /= sum;
		return(_out);
		}
	LFHTEMP
	Tuple<C, TSIZE, Cflag> Tuple<C,TSIZE,Cflag>::normalize(C const & norm) const{
		C sum = data[0];
		Tuple<C, TSIZE, Cflag> _out;
		int i;
		for(i=1;i<TSIZE;i++) sum += data[i];
		C sumi = (ExCo<C>::invert(sum)) * norm;
		for(i=0;i<TSIZE;i++) _out[i] = data[i] * sumi;
		return(_out);
		}


LFHTEMP
template<unsigned int superTSIZE> Tuple<Complex<C>, TSIZE,Cflag> Tuple<C,TSIZE,Cflag>::fourierTransformReal() const{
		//Tuple<Complex<C>, TSIZE,Cflag>& _out;

		//Tuple<Complex<C>, TSIZE,Cflag>& _out;

		//int s

	}





	/*
template <class C, int nbdim>
const SetComparison& HyperPositionProjectCube<C,nbdim>::compare(const HyperPosition<C,nbdim> & other)  const {

	}
	*/
/*
LFHTEMP
template<class I, class O, class MC> Tuple<MC,TSIZE> Tuple<C,TSIZE,Cflag>::operator()(Oper<I,O> const & op) const{
	Tuple<MC,TSIZE> _out;
	int i;
	for(i=0;i<TSIZE;i++) _out = data[i](op);
	return(_out);
}

LFHTEMP
template<class I, class O> Tuple<O,TSIZE> Tuple<C,TSIZE,Cflag>::operator()(Oper< Tuple<I,TSIZE> ,Tuple<O,TSIZE> > const & op){
	Tuple<O,TSIZE> _out;
	int i;
	for(i=0;i<TSIZE;i++) _out = op(data[i]);
	return(_out);
}

LFHTEMP
template<class O> O Tuple<C,TSIZE,Cflag>::operator()(Oper< Tuple<C,TSIZE,Cflag> ,O> const & op) const {return(op(*this));}
	*/

LFHTEMP
Tuple<C,TSIZE,Cflag>::Tuple(C const* const clonefrom){memcpy(data,clonefrom, sizeof(C)*(TSIZE));}

LFHTEMP
template<class O, unsigned int oTSIZE>
Tuple<C,TSIZE,Cflag>::Tuple(Tuple<O,oTSIZE,Cflag> const & other){
	(*this) = other;
	}

LFHTEMP
const Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator-() const{
	Tuple<C,TSIZE,Cflag> tmp;
	int i;
	for(i = TSIZE-1;i>=0;i--) tmp[i] = -data[i];
	return(Tuple<C,TSIZE,Cflag>(tmp));
	}

LFHTEMP Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator=(Tuple<C,TSIZE,Cflag> const & other){
	unsigned int i;
	for(i = TSIZE-1;i < TSIZE ;i--) data[i] = other[i];
	return(*this);
	}
LFHTEMP Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::memmove(Tuple<C,TSIZE,Cflag> & other){
	unsigned int i;
	for(i = TSIZE-1;i < TSIZE ;i--) ExOp::memmove(data[i],other[i]);
	return(*this);
}

LFHTEMP template <class O> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator=(O const & other){
	unsigned int i;
	for(i = TSIZE-1;i < TSIZE ;i--) data[i] = other;
	return(*this);
	}

//LFHTEMP Tuple<C,TSIZE,Cflag>::operator Vector<C,LFHVECTOR_REMOTE> (){return Vector<C,LFHVECTOR_REMOTE>(data,TSIZE);}


LFHTEMP template <class O> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator+=(O const & other){
	int i;
	for(i = TSIZE-1;i>=0;i--) data[i] += other;
	return(*this);
	}
LFHTEMP template <class O> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator-=(O const & other){
	int i;
	for(i = TSIZE-1;i>=0;i--) data[i] -= other;
	return(*this);
	}

LFHTEMP template <class O> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator*=(O const & other){
	int i;
	for(i = TSIZE-1;i>=0;i--) data[i] *= other;
	return(*this);
}
LFHTEMP template <class O> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator/=(O const & other){
	int i;
	for(i = TSIZE-1;i>=0;i--) data[i] /= other;
	return(*this);
	}
LFHTEMP template <class O, unsigned int oTSIZE> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator=(Tuple<O,oTSIZE,Cflag> const & other){
	int i = (oTSIZE > TSIZE) ? TSIZE-1 : oTSIZE-1;
	for(;i>=0;i--) data[i] = other[i];
	return(*this);
	}

	LFHTEMP template <class O>
	Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator+=(Tuple<O,TSIZE,Cflag> const & other){
		int i;
		for(i = TSIZE-1;i>=0;i--) data[i] += other.data[i];
		return(*this);
	}
	LFHTEMP template <class O>
	Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator-=(Tuple<O,TSIZE,Cflag> const & other){
		int i;
		for(i = TSIZE-1;i>=0;i--) data[i] -= other.data[i];
		return(*this);
	}

	LFHTEMP template <class O>
	Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator*=(Tuple<O,TSIZE,Cflag> const & other){
		int i;
		for(i = TSIZE-1;i>=0;i--) data[i] *= other.data[i];
		return(*this);
	}
	LFHTEMP template <class O>
	Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator/=(Tuple<O,TSIZE,Cflag> const & other){
		int i;
		for(i = TSIZE-1;i>=0;i--) data[i] /= other.data[i];
		return(*this);
	}

	LFHTEMP template <class O> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator*=(TMatrix<O,TSIZE,TSIZE> const & other){
		Tuple<C,TSIZE,Cflag> clo(*this);
		unsigned int i,k;
		
		for(i=0;i<TSIZE;i++) data[i] = other.data[i * TSIZE] * clo[0];
		for(k=1;k<TSIZE;k++){
			for(i=0;i<TSIZE;i++) data[i] += other.data[k + i * TSIZE] * clo[k];
		}
		return(*this);
	}
	
	LFHTEMP template <class O> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator/=(TMatrix<O,TSIZE,TSIZE> const & other){
		return(*this);
	}

	

LFHTEMP template <class O> Tuple<C,TSIZE,Cflag> Tuple<C,TSIZE,Cflag>::operator+(KeyElem<unsigned int, O> const & other) const{return( (Tuple<C,TSIZE,Cflag>(*this)) += other );}
LFHTEMP template <class O> Tuple<C,TSIZE,Cflag> Tuple<C,TSIZE,Cflag>::operator-(KeyElem<unsigned int, O> const & other) const{return( (Tuple<C,TSIZE,Cflag>(*this)) -= other );}
LFHTEMP template <class O> Tuple<C,TSIZE,Cflag> Tuple<C,TSIZE,Cflag>::operator*(KeyElem<unsigned int, O> const & other) const{return( (Tuple<C,TSIZE,Cflag>(*this)) *= other );}
LFHTEMP template <class O> Tuple<C,TSIZE,Cflag> Tuple<C,TSIZE,Cflag>::operator/(KeyElem<unsigned int, O> const & other) const{return( (Tuple<C,TSIZE,Cflag>(*this)) /= other );}
LFHTEMP template <class O> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator+=(KeyElem<unsigned int, O> const & other) {data[other.k] += other.d; return(*this);}
LFHTEMP template <class O> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator-=(KeyElem<unsigned int, O> const & other) {data[other.k] -= other.d; return(*this);}
LFHTEMP template <class O> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator*=(KeyElem<unsigned int, O> const & other) {data[other.k] *= other.d; return(*this);}
LFHTEMP template <class O> Tuple<C,TSIZE,Cflag>& Tuple<C,TSIZE,Cflag>::operator/=(KeyElem<unsigned int, O> const & other) {data[other.k] /= other.d; return(*this);}


	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	Tuple<C,TSIZE,Cflag> const &  Tuple<C,TSIZE,Cflag>::genConvolution_Gaussian(double std){
		Tuple<C,TSIZE,Cflag> out;
		int i,j;
		double tmp;
		for(i=0;i<TSIZE;i++){
			tmp = (i -(TSIZE >>1)) /std;
			out.data[i] = (C)(exp(-tmp*tmp/2.0f));
			for(j=0;j<TSIZE*100000000;j++) tmp += j;
			}

		return(out);
		}

			// OPERATOR!!!


	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	template<unsigned int pos>
	void Tuple<C,TSIZE,Cflag>::Selector<pos>::operator()(C & _out, Tuple<C,TSIZE,Cflag> & input) const{
		_out = input[pos];
	}


	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	void Tuple<C,TSIZE,Cflag>::MakeTuple::operator()(Tuple<C,TSIZE,Cflag> & output, C & _in) const{
		for(unsigned int i=0;i<TSIZE;i++)	output[i] = _in;
	}

    template<class C, unsigned int TSIZE, Tuple_flag Cflag>
    void Tuple<C,TSIZE,Cflag>::HouseHolderMultiply(const C * const vec, double denum2, unsigned int length ){
        unsigned int i;
		if ((denum2 == 0.0f)||(!ExCo<double>::isValid(denum2))) return;
        C s = data[0] * vec[0];
        for(i=1;i<length;i++) s+= data[i] * vec[i];
        s = ExOp::mktrju(s) / (-denum2);
        for(i=0;i<length;i++) data[i] += s * vec[i];

    }

	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	template<unsigned int fTSIZE>
	void Tuple<C,TSIZE,Cflag>::Concatenate<fTSIZE>::operator()(Tuple<C,TSIZE,Cflag> & _out, Tuple<C,fTSIZE,Cflag> & input, Tuple<C,TSIZE-fTSIZE,Cflag> & input2) const{
		int i;
		for(i=0;i<fTSIZE;i++) _out.data[i] = input[i];
		for(;i<TSIZE;i++) _out.data[i] = input2[i-fTSIZE];

	}

	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	template<unsigned int fTSIZE>
	void Tuple<C,TSIZE,Cflag>::Deconcatenate<fTSIZE>::operator()(Tuple<C,TSIZE,Cflag> & _out, Tuple<C,fTSIZE,Cflag> & out_2, Tuple<C,TSIZE+fTSIZE,Cflag> & input) const{
		int i;
		for(i=0;i<TSIZE;i++) _out[i] = input.data[i];
		for(;i<TSIZE+ fTSIZE;i++) out_2[i-TSIZE] = input.data[i];

	}

	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	template<unsigned int pos> void Tuple<C,TSIZE,Cflag>::SelectorWrite<pos>::operator()(Tuple<C,TSIZE,Cflag> & _out, C & _input) const{
		_out[pos] = _input;
	}

	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	template<unsigned int TSIZEin, unsigned int pos_in, unsigned int pos_out> void Tuple<C,TSIZE,Cflag>::SelectSelectorWrite<TSIZEin,pos_in,pos_out>::operator()(Tuple<C,TSIZE,Cflag> & _out, Tuple<C,TSIZEin,Cflag> & _input) const{
		_out[pos_out] = _input[pos_in];
	}

	template<class C>
	void Square<C>::operator()(C & _out, C & input) const{
		_out =input*input;
	}
	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	void Tuple<C,TSIZE,Cflag>::L2Norm::operator()(C & _out, Tuple<C,TSIZE,Cflag> & input) const{
		_out = input[0] *input[0];
		int i;
		for(i=1;i<TSIZE;i++) _out += input[1] * input[1];
	}

	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	void Tuple<C,TSIZE,Cflag>::FiniteDifferenceAssign::operator()(Tuple<C,TSIZE,Cflag> &_out) const{
		C tmp = _out[0] - _out[TSIZE-1];
		int i;
		for(i=0;i<TSIZE-1;i++){
			_out[i] = _out[i+1] - _out[i];
			}
			_out[i] = tmp;
		}

	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	void Tuple<C,TSIZE,Cflag>::FiniteDifference::operator()(Tuple<C,TSIZE,Cflag> &_out, Tuple<C,TSIZE,Cflag> & input) const{
		_out[TSIZE-1] = input[0] - input[TSIZE-1];
		int i;
		for(i=0;i<TSIZE-1;i++){
			_out[i] = input[i+1] - input[i];
			}
		}

	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	void Tuple<C,TSIZE,Cflag>::ArgMax::operator()(int &_out, Tuple<C,TSIZE,Cflag> & input) const{
		int b =0;
		C curb = input[0];
		int i;
		for(i=1;i<TSIZE;i++){
			if (input[i] > curb) {curb = input[i]; b= i; }
			}
		_out = b;
		}

	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	void Tuple<C,TSIZE,Cflag>::MassCenter::operator()(double &_out, Tuple<C,TSIZE,Cflag> & input) const{
		double tm = 0.0f;
		double sum = 0.0f;
		double tmp;
		int i;
		for(i=0;i<TSIZE;i++){
			tmp = (double) input[i]*input[i];
			sum += pow(2.0f,0.5f*(i+1-TSIZE)) * tmp;
			tm += tmp;
			}
		_out = sum / tm;
		}

	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	void Tuple<C,TSIZE,Cflag>::MassCenter2::operator()(Tuple<double,3,Cflag> &_out, Tuple<C,TSIZE,Cflag> & input) const{
		double tm = (double) input[0];
		double sum = 0.0f;
		double tmp;
		int i;
		for(i=1;i<TSIZE;i++){
			tmp = (double) input[i];
			sum += i * tmp;
			tm += tmp;
			}
		_out[0] = sum / tm;
		_out[1] = tm;
		_out[2] = sum / tm;
		}

	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	void Tuple<C,TSIZE,Cflag>::PopWeight::operator()(Tuple<C,TSIZE,Cflag> &_out, Tuple<C, TSIZE+1,Cflag> & input) const{
		int i;
		for(i=0;i<TSIZE;i++) _out[i] = input[i] / input[TSIZE];
	}


	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	void Tuple<C,TSIZE,Cflag>::PushWeight::operator()(Tuple<C,TSIZE,Cflag> &_out, Tuple<C, TSIZE-1,Cflag> & input, C & w) const{
		int i;
		_out[TSIZE-1] = w;
		for(i=0;i<TSIZE-1;i++) _out[i] = input[i] * w;
	}

	template<class C, unsigned int TSIZE, Tuple_flag Cflag>
	template<unsigned int winTSIZE>
	void Tuple<C,TSIZE,Cflag>::Convolution<winTSIZE>::operator()(Tuple<C,TSIZE,Cflag> &_out, Tuple<C,TSIZE,Cflag> & input, Tuple<C, winTSIZE,Cflag> & convwind) const{
		int i,j;
		for(i=0;i<TSIZE;i++){

			if (i < (winTSIZE >> 1)) j=(winTSIZE >> 1) - i;
			else j = 0;
			_out[i] = input[i+j-(winTSIZE >> 1)] * convwind[j];
			for(j++;(j<winTSIZE)&&(i+j-(winTSIZE >> 1) < TSIZE);j++) _out[i] += input[i+j-(winTSIZE >> 1)] * convwind[j];
			}
	}

#undef LFHTEMP
#define LFHTEMP template<class C, Tuple_flag Cflag>


LFHTEMP Tuple<C,0u,Cflag>::Tuple(const Tuple<C,0u, Cflag> &other){if (other.tup_size) {tup_size = other.tup_size; data = new C[other.tup_size]; LFH_USE_MEMLEAK_HELP(memleak_helper[(unsigned int)data] =1;) for(unsigned int i=0;i<tup_size;i++) data[i] = other.data[i];}else tup_size = 0u;}
LFHTEMP template<unsigned int OSIZE> Tuple<C,0u,Cflag>::Tuple(const Tuple<C,OSIZE, Cflag> &other){tup_size = OSIZE; data = new C[OSIZE]; LFH_USE_MEMLEAK_HELP(memleak_helper[(unsigned int)data] =1;) for(unsigned int i=0;i<OSIZE;i++) data[i] = other.data[i];}

LFHTEMP Tuple<C,0u,Cflag>::~Tuple(){if (tup_size) {delete[](data); LFH_USE_MEMLEAK_HELP(memleak_helper.erase((unsigned int)data);)}}
LFHTEMP void Tuple<C,0u,Cflag>::setSize(unsigned int s) {
	if (s != tup_size) { if (tup_size) {delete[](data); LFH_USE_MEMLEAK_HELP(memleak_helper.erase((unsigned int)data);)}  data = (s == 0) ? NULL : new C[s]; tup_size=s; LFH_USE_MEMLEAK_HELP(if (data) memleak_helper[(unsigned int)data] =1;)}}



 LFHTEMP void Tuple<C,0u,Cflag>::HouseHolderMultiply(const C * const vec, double denum2, unsigned int length ){
        unsigned int i;
		if ((denum2 == 0.0f)||(!ExCo<double>::isValid(denum2)))  return;
        C s = data[0] * vec[0];
        for(i=1;i<length;i++) s+= data[i] * vec[i];
        s = ExOp::mktrju(s) / (-denum2);
        for(i=0;i<length;i++) data[i] += s * vec[i];
  }


LFHTEMP void Tuple<C,0u,Cflag>::show( FILE *f_out, int l) const{
	switch(l){
		case 1:
		case 0:
			for(unsigned int i=0;i<tup_size;i++) {
				if (i != 0) fprintf(f_out,"\t");
				ExOp::show(data[i],f_out,2);
			}
			if (l==0) fprintf(f_out,"\n");
		break;
		case 2:
			fprintf(f_out,"[");
			for(unsigned int i=0;i<tup_size;i++) {
				if (i != 0) fprintf(f_out,";");
				ExOp::show(data[i],f_out,3);
			}
			fprintf(f_out,"]");
			break;
		default:
			fprintf(f_out,"(");
			for(unsigned int i=0;i<tup_size;i++) {
				if (i != 0) fprintf(f_out,",");
				ExOp::show(data[i],f_out,4);
			}
			fprintf(f_out,")");
			break;
	}
}

LFHTEMP string Tuple<C,0u,Cflag>::type_tostring() const{return string("Tuple<") + ExOp::type_tostring(data[0]) +string(",0u>");}

#undef LFHTEMP
#define LFHTEMP template<class C,unsigned int order>

	LFHTEMP void WeightElem<C,order>::zero(){for(unsigned int i=0;i<order;i++) {w[i] =0.0f; ExOp::toZero(e[i]);} }
	LFHTEMP void WeightElem<C,order>::random(){for(unsigned int i=0;i<order;i++) ExOp::toRand(e[i]); for(unsigned int i=0;i<order;i++) w[i] =1.0f;}
	LFHTEMP void WeightElem<C,order>::one(){ExOp::toOne(e); for(unsigned int i=0;i<order;i++) {w[i] =1.0f;} }
	LFHTEMP void WeightElem<C,order>::save( FILE *f) const {ExOp::save(e,f); ExOp::save(w,f);}
	LFHTEMP void WeightElem<C,order>::load( FILE *f, unsigned int TSIZE) {ExOp::load(e,f);ExOp::load(w,f);}
	LFHTEMP void WeightElem<C,order>::show( FILE *f_out, int l) const{
		switch(l){
			case 0:
				fprintf(f_out,"average:"); ExOp::show(getMean(), f_out,1);
				if (order >1) {fprintf(f_out,"\nvariance:"); ExOp::show(getVar(), f_out,1);}
				if (order >2) {fprintf(f_out,"\nskew:"); ExOp::show(getSkew(), f_out,1);}
				if (order >3) {fprintf(f_out,"\nkurt:"); ExOp::show(getKurt(), f_out,1);}
				fprintf(f_out,"\n");
			break;
			case 1:
				fprintf(f_out,"average:"); ExOp::show(getMean(), f_out,2);
				if (order >1) {fprintf(f_out,"\tvariance:"); ExOp::show(getVar(), f_out,2);}
				if (order >2) {fprintf(f_out,"\tskew:"); ExOp::show(getSkew(), f_out,2);}
				if (order >3) {fprintf(f_out,"\tkurt:"); ExOp::show(getKurt(), f_out,2);}
			break;
		//	case 2:
		//		fprintf(f_out,"[");
		//		for(int i=0;i<TSIZE;i++) {
		//			if (i != 0) fprintf(f_out,";");
		//			ExOp::show(data[i],f_out,3);
		//		}
		//		fprintf(f_out,"]");
		//		break;
		//	default:
		//		fprintf(f_out,"(");
		//		for(int i=0;i<TSIZE;i++) {
		//			if (i != 0) fprintf(f_out,",");
		//			ExOp::show(data[i],f_out,4);
		//		}
		//		fprintf(f_out,")");
		}
	}
LFHTEMP string WeightElem<C,order>::type_tostring() const{char buffer[256]; sprintf(buffer,"%u",order); return string("WeightElem<") + ExOp::type_tostring(e[0]) +string(",")+ ExOp::type_tostring(buffer) + string(">");}


	LFHTEMP WeightElem<C,order> WeightElem<C,order>::operator-(const WeightElem<C,order>& s) const{ WeightElem<C,order> f_out;
		f_out.w = w + s.w;
		for(int i=0;i< order;i++) f_out.e[i] = (i & 1) ? e[i] + s.e[i] : e[i] - s.e[i];
		return(f_out);
	}




	LFHTEMP void WeightElem<C, order>::OpMean::operator()(C & out, WeightElem<C, order> & in) const {out = in.getMean();}



	LFHTEMP void WeightElem<C, order>::OpVar::operator()(C & out, WeightElem<C, 2> & in) const {
		LinkAssert<(order>1) > ass;
		out = in.getVar();
//	printf("%f,%f\n");
	}


	LFHTEMP void WeightElem<C, order>::OpStd::operator()(C & out, WeightElem<C, 2> & in) const { out = ExCo<C>::invintPow(in.getVar(),2);}





#undef LFHTEMP
#define LFHTEMP template<class C, unsigned int sizex, unsigned int sizey>
	LFHTEMP
	TMatrix<C,sizex,sizey>::TMatrix(TMatrix<C, sizex, sizey> const & idata){
		for(unsigned int k=0;k<sizex*sizey;k++) data[k] = idata.data[k];
	}

	LFHTEMP TMatrix<C,sizex,sizey>::TMatrix(DataGrid<C, 2> const & idata){
		Tuple<unsigned int,2> coor;
		for(coor[1] = 0; coor[1] < sizey;coor[1]++)
			for(coor[0] = 0; coor[0] < sizex;coor[0]++) data[coor[0] + coor[1]*sizex] = idata(coor);
	}

	LFHTEMP
	TMatrix<C,sizex,sizey>& TMatrix<C,sizex,sizey>::operator=(TMatrix<C, sizex, sizey> const & idata){
		for(unsigned int k=0;k<sizex*sizey;k++) data[k] = idata.data[k];
		return(*this);
	}





	LFHTEMP
	void TMatrix<C,sizex,sizey>::bidiag(){

		int i,j,k;
		double tau, lead;
		for(i=0;i<sizex;i++){

			for(j=i+1;j<sizex;j++){


			}
		}


	}

	LFHTEMP
	void TMatrix<C,sizex,sizey>::operator()(Tuple<C, sizey> &_out, Tuple<C, sizex> &_in) const{
		_out = (*this) * _in;
	}

	LFHTEMP
	void TMatrix<C,sizex,sizey>::derivative(Tuple<C, sizey> &_out, Tuple<C, sizex > &_in, int in_direct) const{
		int i;
		for(i=0;i<sizey;i++) _out[i] = data[in_direct + i * sizex];

	}

	LFHTEMP
	void TMatrix<C,sizex,sizey>::derivativeMatrix(TMatrix<C, sizex, sizey> & _out, Tuple<C, sizex > &_in) const{
		// HAHA, the derivative Matrix is the Matrix itself, what ever _in is

		_out = (*this);
	}

	LFHTEMP
	void TMatrix<C,sizex,sizey>::set_to_zero(){
		for(unsigned int i=0;i<sizex*sizey;i++) ExOp::toZero(data[i]);
	}

	LFHTEMP
	void TMatrix<C,sizex,sizey>::show(FILE* o, int level) const{
		unsigned int i,j;
		for(j=0;j<sizey;j++){
			for(i=0;i<sizex;i++){
				ExOp::show(data[i + j*sizex],o,1);
				fprintf(o,"%c", i == sizex-1 ? '\n' : '\t');
				}
			}

	}


	LFHTEMP
	TMatrix<C,sizex,sizey> TMatrix<C,sizex,sizey>::inverse() const{
		DataGrid<double, 2> dainv;
		Vector<double> eigen;
		dainv.makeDiagonalizer_ofinverse(*this,eigen);
		TMatrix<double, sizex, sizey> fun(dainv);
		TMatrix<double, sizex, sizey> funt = fun.transpose();
		for(unsigned int i=0;i<sizex;i++) for(unsigned int j=0;j<sizey;j++) fun.data[i+j*sizex] *= eigen[i];
		fun *= funt;
		return(fun);

	}


	LFHTEMP
	TMatrix<C,sizey,sizex> TMatrix<C,sizex,sizey>::transpose() const{ TMatrix<C,sizey,sizex> f_out;
		for(unsigned int i= 0 ; i< sizex;i++)
			for(unsigned int j= 0 ; j< sizey;j++) f_out.data[j+i*sizey] = data[i+j*sizex];
		return(f_out);
	}


	LFHTEMP
	C TMatrix<C,sizex,sizey>::determinant() const{
		DataGrid<double, 2> dainv;
		Vector<double> eigen;
		dainv.makeDiagonalizer_ofinverse(*this,eigen);
		C fout = (C) eigen[0];
		for(unsigned int i=1;i<sizex;i++) fout *= eigen[i];
		return(fout);

	}

	LFHTEMP TMatrix<C,sizex,sizey>& TMatrix<C,sizex,sizey>::toZero(){
		unsigned int i = sizex*sizey;
		for(i--;i != ExCo<unsigned int>::maximum() ; i--) ExOp::toZero(data[i]);
		return(*this);
	}
	LFHTEMP TMatrix<C,sizex,sizey>& TMatrix<C,sizex,sizey>::toOne(){
		unsigned int x,y,z;
		for(y=0,z=0;y<sizey;y++)
			for(x=0;x<sizex;x++,z++)
				if (x == y) ExOp::toOne(data[z]);
				else ExOp::toZero(data[z]);
		return(*this);
	}
	LFHTEMP TMatrix<C,sizex,sizey>& TMatrix<C,sizex,sizey>::toRand(){
		unsigned int i = sizex*sizey;
		for(i--;i != ExCo<unsigned int>::maximum() ; i--) ExOp::toRand(data[i]);
		return(*this);
	}
	LFHTEMP TMatrix<C,sizex,sizey>& TMatrix<C,sizex,sizey>::toRandSymetric(){
		unsigned int x,y,z;
		for(y=0,z=0;y<sizey;y++)
			for(x=0;x<sizex;x++,z++) if (y > x) data[z] = data[y + x * sizex]; else ExOp::toRand(data[z]);
		return(*this);
	}

	LFHTEMP
	template<class O,Tuple_flag FO>
	TMatrix<C,sizex,sizey> TMatrix<C,sizex,sizey>::scale_rows(Tuple<O, sizey, FO> const &fact) const{
		TMatrix<C,sizex,sizey> out;

		unsigned int i,j;
		for(i=0;i<sizex;i++) for(j=0;j<sizey;j++) out.data[i + j * sizex] = data[i + j * sizex] * fact[j];
		return(out);
	}

	LFHTEMP
	template<class O,Tuple_flag FO>
	TMatrix<C,sizex,sizey> TMatrix<C,sizex,sizey>::scale_cols(Tuple<O, sizex, FO> const &fact) const{
		TMatrix<C,sizex,sizey> out;

		unsigned int i,j;
		for(i=0;i<sizex;i++) for(j=0;j<sizey;j++) out.data[i + j * sizex] = data[i + j * sizex] * fact[i];
		return(out);
	}

	LFHTEMP Trianglix<C,sizey> TMatrix<C,sizex,sizey>::operator*(const Trianglix<C,sizex> &input)const{Trianglix<C,sizey> fout;		
		// Tr_OUT = (this) * tr_IN * (this)^T
		TMatrix<C,sizex,sizey> proj;
		unsigned int x,y,k;
		for(y=0;y< sizey;y++)
			for(x=0;x< sizex;x++){
				proj.data[x + y * sizex] = ExOp::mkmult(this->data[x + y * sizex],input.cell(x,x));
				for(k=0;k<x;k++) proj.data[x + y * sizex] += ExOp::mkmult(this->data[k + y * sizex],input.cell(k,x));
				for(k++;k<sizex;k++) proj.data[x + y * sizex] += ExOp::mkmult(this->data[k + y * sizex],ExOp::mktrju(input.cell(k,x)));
			}
		
		proj.show();
		
		for(y=0;y< sizey;y++){
			for(x=0;x<= y;x++){
				fout.cell(y,x) = proj.data[ y* sizex] * ExOp::mktrju(this->data[x * sizex]);
				for(k=1;k<sizex;k++) fout.cell(y,x) += proj.data[k + y* sizex] * ExOp::mktrju(this->data[k + x * sizex]);
			}
		}
		return fout;
	}

    LFHTEMP
    TMatrix<C,sizey,sizex> TMatrix<C,sizex,sizey>::inverse_2() const{ TMatrix<C,sizey,sizex> fout = this->transpose();
        bool dir = (sizex < sizey);
        TMatrix<C,sizex,sizex> HH_X; HH_X.toOne();
        TMatrix<C,sizex,sizex> HH_Y; HH_Y.toOne();
        unsigned int x=0;
        unsigned int y=0;
        unsigned int z;
     //   if (dir) 
            
     //   }
        double tmp;
        C* hh = new C[ (sizex < sizey) ? sizey:sizex];
            double sum2 =0.0f;
			for(z=x+1;z<sizex;z++) {
				hh[z] = fout.data[y+z * sizey];
				sum2 += ExOp::pnorm(hh[z]);
			}

			hh[x] = fout.data[y+x * sizey];
			if ((sum2 != 0.0f)||(fabs(ExOp::pnorm(hh[x]) / sum2) > 1.0f / ExCo<double>::epsilon())) {
			sum2 += ExOp::pnorm(hh[x]);
			tmp = sqrt(sum2);
			hh[x] += ((hh[x] < 0) ? -tmp: tmp);
			tmp = (sum2 + fabs(fout.data[y+x * sizey] *tmp));
            
            }
            y++; 
            dir = !dir;
        delete[](hh);

   //     HH_X.show();
    //    HH_Y.show();

        return fout;
    }

	LFHTEMP
	const TMatrix<C,sizex,sizey>& TMatrix<C,sizex,sizey>::LeftHouseHolderMultiply(const double * const vec, const double& sqrt_den, const int lenght, C* inner, bool hint){
			int inl = (hint) ? sizey + lenght - sizex : sizex; // hint tells that the first collumns are immune to transformation

			if ((sqrt_den == 0.0f)||(!ExCo<double>::isValid(sqrt_den))) return;
		
			double fact = -1.0f / sqrt_den;

			int i,j;
			i=0;

			for(j=0;j< inl;j++) inner[j] = data[ j + i* sizex] * vec[i];

			for(i++;i<lenght;i++){
				for(j=0;j< inl;j++) inner[j] += data[ j + i * sizex] * vec[i];
			}
			for(i=0;i<lenght;i++){
				for(j=0;j< inl;j++) data[j + i * sizex] += inner[j] * (vec[i] * fact);
			}

            return(*this);
		}

	LFHTEMP
	const TMatrix<C,sizex,sizey>& TMatrix<C,sizex,sizey>::leftHouseHolderMultiply(const double * const vec, const int lenght, const double& sqrt_den, bool hint){
			int inl = (hint) ? sizey + lenght - sizex : sizex; // hint tells that the first collumns are immune to transformation
            
			if ((sqrt_den == 0.0f)||(!ExCo<double>::isValid(sqrt_den))) return;
		
			C* min = data + sizex - inl +(sizey - lenght) * sizex;

			double fact = -1.0f / sqrt_den;

			int i,j;
			i=0;
			C* inner = new C[inl];

			for(j=0;j< inl;j++) inner[j] = min[ j + i* sizex] * vec[i];

			for(i++;i<lenght;i++){
				for(j=0;j< inl;j++) inner[j] += min[ j + i * sizex] * vec[i];
			}
			for(i=0;i<lenght;i++){
				for(j=0;j< inl;j++) min[j + i * sizex] += inner[j] * (vec[i] * fact);
			}

			delete[](inner);
            return(*this);
		}

		// hint implies that the outside the region of interes are zeroes only!
	LFHTEMP
	const TMatrix<C,sizex,sizey>& TMatrix<C,sizex,sizey>::rightHouseHolderMultiply(const double * const vec, const int lenght, const double& sqrt_den, bool hint){

			if ((sqrt_den == 0.0f)||(!ExCo<double>::isValid(sqrt_den))) return;
		
			int inl = (hint) ? sizex + lenght - sizey : sizey;
			int minj = sizey - inl;

			C* min = data + sizex - lenght + minj * sizex;

			double fact = -1.0f / sqrt_den;

			int i,j;
			i=0;
			C* inner = new C[inl];

			for(j=0;j< inl;j++) inner[j] = min[ i + j * sizex] * vec[i];

			for(i++;i<lenght;i++){
				for(j=0;j< inl;j++) inner[j] += min[ i +  j * sizex] * vec[i];
			}
			for(i=0;i<lenght;i++){
				for(j=0;j< inl;j++) min[i + j * sizex] += inner[j] * (vec[i] * fact);
			}

			delete[](inner);
            return(*this);
		}

/*

	LFHTEMP template<class O, class A> TMatrix<O,sizex,sizey> TMatrix<C,sizex,sizey>::operator+(A const & other) const{
		TMatrix<O,sizex,sizey> _out; int i;
		for(i=sizex*sizey-1;i!=ExCo<unsigned int>::maximum();i--) _out.data[i] = data[i] + other;
		return(_out);
	}
	LFHTEMP template<class O, class A> TMatrix<O,sizex,sizey> TMatrix<C,sizex,sizey>::operator-(A const & other) const{
		TMatrix<O,sizex,sizey> _out;int i;
		for(i=sizex*sizey-1;i!=ExCo<unsigned int>::maximum();i--) _out.data[i] = data[i] - other;
		return(_out);
	}
	LFHTEMP template<class O, class A> TMatrix<O,sizex,sizey> TMatrix<C,sizex,sizey>::operator*(A const & other) const{
		TMatrix<O,sizex,sizey> _out;int i;
		for(i=sizex*sizey-1;i!=ExCo<unsigned int>::maximum();i--) _out.data[i] = data[i] * other;
		return(_out);
	}

	LFHTEMP template<class O, class A> TMatrix<O,sizex,sizey> TMatrix<C,sizex,sizey>::operator/(A const & other) const{
		TMatrix<O,sizex,sizey> _out;int i;
		for(i=sizex*sizey-1;i!=ExCo<unsigned int>::maximum();i--) _out.data[i] = data[i] / other;
		return(_out);
	}*/




	LFHTEMP
	template<class O>
	TMatrix<C,sizex,sizey>& TMatrix<C,sizex,sizey>::operator*=(TMatrix<O,sizex,sizex> const & other){
		C tmpr[sizex];
		unsigned int i,j,k;

		for(j = 0;j<sizey;j++){
		for(i=0;i<sizex;i++) {
			tmpr[i] = data[i + j*sizex];
			data[i + j*sizex] =tmpr[0] * other.data[i];
		}
		for(k=1;k<sizex;k++){
		for(i=0;i<sizex;i++) data[i + j*sizex] +=tmpr[k] * other.data[i+ k*sizex];
		}
		}
		return(*this);
		}
	LFHTEMP
	template<class O>
	TMatrix<C,sizex,sizey>& TMatrix<C,sizex,sizey>::operator/=(TMatrix<O,sizex,sizex> const & other) {

		return(*this);
		}

	LFHTEMP
	template<class O, unsigned int osize>
	TMatrix<C,osize,sizey> TMatrix<C,sizex,sizey>::operator*(TMatrix<O,osize,sizex> const & other) const{
		TMatrix<C,osize,sizey> _out;
		unsigned int i,j,k;
		for(k=0;k<osize;k++){
			for(j=0;j<sizey;j++){
				_out.data[k + j * osize] = data[j * sizex] * other.data[k];
				for(i=1;i<sizex;i++){
					_out.data[k + j * osize] += data[i + j * sizex] * other.data[k + i * osize];
				}
			}
		}
		return(_out);
	}

	LFHTEMP
	template<class O>
	Tuple<C,sizey> TMatrix<C,sizex,sizey>::operator*(Tuple<O,sizex> const & other) const{
		Tuple<C,sizey> _out;
		unsigned int i,j;
		for(j=0;j<sizey;j++){
			_out.data[j] = data[j * sizex] * other.data[0];
			for(i=1;i<sizex;i++){
				_out.data[j] += data[i + j * sizex] * other.data[i];
			}
		}
		return(_out);
		}


	LFHTEMP
	Tuple<C,sizey> TMatrix<C,sizex,sizey>::operator*(Tuple<C,sizex> const & other) const{
		Tuple<C,sizey> _out;
		unsigned int i,j;
		for(j=0;j<sizey;j++){
			_out.data[j] = data[j * sizex] * other.data[0];
			for(i=1;i<sizex;i++){
				_out.data[j] += data[i + j * sizex] * other.data[i];
			}
		}
		return(_out);
		}
	LFHTEMP
	Continuous<C,sizex,C,sizey>* TMatrix<C,sizex,sizey>::derivative(Tuple<C,sizex> const & where) const{
		return(new TMatrix<C,sizex,sizey>(*this));
		}


//	LFHTEMP
//	PolyThing<C,0> TMatrix<C,sizex,sizey>::makeCharacteristic() const{
		/*
		PolyThing<C,0> _out;
		int i;
		_out.setOrder(sizex);

		_out[0] = data[0];
		for(i=1;i<sizex;i++) _out[0] *= data[i * (sizex+1)];
	*/
//	}


	template<class C, int size>
	Matrianglix<C,size>::Matrianglix(const TMatrix<C,size,size> & what){
		unsigned int i,j,k;
		for(i=0;i<size*size;i++){
			data[i] = what.data[i];
			}

		for(i=0;i<size-1;i++){
			for(j=i+1;j<size;j++){
				data[i + j*size] /= data[i*(size+1)];
				for(k=i+1;k<size;k++){
					data[k + j*size] -= data[i + j*size] * data[k + i*size];
				}
			}
		}
/*
		for(j=1;j<size;j++){



		for(i=0;i<j;i++){
			data[i + j*size] = what.data[i+ j*size];
			for(k=0;k<i;k++) data[i + j*size] -= data[k + j*size] * what.data[i+ k*size];
			data[i + j*size] = data[i + j*size] / data[i + (size+1)];
		}
		for(;i<size;i++){
			data[i + j*size] = what.data[i+ j*size];
			for(k=0;k<j;k++) data[i + j*size] -= data[k + j*size] * what.data[i+ k*size];
			}
		}

	}*/


	}


	template<class C, int size>
	C Matrianglix<C,size>::determinant(){
		C _out = data[0];
		for(unsigned int i=1;i<size;i++) _out *= data[i * (size+1)];
		return(_out);
	}

	template<class C, int size>
	Matrianglix<C,size> Matrianglix<C,size>::operator*(const Matrianglix<C, size> & other){
		return(Matrianglix<C,size>(( ((TMatrix<C,size,size>) (*this)) * ((TMatrix<C,size,size>) (other)) )  ));
	}

	template<class C, int size>
	TMatrix<C,size,size> Matrianglix<C,size>::inverse(){
		int i,j,k;

		Matrianglix<C,size> _tmp;

		for(int i=0;i<size;i++) {
			_tmp.data[i* (size+1)] = ExCo<C>::invert(data[i* (size+1)]); // inverse of D

		}

		// back stubstitution!
		// lower
		for(i=1;i<size;i++){
			for(j=0;j<size-i;j++){
				_tmp.data[i*size+j*(size+1)] = -data[i*size+j*(size+1)];
				_tmp.data[i+j*(size+1)] = -data[i+j*(size+1)];
				for(k=1;k<i;k++) {
					_tmp.data[i*size+j*(size+1)] -= data[k+i*size+j*(size+1)] * _tmp.data[(i-k)*size+j*(size+1)];
					_tmp.data[i+j*(size+1)] -= data[k*size+i+j*(size+1)] * _tmp.data[(i-k)+j*(size+1)];
				}
				_tmp.data[i+j*(size+1)] *= _tmp.data[(j+i)* (size+1)]*_tmp.data[j*(size+1)];
			}

		}
		_tmp.perm = -perm;

	//	_tmp.show(stdout);
		TMatrix<C,size,size> _out;


		for(i=0;i<size;i++){
			_out.data[i * (size+1)] = _tmp.data[i* (size+1)];
			for(k=i+1;k<size;k++) _out.data[i * (size+1)] += _tmp.data[k + i* size] * _tmp.data[i + k* size];
			for(j=i+1;j<size;j++){
				_out.data[i+j  * size] = _tmp.data[i + j* size] * _tmp.data[j + j* size];
				_out.data[j+i  * size] = _tmp.data[j + i* size];
				for(k=j+1;k<size;k++) _out.data[i+j  * size] += _tmp.data[i + k* size] * _tmp.data[k + j* size];
				for(k=j+1;k<size;k++) _out.data[j+i  * size] += _tmp.data[k + i* size] * _tmp.data[j + k* size];
			}
		}




		return(_out);
	}




	template<class C, int size>
	Matrianglix<C,size>::operator TMatrix<C,size,size>() const{
		TMatrix<C,size,size> _out;
		int i,j,k;

		for(i=0;i<size;i++){
			_out.data[perm[i] + i*size] = data[i * (size+1)];
			for(k=1;k<=i;k++) _out.data[perm[i] + i*size] += data[i * (size+1)-k ] * data[i * (size+1)-k*size];

			for(j= i+1;j<size; j++){
				_out.data[perm[i] + j*size] = data[i + j * size] * data[i * (size+1)];
				for(k=1;k<=i;k++) _out.data[perm[i] + j*size]  += data[i-k + j * size] * data[j + (i-k) * size];
				_out.data[perm[j] + i*size] = data[j + i * size];
				for(k=1;k<=i;k++) _out.data[perm[j] + i*size]  += data[j + (i-k) * size] * data[(i-k) + j * size];
			}
		}


		return(_out);
	}


	template<class C, int size>
	void Matrianglix<C,size>::show(FILE* o) const{
		int i,j;
		for(j=0;j<size;j++){
			for(i=0;i<size;i++){
				ExOp::show(data[i + j*size],o,1);
				fprintf(o,"%c", i == size-1 ? '\n' : '\t');
			}
		}
	}
template<class C, int size> string  Matrianglix<C,size>::type_tostring() const{char buffer[256]; sprintf(buffer,"%u",size); return string("Matrianglix<") + ExOp::type_tostring(data[0]) +string(",")+ ExOp::type_tostring(buffer) + string(">");}



	/*
	template<class C, int size>
	const Matrix<C,size,size>& Matrianglix<C,size>::power(double exp){

		Matrix<C,size,size> _out;
		int i,j,k;
		C mdiag[size];
		for(j=0;j<size;j++){
			mdiag[j] = (C) pow((double)diag[size],exp);
		}
		for(j=0;j<size;j++){
			for(i=0;i<size;i++){
				_out.data[i + j*size] = mdiag[0] * (perm.data[0 + j*size] * iperm.data[i + 0*size]).r; // imagenary elems should cancel out! ^^
				for(k=1;k<size;k++){
					_out.data[i + j*size] += mdiag[k] * (perm.data[k + j*size] * iperm.data[i + k*size]).r;
				}
			}
		}

		return(_out);
	}*/


#undef LFHTEMP
#define LFHTEMP template <class C>


LFHTEMP	void Complex<C>::zero() {ExOp::toZero((*this)[0]);ExOp::toZero((*this)[1]);}
LFHTEMP	void Complex<C>::one() {ExOp::toOne((*this)[0]);ExOp::toZero((*this)[1]);}
LFHTEMP	void Complex<C>::random() {ExOp::toRand((*this)[0]);ExOp::toRand((*this)[1]);}


	LFHTEMP	C& Complex<C>::operator[](unsigned int ind) {return data[ind];}
//	LFHTEMP	C Complex<C>::operator[](unsigned int ind) const {return data[ind];}
	LFHTEMP	const C& Complex<C>::operator[](unsigned int ind) const {return data[ind];}


	LFHTEMP	setcomparison Complex<C>::setcmp(const Complex<C> &a) const{
		bool v = data[0] < ExOp::mkzero<C>();
		if (v ^ (a.data[0] < ExOp::mkzero<C>())) return v ? SETCMP_LT : SETCMP_GT;
		else {
			double norm = ExOp::pnorm(data[0]) + ExOp::pnorm(data[1]);
			double onorm = ExOp::pnorm(a.data[0]) + ExOp::pnorm(a.data[1]);
			if (norm == onorm){
				return ExOp::setcmp(data[1], a.data[1]);
			} else return ((v)^(norm > onorm)) ? SETCMP_LT : SETCMP_GT;
		}
	}

	LFHTEMP bool Complex<C>::operator>(const Complex<C> &a)const{return((setcmp(a) & SETCMP_CMP_MASK2) == SETCMP_MASKED_GT);}
	LFHTEMP bool Complex<C>::operator>=(const Complex<C> &a)const{return((setcmp(a) & SETCMP_CMP_MASK) == SETCMP_MASKED_GE);}
	LFHTEMP bool Complex<C>::operator<(const Complex<C> &a)const{return((setcmp(a) & SETCMP_CMP_MASK2) == SETCMP_MASKED_LT);}
	LFHTEMP bool Complex<C>::operator<=(const Complex<C> &a)const{return((setcmp(a) & SETCMP_CMP_MASK) == SETCMP_MASKED_LE);}
	LFHTEMP bool Complex<C>::operator==(const Complex<C> &o)const{return(data[0] == o.data[0]) && (data[1] == o.data[1]);}
	LFHTEMP bool Complex<C>::operator!=(const Complex<C> &o)const{return(data[0] != o.data[0]) || (data[1] != o.data[1]);}




LFHTEMP	template<class A>
Complex<C>& Complex<C>::operator+=(const Complex<A>& other){
	ExOp::toadd((*this)[0],other[0]);
	ExOp::toadd((*this)[1],other[1]);
	return(*this);
	}

LFHTEMP	template<class A>
	Complex<C>& Complex<C>::operator-=(const Complex<A>& other){
		ExOp::tosub((*this)[0],other[0]);
		ExOp::tosub((*this)[1],other[1]);
		return(*this);
	}

LFHTEMP	template<class A>
Complex<C>& Complex<C>::operator*=(const Complex<A>& other){
	C tmp = ExOp::mkmult((*this)[0], other[1]);
	ExOp::tomult((*this)[0], other[0]);
//	ExOp::show((*this)[0]);
	ExOp::tosub((*this)[0], ExOp::mkmult((*this)[1],other[1]));
//	ExOp::show((*this)[0]);
	ExOp::tomult((*this)[1], other[0]);
//	ExOp::show((*this)[1]);
	ExOp::toadd((*this)[1], tmp);
	return(*this);
	}

LFHTEMP	template<class A> Complex<C>& Complex<C>::operator+=(const A& a){(*this)[0] += a;return(*this);}
LFHTEMP	template<class A> Complex<C>& Complex<C>::operator-=(const A& a){(*this)[0] -= a; return(*this);}
LFHTEMP	template<class A> Complex<C>& Complex<C>::operator*=(const A& a){
	(*this)[0] *= a;
	(*this)[1] *= a;
	return(*this);
}
LFHTEMP	template<class A> Complex<C>& Complex<C>::operator/=(const A& a){
	(*this)[0] /= a;
	(*this)[1] /= a;
	return(*this);
}

LFHTEMP	template<class A>
Complex< typename STDRETTYPE2<C,A>::PLUS_TYPE > Complex<C>::operator+(const Complex<A>& other) const{
	return(Complex((*this)[0] + other[0], (*this)[1] + other[1]));
}

LFHTEMP	template<class A>
Complex< typename STDRETTYPE2<C,A>::MINU_TYPE > Complex<C>::operator-(const Complex<A>& other) const{
	return(Complex((*this)[0] - other[0], (*this)[1] - other[1]));
}

LFHTEMP	template<class A> Complex< typename STDRETTYPE2<C,A>::PLUS_TYPE > Complex<C>::operator*(Complex<A> const & other) const{
    return Complex< typename STDRETTYPE2<C,A>::PROD_TYPE >( ((*this)[0] *  other[0])-  ((*this)[1] *  other[1]) , ((*this)[0] *  other[1])+  ((*this)[1] *  other[0]));
    }
LFHTEMP	template<class A> Complex< typename STDRETTYPE2<C,A>::DIVI_TYPE > Complex<C>::operator/(Complex<A> const & other) const{}

LFHTEMP Complex< typename ExCo<C>::NEG_TYPE > Complex<C>::operator-() const{
	return(Complex(-(*this)[0], -(*this)[1]));
}

LFHTEMP double Complex<C>::pnorm() const{ return(ExOp::pnorm(((Tuple<C,2>*)this)->data[0]) + ExOp::pnorm(((Tuple<C,2>*)this)->data[1])); }

LFHTEMP Complex<C> Complex<C>::inverse() const{
	Complex<C> _out;
	C tmpA = ((*this)[0] * (*this)[0]) + ((*this)[1] * (*this)[1]);
	_out[0] =(*this)[0] / tmpA;
	_out[1] =(*this)[1] / tmpA;
	return(Complex<C>(_out));
	}

LFHTEMP void Complex<C>::show(FILE* f, int level) const{
	switch(level){
		case 1:
            ExOp::show((*this)[0],f,level+1);fprintf(f," r\t");ExOp::show((*this)[1],f,level+1);fprintf(f," i");
		break;
		case 0:
            ExOp::show((*this)[0],f,level+1);fprintf(f," r\t");ExOp::show((*this)[1],f,level+1);fprintf(f," i\n");
		break;
		case 2:
            fprintf(f,"[");ExOp::show((*this)[0],f,level+1);fprintf(f," r;");ExOp::show((*this)[1],f,level+1);fprintf(f," i]");
			break;
		default:
            fprintf(f,"(");ExOp::show((*this)[0],f,level+1);fprintf(f," r,");ExOp::show((*this)[1],f,level+1);fprintf(f," i)");
			break;
	}
}
LFHTEMP string Complex<C>::type_tostring() const{return string("Complex<") + ExOp::type_tostring(data[0]) + string(">");}

LFHTEMP	Quaternion<C>::Quaternion() {}

LFHTEMP	const Quaternion<C>& Quaternion<C>::to_normal(const C& _i,const C& _j,const C& _k){
    return(*this);
}

LFHTEMP	const Quaternion<C>& Quaternion<C>::to_normal_and_scale(const C& _i,const C& _j,const C& _k,const C& _s){
    data[0] = sqrt(_s);
    data[1] = 0.0f;
    data[2] = 0.0f;
    data[3] = 0.0f;
    return(*this);
}

LFHTEMP	double Quaternion<C>::pnorm() const{return ExOp::pnorm(data[0]) + ExOp::pnorm(data[1]) +ExOp::pnorm(data[2]) + ExOp::pnorm(data[3]);}

LFHTEMP	const Quaternion<C>& Quaternion<C>::rotateX(double ang){
    double s = sin(ang); double c = cos(ang);
    C tmp = data[0] * c + (data[1] * s);
    data[1] =  data[1] * c + (data[0] * (-s)); data[0] = tmp;
    tmp = data[2] * c + (data[3] * (-s));
    data[3] = data[3] * c + (data[2] * s); data[2] = tmp;
    return(*this);
}
LFHTEMP	const Quaternion<C>& Quaternion<C>::rotateY(double ang){
    double s = sin(ang); double c = cos(ang);
    return(*this);
}
LFHTEMP	const Quaternion<C>& Quaternion<C>::rotateZ(double ang){
    double s = sin(ang); double c = cos(ang);
    return(*this);
}

LFHTEMP	void Quaternion<C>::mk_proj_matrix(TMatrix<C,3,3> &fout)const{
    fout.data[0] = data[0] * data[0];
    C tmp = data[1] * data[1];
    fout.data[4] = fout.data[0] - tmp;
    fout.data[0] += tmp;
    tmp = data[2] * data[2];
    fout.data[8] = fout.data[4] - tmp;
    fout.data[4] += tmp;
    fout.data[0] -= tmp;
    tmp = data[3] * data[3];
    fout.data[8] += tmp;
    fout.data[4] -= tmp;
    fout.data[0] -= tmp;

    fout.data[1] = data[0] * data[3] *2.0f;
    fout.data[5] = data[0] * data[1]*2.0f;
    fout.data[6] = data[0] * data[2]*2.0f;
    tmp = data[1] * data[2]*2.0f;
    fout.data[3] = tmp - fout.data[1];
    fout.data[1] += tmp;


    tmp = data[3] * data[2]*2.0f;
    fout.data[7] = tmp - fout.data[5];
    fout.data[5] += tmp;

    tmp = data[1] * data[3]*2.0f;
    fout.data[2] = tmp - fout.data[6];
    fout.data[6] += tmp;

}
LFHTEMP	void Quaternion<C>::mk_proj_matrix(TMatrix<C,4,4> &fout)const{
    fout.data[0] = data[0] * data[0];
    C tmp = data[1] * data[1];
    fout.data[5] = fout.data[0] - tmp;
    fout.data[0] += tmp;
    tmp = data[2] * data[2];
    fout.data[10] = fout.data[5] - tmp;
    fout.data[5] += tmp;
    fout.data[0] -= tmp;
    tmp = data[3] * data[3];
    fout.data[10] += tmp;
    fout.data[5] -= tmp;
    fout.data[0] -= tmp;

    fout.data[1] = data[0] * data[3]*2.0f;
    fout.data[6] = data[0] * data[1]*2.0f;
    fout.data[8] = data[0] * data[2]*2.0f;
    tmp = data[1] * data[2]*2.0f;
    fout.data[4] = tmp - fout.data[1];
    fout.data[1] += tmp;

    tmp = data[3] * data[2]*2.0f;
    fout.data[9] = tmp - fout.data[6];
    fout.data[6] += tmp;

    tmp = data[1] * data[3]*2.0f;
    fout.data[2] = tmp - fout.data[8];
    fout.data[8] += tmp;

    ExOp::toZero(fout.data[3]);
    ExOp::toZero(fout.data[7]);
    ExOp::toZero(fout.data[11]);
    ExOp::toZero(fout.data[12]);
    ExOp::toZero(fout.data[13]);
    ExOp::toZero(fout.data[14]);
    ExOp::toOne(fout.data[15]);
}
LFHTEMP	Tuple<C,3> Quaternion<C>::mkXvector() const{ Tuple<C,3> fout;

    fout[0] = data[0] * data[0] + data[1] * data[1] - (data[2] * data[2] + data[3] * data[3]);
    fout[1] = ( data[1] * data[2] - data[0] * data[3]) *2.0f;
    fout[2] = ( data[1] * data[3] + data[0] * data[2]) *2.0f;

    return(fout);
}
LFHTEMP	Tuple<C,3> Quaternion<C>::mkYvector() const{Tuple<C,3> fout;
    fout[1] = data[0] * data[0] + data[2] * data[2] - (data[1] * data[1] + data[3] * data[3]);
    fout[2] = ( data[2] * data[3] - data[0] * data[1]) *2.0f;
    fout[0] = ( data[2] * data[1] + data[0] * data[3]) *2.0f;
    return(fout);
}
LFHTEMP	Tuple<C,3> Quaternion<C>::mkZvector() const{Tuple<C,3> fout;
    fout[2] = data[0] * data[0] + data[3] * data[3] - (data[2] * data[2] + data[1] * data[1]);
    fout[0] = ( data[3] * data[1] - data[0] * data[2]) *2.0f;
    fout[1] = ( data[3] * data[2] + data[0] * data[1]) *2.0f;
    return(fout);
}
LFHTEMP	template<class O> void Quaternion<C>::wrXvector(O* fout) const{
    fout[0] = data[0] * data[0] + data[1] * data[1] - (data[2] * data[2] + data[3] * data[3]);
    fout[1] = ( data[1] * data[2] - data[0] * data[3]) *2.0f;
    fout[2] = ( data[1] * data[3] + data[0] * data[2]) *2.0f;
}
LFHTEMP	template<class O> void Quaternion<C>::wrYvector(O* fout) const{
    fout[1] = data[0] * data[0] + data[2] * data[2] - (data[1] * data[1] + data[3] * data[3]);
    fout[2] = ( data[2] * data[3] - data[0] * data[1]) *2.0f;
    fout[0] = ( data[2] * data[1] + data[0] * data[3]) *2.0f;
}
LFHTEMP	template<class O> void Quaternion<C>::wrZvector(O* fout) const{
    fout[2] = data[0] * data[0] + data[3] * data[3] - (data[2] * data[2] + data[1] * data[1]);
    fout[0] = ( data[3] * data[1] - data[0] * data[2]) *2.0f;
    fout[1] = ( data[3] * data[2] + data[0] * data[1]) *2.0f;
}

LFHTEMP    const Quaternion<C>& Quaternion<C>::toUnitQuaternion(const Tuple<double, 3> value){
    double norm = sqrt(value[0] * value[0] + value[1] *value[1] + value[2]*value[2]);
    ExOp::toOne(data[1]); 
    ExOp::toOne(data[2]);
    ExOp::toOne(data[3]);
    
    double ang = norm * M_PI * 0.5f;
    ExOp::toOne(data[0]); data[0] *= cos(ang);
    norm = sin(ang) / norm;
    data[1] *= value[0] * norm;
    data[2] *= value[1] * norm;
    data[3] *= value[2] * norm;
    return(*this);
}
    
// W2 + X2 -Y2 -Z2   2WZ + 2XY     2WY + 2XZ
// 2WZ - 2XY     W2 + Y2 - X2 -Z2   2WX + 2YZ 
// 2WY + 2XZ         2WX - 2YZ   W2 + Z2 - Y2 -X2 

// COS2 + 2SIN2X2 - 1   2SIN2WZ + 2SIN2XY     2WY + 2XZ
// 2COSSINZ - 2SIN2XY     COS2 + 2SIN2Y2 - 1   2WX + 2YZ 
// 2COSSINY + 2SIN2XZ         2WX - 2YZ   COS2 + 2SIN2Z2 - 1

// COS2 + 2SIN2X2 - 1   2SIN2WZ + 2SIN2XY     2WY + 2XZ
// 2SINZ + (-1 + COS)XY     COS2 + 2SIN2Y2 - 1   2WX + 2YZ 
// 2SINY + 2SIN2XZ         2WX - 2YZ   COS2 + 2SIN2Z2 - 1

LFHTEMP	const C&  Quaternion<C>::operator[](unsigned int w) const{return(data[w]);}
LFHTEMP	C& Quaternion<C>::operator[](unsigned int w){return(data[w]);}

LFHTEMP	void Quaternion<C>::zero() {ExOp::toZero((*this)[0]);ExOp::toZero((*this)[1]);ExOp::toZero((*this)[2]);ExOp::toZero((*this)[3]);}
LFHTEMP	void Quaternion<C>::one() {ExOp::toOne((*this)[0]);ExOp::toZero((*this)[1]);ExOp::toZero((*this)[2]);ExOp::toZero((*this)[3]);}
LFHTEMP	void Quaternion<C>::random() {ExOp::toRand((*this)[0]);ExOp::toRand((*this)[1]);ExOp::toRand((*this)[2]);ExOp::toRand((*this)[3]);}

LFHTEMP const Quaternion<C>& Quaternion<C>::inverse() const {
	Quaternion<C> _out;
	C tmpA = ((*this)[0] * (*this)[0]) + ((*this)[1] * (*this)[1])+ ((*this)[2] * (*this)[2])+ ((*this)[3] * (*this)[3]);
	_out[0] =(*this)[0] / tmpA;
	tmpA = -tmpA;
	_out[1] =(*this)[1] / tmpA;
	_out[2] =(*this)[2] / tmpA;
	_out[3] =(*this)[3] / tmpA;
	return(Complex<C>(_out));
}

LFHTEMP	template<class A>
Quaternion<C>& Quaternion<C>::operator+=(Quaternion<A> const & o){
		ExOp::toadd(data[0],o[0]);
		ExOp::toadd(data[1],o[1]);
		ExOp::toadd(data[2],o[2]);
		ExOp::toadd(data[3],o[3]);
		return(*this);
	}

LFHTEMP	template<class A>
Quaternion<C>& Quaternion<C>::operator-=(Quaternion<A> const & o){
		ExOp::tosub(data[0],o[0]);
		ExOp::tosub(data[1],o[1]);
		ExOp::tosub(data[2],o[2]);
		ExOp::tosub(data[3],o[3]);
		return(*this);
	}

LFHTEMP void Quaternion<C>::show(FILE* f, int level) const{
    	switch(level){
		case 1:
		case 0:
            ExOp::show((*this)[0],f,level+1);fprintf(f," r\t");ExOp::show((*this)[1],f,level+1);fprintf(f," i\t");ExOp::show((*this)[2],f,level+1);fprintf(f," j\t");ExOp::show((*this)[3],f,level+1);fprintf(f," k\n");
		break;
		case 2:
            fprintf(f,"[");ExOp::show((*this)[0],f,level+1);fprintf(f," r;");ExOp::show((*this)[1],f,level+1);fprintf(f," i;");ExOp::show((*this)[2],f,level+1);fprintf(f," j;");ExOp::show((*this)[3],f,level+1);fprintf(f," k]");
			break;
		default:
           fprintf(f,"(");ExOp::show((*this)[0],f,level+1);fprintf(f," r,");ExOp::show((*this)[1],f,level+1);fprintf(f," i,");ExOp::show((*this)[2],f,level+1);fprintf(f," j,");ExOp::show((*this)[3],f,level+1);fprintf(f," k)");
			break;
	}
	}
LFHTEMP string Quaternion<C>::type_tostring() const{return string("Quaternion<") + string(">");}

LFHTEMP template<unsigned int SIZE> Trianglix<C, 0u>::Trianglix(const Tuple<C, SIZE>& other): t_size(other.size()) {
    data = new C[totsize()];
    unsigned int i,j,k;
    for(j=0,k=0;j<t_size;j++) for(i=0;i<=j;i++,k++) data[k] = ExOp::mktrju(other[i]) * other[j];     
    }

LFHTEMP Trianglix<C, 0u>::Trianglix(const Tuple<C, 0u>& other): t_size(other.getSize()) {
    data = new C[totsize()];
    unsigned int i,j,k;
    for(j=0,k=0;j<t_size;j++) for(i=0;i<=j;i++,k++) data[k] = ExOp::mktrju(other[i]) * other[j];     
    }

LFHTEMP Trianglix<C, 0u>::Trianglix(const Trianglix<C, 0u>& other): t_size(other.getSize()) {data = (t_size) ? new C[totsize()] : NULL; LFH_USE_MEMLEAK_HELP(if (data) memleak_helper[(unsigned int)data] =2;) unsigned int i,ts; for(i=0,ts=totsize();i<ts;i++) data[i] = other.data[i];}
LFHTEMP template<unsigned int SIZE> Trianglix<C, 0u>::Trianglix(const Trianglix<C, SIZE>& other): t_size(other.getSize()) {data = (t_size) ? new C[totsize()] : NULL; LFH_USE_MEMLEAK_HELP(if (data) memleak_helper[(unsigned int)data] =2;) unsigned int i,ts; for(i=0,ts=totsize();i<ts;i++) data[i] = other.data[i];}


LFHTEMP template<unsigned int SIZE> Trianglix<C, 0u>& Trianglix<C, 0u>::operator=(const Tuple<C, SIZE>& other){
    setSize(SIZE);
    unsigned int i,j,k;
    for(j=0,k=0;j<t_size;j++) for(i=0;i<=j;i++,k++) data[k] = ExOp::mktrju(other[i]) * other[j];     
    return(*this);
}

LFHTEMP Trianglix<C, 0u>& Trianglix<C, 0u>::memmove(Trianglix<C, 0u>& other){
	if (t_size) delete[](data);
	data = other.data;
	t_size = other.t_size;
	other.t_size = 0;
    return(*this);
}

LFHTEMP Trianglix<C, 0u>::~Trianglix(){if (t_size) {delete[](data);LFH_USE_MEMLEAK_HELP(memleak_helper.erase((unsigned int)data);)}}
LFHTEMP void Trianglix<C, 0u>::setSize(unsigned int s){if (s == t_size) return; if (t_size) {delete[](data); LFH_USE_MEMLEAK_HELP(memleak_helper.erase((unsigned int)data);)} t_size =s; data = (t_size == 0) ? NULL : new C[totsize()]; LFH_USE_MEMLEAK_HELP(if (data) memleak_helper[(unsigned int)data] =2;)}



LFHTEMP Matrix<C> Trianglix<C, 0u>::makeMatrix()const{
    Matrix<C> damat;damat.setSizes(t_size,t_size);
    unsigned int i,j,k;
    for(j=0,k=0;j<t_size;j++) {
        for(i=0;i<j;i++,k++) {damat(i,j) = data[k];damat(j,i) = data[k];}
        damat(j,j) = data[k]; k++;
    } 
    return(damat);
    }



LFHTEMP void Trianglix<C, 0u>::HouseHolderMultiply(const C * const vec, double denum2, unsigned int length, C* buf, bool hint ){
    // house holder multiplication, done twice!
    // = T - vec * vec' T - T' vec * vec' + vec * vec' * T *vec * vec'
    unsigned int i,j,k;
    C s;
        if ((denum2 == 0.0f)||(!ExCo<double>::isValid(denum2))) return;
        buf[0] = data[0] * vec[0];
        s = data[0] * ExOp::mktrju(vec[0]) * vec[0];
	
	for(j=1,k=1;j<length;j++,k++) {
        buf[j] = data[k] * vec[0]; buf[0] += ExOp::mktrju(data[k++]) * vec[j];
        for(i=1;i<j;i++) {buf[j] += data[k] * vec[i]; buf[i] += ExOp::mktrju(data[k++]) * vec[j];}
        s += data[k] * ExOp::mktrju(vec[j]) * vec[j] + (buf[j] *  ExOp::mktrju(vec[j]) + ExOp::mktrju(buf[j]) * vec[j]);
        buf[j] += data[k] * vec[j];
	}

        if (!hint){
        for(;j<t_size;j++) {
			buf[j] = data[k++] * vec[0];
			for(i=1;i<length;i++) buf[j] += data[k++] * vec[i];
            k+= j - length+1;
        }
        for(j=0;j<t_size;j++) buf[j] /= denum2;
        }else for(j=0;j<length;j++) buf[j] /= denum2;
        s /= denum2 * denum2; 

        for(j=0,k=0;j<length;j++) {
			for(i=0;i<j;i++) data[k++] += s * vec[j] *ExOp::mktrju(vec[i]) - buf[j]*ExOp::mktrju(vec[i]) - vec[j] * ExOp::mktrju(buf[i]) ;
			data[k++] += s * vec[j] *ExOp::mktrju(vec[j]) - buf[j]*ExOp::mktrju(vec[j]) - vec[j] * ExOp::mktrju(buf[j]);
        }
        if (hint) return;
        for(;j<t_size;j++) {
            for(i=0;i<length;i++) data[k++] -= buf[j]*ExOp::mktrju(vec[i]);
            k+= j - length+1;
        }
}
	
	


LFHTEMP void Trianglix<C, 0u>::offdiagelimination_down(const C &fact,unsigned int col, C* buf){ // multiply row/col "col" by  fact, and adds it in "(col+1)" 
    // = T - vec * vec' T - T' vec * vec' + vec * vec' * T *vec * vec'
    unsigned int i,j,k;

    k = (col *(col+1))/2;
    for(i=0;i<=col;i++) {buf[i] = fact * data[k+i];}
    for(;i<t_size;i++) {k += i; buf[i] = fact * data[k+col];}
    k = ((col+2) *(col+1))/2;
    for(i=0;i<=col;i++) data[k+i] += buf[i];
    k++;
    data[k+col] += buf[col] * ExOp::mktrju(fact) + ExOp::mktrju(buf[i]);
    for(i++;i<t_size;i++) {k += i; data[k+col] += buf[i];}
}

LFHTEMP void Trianglix<C, 0u>::offdiagelimination_up(const C &fact,unsigned int col, C* buf){ // multiply row/col "col" by  fact, and adds it in "(col+1)" 
    // = T - vec * vec' T - T' vec * vec' + vec * vec' * T *vec * vec'
    unsigned int i,j,k;

    k = ((col+2) *(col+1))/2;
    for(i=0;i<=col;i++) buf[i] = fact * data[k+i];
    buf[i] = fact * data[k+i]; 
    k++;
    for(i++;i<t_size;i++) {k += i; buf[i] = fact * data[k+col];}

    k = ((col) *(col+1))/2;
    for(i=0;i<col;i++) data[k+i] += buf[i];
    data[k+col] += ExOp::mktrju(buf[col+1]) * fact + ExOp::mktrju(buf[i]); 
    for(i++;i<t_size;i++) {k += i; data[k+col] += buf[i];}
}

LFHTEMP void Trianglix<C, 0u>::offdiagelimination_up_backroutine(const C &fact,unsigned int col, C* buf){ // multiply row/col "col" by  fact, and adds it in "(col+1)" 
    // = T - vec * vec' T - T' vec * vec' + vec * vec' * T *vec * vec'
    unsigned int i,j,k;

    k = ((col+2) *(col+1))/2;
    i =col+1;
    buf[i] = fact * data[k+i]; 
    k++;
    for(i++;i<t_size;i++) {k += i; buf[i] = fact * data[k+col];}

    k = ((col) *(col+1))/2;
    
    data[k+col] += ExOp::mktrju(buf[col+1]) * fact; 
    for(i=col+1;i<t_size;i++) {k += i; data[k+col] = buf[i];}
}



LFHTEMP void Trianglix<C, 0u>::QR_back(const C &factC,const C &factS,unsigned int col){
    C tmp;
    unsigned int k,i;
    k = (col *(col+1))/2;
    for(i=0;i<col;i++){
    tmp = data[k+i] * factC + data[k+i+col+1] * factS;
    data[k+i+col+1] = data[k+i+col+1] * factC - data[k+i] * factS;
    data[k+i] = tmp;

 //   printf("bepair = %i %i \n",k+i, k+i+col+1);
    }
    // A 

    k += i+col+1;
 //   tmp = data[k-col-1] * factC + data[k] * factS;
 //   C tmp2 = -data[k-col-1] * factS + data[k] * factC;
 //   C tmp3 = data[k] * factC + data[k+1] * factS;
 //   C tmp4 = -data[k] * factS + data[k+1] * factC;
 //   printf("dapox = %i %i %i \n",k-col-1,k, k+1);
//    data[k] = - tmp * factS + tmp3 * factC; // or tmp2 * factC + tmp4 * factS
    tmp = data[k-col-1] * factS* factS + data[k+1] * factC* factC - 2.0f * data[k] * factC * factS;
    data[k] = data[k] * (factC * factC - factS * factS) - (data[k-col-1] - data[k+1]) * factC * factS;
      
    data[k-col-1] += data[k+1] -tmp;
    data[k+1] = tmp;
   
    for(i=col+2;i<t_size;i++){k+=i;
    tmp = data[k] * factC + data[k+1] * factS;
    data[k+1] = data[k+1] * factC - data[k] * factS;
    data[k] = tmp;
 //   printf("afpair = %i %i \n",k, k+1);
    }
    
}
	
LFHTEMP template<unsigned int SIZE> Tuple<C,SIZE> Trianglix<C, 0u>::operator*( const Tuple<C, SIZE>& other) const{ Tuple<C,SIZE> fout;
    unsigned int i,j,k;
    fout[0] = other[0] * data[0];
	for(j=1,k=1;j<SIZE;j++) {
		fout[j] = other[0] * data[k]; 
		fout[0] += other[j] * ExOp::mktrju(data[k]); 
        for(i=1,k++;i<j;i++,k++) {fout[j] += other[i] * data[k]; fout[i] += other[j] * ExOp::mktrju(data[k]);}
        fout[j] += other[j] * data[k++];  
	}	
	return fout;
	}
	
LFHTEMP Tuple<C,0u> Trianglix<C, 0u>::operator*( const Tuple<C,0u>& other) const{ Tuple<C,0u> fout; fout.setSize(other.getSize());
    if (t_size == 0) return fout;
    unsigned int i,j,k;
    fout[0] = other[0] * data[0];
	for(j=1,k=1;j<t_size;j++) {
		fout[j] = other[0] * data[k]; 
		fout[0] += other[j] * ExOp::mktrju(data[k]); 
        for(i=1,k++;i<j;i++,k++) {fout[j] += other[i] * data[k]; fout[i] += other[j] * ExOp::mktrju(data[k]);}
        fout[j] += other[j] * data[k++];  
	}	
	return fout;
	}
	
LFHTEMP template<unsigned int SIZE> Tuple<C,SIZE> Trianglix<C, 0u>::divisionof( const Tuple<C, SIZE>& other) const{ Tuple<C,SIZE> fout;
    unsigned int i,j,k;
    fout[0] = other[0] * data[0];
	for(j=1,k=1;j<SIZE;j++) {
		fout[j] = other[0] * data[k]; 
		fout[0] += other[j] * ExOp::mktrju(data[k]); 
        for(i=1,k++;i<j;i++,k++) {fout[j] += other[i] * data[k]; fout[i] += other[j] * ExOp::mktrju(data[k]);}
        fout[j] += other[j] * data[k++];  
	}	
	return fout;
}

LFHTEMP Tuple<C,0u> Trianglix<C, 0u>::divisionof( const Tuple<C,0u>& other) const{ Tuple<C,0u> fout; fout.setSize(other.getSize());
    if (t_size == 0) return fout;
    unsigned int i,j,k;
    fout[0] = other[0] * data[0];
	for(j=1,k=1;j<t_size;j++) {
		fout[j] = other[0] * data[k]; 
		fout[0] += other[j] * ExOp::mktrju(data[k]); 
        for(i=1,k++;i<j;i++,k++) {fout[j] += other[i] * data[k]; fout[i] += other[j] * ExOp::mktrju(data[k]);}
        fout[j] += other[j] * data[k++];  
	}	
	return fout;
}


LFHTEMP template<unsigned int SIZE> C Trianglix<C, 0u>::Xformed_inner_product( const Tuple<C, SIZE>& other) const{
    C s;ExOp::toZero(s);
    C s2;
    unsigned int i,j,k;
        for(j=0,k=0;j<SIZE;j++) {
        if (j>0) s2 = data[k++] * other[0];
        else ExOp::toZero(s2);
        for(i=1;i<j;i++) s2 += data[k++] * other[i];
        // s += other[j] * ((ExOp::mktrju(s2) + s2) + data[k++] * other[j]); WORKS FOR REAL ONLY!
        s +=  data[k++] * ExOp::mktrju(other[j]) * other[j] + (ExOp::mkrealproj(s2) * ExOp::mkrealproj(other[j]) - ExOp::mkimmaproj(s2) * ExOp::mkimmaproj(other[j])) *2.0f;
        
        }
        
    return s;
    }

LFHTEMP C Trianglix<C, 0u>::Xformed_inner_product( const Tuple<C, 0u>& other) const{
    C s;
    C s2;
    unsigned int i,j,k;
    if (t_size == 0) {ExOp::toZero(s); return s;} 
	s = data[0] * other[0] * ExOp::mktrju(other[0]);
	
    if (ExOp::mktrju(other[0]) * other[0] < 0.0f) { fprintf(stderr, "Should not be negative!: "); ExOp::show(ExOp::mktrju(other[0]) * other[0],stderr); exit(1); }
	for(j=1,k=1;j<t_size;j++) {
        s2 = data[k++] * other[0];
        for(i=1;i<j;i++) s2 += data[k++] * other[i];
        s += 2.0f * other[j] * s2 + data[k++] * other[j]; // WORKS FOR REAL ONLY!
        // s +=  data[k++] * ExOp::mktrju(other[j]) * other[j] + (ExOp::mkrealproj(s2) * ExOp::mkrealproj(other[j]) - ExOp::mkimmaproj(s2) * ExOp::mkimmaproj(other[j])) *2.0f;
		if (ExOp::mktrju(other[j]) * other[j] < 0.0f) { fprintf(stderr,"Should not be negative!: "); ExOp::show(ExOp::mktrju(other[j]) * other[j],stderr); exit(1); }
	}
    return s;	
	}



	
	LFHTEMP C Trianglix<C, 0u>::trace_of_product(const Trianglix<C,0u> &other) const{C fout; ExOp::toZero(fout);
		unsigned int i,j,k;
		for(j=0,k=0;j<t_size;j++) {
			for(i=0;i<j;i++,k++) fout += data[k] * ExOp::mktrju(other.data[k]) + ExOp::mktrju(data[k]) * other.data[k];
			fout += data[k] * other.data[k];
			k++;
		}		
		return fout;
	}
	
LFHTEMP C Trianglix<C, 0u>::trace_of_division(const Trianglix<C,0u> &other) const{C fout; ExOp::toZero(fout);
	unsigned int i,j,k;
	for(j=0,k=0;j<t_size;j++) {
		for(i=0;i<j;i++,k++) fout += data[k] * ExOp::mktrju(other.data[k]) + ExOp::mktrju(data[k]) * other.data[k];
		fout += data[k] * other.data[k];
		k++;
	}		
	return fout;
}

	
LFHTEMP Trianglix<C, 0u> Trianglix<C, 0u>::operator*(const Matrix<C>& other) const{ // does not work with com
    Trianglix<C, 0u> fout; if (other.sizey != t_size) exit(1); fout.setSize(other.sizex);fout.toZero();
    C* s2 = new C[other.sizex];
    unsigned int i,j,k,l;
        for(j=0,k=0;j<t_size;j++) {

        if (j>0) {for(l=0;l<other.sizex;l++) s2[l] = data[k] * other.data[l]; k++;}
        else for(l=0;l<other.sizex;l++) ExOp::toZero(s2[l]);
        for(i=1;i<j;i++,k++) for(l=0;l<other.sizex;l++) s2[l] += data[k] * other.data[l + i * other.sizex];
        for(l=0;l<other.sizex;l++) s2[l] += ExOp::mktrju(s2[l]) +data[k] * other.data[l + j * other.sizex]; // WORKS FOR REAL ONLY!
    //    for(l=0;l<other.sizex;l++) s2[l] = ExOp::mktrju(s2[l]) + s2[l] + data[k] * other.data[l + j * other.sizex]; 
        k++;
            fout.data[0] += s2[0] * other.data[0 + j * other.sizex];
            fout.data[1] += s2[1] * other.data[0 + j * other.sizex] + s2[0] *  other.data[1 + j * other.sizex]; // a*c - b*d
            fout.data[2] += s2[1] * other.data[1 + j * other.sizex];
            fout.data[5] += s2[2] * other.data[2 + j * other.sizex];
            fout.data[6] += s2[3] * other.data[0 + j * other.sizex] + s2[0] *  other.data[3 + j * other.sizex]; // a*c - b*d
            fout.data[7] += s2[3] * other.data[1 + j * other.sizex] + s2[1] *  other.data[3 + j * other.sizex]; // a*c - b*d
        }
    delete[](s2);
    return fout;
}
LFHTEMP const Trianglix<C, 0u>& Trianglix<C, 0u>::operator*=(const Matrix<C>& other){
}

LFHTEMP Trianglix<C, 0u> Trianglix<C, 0u>::inverse_OLD() const{
     
    Trianglix<C, 0u> fout(*this);
   
    if (t_size >2){
    C* buf = new C[t_size*2];
    C* hv = new C[t_size *(t_size+1)/2];
    double* normbuf = new double[t_size-2];
    unsigned int i,j,k;
      for(j=t_size-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mktrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mktrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mktrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }

        fout.data[0] = ExOp::mkinverse(fout.data[0]);
        buf[t_size] = fout.data[0]; 
        for(j=1;j<t_size;j++){
            k = (j*(j+3))/2;
            buf[t_size+j-1] *= fout.data[k-1];
            fout.data[k] = ExOp::mkinverse(fout.data[k] - ExOp::mktrju(fout.data[k-1]) * buf[t_size+j-1]);
            ExOp::toZero(fout.data[k-1]);
            buf[t_size+j] = fout.data[k];
            }


    for(j=t_size-2;j<t_size;j--) fout.offdiagelimination_up( -buf[t_size+j],j,buf);

    for(j=2;j<t_size;j++) fout.HouseHolderMultiply(hv+((j * (j+1))/2), normbuf[j-2], j+1, buf,false);
        
    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (t_size == 2){
        C denum = ExOp::mkinverse(data[0] * data[2] - data[1] * ExOp::mktrju(data[1]));
        fout.data[0] = data[2] * denum;
        fout.data[1] = -data[1] * denum;
        fout.data[2] = data[0] * denum;
    }  else if (t_size == 1) fout.data[0] = ExOp::mkinverse(data[0]);
    return(fout);
	

	
}

LFHTEMP Trianglix<C, 0u> Trianglix<C, 0u>::inverse_MK2() const{
    Trianglix<C, 0u> fout(*this);
    
     Trianglix<C, 0u> backup(*this);
    if (t_size >2){
    C* buf = new C[t_size*2];
    C* hv = new C[t_size *(t_size+1)/2];
    C* bufRQ = new C[t_size*2];
    double normnorm;
    double* normbuf = new double[t_size-2];
    unsigned int i,j,k;
      for(j=t_size-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mktrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mktrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mktrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }
            

        bufRQ[0] = fout.data[0];
        bufRQ[1] = fout.data[1];


        for(j=1;j<t_size;j++){
            k = (j*(j+3))/2;
            
            normnorm = pow(ExOp::pnorm(bufRQ[((j-1)<<1)]) + ExOp::pnorm(fout.data[k-1] ), -0.5f);
            
            if (j > 1) fout.data[k-j-2] = bufRQ[((j-2)<<1) | 1 ] * (1.0f / normnorm);
        //    printf("%e %e %e\n", bufRQ[((j-1)<<1)],fout.data[k-1] , normnorm );
            bufRQ[((j-1)<<1)] *= normnorm;

         //   bufRQ[(j<<1)] = fout.data[k+j+1] * bufRQ[((j-1)<<1)];
          //  bufRQ[(j<<1)] = fout.data[k]*bufRQ[((j-1)<<1)]  ;
            if (j+1 < t_size) bufRQ[(j<<1) | 1] =  bufRQ[((j-1)<<1)] *fout.data[k+j+1];
            bufRQ[(j<<1)] = bufRQ[((j-1)<<1)] * fout.data[k];
            bufRQ[(j<<1)] -= fout.data[k-1] *bufRQ[((j-1)<<1) | 1 ] * normnorm;
            
         //   bufRQ[(j<<1)] -= bufRQ[((j-1)<<1) | 1]  *bufRQ[((j-1)<<1) | 1];

            //printf("partial x=%e y=%e  %e\n", bufRQ[(j<<1)], bufRQ[(j<<1) | 1] , fout.data[k-1] *bufRQ[((j-1)<<1) | 1 ] * normnorm );
            bufRQ[((j-1)<<1)| 1] = fout.data[k-1] * normnorm;

            }
        fout.data[k-1] = bufRQ[((j-2)<<1) | 1 ] * bufRQ[((j-1)<<1)]; 
     //   fout.data[k] =bufRQ[((j-1)<<1)]; 
        backup.show();
        printf("first rotation! %e %e\n", bufRQ[0] ,bufRQ[1] );
        backup.QR_back(bufRQ[0], bufRQ[1], 0);backup.show();
        printf("second rotation! %e %e\n", bufRQ[2] ,bufRQ[3] );
        backup.QR_back(bufRQ[2], bufRQ[3], 1);backup.show();
        printf("third rotation! %e %e\n", bufRQ[4] ,bufRQ[5] );
        backup.QR_back(bufRQ[4], bufRQ[5], 2);backup.show();
        fout.show();
        
        fout.data[0] = ExOp::mkinverse(fout.data[0]);
        buf[t_size] = fout.data[0]; 
        for(j=1;j<t_size;j++){
            k = (j*(j+3))/2;
            buf[t_size+j-1] *= fout.data[k-1];
            fout.data[k] = ExOp::mkinverse(fout.data[k] - ExOp::mktrju(fout.data[k-1]) * buf[t_size+j-1]);
            ExOp::toZero(fout.data[k-1]);
            buf[t_size+j] = fout.data[k];
            }

    for(j=t_size-2;j<t_size;j--) {fout.offdiagelimination_up_backroutine( -buf[t_size+j],j,buf); if (j == t_size-2) fout.show();}
   


    for(j=2;j<t_size;j++) {fout.HouseHolderMultiply(hv+((j * (j+1))/2), normbuf[j-2], j+1, buf,false);} // pts55
        
    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (t_size == 2){
        C denum = ExOp::mkinverse(data[0] * data[2] - data[1] * ExOp::mktrju(data[1]));
        fout.data[0] = data[2] * denum;
        fout.data[1] = -data[1] * denum;
        fout.data[2] = data[0] * denum;
    }  else if (t_size == 1) fout.data[0] = ExOp::mkinverse(data[0]);
    return(fout);
}

LFHTEMP Trianglix<C, 0u> Trianglix<C, 0u>::mkinverse() const{
    Trianglix<C, 0u> fout(*this);
   
    if (t_size >2){
    C* buf = new C[t_size*2];
    C* hv = new C[t_size *(t_size+1)/2];
    double* normbuf = new double[t_size-2];
    C tmp;

    unsigned int i,j,k,l;
      for(j=t_size-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mktrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mktrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mktrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }

        // interative method for the 
        Matrix<double> tmptmp_triangle = Matrix<double>(fout); 

        Tuple<C,0u> tridiago; tridiago.setSize(t_size*2-1);
        Tuple<C,0u> x_tridiago; x_tridiago.setSize(t_size-1);
        Tuple<C,0u> refine_buf; refine_buf.setSize(t_size);
        tridiago[0] = fout.data[0];
        x_tridiago[0] = fout.data[1] / fout.data[0]; 
        j=0;
        while(true){
            k = ((j+1)*(j+4))/2;
            tridiago[(j<<1)| 1] = fout.data[k-1];
            j++;
            tridiago[(j<<1)] = fout.data[k];
            if (j == t_size-1) break;
            x_tridiago[j] = fout.data[k-1] / (fout.data[k] - fout.data[k-1] * x_tridiago[j-1]);
        }


        fout.data[0] = ExOp::mkinverse(fout.data[0]);
        buf[t_size] = fout.data[0]; 
        for(j=1;j<t_size;j++){
            k = (j*(j+3))/2;
            buf[t_size+j-1] *= fout.data[k-1];
            fout.data[k] = ExOp::mkinverse(fout.data[k] - ExOp::mktrju(fout.data[k-1]) * buf[t_size+j-1]);
            ExOp::toZero(fout.data[k-1]);
            buf[t_size+j] = fout.data[k];
            }


    for(j=t_size-2;j<t_size;j--) fout.offdiagelimination_up( -buf[t_size+j],j,buf);
   
    // iterative refinement for tri-diagonal matrices!

 //   k = ((t_size-1) *t_size)/2;
 //   ExOp::show(fout.data[k] * tridiago[0] + fout.data[k+1] * tridiago[1]);
 //   for(i=2,k++;(i>>1) < t_size-1;i+=2, k++){
 //       ExOp::show(fout.data[k-1] * tridiago[i-1] + fout.data[k] * tridiago[i] +fout.data[k+1] * tridiago[i+1]);
 //       }
 //   ExOp::show(fout.data[k-1] * tridiago[i-1] + fout.data[k] * tridiago[i]);


//   printf("before\n");   (Matrix<double>(fout) * tmptmp_triangle).show();
		
    for(l=0;l<t_size;l++){ // makes this n^3 like the hols
    for(j = t_size-1; j > 0; j--){
    k = (j *(j+1))/2;

 //  ExOp::show(fout.data[k] * tridiago[0] + fout.data[k+1] * tridiago[1]);
    refine_buf[0] = -fout.data[k] - ((fout.data[k+1] * tridiago[1]) / tridiago[0]);
    for(i=2,k++;(i>>1) < j ;i+=2, k++){
    //  ExOp::show(fout.data[k-1] * tridiago[i-1] + fout.data[k] * tridiago[i] +fout.data[k+1] * tridiago[i+1]);
        refine_buf[(i>>1)] = ( (fout.data[k-1] + refine_buf[(i>>1)-1]) * tridiago[i-1] + fout.data[k] * tridiago[i] +fout.data[k+1] * tridiago[i+1]  ) / (tridiago[i-1] * x_tridiago[(i>>1)] - tridiago[i]);
	}
    // objective is 1
    tmp = -1.0f;
    if (j < t_size-1)  tmp += tridiago[i+1] * fout.data[k+j+1];
    
  //  ExOp::show( tmp  + fout.data[k-1] * tridiago[i-1] + fout.data[k] * tridiago[i]);
    tmp = ( tmp  + (fout.data[k-1] + refine_buf[(i>>1)-1]) * tridiago[i-1] + fout.data[k] * tridiago[i]) / (tridiago[i-1] * x_tridiago[(i>>1)] - tridiago[i]);
    fout.data[k] += tmp;// printf("%e refine\n",tmp);
    for(i=(j << 1)-2,k--; i != 0xFFFFFFFE;i-=2, k--){
    tmp = refine_buf[(i>>1)] - tmp * x_tridiago[(i>>1)];
    fout.data[k] += tmp;// printf("%e refine\n",tmp);
    }
    }

    fout.data[0] = (1.0f -  fout.data[1] * tridiago[1]) / tridiago[0];
    }

  
/*
    k = ((t_size-1) *t_size)/2;
    ExOp::show(fout.data[0] * tridiago[0] + fout.data[1] * tridiago[1]);
   refine_buf[0] = -fout.data[k] - ((-1.0f + fout.data[k+1] * tridiago[1]) / tridiago[0]);
    j=1;
    for(i=2,k+=j;(i>>1) < t_size-1;i+=2, k+=j){
        ExOp::show(fout.data[k-j] * tridiago[i-1] + fout.data[k] * tridiago[i] +fout.data[k+j+1] * tridiago[i+1]);
        refine_buf[(i>>1)] = ( (fout.data[k-j] + refine_buf[(i>>1)-1]) * tridiago[i-1] + fout.data[k] * tridiago[i] +fout.data[k+j+1] * tridiago[i+1]  ) / (tridiago[i-1] * x_tridiago[(i>>1)] - tridiago[i]);
        
        j++;
        }
    // objective is 1
    ExOp::show(fout.data[k-j] * tridiago[i-1] + fout.data[k] * tridiago[i]);
    tmp = ((fout.data[k-j] + refine_buf[(i>>1)-1]) * tridiago[i-1] + fout.data[k] * tridiago[i]) / (tridiago[i-1] * x_tridiago[(i>>1)] - tridiago[i]);
    fout.data[k] += tmp; printf("%e refine\n",tmp);
    for(i-=2,k-=j; i != 0xFFFFFFFE;i-=2, k-=j){
    tmp = refine_buf[(i>>1)] - tmp * x_tridiago[(i>>1)];
    fout.data[k] += tmp; printf("%e refine\n",tmp);
    }*/
    
 //   printf("after\n");   (Matrix<double>(fout) * tmptmp_triangle).show();


    for(j=2;j<t_size;j++) fout.HouseHolderMultiply(hv+((j * (j+1))/2), normbuf[j-2], j+1, buf,false);
        
    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (t_size == 2){
        C denum = ExOp::mkinverse(data[0] * data[2] - data[1] * ExOp::mktrju(data[1]));
        fout.data[0] = data[2] * denum;
        fout.data[1] = -data[1] * denum;
        fout.data[2] = data[0] * denum;
    }  else if (t_size == 1) fout.data[0] = ExOp::mkinverse(data[0]);
    return(fout);
}
	
LFHTEMP C Trianglix<C, 0u>::maxEigenValue()const{ C daans; return daans;} // TODO!
	
LFHTEMP Tuple<C, 0u> Trianglix<C, 0u>::getEigenValues() const{
Trianglix<C, 0u> fout(*this);
    Tuple<C, 0u> real_fout; real_fout.setSize(t_size);
    C prod;
    if (t_size >2){
    C* buf = new C[t_size*2];
    C* hv = new C[t_size *(t_size+1)/2];
    double* normbuf = new double[t_size-2];
    unsigned int i,j,k;
      for(j=t_size-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mktrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mktrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mktrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }

        real_fout[0] = fout.data[0];
        buf[t_size] = ExOp::mkinverse(fout.data[0]); 
        for(j=1;j<t_size;j++){
            k = (j*(j+3))/2;
            buf[t_size+j-1] *= fout.data[k-1];
            fout.data[k] -= ExOp::mktrju(fout.data[k-1]) * buf[t_size+j-1];
            real_fout[j] = fout.data[k];
            buf[t_size+j] = ExOp::mkinverse(fout.data[k]);
            }
        
    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (t_size == 2) {
        real_fout[0] = fout.data[0] - fout.data[2];
        real_fout[1] = sqrt(real_fout[0] * real_fout[0] + (ExOp::mktrju(fout.data[1]) * fout.data[1] * 4.0f) ) ;  
        real_fout[0] = (fout.data[0] + fout.data[2] - real_fout[1]) *0.5f;
        real_fout[1] += real_fout[0];
    }else if (t_size == 1) real_fout[0] = data[0];
    return(real_fout);
}/*
LFHTEMP C Trianglix<C, 0u>::log_determinant() const{
Trianglix<C, 0u> fout(*this);
   C log_prod;
    if (t_size >2){
    C* buf = new C[t_size*2];
    C* hv = new C[t_size *(t_size+1)/2];
    double* normbuf = new double[t_size-2];
    unsigned int i,j,k;
      for(j=t_size-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mktrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mktrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mktrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }

        log_prod = log(fout.data[0]);
        ExOp::toZero(hv[0]);
        if (fout.data[0] > hv[0]) hv[0] = fout.data[0];
        buf[t_size] = ExOp::mkinverse(fout.data[0]); 
        for(j=1;j<t_size;j++){
            k = (j*(j+3))/2;
            buf[t_size+j-1] *= fout.data[k-1];
            fout.data[k] -= ExOp::mktrju(fout.data[k-1]) * buf[t_size+j-1];
            log_prod  += log(fout.data[k]);
            if (fout.data[k] > hv[0]) hv[0] = fout.data[k];
            buf[t_size+j] = ExOp::mkinverse(fout.data[k]);
            }
        k =0;
        log_prod = (fout.data[k] > 0.0f) ? log(fout.data[k]) : log(hv[0]) + log(ExCo<double>::epsilon());
        for(j=1;j<t_size;j++){k = (j*(j+3))/2;
        log_prod += (fout.data[k] > 0.0f) ? log(fout.data[k]) : log(hv[0]) + log(ExCo<double>::epsilon());
        }
        
    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (t_size == 2) return(log(data[0] *data[2] - data[1] * ExOp::mktrju(data[1])));
    else if (t_size == 1) return(log(data[0]));
    else ExOp::toZero(log_prod);
    return(log_prod);
}*/


LFHTEMP double Trianglix<C, 0u>::bhattacharryya_partial(const Tuple<C, 0u> &dev ,const Trianglix<C, 0u>& other)const{
    double fout;
    // dev (P+P2/2)-1 dev 
    C sum;
    C prod;
    Trianglix<C, 0u> invsc = (*this) + other;
    Tuple<C, 0u> xdev = dev;
    if (xdev.tup_size >2){
    C* buf = new C[xdev.tup_size*2];
    C* hv = new C[xdev.tup_size *(xdev.tup_size+1)/2];
    double* normbuf = new double[xdev.tup_size-2];
    unsigned int i,j,k;
      for(j=xdev.tup_size-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(invsc.data[k]);hv[k] = ExOp::mktrju(invsc.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(invsc.data[k+i]); hv[k+i] = ExOp::mktrju(invsc.data[k+i]);}
            hv[k+i] = ExOp::mktrju(invsc.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            invsc.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
            xdev.HouseHolderMultiply(hv+k, normbuf[j-2], j);
        }

        buf[xdev.tup_size] = ExOp::mkinverse(invsc.data[0]);
        prod =  invsc.data[0];
        

        ExOp::toZero(hv[0]);
        if (invsc.data[0] > hv[0]) hv[0] = invsc.data[0];
        for(j=1;j<xdev.tup_size;j++){
            k = (j*(j+3))/2;
            xdev[j] -= ExOp::mktrju(xdev[j-1]) * buf[xdev.tup_size+j-1]*invsc.data[k-1];
            if (j==1) sum = ExOp::mktrju(xdev[0]) * xdev[0] * buf[xdev.tup_size]; 
            else sum += ExOp::mktrju(xdev[j-1]) * xdev[j-1] * buf[xdev.tup_size+j-1]; 
            invsc.data[k] -= ExOp::mktrju(invsc.data[k-1]) * buf[xdev.tup_size+j-1]*invsc.data[k-1];
            if (invsc.data[k] > hv[0]) hv[0] = invsc.data[k];
            buf[xdev.tup_size+j] = ExOp::mkinverse(invsc.data[k]);
       //     sum += ExOp::mktrju(xdev[j]) * xdev[j] * buf[SIZE+j];
            }
        sum += ExOp::mktrju(xdev[j-1]) * xdev[j-1] * buf[xdev.tup_size+j-1];
        prod = xdev.tup_size * log(0.5f);
        for(j=0;j<xdev.tup_size;j++) {k = (j*(j+3))/2;prod += (invsc.data[k] > 0.0f) ? log(invsc.data[k]) : log(hv[0]) + log(ExCo<double>::epsilon());}
        fout = 0.5f *( ExOp::norm(sum) + prod);

        
    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (xdev.tup_size == 2){
       // TODO

    }  else fout = 0.5f *(ExOp::mktrju(dev[0]) * dev[0] * ExOp::mkinverse(data[0] + other.data[0]) + log(0.5f) + log(data[0] + other.data[0]));


    return(fout);
    }

    
LFHTEMP C Trianglix<C, 0u>::Xformed_inner_product_of_inverse(const Tuple<C, 0u> &dev) const{
    double fout;
	
	
    // dev (P+P2/2)-1 dev 
    C sum;
    C prod;
    Trianglix<C, 0u> invsc = (*this);
    Tuple<C, 0u> xdev = dev;
	if (getSize() < xdev.getSize()) exit(1);
	unsigned int siz = xdev.getSize();
	
    if (siz >2){
		//printf("alloc to die %i\n",siz); fflush(stdout);

    C* buf = new C[siz*2];
    C* hv = new C[siz *(siz+1)/2];
    double* normbuf = new double[siz-2];
		
    unsigned int i,j,k;
	//	printf("%i\t%t\n", invsc.getSize(), xdev.getSize()); fflush(stdout);
		for(j=siz-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(invsc.data[k]);hv[k] = ExOp::mktrju(invsc.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(invsc.data[k+i]); hv[k+i] = ExOp::mktrju(invsc.data[k+i]);}
            hv[k+i] = ExOp::mktrju(invsc.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            invsc.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
            xdev.HouseHolderMultiply(hv+k, normbuf[j-2], j);
        }
		

        buf[xdev.tup_size] = ExOp::mkinverse(invsc.data[0]);
        prod =  invsc.data[0];
        

        ExOp::toZero(hv[0]);
        if (invsc.data[0] > hv[0]) hv[0] = invsc.data[0];
        for(j=1;j<siz;j++){
            k = (j*(j+3))/2;
            xdev[j] -= ExOp::mktrju(xdev[j-1]) * buf[siz+j-1]*invsc.data[k-1];
            if (j==1) sum = ExOp::mktrju(xdev[0]) * xdev[0] * buf[siz]; 
            else sum += ExOp::mktrju(xdev[j-1]) * xdev[j-1] * buf[siz+j-1]; 
            invsc.data[k] -= ExOp::mktrju(invsc.data[k-1]) * buf[siz+j-1]*invsc.data[k-1];
            if (invsc.data[k] > hv[0]) hv[0] = invsc.data[k];
            buf[siz+j] = ExOp::mkinverse(invsc.data[k]);
       //     sum += ExOp::mktrju(xdev[j]) * xdev[j] * buf[SIZE+j];
            }
        sum += ExOp::mktrju(xdev[j-1]) * xdev[j-1] * buf[xdev.tup_size+j-1];
        fout = sum;
		
        
    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (siz == 2){
       fout = ExOp::mktrju(dev[1]) * dev[0] * data[1];
       fout = (ExOp::mktrju(dev[0]) * dev[0]* data[2] + ExOp::mktrju(dev[1]) * dev[1]* data[0] - fout - ExOp::mktrju(fout) ) * ExOp::mkinverse(data[0] * data[2] - data[1] * ExOp::mktrju(data[1]));
    }  else fout = ExOp::mktrju(dev[0]) * dev[0] * ExOp::mkinverse(data[0]);


    return(fout);
    }

    
    
LFHTEMP void Trianglix<C, 0u>::eigen() const{
    Trianglix<C, 0u> fout(*this);
    if (getSize() >2){
    C* buf = new C[getSize()];
    C* hv = new C[getSize() *(getSize()+1)/2];
    double* normbuf = new double[getSize()-2];
    unsigned int i,j,k;
      for(j=getSize()-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mktrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mktrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mktrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf);
        }
        // tri-diago form!

        buf[0] = fout.data[0]; printf("eigen vals! %f", ExOp::mkrealproj(buf[0]) );
        for(j=1;j<getSize();j++){
            k = (j*(j+3))/2; 
            buf[j] = data[k] - data[k-1] * ExOp::mktrju(data[k-1]) / buf[j-1];printf("\t%f", ExOp::mkrealproj(buf[j]));
            }
        printf("\n");
        
        fout.offdiagelimination_down(-fout.data[1] / fout.data[0], 0, buf);


//     for(j=2;j<getSize();j++) fout.HouseHolderMultiply(hv+((j * (j+1))/2), normbuf[j-2], j+1, buf);
        
    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (getSize() == 2){
        C denum = ExOp::mkinverse(data[0] * data[2] - data[1] * ExOp::mktrju(data[1]));
        fout.data[0] = data[2] * denum;
        fout.data[1] = -data[1] * denum;
        fout.data[2] = data[0] * denum;
    }  else if (getSize() == 1) fout.data[0] = ExOp::mkinverse(data[0]);
    return(fout);
}

LFHTEMP void Trianglix<C, 0u>::show(FILE* f, int level)const{
    unsigned int i,ts,k,j; 
    for(i=0,ts=totsize(),j=0,k=2;i<ts;i++) {ExOp::show(data[i],f,level+1); if (i==j) {fprintf(f,"\n"); j+= k++;} else fprintf(f,"\t");}
}
LFHTEMP string  Trianglix<C, 0u>::type_tostring() const{return string("Trianglix<") + ExOp::type_tostring(data[0]) +string(",0>");}


LFHTEMP C Trianglix<C, 0u>::log_determinant() const{

	Trianglix<C, 0u> fout(*this);
    C deter[2];
    C tmp;
    int safe_mag = 0;
    if (getSize() >2){
    C* buf = new C[getSize()*2];
    C* hv = new C[getSize() *(getSize()+1)/2];
    double* normbuf = new double[getSize()-2];
    unsigned int i,j,k;
      for(j=getSize()-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mktrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mktrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mktrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }
    delete[](hv);
    delete[](buf);
    delete[](normbuf);

        deter[0] = fout.data[0];
        ExOp::toOne(deter[1]);

        buf[getSize()] = ExOp::mkinverse(fout.data[0]); 
        for(j=1;j<getSize();j++){
            k = (j*(j+3))/2;
            tmp = deter[((j & 1)^1)] * fout.data[k];
            tmp -= deter[(j & 1)] * (ExOp::mktrju(fout.data[k-1]) * fout.data[k-1]);  
            if (!(ExOp::isValid(tmp))){ // overflow!
               deter[0] *= pow(0.5f, 300.0f); 
               deter[1] *= pow(0.5f, 300.0f); 
                tmp = deter[((j & 1)^1)] * fout.data[k];
                tmp -= deter[(j & 1)] * (ExOp::mktrju(fout.data[k-1]) * fout.data[k-1]);  
                safe_mag++;
            }else if (fabs(tmp) <= pow(0.5f, 300.0f)) { // underflow!
               deter[0] *= pow(2.0f, 300.0f); 
               deter[1] *= pow(2.0f, 300.0f); 
               tmp = deter[((j & 1)^1)] * fout.data[k];
                tmp -= deter[(j & 1)] * (ExOp::mktrju(fout.data[k-1]) * fout.data[k-1]);  
               safe_mag--;
            }
            deter[(j & 1)]  = tmp;

            }
        deter[0] = log(ExOp::norm(deter[(getSize() & 1)^1])); 


       return deter[0] + log(pow(2.0,300.0f)) * safe_mag; 
    }else if (getSize() == 2) return(log(data[0] *data[2] - data[1] * ExOp::mktrju(data[1])));
    else if (getSize() == 1) return(log(data[0]));
	else return(0.0f);
}



#undef LFHTEMP
#define LFHTEMP template <class C, unsigned int SIZE>

LFHTEMP void Trianglix<C, SIZE>::HouseHolderMultiply(const C * const vec, double denum2, unsigned int length, C* buf, bool hint ){
    // house holder multiplication, done twice!
    // = T - vec * vec' T - T' vec * vec' + vec * vec' * T *vec * vec'
	if ((denum2 == 0.0f)||(!ExCo<double>::isValid(denum2))) return;
    unsigned int i,j,k;
    C s;
        
        buf[0] = data[0] * vec[0];
        s = data[0] * ExOp::mktrju(vec[0]) * vec[0];
        for(j=1,k=1;j<length;j++) {
        buf[j] = data[k] * vec[0]; buf[0] += ExOp::mktrju(data[k++]) * vec[j];
        for(i=1;i<j;i++) {buf[j] += data[k] * vec[i]; buf[i] += ExOp::mktrju(data[k++]) * vec[j];}
        s += data[k] * ExOp::mktrju(vec[j]) * vec[j] + (buf[j] *  ExOp::mktrju(vec[j]) + ExOp::mktrju(buf[j]) * vec[j]);
        buf[j] += data[k] * vec[j];
        k++;
        }
        if (!hint){
        for(;j<SIZE;j++) {
        buf[j] = data[k++] * vec[0];
        for(i=1;i<length;i++) buf[j] += data[k++] * vec[i];
            k+= j - length+1;
        }
        for(j=0;j<SIZE;j++) buf[j] /= denum2;
        }else for(j=0;j<length;j++) buf[j] /= denum2;
        s /= denum2 * denum2; 

        for(j=0,k=0;j<length;j++) {
        for(i=0;i<j;i++) data[k++] += s * vec[j] *ExOp::mktrju(vec[i]) - buf[j]*ExOp::mktrju(vec[i]) - vec[j] * ExOp::mktrju(buf[i]) ;
        data[k++] += s * vec[j] *ExOp::mktrju(vec[j]) - buf[j]*ExOp::mktrju(vec[j]) - vec[j] * ExOp::mktrju(buf[j]);
        }
        if (hint) return;
        for(;j<SIZE;j++) {
            for(i=0;i<length;i++) data[k++] -= buf[j]*ExOp::mktrju(vec[i]);
            k+= j - length+1;
        }
}

	LFHTEMP void Trianglix<C, SIZE>::CholeskyStep_up(const C * const vec, unsigned int length, C* buf){
		// = T - vec * vec' T - T' vec * vec' + vec * vec' * T *vec * vec'
		unsigned int i,j,k,l;
		
		
		// computes the ( partial vec * T)
		
		buf[0] = data[0] * vec[0];
		for(i=1,k=1;i<SIZE;i++){
			buf[i] = data[k] * vec[0];
			if (i < length-1)	buf[0] += data[k] *  vec[i];
			else if (i == length -1) buf[0] += data[k] * (vec[length -1] - 1.0f);
			for(j=1,k++;j<i;j++,k++){
				if (j < length-1)	buf[i] += data[k] * vec[j];
				else if (j == length-1) buf[i] += data[k] * (vec[length -1] -1.0f);
				if (i < length-1)	buf[j] += data[k] * vec[i];
				else if (i == length-1) buf[j] += data[k] * (vec[length -1] -1.0f);
			}
			if (i < length-1)	buf[i] += data[k] * vec[i];
			else if (i == length-1) buf[i] += data[k] * (vec[length -1] -1.0f);
			k++;
		}
		
		//for(i=0;i<SIZE;i++) printf("%f\t", buf[i]);
		//printf("\n", buf[i]);
		

		
		k = (length *(length-1))/2;
		for(i=0;i<length-1;i++) data[k++] += 2.0f * buf[i];
		for(l=0;l< length-1;l++) data[k] += buf[l] * vec[l];
		data[k] += 2.0f * buf[i++] + buf[l] * (vec[l] -1.0f);
		for(k+=i+1;i<SIZE;k+=i+1) { data[k] += 2.0f * buf[i++];} 
	}
	
	
	
LFHTEMP void Trianglix<C, SIZE>::offdiagelimination_down(const C &fact,unsigned int col, C* buf){ // multiply row/col "col" by  fact, and adds it in "(col+1)" 
    // = T - vec * vec' T - T' vec * vec' + vec * vec' * T *vec * vec'
    unsigned int i,j,k;

    k = (col *(col+1))/2;
    for(i=0;i<=col;i++) {buf[i] = fact * data[k+i];}
    for(;i<SIZE;i++) {k += i; buf[i] = fact * data[k+col];}
    k = ((col+2) *(col+1))/2;
    for(i=0;i<=col;i++) data[k+i] += buf[i];
    k++;
    data[k+col] += buf[col] * ExOp::mktrju(fact) + ExOp::mktrju(buf[i]);
    for(i++;i<SIZE;i++) {k += i; data[k+col] += buf[i];}
}

LFHTEMP void Trianglix<C, SIZE>::offdiagelimination_up(const C &fact,unsigned int col, C* buf){ // multiply row/col "col" by  fact, and adds it in "(col+1)" 
    // = T - vec * vec' T - T' vec * vec' + vec * vec' * T *vec * vec'
    unsigned int i,j,k;

    k = ((col+2) *(col+1))/2;
    for(i=0;i<=col;i++) buf[i] = fact * data[k+i];
    buf[i] = fact * data[k+i]; 
    k++;
    for(i++;i<SIZE;i++) {k += i; buf[i] = fact * data[k+col];}

    k = ((col) *(col+1))/2;
    for(i=0;i<col;i++) data[k+i] += buf[i];
    data[k+col] += ExOp::mktrju(buf[col+1]) * fact + ExOp::mktrju(buf[i]); 
    for(i++;i<SIZE;i++) {k += i; data[k+col] += buf[i];}
	
}


LFHTEMP Trianglix<C, SIZE>::Trianglix(const Tuple<C, SIZE>& other){
    unsigned int i,j,k;
    for(j=0,k=0;j<SIZE;j++) for(i=0;i<=j;i++,k++) data[k] = ExOp::mktrju(other[i]) * other[j];    
}

LFHTEMP Trianglix<C, SIZE>& Trianglix<C, SIZE>::operator=(const Tuple<C, SIZE>& other){
    unsigned int i,j,k;
    for(j=0,k=0;j<SIZE;j++) for(i=0;i<=j;i++,k++) data[k] = ExOp::mktrju(other[i]) * other[j];
    return(*this);
}

LFHTEMP Tuple<C,SIZE> Trianglix<C, SIZE>::getEigenValues() const{
Trianglix<C, SIZE> fout(*this);
   Tuple<C,SIZE> real_fout;
    if (SIZE >2){
    C* buf = new C[SIZE*2];
    C* hv = new C[SIZE *(SIZE+1)/2];
    double* normbuf = new double[SIZE-2];
    unsigned int i,j,k;
      for(j=SIZE-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mktrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mktrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mktrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }

        real_fout[0] = fout.data[0];
        buf[SIZE] = ExOp::mkinverse(fout.data[0]); 
        for(j=1;j<SIZE;j++){
            k = (j*(j+3))/2;
            buf[SIZE+j-1] *= fout.data[k-1];
            fout.data[k] -= ExOp::mktrju(fout.data[k-1]) * buf[SIZE+j-1];
            real_fout[j] = fout.data[k];
            buf[SIZE+j] = ExOp::mkinverse(fout.data[k]);
            }
        
    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (SIZE == 2) {
    // P(x) = x2 - data[0] - data[2]x data[0] * data[2] - mktr[1] * data[1];
    // discr = (data[0] - data[1])^2 + mktr[1] * data[1]; hence, POSITIVE and real!
    real_fout[0] = fout.data[0] - fout.data[2];
    real_fout[1] = sqrt(real_fout[0] * real_fout[0] + (ExOp::mktrju(fout.data[1]) * fout.data[1] * 4.0f) ) ;  
    real_fout[0] = (fout.data[0] + fout.data[2] - real_fout[1]) *0.5f;
    real_fout[1] += real_fout[0];
    }else if (SIZE == 1) real_fout[0] = data[0];
    return(real_fout);

}

LFHTEMP C Trianglix<C, SIZE>::determinant() const{
Trianglix<C, SIZE> fout(*this);
    C deter[2];
    if (SIZE >2){
    C* buf = new C[SIZE*2];
    C* hv = new C[SIZE *(SIZE+1)/2];
    double* normbuf = new double[SIZE-2];
    unsigned int i,j,k;
      for(j=SIZE-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mktrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mktrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mktrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }

        deter[0] = fout.data[0];
        ExOp::toOne(deter[1]);

        buf[SIZE] = ExOp::mkinverse(fout.data[0]); 
        for(j=1;j<SIZE;j++){
            k = (j*(j+3))/2;
            deter[(j & 1)] *= -(ExOp::mktrju(fout.data[k-1]) * fout.data[k-1]);
            deter[(j & 1)] += deter[((j & 1)^1)] * fout.data[k];
            }
        
    delete[](hv);
    delete[](buf);
    delete[](normbuf);
       return (SIZE & 1) ? -deter[0] : deter[1]; 
    }else if (SIZE == 2) return(log(data[0] *data[2] - data[1] * ExOp::mktrju(data[1])));
    else if (SIZE == 1) return(log(data[0]));

}

LFHTEMP C Trianglix<C, SIZE>::log_determinant() const{
Trianglix<C, SIZE> fout(*this);
    C deter[2];
    C tmp;
    int safe_mag = 0;
    if (SIZE >2){
    C* buf = new C[SIZE*2];
    C* hv = new C[SIZE *(SIZE+1)/2];
    double* normbuf = new double[SIZE-2];
    unsigned int i,j,k;
      for(j=SIZE-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mktrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mktrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mktrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }
    delete[](hv);
    delete[](buf);
    delete[](normbuf);

        deter[0] = fout.data[0];
        ExOp::toOne(deter[1]);

        buf[SIZE] = ExOp::mkinverse(fout.data[0]); 
        for(j=1;j<SIZE;j++){
            k = (j*(j+3))/2;
            tmp = deter[((j & 1)^1)] * fout.data[k];
            tmp -= deter[(j & 1)] * (ExOp::mktrju(fout.data[k-1]) * fout.data[k-1]);  
            if (!(ExOp::isValid(tmp))){ // overflow!
               deter[0] *= pow(0.5f, 300.0f); 
               deter[1] *= pow(0.5f, 300.0f); 
                tmp = deter[((j & 1)^1)] * fout.data[k];
                tmp -= deter[(j & 1)] * (ExOp::mktrju(fout.data[k-1]) * fout.data[k-1]);  
                safe_mag++;
            }else if (fabs(tmp) <= pow(0.5f, 300.0f)) { // underflow!
               deter[0] *= pow(2.0f, 300.0f); 
               deter[1] *= pow(2.0f, 300.0f); 
               tmp = deter[((j & 1)^1)] * fout.data[k];
                tmp -= deter[(j & 1)] * (ExOp::mktrju(fout.data[k-1]) * fout.data[k-1]);  
               safe_mag--;
            }
            deter[(j & 1)]  = tmp;

            }
        deter[0] = log(ExOp::norm(deter[(SIZE & 1)^1])); 


       return deter[0] + log(pow(2.0,300.0f)) * safe_mag; 
    }else if (SIZE == 2) return(log(data[0] *data[2] - data[1] * ExOp::mktrju(data[1])));
    else if (SIZE == 1) return(log(data[0]));
}

	
	LFHTEMP C Trianglix<C, SIZE>::maxEigenValue()const{
		Trianglix<C, SIZE> fout(*this);
		Tuple<C, SIZE> da_vec;
		C da_tmp[2];
		C da_norm;
		C da_dev;
		C da_old_norm;
		double tmptmp;
		if (SIZE >2){
			C* buf = new C[SIZE*2];
			C* hv = new C[SIZE *(SIZE+1)/2];
			double* normbuf = new double[SIZE-2];
			unsigned int i,j,k;
//			fout.show();
			for(j=SIZE-1;j>1;j--){
				k = (j * (j+1))/2;
				normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mktrju(fout.data[k]);
				for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mktrju(fout.data[k+i]);}
//				printf("da pre scale : %e\n", normbuf[j-2]);
				
				hv[k+i] = ExOp::mktrju(fout.data[k+i]);
				tmptmp = 1.0f -sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
				if (ExOp::isValid(tmptmp)) {
					hv[k+i] *=  tmptmp;
					normbuf[j-2] += ExOp::pnorm(hv[k+i]);
					normbuf[j-2] *=0.5f;
				}else{ // hv[k+i] is too small, assume its zero!
//					printf("hello! has zero!\n");
					hv[k+i] = -sqrt(normbuf[j-2]);
				}
				
				ExOp::toZero(hv[k+j]); 
//				printf("da scale : %e\n", normbuf[j-2]);
				fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);//fout.show();
			}
//			fout.show();
			delete[](normbuf);

			for(j=2;j<SIZE;j++){
				k = (j*(j+3))/2;
				fout.data[(j<<1)-1] = fout.data[k-1];
				fout.data[(j<<1)] = fout.data[k];
			}
			
			ExOp::toRand(da_vec);
//			fout.show();
			for(k=0;;k++){
			
			da_tmp[0] = da_vec[0];
			da_vec[0] = da_vec[0] * fout.data[0] + da_vec[1] * ExOp::mktrju(fout.data[1]);
			da_tmp[1] = da_vec[0] * da_old_norm - da_tmp[0];// printf(" %e\t%e\n",da_tmp[0] , da_vec[0] );
			da_dev = da_tmp[1] * da_tmp[1];
 		    da_norm = da_vec[0]* da_vec[0]; 
			for(j=1;j<SIZE-1;j++){
				da_tmp[(j & 1)] = da_vec[j];
				da_vec[j] = da_tmp[((j & 1 ) ^ 1)]  * fout.data[(((j-1) << 1) | 1)]  + da_vec[j] * fout.data[(j<<1)] + da_vec[j+1] * ExOp::mktrju(fout.data[((j<<1) | 1)]);
				da_tmp[((j & 1 ) ^ 1)] = da_tmp[(j & 1)] - da_vec[j]* da_old_norm; // printf(" %e\t%e\n",da_tmp[(j & 1)] , da_vec[j] );
				da_dev += da_tmp[((j & 1 ) ^ 1)] * da_tmp[((j & 1 ) ^ 1)];
				da_norm += da_vec[j]* da_vec[j]; 
			}
				da_tmp[(j & 1)] = da_vec[j];
				da_vec[j] = da_tmp[((j & 1 ) ^ 1)] * fout.data[(((j-1) << 1) | 1)] + da_vec[j] * fout.data[(j<<1)];
			//	printf(" %e\t%e\n",da_tmp[(j & 1)] , da_vec[j] );
				da_tmp[(j & 1)] -= da_vec[j] * da_old_norm;
				da_dev += da_tmp[(j & 1)] * da_tmp[(j & 1)];
				da_norm += da_vec[j]* da_vec[j];
			//	printf("devdev %e\t%e\n",da_norm * ExCo<double>::epsilon(), da_dev );
				if ((da_norm * ExCo<double>::epsilon()  > da_dev)||(k > SIZE*SIZE)) break;
				da_old_norm = ExOp::mkinvintpow(da_norm, -2);
				da_vec *= da_old_norm; // ^ -0.5f
				
			}
			
			delete[](hv);
			delete[](buf);
			return ExOp::mkinvintpow(da_norm, 2); // sqrt
		}else if (SIZE == 2){
			Tuple<C , 2> eigen = getEigenValues();
			return eigen[0] > eigen[1] ? eigen[0] : eigen[1];
		}  else return(fout.data[0]);
		return(da_norm);	
	
	}
	
	
	
	LFHTEMP template<unsigned int SIZ2> C Trianglix<C, SIZE>::Xformed_inner_product( const Tuple<C, SIZ2>& other) const{
		C s;ExOp::toZero(s);
		C s2;
		unsigned int i,j,k;

        for(j=0,k=0;(j<SIZ2) &&  (j<SIZE);j++) {
			if (j>0) s2 = data[k++] * other[0];
			else ExOp::toZero(s2);
			for(i=1;i<j;i++) s2 += data[k++] * other[i];
			// s += other[j] * ((ExOp::mktrju(s2) + s2) + data[k++] * other[j]); WORKS FOR REAL ONLY!
			s +=  data[k++] * ExOp::mktrju(other[j]) * other[j] + (ExOp::mkrealproj(s2) * ExOp::mkrealproj(other[j]) - ExOp::mkimmaproj(s2) * ExOp::mkimmaproj(other[j])) *2.0f;
			
        }
        
		return s;
    }
	
	LFHTEMP C Trianglix<C, SIZE>::Xformed_inner_product( const Tuple<C, 0u>& other) const{
		C s;
		C s2;
		unsigned int i,j,k;
		if (other.size() == 0) {ExOp::toZero(s); return s;} 
		s = data[0] * ExOp::mktrju(other[0]) * other[0];
		if (ExOp::mktrju(other[0]) * other[0] < 0.0f) { fprintf(stderr,"Should not be negative!: "); ExOp::show(ExOp::mktrju(other[0]) * other[0],stderr); exit(1); }
		for(j=1,k=1;j<other.size();j++) {
			s2 = data[k++] * other[0];
			for(i=1;i<j;i++) s2 += data[k++] * other[i];
			s += 2.0f * other[j] * s2 + data[k++] * other[j]; // WORKS FOR REAL ONLY!
			// s +=  data[k++] * ExOp::mktrju(other[j]) * other[j] + (ExOp::mkrealproj(s2) * ExOp::mkrealproj(other[j]) - ExOp::mkimmaproj(s2) * ExOp::mkimmaproj(other[j])) *2.0f;
			if (ExOp::mktrju(other[j]) * other[j] < 0.0f) { fprintf(stderr,"Should not be negative!: "); ExOp::show(ExOp::mktrju(other[j]) * other[j],stderr); exit(1); }
		}
		return s;	
	}
	

/*
LFHTEMP C Trianglix<C, SIZE>::log_determinant() const{
Trianglix<C, SIZE> fout(*this);
   C log_prod;
    if (SIZE >2){
    C* buf = new C[SIZE*2];
    C* hv = new C[SIZE *(SIZE+1)/2];
    double* normbuf = new double[SIZE-2];
    unsigned int i,j,k;
      for(j=SIZE-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mktrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mktrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mktrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]);normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }

        ExOp::toZero(hv[0]);
        if (fout.data[0] > hv[0]) hv[0] = fout.data[0];
        buf[SIZE] = ExOp::mkinverse(fout.data[0]); 
        for(j=1;j<SIZE;j++){
            k = (j*(j+3))/2;
            buf[SIZE+j-1] *= fout.data[k-1];
            fout.data[k] -= ExOp::mktrju(fout.data[k-1]) * buf[SIZE+j-1];
            log_prod  += log(fout.data[k]);
            if (fout.data[k] > hv[0]) hv[0] = fout.data[k];
            buf[SIZE+j] = ExOp::mkinverse(fout.data[k]);
            }
        k =0;
        log_prod = (fout.data[k] > 0.0f) ? log(fout.data[k]) : log(hv[0]) + log(ExCo<double>::epsilon());
        for(j=1;j<SIZE;j++){k = (j*(j+3))/2;
        log_prod += (fout.data[k] > 0.0f) ? log(fout.data[k]) : log(hv[0]) + log(ExCo<double>::epsilon());
        }
        
    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (SIZE == 2) return(log(data[0] *data[2] - data[1] * ExOp::mktrju(data[1])));
    else if (SIZE == 1) return(log(data[0]));
    return(log_prod);

} */

LFHTEMP Trianglix<C, SIZE> Trianglix<C, SIZE>::inverse() const{
    Trianglix<C, SIZE> fout(*this);
   
    if (SIZE >2){
    C* buf = new C[SIZE*2];
    C* hv = new C[SIZE *(SIZE+1)/2];
    double* normbuf = new double[SIZE-2];
    unsigned int i,j,k;
      for(j=SIZE-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mktrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mktrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mktrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]);  normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }

        fout.data[0] = ExOp::mkinverse(fout.data[0]);
        buf[SIZE] = fout.data[0];

        for(j=1;j<SIZE;j++){
            k = (j*(j+3))/2;
            buf[SIZE+j-1] *= fout.data[k-1];
            fout.data[k] = ExOp::mkinverse(fout.data[k] - ExOp::mktrju(fout.data[k-1]) * buf[SIZE+j-1]);
            ExOp::toZero(fout.data[k-1]);
            buf[SIZE+j] = fout.data[k];
            
            }


        for(j=SIZE-2;j<SIZE;j--){
            fout.offdiagelimination_up( -buf[SIZE+j],j,buf);
            }


    for(j=2;j<SIZE;j++) fout.HouseHolderMultiply(hv+((j * (j+1))/2), normbuf[j-2], j+1, buf,false);
        
    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (SIZE == 2){
        C denum = ExOp::mkinverse(data[0] * data[2] - data[1] * ExOp::mktrju(data[1]));
        fout.data[0] = data[2] * denum;
        fout.data[1] = -data[1] * denum;
        fout.data[2] = data[0] * denum;

    }  else fout.data[0] =  ExOp::mkinverse(data[0]);
    return(fout);
}

// fout * fout' = inverse
LFHTEMP TMatrix<C,SIZE,SIZE> Trianglix<C, SIZE>::diagonalizer_of_inverse() const{ TMatrix<C,SIZE,SIZE> f_out;
    ExOp::toOne(f_out);
     
    Trianglix<C, SIZE> fout(*this);
    if (SIZE >1){
    C* buf = new C[SIZE*2];
    C* hv = new C[SIZE *(SIZE+1)/2];
    double* normbuf = new double[SIZE-2];
    unsigned int i,j,k;
      for(j=SIZE-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(fout.data[k]);hv[k] = ExOp::mktrju(fout.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(fout.data[k+i]); hv[k+i] = ExOp::mktrju(fout.data[k+i]);}
            hv[k+i] = ExOp::mktrju(fout.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]);  normbuf[j-2] *=0.5f;
            fout.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
            f_out.LeftHouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
        }

        fout.data[0] = ExOp::mkinverse(fout.data[0]);
        buf[SIZE] = fout.data[0];

        for(j=1;j<SIZE;j++){
            k = (j*(j+3))/2;
             buf[SIZE+j-1] *= fout.data[k-1];
            fout.data[k] = ExOp::mkinverse(fout.data[k] - ExOp::mktrju(fout.data[k-1]) * buf[SIZE+j-1]);
            ExOp::toZero(fout.data[k-1]);
            buf[SIZE+j] = fout.data[k];
            }


        for(j=SIZE-2;j<SIZE;j--){
            fout.offdiagelimination_up( -buf[SIZE+j],j,buf);
            }

        
    delete[](hv);
    delete[](buf);
    delete[](normbuf);       
    }else f_out.data[0] = pow(data[0], -0.5f);
    return(f_out);
    }

LFHTEMP double Trianglix<C, SIZE>::bhattacharryya_partial(const Tuple<C,SIZE> &dev ,const Trianglix<C,SIZE>& other)const{
    double fout;
    // dev (P+P2/2)-1 dev 
    C sum,tmp;
    C deter[2];
    int safe_mag =0;
    Trianglix<C, SIZE> invsc = (*this) + other;
    Tuple<C,SIZE> xdev = dev;
    if (SIZE >1){
    C* buf = new C[SIZE*2];
    C* hv = new C[SIZE *(SIZE+1)/2];
    double* normbuf = new double[SIZE-2];
    unsigned int i,j,k;
      for(j=SIZE-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(invsc.data[k]);hv[k] = ExOp::mktrju(invsc.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(invsc.data[k+i]); hv[k+i] = ExOp::mktrju(invsc.data[k+i]);}
            hv[k+i] = ExOp::mktrju(invsc.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]);  normbuf[j-2] *=0.5f;
            invsc.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
            xdev.HouseHolderMultiply(hv+k, normbuf[j-2], j);
        }

       

        deter[0] = invsc.data[0];
        ExOp::toOne(deter[1]);

         buf[SIZE] = ExOp::mkinverse(invsc.data[0]);
        for(j=1;j<SIZE;j++){
            k = (j*(j+3))/2;

            tmp = deter[((j & 1)^1)] * invsc.data[k];
            tmp -= deter[(j & 1)] * (ExOp::mktrju(invsc.data[k-1]) * invsc.data[k-1]);  
            if (!(ExOp::isValid(tmp))){ // overflow!
               deter[0] *= pow(0.5f, 300.0f); 
               deter[1] *= pow(0.5f, 300.0f); 
                tmp = deter[((j & 1)^1)] * invsc.data[k];
                tmp -= deter[(j & 1)] * (ExOp::mktrju(invsc.data[k-1]) * invsc.data[k-1]);  
                safe_mag++;
            }else if (fabs(tmp) <= pow(0.5f, 300.0f)) { // underflow!
               deter[0] *= pow(2.0f, 300.0f); 
               deter[1] *= pow(2.0f, 300.0f); 
               tmp = deter[((j & 1)^1)] * invsc.data[k];
                tmp -= deter[(j & 1)] * (ExOp::mktrju(invsc.data[k-1]) * invsc.data[k-1]);  
               safe_mag--;
            }
            deter[(j & 1)]  = tmp;
    




            xdev[j] -= ExOp::mktrju(xdev[j-1]) * buf[SIZE+j-1]*invsc.data[k-1];
            if (j==1) sum = ExOp::mktrju(xdev[0]) * xdev[0] * buf[SIZE]; 
            else sum += ExOp::mktrju(xdev[j-1]) * xdev[j-1] * buf[SIZE+j-1]; 
            invsc.data[k] -= ExOp::mktrju(invsc.data[k-1]) * buf[SIZE+j-1]*invsc.data[k-1];
            buf[SIZE+j] = ExOp::mkinverse(invsc.data[k]);
       //     sum += ExOp::mktrju(xdev[j]) * xdev[j] * buf[SIZE+j];
            }
        sum += ExOp::mktrju(xdev[j-1]) * xdev[j-1] * buf[SIZE+j-1];
		
        fout = 0.5f *( 0.5f * ExOp::norm(sum) + log(ExOp::norm(deter[((SIZE & 1)^1)])) + (300.0f * safe_mag - SIZE) *log(2.0f) );


        
    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else fout = 0.5f *(ExOp::mktrju(dev[0]) * dev[0] * ExOp::mkinverse(data[0] + other.data[0]) + log(0.5f) + log(data[0] + other.data[0]));


    return(fout);
    }

LFHTEMP C Trianglix<C, SIZE>::Xformed_inner_product_of_inverse(const Tuple<C, SIZE> &dev) const{
    double fout;
    // dev (P+P2/2)-1 dev 
    C sum;
    C prod;
    Trianglix<C, SIZE> invsc = (*this);
    Tuple<C, SIZE> xdev = dev;
    if (SIZE >2){
    C* buf = new C[SIZE*2];
    C* hv = new C[SIZE *(SIZE+1)/2];
    double* normbuf = new double[SIZE-2];
    unsigned int i,j,k;
      for(j=SIZE-1;j>1;j--){
            k = (j * (j+1))/2;
            normbuf[j-2] = ExOp::pnorm(invsc.data[k]);hv[k] = ExOp::mktrju(invsc.data[k]);
            for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(invsc.data[k+i]); hv[k+i] = ExOp::mktrju(invsc.data[k+i]);}
            hv[k+i] = ExOp::mktrju(invsc.data[k+i]);
            hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
            normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
            invsc.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
            xdev.HouseHolderMultiply(hv+k, normbuf[j-2], j);
        }

        buf[SIZE] = ExOp::mkinverse(invsc.data[0]);
        prod =  invsc.data[0];
        

        ExOp::toZero(hv[0]);
        if (invsc.data[0] > hv[0]) hv[0] = invsc.data[0];
        for(j=1;j<SIZE;j++){
            k = (j*(j+3))/2;
            xdev[j] -= ExOp::mktrju(xdev[j-1]) * buf[SIZE+j-1]*invsc.data[k-1];
            if (j==1) sum = ExOp::mktrju(xdev[0]) * xdev[0] * buf[SIZE]; 
            else sum += ExOp::mktrju(xdev[j-1]) * xdev[j-1] * buf[SIZE+j-1]; 
            invsc.data[k] -= ExOp::mktrju(invsc.data[k-1]) * buf[SIZE+j-1]*invsc.data[k-1];
            if (invsc.data[k] > hv[0]) hv[0] = invsc.data[k];
            buf[SIZE+j] = ExOp::mkinverse(invsc.data[k]);
       //     sum += ExOp::mktrju(xdev[j]) * xdev[j] * buf[SIZE+j];
            }
        sum += ExOp::mktrju(xdev[j-1]) * xdev[j-1] * buf[SIZE+j-1];
        fout = sum;

        
    delete[](hv);
    delete[](buf);
    delete[](normbuf);
    }else if (SIZE == 2){
       fout = ExOp::mktrju(dev[1]) * dev[0] * data[1];
       fout = (ExOp::mktrju(dev[0]) * dev[0]* data[2] + ExOp::mktrju(dev[1]) * dev[1]* data[0] - fout - ExOp::mktrju(fout) ) * ExOp::mkinverse(data[0] * data[2] - data[1] * ExOp::mktrju(data[1]));
    }  else fout = ExOp::mktrju(dev[0]) * dev[0] * ExOp::mkinverse(data[0]);


    return(fout);
    }

LFHTEMP C Trianglix<C, SIZE>::inv_Xformed_inner_product_singularguard(const Tuple<C,SIZE> &dev, double guard_fraction, double *ln_det) const{
        C fout;
		// dev (P+P2/2)-1 dev
        unsigned int count;
		C sum,tmp;
		C deter[2];
        C maxeigen; ExOp::toZero(maxeigen);
		int safe_mag =0;
		Trianglix<C, SIZE> invsc = (*this);
		Tuple<C,SIZE> xdev = dev;
		if (SIZE >1){
			C* buf = new C[SIZE*2];
			C* hv = new C[SIZE *(SIZE+1)/2];
			double* normbuf = new double[SIZE-2];
			unsigned int i,j,k;
			for(j=SIZE-1;j>1;j--){
				k = (j * (j+1))/2;
				normbuf[j-2] = ExOp::pnorm(invsc.data[k]);hv[k] = ExOp::mktrju(invsc.data[k]);
				for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(invsc.data[k+i]); hv[k+i] = ExOp::mktrju(invsc.data[k+i]);}
				hv[k+i] = ExOp::mktrju(invsc.data[k+i]);
				hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
				normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]);  normbuf[j-2] *=0.5f;
				invsc.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
				xdev.HouseHolderMultiply(hv+k, normbuf[j-2], j);
			}
			
			
			
			deter[0] = invsc.data[0];
			ExOp::toOne(deter[1]);
			
			buf[SIZE] = invsc.data[0];
            maxeigen = fabs(buf[SIZE]);
			for(j=1;j<SIZE;j++){
				k = (j*(j+3))/2;
				
				/*
				tmp = deter[((j & 1)^1)] * invsc.data[k];
				tmp -= deter[(j & 1)] * (ExOp::mktrju(invsc.data[k-1]) * invsc.data[k-1]);  
				if (!(ExOp::isValid(tmp))){ // overflow!
					deter[0] *= pow(0.5f, 300.0f); 
					deter[1] *= pow(0.5f, 300.0f); 
					tmp = deter[((j & 1)^1)] * invsc.data[k];
					tmp -= deter[(j & 1)] * (ExOp::mktrju(invsc.data[k-1]) * invsc.data[k-1]);  
					safe_mag++;
				}else if (fabs(tmp) <= pow(0.5f, 300.0f)) { // underflow!
					deter[0] *= pow(2.0f, 300.0f); 
					deter[1] *= pow(2.0f, 300.0f); 
					tmp = deter[((j & 1)^1)] * invsc.data[k];
					tmp -= deter[(j & 1)] * (ExOp::mktrju(invsc.data[k-1]) * invsc.data[k-1]);  
					safe_mag--;
				}
				deter[(j & 1)]  = tmp;*/
								
				xdev[j] -= ExOp::mktrju(xdev[j-1]) * ExOp::mkinverse(buf[SIZE+j-1])*invsc.data[k-1];
				invsc.data[k] -= ExOp::mktrju(invsc.data[k-1]) * buf[SIZE+j-1]*invsc.data[k-1];
				buf[SIZE+j] = invsc.data[k];
                if (fabs(buf[SIZE+j]) > maxeigen) maxeigen = fabs(buf[SIZE+j]); 
            }
			
            deter[0] = 0.0f;
            for(j=0;j<SIZE;j++){
                if (fabs(buf[SIZE+j]) > maxeigen * ExCo<C>::epsilon()){
                    deter[0] += log(fabs(buf[SIZE+j]));
                    count++;
                }
            }
            if (count == 0) count =SIZE;
            if ((guard_fraction == 0.0f) || (count == SIZE)){
				for(j=0;j<SIZE;j++) ExOp::toinverse(buf[SIZE+j]);
				if (ln_det != NULL) (*ln_det) = deter[0] * (((double)SIZE) / ((double)count));
				
			}else{

			for(j=0;j<SIZE;j++){
                if (fabs(buf[SIZE+j]) > maxeigen * ExCo<C>::epsilon()) ExOp::toinverse(buf[SIZE+j]);
                else buf[SIZE+j] = exp( -deter[0] / count) / guard_fraction;
            }       
            if (ln_det != NULL) (*ln_det) = deter[0] * (((double)SIZE) / ((double)count));
				
			}
            fout = ExOp::mktrju(xdev[0]) * xdev[0] * buf[SIZE]; 
            for(j=1;j<SIZE;j++) fout += ExOp::mktrju(xdev[j]) * xdev[j] * buf[SIZE+j];
			
		//	fout = (this_weight + other_weight) * 0.5f *( ExOp::norm(sum) + log(ExOp::norm(deter[((SIZE & 1)^1)])) + log(pow(2.0f,300.0f)) * safe_mag);
			
			
			
			delete[](hv);
			delete[](buf);
			delete[](normbuf);
		}else {
            fout = ExOp::mktrju(dev[0]) * dev[0] * ExOp::mkinverse(data[0]);
            if (ln_det != NULL) *ln_det = log(data[0]);
        }
		
		
		return(fout);
    }

LFHTEMP double Trianglix<C, SIZE>::weighted_bhattacharryya_partial(const Tuple<C,SIZE> &dev ,const Trianglix<C,SIZE>& other, double this_weight, double other_weight)const{
		double fout;
		// dev (P+P2/2)-1 dev 
		C sum,tmp;
		C deter[2];
		int safe_mag =0;
		Trianglix<C, SIZE> invsc = ((*this) * (this_weight / (this_weight + other_weight))) + (other* (other_weight / (this_weight + other_weight))) ;
		Tuple<C,SIZE> xdev = dev;
		if (SIZE >1){
			C* buf = new C[SIZE*2];
			C* hv = new C[SIZE *(SIZE+1)/2];
			double* normbuf = new double[SIZE-2];
			unsigned int i,j,k;
			for(j=SIZE-1;j>1;j--){
				k = (j * (j+1))/2;
				normbuf[j-2] = ExOp::pnorm(invsc.data[k]);hv[k] = ExOp::mktrju(invsc.data[k]);
				for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(invsc.data[k+i]); hv[k+i] = ExOp::mktrju(invsc.data[k+i]);}
				hv[k+i] = ExOp::mktrju(invsc.data[k+i]);
				hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
				normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
				invsc.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
				xdev.HouseHolderMultiply(hv+k, normbuf[j-2], j);
			}
			
			
			
			deter[0] = invsc.data[0];
			ExOp::toOne(deter[1]);
			
			buf[SIZE] = ExOp::mkinverse(invsc.data[0]);
			for(j=1;j<SIZE;j++){
				k = (j*(j+3))/2;
				
				tmp = deter[((j & 1)^1)] * invsc.data[k];
				tmp -= deter[(j & 1)] * (ExOp::mktrju(invsc.data[k-1]) * invsc.data[k-1]);  
				if (!(ExOp::isValid(tmp))){ // overflow!
					deter[0] *= pow(0.5f, 300.0f); 
					deter[1] *= pow(0.5f, 300.0f); 
					tmp = deter[((j & 1)^1)] * invsc.data[k];
					tmp -= deter[(j & 1)] * (ExOp::mktrju(invsc.data[k-1]) * invsc.data[k-1]);  
					safe_mag++;
				}else if (fabs(tmp) <= pow(0.5f, 300.0f)) { // underflow!
					deter[0] *= pow(2.0f, 300.0f); 
					deter[1] *= pow(2.0f, 300.0f); 
					tmp = deter[((j & 1)^1)] * invsc.data[k];
					tmp -= deter[(j & 1)] * (ExOp::mktrju(invsc.data[k-1]) * invsc.data[k-1]);  
					safe_mag--;
				}
				deter[(j & 1)]  = tmp;
				
				
				
				
				
				xdev[j] -= ExOp::mktrju(xdev[j-1]) * buf[SIZE+j-1]*invsc.data[k-1];
				if (j==1) sum = ExOp::mktrju(xdev[0]) * xdev[0] * buf[SIZE]; 
				else sum += ExOp::mktrju(xdev[j-1]) * xdev[j-1] * buf[SIZE+j-1]; 
				invsc.data[k] -= ExOp::mktrju(invsc.data[k-1]) * buf[SIZE+j-1]*invsc.data[k-1];
				buf[SIZE+j] = ExOp::mkinverse(invsc.data[k]);
				//     sum += ExOp::mktrju(xdev[j]) * xdev[j] * buf[SIZE+j];
            }
			sum += ExOp::mktrju(xdev[j-1]) * xdev[j-1] * buf[SIZE+j-1];
			
			fout = (this_weight + other_weight) * 0.5f *( ExOp::norm(sum) + log(ExOp::norm(deter[((SIZE & 1)^1)])) + log(pow(2.0f,300.0f)) * safe_mag);
			
			
			
			delete[](hv);
			delete[](buf);
			delete[](normbuf);
		}else fout = (this_weight + other_weight) * 0.5f *(ExOp::mktrju(dev[0]) * dev[0] * ExOp::mkinverse(data[0] + other.data[0]) + log(0.5f) + log(data[0] + other.data[0]));
		
		
		return(fout);
    }
	LFHTEMP double Trianglix<C, SIZE>::weighted_bhattacharryya_partial_MK2(const Tuple<C,SIZE> &dev ,const Trianglix<C,SIZE>& other, double this_weight, double other_weight)const{
		double fout;
		// dev (P+P2/2)-1 dev 
		C sum,tmp;
		C deter[2];
		int safe_mag =0;
		Trianglix<C, SIZE> invsc = ((*this) * (this_weight / (this_weight + other_weight))) + (other* (other_weight / (this_weight + other_weight))) ;
		Tuple<C,SIZE> xdev = dev;
		if (SIZE >1){
			C* buf = new C[SIZE*2];
			C* hv = new C[SIZE *(SIZE+1)/2];
			double* normbuf = new double[SIZE-2];
			unsigned int i,j,k;
			for(j=SIZE-1;j>1;j--){
				k = (j * (j+1))/2;
				normbuf[j-2] = ExOp::pnorm(invsc.data[k]);hv[k] = ExOp::mktrju(invsc.data[k]);
				for(i=1;i<j-1;i++) {normbuf[j-2] += ExOp::pnorm(invsc.data[k+i]); hv[k+i] = ExOp::mktrju(invsc.data[k+i]);}
				hv[k+i] = ExOp::mktrju(invsc.data[k+i]);
				hv[k+i] *= 1.0f - sqrt(1.0f + normbuf[j-2] / ExOp::pnorm(hv[k+i]));
				normbuf[j-2] += ExOp::pnorm(hv[k+i]);ExOp::toZero(hv[k+j]); normbuf[j-2] *=0.5f;
				invsc.HouseHolderMultiply(hv+k, normbuf[j-2], j+1, buf,true);
				xdev.HouseHolderMultiply(hv+k, normbuf[j-2], j);
			}
			
			
			
			deter[0] = invsc.data[0];
			ExOp::toOne(deter[1]);
			
			buf[SIZE] = ExOp::mkinverse(invsc.data[0]);
			for(j=1;j<SIZE;j++){
				k = (j*(j+3))/2;
				
				tmp = deter[((j & 1)^1)] * invsc.data[k];
				tmp -= deter[(j & 1)] * (ExOp::mktrju(invsc.data[k-1]) * invsc.data[k-1]);  
				if (!(ExOp::isValid(tmp))){ // overflow!
					deter[0] *= pow(0.5f, 300.0f); 
					deter[1] *= pow(0.5f, 300.0f); 
					tmp = deter[((j & 1)^1)] * invsc.data[k];
					tmp -= deter[(j & 1)] * (ExOp::mktrju(invsc.data[k-1]) * invsc.data[k-1]);  
					safe_mag++;
				}else if (fabs(tmp) <= pow(0.5f, 300.0f)) { // underflow!
					deter[0] *= pow(2.0f, 300.0f); 
					deter[1] *= pow(2.0f, 300.0f); 
					tmp = deter[((j & 1)^1)] * invsc.data[k];
					tmp -= deter[(j & 1)] * (ExOp::mktrju(invsc.data[k-1]) * invsc.data[k-1]);  
					safe_mag--;
				}
				deter[(j & 1)]  = tmp;
				
				
				
				
				
				xdev[j] -= ExOp::mktrju(xdev[j-1]) * buf[SIZE+j-1]*invsc.data[k-1];
				if (j==1) sum = ExOp::mktrju(xdev[0]) * xdev[0] * buf[SIZE]; 
				else sum += ExOp::mktrju(xdev[j-1]) * xdev[j-1] * buf[SIZE+j-1]; 
				invsc.data[k] -= ExOp::mktrju(invsc.data[k-1]) * buf[SIZE+j-1]*invsc.data[k-1];
				buf[SIZE+j] = ExOp::mkinverse(invsc.data[k]);
				//     sum += ExOp::mktrju(xdev[j]) * xdev[j] * buf[SIZE+j];
            }
			sum += ExOp::mktrju(xdev[j-1]) * xdev[j-1] * buf[SIZE+j-1];
			
			fout = (this_weight + other_weight) * 0.125f * ExOp::norm(sum) + 0.5f * log(ExOp::norm(deter[((SIZE & 1)^1)])) + 150.0f * log(2.0f) * safe_mag;
			
			
			
			delete[](hv);
			delete[](buf);
			delete[](normbuf);
		}else fout = (this_weight + other_weight) * 0.5f *(ExOp::mktrju(dev[0]) * dev[0] * ExOp::mkinverse(data[0] + other.data[0]) + log(0.5f) + log(data[0] + other.data[0]));
		
		
		return(fout);
    }
	
LFHTEMP void Trianglix<C, SIZE>::show(FILE* f, int level)const{
    unsigned int i,k,j; // fprintf(f,"%u\t%u\n",totsize, sizeof(data) / sizeof(C) );
    for(i=0,j=0,k=2;i<totsize;i++) {ExOp::show(data[i],f,level+1); if (i==j) {fprintf(f,"\n"); j+= k++;} else fprintf(f,"\t");}
}
LFHTEMP string Trianglix<C, SIZE>::type_tostring() const{char buffer[256]; sprintf(buffer,"%u",SIZE); return string("Trianglix<") + ExOp::type_tostring(data[0]) +string(",")+ ExOp::type_tostring(buffer) + string(">");}

	
LFHTEMP double Trianglix<C, SIZE>::WishartLogDensity(const Trianglix<C, SIZE>& var, unsigned int nbsamples) const{
	double fout;
	
	unsigned int i;
	fout = ((double)SIZE * nbsamples) * log(0.5f);
	fout -= ((double)nbsamples) * this->log_determinant();
	fout += ((double)(nbsamples - SIZE - 1)) * var.log_determinant();
	fout -= var.trace_of_division(*this);
	
	fout *=0.5f;
	for(i=0;i<SIZE;i++) fout -= lngamma( ((double)(nbsamples -i)) *0.5f);
	return fout;
}
	
/*
LFHTEMP	template<class A>
const Quaternion<C>& Quaternion<C>::operator*=(Quaternion<A> const & other){
		C tmp[3];
		tmp[0] = (*this)[0];
		tmp[1] = (*this)[1];
		tmp[2] = (*this)[2];
		(*this)[0] = (*this)[0]*other[0] - (*this)[1]*other[1] -(*this)[2]*other[2] -(*this)[3]*other[3];
		(*this)[1] = (*this)[1]*other[0] + tmp[0]*other[1] -(*this)[3]*other[2] +(*this)[2]*other[3];
		(*this)[2] = (*this)[2]*other[0] + (*this)[3]*other[1] +tmp[0]*other[2] -tmp[1]*other[3];
		(*this)[3] = (*this)[3]*other[0] - tmp[2]*other[1] +tmp[1]*other[2] +tmp[0]*other[3];
		return(*this);
	}
	*/

 // end of namespace







