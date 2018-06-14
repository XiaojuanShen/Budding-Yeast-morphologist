/*
 * primitive_tem.h
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

	template<unsigned int order> double cheb_eval(const Tuple<double, order> cs, const double x){
		int j;
		double d  = 0.0;
		double dd = 0.0;
		double y  = 2.0*x;
		double t;
		for(j = order; j>=1; j--) {
			t = d;
			d = y*d - dd + cs[j];
			dd =t;
		}
		d = (0.5 *(y*d + cs[0])) - dd;
		return d;
	}



#undef LFHTEMP
#define LFHTEMP template <class C>


LFHTEMP void HeapTree<C>::insert(const C &what){
	if (hasunsorted){

	int cur = data.getSize();
	data.push_back();
	//push_back(KeyElem<Key, Node>());

	while(true){
	if ((cur == 1)||( data[cur >>1] < what )) break;
	ExOp::memmove(data[cur], data[cur  >>1]);
	cur  >>= 1;
	}
	data[cur] = what;
	}else{
		hasunsorted = true;
		data[0] = what; // quick insert!
	}
	}

LFHTEMP C HeapTree<C>::pop(){
	C swap;
	unsigned int s = data.getSize();
	if (!hasunsorted){
        if (s == 1) {fprintf(stderr,"poping form an empty heap! (required to check first)\n");exit(1);}
		data[0] = data[1];
		swap = *(data.last());
		data.pop_back();s--;
		if (s == 1) return(data[0]);
	}else{
		hasunsorted = false;
		if ((s ==1)||(data[0] < data[1])) {
			return(data[0]);
		}else{
			swap = data[0];
			data[0] = data[1];
		}
        
	}

	unsigned int cur = 2;
	while(true){
		if (cur+1 < s) {
			if (data[cur+1] < data[cur]) cur++;
		}else{
			if (cur+1 !=  s) break;
		}
		if (data[cur] > swap) break;
		data[cur >> 1] = data[cur];
		cur <<=1;
	}
	data[cur >> 1] = swap;

	return(data[0]);
	}

	// pointer specialization!!!
#undef LFHTEMP
#define LFHTEMP template <class C>
	// pointer specialization!!!

	LFHTEMP HeapTree<C*>::HeapTree() : hasunsorted(false){
		this->setSize(1);
	}


	LFHTEMP HeapTree<C*>::HeapTree(void * owner) : Vector< C* >(owner), hasunsorted(false){
		this->setSize(1);
	}

	LFHTEMP bool HeapTree<C*>::isempty(){
		return((!hasunsorted) && (this->size() == 1));
	}

	//template <class Key, class Node>
	//Node& HeapTree<Key,Node>::operator()(const Key & where){ // yet a method to insert

	//	}
	LFHTEMP void HeapTree<C*>::insert(C* what){
		unsigned int cur;
		if (hasunsorted){

			cur = this->size();
			this->push_back();
			//push_back(KeyElem<Key, Node>());

			while(true){
				if ((cur == 1)||( *(*this)[cur >>1] < *what )) break;
				(*this)[cur] = (*this)[cur  >>1];
				if (HeapTreeBackPointerOffset<C>::ans) HeapTreeBackPointerOffset<C>::getBackOffset(* (*this)[cur] ) = cur;
				cur  >>= 1;
			}
		}else{
			hasunsorted = true;
			cur = 0;// quick insert!
		}
		(*this)[cur] = what;

		if (HeapTreeBackPointerOffset<C>::ans) HeapTreeBackPointerOffset<C>::getBackOffset(*what) = cur;

	}


	LFHTEMP C& HeapTree<C*>::pop(){
		C* swap;

		if (!hasunsorted){
			(*this)[0] = (*this)[1];
			swap = *(this->last());
			this->pop_back();
		}else{
			hasunsorted = false;
			if (*(*this)[0] < *(*this)[1]) {
				return(*(*this)[0]);
			}else{
				swap = (*this)[0];
				(*this)[0] = (*this)[1];
			}

		}

		int s = this->size();
		int cur = 1;
		while((cur << 1) < s){
			if (((cur << 1) == s-1)||(*(*this)[cur<<1] < *(*this)[1 | (cur<< 1)])){
				if (*swap < *(*this)[cur<<1]) break;
				(*this)[cur] = (*this)[cur <<1];
				if (HeapTreeBackPointerOffset<C>::ans) HeapTreeBackPointerOffset<C>::getBackOffset(* (*this)[cur] ) = cur;
				cur = cur << 1;
			}else{
				if (*swap < *(*this)[1 | (cur<<1)]) break;
				(*this)[cur] = (*this)[1 | (cur <<1)];
				if (HeapTreeBackPointerOffset<C>::ans) HeapTreeBackPointerOffset<C>::getBackOffset(* (*this)[cur] ) = cur;
				cur = 1 | (cur << 1);
			}
		}

		(*this)[cur] = swap;
		if (HeapTreeBackPointerOffset<C>::ans) HeapTreeBackPointerOffset<C>::getBackOffset(* (*this)[cur] ) = cur;

		return(*((*this)[0]));
	}


	LFHTEMP void HeapTree<C*>::update(unsigned int cur){// Elem was changed! update its position!
		if (cur == 0) return; // unsorted position, ignore!
		C* swap;
		int s = this->size();
		if ((cur & 0xFFFFFFFE)&&(*(*this)[cur >> 1] > *(*this)[cur])){
			// got smaller
			swap = (*this)[cur];
			do{
			(*this)[cur] = (*this)[cur >> 1];
			if (HeapTreeBackPointerOffset<C>::ans) HeapTreeBackPointerOffset<C>::getBackOffset(* (*this)[cur] ) = cur;
			cur = cur >> 1;
			}while((cur & 0xFFFFFFFE)&&(*(*this)[cur >> 1] > *swap));
			(*this)[cur] = swap;
			if (HeapTreeBackPointerOffset<C>::ans) HeapTreeBackPointerOffset<C>::getBackOffset(* (*this)[cur] ) = cur;
			return;
		}
			if ((cur << 1) < s){
				if (((cur << 1) == s-1)||(*(*this)[cur << 1] < *(*this)[1 | (cur << 1)] )){
					if (*(*this)[cur << 1] < *(*this)[cur]){
						swap = (*this)[cur];
						(*this)[cur] = (*this)[cur << 1];
						if (HeapTreeBackPointerOffset<C>::ans) HeapTreeBackPointerOffset<C>::getBackOffset(* (*this)[cur] ) = cur;
						cur = cur << 1;
					}else return;
				}else{
					if (*(*this)[1+(cur << 1)] < *(*this)[cur]){
						swap = (*this)[cur];
						(*this)[cur] = (*this)[1 + (cur << 1)];
						if (HeapTreeBackPointerOffset<C>::ans) HeapTreeBackPointerOffset<C>::getBackOffset(* (*this)[cur] ) = cur;
						cur = 1 | (cur << 1);
					}else return;
				}
			} else return;
			while((cur << 1) < s){
				if (((cur << 1) == s-1)||(*(*this)[cur << 1] < *(*this)[1 | (cur << 1)] )){
					if (*(*this)[cur << 1] < *swap){
						(*this)[cur] = (*this)[cur << 1];
						if (HeapTreeBackPointerOffset<C>::ans) HeapTreeBackPointerOffset<C>::getBackOffset(* (*this)[cur] ) = cur;
						cur = cur << 1;
					}else break;
				}else{
					if (*(*this)[1+(cur << 1)] < *swap){
						(*this)[cur] = (*this)[1 + (cur << 1)];
						if (HeapTreeBackPointerOffset<C>::ans) HeapTreeBackPointerOffset<C>::getBackOffset(* (*this)[cur] ) = cur;
						cur = 1 | (cur << 1);
					}else break;
				}
			}
			(*this)[cur] = swap;
			if (HeapTreeBackPointerOffset<C>::ans) HeapTreeBackPointerOffset<C>::getBackOffset(* (*this)[cur] ) = cur;

	}


#undef LFHTEMP
#define LFHTEMP template <class C,unsigned int TSIZE,Tuple_flag Cflag> template <class O>

LFHTEMP Tuple<typename STDRETTYPE2<C,O>::PLUS_TYPE ,TSIZE,Cflag>
	Tuple<C,TSIZE,Cflag>::operator+(Tuple<O,TSIZE,Cflag> const & other) const{
		Tuple<typename STDRETTYPE2<C,O>::PLUS_TYPE ,TSIZE,Cflag> f_out;
		for(unsigned int i=0;i<TSIZE;i++) f_out[i] = (*this)[i] + other[i];
		return( f_out );
	}
	LFHTEMP Tuple<typename STDRETTYPE2<C,O>::MINU_TYPE ,TSIZE,Cflag>
	Tuple<C,TSIZE,Cflag>::operator-(Tuple<O,TSIZE,Cflag> const & other) const{
		Tuple<typename STDRETTYPE2<C,O>::MINU_TYPE ,TSIZE,Cflag> f_out;
		for(unsigned int i=0;i<TSIZE;i++) f_out[i] = (*this)[i] - other[i];
		return( f_out );
	}
	LFHTEMP Tuple<typename STDRETTYPE2<C,O>::PROD_TYPE ,TSIZE,Cflag>
	Tuple<C,TSIZE,Cflag>::operator*(Tuple<O,TSIZE,Cflag> const & other) const{
		Tuple<typename STDRETTYPE2<C,O>::PROD_TYPE ,TSIZE,Cflag> f_out;
		for(unsigned int i=0;i<TSIZE;i++) f_out[i] = (*this)[i] * other[i];
		return( f_out );
	}
	LFHTEMP Tuple<typename STDRETTYPE2<C,O>::DIVI_TYPE ,TSIZE,Cflag>
	Tuple<C,TSIZE,Cflag>::operator/(Tuple<O,TSIZE,Cflag> const & other) const{
		Tuple<typename STDRETTYPE2<C,O>::DIVI_TYPE ,TSIZE,Cflag> f_out;
		for(unsigned int i=0;i<TSIZE;i++) f_out[i] = (*this)[i] / other[i];
		return( f_out );
	}

LFHTEMP Tuple<typename STDRETTYPE2<C,O>::PLUS_TYPE ,TSIZE,Cflag>
	Tuple<C,TSIZE,Cflag>::operator+(O const & other) const{
		Tuple<typename STDRETTYPE2<C,O>::PLUS_TYPE ,TSIZE,Cflag> f_out;
		for(unsigned int i=0;i<TSIZE;i++) f_out[i] = (*this)[i] + other;
		return( f_out );
	}
LFHTEMP Tuple<typename STDRETTYPE2<C,O>::MINU_TYPE ,TSIZE,Cflag>
	Tuple<C,TSIZE,Cflag>::operator-(O const & other) const{
		Tuple<typename STDRETTYPE2<C,O>::MINU_TYPE ,TSIZE,Cflag> f_out;
		for(unsigned int i=0;i<TSIZE;i++) f_out[i] = (*this)[i] - other;
		return( f_out );
	}
LFHTEMP Tuple<typename STDRETTYPE2<C,O>::PROD_TYPE ,TSIZE,Cflag>
	Tuple<C,TSIZE,Cflag>::operator*(O const & other) const{
		Tuple<typename STDRETTYPE2<C,O>::PROD_TYPE ,TSIZE,Cflag> f_out;
		for(unsigned int i=0;i<TSIZE;i++) f_out[i] = (*this)[i] * other;
		return( f_out );
	}
LFHTEMP Tuple<typename STDRETTYPE2<C,O>::DIVI_TYPE ,TSIZE,Cflag>
	Tuple<C,TSIZE,Cflag>::operator/(O const & other) const{
		Tuple<typename STDRETTYPE2<C,O>::DIVI_TYPE ,TSIZE,Cflag> f_out;
		for(unsigned int i=0;i<TSIZE;i++) f_out[i] = (*this)[i] / other;
		return( f_out );
	}



	template<class TYPE> void comp_FFT_routine(TYPE* data, unsigned char mag, unsigned int mult){
		// assumes size is a power of two, and C can be multiplied by a complex number

		// perorm a partial FFT on mult << mag, mag >=1
		if (mag == 0) return;
		TYPE * icur;
		TYPE * ocur;
		for(icur=data,ocur= data + (mult << (mag - 1)) ;(((unsigned int)(ocur - data))  >> mag) < mult;icur++, ocur++){
			TYPE tmp = (*icur);
			tmp -= (*ocur);
			(*icur) += (*ocur);
			(*ocur) = tmp;
		}


		if (mag == 1) return;

		// TODO


	}

	template<class TYPE> void comp_IFFT_routine(TYPE* data, unsigned char mag, unsigned int mult){
		// assumes size is a power of two, and C can be multiplied by a complex number

		if (mag == 0) return;

		TYPE * icur;
		TYPE * ocur;
		for(icur=data,ocur= data + (mult << (mag - 1)) ;(((unsigned int)(ocur - data))  >> mag) < mult;icur++, ocur++){
			TYPE tmp = (*icur);
			tmp -= (*ocur);
			(*icur) += (*ocur);
			(*ocur) = tmp;
		}

		if (mag == 1) return;

//		complex half = complex(1/((double)size),0.0f);
//		for(icur =0;icur<size;icur++) data[icur] *= half;

	}

template<class TYPE> void pow2_FFT_routine(TYPE* data, unsigned char mag){
		// assumes size is a power of two, and C can be multiplied by a complex number

		unsigned int icur;

        unsigned int step = 1 << (mag-1);

           for(icur=0;icur<step;icur++){
					TYPE tmp = data[icur];
					tmp -=data[(step | icur)];
                    data[ icur] += data[(step | icur)];
                    data[(step | icur)] = tmp;
                }


		if (mag == 1) return;

		complex damultip;
		double tmpdouble;
		unsigned int ang;
		unsigned char stmag;
		for(stmag = mag-2 ;stmag != 0xFF ;stmag--){
            // run
			step= 1 << stmag;
			for(icur=0;(icur >> mag) == 0;icur+= step){
				ang = icur;
		//		printf("%i\t", ang);
				ExCo<unsigned int>::bitreverse(ang);
		//		printf("%i\n", ang);
				if (mag > 31 - stmag) ang <<= mag + stmag- 31;
				else ang >>= 31 - stmag - mag;
                tmpdouble = ((double) M_PI*ang) / (1 << mag);
		//		printf("%i\t%f\n", ang, tmpdouble);
			    damultip = complex(cos(tmpdouble),sin(tmpdouble));
				for(;(icur & step) == 0 ; icur++){
                        data[ (step | icur)] *= damultip;
                        TYPE tmp = data[icur];
                        tmp -= data[ (step | icur)];
						data[icur] += data[(step | icur)];
						data[(step | icur)] = tmp;
				}
			}
		}
	}

	template<class TYPE> void pow2_FFT_routine_zerohalf(TYPE* data, unsigned char mag){ // ignores the 2nd half of data, assumes it is zeros
		// assumes size is a power of two, and C can be multiplied by a complex number

		unsigned int icur;

        unsigned int step = 1 << (mag-1);

		for(icur=0;icur<step;icur++){
			data[(step | icur)] = data[icur];
		}


		if (mag == 1) return;

		complex damultip;
		double tmpdouble;
		unsigned int ang;
		unsigned char stmag;
		for(stmag = mag-2 ;stmag != 0xFF ;stmag--){
            // run
			step= 1 << stmag;
			for(icur=0;(icur >> mag) == 0;icur+= step){
				ang = icur;
				//		printf("%i\t", ang);
				ExCo<unsigned int>::bitreverse(ang);
				//		printf("%i\n", ang);
				if (mag > 31 - stmag) ang <<= mag + stmag- 31;
				else ang >>= 31 - stmag - mag;
                tmpdouble = ((double) M_PI*ang) / (1 << mag);
				//		printf("%i\t%f\n", ang, tmpdouble);
			    damultip = complex(cos(tmpdouble),sin(tmpdouble));
				for(;(icur & step) == 0 ; icur++){
					data[ (step | icur)] *= damultip;
					TYPE tmp = data[icur];
					tmp -= data[ (step | icur)];
					data[icur] += data[(step | icur)];
					data[(step | icur)] = tmp;
				}
			}
		}
	}

template<class TYPE> void pow2_IFFT_routine(TYPE* data, unsigned char mag){

	unsigned int icur;

	unsigned int step = 1 << (mag-1);

	for(icur=0;icur<step;icur++){
		TYPE tmp = data[icur];
		tmp -=data[(step | icur)];
		data[ icur] += data[(step | icur)];
		data[(step | icur)] = tmp;
	}


	if (mag == 1) return;

	complex damultip;
	double tmpdouble;
	unsigned int ang;
	unsigned char stmag;
	for(stmag = mag-2 ;stmag != 0xFF ;stmag--){
		// run
		step= 1 << stmag;
		for(icur=0;(icur >> mag) == 0;icur+= step){
			ang = icur;
			//		printf("%i\t", ang);
			ExCo<unsigned int>::bitreverse(ang);
			//		printf("%i\n", ang);
			if (mag > 31 - stmag) ang <<= mag + stmag- 31;
			else ang >>= 31 - stmag - mag;
			tmpdouble = ((double)-M_PI*ang) / (1 << mag);
			damultip = complex(cos(tmpdouble),sin(tmpdouble));
			for(;(icur & step) == 0 ; icur++){
				data[ (step | icur)] *= damultip;
				TYPE tmp = data[icur];
				tmp -= data[ (step | icur)];
				data[icur] += data[(step | icur)];
				data[(step | icur)] = tmp;
			}
		}
	}

		complex half = complex(1/((double)(1 << mag)),0.0f);
		for(icur =0;(icur >> mag) == 0;icur++) data[icur] *= half;
	}

	 // assumes input is swaped, output does not need swaping
	template<class TYPE> void pow2_IFFT_routine_swapped(TYPE* data, unsigned char mag){
		// assumes size is a power of two, and C can be multiplied by a complex number
        unsigned int step;
		unsigned int icur;
		unsigned int size = 1 << mag;
		unsigned int ang;
		unsigned char stmag;
		complex damultip;
		double tmpdouble;
		 if (mag > 1) {

		 for(step= 1,stmag=0 ; (step << 1) < size ;step=step << 1,stmag++){
		 // run

		 for(icur =0;icur<size;icur+= step){
			 ang = icur;
			 //		printf("%i\t", ang);
			 ExCo<unsigned int>::bitreverse(ang);
			 //		printf("%i\n", ang);
			 if (mag > 31 - stmag) ang <<= mag + stmag- 31;
			 else ang >>= 31 - stmag - mag;
			 tmpdouble = ((double)-M_PI*ang) / (1 << mag);
			 damultip = complex(cos(tmpdouble),sin(tmpdouble));

		 for(;(icur & step) == 0 ; icur++){
	     TYPE tmp = data[icur];
		 tmp -= data[ step | icur];
		 data[icur] += data[step | icur];
		 data[step | icur] = tmp;
			 data[ step | icur] *= damultip;

		 }
		 }


		 }
		 }
		 step = 1 << (mag-1);

		 for(icur=0;icur<step;icur++){
		 TYPE tmp = data[icur] - data[step | icur];
		 data[ icur] += data[step | icur];
		 data[step | icur] = tmp;
		 }

		complex half = complex(1/((double)(1 << mag)),0.0f);
		for(icur =0;(icur >> mag) == 0;icur++) data[icur] *= half;
	}

	// assumes input is swaped, output does not need swaping
	template<class TYPE> void pow2_IFFT_routine_swapped_shift(TYPE* data, unsigned char mag, unsigned int out_size){
		// assumes size is a power of two, and C can be multiplied by a complex number
        int step;
		int icur;
		unsigned int size = 1 << mag;
		unsigned int ang;
		unsigned char stmag;
		complex damultip;
		double tmpdouble;
		if (mag > 1) {

			for(step= 1,stmag=0 ; (step << 1) < size ;step=step << 1,stmag++){
				// run

				for(icur =0;icur<size;icur+= step){
					ang = icur;
					//		printf("%i\t", ang);
					ExCo<unsigned int>::bitreverse(ang);
					//		printf("%i\n", ang);
					if (mag > 31 - stmag) ang <<= mag + stmag- 31;
					else ang >>= 31 - stmag - mag;
					tmpdouble = ((double)-M_PI*ang) / (1 << mag);
					damultip = complex(cos(tmpdouble),sin(tmpdouble));

					for(;(icur & step) == 0 ; icur++){
						TYPE tmp = data[icur];
						tmp -= data[ step | icur];
						data[icur] += data[step | icur];
						data[step | icur] = tmp;
						data[ step | icur] *= damultip;

					}
				}


			}
		}else if (mag == 1) {data[1] += data[0]; data[1] *= 0.5f; return;}
		else return;
		step = 1 << (mag-1);
		complex half = complex(1/((double)(1 << mag)),0.0f);
		for(icur=out_size-1;icur != 0xFFFFFFFF;icur--){
			data[size- out_size + icur ] = (data[icur] + data[step | icur]) * half;
		}

	}


template<class TYPE> void pow2_bitswap_permutation(TYPE* a, unsigned char mag){ // partial
    unsigned int i=0;
    unsigned int j=0;
    unsigned char imag =0; // lead dif
    TYPE swp;
    unsigned char addc =0; // lead dif
    for(imag=0; imag < (mag >> 1); imag++){
            i = 1 << imag;
            j = 1 << (mag - 1 - imag);
            do{
            while(((i >> (mag - 1 - imag)) & 1) == 0){
                swp = a[i]; a[i] = a[j]; a[j] = swp; //  printf("%X\t%X\n", i,j);
                for(addc = imag+1; ; addc++) if (i & (1 << addc)) {i ^= (1 << addc); j ^= (1 << (mag - 1 - addc));} else{i |= (1 << addc); j |= (1 << (mag - 1 - addc)); break;}
            }
            i ^= (1 << (mag - 1 - imag));
            j ^= (1 << imag);
            for(addc = 0;addc < imag; addc++) if (i & (1 << addc)) {i ^= (1 << addc); j ^= (1 << addc); i^= (1 << (mag - 1 - addc)); j ^= (1 << (mag - 1 - addc));} else{i |= (1 << addc);j |= (1 << addc); i |= (1 << (mag - 1 - addc)); j |= (1 << (mag - 1 - addc)); break;}
            }while(addc < imag);
    }
}

	template<class TYPE> void pow2_BLUE_FFT_routine(TYPE* data, Complex<double> &blue, unsigned char mag, unsigned int subsize){
		// assumes size is a power of two, and C can be multiplied by a complex number


	}

	template<class TYPE> void pow2_BLUE_IFFT_routine(TYPE* data, Complex<double> &blue, unsigned char mag, unsigned int subsize){
		// assumes size is a power of two, and C can be multiplied by a complex number

	}



template<class C> void Bluestein::toFT(C* data) const{
	if (mult == 0){// pow2_FFT!
		pow2_FFT_routine(data, pre_pow2); //  1 <= mult <= 7
		pow2_bitswap_permutation(data, pre_pow2); //  1 <= mult <= 7
	}else if (pre_pow2 == 0){ // odd number! D:

		unsigned int i;
		for(i=0;i<mult;i++){
			double ang = i * i * M_PI / mult; //noneg!
			Complex<double> factor(cos(ang),-sin(ang));
			data[i] *= factor;
		}

		for(; (i>> post_pow2) == 0;i++) ExOp::toZero(data[i]);

		pow2_FFT_routine(data, post_pow2);
		for(i=0;(i>> post_pow2) == 0;i++){
			data[i] *= bluewindow[i];
		}

		pow2_IFFT_routine_swapped(data, post_pow2);

		for(i=0;i<mult;i++){
			double ang = i * i * M_PI / mult;//noneg!
			Complex<double> factor(cos(ang),-sin(ang));
			data[i] *= factor;
		}


	}else{
//		comp_FFT_routine(data, pre_pow2, mult);

		unsigned int i;
		for(i=0;i<mult;i++){
			double ang = i * i * M_PI / mult; //noneg!
			Complex<double> factor(cos(ang),-sin(ang));
			data[i] *= factor;
		}

		for(; (i>> post_pow2) == 0;i++) ExOp::toZero(data[i]);

		pow2_FFT_routine(data, post_pow2);
		for(i=0;(i>> post_pow2) == 0;i++){
			data[i] *= bluewindow[i];
		}

		pow2_IFFT_routine_swapped(data, post_pow2);



		for(i=0;i<mult;i++){
			double ang = i * i * M_PI / mult;//noneg!
			Complex<double> factor(cos(ang),-sin(ang));
			data[i] *= factor;
		}
	//	printf("\n");for(i=0;(i>> post_pow2) == 0;i++) ExOp::show(data[i]);
		// A B 0 0 0 0
	}
}

template<class C> void Bluestein::toIFT(C* data) const{
	if (mult == 0){// pow2_FFT!
		pow2_IFFT_routine(data, pre_pow2); //  1 <= mult <= 7
		pow2_bitswap_permutation(data, pre_pow2); //  1 <= mult <= 7
	}else if (pre_pow2 == 0){ // odd window D:

		unsigned int i;
		for(i=0;i<mult;i++){
			double ang = i * i * M_PI / mult;//noneg!
			Complex<double> factor(cos(ang),sin(ang));
			data[i] *= factor;
		}

		for(; (i>> post_pow2) == 0;i++) ExOp::toZero(data[i]);

		pow2_FFT_routine(data, post_pow2);
		for(i=0;(i>> post_pow2)== 0;i++){
			data[i] *= Complex<double>(bluewindow[i][0]/ mult,-bluewindow[i][1] / mult);
		}

		pow2_IFFT_routine_swapped(data, post_pow2);

		for(i=0;i<mult;i++){
			double ang = i * i * M_PI / mult;//noneg!
			Complex<double> factor(cos(ang),sin(ang));
			data[i] *= factor;
		}

	}else{
		comp_IFFT_routine(data, pre_pow2, mult);

		unsigned int i;
		for(i=0;i<mult;i++){
			double ang = i * i * M_PI / mult;//noneg!
			Complex<double> factor(cos(ang),sin(ang));
			data[i] *= factor;
		}

		for(; (i>> post_pow2) == 0;i++) ExOp::toZero(data[i]);

		pow2_FFT_routine(data, post_pow2);
		for(i=0;(i>> post_pow2)== 0;i++){
			data[i] *= Complex<double>(bluewindow[i][0]/ mult,-bluewindow[i][1] / mult);
		}

		pow2_IFFT_routine_swapped(data, post_pow2);

		for(i=0;i<mult;i++){
			double ang = i * i * M_PI / mult;//noneg!
			Complex<double> factor(cos(ang),sin(ang));
			data[i] *= factor;
		}

	}
}


// if N is even:
// F'[k] = \sum_j F[i] * (sin(f i*k) + sin (f * (N-1-i)*k )
// so
// F'[k] = \sum_j F[i] * (2 * sin(f *k (N-1)/2) * cos (f*k* (N-1)/2 - i)

template<class C> void Bluestein::toFT_even(C* data) const{
    toFT(data);

    // need to predend we are combining 2 FFT, but only care about the

    unsigned int i,j;
    for(i= 0; i < tsize/2;i++){
        // f[i] and f[tsize-1-i] needs to be combined
        j = tsize-1-i;
        C tmp = data[i] + data[j];
        data[i] = data[i] - data[j];
        data[j] = tmp;
    }


}

template<class C> void Bluestein::toIFT_even(C* data) const{
    toIFT(data);
}

template<class C, class D> void Bluestein::toFT(D* o_data,const C* i_data) const{
	unsigned int i;
	exit(1);
	if (mult == 0){
		// pow2_FFT!
		for(i=0;i<tsize;i++)  o_data[i] = i_data[i];
		pow2_FFT_routine(o_data, pre_pow2 - 1, 2); //  1 <= mult <= 7
		pow2_bitswap_permutation(o_data, pre_pow2); //  1 <= mult <= 7
	}else{


	}
}

template<class C, class D> void Bluestein::toIFT(D* o_data,const C* i_data) const{
	unsigned int i;
	exit(1);
	if (mult == 0){
		// pow2_FFT!
		for(i=0;i<tsize;i++)  o_data[i] = i_data[i];
		pow2_FFT_routine(o_data, pre_pow2 - 1, 2); //  1 <= mult <= 7
		pow2_bitswap_permutation(o_data, pre_pow2); //  1 <= mult <= 7
	}else{


	}
}


template<class C> void Bluestein::alloc_resample_buffer(C* &iobuf, Bluestein* invblue) const{
}
template<class C> void Bluestein::resample(C* iobuf, Bluestein* invblue) const{
	// apply pow2 transforms
}








	template<int nbstate>
	template <Tuple_flag TF>
	void HMM<nbstate>::init(double transitprob, Tuple<double, nbstate, TF> const & _bound){



		int i,j;
		for(i=0;i<nbstate;i++) {
			boundary[i] = _bound[i];
			for(j=0;j<nbstate;j++){
			transition.data[i+j*nbstate] = _bound[i] * transitprob + ((i == j) ? (1.0f - transitprob) : 0.0f);
		}
		}
	}

	template<int nbstate>
	void HMM<nbstate>::show(FILE* f) const{
		fprintf(f,"boundary:\t");

		int i;
		for(i=0;i<nbstate;i++) fprintf(f,"%e%c", boundary[i], i +1 == nbstate ? '\n': '\t' );
		fprintf(f,"transit:\n");
		transition.show(f);

	}



	template<int nbstate>
	void HMM<nbstate>::swapstates(int a,int b){
		double tmp;
		tmp = boundary[a];
		boundary[a] = boundary[b];
		boundary[b] = tmp;
		int i;
		for(i=0;i< nbstate;i++){
			tmp = transition.data[a + i*nbstate];
			transition.data[a + i*nbstate] = transition.data[b + i*nbstate];
			transition.data[b + i*nbstate] = tmp;
		}
		for(i=0;i< nbstate;i++){
			tmp = transition.data[i + a*nbstate];
			transition.data[i + a*nbstate] = transition.data[i + b*nbstate];
			transition.data[i + b*nbstate] = tmp;
		}
	}

	template<int nbstate>
	template<int lenght>
	void HMM<nbstate>::runHMM(Tuple<Tuple<double, nbstate>, lenght> &likelyhoods, double* weights ){
		runHMM(&(likelyhoods[0]),lenght,weights);
	}

	template<int nbstate>
	void HMM<nbstate>::runHMM(Tuple<double, nbstate>* likelyhoods, int lenght, double* weights){
		Tuple<double, nbstate> buffer[lenght];





	}

	template<int nbstate>
	template< unsigned int nbdim>
	DataGrid< Tuple<double, nbstate>, nbdim> HMM<nbstate>::runHMM(DataGrid< Tuple<double, nbstate>, nbdim> const &like, int dir, bool dolearn){

		DataGrid< Tuple<double, nbstate>, nbdim> _out;

		Tuple<double, nbstate> tmp_out;

		_out.setSizes(like.dims);

		TMatrix<double,nbstate,nbstate> learn;

		int coor[nbdim];
		int z,w;

		TMatrix<double,1,nbstate>* back = new TMatrix<double,1,nbstate>[ like.dims[dir] ];
		TMatrix<double,nbstate,1> tmpb[2];
		int swp=0;
		double tmp;
		int altdir;
		memset(coor,'\0',sizeof(unsigned int)*nbdim);


		TMatrix<double,nbstate,nbstate> tmptransition = transition;

		//	memset(updatetransit,'\0',sizeof(double)*channels*channels); // asume it<s done

		while(true){

			coor[dir] = like.dims[dir];
			for(z=0;z<nbstate;z++)	back[coor[dir]-1].data[z] =1.0f;

			for(coor[dir]--;coor[dir]>0;coor[dir]--) {
				back[coor[dir]-1] = tmptransition * (back[coor[dir]].scale_rows( like(coor) ));
				tmp = back[coor[dir]-1].data[0];
				for(w=1;w<nbstate;w++) tmp += back[coor[dir]-1].data[w];
				if (tmp <= 0.000001f) tmptransition *= 2.0f;
				else if (tmp >= 100000.0f) tmptransition *= 0.5f;


			}


			for(z=0;z<nbstate;z++)	tmpb[swp].data[z] = boundary[z];

			for( ;coor[dir]<like.dims[dir];coor[dir]++) {
				Tuple<double, nbstate> curob = like(coor);
				tmpb[swp^1]= (tmpb[swp] * tmptransition).scale_cols( curob );

				tmp = tmpb[swp^1].data[0];
				for(w=1;w<nbstate;w++) tmp += tmpb[swp^1].data[w];
				if (tmp <= 0.000001f) tmptransition *= 2.0f;
				else if (tmp >= 100000.0f) tmptransition *= 0.5f;

				if ((coor[dir] == 0)||(dolearn == false)){
					for(z=0;z<nbstate;z++) tmp_out[z] = tmpb[swp^1].data[z] * back[coor[dir] ].data[z];
					}else{
						tmp =0;
				for(z=0;z<nbstate;z++) {
					tmp_out[z] = tmpb[swp^1].data[z] * back[coor[dir] ].data[z];


					for(w=0;w<nbstate;w++) {
						learn.data[w + z* nbstate] = tmpb[swp].data[z] * transition.data[w + z * nbstate] * curob[w] * back[coor[dir] ].data[w];
						tmp += learn.data[w + z* nbstate];
					}
				}
						if ((tmp != 0)&&(!isnan(tmp))) {learn *= (1.0f / tmp);

						learned_transition += WeightElem< TMatrix<double,nbstate,nbstate>, 1 >(learn);}
				}
				swp ^=1;
				_out( coor ) =  tmp_out.normalize();
		//		if ((!ExCo<double>::isValid(tmp_out[0]))||(!ExCo<double>::isValid(tmp_out[1]))) printf("%e\t%e\n", tmp_out[0], tmp_out[1]);
		//		if ((!ExCo<double>::isValid(tmp_out.normalize()[0]))||(!ExCo<double>::isValid(tmp_out.normalize()[1]))) printf("%e\t%e\n", tmp_out[0], tmp_out[1]);

			}

			for(altdir=1;altdir<nbdim;altdir++) if (coor[(dir + altdir) % nbdim] >= like.dims[(dir + altdir) % nbdim]-1) coor[(dir + altdir) % nbdim] =0; else {coor[(dir + altdir) % nbdim]++; break;}
			if (altdir == nbdim) break;
		}


		delete[](back);
		return(_out);

	}


	template<int nbstate>
	template<unsigned int nbdim>
	DataGrid< Tuple<double, nbstate>, nbdim> HMM<nbstate>::runHMM_meta(DataGrid< Tuple<double, nbstate>, nbdim> const &like, int dir, double* weights){

		DataGrid< Tuple<double, nbstate>, nbdim> _out;

		DataGrid< Tuple<double, nbstate>, nbdim> backward;

		Tuple<double, nbstate> tmp_out;

		DataGrid< Tuple<double, nbstate>, nbdim-1> forward;
		DataGrid< Tuple<double, nbstate>, nbdim-1> subHMM;

		_out.setSizes(like.dims);
		backward.setSizes(like.dims);

		TMatrix<double,nbstate,nbstate> learn;

		int coor[nbdim];
		int altdim[nbdim-1];
		int altcoor[nbdim-1];
		int i;
		for(i=0;i<nbdim-1;i++){
			altdim[i] = like.dim[(dir +i +1) % nbdim];


		}
		forward.setSizes(altdim);
		subHMM.setSizes(altdim);



		int z,w;

		int swp=0;
		double tmp;
		Tuple<double, nbstate> tmprec;
		int altdir;



		TMatrix<double,nbstate,nbstate> tmptransition = transition;

		//	memset(updatetransit,'\0',sizeof(double)*channels*channels); // asume it<s done
		memset(coor,'\0',sizeof(unsigned int)*nbdim);

		while(true){

			coor[dir] = like.dims[dir];
			for(z=0;z<nbstate;z++)	tmprec[z] =1.0f;

			for(coor[dir]--;coor[dir]>0;coor[dir]--) {
				backward(coor) = tmprec;

				subHMM(altcoor) = tmptransition * (tmprec.scale_rows( like(coor) ));



				tmp = tmprec[0];
				for(w=1;w<nbstate;w++) tmp += tmprec[w];
				if (tmp <= 0.000001f) tmptransition *= 2.0f;
				else if (tmp >= 100000.0f) tmptransition *= 0.5f;
			}
			backward(coor) = tmprec;

			for(altdir=1;altdir<nbdim;altdir++) if (coor[(dir + altdir) % nbdim] >= like.dims[(dir + altdir) % nbdim]-1) coor[(dir + altdir) % nbdim] =0; else {coor[(dir + altdir) % nbdim]++; break;}
			if (altdir == nbdim) break;
		}

		for(i=0;i< like.dims[dir];i++){

		memset(coor,'\0',sizeof(unsigned int)*nbdim);
		memset(altcoor,'\0',sizeof(unsigned int)*nbdim);
		coor[dir] = 0;

		while(true){
		/*
			if (i ==0) for(z=0;z<nbstate;z++) forward(altcoor)[z]
			tmprec = (forward(altcoor) * tmptransition).scale_cols( like(coor) );
			tmp = tmprec[0]
			for(w=1;w<nbstate;w++) tmp += tmprec[w];
			while (tmp <= 0.000001f) {
				tmprec *= 65536.0f;
				tmp *= 65536.0f;
			}
			while (tmp>= 100000.0f) {
				tmprec *= 1.0f /  65536.0f;
				tmp *=  1.0f / 65536.0f;
			}


			forward(altcoor) = tmprec;
			subHMM(altcoor) = tmprec * backward(coor);*/

			for(altdir=1;altdir<nbdim;altdir++) if (coor[(dir + altdir) % nbdim] >= like.dims[(dir + altdir) % nbdim]-1) {coor[(dir + altdir) % nbdim] =0; altcoor[altdir] =0;} else {coor[(dir + altdir) % nbdim]++; altcoor[altdir]++; break;}
			if (altdir == nbdim) break;
		}


		}




		return(_out);

	}




	template<int nbstate>
	void HMM<nbstate>::EMinit(){
		ExOp::toZero(learned_transition);

		TMatrix<double,nbstate,nbstate> dummy;
		int i;
		double a,b;
		ExCo<double>::epsilon(b);
		a = 1.0f - (b*(nbstate-1));
		for(i=0;i<nbstate*nbstate;i++) dummy.data[i] = (0 == (i % (nbstate+1))) ? a: b;
		learned_transition += WeightElem< TMatrix<double,nbstate,nbstate>, 1 >( dummy , 0.01f);
	}
	template<int nbstate>
	void HMM<nbstate>::EMfinit(){

		TMatrix<double,nbstate,nbstate> dummy = learned_transition.getMean();
		if (ExOp::isValid(dummy)){
		Tuple<double,nbstate> sum;
			int z,w;
			for(z=0;z<nbstate;z++) {
				sum[z] = dummy.data[0 + z *  nbstate];
				for(w=1;w<nbstate;w++)  sum[z] += transition.data[w + z *  nbstate];
				if (sum[z] == 0.0f) break;
			}
			if (z == nbstate){
		transition = dummy;
		for(z=0;z<nbstate;z++) {
			sum[z] = transition.data[0 + z *  nbstate];
			for(w=1;w<nbstate;w++)  sum[z] += transition.data[w + z *  nbstate];
			if (sum[z] != 0.0f){
				for(w=0;w<nbstate;w++) transition.data[w + z *  nbstate] /= sum[z];
			}
		}
			boundary = sum.normalize();
			}
		}
	}


	template<int nbstate>
	void HMM<nbstate>::save(FILE* f) const{
		fwrite(&boundary,sizeof(Tuple<double,nbstate> ),1,f);
		fwrite(&transition,sizeof(TMatrix<double,nbstate,nbstate>),1,f);
	}


	template<int nbstate>
	void HMM<nbstate>::load(FILE* f, unsigned int size){
		fread(&boundary,sizeof(Tuple<double,nbstate>),1,f);
		fread(&transition,sizeof(TMatrix<double,nbstate,nbstate>),1,f);
	}


	template<int size>
	Permutation<size>::Permutation(){
		for(int i=0;i<size;i++) map[i] = (char)i;
	}
	template<int size>
	Permutation<size> Permutation<size>::operator-() const{
		Permutation<size> _out;
		int i;
		for(i=0;i<size;i++) _out.map[ (unsigned int)map[i] ] = (char)i;
		return(_out);
	}

	template<int size>
	Permutation<size> Permutation<size>::operator+(const Permutation<size> &other) const{
		Permutation<size> _out;
		int i;
		for(i=0;i<size;i++) _out.map[ i ] =  map[other[i]];
	}

	template<int size>
	int Permutation<size>::operator[](const int &offset) const{
		return((int)map[offset]);
	}

	template<int size>
	template<class C, Tuple_flag Cflag, unsigned int esize>
	Tuple<C, size + esize, Cflag> Permutation<size>::operator+(const Tuple<C, size +  esize, Cflag>& input) const{
		Tuple<C, esize + size, Cflag> _out;

		int i;
			for(i=size + esize-1;i>=size;i--) _out[i] = input[i];
			for(;i>=0;i--) _out[i] = input[map[i]];

		return(_out);
	}



#undef LFHTEMP
#define LFHTEMP template<class C>

LFHTEMP template<class O, class B> void Matrix<C>::alloc_add(Matrix<O>& fout, const Matrix<B>& other)const{
	bool mismatch = false;
	if (other.sizex > sizex){
		mismatch |= true;
		fout.sizex = other.sizex;
	}else {
		mismatch |= (other.sizex<sizex);
		fout.sizex = sizex;
	}

	if (other.sizey > sizey){
		static_warning_handdle << LFH_WARNING_MATRIX_SIZES_MISMATCH;
		fout.sizey = other.sizey;
		}else {
			mismatch |= (other.sizey<sizey);
			fout.sizey = sizey;
		}
		if (mismatch) static_warning_handdle << LFH_WARNING_MATRIX_SIZES_MISMATCH;
		fout.data = new O[fout.sizex*fout.sizey];
	}

LFHTEMP Matrix<C> Matrix<C>::operator+(const Matrix<C>& other)const{Matrix<C> fout; alloc_add(fout,other);
	unsigned int i = sizey*sizex; for(i--;i!=0xFFFFFFFF;i--) fout[i] = data[i] + other.data[i];
	return (fout);
}
LFHTEMP Matrix<C> Matrix<C>::operator-(const Matrix<C>& other)const{Matrix<C> fout; alloc_add(fout,other);
	unsigned int i = sizey*sizex; for(i--;i!=0xFFFFFFFF;i--) fout[i] = data[i] - other.data[i];
	return (fout);
}
LFHTEMP Matrix<C>& Matrix<C>::operator+=(const Matrix<C>& other){
	unsigned int i = sizey*sizex; for(i--;i!=0xFFFFFFFF;i--) data[i] += other.data[i];
	return(*this);
}
LFHTEMP Matrix<C>& Matrix<C>::operator-=(const Matrix<C>& other){
	unsigned int i = sizey*sizex; for(i--;i!=0xFFFFFFFF;i--) data[i] -= other.data[i];
	return(*this);
}

LFHTEMP Matrix<C> Matrix<C>::operator-()const{Matrix<C> fout(sizex,sizey);
	unsigned int i = sizey*sizex; for(i--;i!=0xFFFFFFFF;i--) fout.data[i] = -data[i];
	return (fout);
}

LFHTEMP Matrix<C> Matrix<C>::operator*(const Matrix<C>& other)const{ Matrix<C> fout(other.sizex, sizey);
    C* t = fout.data;
	unsigned int i,j,k,l;
		for(j=0;j<fout.sizey;j++)
            for(i=0;i<fout.sizex;i++,t++){
                (*t) = data[j * sizex] * other.data[i];
                for(k=1;k< sizex;k++) (*t) += data[k + j * sizex] * other.data[i +k * other.sizex];
			}
    return(fout);
}

LFHTEMP void Matrix<C>::InitFromTrianglix(const C* odata){
    unsigned int i,j,k;
        for(j=0,k=0;j<sizex;j++) {
        for(i=0;i<j;i++) {data[i + j * sizex] = odata[k]; data[j + i * sizex] = ExOp::mktrju(odata[k++]);}
        data[j * (sizex+1)] = odata[k++];
        }
    }

LFHTEMP void Matrix<C>::show(FILE *f, int level)const{
    const C* t = data;
    unsigned int i,j;
    if (level == 0) fprintf(f, "Matrix (size = %i,%i):\n",sizex, sizey);
	for(j=0;j<sizey;j++)
        for(i=0;i<sizex;i++,t++){ExOp::show(*t,f,level+2); fprintf(f, "%c", (i+1 == sizex) ? '\n' : '\t');}
    
}

LFHTEMP Matrix<C> Matrix<C>::inverse() const{
}

LFHTEMP double Matrix<C>::pnorm() const{
    C* k = data ;
    double fout =  ExOp::pnorm(*(k++));
    while((k - data) < sizex*sizey) fout += ExOp::pnorm(*(k++));
    return fout;
}
	
template<class C, class V>  GaussScope<C,V>& GaussScope<C,V>::toZero(){w = 0.0f;w2 = 0.0f;ExOp::toZero(mean);ExOp::toZero(var);return(*this);}

// class Trianglix
#undef LFHTEMP
#define LFHTEMP template<class C, unsigned int SIZE>
LFHTEMP Trianglix<C, SIZE> GaussElem< Tuple<C, SIZE>, 0 >::getCovariance()const{
    Trianglix<C, SIZE> tcov = this->getCovariance_biased();
	/*
	C maxEigen = tcov.maxEigenValue(); // scales with
//	printf("Maximum Eigen val: "); ExOp::show(maxEigen);
	maxEigen *= 1.0f / ((w*w/ w2) - 1.0f);
//	printf("Minimum Variance: "); ExOp::show(maxEigen);
	*/
	unsigned int i,j;
	double unbiased_factor =  1.0f / (1.0f - (w2/ (w*w)));
	for(i=0,j=0;j<SIZE;i+=j+1) {tcov.data[i] *= unbiased_factor; j++;}
    return tcov;
    }
	
LFHTEMP Trianglix<C, SIZE> GaussElem< Tuple<C, SIZE>, 0 >::getCovariance_biased()const{Trianglix<C, SIZE> tcov = cov * (1.0f/  w);tcov -= Trianglix<C, SIZE>(mean / w);return tcov;}

LFHTEMP void GaussElem< Tuple<C, SIZE>, 0 >::setCovariance(const Trianglix<C, SIZE>& covar){
    unsigned int i,j,k;
    w2 = 0;
    C tmp;
    for(j=0,k=0;j<SIZE;j++) {
        tmp = ExOp::mktrju(mean[j]) * w;
        for(i=0;i<=j;i++,k++) cov.data[k] = covar.data[k] + tmp * mean[i];
    }
}
	LFHTEMP double GaussElem< Tuple<C, SIZE>, 0 >::likelihoodratio_dist(const GaussElem< Tuple<C, SIZE>, 0 > &other, bool has_weight, bool has_covariance) const{
        
        Tuple<C, SIZE> mdiff = (mean/ w) - (other.mean / other.w) ;
		double fout= 0.0f;
		if (has_covariance){
		
        C t = determinant; ExOp::toZero(t);
		
        if (determinant == t) determinant = getCovariance().log_determinant();
        if (other.determinant == t) other.determinant = other.getCovariance().log_determinant();
		if (has_weight){
				GaussElem< Tuple<C, SIZE>, 0 > temporary = (*this) + other;
				temporary.determinant =  temporary.getCovariance().log_determinant();
				fout = 0.5f * ((w + other.w) * temporary.determinant - w *  determinant - other.w * other.determinant );
			if (fout < 0.0f) {printf("got negative distance! %e\n",fout);this->show();other.show(); exit(1);}
		 return fout;
			}else{
				GaussElem< Tuple<C, SIZE>, 0 > temporary;
				temporary.w = 1.0f;
				temporary.w2 = ((*this).w2 + other.w2)/((*this).w + other.w) ;
				for(unsigned int i=0;i<SIZE;i++) temporary.mean[i] = (mean[i] + other.mean[i])* 0.5f;
				for(unsigned int i=0;i<TEMPLATE_TRIANGLE_NUMBER<SIZE,2>::ans;i++) temporary.cov.data[i] = (cov.data[i] + other.cov.data[i])* 0.5f;
				
				temporary.determinant =  temporary.getCovariance().log_determinant();
				
				return( 0.5f * ( temporary.determinant - 0.5f * ( determinant + other.determinant) ) );
			}
      //  return (tcov.weighted_bhattacharryya_partial(mdiff,ocov, w, other.w) - 0.25f * (w * determinant+ other.w *other.determinant));
        }else{
            Tuple<C, SIZE> tcov = getVar();
            Tuple<C, SIZE> ocov = other.getVar();
            unsigned int i;
			
            if (has_weight) {
                for(i=0;i<SIZE;i++){
                    fout += log(1.0f + mdiff[i] * mdiff[i] * w * other.w / ((w +  other.w)*(w * tcov[i] + other.w * ocov[i]))) + log(w * tcov[i] + other.w *  ocov[i]) - (w * log(tcov[i]) + other.w * log(ocov[i])) / (w + other.w);
                }
				fout -= log(w + other.w) * SIZE;
				fout *= (w + other.w);
            }else{
                for(i=0;i<SIZE;i++){
                    fout += log(1.0f + 0.5f * mdiff[i] * mdiff[i] / (tcov[i] + ocov[i])) + log( tcov[i] +  ocov[i])  +0.5f * ((mdiff[i] * mdiff[i]  / (tcov[i] + ocov[i]))- log(tcov[i]) - log(ocov[i] )  );
                }
				fout -= log(2) * SIZE ; 
            }
            return(0.5f *fout);
        }
	}
	

	LFHTEMP double GaussElem< Tuple<C, SIZE>, 0 >::bhattacharryya_dist(const GaussElem< Tuple<C, SIZE>, 0 > &other, bool has_weight, bool has_covariance) const{
        
        Tuple<C, SIZE> mdiff = (mean/ w) - (other.mean / other.w) ;
		
		if (has_covariance){
        Trianglix<C, SIZE> tcov = getCovariance();
        Trianglix<C, SIZE> ocov = other.getCovariance();
		
        C t = determinant; ExOp::toZero(t);
		
        if (determinant == t) determinant = tcov.log_determinant();
        if (other.determinant == t) other.determinant = ocov.log_determinant();
		
		if (has_weight) return (0.125f * (tcov.weighted_bhattacharryya_partial(mdiff,ocov, w, other.w) - 0.5f * (w * determinant+ other.w *other.determinant)));
		else return tcov.bhattacharryya_partial(mdiff,ocov) - 0.25f * (determinant+ other.determinant);
			
        

        }else{
            Tuple<C, SIZE> tcov = getVar();
            Tuple<C, SIZE> ocov = other.getVar();
            unsigned int i;
            double fout = 0.0f;
        if (has_weight) {
            for(i=0;i<SIZE;i++){
                fout += 0.25f * (w + other.w) * mdiff[i] * mdiff[i] / (w * tcov[i] + other.w * ocov[i]) + log(w * tcov[i] + other.w *  ocov[i]) - log(w + other.w) - ((w * log(tcov[i]) + other.w * log(ocov[i])) / (w + other.w));
            }
			return(0.25f * fout * (w + other.w));
        }else{
            for(i=0;i<SIZE;i++){
                fout += log( tcov[i] +  ocov[i]) - log(2) +0.5f * ((mdiff[i] * mdiff[i]  / (tcov[i] + ocov[i]))- log(tcov[i]) - log(ocov[i] )  );
            }
			return(0.5f * fout);
        }
           
        }
	}
	
	LFHTEMP double GaussElem< Tuple<C, SIZE>, 0 >::Entropy(bool has_covariance) const{
        if (has_covariance) return( 0.5f * (getCovariance_biased().log_determinant() + (1.0f + log(M_PI * 2.0f)) * SIZE));
        else{
            Tuple<C, SIZE> tcov = getVar();
            double fout = (1.0f + log(M_PI * 2.0f)) * SIZE;
            for(unsigned int i=0;i<SIZE;i++) fout += log(tcov[i]);
            return (fout * 0.5f);
        }
    }
	
LFHTEMP double GaussElem< Tuple<C, SIZE>, 0 >::LLikelihood(const Tuple<C, SIZE> &what, bool has_covariance)const{
    double fout;
    Tuple<C, SIZE> dev = what - (mean/ w);
    if (has_covariance){
        C t = determinant; ExOp::toZero(t);
        Trianglix<C, SIZE> cov = getCovariance_biased();
        if (determinant == t){
            determinant = cov.log_determinant();
        }
        fout = (((double)determinant) + log(M_PI * 2.0f) * SIZE);
        fout += cov.Xformed_inner_product_of_inverse(dev);
        fout *= -0.5f;
    }else{
    }
    return fout;
}

LFHTEMP GaussElem< Tuple<C, SIZE>, 0u >& GaussElem< Tuple<C, SIZE>, 0 >::toHouseHolderMultiply(const Tuple<C, SIZE> &vec, const C &denum2){
	
	return(*this);
	}
	

LFHTEMP template<class D, unsigned int O_SIZE> GaussElem<Tuple<C,O_SIZE>,0u> GaussElem<Tuple<C, SIZE>,0u>::operator*(const TMatrix<D,SIZE , O_SIZE> & projection) const{
	GaussElem< Tuple<C, O_SIZE>, 0u> fout;
	fout.w = w;
	fout.w2 = w2;
	fout.mean = projection * mean;
	fout.cov = projection * cov;
	

	return(fout);
	}

LFHTEMP template<unsigned int S_SIZE> GaussElem< Tuple<C, S_SIZE>, 0u > GaussElem< Tuple<C, SIZE>, 0 >::makeSubGauss(const Tuple<unsigned int, S_SIZE>& s_d)const{ GaussElem< Tuple<C, S_SIZE>, 0u > fout;
		fout.w = w;
		fout.w2 = w2;
		unsigned int i;
		for(i=0;i<S_SIZE;i++) fout.mean[i] = mean[s_d[i]];
		unsigned int j,k;
		for(i=0,k=0;i<S_SIZE;i++)
			for(j=0;j<=i;j++,k++) fout.cov.data[k] = cov.cell(s_d[i],s_d[j]);
		return(fout);
	}
	
LFHTEMP template<unsigned int S_SIZE> GaussElem< Tuple<C, S_SIZE>, 0u > GaussElem< Tuple<C, SIZE>, 0 >::makeResidualGauss(const Tuple<unsigned int, S_SIZE>& s_d)const{ GaussElem< Tuple<C, S_SIZE>, 0u > fout;
			fout.w = w;
			fout.w2 = w2;
			unsigned int i,j,k;
			for(i=0;i<S_SIZE;i++) fout.mean[i] = mean[s_d[i]];
			
			Tuple<unsigned int, SIZE> o_d;
	
			// find left over dimentions
			for(j=0,k=0;k<SIZE;k++){
				for(i=0;i < S_SIZE;i++) if (s_d[i] == k) break;
				o_d[ (i == S_SIZE) ? S_SIZE + j++ : k - j ] = k;
			}

			Trianglix<C, SIZE> tmp_cov;
	
			for(i=0,k=0;i<SIZE;i++)
				for(j=0;j<=i;j++,k++) tmp_cov.data[k] = cov.data[k] - mean[i] * mean[j] / w;
	
			// partial Cholesky decomposition!
	
	//tmp_cov.show();
			C ivar;
			for(k=SIZE-1;k>=S_SIZE;k--){
				ivar = ExOp::mkinverse(tmp_cov.cell(o_d[k],o_d[k]));
				if (!ExOp::isValid(ivar)) continue;
				//ExOp::show(ivar);
				for(i=0;i<k;i++)
					for(j=0;j<=i;j++) tmp_cov.cell(o_d[i],o_d[j]) -= ivar * tmp_cov.cell(o_d[i],o_d[k]) * tmp_cov.cell(o_d[j],o_d[k]);
			}
	//tmp_cov.show();
	
			for(i=0,k=0;i<S_SIZE;i++)
				for(j=0;j<=i;j++,k++) fout.cov.data[k] = tmp_cov.cell(s_d[i],s_d[j]) +  (fout.mean[i] * fout.mean[j] / w);
				
				return(fout);
			}
	
	LFHTEMP void GaussElem< Tuple<C, SIZE>, 0 >::setMeanVaronDimention(double mean, double var, unsigned int dim){
		
	}
	#ifdef GNU_SCIENTIFIC_LIBRARY 
	LFHTEMP  double GaussElem< Tuple<C, SIZE>, 0 >::PvalueOf(const GaussElem< Tuple<C, SIZE>, 0u > &x)const{
		unsigned int i,j,k;
		
		// fill permuted trianglix this, unbiased variance
		Trianglix<C, SIZE> mod_this;
		Trianglix<C, SIZE> tmp_cov;
		for(i=0,k=0;i<SIZE;i++)
			for(j=0;j<=i;j++,k++) {
				mod_this.data[k] = (this->cov.data[k] * this->w  - this->mean[i] * this->mean[j]) / (this->w * this->w - this->w2);
				tmp_cov.data[k] = (x.cov.data[k] * this->w - x.mean[i] * x.mean[j]) / (x.w * x.w - x.w2);
			}
		
		/*
		// diagonalize this, multiply that
		C ivar;
		C cho_vec[SIZE];
		C buffer[SIZE];
		for(k=SIZE-1;k>0;k--){
			ivar = ExOp::mkinverse(mod_this.cell(k,k));
			if (!ExOp::isValid(ivar)) continue;
			// double mutiply by 
			
			cho_vec[k] = sqrt(ivar);
			for(i=0;i<k;i++) cho_vec[i] = cho_vec[k] * tmp_cov.cell(i,k);
			mod_this.CholeskyStep_up(cho_vec,k,buffer);
			
			
			//ExOp::show(ivar);
			for(i=0;i<k;i++)
				for(j=0;j<=i;j++) {
					mod_this.cell(i,j) -= ivar * mod_this.cell(i,k) * mod_this.cell(j,k);
				}
		}*/
		
		double l_pval = 0.0f;
		
		for(i=0;i<SIZE;i++) {k = (i * (i+3)) / 2; l_pval += log(Pvalue_Gamma_Ptail(tmp_cov.data[k] , 0.5f * this->w,  mod_this.data[k] * 2.0f));}
		
		return Pvalue_SumLogPvalues_Ptail(l_pval, SIZE);
	}

	LFHTEMP template<unsigned int S_SIZE> double GaussElem< Tuple<C, SIZE>, 0 >::PvalueOfResidual(const GaussElem< Tuple<C, SIZE>, 0u > &x , const Tuple<unsigned int, S_SIZE>& s_d) const{
		unsigned int i,j,k;
		
		Tuple<unsigned int, SIZE> o_d;
		
		// find left over dimentions
		for(j=0,k=0;k<SIZE;k++){
			for(i=0;i < S_SIZE;i++) if (s_d[i] == k) break;
			o_d[ (i == S_SIZE) ? S_SIZE + j++ : k - j ] = k;
		}
		
		// fill permuted trianglix this, biased variance
		Trianglix<C, SIZE> tmp_cov;
		for(i=0,k=0;i<SIZE;i++)
			for(j=0;j<=i;j++) tmp_cov.data[k++] = (x.cov.cell(o_d[i],o_d[j]) * this->w - x.mean[o_d[i]] * x.mean[o_d[j]]) / (x.w * x.w - x.w2);

		// fill permuted trianglix this, unbiased variance
		Trianglix<C, SIZE> mod_this;
		for(i=0,k=0;i<SIZE;i++)
			for(j=0;j<=i;j++) mod_this.data[k++] = (this->cov.cell(o_d[i],o_d[j]) * this->w  - this->mean[o_d[i]] * this->mean[o_d[j]]) / (this->w * this->w - this->w2);
		
		C ivar;
		C cho_vec[SIZE];
		C buffer[SIZE];
		
		for(k=SIZE-1;k>=S_SIZE;k--){
			ivar = ExOp::mkinverse(tmp_cov.cell(k,k));
			if (!ExOp::isValid(ivar)) continue;
			// double mutiply by 
			
			cho_vec[k] = sqrt(ivar);
			for(i=0;i<k;i++) cho_vec[i] = cho_vec[k] * tmp_cov.cell(i,k);
			mod_this.CholeskyStep_up(cho_vec,k,buffer);
			
			
			//ExOp::show(ivar);
			for(i=0;i<k;i++)
				for(j=0;j<=i;j++) {
					tmp_cov.cell(i,j) -= ivar * tmp_cov.cell(i,k) * tmp_cov.cell(j,k);
					}
		}
		
		return 0.0f;
	
	}
	
	LFHTEMP double GaussElem< Tuple<C, SIZE>, 0 >::Pvalue_Hotelling(const GaussElem< Tuple<C, SIZE>, 0u > &other, bool LogPvalue) const{
		// Assumes that all weights sumed were 1 (or equal), so that n = w*w/w2
		
		double n[2];
		n[0] = w*w/w2;
		n[1] = other.w*other.w/other.w2;
		if ((mean.getSize() == 0)||(w == 0)||(other.w == 0)) (LogPvalue) ? 0.0f : 1.0f; // invalid, or no data for one of the two; 
		
		if (n[0] + n[1] <= 1.0f + SIZE){
			// invalid, too little samples to get any pvalue
			return (LogPvalue) ? 0.0f : 1.0f;
		}else{
			Trianglix<C,SIZE> tcov = (cov - Trianglix<C,SIZE>(mean * pow(w,-0.5f) )) * (n[0] / w2) + (other.cov - Trianglix<C,SIZE>(other.mean * pow(other.w,-0.5f))) * (n[1] / other.w2);
			Tuple<C,SIZE> dif = (mean / w) - (other.mean / other.w);
			C inner = tcov.Xformed_inner_product_of_inverse(dif);
			
			double t = inner * n[0] * n[1] * (n[0] + n[1] - SIZE -1) / (SIZE * (n[0] + n[1])) ;
			
			n[1] = n[0] + n[1] - SIZE - 1;
			
			if (t < 0.0f) {printf("Pvalue_Hotelling has failed (singular matrix?) ...\n");return (LogPvalue) ? 0.0f : 1.0f;}
			return (LogPvalue) ? LogPvalue_Fdistrib_Ptail(t,SIZE, n[1]) : Pvalue_Fdistrib_Ptail(t,SIZE, n[1] );
		}
	}

#endif

LFHTEMP GaussElem< Tuple<C, SIZE>, 0 >& GaussElem< Tuple<C, SIZE>, 0 >::addGaussElem_free_mean(const GaussElem< Tuple<C, SIZE>, TARG2 > &other){
	cov -= Trianglix<C, SIZE>( mean * pow(w, -0.5f));
	w+=other.w; w2+=other.w2; mean += other.mean; cov += other.getCovariance_biased() * other.w;
	cov += Trianglix<C, SIZE>( mean * pow(w, -0.5f));	
	ExOp::toZero(determinant); return(*this);
}


	LFHTEMP void GaussElem< Tuple<C, SIZE>, 0 >::load(FILE *f, unsigned int size){
		fread(&w,sizeof(double),1,f);
		fread(&w2,sizeof(double),1,f);
		ExOp::load(determinant,f);
		ExOp::load(mean,f);
		ExOp::load(cov,f);
	}
	LFHTEMP void GaussElem< Tuple<C, SIZE>, 0 >::save(FILE *f) const{
		fwrite(&w,sizeof(double),1,f);
		fwrite(&w2,sizeof(double),1,f);
		C t = determinant; ExOp::toZero(t);
		if (determinant == t) {
			Trianglix<C, SIZE> tcov = getCovariance();
			determinant = tcov.log_determinant();
		}
		ExOp::save(determinant,f);
		ExOp::save(mean,f);
		ExOp::save(cov,f);
	}
	
#undef LFHTEMP
#define LFHTEMP template<class C>
	LFHTEMP Trianglix<C,0u> GaussElem< Tuple<C,0u>, 0 >::getCovariance()const{
		Trianglix<C,0u> tcov = this->getCovariance_biased();
		unsigned int i,j;
		double unbiased_factor =  1.0f / (1.0f - (w2/ (w*w)));
		for(i=0,j=0;j<mean.getSize();i+=j+1) {tcov.data[i] *= unbiased_factor; j++;}
		return tcov;
    }
	
	LFHTEMP Trianglix<C,0u> GaussElem< Tuple<C,0u>, 0 >::getCovariance_biased()const{Trianglix<C,0u> tcov = cov * (1.0f/  w);tcov -= Trianglix<C,0u>(mean / w);return tcov;}
	
	LFHTEMP void GaussElem< Tuple<C,0u>, 0 >::setCovariance(const Trianglix<C,0u>& covar){
		unsigned int i,j,k;
		w2 = 0;
		C tmp;
		for(j=0,k=0;j<mean.getSize();j++) {
			tmp = ExOp::mktrju(mean[j]) * w;
			for(i=0;i<=j;i++,k++) cov.data[k] = covar.data[k] + tmp * mean[i];
		}
	}
	LFHTEMP double GaussElem< Tuple<C,0u>, 0 >::likelihoodratio_dist(const GaussElem< Tuple<C,0u>, 0 > &other, bool has_weight, bool has_covariance) const{
        
        Tuple<C,0u> mdiff = (mean/ w) - (other.mean / other.w) ;
		double fout= 0.0f;
		if (has_covariance){
			
			C t = determinant; ExOp::toZero(t);
			
			if (determinant == t) determinant = getCovariance().log_determinant();
			if (other.determinant == t) other.determinant = other.getCovariance().log_determinant();
			if (has_weight){
				GaussElem< Tuple<C,0u>, 0 > temporary = (*this) + other;
				temporary.determinant =  temporary.getCovariance().log_determinant();
				fout = 0.5f * ((w + other.w) * temporary.determinant - w *  determinant - other.w * other.determinant );
				if (fout < 0.0f) {printf("got negative distance! %e\n",fout);this->show();other.show(); exit(1);}
				return fout;
			}else{
				GaussElem< Tuple<C,0u>, 0 > temporary;
				temporary.w = 1.0f;
				temporary.w2 = ((*this).w2 + other.w2)/((*this).w + other.w) ;
				for(unsigned int i=0;i<mean.getSize();i++) temporary.mean[i] = (mean[i] + other.mean[i])* 0.5f;
				for(unsigned int i=0;i< (1+mean.getSize())*(mean.getSize()) /2;i++) temporary.cov.data[i] = (cov.data[i] + other.cov.data[i])* 0.5f;
				
				temporary.determinant =  temporary.getCovariance().log_determinant();
				
				return( 0.5f * ( temporary.determinant - 0.5f * ( determinant + other.determinant) ) );
			}
        }else{
            Tuple<C,0u> tcov = getVar();
            Tuple<C,0u> ocov = other.getVar();
            unsigned int i;
			
            if (has_weight) {
                for(i=0;i<mean.getSize();i++){
                    fout += log(1.0f + mdiff[i] * mdiff[i] * w * other.w / ((w +  other.w)*(w * tcov[i] + other.w * ocov[i]))) + log(w * tcov[i] + other.w *  ocov[i]) - (w * log(tcov[i]) + other.w * log(ocov[i])) / (w + other.w);
                }
				fout -= log(w + other.w) * mean.getSize();
				fout *= (w + other.w);
            }else{
                for(i=0;i<mean.getSize();i++){
                    fout += log(1.0f + 0.5f * mdiff[i] * mdiff[i] / (tcov[i] + ocov[i])) + log( tcov[i] +  ocov[i])  +0.5f * ((mdiff[i] * mdiff[i]  / (tcov[i] + ocov[i]))- log(tcov[i]) - log(ocov[i] )  );
                }
				fout -= log(2) * mean.getSize() ; 
            }
            return(0.5f *fout);
        }
	}
	
	
	LFHTEMP double GaussElem< Tuple<C,0u>, 0 >::bhattacharryya_dist(const GaussElem< Tuple<C,0u>, 0 > &other, bool has_weight, bool has_covariance) const{
        
        Tuple<C,0u> mdiff = (mean/ w) - (other.mean / other.w) ;
		
		if (has_covariance){
			Trianglix<C,0u> tcov = getCovariance();
			Trianglix<C,0u> ocov = other.getCovariance();
			
			C t = determinant; ExOp::toZero(t);
			
			if (determinant == t) determinant = tcov.log_determinant();
			if (other.determinant == t) other.determinant = ocov.log_determinant();
			
			//if (has_weight) return (0.125f * (tcov.weighted_bhattacharryya_partial(mdiff,ocov, w, other.w) - 0.5f * (w * determinant+ other.w *other.determinant)));
			//else
			return tcov.bhattacharryya_partial(mdiff,ocov) - 0.25f * (determinant+ other.determinant);
			
			
			
        }else{
            Tuple<C,0u> tcov = getVar();
            Tuple<C,0u> ocov = other.getVar();
            unsigned int i;
            double fout = 0.0f;
			if (has_weight) {
				for(i=0;i<mean.getSize();i++){
					fout += 0.25f * (w + other.w) * mdiff[i] * mdiff[i] / (w * tcov[i] + other.w * ocov[i]) + log(w * tcov[i] + other.w *  ocov[i]) - log(w + other.w) - ((w * log(tcov[i]) + other.w * log(ocov[i])) / (w + other.w));
				}
				return(0.25f * fout * (w + other.w));
			}else{
				for(i=0;i<mean.getSize();i++){
					fout += log( tcov[i] +  ocov[i]) - log(2) +0.5f * ((mdiff[i] * mdiff[i]  / (tcov[i] + ocov[i]))- log(tcov[i]) - log(ocov[i] )  );
				}
				return(0.5f * fout);
			}
			
        }
	}
	
	LFHTEMP double GaussElem< Tuple<C,0u>, 0 >::Entropy(bool has_covariance) const{
        if (has_covariance){

			return( 0.5f * (getCovariance_biased().log_determinant() + (1.0f + log(M_PI * 2.0f)) * mean.getSize()));
        }else{
            Tuple<C,0u> tcov = getVar();
            double fout = (1.0f + log(M_PI * 2.0f)) * mean.getSize();
            for(unsigned int i=0;i<mean.getSize();i++) fout += log(tcov[i]);
            return (fout * 0.5f);
        }
    }
	
	LFHTEMP GaussElem< Tuple<C,0u>, 0u >& GaussElem< Tuple<C,0u>, 0 >::toHouseHolderMultiply(const Tuple<C,0u> &vec, const C &denum2){
		
		return(*this);
	}
	
	LFHTEMP template<unsigned int S_SIZE> GaussElem< Tuple<C, S_SIZE>, 0u > GaussElem< Tuple<C,0u>, 0 >::makeSubGauss(const Tuple<unsigned int, S_SIZE>& s_d)const{ GaussElem< Tuple<C, S_SIZE>, 0u > fout;
		fout.w = w;
		fout.w2 = w2;
		unsigned int i;
		for(i=0;i<S_SIZE;i++) fout.mean[i] = mean[s_d[i]];
		unsigned int j,k;
		for(i=0,k=0;i<S_SIZE;i++)
			for(j=0;j<=i;j++,k++) fout.cov.data[k] = cov.cell(s_d[i],s_d[j]);
		return(fout);
	}
	
	LFHTEMP template<unsigned int S_SIZE> GaussElem< Tuple<C, S_SIZE>, 0u > GaussElem< Tuple<C,0u>, 0 >::makeResidualGauss(const Tuple<unsigned int, S_SIZE>& s_d)const{ GaussElem< Tuple<C, S_SIZE>, 0u > fout;
		fout.w = w;
		fout.w2 = w2;
		unsigned int i,j,k;
		for(i=0;i<S_SIZE;i++) fout.mean[i] = mean[s_d[i]];
		
		Tuple<unsigned int,0u> o_d;
		
		// find left over dimentions
		for(j=0,k=0;k<mean.getSize();k++){
			for(i=0;i < S_SIZE;i++) if (s_d[i] == k) break;
			o_d[ (i == S_SIZE) ? S_SIZE + j++ : k - j ] = k;
		}
		
		Trianglix<C,0u> tmp_cov;
		
		for(i=0,k=0;i<mean.getSize();i++)
			for(j=0;j<=i;j++,k++) tmp_cov.data[k] = cov.data[k] - mean[i] * mean[j] / w;
		
		// partial Cholesky decomposition!
		
		//tmp_cov.show();
		C ivar;
		for(k=mean.getSize()-1;k>=S_SIZE;k--){
			ivar = ExOp::mkinverse(tmp_cov.cell(o_d[k],o_d[k]));
			if (!ExOp::isValid(ivar)) continue;
			//ExOp::show(ivar);
			for(i=0;i<k;i++)
				for(j=0;j<=i;j++) tmp_cov.cell(o_d[i],o_d[j]) -= ivar * tmp_cov.cell(o_d[i],o_d[k]) * tmp_cov.cell(o_d[j],o_d[k]);
		}
		//tmp_cov.show();
		
		for(i=0,k=0;i<S_SIZE;i++)
			for(j=0;j<=i;j++,k++) fout.cov.data[k] = tmp_cov.cell(s_d[i],s_d[j]) +  (fout.mean[i] * fout.mean[j] / w);
		
		return(fout);
	}
	
	LFHTEMP GaussElem< Tuple<C, 0u> > GaussElem< Tuple<C,0u>, 0 >::mkSubGaussElem(const Tuple<bool> filter)const{ GaussElem< Tuple<C, 0u>, 0u > fout;
		unsigned int i,j,k,l;
		fout.w = w;
		fout.w2 = w2;
		for(j=0,i=0;i< filter.getSize();i++) if (filter[i]) j++;
		fout.setSize(j);
		for(j=0,i=0;i< filter.getSize();i++) if (filter[i]) {
			fout.mean[j++] = mean[i];
			for(l=0,k=0;k <= i;k++) if (filter[k]) {
				fout.cov.cell(j,l) = fout.cov.cell(i,k); l++;
			}
		}
		return(fout);
	}

	LFHTEMP  double GaussElem< Tuple<C,0u>, 0 >::PvalueOf(const GaussElem< Tuple<C,0u>, 0u > &x)const{
		unsigned int i,j,k;
		
		// fill permuted trianglix this, unbiased variance
		Trianglix<C,0u> mod_this;
		Trianglix<C,0u> tmp_cov;
		for(i=0,k=0;i<mean.getSize();i++)
			for(j=0;j<=i;j++,k++) {
				mod_this.data[k] = (this->cov.data[k] * this->w  - this->mean[i] * this->mean[j]) / (this->w * this->w - this->w2);
				tmp_cov.data[k] = (x.cov.data[k] * this->w - x.mean[i] * x.mean[j]) / (x.w * x.w - x.w2);
			}
		
	
		double l_pval = 0.0f;
		
		for(i=0;i<mean.getSize();i++) {k = (i * (i+3)) / 2; + l_pval += log(Pvalue_Gamma_Ptail(tmp_cov.data[k] , 0.5f * this->w,  2.0f * mod_this.data[k] ));}
		
		return Pvalue_SumLogPvalues_Ptail(l_pval, mean.getSize());
	}
	
	LFHTEMP template<unsigned int S_SIZE> double GaussElem< Tuple<C,0u>, 0 >::PvalueOfResidual(const GaussElem< Tuple<C,0u>, 0u > &x , const Tuple<unsigned int, S_SIZE>& s_d) const{
		unsigned int i,j,k;
		
		Tuple<unsigned int,0u> o_d;
		
		// find left over dimentions
		for(j=0,k=0;k<mean.getSize();k++){
			for(i=0;i < S_SIZE;i++) if (s_d[i] == k) break;
			o_d[ (i == S_SIZE) ? S_SIZE + j++ : k - j ] = k;
		}
		
		// fill permuted trianglix this, biased variance
		Trianglix<C,0u> tmp_cov;
		for(i=0,k=0;i<mean.getSize();i++)
			for(j=0;j<=i;j++) tmp_cov.data[k++] = (x.cov.cell(o_d[i],o_d[j]) * this->w - x.mean[o_d[i]] * x.mean[o_d[j]]) / (x.w * x.w - x.w2);
		
		// fill permuted trianglix this, unbiased variance
		Trianglix<C,0u> mod_this;
		for(i=0,k=0;i<mean.getSize();i++)
			for(j=0;j<=i;j++) mod_this.data[k++] = (this->cov.cell(o_d[i],o_d[j]) * this->w  - this->mean[o_d[i]] * this->mean[o_d[j]]) / (this->w * this->w - this->w2);
		
		C ivar;
		C* cho_vec =new C[mean.getSize()];
		C* buffer =new C[mean.getSize()];
		
		for(k=mean.getSize()-1;k>=S_SIZE;k--){
			ivar = ExOp::mkinverse(tmp_cov.cell(k,k));
			if (!ExOp::isValid(ivar)) continue;
			// double mutiply by 
			
			cho_vec[k] = sqrt(ivar);
			for(i=0;i<k;i++) cho_vec[i] = cho_vec[k] * tmp_cov.cell(i,k);
			mod_this.CholeskyStep_up(cho_vec,k,buffer);
			
			
			//ExOp::show(ivar);
			for(i=0;i<k;i++)
				for(j=0;j<=i;j++) {
					tmp_cov.cell(i,j) -= ivar * tmp_cov.cell(i,k) * tmp_cov.cell(j,k);
				}
		}
		delete[](cho_vec);
		delete[](buffer);
		
		return 0.0f;
		
	}
#ifdef GNU_SCIENTIFIC_LIBRARY
LFHTEMP double GaussElem< Tuple<C,0u>, 0 >::Pvalue_Hotelling(const GaussElem< Tuple<C,0u>, 0u > &other, bool LogPvalue) const{
	// Assumes that all weights sumed were 1 (or equal), so that n = w*w/w2
	
	double n[2];
	n[0] = w*w/w2;
	n[1] = other.w*other.w/other.w2;
	if ((mean.getSize() != other.mean.getSize())||(mean.getSize() == 0)||(w == 0)||(other.w == 0)) return 1.0f; // invalid, or no data for one of the two; 
	
	if (n[0] + n[1] <= 1.0f + mean.getSize()){
		// invalid, too little samples to get any pvalue
		return (LogPvalue) ? 1.0f : 0.0f;
	}else{
		Trianglix<C,0u> tcov = (cov - Trianglix<C,0u>(mean * pow(w,-0.5f) )) * (n[0] / w2) + (other.cov - Trianglix<C,0u>(other.mean * pow(other.w,-0.5f))) * (n[1] / other.w2);
		Tuple<C,0u> dif = (mean / w) - (other.mean / other.w);
		C inner = tcov.Xformed_inner_product_of_inverse(dif);
		
		double t = inner * n[0] * n[1] * (n[0] + n[1] - mean.getSize() -1) / (mean.getSize() * (n[0] + n[1])) ;
		
		n[1] = n[0] + n[1] - mean.getSize() - 1;
		
		if (t < 0.0f) {printf("Pvalue_Hotelling has failed (singular matrix?) ...\n");return (LogPvalue) ? 0.0f : 1.0f;}
		return (LogPvalue) ? LogPvalue_Fdistrib_Ptail(t,mean.getSize(), n[1]) : Pvalue_Fdistrib_Ptail(t,mean.getSize(), n[1] );
	}
}

LFHTEMP double GaussElem< Tuple<C,0u>, 0 >::Pvalue_LRTest_UnequalMean(const GaussElem< Tuple<C, 0u>, 0u > &other, const Tuple<bool>* channel_filter, bool LogPvalue) const{ // Ho = mu_a = mu_b V_a = V_b   H1 = mu_a != mu_b V_a = V_b
	unsigned int d;
	if ((w == 0)||(other.w == 0)||((d = mean.getSize()) != other.mean.getSize())) return (LogPvalue) ? 0.0f : 1.0f;
	const GaussElem<Tuple<C,0u>, 0 > * targ[2];
	if (channel_filter){
		targ[0] = new GaussElem<Tuple<C,0u>, 0 >(this->mkSubGaussElem(*channel_filter));
		targ[1] = new GaussElem<Tuple<C,0u>, 0 >(other.mkSubGaussElem(*channel_filter));
	}else{
		targ[0] = this;
		targ[1] = &other;
	}
	
	
	Trianglix<C, 0u> mid_covar = ((*this).cov + other.cov - Trianglix<C>(this->mean * pow(this->w, -0.5f)) - Trianglix<C>(other.mean * pow(other.w, -0.5f))) / (w + other.w);
	
	double middleLL = (w + other.w) * (d * log(2 * M_PI) + mid_covar.log_determinant());	
	
	Trianglix<C, 0u> union_inv = mid_covar.mkinverse();
	
	middleLL += (cov - Trianglix<C,0u>(mean * pow(w,-0.5f))).trace_of_product(union_inv);
	middleLL += (other.cov - Trianglix<C,0u>(other.mean * pow(other.w,-0.5f))).trace_of_product(union_inv);
	
	double ratio = -0.5f * middleLL - ((*this) + other).LLikelihood();
	
	if (channel_filter) {delete(targ[0]);delete(targ[1]);}
	if (ratio >= 0.0f) return (LogPvalue) ? LogPvalue_chisquarre_Ptail(2.0f*ratio, 0.5f * d) : Pvalue_chisquarre_Ptail(2.0f*ratio, 0.5f * d);
	else return (LogPvalue) ? 0.0f : 1.0f;
}

LFHTEMP double GaussElem< Tuple<C,0u>, 0 >::Pvalue_LRTest_UnequalMeanandVariance(const GaussElem< Tuple<C, 0u>, 0u > &other, const Tuple<bool>* channel_filter, bool LogPvalue) const{ // Ho = mu_a = mu_b V_a = V_b   H1 = mu_a != mu_b V_a != V_b
	unsigned int d;
	if ((w == 0)||(other.w == 0)) return (LogPvalue) ? 0.0f : 1.0f;
	const GaussElem<Tuple<C,0u>, 0 > * targ[2];
	if (channel_filter){
		targ[0] = new GaussElem<Tuple<C,0u>, 0 >(this->mkSubGaussElem(*channel_filter));
		targ[1] = new GaussElem<Tuple<C,0u>, 0 >(other.mkSubGaussElem(*channel_filter));
	}else{
		targ[0] = this;
		targ[1] = &other;
	} d = targ[0]->mean.getSize();
	double ratio = targ[0]->LLikelihood() + targ[1]->LLikelihood() - ((*targ[0]) + (*targ[1])).LLikelihood();
	if (channel_filter) {delete(targ[0]);delete(targ[1]);}
	if (ratio >= 0.0f) return (LogPvalue) ? LogPvalue_chisquarre_Ptail(2.0f*ratio, 0.5f * ((d+3)*d)) : Pvalue_chisquarre_Ptail(2.0f*ratio, 0.5f * ((d+3)*d));
	else return (LogPvalue) ? 0.0f : 1.0f;
}

LFHTEMP double GaussElem< Tuple<C,0u>, 0 >::Pvalue_LRTest_UnequalVariance(const GaussElem< Tuple<C, 0u>, 0u > &other, const Tuple<bool>* channel_filter, bool LogPvalue) const{ // Ho = mu_a != mu_b V_a = V_b   H1 = mu_a != mu_b V_a != V_b
	unsigned int d;
	if ((w == 0)||(other.w == 0)) return (LogPvalue) ? 0.0f : 1.0f;
	const GaussElem<Tuple<C,0u>, 0 > * targ[2];
	if (channel_filter){
		targ[0] = new GaussElem<Tuple<C,0u>, 0 >(this->mkSubGaussElem(*channel_filter));
		targ[1] = new GaussElem<Tuple<C,0u>, 0 >(other.mkSubGaussElem(*channel_filter));
	}else{
		targ[0] = this;
		targ[1] = &other;
	} d = targ[0]->mean.getSize();
	// H0 likelihood, unequal mean allowed, which yeild Max Likelihood when matching the sample mean, independently of the fixed Covariance Matrix.

	Trianglix<C, 0u> mid_covar = (targ[0]->cov + targ[1]->cov - Trianglix<C>(targ[0]->mean * pow(targ[0]->w, -0.5f)) - Trianglix<C>(targ[1]->mean * pow(targ[1]->w, -0.5f))) / (targ[0]->w + targ[1]->w);
	
	Trianglix<C, 0u> union_inv = mid_covar.mkinverse();
	double middleLL =  (targ[0]->w + targ[1]->w) * (d * log(2 * M_PI) + mid_covar.log_determinant());
	
	middleLL += (targ[0]->cov - Trianglix<C,0u>(targ[0]->mean * pow(targ[0]->w,-0.5f))).trace_of_product(union_inv);
	middleLL += (targ[1]->cov - Trianglix<C,0u>(targ[1]->mean * pow(targ[1]->w,-0.5f))).trace_of_product(union_inv);
	
	double ratio = targ[0]->LLikelihood() + targ[1]->LLikelihood() + 0.5f * middleLL;
	if (channel_filter) {delete(targ[0]);delete(targ[1]);}
	if (ratio >= 0.0f) return (LogPvalue) ? LogPvalue_chisquarre_Ptail(2.0f*ratio, 0.5f * ((d+1)*d)) : Pvalue_chisquarre_Ptail(2.0f*ratio, 0.5f * ((d+1)*d));
	else return (LogPvalue) ? 0.0f : 1.0f;
}
#endif
LFHTEMP double GaussElem< Tuple<C,0u>, 0 >::getLogFoldinVariance_inMeanMajorAxis(const GaussElem< Tuple<C, 0u>, 0u > &x) const{
	Tuple<C> meandif = this->getMean() - x.getMean();
	meandif /= ExOp::norm(meandif);
	return 0.5f * (ExOp::lognorm(this->getCovariance_biased().Xformed_inner_product(meandif)) - ExOp::lognorm(x.getCovariance_biased().Xformed_inner_product(meandif)));
}

LFHTEMP double GaussElem< Tuple<C,0u>, 0 >::getLogFoldinVariance_Determinant(const GaussElem< Tuple<C, 0u>, 0u > &x) const{
	double fout = (ExOp::norm(this->getCovariance_biased().log_determinant()) - ExOp::norm(x.getCovariance_biased().log_determinant())) * (0.5f / this->getSize());
	if (ExOp::isValid(fout)) return fout;
	Trianglix<double> tmptmp = this->getCovariance_biased();
	printf("\t%e\t", tmptmp.data[0] *tmptmp.data[2] - tmptmp.data[1] * tmptmp.data[1]);
	printf("error: %f\t%f\t%i\t%i\n", this->getCovariance_biased().log_determinant(), x.getCovariance_biased().log_determinant(), this->getSize(), x.getSize());
	return 0.0f;
}

LFHTEMP double GaussElem< Tuple<C,0u>, 0 >::getLogFoldinVariance_Diagonal(const GaussElem< Tuple<C, 0u>, 0u > &x) const{
	double fout = 0;
	Trianglix<double> tmptmp = this->getCovariance_biased();
	for(int i =0; i< tmptmp.getSize();i++) fout += log(tmptmp.cell(i,i));
	tmptmp = x.getCovariance_biased();
	for(int i =0; i< tmptmp.getSize();i++) fout -= log(tmptmp.cell(i,i));
	if (ExOp::isValid(fout)) return (fout * (0.5f / this->getSize()));
	return 0.0f;
}

LFHTEMP double GaussElem< Tuple<C,0u>, 0 >::getDifferenceinMean_inVarianceMajorAxis(const GaussElem< Tuple<C, 0u>, 0u > &x) const{
	
}


LFHTEMP GaussElem< Tuple<C, 0u>, 0 >& GaussElem< Tuple<C,0u>, 0 >::addGaussElem_free_mean(const GaussElem< Tuple<C, 0u>, TARG2 > &other){
	if (w > 0.0f) cov -= Trianglix<C, 0u>( mean * pow(w, -0.5f));
	w+=other.w; w2+=other.w2; mean += other.mean; cov += other.getCovariance_biased() * other.w;
	if (w > 0.0f) cov += Trianglix<C, 0u>( mean * pow(w, -0.5f));	
	ExOp::toZero(determinant); return(*this);
}

/*
LFHTEMP double GaussElem< Tuple<C,0u>, 0 >::Pvalue_LRTest_UnequalVariance_Masked(const GaussElem< Tuple<C, 0u>, 0u > &other, const Trianglix<bool> &mask, bool LogPvalue) const{ // Ho = mu_a != mu_b V_a = V_b   H1 = mu_a != mu_b V_a != V_b
	unsigned int d;
	if ((w == 0)||(other.w == 0)||((d = mean.getSize()) != other.mean.getSize())) return (LogPvalue) ? 0.0f : 1.0f;
	GaussElem<Tuple<C,0u>, 0 > * targ[2];
	if (channel_filter){
		targ[0] = new GaussElem<Tuple<C,0u>, 0 >(this->mkSubGaussElem(*channel_filter));
		targ[1] = new GaussElem<Tuple<C,0u>, 0 >(other.mkSubGaussElem(*channel_filter));
	}else{
		targ[0] = this;
		targ[1] = &other;
	} d = targ[0]->mean.getSize();	
	Trianglix<C, 0u> union_cov = ((*this) + other).getCovariance_biased();
	Trianglix<C, 0u> union_inv = union_cov.inverse();
	double middleLL =  (w + other.w) * (d * log(2 * M_PI) + union_cov.log_determinant());
	
	Trianglix<C,0u> covA = cov - Trianglix<C,0u>(mean * pow(w,-0.5f));
	Trianglix<C,0u> covB = other.cov - Trianglix<C,0u>(other.mean * pow(other.w,-0.5f));
	
	middleLL += covA.trace_of_product(union_inv);
	middleLL += covB.trace_of_product(union_inv);

	double ratio = 0.5f * middleLL;
	if (ratio >= 0.0f) return (LogPvalue) ? LogPvalue_chisquarre_Ptail(2.0f*ratio, 0.5f * ((d+1)*d)) : Pvalue_chisquarre_Ptail(2.0f*ratio, 0.5f * ((d+1)*d));
	else return (LogPvalue) ? 0.0f : 1.0f;
}*/


	LFHTEMP void GaussElem< Tuple<C,0u>, 0 >::load(FILE *f, unsigned int size){
		fread(&w,sizeof(double),1,f);
		fread(&w2,sizeof(double),1,f);
		ExOp::load(determinant,f);
		ExOp::load(mean,f);
		ExOp::load(cov,f);
	}
	
	LFHTEMP void GaussElem< Tuple<C,0u>, 0 >::save(FILE *f) const{
		fwrite(&w,sizeof(double),1,f);
		fwrite(&w2,sizeof(double),1,f);
		C t = determinant; ExOp::toZero(t);
		if (determinant == t) {
			determinant = getCovariance().log_determinant();
		}
		ExOp::save(determinant,f);
		ExOp::save(mean,f);
		ExOp::save(cov,f);
	}

	template <class C> const C& getValue(int flagindex){ C tmp; exit(1); return(tmp);}


	template <int buffersize>
	SuperString<buffersize>::SuperString(){
		buffer[0] = '\0';
	}

	template <int buffersize>
	char* SuperString<buffersize>::operator()(){
		return(buffer);
	}

	template <int buffersize>
	void SuperString<buffersize>::setChunk(int chunkId, char* string, int stringlength){
		int i = chunks.size();
		while(i <= chunkId ) {
			if (i == 0) chunks.push_back(0);
			else chunks.push_back(chunks[i-1]);
			i++;
		}

		int j = (stringlength == 0) ? strlen(string): stringlength;

		int d = ((chunkId == 0) ? 0 : chunks[chunkId-1]);
		int sh = j - chunks[chunkId] + d;
		int e = chunks[chunkId];

		int c = chunks[i-1];

		if (sh > 0){
			for(;c>=e;c--) buffer[c+sh] = buffer[c];
			for(c = chunkId;c<i;c++) chunks[c] += sh;
			memcpy(buffer+d, string, sizeof(char)*j);
		}else if (sh != 0) {
			memcpy(buffer+d, string, sizeof(char)*j);
			for(;c>=e;e++) buffer[e+sh] = buffer[e];
			for(c = chunkId;c<i;c++) chunks[c] += sh;
		} else memcpy(buffer+d, string, sizeof(char)*j);
	}

	template <int buffersize>
	char* SuperString<buffersize>::operator[](int chunkId){
		return((chunkId == 0) ? buffer : buffer + chunks[chunkId-1]);
	}


	template <class C, unsigned int channels, Tuple_flag Cflag>
	void TiffFile::fetch(vector<DataGrid<Tuple<C,channels,Cflag>,2> > & _out){
		char buffer[65536];
		unsigned short nbflag;
		int i;
		int typesize;
		DataGrid<Tuple<C,channels>,2> *cur;
		unsigned int size[2];
		unsigned int nbchannels;
		unsigned int stripof, stripro, stripco;
		unsigned int required;
		int tmp_iscompressed;
		while(gotoNext()){
			required =0; nbchannels= 1;
			for(i=0;i<curflaglist.getSize()/12;i++){
				switch (*(unsigned short*)&(curflaglist[i*12])){
					case 256: size[0] = getValue<unsigned int>(i); required |= 1; break;//memcpy(buffer+512,"ImageWidth",sizeof(char)*11);break;
					case 257: size[1] = getValue<unsigned int>(i); required |= 2;break;//memcpy(buffer+512,"ImageHeight",sizeof(char)*12);break;
					case TIFFFLAG_Compression: tmp_iscompressed = getValue<int>(i); if ((tmp_iscompressed != 0)&&(tmp_iscompressed != 1)) {fprintf(stderr,"Tiff file is saved in a compressed format, it needs to be uncompressed.\n"); exit(1);} break; 
					case 258: break;
					case 273: stripof = getValue<unsigned int>(i); required |= 8;break;//memcpy(buffer+512,"SkipOffsets",sizeof(char)*12);break;
					case 277: nbchannels= getValue<unsigned int>(i); break;//memcpy(buffer+512,"SamplesPerPixel",sizeof(char)*16);break;
					case 278: stripro= getValue<unsigned int>(i); required |= 16;break;//memcpy(buffer+512,"RowsPerStrip",sizeof(char)*13);break;
					case 279: stripco= getValue<unsigned int>(i); required |= 32;break;//memcpy(buffer+512,"StripByteCounts",sizeof(char)*16);break;
				}
			}
			if (required == 59){
				_out.push_back(DataGrid<Tuple<C,channels>,2>());
				cur = &(_out[_out.getSize()-1]);
				cur->setSizes(size);
				fseek(f,stripof,SEEK_SET);
			}
		}
	}

	template <class C>
	bool TiffFile::fetch(DataGrid<C,3>& f_out, char * type_out){
//		char buffer[65536];
		unsigned int i,j,k;
		typedef class DataGrid<C,3>::KeyIterator itetype;
		itetype ite = f_out.getKeyIterator();
		unsigned int size[3];
		unsigned int stripro;

		vector<unsigned int> stripof, stripco;
		unsigned int required;
		unsigned int unitsize;
		unsigned int type =1;
		unsigned short tmpflag;
		int tmp_iscompressed;
		if (!gotoNext()){
			return(false);
		}else{
			required =0; size[0] = 1;
			for(i=0;i<curflaglist.getSize()/12;i++){
				tmpflag = *(unsigned short*)&(curflaglist[i*12]);
				if (inv) ExOp::bytereverse(tmpflag);
				switch (tmpflag){
					case 256: size[1] = getValue<unsigned int>(i); required |= 1; break;//memcpy(buffer+512,"ImageWidth",sizeof(char)*11);break;
					case 257: size[2] = getValue<unsigned int>(i); required |= 2;break;//memcpy(buffer+512,"ImageHeight",sizeof(char)*12);break;
					case TIFFFLAG_Compression: tmp_iscompressed = getValue<int>(i); if ((tmp_iscompressed != 0)&&(tmp_iscompressed != 1)) {fprintf(stderr,"Tiff file is saved in a compressed format, it needs to be uncompressed.\n"); exit(1);} break; 
					case 258: unitsize = getValue<unsigned int>(i); break;
					case 273: stripof = getValue< vector<unsigned int> >(i); required |= 8;break;//memcpy(buffer+512,"SkipOffsets",sizeof(char)*12);break;
					case 277: size[0] = getValue<unsigned int>(i); break;//memcpy(buffer+512,"SamplesPerPixel",sizeof(char)*16);break;
					case 278: stripro= getValue<unsigned int>(i); required |= 16;break;//memcpy(buffer+512,"RowsPerStrip",sizeof(char)*13);break;
					case 279: stripco= getValue< vector<unsigned int> >(i); required |= 32;break;//memcpy(buffer+512,"StripByteCounts",sizeof(char)*16);break;
					case 339: type = getValue<unsigned int>(i);	break; //SampleFormat
				}
			}
				if (required == 59){
					f_out.setSizes(size);
				//	printf("%i,,,%i,,,%i\n",stripro, size[2] ,stripco.getSize());
					fseek(f,stripof[0],SEEK_SET); k=j=i=0;
                //    printf("datype %i! %i\n",type,unitsize);
					switch(type){
						case 1: // unsigned int
                          //  printf("ualllllo! %i\n",unitsize);
							if (unitsize == sizeof(unsigned char)*8){
								if (type_out) type_out[0] = 'C';
								//printf("RecognizedQ unsigned char!\n");
								{ //:
									unsigned char* buffer = new unsigned char[stripro * size[1]* size[0]];
									fread(buffer,sizeof(unsigned char), stripro * size[1]* size[0],f);

								if (ite.first()) do{
									f_out.data[j++] = (C) buffer[k++];
									if (k >= stripro * size[1]* size[0]){
										i++;
										fseek(f,stripof[i],SEEK_SET);
										fread(buffer,sizeof(unsigned char), stripro * size[1]* size[0],f);
										k=0;
									}
								} while(ite.next());
									delete[](buffer);

								} //:
							}else if (unitsize == sizeof(unsigned short)*8){
								{ //:
									//printf("RecognizedQ unsigned short!\n");
									if (type_out) type_out[0] = 'S';
									unsigned short* buffer = new unsigned short[stripro * size[1]* size[0]];
									fread(buffer,sizeof(unsigned short), stripro * size[1]* size[0],f);
									if (ite.first()) do{
										if (inv) ExOp::bytereverse(buffer[k]);
										f_out.data[j++] = (C) buffer[k++];
										if (k >= stripro * size[1]* size[0]){
											i++;
											fseek(f,stripof[i],SEEK_SET);
											fread(buffer,sizeof(unsigned short), stripro * size[1]* size[0],f);
											k=0;
										}
									} while(ite.next());
									delete[](buffer);
								} //:
							}else if (unitsize == sizeof(unsigned int)*8){
								{ //:
									if (type_out) type_out[0] = 'I'; //printf("helllllo!\n");
									//printf("RecognizedQ unsigned int!\n");
									unsigned int* buffer = new unsigned int[stripro * size[1]* size[0]];
									fread(buffer,sizeof(unsigned int), stripro * size[1]* size[0],f);

									if (ite.first()) do{
										if (inv) ExOp::bytereverse(buffer[k]);
										f_out.data[j++] = (C) buffer[k++];
										if (k >= stripro * size[1]* size[0]){
											i++;
											fseek(f,stripof[i],SEEK_SET);
											fread(buffer,sizeof(unsigned int), stripro * size[1]* size[0],f);
											k=0;
										}
									} while(ite.next());
									delete[](buffer);

								} //:
							}
						break;
						case 2: // signed int
                            // printf("salllllo! %i\n",unitsize);
							if (unitsize == sizeof(unsigned char)*8){
								if (type_out) type_out[0] = 'c';
								//printf("RecognizedQ unsigned char!\n");
								{ //:
									char* buffer = new char[stripro * size[1]* size[0]];
									fread(buffer,sizeof(char), stripro * size[1]* size[0],f);
									
									if (ite.first()) do{
										f_out.data[j++] = (C) buffer[k++];
										if (k >= stripro * size[1]* size[0]){
											i++;
											fseek(f,stripof[i],SEEK_SET);
											fread(buffer,sizeof(char), stripro * size[1]* size[0],f);
											k=0;
										}
									} while(ite.next());
									delete[](buffer);
									
								} //:
							}else if (unitsize == sizeof(unsigned short)*8){
								{ //:
									//printf("RecognizedQ unsigned short!\n");
									if (type_out) type_out[0] = 's';
									short* buffer = new short[stripro * size[1]* size[0]];
									fread(buffer,sizeof(short), stripro * size[1]* size[0],f);
									
									if (ite.first()) do{
										if (inv) ExOp::bytereverse(buffer[k]);
										f_out.data[j++] = (C) buffer[k++];
										if (k >= stripro * size[1]* size[0]){
											i++;
											fseek(f,stripof[i],SEEK_SET);
											fread(buffer,sizeof(short), stripro * size[1]* size[0],f);
											k=0;
										}
									} while(ite.next());
									delete[](buffer);
								} //:
							}else if (unitsize == sizeof(unsigned int)*8){
								{ //:
									if (type_out) type_out[0] = 'i'; //printf("alllllo!\n");
									//printf("RecognizedQ unsigned int!\n");
									int* buffer = new int[stripro * size[1]* size[0]];
									fread(buffer,sizeof(int), stripro * size[1]* size[0],f);
									
									if (ite.first()) do{
										if (inv) ExOp::bytereverse(buffer[k]);
										f_out.data[j++] = (C) buffer[k++];
										if (k >= stripro * size[1]* size[0]){
											i++;
											fseek(f,stripof[i],SEEK_SET);
											fread(buffer,sizeof(int), stripro * size[1]* size[0],f);
											k=0;
										}
									} while(ite.next());
									delete[](buffer);
								} //:
							}
						break;
						case 3: // float
							if (unitsize == sizeof(float)*8){
								{ //:
									if (type_out) type_out[0] = 'f';
									//printf("RecognizedQ float!\n");
									float* buffer = new float[stripro * size[1]* size[0]];
									fread(buffer,sizeof(float), stripro * size[1]* size[0],f);

									if (ite.first()) do{
										f_out.data[j++] = (C) buffer[k++];
										if (k >= stripro * size[1]* size[0]){
											i++;
										//	printf("%i!!!%i!!!%i\n",stripro, size[2] ,stripco.getSize());
											fseek(f,stripof[i],SEEK_SET);
											fread(buffer,sizeof(float), stripro * size[1]* size[0],f);
											k=0;
										}
									} while(ite.next());
									delete[](buffer);

								} //:
							}else if (unitsize == sizeof(double)*8){
								{ //:
									if (type_out) type_out[0] = 'd';
									//printf("RecognizedQ double!\n");
									double* buffer = new double[stripro * size[1]* size[0]];
									fread(buffer,sizeof(double), stripro * size[1]* size[0],f);

									if (ite.first()) do{
										f_out.data[j++] = (C) buffer[k++];
										if (k >= stripro * size[1]){
											i++;
											fseek(f,stripof[i],SEEK_SET);
											fread(buffer,sizeof(double), stripro * size[1]* size[0],f);
											k=0;
										}
									} while(ite.next());
									delete[](buffer);
								} //:
							}
						break;
					}
				}
		}
		return(true);
	}

	template <class C>
	bool TiffFile::fetch(DataGrid<C,2>& f_out, const unsigned int channel){
		//		char buffer[65536];
		unsigned int i,j,k;


		typedef class DataGrid<C,2>::KeyIterator itetype;
		itetype ite = f_out.getKeyIterator();


		unsigned int size[3];
		unsigned int stripro;

		vector<unsigned int> stripof, stripco;
		unsigned int required;
		unsigned int unitsize;
		unsigned int type =1;
		unsigned short tmpflag;
		int tmp_iscompressed;
		if (!gotoNext()){
			return(false);
		}else{
			required =0; size[0] = 1;
			for(i=0;i<curflaglist.getSize()/12;i++){
				//	printf("flagtype %i\n", *(unsigned short*)&(curflaglist[i*12])); fflush(stdout);
				tmpflag = *(unsigned short*)&(curflaglist[i*12]);
				if (inv) ExOp::bytereverse(tmpflag);

				switch (tmpflag){
					case 256: size[1] = getValue<unsigned int>(i); required |= 1; break;//memcpy(buffer+512,"ImageWidth",sizeof(char)*11);break;
					case 257: size[2] = getValue<unsigned int>(i); required |= 2;break;//memcpy(buffer+512,"ImageHeight",sizeof(char)*12);break;
					case TIFFFLAG_Compression: tmp_iscompressed = getValue<int>(i); if ((tmp_iscompressed != 0)&&(tmp_iscompressed != 1)) {fprintf(stderr,"Tiff file is saved in a compressed format, it needs to be uncompressed.\n"); exit(1);} break; 
					case 258: unitsize = getValue<unsigned int>(i); break;
					case 273: stripof = getValue< vector<unsigned int> >(i); required |= 8;break;//memcpy(buffer+512,"SkipOffsets",sizeof(char)*12);break;
					case 277: size[0] = getValue<unsigned int>(i); break;//memcpy(buffer+512,"SamplesPerPixel",sizeof(char)*16);break;
					case 278: stripro= getValue<unsigned int>(i); required |= 16;break;//memcpy(buffer+512,"RowsPerStrip",sizeof(char)*13);break;
					case 279: stripco= getValue< vector<unsigned int> >(i); required |= 32;break;//memcpy(buffer+512,"StripByteCounts",sizeof(char)*16);break;
					case 339: type = getValue<unsigned int>(i);	break; //SampleFormat
				}
			}
			if (required == 59){
				f_out.setSizes(size+1);
				//	printf("%i,,,%i,,,%i\n",stripro, size[2] ,stripco.getSize());
				fseek(f,stripof[0],SEEK_SET); k=j=i=0;
				switch(type){
					case 1: // unsigned int
						if (unitsize == sizeof(unsigned char)*8){
							//printf("RecognizedQ unsigned char!\n");
							{ //:
								unsigned char* buffer = new unsigned char[stripro * size[1]];
								fread(buffer,sizeof(unsigned char), stripro * size[1],f);

								if (ite.first()) do{
									f_out.data[j++] = (C) buffer[k++];
									if (k >= stripro * size[1]* size[0]){
										i++;
										fseek(f,stripof[i],SEEK_SET);
										if (size[0] == 1) fread(buffer,sizeof(unsigned char), stripro * size[1],f);
										else{

										}
										k=0;
									}
								} while(ite.next());
								delete[](buffer);

							} //:
						}else if (unitsize == sizeof(unsigned short)*8){
							{ //:
								//printf("RecognizedQ unsigned short!\n");

								unsigned short* buffer = new unsigned short[stripro * size[1]];
								fread(buffer,sizeof(unsigned short), stripro * size[1],f);

								if (ite.first()) do{
									f_out.data[j++] = (C) buffer[k++];
									if (k >= stripro * size[1]* size[0]){
										i++;
										fseek(f,stripof[i],SEEK_SET);
										if (size[0] == 1) fread(buffer,sizeof(unsigned short), stripro * size[1],f);
										else{// TODO!

										}
										k=0;
									}
								} while(ite.next());
								delete[](buffer);
							} //:
						}else if (unitsize == sizeof(unsigned int)*8){
							{ //:
								//printf("RecognizedQ unsigned int!\n");
								unsigned int* buffer = new unsigned int[stripro * size[1]];
								fread(buffer,sizeof(unsigned int), stripro * size[1],f);

								if (ite.first()) do{
									f_out.data[j++] = (C) buffer[k++];
									if (k >= stripro * size[1]* size[0]){
										i++;
										fseek(f,stripof[i],SEEK_SET);
										if (size[0] == 1) fread(buffer,sizeof(unsigned int), stripro * size[1],f);
										else{ // TODO!
										}
										k=0;
									}
								} while(ite.next());
								delete[](buffer);

							} //:
						}else if (unitsize == 64){
						}
						break;
					case 2: // signed int
						break;
					case 3: // float
						if (unitsize == sizeof(float)*8){
							{ //:
								//printf("RecognizedQ float!\n");
								float* buffer = new float[stripro * size[1]];
								fread(buffer,sizeof(float), stripro * size[1],f);

								if (ite.first()) do{
									f_out.data[j++] = (C) buffer[k++];
									if (k >= stripro * size[1]* size[0]){
										i++;
										//	printf("%i!!!%i!!!%i\n",stripro, size[2] ,stripco.getSize());
										fseek(f,stripof[i],SEEK_SET);
										if (size[0] == 1) fread(buffer,sizeof(float), stripro * size[1],f);
										else{

										}
										k=0;
									}
								} while(ite.next());
								delete[](buffer);

							} //:
						}else if (unitsize == sizeof(double)*8){
							{ //:
								//printf("RecognizedQ double!\n");
								double* buffer = new double[stripro * size[1]];
								fread(buffer,sizeof(double), stripro * size[1],f);

								if (ite.first()) do{
									f_out.data[j++] = (C) buffer[k++];
									if (k >= stripro * size[1]){
										i++;
										fseek(f,stripof[i],SEEK_SET);
										if (size[0] == 1) fread(buffer,sizeof(double), stripro * size[1],f);
										else{
										}
										k=0;
									}
								} while(ite.next());
								delete[](buffer);
							} //:
						}
						break;
				}
			}
		}
		return(true);
	}

	template <class C>
	bool TiffFile::fetch_grayscale(DataGrid<C,2>& f_out){
		//		char buffer[65536];
		unsigned int i,j,k;


		typedef class DataGrid<C,2>::KeyIterator itetype;
		itetype ite = f_out.getKeyIterator();


		unsigned int size[3];
		unsigned int stripro;

		vector<unsigned int> stripof, stripco;
		unsigned int required;
		unsigned int unitsize;
		unsigned int type =1;
		unsigned short tmpflag;
		int tmp_iscompressed;
		double value;
		if (!gotoNext()){
			return(false);
		}else{
			required =0; size[0] =1;
			for(i=0;i<curflaglist.getSize()/12;i++){
				//	printf("flagtype %i\n", *(unsigned short*)&(curflaglist[i*12])); fflush(stdout);
				tmpflag = *(unsigned short*)&(curflaglist[i*12]);
				if (inv) ExOp::bytereverse(tmpflag);

				switch (tmpflag){
					case 256: size[1] = getValue<unsigned int>(i); required |= 1; break;//memcpy(buffer+512,"ImageWidth",sizeof(char)*11);break;
					case 257: size[2] = getValue<unsigned int>(i); required |= 2;break;//memcpy(buffer+512,"ImageHeight",sizeof(char)*12);break;
					case TIFFFLAG_Compression: tmp_iscompressed = getValue<int>(i); if ((tmp_iscompressed != 0)&&(tmp_iscompressed != 1)) {fprintf(stderr,"Tiff file is saved in a compressed format, it needs to be uncompressed.\n"); exit(1);} break; 
					case 258: unitsize = getValue<unsigned int>(i); break;
					case 273: stripof = getValue< vector<unsigned int> >(i); required |= 8;break;//memcpy(buffer+512,"SkipOffsets",sizeof(char)*12);break;
					case 277: size[0] = getValue<unsigned int>(i); break;//memcpy(buffer+512,"SamplesPerPixel",sizeof(char)*16);break;
					case 278: stripro= getValue<unsigned int>(i); required |= 16;break;//memcpy(buffer+512,"RowsPerStrip",sizeof(char)*13);break;
					case 279: stripco= getValue< vector<unsigned int> >(i); required |= 32;break;//memcpy(buffer+512,"StripByteCounts",sizeof(char)*16);break;
					case 339: type = getValue<unsigned int>(i);	break; //SampleFormat
				}
			}
			if (required == 59){
				f_out.setSizes(size+1);
				//	printf("%i,,,%i,,,%i\n",stripro, size[2] ,stripco.getSize());
				fseek(f,stripof[0],SEEK_SET); k=j=i=0;
				switch(type){
					case 1: // unsigned intTIFFFLAG_Compression
						if (unitsize == sizeof(unsigned char)*8){
							//printf("RecognizedQ unsigned char!\n");
							{ //:
								unsigned char* buffer = new unsigned char[stripro * size[1]* size[0]];
								fread(buffer,sizeof(unsigned char), stripro * size[1]* size[0],f);

								if (ite.first()) do{
								    if (size[0] ==3){
								        value = 0.30f * pow(buffer[k++], 2.2f);
								        value += 0.59f * pow(buffer[k++], 2.2f);
								        value += 0.11f * pow(buffer[k++], 2.2f);
                                        f_out.data[j++] = (C) pow(value, 1.0f / 2.2f);
								    }else {
                                        f_out.data[j++] = buffer[k]; k+= size[0];
								    }
									if (k >= stripro * size[1]* size[0]){
										i++;
										fseek(f,stripof[i],SEEK_SET);
										fread(buffer,sizeof(unsigned char), stripro * size[1]* size[0],f);
										k=0;
									}
								} while(ite.next());
								delete[](buffer);

							} //:
						}else if (unitsize == sizeof(unsigned short)*8){
							{ //:
								//printf("RecognizedQ unsigned short!\n");

								unsigned short* buffer = new unsigned short[stripro * size[1]* size[0]];
								fread(buffer,sizeof(unsigned short), stripro * size[1],f);

								if (ite.first()) do{
									f_out.data[j++] = (C) buffer[k++];
									if (k >= stripro * size[1]* size[0]){
										i++;
										fseek(f,stripof[i],SEEK_SET);
										if (size[0] == 1) fread(buffer,sizeof(unsigned short), stripro * size[1],f);
										else{// TODO!

										}
										k=0;
									}
								} while(ite.next());
								delete[](buffer);
							} //:
						}else if (unitsize == sizeof(unsigned int)*8){
							{ //:
								//printf("RecognizedQ unsigned int!\n");
								unsigned int* buffer = new unsigned int[stripro * size[1]];
								fread(buffer,sizeof(unsigned int), stripro * size[1],f);

								if (ite.first()) do{
									f_out.data[j++] = (C) buffer[k++];
									if (k >= stripro * size[1]* size[0]){
										i++;
										fseek(f,stripof[i],SEEK_SET);
										if (size[0] == 1) fread(buffer,sizeof(unsigned int), stripro * size[1],f);
										else{ // TODO!
										}
										k=0;
									}
								} while(ite.next());
								delete[](buffer);

							} //:
						}else if (unitsize == 64){
						}
						break;
					case 2: // signed int
						break;
					case 3: // float
						if (unitsize == sizeof(float)*8){
							{ //:
								//printf("RecognizedQ float!\n");
								float* buffer = new float[stripro * size[1]];
								fread(buffer,sizeof(float), stripro * size[1],f);

								if (ite.first()) do{
									f_out.data[j++] = (C) buffer[k++];
									if (k >= stripro * size[1]* size[0]){
										i++;
										//	printf("%i!!!%i!!!%i\n",stripro, size[2] ,stripco.getSize());
										fseek(f,stripof[i],SEEK_SET);
										if (size[0] == 1) fread(buffer,sizeof(float), stripro * size[1],f);
										else{

										}
										k=0;
									}
								} while(ite.next());
								delete[](buffer);

							} //:
						}else if (unitsize == sizeof(double)*8){
							{ //:
								//printf("RecognizedQ double!\n");
								double* buffer = new double[stripro * size[1]];
								fread(buffer,sizeof(double), stripro * size[1],f);

								if (ite.first()) do{
									f_out.data[j++] = (C) buffer[k++];
									if (k >= stripro * size[1]){
										i++;
										fseek(f,stripof[i],SEEK_SET);
										if (size[0] == 1) fread(buffer,sizeof(double), stripro * size[1],f);
										else{
										}
										k=0;
									}
								} while(ite.next());
								delete[](buffer);
							} //:
						}
						break;
				}
			}
		}
		return(true);
	}

	template<class TYPE> void TiffFile::writeDescriptor(TYPE min, TYPE max, char * &p){
		if (isTypeEqual<TYPE,unsigned char>::ans){
			*(short*)(p) = 339;
			*(short*)(p+2) = 3;
			*(int*)(p+4) = 1;
			*(short*)(p+8) = 1; //format
			p+=12;

			*(short*)(p) = 340;
			*(short*)(p+2) = 4;
			*(int*)(p+4) = 1;
			*(int*)(p+8) = 0; //min value
			p+=12;

			*(short*)(p) = 341;
			*(short*)(p+2) = 4;
			*(int*)(p+4) = 1;
			*(unsigned int*)(p+8) = 255; //max value

			p+=12;
		}else if (isTypeEqual<TYPE,char>::ans){
			*(short*)(p) = 339;
			*(short*)(p+2) = 3;
			*(int*)(p+4) = 1;
			*(short*)(p+8) = 2; //format
			p+=12;

			*(short*)(p) = 340;
			*(short*)(p+2) = 4;
			*(int*)(p+4) = 1;
			*(int*)(p+8) = -128; //min value
			p+=12;

			*(short*)(p) = 341;
			*(short*)(p+2) = 4;
			*(int*)(p+4) = 1;
			*(unsigned int*)(p+8) = 127; //max value

			p+=12;
		}else if (isTypeEqual<TYPE,unsigned short>::ans){
			*(short*)(p) = 339;
			*(short*)(p+2) = 3;
			*(int*)(p+4) = 1;
			*(short*)(p+8) = 1; //format
			p+=12;

			*(short*)(p) = 340;
			*(short*)(p+2) = 4;
			*(int*)(p+4) = 1;
			*(int*)(p+8) = 0; //min value
			p+=12;

			*(short*)(p) = 341;
			*(short*)(p+2) = 4;
			*(int*)(p+4) = 1;
			*(unsigned int*)(p+8) = 65535; //max value
			p+=12;
		}else if (isTypeEqual<TYPE,short>::ans){
			*(short*)(p) = 339;
			*(short*)(p+2) = 3;
			*(int*)(p+4) = 1;
			*(short*)(p+8) = 2; //format
			p+=12;

			*(short*)(p) = 340;
			*(short*)(p+2) = 4;
			*(int*)(p+4) = 1;
			*(int*)(p+8) = -32768; //min value
			p+=12;

			*(short*)(p) = 341;
			*(short*)(p+2) = 4;
			*(int*)(p+4) = 1;
			*(unsigned int*)(p+8) = 32767; //max value
			p+=12;
		}else if (isTypeEqual<TYPE,unsigned int>::ans){
			*(short*)(p) = 339;
			*(short*)(p+2) = 3;
			*(int*)(p+4) = 1;
			*(short*)(p+8) = 1; //format
			p+=12;

			*(short*)(p) = 340;
			*(short*)(p+2) = 4;
			*(int*)(p+4) = 1;
			*(int*)(p+8) = 0; //min value
			p+=12;

			*(short*)(p) = 341;
			*(short*)(p+2) = 4;
			*(int*)(p+4) = 1;
			*(unsigned int*)(p+8) = 0xFFFFFFFF; //max value
			p+=12;
		}else if (isTypeEqual<TYPE,int>::ans){
			*(short*)(p) = 339;
			*(short*)(p+2) = 3;
			*(int*)(p+4) = 1;
			*(short*)(p+8) = 2; //format
			p+=12;

			*(short*)(p) = 340;
			*(short*)(p+2) = 4;
			*(int*)(p+4) = 1;
			*(int*)(p+8) = 0x80000000; //min value
			p+=12;

			*(short*)(p) = 341;
			*(short*)(p+2) = 4;
			*(int*)(p+4) = 1;
			*(int*)(p+8) = 0x7FFFFFFF; //max value
			p+=12;
		}else{
			*(short*)(p) = 339;
			*(short*)(p+2) = 3;
			*(int*)(p+4) = 1;
			*(short*)(p+8) = 3; //format
			p+=12;

			*(short*)(p) = 340;
			*(short*)(p+2) = 4;
			*(int*)(p+4) = 1;
			*(int*)(p+8) = 0; //min value
			p+=12;

			*(short*)(p) = 341;
			*(short*)(p+2) = 4;
			*(int*)(p+4) = 1;
			*(int*)(p+8) = 1; //max value
			p+=12;
		}
	}


	template <class C, class TYPE>
	bool TiffFile::put(DataGrid<C,3>& f_out, TYPE min, TYPE max,bool updateRange){

	//	printf("%i\t%i\t%i\n", f_out.dims[0],f_out.dims[1],f_out.dims[2]);
		if (endfile_pos == 0){ // nothing was written so far
			endfile_pos = curfp_pos + 4;
		}

		fseek(f,curfp_pos,SEEK_SET);
		fwrite(&endfile_pos,sizeof(unsigned int),1, f);
		fseek(f,endfile_pos,SEEK_SET);
		unsigned int i;

		char buffer[65536];
		char *p = buffer;
		int nbflag=17;
		int extra=0;
		Tuple<unsigned int,3> coor;
			*(short*)(p) = nbflag;
			p+=2;

			*(short*)(p) = 254;
			*(short*)(p+2) = 4;
			*(int*)(p+4) = 1;
			*(int*)(p+8) = 0;
			p+=12;

			*(short*)(p) = 256;
			*(short*)(p+2) = 4;
			*(int*)(p+4) = 1;
			*(int*)(p+8) = f_out.dims[1];
			p+=12;

			*(short*)(p) = 257;
			*(short*)(p+2) = 4;
			*(int*)(p+4) = 1;
			*(int*)(p+8) = f_out.dims[2];
			p+=12;

			*(short*)(p) = 258;
			*(short*)(p+2) = 3;
			*(int*)(p+4) = f_out.dims[0];
			if (f_out.dims[0] < 3){
				*(short*)(p+8) = 8* sizeof(TYPE); //bitpersample
				*(short*)(p+10) = 8* sizeof(TYPE); //bitpersample
			}else{
				*(int*)(p+8) = endfile_pos + 6 + 12 * nbflag + extra;
				extra += sizeof(short) * f_out.dims[0];
			}
			p+=12;

			*(short*)(p) = 259;
			*(short*)(p+2) = 3;
			*(int*)(p+4) = 1;
			*(short*)(p+8) = 1; //compression
			p+=12;

			*(short*)(p) = 262;
			*(short*)(p+2) = 3;
			*(int*)(p+4) = 1;
			switch(f_out.dims[0]){//photothing
				case 1: *(short*)(p+8) = 1; break;
				case 2: *(short*)(p+8) = 1; break;
				default: *(short*)(p+8) = 2; break;
			}
			p+=12;

			*(short*)(p) = 282;
			*(short*)(p+2) = 5;
			*(int*)(p+4) = 1;
			*(int*)(p+8) = endfile_pos + 6 + 12 * nbflag + extra;
			extra += sizeof(int)*2;
			p+=12;

			*(short*)(p) = 283;
			*(short*)(p+2) = 5;
			*(int*)(p+4) = 1;
			*(int*)(p+8) = endfile_pos + 6 + 12 * nbflag + extra;
			extra += sizeof(int)*2;
			p+=12;

			*(short*)(p) = 273;
			*(short*)(p+2) = 4;
			*(int*)(p+4) = 1;
			*(int*)(p+8) = endfile_pos + 6 + 12 * nbflag + extra; //data start
			p+=12;

			*(short*)(p) = 277;
			*(short*)(p+2) = 3;
			*(int*)(p+4) = 1;
			*(short*)(p+8) = f_out.dims[0];
			p+=12;

			*(short*)(p) = 278;
			*(short*)(p+2) = 4;
			*(int*)(p+4) = 1;
			*(int*)(p+8) = f_out.dims[2]; //Rowperstrip
			p+=12;

			*(short*)(p) = 279;
			*(short*)(p+2) = 4;
			*(int*)(p+4) = 1;
			*(int*)(p+8) = f_out.dims[0]*f_out.dims[1]*f_out.dims[2] * sizeof(TYPE); //NBbyte
			p+=12;


			*(short*)(p) = 284;
			*(short*)(p+2) = 3;
			*(int*)(p+4) = 1;
			*(short*)(p+8) = 1; //planar_rep
			p+=12;

			*(short*)(p) = 296;
			*(short*)(p+2) = 3;
			*(int*)(p+4) = 1;
			*(short*)(p+8) = 0; //unit
			p+=12;
/*
			*(short*)(p) = 338; // extra sample meaning
			*(short*)(p+2) = 3;
			*(int*)(p+4) = 1;
			*(short*)(p+8) = 0;
			p+=12;*/



			writeDescriptor(min, max, p);

			*(int*)(p) = 0;p+= 4;

			curfp_pos = endfile_pos + 2 + 12*nbflag;

			endfile_pos += 6 + 12*nbflag + extra +  f_out.dims[0]*f_out.dims[1]*f_out.dims[2] * sizeof(TYPE);


			if (f_out.dims[0] < 3){
				*(int*)(p) = 1;
				*(int*)(p+4) = 1;
				*(int*)(p+8) = 1;
				*(int*)(p+12) = 1;
			}else{

				for(i=0;i<f_out.dims[0];i++) {*(short*)p = 8 * sizeof(TYPE); p+= 2;}
				*(int*)(p) = 1;
				*(int*)(p+4) = 1;
				*(int*)(p+8) = 1;
				*(int*)(p+2) = 1;
			}
			fwrite(buffer,sizeof(char), 6 + 12*nbflag + extra,f);


		if (isTypeEqual<TYPE,unsigned char>::ans){

			{//:

				unsigned char* buffer = new unsigned char[f_out.dims[0] * f_out.dims[1] ];
				for(coor[2] = 0; coor[2]< f_out.dims[2];coor[2]++){
					i=0;
					for(coor[1] = 0; coor[1]< f_out.dims[1];coor[1]++)
						for(coor[0] = 0; coor[0]< f_out.dims[0];coor[0]++,i++){
							buffer[i] =  (unsigned char) f_out(coor);
						}
					fwrite(buffer,sizeof(unsigned char), f_out.dims[0] * f_out.dims[1],f);
				}
				delete[](buffer);

			}//:

		} else if (isTypeEqual<TYPE,char>::ans){

			{//:

				char* buffer = new char[f_out.dims[0] * f_out.dims[1] ];
				for(coor[2] = 0; coor[2]< f_out.dims[2];coor[2]++){
					i=0;
					for(coor[1] = 0; coor[1]< f_out.dims[1];coor[1]++)
						for(coor[0] = 0; coor[0]< f_out.dims[0];coor[0]++,i++){
							buffer[i] =  (char) f_out(coor);
						}
					fwrite(buffer,sizeof(char), f_out.dims[0] * f_out.dims[1],f);
				}
				delete[](buffer);

			}//:

		} else if (isTypeEqual<TYPE,unsigned short>::ans){

			{//:

				unsigned short* buffer = new unsigned short[f_out.dims[0] * f_out.dims[1] ];
				for(coor[2] = 0; coor[2]< f_out.dims[2];coor[2]++){
					i=0;
					for(coor[1] = 0; coor[1]< f_out.dims[1];coor[1]++)
						for(coor[0] = 0; coor[0]< f_out.dims[0];coor[0]++,i++){
							buffer[i] = (unsigned short) f_out(coor);
						}
					fwrite(buffer,sizeof(unsigned short), f_out.dims[0] * f_out.dims[1],f);
				}
				delete[](buffer);

			}//:

		} else if (isTypeEqual<TYPE,short>::ans){

			{//:

				short* buffer = new short[f_out.dims[0] * f_out.dims[1] ];
				for(coor[2] = 0; coor[2]< f_out.dims[2];coor[2]++){
					i=0;
					for(coor[1] = 0; coor[1]< f_out.dims[1];coor[1]++)
						for(coor[0] = 0; coor[0]< f_out.dims[0];coor[0]++,i++){
							buffer[i] = (short) f_out(coor);
						}
					fwrite(buffer,sizeof(short), f_out.dims[0] * f_out.dims[1],f);
				}
				delete[](buffer);

			}//:

		} else if (isTypeEqual<TYPE,unsigned int>::ans){

			{//:

				unsigned int* buffer = new unsigned int[f_out.dims[0] * f_out.dims[1] ];
				for(coor[2] = 0; coor[2]< f_out.dims[2];coor[2]++){
					i=0;
					for(coor[1] = 0; coor[1]< f_out.dims[1];coor[1]++)
						for(coor[0] = 0; coor[0]< f_out.dims[0];coor[0]++,i++){
							buffer[i] = (unsigned int) f_out(coor);
						}
					fwrite(buffer,sizeof(unsigned int), f_out.dims[0] * f_out.dims[1],f);
				}
				delete[](buffer);

			}//:

		} else if (isTypeEqual<TYPE,int>::ans){

			{//:

				int* buffer = new int[f_out.dims[0] * f_out.dims[1] ];
				for(coor[2] = 0; coor[2]< f_out.dims[2];coor[2]++){
					i=0;
					for(coor[1] = 0; coor[1]< f_out.dims[1];coor[1]++)
						for(coor[0] = 0; coor[0]< f_out.dims[0];coor[0]++,i++){
							buffer[i] = (int) f_out(coor);
						}
					fwrite(buffer,sizeof(int), f_out.dims[0] * f_out.dims[1],f);
				}
				delete[](buffer);

			}//:

		} else if (sizeof(TYPE) == sizeof(float)){
				{ //:
					float* buffer = new float[f_out.dims[0] * f_out.dims[1] ];
					for(coor[2] = 0; coor[2]< f_out.dims[2];coor[2]++){
						i=0;
						for(coor[1] = 0; coor[1]< f_out.dims[1];coor[1]++)
							for(coor[0] = 0; coor[0]< f_out.dims[0];coor[0]++,i++){
								buffer[i] = (float) f_out(coor);
							}
						fwrite(buffer,sizeof(float), f_out.dims[0] * f_out.dims[1],f);
					}
					delete[](buffer);

				} //:
		}else if (isTypeEqual<TYPE,double>::ans){
				{ //:
					double* buffer = new double[f_out.dims[0] * f_out.dims[1] ];
					for(coor[2] = 0; coor[2]< f_out.dims[2];coor[2]++){
						i=0;
						for(coor[1] = 0; coor[1]< f_out.dims[1];coor[1]++)
							for(coor[0] = 0; coor[0]< f_out.dims[0];coor[0]++,i++){
								buffer[i] = (double) f_out(coor);
							}
						fwrite(buffer,sizeof(double), f_out.dims[0] * f_out.dims[1],f);
					}
					delete[](buffer);
				} //:
		}



		//	fwrite(list[imc]->data,sizeof(float), list[imc]->sizex* list[imc]->sizey* list[imc]->channels,out);
		return(true);
	}

	template <class C, class TYPE>
	bool TiffFile::put(DataGrid<C,2>& f_out, TYPE min, TYPE max,bool updateRange){

		//	printf("%i\t%i\t%i\n", f_out.dims[0],f_out.dims[1],f_out.dims[2]);
		if (endfile_pos == 0){ // nothing was written so far
			endfile_pos = curfp_pos + 4;
		}

		fseek(f,curfp_pos,SEEK_SET);
		fwrite(&endfile_pos,sizeof(unsigned int),1, f);
		fseek(f,endfile_pos,SEEK_SET);
		int i;
//		int type;

		char buffer[65536];
		char *p = buffer;
		int nbflag=18;
		int extra=0;
		Tuple<unsigned int,2> coor;
		*(short*)(p) = nbflag;
		p+=2;

		*(short*)(p) = 254;
		*(short*)(p+2) = 4;
		*(int*)(p+4) = 1;
		*(int*)(p+8) = 0;
		p+=12;

		*(short*)(p) = 256;
		*(short*)(p+2) = 4;
		*(int*)(p+4) = 1;
		*(int*)(p+8) = f_out.dims[0];
		p+=12;

		*(short*)(p) = 257;
		*(short*)(p+2) = 4;
		*(int*)(p+4) = 1;
		*(int*)(p+8) = f_out.dims[1];
		p+=12;

		*(short*)(p) = 258;
		*(short*)(p+2) = 3;
		*(int*)(p+4) = 1;
		*(short*)(p+8) = 8* sizeof(TYPE); //bitpersample
		*(short*)(p+10) = 8* sizeof(TYPE); //bitpersample



		p+=12;

		*(short*)(p) = 259;
		*(short*)(p+2) = 3;
		*(int*)(p+4) = 1;
		*(short*)(p+8) = 1; //compression
		p+=12;

		*(short*)(p) = 262;
		*(short*)(p+2) = 3;
		*(int*)(p+4) = 1;
		switch(f_out.dims[0]){//photothing
			case 1: *(short*)(p+8) = 1; break;
			case 2: *(short*)(p+8) = 1; break;
			default: *(short*)(p+8) = 2; break;
		}
		p+=12;

		*(short*)(p) = 282;
		*(short*)(p+2) = 5;
		*(int*)(p+4) = 1;
		*(int*)(p+8) = endfile_pos + 6 + 12 * nbflag + extra;
		extra += sizeof(int)*2;
		p+=12;

		*(short*)(p) = 283;
		*(short*)(p+2) = 5;
		*(int*)(p+4) = 1;
		*(int*)(p+8) = endfile_pos + 6 + 12 * nbflag + extra;
		extra += sizeof(int)*2;
		p+=12;

		*(short*)(p) = 273;
		*(short*)(p+2) = 4;
		*(int*)(p+4) = 1;
		*(int*)(p+8) = endfile_pos + 6 + 12 * nbflag + extra; //data start
		p+=12;

		*(short*)(p) = 277;
		*(short*)(p+2) = 3;
		*(int*)(p+4) = 1;
		*(short*)(p+8) = 1;
		p+=12;

		*(short*)(p) = 278;
		*(short*)(p+2) = 4;
		*(int*)(p+4) = 1;
		*(int*)(p+8) = f_out.dims[1]; //Rowperstrip
		p+=12;

		*(short*)(p) = 279;
		*(short*)(p+2) = 4;
		*(int*)(p+4) = 1;
		*(int*)(p+8) = f_out.dims[0]*f_out.dims[1] * sizeof(TYPE); //NBbyte
		p+=12;


		*(short*)(p) = 284;
		*(short*)(p+2) = 3;
		*(int*)(p+4) = 1;
		*(short*)(p+8) = 1; //planar_rep
		p+=12;

		*(short*)(p) = 296;
		*(short*)(p+2) = 3;
		*(int*)(p+4) = 1;
		*(short*)(p+8) = 0; //unit
		p+=12;

		*(short*)(p) = 338; // extra sample meaning
		*(short*)(p+2) = 3;
		*(int*)(p+4) = 1;
		*(short*)(p+8) = 0;
		p+=12;



		writeDescriptor(min, max, p);


		*(int*)(p) = 0;p+= 4;

		curfp_pos = endfile_pos + 2 + 12*nbflag;

		endfile_pos += 6 + 12*nbflag + extra +  f_out.dims[0]*f_out.dims[1] * sizeof(TYPE);


		*(int*)(p) = 1;
		*(int*)(p+4) = 1;
		*(int*)(p+8) = 1;
		*(int*)(p+12) = 1;
 
		fwrite(buffer,sizeof(char), 6 + 12*nbflag + extra,f);

		if (isTypeEqual<TYPE,unsigned char>::ans){

			{//:
				unsigned char* buffer = new unsigned char[f_out.dims[0]];
				for(coor[1] = 0; coor[1]< f_out.dims[1];coor[1]++){
					i=0;
					for(coor[0] = 0; coor[0]< f_out.dims[0];coor[0]++) buffer[i++] =  (unsigned char) f_out(coor);

					fwrite(buffer,sizeof(unsigned char), f_out.dims[0],f);
				}
				delete[](buffer);

			}//:

		} else if (isTypeEqual<TYPE,char>::ans){

			{//:

				char* buffer = new char[f_out.dims[0]];
				for(coor[1] = 0; coor[1]< f_out.dims[1];coor[1]++){
					i=0;
					for(coor[0] = 0; coor[0]< f_out.dims[0];coor[0]++) buffer[i++] =  (char) f_out(coor);

					fwrite(buffer,sizeof(char), f_out.dims[0],f);
				}
				delete[](buffer);

			}//:

		} else if (isTypeEqual<TYPE,unsigned short>::ans){

			{//:

				unsigned short* buffer = new unsigned short[f_out.dims[0] ];
				for(coor[1] = 0; coor[1]< f_out.dims[1];coor[1]++){
					i=0;
					for(coor[0] = 0; coor[0]< f_out.dims[0];coor[0]++) buffer[i++] = (unsigned short) f_out(coor);

					fwrite(buffer,sizeof(unsigned short), f_out.dims[0] ,f);
				}
				delete[](buffer);

			}//:

		} else if (isTypeEqual<TYPE,short>::ans){

			{//:

				short* buffer = new short[f_out.dims[0]];
				for(coor[1] = 0; coor[1]< f_out.dims[1];coor[1]++){
					i=0;
					for(coor[0] = 0; coor[0]< f_out.dims[0];coor[0]++) buffer[i++] =  (short) f_out(coor);

					fwrite(buffer,sizeof(short), f_out.dims[0],f);
				}
				delete[](buffer);

			}//:

		} else if (isTypeEqual<TYPE,unsigned int>::ans){

			{//:

				unsigned int* buffer = new unsigned int[f_out.dims[0] ];
				for(coor[1] = 0; coor[1]< f_out.dims[1];coor[1]++){
					i=0;
					for(coor[0] = 0; coor[0]< f_out.dims[0];coor[0]++)	buffer[i++] = (unsigned int) f_out(coor);
					fwrite(buffer,sizeof(unsigned int), f_out.dims[0] ,f);
				}
				delete[](buffer);

			}//:

		} else if (isTypeEqual<TYPE, int>::ans){

			{//:

				int* buffer = new int[f_out.dims[0] ];
				for(coor[1] = 0; coor[1]< f_out.dims[1];coor[1]++){
					i=0;
					for(coor[0] = 0; coor[0]< f_out.dims[0];coor[0]++)	buffer[i++] = (int) f_out(coor);
					fwrite(buffer,sizeof(int), f_out.dims[0] ,f);
				}
				delete[](buffer);

			}//:

		} else if (sizeof(TYPE) == sizeof(float)){
			{ //:
				float* buffer = new float[f_out.dims[0] ];
				for(coor[1] = 0; coor[1]< f_out.dims[1];coor[1]++){
					i=0;
					for(coor[0] = 0; coor[0]< f_out.dims[0];coor[0]++)	buffer[i++] = (float) f_out(coor);
					fwrite(buffer,sizeof(float), f_out.dims[0] ,f);
				}
				delete[](buffer);

			} //:
		}else if (isTypeEqual<TYPE,double>::ans){
			{ //:
				double* buffer = new double[f_out.dims[0]  ];
				for(coor[1] = 0; coor[1]< f_out.dims[1];coor[1]++){
					i=0;
					for(coor[0] = 0; coor[0]< f_out.dims[0];coor[0]++)	buffer[i++] = (double) f_out(coor);

					fwrite(buffer,sizeof(double), f_out.dims[0] ,f);
				}
				delete[](buffer);
			} //:
		}


		//	fwrite(list[imc]->data,sizeof(float), list[imc]->sizex* list[imc]->sizey* list[imc]->channels,out);
		return(true);
	}

/*
	template <class C>
	bool TiffFile::put(const DataGrid<C,3>& f_out, Vector< char* > &options){

		int nbflag=18;
		unsigned int i;
		unsigned int cur_str_len;
		for(i=0;i<options.getSize();i++){
			switch(options[i][0]){
				case 'D': //desc
					nbflag += 1;
				break;
				case 'R': // rngx/rngy
					nbflag += 2;
					// XPosition 286  RATIONAL
					// YPosition 287  RATIONAL
					// XResolution  282 RATIONAL
					// YResolution  283 RATIONAL
				break;
			}
		}

		//	printf("%i\t%i\t%i\n", f_out.dims[0],f_out.dims[1],f_out.dims[2]);
		if (endfile_pos == 0){ // nothing was written so far
			endfile_pos = curfp_pos + 4;
		}

		fseek(f,curfp_pos,SEEK_SET);
		fwrite(&endfile_pos,sizeof(unsigned int),1, f);
		fseek(f,endfile_pos,SEEK_SET);

		//		int type;

		char buffer[65536];
		char *p = buffer;

		int extra=0;

		// write mandatory flags

		Tuple<unsigned int,3> coor;
		*(short*)(p) = nbflag;
		p+=2;

		*(short*)(p) = 254;
		*(short*)(p+2) = 4;
		*(int*)(p+4) = 1;
		*(int*)(p+8) = 0;
		p+=12;

		*(short*)(p) = 256;
		*(short*)(p+2) = 4;
		*(int*)(p+4) = 1;
		*(int*)(p+8) = f_out.dims[1];
		p+=12;

		*(short*)(p) = 257;
		*(short*)(p+2) = 4;
		*(int*)(p+4) = 1;
		*(int*)(p+8) = f_out.dims[2];
		p+=12;

		*(short*)(p) = 258;
		*(short*)(p+2) = 3;
		*(int*)(p+4) = f_out.dims[0];
		if (f_out.dims[0] < 3){
			*(short*)(p+8) = 8* sizeof(C); //bitpersample
			*(short*)(p+10) = 8* sizeof(C); //bitpersample
		}else{
			*(int*)(p+8) = endfile_pos + 6 + 12 * nbflag + extra;
			extra += sizeof(short) * f_out.dims[0];
		}



		p+=12;

		*(short*)(p) = 259;
		*(short*)(p+2) = 3;
		*(int*)(p+4) = 1;
		*(short*)(p+8) = 1; //compression
		p+=12;

		*(short*)(p) = 262;
		*(short*)(p+2) = 3;
		*(int*)(p+4) = 1;
		switch(f_out.dims[0]){//photothing
			case 1: *(short*)(p+8) = 1; break;
			case 2: *(short*)(p+8) = 1; break;
			default: *(short*)(p+8) = 2; break;
		}
		p+=12;

		*(short*)(p) = 282;
		*(short*)(p+2) = 5;
		*(int*)(p+4) = 1;
		*(int*)(p+8) = endfile_pos + 6 + 12 * nbflag + extra;
		extra += sizeof(int)*2;
		p+=12;

		*(short*)(p) = 283;
		*(short*)(p+2) = 5;
		*(int*)(p+4) = 1;
		*(int*)(p+8) = endfile_pos + 6 + 12 * nbflag + extra;
		extra += sizeof(int)*2;
		p+=12;

		*(short*)(p) = 277;
		*(short*)(p+2) = 3;
		*(int*)(p+4) = 1;
		*(short*)(p+8) = f_out.dims[0];
		p+=12;

		*(short*)(p) = 278;
		*(short*)(p+2) = 4;
		*(int*)(p+4) = 1;
		*(int*)(p+8) = f_out.dims[2]; //Rowperstrip
		p+=12;

		*(short*)(p) = 279;
		*(short*)(p+2) = 4;
		*(int*)(p+4) = 1;
		*(int*)(p+8) = f_out.dims[1]*f_out.dims[2] *f_out.dims[0] * sizeof(C); //NBbyte
		p+=12;


		*(short*)(p) = 284;
		*(short*)(p+2) = 3;
		*(int*)(p+4) = 1;
		*(short*)(p+8) = 1; //planar_rep
		p+=12;

		*(short*)(p) = 296;
		*(short*)(p+2) = 3;
		*(int*)(p+4) = 1;
		*(short*)(p+8) = 0; //unit
		p+=12;

		*(short*)(p) = 338; // extra sample meaning
		*(short*)(p+2) = 3;
		*(int*)(p+4) = 1;
		*(short*)(p+8) = 0;
		p+=12;


		if (f_out.dims[0] < 3){
			*(int*)(p) = 1;
			*(int*)(p+4) = 1;
			*(int*)(p+8) = 1;
			*(int*)(p+12) = 1;
		}else{

			for(i=0;i<f_out.dims[0];i++) {*(short*)p = 8 * sizeof(C); p+= 2;}
			*(int*)(p) = 1;
			*(int*)(p+4) = 1;
			*(int*)(p+8) = 1;
			*(int*)(p+12) = 1;
		}

		fwrite(buffer,sizeof(char), 6 + 12*nbflag + extra,f);

		for(i=0;i<options.getSize();i++){
			switch(options[i][0]){
				case 'D': //desc
					*(short*)(p) = 338; // extra sample meaning
					*(short*)(p+2) = 2;
					*(int*)(p+4) = strlen(options[i])-3;
					*(short*)(p+8) = endfile_pos + 6 + 12 * nbflag + extra;
					extra += sizeof(char)*(*(int*)(p+4));
					p+=12;
					break;
				case 'R': // rngx/rngy
					// XPosition 286  RATIONAL
					// YPosition 287  RATIONAL
					*(short*)(p) = (options[i][3] == 'X') ? 286 : 287;
					*(short*)(p+2) = 12;
					*(int*)(p+4) = 1;
					*(short*)(p+8) = endfile_pos + 6 + 12 * nbflag + extra;
					extra += sizeof(double);
					p+=12;

					// XResolution  282 RATIONAL
					// YResolution  283 RATIONAL
					*(short*)(p) = (options[i][3] == 'X') ? 282 : 283;
					*(short*)(p+2) = 12;
					*(int*)(p+4) = 1;
					*(short*)(p+8) = endfile_pos + 6 + 12 * nbflag + extra;
					extra += sizeof(double);
					p+=12;
					break;
			}
		}


		*(short*)(p) = 273;
		*(short*)(p+2) = 4;
		*(int*)(p+4) = 1;
		*(int*)(p+8) = endfile_pos + 6 + 12 * nbflag + extra; //data start
		p+=12;



		C min = ExCo<C>::maximum();
		C max = ExCo<C>::minimum();
		for(coor[2] = 0; coor[2]< f_out.dims[2];coor[2]++){
			for(coor[1] = 0; coor[1]< f_out.dims[1];coor[1]++)
				for(coor[0] = 0; coor[0]< f_out.dims[0];coor[0]++,i++){
					if (ExOp::isValid(f_out(coor))){
						if (min > f_out(coor)) min = f_out(coor);
						if (max < f_out(coor)) max = f_out(coor);
					}
			}
		}
		if (min > max){ExOp::toZero(min);ExOp::toZero(max);}

		writeDescriptor(min, max, p);


		*(int*)(p) = 0;p+= 4;

		curfp_pos = endfile_pos + 2 + 12*nbflag;

		endfile_pos += 6 + 12*nbflag + extra +  f_out.dims[0]*f_out.dims[1] *f_out.dims[2] * sizeof(C);




		p+= 16;

		for(i=0;i<options.getSize();i++){
			switch(options[i][0]){
				case 'D': //desc
					cur_str_len = strlen(options[i])-3;
					memcpy(p,options[i]+4, cur_str_len);
					p += cur_str_len;
					break;
				case 'R': // rngx/rngy
					// XPosition 286  RATIONAL
					// YPosition 287  RATIONAL

					*(double*)p = ((double*)(options[i]+4))[0];
					p += sizeof(double);
					*(double*)p = ((double*)(options[i]+4))[1];
					p += sizeof(double);

					// XResolution  282 RATIONAL
					// YResolution  283 RATIONAL

					break;
			}
		}



		fwrite(buffer,sizeof(char), 6 + 12*nbflag + extra,f);

		C* loc_buffer = new C[f_out.dims[0] * f_out.dims[1] ];
		for(coor[2] = 0; coor[2]< f_out.dims[2];coor[2]++){
			i=0;
			for(coor[1] = 0; coor[1]< f_out.dims[1];coor[1]++)
				for(coor[0] = 0; coor[0]< f_out.dims[0];coor[0]++,i++){
					loc_buffer[i] = f_out(coor);
				}
			fwrite(loc_buffer,sizeof(C), f_out.dims[0] * f_out.dims[1],f);
		}
		delete[](loc_buffer);

		//	fwrite(list[imc]->data,sizeof(float), list[imc]->sizex* list[imc]->sizey* list[imc]->channels,out);
		return(true);
	}*/

#undef LFHTEMP
#define LFHTEMP template<class C, int nbdata>


	LFHTEMP
	void Forest<C,nbdata>::cluster(Vector<C> &data, double (*metric)(const C&, const C&), C (*merge)(const C&, const C&)){
		int s = data.getSize();
		this->setSize(s*2 -1);
		HeapTree< KeyElem<double, Tuple<int,2> > > pairs;

		unsigned int i;
		Tuple<unsigned int,2> coor;

		for(i=0;i<s;i++){
			(*this)[i+s-1].second = data[i];
			clear(i+s-1);
			coor[1] = i+s-1;
			for(coor[0]=s-1;coor[0]<coor[1];coor[0]++){
				pairs.insert( KeyElem<double, Tuple<int,2> >(metric((*this)[coor[0]].second,(*this)[coor[1]].second),coor));
			}
		}
		KeyElem<double, Tuple<int,2> > cur;
		i = s -2;

		while(!pairs.isempty()){
			cur = pairs.pop();
			clear(i);
			if (hasParent(cur.d[0])||hasParent(cur.d[1])) continue;
			if ((*this)[ cur.d[0] ].second < (*this)[ cur.d[1] ].second) {
				makeLeftof(cur.d[0],i);
				makeRightof(cur.d[1],i);

			}else{
				makeLeftof(cur.d[1],i);
				makeRightof(cur.d[0],i);
			}
			if (merge == NULL) (*this)[i].second = (*this)[ cur.d[0] ].second + (*this)[ cur.d[1] ].second;
			else (*this)[i].second = merge((*this)[ cur.d[0] ].second, (*this)[ cur.d[1] ].second );
			coor[0] = i;
			for(coor[1]=2*s-2;coor[0]<coor[1];coor[1]--){
				if (!hasParent(coor[1])){
					pairs.insert( KeyElem<double, Tuple<int,2> >(metric((*this)[coor[0]].second,(*this)[coor[1]].second), coor));
				}
			}
			i--;
		}


	}

	LFHTEMP
	void Forest<C,nbdata>::cluster_singlelink(Vector<C> &data, double (*metric)(const C&, const C&)){
		int s = data.getSize();
		this->setSize(s*2 -1);
		HeapTree< KeyElem<double, Tuple<int,2> > > pairs;

		unsigned int i;
		Tuple<unsigned int,2> coor;

		for(i=0;i<s;i++){
			(*this)[i+s-1].second = data[i];
			clear(i+s-1);
			coor[1] = i+s-1;
			for(coor[0]=s-1;coor[0]<coor[1];coor[0]++){
				pairs.insert( KeyElem<double, Tuple<int,2> >(metric((*this)[coor[0]].second,(*this)[coor[1]].second),coor));
			}
		}
		KeyElem<double, Tuple<int,2> > cur;
		i = s -2;
		unsigned int pa;
		unsigned int pb;





		while(!pairs.isempty()){
			cur = pairs.pop();
			clear(i);
			pa = cur.d[0];
			while(hasParent(pa)) pa = getParent(pa);
			pb = cur.d[1];
			while(hasParent(pb)) pb = getParent(pb);
			if (pa == pb) continue;
			//			if ((*this)[ pa ].second < (*this)[ pb ].second) {
			makeLeftof(pa,i);
			makeRightof(pb,i);
			//			}else{
			//				makeLeftof(pa,i);
			//				makeRightof(pb,i);
			//			}
			i--;
		}


	}

	LFHTEMP
	void Forest<C,nbdata>::cluster_completelink(Vector<C> &data, double (*metric)(const C&, const C&)){
		int s = data.getSize();
		this->setSize(s*2 -1);
		HeapTree< KeyElem<double, Tuple<int,2> > > pairs;


		unsigned int i;
		Tuple<unsigned int,2> coor;

		for(i=0;i<s;i++){
			(*this)[i+s-1].second = data[i];
			clear(i+s-1);
			coor[1] = i+s-1;
			for(coor[0]=s-1;coor[0]<coor[1];coor[0]++){
				pairs.insert( KeyElem<double, Tuple<int,2> >(metric((*this)[coor[0]].second,(*this)[coor[1]].second),coor));
			}
		}
		KeyElem<double, Tuple<int,2> > cur;

		unsigned int pa;
		unsigned int pb;
		unsigned int j;
		pair<unsigned int, unsigned int*>* counts = new pair<unsigned int, unsigned int*>[s-1];
		for(i=0;i<s-1;i++) counts[i].second = NULL;
		i = s -2;
		while(!pairs.isempty()){
			cur = pairs.pop();
			pa = cur.d[0];
			while(hasParent(pa)) pa = getParent(pa);
			pb = cur.d[1];
			while(hasParent(pb)) pb = getParent(pb);

			if (pa == pb) continue;
			if (pa > pb){
				j =pa;
				pa = pb;
				pb=j;
			}
			//		printf("hello! %i,%i\n", pa, pb);fflush(stdout);
			if (pa < s-1){
				counts[pa].second[pb-1]++;
				if (pb >= s-1){
					if (counts[pa].second[pb-1] == counts[pa].first){
						//				printf("didultramerge! %i,%i ->%i\n", pa, pb,i);fflush(stdout);
						clear(i);
						makeLeftof(pa,i);
						makeRightof(pb,i);

						counts[i].first = counts[pa].first + 1;
						counts[i].second = counts[pa].second; counts[pa].second = NULL;
						for(j=i+1;j< s-1;j++){
							if (counts[j].second != NULL){
								counts[i].second[j-1] += counts[j].second[pb-1] + ((pa < j) ? 0 :  counts[j].second[pa-1]);
							}
						}
						i--;
					}
				}else{
					if (counts[pa].second[pb-1] == counts[pa].first * counts[pb].first){
						//				printf("didstdomegamerge! %i,%i ->%i\n", pa, pb,i);fflush(stdout);
						clear(i);
						if (counts[pa].first > counts[pb].first){
							makeLeftof(pa,i);
							makeRightof(pb,i);
						}else{
							makeLeftof(pb,i);
							makeRightof(pa,i);
						}
						counts[i].first = counts[pa].first + counts[pb].first;
						for(j = pb-1;j< 2*s-2;j++){
							counts[pa].second[j] += counts[pb].second[j];
						}
						counts[i].second = counts[pa].second; counts[pa].second = NULL;
						delete[](counts[pb].second); counts[pb].second = NULL;
						for(j=i+1;j< s-1;j++){
							if (counts[j].second != NULL){
								counts[i].second[j-1] += (pa < j) ? 0 :  counts[j].second[pa-1];
								counts[i].second[j-1] += (pb < j) ? 0 :  counts[j].second[pb-1];

							}
						}
						i--;

					}
				}
			}else{
				//			printf("didstdmerge! %i,%i ->%i\n", pa, pb,i);fflush(stdout);
				clear(i);
				makeLeftof(pa,i);
				makeRightof(pb,i);
				counts[i].first =2;
				counts[i].second = new unsigned int[2*s-2];
				memset(counts[i].second,'\0',sizeof(unsigned int)*(2*s-2));
				for(j=i+1;j< s-1;j++){
					if (counts[j].second != NULL){
						counts[i].second[j-1] = counts[j].second[pa-1] + counts[j].second[pb-1];
					}
				}
				i--;
			}
		}


		delete[](counts[0].second);
		delete[](counts);

	}

	/*
	LFHTEMP template<class D>
	void Forest<C,nbdata>::cluster(Vector<D> &data, C (*metric)(const D&, const D&), D (*merge)(const D&, const D&)){
		int s = data.getSize();
		this->setSize(s*2 -1);
		HeapTree< KeyElem<double, Tuple<int,2> > > pairs;

		unsigned int i;
		Tuple<unsigned int,2> coor;

		for(i=0;i<s;i++){
			(*this)[i+s-1].second = data[i];
			clear(i+s-1);
			coor[1] = i+s-1;
			for(coor[0]=s-1;coor[0]<coor[1];coor[0]++){
				pairs.insert( KeyElem<double, Tuple<int,2> >(metric((*this)[coor[0]].second,(*this)[coor[1]].second),coor));
			}
		}
		KeyElem<double, Tuple<int,2> > cur;
		i = s -2;

		while(!pairs.isempty()){
			cur = pairs.pop();
			clear(i);
			if (hasParent(cur.d[0])||hasParent(cur.d[1])) continue;
			if ((*this)[ cur.d[0] ].second < (*this)[ cur.d[1] ].second) {
				makeLeftof(cur.d[0],i);
				makeRightof(cur.d[1],i);

			}else{
				makeLeftof(cur.d[1],i);
				makeRightof(cur.d[0],i);
			}
			if (merge == NULL) (*this)[i].second = (*this)[ cur.d[0] ].second + (*this)[ cur.d[1] ].second;
			else (*this)[i].second = merge((*this)[ cur.d[0] ].second, (*this)[ cur.d[1] ].second );
			coor[0] = i;
			for(coor[1]=2*s-2;coor[0]<coor[1];coor[1]--){
				if (!hasParent(coor[1])){
					pairs.insert( KeyElem<double, Tuple<int,2> >(metric((*this)[coor[0]].second,(*this)[coor[1]].second), coor));
				}
			}
			i--;
		}
	}*/

	LFHTEMP template<class D>
	void Forest<C,nbdata>::cluster_singlelink(Vector<D> &data, C (*metric)(const D&, const D&)){
		nbroots =1;
		int s = data.getSize();
		this->setSize(s*2 -1);
		HeapTree< KeyElem<C, Tuple<int,2> > > pairs;

		unsigned int i;
		Tuple<unsigned int,2> coor;

		for(i=0;i<s;i++){
			clear(i+s-1);
			coor[1] = i+s-1;
			for(coor[0]=s-1;coor[0]<coor[1];coor[0]++){
				pairs.insert( KeyElem<C, Tuple<int,2> >(metric( data[coor[0]-s+1],data[coor[1]-s+1]),coor));
			}
		}
		KeyElem<C, Tuple<int,2> > cur;
		i = s -2;
		unsigned int pa;
		unsigned int pb;


		while(!pairs.isempty()){
			cur = pairs.pop();
			clear(i);
			pa = cur.d[0];
			while(hasParent(pa)) pa = getParent(pa);
			pb = cur.d[1];
			while(hasParent(pb)) pb = getParent(pb);
			if (pa == pb) continue;
				(*this)[i].second = cur.first();
				makeLeftof(pa,i);
				makeRightof(pb,i);
			i--;
		}


	}

	LFHTEMP template<class D>
	void Forest<C,nbdata>::cluster_completelink(Vector<D> &data, C (*metric)(const D&, const D&)){
		nbroots =1;
		int s = data.getSize();
		this->setSize(s*2 -1);
		HeapTree< KeyElem<C, Tuple<int,2> > > pairs;


		unsigned int i;
		Tuple<unsigned int,2> coor;

		for(i=0;i<s;i++){

			clear(i+s-1);
			coor[1] = i+s-1;
			for(coor[0]=s-1;coor[0]<coor[1];coor[0]++){
				pairs.insert( KeyElem<C, Tuple<int,2> >(metric( data[coor[0]-s+1],data[coor[1]-s+1]),coor));
			}
		}
		KeyElem<C, Tuple<int,2> > cur;

		unsigned int pa;
		unsigned int pb;
		unsigned int j;
		pair<unsigned int, unsigned int*>* counts = new pair<unsigned int, unsigned int*>[s-1];
		for(i=0;i<s-1;i++) counts[i].second = NULL;
		i = s -2;
		while(!pairs.isempty()){
			cur = pairs.pop();
			pa = cur.d[0];
			while(hasParent(pa)) pa = getParent(pa);
			pb = cur.d[1];
			while(hasParent(pb)) pb = getParent(pb);

			if (pa == pb) continue;
			if (pa > pb){
				j =pa;
				pa = pb;
				pb=j;
			}
	//		printf("hello! %i,%i\n", pa, pb);fflush(stdout);
			if (pa < s-1){
				counts[pa].second[pb-1]++;
				if (pb >= s-1){
				if (counts[pa].second[pb-1] == counts[pa].first){
	//				printf("didultramerge! %i,%i ->%i\n", pa, pb,i);fflush(stdout);
					clear(i);
					(*this)[i].second = cur.k;
					makeLeftof(pa,i);
					makeRightof(pb,i);

					counts[i].first = counts[pa].first + 1;
					counts[i].second = counts[pa].second; counts[pa].second = NULL;
					for(j=i+1;j< s-1;j++){
						if (counts[j].second != NULL){
							counts[i].second[j-1] += counts[j].second[pb-1] + ((pa < j) ? 0 :  counts[j].second[pa-1]);
						}
					}
					i--;
				}
			}else{
				if (counts[pa].second[pb-1] == counts[pa].first * counts[pb].first){
	//				printf("didstdomegamerge! %i,%i ->%i\n", pa, pb,i);fflush(stdout);
					clear(i);
					(*this)[i].second = cur.k;
					if (counts[pa].first > counts[pb].first){
						makeLeftof(pa,i);
						makeRightof(pb,i);
					}else{
						makeLeftof(pb,i);
						makeRightof(pa,i);
					}
					counts[i].first = counts[pa].first + counts[pb].first;
					for(j = pb-1;j< 2*s-2;j++){
						counts[pa].second[j] += counts[pb].second[j];
					}
					counts[i].second = counts[pa].second; counts[pa].second = NULL;
					delete[](counts[pb].second); counts[pb].second = NULL;
					for(j=i+1;j< s-1;j++){
						if (counts[j].second != NULL){
							counts[i].second[j-1] += (pa < j) ? 0 :  counts[j].second[pa-1];
							counts[i].second[j-1] += (pb < j) ? 0 :  counts[j].second[pb-1];

						}
					}
					i--;

				}
			}
			}else{
	//			printf("didstdmerge! %i,%i ->%i\n", pa, pb,i);fflush(stdout);
				clear(i);
				(*this)[i].second = cur.k;
				makeLeftof(pa,i);
				makeRightof(pb,i);
				counts[i].first =2;
				counts[i].second = new unsigned int[2*s-2];
				memset(counts[i].second,'\0',sizeof(unsigned int)*(2*s-2));
				for(j=i+1;j< s-1;j++){
					if (counts[j].second != NULL){
						counts[i].second[j-1] = counts[j].second[pa-1] + counts[j].second[pb-1];
					}
				}
				i--;
			}
		}


		delete[](counts[0].second);
		delete[](counts);

	}

	LFHTEMP template<class D>
	void Forest<C,nbdata>::cluster_majoritylink(Vector<D> &data, C (*metric)(const D&, const D&)){
		nbroots =1;
		int s = data.getSize();
		this->setSize(s*2 -1);
		HeapTree< KeyElem<C, Tuple<int,2> > > pairs;


		unsigned int i;
		Tuple<unsigned int,2> coor;

		for(i=0;i<s;i++){

			clear(i+s-1);
			coor[1] = i+s-1;
			for(coor[0]=s-1;coor[0]<coor[1];coor[0]++){
				pairs.insert( KeyElem<C, Tuple<int,2> >(metric( data[coor[0]-s+1],data[coor[1]-s+1]),coor));
			}
		}
		KeyElem<C, Tuple<int,2> > cur;

		unsigned int pa;
		unsigned int pb;
		unsigned int j;
		pair<unsigned int, unsigned int*>* counts = new pair<unsigned int, unsigned int*>[s-1];
		for(i=0;i<s-1;i++) counts[i].second = NULL;
		i = s -2;
		while(!pairs.isempty()){
			cur = pairs.pop();
			pa = cur.d[0];
			while(hasParent(pa)) pa = getParent(pa);
			pb = cur.d[1];
			while(hasParent(pb)) pb = getParent(pb);

			if (pa == pb) continue;
			if (pa > pb){
				j =pa;
				pa = pb;
				pb=j;
			}
			//		printf("hello! %i,%i\n", pa, pb);fflush(stdout);
			if (pa < s-1){
				counts[pa].second[pb-1]++;
				if (pb >= s-1){
					if (counts[pa].second[pb-1] > counts[pa].first / 2){
						//				printf("didultramerge! %i,%i ->%i\n", pa, pb,i);fflush(stdout);
						clear(i);
						(*this)[i].second = cur.k;
						makeLeftof(pa,i);
						makeRightof(pb,i);

						counts[i].first = counts[pa].first + 1;
						counts[i].second = counts[pa].second; counts[pa].second = NULL;
						for(j=i+1;j< s-1;j++){
							if (counts[j].second != NULL){
								counts[i].second[j-1] += counts[j].second[pb-1] + ((pa < j) ? 0 :  counts[j].second[pa-1]);
							}
						}
						i--;
					}
				}else{
					if (counts[pa].second[pb-1] > counts[pa].first * counts[pb].first / 2){
						//				printf("didstdomegamerge! %i,%i ->%i\n", pa, pb,i);fflush(stdout);
						clear(i);
						(*this)[i].second = cur.k;
						if (counts[pa].first > counts[pb].first){
							makeLeftof(pa,i);
							makeRightof(pb,i);
						}else{
							makeLeftof(pb,i);
							makeRightof(pa,i);
						}
						counts[i].first = counts[pa].first + counts[pb].first;
						for(j = pb-1;j< 2*s-2;j++){
							counts[pa].second[j] += counts[pb].second[j];
						}
						counts[i].second = counts[pa].second; counts[pa].second = NULL;
						delete[](counts[pb].second); counts[pb].second = NULL;
						for(j=i+1;j< s-1;j++){
							if (counts[j].second != NULL){
								counts[i].second[j-1] += (pa < j) ? 0 :  counts[j].second[pa-1];
								counts[i].second[j-1] += (pb < j) ? 0 :  counts[j].second[pb-1];

							}
						}
						i--;

					}
				}
			}else{
				//			printf("didstdmerge! %i,%i ->%i\n", pa, pb,i);fflush(stdout);
				clear(i);
				(*this)[i].second = cur.k;
				makeLeftof(pa,i);
				makeRightof(pb,i);
				counts[i].first =2;
				counts[i].second = new unsigned int[2*s-2];
				memset(counts[i].second,'\0',sizeof(unsigned int)*(2*s-2));
				for(j=i+1;j< s-1;j++){
					if (counts[j].second != NULL){
						counts[i].second[j-1] = counts[j].second[pa-1] + counts[j].second[pb-1];
					}
				}
				i--;
			}
		}


		delete[](counts[0].second);
		delete[](counts);

	}

	LFHTEMP template<class D, unsigned int TSIZE>
	void Forest<C,nbdata>::cluster_bhattacharryya(const Vector< Tuple<D, TSIZE > > &data, C (*report)(const GaussElem< Tuple<D, TSIZE > > &, const GaussElem< Tuple<D, TSIZE > > &)){
        Tuple<D, TSIZE > tmp,mean; 
        unsigned int totsize = data.getSize();
		
		this->setSize(totsize*2 -1);

		HeapTree< KeyElem<C, Tuple<int,2> > > pairs;

        HeapTree< KeyElem<double, unsigned int > >* loc_merge = new HeapTree< KeyElem<double, unsigned int > >[totsize];
        GaussElem< Tuple<D, TSIZE > >* stats = new GaussElem< Tuple<D, TSIZE > >[totsize];
        HeapTree< KeyElem<double, unsigned int > > merge;
        GaussElem< Tuple<D, TSIZE > > prior;
        unsigned int* links = new unsigned int[totsize*2 -1];
        unsigned int* tom = new unsigned int[totsize];

        
		unsigned int i;
		KeyElem<double, unsigned int > toins;

        printf("Initializing Clustering: "); fflush(stdout);
        stats[0] = GaussElem< Tuple<D, TSIZE > >(data[0]);
        prior = stats[0];
        for(i=1;i<totsize;i++){
            stats[i] = GaussElem< Tuple<D, TSIZE > >(data[i]); 
            prior += stats[i];
            }

    //    prior.setWeight(0.1f);
        /*
        
        printf("fact %e\n", 0.00000000000625f * pow(prior.w, -2.0f / (4+TSIZE) )); // 16 fact => http://www.math.siu.edu/mugdadi/hellinger.pdf
        
     
        
        diagonalizer = prior.getCovariance().diagonalizer_of_inverse();

        tmp = data[0] - mean;
        tmp *= diagonalizer;
        stats[0] = GaussElem< Tuple<D, TSIZE > >(tmp);
        prior = stats[0];
        for(i=1;i<totsize;i++){
            tmp = data[i] - mean;
            tmp *= diagonalizer;
            stats[i] = GaussElem< Tuple<D, TSIZE > >(data[i]); 
            prior += stats[i];
            }
        prior.show();
        
        */
		/*
       
        mean = prior.getMean(); printf("dascale %e\n",prior_covar.log_determinant());
        double scale = exp(-0.5f * (prior_covar.log_determinant() / TSIZE));
        prior_covar *= pow(scale,2.0f); 
        prior_covar *= pow(0.25f * prior.w / (2u+TSIZE), 2.0f / (4u+TSIZE) );
		 */
		double scale = pow(0.25f * prior.w / (2u+TSIZE), 2.0f / (4u+TSIZE) );
		//	prior.setWeight(  (scale < 0.5f) ?  scale / (1.0f - scale) : 1.0f);
		Trianglix<D, TSIZE > prior_covar; ExOp::toOne(prior_covar);
	//	prior_covar *= 0.01f;
		
		for(i=totsize-1;i<totsize*2-1;i++){
			links[i] = i+ 1- totsize;
     //       stats[links[i]] = GaussElem< Tuple<D, TSIZE > >((data[links[i]] - mean) *= scale);
            stats[links[i]].setCovariance(prior_covar);
		//	stats[links[i]] += prior;
        //    printf("da det = %e  is  %e ?\n",  stats[links[i]].getCovariance().log_determinant(),prior_covar.log_determinant());fflush(stdout);
            tom[links[i]] = i;
            for(toins.d=totsize-1;toins.d<i;toins.d++){
                toins.k = stats[links[i]].dist(stats[links[toins.d]]); 
				if (!(ExOp::isValid(toins.k))) ExOp::toMax(toins.k);
                loc_merge[links[i]].insert(toins);
			}

            if (false == loc_merge[links[i]].isempty()){
                toins = loc_merge[links[i]].top();
                printf("partial best : (%i,%i) with %e (%.1f%c done!)\n",i, toins.d, toins.k,  ((50.0f*(i-totsize+1))*(i-totsize+2)) / (double)(totsize * (totsize+1)),'%');
                toins.d= i;
                merge.insert(toins);
                }
		}
       
        printf("(DONE)\nClustering: "); fflush(stdout);
        nbroots =totsize-2;
        while(!merge.isempty()){
            toins = merge.pop();
            i = toins.d;
            if (links[i] == 0xFFFFFFFF) continue;
 //printf("popinmerg %i, %i \n",i,links[i]);fflush(stdout);
            toins = loc_merge[links[i]].pop();
   //          printf("popinmerg\n");fflush(stdout);
            if (links[toins.d] == 0xFFFFFFFF){
                while(!loc_merge[links[i]].isempty()){
                    toins = loc_merge[links[i]].top();
             //       printf("popsinmerg\n");fflush(stdout);
                    if (links[toins.d] == 0xFFFFFFFF) loc_merge[links[i]].pop();
                   
                    else {toins.d= i;merge.insert(toins); break;}
 //printf("popsinmerg\n");fflush(stdout);
                }
            }else{
          //      printf("damerge!\n");
			//	printf("%f\t%f\n", stats[links[i]].w,stats[links[toins.d]].w);
                if (stats[links[i]].w > stats[links[toins.d]].w){
                	makeLeftof(i, nbroots);
                    makeRightof(toins.d, nbroots);
                }else{
                	makeLeftof(toins.d, nbroots);
                    makeRightof(i, nbroots);
                }
                printf("%i <- %i,%i with %e (%.1f%c done!)\n",nbroots,i, toins.d, toins.k,  100.0f - ((50.0f*(nbroots+1))*(nbroots)) / (double)(totsize * (totsize+1)),'%');
                (*this)[nbroots].second = report( stats[links[i]], stats[links[toins.d]]);
                stats[links[i]] += stats[links[toins.d]];
                links[nbroots] = links[i];
                loc_merge[links[toins.d]].clear();
                loc_merge[links[i]].clear();
                tom[links[toins.d]] = 0xFFFFFFFF;
                tom[links[i]] = nbroots;
                links[i] = 0xFFFFFFFF;
                links[toins.d] = 0xFFFFFFFF;


                for(i=0;i<totsize;i++){
                    if (tom[i] == 0xFFFFFFFF) continue;
                    if (tom[i] == nbroots) continue;
                    toins.d = tom[i]; // printf("%i is avail\n", tom[i]);
                    toins.k = stats[links[nbroots]].dist(stats[links[toins.d]]); 
					if (!(ExOp::isValid(toins.k))) ExOp::toMax(toins.k);
                    loc_merge[links[nbroots]].insert(toins);
                }

                if (!loc_merge[links[nbroots]].isempty()){
                    toins = loc_merge[links[nbroots]].top();
                    toins.d= nbroots;
                    merge.insert(toins);
                }
                nbroots--;
            }
        }
        
        delete[](loc_merge);
        delete[](stats);
        delete[](tom);
        delete[](links);
        nbroots =1;

        printf("(DONE)\n"); fflush(stdout);

	//	delete[](counts[0].second);
	//	delete[](counts);

	}
	
	LFHTEMP template<class D, unsigned int TSIZE>
	void Forest<C,nbdata>::cluster_bhattacharryya_simple(const Vector< Tuple<WeightElem<D,2> , TSIZE > > &data, C (*report)(const GaussElem< Tuple<D, TSIZE > > &, const GaussElem< Tuple<D, TSIZE > > &)){
        Tuple<D, TSIZE > tmp,mean; 		
        unsigned int totsize = data.getSize();
		
		this->setSize(totsize*2 -1);
		
		HeapTree< KeyElem<C, Tuple<int,2> > > pairs;
		
        HeapTree< KeyElem<double, unsigned int > >* loc_merge = new HeapTree< KeyElem<double, unsigned int > >[totsize];
        GaussElem< Tuple<D, TSIZE > >* stats = new GaussElem< Tuple<D, TSIZE > >[totsize];
        HeapTree< KeyElem<double, unsigned int > > merge;
        GaussElem< Tuple<D, TSIZE > > prior;
        unsigned int* links = new unsigned int[totsize*2 -1];
        unsigned int* tom = new unsigned int[totsize];
		
        
		unsigned int i;
		KeyElem<double, unsigned int > toins;
		
        printf("Initializing Clustering: "); fflush(stdout);
        stats[0] = GaussElem< Tuple<D, TSIZE > >(data[0]);
        prior = stats[0];
        for(i=1;i<totsize;i++){
            stats[i] = GaussElem< Tuple<D, TSIZE > >(data[i]); 
            prior += stats[i];
			
		}
		ExOp::show(prior);
		//    prior.setWeight(0.000001f);
		/*
		 mean = prior.getMean(); printf("dascale %e\n",prior_covar.log_determinant());
		 double scale = exp(-0.5f * (prior_covar.log_determinant() / TSIZE));
		 prior_covar *= pow(scale,2.0f); 
		 prior_covar *= pow(0.25f * prior.w / (2u+TSIZE), 2.0f / (4u+TSIZE) );
		 
		 double scale = pow(0.25f * prior.w / (2u+TSIZE), 2.0f / (4u+TSIZE) );
		 //	prior.setWeight(  (scale < 0.5f) ?  scale / (1.0f - scale) : 1.0f);
		 Trianglix<D, TSIZE > prior_covar; ExOp::toOne(prior_covar);*/
		//	prior_covar *= 0.01f;
		
		for(i=totsize-1;i<totsize*2-1;i++){
			links[i] = i+ 1- totsize;
			//       stats[links[i]] = GaussElem< Tuple<D, TSIZE > >((data[links[i]] - mean) *= scale);
			//     stats[links[i]].setCovariance(prior_covar);
			//	stats[links[i]] += prior;
			//    printf("da det = %e  is  %e ?\n",  stats[links[i]].getCovariance().log_determinant(),prior_covar.log_determinant());fflush(stdout);
            tom[links[i]] = i;
            for(toins.d=totsize-1;toins.d<i;toins.d++){
				
                toins.k = stats[links[i]].bhattacharryya_dist(stats[links[toins.d]]);
				//	printf("%e\t%i\n", toins.k,toins.d);
				if (!(ExOp::isValid(toins.k))) ExOp::toMax(toins.k);
                loc_merge[links[i]].insert(toins);
			}
			
            if (false == loc_merge[links[i]].isempty()){
                toins = loc_merge[links[i]].top();
                printf("partial best : (%i,%i) with %e (%.1f%c done!)\n",i, toins.d, toins.k,  ((50.0f*(i-totsize+1))*(i-totsize+2)) / (double)(totsize * (totsize+1)),'%');
                toins.d= i;
                merge.insert(toins);
			}
		}
		
        printf("(DONE)\nClustering: "); fflush(stdout);
        nbroots =totsize-2;
        while(!merge.isempty()){
            toins = merge.pop();
            i = toins.d;
            if (links[i] == 0xFFFFFFFF) continue;
			//printf("popinmerg %i, %i \n",i,links[i]);fflush(stdout);
            toins = loc_merge[links[i]].pop();
			//          printf("popinmerg\n");fflush(stdout);
            if (links[toins.d] == 0xFFFFFFFF){
                while(!loc_merge[links[i]].isempty()){
                    toins = loc_merge[links[i]].top();
					//       printf("popsinmerg\n");fflush(stdout);
                    if (links[toins.d] == 0xFFFFFFFF) loc_merge[links[i]].pop();
					
                    else {toins.d= i;merge.insert(toins); break;}
					//printf("popsinmerg\n");fflush(stdout);
                }
            }else{
				//      printf("damerge!\n");
				//	printf("%f\t%f\n", stats[links[i]].w,stats[links[toins.d]].w);
                if (stats[links[i]].w > stats[links[toins.d]].w){
                	makeLeftof(i, nbroots);
                    makeRightof(toins.d, nbroots);
                }else{
                	makeLeftof(toins.d, nbroots);
                    makeRightof(i, nbroots);
                }
                printf("%i <- %i,%i with %e (%.1f%c done!)\n",nbroots,i, toins.d, toins.k,  100.0f - ((50.0f*(nbroots+1))*(nbroots)) / (double)(totsize * (totsize+1)),'%');
                (*this)[nbroots].second = report( stats[links[i]], stats[links[toins.d]]);
				
                stats[links[i]] += stats[links[toins.d]];
                links[nbroots] = links[i];
                loc_merge[links[toins.d]].clear();
                loc_merge[links[i]].clear();
                tom[links[toins.d]] = 0xFFFFFFFF;
                tom[links[i]] = nbroots;
                links[i] = 0xFFFFFFFF;
                links[toins.d] = 0xFFFFFFFF;
				
				
                for(i=0;i<totsize;i++){
                    if (tom[i] == 0xFFFFFFFF) continue;
                    if (tom[i] == nbroots) continue;
                    toins.d = tom[i]; // printf("%i is avail\n", tom[i]);
                    toins.k = stats[links[nbroots]].bhattacharryya_dist(stats[links[toins.d]]); 
					if (!(ExOp::isValid(toins.k))) ExOp::toMax(toins.k);
                    loc_merge[links[nbroots]].insert(toins);
                }
				
                if (!loc_merge[links[nbroots]].isempty()){
                    toins = loc_merge[links[nbroots]].top();
                    toins.d= nbroots;
                    merge.insert(toins);
                }
                nbroots--;
            }
        }
        
        delete[](loc_merge);
        delete[](stats);
        delete[](tom);
        delete[](links);
        nbroots =1;
		
        printf("(DONE)\n"); fflush(stdout);
		
		//	delete[](counts[0].second);
		//	delete[](counts);
		
	}

	LFHTEMP template<class D, unsigned int TSIZE>
	void Forest<C,nbdata>::cluster_likelihood_ratio(const Vector< Tuple<WeightElem<D,2> , TSIZE > > &data, C (*report)(const GaussElem< Tuple<D, TSIZE > > &, const GaussElem< Tuple<D, TSIZE > > &)){
        Tuple<D, TSIZE > tmp,mean; 		
        unsigned int totsize = data.getSize();
		
		this->setSize(totsize*2 -1);
		
		HeapTree< KeyElem<C, Tuple<int,2> > > pairs;
		
        HeapTree< KeyElem<double, unsigned int > >* loc_merge = new HeapTree< KeyElem<double, unsigned int > >[totsize];
        GaussElem< Tuple<D, TSIZE > >* stats = new GaussElem< Tuple<D, TSIZE > >[totsize];
        HeapTree< KeyElem<double, unsigned int > > merge;
        GaussElem< Tuple<D, TSIZE > > prior;
        unsigned int* links = new unsigned int[totsize*2 -1];
        unsigned int* tom = new unsigned int[totsize];
		
        
		unsigned int i;
		KeyElem<double, unsigned int > toins;
		
        printf("Initializing Clustering: "); fflush(stdout);
        stats[0] = GaussElem< Tuple<D, TSIZE > >(data[0]);
        prior = stats[0];
        for(i=1;i<totsize;i++){
            stats[i] = GaussElem< Tuple<D, TSIZE > >(data[i]); 
            prior += stats[i];
			
		}
		ExOp::show(prior);
		//    prior.setWeight(0.000001f);
		/*
		 mean = prior.getMean(); printf("dascale %e\n",prior_covar.log_determinant());
		 double scale = exp(-0.5f * (prior_covar.log_determinant() / TSIZE));
		 prior_covar *= pow(scale,2.0f); 
		 prior_covar *= pow(0.25f * prior.w / (2u+TSIZE), 2.0f / (4u+TSIZE) );
		 
		 double scale = pow(0.25f * prior.w / (2u+TSIZE), 2.0f / (4u+TSIZE) );
		 //	prior.setWeight(  (scale < 0.5f) ?  scale / (1.0f - scale) : 1.0f);
		 Trianglix<D, TSIZE > prior_covar; ExOp::toOne(prior_covar);*/
		//	prior_covar *= 0.01f;
		
		for(i=totsize-1;i<totsize*2-1;i++){
			links[i] = i+ 1- totsize;
			//       stats[links[i]] = GaussElem< Tuple<D, TSIZE > >((data[links[i]] - mean) *= scale);
			//     stats[links[i]].setCovariance(prior_covar);
			//	stats[links[i]] += prior;
			//    printf("da det = %e  is  %e ?\n",  stats[links[i]].getCovariance().log_determinant(),prior_covar.log_determinant());fflush(stdout);
            tom[links[i]] = i;
            for(toins.d=totsize-1;toins.d<i;toins.d++){
				
                toins.k = stats[links[i]].likelihoodratio_dist(stats[links[toins.d]]);
				//	printf("%e\t%i\n", toins.k,toins.d);
				if (!(ExOp::isValid(toins.k))) ExOp::toMax(toins.k);
                loc_merge[links[i]].insert(toins);
			}
			
            if (false == loc_merge[links[i]].isempty()){
                toins = loc_merge[links[i]].top();
                printf("partial best : (%i,%i) with %e (%.1f%c done!)\n",i, toins.d, toins.k,  ((50.0f*(i-totsize+1))*(i-totsize+2)) / (double)(totsize * (totsize+1)),'%');
                toins.d= i;
                merge.insert(toins);
			}
		}
		
        printf("(DONE)\nClustering: "); fflush(stdout);
        nbroots =totsize-2;
        while(!merge.isempty()){
            toins = merge.pop();
            i = toins.d;
            if (links[i] == 0xFFFFFFFF) continue;
			//printf("popinmerg %i, %i \n",i,links[i]);fflush(stdout);
            toins = loc_merge[links[i]].pop();
			//          printf("popinmerg\n");fflush(stdout);
            if (links[toins.d] == 0xFFFFFFFF){
                while(!loc_merge[links[i]].isempty()){
                    toins = loc_merge[links[i]].top();
					//       printf("popsinmerg\n");fflush(stdout);
                    if (links[toins.d] == 0xFFFFFFFF) loc_merge[links[i]].pop();
					
                    else {toins.d= i;merge.insert(toins); break;}
					//printf("popsinmerg\n");fflush(stdout);
                }
            }else{
				//      printf("damerge!\n");
				//	printf("%f\t%f\n", stats[links[i]].w,stats[links[toins.d]].w);
                if (stats[links[i]].w > stats[links[toins.d]].w){
                	makeLeftof(i, nbroots);
                    makeRightof(toins.d, nbroots);
                }else{
                	makeLeftof(toins.d, nbroots);
                    makeRightof(i, nbroots);
                }
                printf("%i <- %i,%i with %e (%.1f%c done!)\n",nbroots,i, toins.d, toins.k,  100.0f - ((50.0f*(nbroots+1))*(nbroots)) / (double)(totsize * (totsize+1)),'%');
                (*this)[nbroots].second = toins.k; //report( stats[links[i]], stats[links[toins.d]]);
				
                stats[links[i]] += stats[links[toins.d]];
                links[nbroots] = links[i];
                loc_merge[links[toins.d]].clear();
                loc_merge[links[i]].clear();
                tom[links[toins.d]] = 0xFFFFFFFF;
                tom[links[i]] = nbroots;
                links[i] = 0xFFFFFFFF;
                links[toins.d] = 0xFFFFFFFF;
				
				
                for(i=0;i<totsize;i++){
                    if (tom[i] == 0xFFFFFFFF) continue;
                    if (tom[i] == nbroots) continue;
                    toins.d = tom[i]; // printf("%i is avail\n", tom[i]);
                    toins.k = stats[links[nbroots]].likelihoodratio_dist(stats[links[toins.d]]); 
					if (!(ExOp::isValid(toins.k))) ExOp::toMax(toins.k);
                    loc_merge[links[nbroots]].insert(toins);
                }
				
                if (!loc_merge[links[nbroots]].isempty()){
                    toins = loc_merge[links[nbroots]].top();
                    toins.d= nbroots;
                    merge.insert(toins);
                }
                nbroots--;
            }
        }
        
        delete[](loc_merge);
        delete[](stats);
        delete[](tom);
        delete[](links);
        nbroots =1;
		
        printf("(DONE)\n"); fflush(stdout);
		
		//	delete[](counts[0].second);
		//	delete[](counts);
		
	}

	LFHTEMP template<class D, unsigned int TSIZE>
	void Forest<C,nbdata>::cluster_likelihood_ratio(const Vector< GaussElem< Tuple<D, TSIZE > > > &data, C (*report)(const GaussElem< Tuple<D, TSIZE > > &, const GaussElem< Tuple<D, TSIZE > > &)){
        Tuple<D, TSIZE > tmp,mean;		
        unsigned int totsize = data.getSize();
		
		this->setSize(totsize*2 -1);
		
		HeapTree< KeyElem<C, Tuple<int,2> > > pairs;
		
        HeapTree< KeyElem<double, unsigned int > >* loc_merge = new HeapTree< KeyElem<double, unsigned int > >[totsize];
        GaussElem< Tuple<D, TSIZE > >* stats = new GaussElem< Tuple<D, TSIZE > >[totsize];
        HeapTree< KeyElem<double, unsigned int > > merge;

        unsigned int* links = new unsigned int[totsize*2 -1];
        unsigned int* tom = new unsigned int[totsize];
		
        
		unsigned int i;
		KeyElem<double, unsigned int > toins;
		
        printf("Initializing Clustering: "); fflush(stdout);
        for(i=0;i<totsize;i++) stats[i] = data[i]; 
		
		
		for(i=totsize-1;i<totsize*2-1;i++){
			links[i] = i+ 1- totsize;
            tom[links[i]] = i;
            for(toins.d=totsize-1;toins.d<i;toins.d++){
				
                toins.k = stats[links[i]].likelihoodratio_dist(stats[links[toins.d]]);
				if (!(ExOp::isValid(toins.k))) ExOp::toMax(toins.k);
                loc_merge[links[i]].insert(toins);
			}
			
            if (false == loc_merge[links[i]].isempty()){
                toins = loc_merge[links[i]].top();
                printf("partial best : (%i,%i) with %e (%.1f%c done!)\n",i, toins.d, toins.k,  ((50.0f*(i-totsize+1))*(i-totsize+2)) / (double)(totsize * (totsize+1)),'%'); if ((i & 7) ==0)fflush(stdout);
                toins.d= i;
                merge.insert(toins);
			}
		}
		
        printf("(DONE)\nClustering: "); fflush(stdout);
        nbroots =totsize-2;
        while(!merge.isempty()){
            toins = merge.pop();
            i = toins.d;
            if (links[i] == 0xFFFFFFFF) continue;
			//printf("popinmerg %i, %i \n",i,links[i]);fflush(stdout);
            toins = loc_merge[links[i]].pop();
			//          printf("popinmerg\n");fflush(stdout);
            if (links[toins.d] == 0xFFFFFFFF){
                while(!loc_merge[links[i]].isempty()){
                    toins = loc_merge[links[i]].top();
					//       printf("popsinmerg\n");fflush(stdout);
                    if (links[toins.d] == 0xFFFFFFFF) loc_merge[links[i]].pop();
					
                    else {toins.d= i;merge.insert(toins); break;}
					//printf("popsinmerg\n");fflush(stdout);
                }
            }else{
				//      printf("damerge!\n");
				//	printf("%f\t%f\n", stats[links[i]].w,stats[links[toins.d]].w);
                if (stats[links[i]].w > stats[links[toins.d]].w){
                	makeLeftof(i, nbroots);
                    makeRightof(toins.d, nbroots);
                }else{
                	makeLeftof(toins.d, nbroots);
                    makeRightof(i, nbroots);
                }
                printf("%i <- %i,%i with %e (%.1f%c done!)\n",nbroots,i, toins.d, toins.k,  100.0f - ((50.0f*(nbroots+1))*(nbroots)) / (double)(totsize * (totsize+1)),'%');if ((nbroots & 7) ==0)fflush(stdout);
                (*this)[nbroots].second = toins.k ; //report( stats[links[i]], stats[links[toins.d]]);

                stats[links[i]] += stats[links[toins.d]];
                links[nbroots] = links[i];
                loc_merge[links[toins.d]].clear();
                loc_merge[links[i]].clear();
                tom[links[toins.d]] = 0xFFFFFFFF;
                tom[links[i]] = nbroots;
                links[i] = 0xFFFFFFFF;
                links[toins.d] = 0xFFFFFFFF;
				
				
                for(i=0;i<totsize;i++){
                    if (tom[i] == 0xFFFFFFFF) continue;
                    if (tom[i] == nbroots) continue;
                    toins.d = tom[i]; // printf("%i is avail\n", tom[i]);
                    toins.k = stats[links[nbroots]].likelihoodratio_dist(stats[links[toins.d]]); 
					if (!(ExOp::isValid(toins.k))) ExOp::toMax(toins.k);
                    loc_merge[links[nbroots]].insert(toins);
                }
				
                if (!loc_merge[links[nbroots]].isempty()){
                    toins = loc_merge[links[nbroots]].top();
                    toins.d= nbroots;
                    merge.insert(toins);
                }
                nbroots--;
            }
        }
        
        delete[](loc_merge);
        delete[](stats);
        delete[](tom);
        delete[](links);
        nbroots =1;
		
        printf("(DONE)\n"); fflush(stdout);
		
		//	delete[](counts[0].second);
		//	delete[](counts);
		
	}
	
LFHTEMP template<class D>
void Forest<C,nbdata>::cluster_likelihood_ratio(const Vector< GaussElem< Tuple<D, 0u > > > &data, C (*report)(const GaussElem< Tuple<D, 0u > > &, const GaussElem< Tuple<D, 0u > > &)){
	unsigned int totsize = data.getSize();
	if (totsize == 0) return;
	this->setSize(totsize*2 -1);
	
	HeapTree< KeyElem<C, Tuple<int,2> > > pairs;
	
	HeapTree< KeyElem<double, unsigned int > >* loc_merge = new HeapTree< KeyElem<double, unsigned int > >[totsize];
	GaussElem< Tuple<D, 0u > >* stats = new GaussElem< Tuple<D, 0u > >[totsize];
	HeapTree< KeyElem<double, unsigned int > > merge;
	unsigned int* links = new unsigned int[totsize*2 -1];
	unsigned int* tom = new unsigned int[totsize];
	
	
	unsigned int i;
	KeyElem<double, unsigned int > toins;
	
	printf("Initializing Clustering: "); fflush(stdout);
	for(i=0;i<totsize;i++){
		stats[i] = data[i]; 
	}
	
	for(i=totsize-1;i<totsize*2-1;i++){
		links[i] = i+ 1- totsize;
		tom[links[i]] = i;
		for(toins.d=totsize-1;toins.d<i;toins.d++){
			toins.k = stats[links[i]].likelihoodratio_dist(stats[links[toins.d]]);
			if (!(ExOp::isValid(toins.k))) ExOp::toMax(toins.k);
			loc_merge[links[i]].insert(toins);
		}
		
		if (false == loc_merge[links[i]].isempty()){
			toins = loc_merge[links[i]].top();
			printf("partial best : (%i,%i) with %e (%.1f%c done!)\n",i, toins.d, toins.k,  ((50.0f*(i-totsize+1))*(i-totsize+2)) / (double)(totsize * (totsize+1)),'%'); if ((i & 7) ==0)fflush(stdout);
			toins.d= i;
			merge.insert(toins);
		}
	}
	
	printf("(DONE)\nClustering: "); fflush(stdout);
	nbroots =totsize-2;
	while(!merge.isempty()){
		toins = merge.pop();
		i = toins.d;
		if (links[i] == 0xFFFFFFFF) continue;
		//printf("popinmerg %i, %i \n",i,links[i]);fflush(stdout);
		toins = loc_merge[links[i]].pop();
		//          printf("popinmerg\n");fflush(stdout);
		if (links[toins.d] == 0xFFFFFFFF){
			while(!loc_merge[links[i]].isempty()){
				toins = loc_merge[links[i]].top();
				//       printf("popsinmerg\n");fflush(stdout);
				if (links[toins.d] == 0xFFFFFFFF) loc_merge[links[i]].pop();
				
				else {toins.d= i;merge.insert(toins); break;}
				//printf("popsinmerg\n");fflush(stdout);
			}
		}else{
			//      printf("damerge!\n");
			//	printf("%f\t%f\n", stats[links[i]].w,stats[links[toins.d]].w);
			if (stats[links[i]].w > stats[links[toins.d]].w){
				makeLeftof(i, nbroots);
				makeRightof(toins.d, nbroots);
			}else{
				makeLeftof(toins.d, nbroots);
				makeRightof(i, nbroots);
			}
			printf("%i <- %i,%i with %e (%.1f%c done!)\n",nbroots,i, toins.d, toins.k,  100.0f - ((50.0f*(nbroots+1))*(nbroots)) / (double)(totsize * (totsize+1)),'%');if ((nbroots & 7) ==0)fflush(stdout);
			(*this)[nbroots].second = toins.k ; //report( stats[links[i]], stats[links[toins.d]]);
			
			stats[links[i]] += stats[links[toins.d]];
			links[nbroots] = links[i];
			loc_merge[links[toins.d]].clear();
			loc_merge[links[i]].clear();
			tom[links[toins.d]] = 0xFFFFFFFF;
			tom[links[i]] = nbroots;
			links[i] = 0xFFFFFFFF;
			links[toins.d] = 0xFFFFFFFF;
			
			
			for(i=0;i<totsize;i++){
				if (tom[i] == 0xFFFFFFFF) continue;
				if (tom[i] == nbroots) continue;
				toins.d = tom[i]; // printf("%i is avail\n", tom[i]);
				toins.k = stats[links[nbroots]].likelihoodratio_dist(stats[links[toins.d]]); 
				if (!(ExOp::isValid(toins.k))) ExOp::toMax(toins.k);
				loc_merge[links[nbroots]].insert(toins);
			}
			
			if (!loc_merge[links[nbroots]].isempty()){
				toins = loc_merge[links[nbroots]].top();
				toins.d= nbroots;
				merge.insert(toins);
			}
			nbroots--;
		}
	}
	
	delete[](loc_merge);
	delete[](stats);
	delete[](tom);
	delete[](links);
	nbroots =1;
	
	printf("(DONE)\n"); fflush(stdout);
	
	//	delete[](counts[0].second);
	//	delete[](counts);
	
}

	LFHTEMP template<class D, unsigned int TSIZE>
	void Forest<C,nbdata>::cluster_likelihood(const Vector< Tuple<WeightElem<D,2> , TSIZE > > &data, C (*report)(const GaussElem< Tuple<D, TSIZE > > &, const GaussElem< Tuple<D, TSIZE > > &)){
        Tuple<D, TSIZE > tmp,mean;		
        unsigned int totsize = data.getSize();
		
		this->setSize(totsize*2 -1);
		
		HeapTree< KeyElem<C, Tuple<int,2> > > pairs;
		
        HeapTree< KeyElem<double, unsigned int > >* loc_merge = new HeapTree< KeyElem<double, unsigned int > >[totsize];
        GaussElem< Tuple<D, TSIZE > >* stats = new GaussElem< Tuple<D, TSIZE > >[totsize];
        HeapTree< KeyElem<double, unsigned int > > merge;
        GaussElem< Tuple<D, TSIZE > > prior;
        unsigned int* links = new unsigned int[totsize*2 -1];
        unsigned int* tom = new unsigned int[totsize];
		
        
		unsigned int i;
		KeyElem<double, unsigned int > toins;
		
        printf("Initializing Clustering: "); fflush(stdout);
        stats[0] = GaussElem< Tuple<D, TSIZE > >(data[0]);
        prior = stats[0];
        for(i=1;i<totsize;i++){
            stats[i] = GaussElem< Tuple<D, TSIZE > >(data[i]); 
            prior += stats[i];
			
		}
		ExOp::show(prior);
		//    prior.setWeight(0.000001f);
		 /*
		 mean = prior.getMean(); printf("dascale %e\n",prior_covar.log_determinant());
		 double scale = exp(-0.5f * (prior_covar.log_determinant() / TSIZE));
		 prior_covar *= pow(scale,2.0f); 
		 prior_covar *= pow(0.25f * prior.w / (2u+TSIZE), 2.0f / (4u+TSIZE) );
		 
		double scale = pow(0.25f * prior.w / (2u+TSIZE), 2.0f / (4u+TSIZE) );
		//	prior.setWeight(  (scale < 0.5f) ?  scale / (1.0f - scale) : 1.0f);
		Trianglix<D, TSIZE > prior_covar; ExOp::toOne(prior_covar);*/
		//	prior_covar *= 0.01f;
		
		for(i=totsize-1;i<totsize*2-1;i++){
			links[i] = i+ 1- totsize;
			//       stats[links[i]] = GaussElem< Tuple<D, TSIZE > >((data[links[i]] - mean) *= scale);
       //     stats[links[i]].setCovariance(prior_covar);
			//	stats[links[i]] += prior;
			//    printf("da det = %e  is  %e ?\n",  stats[links[i]].getCovariance().log_determinant(),prior_covar.log_determinant());fflush(stdout);
            tom[links[i]] = i;
            for(toins.d=totsize-1;toins.d<i;toins.d++){
				
                toins.k = stats[links[i]].dist(stats[links[toins.d]]);
			//	printf("%e\t%i\n", toins.k,toins.d);
				if (!(ExOp::isValid(toins.k))) ExOp::toMax(toins.k);
                loc_merge[links[i]].insert(toins);
			}
			
            if (false == loc_merge[links[i]].isempty()){
                toins = loc_merge[links[i]].top();
                printf("partial best : (%i,%i) with %e (%.1f%c done!)\n",i, toins.d, toins.k,  ((50.0f*(i-totsize+1))*(i-totsize+2)) / (double)(totsize * (totsize+1)),'%');
                toins.d= i;
                merge.insert(toins);
			}
		}
		
        printf("(DONE)\nClustering: "); fflush(stdout);
        nbroots =totsize-2;
        while(!merge.isempty()){
            toins = merge.pop();
            i = toins.d;
            if (links[i] == 0xFFFFFFFF) continue;
			//printf("popinmerg %i, %i \n",i,links[i]);fflush(stdout);
            toins = loc_merge[links[i]].pop();
			//          printf("popinmerg\n");fflush(stdout);
            if (links[toins.d] == 0xFFFFFFFF){
                while(!loc_merge[links[i]].isempty()){
                    toins = loc_merge[links[i]].top();
					//       printf("popsinmerg\n");fflush(stdout);
                    if (links[toins.d] == 0xFFFFFFFF) loc_merge[links[i]].pop();
					
                    else {toins.d= i;merge.insert(toins); break;}
					//printf("popsinmerg\n");fflush(stdout);
                }
            }else{
				//      printf("damerge!\n");
				//	printf("%f\t%f\n", stats[links[i]].w,stats[links[toins.d]].w);
                if (stats[links[i]].w > stats[links[toins.d]].w){
                	makeLeftof(i, nbroots);
                    makeRightof(toins.d, nbroots);
                }else{
                	makeLeftof(toins.d, nbroots);
                    makeRightof(i, nbroots);
                }
                printf("%i <- %i,%i with %e (%.1f%c done!)\n",nbroots,i, toins.d, toins.k,  100.0f - ((50.0f*(nbroots+1))*(nbroots)) / (double)(totsize * (totsize+1)),'%');
                (*this)[nbroots].second = toins.k; //report( stats[links[i]], stats[links[toins.d]]);
				
                stats[links[i]] += stats[links[toins.d]];
                links[nbroots] = links[i];
                loc_merge[links[toins.d]].clear();
                loc_merge[links[i]].clear();
                tom[links[toins.d]] = 0xFFFFFFFF;
                tom[links[i]] = nbroots;
                links[i] = 0xFFFFFFFF;
                links[toins.d] = 0xFFFFFFFF;
				
				
                for(i=0;i<totsize;i++){
                    if (tom[i] == 0xFFFFFFFF) continue;
                    if (tom[i] == nbroots) continue;
                    toins.d = tom[i]; // printf("%i is avail\n", tom[i]);
                    toins.k = stats[links[nbroots]].dist(stats[links[toins.d]]); 
					if (!(ExOp::isValid(toins.k))) ExOp::toMax(toins.k);
                    loc_merge[links[nbroots]].insert(toins);
                }
				
                if (!loc_merge[links[nbroots]].isempty()){
                    toins = loc_merge[links[nbroots]].top();
                    toins.d= nbroots;
                    merge.insert(toins);
                }
                nbroots--;
            }
        }
        
        delete[](loc_merge);
        delete[](stats);
        delete[](tom);
        delete[](links);
        nbroots =1;
		
        printf("(DONE)\n"); fflush(stdout);
		
		//	delete[](counts[0].second);
		//	delete[](counts);
		
	}

	LFHTEMP template<class D, unsigned int TSIZE>
	void Forest<C,nbdata>::cluster_bhattacharryya_complex(const Vector< GaussElem< Tuple<D, TSIZE > > > &data, C (*report)(const GaussElem< Tuple<D, TSIZE > > &, const GaussElem< Tuple<D, TSIZE > > &)){
        Tuple<D, TSIZE > tmp,mean;		
        unsigned int totsize = data.getSize();
		
		this->setSize(totsize*2 -1);
		
		HeapTree< KeyElem<C, Tuple<int,2> > > pairs;
		
        HeapTree< KeyElem<double, unsigned int > >* loc_merge = new HeapTree< KeyElem<double, unsigned int > >[totsize];
        GaussElem< Tuple<D, TSIZE > >* stats = new GaussElem< Tuple<D, TSIZE > >[totsize];
        HeapTree< KeyElem<double, unsigned int > > merge;
     //   GaussElem< Tuple<D, TSIZE > > prior;
        unsigned int* links = new unsigned int[totsize*2 -1];
        unsigned int* tom = new unsigned int[totsize];
		
        
		unsigned int i;
		KeyElem<double, unsigned int > toins;
		
        printf("Initializing Clustering: "); fflush(stdout);

        for(i=0;i<totsize;i++){stats[i] = data[i];}
		
		for(i=totsize-1;i<totsize*2-1;i++){
			links[i] = i+ 1- totsize;
			//       stats[links[i]] = GaussElem< Tuple<D, TSIZE > >((data[links[i]] - mean) *= scale);
			//     stats[links[i]].setCovariance(prior_covar);
			//	stats[links[i]] += prior;
			//    printf("da det = %e  is  %e ?\n",  stats[links[i]].getCovariance().log_determinant(),prior_covar.log_determinant());fflush(stdout);
            tom[links[i]] = i;
            for(toins.d=totsize-1;toins.d<i;toins.d++){
				
                toins.k = stats[links[i]].bhattacharryya_dist(stats[links[toins.d]]);
				//	printf("%e\t%i\n", toins.k,toins.d);
				if (!(ExOp::isValid(toins.k))) ExOp::toMax(toins.k);
                loc_merge[links[i]].insert(toins);
			}
			
            if (false == loc_merge[links[i]].isempty()){
                toins = loc_merge[links[i]].top();
                printf("partial best : (%i,%i) with %e (%.1f%c done!)\n",i, toins.d, toins.k,  ((50.0f*(i-totsize+1))*(i-totsize+2)) / (double)(totsize * (totsize+1)),'%');
                toins.d= i;
                merge.insert(toins);
			}
		}
		
        printf("(DONE)\nClustering: "); fflush(stdout);
        nbroots =totsize-2;
        while(!merge.isempty()){
            toins = merge.pop();
            i = toins.d;
            if (links[i] == 0xFFFFFFFF) continue;
			//printf("popinmerg %i, %i \n",i,links[i]);fflush(stdout);
            toins = loc_merge[links[i]].pop();
			//          printf("popinmerg\n");fflush(stdout);
            if (links[toins.d] == 0xFFFFFFFF){
                while(!loc_merge[links[i]].isempty()){
                    toins = loc_merge[links[i]].top();
					//       printf("popsinmerg\n");fflush(stdout);
                    if (links[toins.d] == 0xFFFFFFFF) loc_merge[links[i]].pop();
					
                    else {toins.d= i;merge.insert(toins); break;}
					//printf("popsinmerg\n");fflush(stdout);
                }
            }else{
				//      printf("damerge!\n");
				//	printf("%f\t%f\n", stats[links[i]].w,stats[links[toins.d]].w);
                if (stats[links[i]].w > stats[links[toins.d]].w){
                	makeLeftof(i, nbroots);
                    makeRightof(toins.d, nbroots);
                }else{
                	makeLeftof(toins.d, nbroots);
                    makeRightof(i, nbroots);
                }
                printf("%i <- %i,%i with %e (%.1f%c done!)\n",nbroots,i, toins.d, toins.k,  100.0f - ((50.0f*(nbroots+1))*(nbroots)) / (double)(totsize * (totsize+1)),'%');
                (*this)[nbroots].second = report( stats[links[i]], stats[links[toins.d]]);
				
                stats[links[i]] += stats[links[toins.d]];
                links[nbroots] = links[i];
                loc_merge[links[toins.d]].clear();
                loc_merge[links[i]].clear();
                tom[links[toins.d]] = 0xFFFFFFFF;
                tom[links[i]] = nbroots;
                links[i] = 0xFFFFFFFF;
                links[toins.d] = 0xFFFFFFFF;
				
				
                for(i=0;i<totsize;i++){
                    if (tom[i] == 0xFFFFFFFF) continue;
                    if (tom[i] == nbroots) continue;
                    toins.d = tom[i]; // printf("%i is avail\n", tom[i]);
                    toins.k = stats[links[nbroots]].bhattacharryya_dist(stats[links[toins.d]]); 
					if (!(ExOp::isValid(toins.k))) ExOp::toMax(toins.k);
                    loc_merge[links[nbroots]].insert(toins);
                }
				
                if (!loc_merge[links[nbroots]].isempty()){
                    toins = loc_merge[links[nbroots]].top();
                    toins.d= nbroots;
                    merge.insert(toins);
                }
                nbroots--;
            }
        }
        
        delete[](loc_merge);
        delete[](stats);
        delete[](tom);
        delete[](links);
        nbroots =1;
		
        printf("(DONE)\n"); fflush(stdout);
		
		//	delete[](counts[0].second);
		//	delete[](counts);
		
	}

	
	LFHTEMP template<class D>
	void Forest<C,nbdata>::partition_squarreroot(const Vector<D> &data, C (*metric)(const D&, const D&)){
		int s = data.getSize();
		this->setSize(s*2 -1);
		HeapTree< KeyElem<C, Tuple<int,2> > > pairs;

		// partition into groups, which are at maximum sqrt(n) in size

		int maxsize = 1+ (int) sqrt((double)s);
		unsigned int i;
		Tuple<unsigned int,2> coor;

		for(i=0;i<s;i++){

			clear(i+s-1);
			coor[1] = i+s-1;
			for(coor[0]=s-1;coor[0]<coor[1];coor[0]++){
				pairs.insert( KeyElem<C, Tuple<int,2> >(metric( data[coor[0]-s+1],data[coor[1]-s+1]),coor));
			}
		}
		KeyElem<C, Tuple<int,2> > cur;

		unsigned int pa;
		unsigned int pb;
		unsigned int j;
		pair<unsigned int, unsigned int*>* counts = new pair<unsigned int, unsigned int*>[s-1];
		for(i=0;i<s-1;i++) counts[i].second = NULL;
		i = s -2;
		while(!pairs.isempty()){
			cur = pairs.pop();
			pa = cur.d[0];
			while(hasParent(pa)) pa = getParent(pa);
			pb = cur.d[1];
			while(hasParent(pb)) pb = getParent(pb);

			if (pa == pb) continue;
			if (pa > pb){
				j =pa;
				pa = pb;
				pb=j;
			}
			//		printf("hello! %i,%i\n", pa, pb);fflush(stdout);
			if (pa < s-1){
				counts[pa].second[pb-1]++;
				if (pb >= s-1){
					if ((counts[pa].second[pb-1] > counts[pa].first / 2)&&(counts[pa].first < maxsize)){
						//				printf("didultramerge! %i,%i ->%i\n", pa, pb,i);fflush(stdout);

						clear(i);
						(*this)[i].second = cur.k;
						makeLeftof(pa,i);
						makeRightof(pb,i);

						counts[i].first = counts[pa].first + 1;
						counts[i].second = counts[pa].second; counts[pa].second = NULL;
						for(j=i+1;j< s-1;j++){
							if (counts[j].second != NULL){
								counts[i].second[j-1] += counts[j].second[pb-1] + ((pa < j) ? 0 :  counts[j].second[pa-1]);
							}
						}
						i--;
					}
				}else{
					if ((counts[pa].second[pb-1] > counts[pa].first * counts[pb].first / 2)&&(counts[pa].first + counts[pb].first <= maxsize)){
						//				printf("didstdomegamerge! %i,%i ->%i\n", pa, pb,i);fflush(stdout);
						clear(i);
						(*this)[i].second = cur.k;
						if (counts[pa].first > counts[pb].first){
							makeLeftof(pa,i);
							makeRightof(pb,i);
						}else{
							makeLeftof(pb,i);
							makeRightof(pa,i);
						}
						counts[i].first = counts[pa].first + counts[pb].first;
						for(j = pb-1;j< 2*s-2;j++){
							counts[pa].second[j] += counts[pb].second[j];
						}
						counts[i].second = counts[pa].second; counts[pa].second = NULL;
						delete[](counts[pb].second); counts[pb].second = NULL;
						for(j=i+1;j< s-1;j++){
							if (counts[j].second != NULL){
								counts[i].second[j-1] += (pa < j) ? 0 :  counts[j].second[pa-1];
								counts[i].second[j-1] += (pb < j) ? 0 :  counts[j].second[pb-1];

							}
						}
						i--;

					}
				}
			}else{
				//			printf("didstdmerge! %i,%i ->%i\n", pa, pb,i);fflush(stdout);
				clear(i);
				(*this)[i].second = cur.k;
				makeLeftof(pa,i);
				makeRightof(pb,i);
				counts[i].first =2;
				counts[i].second = new unsigned int[2*s-2];
				memset(counts[i].second,'\0',sizeof(unsigned int)*(2*s-2));
				for(j=i+1;j< s-1;j++){
					if (counts[j].second != NULL){
						counts[i].second[j-1] = counts[j].second[pa-1] + counts[j].second[pb-1];
					}
				}
				i--;
			}
		}
		nbroots=0;
		for(j=i+1;j< 2*s-1;j++){
			if (!hasParent(j)){
				moveNode(j,nbroots);
				nbroots++;
				}


		}





		delete[](counts[0].second);
		delete[](counts);

	}

	LFHTEMP template<class D, class COMP>
	void Forest<C,nbdata>::partition_squarreroot_nbcluster(const Vector<D> &data, COMP (*metric)(const D&, const D&), int nbgroups){
		int s = data.getSize();
		this->setSize(s*2 -1);
		HeapTree< KeyElem<C, Tuple<int,2> > > pairs;

		// partition into groups, which are at maximum sqrt(n) in size

		int maxsize = 1+ (int) sqrt((double)s);
		unsigned int i;
		Tuple<unsigned int,2> coor;

		for(i=0;i<s;i++){

			if (isTypeEqual<unsigned int, C>::ans){
				(*this)[i+s-1].second = i;
				}
			clear(i+s-1);
			coor[1] = i+s-1;
			for(coor[0]=s-1;coor[0]<coor[1];coor[0]++){
				pairs.insert( KeyElem<C, Tuple<int,2> >(metric( data[coor[0]-s+1],data[coor[1]-s+1]),coor));
			}
		}
		KeyElem<C, Tuple<int,2> > cur;

		unsigned int pa;
		unsigned int pb;
		unsigned int j;
		pair<unsigned int, unsigned int*>* counts = new pair<unsigned int, unsigned int*>[s-1];
		for(i=0;i<s-1;i++) counts[i].second = NULL;
		i = s -2;
		while(!pairs.isempty()){
			cur = pairs.pop();
			pa = cur.d[0];
			while(hasParent(pa)) pa = getParent(pa);
			pb = cur.d[1];
			while(hasParent(pb)) pb = getParent(pb);

			if (pa == pb) continue;
			if (pa > pb){
				j =pa;
				pa = pb;
				pb=j;
			}
			//		printf("hello! %i,%i\n", pa, pb);fflush(stdout);
			if (pa < s-1){
				counts[pa].second[pb-1]++;
				if (pb >= s-1){
					if ((counts[pa].second[pb-1] > counts[pa].first / 2)&&(counts[pa].first < maxsize)){
						//				printf("didultramerge! %i,%i ->%i\n", pa, pb,i);fflush(stdout);

						clear(i);

						if (isTypeEqual<COMP, C>::ans) (*this)[i].second = cur.k;
						makeLeftof(pa,i);
						makeRightof(pb,i);

						counts[i].first = counts[pa].first + 1;
						counts[i].second = counts[pa].second; counts[pa].second = NULL;
						for(j=i+1;j< s-1;j++){
							if (counts[j].second != NULL){
								counts[i].second[j-1] += counts[j].second[pb-1] + ((pa < j) ? 0 :  counts[j].second[pa-1]);
							}
						}
						i--;
						if (i+1 < nbgroups) break;
					}
				}else{
					if ((counts[pa].second[pb-1] > counts[pa].first * counts[pb].first / 2)&&(counts[pa].first + counts[pb].first <= maxsize)){
						//				printf("didstdomegamerge! %i,%i ->%i\n", pa, pb,i);fflush(stdout);
						clear(i);
						if (isTypeEqual<COMP, C>::ans) (*this)[i].second = cur.k;
						if (counts[pa].first > counts[pb].first){
							makeLeftof(pa,i);
							makeRightof(pb,i);
						}else{
							makeLeftof(pb,i);
							makeRightof(pa,i);
						}
						counts[i].first = counts[pa].first + counts[pb].first;
						for(j = pb-1;j< 2*s-2;j++){
							counts[pa].second[j] += counts[pb].second[j];
						}
						counts[i].second = counts[pa].second; counts[pa].second = NULL;
						delete[](counts[pb].second); counts[pb].second = NULL;
						for(j=i+1;j< s-1;j++){
							if (counts[j].second != NULL){
								counts[i].second[j-1] += (pa < j) ? 0 :  counts[j].second[pa-1];
								counts[i].second[j-1] += (pb < j) ? 0 :  counts[j].second[pb-1];

							}
						}
						i--;
						if (i+1 < nbgroups) break;
					}
				}
			}else{
				//			printf("didstdmerge! %i,%i ->%i\n", pa, pb,i);fflush(stdout);
				clear(i);

				(*this)[i].second = cur.k;
				makeLeftof(pa,i);
				makeRightof(pb,i);
				counts[i].first =2;
				counts[i].second = new unsigned int[2*s-2];
				memset(counts[i].second,'\0',sizeof(unsigned int)*(2*s-2));
				for(j=i+1;j< s-1;j++){
					if (counts[j].second != NULL){
						counts[i].second[j-1] = counts[j].second[pa-1] + counts[j].second[pb-1];
					}
				}
				i--;
			}
		}
		nbroots=0;
		for(j=i+1;j< 2*s-1;j++){
			if (!hasParent(j)){
				if (j != nbroots) moveNode(j,nbroots);
				nbroots++;
			}


		}





		delete[](counts[0].second);
		delete[](counts);

	}





	LFHTEMP
	void Forest<C,nbdata>::saveGTRfile(const char* const path, const char* const  outer,double (*metric)(const C&, const C&)) const{
		LinkAssert< (nbdata >1) > ass;
		unsigned int s = this->size();
		FILE *f = fopen(path,"w+");
		unsigned int nbinner = (s+1)>>1;
		unsigned int i;
		for(i=nbinner-2;i != 0xFFFFFFFF;i--){
			fprintf(f,"NODE%iX\t",i);
			if ((*this)[i].first[0] >= nbinner-1) fprintf(f,"%s%iX\t",  outer,(*this)[i].first[0]-nbinner+2);
			else fprintf(f,"NODE%iX\t",  (*this)[i].first[0]);
			if ((*this)[i].first[1] >= nbinner-1) fprintf(f,"%s%iX\t",  outer,(*this)[i].first[1]-nbinner+2);
			else fprintf(f,"NODE%iX\t",  (*this)[i].first[1]);
			if (metric == NULL) fprintf(f,"0.0\n",i);
			else fprintf(f,"%f\n",  exp(-0.25f * metric((*this)[(*this)[i].first[0]].second,(*this)[(*this)[i].first[1]].second)));

		}
        fclose(f);
	}
	LFHTEMP
	void Forest<C,nbdata>::saveGTRfile(const char* const path, const char* const  outer,double (*xform)(const C&))const{
		LinkAssert< (nbdata >1) > ass;
		unsigned int s = this->size();
		FILE *f = fopen(path,"w+");
		unsigned int nbinner = (s+1)>>1;
		unsigned int i;
		for(i=nbinner-2;i != 0xFFFFFFFF;i--){
			fprintf(f,"NODE%iX\t",i);
			if ((*this)[i].first[1] >= nbinner-1) fprintf(f,"%s%iX\t",  outer,(*this)[i].first[1]-nbinner+2);
			else fprintf(f,"NODE%iX\t",  (*this)[i].first[1]);
			if ((*this)[i].first[0] >= nbinner-1) fprintf(f,"%s%iX\t",  outer,(*this)[i].first[0]-nbinner+2);
			else fprintf(f,"NODE%iX\t",  (*this)[i].first[0]);
			if (xform == NULL) fprintf(f,"%f\n",(double) (*this)[i].second );
			else fprintf(f,"%f\n",  xform((*this)[i].second));

		}
        fclose(f);
	}

	LFHTEMP
    template<class DATA> void Forest<C,nbdata>::saveCDTfile(FILE* f, Vector<char*> &header , Vector<DATA> &data) const{
		LinkAssert< (nbdata >1) > ass;
		unsigned int s = this->size();
		unsigned int nbinner = (s+1)>>1;
		unsigned int i,j;

        stack<unsigned int> path;

        path.push(0u);
        while(!path.empty()){
			i = path.top(); path.pop();
			if (i < nbinner-1){
            path.push((*this)[i].first[1]);
            path.push((*this)[i].first[0]);
			}else {fprintf(f,"%s\t", header[i-nbinner+1]);ExOp::show(data[i-nbinner+1],f,1);fprintf(f,"\n");}
			
			
            }
	}
	
	LFHTEMP
    template<class DATA,  unsigned int SIZE, unsigned int ORDER> void Forest<C,nbdata>::saveCDTfile_W(FILE* f, Vector<char*> &header , Vector< Tuple<WeightElem<DATA, ORDER>, SIZE >  > &data) const{
		LinkAssert< (nbdata >1) > ass;
		unsigned int s = this->size();
		unsigned int nbinner = (s+1)>>1;
		unsigned int i,j;
		
        stack<unsigned int> path;
		
        path.push(0u);
        while(!path.empty()){
			i = path.top(); path.pop();
			if (i < nbinner-1){
				path.push((*this)[i].first[1]);
				path.push((*this)[i].first[0]);
			}else {
				fprintf(f,"%s", header[i-nbinner+1]);
				for(j=0;j<SIZE;j++) {fprintf(f,"\t");ExOp::show(data[i-nbinner+1][j].getMean(),f,1);}
				for(j=0;j<SIZE;j++) {fprintf(f,"\t"); ExOp::show(log(data[i-nbinner+1][j].getVar()),f,1);}
				fprintf(f,"\n");
			}
			
			
		}
	}
	
	LFHTEMP
    template<class DATA,  unsigned int SIZE> void Forest<C,nbdata>::saveCDTfile_W(FILE* f, Vector<char*> &header , Vector< GaussElem< Tuple<DATA, SIZE > > > &data) const{
		LinkAssert< (nbdata >1) > ass;
		unsigned int s = this->size();
		unsigned int nbinner = (s+1)>>1;
		unsigned int i,j;
		
        stack<unsigned int> path;
		
		Tuple<DATA, SIZE > tmean;
		
        path.push(0u);
        while(!path.empty()){
			i = path.top(); path.pop();
			if (i < nbinner-1){
				path.push((*this)[i].first[1]);
				path.push((*this)[i].first[0]);
			}else {
				fprintf(f,"%s", header[i-nbinner+1]);
				tmean = data[i-nbinner+1].getMean();
				for(j=0;j<SIZE;j++) {fprintf(f,"\t");ExOp::show(tmean[j],f,1);}
				tmean = data[i-nbinner+1].getVar();
				for(j=0;j<SIZE;j++) {fprintf(f,"\t"); ExOp::show(log(tmean[j]),f,1);}
				fprintf(f,"\n");
			}
			
			
		}
	}
	

	LFHTEMP
	void saveGTRfile(char* path, char* outer){
	}


	LFHTEMP
	void Forest<C,nbdata>::makeLeftof(unsigned int left, unsigned int par){
		if (nbdata & 1) (*this)[left].first[nbdata-1] = par;
		if (nbdata > 1) (*this)[par].first[0] = left;
	}

	LFHTEMP
	void Forest<C,nbdata>::makeRightof(unsigned int right, unsigned int par){
		if (nbdata & 1) (*this)[right].first[nbdata-1] = par;
		if (nbdata > 1) (*this)[par].first[1] = right;
	}

	LFHTEMP
	void Forest<C,nbdata>::clear(unsigned int node){
		int i;
		for(i=0;i<nbdata;i++) (*this)[node].first[i] = node;
	}
	LFHTEMP
	bool Forest<C,nbdata>::isclear(unsigned int node){ return (*this)[node].first[0] == node;}

	LFHTEMP
	bool Forest<C,nbdata>::hasParent(unsigned int node){
		LinkAssert<(nbdata & 1) > ass;
		return (*this)[node].first[nbdata-1] != node;
	}
	LFHTEMP
	bool Forest<C,nbdata>::hasLeft(unsigned int node){
		LinkAssert< (nbdata > 1) > ass;
		return((*this)[node].first[0] != node);
	}
	LFHTEMP
	bool Forest<C,nbdata>::hasRight(unsigned int node){
		LinkAssert< (nbdata > 1) > ass;
		return((*this)[node].first[1] != node);
	}
	LFHTEMP
	unsigned int Forest<C,nbdata>::getParent(unsigned int node)const{
		LinkAssert< (nbdata & 1) > ass;
		return((*this)[node].first[nbdata-1]);
	}
	LFHTEMP
	unsigned int Forest<C,nbdata>::getLeft(unsigned int node)const{
		LinkAssert< (nbdata > 1) > ass;
		return((*this)[node].first[0]);
	}
	LFHTEMP
	unsigned int Forest<C,nbdata>::getRight(unsigned int node)const{
		LinkAssert< (nbdata > 1) > ass;
		return((*this)[node].first[1]);

	}

	LFHTEMP
	void Forest<C,nbdata>::moveNode(unsigned int node,unsigned int tar){
		LinkAssert< (nbdata & 1)&&(nbdata > 2 ) > ass;
		clear(tar);
		unsigned int tmp = getParent(node);
		if (tmp != node) {
			if (getLeft(tmp) != node) makeRightof(tar , tmp);
			else makeLeftof(tar,tmp);
			}
		tmp = getLeft(node);
		if (tmp != node) makeLeftof(tmp , tar);
		tmp = getRight(node);
		if (tmp != node) makeRightof(tmp , tar);
		(*this)[tar].second = (*this)[node].second;
		}


#undef LFHTEMP
#define LFHTEMP template<class IDTYPE>

	LFHTEMP	SerialStore<IDTYPE>::SerialStore(const char *f_path): wasted(0){
		int i;
		if (f_path == NULL) path = NULL;
			else{
		f = tmpfile();
		if (f == NULL) {fprintf(stderr,"Cant make Temporary file %s!\n", path); exit(1);}

		FILE* g = fopen(f_path, "rb+");
		if (g != NULL) {
		KeyElem<unsigned int , unsigned int> off;
		unsigned int offset;
		IDTYPE idid;
		fread(&offset,sizeof(unsigned int),1,g);
			i = offset / (sizeof(IDTYPE) + sizeof(unsigned int)); // position of the first data chunk, tell nb elements
	//		printf("%i entries!\n", offset);
			off.k = offset;
			for(i--;i>0;i--){
				fread(&idid,sizeof(IDTYPE),1,g);
				fread(&offset,sizeof(unsigned int),1,g);
				off.d = offset - off.k;
				hash_scope[idid] = off;
				off.k = offset;
			}
			fread(&idid,sizeof(IDTYPE),1,g);
			fseek(g, 0, SEEK_END);
			size = ftell(g);
			off.d = ftell(g) - off.k;
			hash_scope[idid] = off;
			char buffer[65536];
			fseek(g, 0, SEEK_SET);
			for(i= (size>>16);i>0;i--){
				fread(buffer, sizeof(char), 65536, g);
				fwrite(buffer, sizeof(char), 65536, f);
			}
			if (size & 65535) {
				fread(buffer, sizeof(char), size & 65535, g);
				fwrite(buffer, sizeof(char), size & 65535, f);
			}
			fclose(g);
		} else size =0;
		i = strlen(f_path)+1;
		path = new char[i];
		memcpy(path,f_path,i);
		}
	}

	LFHTEMP	SerialStore<IDTYPE>::~SerialStore(){
		if (path != NULL) {
		flush();
		fclose(f);delete[](path);
		}
	}


	LFHTEMP	void SerialStore<IDTYPE>::compress(char* buffer){

		Vector< KeyElem<unsigned int, pair<IDTYPE, unsigned int> > > chlist;
		KeyElem<unsigned int, pair<IDTYPE, unsigned int> > ttt;

		typename map<IDTYPE, KeyElem<unsigned int, unsigned int> >::iterator ite2;
		for(ite2 =  hash_scope.begin(); ite2 != hash_scope.end();ite2++) {
			ttt.k = ite2->second.k;
			ttt.d.first = ite2->first;
			ttt.d.second = ite2->second.d;
			chlist.push_back(ttt);
		}
		chlist.sort();
		unsigned int off = chlist.getSize() * (sizeof(unsigned int) + sizeof(IDTYPE));
		unsigned int i=0;

		// step 1: find shifting center
		while ((i< chlist.getSize())&&(off > chlist[i].k)) off += chlist[i++].d.second;
		unsigned int j = i-1;
		unsigned int l;
		unsigned int a,b;

		// step 2:  moving up chunks!
		while(j != 0xFFFFFFFF){
			l = chlist[j].d.second;
			a = chlist[j].k+l;
			b = off;
			off -= chlist[j].d.second;
			if (a != b) {
			while(l >= 65536) {
				a -= 65536;
				b -= 65536;
				fseek(f, a,SEEK_SET);
				fread(buffer,sizeof(char),65536,f);
				fseek(f, b,SEEK_SET);
				fwrite(buffer,sizeof(char),65536,f);
				l -= 65536;
			}
			a -= l;
			b -= l;
			fseek(f, a,SEEK_SET);
			fread(buffer,sizeof(char),l,f);
			fseek(f, b,SEEK_SET);
			fwrite(buffer,sizeof(char),l,f);
			}
			chlist[j].k = b;
			j--;
		}

		// step 3:  moving down chunks!

		off = (i == 0) ? chlist.getSize() * (sizeof(unsigned int) + sizeof(IDTYPE)) : chlist[i-1].k + chlist[i-1].d.second;
		j =i;
		while(j < chlist.getSize()){
			l = chlist[j].d.second;
			a = chlist[j].k;
			b = off;
			chlist[j].k = b;
			if (a != b) {
			while(l >= 65536) {
				fseek(f, a,SEEK_SET);
				fread(buffer,sizeof(char),65536,f);
				fseek(f, b,SEEK_SET);
				fwrite(buffer,sizeof(char),65536,f);
				l -= 65536;
				a += 65536;
				b += 65536;
			}
			fseek(f, a,SEEK_SET);
			fread(buffer,sizeof(char),l,f);
			fseek(f, b,SEEK_SET);
			fwrite(buffer,sizeof(char),l,f);
			}
			off += chlist[j].d.second;
			j++;
		}

		// step 4:  write header!

		fseek(f, 0,SEEK_SET);


		KeyElem<unsigned int, unsigned int> tmpnewin;

		for(j=0;j<chlist.getSize();j++){
			fwrite(&(chlist[j].k),sizeof(unsigned int),1,f);
			fwrite(&(chlist[j].d.first),sizeof(IDTYPE),1,f);
			tmpnewin.k = chlist[j].k;
			tmpnewin.d = chlist[j].d.second;
			hash_scope[chlist[j].d.first] = tmpnewin;
		}
		j = chlist.getSize()-1;
		size = (j+1 != 0) ? chlist[j].k + chlist[j].d.second : 0;
	}

	LFHTEMP	void SerialStore<IDTYPE>::flush(char* alt_path){
		char buffer[65536];
		FILE *g = fopen(alt_path ? alt_path : path, "wb+");
		compress(buffer);
		if (g != NULL){
			fseek(f, 0,SEEK_SET);
			for(int l=(size >> 16);l>0;l--) {
				fread(buffer,sizeof(char),65536,f);
				fwrite(buffer,sizeof(char),65536,g);
			}
			if ((size &65535) != 0 ) {
				fread(buffer,sizeof(char),(size &65535),f);
				fwrite(buffer,sizeof(char),(size &65535),g);
			}
		fclose(g);
		}
	}

	LFHTEMP	bool SerialStore<IDTYPE>::has(IDTYPE whatyp) const{
		return (hash_scope.find(whatyp) != hash_scope.end());
	}

	LFHTEMP	template<class T> void SerialStore<IDTYPE>::load(IDTYPE whatyp, T & what){

		typename map<IDTYPE, KeyElem<unsigned int, unsigned int> >::iterator ite2 = hash_scope.find(whatyp);
		if (ite2 == hash_scope.end()) {
		}else{
			fseek(f,ite2->second.k,SEEK_SET);
			ExOp::load(what,f, ite2->second.d);
		}
	}

	LFHTEMP	template<class T> void SerialStore<IDTYPE>::save(IDTYPE whatyp, const T & what){
		typename map<IDTYPE, KeyElem<unsigned int, unsigned int> >::iterator ite2 = hash_scope.find(whatyp);
		KeyElem<unsigned int, unsigned int> ocoo;
		fseek(f, size, SEEK_SET);
		ocoo.k = ftell(f);
		ExOp::save(what, f);
		size = ftell(f);
		ocoo.d = size - ocoo.k;
		if (ite2 != hash_scope.end()) {
		hash_scope[whatyp] = ocoo;
		} else hash_scope[whatyp] = ocoo;
	}






	LFHTEMP	unsigned int SerialStore<IDTYPE>::itemSize(IDTYPE whatyp) const{
		typename map<IDTYPE, KeyElem<unsigned int, unsigned int> >::const_iterator ite2 = hash_scope.find(whatyp);
		return (ite2 == hash_scope.end()) ? 0: ite2->second.d;
	}

	LFHTEMP	template<class T> void SerialStore<IDTYPE>::load_arrayitem(IDTYPE whatyp, T & what,unsigned int offset, unsigned int unitsize) const{
		typename map<IDTYPE, KeyElem<unsigned int, unsigned int> >::const_iterator ite2 = hash_scope.find(whatyp);
		
        if (ite2 == hash_scope.end()) {
		printf("Did not find Entry!\n");
        }else{ printf("SEEEK %i find Entry!\n",ite2->second.k + unitsize * offset );
			fseek(f,ite2->second.k + unitsize * offset ,SEEK_SET);
			ExOp::load(what,f, sizeof(T));
		}
	}


	LFHTEMP	template<class T> void SerialStore<IDTYPE>::save_arrayitem(IDTYPE whatyp, const T & what, unsigned int offset, unsigned int unitsize){
		typename map<IDTYPE, KeyElem<unsigned int, unsigned int> >::iterator ite2 = hash_scope.find(whatyp);
		if (ite2 == hash_scope.end()) {
		}else{
		    if ((offset+1) * unitsize > ite2->second.d){
                // illegal write! crashing!
                exit(1);
		    }
			fseek(f,ite2->second.k + unitsize * offset ,SEEK_SET);
			ExOp::save(what,f);
		}
	}


	LFHTEMP void SerialStore<IDTYPE>::save_reservearray(IDTYPE whatyp, unsigned int arr_size){
		typename map<IDTYPE, KeyElem<unsigned int, unsigned int> >::iterator ite2 = hash_scope.find(whatyp);
		KeyElem<unsigned int, unsigned int> ocoo;
		fseek(f, size, SEEK_SET);
		ocoo.k = ftell(f);
        size += arr_size;
		ocoo.d = arr_size;
		if (ite2 != hash_scope.end()) {
		hash_scope[whatyp] = ocoo;
		} else hash_scope[whatyp] = ocoo;
	}

 LFHTEMP FILE* SerialStore<IDTYPE>::getItemHandle(IDTYPE whatyp) const{
		typename map<IDTYPE, KeyElem<unsigned int, unsigned int> >::const_iterator ite2 = hash_scope.find(whatyp);
		if (ite2 == hash_scope.end()) return(NULL);
		else{
            fseek(f,ite2->second.k ,SEEK_SET);
            return(f);
		}
}

	LFHTEMP	void SerialStore<IDTYPE>::show(FILE* f_f) const{

		Vector< KeyElem<unsigned int, pair<IDTYPE, unsigned int> > > chlist;

		KeyElem<unsigned int, pair<IDTYPE, unsigned int> > ttt;

		typename map< IDTYPE, KeyElem<unsigned int, unsigned int> >::const_iterator ite2 =  hash_scope.begin();
		while(ite2 != hash_scope.end()) {
			ttt.k = ite2->second.k;
			ttt.d.first = ite2->first;
			ttt.d.second = ite2->second.d;
			chlist.push_back(ttt); ite2++;
		}

		chlist.sort();


		

		for(unsigned int i=0;i<chlist.getSize();i++){
			fprintf(f_f,"(type=%i,start=%i,size=%i)\n", (unsigned int) chlist[i].d.first, chlist[i].k,chlist[i].d.second);
		}
	}

#undef LFHTEMP
#define LFHTEMP template<class C, unsigned int NBCHAN>


LFHTEMP void CharArea<C,NBCHAN>::setDefaultStyle(){
	ExOp::toMax(color_bg); ExOp::toZero(color_axe); ExOp::toZero(color_axetext);ExOp::toMax(color_min);	ExOp::toZero(color_max);
	}
LFHTEMP void CharArea<C,NBCHAN>::initialize(const Tuple<unsigned int,2> dims){
	Tuple<unsigned int, 3> coors;
	coors[0] = 3;
	coors[1] = dims[0];
	coors[2] = dims[1];
	image.setSizes(coors);
	inner_rect[0] =0;
	inner_rect[1] =0;
	inner_rect[2] =dims[0];
	inner_rect[3] =dims[1];
	ExOp::toZero(image);
}

LFHTEMP void CharArea<C,NBCHAN>::setAxes(double Xmin,double Ymin, double Xmax, double Ymax, bool Xlog, bool Ylog){
	inner_rect[0] += 14;
	inner_rect[2] -= 14;
	inner_rect[3] -= 14;
	Tuple<unsigned int,2> line[2];
	
	line[0][0] = 13;
	line[0][1] = inner_rect[3];
	line[1][0] = 13;
	line[1][1] = 0;
	Tuple<unsigned int, 3> coors;
	
	
	for(coors[2] =0; coors[2] < image.dims[2]; coors[2]++) for(coors[1] =0; coors[1] < image.dims[1]; coors[1]++) for(coors[0] =0; coors[0] < image.dims[0]; coors[0]++) image(coors) = color_bg[coors[0]];



//	printf("%i\t%i\n", line[0][0],line[1][0]);
	image.drawLine(line[0],line[1], color_axe);
	line[1][0] = inner_rect[0]+inner_rect[2] -1;
	line[1][1] = inner_rect[3];
//	printf("%i\t%i\n", line[0][0],line[1][0]);
	image.drawLine(line[0],line[1], color_axe);

	unsigned int nbbar = 3;
 	
	double step = (Xmax - Xmin) / (nbbar-1);
	double tmin = Xmin;
	double tmax = Xmax;
	for(coors[2] =0; coors[2] <2;coors[2]++){
	step = (tmax  - tmin ) / nbbar;
	
	int lstep = floor(log(step) / log(10)); //  power of 10!
	double lead = step * exp((-log(10))*lstep);
	double s = floor(lead);
		
	if (s < 3.0f){
		if (((s * (nbbar-1)) * exp(log(10)*lstep)) < (tmax -tmin) *0.75f) s += 0.5f;  
		
	}else if (s == 7.0f) s = 8.0f;
	else if (s == 9.0f) s = 10.0f;

	// lead is the maximal step. imposible to have 5 or more otherwise
	// step is allowed to change from 1 to 1 - 1/nbbar
	
	double f = ceil(tmin * exp((-log(10))*lstep));
		unsigned int i;	
	if (((tmax > 0.0f)&&(tmin < 0.0f))||(tmin == 0.0f)||(tmax == 0.0f)) { // forces to pass by 0.0
		
		for(i=0;((f+ s*i) * exp(log(10)*lstep)) < tmax ;i++) {
			if ((i==0)||(fabs(f + s * i) < fabs(lead))) lead = f + s * i;
		}
		f -= lead;
		if (((f-s) * exp(log(10)*lstep)) > tmin) f -= s;
		else if ((f * exp(log(10)*lstep)) < tmin) f += s;
		
	}
		
	
		
	

	
	char buffer[65536];
	unsigned int str_size;
	for(i=0;(f * exp(log(10)*lstep)) < tmax ;i++) {
		sprintf(buffer,"%.1e", f);
		str_size=0;
		for(coors[1]=0; buffer[coors[1]] != '\0'; coors[1]++) if (buffer[coors[1]] != 'e') str_size += 8;
	
		if (coors[2] == 0){
		line[1][0] = line[0][0] = inner_rect[0] + floor((((double)inner_rect[2]) * ((f * exp(log(10)*lstep)) - tmin)  / (tmax - tmin)));
		line[0][1] = inner_rect[3];
		line[1][1] = inner_rect[3]+ 2;
			
		image.drawLine(line[0],line[1], color_axe);
		line[0][0] -= (str_size >> 1);
		line[0][1] = image.dims[2]-1;
		if (line[0][0] - inner_rect[0] >= (inner_rect[0] << 1) +inner_rect[2]) line[0][0] = inner_rect[0];
		else if (line[0][0] + str_size >= inner_rect[0]+inner_rect[2]) line[0][0] = inner_rect[0]+inner_rect[2] - str_size;
			for(coors[0]=0;buffer[coors[0]] != '\0';coors[0]++){if (buffer[coors[0]] != 'e') {image.drawChar(buffer[coors[0]],  line[0], color_axetext); line[0][0] += 8;}}
		}else{
			line[0][0] = inner_rect[0]-1;
			line[1][0] = inner_rect[0]-3;
			
			line[1][1] = line[0][1] = floor(((double)inner_rect[3]) * (tmax - (f * exp(log(10)*lstep)))  / (tmax - tmin) );

			image.drawLine(line[0],line[1], color_axe);
			line[0][1] -= (str_size >> 1);
			line[0][0] = 0;			
			if (line[0][1] >= inner_rect[1]+inner_rect[3]) line[0][1] = 0;
			else if (line[0][1] + str_size >= inner_rect[1]+inner_rect[3]) line[0][1] = inner_rect[1]+inner_rect[3]-1 - str_size;
			for(coors[0]=0;buffer[coors[0]] != '\0';coors[0]++){if (buffer[coors[0]] != 'e') {image.drawDownChar(buffer[coors[0]],  line[0], color_axetext); line[0][1] += 8;}}		
		}
		f+= s;
	}
		tmin = Ymin;
		tmax = Ymax;
	}
	
}

LFHTEMP template<class D, class E> void CharArea<C,NBCHAN>::drawFrame(const DataGrid<D, 2> &fr,const E &min,const E &max){
	Tuple<unsigned int, 3> coors;
	Tuple<unsigned int, 2> coor;
	Tuple<C, NBCHAN> dev = color_max - color_min;
	for(coor[1] =0; coor[1] < fr.dims[1]; coor[1]++){ coors[2] = inner_rect[3] - coor[1] - 1;
		for(coor[0] =0; coor[0] < fr.dims[0]; coor[0]++){ coors[1] = coor[0] + inner_rect[0];
			if (fr(coor) <= min) {for(coors[0]=0;coors[0]<NBCHAN;coors[0]++) image(coors) = color_min[coors[0]];}
			else if (fr(coor) >= max) {for(coors[0]=0;coors[0]<NBCHAN;coors[0]++) image(coors) = color_max[coors[0]];}
			else {for(coors[0]=0;coors[0]<NBCHAN;coors[0]++) image(coors) = color_min[coors[0]] + (((fr(coor) - min) * dev[coors[0]]) / (max - min));}
		}
	}
	
}

LFHTEMP void CharArea<C,NBCHAN>::overwriteFrame(const DataGrid<C, 3> &fr){
	Tuple<unsigned int, 3> coors;
	Tuple<unsigned int, 3> coor;
	Tuple<C, NBCHAN> dev = color_max - color_min;
	for(coor[2] =0; coor[2] < fr.dims[2]; coor[2]++){ coors[2] = inner_rect[3] - coor[2] - 1;
		for(coor[1] =0; coor[1] < fr.dims[1]; coor[1]++){ coors[1] = coor[1] + inner_rect[0];
			for(coors[0]=0,coor[0]=0;coors[0]<NBCHAN;coors[0]++,coor[0]++) image(coors) = fr(coor);
		}
	}
	
}



