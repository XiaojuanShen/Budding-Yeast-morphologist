/*
 * Vector.hpp
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

#undef LFHTEMP
#define LFHTEMP template <class C>
	
	
	LFHTEMP Vector<C>::Vector(void* owner): asize(0), darray((C*)owner){
		//printf("vec assign! %i\n", (int)owner);
	}
	
	LFHTEMP Vector<C>::Vector(C* _data, int _size){
		unsigned int i;
		setSize(_size);
		for(i=0;i<_size;i++) darray[i] = _data[i];
	}
	
	
	
	LFHTEMP void Vector<C>::show(FILE* f_out, int l) const{
		unsigned int s = getSize();
		switch(l){
			case 1:
			case 0:
				for(unsigned int i=0;i<s;i++) {
					if (i != 0) fprintf(f_out,"\t");
					ExOp::show(darray[i],f_out,2);
				}
				if (l==0) fprintf(f_out,"\n");
				break;
			case 2:
				fprintf(f_out,"[");
				for(unsigned int i=0;i<s;i++) {
					if (i != 0) fprintf(f_out,";");
					ExOp::show(darray[i],f_out,3);
				}
				fprintf(f_out,"]");
				break;
			default:
				fprintf(f_out,"(");
				for(unsigned int i=0;i<s;i++) {
					if (i != 0) fprintf(f_out,",");
					ExOp::show(darray[i],f_out,4);
				}
				fprintf(f_out,")");
				break;
		}
	}
	
	
	
	LFHTEMP void Vector<C>::save(FILE* f) const {
		unsigned int s = getSize();
		fwrite(&s,sizeof(unsigned int), 1, f);
		if (ExCo<C>::IsPOD){
			fwrite(darray,sizeof(C),getSize(),f);
		}else{
			unsigned int s = getSize(); for(unsigned int i=0;i<s;i++) ExOp::save(darray[i],f);
		}
	}
	LFHTEMP void Vector<C>::load(FILE* f, unsigned int ch_size) {
		unsigned int s;
		fread(&s,sizeof(unsigned int), 1, f);
		this->setSize(s);
		
		if (ExCo<C>::IsPOD){
			//s = ch_size / sizeof(C);
			//this->setSize(s);
			fread(darray,sizeof(C),s,f);
		}else{
			//s = ((unsigned int)ftell(f)) + ch_size;
			for(unsigned int i=0;i<s;i++) ExOp::load(darray[i],f);
		}
	}
	
	
	
	
	LFHTEMP template <class OC>
	Vector<C>::Vector(const Vector<OC> & v): asize(0), darray(NULL){
		(*this) = v;
	}
	
	LFHTEMP template<class O, unsigned int Osize> Vector<C>::Vector(const Tuple<O,Osize> &from) : darray(NULL){
		//	printf("tuble?vec assign!\n");
		setSize(Osize);
		for(unsigned int i=0;i<Osize;i++) darray[i] = (C)from[i];
	}
	LFHTEMP Vector<C>::~Vector(){
		
		if (asize != 0) delete[](darray);
		
	}
	
	
	LFHTEMP template<class D> Vector<C>::Vector(Vector<D> const & other,void* owner) : asize(0){
		setSize(other.getSize());
		for(unsigned int i=0;i<(unsigned int)(asize & 0x7FFFFFFF);i++) darray[i] = (C) other.darray[i];
	//	if (ExCo<C>::NeedsAddLink) LinkMem.setOwner(owner, darray, darray + sizeof(C) * (*this).allocatedSize() )(darray);
	}
	
	/*LFHTEMP void Vector<C>::setLinkMemory(void* new_mem){
		LinkAssert< ExCo<C>::NeedsAddLink == false > ass;
		
		if (asize == 0) darray = new_mem;
		else{
			LinkMem.setOwner(new_mem, darray, darray + sizeof(C) * (*this).allocatedSize() )(darray);
		}
		
	}*/
	
	
	LFHTEMP unsigned int Vector<C>::allocatedSize() const{
		
		unsigned int tmp = (asize-1) & 0x7FFFFFFF;
		if (tmp == 0) return(0);
		else tmp--;
		tmp |= tmp >> 8;
		tmp |= tmp >> 4;
		tmp |= tmp >> 2;
		tmp |= tmp >> 1;
		tmp++;
		return  (asize & 0x80000000) ? (tmp << 1) : tmp;
	}
	
	
	
	LFHTEMP C&  Vector<C>::push_back(){
		int rsize = (asize & 0x7FFFFFFF);
        C* swp;
        int i;
		if (((asize + 1)^(asize))> rsize){
			if (asize & 0x80000000) asize &= 0x7FFFFFFf;
			else{
					// needs up alloc
					swp = new C[((rsize+1) << 1)];
                    if (swp == NULL) {fprintf(stderr,"could not allocate %i bytes of memory!\n", ((rsize+1) << 1) * sizeof(C) ); exit(1);}
                   // printf("triedand got %i!\n", swp ); fflush(stdout);
					
					if (rsize> 0){
						for(i=0;i<rsize;i++){
							ExOp::memmove(swp[i], darray[i]);
                            //swp[i]= darray[i];
						}
						delete[](darray);
					}
				//	if (ExCo<C>::NeedsAddLink){
				//		if (rsize == 0) LinkMem.addOwner(darray,swp,swp +1);
				//		else LinkMem.moveOwner(darray, swp, swp+ ((rsize+1) << 1));
				//	}
					darray = swp;

			}
		}
        return(darray[((asize++) & 0x7FFFFFFF)]);
	}
	
	LFHTEMP void Vector<C>::pop_back(){
		int rsize = (asize & 0x7FFFFFFF);
		if (((asize)^(asize-1))> rsize-1){
			if (rsize-1 == 0){
				delete[](darray);
				//if (ExCo<C>::NeedsAddLink){
				//	darray =(C*) LinkMem.removeOwner(darray);
				//}else 
                darray = NULL;
				asize =0;
				return;
			}else{
				if ((asize & 0x80000000) == 0) asize |= 0x80000000;
				else{
					// needs down alloc
					{//:
						C* swp = new C[(rsize) << 1];
						int i;
						for(i=0;i<rsize;i++){
							ExOp::memmove(swp[i], darray[i]);
						}
						delete[](darray);
						//if (ExCo<C>::NeedsAddLink) LinkMem.moveOwner(darray, swp, swp+ ((rsize) << 1));
						darray = swp;
					}//:
				}
			}
		}
		asize--;
	}
	
	LFHTEMP void Vector<C>::pop_swap(unsigned int w){
		if (w+1 != (asize & 0x7FFFFFFF) ) ExOp::memmove(darray[w], darray[(asize & 0x7FFFFFFF)-1]);
		int rsize = (asize & 0x7FFFFFFF);
		if (((asize)^(asize-1))> rsize-1){
			if (rsize-1 == 0){
				delete[](darray);
				// if (ExCo<C>::NeedsAddLink) darray =(C*) LinkMem.removeOwner(darray); else 
                darray = NULL;
				asize =0;
				return;
			}else{
				if ((asize & 0x80000000) == 0) asize |= 0x80000000;
				else{
					// needs down alloc
					{//:
						C* swp = new C[(rsize) << 1];
						int i;
						for(i=0;i<rsize;i++){
							ExOp::memmove(swp[i], darray[i]);
						}
						delete[](darray);
						//if (ExCo<C>::NeedsAddLink) LinkMem.moveOwner(darray, swp, swp+ ((rsize) << 1));
						darray = swp;
					}//:
				}
				
				
			}
		}
		asize--;
	}
	
	LFHTEMP void Vector<C>::setSize(unsigned int nsize) {
		if (nsize == (unsigned int)(asize &0x7FFFFFFF)) return;
		if (asize != 0) delete[](darray);
        asize = (int)nsize;
        if (nsize == 0) return;
		unsigned int i =1;
		while(i< nsize) i = (i<< 1);
		darray = new C[i];
//printf("allocating (%x) (owner=%x, size=%i)\n",(int)darray, (int)this,i);fflush(stdout);
		//	printf("allocating (%x) (owner=%x)\n",(int)darray, (int)this);fflush(stdout);
		
	}

	LFHTEMP void Vector<C>::DownSize(unsigned int nsize) {
		if (nsize == (unsigned int)(asize &0x7FFFFFFF)) return;

        if (nsize == 0) return;
		unsigned int i =1;
		while(i< nsize) i = (i<< 1);
        C* tswap;
        if (i < (asize & 0x7FFFFFFF)){ 
            tswap = new C[i];
            asize = (int)nsize;
            for(i=0;i< nsize;i++) ExOp::memmove(tswap[i], darray[i]);
            delete[](darray);
            darray = tswap;
        }else{ // no resize
            if ((asize & 0x7FFFFFFF) < nsize) {fprintf(stderr,"Illegal Vector DownSize!"); exit(1);}
        }

        //printf("allocating (%x) (owner=%x, size=%i)\n",(int)darray, (int)this,i);fflush(stdout);
		//	printf("allocating (%x) (owner=%x)\n",(int)darray, (int)this);fflush(stdout);
		
	}

	
	LFHTEMP unsigned int Vector<C>::size() const {return (unsigned int)(asize & 0x7FFFFFFF);}
	LFHTEMP unsigned int Vector<C>::getSize() const {return (unsigned int)(asize & 0x7FFFFFFF);}
	
	LFHTEMP C* Vector<C>::begin() const{
		return(darray);
	}
	LFHTEMP C* Vector<C>::end() const{
		return(darray + (asize & 0x7FFFFFFF));
	}
	LFHTEMP C* Vector<C>::last() const{
		return(darray + (asize & 0x7FFFFFFF) -1);
	}
	
	LFHTEMP C& Vector<C>::operator[](int const which){
		return(darray[which]);
	}
	
	
	LFHTEMP
	void Vector<C>::clear(){
		if (asize != 0) delete[](darray);
		asize =0;
	}
	
	//LFHTEMP
	//	template<int flag2>	Vector<C>::operator PolyThing<C2> () const{
	//		return(PolyThing<C>::interpolate(darray, asize & 0x7FFFFFFF) );
	//	}
	
	
LFHTEMP	Vector<C>& Vector<C>::operator=(const Vector<C> & v){
		setSize(v.getSize());
		for(unsigned int i=0;i< (unsigned int)(asize &0x7FFFFFFF) ;i++) darray[i] = v.darray[i];
		return(*this);
	}
	
LFHTEMP Vector<C>& Vector<C>::memmove(Vector<C> & v){
		asize = v.asize;
		v.asize = 0;
		darray = v.darray;
		return(*this);
	}
    
	LFHTEMP
	template<class OC>
	Vector<C>& Vector<C>::operator=(const Vector<OC> & v){
		//	printf("vector assignement\n");
		setSize(v.getSize());
		for(unsigned int i=0;i<(unsigned int)(asize & 0x7FFFFFFF);i++) darray[i] = (C) v.darray[i];
		return(*this);
	}
	
	LFHTEMP
	template<class D> Vector<C>& Vector<C>::operator+=(const Vector<D> & v){
		for(unsigned int i=0;i<(unsigned int)(asize & 0x7FFFFFFF);i++) darray[i] += v[i];
		return(*this);
	}
	
	LFHTEMP
	template<class D> Vector<C>& Vector<C>::operator+=(const D & v){
		for(unsigned int i=0;i<(unsigned int)(asize & 0x7FFFFFFF);i++) darray[i] += v;
		return(*this);
	}
	LFHTEMP
	template<class D> Vector<C>& Vector<C>::operator-=(const Vector<D> & v){
		for(unsigned int i=0;i<(unsigned int)(asize & 0x7FFFFFFF);i++) darray[i] -= v[i];
		return(*this);
	}
	LFHTEMP
	template<class D> Vector<C>& Vector<C>::operator-=(const D & v){
		for(unsigned int i=0;i<(unsigned int)(asize & 0x7FFFFFFF);i++) darray[i] -= v;
		return(*this);
	}
	LFHTEMP
	template<class D> Vector<C>& Vector<C>::operator*=(const Vector<D> & v){
		for(unsigned int i=0;i<(unsigned int)(asize & 0x7FFFFFFF);i++) darray[i] *= v[i];
		return(*this);
	}
	LFHTEMP
	template<class D> Vector<C>& Vector<C>::operator*=(const D & v){
		for(unsigned int i=0;i<(unsigned int)(asize & 0x7FFFFFFF);i++) darray[i] *= v;
		return(*this);
	}
	LFHTEMP
	template<class D> Vector<C>& Vector<C>::operator/=(const Vector<D> & v){
		for(unsigned int i=0;i<(unsigned int)(asize & 0x7FFFFFFF);i++) darray[i] /= v[i];
		return(*this);
	}
	LFHTEMP
	template<class D> Vector<C>& Vector<C>::operator/=(const D & v){
		for(unsigned int i=0;i<(unsigned int)(asize & 0x7FFFFFFF);i++) darray[i] /= v;
		return(*this);
	}
	
		LFHTEMP void Vector<C>::sort(){
		C swap;
		int i;
		int j;
		stack<unsigned int> maxi;
		if ( (asize& 0x7FFFFFFF) == 0) return;
		maxi.push(asize & 0x7FFFFFFF);
		int mini = 0;
		while(!(maxi.empty())){
			//		printf("doing (%i, %i)!\n", mini, maxi.top()-1);
			switch(maxi.top() - mini){
				case 0:
				case 1:
					mini = maxi.top()+1; maxi.pop();
					break;
				case 2:
					if (ExOp::isGT(darray[mini], darray[mini+1])){
						ExOp::memmove(swap, darray[mini+1]);
						ExOp::memmove(darray[mini+1], darray[mini]);
						ExOp::memmove(darray[mini], swap);
					}
					mini += 3;  maxi.pop();
					break;
				case 3:
					if (ExOp::isLT(darray[mini], darray[mini+2])){
						if (darray[mini+1] < darray[mini]){
							ExOp::memmove(swap,darray[mini+1]);
							ExOp::memmove(darray[mini+1], darray[mini]);
							ExOp::memmove(darray[mini], swap);
						}else if (ExOp::isLT(darray[mini+2] , darray[mini+1])){
							ExOp::memmove(swap,darray[mini+2]);
							ExOp::memmove(darray[mini+2],darray[mini+1]);
							ExOp::memmove(darray[mini+1], swap);
						}
					}else{
						ExOp::memmove(swap, darray[mini+2]);
						if (ExOp::isLT(darray[mini+1] , darray[mini+2])){
							ExOp::memmove(darray[mini+2], darray[mini]);
							ExOp::memmove(darray[mini], darray[mini+1]);
							ExOp::memmove(darray[mini+1], swap);
						}else{
							if (ExOp::isLT(darray[mini] , darray[mini+1])){
								ExOp::memmove(darray[mini+2], darray[mini+1]);
								ExOp::memmove(darray[mini+1], darray[mini]);
							}else{
								ExOp::memmove(darray[mini+2], darray[mini]);
								
							}
							ExOp::memmove(darray[mini], swap);
						}
					}
					mini += 4;  maxi.pop();
					break;
				default:
					j = maxi.top()-1;
					i = (mini + j) >> 1; // printf("%i\t%i\n", i, j);
					switch( ((ExOp::isLT(darray[mini], darray[i])) ? 1 : 0) | ((darray[j-1] < darray[j]) ? 2 : 0) ){
						case 0:
							if (ExOp::isLT(darray[mini] , darray[j-1])){
								// daray[i] < darray[mini] < darray[j-1]  ?  darray[j]
								ExOp::memmove(swap, darray[mini]);
								ExOp::memmove(darray[mini], darray[i]);
								ExOp::memmove(darray[i], darray[j]);
								ExOp::memmove(darray[j], darray[j-1]);
							}else{
								// daray[j] < darray[j-1] < darray[mini]  ?  darray[i]
								ExOp::memmove(swap, darray[mini]);
								ExOp::memmove(darray[mini], darray[j]);
								ExOp::memmove(darray[j], swap);
								ExOp::memmove(swap, darray[j-1]);
							}
							break;
						case 1:
							
							if (ExOp::isLT(darray[i] ,darray[j-1])){
								// daray[mini] < darray[i] < darray[j-1]  ?  darray[j]
								ExOp::memmove(swap, darray[i]);
								ExOp::memmove(darray[i], darray[j]);
								ExOp::memmove(darray[j], darray[j-1]);
							}else{
								// daray[j] < darray[j-1] < darray[i]  ?  darray[mini]
								ExOp::memmove(swap, darray[i]);
								ExOp::memmove(darray[i], darray[mini]);
								ExOp::memmove(darray[mini], darray[j]);
								ExOp::memmove(darray[j], swap);
								ExOp::memmove(swap, darray[j-1]);
							}
							break;
						case 2:
							if (ExOp::isLT(darray[mini] , darray[j])){
								// daray[i] < darray[mini] < darray[j]  ?  darray[j-1]
								ExOp::memmove(swap, darray[mini]);
								ExOp::memmove(darray[mini], darray[i]);
								ExOp::memmove(darray[i], darray[j-1]);
							}else{
								// daray[j-1] < darray[j] < darray[mini]  ?  darray[i]
								ExOp::memmove(swap, darray[j]);
								ExOp::memmove(darray[j], darray[mini]);
								ExOp::memmove(darray[mini], darray[j-1]);
							}
							break;
						case 3:
							if (ExOp::isLT(darray[i] , darray[j])){
								// daray[mini] < darray[i] < darray[j]  ?  darray[j-1]
								ExOp::memmove(swap, darray[i]);
								ExOp::memmove(darray[i], darray[j-1]);
							}else{
								// daray[j-1] < darray[j] < darray[i]  ?  darray[mini]
								swap = darray[j];
								ExOp::memmove(darray[j], darray[i]);
								ExOp::memmove(darray[i], darray[mini]);
								ExOp::memmove(darray[mini], darray[j-1]);
							}
							break;
					}
					//		printf("pivot! %i!\n", swap);
					i = mini+1;
					j--;
					//		printf("pospos! %i,%i!\n", i,j);
					while(true){
						while(ExOp::isLE(darray[i] , swap)) {
							if ((++i) == j) break;
							if (darray[i] >= swap) break;
							if ((++i) == j) break;
						}
						//			printf("%i,%i mint\n", i,j);
						if (i >= j) break;
						ExOp::memmove(darray[j], darray[i]);
						while(ExOp::isGE(darray[j], swap)) {
							if (i == (--j)) break;
							if (darray[j] <= swap) break;
							if (i == (--j)) break;
						}
						//			printf("%i,%i maxt\n", i,j);
						if (i >= j) break;
						ExOp::memmove(darray[i], darray[j]);
					}
					ExOp::memmove(darray[i], swap);
					maxi.push(i);
					
					
			}
		}
	}

	
	LFHTEMP void Vector<C>::sort_unique(){
        if ((asize & 0x7FFFFFFF) == 0) return;
        this->sort();
		C *hare, *turtle;
        hare = darray;
         for(hare = darray; hare - darray < (asize & 0x7FFFFFFF) ; hare++){
            if (ExOp::isEQ(hare[0],hare[1])) break;
        }       
        if (hare - darray < (asize & 0x7FFFFFFF)){
            turtle = hare; hare++;
            for(hare++; hare - darray < (asize & 0x7FFFFFFF) ; hare++){
                if (ExOp::isNQ(*turtle,*hare)) {
                    turtle++;ExOp::memmove(*turtle, *hare);
                }
            }
            turtle++;
            this->DownSize(turtle - darray);
        } // else all were unique
	}

LFHTEMP void Vector<C>::reverse(){
	C swap;
	C* i;
	C* j;
	for(i=darray,j=darray + ((asize & 0x7FFFFFFF)-1); i<j; i--,j++){
		ExOp::memmove(swap,*i);
		ExOp::memmove(*i,*j);
		ExOp::memmove(*j,swap);	
	}
}

	LFHTEMP bool Vector<C>::issorted(){
		int i;
		for(i=(asize & 0x7FFFFFFF)-2;i>=0;i--) if (darray[i] > darray[i+1]) break;
		return(i < 0);
	}

LFHTEMP void Vector<C>::sort_decr(){
	C swap;
	int i;
	int j;
	stack<unsigned int> maxi;
	if ( (asize& 0x7FFFFFFF) == 0) return;
	maxi.push(asize & 0x7FFFFFFF);
	int mini = 0;
	while(!(maxi.empty())){
		//		printf("doing (%i, %i)!\n", mini, maxi.top()-1);
		switch(maxi.top() - mini){
			case 0:
			case 1:
				mini = maxi.top()+1; maxi.pop();
				break;
			case 2:
				if (ExOp::isLT(darray[mini], darray[mini+1])){
					ExOp::memmove(swap, darray[mini+1]);
					ExOp::memmove(darray[mini+1], darray[mini]);
					ExOp::memmove(darray[mini], swap);
				}
				mini += 3;  maxi.pop();
				break;
			case 3:
				if (ExOp::isGT(darray[mini], darray[mini+2])){
					if (darray[mini+1] < darray[mini]){
						ExOp::memmove(swap,darray[mini+1]);
						ExOp::memmove(darray[mini+1], darray[mini]);
						ExOp::memmove(darray[mini], swap);
					}else if (ExOp::isGT(darray[mini+2] , darray[mini+1])){
						ExOp::memmove(swap,darray[mini+2]);
						ExOp::memmove(darray[mini+2],darray[mini+1]);
						ExOp::memmove(darray[mini+1], swap);
					}
				}else{
					ExOp::memmove(swap, darray[mini+2]);
					if (ExOp::isGT(darray[mini+1] , darray[mini+2])){
						ExOp::memmove(darray[mini+2], darray[mini]);
						ExOp::memmove(darray[mini], darray[mini+1]);
						ExOp::memmove(darray[mini+1], swap);
					}else{
						if (ExOp::isGT(darray[mini] , darray[mini+1])){
							ExOp::memmove(darray[mini+2], darray[mini+1]);
							ExOp::memmove(darray[mini+1], darray[mini]);
						}else{
							ExOp::memmove(darray[mini+2], darray[mini]);
							
						}
						ExOp::memmove(darray[mini], swap);
					}
				}
				mini += 4;  maxi.pop();
				break;
			default:
				j = maxi.top()-1;
				i = (mini + j) >> 1; // printf("%i\t%i\n", i, j);
				switch( ((ExOp::isGT(darray[mini], darray[i])) ? 1 : 0) | ((darray[j-1] < darray[j]) ? 2 : 0) ){
					case 0:
						if (ExOp::isGT(darray[mini] , darray[j-1])){
							// daray[i] < darray[mini] < darray[j-1]  ?  darray[j]
							ExOp::memmove(swap, darray[mini]);
							ExOp::memmove(darray[mini], darray[i]);
							ExOp::memmove(darray[i], darray[j]);
							ExOp::memmove(darray[j], darray[j-1]);
						}else{
							// daray[j] < darray[j-1] < darray[mini]  ?  darray[i]
							ExOp::memmove(swap, darray[mini]);
							ExOp::memmove(darray[mini], darray[j]);
							ExOp::memmove(darray[j], swap);
							ExOp::memmove(swap, darray[j-1]);
						}
						break;
					case 1:
						
						if (ExOp::isGT(darray[i] ,darray[j-1])){
							// daray[mini] < darray[i] < darray[j-1]  ?  darray[j]
							ExOp::memmove(swap, darray[i]);
							ExOp::memmove(darray[i], darray[j]);
							ExOp::memmove(darray[j], darray[j-1]);
						}else{
							// daray[j] < darray[j-1] < darray[i]  ?  darray[mini]
							ExOp::memmove(swap, darray[i]);
							ExOp::memmove(darray[i], darray[mini]);
							ExOp::memmove(darray[mini], darray[j]);
							ExOp::memmove(darray[j], swap);
							ExOp::memmove(swap, darray[j-1]);
						}
						break;
					case 2:
						if (ExOp::isGT(darray[mini] , darray[j])){
							// daray[i] < darray[mini] < darray[j]  ?  darray[j-1]
							ExOp::memmove(swap, darray[mini]);
							ExOp::memmove(darray[mini], darray[i]);
							ExOp::memmove(darray[i], darray[j-1]);
						}else{
							// daray[j-1] < darray[j] < darray[mini]  ?  darray[i]
							ExOp::memmove(swap, darray[j]);
							ExOp::memmove(darray[j], darray[mini]);
							ExOp::memmove(darray[mini], darray[j-1]);
						}
						break;
					case 3:
						if (ExOp::isGT(darray[i] , darray[j])){
							// daray[mini] < darray[i] < darray[j]  ?  darray[j-1]
							ExOp::memmove(swap, darray[i]);
							ExOp::memmove(darray[i], darray[j-1]);
						}else{
							// daray[j-1] < darray[j] < darray[i]  ?  darray[mini]
							swap = darray[j];
							ExOp::memmove(darray[j], darray[i]);
							ExOp::memmove(darray[i], darray[mini]);
							ExOp::memmove(darray[mini], darray[j-1]);
						}
						break;
				}
				//		printf("pivot! %i!\n", swap);
				i = mini+1;
				j--;
				//		printf("pospos! %i,%i!\n", i,j);
				while(true){
					while(ExOp::isGE(darray[i] , swap)) {
						if ((++i) == j) break;
						if (darray[i] >= swap) break;
						if ((++i) == j) break;
					}
					//			printf("%i,%i mint\n", i,j);
					if (i >= j) break;
					ExOp::memmove(darray[j], darray[i]);
					while(ExOp::isLE(darray[j], swap)) {
						if (i == (--j)) break;
						if (darray[j] <= swap) break;
						if (i == (--j)) break;
					}
					//			printf("%i,%i maxt\n", i,j);
					if (i >= j) break;
					ExOp::memmove(darray[i], darray[j]);
				}
				ExOp::memmove(darray[i], swap);
				maxi.push(i);
				
				
		}
	}
}
LFHTEMP bool Vector<C>::issorted_decr(){
	int i;
	for(i=(asize & 0x7FFFFFFF)-2;i>=0;i--) if (darray[i] > darray[i+1]) break;
	return(i < 0);
}

	
	LFHTEMP void Vector<C>::random_permute(){
		unsigned int i;
		unsigned int j = asize & 0x7FFFFFFF;
		C swap;
		if (j < 2) return;
		do{
			i = rand() % (j--);
			if (i == j) continue;
			ExOp::memmove(swap, darray[j]);
			ExOp::memmove(darray[j], darray[i]);
			ExOp::memmove(darray[i], swap);
		}while (j > 1);
	}
	
	LFHTEMP Vector<complex> Vector<C>::bluesteinWindow(unsigned int size){
		unsigned int t = (size >> 2);
		int i;
		for(i=0;t != 0;i++) t >>=1;
		
        printf("%i is order\n",i); fflush(stdout);
		
		
        return(Vector<complex>());
		switch(i){
			case 1: return(Vector<complex>(Tuple<complex,16>::bluesteinWindow(size)));
			case 2: return(Vector<complex>(Tuple<complex,32>::bluesteinWindow(size)));
			case 3: return(Vector<complex>(Tuple<complex,64>::bluesteinWindow(size)));
			case 4: return(Vector<complex>(Tuple<complex,128>::bluesteinWindow(size)));
			case 5: return(Vector<complex>(Tuple<complex,256>::bluesteinWindow(size)));
			case 6: return(Vector<complex>(Tuple<complex,512>::bluesteinWindow(size)));
			case 7: return(Vector<complex>(Tuple<complex,1024>::bluesteinWindow(size)));
			case 8: return(Vector<complex>(Tuple<complex,2048>::bluesteinWindow(size)));
			case 9: return(Vector<complex>(Tuple<complex,4096>::bluesteinWindow(size)));
			case 10: return(Vector<complex>(Tuple<complex,8192>::bluesteinWindow(size)));
			case 11: return(Vector<complex>(Tuple<complex,16384>::bluesteinWindow(size)));
			case 12: return(Vector<complex>(Tuple<complex,32768>::bluesteinWindow(size)));
			case 13: return(Vector<complex>(Tuple<complex,65536>::bluesteinWindow(size)));
			default:
				return(Vector<complex>(Tuple<complex,8>::bluesteinWindow(size)));
				
		}
		
	}
	/*
	 LFHTEMP void Vector<C>::bluesteinWindow(complex *& fout, unsigned int size){
	 unsigned int t = (size >> 2);
	 int i;
	 for(i=0;t != 0;i++) t >>=1;
	 printf("%i is order\n",i); fflush(stdout);
	 switch(i){
	 case 1: Tuple<complex,16>::bluesteinWindow(fout,size);
	 case 2: Tuple<complex,32>::bluesteinWindow(fout,size);
	 case 3:Tuple<complex,64>::bluesteinWindow(fout,size);
	 case 4: Tuple<complex,128>::bluesteinWindow(fout,size);
	 case 5: Tuple<complex,256>::bluesteinWindow(fout,size);
	 case 6: Tuple<complex,512>::bluesteinWindow(fout,size);
	 case 7: Tuple<complex,1024>::bluesteinWindow(fout,size);
	 case 8: Tuple<complex,2048>::bluesteinWindow(fout,size);
	 case 9:Tuple<complex,4096>::bluesteinWindow(fout,size);
	 case 10: Tuple<complex,8192>::bluesteinWindow(fout,size);
	 case 11: Tuple<complex,16384>::bluesteinWindow(fout,size);
	 case 12: Tuple<complex,32768>::bluesteinWindow(fout,size);
	 case 13: Tuple<complex,65536>::bluesteinWindow(fout,size);
	 default:
	 Tuple<complex,8>::bluesteinWindow(fout,size);
	 }
	 }*/
	
	LFHTEMP Vector<C> Vector<C>::fourierTransform(const Vector<complex>& bluewindow) const{
		unsigned int t = (bluewindow.getSize())>>4;
		int i;
		for(i=0;t != 0;i++) t >>=1;
		
		switch(i){
			case 1: return(fourierTransform_routine<16>(bluewindow));
			case 2: return(fourierTransform_routine<32>(bluewindow));
			case 3: return(fourierTransform_routine<64>(bluewindow));
			case 4: return(fourierTransform_routine<128>(bluewindow));
			case 5: return(fourierTransform_routine<256>(bluewindow));
			case 6: return(fourierTransform_routine<512>(bluewindow));
			case 7: return(fourierTransform_routine<1024>(bluewindow));
			case 8: return(fourierTransform_routine<2048>(bluewindow));
			case 9: return(fourierTransform_routine<4096>(bluewindow));
			case 10: return(fourierTransform_routine<8192>(bluewindow));
			case 11: return(fourierTransform_routine<16384>(bluewindow));
			case 12: return(fourierTransform_routine<32768>(bluewindow));
			case 13: return(fourierTransform_routine<65536>(bluewindow));
			default:
				return(fourierTransform_routine<8>(bluewindow));
				
		}
		
	}
	LFHTEMP Vector<C> Vector<C>::invfourierTransform(const Vector<complex>& bluewindow) const{
		unsigned int t = (bluewindow.getSize())>>4;
		int i;
		for(i=0;t != 0;i++) t >>=1;
		
		switch(i){
			case 1: return(invfourierTransform_routine<16>(bluewindow));
			case 2: return(invfourierTransform_routine<32>(bluewindow));
			case 3: return(invfourierTransform_routine<64>(bluewindow));
			case 4: return(invfourierTransform_routine<128>(bluewindow));
			case 5: return(invfourierTransform_routine<256>(bluewindow));
			case 6: return(invfourierTransform_routine<512>(bluewindow));
			case 7: return(invfourierTransform_routine<1024>(bluewindow));
			case 8: return(invfourierTransform_routine<2048>(bluewindow));
			case 9: return(invfourierTransform_routine<4096>(bluewindow));
			case 10: return(invfourierTransform_routine<8192>(bluewindow));
			case 11: return(invfourierTransform_routine<16384>(bluewindow));
			case 12: return(invfourierTransform_routine<32768>(bluewindow));
			case 13: return(invfourierTransform_routine<65536>(bluewindow));
			default:
				return(invfourierTransform_routine<8>(bluewindow));
				
		}
	}
	
	
	
	LFHTEMP template<int supersize> Vector<C> Vector<C>::fourierTransform_routine(const Vector<complex>& bluewindow) const{
		
		Tuple<C, supersize> tmp;
		int i;
		int j=0;
		int s;
		int vsize = getSize();
		for(i=0;i<vsize;i++){
			double ang = -i * i * M_PI / vsize;
			complex factor = complex(cos(ang),sin(ang));
			tmp[j] = darray[i];
			tmp[j] *= factor;
			for(s = supersize >> 1;(j & s);s >>=1) j ^= s;
			j |= s;
		}
		for(;i<supersize;i++){
			tmp[j] = ExCo<C>::zero();
			for(s = supersize >> 1;(j & s);s >>=1) j ^= s;
			j |= s;
		}
		
		//				for(i=0;i<supersize;i++) printf("1:%i: %f\t%f\n",i,tmp[i][0],tmp[i][1]);
		tmp.fourierTransform_routine();
		//				for(i=0;i<supersize;i++) printf("1:%i: %f\t%f\n",i,tmp[i][0],tmp[i][1]);
		
		Tuple<C, supersize> tmp2;
		//		Tuple<C, supersize,Cflag> tmpb = bluesteinWindow<supersize>();
		
		//		for(i=0;i<supersize;i++) printf("b:%i: %f\t%f\n",i,tmpb[i][0],tmpb[i][1]);
		
		for(i=0;i<supersize;i++){
			tmp2[j]	= tmp[i];
			tmp2[j]	*= bluewindow[i];
			for(s = supersize >> 1;(j & s);s >>=1) j ^= s;
			j |= s;
		}
		
		//		for(i=0;i<supersize;i++) printf("2:%i: %f\t%f\n",i,tmp2[i][0],tmp2[i][1]);
		tmp2.invfourierTransform_routine();
		//		for(i=0;i<supersize;i++) printf("2:%i: %f\t%f\n",i,tmp2[i][0],tmp2[i][1]);
		
		Vector<C> _out;
		_out.setSize(vsize);
		for(i=0;i<vsize;i++) {
			double ang = -i * i * M_PI / vsize;
			complex factor = complex(cos(ang),sin(ang));
			_out[i] = tmp2[i];
			_out[i] *= factor;
			
		}
		return(_out);
	}
	
	
	LFHTEMP template<int supersize> Vector<C> Vector<C>::invfourierTransform_routine(const Vector<complex>& bluewindow) const{
		
		// implement Bluestein algorithm!
		
		// B windows : 2 {1.0+1.0i, 1.0-1.0i }
		// B windows : 4 {0.0+0.0i, 2.0-0.0i,(1 + i) * sqrt(2)/2, (1 + i) * sqrt(2)/2}
		
		
		Tuple<C, supersize> tmp;
		int i;
		int j=0;
		int s;
		int vsize = getSize();
		
		for(i=0;i<vsize;i++){
			double ang = i * i * M_PI / vsize;
			complex factor = complex(cos(ang),sin(ang));
			tmp[j] = darray[i];
			tmp[j] *= factor;
			for(s = supersize >> 1;(j & s);s >>=1) j ^= s;
			j |= s;
		}
		for(;i<supersize;i++){
			tmp[j] = ExCo<C>::zero();
			for(s = supersize >> 1;(j & s);s >>=1) j ^= s;
			j |= s;
		}
		
		//				for(i=0;i<supersize;i++) printf("1:%i: %f\t%f\n",i,tmp[i][0],tmp[i][1]);
		tmp.fourierTransform_routine();
		//				for(i=0;i<supersize;i++) printf("1:%i: %f\t%f\n",i,tmp[i][0],tmp[i][1]);
		
		Tuple<C, supersize> tmp2;
		//		Tuple<C, supersize,Cflag> tmpb = bluesteinWindow<supersize>();
		
		//		for(i=0;i<supersize;i++) printf("b:%i: %f\t%f\n",i,tmpb[i][0],tmpb[i][1]);
		
		for(i=0;i<supersize;i++){
			tmp2[j]	= tmp[i];
			complex factor =bluewindow[i];
			factor[1] = -factor[1];
			tmp2[j]	*= factor;
			for(s = supersize >> 1;(j & s);s >>=1) j ^= s;
			j |= s;
		}
		
		//		for(i=0;i<supersize;i++) printf("2:%i: %f\t%f\n",i,tmp2[i][0],tmp2[i][1]);
		tmp2.invfourierTransform_routine();
		//		for(i=0;i<supersize;i++) printf("2:%i: %f\t%f\n",i,tmp2[i][0],tmp2[i][1]);
		
		Vector<C> _out;
		_out.setSize(vsize);
		for(i=0;i<vsize;i++) {
			double ang = i * i * M_PI / vsize;
			complex factor = complex(cos(ang)/vsize,sin(ang)/vsize);
			_out[i] = tmp2[i];
			_out[i] *= factor;
			
		}
		return(_out);
		
	}
	
	/*
	 LFHTEMP void Vector<C>::fourierTransform_semiroutine(const Vector<complex>&  bluewindow) {
	 
	 unsigned int supersize = getSize();
	 int i;
	 int j=0;
	 int s;
	 int vsize = getSize();
	 for(i=0;i<vsize;i++){
	 double ang = -i * i * M_PI / vsize;
	 complex factor = complex(cos(ang),sin(ang));
	 tmp[j] = darray[i];
	 tmp[j] *= factor;
	 for(s = supersize >> 1;(j & s);s >>=1) j ^= s;
	 j |= s;
	 }
	 for(;i<supersize;i++){
	 tmp[j] = ExCo<C>::zero();
	 for(s = supersize >> 1;(j & s);s >>=1) j ^= s;
	 j |= s;
	 }
	 
	 
	 int step;
	 int cur;
	 int icur;
	 complex multip[size];
	 multip[0] = complex(1.0f,0.0f);
	 
	 for(step=1;step<size;step<<=1){
	 for(icur=1;icur<(step<<1);icur++){
	 double tmpdouble =  ((double)icur * M_PI) / step;
	 multip[icur] = complex(cos(tmpdouble),sin(tmpdouble));
	 }
	 for(cur =0;cur<size;cur+= (step<< 1)){
	 for(icur=0;icur<step;icur++){
	 //					printf("%f\t%f\t",((complex)data[cur | step | icur])[0],((complex)data[cur | step | icur])[1]);
	 C tmp = data[cur | step | icur];
	 tmp *= multip[step | icur];
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
	 
	 
	 
	 //				for(i=0;i<supersize;i++) printf("1:%i: %f\t%f\n",i,tmp[i][0],tmp[i][1]);
	 //	tmp.fourierTransform_routine();
	 //				for(i=0;i<supersize;i++) printf("1:%i: %f\t%f\n",i,tmp[i][0],tmp[i][1]);
	 
	 Tuple<C, supersize> tmp2;
	 //		Tuple<C, supersize,Cflag> tmpb = bluesteinWindow<supersize>();
	 
	 //		for(i=0;i<supersize;i++) printf("b:%i: %f\t%f\n",i,tmpb[i][0],tmpb[i][1]);
	 
	 for(i=0;i<supersize;i++){
	 tmp2[j]	= tmp[i];
	 tmp2[j]	*= bluewindow[i];
	 for(s = supersize >> 1;(j & s);s >>=1) j ^= s;
	 j |= s;
	 }
	 
	 //		for(i=0;i<supersize;i++) printf("2:%i: %f\t%f\n",i,tmp2[i][0],tmp2[i][1]);
	 //	tmp2.invfourierTransform_routine();
	 //		for(i=0;i<supersize;i++) printf("2:%i: %f\t%f\n",i,tmp2[i][0],tmp2[i][1]);
	 
	 Vector<C> _out;
	 _out.setSize(vsize);
	 for(i=0;i<vsize;i++) {
	 double ang = -i * i * M_PI / vsize;
	 complex factor = complex(cos(ang),sin(ang));
	 _out[i] = tmp2[i];
	 _out[i] *= factor;
	 
	 }
	 
	 
	 
	 }
	 
	 
	 LFHTEMP void Vector<C>::invfourierTransform_routine(Vector<C> & fout, const complex * bluewindow, unsigned char mag) const{
	 
	 // implement Bluestein algorithm!
	 
	 // B windows : 2 {1.0+1.0i, 1.0-1.0i }
	 // B windows : 4 {0.0+0.0i, 2.0-0.0i,(1 + i) * sqrt(2)/2, (1 + i) * sqrt(2)/2}
	 
	 
	 Tuple<C, supersize> tmp;
	 int i;
	 int j=0;
	 int s;
	 int vsize = getSize();
	 
	 for(i=0;i<vsize;i++){
	 double ang = i * i * M_PI / vsize;
	 complex factor = complex(cos(ang),sin(ang));
	 tmp[j] = darray[i];
	 tmp[j] *= factor;
	 for(s = supersize >> 1;(j & s);s >>=1) j ^= s;
	 j |= s;
	 }
	 for(;i<supersize;i++){
	 tmp[j] = ExCo<C>::zero();
	 for(s = supersize >> 1;(j & s);s >>=1) j ^= s;
	 j |= s;
	 }
	 
	 //				for(i=0;i<supersize;i++) printf("1:%i: %f\t%f\n",i,tmp[i][0],tmp[i][1]);
	 tmp.fourierTransform_routine();
	 //				for(i=0;i<supersize;i++) printf("1:%i: %f\t%f\n",i,tmp[i][0],tmp[i][1]);
	 
	 Tuple<C, supersize> tmp2;
	 //		Tuple<C, supersize,Cflag> tmpb = bluesteinWindow<supersize>();
	 
	 //		for(i=0;i<supersize;i++) printf("b:%i: %f\t%f\n",i,tmpb[i][0],tmpb[i][1]);
	 
	 for(i=0;i<supersize;i++){
	 tmp2[j]	= tmp[i];
	 complex factor =bluewindow[i];
	 factor[1] = -factor[1];
	 tmp2[j]	*= factor;
	 for(s = supersize >> 1;(j & s);s >>=1) j ^= s;
	 j |= s;
	 }
	 
	 //		for(i=0;i<supersize;i++) printf("2:%i: %f\t%f\n",i,tmp2[i][0],tmp2[i][1]);
	 tmp2.invfourierTransform_routine();
	 //		for(i=0;i<supersize;i++) printf("2:%i: %f\t%f\n",i,tmp2[i][0],tmp2[i][1]);
	 
	 Vector<C> _out;
	 _out.setSize(vsize);
	 for(i=0;i<vsize;i++) {
	 double ang = i * i * M_PI / vsize;
	 complex factor = complex(cos(ang)/vsize,sin(ang)/vsize);
	 _out[i] = tmp2[i];
	 _out[i] *= factor;
	 
	 }
	 return(_out);
	 
	 }*/
	
	
	// uses interval queries D must be comparable to C !
	LFHTEMP template<class D>
	void Vector<C>::getIntersection(const IntervalSet<D> &query, Vector<C> &out) const{
		int node = asize & 0x7FFFFFF;
		stack<unsigned int> indexes;
		int imin = 0;
		int imax = node -1;
		int imid;
		SetComparison cmp, cmp2;
		cmp = query.compareInterval((*this)[imin],(*this)[imax]);
		if (cmp.areMonotonicDisjoint()) return;
		while(true){
			while(true){
				//printf("%i,%i\n",imin,imax); fflush(stdout);
				if (imax - imin  < 2){
					for(imid = imin; imid <= imax;imid++){
						cmp = query.compare((*this)[imid]);
						if (cmp.contains()) out.push_back((*this)[imid]);
					}
					break;
				}else{
					imid = (imin + imax) >> 1;
					cmp = query.compareInterval((*this)[imin], (*this)[imid]);
					cmp2 = query.compareInterval((*this)[imid+1], (*this)[imax]);
					if (cmp.areDisjoint()){
						if (cmp2.areDisjoint()) break;
						else{
							imin = imid+1;
						}
					}else{
						if (!(cmp2.areDisjoint())){
							indexes.push(imid+1);
							indexes.push(imax);
						}
						imax = imid;
					}
				}
			}
			
			if (indexes.empty()) break;
			imax = indexes.top();indexes.pop();
			imin = indexes.top();indexes.pop();
		}
		
	}
	
	/*
	 template <class C,int flag>
	 C Vector<C>::operator[](int const which) const{
	 return(darray[which]);
	 }*/
	
	LFHTEMP
	C const & Vector<C>::operator[](int const which) const{
		return(darray[which]);
	}
	
	LFHTEMP
	template<class I> void Vector<C>::operator() (Oper1<I> const & op){ // not a match
		int i = asize -1;
		for(;i>=0;i--) darray[i](op);
	}
	LFHTEMP
	void Vector<C>::operator() (Oper1< C> const & op){ // match
		int i = asize -1;
		for(;i>=0;i--) op(darray[i], darray[i]);
	}
	
	LFHTEMP
	template<class A_1, class A_2, class C_I> void Vector<C>::operator() (Oper2<A_1,A_2> const & op, const Vector<C_I> & a_2 ){ // not a match
		//setSize(_in.getSize());
		int i = asize -1;
		for(;i>=0;i--) (darray[i])(op,a_2.darray[i]);
	}
	
	LFHTEMP
	template<class C_I> void Vector<C>::operator() (Oper2<C,C_I> const & op, const Vector<C_I> & a_2){ // match
		//setSize(_in.getSize());
		int i = asize -1;
		for(;i>=0;i--) op(darray[i], a_2.darray[i]);
	}
	
	LFHTEMP
	template<class A_1, class A_2, class C_I> void Vector<C>::operator() (Oper2<A_1,A_2> const & op, Vector<C_I> & a_2 ){ // not a match
		//setSize(_in.getSize());
		int i = asize -1;
		for(;i>=0;i--) (darray[i])(op,a_2.darray[i]);
	}
	
	LFHTEMP
	template<class C_I> void Vector<C>::operator() (Oper2<C,C_I> const & op, Vector<C_I> & a_2){ // match
		//setSize(_in.getSize());
		int i = asize -1;
		for(;i>=0;i--) op(darray[i], a_2.darray[i]);
	}
	
	LFHTEMP
	template<class A_1, class A_2, class A_3, class C_2, class C_3> void Vector<C>::operator() (Oper3<A_1,A_2,A_3> const & op, Vector<C_2> const & a_2, Vector<C_3> const & a_3 ){ // not a match
		int i = asize -1;
		for(;i>=0;i--) (darray[i])(op,a_2.darray[i],a_3.darray[i]);
	}
	
	LFHTEMP
	template<class C_2, class C_3> void Vector<C>::operator() (Oper3<C,C_2,C_3> const & op, Vector<C_2> const & a_2, Vector<C_3> const & a_3 ){ // not a match
		int i = asize -1;
		for(;i>=0;i--) op(darray[i],a_2.darray[i],a_3.darray[i]);
	}
	
	LFHTEMP
	template<class A_1, class A_2, class A_3, class C_2, class C_3> void Vector<C>::operator() (Oper3<A_1,A_2,A_3> const & op, Vector<C_2> & a_2, Vector<C_3> const & a_3 ){ // not a match
		int i = asize -1;
		for(;i>=0;i--) (darray[i])(op,a_2.darray[i],a_3.darray[i]);
	}
	
	LFHTEMP
	template<class C_2, class C_3> void Vector<C>::operator() (Oper3<C,C_2,C_3> const & op, Vector<C_2> & a_2, Vector<C_3> const & a_3 ){ // not a match
		int i = asize -1;
		for(;i>=0;i--) op(darray[i],a_2.darray[i],a_3.darray[i]);
	}
	
	LFHTEMP
	template<class A_1, class A_2, class A_3, class C_2, class C_3> void Vector<C>::operator() (Oper3<A_1,A_2,A_3> const & op, Vector<C_2> & a_2, Vector<C_3> & a_3 ){ // not a match
		int i = asize -1;
		for(;i>=0;i--) (darray[i])(op,a_2.darray[i],a_3.darray[i]);
	}
	
	LFHTEMP
	template<class C_2, class C_3> void Vector<C>::operator() (Oper3<C,C_2,C_3> const & op, Vector<C_2> & a_2, Vector<C_3> & a_3 ){ // not a match
		int i = asize -1;
		for(;i>=0;i--) op(darray[i],a_2.darray[i],a_3.darray[i]);
	}
	
	template <class C>
	angle<C>::angle(C const & value) : ang(value){}
	template <class C>
	angle<C>::operator double() const{
		return( (double)ang);
	}
	
	template <class C>
	angle<C>::operator complex() const{
		double a = (double)(*this);
		return(complex(cos(a),sin(a)));
	}
	template <class C>
	complex angle<C>::getcomplex() const{
		double a = (double)(*this);
		return(complex(cos(a),sin(a)));
	}
	
	
	LFHTEMP
	void Vector<C>::GP_Covariance(DataGrid<double, 2> &f_out, double (*metric)(const C&, const C&), double noise_variance) const{
		
		unsigned int coor[2];
		coor[1] =coor[0] = getSize();
		f_out.setSizes(coor);
		
		double tmp;
		
		typedef class DataGrid<double,2>::KeyIterator itetype;
		itetype ite = f_out.getKeyIterator();
		
		if (ite.first()) do{
			
			if (ite()[0] == ite()[1]) {
				f_out(ite()) = 1.0f + noise_variance;
			}else{
				tmp = metric((*this)[ite()[0]], (*this)[ite()[1]]);
				f_out(ite()) = exp(-0.5f * tmp * tmp);
			}
			
		} while(ite.next());
		
		
		
	}
	
	LFHTEMP  void Vector<C>::GP_Covariance_cross(DataGrid<double, 2> &f_out, double (*metric)(const C&, const C&), Vector<C>& query) const{
		
		unsigned int coor[2];
		coor[0] = getSize();
		coor[1] = query.getSize();
		
		f_out.setSizes(coor);
		
		double tmp;
		
		typedef class DataGrid<double,2>::KeyIterator itetype;
		itetype ite = f_out.getKeyIterator();
		
		if (ite.first()) do{
			tmp = metric((*this)[ite()[0]], query[ite()[1]]);
			f_out(ite()) = exp(-0.5f * tmp * tmp);
		} while(ite.next());
		
		
		
	}
	
	LFHTEMP
	void Vector<C>::GP_fromCorrel(DataGrid<double, 2> &f_out, double (*metric)(const C&, const C&)) const{
		
		unsigned int coor[2];
		coor[1] =coor[0] = getSize();
		f_out.setSizes(coor);
		
		double tmp;
		
		typedef class DataGrid<double,2>::KeyIterator itetype;
		itetype ite = f_out.getKeyIterator();
		
		if (ite.first()) do{
			
			if (ite()[0] == ite()[1]) {
				f_out(ite()) = metric((*this)[ite()[0]], (*this)[ite()[1]]);
			}else{
				f_out(ite()) = metric((*this)[ite()[0]], (*this)[ite()[1]]);
			}
			
		} while(ite.next());
		
		
		
	}
	
	LFHTEMP void Vector<C>::GP_fromCorrel_cross(DataGrid<double, 2> &f_out, Vector<C>& query, double (*correl)(const C&, const C&)) const{
		
		unsigned int coor[2];
		coor[0] = getSize();
		coor[1] = query.getSize();
		
		f_out.setSizes(coor);
		
		double tmp;
		
		typedef class DataGrid<double,2>::KeyIterator itetype;
		itetype ite = f_out.getKeyIterator();
		
		if (ite.first()) do{
			tmp = metric((*this)[ite()[0]], query[ite()[1]]);
			f_out(ite()) = exp(-0.5f * tmp * tmp);
		} while(ite.next());
		
		
		
	}
	
	LFHTEMP void Vector<C>::GP_fromCorrel_complete(DataGrid<double, 2> &f_out, Vector<C>& query, double (*corr_metric)(const C&, const C&)) const{
		
		unsigned int coor[2];
		coor[0] = getSize();
		coor[1] = getSize() + query.getSize();
		
		f_out.setSizes(coor);
		int s = getSize();
		
		typedef class DataGrid<double,2>::KeyIterator itetype;
		itetype ite = f_out.getKeyIterator();
		
		if (ite.first()) do{
			if (ite()[1] < s){
				f_out(ite()) = corr_metric((*this)[ite()[0]], (*this)[ite()[1]]);
			}else{
				f_out(ite()) = corr_metric((*this)[ite()[0]], query[ite()[1] - s]);
			}
		} while(ite.next());
		
	}
	
    LFHTEMP Matrix<C> Vector<C>::mkouterprod() const{unsigned int s = asize & 0x7FFFFFFF; Matrix<C> fout(s, s);
        unsigned i,j;
        C* cur = fout.data;
        for(j=0;j<s;j++) for(i=0;i<s;i++) *(cur++) = darray[i] * darray[j];
        return(fout);
    }




