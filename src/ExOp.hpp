/*
 * ExOp.hpp
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

#pragma once

#undef LFHTEMP	
#define LFHTEMP template<class A>

LFHTEMP bool ExCo<vector<A>,0u >::isValid(const vector<A> &a){unsigned int i; for(i=0;i< a.size();i++) if (! ExOp::isValid(a[i])) break; return i == a.size();}
LFHTEMP bool ExCo<vector<A>,0u >::isZero(const vector<A> &a){unsigned int i; for(i=0;i< a.size();i++) if (! ExOp::isZero(a[i])) break; return i == a.size();}
LFHTEMP bool ExCo<vector<A>,0u >::isOne(const vector<A> &a){unsigned int i; for(i=0;i< a.size();i++) if (! ExOp::isOne(a[i])) break; return i == a.size();}




#undef LFHTEMP
#define LFHTEMP template<class A, unsigned int SIZE>

LFHTEMP void ExCo<A[SIZE],0u>::zero( A (&a)[SIZE]){	for(unsigned int i=0;i<SIZE;i++) ExOp::toZero(a[i]);}
LFHTEMP void ExCo<A[SIZE],0u>::random( A (&a)[SIZE]){	for(unsigned int i=0;i<SIZE;i++) ExOp::toRand(a[i]);}
LFHTEMP void ExCo<A[SIZE],0u>::one( A (&a)[SIZE]){	for(unsigned int i=0;i<SIZE;i++) ExOp::toOne(a[i]);}
LFHTEMP void ExCo<A[SIZE],0u>::minimum( A (&a)[SIZE]){	for(unsigned int i=0;i<SIZE;i++) ExOp::toMin(a[i]);}
LFHTEMP void ExCo<A[SIZE],0u>::maximum( A (&a)[SIZE]){	for(unsigned int i=0;i<SIZE;i++) ExOp::toMax(a[i]);}
LFHTEMP void ExCo<A[SIZE],0u>::show(const A (&a)[SIZE], FILE* f_out, int l){
	switch(l){
		case 0:
			for(unsigned int i=0;i<SIZE;i++) {
				ExOp::show(a[i],f_out,1);
				fprintf(f_out,"\n");
			}
			fprintf(f_out,"\n");
			break;
		case 1:
			for(unsigned int i=0;i<SIZE;i++) {
				if (i != 0) fprintf(f_out,"\t");
				ExOp::show(a[i],f_out,2);
			}
			if (l==0) fprintf(f_out,"\n");
			break;
		case 2:
			fprintf(f_out,"[");
			for(unsigned int i=0;i<SIZE;i++) {
				if (i != 0) fprintf(f_out,";");
				ExOp::show(a[i],f_out,3);
			}
			fprintf(f_out,"]");
			break;
		default:
			fprintf(f_out,"(");
			for(unsigned int i=0;i<SIZE;i++) {
				if (i != 0) fprintf(f_out,",");
				ExOp::show(a[i],f_out,4);
			}
			fprintf(f_out,")");
			break;
	}
}


#undef LFHTEMP	
#define LFHTEMP template<class A>

// const mappings

LFHTEMP void ExCo<A*, 0u>::show(const A * & what, FILE *f,int level) {ExOp::show(*what,f, level);}

#undef LFHTEMP
#define LFHTEMP template<class A, unsigned int flag>

// assignments		


LFHTEMP double ExCo<A,flag>::pdist(const A& a, const A& b){return ExOp::pnorm(a-b);}
LFHTEMP void ExCo<A,flag>::tointpow(A& what, const int pow) {
		int modp;
		if (pow >1){
		 modp = pow;
		A cp = what;
		while(modp > 1){
			ExOp::tosquare(what);
			if (modp & 1) ExOp::tomult(what, cp);
			modp = modp >> 1;
			}
		}else if (pow < 0){
			ExOp::toinverse(what);
			if (pow < -1){
			 modp = -pow;
			 A cp = what;
			 while(modp > 1){
			ExOp::tosquare(what);
			if (modp & 1) ExOp::tomult(what, cp);
			modp = modp >> 1;
			}
			}
			}else if (pow == 0) ExOp::toOne(what);
		}
LFHTEMP A ExCo<A,flag>::mkintpow(const A& what, const int pow){
		int modp;
		if (pow >1){
		 modp = pow;
		A cp = what;
		A fout = what;
		while(modp > 1){
			ExOp::tosquare(fout);
			if (modp & 1) ExOp::tomult(fout, cp);
			modp = modp >> 1;
			}
			return(fout);
		}else if (pow < 0){
			A fout = mkinverse(what);

			if (pow < -1){
			 modp = -pow;
			 A cp = fout;
			 while(modp > 1){
			ExOp::tosquare(fout);
			if (modp & 1) ExOp::tomult(fout, cp);
			modp = modp >> 1;
			}
			}
			return(fout);
			}else if (pow == 0) return ExOp::mkone<A>();
			else return what;
	}

// builds

LFHTEMP A ExCo<A,flag>::mkinvintpow(const A& what, const int pow){return ExOp::mkintpow(ExOp::mkinverse(what), pow);}
LFHTEMP A ExCo<A,flag>::mknegative(const A& a){return ExOp::mksub(ExOp::mkzero<A>(), a);}
LFHTEMP A ExCo<A,flag>::mkinverse(const A& a){return ExOp::mkdivi(ExOp::mkone<A>(), a);}
LFHTEMP A ExCo<A,flag>::mksquare(const A& a){return ExOp::mkmult(a,a);}


	// assignments

LFHTEMP bool ExCo<A,flag>::isGT(const A &a,const A &b){return((ExOp::setcmp(a,b) & SETCMP_CMP_MASK2) == SETCMP_MASKED_GT);}
LFHTEMP bool ExCo<A,flag>::isGE(const A &a,const A &b){return((ExOp::setcmp(a,b) & SETCMP_CMP_MASK2) != SETCMP_MASKED_LT);}
LFHTEMP bool ExCo<A,flag>::isLT(const A &a,const A &b){return((ExOp::setcmp(a,b) & SETCMP_CMP_MASK2) == SETCMP_MASKED_LT);}
LFHTEMP bool ExCo<A,flag>::isLE(const A &a,const A &b){return((ExOp::setcmp(a,b) & SETCMP_CMP_MASK2) != SETCMP_MASKED_GT);}
LFHTEMP bool ExCo<A,flag>::isEQ(const A &a,const A &b){return((ExOp::setcmp(a,b) & SETCMP_DISJOINT) == SETCMP_EQUAL);}
LFHTEMP bool ExCo<A,flag>::isNQ(const A &a,const A &b){return((ExOp::setcmp(a,b) & SETCMP_DISJOINT) == SETCMP_DISJOINT);}



	// const mappings
#undef LFHTEMP
#define LFHTEMP template<class A>

LFHTEMP void ExCo<const A,0u>::show(const A &a, FILE* fout, int level){ExCo<A>::show(a,fout,level);}
LFHTEMP bool ExCo<const A,0u>::isValid(const A& a){return(ExCo<A>::isValid(a));}



#undef LFHTEMP
#define LFHTEMP template<class A, class B>


LFHTEMP bool ExCo<pair<A,B> ,0u>::isValid(const pair<A,B> &a){ return ExOp::isValid(a.first) &&  ExOp::isValid(a.second) ;}

LFHTEMP bool ExCo<pair<A,B> ,0u>::isZero(const pair<A,B> &a){return ExOp::isZero(a.first) && ExOp::isZero(a.second) ;}
LFHTEMP bool ExCo<pair<A,B> ,0u>::isOne(const pair<A,B> &a){return ExOp::isOne(a.first) && ExOp::isOne(a.second) ;}

LFHTEMP void ExCo<pair<A,B> ,0u>::zero(pair<A,B> &a){ExOp::toZero(a.first); ExOp::toZero(a.second);}
LFHTEMP void ExCo<pair<A,B>,0u >::one(pair<A,B> &a){ExOp::toOne(a.first); ExOp::toOne(a.second);}
LFHTEMP void ExCo<pair<A,B>,0u >::random(pair<A,B> &a){ExOp::toRand(a.first); ExOp::toRand(a.second);}

LFHTEMP void ExCo<pair<A,B> ,0u>::show(const pair<A,B> &a, FILE* f_out, int l){
		switch(l){
		case 0:
			ExOp::show(a.first,f_out,1);
			fprintf(f_out,"\n");
			ExOp::show(a.second,f_out,1);
			fprintf(f_out,"\n\n");
			break;
		case 1:
			ExOp::show(a.first,f_out,2);
			fprintf(f_out,"\t");
			ExOp::show(a.second,f_out,2);
			fprintf(f_out,"\n");
			break;
		case 2:
			fprintf(f_out,"[");
			ExOp::show(a.first,f_out,3);
			fprintf(f_out,";");
			ExOp::show(a.second,f_out,3);
			fprintf(f_out,"]");
			break;
		default:
			fprintf(f_out,"(");
			ExOp::show(a.first,f_out,4);
			fprintf(f_out,",");
			ExOp::show(a.second,f_out,4);
			fprintf(f_out,")");
			break;
	}

	}

LFHTEMP string ExCo<pair<A,B> ,0u>::type_tostring(const pair<A,B> &a){return string("pair<") + ExOp::type_tostring(a.first) +string(",")+ ExOp::type_tostring(a.second) + string(">");}

	// ExFn specializations


	template<class F, class A, class B, unsigned int dims> void ExFn<false>::compose(F &func, DataGrid<A,dims> &a, const DataGrid<B,dims> &b)	{
			a.setSizes(b.dims);
			typename DataGrid<A,dims>::KeyIterator ite = a.getKeyIterator();
			if (ite.first()) do{
				ExOp::comp(func,(A&) a(ite()),(const B&) b(ite()));
				} while(ite.next());
			}

	template<class F, class A, class B, unsigned int DA, unsigned int DB> inline static void compose(F &func, DataGrid<A,DA> &a, const DataGrid<B,DB> &b){
		 if (DA < DB){
			 a.setSizes(b.dims + DB-DA);
			 			typename DataGrid<A,DA>::KeyIterator ite = a.getKeyIterator();
			if (ite.first()) do{
				ExOp::comp(func,(A&)a(ite()),(B&)b(ite()));
				} while(ite.next());
			 }else{
		 }
	}


