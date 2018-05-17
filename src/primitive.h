/*
 * primitive.h
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

// optionnal external dependency, which uses a modified file "GSLfunc.hpp"
// #define GNU_SCIENTIFIC_LIBRARY


#ifndef _defined_Primitive
#define _defined_Primitive

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stddef.h>

#include <string.h>

#include <time.h>
#include <math.h>
#include <map>

#include <list>
#include <queue>
#include <stack>
#include <vector>
#include <iostream>

#include <float.h>

using namespace std;

#define LFH_WARNINGS true
#define LFH_STR(tok) #tok

#define STRINGIFY(A) #A


//#define LFH_USE_MEMLEAK_HELP(CMD)   CMD
#define LFH_USE_MEMLEAK_HELP(CMD)   


// Target Table
//                     One client                 Many clients             Dependent Client (needs simutaneous destory )
//             Immovible        Ammovible  Immovible      Ammovible        Immovible      Ammovible
// Permanent   pointer        double-ptr     ptr          alias            ptr            alias/listener
// Temporary   double-ptr     double-ptr   alias          alias         alias/listener    alias/listener
//
#define LFH_DEFAULT_PORT "27015"

#define LFH_GOLD
#define LFH_FAULTY


#ifndef M_EULER
#define M_EULER    0.57721566490153286060651209008      /* Euler constant */
#endif

#ifndef M_SQRTPI
#define M_SQRTPI   1.77245385090551602729816748334      /* sqrt(pi) */
#endif

#define LFHDebug_RETURN_void(Test) if (!(Test)) return
#define LFHDebug_RETURN(Test, out) if (!(Test)) return(out)

#define LFHDebug(Test, String) if (Test) fprintf(stderr, String)
#define LFHDebug_LOOP(Test) while(true)

// 8bit guarrantied
#define schar char
#define uchar unsigned char

// 16bit guarrantied
#define sshort short
#define ushort unsigned short

// 32bit guarrantied
#define sint int
#define uint unsigned int

// 64bit guarrantied
#define slong int
#define ulong unsigned int

#define LFHPRIMITIVE_ANGCONVERT 0.000095873799242852576857380474343247f

//#define myexit(text) exit(MessageBox(NULL,text,"ERROR",MB_OK|MB_ICONSTOP))
#define myexit(text) exit(1 | fprintf(stderr,text))


#define mybreak(text) MessageBox(NULL,text,"ERROR",MB_OK|MB_ICONSTOP)

#ifndef M_PI_2
#define M_PI_2 ((double)1.5707963267948966192313216916398)
#endif

#ifndef M_PI
#define M_PI ((double)3.1415926535897932384626433832795)
#endif

#ifndef M_SQRT_1_2
#define M_SQRT_1_2 ((double)0.70710678118654752440084436210485)
#endif

#ifndef M_FLOAT_NAN
#define M_FLOAT_NAN ((float)0.0f/ 0.0f)
#endif

#define M_2PI (6.283185307179586476925286766559)
#define AngConvert 0.000095873799242852576857380474343247

#define LFHCLONE(PoInTeR, ClAsSnAmE) ((PoInTeR == NULL) ? NULL : (ClAsSnAmE)(PoInTeR->clone()) )

#define LFHDECL_DESCTRUCTOR(ClAsSnAmE) virtual ~ClAsSnAmE(); ClAsSnAmE(const ClAsSnAmE & other); ClAsSnAmE& operator=(const ClAsSnAmE & other);
#define LFHDECL_VIRTUALCLONE(ClAsSnAmE) virtual void* clone()=0;
#define LFHDECL_CLONE(ClAsSnAmE) virtual void* clone() {return (void*) new ClAsSnAmE(*this) ;}

#define LFHDECL_COPYCON_BODY(nAmEsPaCe , ClAsSnAmE) nAmEsPaCe::ClAsSnAmE ::ClAsSnAmE(const nAmEsPaCe::ClAsSnAmE & other)
#define LFHDECL_COPYCON_BODY_AUTO(nAmEsPaCe , ClAsSnAmE) nAmEsPaCe::ClAsSnAmE ::ClAsSnAmE(const nAmEsPaCe::ClAsSnAmE & other)
#define LFHDECL_ASSIGN_BODY(nAmEsPaCe , ClAsSnAmE , BoDy) nAmEsPaCe::ClAsSnAmE & nAmEsPaCe::ClAsSnAmE ::operator=(const nAmEsPaCe::ClAsSnAmE & other)
#define LFHDECL_CLONEINIT(MeMbEr , TyPe) MeMbEr = (TyPe)other.MeMbEr->clone()
#define LFHDECL_CLONEASSIGN(MeMbEr , TyPe) delete(MeMbEr); MeMbEr = (TyPe)other.MeMbEr->clone()

#define LFHDECL_OPER1(ClAsS_A) public Oper1<ClAsS_A> { public: void operator()(ClAsS_A &) const;}
#define LFHDECL_OPER2(ClAsS_A, ClAsS_B) public Oper2<ClAsS_A, ClAsS_B> { public: void operator()(ClAsS_A &, ClAsS_B &) const;}
#define LFHDECL_OPER3(ClAsS_A, ClAsS_B, ClAsS_C) public Oper3<ClAsS_A, ClAsS_B, ClAsS_C> { public: void operator()(ClAsS_A &, ClAsS_B &, ClAsS_C &) const;}

// more of a reminder:
#define LFHCONCAT2E(a,b) a##b

#define LFHCONCAT2(a,b) a,b
#define LFHCONCAT3(a,b,c) a,b,c
#define LFHDECL_TRIVIAL_OPERATOR(oPeRaToR, oUtPuT) template<class TrIvIaLaRgUmEnT> oUtPuT& operator oPeRaToR(TrIvIaLaRgUmEnT const &other) {return( ((oUtPuT(*this)) oPeRaToR##= other) );}

#define LFHALIVE(text) fprintf(stdout, "Debug flag: %s on %s::line#%d  (time=%i)\n", text, __FILE__, __LINE__, clock());fflush(stdout)
#define LFH_STOPPOPUP(a) printf("%c", printf(a) == -1 ? '\t' : '\n')
#define LFH_ENDPOPUP(a) exit(printf("%c", printf(a) == -1 ? '\t' : '\n'))

#define LFH_address unsigned int
#define myoffsetof(A,B,C) (int(&(((A*)0)->B)))

#define LFH_ALIVE printf("%i th line in %s reached within func %s\n", __LINE__ , __FILE__ ,  __func__) | fflush(stdout) 
#define LFH_FUNCTION_MISSING printf("function  %s in %s is not done yet!\n", __func__ , __FILE__);fflush(stdout);exit(1)

#ifdef GNU_SCIENTIFIC_LIBRARY
struct gsl_sf_result_struct {
	double val;
	double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
	double val;
	double err;
	int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;
#endif


namespace LFHPrimitive{

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Meta-programmation section
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	template <bool A> class YESNO{	public:	enum {ans = A};	};

	enum FLAGFORTYPE{
		FLAGFORTYPE_NULL =0,
		FLAGFORTYPE_CONST =1,
		FLAGFORTYPE_REF =2,
		FLAGFORTYPE_CONSTREF =3,
		};

	template< class A, FLAGFORTYPE isconst>
	class TYPEMANIP{public:typedef A TYPE;};

	template<class A> class TYPEMANIP<A, FLAGFORTYPE_CONST>{public: typedef const A TYPE;};
	template<class A> class TYPEMANIP<A, FLAGFORTYPE_REF>{public:typedef A& TYPE;};
	template<class A> class TYPEMANIP<A, FLAGFORTYPE_CONSTREF>{	public:	typedef const A& TYPE;	};

    template<class A> class TYPEMANIP<const A, FLAGFORTYPE_NULL>{public: typedef A TYPE;};
	template<class A> class TYPEMANIP<const A, FLAGFORTYPE_REF>{public:typedef A& TYPE;};
	template<class A> class TYPEMANIP<const A, FLAGFORTYPE_CONSTREF>{	public:	typedef const A& TYPE;	};

    template<class A> class TYPEMANIP<A&, FLAGFORTYPE_NULL>{public: typedef A TYPE;};
	template<class A> class TYPEMANIP<A&, FLAGFORTYPE_CONST>{public: typedef const A TYPE;};
	template<class A> class TYPEMANIP<A&, FLAGFORTYPE_CONSTREF>{	public:	typedef const A& TYPE;	};

    template<class A> class TYPEMANIP<const A&, FLAGFORTYPE_NULL>{public: typedef A TYPE;};
	template<class A> class TYPEMANIP<const A&, FLAGFORTYPE_CONST>{public: typedef const A TYPE;};
	template<class A> class TYPEMANIP<const A&, FLAGFORTYPE_REF>{public:typedef A& TYPE;};

	extern unsigned int def_font[];
    enum metaop_type{
        METAOP_ZERO,
        METAOP_AND,
        METAOP_OR,
        METAOP_XOR,
        METAOP_PLUS,
        METAOP_MINU,
        METAOP_PROD,
        METAOP_DIVI,
        METAOP_MOD,
    };

    template<metaop_type METAOP, class A, class B, class TYPE = int>
    class MetaOperation{
    public:
        static const TYPE val=0;
    };

    template<class A, class B, class TYPE> class MetaOperation<METAOP_AND,A,B,TYPE>{ public: static const TYPE val=(A::val) & (B::val); };
    template<class A, class B, class TYPE> class MetaOperation<METAOP_OR,A,B,TYPE>{ public: static const TYPE val=(A::val) | (B::val); };
    template<class A, class B, class TYPE> class MetaOperation<METAOP_XOR,A,B,TYPE>{ public: static const TYPE val=(A::val) ^ (B::val); };
    template<class A, class B, class TYPE> class MetaOperation<METAOP_PLUS,A,B,TYPE>{ public: static const TYPE val=(A::val) + (B::val); };
    template<class A, class B, class TYPE> class MetaOperation<METAOP_MINU,A,B,TYPE>{ public: static const TYPE val=(A::val) - (B::val); };
    template<class A, class B, class TYPE> class MetaOperation<METAOP_PROD,A,B,TYPE>{ public: static const TYPE val=(A::val) * (B::val); };
    template<class A, class B, class TYPE> class MetaOperation<METAOP_DIVI,A,B,TYPE>{ public: static const TYPE val=(A::val) / (B::val); };
    template<class A, class B, class TYPE> class MetaOperation<METAOP_MOD,A,B,TYPE>{ public: static const TYPE val=(A::val) % (B::val); };

	// meant for small fatorial only! 12! at most
	template <int base, int nbfact = (base< 13) ? base : 0> class TEMPLATE_FACTORIAL{
public:
	enum {ans = base * TEMPLATE_FACTORIAL<base-1, nbfact-1>::ans};
};

template<unsigned int which> class TEMPLATE_PRIME{
public:	
    static const unsigned int val = 0x7FFFFFFF;
};
template< >  class TEMPLATE_PRIME<1u>{public:	static const unsigned int val = 2;};
template< >  class TEMPLATE_PRIME<2u>{public:	static const unsigned int val = 3;};
template< >  class TEMPLATE_PRIME<3u>{public:	static const unsigned int val = 5;};
template< >  class TEMPLATE_PRIME<4u>{public:	static const unsigned int val = 7;};
template< >  class TEMPLATE_PRIME<5u>{public:	static const unsigned int val = 11;};
template< >  class TEMPLATE_PRIME<6u>{public:	static const unsigned int val = 13;};
template< >  class TEMPLATE_PRIME<7u>{public:	static const unsigned int val = 17;};
template< >  class TEMPLATE_PRIME<8u>{public:	static const unsigned int val = 19;};
template< >  class TEMPLATE_PRIME<9u>{public:	static const unsigned int val = 23;};
template< >  class TEMPLATE_PRIME<10u>{public:	static const unsigned int val = 29;};
template< >  class TEMPLATE_PRIME<11u>{public:	static const unsigned int val = 31;};
template< >  class TEMPLATE_PRIME<12u>{public:	static const unsigned int val = 37;};

template <int base> class TEMPLATE_FACTORIAL<base,0>{
public:
	enum {ans = 1};
};

template<unsigned int which>
class StConstList{};
template< > class StConstList<0u>{public: static const unsigned int prime = 2;};
template< > class StConstList<1u>{public: static const unsigned int prime = 3;};
template< > class StConstList<2u>{public: static const unsigned int prime = 5;};
template< > class StConstList<3u>{public: static const unsigned int prime = 7;};
template< > class StConstList<4u>{public: static const unsigned int prime = 11;};
template< > class StConstList<5u>{public: static const unsigned int prime = 13;};
template< > class StConstList<6u>{public: static const unsigned int prime = 17;};
template< > class StConstList<7u>{public: static const unsigned int prime = 19;};
template< > class StConstList<8u>{public: static const unsigned int prime = 23;	typedef char SINT;	typedef unsigned char UINT;};
template< > class StConstList<9u>{public: static const unsigned int prime = 29;};
template< > class StConstList<10u>{public: static const unsigned int prime = 31;};
template< > class StConstList<11u>{public: static const unsigned int prime = 37;};
template< > class StConstList<12u>{public: static const unsigned int prime = 43;};
template< > class StConstList<16u>{public: static const unsigned int prime = 2;	typedef short SINT;	typedef unsigned short UINT;};
template< > class StConstList<32u>{public: static const unsigned int prime = 2;	typedef int SINT;	typedef unsigned int UINT;};
template< > class StConstList<64>{public: static const unsigned int prime = 2;typedef long int SINT; typedef unsigned long int UINT;};


template <unsigned int val, unsigned int base, unsigned int ite = 1u, unsigned int div = 1u> class TEMPLATE_DIVIDE_BY_MAX_DIVISOR{
public:
static const unsigned int ans = ( TEMPLATE_PRIME<ite>::val > base ? val :
			 ((((val % TEMPLATE_PRIME<ite>::val) == 0u)&&(TEMPLATE_PRIME<ite>::val <= base / div)) ?
			  TEMPLATE_DIVIDE_BY_MAX_DIVISOR<val / TEMPLATE_PRIME<ite>::val, base, ite, div * TEMPLATE_PRIME<ite>::val >::ans   :
			  TEMPLATE_DIVIDE_BY_MAX_DIVISOR<val, base, ite+1, 1u>::ans)  );
};

template <unsigned int val, unsigned int base, unsigned int div> class TEMPLATE_DIVIDE_BY_MAX_DIVISOR<val, base,11u,div>{
public:	static const unsigned int ans = val;};

template <unsigned int base, unsigned int ite, unsigned int div> class TEMPLATE_DIVIDE_BY_MAX_DIVISOR<0u, base,ite,div>{
public:	static const unsigned int ans = 0;};

template <unsigned int lenght,unsigned int nbdims> class TEMPLATE_TRIANGLE_NUMBER{
public:	enum {ans = 0};};

template <unsigned int lenght> class TEMPLATE_TRIANGLE_NUMBER<lenght, 1u>{	public:	static const unsigned int ans = lenght;	};
template <unsigned int lenght> class TEMPLATE_TRIANGLE_NUMBER<lenght, 2u>{	public: static const unsigned int ans = (TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght,2u>::ans * TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght+1,2u>::ans);};
template <unsigned int lenght> class TEMPLATE_TRIANGLE_NUMBER<lenght, 3u>{	public: static const unsigned int ans = (((lenght % 2) == 0 ? 2 : 1) * TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght,3u>::ans * TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght+1,3u>::ans*TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght+2,3u>::ans);	};
template <unsigned int lenght> class TEMPLATE_TRIANGLE_NUMBER<lenght, 4u>{	public: static const unsigned int ans = (((lenght % 3) == 0 ? 3 : 1) * TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght,4u>::ans * TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght+1,4u>::ans*TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght+2,4u>::ans*TEMPLATE_DIVIDE_BY_MAX_DIVISOR<lenght+3,4u>::ans);	};
template <unsigned int lenght> class TEMPLATE_TRIANGLE_NUMBER<lenght, 5u>{	public: static const unsigned int ans = (lenght * (lenght +1)*(lenght +2)*(lenght +3)*(lenght +4)/120);	};


template <int base, int exp> class TEMPLATE_INT_POWER{
public: enum {ans = base * TEMPLATE_INT_POWER<base, exp - 1>::ans};};

template <int base> class TEMPLATE_INT_POWER<base,0>{
public:	enum {ans = 1};};

template <int val, class C, class D> class MT_IFTYPE{public:typedef C TYPE;};
template <class C, class D> class MT_IFTYPE<0,C,D>{public: typedef D TYPE;};

template <class C> class TEMP_IS_CONST{ public: static const bool ans = false;  };
template <class C> class TEMP_IS_CONST<const C>{ public: static const bool ans = true;  };

template <class C> class TEMP_IS_REF{ public: static const bool ans = false; };
template <class C> class TEMP_IS_REF<C&>{ public: static const bool ans = true;  };

class ExOp;

// convention
template<char C> class LTRType{typedef void TYPE;};
template< > class LTRType<'c'>{typedef char TYPE;};
template< > class LTRType<'C'>{typedef unsigned char TYPE;};
template< > class LTRType<'s'>{typedef StConstList<16>::SINT TYPE;};
template< > class LTRType<'S'>{typedef StConstList<16>::UINT TYPE;};
template< > class LTRType<'l'>{typedef StConstList<64>::SINT TYPE;};
template< > class LTRType<'L'>{typedef StConstList<64>::UINT TYPE;};
template< > class LTRType<'i'>{typedef StConstList<32>::SINT TYPE;};
template< > class LTRType<'I'>{typedef StConstList<32>::UINT TYPE;};
template< > class LTRType<'f'>{typedef float TYPE;};
template< > class LTRType<'d'>{typedef double TYPE;};

template<class A> class IsLFHPrimitive{public: enum{ans = 0};};

template <class A, class B> class isTypeEquivalent {enum {ans = false };

};

template<class lhs, class rhs> struct isTypeEqual {enum {ans = false }; };
template<class lhs> struct isTypeEqual<lhs, lhs> {enum { ans = true };  };

template<class A, class B> class STDRETTYPE2;

	template<class A, class B>
	class STDRETTYPE2{
		public:
		typedef A DEF_TYPE;
		typedef A PLUS_TYPE; // operator+
		typedef A MINU_TYPE; // operator-
		typedef A PROD_TYPE; // operator*
		typedef A DIVI_TYPE; // operator/
		typedef A OP2_TYPE; // operator()

        typedef A RMUL_TYPE; // mult_rigth  a *bH
		typedef A LMUL_TYPE; // mult_left   bH * a
        typedef A RDIV_TYPE; // mult_rigth  a / bH
		typedef A LDIV_TYPE; // mult_left   bH \ a
	};



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// End of Meta-programattion
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Class declarations
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	class ProgressBarPrint;
	class Function;

	template<class KEY, class DATA, bool isMax = true> class ExtremumScope;

    class emptyclass;

	template<typename A> class EnumBox;


	template<class A> class classarg : public A{};
	template<bool> class WarH; // Handdles warnings!
	template<bool> class LinkAssert; // if false, gets and error

	class Weight; // simple double, may have special meaning

	template<class C, unsigned int nbstates> class Classifier;

	template<int charsize, class T> class Anything;

	template<int buffersize> class SuperString;



	class Bluestein;



	class TiffFile;
	template<class C> class TiffFileTypeInterpretation;

	//class PriorityQueue;

	template<class C> class Complex;
	template<class C> class Quaternion;



	template<class A> class Oper1;
	template<class A, class B> class Oper2;
	template<class A, class B, FLAGFORTYPE AT, FLAGFORTYPE BT> class OperMK2;
	template<class A, class B, class C> class Oper3;

	template<class C, bool hasone> class mantissa;

	template<class C> class angle;

	template<class I, int ISIZE, class O, int OSIZE> class Continuous;

	enum Tuple_flag{
		TUPLE_FLAG_NULL =0, // default
		TUPLE_FLAG_COMPLEX,
		TUPLE_FLAG_HYPERPOSITION,
		TUPLE_FLAG_WEIGHT_ORDER_0= 256,
		TUPLE_FLAG_WEIGHT_ORDER_1= 257,
		TUPLE_FLAG_WEIGHT_ORDER_2= 258,
		TUPLE_FLAG_WEIGHT_ORDER_3= 259,
		TUPLE_FLAG_WEIGHT_ORDER_4= 260,
		TUPLE_FLAG_WEIGHT_ORDER_5= 261,
		TUPLE_FLAG_WEIGHT_ORDER_6= 262,
		TUPLE_FLAG_WEIGHT_ORDER_7= 263,
        TUPLE_FLAG_REMOTE_MEMORY=1024
	};

	template<class C, unsigned int size =0u, Tuple_flag flag= TUPLE_FLAG_NULL> class Tuple;
	template<class C, unsigned int order> class WeightElem_baseclass;
	template<class C, unsigned int order> class WeightElem;
	template<class C, unsigned int flag = 0u> class GaussElem; // order 2 with "covariances"

	template<class C, class V> class GaussScope;
		template <class C, int size, Tuple_flag Cflag> class PartitionTupleComparator;


	template<class C, unsigned int sizex, unsigned int sizey> class TMatrix;
	template<class C> class Matrix;

	template<class C, int size> class Matrianglix; // squarre Matrix stored in triangular format!

enum LFHVectorflag{
	LFHVECTOR_NORMAL=0,
	LFHVECTOR_REMOTE=1, // data is not owned! cant change size or deallocate!
	LFHVECTOR_OWNED_LINKMEMORY=2, // data is owned and the Vector must maintain a scope access for its elements, scope is owned so it must be deleted uppon destruction
	LFHVECTOR_LINKMEMORY=3 // data is owned and the Vector must maintain a scope access for its elements
};
	template<class C> class Vector;
	template <class C> class HeapTree;

	// abstract class, read_only
	template<class C, unsigned int nbdim> class ConstGrid;
	template<class C, unsigned int nbdim = 2u> class DataGrid;
	template<class C, class D, class FUNC, unsigned int nbdim> class FuncGrid;

	class SetComparison;
	template<class C> class IntervalSet;
	template<class C> class Interval;

	template<unsigned int nbdim> class GaussianDistribution;
	template<class C, int nbrel> class Forest; // 1 = parent, 2 = LR, 3= LRP, 4 = LRBA, 5 = LRBAP


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// End of Class declarations
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


// funky functions - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//	static void showDoubleMap(char* path, double* data, int resolution, int sizex, int sizey, int nbchannels);
//	static double quadsolve(double* coefs); // return smallest positive root, returns a negative value if none exists
//  static double orientedquadsolve(double* coefs); // return smallest positive root, where the derivative is negative, returns a negative value if none exists
//  static double negorientedquadsolve(double* coefs); // return smallest positive root, where the derivative is negative, returns a negative value if none exists
//	static void crossProduct(double * const out,double const * const vec,double const * const vec2);

	template <class A,class B> class AbstractKeyIterator{
		public:
		A curkey;
		const B& target;
		AbstractKeyIterator(const B& _target) : target(_target){}
		operator const A& () const{return(curkey);}
		const A& operator()()const{return(curkey);}
		virtual bool first()=0;
		virtual bool last()=0;
		virtual bool next()=0;
		virtual bool prev()=0;
	};

	template <class A,class B> class Iterator{
	public:
		Iterator(){}
		A* first(B & ob){init(ob); return(next(ob));}
		virtual void init( B &)=0;
		virtual A* next( B &)=0;
	};

	enum functionname{
		FUNCNAME_NULL =0,
		FUNCNAME_SQUARRE,
		FUNCNAME_INVERSE,
		FUNCNAME_SIN,
		FUNCNAME_COS,
		FUNCNAME_LNGAMMA,
		FUNCNAME_POLYGAMMA_0,
		FUNCNAME_SINH,
		FUNCNAME_COSH,
		FUNCNAME_BESSEL_J0,
		FUNCNAME_BESSEL_J1,
		FUNCNAME_BESSEL_I0,
		FUNCNAME_BESSEL_I1
	};

	template <class K, class D> class KeyElem;

    template <class Key> class defaultHashFnc;
    template <class Key, class Data, class HashFnc = defaultHashFnc<Key> > class myHashmap;
    template <class Key, class BackKey, class Data = emptyclass, class HashFnc = defaultHashFnc<Key> , class BackHashFnc = defaultHashFnc<BackKey> > class dualHashmap;

#define ExCoMeMdEcLaRe(ClAsS_A) void save(FILE* f) const; void load(FILE* f, unsigned int arr_size = 0); void show(FILE* f = stdout, int level = 0) const; string type_tostring()const; void zero();  void one(); void random();

LFH_USE_MEMLEAK_HELP(LFHCONCAT2(extern myHashmap<unsigned int, unsigned int> memleak_helper;))

enum setcomparison{
	SETCMP_EQUAL=0,
	SETCMP_MASKED_GE=1, // true (min(X) >= max(Y))
	SETCMP_MASKED_LE=2, // true (max(X) <= min(Y))
	SETCMP_CMP_MASK=3, // or is undefined!
	SETCMP_DISJOINT=4, // true if A x A y, x != y
	SETCMP_MASKED_GT=5, // true (min(X) >= max(Y))
	SETCMP_MASKED_LT=6, // true (max(X) <= min(Y))
	SETCMP_CMP_MASK2=7, // or is undefined!

	SETCMP_GE= 24 +32 +128 +1, // is there is 1 elem,
	SETCMP_LE= 24 +64 +256 +2,

	SETCMP_GT=24 +32 +128 +5, // true if A x A y, x != y
	SETCMP_LT=24 +64 +256 +6, // true if A x A y, x != y

	SETCMP_IS_NOT_SUBSET=8, // true if E x st A y, x != y
	SETCMP_IS_NOT_SUPERSET=16, // true if E y st A x, x != y

	SETCMP_MAX_BIGGER=32, // true if max(X) >= max(Y)
	SETCMP_MAX_SMALLER=64, // true if max(X) <= max(Y)
	SETCMP_MIN_BIGGER=128, // true if min(X) >= min(Y)
	SETCMP_MIN_SMALLER=256, // true if min(X) <= min(Y)

	SETCMP_SUP_BOUNDED= 32 | 256,
	SETCMP_SUB_BOUNDED= 64 | 128,
};

template<class A,unsigned int flag = 0> class ExCo; // (std::is_enum<A>::value) ? 1 : 0


template<class C, class V>
class GaussScope{ // the weight of vector under a weight
public:
	double w,w2;
    C mean;
	V var;

    GaussScope(): w(0),w2(0){}
    GaussScope(const C &pos, const V &_var, double _w = 1.0f, double _w2 = 0.0f):w(_w),w2(_w2){mean = pos * _w;var = var * _w; if (w2 == 0.0f) w2 = _w*_w; }

	GaussScope<C,V> operator+(const GaussScope<C,V>& other)const {return GaussScope<C,V>(mean+ other.mean, var+ other.var, w+other.w,w2+other.w2);}
	GaussScope<C,V> operator-(const GaussScope<C,V>& other)const {return GaussScope<C,V>(mean- other.mean, var+ other.var, w+other.w,w2+other.w2);}
	GaussScope<C,V>& operator+=(const GaussScope<C,V>& other){
		w += other.w;
        w2 += other.w2;
		mean += other.mean;
        var += other.var;
		return(*this);
	}

	GaussScope<C,V>& operator-=(const GaussScope<C,V>& other){ // TODO
		w += other.w;
        w2 += other.w2;
		mean += other.mean;
        var += other.var;
		return(*this);
	}

	C getMean() const{return( mean / w);}
	GaussScope<C,V>& toZero();
};

#include "primitive_exop.hpp"

template <class A, class B> struct IteratorScope : public Iterator<A,B> {
	enum {valid = false};
};

	template<class KEY, class DATA, bool isMax>
	class ExtremumScope{
	public:
		KEY best_key;
		DATA best;
		inline ExtremumScope<KEY,DATA,isMax>& toZero();
		inline void init(const KEY& key, const DATA& data);
		inline void regist(const KEY& key, const DATA& data);
	};

	template<class KEY, class DATA>
	class ExtremumScope<KEY,DATA,false>{
	public:
		KEY best_key;
		DATA best;
		inline ExtremumScope<KEY,DATA,false>&  toZero();
		inline void init(const KEY& key, const DATA& data);
		inline void regist(const KEY& key, const DATA& data);
	};


	class SetComparison{
	public:
		int comp;
		SetComparison();
		SetComparison(int v);

		template<class O> SetComparison(const O& , const O& , setcomparison mask = SETCMP_CMP_MASK); // MASK FOR DESIRED INFO

		bool areMonotonicDisjoint() {return((comp & 3) != 0);}
		bool areDisjoint() {return((comp & 16) != 0);}
		bool contains() {return((comp & 12) == 4);}
		bool iscontained() {return((comp & 12) == 8);}


		int operator&(const int& mask){ comp &= mask; return(comp);}
		int operator&(const setcomparison& mask){ comp &= (int)mask;return(comp);}


		void show(FILE* f= stdout) const{
			switch(comp){
				case 0: fprintf(f,"Set A = B\n"); break;
				case 1: fprintf(f,"Set A contains B\n"); break;
				case 2: fprintf(f,"Set B contains A\n"); break;
				case 3: fprintf(f,"Set are disjoint\n"); break;
				case 7: fprintf(f,"Set are disjoint and A>=B\n"); break;
				case 11: fprintf(f,"Set are disjoint and B>=A\n"); break;
			}
		};
	};

	template <class Node, int base_resol, int incr_resol, int nb_resol, class intclass, int nbdim, int flag> class MultiResolHashTable;




//	static void loadBMP(const char* const path,DataGrid<Tuple<unsigned char, 3>,2> &im);
//	static void saveBMP(const char* const path,const DataGrid<Tuple<unsigned char, 3>,2> &im);

	template<class prec, unsigned int nbchan, Tuple_flag flag> static void loadWAV(const char* path,DataGrid<Tuple<prec, nbchan>,1> &im);
	template<class prec, unsigned int nbchan, Tuple_flag flag> static void saveWAV(const char* path,DataGrid<Tuple<prec, nbchan>,1> &im);

 //   template<class TYPE> static void FFT_pow2(TYPE* array, unsigned char mag, bool inverse = false);

	template<class TYPE> void pow2_BLUE_FFT_routine(TYPE* data, Complex<double>* blue, unsigned char mag, unsigned int subsize); //  1 <= mult <= 7
	template<class TYPE> void pow2_BLUE_IFFT_routine(TYPE* data, Complex<double>* blue, unsigned char mag, unsigned int subsize); //  1 <= mult <= 7


	char* cloneString(const char* const what);
	template<class TYPE> void pow2_FFT_routine(TYPE* data, unsigned char mag); //output may need swaping
	template<class TYPE> void pow2_IFFT_routine(TYPE* data, unsigned char mag); //output may need swaping
	template<class TYPE> void pow2_IFFT_routine_swapped(TYPE* data, unsigned char mag); // assumes input is swaped, output does not need swaping
	template<class TYPE> void pow2_IFFT_routine_swapped_shift(TYPE* data, unsigned char mag, unsigned int size); // assumes input is swaped, output does not need swaping, only writes the size value, shifted to fit at the end of the array
	template<class TYPE> void comp_FFT_routine(TYPE* data, unsigned char mag, unsigned int mult);
    template<class TYPE> void comp_IFFT_routine(TYPE* data, unsigned char mag, unsigned int mult);
    template<class TYPE> static void pow2_bitswap_permutation(TYPE* array, unsigned char mag);

	// special fouriers
	template<class TYPE> void pow2_FFT_routine_zerohalf(TYPE* data, unsigned char mag); // ignores the 2nd half, assumes they are zeroes



	template<class prec, int nbchan, Tuple_flag flag> static void rescale(DataGrid<Tuple<prec, nbchan>,1> &out,DataGrid<Tuple<prec, nbchan>,1> &sour, double freq);
	template<class prec, int nbchan, Tuple_flag flag> static void lowpassfilter(DataGrid<Tuple<prec, nbchan>,1> &out,DataGrid<Tuple<prec, nbchan>,1> &sour, double freq);
	template<class prec, int nbchan, Tuple_flag flag> static void highpassfilter(DataGrid<Tuple<prec, nbchan>,1> &out,DataGrid<Tuple<prec, nbchan>,1> &sour, double freq);
	template<class prec, int nbchan, Tuple_flag flag> static void lowandhighpassfilter(DataGrid<Tuple<prec, nbchan>,1> &out,DataGrid<Tuple<prec, nbchan>,1> &sour, double low, double high);

	template<class prec, int nbchan, Tuple_flag flag> static void fourrierfilter(DataGrid<Tuple<prec, nbchan>,1> &out,DataGrid<Tuple<prec, nbchan>,1> &sour, double freq);



	template<class prec, int nbchan, Tuple_flag flag> static void extractFreqfromWAV(DataGrid<Tuple<prec, nbchan>,1> &out,DataGrid<Tuple<prec, nbchan>,1> &sour, double freq);


	class ArgumentParser;
	template<int TEMPVAL>
	class templateInt{
		static const int value = TEMPVAL;
	};

	extern WarH<LFH_WARNINGS> static_warning_handdle;
	extern map< void* , pair< unsigned int, unsigned int > > man_in_the_middle;
  //  extern AliasBank static_alias_bank;


	enum ODEfunctions{
	LFHP_ODE_NONE,
	LFHP_ODE_2ND_COS, // k*f(x) * cos = f"(x);
	LFHP_ODE_2ND_SIN, // k*f(x) * sin = f"(x);
	LFHP_ODE_2ND_CYLINDER //  f"(x) = k*( sin(f(x)) * u * cos(f(x)));
	};


	enum LFHwarnings{
		LFH_WARNING_NONE =0,
		LFH_WARNING_NAN =1,
		LFH_WARNING_MAXITE =2,
		LFH_WARNING_UNEXPECTED_INPUT=3,
		LFH_WARNING_MATRIX_SIZES_MISMATCH=4,
		LFH_ERROR_CANT_OPEN_FILE=5,
		LFH_ERROR_CANT_ALLOCATE_MEMORY=6
	};
template<bool a>
class WarH{
public:

	int makeWarning(char*)const{return(0);}
	void report(char*)const{}
	void Warning(int) const{}
	void operator<<(LFHwarnings what)const{}
	void operator()(bool test, LFHwarnings what){}
};

template< >
class WarH<true>{
	int lastw;
	int wcount;
	map<int, char*> warns;
	void reportlast(){
		if (lastw == 0) return;
		if (wcount > 1) fprintf(stderr,"%ix:", wcount);

		if (lastw < 10){
			switch((LFHwarnings)lastw){
				case LFH_WARNING_NAN: fprintf(stderr,"NaN obtained!\n"); break;
				case LFH_WARNING_MAXITE: fprintf(stderr,"Max Iteration Reached!\n"); break;
				case LFH_WARNING_UNEXPECTED_INPUT:  fprintf(stderr,"Function received unexpected input!\n"); break;
				case LFH_ERROR_CANT_OPEN_FILE:  fprintf(stderr,"Could not open file! (critical error)\n"); break;
				case LFH_ERROR_CANT_ALLOCATE_MEMORY:  myexit("Could allocate memory! (critical error)\n");
				default:  fprintf(stderr,"Unknown Warning used!\n");
			}
		}else fprintf(stderr,warns[lastw]);
		fflush(stderr);
	}
public:
	WarH(): lastw(0){

	}
	~WarH(){ reportlast();
	}
	int makeWarning(char*)const{
		return(0);
	}
	void report(char* arg)const{fprintf(stderr,arg);}
	void operator<<(LFHwarnings what){this->Warning((int) what);}

	void operator()(bool test, LFHwarnings what){if (test) this->Warning((int) what);}
	void Warning(int w){
		if (lastw == w) wcount++;
		else{
			reportlast();
			lastw = w;
			wcount++;
		}
	}

};

	template<bool a> class LinkAssert{public: char junk[ a ? 1 : -1]; };

	// Polyvalent Reader Writer, which call convertions from what is read into the desired output
	class WiseIOUnit{
	public:
		template<class T> static T getFrom(const char* path);
	};

	// event type, entity array

	template<class A>
	class Oper1{
	public:
		static const bool NBARG = 1;
		static const bool NBCONSTARG = 0;
		typedef A& ARG1;
		typedef void RETT;
		virtual void operator()(A &) const =0;

		template <class CA> void operator()(CA & ca) const{
			if (IteratorScope< A , CA >::valid) {
			IteratorScope< A , CA > ite;
			A* pt;
			for( pt = ite.first(ca);pt != NULL;pt = ite.next(ca)) (*this)(*pt);
			} else fprintf(stderr,"No valid Match for operator!\n");
		}

		template <class CA> void apply(CA & ca) const{
			if (IteratorScope< A , CA >::valid) {
				IteratorScope< A , CA > ite;
				A* pt;
				for( pt = ite.first(ca);pt != NULL;pt = ite.next(ca)) (*this)(*pt);
			} else fprintf(stderr,"No valid Match for operator!\n");
		}
	};





	template<class A, class B>
	class Oper2{
	public:
		static const bool IsPOD = false;
		static const bool NeedsAddLink = false; // containers needs to update addresses in link registers

		typedef A ARG1;
		typedef B ARG2;



		virtual void operator()(A& , B& ) const =0;
		template <class CA,class CB> void operator()(CA & ca, CB & cb) const{ ca(*this,cb);}
		template <class CA,class CB> void operator()(CA & ca, const CB & cb) const{ ca(*this,cb);}


		template <class CA> void apply(CA & ca , const B & b) const{
			if (IteratorScope< A , CA >::valid) {
				IteratorScope< A , CA > ite;
				A* pt;
				for( pt = ite.first(ca);pt != NULL;pt = ite.next(ca)) (*this)(*pt, b);
			} else fprintf(stderr,"No valid Match for operator!\n");
		}

		template <class CA,class CB> void apply(CA & ca , CB & cb) const{ // assumes containners are of the same size!
			if ((IteratorScope< A , CA >::valid)&&(IteratorScope< B , CB >::valid)) {
				IteratorScope< A , CA > itea;
				IteratorScope< B , CB > iteb;
				A* pa; B* pb;
				for( pa = itea.first(ca) , pb = iteb.first(cb) ;pa != NULL;pa = itea.next(ca),pb = iteb.next(cb)) (*this)(*pa, *pb);
			} else fprintf(stderr,"No valid Match for operator!\n");
		}
		template <class CA,class CB> void apply(CA & ca , const CB & cb) const{ // assumes containners are of the same size!
			if ((IteratorScope< A , CA >::valid)&&(IteratorScope< B , CB >::valid)) {
				IteratorScope< A , CA > itea;
				IteratorScope< B , CB > iteb;
				A* pa; B* pb;
				for( pa = itea.first(ca) , pb = iteb.first(cb) ;pa != NULL;pa = itea.next(ca),pb = iteb.next(cb)) (*this)(*pa, *pb);
			} else fprintf(stderr,"No valid Match for operator!\n");
		}
	};

	template<class A, class B, class C>
	class Oper3{
	public:
		virtual void operator()(A &, B &, C &) const =0;
		template <class CA,class CB,class CC> void operator()(CA & ca, CB & cb, CC & cc) const{ ca(*this,cb,cc);}
		template <class CA,class CB,class CC> void operator()(CA & ca, CB & cb, const CC & cc) const{ ca(*this,cb,cc);}
		template <class CA,class CB,class CC> void operator()(CA & ca, const CB & cb, const CC & cc) const{ ca(*this,cb,cc);}
	};

	template<int order> double cheb_eval(const Tuple<double, order> cs, const double x);

	// functions!

	
LFH_GOLD	double lngamma(double x);
#ifdef GNU_SCIENTIFIC_LIBRARY
double log_1plusx_mx_e(const double x); // Log(1 + x) - x
LFH_GOLD	double logdensity_chisquarre(double x, double free);
LFH_GOLD	double Pvalue_chisquarre_Ptail(double x, double free);
LFH_GOLD	double Pvalue_chisquarre_Ntail(double x, double free);
LFH_GOLD	double Pvalue_Gamma_Ptail(double x, double k, double theta);
LFH_GOLD	double Pvalue_Gamma_Ntail(double x, double k, double theta);
LFH_GOLD	double Pvalue_GammaRate_Ptail(double x, double alpha, double beta);
LFH_GOLD	double Pvalue_GammaRate_Ntail(double x, double alpha, double beta);
LFH_GOLD	double Pvalue_SumLogPvalues_Ptail(double sum, int nbpvals);
LFH_GOLD	double Pvalue_Beta_Ntail(double x, double a, double b);
LFH_GOLD	double Pvalue_Beta_Ptail(double x, double a, double b);
LFH_GOLD	double Pvalue_Fdistrib_Ntail(double x, double a, double b);
LFH_GOLD	double Pvalue_Fdistrib_Ptail(double x, double a, double b);

LFH_GOLD	double LogPvalue_chisquarre_Ptail(double x, double free);
LFH_GOLD	double LogPvalue_chisquarre_Ntail(double x, double free);
LFH_GOLD	double LogPvalue_Gamma_Ptail(double x, double k, double theta);
LFH_GOLD	double LogPvalue_Gamma_Ntail(double x, double k, double theta);
LFH_GOLD	double LogPvalue_GammaRate_Ptail(double x, double alpha, double beta);
LFH_GOLD	double LogPvalue_GammaRate_Ntail(double x, double alpha, double beta);
LFH_GOLD	double LogPvalue_SumLogPvalues_Ptail(double sum, int nbpvals);
LFH_GOLD	double LogPvalue_Beta_Ntail(double x, double a, double b);
LFH_GOLD	double LogPvalue_Beta_Ptail(double x, double a, double b);
LFH_GOLD	double LogPvalue_Fdistrib_Ntail(double x, double a, double b);
LFH_GOLD	double LogPvalue_Fdistrib_Ptail(double x, double a, double b);

double BesselJ0(double);
double BesselJ1(double);
double BesselI0(double);
double BesselI1(double);

#endif

	double gammastar(double);
	double incgamma_dom(double, double); // D(a,x) := x^a e^(-x) / Gamma(a+1), domminant part
	double incgamma_frac(double, double); // CDP looking incomplete gamma
	double polygamma0(double);
	double sinc(double);
	double sinc_d(double);

	void addRuler(DataGrid<unsigned char, 3> &RGBimage, double bounds[4], int nbsep[2]);


	double sampleGaussian();

	double distanceToEllipse(double squared_dist1, double squared_dist2, double squared_width, double squarred_focal_dist);
	double CubicRealRoot(double*, bool want_middle); // smallest root returned otherwise
	double CubicInterpolationRoot(double* , bool is_monotonic);

    void RGB2Dvec(unsigned char *, double *); // reads 2 double, write 3 chars

enum struct_property{
	STRUCTPROP_IS_VALID=1, // is special values exists, then this checks for them
	STRUCTPROP_HAS_VALID_INVERSE=2, // check if 1.0f / X is a special value!
	STRUCTPROP_CAN_COMMUTE=4, // check is for all Y, YX = XY
};



void GPParamSearch(const Vector<KeyElem<double, double> >&, double &mean, double &var_signal, double &var_noise, double &scale);

// abstract class

	// treats A as an enum
	template<typename A>
	class EnumBox{
		public:
		static const bool IsPOD = true;
		static const bool NeedsAddLink = false; // containers needs to update addresses in link registers


		A data;
		EnumBox():data((A)0){}
		EnumBox(const A& f_init):data((A)f_init){}

		bool operator>(const EnumBox<A> &o)const { return data > o.data;}
		bool operator>=(const EnumBox<A> &o)const{ return data >= o.data;}
		bool operator<(const EnumBox<A> &o)const{ return data < o.data;}
		bool operator<=(const EnumBox<A> &o)const{ return data <= o.data;}
		bool operator==(const EnumBox<A> &o)const{ return data == o.data;}
		bool operator!=(const EnumBox<A> &o)const{ return data != o.data;}


		};


template <class C, unsigned int DIM>
class ConstGrid{
	public:
		virtual C operator()(const Tuple<unsigned int, DIM> &coor) const =0;
		virtual void getDims(Tuple<unsigned int, DIM> &o_dims) const =0;
	};


// list of the member function, any struc in privitive structures needs!








#define HAS_MEMBER_FUNC(func, name)                                        \
template<typename T,typename Sign>                                 \
struct name {                                                       \
typedef char yes[1];                                            \
typedef char no [2];                                            \
template <typename U, U> struct type_check;                     \
template <typename _1> static yes &chk(type_check<Sign, &_1::func> *); \
template <typename   > static no  &chk(...);                    \
static bool const ans = sizeof(chk<T>(0)) == sizeof(yes);     \
}

HAS_MEMBER_FUNC(memmove, Exlisten_memmove);
HAS_MEMBER_FUNC(getWeight, Exlisten_getWeight);
HAS_MEMBER_FUNC(show, Exlisten_show);

HAS_MEMBER_FUNC(load, Exlisten_load);
HAS_MEMBER_FUNC(save, Exlisten_save);

HAS_MEMBER_FUNC(toZero, Exlisten_toZero);
HAS_MEMBER_FUNC(toRand, Exlisten_toRand);
HAS_MEMBER_FUNC(toOne, Exlisten_toOne);
HAS_MEMBER_FUNC(toMin, Exlisten_toMin);
HAS_MEMBER_FUNC(toMax, Exlisten_toMax);

HAS_MEMBER_FUNC(isZero, Exlisten_isZero);
HAS_MEMBER_FUNC(isOne, Exlisten_isOne);

HAS_MEMBER_FUNC(lognorm, Exlisten_lognorm);
HAS_MEMBER_FUNC(pnorm, Exlisten_pnorm);
HAS_MEMBER_FUNC(pdist, Exlisten_pdist);

HAS_MEMBER_FUNC(tointpow, Exlisten_tointpow);
HAS_MEMBER_FUNC(tosquare, Exlisten_tosquare);
HAS_MEMBER_FUNC(toinverse, Exlisten_toinverse);
HAS_MEMBER_FUNC(tonegative, Exlisten_tonegative);

HAS_MEMBER_FUNC(operator+=, Exlisten_toadd);
HAS_MEMBER_FUNC(operator-=, Exlisten_tosub);
HAS_MEMBER_FUNC(operator*=, Exlisten_tomult);
HAS_MEMBER_FUNC(operator/=, Exlisten_todivi);

HAS_MEMBER_FUNC(operator+, Exlisten_mkadd);
HAS_MEMBER_FUNC(operator-, Exlisten_mksub);
HAS_MEMBER_FUNC(operator*, Exlisten_mkmult);
HAS_MEMBER_FUNC(operator/, Exlisten_mkdivi);

HAS_MEMBER_FUNC(operator<, Exlisten_lt);
HAS_MEMBER_FUNC(operator<=, Exlisten_le);
HAS_MEMBER_FUNC(operator>, Exlisten_gt);
HAS_MEMBER_FUNC(operator>=, Exlisten_ge);
HAS_MEMBER_FUNC(operator==, Exlisten_eq);
HAS_MEMBER_FUNC(operator!=, Exlisten_nq);

HAS_MEMBER_FUNC(mkrealproj, Exlisten_mkrealproj);
HAS_MEMBER_FUNC(mkimmaproj, Exlisten_mkimmaproj);
HAS_MEMBER_FUNC(mkjmmaproj, Exlisten_mkjmmaproj);
HAS_MEMBER_FUNC(mkkmmaproj, Exlisten_mkkmmaproj);


HAS_MEMBER_FUNC(mksquare, Exlisten_mksquare);
HAS_MEMBER_FUNC(mkinverse, Exlisten_mkinverse);
HAS_MEMBER_FUNC(mkintpow, Exlisten_mkintpow);
HAS_MEMBER_FUNC(mkinvintpow, Exlisten_mkinvintpow);
HAS_MEMBER_FUNC(mktrju, Exlisten_mktrju); // conjugate transpose
HAS_MEMBER_FUNC(mkbmul, Exlisten_mkbmul); // left multiplication
HAS_MEMBER_FUNC(tobmul, Exlisten_tobmul); // left multiplication assignment

HAS_MEMBER_FUNC(setcmp, Exlisten_setcmp);
HAS_MEMBER_FUNC(operator(), Exlisten_mainoper);
HAS_MEMBER_FUNC(isValid, Exlisten_isValid);

HAS_MEMBER_FUNC(wrFirstIterator, Exlisten_wrFirstIterator);
HAS_MEMBER_FUNC(wrLastIterator, Exlisten_wrLastIterator);
HAS_MEMBER_FUNC(wrNextIterator, Exlisten_wrNextIterator);
HAS_MEMBER_FUNC(wrPrevIterator, Exlisten_wrPrevIterator);
HAS_MEMBER_FUNC(wrEndIterator, Exlisten_wrEndIterator);

HAS_MEMBER_FUNC(mkouterprod, Exlisten_mkouterprod); // conjugate transpose
HAS_MEMBER_FUNC(mkvectorization, Exlisten_mkvectorization); // conjugate transpose

HAS_MEMBER_FUNC(type_dimentions, Exlisten_type_dimentions);

HAS_MEMBER_FUNC(mkgaussstat, Exlisten_mkgaussstat);

HAS_MEMBER_FUNC(type_tostring, Exlisten_type_tostring);

enum exo{
    EXOP_zero,
};

template<exo EXO, class NEXT>
class ExCm{

};

template<bool hasfunc>
class ExFn{// ExFn<true>
	public:
		double dummy();

		template<class A> inline static bool isZero(const A& what) {return what.isZero();}
		template<class A> inline static bool isOne(const A& what) {return what.isOne();}
		template<class A> inline static A& toZero(A& what) {return what.toZero();}
		template<class A> inline static A& toOne(A& what) {return what.toOne();}
		template<class A> inline static A& toRand(A& what) {return what.toRand();}
		template<class A> inline static A& toMin(A& what) {return what.toMin();}
		template<class A> inline static A& toMax(A& what) {return what.toMax();}

		// builders

		template<class A> inline static A mknegative(const A& a) {return -a;}
		template<class A> inline static A mkinverse(const A& a) {return a.mkinverse();}

		template<class A> inline static A mkintpow(const A& a, const int pow) {return a.mkintpow(pow);}
		template<class A> inline static A mkinvintpow(const A& a, const int pow) {return a.mkinvintpow(pow);}

		template<class A> inline static A mksquare(const A& a) {return a.mksquare();}
        template<class A> inline static typename ExCo<A>::TRJU_TYPE mktrju(const A& a){return a.mktrju();}

		// assigns

		template<class A> inline static A& tonegative(A& a) {return a.tonegative();}
		template<class A> inline static A& toinverse(A& a) {return a.toinverse();}

		template<class A> inline static void tointpow(A& what, const int pow) {what.tointpow(pow);}
		template<class A> inline static A& tosquare(A& a) {return ExCo<A>::tosquare(a);}


	template<class A> inline static	A& toadd(A& a, const A& b) {return (a += b);}
	template<class A> inline static A& tosub(A& a, const A& b) {return (a -= b);}
	template<class A> inline static A& tomult(A& a, const A& b) {return (a *= b);}
	template<class A> inline static A& todivi(A& a, const A& b) {return (a /= b);}
	template<class A, class B> inline static A& toadd(A& a, const B& b) {return (a += b);}
	template<class A, class B> inline static A& tosub(A& a, const B& b) {return (a -= b);}
	template<class A, class B> inline static A& tomult(A& a, const B& b) {return (a *= b);}
	template<class A, class B> inline static A& todivi(A& a, const B& b) {return (a /= b);}

	template<class A> inline static	A mkadd(const A& a, const A& b) {return a + b;}
	template<class A> inline static A mksub(const A& a, const A& b) {return a - b;}
	template<class A> inline static A mkmult(const A& a, const A& b) {return a * b;}
	template<class A> inline static A mkdivi(const A& a, const A& b) {return a / b;}
	template<class A, class B> inline static typename STDRETTYPE2<A,B>::PLUS_TYPE mkadd(const A& a, const B& b) {return a + b;}
	template<class A, class B> inline static typename STDRETTYPE2<A,B>::MINU_TYPE mksub(const A& a, const B& b) {return a - b;}
	template<class A, class B> inline static typename STDRETTYPE2<A,B>::PROD_TYPE mkmult(const A& a, const B& b) {return a * b;}
	template<class A, class B> inline static typename STDRETTYPE2<A,B>::DIVI_TYPE mkdivi(const A& a, const B& b) {return a / b;}
	
	template<class A, class B> inline static typename STDRETTYPE2<A,B>::PROD_TYPE mkmult_trju(const A& a, const B& b) {return b * a;}

	template<class A, class B> inline static A mkadd_alt(const A& a, const B& b) {A ca=a; return ca += b;}
	template<class A, class B> inline static A mksub_alt(const A& a, const B& b) {A ca=a; return ca -= b;}
	template<class A, class B> inline static A mkmult_alt(const A& a, const B& b) {A ca=a; return ca *= b;}
	template<class A, class B> inline static A mkdivi_alt(const A& a, const B& b) {A ca=a; return ca /= b;}
	
	
	

	template<class A> inline static typename ExCo<A>::REAL_TYPE mkrealproj(const A& a){return a.mkrealproj();}
	template<class A> inline static typename ExCo<A>::REAL_TYPE mkimmaproj(const A& a){return a.mkimmaproj();}
	template<class A> inline static typename ExCo<A>::REAL_TYPE mkjmmaproj(const A& a){return a.mkjmmaproj();}
	template<class A> inline static typename ExCo<A>::REAL_TYPE mkkmmaproj(const A& a){return a.mkkmmaproj();}

    template<class A> inline static typename ExCo<A>::GAUS_TYPE mkgaussstat(const A& a,double & weight){return a.mkgaussstat( weight);}

	template<class A, class B> inline static const A& tobmul(A& a, const B& b) {return (a.tobmul(b));}
	template<class A, class B> inline static A mkbmul(const A& a, const B& b) {return (a.mkbmul(b));}
	
	
	template<class A> inline static typename ExCo<A>::OUTER_TYPE mkouterprod(const A& a, const A& b){return a.mkouterprod(b);}
	

    template<class A, class B> inline static bool comp_lt(const A& a, const B& b) {return a < b;}
    template<class A, class B> inline static bool comp_le(const A& a, const B& b) {return a <= b;}
    template<class A, class B> inline static bool comp_gt(const A& a, const B& b) {return a > b;}
    template<class A, class B> inline static bool comp_ge(const A& a, const B& b) {return a >= b;}
    template<class A, class B> inline static bool comp_eq(const A& a, const B& b) {return a == b;}
    template<class A, class B> inline static bool comp_nq(const A& a, const B& b) {return a != b;}

    template<class A, class B> inline static bool comp_lti(const A& a, const B& b) {return b > a;}
    template<class A, class B> inline static bool comp_lei(const A& a, const B& b) {return b >= a;}
    template<class A, class B> inline static bool comp_gti(const A& a, const B& b) {return b < a;}
    template<class A, class B> inline static bool comp_gei(const A& a, const B& b) {return b <= a;}
    template<class A, class B> inline static bool comp_eqi(const A& a, const B& b) {return b == a;}
    template<class A, class B> inline static bool comp_nqi(const A& a, const B& b) {return b != a;}



		template<class F, class I> static typename ExCo<F>::template RETT<I>::TYPE comp3(F f, I f_in) {return f(f_in);}


		template<class A> inline static setcomparison setcmp(const A &a, const A &b) {return a.setcmp(b);}

		template<class A> inline static void save(const A& what, FILE *f) {what.save(f);}
		template<class A> inline static void load(A& what, FILE *f, unsigned int lenght) {what.load(f, lenght);}
		template<class A> inline static void save_ISA_pod(const A& what, FILE *f) {fwrite(&what,sizeof(A),1,f);}
		template<class A> inline static void load_ISA_pod(A& what, FILE *f, unsigned int lenght) { fread(&what,sizeof(A),1,f);}
	
		template<class A> inline static A& memmove(A& a, A& o){return a.memmove(o);}
		template<class A> inline static double getWeight(const A& a){return a.getWeight; }


		template<class A> inline static double pdist(const A& a, const A& b){return a.pdist(b);}
		template<class A> inline static double pnorm(const A& a){return a.pnorm(); }
		template<class A> inline static double lognorm(const A& a){return a.lognorm(); }
		template<class A> inline static void show(const A& a, FILE*f , int level){return a.show(f,level);}
        template<class A> inline static string type_tostring(const A& a){return  a.type_tostring();}


		template<class F,exo EXO, class A> inline static void compose(const ExCm<EXO,F> &func, A &a);

		template<class F, class A> inline static void compose(F &func, A &a){return func(a);}
		template<class F, class A> inline static void compose(F &func, const A &a){return func(a);}
		template<class F, class A> inline static void compose(const F &func, A &a){return func(a);}

		template<class F, class A, class B> inline static void compose(F &func, A &a, B &b){return func(a,b);}
		template<class F, class A, class B> inline static void compose(F &func, A &a, const B &b){return func(a,b);}
		template<class F, class A, class B> inline static void compose(F &func, const A &a, const B &b){return func(a,b);}
        template<class F, class A, class B> inline static void compose(const F &func, A &a, const B &b){return func(a,b);}
        template<class F, class A, class B> inline static void compose(const F &func, A &a, B &b){return func(a,b);}

		template<class A> inline static bool isValid(const A &a){return a.isValid();}
		template<class A> inline static bool isValid_cmpchk(const A &a){return ExCo<A>::isValid(a);}

	};

template< >
class ExFn<false>{
	public:
		
        template<class A> inline static bool isZero(const A& what) {return ExCo<A>::isZero(what);}
		template<class A> inline static bool isOne(const A& what) {return ExCo<A>::isOne(what);}
        template<class A> inline static A& toZero(A& what) {ExCo<A>::zero(what); return(what);}
		template<class A> inline static A& toOne(A& what) {ExCo<A>::one(what); return what;}
		template<class A> inline static A& toRand(A& what) {ExCo<A>::random(what); return(what);}
		template<class A> inline static A& toMin(A& what) {ExCo<A>::minimum(what);return what;}
		template<class A> inline static A& toMax(A& what) {ExCo<A>::maximum(what);return what;}



		// fundamental builders

		// builders

		template<class A> inline static A mknegative(const A& a) {return ExCo<A>::mknegative(a);}
		template<class A> inline static A mkinverse(const A& a) {return ExCo<A>::mkinverse(a);}
		template<class A> inline static A mksquare(const A& a) {return ExCo<A>::mksquare(a);}
		template<class A> inline static A mkintpow(const A& what, const int pow) {return ExCo<A>::mkintpow(what,pow);}
		template<class A> inline static A mkinvintpow(const A& what, const int pow) {return ExCo<A>::mkinvintpow(what,pow);}
        template<class A> inline static typename ExCo<A>::TRJU_TYPE mktrju(const A& a){return A(a);}

		// assigns

		template<class A> inline static A& tonegative(A& a) {return ExCo<A>::tonegative(a);}
		template<class A> inline static A& toinverse(A& a) {return ExCo<A>::toinverse(a);}
		template<class A> inline static A& tosquare(A& a) {return ExCo<A>::tosquare(a);}
		template<class A> inline static void tointpow(A& what, const int pow) {ExCo<A>::tointpow(what,pow);}

		template<class A> inline static	A& toadd(A& a, const A& b) {return ExCo<A>::toadd(a,b);}
		template<class A> inline static A& tosub(A& a, const A& b) {return ExCo<A>::tosub(a,b);}
		template<class A> inline static A& tomult(A& a, const A& b) {return ExCo<A>::tomult(a,b);}
		template<class A> inline static A& todivi(A& a, const A& b) {return ExCo<A>::todivi(a,b);}
		template<class A,class B> inline static	A& toadd(A& a, const B& b) {return ExCo<A>::toadd(a,b);}
		template<class A,class B> inline static A& tosub(A& a, const B& b) {return ExCo<A>::tosub(a,b);}
		template<class A,class B> inline static A& tomult(A& a, const B& b) {return ExCo<A>::tomult(a,b);}
		template<class A,class B> inline static A& todivi(A& a, const B& b) {return ExCo<A>::todivi(a,b);}


	//	template<class A, class B> inline static typename STDRETTYPE2<A,B>::PLUS_TYPE mkadd(const A& a, const B& b) {return ExCo<A>::mkadd(a,b);}
	//	template<class A, class B> inline static typename STDRETTYPE2<A,B>::MINU_TYPE mksub(const A& a, const B& b) {return ExCo<A>::mksub(a,b);}
	//	template<class A, class B> inline static typename STDRETTYPE2<A,B>::PROD_TYPE mkmult(const A& a, const B& b) {return ExCo<A>::mkmult(a,b); }
	//	template<class A, class B> inline static typename STDRETTYPE2<A,B>::DIVI_TYPE mkdivi(const A& a, const B& b) {return ExCo<A>::mkdivi(a,b);}
		template<class A> inline static	A mkadd(const A& a, const A& b)  {return ExCo<A>::mkadd(a,b);}
		template<class A> inline static A mksub(const A& a, const A& b)  {return ExCo<A>::mksub(a,b);}
	
//		template<class A> inline static A mkmult(const A& a, const A& b) {return ExFn< Exlisten_mktrju<decltype(b * a), decltype(b * a)& (ExCo<decltype(b * a)>:  >::mkmult_trju(a,b);}
	
	//	template<class A, class B> inline static STDRETTYPE2<B,A>::PROD_TYPE mkmult(const A& a, const B& b) {return ExFn< Exlisten_mkmult<B, typename STDRETTYPE2<B,A>::PROD_TYPE (ExCo<B>::SAFETYPE::*)(const A& a)const>  >::mkmult_mktrju(a,b);}
		template<class A> inline static A mkmult(const A& a, const A& b) {return ExCo<A>::mkmult(a,b);}
	
	
	
		template<class A> inline static A mkdivi(const A& a, const A& b) {return ExCo<A>::mkdivi(a,b);}
		template<class A, class B> inline static A mkadd(const A& a, const B& b)  {return ExFn< Exlisten_toadd< A , A& (ExCo<A>::SAFETYPE::*)(const B&) >::ans >::mkadd_alt(a,b);}
		template<class A, class B> inline static A mksub(const A& a, const B& b)  {return ExFn< Exlisten_tosub< A , A& (ExCo<A>::SAFETYPE::*)(const B&) >::ans >::mksub_alt(a,b);}
		template<class A, class B> inline static A mkmult(const A& a, const B& b) {return ExFn< Exlisten_tomult< A , A& (ExCo<A>::SAFETYPE::*)(const B&) >::ans >::mkmult_alt(a,b);}
		template<class A, class B> inline static A mkdivi(const A& a, const B& b) {return ExFn< Exlisten_todivi< A , A& (ExCo<A>::SAFETYPE::*)(const B&) >::ans >::mkdivi_alt(a,b);}

		template<class A> inline static typename ExCo<A>::OUTER_TYPE mkouterprod(const A& a, const A& b){return ExFn< Exlisten_mkouterprod< A , typename ExCo<A>::OUTER_TYPE (ExCo<A>::SAFETYPE::*)(const A&) const >::ans  >::mkouterprod(a,b);}
	

		template<class A, class B> inline static const A& tobmul(A& a, const B& b);
		template<class A, class B> inline static A mkbmul(const A& a, const B& b);

		template<class A> inline static setcomparison setcmp(const A &a, const A &b) {return ExCo<A>::setcmp(a,b);}

	template<class A> inline static void save(const A& what, FILE *f) {ExFn< ExCo<A>::IsPOD >::save_ISA_pod(what,f);}
		template<class A> inline static void load(A& what, FILE *f, unsigned int lenght) { ExFn< ExCo<A>::IsPOD >::load_ISA_pod(what,f,lenght);}
		template<class A> inline static void save_ISA_pod(const A& what, FILE *f) {ExCo<A>::save(what,f);}
		template<class A> inline static void load_ISA_pod(A& what, FILE *f, unsigned int lenght) { ExCo<A>::load(what,f,lenght);}
		template<class A> inline static A& memmove(A& a, A& o){return (a = o);}
		template<class A> inline static double getWeight(const A& a){return 1.0f; }

		template<class A> inline static double pdist(const A& a, const A& b){return ExCo<A>::pdist(a,b);}
		template<class A> inline static double pnorm(const A& a){return ExCo<A>::pnorm(a);}
		template<class A> inline static double lognorm(const A& a){return ExCo<A>::lognorm(a);}

		template<class A> inline static void show(const A& a, FILE*f , int level){return ExCo<A>::show(a,f,level);}
        template<class A> inline static string type_tostring(const A& a){return ExCo<A>::type_tostring(a);}
		template<class F, class A, class B, unsigned int dims> inline static void compose(F &func, DataGrid<A,dims> &a, const DataGrid<B,dims> &b);

		template<class F, class A, class B, int DA, int DB> inline static void compose(F &func, DataGrid<A,DA> &a, const DataGrid<B,DB> &b);


        template<class A, class B> inline static bool comp_lt(const A& a, const B& b) {return( ExFn< Exlisten_gt<B, bool (ExCo<B>::SAFETYPE::*)(const A&)const >::ans >::comp_lti(a,b));}
        template<class A, class B> inline static bool comp_le(const A& a, const B& b) {return( ExFn< Exlisten_ge<B, bool (ExCo<B>::SAFETYPE::*)(const A&)const >::ans >::comp_lei(a,b));}
        template<class A, class B> inline static bool comp_gt(const A& a, const B& b) {return( ExFn< Exlisten_lt<B, bool (ExCo<B>::SAFETYPE::*)(const A&)const >::ans >::comp_gti(a,b));}
        template<class A, class B> inline static bool comp_ge(const A& a, const B& b) {return( ExFn< Exlisten_le<B, bool (ExCo<B>::SAFETYPE::*)(const A&)const >::ans >::comp_gei(a,b));}
        template<class A, class B> inline static bool comp_eq(const A& a, const B& b) {return( ExFn< Exlisten_eq<B, bool (ExCo<B>::SAFETYPE::*)(const A&)const >::ans >::comp_eqi(a,b));}
        template<class A, class B> inline static bool comp_nq(const A& a, const B& b) {return( ExFn< Exlisten_nq<B, bool (ExCo<B>::SAFETYPE::*)(const A&)const >::ans >::comp_nqi(a,b));}

        template<class A, class B> inline static bool comp_lti(const A& a, const B& b) {return b > a;}
        template<class A, class B> inline static bool comp_lei(const A& a, const B& b) {return b >= a;}
        template<class A, class B> inline static bool comp_gti(const A& a, const B& b) {return b < a;}
        template<class A, class B> inline static bool comp_gei(const A& a, const B& b) {return b <= a;}
        template<class A, class B> inline static bool comp_eqi(const A& a, const B& b) {return b == a;}
        template<class A, class B> inline static bool comp_nqi(const A& a, const B& b) {return b != a;}

		template<class A> inline static bool isValid(const A &a){return ExFn< Exlisten_isValid< ExCo<A> , bool (*)(const A&) >::ans >::isValid_cmpchk(a);}
		template<class A> inline static bool isValid_cmpchk(const A &a){return true;}

		template<class F, class I, unsigned int SIZE> static Tuple< typename ExCo<F>::template RETT<I>::TYPE , SIZE > comp3(F f, Tuple<I, SIZE> f_in);
        template<class A> inline static typename ExCo<A>::GAUS_TYPE mkgaussstat(const A& a,double & weight){return ExCo<A>::mkgaussstat(a,weight);}

        template<class A> inline static typename ExCo<A>::REAL_TYPE mkrealproj(const A& a){return ExCo<A>::mkrealproj(a);}
        template<class A> inline static typename ExCo<A>::REAL_TYPE mkimmaproj(const A& a){return ExCo<A>::mkimmaproj(a);}
        template<class A> inline static typename ExCo<A>::REAL_TYPE mkjmmaproj(const A& a){return ExCo<A>::mkjmmaproj(a);}
        template<class A> inline static typename ExCo<A>::REAL_TYPE mkkmmaproj(const A& a){return ExCo<A>::mkkmmaproj(a);}
	};

	class ExOp{
	public:
		// class Properties

		template<class A> static void bitreverse(A& a){ExCo<A>::bitreverse(a);}
		template<class A> static void bytereverse(A& a){ExCo<A>::bytereverse(a);}


		//template<class A> static void next(A& a){ExOp_body<A>::next(a);}



		template<class A> static unsigned char upperbound_pow_of_2(const A& a){return ExCo<A>::upperbound_pow_of_2(a); }

		//template<class F> static typename ExCo<F>::DEFRETT comp2(F f, typename ExCo<F>:: f_in) {return f(f_in);}


		//		template<class F, class I> static auto comp3(F f, I f_in) -> int {return ExFn< Exlisten_mainoper< typename ExCo<F>::SAFETYPE, typename ExCo<F>::template RETT<I>::TYPE (ExCo<F>::SAFETYPE::*)(I b) >::ans  >::comp3(f,f_in);}
				template<class F, class I> static typename ExCo<F>::template RETT<I>::TYPE comp3(F f, I f_in) {return ExFn< Exlisten_mainoper< typename ExCo<F>::SAFETYPE, typename ExCo<F>::template RETT<I>::TYPE (ExCo<F>::SAFETYPE::*)(I b) >::ans  >::comp3(f,f_in);}

		//		template<class F> static int comp3(F f, double f_in) {return ExFn< true  >::comp3(f,f_in);}

		template<class F, class I> static typename ExCo<F>::template RETT<I>::TYPE comp2(F f, I f_in) {return f(f_in);}

		template<class O, class I> static void comp(O& o, void (*fun)(O&, const I&), const I & in){fun(o, in);}
		template<class O, class I> static void comp(O& o, O (*fun)(const I&), const I& in){o = fun(in);}
		template<class O, class I> static void comp(O& o, void (I::*fun)(O&), const I& in){in.fun(o);}
		template<class O, class I> static void comp(O& o, O (I::*fun)(), const I& in){o = in.fun();}
		template<class O, class I, class F> static void comp(O& o, void (F::*fun)(O&, const I&), const I& in){fun(o, in);}
		template<class O, class I, class F> static void comp(O& o, O (F::*fun)(const I&), const I& in){o = fun(in);}



		template<class O, class I, class C, class D> static void comp(C& c, void (*fun)(O&, const I&), const D& in){
			I itmp = (I)in;
			O tmp;
			fun(tmp, itmp);
			c = (C) tmp;
			}
		template<class O, class I, class F, class C, class D> static void comp(C& c, void (F::*fun)(O&, const I&), const D& in){
			I itmp = (I)in;
			O tmp;
			fun(tmp, itmp);
			c = (C) tmp;
			}

		template<class O, class I, int size> static void comp(Tuple<O,size> &c, void (*fun)(O&, const I&), const Tuple<I, size>& in){
			for(unsigned int i =0; i < size;i++) fun(c[i], in[i]);
			}
		template<class O, class I, class F, int size> static void comp(Tuple<O,size> &c, void (F::*fun)(O&, const I&), const Tuple<I, size>& in){
			for(unsigned int i =0; i < size;i++) fun(c[i], in[i]);
			}

		template<class O, class I, unsigned int nbdim> static void comp(DataGrid<O,nbdim> &c, void (*fun)(O&, const I&), const DataGrid<I, nbdim>& in){
			c.setSizes(in.dims);
			for(unsigned int i = in.totsize(); i != 0xFFFFFFFF;i--) fun(c.data[i], in.data[i]);
			}
		template<class O, class I, class F, unsigned int nbdim> static void comp(DataGrid<O,nbdim> &c, void (F::*fun)(O&, const I&), const DataGrid<I,nbdim>& in){
			c.setSizes(in.dims);
			for(unsigned int i = in.totsize(); i != 0xFFFFFFFF;i--) fun(c.data[i], in.data[i]);
			}



	//	template<class A> static ExCo<A>::TYPE_norm norm_cmp_set(const A&);
	//	template<class A> static void norm_cmp_add(const A&, classarg<ExCo<A>::TYPE>& _out);
	//	template<class A> static double norm_cmp_extract(const A&);

	/*	template<unsigned int query,class A> static bool hasProperties(const A& a){
			if ((query & STRUCTPROP_IS_VALID)&&(!ExOp_flagbody<STRUCTPROP_IS_VALID,A>::hasProperties(a)) )return(false);
			if ((query & STRUCTPROP_HAS_VALID_INVERSE)&&(!ExOp_flagbody<STRUCTPROP_IS_VALID,A>::hasProperties(a)) )return(false);
			if ((query & STRUCTPROP_CAN_COMMUTE)&&(!ExOp_flagbody<STRUCTPROP_CAN_COMMUTE,A>::hasProperties(a)) )return(false);
			return(true);
		}*/


	//	template<class A> static bool isValid(const A& a){return( ExCo<A>::isValid(a));}

		template<class A> static bool isValid(const A& a){return ExFn< Exlisten_isValid<A, bool (ExCo<A>::SAFETYPE::*)() const >::ans >::isValid(a);}

		template<class A> static double getWeight(const A& a){return(ExCo<A>::getWeight(a));}

		template<class A> static void move(const A& a,const A& b){
			if (ExCo<A>::isMobile::ans) a =b;
			else a.moveFrom(b);
		}

		template<class A> static A* new_class(const A & a) {return (new A(a));} // sementic goody!


		//template<class A> static void showonline(const A& val,FILE* out){ExCo<A>::showonline(val,out);}
		//template<class A> static void show(const A& val){ExCo<A>::show(val,stdout);}
		//template<class A> static void showonline(const A& val){ExCo<A>::showonline(val,stdout);}


		template<class A> static SetComparison compare(const A& a, const A& b) {return( SetComparison( ((a < b) ?  30 | 64 | 256 :  0 ) | ((a > b) ? 29 | 32 | 128 :  0) )); }

		template<class A> inline static double dist(const A& a, const A& b){return sqrt(ExFn< Exlisten_pdist<A, double (ExCo<A>::SAFETYPE::*)(const A&)const >::ans >::pdist(a,b));}
		template<class A> inline static double pnorm(const A& a){return ExFn< Exlisten_pnorm<A, double (ExCo<A>::SAFETYPE::*)()const >::ans >::pnorm(a);}
		template<class A> inline static double norm(const A& a){return sqrt(ExFn< Exlisten_pnorm<A, double (ExCo<A>::SAFETYPE::*)()const >::ans >::pnorm(a));}
		template<class A> inline static double lognorm(const A& a){return sqrt(ExFn< Exlisten_lognorm<A, double (ExCo<A>::SAFETYPE::*)()const >::ans >::lognorm(a));}
		
		template<class A> inline static bool isZero(const A& a){return ExFn< Exlisten_isZero<A, bool (ExCo<A>::SAFETYPE::*)()const >::ans >::isZero(a);}
		template<class A> inline static bool isOne(const A& a){return ExFn< Exlisten_isOne<A, bool (ExCo<A>::SAFETYPE::*)()const >::ans >::isOne(a);}
		template<class A> inline static A& toZero(A& a){return ExFn< Exlisten_toZero<A, A& (ExCo<A>::SAFETYPE::*)() >::ans >::toZero(a);}
		template<class A> inline static A& toRand(A& a){return ExFn< Exlisten_toRand<A, A& (ExCo<A>::SAFETYPE::*)() >::ans >::toRand(a);}
		template<class A> inline static A& toOne(A& a){return ExFn< Exlisten_toOne<A, A& (ExCo<A>::SAFETYPE::*)() >::ans >::toOne(a);}
		template<class A> inline static A& toMin(A& a){return ExFn< Exlisten_toMin<A, A& (ExCo<A>::SAFETYPE::*)() >::ans >::toMin(a);}
		template<class A> inline static A& toMax(A& a){return ExFn< Exlisten_toMax<A, A& (ExCo<A>::SAFETYPE::*)() >::ans >::toMax(a);}

		template<class A> inline static A mkone(){A fout; ExOp::toOne(fout); return fout;}
		template<class A> inline static A mkzero(){A fout; ExOp::toZero(fout); return fout;}

		template<class A> static A& tonegative(A& a){return ExFn< Exlisten_tonegative<A, A& (ExCo<A>::SAFETYPE::*)() >::ans >::tonegative(a);}
		template<class A> static A& toinverse(A& a){return ExFn< Exlisten_toinverse<A, A& (ExCo<A>::SAFETYPE::*)() >::ans >::toinverse(a);}
		template<class A> static A& tosquare(A& a){return ExFn< Exlisten_tosquare<A, A& (ExCo<A>::SAFETYPE::*)() >::ans >::tosquare(a);}
		template<class A> static void tointpow(A& a, const int pow){ExFn< Exlisten_tointpow<A, void (ExCo<A>::SAFETYPE::*)(const int) >::ans >::tointpow(a,pow);}

		template<class A> inline static A& toadd(A& a, const A& b){return ExFn< Exlisten_toadd<A, const A& (ExCo<A>::SAFETYPE::*)(const A& b) >::ans >::toadd(a,b);}
		template<class A> inline static A& tosub(A& a, const A& b){return ExFn< Exlisten_tosub<A, const A& (ExCo<A>::SAFETYPE::*)(const A& b) >::ans >::tosub(a,b);}
		template<class A> inline static A& tomult(A& a, const A& b){return ExFn< Exlisten_tomult<A, const A& (ExCo<A>::SAFETYPE::*)(const A& b) >::ans >::tomult(a,b);}
		template<class A> inline static A& todivi(A& a, const A& b){return ExFn< Exlisten_todivi<A, const A& (ExCo<A>::SAFETYPE::*)(const A& b) >::ans >::todivi(a,b);}
		template<class A, class B> inline static A& toadd(A& a, const B& b){return ExFn< Exlisten_toadd<A, const A& (ExCo<A>::SAFETYPE::*)(const B& b) >::ans >::toadd(a,b);}
		template<class A, class B> inline static A& tosub(A& a, const B& b){return ExFn< Exlisten_tosub<A, const A& (ExCo<A>::SAFETYPE::*)(const B& b) >::ans >::tosub(a,b);}
		template<class A, class B> inline static A& tomult(A& a, const B& b){return ExFn< Exlisten_tomult<A, const A& (ExCo<A>::SAFETYPE::*)(const B& b) >::ans >::tomult(a,b);}
		template<class A, class B> inline static A& todivi(A& a, const B& b){return ExFn< Exlisten_todivi<A, const A& (ExCo<A>::SAFETYPE::*)(const B& b) >::ans >::todivi(a,b);}

		template<class A, class B> inline static const A& tobmul(A& a, const B& b){return ExFn< Exlisten_tobmul<A, const A& (ExCo<A>::SAFETYPE::*)(const B& b) >::ans >::tobmul(a,b);}
//		template<class A, class B> inline static const A& toRmult(A& a, const B& b){return ExFn< Exlisten_toRmult<A, void (ExCo<A>::SAFETYPE::*)(const B& b) >::ans >::toRmult(a,b);}

		template<class A, class B> inline static typename ExCo<A>::LMUL_TYPE tobmul(const A& a, const B& b){return ExFn< Exlisten_mkbmul<A, void (ExCo<A>::SAFETYPE::*)(const B& b)const >::ans >::toLmult(a,b);}
//		template<class A> inline static typename ExCo<A>::INNER_TYPE mkRmult(const A& a, const A& b){return ExFn< Exlisten_mkrmul<A, void (ExCo<A>::SAFETYPE::*)(const A& b)const >::ans >::toRmult(a,b);}
        template<class A> inline static typename ExCo<A>::TRJU_TYPE mktrju(const A& a){return ExFn< Exlisten_mktrju<A, typename ExCo<A>::TRJU_TYPE (ExCo<A>::SAFETYPE::*)()const  >::ans >::mktrju(a);}

		template<class A> inline static A mknegative(const A& a){return ExFn< Exlisten_mksub<A, A (ExCo<A>::SAFETYPE::*)()const >::ans >::mknegative(a);}
		template<class A> inline static A mkinverse(const A& a){return ExFn< Exlisten_mkinverse<A, A (ExCo<A>::SAFETYPE::*)()const >::ans >::mkinverse(a);}
		template<class A> inline static A mksquare(const A& a){return ExFn< Exlisten_mksquare<A, A (ExCo<A>::SAFETYPE::*)()const  >::ans >::mksquare(a);}

		template<class A> inline static typename ExCo<A>::REAL_TYPE mkrealproj(const A& a){return ExFn< Exlisten_mkrealproj<A, typename ExCo<A>::REAL_TYPE (ExCo<A>::SAFETYPE::*)()const  >::ans >::mkrealproj(a);}
		template<class A> inline static typename ExCo<A>::REAL_TYPE mkimmaproj(const A& a){return ExFn< Exlisten_mkimmaproj<A, typename ExCo<A>::REAL_TYPE (ExCo<A>::SAFETYPE::*)()const  >::ans >::mkimmaproj(a);}
		template<class A> inline static typename ExCo<A>::REAL_TYPE mkjmmaproj(const A& a){return ExFn< Exlisten_mkjmmaproj<A, typename ExCo<A>::REAL_TYPE (ExCo<A>::SAFETYPE::*)()const  >::ans >::mkjmmaproj(a);}
		template<class A> inline static typename ExCo<A>::REAL_TYPE mkkmmaproj(const A& a){return ExFn< Exlisten_mkkmmaproj<A, typename ExCo<A>::REAL_TYPE (ExCo<A>::SAFETYPE::*)()const  >::ans >::mkkmmaproj(a);}
        template<class A> inline static typename ExCo<A>::GAUS_TYPE mkgaussstat(const A& a,double weight = 1.0f){return ExFn< Exlisten_mkgaussstat<A, typename ExCo<A>::GAUS_TYPE (ExCo<A>::SAFETYPE::*)(double&)const  >::ans >::mkgaussstat(a,weight);}


		template<class A, class B> inline static typename STDRETTYPE2<A,B>::PLUS_TYPE mkadd(const A& a, const B& b){return ExFn< Exlisten_mkadd<A, typename STDRETTYPE2<A,B>::PLUS_TYPE (ExCo<A>::SAFETYPE::*)(const B& b)const >::ans >::mkadd(a,b);}
		template<class A, class B> inline static typename STDRETTYPE2<A,B>::MINU_TYPE mksub(const A& a, const B& b){return ExFn< Exlisten_mksub<A, typename STDRETTYPE2<A,B>::MINU_TYPE (ExCo<A>::SAFETYPE::*)(const B& b)const >::ans >::mksub(a,b);}
		template<class A, class B> inline static typename STDRETTYPE2<A,B>::PROD_TYPE mkmult(const A& a, const B& b){return ExFn< Exlisten_mkmult<A, typename STDRETTYPE2<A,B>::PROD_TYPE (ExCo<A>::SAFETYPE::*)(const B& b)const >::ans >::mkmult(a,b);}
		template<class A, class B> inline static typename STDRETTYPE2<A,B>::DIVI_TYPE mkdivi(const A& a, const B& b){return ExFn< Exlisten_mkdivi<A, typename STDRETTYPE2<A,B>::DIVI_TYPE (ExCo<A>::SAFETYPE::*)(const B& b)const >::ans >::mkdivi(a,b);}

		template<class A> inline static A mkintpow(const A& a, const int pow){return ExFn< Exlisten_mkintpow<A, A (ExCo<A>::SAFETYPE::*)(const int)const >::ans >::mkintpow(a,pow);}
		template<class A> inline static A mkinvintpow(const A& a, const int pow){return ExFn< Exlisten_mkinvintpow<A, A (ExCo<A>::SAFETYPE::*)(const int)const >::ans >::mkinvintpow(a,pow);}
		template<class A> inline static setcomparison setcmp(const A& a, const A& b) {return ExFn< Exlisten_setcmp<A, void (ExCo<A>::SAFETYPE::*)(const A&)const >::ans >::setcmp(a,b);}

// comparisons

        template<class A, class B> inline static bool isLT(const A &a, const B &b){ return( ExFn< Exlisten_lt<A, bool (ExCo<A>::SAFETYPE::*)(const B&)const >::ans >::comp_lt(a,b));}
        template<class A, class B> inline static bool isGT(const A &a, const B &b){ return( ExFn< Exlisten_gt<A, bool (ExCo<A>::SAFETYPE::*)(const B&)const >::ans >::comp_gt(a,b));}
        template<class A, class B> inline static bool isLE(const A &a, const B &b){ return( ExFn< Exlisten_le<A, bool (ExCo<A>::SAFETYPE::*)(const B&)const >::ans >::comp_le(a,b));}
        template<class A, class B> inline static bool isGE(const A &a, const B &b){ return( ExFn< Exlisten_ge<A, bool (ExCo<A>::SAFETYPE::*)(const B&)const >::ans >::comp_ge(a,b));}
        template<class A, class B> inline static bool isEQ(const A &a, const B &b){ return( ExFn< Exlisten_eq<A, bool (ExCo<A>::SAFETYPE::*)(const B&)const >::ans >::comp_eq(a,b));}
        template<class A, class B> inline static bool isNQ(const A &a, const B &b){ return( ExFn< Exlisten_nq<A, bool (ExCo<A>::SAFETYPE::*)(const B&)const >::ans >::comp_nq(a,b));}



		// same as operator ==, but the source looses the ownership of its pointers
		template<class A> inline static A& memmove(A& a, A& o){return ExFn< Exlisten_memmove<A, A& (ExCo<A>::SAFETYPE::*)( A &) >::ans >::memmove(a,o);}

		template<class A> inline static double getWeight(A& a){return  ExFn< Exlisten_getWeight<A, double (ExCo<A>::SAFETYPE::*)()const >::ans >::getWeight(a);}
		template<class A> inline static void save(const A& w, FILE *f) {ExFn< Exlisten_save<A, void ( ExCo<A>::SAFETYPE::*)(FILE*) const >::ans >::save(w,f);}
		template<class A> inline static void load(A& w, FILE *f, unsigned int lenght = 0xFFFFFFFF) {ExFn< Exlisten_load<A, void ( ExCo<A>::SAFETYPE::*)(FILE*, unsigned int)>::ans >::load(w,f,lenght);}


		// THE SHOW COMMAND! level 1: '\n' forbidden, level 2: '\t' forbidden, level 3: '[;]' forbidden, level 4: '(,)' forbidden,
		template<class A> inline static void show(const A& a, FILE* f_out = stdout, int level = 0){return  ExFn< Exlisten_show<A, void (ExCo<A>::SAFETYPE::*)(FILE* f, int level)const >::ans >::show(a,f_out,level);}
        template<class A> inline static string type_tostring(const A& a){return  ExFn< Exlisten_type_tostring<A, string (ExCo<A>::SAFETYPE::*)() const >::ans >::type_tostring(a);}

		// 1 arg operator
		template<class F, class A> inline static void comp(F &func, A &a){
			ExFn< (Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(A&)>::ans)
				||(Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(const A&)>::ans)
				||(Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(A&) const>::ans)
			 >::compose(func,a);}
		template<class F, class A> inline static void comp(const F &func, A &a){ExFn< (Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(A&) const>::ans)>::compose(func,a);}
		template<class F, class A> inline static void comp(F &func, const A &a){ExFn< (Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(const A&)>::ans)>::compose(func,a);}

		// 2  arg operator
		template<class F, class A, class B> inline static void comp(F &func, A &a, B &b){
			ExFn< (Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(A& , B&)>::ans)
				||(Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(A& , const B&)>::ans)
				||(Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(const A& , const B&)>::ans)
				||(Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(A& , B&) const>::ans)
				||(Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(A& , const B&) const>::ans)
			 >::compose(func,a,b);}

		template<class F, class A, class B> inline static void comp(const F &func, A &a, B &b){
			ExFn< (Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(A& , B&) const >::ans)
			    ||(Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(A& , const B&) const>::ans) >::compose(func,a,b);}

		template<class F, class A, class B> inline static void comp(F &func, A &a, const B &b){
			ExFn< (Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(A& , const B&)>::ans)
			    ||(Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(A& , const B&) const>::ans)
			    ||(Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(const A& , const B&) >::ans) >::compose(func,a,b);}

		template<class F, class A, class B> inline static void comp(const F &func, A &a, const B &b){
			ExFn< Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(A& , const B&) const >::ans >::compose(func,a,b);}

		template<class F, class A, class B> inline static void comp(F &func, const A &a, const B &b){
			ExFn< Exlisten_mainoper<F, void (ExCo<F>::SAFETYPE::*)(const A& , const B&)>::ans >::compose(func,a,b);}

						// a * b^H
	//	template<class A> inline static void mkcomult(const A& a,const A& b)



	};


	template<int size>
	class ExOpa{
	public:
		template<class A> static Tuple<A, size> intPows(const A& v){
			Tuple<A, size> _out;
			_out[0] = v;
			int mp =1;
			int i;
			for(i=1;i<size;i++){
				_out[i] = _out[mp -1] * _out[i- mp];
				if (mp<<1 == i+1) mp <<=1;
			}
			return(_out);
		}
	};




template<class NODE, int nbrel>
class Forest : public Vector< pair<Tuple<unsigned int, nbrel>, NODE > >{
public:
	unsigned int nbroots; // number of roots
	void cluster(Vector<NODE> &data, double (*metric)(const NODE&, const NODE&), NODE (*merge)(const NODE&, const NODE&) = NULL);
	void cluster_singlelink(Vector<NODE> &data, double (*metric)(const NODE&, const NODE&));
	void cluster_completelink(Vector<NODE> &data, double (*metric)(const NODE&, const NODE&));

    template<class DATA, unsigned int TSIZE > void cluster_bhattacharryya( const Vector< Tuple<DATA, TSIZE > > &data, NODE (*report)(const GaussElem< Tuple<DATA, TSIZE > > &,const GaussElem< Tuple<DATA, TSIZE > > &));
    template<class DATA, unsigned int TSIZE > void cluster_likelihood_ratio( const Vector< Tuple< WeightElem<DATA, 2> , TSIZE > > &data, NODE (*report)(const GaussElem< Tuple<DATA, TSIZE > > &,const GaussElem< Tuple<DATA, TSIZE > > &));

    template<class DATA, unsigned int TSIZE > void cluster_likelihood_ratio( const Vector< GaussElem< Tuple<DATA, TSIZE > > > &data, NODE (*report)(const GaussElem< Tuple<DATA, TSIZE > > &,const GaussElem< Tuple<DATA, TSIZE > > &));
    template<class DATA> void cluster_likelihood_ratio( const Vector< GaussElem< Tuple<DATA, 0u > > > &data, NODE (*report)(const GaussElem< Tuple<DATA, 0u > > &,const GaussElem< Tuple<DATA, 0u > > &));

	template<class DATA, unsigned int TSIZE > void cluster_bhattacharryya_simple( const Vector< Tuple< WeightElem<DATA, 2> , TSIZE > > &data, NODE (*report)(const GaussElem< Tuple<DATA, TSIZE > > &,const GaussElem< Tuple<DATA, TSIZE > > &));
    template<class DATA, unsigned int TSIZE > void cluster_bhattacharryya_complex( const Vector< GaussElem< Tuple<DATA, TSIZE > > > &data, NODE (*report)(const GaussElem< Tuple<DATA, TSIZE > > &,const GaussElem< Tuple<DATA, TSIZE > > &));

	template<class DATA, unsigned int TSIZE > void cluster_likelihood( const Vector< Tuple< WeightElem<DATA, 2> , TSIZE > > &data, NODE (*report)(const GaussElem< Tuple<DATA, TSIZE > > &,const GaussElem< Tuple<DATA, TSIZE > > &));

	
	template<class DATA> void cluster(Vector<DATA> &data, NODE (*metric)(const DATA&, const DATA&), DATA (*merge)(const DATA&, const DATA&) = NULL);
	template<class DATA> void cluster_singlelink(Vector<DATA> &data, NODE (*metric)(const DATA&, const DATA&));
	template<class DATA> void cluster_completelink(Vector<DATA> &data, NODE (*metric)(const DATA&, const DATA&));

	template<class DATA> void cluster_majoritylink(Vector<DATA> &data, NODE (*metric)(const DATA&, const DATA&));
	template<class DATA> void partition_squarreroot(const Vector<DATA> &data, NODE (*metric)(const DATA&, const DATA&));
	template<class DATA, class COMP> void partition_squarreroot_nbcluster(const Vector<DATA> &data, COMP (*metric)(const DATA&, const DATA&), int nbgroups);

	void saveGTRfile(const char* const path, const char* const  outer, double (*metric)(const NODE&, const NODE&)) const;
	void saveGTRfile(const char* const path, const char* const  outer, double (*xform)(const NODE&) = NULL) const; // simple convertion otherwise
    template<class DATA> void saveCDTfile(FILE* f, Vector<char*> &header , Vector<DATA> &data) const;
	template<class DATA, unsigned int SIZE, unsigned int ORDER> void saveCDTfile_W(FILE* f, Vector<char*> &header , Vector< Tuple<WeightElem<DATA, ORDER>, SIZE >  > &data) const;
	template<class DATA, unsigned int SIZE> void saveCDTfile_W(FILE* f, Vector<char*> &header , Vector< GaussElem< Tuple<DATA, SIZE > > > &data) const;

	void makeLeftof(unsigned int left, unsigned int par);
	void makeRightof(unsigned int left, unsigned int par);
	bool hasParent(unsigned int node);
	bool hasLeft(unsigned int node);
	bool hasRight(unsigned int node);
	unsigned int getParent(unsigned int node)const;
	unsigned int getLeft(unsigned int node)const;
	unsigned int getRight(unsigned int node)const;
	void clear(unsigned int node);
	bool isclear(unsigned int node);

	void moveNode(unsigned int Node, unsigned int Target);


};

template<class C>
class Complex{
public:
		typedef C REAL_TYPE;
		typedef Complex<C> COMPLEX_TYPE;
		typedef Quaternion<C> QUATERNION_TYPE;
        typedef Complex<C> VECTOR_TYPE;
        typedef YESNO< ExCo<C>::IS_COMMUTATIVE::ans > IS_COMMUTATIVE; 
        Complex<C> mkVectorization()const{return(*this);}

		C data[2];
		// ExOp section
		static const bool IsPOD = ExCo<C>::IsPOD;
		// ExOp section
		static const bool NeedsAddLink = false; // containers needs to update addresses in link registers


		ExCoMeMdEcLaRe( Complex<C> );

		Complex(){}
		Complex(C const & value) {data[0] = value; data[1] = ExCo<C>::zero();}
		Complex(C const & real,C const & im) {data[0] = real; data[1] = im;}
		C& operator[](unsigned int);
		const C& operator[](unsigned int) const;

		Complex<C> inverse() const;
        Complex<C>& toOne(){ExOp::toOne(data[0]);ExOp::toZero(data[1]); return(*this);}
        Complex<C>& toZero(){ExOp::toZero(data[0]);ExOp::toZero(data[1]);return(*this);}
        Complex<C>& toRand(){ExOp::toRand(data[0]);ExOp::toRand(data[1]);return(*this);}
        Complex<C> mkinverse()const{return (Complex<C>(data[0],-data[1]) /= (data[0]*data[0] + data[1]*data[1]));}

		double pnorm() const;
		double norm() const {return sqrt(pnorm());}

		double sign() const;
		template<class A> Complex<C>& operator+=(const Complex<A>&);
		template<class A> Complex<C>& operator-=(const Complex<A>&);
		template<class A> Complex<C>& operator*=(const Complex<A>&);
		template<class A> Complex<C>& operator/=(const Complex<A>& other){return((*this) *= (other.inverse()));}

		template<class A> Complex<C>& operator+=(const A&);
		template<class A> Complex<C>& operator-=(const A&);
		template<class A> Complex<C>& operator*=(const A&);
		template<class A> Complex<C>& operator/=(const A&);


		template<class A> Complex< typename STDRETTYPE2<C,A>::PLUS_TYPE > operator+(const Complex<A>&) const;
		template<class A> Complex< typename STDRETTYPE2<C,A>::MINU_TYPE > operator-(const Complex<A>&) const;
		template<class A> Complex< typename STDRETTYPE2<C,A>::PLUS_TYPE > operator*(Complex<A> const & other) const;
		template<class A> Complex< typename STDRETTYPE2<C,A>::DIVI_TYPE > operator/(Complex<A> const & other) const;

		Complex< typename ExCo<C>::NEG_TYPE > operator-() const;

		setcomparison setcmp(const Complex<C> &) const;

		bool operator>(const Complex<C> &)const;
		bool operator>=(const Complex<C> &)const;
		bool operator<(const Complex<C> &)const;
		bool operator<=(const Complex<C> &)const;
		bool operator==(const Complex<C> &)const;
		bool operator!=(const Complex<C> &)const;


		Complex<C> mktrju() const {return Complex<C>((*this)[0],ExOp::mknegative((*this)[1]));}


		inline C mkrealproj()const{return data[0];}
		inline C mkimmaproj()const{return data[1];}
		inline C mkjmmaproj()const{return data[0];}
		inline C mkkmmaproj()const{return data[1];}
	};
	typedef Complex<double> complex;

	template<class C>
	class Quaternion{
	public:
		C data[4];
		Quaternion();
		typedef C REAL_TYPE;
		typedef Quaternion<C> COMPLEX_TYPE;
		typedef Quaternion<C> QUATERNION_TYPE;
        typedef Quaternion<C> VECTOR_TYPE;
        typedef YESNO<false> IS_COMMUTATIVE;
        Quaternion<C> mkVectorization()const{return(*this);}
        ExCoMeMdEcLaRe( Quaternion  <C> );


        const Quaternion<C>& to_normal(const C& _i,const C& _j,const C& _k);
        const Quaternion<C>& to_normal_and_scale(const C& _i,const C& _j,const C& _k,const C& _s);
		void mk_proj_matrix(TMatrix<C,3,3> &)const;
		void mk_proj_matrix(TMatrix<C,4,4> &)const;


        const Quaternion<C>& rotateX(double);
        const Quaternion<C>& rotateY(double);
        const Quaternion<C>& rotateZ(double);

        Tuple<C,3u> mkXvector() const;
        Tuple<C,3u> mkYvector() const;
        Tuple<C,3u> mkZvector() const;
        template<class O> void wrXvector(O*) const;
        template<class O> void wrYvector(O*) const;
        template<class O> void wrZvector(O*) const;

		double pnorm() const;
		double norm() const {return sqrt(pnorm());}

		const Quaternion<C>& inverse() const;

        Quaternion<C>& toOne(){ExOp::toOne(data[0]);ExOp::toZero(data[1]);ExOp::toZero(data[2]); ExOp::toZero(data[3]); return(*this);}
        Quaternion<C>& toZero(){ExOp::toZero(data[0]);ExOp::toZero(data[1]);ExOp::toZero(data[2]);ExOp::toZero(data[3]);return(*this);}
        Quaternion<C>& toRand(){ExOp::toRand(data[0]);ExOp::toRand(data[1]);ExOp::toRand(data[2]);ExOp::toRand(data[3]);return(*this);}

		const C& operator[](unsigned int) const;
		C& operator[](unsigned int);

        template<class A> Quaternion<C>& operator+=(const A& other){data[0] += other;return(*this);}
        template<class A> Quaternion<C>& operator-=(const A& other){data[0] -= other;return(*this);}
        template<class A> Quaternion<C>& operator*=(const A& other){data[0] *= other;data[1] *= other;data[2] *= other;data[3] *= other; return(*this);}
        template<class A> Quaternion<C>& operator/=(const A& other){data[0] /= other;data[1] /= other;data[2] /= other;data[3] /= other; return(*this);}

		template<class A> Quaternion<C>& operator+=(const Quaternion<A>& other);
		template<class A> Quaternion<C>& operator-=(const Quaternion<A>& other);

	//	template<class A> const Quaternion<C>& operator*=(const Quaternion<A>& other);
	//	template<class A> const Quaternion<C>& operator/=(const Quaternion<A>& other){return((*this) *= (other.inverse()));}

		// has NO multiplication

		inline C mkrealproj()const{return data[0];}
		inline C mkimmaproj()const{return data[1];}
		inline C mkjmmaproj()const{return data[2];}
		inline C mkkmmaproj()const{return data[3];}

        const Quaternion<C>& toUnitQuaternion(const Tuple<double, 3> value);
	};

template<int size>
class Permutation{
public:
	unsigned char map[size];
	Permutation();
	int operator[](const int &offset) const;
	Permutation<size> operator-() const;
	Permutation<size> operator+(const Permutation<size> &) const;

	template<class C, Tuple_flag Cflag, unsigned int esize> Tuple<C, size+esize, Cflag> operator+(const Tuple<C, size+esize, Cflag>&) const;
};


	// or mantissa [0,1] range

	template<class C = int>
	class angle{

	public:
	C ang;

	angle();
	angle(C const & value);
	operator double() const;
	operator complex() const;



	double real(){return(cos((double)(*this)));}
	double imag(){return(sin((double)(*this)));}
	complex getcomplex() const;
};


template<class C, class B>
class KeyElem{
public:
	static const bool IsPOD = false;
	static const bool NeedsAddLink = ExCo<B>::NeedsAddLink;
	C k;
	B d;

	KeyElem();
	KeyElem(const C& _k, const B& _n);

    KeyElem<C,B>& toRand(){ExOp::toRand(k);ExOp::toRand(d); return(*this);}
    KeyElem<C,B>& toZero(){ExOp::toZero(k);ExOp::toZero(d); return(*this);}

	KeyElem<C,B>& operator=(const KeyElem<C,B> &other){k = other.k; d = other.d; return(*this);}
    KeyElem<C,B>& memmove(KeyElem<C,B> &other){ExOp::memmove(k,other.k); ExOp::memmove(d,other.d); return(*this);}

	bool isValid() const;
	
	template<class A> bool operator>(const KeyElem<C,A> &other) const {return( (k > other.k) || ((k == other.k) && (d > other.d)) );}
	template<class A> bool operator<(const KeyElem<C,A> &other) const {return( (k < other.k) || ((k == other.k) && (d < other.d)));}
	template<class A> bool operator>=(const KeyElem<C,A> &other) const {return( k > other.k) || ((k == other.k) && (d >= other.d));}
	template<class A> bool operator<=(const KeyElem<C,A> &other) const {return( (k < other.k) || ((k == other.k) && (d <= other.d)));}
	template<class A> bool operator==(const KeyElem<C,A> &other) const {return( (k == other.k) && (d >= other.d));}
	template<class A> bool operator!=(const KeyElem<C,A> &other) const {return( (k != other.k) || (d != other.d));}
    bool operator>(const C &other) const {return( k > other);}
	bool operator<(const C &other) const {return( k < other);}
	bool operator>=(const C &other) const {return( k >= other);}
	bool operator<=(const C &other) const {return( k <= other);}
	bool operator==(const C &other) const {return( k == other);}
	bool operator!=(const C &other) const {return( k != other);}
    template<class A> bool operator>(const A &other) const {return( ExOp::isGT(k, other));}
	template<class A> bool operator<(const A &other) const {return( ExOp::isLT(k, other));}
	template<class A> bool operator>=(const A &other) const {return( ExOp::isGE(k, other));}
	template<class A> bool operator<=(const A &other) const {return( ExOp::isLE(k, other));}
	template<class A> bool operator==(const A &other) const {return( ExOp::isEQ(k, other));}
	template<class A> bool operator!=(const A &other) const {return( ExOp::isNQ(k, other));}

    template<class A,class D> bool operator>(const KeyElem<D,A> &other) const {return( ExOp::isGT(k, other.k));}
	template<class A,class D> bool operator<(const KeyElem<D,A> &other) const {return( ExOp::isLT(k, other.k));}
	template<class A,class D> bool operator>=(const KeyElem<D,A> &other) const {return( ExOp::isGE(k, other.k));}
	template<class A,class D> bool operator<=(const KeyElem<D,A> &other) const {return( ExOp::isLE(k, other.k));}
	template<class A,class D> bool operator==(const KeyElem<D,A> &other) const {return( ExOp::isEQ(k, other.k));}
	template<class A,class D> bool operator!=(const KeyElem<D,A> &other) const {return( ExOp::isNQ(k, other.k));}



	void show(FILE* out = stdout, int level=0) const{
		switch(level){
		case 0:
			fprintf(out,"Key: "); ExOp::show(k,out, 1);
			fprintf(out,"\nData: "); ExOp::show(d,out, 1);
			fprintf(out,"\n");
			break;
		case 1:
			fprintf(out,"Key: ["); ExOp::show(k,out, 2);
			fprintf(out,"]\tData: ["); ExOp::show(d,out, 2);
			fprintf(out,"]");
			break;
		case 2:
			fprintf(out,"K{"); ExOp::show(k,out, 3);
			fprintf(out,"}#D{"); ExOp::show(d,out, 3);
			fprintf(out,"}");
			break;
		case 3:
			fprintf(out,"("); ExOp::show(k,out, 4);
			fprintf(out,");("); ExOp::show(d,out, 4);
			fprintf(out,")");
			break;
		}
	}
    void save(FILE*f) const{ ExOp::save(k,f);ExOp::save(d,f);}
    void load(FILE*f, unsigned int l){ ExOp::load(k,f);ExOp::load(d,f);}
    string type_tostring()const{return string("KeyElem<") + ExOp::type_tostring(k) + string(",") + ExOp::type_tostring(d) + string(">");}
};

template<int buffersize>
class SuperString{
public:
	char buffer[buffersize];
	vector<int> chunks;

	SuperString();

	char* operator[](int chunkId);

	void setChunk(int chunkId, char* string, int stringlength =0);

	char* operator()();
};

enum tiffflag{
    TIFFFLAG_Compression = 259, 
    TIFFFLAG_PHOTOMETRICINTERPRETATION = 262, // NEEDED , short
    TIFFFLAG_DOCUMENTNAME= 269, // ASCII
    TIFFFLAG_PAGENAME = 285, // ASCII
    TIFFFLAG_XPOSITION = 286, // rationnal
    TIFFFLAG_YPOSITION = 287, // rationnal
    TIFFFLAG_PAGENUMBER = 297 // short
    };


template<int charsize, class T = void*>
class Anything{
public:
	union{
		char data[charsize];
		T* ptr;
	};
	T* operator->(){return(sizeof(T) > charsize ? ptr : (T*)data);}
};


template<class C,unsigned int size, Tuple_flag Cflag> class IsLFHPrimitive< Tuple<C, size,Cflag> > {public: enum{ans = true};};



template<class A, class B, unsigned int size, Tuple_flag Cflag>
class STDRETTYPE2<Tuple<A,size,Cflag> , Tuple<B,size,Cflag > >{
public:
	typedef Tuple<typename STDRETTYPE2<A,B>::DEF_TYPE, size, Cflag> DEF_TYPE; // default
	typedef Tuple<typename STDRETTYPE2<A,B>::PLUS_TYPE, size, Cflag> PLUS_TYPE; // operator+
	typedef Tuple<typename STDRETTYPE2<A,B>::MINU_TYPE, size, Cflag> MINU_TYPE; // operator-
	typedef Tuple<typename STDRETTYPE2<A,B>::PROD_TYPE, size, Cflag> PROD_TYPE; // operator*
	typedef Tuple<typename STDRETTYPE2<A,B>::DIVI_TYPE, size, Cflag> DIVI_TYPE; // operator/
	typedef Tuple<typename STDRETTYPE2<A,B>::OP2_TYPE, size, Cflag> OP2_TYPE; // operator()
};

template<class A, class B, unsigned int size, Tuple_flag Cflag>
class STDRETTYPE2<Tuple<A,size,Cflag> , B >{
public:
	typedef Tuple<typename STDRETTYPE2<A,B>::DEF_TYPE, size, Cflag> DEF_TYPE; // default
	typedef Tuple<typename STDRETTYPE2<A,B>::PLUS_TYPE, size, Cflag> PLUS_TYPE; // operator+
	typedef Tuple<typename STDRETTYPE2<A,B>::MINU_TYPE, size, Cflag> MINU_TYPE; // operator-
	typedef Tuple<typename STDRETTYPE2<A,B>::PROD_TYPE, size, Cflag> PROD_TYPE; // operator*
	typedef Tuple<typename STDRETTYPE2<A,B>::DIVI_TYPE, size, Cflag> DIVI_TYPE; // operator/
	typedef Tuple<typename STDRETTYPE2<A,B>::OP2_TYPE, size, Cflag> OP2_TYPE; // operator()
};

// for matrix operation, these are treated as a diagonal matrix (not a colunm/row vector, use TMatrix for that)
template<class C,unsigned int TSIZE, Tuple_flag Cflag>
class Tuple{ //: public Iterable{
public:
	C data[TSIZE];

    unsigned int getSize() const{return TSIZE;}

	typedef Tuple< typename ExCo<C>::REAL_TYPE, TSIZE> REAL_TYPE;
	typedef Tuple< typename ExCo<C>::COMPLEX_TYPE, TSIZE> COMPLEX_TYPE;
    typedef GaussScope< Tuple<C,TSIZE> , typename MT_IFTYPE<ExCo<C>::IS_COMMUTATIVE::ans, Tuple<C, TEMPLATE_TRIANGLE_NUMBER<TSIZE, 2>::ans >,TMatrix<C, TSIZE ,TSIZE> >::TYPE > GAUS_TYPE;
	typedef YESNO< ExCo<C>::IS_COMMUTATIVE::ans > IS_COMMUTATIVE; 

	// ExOp Section:
	static const bool IsPOD = ExCo<C>::IsPOD;
	static const bool NeedsAddLink = ExCo<C>::NeedsAddLink;

    typedef DataGrid<C,2> TRJU_TYPE;
    typedef DataGrid<C,2> LMUL_TYPE;
    typedef C INNER_TYPE;
	
	typedef TMatrix< C,TSIZE,TSIZE > OUTER_TYPE;
	
	typedef unsigned int ITERATOR_TYPE;
    template <class S> class SUBS_INNER{public: typedef Tuple<S, TSIZE, Cflag> TYPE;};

	//typedef YESNO< false > IsPOD;

	ExCoMeMdEcLaRe( LFHCONCAT3(Tuple<C, TSIZE, Cflag>) )

	// ExOp Section end

	Tuple(){}
//	Tuple(C const & val); // update all using single val
	Tuple(Tuple<C,TSIZE> const & clonefrom){for(unsigned int i=0;i<TSIZE;i++) data[i] =clonefrom.data[i];}
	template<class O> Tuple(Tuple<O,TSIZE> const & clonefrom){unsigned int i; for(i=0;i<TSIZE;i++) data[i] = C(clonefrom.data[i]);}

	operator const C*()const{return data;}
	operator C*(){return data;}

	Tuple(C const * const clonefrom);
	//operator Vector<C,LFHVECTOR_REMOTE> ();

	template<class O, unsigned int oTSIZE> Tuple(Tuple<O,oTSIZE,Cflag> const & other);
	int getTSIZE() const{ return(TSIZE);}

	void fourierTransform_routine();
	void invfourierTransform_routine();

    void fourierTransform_routine2();
	void invfourierTransform_routine2();
	void pow2_bitswap_permutation();

    void HouseHolderMultiply(const C * const vec, double denum2, unsigned int length );

   static void fourierTransform_routine(complex*);
	static void invfourierTransform_routine(complex*);

	static Tuple<complex, TSIZE,Cflag> bluesteinWindow(int wTSIZE);
    static void bluesteinWindow(complex*&,int wTSIZE);

	template<unsigned int superTSIZE> Tuple<C, TSIZE,Cflag> fourierTransform(const Tuple<complex, superTSIZE,Cflag>& bluewindow) const;
	template<unsigned int superTSIZE> Tuple<C, TSIZE,Cflag> invfourierTransform(const Tuple<complex, superTSIZE,Cflag>& bluewindow) const;
	template<unsigned int superTSIZE> Tuple<Complex<C>, TSIZE,Cflag> fourierTransformReal() const; // C cannot be multiplied by complex

	Tuple<C,TSIZE,Cflag>& toZero(){for(unsigned int i=0;i<TSIZE;i++) ExOp::toZero(data[i]);return(*this);}
    Tuple<C,TSIZE,Cflag>& toOne(){for(unsigned int i=0;i<TSIZE;i++) ExOp::toOne(data[i]);return(*this);}
    Tuple<C,TSIZE,Cflag>& toRand(){for(unsigned int i=0;i<TSIZE;i++) ExOp::toRand(data[i]);return(*this);}
    Tuple<C,TSIZE,Cflag>& toMin(){for(unsigned int i=0;i<TSIZE;i++) ExOp::toMin(data[i]);return(*this);}
    Tuple<C,TSIZE,Cflag>& toMax(){for(unsigned int i=0;i<TSIZE;i++) ExOp::toMax(data[i]);return(*this);}

    OUTER_TYPE mkouterprod(const Tuple<C, TSIZE> &)const;
	
	Tuple<C, TSIZE, Cflag> normalize() const;
	Tuple<C, TSIZE, Cflag> normalize(C const & norm) const;

	C max() const;
	C min() const;

	bool isValid() const;

	inline Tuple< typename ExCo<C>::REAL_TYPE, TSIZE> mkrealproj()const{Tuple< typename ExCo<C>::REAL_TYPE, TSIZE> fout; for(unsigned int i=0;i<TSIZE;i++) fout[i] = ExOp::mkrealproj(data[i]); return fout;}
	inline Tuple< typename ExCo<C>::REAL_TYPE, TSIZE> mkimmaproj()const{Tuple< typename ExCo<C>::REAL_TYPE, TSIZE> fout; for(unsigned int i=0;i<TSIZE;i++) fout[i] = ExOp::mkimmaproj(data[i]); return fout;}
	inline Tuple< typename ExCo<C>::REAL_TYPE, TSIZE> mkjmmaproj()const{Tuple< typename ExCo<C>::REAL_TYPE, TSIZE> fout; for(unsigned int i=0;i<TSIZE;i++) fout[i] = ExOp::mkjmmaproj(data[i]); return fout;}
	inline Tuple< typename ExCo<C>::REAL_TYPE, TSIZE> mkkmmaproj()const{Tuple< typename ExCo<C>::REAL_TYPE, TSIZE> fout; for(unsigned int i=0;i<TSIZE;i++) fout[i] = ExOp::mkkmmaproj(data[i]); return fout;}

	Tuple<C,TSIZE,Cflag>& operator=(Tuple<C,TSIZE,Cflag> const & other);
	template<class O> Tuple<C,TSIZE,Cflag>& operator=(O const & other);
    Tuple<C,TSIZE,Cflag>& memmove(Tuple<C,TSIZE,Cflag> & other);

//	template<class O> const Tuple<C,TSIZE,Cflag>& operator=(Tuple<O,TSIZE,Cflag> const & other);
	template<class O, unsigned int  OTSIZE> Tuple<C,TSIZE,Cflag>& operator=(Tuple<O,OTSIZE,Cflag> const & other);

#undef LFHTEMP
#define LFHTEMP template<class O,unsigned int oTSIZE, Tuple_flag Oflag>
	LFHTEMP char compare(const Tuple<O,oTSIZE, Oflag>& other) const;
	LFHTEMP bool operator>(const Tuple<O,oTSIZE, Oflag>& other) const{return(compare(other) == 1);}
	LFHTEMP bool operator<(const Tuple<O,oTSIZE, Oflag>& other) const{return(compare(other) == 2);}
	LFHTEMP bool operator>=(const Tuple<O,oTSIZE, Oflag>& other) const{return(compare(other) != 2);}
	LFHTEMP bool operator<=(const Tuple<O,oTSIZE, Oflag>& other) const{return(compare(other) != 1);}
	LFHTEMP bool operator==(const Tuple<O,oTSIZE, Oflag>& other) const{return(compare(other) == 0);}
	LFHTEMP bool operator!=(const Tuple<O,oTSIZE, Oflag>& other) const{return(compare(other) != 0);}

	template<class A_1> void operator() (Oper1<A_1> const & op); // not a match
	void operator() (Oper1<C> const & op); // match
	const Tuple<C,TSIZE,Cflag>& operator-() const;

	template<class A_1, class A_2, class C_2, unsigned int TSIZE_2> void operator() (Oper2<A_1,A_2> const & op, Tuple<C_2, TSIZE_2,Cflag> const & ); // not a match
	template<class C_2, unsigned int TSIZE_2> void operator() (Oper2<C,C_2> const & op, Tuple<C_2, TSIZE_2,Cflag> const & ); // match
	template<class A_1, class A_2, class C_2, unsigned int TSIZE_2> void operator() (Oper2<A_1,A_2> const & op, Tuple<C_2, TSIZE_2,Cflag> & ); // not a match
	template<class C_2, unsigned int TSIZE_2> void operator() (Oper2<C,C_2> const & op, Tuple<C_2, TSIZE_2,Cflag> & ); // match
#undef LFHTEMP
#define LFHTEMP template<class A_1, class A_2, class A_3, class C_2, unsigned int TSIZE_2, class C_3, unsigned int TSIZE_3>
	LFHTEMP	void operator() (Oper3<A_1,A_2,A_3> const & op, Tuple<C_2, TSIZE_2,Cflag> const &, Tuple<C_3, TSIZE_3,Cflag> const & ); // not a match
	LFHTEMP	void operator() (Oper3<A_1,A_2,A_3> const & op, Tuple<C_2, TSIZE_2,Cflag> &, Tuple<C_3, TSIZE_3,Cflag> const & ); // not a match
	LFHTEMP	void operator() (Oper3<A_1,A_2,A_3> const & op, Tuple<C_2, TSIZE_2,Cflag> &, Tuple<C_3, TSIZE_3,Cflag> & ); // not a match
#undef LFHTEMP
#define LFHTEMP template<class C_2, unsigned int TSIZE_2, class C_3, unsigned int TSIZE_3>
	LFHTEMP void operator() (Oper3<C,C_2,C_3> const & op, Tuple<C_2, TSIZE_2,Cflag> const &, Tuple<C_3, TSIZE_3,Cflag> const & ); // match
	LFHTEMP void operator() (Oper3<C,C_2,C_3> const & op, Tuple<C_2, TSIZE_2,Cflag> &, Tuple<C_3, TSIZE_3,Cflag> const & ); // match
	LFHTEMP void operator() (Oper3<C,C_2,C_3> const & op, Tuple<C_2, TSIZE_2,Cflag> &, Tuple<C_3, TSIZE_3,Cflag> & ); // match

#undef LFHTEMP
#define LFHTEMP template<class O>

	LFHTEMP Tuple<C,TSIZE,Cflag>& operator+=(Tuple<O,TSIZE,Cflag> const & other);
	LFHTEMP Tuple<C,TSIZE,Cflag>& operator-=(Tuple<O,TSIZE,Cflag> const & other);
	LFHTEMP Tuple<C,TSIZE,Cflag>& operator*=(Tuple<O,TSIZE,Cflag> const & other);
	LFHTEMP Tuple<C,TSIZE,Cflag>& operator/=(Tuple<O,TSIZE,Cflag> const & other);
	LFHTEMP Tuple< typename STDRETTYPE2<C,O>::PLUS_TYPE ,TSIZE,Cflag> operator+(Tuple<O,TSIZE,Cflag> const & other) const;
	LFHTEMP Tuple< typename STDRETTYPE2<C,O>::MINU_TYPE ,TSIZE,Cflag> operator-(Tuple<O,TSIZE,Cflag> const & other) const;
	LFHTEMP Tuple< typename STDRETTYPE2<C,O>::PROD_TYPE ,TSIZE,Cflag> operator*(Tuple<O,TSIZE,Cflag> const & other) const;
	LFHTEMP Tuple< typename STDRETTYPE2<C,O>::DIVI_TYPE ,TSIZE,Cflag> operator/(Tuple<O,TSIZE,Cflag> const & other) const;
	
	// expending to all values
	LFHTEMP Tuple<C,TSIZE,Cflag>& operator+=(O const & other);
	LFHTEMP Tuple<C,TSIZE,Cflag>& operator-=(O const & other);
	LFHTEMP Tuple<C,TSIZE,Cflag>& operator*=(O const & other);
	LFHTEMP Tuple<C,TSIZE,Cflag>& operator/=(O const & other);
	LFHTEMP	Tuple< typename STDRETTYPE2<C,O>::PLUS_TYPE ,TSIZE,Cflag> operator+(O const & other) const; //
	LFHTEMP Tuple< typename STDRETTYPE2<C,O>::MINU_TYPE ,TSIZE,Cflag> operator-(O const & other) const; //
	LFHTEMP Tuple< typename STDRETTYPE2<C,O>::PROD_TYPE ,TSIZE,Cflag> operator*(O const & other) const; //
	LFHTEMP Tuple< typename STDRETTYPE2<C,O>::DIVI_TYPE ,TSIZE,Cflag> operator/(O const & other) const; //

	LFHTEMP Tuple<C,TSIZE,Cflag>& operator*=(TMatrix<O,TSIZE,TSIZE> const & other);
	LFHTEMP Tuple<C,TSIZE,Cflag>& operator/=(TMatrix<O,TSIZE,TSIZE> const & other);
	
	LFHTEMP	Tuple<C,TSIZE,Cflag> operator+(KeyElem<unsigned int, O> const & other) const;
	LFHTEMP	Tuple<C,TSIZE,Cflag> operator-(KeyElem<unsigned int, O> const & other) const;
	LFHTEMP	Tuple<C,TSIZE,Cflag> operator*(KeyElem<unsigned int, O> const & other) const;
	LFHTEMP	Tuple<C,TSIZE,Cflag> operator/(KeyElem<unsigned int, O> const & other) const;
	LFHTEMP	Tuple<C,TSIZE,Cflag>& operator+=(KeyElem<unsigned int, O> const & other);
	LFHTEMP	Tuple<C,TSIZE,Cflag>& operator-=(KeyElem<unsigned int, O> const & other);
	LFHTEMP	Tuple<C,TSIZE,Cflag>& operator*=(KeyElem<unsigned int, O> const & other);
	LFHTEMP	Tuple<C,TSIZE,Cflag>& operator/=(KeyElem<unsigned int, O> const & other);

	//template<unsigned int PSIZE> Tuple<PolyThing<C,PSIZE>,TSIZE,Cflag> operator*(const PolyThing<C,PSIZE>&) const;
	
	
	C& operator[](int const pos);
	const C& operator[](int const pos) const;

	template<class D> operator Tuple<D,TSIZE> () const;
	double weight();

	inline static const Tuple<C,TSIZE,Cflag>& genConvolution_Gaussian(double std);


	// operators

	template<unsigned int pos> class Selector : LFHDECL_OPER2(C,LFHCONCAT3(Tuple<C, TSIZE,Cflag>));

	template<unsigned int firsTSIZE> class Concatenate : LFHDECL_OPER3(LFHCONCAT3(Tuple<C, TSIZE,Cflag>),LFHCONCAT3(Tuple<C,firsTSIZE,Cflag>), LFHCONCAT3(Tuple<C, TSIZE - firsTSIZE,Cflag>));
	template<unsigned int firsTSIZE> class Deconcatenate : LFHDECL_OPER3(LFHCONCAT3(Tuple<C, TSIZE,Cflag>),LFHCONCAT3(Tuple<C,firsTSIZE,Cflag>), LFHCONCAT3(Tuple<C, TSIZE + firsTSIZE,Cflag>));


	template<unsigned int pos> class SelectorWrite : LFHDECL_OPER2(LFHCONCAT3(Tuple<C, TSIZE,Cflag>),C);
	template<unsigned int TSIZEin, unsigned int pos_in, unsigned int pos_out> class SelectSelectorWrite : LFHDECL_OPER2(LFHCONCAT3(Tuple<C, TSIZE,Cflag>),LFHCONCAT3(Tuple<C, TSIZEin,Cflag>));


//	class FiniteDifference: LFHDECL_OPER2(LFHCONCAT(Tuple<C, TSIZE>),LFHCONCAT(Tuple<C, TSIZE>));
	class FiniteDifference:  public Oper2< Tuple<C, TSIZE,Cflag> , Tuple<C, TSIZE,Cflag> > { public: void operator()(Tuple<C, TSIZE,Cflag>  &, Tuple<C, TSIZE,Cflag> &) const;};
	class FiniteDifferenceAssign : LFHDECL_OPER1(LFHCONCAT3(Tuple<C, TSIZE,Cflag>));
	class L2Norm : LFHDECL_OPER2(C,LFHCONCAT3(Tuple<C, TSIZE,Cflag>));

	class ArgMax : LFHDECL_OPER2(int,LFHCONCAT3(Tuple<C, TSIZE,Cflag>));

	class MassCenter : LFHDECL_OPER2(double,LFHCONCAT3(Tuple<C, TSIZE,Cflag>));
	class MassCenter2 : LFHDECL_OPER2(LFHCONCAT3(Tuple<double, 3,Cflag>),LFHCONCAT3(Tuple<C, TSIZE,Cflag>));


	template<unsigned int winTSIZE>	class Convolution : LFHDECL_OPER3(LFHCONCAT3(Tuple<C, TSIZE,Cflag>),LFHCONCAT3(Tuple<C, TSIZE,Cflag>),LFHCONCAT3(Tuple<C, winTSIZE,Cflag>));

	class MakeTuple : LFHDECL_OPER2(LFHCONCAT3(Tuple<C, TSIZE,Cflag>), C);
	class PopWeight : LFHDECL_OPER2(LFHCONCAT3(Tuple<C, TSIZE,Cflag>),LFHCONCAT3(Tuple<C, TSIZE+1,Cflag>));
	class PushWeight : LFHDECL_OPER3(LFHCONCAT3(Tuple<C, TSIZE,Cflag>),LFHCONCAT3(Tuple<C, TSIZE-1,Cflag>), C);

	class VarNormalFiniteDifference:  LFHDECL_OPER2(LFHCONCAT3(Tuple<C, TSIZE,Cflag>),LFHCONCAT3(Tuple<C, TSIZE*2,Cflag>));

	C dotProduct(const Tuple<C,TSIZE, Cflag>& other)const{
		C _out = data[0] * other[0];
		int i;
		for(i=1;i<TSIZE;i++) _out += data[i] * other[i];
		return(_out);
	}

	template<unsigned int order> Tuple<C,  TEMPLATE_TRIANGLE_NUMBER<TSIZE, order+1>::ans  > genProducts(){
		C partial[order-1];
		Tuple<C,  TEMPLATE_TRIANGLE_NUMBER<TSIZE,order+1>::ans  > _out;
		int coor[order];
		int cor;
		int j;
		for(j=0;j<TSIZE;j++) _out[j] = data[j];
		int k = TSIZE;
		for(cor=2;cor<=order;cor++){
			for(j=0;j<cor-1;j++) {
				coor[j] = j;
				partial[j] = (j ==0) ?  data[0] : data[j] * partial[j-1];
			}
			coor[j] = j;
			do{
				_out[k] = data[k] * partial[cor-2];k++;

				j = cor-1;
				if (coor[j] = TSIZE-1) {
					for(j--; j>=0 ;j--) if (coor[j] != coor[j+1] -1) break;
					if (j < 0) break;
					coor[j]++;
					partial[j] = (j ==0) ?  data[0] : data[j] * partial[j-1];
					for(j++;j<cor-1;j++) {
						coor[j] = coor[j-1]+1;
						partial[j] = data[ coor[j] ] * partial[j-1];
					}
				}
			}while(true);
		}
		return(_out);
	}
	double pnorm() const;
	double norm() const;

    GaussScope< Tuple<C,TSIZE> , typename MT_IFTYPE<ExCo<C>::IS_COMMUTATIVE::ans, Tuple<C, TEMPLATE_TRIANGLE_NUMBER<TSIZE, 2>::ans> , TMatrix<C,TSIZE,TSIZE> >::TYPE > mkgaussstat(double &w) const;
    
    static double overlap(const GaussScope< Tuple<C,TSIZE> , typename MT_IFTYPE<ExCo<C>::IS_COMMUTATIVE::ans, Tuple<C, TEMPLATE_TRIANGLE_NUMBER<TSIZE, 2>::ans> , TMatrix<C,TSIZE,TSIZE> >::TYPE >&a, const GaussScope< Tuple<C,TSIZE> , typename MT_IFTYPE<ExCo<C>::IS_COMMUTATIVE::ans, Tuple<C, TEMPLATE_TRIANGLE_NUMBER<TSIZE, 2>::ans> , TMatrix<C,TSIZE,TSIZE> >::TYPE >&b);
};

template<class C, Tuple_flag Cflag>
class Tuple<C,0u, Cflag>{
public:
    unsigned int tup_size;
    C* data;
	static const bool IsPOD = false;

    Tuple(): tup_size(0){} // size == 0 iif data is unitialized
    Tuple(const Tuple<C,0u, Cflag> &other);
    template<unsigned int OSIZE> Tuple(const Tuple<C,OSIZE, Cflag> &other);

    ~Tuple();
    Tuple<C,0u, Cflag>& operator=(const Tuple<C,0u, Cflag> &other){clear();if (other.tup_size) {this->setSize(other.tup_size); for(unsigned int i=0;i<tup_size;i++) data[i] = other.data[i];}else tup_size = 0u; return(*this);}
    template<unsigned int OSIZE> Tuple<C,0u, Cflag>& operator=(const Tuple<C,OSIZE, Cflag> &other){this->setSize(OSIZE); for(unsigned int i=0;i<OSIZE;i++) data[i] = other.data[i]; return(*this);}
    void clear(){if (tup_size) delete[](data); tup_size = 0u;}

    unsigned int getSize() const{return tup_size;}
    void setSize(unsigned int s);
	
	bool isValid()const{for(unsigned int i =0; i < tup_size;i++) if (!ExOp::isValid(data[i])) return false; return true;}
	
    Tuple<C,0u, Cflag>& toZero(){for(unsigned int i=0;i<tup_size;i++) ExOp::toZero(data[i]); return(*this);}
    Tuple<C,0u, Cflag>& toOne(){for(unsigned int i=0;i<tup_size;i++) ExOp::toOne(data[i]); return(*this);}
    Tuple<C,0u, Cflag>& toRand(){for(unsigned int i=0;i<tup_size;i++) ExOp::toRand(data[i]); return(*this);}

    Tuple<C,0u, Cflag>& operator+=(const Tuple<C,0u, Cflag> &other){unsigned int j = tup_size < other.tup_size ? tup_size : other.tup_size; for(unsigned int i=0;i<j;i++) data[i] += other.data[i]; return(*this);}
    Tuple<C,0u, Cflag>& operator-=(const Tuple<C,0u, Cflag> &other){unsigned int j = tup_size < other.tup_size ? tup_size : other.tup_size; for(unsigned int i=0;i<j;i++) data[i] -= other.data[i]; return(*this);}
    Tuple<C,0u, Cflag>& operator*=(const Tuple<C,0u, Cflag> &other){unsigned int j = tup_size < other.tup_size ? tup_size : other.tup_size; for(unsigned int i=0;i<j;i++) data[i] *= other.data[i]; return(*this);}
    Tuple<C,0u, Cflag>& operator/=(const Tuple<C,0u, Cflag> &other){unsigned int j = tup_size < other.tup_size ? tup_size : other.tup_size; for(unsigned int i=0;i<j;i++) data[i] /= other.data[i]; return(*this);}
    template<unsigned int OSIZE> Tuple<C,0u, Cflag>& operator+=(const Tuple<C,OSIZE, Cflag> &other){unsigned int j = tup_size < OSIZE ? tup_size : OSIZE; for(unsigned int i=0;i<j;i++) data[i] += other.data[i]; return(*this);}
    template<unsigned int OSIZE> Tuple<C,0u, Cflag>& operator-=(const Tuple<C,OSIZE, Cflag> &other){unsigned int j = tup_size < OSIZE ? tup_size : OSIZE; for(unsigned int i=0;i<j;i++) data[i] -= other.data[i]; return(*this);}
    template<unsigned int OSIZE> Tuple<C,0u, Cflag>& operator*=(const Tuple<C,OSIZE, Cflag> &other){unsigned int j = tup_size < OSIZE ? tup_size : OSIZE; for(unsigned int i=0;i<j;i++) data[i] *= other.data[i]; return(*this);}
    template<unsigned int OSIZE> Tuple<C,0u, Cflag>& operator/=(const Tuple<C,OSIZE, Cflag> &other){unsigned int j = tup_size < OSIZE ? tup_size : OSIZE; for(unsigned int i=0;i<j;i++) data[i] /= other.data[i]; return(*this);}
    template<class O> Tuple<C,0u, Cflag>& operator+=(const O &other){for(unsigned int i=0;i<tup_size;i++) data[i] += other; return(*this);}
    template<class O> Tuple<C,0u, Cflag>& operator-=(const O &other){for(unsigned int i=0;i<tup_size;i++) data[i] -= other; return(*this);}
    template<class O> Tuple<C,0u, Cflag>& operator*=(const O &other){for(unsigned int i=0;i<tup_size;i++) data[i] *= other; return(*this);}
    template<class O> Tuple<C,0u, Cflag>& operator/=(const O &other){for(unsigned int i=0;i<tup_size;i++) data[i] /= other; return(*this);}

    Tuple<C,0u, Cflag> operator+(const Tuple<C,0u, Cflag> &other)const{return((Tuple<C,0u, Cflag>(*this)) += other);}
    Tuple<C,0u, Cflag> operator-(const Tuple<C,0u, Cflag> &other)const{return((Tuple<C,0u, Cflag>(*this)) -= other);}
    Tuple<C,0u, Cflag> operator*(const Tuple<C,0u, Cflag> &other)const{return((Tuple<C,0u, Cflag>(*this)) *= other);}
    Tuple<C,0u, Cflag> operator/(const Tuple<C,0u, Cflag> &other)const{return((Tuple<C,0u, Cflag>(*this)) /= other);}
    template<unsigned int OSIZE> Tuple<C,0u, Cflag> operator+(const Tuple<C,OSIZE, Cflag> &other)const{return((Tuple<C,0u, Cflag>(*this)) += other);}
    template<unsigned int OSIZE> Tuple<C,0u, Cflag> operator-(const Tuple<C,OSIZE, Cflag> &other)const{return((Tuple<C,0u, Cflag>(*this)) -= other);}
    template<unsigned int OSIZE> Tuple<C,0u, Cflag> operator*(const Tuple<C,OSIZE, Cflag> &other)const{return((Tuple<C,0u, Cflag>(*this)) *= other);}
    template<unsigned int OSIZE> Tuple<C,0u, Cflag> operator/(const Tuple<C,OSIZE, Cflag> &other)const{return((Tuple<C,0u, Cflag>(*this)) /= other);}
    template<class O> Tuple<C,0u, Cflag> operator+(const O &other)const{return((Tuple<C,0u, Cflag>(*this)) += other);}
    template<class O> Tuple<C,0u, Cflag> operator-(const O &other)const{return((Tuple<C,0u, Cflag>(*this)) -= other);}
    template<class O> Tuple<C,0u, Cflag> operator*(const O &other)const{return((Tuple<C,0u, Cflag>(*this)) *= other);}
    template<class O> Tuple<C,0u, Cflag> operator/(const O &other)const{return((Tuple<C,0u, Cflag>(*this)) /= other);}

    void HouseHolderMultiply(const C * const vec, double denum2, unsigned int length );

    const C& operator[](unsigned int index)const{return data[index];}
    C& operator[](unsigned int index){return data[index];}

	Tuple<C,0u, Cflag>& memmove(Tuple<C,0u, Cflag>& source){clear(); data = source.data; tup_size = source.tup_size; source.tup_size =0;}
	
	void show(FILE * f = stdout, int level = 0) const;
    string type_tostring()const;
};


template<class C, int size> struct IteratorScope< C, Tuple<C,size> > : public Iterator< C, Tuple<C,size> > {
	enum {valid = true};
	unsigned int i;
	IteratorScope(): i(0){}
	void init(Tuple<C,size> & ob) {i=0;}
	C* next(Tuple<C,size> & ob){ i++; return((i <= size)? ob.data + i - 1 : NULL);}
};

template<class C, class D, int size> struct IteratorScope<D, Tuple<C,size> > : public Iterator<D, Tuple<C,size> >{
	enum {valid = true};
	unsigned int i;
	IteratorScope<D, C> j;
	IteratorScope(Tuple<C,size> & ob){}
	void init(Tuple<C,size> & ob) {i=0; j.init();}

	C* next(Tuple<C,size> & ob){
		D* _out = j.next(ob.data[i]);
		while (_out == NULL){
			i++;
		if (i == size) return(NULL);
		j.init();
		_out = j.next(ob.data[i]);
		}
		return(_out);
		}
};

template<class C, int size,Tuple_flag Cflag> struct isTypeEquivalent< C*, Tuple<C,size, Cflag> > {enum {ans = true }; };
template<class C, int size,Tuple_flag Cflag> struct isTypeEquivalent< Tuple<C,size, Cflag> , C*> {enum {ans = true }; };

	class ProgressBarPrint{
		public:
		unsigned int state;
		unsigned int lenght;
		unsigned int lasttime;
		void start(const char*);
		void update(double fraction);
		void finish();
	};

class ArgumentParser{
public:
	void readList(char* const, Vector<unsigned int> & _out);
	void readList(char* const, vector<unsigned int> & _out);
	int operator()(int argv, char * const * args);
	virtual void nbaddtoken(char const * const token, int& min, int& max)=0;
	virtual void store(char* const * token, int nbtoken)=0; // return nb token used
	virtual int defstore(char* const * token, int nbtoken)=0; // return nb token used
	virtual void help() {printf("help unavailable!\n");}
};

template<int nbstate>
class HMM{
public:
	static const bool IsPOD = false;
	static const bool NeedsAddLink = false;

	//typedef LFHPrimitive::YESNO<true> IsPOD;
	TMatrix<double,nbstate,nbstate> transition;
	Tuple<double,nbstate> boundary;
	WeightElem< TMatrix<double,nbstate,nbstate>, 1 > learned_transition; // EM step buffer;


	template <Tuple_flag TF> void init(double transitrate, Tuple<double, nbstate, TF> const &);

	void EMinit();
	template<int lenght> void runHMM(Tuple<Tuple<double, nbstate>, lenght> &likelyhoods, double* weights = NULL);
	void runHMM(Tuple<double, nbstate>* likelyhoods, int lenght, double* weights = NULL);
	template< unsigned int nbdimss> DataGrid< Tuple<double, nbstate>, nbdimss> runHMM(DataGrid< Tuple<double, nbstate>, nbdimss> const &likelyhoods, int direction, bool learn = false);
	template< unsigned int nbdimss> DataGrid< Tuple<double, nbstate>, nbdimss> runHMM_meta(DataGrid< Tuple<double, nbstate>, nbdimss> const &likelyhoods, int direction, double* weights = NULL);

	void swapstates(int a,int b);

	void EMfinit();

	void show(FILE* f= stdout)const;
	void save(FILE* f) const;
	void load(FILE* f, unsigned int size = 0);
};

template<class A,class B>
class Convert : public Oper2<A , B>{
	public:
	void operator()(A  &, B &) const;
	void operator()(A  &, const B &) const;
	};

template<class A>
class getNorm : public Oper2<double , A>{
public:
	void operator()(double &, A &) const;
};

// if any negative, its zero
template<class A>
class getPositiveNorm : public Oper2<double , A>{
public:
	void operator()(double &, A &) const;
};

template<class C>
class Square : LFHDECL_OPER2(C,C);

// Gradient Search scope stores the last given derivative of a Function on which gradient Ascent or Descent is Performed in order to check an bound on the non-linearity of the function.
class GradientSearchScope{
	int size;
	double* lastderiv;
	double lastvalue;
	double expdiff;
	double log_alpha;
public:
	GradientSearchScope(): lastderiv(NULL){}
	GradientSearchScope(const GradientSearchScope&);
	~GradientSearchScope() {delete[](lastderiv);}
	GradientSearchScope& operator=(const GradientSearchScope&);

	double operator()() const;

	void init(double initalalpha, int in_size);
	void registAscent(const double &value, double const * const deriv, double mod = 1.0f);
	void registDescent(const double &value, double const * const deriv, double mod = 1.0f);

	double updateAscent(const double &value, double * guess,  double const * const deriv); // returns a safe bound on the parameter changes;
	double updateDescent(const double &value, double * guess , double const * const deriv); // returns a safe bound on the parameter changes;

	void punish(double log_mag);

};

class CurvatureSearchScope{
	unsigned int size;
	double* lastderiv;
	double lastvalue;
	double expdiff;
	double* curvature;
	double log_alpha; // temp... to remove
    unsigned int shrinkcount;

public:

	CurvatureSearchScope(): lastderiv(NULL), curvature(NULL){}
	CurvatureSearchScope(const CurvatureSearchScope&);
	~CurvatureSearchScope() {delete[](lastderiv);}
	CurvatureSearchScope& operator=(const CurvatureSearchScope&);

	void init(double initalalpha, unsigned int in_size);

	double updateAscent(const double &value, double * guess,  double const * const deriv); // returns a safe bound on the parameter changes;
	double updateDescent(const double &value, double * guess , double const * const deriv); // returns a safe bound on the parameter changes;

    void show(FILE* f = stdout, int level=0) const{fprintf(f,"Curvature scope in %i dimentions, last value %e\t expected difference %e, surrent log-alpha %e\n",size,lastvalue,expdiff,log_alpha);}
};

template <class C>
class Vector{
	template<int supersize> Vector<C> fourierTransform_routine(const Vector<complex>& bluewindow) const;
	template<int supersize> Vector<C> invfourierTransform_routine(const Vector<complex>& bluewindow) const;
    void fourierTransform_routine(Vector<C> & fout, const complex * bluewindow, unsigned char mag) const;
    void invfourierTransform_routine(Vector<C> & fout,const complex * bluewindow, unsigned char mag) const;

public:
int asize;
C* darray;

static const bool IsPOD = false;
static const bool NeedsAddLink = ExCo<C>::NeedsAddLink; // containers needs to update addresses in link registers//typedef YESNO<true> IsPOD;

typedef DataGrid<C,2> TRJU_TYPE;
typedef DataGrid<C,2> LMUL_TYPE;
typedef C INNER_TYPE;
typedef unsigned int ITERATOR_TYPE;

typedef Matrix<C> OUTER_TYPE;
	
ExCoMeMdEcLaRe( Vector<C> )

Vector(void* owner = NULL);
Vector(C* data, int size);
Vector(const Vector<C> & v){setSize(v.getSize()); for(unsigned int i=0;i<v.getSize();i++) darray[i] = v[i];}
template<class OC> Vector(const Vector<OC> & v);

~Vector();

template<class O, unsigned int Osize> Vector(const Tuple<O,Osize> &);

void setLinkMemory(void* _newmem);
void* getLinkMemory() const;

Vector<C>& toZero(){unsigned int s = getSize();for(unsigned int i=0;i<s;i++) ExOp::toZero(darray[i]); return(*this);}
Vector<C>& toOne(){unsigned int s = getSize();for(unsigned int i=0;i<s;i++) ExOp::toOne(darray[i]); return(*this);}
Vector<C>& toRand(){unsigned int s = getSize();for(unsigned int i=0;i<s;i++) ExOp::toRand(darray[i]); return(*this);}
Vector<C>& toMin(){unsigned int s = getSize();for(unsigned int i=0;i<s;i++) ExOp::toMin(darray[i]); return(*this);}
Vector<C>& toMax(){unsigned int s = getSize();for(unsigned int i=0;i<s;i++) ExOp::toMax(darray[i]); return(*this);}

unsigned int allocatedSize() const;

Vector<C>& operator=(const Vector<C> & v);
template<class OC> Vector<C>& operator=(const Vector<OC> & v);

Vector<C>& memmove(Vector<C>& source);


template<class D> Vector(Vector<D> const & other,void* owner = NULL);

template<class A_1> void operator() (Oper1<A_1> const & op); // not a match
void operator() (Oper1< C> const & op); // match

template<class A_1, class A_2, class C_I> void operator() (Oper2<A_1,A_2> const & op, Vector<C_I> const & _in ); // not a match
template<class C_I> void operator() (Oper2<C,C_I> const & op, Vector<C_I> const & _in); // match
template<class A_1, class A_2, class C_I> void operator() (Oper2<A_1,A_2> const & op, Vector<C_I> & _in ); // not a match
template<class C_I> void operator() (Oper2<C,C_I> const & op, Vector<C_I>  & _in); // match

template<class A_1, class A_2, class A_3, class C_2, class C_3> void operator() (Oper3<A_1,A_2,A_3> const & op, Vector<C_2> const &, Vector<C_3> const &); // not a match
template<class C_2, class C_3> void operator() (Oper3<C,C_2,C_3> const & op, Vector<C_2> const &, Vector<C_3> const &); // match
template<class A_1, class A_2, class A_3, class C_2, class C_3> void operator() (Oper3<A_1,A_2,A_3> const & op, Vector<C_2>  &, Vector<C_3> const &); // not a match
template<class C_2, class C_3> void operator() (Oper3<C,C_2,C_3> const & op, Vector<C_2> &, Vector<C_3> const &); // match
template<class A_1, class A_2, class A_3, class C_2, class C_3> void operator() (Oper3<A_1,A_2,A_3> const & op, Vector<C_2>  &, Vector<C_3>  &); // not a match
template<class C_2, class C_3> void operator() (Oper3<C,C_2,C_3> const & op, Vector<C_2>  &, Vector<C_3> &); // match


C& push_back(const C& entry) {return (push_back() = entry);}
C& push_back();

void pop_back();
void pop_swap(unsigned int a);

//	void up_alloc();
//	void down_alloc();

void clear();
unsigned int getSize() const;unsigned int size() const;
void setSize(unsigned int nsize); 
void DownSize(unsigned int nsize); // conserve the first elems;
C* begin() const;
C* end() const;
C* last() const;

C & operator[](int const which);
const C & operator[](int const which)const;

//template<int flag2> operator PolyThing<C2> () const;
template<int size> operator Tuple<C,size> () const;

LFHDECL_TRIVIAL_OPERATOR(+,Vector<C>)
LFHDECL_TRIVIAL_OPERATOR(-,Vector<C>)
LFHDECL_TRIVIAL_OPERATOR(*,Vector<C>)
LFHDECL_TRIVIAL_OPERATOR(/,Vector<C>)
template<class D> Vector<C>& operator+=(const Vector<D> & v);
template<class D> Vector<C>& operator+=(const D & v);
template<class D> Vector<C>& operator-=(const Vector<D> & v);
template<class D> Vector<C>& operator-=(const D & v);
template<class D> Vector<C>& operator*=(const Vector<D> & v);
template<class D> Vector<C>& operator*=(const D & v);
template<class D> Vector<C>& operator/=(const Vector<D> & v);
template<class D> Vector<C>& operator/=(const D & v);

LFH_FAULTY Vector<unsigned int> mksortindex() const; // makes a vector of indexes which access elements in a increasing order
	
void sort(); void sort_unique();
bool issorted();
void sort_decr(); // Decreasing Order
bool issorted_decr(); // Decreasing Order
void reverse();
void random_permute();
template<class D> void getIntersection(const IntervalSet<D> & query , Vector<C> &out) const;

	static Vector<complex> bluesteinWindow(unsigned int size);
	static Vector<complex> bluesteinWindow_new(unsigned int size);
//	static void bluesteinWindow(complex *&, unsigned int size);
	Vector<C> fourierTransform(const Vector<complex>& bluewindow) const;
	Vector<C> invfourierTransform(const Vector<complex>& bluewindow) const;

	void fourierTransform_semiroutine(const Vector<complex>& bluewindow);
	void invfourierTransform_semiroutine(const Vector<complex>& bluewindow);

Matrix<C> mkouterprod() const;

void GP_Covariance(DataGrid<double, 2> &f_out, double (*metric)(const C&, const C&), double noise_variance =0.0f) const;
void GP_Covariance_cross(DataGrid<double, 2> &f_out, double (*metric)(const C&, const C&), Vector<C>& query) const;
void GP_fromCorrel(DataGrid<double, 2> &f_out, double (*correl)(const C&, const C&)) const;
void GP_fromCorrel_cross(DataGrid<double, 2> &f_out, Vector<C>& query, double (*correl)(const C&, const C&)) const;
void GP_fromCorrel_complete(DataGrid<double, 2> &f_out, Vector<C>& query, double (*correl)(const C&, const C&)) const;

};

// class which stores structured data
// data may be owned or not
// data may have a addresslink
template <class C, unsigned int nbdim>
class DataGrid : public ConstGrid<C,nbdim>{
	public:
	static const unsigned int type_dimentions = nbdim;
	typedef C INNER_TYPE;
	typedef typename MT_IFTYPE<nbdim-1, Tuple<unsigned int, nbdim> , unsigned int >::TYPE ITERATOR_TYPE;
	static const bool IsPOD = false;
	static const bool NeedsAddLink = false; // containers needs to update addresses in link registers
    template <class S> class SUBS_INNER{public: typedef DataGrid<S, nbdim> type;};

	class KeyIterator : public AbstractKeyIterator< Tuple<unsigned int, nbdim> , DataGrid<C, nbdim> >{
	public:
		KeyIterator(const DataGrid<C, nbdim>& tar) :AbstractKeyIterator< Tuple<unsigned int, nbdim> , DataGrid<C, nbdim> >(tar){}
		bool first(){ExOp::toZero((*this).curkey); return((*this).target.data != NULL);}
		bool last(){for(unsigned int i=0;i<nbdim;i++) (*this).curkey[i] = (*this).target.dims[i]-1;return((*(Vector<C>*)this).getSize() != 0);}
		bool next(){
			unsigned int dir=0;
			for(dir=0;dir< nbdim;dir++) if ((*this).curkey[dir] == (*this).target.dims[dir]-1) (*this).curkey[dir] =0; else {(*this).curkey[dir]++; break;}
			return(dir < nbdim);
		}
		bool prev(){
			unsigned int dir=0;
			for(dir=0;dir< nbdim;dir++) if ((*this).curkey[dir] == 0) (*this).curkey[dir] =(*this).target.dims[dir]-1; else {(*this).curkey[dir]--; break;}
			return(dir < nbdim);
		}
		Tuple<unsigned int, nbdim+1> padwith(const unsigned int pad) const{
			Tuple<unsigned int, nbdim+1> fout;
			memcpy(&(fout[1]), (*this).curkey.data, sizeof(unsigned int) * nbdim);
			fout[0] = pad;
			return(fout);
		}
		void write_to(unsigned int* target) const {memcpy(target, (*this).curkey.data, sizeof(unsigned int) * nbdim);}
	};
	
	
    ExCoMeMdEcLaRe( LFHCONCAT2(DataGrid<C,nbdim>) );
		C* data;
		unsigned int dims[nbdim];
		DataGrid(void* owner = NULL);
		DataGrid(char* path);
		DataGrid(unsigned int* dims, void* owner = NULL);
		DataGrid(const DataGrid<C,nbdim>&);
	template<unsigned int fsize> DataGrid(const DataGrid<Tuple<C,fsize> ,nbdim-1>&);
	template<unsigned int fsize> DataGrid(const DataGrid< C[fsize]  ,nbdim-1>&);
	template<unsigned int fsize> DataGrid(const DataGrid< DataGrid<C,fsize> , nbdim-fsize >&);
	~DataGrid();

    DataGrid<C, nbdim>& toZero();
    DataGrid<C, nbdim>& toOne();
    DataGrid<C, nbdim>& toRand();
    DataGrid<C,nbdim>& memmove(DataGrid<C,nbdim>& source);
	DataGrid<C,nbdim-1> makeSlice(int which_dim, int coor)const;
	DataGrid<C,nbdim+1> fromSlice(int which_dim) const;
    DataGrid<C,nbdim> mkconcatenate(const DataGrid<C,nbdim> &, unsigned int direction)const;
	DataGrid<C,nbdim> selectedSlices(unsigned int direction, Vector<unsigned int> coor); // TODO
    template<unsigned int fsize> DataGrid<Tuple<C,fsize>,nbdim-1> mkTupleGrid(const Tuple<unsigned int,fsize> &which) const;

	
	

	KeyIterator getKeyIterator() const{return(KeyIterator(*this));}

//operators
DataGrid<C,nbdim>& operator=(const DataGrid<C,nbdim>& other);
template<class O> DataGrid<C,nbdim>& operator=(const DataGrid<O,nbdim>& other);

DataGrid<C,nbdim> mkOrthoFlip(unsigned int which_dim_source, unsigned int which_dim_target) const;
const DataGrid<C,nbdim>& toOrthoFlip(unsigned int which_dim_source, unsigned int which_dim_target){DataGrid<C,nbdim> tmp = this->mkOrthoFlip(which_dim_source, which_dim_target);(*this) = tmp; return (*this);}
DataGrid<C,nbdim> mkDimFlip(unsigned int which_dim) const;
const DataGrid<C,nbdim>& toDimFlip(unsigned int which_dim);

unsigned int totsize() const; // total size of array

C& operator[](unsigned int);
C& operator[](int);
C operator[](unsigned int) const;
C operator[](int) const;

C& operator()(unsigned int*);
C& operator()(int*);
C operator()(unsigned int*) const;
C operator()(int*) const;

typename MT_IFTYPE<nbdim-1, DataGrid<C,nbdim-1>,C>::TYPE operator()(unsigned int) const;
typename MT_IFTYPE<nbdim-2, DataGrid<C,nbdim-2>,C>::TYPE operator()(unsigned int,unsigned int) const;
typename MT_IFTYPE<nbdim-3, DataGrid<C,nbdim-3>,C>::TYPE operator()(unsigned int,unsigned int,unsigned int) const;

C linearInpterpolation(const Tuple<double,nbdim>&) const;

C& operator()(const Tuple<unsigned int,nbdim>& );
C& operator()(const Tuple<int,nbdim>&);
C operator()(const Tuple<unsigned int,nbdim>&) const; // constgrid func
C operator()(const Tuple<int,nbdim>&) const; // constgrid func
void getDims(Tuple<unsigned int, nbdim> &o_dims) const; // constgrid func

void drawLine(const Tuple<unsigned int,nbdim> &a, const Tuple<unsigned int,nbdim> &b, const C&);
template<unsigned int size> void drawLine(const Tuple<unsigned int,nbdim-1> &a, const Tuple<unsigned int,nbdim-1> &b, const Tuple<C,size> &);

void drawChar(unsigned char which, const Tuple<unsigned int,nbdim> &where, const C &what);
void drawCharFlip(unsigned char which, const Tuple<unsigned int,nbdim> &where, const C &what);
void drawUpChar(unsigned char which, const Tuple<unsigned int,nbdim> &where, const C &what);
void drawDownChar(unsigned char which, const Tuple<unsigned int,nbdim> &where, const C &what);
template<unsigned int size> void drawChar(unsigned char which, const Tuple<unsigned int,nbdim-1> &where, const Tuple<C,size> &what);
template<unsigned int size> void drawCharFlip(unsigned char which, const Tuple<unsigned int,nbdim-1> &where, const Tuple<C,size> &what);
template<unsigned int size> void drawUpChar(unsigned char which, const Tuple<unsigned int,nbdim-1> &where, const Tuple<C,size> &what);
template<unsigned int size> void drawDownChar(unsigned char which, const Tuple<unsigned int,nbdim-1> &where, const Tuple<C,size> &what);

Tuple<unsigned int, nbdim> addressof(C const * const item) const;

template <int ssize> C* operator()(const Tuple<unsigned int,ssize>&);
template <int ssize> C* operator()(const Tuple<int,ssize>&);
template <int ssize> const C * const operator()(const Tuple<unsigned int,ssize>&) const;
template <int ssize> const C * const operator()(const Tuple<int,ssize>&) const;

void setRow(const C* const,const Tuple<unsigned int, nbdim>& coor, const int direction);

template<class O> DataGrid<C,nbdim>& operator+=(KeyElem< Tuple<double,nbdim>, O > const & other);

template<class O> DataGrid<C,nbdim>& operator+=(DataGrid<O,nbdim> const & other);
template<class O> DataGrid<C,nbdim>& operator-=(DataGrid<O,nbdim>  const & other);
template<class O> DataGrid<C,nbdim>& operator*=(DataGrid<O,nbdim>  const & other);
template<class O> DataGrid<C,nbdim>& operator/=(DataGrid<O,nbdim>  const & other);
template<class O> DataGrid<C,nbdim>& operator+=(O const & other);
template<class O> DataGrid<C,nbdim>& operator-=(O const & other);
template<class O> DataGrid<C,nbdim>& operator*=(O const & other);
template<class O> DataGrid<C,nbdim>& operator/=(O const & other);

//////////////////////////////
// custom constructors
DataGrid<C,nbdim> Crop(const Tuple<unsigned int, nbdim> &min, const Tuple<unsigned int, nbdim> &max) const;



DataGrid<C,2> pseudoInverse(DataGrid<C,nbdim>*L = NULL, DataGrid<C,nbdim>*R = NULL) const;// return pseudo inverse, if factors R and/or L are provided, returns L * ((this)^-1) * R
DataGrid<C,nbdim> Inverse() const; // assumes the Matrix is inversible!!!

void initeye();

void update_pseudoInverse(const DataGrid<C,nbdim>& target);
void update_Inverse(const DataGrid<C,nbdim>& target);// assumes the Matrix is still inversible!!!

		void makeFreqConvolutionMap(int size, double low, double high);

		template <class D> DataGrid(DataGrid<D,nbdim> const & other);

	void setSizes(unsigned int const * const dims);
	unsigned int const * const getSizes() const;
	void setSizes(const Tuple<unsigned int,nbdim> &dims);
//	Tuple<unsigned int,nbdim> getSizes() const;

    DataGrid< Tuple<C, nbdim> ,nbdim> makeGradient() const;
    DataGrid<double,nbdim> makeGradientNorm() const;

	
    DataGrid<unsigned int,nbdim> makeWatershedSortIndex(unsigned int &nbpix, bool min_gradient_bassin = true, bool tree_skip_save = true) const;
	DataGrid<unsigned int, nbdim> makeSlicesIndexes(const unsigned int nbslices) const{Vector<double> daslices; for(unsigned int i=0;i < nbslices-1;i++) daslices.push_back(((double)i+1)/nbslices); return makeSlicesIndexes(daslices);}
	DataGrid<unsigned int, nbdim> makeSlicesIndexes(const Vector<double> &fractions)const;

		template<class D> void matrixleftmultiply(const ConstGrid<D,2>& R);
		void settomatrixMultiply(const DataGrid<C,nbdim>& L, const DataGrid<C,nbdim>& R);
		void settomatrixMultiplyofTransposed(const DataGrid<C,nbdim>& L, const DataGrid<C,nbdim>& Transposed_R);

		void leftHouseHolderMultiply_lame(const double * const vec, int lenght,bool hint);
		void rightHouseHolderMultiply_lame(const double * const vec, int lenght,bool hint);
		void leftHouseHolderMultiply_lame(const double * const vec, int lenght, const double& sqrt_den,bool hint);
		void rightHouseHolderMultiply_lame(const double * const vec, int lenght, const double& sqrt_den,bool hint);
		void leftHouseHolderMultiply(const double * const vec, int lenght, const double& sqrt_den,bool hint);
		void rightHouseHolderMultiply(const double * const vec, int lenght, const double& sqrt_den,bool hint);
		void dualHouseHolderMultiply(const double * const vec, int lenght, const double& sqrt_den,bool hint); // multiply on both sides, assume squarre symetric Matrix!

		DataGrid<C,nbdim>& convolvecircle_zeroborder(double radius);

		const DataGrid<C,nbdim>& blur(const GaussianDistribution<nbdim> &);
		const DataGrid<C,nbdim>& blur_zeroborder(const GaussianDistribution<nbdim> &);
		const DataGrid<C,nbdim>& blur_crude(const GaussianDistribution<nbdim> &);

		template<class D> void leftHouseHolderMultiply(const D * const vec, int lenght, const D& sqrt_den,bool hint);
		template<class D> void rightHouseHolderMultiply(const D * const vec, int lenght, const D& sqrt_den,bool hint);

		// row columns operations

		void leftoffdiagelimination(double factor, int from, int to);// add collumns (top) to // add collumns
		void rightoffdiagelimination(double factor, int from, int to);

		void positivedefinite_inverse_and_involution(DataGrid<C, nbdim> &f_out, const DataGrid<C, nbdim> &invol);
		void solveLinearSystem(const C* const target, C * _out);
		void LinearTransform(const C* const target, C * _out);

		// this computes (R*M*(R^-1))^-1, where M is a large Matrix, and R a reduction
		void downsampled_Matrix(const DataGrid<double, nbdim> &f_inner,const DataGrid<double, nbdim> &f_reduction);

		void Matrix_factor_out_from_symetric_Matrix(DataGrid<C, nbdim> &io_symmetric_Matrix);

	void makeDiagonalizer(const ConstGrid<double,2> &symmetric_Matrix, Vector<double> &eigenval);
	void makeDiagonalizer_ofinverse(const ConstGrid<double,2> &symmetric_Matrix, Vector<double> &eigenval);

	template<unsigned int size> void makeDiagonalizer(const TMatrix<double,size,size> &symmetric_Matrix, Vector<double> &eigenval);
	template<unsigned int size> void makeDiagonalizer_ofinverse(const TMatrix<double,size,size> &symmetric_Matrix, Vector<double> &eigenval);

		void solve_sym_and_project(DataGrid<double, 2> &project, const ConstGrid<double,2> &symmetric_Matrix); // TODO // computes : project dot symm_Matrix^-1 dot (this)

		void makeInverse(const ConstGrid<double,2> &symmetric_Matrix);

		void LinearSystem_fixpoint_solve(const DataGrid<C, nbdim>&target, double (*weights)(int a, int b));
		void LinearSystem_residual(const DataGrid<C, nbdim>&guess, const DataGrid<C, nbdim>&target, const ConstGrid<double,2>&);

		DataGrid<int, nbdim> directNeightbor_Count_Compare(SetComparison (*query)(const C&, const C&), bool equalitycount = false);
		DataGrid<int, nbdim> indirectNeightbor_Count_Compare(SetComparison (*query)(const C&, const C&), bool equalitycount = false);
		DataGrid<int, nbdim> knightNeightbor_Count_Compare(SetComparison (*query)(const C&, const C&), bool equalitycount = false);

		vector< Tuple<unsigned int,nbdim > > getCoor_Mathching(const C&)const;

		unsigned int get_directNeightbor(const Tuple<unsigned int , nbdim> & query, Tuple< Tuple<unsigned int , nbdim>,  nbdim * 2> &_out )const;
		unsigned int get_indirectNeightbor(const Tuple<unsigned int , nbdim> & query, Tuple< Tuple<unsigned int , nbdim>,  (unsigned int)(TEMPLATE_INT_POWER<3,nbdim>::ans) -1u > &_out )const;
		unsigned int get_knightNeightbor(const Tuple<unsigned int , nbdim> & query, Tuple< Tuple<unsigned int , nbdim>, nbdim * (nbdim -1) * 4  > &_out )const;

		unsigned int getDimSize(unsigned int dim);
		void getDims(unsigned int* out_dim);
		void exponentialBlur(double std_dev, bool circular);
		template<unsigned int size, Tuple_flag Cflag> DataGrid<Tuple<C, size,Cflag>, nbdim> exponentialBlur(Tuple<double,size,Cflag> apertures, bool circular);
		template<unsigned int size, Tuple_flag Cflag> DataGrid<Tuple<C, size,Cflag>, nbdim> crudeGaussianBlur(Tuple<double,size,Cflag> apertures, bool circular);
		DataGrid<C, nbdim> crudeGaussianBlur(double aperture, bool circular = false);

		DataGrid<C,nbdim>& toresize(unsigned int* newdims);
		DataGrid<C,nbdim> resize(unsigned int* newdims) const;
		DataGrid<C,nbdim> resize_dir(unsigned int dims, unsigned int size) const;

		DataGrid<C,nbdim>& toresize_crude(unsigned int* newdims);
		DataGrid<C,nbdim> resize_crude(unsigned int* newdims) const;
		DataGrid<C,nbdim> resize_dir_crude(unsigned int dims, unsigned int size) const;
	

DataGrid<double ,nbdim>  ExpectedDistanceToOther(double exit_prob) const;


template<class ARR_double> void ExpectedDistanceToOther(const DataGrid< ARR_double ,nbdim> &source, const DataGrid<double, 2>& transition);
template<class ARR_double> void ExpectedDistanceToOther_instate(const DataGrid< ARR_double ,nbdim> &source,const DataGrid<double, 2>& transition, int state);

void ExpectedDistanceToOther_I(const DataGrid< double ,nbdim> &source, const DataGrid<double, 2>& transition);
void ExpectedDistanceToOther_instate_I(const DataGrid< double ,nbdim> &source,const DataGrid<double, 2>& transition, int state);

void ExpectedDistanceToOther(const DataGrid< double ,nbdim+1> &source, const DataGrid<double, 2>& transition, bool state_zero_special = false);
void ExpectedDistanceToOther_instate(const DataGrid< double ,nbdim+1> &source,const DataGrid<double, 2>& transition, int state);


Vector< Tuple<unsigned int, nbdim> > get_ordered_Data_Coors() const;
Vector< KeyElem<C, Tuple<unsigned int, nbdim> > > get_ordered_Data_Coors_Tagged() const;

DataGrid< Tuple<unsigned int, nbdim> ,nbdim> LocalMaxMap() const;
DataGrid< Tuple<unsigned int, nbdim> ,nbdim> LocalMinMap() const;

DataGrid< unsigned int ,nbdim> SegementClimbToMax(unsigned int &nbsegs, bool dist_penal = false, bool knight = false , bool eq_join = true, bool mark_da_maxima = false, const DataGrid< bool ,nbdim>* ignore_filter = NULL) const;
DataGrid< unsigned int ,nbdim> SegementClimbToMaxMax(unsigned int &nbsegs, bool eq_join = true) const;
DataGrid< unsigned int ,nbdim> SegementClimbToMin(unsigned int &nbsegs, bool dist_penal = false, bool knight = false , bool eq_join = true, bool mark_da_maxima = false, const DataGrid< bool ,nbdim>* ignore_filter = NULL) const;

DataGrid<C,nbdim> FFtransform_old() const;
DataGrid<C,nbdim> invFFtransform_old() const;
DataGrid<typename ExCo<C>::COMPLEX_TYPE ,nbdim> FFtransform() const;
DataGrid<typename ExCo<C>::COMPLEX_TYPE ,nbdim> invFFtransform() const;

	DataGrid<typename ExCo<C>::COMPLEX_TYPE ,nbdim> FFtransform_2(Bluestein* blue = NULL, typename ExCo<C>::COMPLEX_TYPE * buffer= NULL) const;
	DataGrid<typename ExCo<C>::COMPLEX_TYPE ,nbdim> invFFtransform_2(Bluestein* blue = NULL, typename ExCo<C>::COMPLEX_TYPE * buffer= NULL) const;

    // emulate a transform on a twice bigger array, where F[x] = F[-x] so that F'[w] = F'[-w], hence there is no need to have the
    DataGrid<typename ExCo<C>::COMPLEX_TYPE ,nbdim> FFtransform_even(Bluestein* blue = NULL, typename ExCo<C>::COMPLEX_TYPE * buffer= NULL) const;
	DataGrid<typename ExCo<C>::COMPLEX_TYPE ,nbdim> invFFtransform_even(Bluestein* blue = NULL, typename ExCo<C>::COMPLEX_TYPE * buffer= NULL) const;

	DataGrid<typename ExCo<C>::COMPLEX_TYPE ,nbdim> FFtransform_zeroborder(Bluestein* blue = NULL, typename ExCo<C>::COMPLEX_TYPE * buffer= NULL) const;
	DataGrid<typename ExCo<C>::COMPLEX_TYPE ,nbdim> invFFtransform_zeroborder(Bluestein* blue = NULL, typename ExCo<C>::COMPLEX_TYPE * buffer= NULL) const;

void applycircleFilter(double minx, double maxx);
DataGrid<C,nbdim> circleFilter(double minx, double maxx);

DataGrid<unsigned int, nbdim> LabelConnected( bool (*isToLabel)(const C&), unsigned int *out_nblabels = NULL, Vector< Tuple<unsigned int, nbdim*2> >* out_boundingrect= NULL);

template<class A_1, class A_2, class C_2, unsigned int size_2> void operator() (Oper2<A_1,A_2> const & op, DataGrid<C_2,size_2> const & _in ); // not a match
template<class C_2, unsigned int size_2> void operator() (Oper2<C,C_2> const & op, DataGrid<C_2,size_2> const & _in); // match
template<class A_1, class A_2, class C_2, unsigned int size_2> void operator() (Oper2<A_1,A_2> const & op, DataGrid<C_2,size_2> & _in ); // not a match
template<class C_2, unsigned int size_2> void operator() (Oper2<C,C_2> const & op, DataGrid<C_2,size_2> & _in); // match

template<class A_1, class A_2, class A_3, class C_2, unsigned int nbdim_2, class C_3, unsigned int nbdim_3> void operator() (Oper3<A_1,A_2,A_3> const & op, DataGrid<C_2, nbdim_2>  &, DataGrid<C_3, nbdim_3>  &); // not a match
template<class C_2, unsigned int nbdim_2, class C_3, unsigned int nbdim_3> void operator() (Oper3<C,C_2,C_3> const & op, DataGrid<C_2, nbdim_2>  &, DataGrid<C_3, nbdim_3> &); // match

class FiniteDifference:  LFHDECL_OPER2(LFHCONCAT2(DataGrid<C, nbdim>),LFHCONCAT2(DataGrid<C, nbdim>));
class FiniteSum:  LFHDECL_OPER2(LFHCONCAT2(DataGrid<C, nbdim>),LFHCONCAT2(DataGrid<C, nbdim>));

class Convolution : LFHDECL_OPER3(LFHCONCAT2(DataGrid<C, nbdim>),LFHCONCAT2(DataGrid<C, nbdim>), LFHCONCAT2(DataGrid<double, nbdim>));
template<class D> class RankMap:  LFHDECL_OPER2(LFHCONCAT2(DataGrid<C, nbdim>),LFHCONCAT2(DataGrid<D, nbdim>));

}; // end of class DataGrid

template<class C, unsigned int SIZE = 0u>
class Trianglix{
    void HouseHolderMultiply(const C * const vec, double denum2, unsigned int  lenght, C * buffer, bool hint); // * (I + vt^H)
    void offdiagelimination_down(const C &fact,unsigned int col, C * buffer);
    void offdiagelimination_up(const C &fact,unsigned int col, C * buffer);
public:
    static const unsigned int totsize = TEMPLATE_TRIANGLE_NUMBER<SIZE,2>::ans;
	static const bool IsPOD = ExCo<C>::IsPOD;
    C data[TEMPLATE_TRIANGLE_NUMBER<SIZE,2>::ans];
    Trianglix(){}
	Trianglix(const Tuple<C, SIZE>& other);
    Trianglix(const Trianglix<C, SIZE> &other){for(unsigned int i=0;i<totsize;i++) data[i] = other.data[i];}
    unsigned int getSize()const{return SIZE;}

	void CholeskyStep_up(const C * const vec, unsigned int length, C* buf); // L * T *L' 

bool isValid()const{for(unsigned int i =0; i < totsize;i++) if (!ExOp::isValid(data[i])) return false; return true;}
 
	C maxEigenValue() const; // finds the maximum eigen value, uses the "Power iteration" method

operator TMatrix<C, SIZE, SIZE> () const{ TMatrix<C, SIZE, SIZE> fout; unsigned int i,j,k; for(i=0,j=0,k=2;i<totsize;i++) if (i == j) { fout.data[(k-2) * (SIZE +1)] = data[i]; j += k++;}else {fout.data[(k-2) + (k-2-j+i) * SIZE] = data[i]; fout.data[(k-2-j+i) + (k-2) * SIZE] = ExOp::mktrju(data[i]);} return fout; }

//	Trianglix<C, SIZE>& makeNonSingular(double factor); // assumes all eigen values are positive, then makes the eigen value converge to the maximum eigne value, depending on factor (0 does nothing, 1 all eigen values matches (diagonal matrix)


	double WishartLogDensity(const Trianglix<C, SIZE>&, unsigned int nbsamples) const;



    template<class O> Trianglix<C, SIZE>& operator=( const Trianglix<O,SIZE>& other){ for(unsigned int i=0;i<totsize;i++) data[i] = other.data[i]; return(*this);}
    template<class O> Trianglix<C, SIZE>& operator+=( const Trianglix<O,SIZE>& other){ for(unsigned int i=0;i<totsize;i++) data[i] += other.data[i]; return(*this);}
    template<class O> Trianglix<C, SIZE>& operator-=( const Trianglix<O,SIZE>& other){ for(unsigned int i=0;i<totsize;i++) data[i] -= other.data[i]; return(*this);}
	Trianglix<C, SIZE>& operator*=( const double& other) {for(unsigned int i=0;i<totsize;i++) data[i] *= other; return(*this);}
	Trianglix<C, SIZE>& operator/=( const double& other) {for(unsigned int i=0;i<totsize;i++) data[i] /= other; return(*this);}
    

	template<class O> Trianglix<C, SIZE> operator/( const O& other)const{return (Trianglix<C, SIZE>(*this) /= other);}
    template<class O> Trianglix<C, SIZE> operator*( const O& other)const{return (Trianglix<C, SIZE>(*this) *= other);}
    template<class O> Trianglix<C, SIZE> operator+( const Trianglix<O, SIZE>& other)const{ return (Trianglix<C, SIZE>(*this) += other);}
    template<class O> Trianglix<C, SIZE> operator-( const Trianglix<O, SIZE>& other)const{ return (Trianglix<C, SIZE>(*this) -= other);}

    Trianglix<C, SIZE>& toZero(){ unsigned int i; for(i=0;i<totsize;i++) ExOp::toZero(data[i]); return(*this);}
    Trianglix<C, SIZE>& toOne(){ unsigned int i,j,k; for(i=0,j=0,k=2;i<totsize;i++) if (i == j) {ExOp::toOne(data[i]); j += k++;}else ExOp::toZero(data[i]);}
    Trianglix<C, SIZE>& toRand(){ unsigned int i,j,k; for(i=0,j=0,k=2;i<totsize;i++) if (i == j) {ExOp::toRand(data[i]); data[i] = ExOp::mkrealproj(data[i]); j += k++;}else ExOp::toRand(data[i]);}

    Trianglix<C,SIZE>& operator=(const Tuple<C, SIZE>& other);
    Trianglix<C, SIZE> inverse() const;


	C trace_of_division(const Trianglix<C,SIZE> &divisor) const;


	C& cell(unsigned int x, unsigned int y){return data[ (x>= y) ? y + ((x * (x+1)) /2) : x + ((y * (y+1)) /2)];}
	C cell(unsigned int x, unsigned int y) const {return data[ (x>= y) ? y + ((x * (x+1)) /2) : x + ((y * (y+1)) /2)];}

	Tuple<C,SIZE> getEigenValues()const;
	Tuple<C,SIZE> getDiagonal()const{Tuple<C,SIZE> fout;  unsigned int i,j; for(i=0,j=0;j<SIZE;i+=j+1) fout[j++] = data[i]; return(fout);}
    C determinant() const;
    C log_determinant() const;

  //  C determinant_singularguard() const; // cannot be zero but if all eigneval are 0, forces eigenval to the geometric mean of not-zero eignevals
  //  C log_determinant_singularguard() const; 
    C Xformed_inner_product_of_inverse(const Tuple<C, SIZE>& other) const;
    C inv_Xformed_inner_product_singularguard( const Tuple<C, SIZE>& other, double guard_fraction = 1.0f, double *log_det = NULL) const;


    TMatrix<C,SIZE,SIZE> diagonalizer_of_inverse() const;

    double bhattacharryya(const Tuple<C,SIZE>&dev ,const Trianglix<C,SIZE>& other)const{return bhattacharryya_partial(dev,other) - 0.25f * (log(determinant())+ log(other.determinant())); }
	double bhattacharryya_partial(const Tuple<C,SIZE>& ,const Trianglix<C,SIZE>&)const;

	template< unsigned int SIZE2> C Xformed_inner_product( const Tuple<C, SIZE2>& other) const;
	C Xformed_inner_product( const Tuple<C, 0u>& other) const;
	
	double weighted_bhattacharryya(const Tuple<C,SIZE>&dev ,const Trianglix<C,SIZE>& other, double this_weight, double other_weight)const{return weighted_battacharya_partial(dev,other,this_weight,other_weight) - 0.5f * (this_weight*log(determinant())+ other_weight*log(other.determinant())); }
	double weighted_bhattacharryya_partial(const Tuple<C,SIZE>& ,const Trianglix<C,SIZE>&, double this_weight, double other_weight)const;
	double weighted_bhattacharryya_partial_MK2(const Tuple<C,SIZE>& ,const Trianglix<C,SIZE>&, double this_weight, double other_weight)const;




    void show(FILE* f = stdout, int level= 0)const;
    string type_tostring()const;

void save(FILE* f)const{fwrite(data,sizeof(C),TEMPLATE_TRIANGLE_NUMBER<SIZE,2>::ans,f);}
void load(FILE* f,unsigned int lenght =0) {fread(data,sizeof(C),TEMPLATE_TRIANGLE_NUMBER<SIZE,2>::ans,f);}

};

template<class C>
class Trianglix<C, 0u>{ // square, symetric matrix (anti-symmetric in imaginary), has special mulitiplation to ensure symmetry (multiply other to the left and right at the same time)
    void HouseHolderMultiply(const C * const vec, double denum2, unsigned int  lenght, C * buffer, bool hint); // * (I + vt^H)
    void offdiagelimination_down(const C &fact,unsigned int col, C * buffer);
    void offdiagelimination_up(const C &fact,unsigned int col, C * buffer);
    void offdiagelimination_up_backroutine(const C &fact,unsigned int col, C * buffer); // assumes that the row of interest is triagonal, so some values are assumed to be zero

    public:
	static const bool IsPOD = ExCo<C>::IsPOD;
    void QR_back(const C &factC,const C &factS,unsigned int col);

    typedef YESNO< ExCo<C>::IS_COMMUTATIVE::ans > IS_COMMUTATIVE;
    C* data;
    unsigned int t_size;
    Trianglix(): t_size(0){}

	C maxEigenValue() const; // finds the maximum eigen value, uses the "Power iteration" method
    unsigned int getSize()const{return t_size;}

    Trianglix(const Tuple<C, 0u>& other); 
    template<unsigned int SIZE> Trianglix(const Tuple<C, SIZE>& other);

    Trianglix(const Trianglix<C, 0u>& other);
    template<unsigned int SIZE> Trianglix(const Trianglix<C, SIZE>& other);

    ~Trianglix();
	Trianglix<C, 0u>& operator=(const Trianglix<C, 0u>&other){setSize(other.t_size); unsigned int i,ts; for(i=0,ts=totsize();i<ts;i++) data[i] = other.data[i]; return(*this); }
    
	bool isValid()const{if (t_size == 0) return(true); if (data == NULL) return(false); for(unsigned int i =0; i < totsize();i++) if (!ExOp::isValid(data[i])) return false; return true;  }
	
	template<unsigned int SIZE> Trianglix<C, 0u>& operator=(const Trianglix<C,SIZE>&other){setSize(SIZE); unsigned int i,ts; for(i=0,ts=totsize();i<ts;i++) data[i] = other.data[i]; return(*this); }
    Trianglix<C, 0u>& memmove(Trianglix<C, 0u>& other);
	
    unsigned int totsize()const{return (t_size & 1) ? t_size * ((t_size>>1)+1) : (t_size>>1) * (t_size+1);}
    void setSize(unsigned int s);

	template<class O> Trianglix<C, 0u> operator+( const O& other)const{ return (Trianglix<C, 0u>(*this) += other);}
    template<class O> Trianglix<C, 0u> operator-( const O& other)const{ return (Trianglix<C, 0u>(*this) -= other);}
	template<class O> Trianglix<C, 0u> operator*( const O& other)const{ return (Trianglix<C, 0u>(*this) *= other);}
    template<class O> Trianglix<C, 0u> operator/( const O& other)const{ return (Trianglix<C, 0u>(*this) /= other);}
	
    template<class O> Trianglix<C, 0u> operator+( const Trianglix<O, 0u>& other)const{ return (Trianglix<C, 0u>(*this) += other);}
    template<class O> Trianglix<C, 0u> operator-( const Trianglix<O, 0u>& other)const{ return (Trianglix<C, 0u>(*this) -= other);}
    template<class O> Trianglix<C, 0u>& operator+=( const Trianglix<O, 0u>& other){unsigned int i,ts; for(i=0,ts=totsize();i<ts;i++) data[i] += other.data[i]; return(*this);}
    template<class O> Trianglix<C, 0u>& operator-=( const Trianglix<O, 0u>& other){unsigned int i,ts; for(i=0,ts=totsize();i<ts;i++) data[i] -= other.data[i]; return(*this);}

    Trianglix<C, 0u>& operator*=( const double& other){unsigned int i,ts; for(i=0,ts=totsize();i<ts;i++) data[i] *= other; return(*this);}
    Trianglix<C, 0u>& operator/=( const double& other){unsigned int i,ts; for(i=0,ts=totsize();i<ts;i++) data[i] /= other; return(*this);}

    Trianglix<C, 0u>& toZero(){ unsigned int i,ts; ts = totsize(); for(i=0;i<ts;i++) ExOp::toZero(data[i]); return(*this);}
    Trianglix<C, 0u>& toOne(){ unsigned int i,j,k,ts; ts = totsize(); for(i=0,j=0,k=2;i<ts;i++) if (i == j) {ExOp::toOne(data[i]); j += k++;}else ExOp::toZero(data[i]);}
    Trianglix<C, 0u>& toRand(){ unsigned int i,j,k,ts; ts = totsize(); for(i=0,j=0,k=2;i<ts;i++) if (i == j) {ExOp::toRand(data[i]); data[i] = ExOp::mkrealproj(data[i]);j += k++;}else ExOp::toRand(data[i]);}
    C operator()(unsigned int x, unsigned int y)const{return( (x < y) ? data[((y * (y+1)) / 2) + x]: data[((x * (x+1)) / 2) + y]);}
    void show(FILE* f = stdout, int level= 0)const;
    string type_tostring()const;
    
    template< unsigned int SIZE> C Xformed_inner_product( const Tuple<C, SIZE>& other) const;
    C Xformed_inner_product( const Tuple<C, 0u>& other) const;
	
    C Xformed_inner_product_of_inverse(const Tuple<C, 0u>& other) const;
    
    C gaussLL(const Tuple<C, 0u>& other) const;

	template< unsigned int SIZE> Tuple<C, SIZE> operator*( const Tuple<C, SIZE>& other) const;
    Tuple<C, 0u> operator*( const Tuple<C, 0u>& other) const;
	
	template< unsigned int SIZE> Tuple<C, SIZE> divisionof( const Tuple<C, SIZE>& other) const;
    Tuple<C, 0u> divisionof( const Tuple<C, 0u>& other) const;
	
	
	
	
    Trianglix<C, 0u> Xformed_outer_product( const Tuple<C, 0u>& o) const{ return(Trianglix( (*this) * o ) ); } // Kyy^TK^T
    Trianglix<C, 0u> Xformed_outer_product_of_inverse( const Tuple<C, 0u>& o) const{ return(Trianglix( this->divisionof(o) ) );} // K^{-1}yy^T(K^{-1})^T

    Trianglix<C, 0u> operator*(const Matrix<C>& other) const; // does not work... yet
    const Trianglix<C, 0u>& operator*=(const Matrix<C>& other); // does not work... yet

    template<unsigned int SIZE> Trianglix<C>& operator=(const Tuple<C, SIZE>& other);

	C trace_of_product(const Trianglix<C,0u> &other) const;
	C trace_of_division(const Trianglix<C,0u> &divisor) const;
    
    Matrix<C> makeMatrix()const;

    Tuple<C,0u> getEigenValues() const;
	Tuple<C,0u> getDiagonal()const{Tuple<C,0u> fout; fout.setSize(t_size); unsigned int i,j; for(i=0,j=0;j<t_size;i+=j+1) fout[j++] = data[i]; return(fout);}

    C determinant() const{if (t_size == 0) myexit("determinant of a size 0 matrix, check for that!\n"); Tuple<C, 0u> eigens = getEigenValues(); for(unsigned i=1;i< t_size;i++) eigens[0] *= eigens[i]; return(eigens[0]);}
    C log_determinant() const;
	C trace()const{if (t_size == 0) myexit("determinant of a size 0 matrix, check for that!\n"); C fout = data[0]; unsigned int i,j; for(i=2,j=1;j<t_size;i+=j+1) {fout += data[i];j++;} return(fout); }
	
    double bhattacharryya(const Tuple<C, 0u>&dev ,const Trianglix<C, 0u>& other)const{return battacharya_partial(dev,other) - 0.25f * (log(determinant())+ log(other.determinant())); }
    double bhattacharryya_partial(const Tuple<C, 0u>& ,const Trianglix<C, 0u>&)const;
	
	
    Trianglix<C, 0u> inverse_OLD() const;
	Trianglix<C, 0u>& toinverse() {Trianglix<C, 0u> tmp = this->mkinverse(); return ExOp::memmove(*this,tmp); }

    Trianglix<C, 0u> inverse_MK2() const;
    Trianglix<C, 0u> mkinverse() const;
	
	C& cell(unsigned int x, unsigned int y){return data[ (x>= y) ? y + ((x * (x+1)) /2) : x + ((y * (y+1)) /2)];}
	C cell(unsigned int x, unsigned int y) const {return data[ (x>= y) ? y + ((x * (x+1)) /2) : x + ((y * (y+1)) /2)];}

    void eigen() const;
};

template<class C>
class Matrix{ // same encoding, different interpretation
	template<class O, class B> void alloc_add(Matrix<O>& fout, const Matrix<B>& other)const; // assumes fout.data == NULL
	template<class O, class B> void mkmult(Matrix<O>& fout, const B& other, Tuple< char, 0> *)const;
	template<class O, class B> void mkmult(Matrix<O>& fout, const B& other, Tuple< char, 1> *)const;
	template<class O, class B> void mkmult(Matrix<O>& fout, const B& other, Tuple< char, 2> *)const;
    void InitFromTrianglix(const C* data);
public:
    typedef YESNO< false > IS_COMMUTATIVE; 
	unsigned int sizex;
	unsigned int sizey;
	C* data;
	Matrix(): data(NULL){}
	Matrix(unsigned int _sizex, unsigned int _sizey): sizex(_sizex), sizey(_sizey),data(new C[_sizex*_sizey]){}
	Matrix(const DataGrid<C,2>& other):sizex(other.dims[0]), sizey(other.dims[1]) {data = new C[sizex*sizey];}

    Matrix(const Trianglix<C,0u>& other):sizex(other.t_size), sizey(other.t_size) {data = new C[sizex*sizey];InitFromTrianglix(other.data);}
    template <unsigned int SIZE> Matrix(const Trianglix<C,SIZE>& other):sizex(SIZE), sizey(SIZE) {data = new C[sizex*sizey];InitFromTrianglix(other.data);}
    Matrix<C>& operator=(const Matrix<C>& other) {delete[](data); sizex = other.sizex; sizey = other.sizey; data = new C[sizex*sizey];for(unsigned int i=0;i<sizex*sizey;i++) data[i] = other.data[i]; }
    Matrix<C>& operator=(const Trianglix<C,0u>& other) {delete[](data); sizex = other.size; sizey = other.size; data = new C[sizex*sizey];InitFromTrianglix(other.data);}
    template <unsigned int SIZE> Matrix<C>& operator=(const Trianglix<C,SIZE>& other) {delete[](data); sizex = SIZE; sizey = SIZE; data = new C[sizex*sizey];InitFromTrianglix(other.data);}

	~Matrix(){delete[](data);}
    void setSizes(unsigned int x, unsigned int y){delete[](data); sizex = x; sizey=y; data = new C[sizex*sizey];}

	Matrix<C> operator+(const Matrix<C>& other)const;
	Matrix<C> operator-(const Matrix<C>& other)const;
	Matrix<C> operator-()const;
	Matrix<C>& operator+=(const Matrix<C>& other);
	Matrix<C>& operator-=(const Matrix<C>& other);
    Matrix<C>& toZero(){for(unsigned int i=0;i<sizex*sizey;i++) ExOp::toZero(data[i]); return(*this); }
    Matrix<C>& toRand(){for(unsigned int i=0;i<sizex*sizey;i++) ExOp::toRand(data[i]); return(*this); }


    Matrix<C> operator*(const Matrix<C>& other)const;

	Matrix<C> inverse() const;// TODO!

    double pnorm() const;
    double norm() const{return(sqrt(pnorm));}

    void show(FILE *f=stdout, int level=0)const;
};

template<class C, class D, class FUNC, unsigned int DIM>
class FuncGrid: public ConstGrid<C,DIM>{
	unsigned int start;
	unsigned int rangesize;
	public:
	Vector<D> data;
	FUNC* cl; // must implement one of:
	/*
	//	C operator()(const D*, const Tuple<unsigned int, DIM>& );
	//	C operator()(const D,const D,const D,... );
		C operator()(const Tuple<D,DIM> );
	*/

	FuncGrid(): start(0xFFFFFFFF) {}
	C operator()(const Tuple<unsigned int, DIM>&coor)const {
		if (start == 0xFFFFFFFF) return (*cl)(data.darray, coor);
		Tuple<unsigned int, DIM> ncoor = coor;
		unsigned int i;
		for(i=0;i<DIM;i++) ncoor[i] += start;
		return (*cl)(data.darray, ncoor);
		}
	void specialize_range(int i_start, int i_rangesize){ start = i_start; rangesize = i_rangesize;}
	void specialize_range_clear(){ start = 0xFFFFFFFF;}

	void getDims(Tuple<unsigned int, DIM> &o_dims) const{
		unsigned int dim = (start == 0xFFFFFFFF) ? data.getSize() : rangesize;
		unsigned int i;for(i=0;i<DIM;i++) o_dims[i] = dim;
		}
	void show(FILE* f = stdout) const{ // ignore shifts
		int ds = data.getSize();
		if (DIM == 1){
			Tuple<unsigned int, 1> ite1;
			for(ite1[0]=0;ite1[0]<ds;ite1[0]++) {
				ExOp::show( (*cl)(data.darray, ite1) , f,1); fprintf(f,"%c", ite1[0] == ds-1 ? '\n' : '\t');
			}
		}else if (DIM == 2){
			Tuple<unsigned int, 2> ite2;
			for(ite2[1]=0;ite2[1]<ds;ite2[1]++)
			for(ite2[0]=0;ite2[0]<ds;ite2[0]++) {
				ExOp::show((*cl)(data.darray, ite2) , f,1); fprintf(f,"%c", ite2[0] == ds-1 ? '\n' : '\t');
			}


		}
		}


		void CholeskyDecom(C* &f_out){ // to do, maybe
			int k,l;
			Tuple<unsigned int, DIM> coor;
			if (DIM == 2){
				k = data.getSize();
				f_out = new C[(k* (k+1)) /2];
				k=0;
				for(coor[1]=0;coor[1]<data.getSize();coor[1]++){
					for(coor[0]=0;coor[0]<coor[1];coor[0]++){
						f_out[k] = (*this)(coor);
				//		for(l = coor[0]-1; l>=0;l--) f_out[k] -= f_out[k+l-coor[0]] * f_out[k+l-coor[0]] * f_out[(l * (l+1)/2)-1];
						f_out[k] /= f_out[(coor[0] * (coor[0]+3))/2];

						k++;
					}
					f_out[k] = (*this)(coor);
					for(l = coor[0]-1; l>=0;l--) f_out[k] -= f_out[k+l-coor[0]] * f_out[k+l-coor[0]] * f_out[(l*(l+3))/2];
					k++;
				}
			}
		}
	};




class TiffFile{
	template<class TYPE> void writeDescriptor(TYPE min, TYPE max,char * &p);
public:
	FILE* f;
	bool inv;
	Vector<unsigned char> curflaglist;
	TiffFile(const char*, bool writeonly = false);
	~TiffFile();
	unsigned int curfp_pos;
	unsigned int endfile_pos;
	bool gotoNext();
	template <class C, unsigned int channels, Tuple_flag Cflag> void fetch(vector< DataGrid<Tuple<C,channels, Cflag>,2> > &out);
	template <class C> bool fetch(DataGrid<C,3u>& f_out, char * imageType = NULL);
	template <class C> bool fetch(DataGrid<C,2u>& f_out, const unsigned int channel = 0);
	template <class C> bool fetch_grayscale(DataGrid<C,2>& f_out);
	
	
	template <class C, class TYPE> bool put(DataGrid<C,3>& f_out, TYPE min, TYPE max, bool updateRange =false);
	template <class C, class TYPE> bool put(DataGrid<C,2>& f_out, TYPE min, TYPE max, bool updateRange =false);
	
    void savePLYasTIFF(const char * path);
	
	void addopt_Description(Vector< char* > & opt, char const * const descript) const;
	void addopt_Xscale(Vector< char* > & opt, double min_x, double max_x) const;
	void addopt_Yscale(Vector< char* > & opt, double min_y, double max_y) const;
	//	template <class TYPE>void addopt_Xscale(Vector< char* > & opt, TYPE min_x = ExCo<TYPE>::zero(), TYPE max_x = ExCo<TYPE>::zero() ) const;
	
	//	template <class C> bool put(const DataGrid<C,3>& f_out, Vector< char* > &options);
	
	
	//	template <class C, int channels, int flag, Tuple_flag Cflag> void put(vector< DataGrid<Tuple<C,channels, Cflag>,1> > &out);
	int flagType(int flagindex);
	
    class WriteScope{
	public:
        unsigned int nbflags;
        //Vector<char> extra_data;
    };
    WriteScope* current;
    void startFrameWrite();
    void endFrameWrite();
	
	template <class C> C getValue(int flagindex);
};

template< > vector<unsigned int> TiffFile::getValue< vector<unsigned int> >(int flagindex);
template< > unsigned int TiffFile::getValue<unsigned int>(int flagindex);
template< > int TiffFile::getValue<int>(int flagindex);

template<class C, int nbc>
class imageIO{
	public:
	static void importBMP(char* path, DataGrid<Tuple<C, nbc>,2> &out);
	static void exportBMP(char* path, DataGrid<Tuple<C, nbc>,2> &out);
};

template<class O, class I, int _out_dim, int _in_dim>
class ContinuousFunction : public Oper2< Tuple<O, _out_dim> , Tuple<I, _in_dim > >{
public:
	virtual void operator()(Tuple<O, _out_dim> &, Tuple<I, _in_dim > &) const =0;
	virtual void derivative(Tuple<O, _out_dim> &, Tuple<I, _in_dim > &, int in_direct) const;
	virtual void derivativeMatrix(TMatrix<O, _in_dim, _out_dim> &, Tuple<I, _in_dim > &) const;
	void derivativeMatrix_default(TMatrix<O, _in_dim, _out_dim> &, Tuple<I, _in_dim > &) const;
	void derivative_from_difference(Tuple<O, _out_dim> &, Tuple<I, _in_dim > &, int in_direct) const;
	void newtonStep(Tuple<I, _in_dim > &);
};


template<class C, unsigned int sizex, unsigned int sizey = sizex>
class TMatrix : public ConstGrid<C,2> {
	const TMatrix<C,sizex,sizey>& leftHouseHolderMultiply(const double * const vec, const int lenght, const double& sqrt_den, bool hint); // bottom of matrix
	const TMatrix<C,sizex,sizey>& rightHouseHolderMultiply(const double * const vec, const int lenght, const double& sqrt_den, bool hint);
public:
	const TMatrix<C,sizex,sizey>& LeftHouseHolderMultiply(const double * const vec, const double& sqrt_den, const int lenght, C* inner, bool hint); // top of matrix
    typedef YESNO<false> IS_COMMUTATIVE; 
	// ExOp section
	static const bool IsPOD = ExCo<C>::IsPOD;
	// ExOp section
	static const bool NeedsAddLink = ExCo<C>::NeedsAddLink;
	C data[sizex*sizey];
	bool fake_inverse;
	TMatrix(){}
	TMatrix(C* idata){memcpy(data,idata,sizeof(C)*sizex*sizey);}
	TMatrix(TMatrix<C, sizex, sizey> const & idata);
	TMatrix(DataGrid<C, 2> const & idata);
	void zero(){unsigned int i = sizex*sizey;for(i--;i != ExCo<unsigned int>::maximum() ; i--) ExOp::toZero(data[i]);}
	TMatrix<C,sizex,sizey>& toZero();
	TMatrix<C,sizex,sizey>& toOne();
	TMatrix<C,sizex,sizey>& toRand();
    TMatrix<C,sizex,sizey>& toRandSymetric();

	Trianglix<C,sizey> operator*(const Trianglix<C,sizex> &input) const;	
		
	C operator()(const Tuple<unsigned int, 2> &coor) const{return data[coor[0] + coor[1] * sizex];}
	void getDims(Tuple<unsigned int, 2> &o_dims) const{o_dims[0] = sizex;o_dims[1] = sizey; }

	void bidiag(); // test

	bool isValid() const{
		int i;
		for(i=0;i<sizex*sizey;i++) if (!ExOp::isValid(data[i])) break;
		return (i == sizex*sizey);
	}

//	template<class O> Tuple<C,size>& operator+=(Tuple<O,size> const & other);
//	template<class O> Tuple<C,size>& operator-=(Tuple<O,size> const & other);
//	template<class O> Tuple<C,size>& operator*=(Tuple<O,size> const & other);
//	template<class O> Tuple<C,size>& operator/=(Tuple<O,size> const & other);

	TMatrix<C,sizex,sizey>& operator=(TMatrix<C, sizex, sizey> const & idata);

	void operator()(Tuple<C, sizey> &, Tuple<C, sizex> &) const;
	void derivative(Tuple<C, sizey> &, Tuple<C, sizex > &, int in_direct) const;
	void derivativeMatrix(TMatrix<C, sizex, sizey> &, Tuple<C, sizex > &) const;

	template<class O> TMatrix<C,sizex,sizey>& operator+=(TMatrix<O,sizex,sizey> const & other){for(unsigned int i=sizex*sizey-1;i!=ExCo<unsigned int>::maximum();i--) data[i] += other.data[i];return(*this);}
	template<class O> TMatrix<C,sizex,sizey>& operator-=(TMatrix<O,sizex,sizey> const & other){for(unsigned int i=sizex*sizey-1;i!=ExCo<unsigned int>::maximum();i--) data[i] -= other.data[i];return(*this);}
	template<class O> TMatrix<C,sizex,sizey>& operator*=(TMatrix<O,sizex,sizex> const & other);
	template<class O> TMatrix<C,sizex,sizey>& operator/=(TMatrix<O,sizex,sizex> const & other);




	template<class O, unsigned int osize> TMatrix<C,osize,sizey> operator*(TMatrix<O,osize,sizex> const & other) const;
	template<class O> Tuple<C,sizey> operator*(Tuple<O,sizex> const & other) const;
	Tuple<C,sizey> operator*(Tuple<C,sizex> const & other) const;


	// template<class O, class A, Tuple_flag FA, Tuple_flag FO> Tuple<O, sizex,FO> operator*(Tuple<A,sizey, FA> const & in);


template<class O> TMatrix<C,sizex,sizey>& operator+=(O const & other){for(unsigned int i=sizex*sizey-1;i!=ExCo<unsigned int>::maximum();i--) data[i] += other;return(*this);}
template<class O> TMatrix<C,sizex,sizey>& operator-=(O const & other){for(unsigned int i=sizex*sizey-1;i!=ExCo<unsigned int>::maximum();i--) data[i] -= other;return(*this);}
template<class O> TMatrix<C,sizex,sizey>& operator*=(O const & other){for(unsigned int i=sizex*sizey-1;i!=ExCo<unsigned int>::maximum();i--) data[i] *= other;return(*this);}
template<class O> TMatrix<C,sizex,sizey>& operator/=(O const & other){for(unsigned int i=sizex*sizey-1;i!=ExCo<unsigned int>::maximum();i--) data[i] /= other;return(*this);}

template<class A> TMatrix<C,sizex,sizey> operator+(A const & other) const {return(TMatrix<C,sizex,sizey>(*this) += other);}
template<class A> TMatrix<C,sizex,sizey> operator-(A const & other) const {return(TMatrix<C,sizex,sizey>(*this) -= other);}
template<class A> TMatrix<C,sizex,sizey> operator*(A const & other) const {return(TMatrix<C,sizex,sizey>(*this) *= other);}
template<class A> TMatrix<C,sizex,sizey> operator/(A const & other) const {return(TMatrix<C,sizex,sizey>(*this) /= other);}
/*
template<class O, class A> TMatrix<O,sizex,sizey> operator+(A const & other) const;
template<class O, class A> TMatrix<O,sizex,sizey> operator-(A const & other) const;
template<class O, class A> TMatrix<O,sizex,sizey> operator*(A const & other) const;
template<class O, class A> TMatrix<O,sizex,sizey> operator/(A const & other) const;*/

	template<class O,Tuple_flag FO> TMatrix<C,sizex,sizey> scale_rows(Tuple<O, sizey, FO> const &) const;
	template<class O,Tuple_flag FO> TMatrix<C,sizex,sizey> scale_cols(Tuple<O, sizex, FO> const &) const;


//template<class O, unsigned int osize> TMatrix<C,sizex,sizey> operator-(O const & other) const{return( (TMatrix(*this)) -= other );}

//	template<class O> TMatrix<C,sizex,sizey> operator*(O const & other) const{return( (TMatrix(*this)) *= other );}
//	template<class O> TMatrix<C,sizex,sizey> operator/(O const & other) const{return( (TMatrix(*this)) /= other );}


	TMatrix<C,sizex,sizey> inverse() const;
	TMatrix<C,sizey,sizex> transpose() const;

//	PolyThing<C,0> makeCharacteristic() const;
	Continuous<C,sizex,C,sizey>* derivative(Tuple<C,sizex> const & where) const; // a TMatrix is it's own derivative!
	void show(FILE* o = stdout, int level=0) const;

    TMatrix<C,sizey,sizex> inverse_2() const;

	C determinant() const;
	void set_to_zero();
};



template<class C, int sizex, int sizey> struct isTypeEquivalent< C*, TMatrix<C,sizex,sizey> > {enum {ans = true }; };
template<class C, int sizex, int sizey> struct isTypeEquivalent< TMatrix<C,sizex,sizey> , C*> {enum {ans = true }; };
template<class C, int sizex, int sizey, int size> struct isTypeEquivalent< Tuple<C,size> ,  TMatrix<C,sizex,sizey> > {enum {ans = true }; };
template<class C, int sizex, int sizey, int size> struct isTypeEquivalent< TMatrix<C,sizex,sizey> , Tuple<C,size> > {enum {ans = true }; };

// LDU decomposed TMatrix, diagonal of L and U is trivial(and not stored)
template<class C, int size>
class Matrianglix{
public:
	C data[size*size];
	Permutation<size> perm; // stores permutation Matrix!

	Matrianglix(){}
	Matrianglix(const TMatrix<C,size,size> & what);
//
	Matrianglix<C,size> operator*(const Matrianglix<C, size> & other);

	TMatrix<C,size,size> inverse();


    string type_tostring()const;
	C determinant();
	operator TMatrix<C,size,size>() const;
//	const TMatrix<C,size,size>& power(double exp);
//  const and io functions
	void show(FILE* o= stdout) const;

};

// heapTree, may contain an unsorted node at position 0!
template <class C>
class HeapTree{
public:
    Vector<C> data;
	bool hasunsorted;
	HeapTree() : hasunsorted(false){data.setSize(1);}
	HeapTree(void* owner): data(owner), hasunsorted(false){data.setSize(1);}
//
	bool isempty() const{return((!hasunsorted) && (data.getSize() == 1));}
    void clear(){data.setSize(1);hasunsorted = false;}
	void insert(const C &what);
	C top() const{return( data[ (  (hasunsorted)&&((data.getSize() == 1)||(data[0] < data[1])) ? 0 : 1) ]);}
	C pop();
};



template <class C>
class HeapTree<C*> : public Vector<C*>{
public:
bool hasunsorted;
HeapTree();
HeapTree(void* owner);

//	Node& operator()(const Key & where);
bool isempty();
void insert(C* what);
C& pop();
void update(unsigned int offset); // Elem was changed! update its position!
C& top()const{return( *(*this)[ (  (hasunsorted)&&(((*this).getSize() == 1)||((*this)[0] < (*this)[1])) ? 0 : 1) ]);}
//	KeyElem<Key, Node> pop();
//	KeyElem<Key, Node> ipop(Node what, Key where); // insert and pop!

};

template <class D>
class defaultHeapBackPointer{
public:
	unsigned int backptr;
	D data;
	static const bool IsPOD = false;
	static const bool NeedsAddLink = false;
	
	defaultHeapBackPointer<D>& operator=(const defaultHeapBackPointer<D> &other){data = other.data; return(*this);}
	bool operator>(const defaultHeapBackPointer<D> &other) const {return( data > other.data);}
	bool operator<(const defaultHeapBackPointer<D> &other) const {return( data < other.data);}
	bool operator>=(const defaultHeapBackPointer<D> &other) const {return( data >= other.data);}
	bool operator<=(const defaultHeapBackPointer<D> &other) const {return( data <= other.data);}
	bool operator==(const defaultHeapBackPointer<D> &other) const {return( data == other.data);}
	bool operator!=(const defaultHeapBackPointer<D> &other) const {return( data != other.data);}
};

template<class C> class HeapTreeBackPointerOffset{
public:
	enum{ans = 0}; // put non-zero if non trivial
	// static int& getBackPtr(C&)
};


template<class D>
class HeapTreeBackPointerOffset< defaultHeapBackPointer<D> >{
public:
	enum{ans = 1};
	static unsigned int& getBackOffset(defaultHeapBackPointer<D>& targ){return(targ.backptr);}
};

class Weight{ // double with special meaning
public:
	double w;
	Weight(){}
	Weight(const double&v):w(v){}
	const double& operator()()const{return(w);}
};

// base struct for recofering mean variance and higher momments of arbritrary weighted objects
template<class C,unsigned int order =1>
class WeightElem{
public:

	typedef WeightElem< typename ExCo<C>::COMPLEX_TYPE, order> COMPLEX_TYPE;
	typedef WeightElem< typename ExCo<C>::REAL_TYPE, order> REAL_TYPE;

	// ExOp Section:
	static const bool IsPOD = ExCo<C>::IsPOD;
	static const bool NeedsAddLink = false; // containers needs to update addresses in link registers

	ExCoMeMdEcLaRe( LFHCONCAT2(WeightElem<C, order>) )

	// ExOp Section end

	Tuple<double, order> w;
	Tuple<C, order> e;
	WeightElem(){}
	WeightElem(const C & ob, double weight =1.0f){

		w[0] =weight;
		for(unsigned int i=1;i<order;i++) w[i] = w[i-1] * weight;
		e[0] = ob * weight;
		for(unsigned int i=1;i<order;i++) e[i] = e[i-1] * ob;

	}
	WeightElem(const WeightElem<C,order>& clonefrom) : w(clonefrom.w), e(clonefrom.e) {}
	template<class O> WeightElem(const WeightElem<O,order>& clonefrom) : w(clonefrom.w), e(clonefrom.e) {}
	WeightElem(const Tuple<C, order>& _e,const Tuple<double, order>& _w) : w(_w), e(_e) {}

	operator C() const{return( (e[0] / w[0]) );}

	WeightElem<C,order>& toZero(){ExOp::toZero(w); ExOp::toZero(e); return(*this);}
	WeightElem<C,order>& toRand(){ExOp::toRand(e[0]); ExOp::toOne(w); for(unsigned int i=1;i<order;i++) e[i] = e[i-1] * e[0];return(*this);}
	WeightElem<C,order>& toOne(){ExOp::toOne(e); ExOp::toOne(w); return(*this);}

	WeightElem<C,order> operator+(const C& s) const {WeightElem<C,order> _out = *this; return(_out += s);}
	WeightElem<C,order> operator-(const C& s) const {WeightElem<C,order> _out = *this; return(_out -= s);}
	WeightElem<C,order> operator*(const C& s) const {WeightElem<C,order> _out = *this; return(_out *= s);}
	WeightElem<C,order> operator/(const C& s) const {WeightElem<C,order> _out = *this; return(_out /= s);}


inline WeightElem< typename ExCo<C>::REAL_TYPE, order> mkrealproj()const{return WeightElem< typename ExCo<C>::REAL_TYPE, order>(ExOp::mkrealproj(e),w);}
inline WeightElem< typename ExCo<C>::REAL_TYPE, order> mkimmaproj()const{return WeightElem< typename ExCo<C>::REAL_TYPE, order>(ExOp::mkimmaproj(e),w);}
inline WeightElem< typename ExCo<C>::REAL_TYPE, order> mkjmmaproj()const{return WeightElem< typename ExCo<C>::REAL_TYPE, order>(ExOp::mkjmmaproj(e),w);}
inline WeightElem< typename ExCo<C>::REAL_TYPE, order> mkkmmaproj()const{return WeightElem< typename ExCo<C>::REAL_TYPE, order>(ExOp::mkkmmaproj(e),w);}
WeightElem<C,order>& operator=(const WeightElem<C,order>& o){w = o.w;e = o.e;return(*this);}
	WeightElem<C,order>& operator+=(const C& s){
		Tuple<C,order> buf = ExOpa<order>::intPows(s);

		if (order >2){
			e[1] += (buf[0] * e[1] + buf[1] * e[0]) * 3.0f + buf[2] * w[0];
		}
		if (order >1){
			e[1] += buf[0] * e[0] * 2.0f + buf[1] * w[0];
		}
		e[0] += s * w[0];
		return(*this);
	}

	WeightElem<C,order>& operator-=(const C& s){
		Tuple<C,order> buf = ExOpa<order>::intPows(s);
		if (order >2){
			e[1] += (buf[0] * e[1] - buf[1] * e[0]) * -3.0f - buf[2] * w[0];
		}
		if (order >1){
			e[1] += (s * e[0]) * -2.0f + ExCo<C>::intPow(s,2) * w[0];
		}
		e[0] -= s * w[0];
		return(*this);
	}

	WeightElem<C,order>& operator*=(const C& s){
		Tuple<C,order> buf = ExOpa<order>::intPows(s);
		if (order >2) e[2] *= buf[2];
		if (order >1) e[1] *= buf[1];
		e[0] *= s;
		return(*this);
	}


	WeightElem<C,order>& operator+=(const WeightElem<C,order>& o ){w += o.w;  e += o.e;return(*this);}
	WeightElem<C,order>& operator-=(const WeightElem<C,order>& o ){w += o.w;  e -= o.e;return(*this);}
	WeightElem<C,order>& operator*=(const Weight& _w){
		double o = _w();
		double f = o;
		e *= o;
		w[0] *=o;
		for(unsigned int i=1;i<order;i++) {f *=o; w[i] *= f;}
		return(*this);
		}
	WeightElem<C,order>& operator/=(const Weight& _w ){
		double o = 1.0f / _w();
		double f = o;
		e *= o;
		w[0] *=o;
		for(unsigned int i=1;i<order;i++) {f *=o; w[i] *= f;}
		return(*(WeightElem<C,order>*)this);
		}
	WeightElem<C,order> operator+(const WeightElem<C,order>& o )const {return(WeightElem<C,order>(e +o.e,w +o.w)); }
	WeightElem<C,order> operator-(const WeightElem<C,order>& o )const;
	template<class O> WeightElem<C,order> operator+(const O&o)const{WeightElem<C,order> f_out = *this; return (f_out += o);}
	template<class O> WeightElem<C,order> operator-(const O&o)const{WeightElem<C,order> f_out = *this; return (f_out -= o);}
	template<class O> WeightElem<C,order> operator*(const O&o)const{WeightElem<C,order> f_out = *this; return (f_out *= o);}
	template<class O> WeightElem<C,order> operator/(const O&o)const{WeightElem<C,order> f_out = *this; return (f_out /= o);}

void setWeight(double n_w){e *= (n_w / w[0]); if (order >2) w[1] *= (n_w / w[0]); w[0] = n_w; }

	C getMean() const {	return( (e[0] / w[0]) );}
	C getVar() const {return( (this->e[1]*this->w[0] - this->e[0] * this->e[0]) * (1.0f / (this->w[0]*this->w[0] - this->w[1])) );}
	C getVar_biaised() const {return( (this->e[1]*this->w[0] - this->e[0] * this->e[0]) * (1.0f / (this->w[0]*this->w[0])) );}

	// adds a weight of this->w[1] / this->w[0]*this->w[0] for a given variance
//	C getVar_underPrior(C prior) const {return(  (this->w[1] >= this->w[0]*this->w[0]) ? prior : (this->e[1]*this->w[0] - this->e[0] * this->e[0] + prior * this->w[1] ) * (1.0f / (this->w[0]*this->w[0]))    );}
	// adds a weight of (this->w[1] / this->w[0]*this->w[0])^2 for a given variance
    C getVar_underPrior(C prior) const {double ratio = this->w[1] / (this->w[0]*this->w[0]); return(  (ratio > 0.5f)?  prior :  (this->e[1]*this->w[0] - this->e[0] * this->e[0]) * (1.0f / (this->w[0]*this->w[0] - this->w[1]))) ;}

	bool operator<(const WeightElem<C,order> &other){
		return( w[0] < other.w[0] );
	}
	C getSkew() const {
//		LinkAssert< (order>2) > ass;
		double w2 = this->w[0]*this->w[0];
		/*
		 printf("%f\t%f\t%f\n",this->e[3],this->w[0],this->e[3]/this->w[0]);
		 printf("%f\t%f\t%f\n",this->e[4],this->w[1],this->e[4]/this->w[1]);
		 printf("%f\t%f\t%f\n",this->e[5],this->w[2],this->e[5]/this->w[2]);

		 printf("%f\t%f\t%f\n",- (this->e[0]*this->e[2] - this->e[5])*3.0f,this->w[0]*this->w[1] - this->w[2],- (this->e[0]*this->e[2] - this->e[5])*3.0f/(this->w[0]*this->w[1] - this->w[2]));
		 printf("%f\t%f\t%f\n",- (this->e[0]*this->e[1] - this->e[4])*3.0f,w2 - this->w[1],- (this->e[0]*this->e[1] - this->e[4])*3.0f/(w2 - this->w[1]));


		 printf("%f\t%f\t%f\n",(this->e[0]*(this->e[0]*this->e[0] - this->e[2]*3.0f) + this->e[5] * 2.0f )*2.0f , (this->w[0]*(w2 -3.0f*this->w[1]) + 2.0f* this->w[2]),(this->e[0]*(this->e[0]*this->e[0] - this->e[2]*3.0f) + this->e[5] * 2.0f )*2.0f / (this->w[0]*(w2 -3.0f*this->w[1]) + 2.0f* this->w[2]));

		 printf("%f\t%f\t%f\n",this->w[0]*this->w[0]*this->e[3] -3*this->w[0]*this->e[0]*this->e[1] +2*this->e[0]*this->e[0]*this->e[0], (this->w[0]*(w2 -3.0f*this->w[1]) + 2.0f* this->w[2]), (this->w[0]*this->w[0]*this->e[3] -3*this->w[0]*this->e[0]*this->e[1] +2*this->e[0]*this->e[0]*this->e[0])/ ((this->w[0]*(w2 -3.0f*this->w[1]) + 2.0f* this->w[2])));
		 */
		/*
		 return( this->e[3] * (1 / this->w[0])
		 - (this->e[0]*this->e[1] - this->e[4])*(3.0f / ( w2 - this->w[1]))
		 + (this->e[0]*(this->e[0]*this->e[0] - this->e[2]*2.0f) + this->e[5] * 1.0f )*(2.0f/(this->w[0]*(w2 -2.0f*this->w[1]) + 1.0f* this->w[2]) )
		 );

		 */
		//	 return((w2*this->e[2] -3*this->w[0]*this->e[0]*this->e[1] +2*this->e[0]*this->e[0]*this->e[0])
		//	 / ((this->w[0]*(w2 -3.0f*this->w[1]) + 2.0f* this->w[2])));
		//		return(((this->e[2]* w2) +this->e[0]*this->e[1]*(-3.0f*this->w[0]) +this->e[0]*this->e[0]*this->e[0]*2.0f) *
		//			   (1.0f / ((this->w[0]*(w2 -3.0f*this->w[1]) + 2.0f* this->w[2]))));

				return((this->e[2]* (1.0f/this->w[0])) +this->e[0]*this->e[1]*(-3.0f/w2) +this->e[0]*this->e[0]*this->e[0]*(2.0f / (w2 * this->w[0]))) ;
		//		return( (this->e[2] +this->e[0] * (
		//				((this->e[1]*this->w[0] - this->e[0] * this->e[0]) * (-3.0f / (this->w[0]*this->w[0] - this->w[1]))) + this->e[0]*this->e[0]*(2.0f /(this->w[0]*this->w[0]))
		//			   ) )  * (1.0f/this->w[0])) ;
	}
	C getVar_scaleinv() const {	double v = getVar();if (v == 0.0f) return(0.0f);return sqrt(v) / getMean();}
	C getSkew_scaleinv() const {double v = getVar();if (v == 0.0f) return(0.0f);return getSkew() * pow(v, -1.5f);}
	C getKurt_scaleinv() const {double v = getVar();if (v == 0.0f) return(0.0f); return getKurt() * pow(v, -2.0f);}

	C getKurt() const {
//		LinkAssert< (order>3) > ass;
		//$K = \frac{V_4}{C_1}
		//- \frac{4V_1V_3 + 3V_2^2 - 7W_4}{C_1^2 - C_2}
		// + 12*\frac{V_1^2V_2 - 2V_1W_3 - V_2W_2 + 2W_4}{C_1^3- 3C_1C_2+ 2C_3}
		// - 6*\frac{V_1^4 - 6V_1^2W_2 +8V_1W_3 + 3W_2^2 -6W_4}{C_1^4 - 6C_1^2C_2 + 8C_1C_3 +3C_2^2 -6C_4}$

		/*
		 return( this->e[6] * (1 / this->w[0])
		 - (this->e[0]*this->e[3]*4.0f + this->e[1]*this->e[1]*3.0f + this->e[7]* -7.0f)*(1.0f / ( w2 - this->w[1]))
		 + ( this->e[1] * (e2 - this->e[2]) + this->e[0]* this->e[4] * -2.0f + this->e[8] * 2.0f )
		 *(12.0f/(this->w[0]*(w2 -3.0f*this->w[1]) + 2.0f* this->w[2]) )
		 - ( e2*(e2 + this->e[2]*-6.0f) + this->e[0] * this->e[5] *8.0f + this->e[2] *this->e[2] *3.0f + this->e[9] * -6.0f )
		 *(6.0f/(w2*w2 + 8.0f*this->w[0]*this->w[2] - 3.0f*this->w[1]*(2.0f*w2 - this->w[1]) - 6*this->w[3]))
		 );

		 */
		double w2 = this->w[0]*this->w[0];

		C e2 = this->e[0]*this->e[0] * (1.0f / w2) ;
		/*
		double den = (1.0f /(w2*w2 + 8.0f*this->w[0]*this->w[2] - 3.0f*this->w[1]*(2.0f*w2 - this->w[1]) - 6*this->w[3]) );
		return( ((this->e[3]*(this->w[0]) + this->e[2]*this->e[0]*(-4.0f) + this->e[1]*this->e[1]*(-3.0f))*w2 +this->e[1]*e2*(12.0f*this->w[0]) +e2*e2*(-6.0f) ) *den);*/
		
	//	return( this->e[3]*(1.0f / this->w[0]) + this->e[2]*this->e[0]*(-4.0f / w2) + this->e[1]*this->e[1]*(-3.0f / w2) +this->e[1]*e2*(12.0f /this->w[0]) +e2*e2*(-6.0f) );
		return( this->e[3]*(1.0f / this->w[0]) + this->e[2]*this->e[0]*(-4.0f / w2) +e2* (this->e[1]*(6.0f /this->w[0]) +e2*(-3.0f)) );
		//		 return( (w2*this->w[0]*this->e[3] -4*w2*this->e[2]*this->e[0] +6*this->w[0]*this->e[1]*e2 -3*e2*e2 )
		//		 /(w2*w2 + 8.0f*this->w[0]*this->w[2] - 3.0f*this->w[1]*(2.0f*w2 - this->w[1]) - 6*this->w[3]) );
	}

	C getSecondMomment(){ // E[(X - mu)^2] , biased variance
		return (this->e[1]  - this->e[0] * this->e[0] * (1.0f / this->w[0])) * (1.0f / this->w[0]);
	}

	C getFouthMomment(){ // E[(X - mu)^4]
		C e2 = this->e[0] * (1.0f / this->w[0]);
		return (this->e[3] - e2 * (this->e[2] * 4.0f - e2 * (this->e[1] * 6.0f - this->e[0] * e2 * 3.0f))) * (1.0f / this->w[0]);
	}


	C operator[](int i){
		if (i>order) return(ExCo<C>::zero());
		switch(i){
			case 1: return(this->getMean());
			case 2: return(this->getVar());
			case 3: return(this->getSkew());
			case 4: return(this->getKurt());
			default: return(this->getMean());
		}
	}

	class OpMean : LFHDECL_OPER2(C,LFHCONCAT2(WeightElem<C, order>));
	class OpMom : LFHDECL_OPER2(LFHCONCAT2(Tuple<C, order>),LFHCONCAT2(WeightElem<C, order>));
	class OpVar : LFHDECL_OPER2(C,LFHCONCAT2(WeightElem<C, 2>));
	class OpStd : LFHDECL_OPER2(C,LFHCONCAT2(WeightElem<C, 2>));


	double getHDist(const WeightElem<C,order> &other){return sqrt(1.0f - exp(-getBattDist(other)));}
	double getBattDist(const WeightElem<C,order> &other) const {
	LinkAssert< (order>1) > ass;
		C diff = (other.e[0] * this->w[0] - this->e[0] * other.w[0]);
		double d = ExOp::pnorm(diff) / (w[0] * other.w[0]);
		double v1 = ExOp::norm(this->getVar());
		double v2 = ExOp::norm(other.getVar());

		d /= 0.5f * this->w[0] * other.w[0] * (v1+v2);
		return((d - log (v1) -log(v2)) * 0.5f + log(v1+v2) - log(2));
	}


	void show(FILE* out = stdout, int level=0){
		switch(level){
			case 0:
				fprintf(out,"Weight: %f\n", w[0]);
				fprintf(out,"Mean: ", w[0]); ExOp::show(getMean(),out, 1);
				if (order>1) {fprintf(out,"\nVar: ", w[0]); ExOp::show(getVar(),out, 1);}
				fprintf(out,"\n");
			break;
			case 1:

				break;
		}
	}
	double getWeight() const{return w[0];}
};

#ifdef GNU_SCIENTIFIC_LIBRARY
class LogPvalSum{
public:
	unsigned int nb_independent;
	double value;
	
	LogPvalSum& toZero(){nb_independent =0; value=0.0f; return(*this);}
	void operator+=(const double val){if (ExOp::isValid(val)) {nb_independent++; value+= val;} }
	double operator()() const{ return (nb_independent > 1) ?  LogPvalue_SumLogPvalues_Ptail(value, nb_independent ) : value; }
	void compress(){if (nb_independent > 1) {value = LogPvalue_SumLogPvalues_Ptail(value, nb_independent ); nb_independent =1;} }
};
#endif

template<class C, unsigned int flag>
class GaussElem{
public:
    GaussElem();
    // Tuple<C,SIZE> getMean() const; inherited 
};

template<class C, unsigned int SIZE>
//class GaussElem< Tuple<C, SIZE>, (ExCo<C>::IS_COMMUTATIVE::ans == 0) ? 0u : 1u >{ // COMMUTATIVE!
class GaussElem< Tuple<C, SIZE>, 0u >{ // COMMUTATIVE!
public:
    double w,w2;
    Tuple<C, SIZE> mean;
    Trianglix<C, SIZE> cov;
    mutable C determinant; // set when needed
	static const int TARG2 =0;
	
	static const bool IsPOD = Tuple<C, SIZE>::IsPOD && Trianglix<C, SIZE>::IsPOD;
	
	
    GaussElem(){ExOp::toZero(determinant);}
	
    GaussElem(const Tuple<WeightElem<C,2>, SIZE> &o, double _w=1.0f): w(_w), w2(0.0f){unsigned int i,j,k; for(i=0,k=0;i< SIZE;i++) {mean[i] = o[i].getMean(); for(j=0;j< i;j++) {cov.data[k++] =mean[i] * mean[j];} cov.data[k++] = o[i].getVar() + mean[i] * mean[i];}
		ExOp::toZero(determinant);
		}
    GaussElem(const Tuple<C, SIZE> &o, double _w=1.0f): w(_w), w2(_w*_w), mean(o * w), cov(o * sqrt(w)){ExOp::toZero(determinant); }
    GaussElem(const Tuple<C, SIZE> &_mean, const Trianglix<C, SIZE>& _cov, double _w, double _w2): w(_w), w2(_w2), mean(_mean), cov(_cov){ExOp::toZero(determinant);}

    // ExOp
    double getWeight()const{return w;}
	double getN()const{return (w == 0) ? 0.0f : w*w/w2;}
	inline unsigned int getSize() const{return SIZE;}
	void setSize(unsigned int size){mean.setSize(size); cov.setSize(size); }
 //   void show(FILE* f = stdout, int level =0)const;
	
	GaussElem< Tuple<C, SIZE>, 0 >& operator+=(const Tuple<C, SIZE>& shift){unsigned int i,j,k; for(j=0,k=0;j<SIZE;j++) for(i=0;i<=j;i++,k++) cov.data[k] += (((ExOp::mktrju(shift[i]) * shift[j]) * w) + ExOp::mktrju(shift[i]) * mean[j]+ ExOp::mktrju(mean[i]) * shift[j]); mean += shift * w; return(*this);}
	GaussElem< Tuple<C, SIZE>, 0 >& operator-=(const Tuple<C, SIZE>& shift){unsigned int i,j,k; for(j=0,k=0;j<SIZE;j++) for(i=0;i<=j;i++,k++) cov.data[k] += (((ExOp::mktrju(shift[i]) * shift[j]) * w) - ExOp::mktrju(shift[i]) * mean[j]- ExOp::mktrju(mean[i]) * shift[j]); mean -= shift * w; return(*this);}


		
	Tuple<C,SIZE> getMean() const {	return( mean / w );}
    Tuple<C,SIZE> getVar() const {Tuple<C,SIZE> fout;  double fact = (1.0f / (w*w - w2)); unsigned int i,j; for(i=0,j=0;j<SIZE;i+=j+1) {fout[j] = ((cov.data[i]*w) - (mean[j] * mean[j])); j++;} fout *= fact; return(fout);}
	Tuple<C,SIZE> getVar_biaised() const {Tuple<C,SIZE> fout;  double fact = (1.0f / (w*w)); unsigned int i,j; for(i=0,j=0;j<SIZE;i+=j+1) {fout[j] = ((cov.data[i]*w) - (mean[j] * mean[j])); j++;} fout *= fact; return(fout);}
	
    GaussElem< Tuple<C, SIZE>, 0 > & toZero(){w=0.0f;w2=0.0f;ExOp::toZero(mean); ExOp::toZero(cov);ExOp::toZero(determinant);return(*this);}
    GaussElem< Tuple<C, SIZE>, 0 > & toRand(){w=1.0f;w2=1.0f;ExOp::toRand(mean); cov = Trianglix<C, SIZE>(mean);ExOp::toZero(determinant); return(*this);}
    GaussElem< Tuple<C, SIZE>, 0 > & operator=(const GaussElem< Tuple<C, SIZE>,TARG2 > &other){mean = other.mean;cov = other.cov;w = other.w;w2 = other.w2; determinant = other.determinant;return(*this);}

	GaussElem< Tuple<C, SIZE>, TARG2 >& operator+=(const GaussElem< Tuple<C, SIZE>,TARG2 > &other){w+=other.w; w2+=other.w2; mean += other.mean; cov += other.cov; ExOp::toZero(determinant); return(*this);}
    GaussElem< Tuple<C, SIZE>, TARG2 > operator+(const GaussElem< Tuple<C, SIZE>, TARG2 > &other)const {return (GaussElem< Tuple<C, SIZE>, TARG2 >(mean +other.mean,  cov + other.cov, w+other.w, w2+other.w2));}

	// free mean addition
	GaussElem< Tuple<C, SIZE>, 0 >& addGaussElem_free_mean(const GaussElem< Tuple<C, SIZE>, TARG2 > &other);
	
	// REVERT SCOPE!
	GaussElem< Tuple<C, SIZE>, TARG2 >& operator-=(const GaussElem< Tuple<C, SIZE>,TARG2 > &other){w-=other.w; w2-=other.w2; mean -= other.mean; cov -= other.cov; ExOp::toZero(determinant); return(*this);}

	GaussElem< Tuple<C, SIZE>, 0u >& operator*=(const C& value){ mean *= value; cov *= ExOp::mkintpow(value,2); return (*this);}
	
	GaussElem< Tuple<C, SIZE>, 0u >& operator*=(const Weight& inc_w){w *= inc_w();w2 *= inc_w()*inc_w(); mean *= inc_w(); cov *= inc_w(); return (*this);}
	GaussElem< Tuple<C, SIZE>, 0u > operator*(const Weight& inc_w) const {return GaussElem< Tuple<C, SIZE>, 0u >(mean * inc_w(), cov * inc_w(), w * inc_w(),w2 * inc_w()*inc_w());}
	GaussElem< Tuple<C, SIZE>, 0u >& operator/=(const Weight& inc_w){w /= inc_w();w2 /= inc_w()*inc_w(); mean /= inc_w(); cov /= inc_w(); return (*this);}
	GaussElem< Tuple<C, SIZE>, 0u > operator/(const Weight& inc_w) const {return GaussElem< Tuple<C, SIZE>, 0u >(mean / inc_w(), cov / inc_w(), w / inc_w(),w2 / (inc_w()*inc_w()));}
	
	GaussElem< Tuple<C, SIZE>, 0u >& toHouseHolderMultiply(const Tuple<C, SIZE> &vec, const C &denum2);
	GaussElem< Tuple<C, SIZE>, 0u > mkHouseHolderMultiply(const Tuple<C, SIZE> &vec, const C &denum2) const{GaussElem< Tuple<C, SIZE>, 0u > fout(*this); fout.toHouseHolderMultiply(vec,denum2); return(fout);}
	
	
LFH_GOLD double likelihoodratio_dist(const GaussElem< Tuple<C, SIZE>, TARG2 > &other, bool has_weight = true, bool has_covariance = true) const;
    double bhattacharryya_dist(const GaussElem< Tuple<C, SIZE>, TARG2 > &other, bool has_weight = true, bool has_covariance = true) const; // weighted does not work
	
    
LFH_GOLD double LLikelihood(bool has_covariance =true) const{return (w == 0.0f)? 0.0f : Entropy(has_covariance) * -w;}
LFH_GOLD double Entropy(bool has_covariance =true ) const;
	
    double LLikelihood(const Tuple<C, SIZE> &what, bool has_covariance = true)const;
    void setWeight(double _w){ mean *= (_w / w); cov *= (_w / w); w2 *= (_w / w); w = _w;ExOp::toZero(determinant);}

    Trianglix<C, SIZE> getCovariance()const;
	Trianglix<C, SIZE> getCovariance_biased()const;

	template<class D, unsigned int O_SIZE> GaussElem< Tuple<C, O_SIZE>, 0u > operator*(const TMatrix<D,SIZE , O_SIZE> &) const;
	
	template<class D> GaussElem< Tuple<C, SIZE>, 0u > operator*(const Matrix<D> &);
	
    void setCovariance(const Trianglix<C, SIZE>& covar);
    void setCovariance(const GaussElem< Tuple<C, SIZE>, 0 >& other){setCovariance(other.getCovariance());}

//	C getVar() const {return( (this->e[1]*this->w[0] - this->e[0] * this->e[0]) * (1.0f / (this->w[0]*this->w[0] - this->w[1])) );}
    void show(FILE *f = stdout, int level=0)const{fprintf(f,"mean:\t"); ExOp::show(mean / w,f,1); fprintf(f,"\nCovar\n"); ExOp::show(getCovariance(),f,0);  }
	
	template<unsigned int S_SIZE> GaussElem< Tuple<C, S_SIZE>, 0u > makeSubGauss(const Tuple<unsigned int, S_SIZE>& selected_dims)const;
	template<unsigned int S_SIZE> GaussElem< Tuple<C, S_SIZE>, 0u > makeResidualGauss(const Tuple<unsigned int, S_SIZE>& selected_dims)const;
	
#ifdef GNU_SCIENTIFIC_LIBRARY 
	double PvalueOf(const GaussElem< Tuple<C, SIZE>, 0u > &x)const;
	LFH_GOLD double Pvalue_Hotelling(const GaussElem< Tuple<C, SIZE>, 0u > &x, bool LogPvalue = false) const;
	template<unsigned int S_SIZE> double PvalueOfResidual(const GaussElem< Tuple<C, SIZE>, 0u > &x , const Tuple<unsigned int, S_SIZE>& selected_dims)const; // TODO!
#endif
	void setMeanVaronDimention(double mean, double var, unsigned int dim);
	
	
	
	
	void load(FILE *f, unsigned int size=0);
	void save(FILE *f) const;
	
};

template<class C>
class GaussElem< Tuple<C, 0u>, 0u >{ // COMMUTATIVE!
public:
    double w,w2;
    Tuple<C, 0u> mean;
    Trianglix<C, 0u> cov;
    mutable C determinant; // set when needed
	static const int TARG2 =0;
    GaussElem(){ExOp::toZero(determinant);}
	
    GaussElem(const Tuple<WeightElem<C,2>, 0u> &o, double _w=1.0f): w(_w), w2(0.0f){unsigned int i,j,k; for(i=0,k=0;i< mean.getSize();i++) {mean[i] = o[i].getMean(); for(j=0;j< i;j++) {cov.data[k++] =mean[i] * mean[j];} cov.data[k++] = o[i].getVar() + mean[i] * mean[i];}
		ExOp::toZero(determinant);
	}
    GaussElem(const Tuple<C, 0u> &o, double _w=1.0f): w(_w), w2(_w*_w), mean(o * w), cov(o * sqrt(w)){ExOp::toZero(determinant); }
    GaussElem(const Tuple<C, 0u> &_mean, const Trianglix<C, 0u>& _cov, double _w, double _w2): w(_w), w2(_w2), mean(_mean), cov(_cov){ExOp::toZero(determinant);}
	
    // ExOp
    double getWeight()const{return w;}
	double getN()const{return (w == 0) ? 0.0f : w*w/w2;}
	
	void setSize(unsigned int size){mean.setSize(size); cov.setSize(size); }
	inline unsigned int getSize() const{return mean.getSize();}
	//   void show(FILE* f = stdout, int level =0)const;
	Tuple<C,0u> getMean() const {	return( mean / w );}
    Tuple<C,0u> getVar() const {Tuple<C,0u> fout;  fout.setSize(mean.getSize()); double fact = (1.0f / (w*w - w2)); unsigned int i,j; for(i=0,j=0;j<mean.getSize();i+=j+1) {fout[j] = ((cov.data[i]*w) - (mean[j] * mean[j])); j++;} fout *= fact; return(fout);}
	Tuple<C,0u> getVar_biaised() const {Tuple<C,0u> fout; fout.setSize(mean.getSize()); double fact = (1.0f / (w*w)); unsigned int i,j; for(i=0,j=0;j<mean.getSize();i+=j+1) {fout[j] = ((cov.data[i]*w) - (mean[j] * mean[j])); j++;} fout *= fact; return(fout);}
	
    GaussElem< Tuple<C, 0u>, 0 > & toZero(){w=0.0f;w2=0.0f;ExOp::toZero(mean); ExOp::toZero(cov);ExOp::toZero(determinant);return(*this);}
    GaussElem< Tuple<C, 0u>, 0 > & toRand(){w=1.0f;w2=1.0f;ExOp::toRand(mean); cov = Trianglix<C, 0u>(mean);ExOp::toZero(determinant); return(*this);}
    GaussElem< Tuple<C, 0u>, 0 > & operator=(const GaussElem< Tuple<C, 0u>,TARG2 > &other){mean = other.mean;cov = other.cov;w = other.w;w2 = other.w2; determinant = other.determinant;return(*this);}
	
	GaussElem< Tuple<C, 0u>, TARG2 >& operator+=(const GaussElem< Tuple<C, 0u>,TARG2 > &other){w+=other.w; w2+=other.w2; mean += other.mean; cov += other.cov; ExOp::toZero(determinant); return(*this);}
    GaussElem< Tuple<C, 0u>, TARG2 > operator+(const GaussElem< Tuple<C, 0u>, TARG2 > &other)const {return (GaussElem< Tuple<C, 0u>, TARG2 >(mean +other.mean,  cov + other.cov, w+other.w, w2+other.w2));}
	// free mean addition
	GaussElem< Tuple<C, 0u>, 0 >& addGaussElem_free_mean(const GaussElem< Tuple<C, 0u>, TARG2 > &other);
	
	// REVERT SCOPE!
	GaussElem< Tuple<C, 0u>, TARG2 >& operator-=(const GaussElem< Tuple<C, 0u>,TARG2 > &other){w-=other.w; w2-=other.w2; mean -= other.mean; cov -= other.cov; ExOp::toZero(determinant); return(*this);}
	
	GaussElem< Tuple<C, 0u>, 0u >& operator*=(const C& value){ mean *= value; cov *= ExOp::mkintpow(value,2); return (*this);}
	
	GaussElem< Tuple<C, 0u>, 0u >& operator*=(const Weight& inc_w){w *= inc_w();w2 *= inc_w()*inc_w(); mean *= inc_w(); cov *= inc_w(); return (*this);}
	GaussElem< Tuple<C, 0u>, 0u > operator*(const Weight& inc_w) const {return GaussElem< Tuple<C, 0u>, 0u >(mean * inc_w(), cov * inc_w(), w * inc_w(),w2 * inc_w()*inc_w());}
	GaussElem< Tuple<C, 0u>, 0u >& operator/=(const Weight& inc_w){w /= inc_w();w2 /= inc_w()*inc_w(); mean /= inc_w(); cov /= inc_w(); return (*this);}
	GaussElem< Tuple<C, 0u>, 0u > operator/(const Weight& inc_w) const {return GaussElem< Tuple<C, 0u>, 0u >(mean / inc_w(), cov / inc_w(), w / inc_w(),w2 /(inc_w()*inc_w()));}
	
	GaussElem< Tuple<C, 0u>, 0u >& toHouseHolderMultiply(const Tuple<C, 0u> &vec, const C &denum2);
	GaussElem< Tuple<C, 0u>, 0u > mkHouseHolderMultiply(const Tuple<C, 0u> &vec, const C &denum2) const{GaussElem< Tuple<C, 0u>, 0u > fout(*this); fout.toHouseHolderMultiply(vec,denum2); return(fout);}
	
	
	LFH_GOLD double likelihoodratio_dist(const GaussElem< Tuple<C, 0u>, TARG2 > &other, bool has_weight = true, bool has_covariance = true) const;
    double bhattacharryya_dist(const GaussElem< Tuple<C, 0u>, TARG2 > &other, bool has_weight = true, bool has_covariance = true) const;
	
    
	LFH_GOLD double LLikelihood(bool has_covariance =true ) const{return (w == 0.0f)? 0.0f : Entropy(has_covariance) * -w;}
	LFH_GOLD double Entropy(bool has_covariance =true ) const;
	
    void setWeight(double _w){ mean *= (_w / w); cov *= (_w / w); w2 *= (_w / w); w = _w;ExOp::toZero(determinant);}
	
    Trianglix<C, 0u> getCovariance()const;
	Trianglix<C, 0u> getCovariance_biased()const;
	
	template<class D> GaussElem< Tuple<C, 0u>, 0u > operator*(const Matrix<D> &);
	
    void setCovariance(const Trianglix<C, 0u>& covar);
    void setCovariance(const GaussElem< Tuple<C, 0u>, 0 >& other){setCovariance(other.getCovariance());}
	
	//	C getVar() const {return( (this->e[1]*this->w[0] - this->e[0] * this->e[0]) * (1.0f / (this->w[0]*this->w[0] - this->w[1])) );}
    void show(FILE *f = stdout, int level=0)const{fprintf(f,"mean:\t"); ExOp::show(mean / w,f,1); fprintf(f,"\nCovar\n"); ExOp::show(getCovariance(),f,0);  }
	
	template<unsigned int S_SIZE> GaussElem< Tuple<C, S_SIZE>, 0u > makeSubGauss(const Tuple<unsigned int, S_SIZE>& selected_dims)const;
	template<unsigned int S_SIZE> GaussElem< Tuple<C, S_SIZE>, 0u > makeResidualGauss(const Tuple<unsigned int, S_SIZE>& selected_dims)const;
	
	double PvalueOf(const GaussElem< Tuple<C, 0u>, 0u > &x)const;
	template<unsigned int S_SIZE> double PvalueOfResidual(const GaussElem< Tuple<C, 0u>, 0u > &x , const Tuple<unsigned int, S_SIZE>& selected_dims)const; // TODO!



	double getLogFoldinVariance_inMeanMajorAxis(const GaussElem< Tuple<C, 0u>, 0u > &x) const;
	double getLogFoldinVariance_Determinant(const GaussElem< Tuple<C, 0u>, 0u > &x) const;
	double getLogFoldinVariance_Diagonal(const GaussElem< Tuple<C, 0u>, 0u > &x) const;
	double getDifferenceinMean_inVarianceMajorAxis(const GaussElem< Tuple<C, 0u>, 0u > &x) const;
		#ifdef GNU_SCIENTIFIC_LIBRARY
	
	LFH_GOLD double Pvalue_Hotelling(const GaussElem< Tuple<C, 0u>, 0u > &x, bool LogPvalue = false) const;
	double Pvalue_LRTest_UnequalMean(const GaussElem< Tuple<C, 0u>, 0u > &x, const Tuple<bool> *channel_filter = NULL,  bool LogPvalue = false) const; // Ho = mu_a = mu_b V_a = V_b   H1 = mu_a != mu_b V_a = V_b
	double Pvalue_LRTest_UnequalMeanandVariance(const GaussElem< Tuple<C, 0u>, 0u > &x, const Tuple<bool> *channel_filter = NULL, bool LogPvalue = false) const; // Ho = mu_a = mu_b V_a = V_b   H1 = mu_a != mu_b V_a != V_b
	double Pvalue_LRTest_UnequalVariance(const GaussElem< Tuple<C, 0u>, 0u > &x, const Tuple<bool> *channel_filter = NULL, bool LogPvalue = false) const; // Ho = mu_a != mu_b V_a = V_b   H1 = mu_a != mu_b V_a != V_b
#endif
//	double Pvalue_LRTest_UnequalVariance_Masked(const GaussElem< Tuple<C, 0u>, 0u > &x, const Trianglix<bool> &mask, bool LogPvalue = false) const; // Ho = mu_a != mu_b V_a = V_b   H1 = mu_a != mu_b V_a != V_b

	
	void load(FILE *f, unsigned int size=0);
	void save(FILE *f) const;
	
	GaussElem<Tuple<C, 0u> > mkSubGaussElem(const Tuple<bool> filter)const;
	
};

template<class C, unsigned int SIZE>
//class GaussElem< Tuple<C, SIZE>, (ExCo<C>::IS_COMMUTATIVE::ans == 0) ? 1u : 0u >{ // NON-COMMUTATIVE!
class GaussElem< Tuple<C, SIZE>, 1u >{ // NON-COMMUTATIVE!
public:
    double w,w2;
    Tuple<C, SIZE> mean;
	static const int TARG2 =1; 

   // void show(FILE* f = stdout, int level =0)const;
	Tuple<C,SIZE> getMean() const {	return( (C)(mean / w) );}
//	C getVar() const {return( (this->e[1]*this->w[0] - this->e[0] * this->e[0]) * (1.0f / (this->w[0]*this->w[0] - this->w[1])) );}
    void show(FILE *f = stdout, int level=0)const{fprintf(f,"dont commute!\n");}
};
// Assumes C are commutative!

template<class Key>
class defaultHashFnc{
    public:
       static inline unsigned int makeSeed(const Key &inv) {  return (unsigned int) inv;  } 
};

template<class Key, class Key2>
class defaultHashFnc< pair< Key, Key2>  >{
    public:
       static inline unsigned int makeSeed(const pair< Key, Key2> &in){ return defaultHashFnc<Key>::makeSeed(in.first) ^ defaultHashFnc<Key2>::makeSeed(in.second) ;  }
};

template<class Key, class Key2>
class defaultHashFnc< KeyElem< Key, Key2>  >{
    public:
       static inline unsigned int makeSeed(const KeyElem< Key, Key2> &in){ return defaultHashFnc<Key>::makeSeed(in.first) ^ defaultHashFnc<Key2>::makeSeed(in.second) ;  }
};

template< > class defaultHashFnc<string>
{
public: static inline unsigned int makeSeed(const string __s){   unsigned long __h = 0; 
const char* tmp = __s.c_str();
for ( ; *tmp; ++tmp)
__h = 5*__h + *tmp;
return size_t(__h); }
};


template<unsigned int SIZE > class defaultHashFnc<Tuple<char,SIZE> >
{
public: static inline unsigned int makeSeed(const Tuple<char,4> __s) {   unsigned long __h = 0; 
	const char* tmp = & __s[0];
	for ( ; (unsigned int)(tmp - &__s[0]) < SIZE ; ++tmp)
		__h = 5*__h + *tmp;
return size_t(__h); }
};

template< > class defaultHashFnc<char*>
{
  public: static inline unsigned int makeSeed(const char* __s) {   unsigned long __h = 0; 
  for ( ; *__s; ++__s)
    __h = 5*__h + *__s;
  
  return size_t(__h); }
};

template< > class defaultHashFnc<const char*>
{
  public: static inline unsigned int makeSeed(const char* __s) {   unsigned long __h = 0; 
  for ( ; *__s; ++__s)
    __h = 5*__h + *__s;
  
  return size_t(__h); }
};

template< > class defaultHashFnc<char> { public: static inline unsigned int  makeSeed(char __x) { return __x; }};
template< > class defaultHashFnc<unsigned char> { public: static inline unsigned int makeSeed(unsigned char __x) { return __x; }};
template< > class defaultHashFnc<signed char> { public: static inline unsigned int makeSeed(unsigned char __x)  { return __x; }};
template< > class defaultHashFnc<short> { public: static inline unsigned int makeSeed(short __x)  { return __x; }};
template< > class defaultHashFnc<unsigned short> { public: static inline unsigned int makeSeed(unsigned short __x)  { return __x; }};
template< > class defaultHashFnc<int> { public: static inline unsigned int makeSeed(int __x)  { return __x; }};
template< > class defaultHashFnc<unsigned int> { public:  static inline unsigned int makeSeed(unsigned int __x)  { return __x; }};
template< > class defaultHashFnc<long> { public: static inline unsigned int makeSeed(long __x)  { return __x; }};
template< > class defaultHashFnc<unsigned long> { public: static inline unsigned int makeSeed(unsigned long __x)  { return __x; }};
template<class C> class defaultHashFnc<C*> { public: static inline unsigned int makeSeed(C* inv) {
	unsigned long __h = 0; 
	const char* tmp = (const char*)  &inv;
	for ( ; (unsigned int)(tmp - ((const char*)&inv)) < sizeof(C*) ; ++tmp)
		__h = 5*__h + *tmp;
	return size_t(__h); 
}};  // from address

template<class Key, class Data, class HashFnc> 
class myHashmap{
    unsigned int hashpos(unsigned int seed) const;
    void swap_and_pop(unsigned int ite);
    public:
    typedef unsigned int iterator;
    Vector< pair< KeyElem<Key, Data> , unsigned int > > heap;

    unsigned int* hash;

    unsigned char hash_mag;

    myHashmap(): hash_mag(0){}
    Data operator[](const Key &) const;
    Data& operator[](const Key &);
    unsigned int find(const Key &) const;

    Data& deref(unsigned int ite){return heap[ite].first.d;}
    const Data& deref(unsigned int ite) const {return heap[ite].first.d;}

    void erase(const Key &);
    void erase_from_iterator(unsigned int ite);
    inline unsigned int getSize() const{return heap.getSize();}
    void clear(){heap.clear(); if (hash_mag != 0) delete[](hash); hash_mag = 0; }

    void rehash(unsigned char _new_mag);

    void show(FILE*f = stdout, int level=0)const ;
};

template<class Key, class BackKey, class Data, class HashFnc, class BackHashFnc> 
class dualHashmap{
    unsigned int hashpos(unsigned int seed) const{return (seed & ((0xFFFFFFFF << hash_mag)^ 0xFFFFFFFF));}
    void swap_and_pop(unsigned int ite);
    public:
    typedef unsigned int iterator;
    Vector< pair< KeyElem< pair<Key, BackKey > , Data> ,   pair<unsigned int,unsigned int> > > heap;

    unsigned int* hash;
    unsigned char hash_mag;

    dualHashmap(): hash_mag(0){}
    void set(const Key &, const BackKey &); 
 //   void set(const Key &, const BackKey &, const Data&);

    unsigned int find(const Key &) const;
    unsigned int back_find(const BackKey &) const;

    Data& operator[](const Key &);
    Data& dualHash(const BackKey &);

    Data& deref(unsigned int ite){return heap[ite].first.d;}
    BackKey back_deref(unsigned int ite) const {return heap[ite].first.k.second;}
    Key front_deref(unsigned int ite) const {return heap[ite].first.k.first;}

    void erase(const Key &k){erase_from_iterator(find(k));}
    void back_erase(const BackKey &k){erase_from_iterator(back_find(k));}

    void erase_from_iterator(unsigned int ite);

    Key frontKey(const BackKey& b) const {return heap[back_find(b)].first.k.first;}
    BackKey backKey(const Key& k) const {return heap[find(k)].first.k.second;}

    void rehash(unsigned char _new_mag);

    void show(FILE*f = stdout, int level=0)const ;

    void test(FILE*f = stdout);
};



class TableReader{ // Table read *must* have a header row, and '+' '-' '*' '/' '(' ')' are treated as special characters, no collumn name should start with a digit too 
	unsigned int* read_info;
	unsigned int nbread;
	Vector<unsigned int> oper;
public:
	FILE *f;
	unsigned int nbcols;
	
	
	string* current_row;
	string* col_name;
	TableReader(const char*, const Vector<const char*> &collumns); // collumns that are integer representation (only digits) are used as collumn offsets
	~TableReader();
	const char* operator[](const unsigned int i ) const{return current_row[i].c_str();}
	bool nextRow();
};

// store serial elements, assumes that store elements
template<class IDTYPE>
class SerialStore{
	void compress(char* buffer);
public:
	map<IDTYPE, KeyElem<unsigned int, unsigned int> > hash_scope;
	FILE *f;
	char *path;
	unsigned int wasted;
	unsigned int size;
	
	SerialStore(const char* f_path);
	~SerialStore();
	
	template<class T> void load(IDTYPE, T &);
	template<class T> void save(IDTYPE, const T &);
	
	
    unsigned int itemSize(IDTYPE) const;
	template<class T> void load_arrayitem(IDTYPE, T &, unsigned int posis, unsigned int unitsize = sizeof(T)) const;
	template<class T> void save_arrayitem(IDTYPE, const T &, unsigned int posis, unsigned int unitsize = sizeof(T));
	void save_reservearray(IDTYPE, unsigned int size);
	
    FILE* getItemHandle(IDTYPE) const; 
	
	
	void flush(char* newpath = NULL); //
	bool has(IDTYPE) const;
	
	void show(FILE* f = stdout) const;
};


// Scope for batch Fourier Transforms, on arrays of ANY size
class Bluestein{
public:
	unsigned char pre_pow2;
	unsigned int mult;
	unsigned char post_pow2;
	Vector< Complex<double> > bluewindow;
	unsigned int tsize;

	Bluestein(){}
	Bluestein(unsigned int size){setSize(size);}
	void setSize(unsigned int n_size);

	unsigned int getBufferSize()const;

	template<class C> void toFT(C* data) const;
	template<class C> void toIFT(C* data) const;

	template<class C> void toFT_even(C* data) const;
	template<class C> void toIFT_even(C* data) const;

	template<class C, class D> void toFT(D* o_data,const C* i_data) const;
	template<class C, class D> void toIFT(D* o_data,const C* i_data) const;

	template<class C> void alloc_resample_buffer(C* &iobuf, Bluestein* invblue) const;
	template<class C> void resample(C* resample_buffer, Bluestein* invblue) const;

	void show(FILE* f = stdout, int level =0)const;
};




template<class C = unsigned char,unsigned int NBCHAN= 3u>
class CharArea{
	public:
	DataGrid<C ,3> image;
	Tuple<unsigned int , 4> inner_rect;
	Tuple<C, NBCHAN> color_bg;
	Tuple<C, NBCHAN> color_axe;
	Tuple<C, NBCHAN> color_axetext;

	Tuple<C, NBCHAN> color_min;
	Tuple<C, NBCHAN> color_max;

	void setDefaultStyle();
	void initialize(const Tuple<unsigned int,2> dims);
	void setAxes(double Xmin,double Ymin, double Xmax, double Ymax, bool Xlog= false, bool Ylog = false);

	template<class D, class E> void drawFrame(const DataGrid<D, 2> &fr,const E &min,const E &max);

	void overwriteFrame(const DataGrid<C, 3> &fr);

//	template<class D> DataGrid<C ,3> makeChart(const DataGrid<D,2> &image)const;

};

// GSL functions

#ifdef GNU_SCIENTIFIC_LIBRARY
int gsl_sf_lnfact_e(const unsigned int n, gsl_sf_result * result);
#endif

#include "HashMap.hpp"
#include "Vector.hpp"
#include "DataGrid.hpp"
#include "ExOp.hpp"


#include "primitive_basic.h"
#include "primitive_stats.h"
#include "primitive_tem.h"
} // end of namespace LFHPrimitive

#endif


