/*
 * primitive_exop.hpp
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


 // namespace LFHPrimitive{


#define LFH_typename_check(TyPeNaMe) \
template <class C>  int doexist_##TyPeNaMe( typename C::TyPeNaMe const * ); \
template <typename T>   char doexist_##TyPeNaMe( ... ); \
template<int hasfunc, class C> class isdef_##TyPeNaMe{public:typedef typename C::TyPeNaMe TYPE;}; \
template<class C> class isdef_##TyPeNaMe<sizeof(char), C>{public:typedef C TYPE; }

#define LFH_typename_check_dep(TyPeNaMe, WoUlDbEtYpE) \
template <class C>  int doexist_##TyPeNaMe( typename C::TyPeNaMe const * ); \
template <typename T>   char doexist_##TyPeNaMe( ... ); \
template<int hasfunc, class C> class isdef_##TyPeNaMe{public:typedef typename C::TyPeNaMe TYPE;}; \
template<class C> class isdef_##TyPeNaMe<sizeof(char), C>{public:typedef WoUlDbEtYpE TYPE; }

LFH_typename_check(DEF_TYPE);
LFH_typename_check(NEG_TYPE);
LFH_typename_check(TRJU_TYPE);
LFH_typename_check(LMUL_TYPE);
LFH_typename_check(INNER_TYPE);
LFH_typename_check(OUTER_TYPE);
LFH_typename_check(VECTOR_TYPE);
LFH_typename_check(DERIVATIVE_TYPE);

LFH_typename_check(REAL_TYPE);
LFH_typename_check_dep(COMPLEX_TYPE, Complex<C>);
LFH_typename_check_dep(QUATERNION_TYPE, Quaternion<C>);
LFH_typename_check_dep(ITERATOR_TYPE, void);
LFH_typename_check_dep(GAUS_TYPE, LFHCONCAT2(WeightElem<C,2>) );

LFH_typename_check_dep(IS_COMMUTATIVE, YESNO<false> );
LFH_typename_check_dep(IS_OWNED, YESNO<false> );
LFH_typename_check_dep(IS_POD, YESNO<false> );

	template<class A, unsigned int flag>
	class ExCo{
	public:

		typedef A TYPE;
        typedef typename isdef_DEF_TYPE<sizeof(doexist_DEF_TYPE<A>(0)),A>::TYPE DEF_TYPE;
		typedef typename isdef_NEG_TYPE<sizeof(doexist_NEG_TYPE<A>(0)),A>::TYPE NEG_TYPE;
		typedef typename isdef_TRJU_TYPE<sizeof(doexist_TRJU_TYPE<A>(0)),A>::TYPE TRJU_TYPE;
        typedef typename isdef_LMUL_TYPE<sizeof(doexist_LMUL_TYPE<A>(0)),A>::TYPE LMUL_TYPE;
		
		typedef typename isdef_INNER_TYPE<sizeof(doexist_INNER_TYPE<A>(0)),A>::TYPE INNER_TYPE;
        typedef typename isdef_VECTOR_TYPE<sizeof(doexist_VECTOR_TYPE<A>(0)),A>::TYPE VECTOR_TYPE;
		typedef typename isdef_DERIVATIVE_TYPE<sizeof(doexist_DERIVATIVE_TYPE<A>(0)),A>::TYPE DERIVATIVE_TYPE;
		
		typedef typename isdef_OUTER_TYPE<sizeof(doexist_OUTER_TYPE<A>(0)),A>::TYPE OUTER_TYPE; // outer product result, (increase in the number of dimentions!)
		

		typedef typename isdef_REAL_TYPE<sizeof(doexist_REAL_TYPE<A>(0)),A>::TYPE REAL_TYPE;
		typedef typename isdef_COMPLEX_TYPE<sizeof(doexist_COMPLEX_TYPE<A>(0)),A>::TYPE COMPLEX_TYPE;
		typedef typename isdef_QUATERNION_TYPE<sizeof(doexist_QUATERNION_TYPE<A>(0)),A>::TYPE QUATERNION_TYPE;
        typedef typename isdef_GAUS_TYPE<sizeof(doexist_GAUS_TYPE<A>(0)),A>::TYPE GAUS_TYPE;
		
        typedef typename isdef_IS_COMMUTATIVE<sizeof(doexist_IS_COMMUTATIVE<A>(0)),A>::TYPE IS_COMMUTATIVE;
        typedef typename isdef_IS_OWNED<sizeof(doexist_IS_OWNED<A>(0)),A>::TYPE IS_OWNED;
        typedef typename isdef_IS_POD<sizeof(doexist_IS_POD<A>(0)),A>::TYPE IS_POD;

		// typedef YESNO< A::IsPOD::ans > IsPOD;
		static const bool IsPOD = IS_POD::ans ;
		static const bool NeedsAddLink = IS_OWNED::ans; // containers needs to update addresses in link registers

		
		
		typedef A SAFETYPE;
		typedef unsigned char MAGNITUDE_TYPE;

		static A one(){A a; a.toOne(); return(a);}
		static A zero(){A a; a.toZero(); return(a);}
		static A random(){A a; a.toRand(); return(a);}
		static A intPow(const A& val, int power){
			switch(power){
				case -1: return(invert(val));
				case 0: return(one());
				case 1: return(val);
				case 2: return(val*val);
				case 3: return(val*val*val);

			}
		}

		static A invert(const A& val){ return(val);} // needed! if possible
		static A floPow(const A& val, float power){ return(floPow(val,(double)power));} // needed! if possible
		static A floPow(const A& val, double power){ return(val);} // needed! if possible
		static double sign(const A& val);
		template<int v> static A intPow(const  A& what) {return(intPow(v));}
		static double arctan(const A& r,const A& i);

		static void show(const A& val, FILE* out = stdout, int level = 0){val.show(out,level);}

		static void hasTransferFunction(){return(false);}// has function which move ressources, and leaves the origin ready for destruction

		static double pdist(const A& a, const A& b);
		static double norm(A& a){return ExCo<A>::norm(a);}
		static double norm_squarred(A& a){return ExCo<A>::norm_squarred(a);}

		static void save(const A& what, FILE *f) {
			//	if (hasfunc)what.save(f);
			//	else
			fwrite(&what,sizeof(A),1,f);
		} //
		static void load(A& what, FILE *f, unsigned int lenght){
			//		if (hasfunc) what.load(f, lenght);
			//	else
			fread(&what,sizeof(A),1,f);
		}

		// build missing operations from primitives!

		// builds

		inline static A mknegative(const A& a);
		inline static A mkinverse(const A& a);
		inline static A mksquare(const A& a);
		inline static A mkintpow(const A& what, const int pow);
		inline static A mkinvintpow(const A& what, const int pow); // uses mkinverse and mkintpow
		
		
		// assignments

		inline static A& toinverse(A& what);
		inline static A& tonegative(A& what);
		inline static void tomult(A& a, const A& b);
		inline static A& tosquare(A& a);

		inline static A mkadd(const A& a, const A&b);
		inline static A mksub(const A& a, const A&b);
		inline static A mkmult(const A& a, const A&b);
		inline static A mkdivi(const A& a, const A&b);
		
		// compare

		inline static bool isGT(const A &,const A &);
		inline static bool isGE(const A &,const A &);
		inline static bool isLT(const A &,const A &);
		inline static bool isLE(const A &,const A &);
		inline static bool isEQ(const A &,const A &);
		inline static bool isNQ(const A &,const A &);


		static void tointpow(A& what, const int pow);
		//	static double pdist(const A& a, const A& b){return ExOp::pnorm(a-b);}

	};

// dead-end for object manipulation, cant tell the array size at this point
template<class A>
class ExCo<A*, 0u>{
public:
	//typedef YESNO<false> IsPOD;
	static const bool IsPOD = false;
	static const bool NeedsAddLink = false;
	typedef ExCo<A*> SAFETYPE;
	typedef ExCo<A*> MAGNITUDE_TYPE;

	typedef ExCo<A*> COMPLEX_TYPE;
	typedef ExCo<A*> QUATERNION_TYPE;

	static void random(A * & a){a = NULL;}
	static void one(A * & a){a = NULL;}
	static void zero(A * & a){a =NULL;}
	static A* one(){return( NULL);}
	static A* zero(){return( NULL);}
	static A* random(){return(NULL);}

	static void save(A * & what, FILE *f) { what->save(f);}
	static void load(const A * & what, FILE *f, unsigned int lenght =0) { what->load(f);}
	static void show(const A * & what, FILE *f,int level);
	
    static A* & memmove(A * & tar, A * & scr){tar = scr; scr = NULL; return tar;}

	
};

// dead-end for object manipulation, cant tell the array size at this point
template< >
class ExCo<char*, 0u>{
public:
	//typedef YESNO<false> IsPOD;
	static const bool IsPOD = false;
	static const bool NeedsAddLink = false;
	typedef ExCo<char*> SAFETYPE;
	typedef char MAGNITUDE_TYPE;
	
	typedef char* COMPLEX_TYPE;
	typedef char* QUATERNION_TYPE;
	typedef void* DERIVATIVE_TYPE;
	
	static void random(char * & a){a = NULL;}
	static void one(char * & a){a = NULL;}
	static void zero(char * & a){a =NULL;}
	static char* one(){return( NULL);}
	static char* zero(){return( NULL);}
	static char* random(){return(NULL);}
	
	static void save(char* const & what, FILE *f) { fwrite(what,1,strlen(what)+1,f);}
	static void load(char* & what, FILE *f, unsigned int lenght =0) { char buffer[65536]; char* tmp = buffer; do {fread(tmp,sizeof(unsigned char), 1,f); } while(*(tmp++) != '\0'); what = new char[(unsigned int)(tmp - buffer)]; memcpy(what, buffer, sizeof(char) * (unsigned int)(tmp - buffer));}

	static void show(char* const & what, FILE *f,int level) { if (level == 0) fprintf(f,"%s\n",what); else fprintf(f,"%s",what);}
	
    static char* & memmove(char* & tar, char* & scr){tar = scr; scr = NULL; return tar;}
	
};


template<class A, unsigned int SIZE>
	class ExCo<A[SIZE], 0u>{
	public:
		typedef typename ExCo<A>::SAFETYPE SAFETYPE;
		typedef typename ExCo<A>::MAGNITUDE_TYPE MAGNITUDE_TYPE;


        typedef typename isdef_DEF_TYPE<sizeof(doexist_DEF_TYPE<A>(0)),A>::TYPE DEF_TYPE[SIZE];
		typedef typename isdef_NEG_TYPE<sizeof(doexist_NEG_TYPE<A>(0)),A>::TYPE NEG_TYPE[SIZE];
        typedef DataGrid<A,2> TRJU_TYPE;
        typedef DataGrid<A,2> LMUL_TYPE;
	//	typedef typename ExCo<A>::DERIVATIVE_TYPE DERIVATIVE_TYPE[SIZE];
        typedef A INNER_TYPE;

		typedef typename ExCo<A>::COMPLEX_TYPE* COMPLEX_TYPE;
		typedef typename ExCo<A>::QUATERNION_TYPE* QUATERNION_TYPE;

        template <class S> class SUBS_INNER{public: typedef Tuple<S, SIZE> TYPE;};

		static const bool IsPOD = ExCo<A>::IsPOD;
		static const bool NeedsAddLink = ExCo<A>::NeedsAddLink;


		static void zero(A (&a)[SIZE]);
		static void random(A (&a)[SIZE]);
        static void one(A (&a)[SIZE]);
		static void minimum(A (&a)[SIZE]);
        static void maximum(A (&a)[SIZE]);

		static void show(const A (&a)[SIZE], FILE* out=stdout, int level=0);

		static bool isValid(A (&a)[SIZE]);

		static const Tuple<A, SIZE>& toadd(A (&a)[SIZE], const A (&b)[SIZE]);
		static const Tuple<A, SIZE>& tosub(A (&a)[SIZE], const A (&b)[SIZE]);
		static const Tuple<A, SIZE>& tomult(A (&a)[SIZE], const A (&b)[SIZE]);
		static const Tuple<A, SIZE>& todivi(A (&a)[SIZE], const A (&b)[SIZE]);
        

		
		/*
		 static const bool IsPOD = true;
		 static const bool hasSpeMemMove = false;


		 static double zero(){return(0.0f);}
		 static void zero(double &a){a =0.0f;}
		 static double one(){return(1.0f);}
		 static void one(double &a){a =1.0f;}
		 static double random(){return((((unsigned int)rand()) ^ ((unsigned int)(rand() << 12)) ^ ((unsigned int)(rand() << 24))) * pow(0.5f, 32));  };
		 static void random(double &a){a = (((unsigned int)rand()) ^ ((unsigned int)(rand() << 12)) ^ ((unsigned int)(rand() << 24))) * pow(0.5f, 32);  };

		 static bool isValid(const double& a){
		 return(a + (-a) == 0);
		 }
		 static bool isnegative(const double& a){
		 return(a< 0.0f);
		 }
		 static double max(){return(DBL_MAX);}
		 static double min(){return(DBL_MIN);}

		 static double invert(const double& val){return(1.0f / val);}
		 static double intPow(const double& val, int power){return(pow(val, (double) power));}
		 static double floPow(const double& val, double power){return(pow(val, power));}
		 static double sign(const double& val){return(val >= 0 ? 1.0f : -1.0f);}
		 static double invintPow(const double& val, int power){return(pow(val, 1.0f / ((double)power)));}

		 static double arctan(const double& r,const double& i){return(atan2(r,i));}

		 static double doubleratio(const double& num,const double& den){return(num/den);}



		 static void show(const double &val, FILE* out=stdout, int level=0){fprintf(out, "%e", val); if (level == 0) fprintf(out, "\n"); }

		 static void save(const double& what, FILE *f) { fwrite(&what,sizeof(double), 1,f);}
		 static void load(double & what, FILE *f) { fread(&what,sizeof(double), 1,f);}
		 static void memmove(double& a,double& o){a=o;}

		 static double norm(const double& a){return fabs(a);}
		 static double pnorm(const double& a){return a*a ;}
		 */
	};


	template< class R, class A>
	class ExCo< R (*)(A), 0u >{ // function that returns A with arg B
		public:
		static const bool IsPOD = false;
		static const bool NeedsAddLink = false; // containers needs to update addresses in link registers

		typedef void COMPLEX_TYPE;
		typedef void QUATERNION_TYPE;

		typedef ExCo< R (*)(A) > SAFETYPE;

		template< class C > class RETT{ // is a simple function: the return type is always the same for any input
		public:
			typedef R TYPE;
			C junkdata;
		};

		R operator()(A a); // fake declaration! (used to tell if input arguments are legal for function

	//	R operator()(A&){} // fake declaration!
	//	R operator()(const A&){} // fake declaration!
	//	R operator()(A) const{} // fake declaration!
	//	R operator()(A&) const{} // fake declaration!
	//	R operator()(const A&) const{} // fake declaration!

	};




	template<class A, class B>
	class ExCo<pair<A,B>, 0u >{
	public:
		typedef typename ExCo<A>::SAFETYPE SAFETYPE;
		typedef pair< typename ExCo<A>::MAGNITUDE_TYPE, typename ExCo<B>::MAGNITUDE_TYPE> MAGNITUDE_TYPE;
		static const bool NeedsAddLink = ExCo<A>::NeedsAddLink ||  ExCo<B>::NeedsAddLink ; // containers needs to update addresses in link registers
		typedef pair< typename ExCo<A>::COMPLEX_TYPE, typename ExCo<B>::COMPLEX_TYPE>  COMPLEX_TYPE;
		typedef pair< typename ExCo<A>::QUATERNION_TYPE, typename ExCo<B>::QUATERNION_TYPE> QUATERNION_TYPE;


        inline static bool isZero(const pair<A,B> &a);
        inline static bool isOne(const pair<A,B> &a);

        static void zero(pair<A,B> &a);
        static void one(pair<A,B> &a);
        static void random(pair<A,B> &a);
		static bool isValid(const pair<A,B> &a);

		static const bool IsPOD = true;
		static void show(const pair<A,B> &a, FILE* out=stdout, int level=0);
        static string type_tostring(const pair<A,B> &a);

	};

template<class A>
class ExCo< vector<A>, 0u >{
public:
	typedef typename ExCo<A>::SAFETYPE SAFETYPE;
	typedef vector< typename ExCo<A>::MAGNITUDE_TYPE> MAGNITUDE_TYPE;
	static const bool NeedsAddLink = ExCo<A>::NeedsAddLink ; // containers needs to update addresses in link registers
	typedef vector< typename ExCo<A>::COMPLEX_TYPE>  COMPLEX_TYPE;
	typedef vector< typename ExCo<A>::QUATERNION_TYPE> QUATERNION_TYPE;
	static const bool IsPOD = false;	
	
	inline static bool isZero(const vector<A> &a);
	inline static bool isOne(const vector<A> &a);
	inline static bool isValid(const vector<A> &a);
	
/*	static void zero(vector<A> &a);
	static void one(vector<A> &a);
	static void random(vector<A> &a);

	static void show(const vector<A> &a, FILE* out=stdout, int level=0);
	static string type_tostring(const vector<A> &a);*/
};

template<class A>
class ExCo<const A, 0u>{
public:
	typedef typename ExCo<A>::SAFETYPE SAFETYPE;
	typedef typename ExCo<A>::MAGNITUDE_TYPE MAGNITUDE_TYPE;

	typedef typename ExCo<A>::COMPLEX_TYPE  COMPLEX_TYPE;
	typedef typename ExCo<A>::QUATERNION_TYPE QUATERNION_TYPE;

	static const bool NeedsAddLink = ExCo<A>::NeedsAddLink; // containers needs to update addresses in link registers
	static const bool IsPOD = true;


	static inline void show(const A &a, FILE* out=stdout, int level=0);
	static inline bool isValid(const A& a);

};


#define PRIMTYPE_CONSTRUCT_FUNCTION(cLaSs) \
static inline cLaSs maximum(){cLaSs a; maximum(a); return(a);} \
static inline cLaSs minimum(){cLaSs a; minimum(a); return(a);} \
static inline cLaSs zero(){cLaSs a; zero(a); return(a);} \
static inline cLaSs one(){cLaSs a; one(a); return(a);} \
static inline cLaSs random(){cLaSs a; random(a); return(a);} \
static inline cLaSs epsilon(){cLaSs a; epsilon(a); return(a);} /* minimum value that changes the representation of 1 when added */ \
static inline cLaSs delta(){cLaSs a; delta(a); return(a);} /* minimum value that is larger than 0 */ \
typedef cLaSs DEF_TYPE; \
typedef cLaSs NEG_TYPE; \
typedef cLaSs TRJU_TYPE; \
typedef cLaSs OUTER_TYPE; \
typedef cLaSs LMUL_TYPE; \
typedef cLaSs REAL_TYPE; \
typedef Complex<cLaSs>  COMPLEX_TYPE; \
typedef Quaternion<cLaSs> QUATERNION_TYPE; \
typedef GaussScope<double, double> GAUS_TYPE; \
typedef YESNO<true> IS_COMMUTATIVE; \
typedef YESNO<false> IS_OWNED; \
typedef YESNO<true> IS_POD; \
typedef cLaSs VECTOR_TYPE;\
static inline GaussScope<double, double> mkgaussstat(const cLaSs &a, double weight = 1.0f){return (GaussScope<double, double>((double)a, ((double)a)*((double)a),weight) );} \
static inline cLaSs mkrealproj(const cLaSs& a){return a;}\
static inline cLaSs mkimmaproj(const cLaSs& a){return zero();}\
static inline cLaSs mkjmmaproj(const cLaSs& a){return zero();}\
static inline cLaSs mkkmmaproj(const cLaSs& a){return zero();}

#define SIMPLE_FUNCTIONS(cLaSs)                                        \
static inline cLaSs mkadd(const cLaSs &a, const cLaSs &b){return a+b;} \
static inline cLaSs mksub(const cLaSs &a, const cLaSs &b){return a-b;} \
static inline cLaSs mkmult(const cLaSs &a, const cLaSs &b){return a*b;} \
static inline cLaSs mkdivi(const cLaSs &a, const cLaSs &b){return a/b;} \
static inline cLaSs mknegative(const cLaSs &a) {return -a;}\
static inline cLaSs& toadd(cLaSs &a, const cLaSs &b){return (a +=b);} \
static inline cLaSs& tosub(cLaSs &a, const cLaSs &b){return (a -=b);} \
static inline cLaSs& tomult(cLaSs &a, const cLaSs &b){return (a *=b);} \
static inline cLaSs& todivi(cLaSs &a, const cLaSs &b){return (a /=b);} \
static inline setcomparison setcmp(const cLaSs &a, const cLaSs &b) {return ((a > b) ? SETCMP_GT : ((a == b) ? SETCMP_EQUAL : SETCMP_LT));} \
static inline bool isGT(const cLaSs &a, const cLaSs &b){return a > b;} \
static inline bool isGE(const cLaSs &a, const cLaSs &b){return a >= b;}\
static inline bool isLT(const cLaSs &a, const cLaSs &b){return a < b;} \
static inline bool isLE(const cLaSs &a, const cLaSs &b){return a <= b;}\
static inline bool isEQ(const cLaSs &a, const cLaSs &b){return a == b;}\
static inline bool isNQ(const cLaSs &a, const cLaSs &b){return a != b;}



	template< >
	class ExCo<double, 0u>{
	public:
		typedef short MAGNITUDE_TYPE;
		typedef void DERIVATIVE_TYPE;
		static const bool IsPOD = true;
		static const bool hasSpeMemMove = false;

		static const bool NeedsAddLink = false; // containers needs to update addresses in link registers
		typedef ExOp SAFETYPE;

		static bool isZero(const double &a){return a == 0.0f;}
		static bool isOne(const double &a){return a == 1.0f;}
		static void zero(double &a){a =0.0f;}
		static void one(double &a){a =1.0f;}
		static void random(double &a){a = (((unsigned int)rand()) ^ ((unsigned int)(rand() << 12)) ^ ((unsigned int)(rand() << 24))) * pow(0.5f, 32);  };
		static void epsilon(double &a){a =DBL_EPSILON;} // minimum value that changes the representation of 1 when added
		static void delta(double &a){a = DBL_MIN;}  // minimum value that is larger than 0
		
		static void maximum(double &a){a = DBL_MAX;}
		static void minimum(double &a){a = -DBL_MAX;}

		PRIMTYPE_CONSTRUCT_FUNCTION(double)

		inline static double& toinverse(double &a){return (a = 1.0f / a);}
		inline static double& tonegative(double &a){return (a = -a);}
		inline static double& tosquarre(double &a){ return(a *= a); }
		
		inline static double mkinverse(const double &a){return(1.0f / a);}
		inline static double mkintpow(const double& what, const int power){ return( pow(what, (double) power)); }
		inline static double mkinvintpow(const double& what, const int power){ return( pow(what, 1.0f / ((double) power))); }

		SIMPLE_FUNCTIONS(double)

		static bool isValid(const double& a){
			return(a + (-a) == 0);
		}
		static bool isnegative(const double& a){
			return(a< 0.0f);
		}
		static double max(){return(DBL_MAX);}
		static double min(){return(DBL_MIN);}

		static double invert(const double& val){return(1.0f / val);}
		static double intPow(const double& val, int power){return(pow(val, (double) power));}
		static double floPow(const double& val, double power){return(pow(val, power));}
		static double sign(const double& val){return(val >= 0 ? 1.0f : -1.0f);}
		static double invintPow(const double& val, int power){return(pow(val, 1.0f / ((double)power)));}

		static double arctan(const double& r,const double& i){return(atan2(r,i));}

		static double doubleratio(const double& num,const double& den){return(num/den);}



		static void show(const double &val, FILE* out=stdout, int level=0){fprintf(out, "%e", val); if (level == 0) fprintf(out, "\n"); }
        static string type_tostring(const double& a){return string("double");}

		static void save(const double& what, FILE *f) { fwrite(&what,sizeof(double), 1,f);}
		static void load(double & what, FILE *f, unsigned int lenght=0) { fread(&what,sizeof(double), 1,f);}
		static void memmove(double& a,double& o){a=o;}

		static double lognorm(const double& a){return log(fabs(a));}
		static double norm(const double& a){return fabs(a);}
		static double pnorm(const double& a){return a*a;}


		template<unsigned int SIZE> static unsigned int getHyperDirection(const double (&a)[SIZE], const double (&b)[SIZE]){ // direction for hyper-ordering of arrays
			unsigned int best =0;

			typename StConstList<16>::UINT expb, expc;

			typename StConstList<64>::UINT xorb = (*((typename StConstList<64>::UINT *) &a)) ^ (*((typename StConstList<64>::UINT *) &b)); // xor and make unsigned
				// if exponent are different, the max eponent wins
			    // if exponent are the same


			if ( ((typename StConstList<16>::UINT *)&xorb)[3] & 0xFFF0){
				if ( ((typename StConstList<16>::UINT *)&xorb)[3] & 0x8000) return(0);
				// most significant is max exponent

				((typename StConstList<16>::UINT *)&xorb)[3] = 0x0010; // most-significant turned on!
				expb = (((typename StConstList<16>::UINT *)&a)[3] > ((typename StConstList<16>::UINT *)&b)[3] ? ((typename StConstList<16>::UINT *)&a)[3] : ((typename StConstList<16>::UINT *)&b)[3]) >> 2;
				// gets greater exponent

			}else{
				if (xorb != 0) expb = (((typename StConstList<16>::UINT *)&a)[3]) >> 2;
				else  expb = 0;
				// gets exponent (same for both numbers)
			}




			for(unsigned int cur=1;cur<SIZE;cur++) {
			typename StConstList<64>::UINT xorc = (*((typename StConstList<64>::UINT *) &(a + cur))) ^ (*((typename StConstList<64>::UINT *) &(b + cur))); // xor and make unsigned
				if ( ((typename StConstList<16>::UINT *)&xorc)[3] & 0xFFF0){
					if ( ((typename StConstList<16>::UINT *)&xorc)[3] & 0x8000) return(cur);

					expc = (((typename StConstList<16>::UINT *)&(a+cur))[3] > ((typename StConstList<16>::UINT *)&(b+cur))[3] ? ((typename StConstList<16>::UINT *)&(a+cur))[3] : ((typename StConstList<16>::UINT *)&(b+cur))[3]) >> 2;
					// gets greater exponent
					((typename StConstList<16>::UINT *)&xorc)[3] = 0x0010;


				}else expc = (((typename StConstList<16>::UINT *)& (a +cur) )[3]) >> 2;

				if (expc > expb){
					if (expc > expb + 52){
						best = cur; expb = expc; xorb = xorc;
					}else{
						xorb >>=  expc - expb;
						if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
						expb = expc;
					}
				}else{
					if (expb <= expc + 52){  // guarrantied fail otherwise
						xorc >>=  expb - expc;
						if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
					}
				}

			}

			return(best);
		}



	};

template< >
class ExCo<float>{
public:
	typedef short MAGNITUDE_TYPE;
	typedef void DERIVATIVE_TYPE;
	
	static const bool IsPOD = true;
	static const bool hasSpeMemMove = false;

	static const bool NeedsAddLink = false; // containers needs to update addresses in link registers
	typedef ExOp SAFETYPE;

	static bool isZero(const float &a){return a == 0.0f;}
	static bool isOne(const float &a){return a == 1.0f;}
	static void zero(float &a){a =0.0f;}
	static void one(float &a){a =1.0f;}
	static void random(float &a){a = (((unsigned int)rand()) ^ ((unsigned int)(rand() << 12)) ^ ((unsigned int)(rand() << 24))) * pow(0.5f, 32);  };
	static void epsilon(float &a){a =FLT_EPSILON;}
	static void delta(float &a){a = FLT_MIN;} 
	static void maximum(float &a){a = FLT_MAX;}
	static void minimum(float &a){a = -FLT_MAX;}

	PRIMTYPE_CONSTRUCT_FUNCTION(float)

	static float mkinverse(const float &a){return(1.0f / a);}
	static float mkintpow(const float& what, const int power){ return( pow(what, (double) power)); }
	static float mkinvintpow(const float& what, const int power){ return( pow(what, 1.0f / ((double) power))); }



	SIMPLE_FUNCTIONS(float)

	static bool isValid(const float& a){
		return(a + (-a) == 0);
	}
	static bool isnegative(const float& a){
		return(a< 0.0f);
	}
	static float max(){return(DBL_MAX);}
	static float min(){return(DBL_MIN);}

	static float invert(const float& val){return(1.0f / val);}
	static float intPow(const float& val, int power){return(pow(val, (float) power));}
	static float floPow(const float& val, float power){return(pow(val, power));}
	static float sign(const float& val){return(val >= 0 ? 1.0f : -1.0f);}
	static float invintPow(const float& val, int power){return(pow(val, 1.0f / ((float)power)));}
	static float arctan(const float& r,const float& i){return(atan2(r,i));}
	static float floatratio(const float& num,const float& den){return(num/den);}
	static void show(const float &val, FILE* out=stdout, int level=0){fprintf(out, "%e", val); if (level == 0) fprintf(out, "\n"); }
    static string type_tostring(const float& a){return string("float");}
	static void save(const float& what, FILE *f) { fwrite(&what,sizeof(float), 1,f);}
	static void load(float & what, FILE *f, unsigned int lenght=0) { fread(&what,sizeof(float), 1,f);}
	static void memmove(float& a,float& o){a=o;}

	static double lognorm(const float& a){return log(fabs(a));}
	static double norm(const float& a){return fabs(a);}
	static double pnorm(const float& a){return a*a ;}


	template<unsigned int SIZE> static unsigned int getHyperDirection(const float (&a)[SIZE], const float (&b)[SIZE]){ // direction for hyper-ordering of arrays
		unsigned int best =0;

		typename StConstList<16>::UINT expb, expc;

		typename StConstList<64>::UINT xorb = (*((typename StConstList<64>::UINT *) &a)) ^ (*((typename StConstList<64>::UINT *) &b)); // xor and make unsigned
		// if exponent are different, the max eponent wins
		// if exponent are the same


		if ( ((typename StConstList<16>::UINT *)&xorb)[3] & 0xFFF0){
			if ( ((typename StConstList<16>::UINT *)&xorb)[3] & 0x8000) return(0);
			// most significant is max exponent

			((typename StConstList<16>::UINT *)&xorb)[3] = 0x0010; // most-significant turned on!
			expb = (((typename StConstList<16>::UINT *)&a)[3] > ((typename StConstList<16>::UINT *)&b)[3] ? ((typename StConstList<16>::UINT *)&a)[3] : ((typename StConstList<16>::UINT *)&b)[3]) >> 2;
			// gets greater exponent

		}else{
			if (xorb != 0) expb = (((typename StConstList<16>::UINT *)&a)[3]) >> 2;
			else  expb = 0;
			// gets exponent (same for both numbers)
		}




		for(unsigned int cur=1;cur<SIZE;cur++) {
			typename StConstList<64>::UINT xorc = (*((typename StConstList<64>::UINT *) &(a + cur))) ^ (*((typename StConstList<64>::UINT *) &(b + cur))); // xor and make unsigned
			if ( ((typename StConstList<16>::UINT *)&xorc)[3] & 0xFFF0){
				if ( ((typename StConstList<16>::UINT *)&xorc)[3] & 0x8000) return(cur);

				expc = (((typename StConstList<16>::UINT *)&(a+cur))[3] > ((typename StConstList<16>::UINT *)&(b+cur))[3] ? ((typename StConstList<16>::UINT *)&(a+cur))[3] : ((typename StConstList<16>::UINT *)&(b+cur))[3]) >> 2;
				// gets greater exponent
				((typename StConstList<16>::UINT *)&xorc)[3] = 0x0010;


			}else expc = (((typename StConstList<16>::UINT *)& (a +cur) )[3]) >> 2;

			if (expc > expb){
				if (expc > expb + 52){
					best = cur; expb = expc; xorb = xorc;
				}else{
					xorb >>=  expc - expb;
					if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
					expb = expc;
				}
			}else{
				if (expb <= expc + 52){  // guarrantied fail otherwise
					xorc >>=  expb - expc;
					if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
				}
			}

		}

		return(best);
	}



};

	template< >
	class ExCo<char>{
	public:
		typedef unsigned char MAGNITUDE_TYPE;
		typedef void DERIVATIVE_TYPE;
		static bool const IsPOD = true;
		static bool const hasSpeMemMove = false;
		static const bool NeedsAddLink = false; // containers needs to update addresses in link registers
		typedef ExOp SAFETYPE;


        static bool isZero(const char &a){return a == '\0';}
        static bool isOne(const char &a){return a == 1;}

		static void zero(char & a){ a = '\0' ;}
		static void one(char & a){ a = 1 ;}
		static void random(char & a){ a = (rand() & 255) ;}
		static void epsilon(char & a){ a = 1;}
		static void delta(char & a){ a = 1;}
		static void minimum(char & a){ a = -128;}
		static void maximum(char & a){ a = 127;}


		PRIMTYPE_CONSTRUCT_FUNCTION(char)


		static double sign(const int& val){return(val >= 0 ? 1.0f : -1.0f);}

		static void show(const char &val, FILE* out=stdout, int level=0){fprintf(out, "%c", val); if (level == 0) fprintf(out, "\n"); }
        static string type_tostring(const char& a){return string("char");}

		static void save(const char& what, FILE *f) { fwrite(&what,sizeof(char), 1,f);}
		static void load(char& what, FILE *f, unsigned int lenght=0) { fread(&what,sizeof(char), 1,f);}
		static void memmove(char& a ,char& o){a=o;}
		static void bitreverse(char& a){
			a = ((a << 4 ) & 0xF0) | ((a >> 4)& 0x0F);
			a = ((a << 2 ) & 0xCC) | ((a >> 2)& 0x33);
			a = ((a << 1 ) & 0xAA) | ((a >> 1)& 0x55);
		}
		static void bytereverse(char& a){}

		template<unsigned int SIZE> static unsigned int getHyperDirection(const char (&a)[SIZE], const char (&b)[SIZE]){ // direction for hyper-ordering of arrays
			unsigned char xorb,xorc;
			*(char*)(&xorb) = a[0] ^ b[0]; // xor and make unsigned
			unsigned int best =0;
			for(unsigned int cur=1;cur<SIZE;cur++) {
				*(char*)(&xorc) = a[cur] ^ b[cur]; // xor and make unsigned
				if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
			}
			return(best);
		}

	};


	template< >
	class ExCo<unsigned char>{
	public:
		typedef unsigned char MAGNITUDE_TYPE;
		typedef void DERIVATIVE_TYPE;
		static const bool IsPOD = true;
		//	typedef YESNO<true> IsPOD;
		static bool const hasSpeMemMove = false;
		static const bool NeedsAddLink = false; // containers needs to update addresses in link registers
		typedef ExOp SAFETYPE;

        static bool isZero(const unsigned char &a){return a == '\0';}
        static bool isOne(const unsigned char &a){return a == 1;}

		static void zero(unsigned char & a){ a = 0 ;};
		static void one(unsigned char & a){ a = 1 ;};
		static void random(unsigned char & a){ a =  (rand() & 255) ;};
		static void epsilon(unsigned char & a){ a = 1;}
		static void delta(unsigned char & a){ a = 1;}
		static void minimum(unsigned char & a){ a = 0;}
		static void maximum(unsigned char & a){ a = 255;}

		PRIMTYPE_CONSTRUCT_FUNCTION(unsigned char)

		static void show(const unsigned char &val, FILE* out=stdout, int level=0){fprintf(out, "%c", val); if (level == 0) fprintf(out, "\n"); }
        static string type_tostring(const unsigned char & a){return string("unsigned char");}
		static void save(const unsigned char& what, FILE *f) { fwrite(&what,sizeof(unsigned char), 1,f);}
		static void load( unsigned char& what, FILE *f, unsigned int lenght=0) { fread(&what,sizeof(unsigned char), 1,f);}


		static void memmove(unsigned char& a ,unsigned char& o){a=o;}
		static unsigned char upperbound_pow_of_2(const unsigned char &i){
			if (i <= 16){
				if (i <= 4){
					if (i > 2) return(2);
					else if (i < 2) return(0);
					else return(0);
				}else return( (i <=8) ? 3 : 4 );
			}else return ((i <= 64) ? ((i <=32) ? 5 : 6) : ((i <=128) ? 7 : 8 ));
		}
		static void bitreverse(unsigned char& a){
			a = ((a << 4 ) & 0xF0) | ((a >> 4)& 0x0F);
			a = ((a << 2 ) & 0xCC) | ((a >> 2)& 0x33);
			a = ((a << 1 ) & 0xAA) | ((a >> 1)& 0x55);
		}
		static void bytereverse(unsigned char& a){}
		template<unsigned int SIZE> static unsigned int getHyperDirection(const unsigned char (&a)[SIZE], const unsigned char (&b)[SIZE]){ // direction for hyper-ordering of arrays
			unsigned char xorb = a[0] ^ b[0]; // xor and make unsigned
			unsigned int best =0;
			for(unsigned int cur=1;cur<SIZE;cur++) {
				unsigned char xorc = a[cur] ^ b[cur]; // xor and make unsigned
				if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
			}
			return(best);
		}

	};

	template< >
	class ExCo<int>{
	public:
		typedef unsigned char MAGNITUDE_TYPE;
		typedef void DERIVATIVE_TYPE;
		static const bool IsPOD = true;
		static const bool NeedsAddLink = false; // containers needs to update addresses in link registers
		//	typedef YESNO<true> IsPOD;
		typedef ExOp SAFETYPE;
		static bool const hasSpeMemMove = false;

		
		PRIMTYPE_CONSTRUCT_FUNCTION(int)
		
        static bool isZero(const int &a){return a == 0;}
        static bool isOne(const int &a){return a == 1;}

		static void zero(int & a){a = 0;}
		static void one(int & a){a = 1;}
		static void random(int & a ){ a = rand() ^ (rand() << 12) ^ (rand() << 24);};

		static double sign(const int& val){return(val >= 0 ? 1.0f : -1.0f);}
		static void show(const int &val, FILE* out=stdout, int level=0){fprintf(out, "%i", val); if (level == 0) fprintf(out, "\n"); }
        static string type_tostring(const int & a){return string("int");}

		static void minimum(int& a) {a=0x80000000;}
		static void maximum(int& a) {a=0x7FFFFFFF;}
		static void epsilon(int& a) {a=1;}
		static void delta(int& a) {a=1;}
		
		
		static void save(const int& what, FILE *f) { fwrite(&what,sizeof(int), 1,f);}
		static void load(int& what, FILE *f, unsigned int lenght=0) { fread(&what,sizeof(int), 1,f);}

		static void memmove(int& a ,int& o){a=o;}
		static void bitreverse(int& a){
			a = ((a << 16 ) & 0xFFFF0000) | ((a >> 16)& 0x0000FFFF);
			a = ((a << 8 ) & 0xFF00FF00) | ((a >> 8)& 0x00FF00FF);
			a = ((a << 4 ) & 0xF0F0F0F0) | ((a >> 4)& 0x0F0F0F0F);
			a = ((a << 2 ) & 0xCCCCCCCC) | ((a >> 2)& 0x33333333);
			a = ((a << 1 ) & 0xAAAAAAAA) | ((a >> 1)& 0x55555555);
		}
		static void bytereverse(int& a){
			a = ((a << 16 ) & 0xFFFF0000) | ((a >> 16)& 0x0000FFFF);
			a = ((a << 8 ) & 0xFF00FF00) | ((a >> 8)& 0x00FF00FF);
		}
		template<unsigned int SIZE> static unsigned int getHyperDirection(const int (&a)[SIZE], const int (&b)[SIZE]){ // direction for hyper-ordering of arrays
			unsigned int xorb,xorc;
			*(unsigned int*)(&xorb) = a[0] ^ b[0]; // xor and make unsigned
			unsigned int best =0;
			for(unsigned int cur=1;cur<SIZE;cur++) {
				*(unsigned int*)(&xorc) = a[cur] ^ b[cur]; // xor and make unsigned
				if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
			}
			return(best);
		}

	};


	template< >
	class ExCo<unsigned int>{
	public:
		typedef unsigned char MAGNITUDE_TYPE;
typedef void DERIVATIVE_TYPE;
		static const bool IsPOD = true;
		//	typedef YESNO<true> IsPOD;
		static bool const hasSpeMemMove = false;
		static const bool NeedsAddLink = false; // containers needs to update addresses in link registers
		typedef ExOp SAFETYPE;

		PRIMTYPE_CONSTRUCT_FUNCTION(unsigned int)

		
        static bool isZero(const unsigned int &a){return a == 0;}
        static bool isOne(const unsigned int &a){return a == 1;}

		static void zero(unsigned int & a){a = 0;}
		static void one(unsigned int & a){a = 1;}
		static void random(unsigned int & a ){ a = rand() ^ (rand() << 12) ^ (rand() << 24);};

		static void minimum(unsigned int& a) {a=0;}
		static void maximum(unsigned int& a) {a=0xFFFFFFFF;}
		static void epsilon(unsigned int& a) {a=1;}
		static void delta(unsigned int& a) {a=1;}

		static void show(const unsigned int&val, FILE* out=stdout, int level=0){fprintf(out, "%u", val); if (level == 0) fprintf(out, "\n"); }
        static string type_tostring(const unsigned int & a){return string("unsigned int");}

		static void save(const unsigned int& what, FILE *f) { fwrite(&what,sizeof(unsigned int), 1,f);}
		static void load( unsigned int& what, FILE *f, unsigned int lenght=0) { fread(&what,sizeof(unsigned int), 1,f);}
		static void memmove(unsigned int& a ,unsigned int& o){a=o;}


		static unsigned char upperbound_pow_of_2(const unsigned int &i){
			if (i <= 65536){
				if (i <= 256){
					if (i <= 16){
						if (i <= 4){
							if (i > 2) return(2);
							else if (i < 2) return(0);
							else return(0);
						}else return( (i <=8) ? 3 : 4 );
					}else return((i <= 64) ? ((i <=32) ? 5 : 6) : ((i <=128) ? 7 : 8 ));
				}else return( (i<= 0x1000) ? ((i <= 0x0400) ? ((i <=0x200) ? 9 : 10) : ((i <=0x0800) ? 11 : 12 )) : ((i <= 0x4000) ? ((i <=0x2000) ? 13 : 14) : ((i <=0x8000) ? 15 : 16 )) );
			}else return( (i<= 0x01000000) ? ((i<= 0x00100000) ? ((i <= 0x00040000) ? ((i <=0x00020000) ? 17 : 18) : ((i <=0x00080000) ? 19 : 20 )) : ((i <= 0x00400000) ? ((i <=0x00200000) ? 21 : 22) : ((i <=0x00800000) ? 23 : 24 )) ) : ((i<= 0x10000000) ? ((i <= 0x04000000) ? ((i <=0x02000000) ? 25 : 26) : ((i <=0x08000000) ? 27 : 28 )) : ((i <= 0x40000000) ? ((i <=0x20000000) ? 29 : 30) : ((i <=0x80000000) ? 31 : 32 )) ));
		}


		static void bitreverse(unsigned int& a){
			a = ((a << 16 ) & 0xFFFF0000) | ((a >> 16)& 0x0000FFFF);
			a = ((a << 8 ) & 0xFF00FF00) | ((a >> 8)& 0x00FF00FF);
			a = ((a << 4 ) & 0xF0F0F0F0) | ((a >> 4)& 0x0F0F0F0F);
			a = ((a << 2 ) & 0xCCCCCCCC) | ((a >> 2)& 0x33333333);
			a = ((a << 1 ) & 0xAAAAAAAA) | ((a >> 1)& 0x55555555);
		}
		static void bytereverse(unsigned int& a){
			a = ((a << 16 ) & 0xFFFF0000) | ((a >> 16)& 0x0000FFFF);
			a = ((a << 8 ) & 0xFF00FF00) | ((a >> 8)& 0x00FF00FF);
		}
		template<unsigned int SIZE> static unsigned int getHyperDirection(const unsigned int (&a)[SIZE], const unsigned int (&b)[SIZE]){ // direction for hyper-ordering of arrays
			unsigned int xorb = a[0] ^ b[0]; // xor and make unsigned
			unsigned int best =0;
			for(unsigned int cur=1;cur<SIZE;cur++) {
				unsigned int xorc = a[cur] ^ b[cur]; // xor and make unsigned
				if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
			}
			return(best);
		}

	};
	template< >
	class ExCo<short>{
	public:
		typedef unsigned char MAGNITUDE_TYPE;
		typedef void DERIVATIVE_TYPE;
		static const bool IsPOD = true;
		//	typedef YESNO<true> IsPOD;
		static bool const hasSpeMemMove = false;
		static const bool NeedsAddLink = false; // containers needs to update addresses in link registers
		typedef ExOp SAFETYPE;


        PRIMTYPE_CONSTRUCT_FUNCTION(short)

        static bool isZero(const short &a){return a == 0;}
        static bool isOne(const short &a){return a == 1;}
		static void zero(short & a){a = 0;}
		static void one(short & a){a = 1;}
		static void random(short & a ){ a = rand() ^ (rand() << 12);};

		static void minimum(short& a) {a=-32768;}
		static void maximum(short& a) {a=32767;}
        static void epsilon(short & a){a = 1;}
		static void delta(short & a){a = 1;}

        static double sign(const short& val){return(0.0f);}
		static void show(const short &val, FILE* out=stdout, int level=0){fprintf(out, "%i", val); if (level == 0) fprintf(out, "\n"); }
        static string type_tostring(const short & a){return string("short");}


		static void save(const short& what, FILE *f) { fwrite(&what,sizeof(short), 1,f);}
		static void load(short& what, FILE *f, unsigned int lenght=0) { fread(&what,sizeof(short), 1,f);}

		static void memmove(short& a ,short& o){a=o;}
		static void bitreverse(short& a){
			a = ((a << 8 ) & 0xFF00) | ((a >> 8)& 0x00FF);
			a = ((a << 4 ) & 0xF0F0) | ((a >> 4)& 0x0F0F);
			a = ((a << 2 ) & 0xCCCC) | ((a >> 2)& 0x3333);
			a = ((a << 1 ) & 0xAAAA) | ((a >> 1)& 0x5555);
		}
		static void bytereverse(short& a){
			a = ((a << 8 ) & 0xFF00) | ((a >> 8)& 0x00FF);
		}
		template<unsigned int SIZE> static unsigned int getHyperDirection(const short (&a)[SIZE], const short (&b)[SIZE]){ // direction for hyper-ordering of arrays
			unsigned short xorb,xorc;
			*(unsigned short*)(&xorb) = a[0] ^ b[0]; // xor and make unsigned
			unsigned int best =0;
			for(unsigned int cur=1;cur<SIZE;cur++) {
				*(unsigned short*)(&xorc) = a[cur] ^ b[cur]; // xor and make unsigned
				if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
			}
			return(best);
		}
	};


	template< >
	class ExCo<unsigned short>{
	public:
		typedef unsigned char MAGNITUDE_TYPE;
		static const bool IsPOD = true;
		typedef void DERIVATIVE_TYPE;
		//	typedef YESNO<true> IsPOD;

		static bool const hasSpeMemMove = false;
		static const bool NeedsAddLink = false; // containers needs to update addresses in link registers
		typedef ExOp SAFETYPE;

        PRIMTYPE_CONSTRUCT_FUNCTION(unsigned short)

        static bool isZero(const unsigned short &a){return a == 0;}
        static bool isOne(const unsigned short &a){return a == 1;}

		static void zero(unsigned short & a){a = 0;}
		static void one(unsigned short & a){a = 1;}
		static void random(unsigned short & a ){ a = rand() ^ (rand() << 12);};
        static void epsilon(unsigned short & a){a = 1;}
		static void delta(unsigned short & a){a = 1;}
		static void minimum(unsigned short& a) {a=0;}
		static void maximum(unsigned short& a) {a=65535;}

		static void show(const unsigned short &val, FILE* out=stdout, int level=0){fprintf(out, "%u", val); if (level == 0) fprintf(out, "\n"); }
        static string type_tostring(const unsigned short & a){return string("unsigned short");}

		static void save(const unsigned short& what, FILE *f) { fwrite(&what,sizeof(unsigned short), 1,f);}
		static void load( unsigned short& what, FILE *f, unsigned int lenght=0) { fread(&what,sizeof(unsigned short), 1,f);}
		static void memmove(unsigned short& a ,unsigned short& o){a=o;}

		static unsigned char upperbound_pow_of_2(const unsigned short &i){
			if (i <= 256){
				if (i <= 16){
					if (i <= 4){
						if (i > 2) return(2);
						else if (i < 2) return(0);
						else return(0);
					}else return( (i <=8) ? 3 : 4 );
				}else return((i <= 64) ? ((i <=32) ? 5 : 6) : ((i <=128) ? 7 : 8 ));
			}else return( (i<= 0x1000) ? ((i <= 0x0400) ? ((i <=0x200) ? 9 : 10) : ((i <=0x0800) ? 11 : 12 )) : ((i <= 0x4000) ? ((i <=0x2000) ? 13 : 14) : ((i <=0x8000) ? 15 : 16 )) );
		}
		static void bitreverse(unsigned short& a){
			a = ((a << 8 ) & 0xFF00) | ((a >> 8)& 0x00FF);
			a = ((a << 4 ) & 0xF0F0) | ((a >> 4)& 0x0F0F);
			a = ((a << 2 ) & 0xCCCC) | ((a >> 2)& 0x3333);
			a = ((a << 1 ) & 0xAAAA) | ((a >> 1)& 0x5555);
		}
		static void bytereverse(unsigned short& a){
			a = ((a << 8 ) & 0xFF00) | ((a >> 8)& 0x00FF);
		}
		template<unsigned int SIZE> static unsigned int getHyperDirection(const unsigned short (&a)[SIZE], const unsigned short (&b)[SIZE]){ // direction for hyper-ordering of arrays
			unsigned short xorb = a[0] ^ b[0]; // xor and make unsigned
			unsigned int best =0;
			for(unsigned int cur=1;cur<SIZE;cur++) {
				unsigned short xorc = a[cur] ^ b[cur]; // xor and make unsigned
				if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
			}
			return(best);
		}
	};

	template< >
	class ExCo<string>{
	public:
		typedef unsigned char MAGNITUDE_TYPE;
		static const bool IsPOD = false;
		static const bool NeedsAddLink = false;
		typedef ExCo<string> SAFETYPE;
		typedef string REAL_TYPE;
		typedef string COMPLEX_TYPE;
		typedef string DEF_TYPE;
		typedef string NEG_TYPE;
		typedef string TRJU_TYPE;
        typedef string LMUL_TYPE;
		typedef string INNER_TYPE;
        typedef string VECTOR_TYPE;
		typedef string OUTER_TYPE; // outer product result, (increase in the number of dimentions!)
		
		typedef string QUATERNION_TYPE;
		typedef string GAUS_TYPE;
		
        typedef YESNO<false> IS_COMMUTATIVE;
        typedef YESNO<false> IS_OWNED;
        typedef YESNO<false> IS_POD;
		
		static string& memmove(string & tar, string & scr){tar.swap(scr); scr.clear(); return tar;}
		
		static void show(const string& what, FILE *f,int level) { if (level == 0) fprintf(f,"%s\n",what.c_str()); else fprintf(f,"%s",what.c_str());}
		
		static void save(const string& what, FILE *f) {unsigned int s = what.length(); fwrite(what.c_str(),sizeof(char), s+1,f); printf("da super write!\n");}
		static void load( string& what, FILE *f, unsigned int lenght=0) { char buffer[65536]; char* tmp = buffer; do {fread(tmp,sizeof(unsigned char), 1,f); } while(*(tmp++) != '\0'); what = string(buffer);}
	};

	/*
	 template< >
	 class ExCo<char*>{
	 public:
	 typedef YESNO<false> IsPOD;


	 static void save(const char*& what, FILE *f) { fwrite(&what,sizeof(int), 1,f);}
	 static void load(char*& what, FILE *f,) { fread(&what,sizeof(int), 1,f);}
	 };*/

	/*
	 template<class C, int order>
	 class ExCo< WeightElem<C, order> >{
	 public:
	 static WeightElem<C, order> zero(){
	 int i;
	 WeightElem<C, order> out;

	 for(i=0;i<order;i++) {out.w[i] = 0.0f; out.e[i] = ExCo<C>::zero();}
	 return(out);
	 }
	 static WeightElem<C, order> random(){
	 int i;
	 WeightElem<C, order> out;
	 for(i=0;i<order;i++) out.w[i] = 1.0f;
	 out.e[0] = ExCo<C>::rand();
	 for(i=1;i<order;i++) out.e[i] = out.e[i-1] * out.e[0];
	 return(out);
	 }

	 static void save(const  WeightElem<C, order>& what, FILE *f) { fwrite(&what.w,sizeof(double), order,f);}
	 static void load( WeightElem<C, order>& what, FILE *f,) { fread(&what.w,sizeof(double), order,f);}
	 };

	 template<class C, int order, Tuple_flag tflag>
	 class ExCo<Tuple<C, order, tflag> >{
	 public:

	 typedef YESNO< ExCo<C>::IsPOD::ans > IsPOD;

	 static Tuple<C, order, tflag> zero(){
	 int i;
	 Tuple<C, order, tflag> f_out;
	 for(i=0;i<order;i++) f_out[i] = ExCo<C>::zero();
	 return(f_out);
	 }
	 static Tuple<C, order, tflag> rand(){
	 int i;
	 Tuple<C, order, tflag> f_out;
	 for(i=0;i<order;i++) f_out[i] = ExCo<C>::rand();
	 return(f_out);
	 }

	 static double isValid(const Tuple<C, order, tflag> & a){return(a.isValid());}

	 static void show(const Tuple<C, order, tflag> &val, FILE* out){
	 val.show(out);
	 }

	 static void save(const  Tuple<C, order, tflag>& what, FILE *f) {
	 for(int i=0;i<order;i++) ExOp::save(what[i],f);
	 }
	 static void load( Tuple<C, order, tflag>& what, FILE *f) {
	 for(int i=0;i<order;i++) ExOp::load(what[i],f);
	 }
	 };

	 template<class C, int size, LFHVectorflag flag>
	 class ExCo< DataGrid<C,size> >{
	 public:
	 typedef YESNO< false > IsPOD;

	 static void zero(DataGrid<C,size> & targ){
	 int i = ((Vector<C>&)targ).size();
	 for(i--;i>=0;i--) ((Vector<C>&)targ).darray[i] = ExCo<C>::zero();
	 }
	 static void random( DataGrid<C,size> & targ){
	 int i = ((Vector<C>&)targ).size();
	 for(i--;i>=0;i--) ((Vector<C>&)targ).darray[i] = ExCo<C>::random();
	 }

	 static void hasTransferFunction(){return(true);} // has function which move ressources, and leaves the origin ready for destruction

	 };

	 template<class C>
	 class ExCo< Complex<C> >{
	 public:
	 typedef YESNO< ExCo<C>::IsPOD::ans > IsPOD;

	 static Complex<C> zero(){
	 Complex<C> out;
	 out.data[0] = ExCo<C>::zero();
	 out.data[1] = ExCo<C>::zero();
	 return(out);
	 }
	 static C norm(const Complex<C>& what ){ return(what.norm());}
	 static double sign(const Complex<C>& val){
	 return(cos(ExCo<C>::arctan(val.data[1],val.data[0])));
	 }




	 };*/






// } // end of namespace
