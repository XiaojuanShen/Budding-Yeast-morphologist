/*
 * Madstructs.h
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

#define REGRESS_POWER (8)
#define IMAGE_PREC unsigned char

#define MY_MIN(a,b,c,d) ( (c = a) < (d = b) ? c : d)
#define MY_MAX(a,b,c,d) ( (c = a) > (d = b) ? c : d)
#define MY_ERF_P(a) ( (a < 4.405459173f) ? 1.0f / 4294967296.0f : ((a > 10.0f)? 4294967295.0f / 4294967296.0f : (1.0f + erf(a))/2.0f) )
#define MY_LOGISTIC_P(a) ( (fabs(a) < 9.632959861146281f) ? 1 / (1 + exp(-a)) : ((a< 0) ?  1.0f / 4294967296.0f : 4294967295.0f / 4294967296.0f) )

#define MY_LOGOF_LOGISTIC_P(a) ( (fabs(a) < 9.632959861146281f) ? -log(1.0f + exp(-a)) : ((a< 0) ?  a : exp(-a) ) )

#define MY_LOGISTIC_DP(a) ( (fabs(a) < 20.0f) ? 0.5f / (1+ cosh(a)) : 0.0f)

//#define MY_LOGOFLOGISTIC_P(a) ( (fabs(a) < 9.632959861146281f) ? 1 / (1 + exp(-a)) : ((a< 0) ?  1.0f / 4294967296.0f : 4294967295.0f / 4294967296.0f) )

#define MY_SQRT_LOGISTIC(a) ((1 + (a)/sqrt((a)*(a) +4.0f))/2.0f)
#define MY_SQRT_LOGISTIC_INV(a) ((1 + (a)/sqrt((a)*(a) +4.0f))/2.0f)

#define MY_SHOW(a) a

#define TEMPLATE_IMAGE_TYPE unsigned char

using namespace std;
#include "primitive.h"
using namespace LFHPrimitive;

enum cellstuffenum{
	CELLSTUFF_NULL=0,
	CELLSTUFF_FOREBACK_FRACTION=1, // Tuple<double,2>
	CELLSTUFF_TAGSTRING=2, // Tuple<char,32>
	
	CELLSTUFF_GRBACK1DIST=240, // GaussianDistribution<2>
	CELLSTUFF_GRBACK2DIST=241, // GaussianDistribution<2>
	CELLSTUFF_GRBACK3DIST=242, // GaussianDistribution<2>
	CELLSTUFF_GRBACK4DIST=243, // GaussianDistribution<2>
	CELLSTUFF_GRFORE1DIST=244, // GaussianDistribution<2>
	CELLSTUFF_GRFORE2DIST=245, // GaussianDistribution<2>
	CELLSTUFF_GRFORE3DIST=246, // GaussianDistribution<2>
	CELLSTUFF_GRFORE4DIST=247, // GaussianDistribution<2>
	CELLSTUFF_GRTRANS1=248, //HMM<2>
	CELLSTUFF_GRTRANS2=249, //HMM<2>
	CELLSTUFF_GRTRANS3=250, //HMM<2>
	CELLSTUFF_GRTRANS4=251, //HMM<2>
	CELLSTUFF_GRFRAC1=252, //Tuple<double,4> // nb back, nb front, nb outlier,  unknown density
	CELLSTUFF_GRFRAC2=253,
	CELLSTUFF_GRFRAC3=254, 
	CELLSTUFF_GRFRAC4=255,
	
	CELLSTUFF_BACK1DIST=256, // GaussianDistribution<2>
	CELLSTUFF_BACK2DIST=257, // GaussianDistribution<2>
	CELLSTUFF_BACK3DIST=258, // GaussianDistribution<2>
	CELLSTUFF_BACK4DIST=259, // GaussianDistribution<2>
	CELLSTUFF_FORE1DIST=260, // GaussianDistribution<2>
	CELLSTUFF_FORE2DIST=261, // GaussianDistribution<2>
	CELLSTUFF_FORE3DIST=262, // GaussianDistribution<2>
	CELLSTUFF_FORE4DIST=263, // GaussianDistribution<2>
	CELLSTUFF_TRANS1=264, //HMM<2>
	CELLSTUFF_TRANS2=265, //HMM<2>
	CELLSTUFF_TRANS3=266, //HMM<2>
	CELLSTUFF_TRANS4=267, //HMM<2>
	CELLSTUFF_FRAC1=268, //Tuple<double,4> // nb back, nb front, nb outlier,  unknown density
	CELLSTUFF_FRAC2=269,
	CELLSTUFF_FRAC3=270, 
	CELLSTUFF_FRAC4=271,
	
	CELLSTUFF_GRBACK1DISTb=340, // GaussianDistribution<2>
	CELLSTUFF_GRBACK2DISTb=341, // GaussianDistribution<2>
	CELLSTUFF_GRBACK3DISTb=342, // GaussianDistribution<2>
	CELLSTUFF_GRBACK4DISTb=343, // GaussianDistribution<2>
	CELLSTUFF_GRFORE1DISTb=344, // GaussianDistribution<2>
	CELLSTUFF_GRFORE2DISTb=345, // GaussianDistribution<2>
	CELLSTUFF_GRFORE3DISTb=346, // GaussianDistribution<2>
	CELLSTUFF_GRFORE4DISTb=347, // GaussianDistribution<2>
	CELLSTUFF_GRTRANS1b=348, //HMM<2>
	CELLSTUFF_GRTRANS2b=349, //HMM<2>
	CELLSTUFF_GRTRANS3b=350, //HMM<2>
	CELLSTUFF_GRTRANS4b=351, //HMM<2>
	CELLSTUFF_GRFRAC1b=352, //Tuple<double,4> // nb back, nb front, nb outlier,  unknown density
	CELLSTUFF_GRFRAC2b=353,
	CELLSTUFF_GRFRAC3b=354, 
	CELLSTUFF_GRFRAC4b=355,
	
	CELLSTUFF_BACK1DISTb=356, // GaussianDistribution<2>
	CELLSTUFF_BACK2DISTb=357, // GaussianDistribution<2>
	CELLSTUFF_BACK3DISTb=358, // GaussianDistribution<2>
	CELLSTUFF_BACK4DISTb=359, // GaussianDistribution<2>
	CELLSTUFF_FORE1DISTb=360, // GaussianDistribution<2>
	CELLSTUFF_FORE2DISTb=361, // GaussianDistribution<2>
	CELLSTUFF_FORE3DISTb=362, // GaussianDistribution<2>
	CELLSTUFF_FORE4DISTb=363, // GaussianDistribution<2>
	CELLSTUFF_TRANS1b=364, //HMM<2>
	CELLSTUFF_TRANS2b=365, //HMM<2>
	CELLSTUFF_TRANS3b=366, //HMM<2>
	CELLSTUFF_TRANS4b=367, //HMM<2>
	CELLSTUFF_FRAC1b=368, //Tuple<double,4> // nb back, nb front, nb outlier,  unknown density
	CELLSTUFF_FRAC2b=369,
	CELLSTUFF_FRAC3b=370, 
	CELLSTUFF_FRAC4b=371,
	
};


namespace Madstructs {
	

	template<class T> class Image;
	template<class T> class ImageArray;
	class WImage;
	template<int islog> class Classifier;

	
	class Evaluatable;
	class Transformed;
	class ConcatOperator;
	class SubPartOperator;
	class GenericCellLoader;
	
	class TmpScope;class TmpScope_V2;

	class CellPose;
	class CellCover;
	class MultiCover;
	
	class SuperString;
	class StringList;

	class Table;
	
	class CellStageKeyFrame;

	
	union Anything;
	template<class T, int SIZE> class Chunk;
	
	class CellRecord;
	class CellRecord_Twomodes;
	class CellRecord_Twomodes2;
	class CellRecord_Header;
	class CellRecord_Header_Adapted;
	class CellRecord_Header_Reborn;
	class CellRecord_Header_Reborn_AGAIN;
	char* cloneString(char*);
	double lngamma(double);
	double polygamma0(double);
	void memecure(char* sequence);
	
	void solveCubic(double* coef, double *zeros); // 4 double, 6 double
	void solve3linear(double* m, double* outval); //12 double, 3 double
	void solvelinear(int size, double* matrix, double* outval);
	void solvegeteigen(int size, double* matrix, double* outval);
	double determinant(double* matrix, int size, int esize = 0);
	double* makeinverse(double* matrix, int size);
	double L2norm(double x, double y);
	
	char bitwiseswap(char r);
	short bitwiseswap(short r);
	int bitwiseswap(int r);
	void fastaout(FILE* out, char* sequence, int chunksize);
	
	double lncumhypergeoinv(int i, int a, int b, int n);
	double lncumhypergeo(int i, int a, int b, int n);
	
	int bimodalAnalaticalRoutine(double* rawmomments, double* in_std, double* out_means, double &out_mix);
	void bimodalsolverold(double* rawmomments, double* out_means, double* out_std, double &out_mix);
	void bimodalsolversimple(double* rawmomments, double* out_means, double* out_std, double &out_mix);
	
	void bimodalsolver(double* rawmomments, double &mean, double &var, double &mix, double &dif_mean, double &dif_var);
	
	void bimodalsolverGamma(double* rawmomments, double *mean, double *var, double &mix);
	
	void bimodalmommentsEMGamma(double* rawmomments, double sqrtshift, double *mean, double *var, double &mix);
	
	void polysolver(double* polynonimal, int order, LFHPrimitive::complex *outval);
	
	double polyeval(double* polynomial, int order, double x);
	void polyrootdivide(double* polynomial, int order, double root, double* out_poly);
	void polyderivative(double* polynomial, int order, double* out_poly);
	void polyprintf(double* polynomial, int order);
	double polybisection(double* polynomial, int order, double a, double b);
	void polysolverreal(double* polynonimal, int order, vector<double> &out);
	double sampleGaussian();
	double sampleGamma(double k, double t);
	
	void HiddenmapMake(vector<Madstructs::Image<float>* > & hiddenmap, vector<Madstructs::Image<unsigned short>* > &ima, vector<Madstructs::Image<float>* > &seg,  bool EMUP);
	
class Evaluatable{
public:
	virtual void eval(double* coor,double* output)=0;
	virtual void evalDerivatives(double* coor,double* output)=0;
	virtual int inputsize()=0;
	virtual int outputsize()=0;
};

	class GenericCellLoader{
	public:
		
		// Arguments
		bool no_daugther,no_clumped;
		double min_size;
		double max_size;
		Vector<const char *> collumn_names;
		Vector<double> factor;
		Vector<double> minus;
		Vector<unsigned int> take_log;
		Vector<unsigned int> filter_for_lone;
		
		myHashmap<string, unsigned int> protein_class_list;
		Vector<unsigned int> protein_class;
		Vector< string > protein_classes_names;
		
		bool include_area_in_data;
		float no_confidence;
		bool no_daughter_label; // they become 'c'
		bool no_lone_label; // they become 'c'
		mutable ProgressBarPrint daprogbar; 
		// Output
		Vector<KeyElem< Tuple<double, 4 > , Tuple<double, 0u > > > mother_buds;
		Vector<KeyElem< Tuple<double, 3 > , Tuple<double, 0u > > > daughters;
		Vector<KeyElem< Tuple<double, 3 > , Tuple<double, 0u > > > clumped;
		Vector<KeyElem< Tuple<double, 3 > , Tuple<double, 0u > > > loners;
		
		Vector<string> prot_orf;
		Vector<string> prot_name;
		
		Tuple<Trianglix<double> > bandwidth;
		Vector<double> bandwidth_scale;
		
		Vector<Tuple<unsigned int, 8u> > offsets;
		
		Vector<Tuple<GaussElem<Tuple<double> ,0u> > > cell_stage_bins;
		Vector<KeyElem<string, unsigned int> > bin_names;
		
		GenericCellLoader();
		void loadProteinClassLabel(const char* proteinclass_filepath);
		void loadAllCells(const char* pathlist, const char* prefix, const char* suffix);
		void loadAllCells_into_bins(const char* pathlist, const char* prefix, const char* suffix); // uses far less memory
		
		void loadAllCells_into_bins_tableread(const char* pathlist, const char* prefix, const char* suffix); // uses far less memory

		void reportN_in_bins() const;
		
		void MLGP() const;
		void MLClustering(const char* path) const;
		
	#ifdef GNU_SCIENTIFIC_LIBRARY
		Vector< KeyElem<double, Tuple<unsigned int, 2u> > >  LRTestAll(int which_test =0, int minimum_nbcell =1, bool MB_label = true, bool C_label = false,bool D_label = false,bool L_label = false) const; // NO LOESS, Confidence Weighting
		
		
		Vector< KeyElem<double, Tuple<unsigned int, 2u> > >  LRTestAllpairs(int which_test =0, int minimum_nbcell =1, bool MB_label = true, bool C_label = false,bool D_label = false,bool L_label = false) const; // NO LOESS, Confidence Weighting
		Vector< KeyElem<double, Tuple<unsigned int, 2u> > >  HottelingTestAllpairs(double confidence_thr,int minimum_nbcell =1, bool MB_label = true, bool C_label = false,bool D_label = false,bool L_label = false) const; // NO LOESS, Confidence Weighting
		
		Vector< KeyElem<double, Tuple<unsigned int, 2u> > >  LRTestAllpairs_frombins(int which_test =0, int minimum_nbcell =1, const Tuple<bool> *filter = NULL) const; // NO LOESS, Confidence Weighting
		Vector< KeyElem<double, Tuple<unsigned int, 2u> > >  HottelingTestAllpairs_frombins(int minimum_nbcell =1, const Tuple<bool> *filter =NULL) const; // NO LOESS, Confidence Weighting
		Vector< Tuple<double, 8u> >  LRTestAll_frombins(int minimum_nbcell, const Tuple<bool> *bin_filter =NULL, const Tuple<bool> *channel_filter =NULL, Vector<Tuple< Tuple< double > > > *local = NULL) const; // NO LOESS, Confidence Weighting

	//	Vector< Tuple<double, 8u> >  LRTestAll_frombins_median(int minimum_nbcell, const Tuple<bool> *bin_filter =NULL, const Tuple<bool> *channel_filter =NULL, Vector<Tuple< Tuple< double > > > *local = NULL) const; // NO LOESS, Confidence Weighting

		void find_bin_bandwidth_lowess(double fraction, bool class_partitionned = false);
	//	void find_bin_bandwidth_undiscrim(double signif_threshold);
#endif
		
		Vector<Tuple< Tuple< double > > > getMean() const;
		
		void two_ended_variance_analysis() const;

		void printCollumnName(FILE * f) const;
	};

	double cleveratof(const char* buf);
	double keyframeCorr(const Madstructs::CellStageKeyFrame &a, const Madstructs::CellStageKeyFrame &b);
	double keyframeNois(const Madstructs::CellStageKeyFrame &a);

	extern unsigned short font[];
	



class Transformed : public Evaluatable{
public:
	Madstructs::Evaluatable* ev; // not owned
	
	double* input_matrix;
	double* output_matrix;
	bool intrivial;
	bool outtrivial;
	
	int insize;
	int outsize;
	

	Transformed(Madstructs::Evaluatable* target);
	Transformed(int ins, int outs);
	~Transformed();
	void eval(double* coor,double* output);
	void evalDerivatives(double* coor,double* output);
	int inputsize();
	int outputsize();
	
	void invert(); // makes input be output and viseversa
	
	static Madstructs::Transformed* twopointsTransform(int insize, int outsize, double* sour, double* sink); // coordinates of two (distinc) point in the input and output, which defines a mapping
	void rescale(int sizex, int sizey, double mean); // works if insize ==2
	void renormalize(int sizex, int sizey, double mean, double std); // works if insize ==2
};



class ConcatOperator : public Madstructs::Evaluatable{
public:
	Madstructs::Evaluatable* para; // not owned
	Madstructs::Evaluatable* parb; // not owned
	ConcatOperator(Madstructs::Evaluatable* _para = NULL,Madstructs::Evaluatable* _parb = NULL);
	void eval(double* coor,double* output);
	void evalDerivatives(double* coor,double* output);
	int inputsize();
	int outputsize();
};

class SubPartOperator : public Madstructs::Evaluatable{
public:
	Madstructs::Evaluatable* para; // not owned
	int start;
	int nbchannels;
	SubPartOperator(int start =0, int nbchannels =1, Madstructs::Evaluatable* what = NULL);
	void eval(double* coor,double* output);
	void evalDerivatives(double* coor,double* output);
	int inputsize();
	int outputsize();
};



// flag 1 bit for log transform, 2 bit for distance weight channel
template<int flag =0>
class Classifier{
public:
	int nbstates;
	int inputsize;
	double* transit;
	double* boundspr;
	double* counts;
	double* means;
	double* stds;
	Classifier(int,int);
	template<class T> Madstructs::Image<float>* initializeOutput(vector<Madstructs::Image<T>*> &input);
	void applyConstraints();
	template<class T> void applyClassifier(vector<Madstructs::Image<T>*> &input,Madstructs::Image<float>*);
	template<class T> void updateClassifier(vector<Madstructs::Image<T>*> &input,Madstructs::Image<float>*, Madstructs::Image<float>* filter = NULL);
	template<class T> void applyClassifier(vector<Madstructs::Image<T>*> &input,Madstructs::Image<float>*, int rect[4]);
	template<class T> void updateClassifier(vector<Madstructs::Image<T>*> &input,Madstructs::Image<float>*, Madstructs::Image<float>* filter, int rect[4], int orifilter);
	template<class T> void updateClassifier(vector<Madstructs::Image<T>*> &input,Madstructs::Image<float>*, Madstructs::Image<float>* filter, int rect[4], double* expect_transit);
	void showparam();
	double getPriorMean(int inp,int channel);
	double getPriorStd(int inp,int state);


	void permState(int a, int b);
};


class CellPose{
public:
	double center[2];
	double eccentric[2];// eccentricity vector
	double width;
	double error[2];
	double area;
	
	bool isInside(int x, int y);
	//double distToEdge(int x, int y);
	double distsum(double x, double y);
	double cmpArea() const;
	
	void initFromCummul(double* data); // init from s,sx,sy,sxx,syy,syx
	
	void show(FILE*f = stdout,int level =0){fprintf(f, "center(%f,%f)\teccentric(%f,%f)\twidth(%f)\tarea(%f)\n", center[0], center[1], eccentric[0], eccentric[1], width, area);}
};


class CellCover{
public:
//	int plate;
	int rect[4];
	vector<Madstructs::CellPose> cover;
	bool findUnclassifCellpt(Madstructs::Image<double>* source, Madstructs::Image<double>* label ,double* out_coor, double* out_width);
	bool findUnclassifCellpt_2(Madstructs::Image<double>* source, Madstructs::Image<double>* label ,double* out_coor, double* out_width);
	
	void cellLabel(Madstructs::Image<double>* source, Madstructs::Image<double>* label);
	void cellDistanceLabel(Madstructs::Image<double>* source, Madstructs::Image<double>* label);
	
	void cellCompleteLabel(Madstructs::Image<double>* source, Madstructs::Image<double>* label, LFHPrimitive::GaussianDistribution<1>*);
	
	void computeerror(Madstructs::Image<double>* source);
	
	void cellUpdate(Madstructs::Image<double>* source, Madstructs::Image<double>* target, double alpha);
	void cellUpdate2(Madstructs::Image<double>* source, Madstructs::Image<double>* target,double expc, double alpha, double constrex);
	void cellUpdate3(Madstructs::Image<double>* source, Madstructs::Image<double>* target, double alpha, bool issearch = true);
	void cellUpdate_EM(Madstructs::Image<double>* source, Madstructs::Image<double>* target, double alpha, bool issearch = true);
	void cellUpdate4(Madstructs::Image<double>* source, Madstructs::Image<double>* target, double alpha);
	
	void cellUpdate6(Madstructs::Image<double>* source, Madstructs::Image<double>* target, double alpha);
	
	void cellIntUpdate(Madstructs::Image<double>* source, Madstructs::Image<double>* target, LFHPrimitive::GaussianDistribution<1>* );
	
	void cellEllipseFit(Madstructs::Image<double>* source, Madstructs::Image<double>* target);
	void cellUpdateEllipse(Madstructs::Image<double>* source, Madstructs::Image<double>* target, bool ellipseMode);
	void cellUpdateEllipse_EM(Madstructs::Image<double>* source, Madstructs::Image<double>* target, bool ellipseMode);
	
	void cellGenPrattFit(Madstructs::Image<double>* source, Madstructs::Image<double>* target, bool lowconf);
	
	
	void cellUpdateEllipse_EM_steps(Madstructs::Image<double>* source, Madstructs::Image<double>* target, bool ellipseMode);
	
	void addframe(vector<Madstructs::Image<unsigned char>*> &_out, Madstructs::Image<double> *label);
	
	void addframe_errormap(vector<Madstructs::Image<unsigned char>*> &_out, Madstructs::Image<double> *label);
	void addframe_distancemap(vector<Madstructs::Image<unsigned char>*> &_out, Madstructs::Image<double> *label);
	
	void cellError(Madstructs::Image<double>* source, Madstructs::Image<double>* error);
	
	LFHPrimitive::WeightElem<double, 2>* computeErrorunderHeuristic(const LFHPrimitive::DataGrid< LFHPrimitive::Tuple<double,1>, 2 > &source, LFHPrimitive::DataGrid< int ,2> incmessage) const;

//	LFHPrimitive::WeightElem<double, 2> computeErrorunderHeuristic(const LFHPrimitive::DataGrid< LFHPrimitive::Tuple<double,1>, 2 > &source, LFHPrimitive::Tuple<int,2> pt, LFHPrimitive::DataGrid< int ,2> incmessage, LFHPrimitive::DataGrid< double ,2> &wei) const;
	void findCellCoverHeuristic(Madstructs::Image<double>* source, double epsilon, LFHPrimitive::Classifier<LFHPrimitive::Tuple<double,1>,2>*);
	void findCellCoverHeuristic_Sep2011(Madstructs::Image<double>* source, double epsilon, LFHPrimitive::Classifier<LFHPrimitive::Tuple<double,1>,2>*);
	
	void findCellCover(Madstructs::Image<double>* source);
	void findCellCover_2(Madstructs::Image<double>* source, const LFHPrimitive::GaussianDistribution<1>& prbg, const LFHPrimitive::GaussianDistribution<1>& prcell,bool movie);
	void findCellCover_hiddenmap(Madstructs::Image<double>* source, const LFHPrimitive::GaussianDistribution<1>& prbg, const LFHPrimitive::GaussianDistribution<1>& prcell,bool movie , vector<Madstructs::Image<float>* > *PMHidden_Map, bool noem = false);
	void findCellCover_hiddenmap2(Madstructs::Image<double>* source, bool movie , vector<Madstructs::Image<float>* > *PMHidden_Map, bool noem = false);
	
	void findCellCoverAlgerabric(Madstructs::Image<double>* source, double contour_min, double contour_max);
	void findCellCoverAlgerabricHeuristic(Madstructs::Image<double>* source, double contour_min, double contour_max,double epsilon,  unsigned int sub_problems_radii, unsigned int hint_max_width =0);
	
	void findCellCover_hiddenmap_upgr(Madstructs::Image<double>* source, bool movie , vector<Madstructs::Image<float>* > *PMHidden_Map, bool noem = false); // (not working)
	void findCellCover_HMM(Madstructs::Image<double>* source); // (not working)
	
	
	void findCellCoverHeuristic_agglo(Madstructs::Image<double>* source); // (not working)
	void findCellCoverHeuristic_agglo2(Madstructs::Image<double>* source, double min_area); // (not working)
	
	void findCellCover_contour(Madstructs::Image<double>* source, double perim_dist);
	void update_from_contour(Madstructs::Image<double>* source, double perim_dist);

	
	
//	void findCellCover_distance(Madstructs::Image<double>* source);
	
	
//	double EllipseClusterLabel(LFHPrimitive::DataGrid<LFHPrimitive::Tuple<double,3>,2 >& source, LFHPrimitive::DataGrid<LFHPrimitive::Tuple<int,2>,2 >& label, LFHPrimitive::Tuple<int , 2> &worstpt);
//	void findCellCover(LFHPrimitive::DataGrid<LFHPrimitive::Tuple<double,3>,2 >& source, double minradii);
	
	
	
	
	
	void save(char* path);
	void load(char* path);
	
	bool isInside(int k);
	
	
	
	Madstructs::Image<double>* makeVectField();
	Madstructs::Image<double>* makeVectField(double zval);
	void fillVectField(Madstructs::Image<double> *vecf);
	
	Madstructs::Transformed* genTransformation(int modelsize);
	Madstructs::Transformed* genTransformation(int modelsize, int which);
	
	void getStage(double* _out, int* otherid);
	void getStageType(char* _out);

    Madstructs::Transformed* genTransposer();
	
	void showLinks(FILE* out);
	void showTripplets(FILE* out);
	void showLinksErdist(FILE* out);
	
	bool isValid();
};

class MultiCover{
public:
	vector<Madstructs::CellCover> group;
	
	
	
	void save(char* path);
	void load(char* path);
	
	void save(FILE* stre);
	void load(FILE* stre);
	
	template<class C> void drawDarkCircles(Madstructs::Image<C>* where);
	template<class C> void drawLabel(Madstructs::Image<C>* where);
	void drawCircles(Madstructs::Image<unsigned char>* where, bool ignore_conf);
	void evalCircles(Madstructs::Image<unsigned char>* where);
	
	template<class T> void drawCircles(Madstructs::Image<T>* where, int channel, double value, bool ignoreconf);
	template<class T> void drawACircle(Madstructs::Image<T>* where, int channel, double value, int i, int j, bool ignoreconf);
	template<class T> void drawCross(Madstructs::Image<T>* where, int channel, double value);
	
	void evalfit(FILE*f,Madstructs::Image<float>* segmented);
	void filterbyerror(double threshold);
	void filterbydist(double threshold);
	bool getNextValid(int &i);
	void flush();
	
	void show(FILE *where);
	void showtable(FILE *where);
};



class SuperString{
public:
	char buffer[256];
	vector<int> chunks;
	int startlength;
	SuperString(char* start);
	void setEnd(char* what);
	void setChar(int pos, char value);
	
	void push(char* string);
	void pop();
	char* getBuffer();
	void digitsWrite(int value, int start , int prec);
	
};

class StringList{
public:
	char buffer[256];
	FILE *f;
	StringList(char* path);
	~StringList();
	bool getNext();
	void reset();
};
 
union Anything { // BEWARE! does not delete data pointer to
	void* p;
	double d;
	int i;
	short s;
	char sc;
	char* c;
	float f;
	Anything();
	Anything(char* p);
	Anything(float p);
	Anything(double p);
	Anything(int p);
};

class Table{
public:
	int nbrows;
	int nbcols;
	vector<Madstructs::Anything> t;
	char *rowformat;
	
	Table(): rowformat(NULL){}
	Table(int _nbcols);
	Table(const char* filepath); // tabular format assumed
	Table(const char* filepath, const char* predicted_format); // no format assumed
	~Table(){clear();}
	void wiseload(const char* filepath);
	Madstructs::Anything getValue(int x, int y);
	void defineCol(int colId,char dformat[2], Madstructs::Anything label);
	void newRow();	
	void setValue(int x, int y, Madstructs::Anything what);
	
	char* colname(int col)const;
	void writecoltype(char *target, int col){target[0] =rowformat[col*2];target[1] =rowformat[col*2+1];target[2] ='\0'; }
	bool getValueFromTabular(char* tmp, FILE* f, char sep);
	
	Madstructs::Anything& operator()(int x, int y);
	
	int findValue(Madstructs::Anything val, int col) const;
	
	int findCol(const char* colname)const;
	void genHashMap(map<string, int> &_out, int col);
	
	void printHTML(char* format, vector<Madstructs::Anything>::iterator &i, FILE* f);
	void printTabular(char* format, vector<Madstructs::Anything>::iterator &i, FILE* f);
	void makeHTML(const char* path);
	void makeTabular(char* path);
//	void show(FILE* f) const;
	void clear();
};

class CellRecord{
public:
	char name[16];
	int groupId;
	int groupSize;
	int cellsubId;
	double cellstage;
	double intensity[6];
	double width;
	double error[2];
	double factor;
	unsigned short intensitycount[256];
};

class CellRecord_Twomodes{
public:
	char name[32];
	int groupId;
	int groupSize;
	int cellsubId;
	int otherId;
	double cellstage;
	double intensity[5]; // pixcount, mean, var ... 
	double HMMmatrix[4];
	double lowmode_intensity[5]; // pixcount, mean var
	double highmode_intensity[5]; // pixcount, mean var
	int didconverge;
	double width;
	double error[2];
	double factor;
	unsigned short intensitycount[256];
	void show(FILE* out);
};

class CellRecord_Twomodes2{
public:
	char name[32];
	int groupId;
	int groupSize;
	int cellsubId;
	int quad;

	double cellstage;
	double red_intensity[5]; // pixcount, mean, var ... 
	double intensity[5]; // pixcount, mean, var ... 
	double HMMmatrix[4];
	double lowmode_intensity[5]; // pixcount, mean var
	double highmode_intensity[5]; // pixcount, mean var
	int didconverge;
	double width;
	double error[2];
	double factor;
	unsigned short intensitycount[256];
	void show(FILE* out);
};

class CellRecord_Header{
public:
	char name[32];
	int groupId;
	int groupSize;
	int cellsubId;
	int otherId;
	int quad;
	double cellstage;
	int didconverge;
	double width;
	double error[2];
	double factor;
	unsigned short intensitycount[256];
	void show(FILE* out);
};

class CellRecord_Header_Adapted{
public:
	char name[32];
	int groupId;
	int groupSize;
	int cellsubId;
	int otherId;
	int quad;
	double cellstage;
	double width;
	double error[2];
	double area;
	void show(FILE* out);
};

class CellRecord_Header_Reborn{
public:
	int cellID;
	int otherID;
	int groupID;
	unsigned char groupSize;
	unsigned char groupSize_before;
	char subcellID;
	char subcellOtherID;
	double width;
	double raw_red_deviation;
	double error_red;
	double raw_error_dist;
	double error_dist;
	double error_ellipse_dist;
	double area;
	double delta_MD_distance;
	double center_x;
	double center_y;
	double excentric_x;
	double excentric_y;
	void show(FILE* out);
};

class CellRecord_Header_Reborn_AGAIN{
public:
	int frameID;
	unsigned char groupSize;
	int cellID;

	int otherID; 
	double other_dist; 
	double other_major; 
	double other_area; 

	double other_dist_err; 
	double other_density; 
	double other_red_dev; 
	double other_contour; 
	double other_ramanujan;

	
	char celltype;
	
	double raw_red_deviation;
	double error_red;
	double error_dist;
	double density_error;
	
	double area;
	double edgedist_raw_error;
	double edgedist_radius;
	double delta_MD_distance;
	double fold_red_err;
	double red_green_correlation;
	
	double cell_prob;
	double relation_error;
	double cell_error;
	double relation_prob;
	double cell_prob_factor;
	double relation_prob_factor;
	
	double contour; 
	double ramanujan;

	double mc_x;
	double mc_y;
	double mcdist;
	double budneck_x;
	double budneck_y;
	double neckdist;
	
	// ellipse parameters
	double width;
	double center_x;
	double center_y;
	double excentric_x;
	double excentric_y;
	
	double g_area;
	double g_width;
	double g_center_x;
	double g_center_y;
	double g_excentric_x;
	double g_excentric_y;
	
	void show(FILE* out);
};

class CellStageKeyFrame{
public:
	static const bool IsPOD = false;
	static const bool IsComplex = false;
	static const bool NeedsAddLink = false;
	
	char type;
	double MD_stage;
	double size_stage;
	double MDprob; // conditionnal probability on the fact both cells are not artefacts
	double CellProb;
	double noise;
	
	void show(FILE* f, unsigned int level) const;
	double getWeight() const{return CellProb;}
	double getNorm(const Madstructs::CellStageKeyFrame &other) const;
};



class TmpScope{
public:
	static const bool IsPOD = true;
	static const bool IsComplex = false;
	static const bool NeedsAddLink = false; // containers needs to update addresses in link registers
	
	unsigned int frameID;

	unsigned int backID;
	unsigned int otherID;
	double fold_red_err;
	
	LFHPrimitive::WeightElem<double, 2> distance_error;

	LFHPrimitive::WeightElem<double, 4> green;
	LFHPrimitive::WeightElem<double, 4> red;
	LFHPrimitive::WeightElem<double, 1> red_green;
	LFHPrimitive::WeightElem<double, 1> red_intensity;

	LFHPrimitive::WeightElem<double, 4> intensity;

	LFHPrimitive::WeightElem<double, 4> bn_dist;
	LFHPrimitive::WeightElem<double, 4> ed_dist;
	LFHPrimitive::WeightElem<double, 4> sd_dist;
	LFHPrimitive::WeightElem<double, 4> cm_dist;
	LFHPrimitive::WeightElem<double, 4> cn_dist;

	
	LFHPrimitive::WeightElem<double, 4> r_bn_dist;
	LFHPrimitive::WeightElem<double, 4> r_ed_dist;
	LFHPrimitive::WeightElem<double, 4> r_sd_dist;
	LFHPrimitive::WeightElem<double, 4> r_cn_dist;
	LFHPrimitive::WeightElem<double, 4> r_gcm_dist;
	LFHPrimitive::WeightElem<double, 4> r_cm_dist;
	LFHPrimitive::WeightElem<double, 4> rg_sd_dist;
	
	

	double cell_prob;
	double cell_prob_factor;
	
	unsigned int contour;
	double ramanujan;
	double density_error;
	double dist_error;
	
	unsigned int neighbor[32]; 
	unsigned int neighbor_size;
	
	char type;
	double ellipic[6];
	double neck[2];
	double centerofmass[5];
	TmpScope(): neighbor_size(0){memset(ellipic,'\0',sizeof(double) *6); memset(centerofmass,'\0',sizeof(double) *5); LFHPrimitive::ExOp::toZero(red); LFHPrimitive::ExOp::toZero(green); 
	LFHPrimitive::ExOp::toZero( bn_dist); LFHPrimitive::ExOp::toZero(ed_dist); LFHPrimitive::ExOp::toZero(sd_dist);LFHPrimitive::ExOp::toZero(cm_dist);LFHPrimitive::ExOp::toZero(cn_dist);
	LFHPrimitive::ExOp::toZero( r_bn_dist); LFHPrimitive::ExOp::toZero(r_ed_dist); LFHPrimitive::ExOp::toZero(r_sd_dist);LFHPrimitive::ExOp::toZero(r_cm_dist);	LFHPrimitive::ExOp::toZero(r_cn_dist); LFHPrimitive::ExOp::toZero(r_gcm_dist);
		LFHPrimitive::ExOp::toZero(rg_sd_dist);
		LFHPrimitive::ExOp::toZero(red_green);
		LFHPrimitive::ExOp::toZero(red_intensity);
		LFHPrimitive::ExOp::toZero(intensity);
	}
};

class TmpScope_V2{
public:
	static const bool IsPOD = true;
	static const bool IsComplex = false;
	static const bool NeedsAddLink = false; // containers needs to update addresses in link registers
	
	unsigned int frameID;
	unsigned int backID;
	unsigned int otherID;
	double fold_red_err;
	
	double cell_prob;
	double cell_prob_factor;
	unsigned int contour;
	double ramanujan;
	double density_error;
	double dist_error;
	unsigned int neighbor[32]; 
	unsigned int neighbor_size;
	char type;
	double ellipic[6];
	double neck[2];
	double centerofmass[6];
	LFHPrimitive::WeightElem<double, 2> distance_error;

	LFHPrimitive::WeightElem<double, 4> green;
	LFHPrimitive::WeightElem<double, 4> red;
	
	double crossdist[10];
	
	LFHPrimitive::WeightElem<double, 2> red_ell_dist[15];
	LFHPrimitive::WeightElem<double, 2> green_ell_dist[15];
	LFHPrimitive::WeightElem<double, 2> area_ell_dist[15];
	
	TmpScope_V2(): neighbor_size(0){LFHPrimitive::ExOp::toZero(ellipic); LFHPrimitive::ExOp::toZero(centerofmass);LFHPrimitive::ExOp::toZero(distance_error);
		LFHPrimitive::ExOp::toZero(red); LFHPrimitive::ExOp::toZero(green); 
		LFHPrimitive::ExOp::toZero(red_ell_dist);
		LFHPrimitive::ExOp::toZero(green_ell_dist);
		LFHPrimitive::ExOp::toZero(area_ell_dist);
	}
	
	
	
	
};


	template<class T>
	class Image : public Madstructs::Evaluatable{
	public:
		int sizex;
		int sizey;
		int bits;
		int channels;
		
		T* data;
		Image();
		~Image();
		
		template<class D, unsigned int NBCHAN> Image( LFHPrimitive::DataGrid<LFHPrimitive::Tuple<D, NBCHAN>,2> & other);
		template<class D> Image( LFHPrimitive::DataGrid<D,2> & other);
		
		static void LoadTiffImage(vector<Madstructs::Image<T>* > &out, const char* path);
		static void LoadTiffImageBin(vector<Madstructs::Image<unsigned char>* > &out, const char* path);
		static void SaveTiffImage(vector<Madstructs::Image<T>* > &out, const char* path, int ignore = 0, int nbchannels = 0);
		
		static Image<T>* makeImageFromFunction(int sizex, int sizey, Madstructs::Evaluatable* target);
		static Image<T>* makeImageFromFunction(int x, int y, int sizex, int sizey, Madstructs::Evaluatable* target);
		
		
		static void addWeigthedPixel(vector<Madstructs::Image<T>* > im, double x, double y, double z, double* pix, double weigth =1.0f); // 3d
		static void propagablur(vector<Madstructs::Image<T>* >);
		static void propagablurexp(vector<Madstructs::Image<T>* > &im);
		static void weightrescale(vector<Madstructs::Image<T>* > &im);
		void loadRawImage(char* path);
		void saveRawImage(char* path);
		
		template<class D, unsigned int NBCHAN> operator LFHPrimitive::DataGrid<LFHPrimitive::Tuple<D, NBCHAN>,2>();
		template<class D, unsigned int NBCHAN> LFHPrimitive::DataGrid<LFHPrimitive::Tuple<D, NBCHAN>,2> makedatagrid();
		
		//	template<class D> operator LFHPrimitive::DataGrid<D,3>();
		//	template<class D> LFHPrimitive::DataGrid<D,3> makedatagrid();
		
		Madstructs::Image<double>* makeMorphingInput();
		
		// Evaluatable functions
		void eval(double* coor,double* output);
		void evalDerivatives(double* coor,double* output);
		void rescale(int x, int y);
		
		
		int inputsize();
		int outputsize();
		// Evaluatable functions
		
		void allocateBuffer();
		
		
		template <class S> void allocateBuffer(Madstructs::Image<S>* Toclone);
		template <class S> void allocateBuffer(LFHPrimitive::DataGrid<S,2> &Toclone, int _nbchannels);
		void makeEmptyImage(int n_sizex, int n_sizey, int n_channels, int n_bits);
		void initBlack();
		void LoadBmpImage(char* path);
		void SaveBmpImage(char* path);
		
		
		void getPixel(int x, int y, double *outp);
		void setPixel(int x, int y, double *outp);
		void addPixel(int x, int y, double *outp, double factor =1);
		void scaleImage(double factor);
		void blur(double radius);
		void scaleweight(float fact);
		
		void getMaxPixel(double *out);
		void shiftColor(int x, int y, double factor);
		
		void makeImageFromChannels(vector<Madstructs::Image<T>*> &, vector<int>&);
		
		void makeImageFromChannels(int sizex, int sizey,vector<Madstructs::Evaluatable*> &, vector<int> &);
		
		template<class S> void addImage(Madstructs::Image<S>* );
		
		void computeColorDistribution(double center[2], double radius, double* meanvar, double& sampleweigth, double *filter);
		double computeWhiteness(double center[2], double radius);
		// finds the best coordinates which fits a grid of spots of the provided number of row and columns
		double findGrid_evalroutine(double *coords, int col, int row,double size);
		double findGrid_evalroutine2(double *coords, int col, int row,double size,double size2,double size3, double* filter);
		double findGrid_evalroutine3(double *coords, int col, int row,double size2);

		
		double* makecolorDistanceMap();
		double* makecolorLinearityMap();
		void markbymap(double* map, double threshold);
		void applyLinearityFilter(double radius);
		

		
		void fillHistogram(vector<int*> counts);
		
		template<class S> void copyFrom(Madstructs::Image<S>* other);
		
		void GradientFilter(Madstructs::Image<unsigned char>* out);
		void GradientFilter(Madstructs::Image<unsigned short>* out);
		void BFpropagation(double *transit, double *bounds, int dir, Madstructs::Image<float>* out);
		void BFpropagation(double *transit, double *bounds, int dir, Madstructs::Image<float>* out, int rect[4]);
		void BFpropagation(double *transit, double *bounds, int dir, Madstructs::Image<float>* out, int rect[4], double* updatetransit);
		
		//	void WavyBFpropagation(double *transit, double *bounds, int dir, Madstructs::Image<float>* out);
		
		void computeDistances(Madstructs::Image<float>* out, Madstructs::Image<float>* rprob);
		void groupPixel(Madstructs::Image<unsigned char>* origred,Madstructs::Image<unsigned char>* origreen, char* name);
		void groupCovPixel(Madstructs::Image<unsigned char>* origred,Madstructs::Image<unsigned char>* origreen, char* name);
		void groupCovPixel(Madstructs::Image<unsigned char>* origred,Madstructs::Image<unsigned char>* origreen, FILE* name);
		void groupCovPixel(Madstructs::Image<unsigned short>* origred,Madstructs::Image<unsigned short>* origreen, FILE* where);
		void groupCovPixel(Madstructs::Image<float>* origred, FILE* where);
		
		void drawLetter(char what, int x, int y,int c);
		void drawLine(int x, int y, int x2, int y2,int c);
		
		template<class S> void groupCovPixel(Madstructs::Image<S>* origred,Madstructs::Image<S>* origreen, int plateId, Madstructs::MultiCover &_out);
		
		void initMorphMap(Madstructs::Evaluatable* source, Madstructs::Evaluatable* sink);
		void morphingStepSym(Madstructs::Evaluatable* source, Madstructs::Evaluatable* sink);
		void morphingStep(Madstructs::Evaluatable* source, Madstructs::Evaluatable* sink, double magnitude);
		void probmorphingStep(Madstructs::Evaluatable* source, Madstructs::Evaluatable* sink, double magnitude);
		double evalMap(Madstructs::Evaluatable* source, Madstructs::Evaluatable* sink, double &var);
	//	void writeTransformedPixels(Madstructs::Evaluatable* source, Madstructs::Transformed* xform, Madstructs::ImageArray<double>* target, double weight =1.0f);
		void writeTransformedPixels(Madstructs::Evaluatable* source, Madstructs::Transformed* xform, Madstructs::Image<double>* target, double weight =1.0f);
		double computeArea(Madstructs::Evaluatable* source);
		void vectorfieldRescale(double factor);
		void genStdDevs(int start, int nbchannels);
		
		void pixelDistribution(int rect[4], double* _out);
		
		
		Image<double>* exponentialblur_x(double coef,bool hasweights);
		Image<double>* exponentialblur_y(double coef,bool hasweights);
		Image<double>* exponentialGradientNorm(double coef,bool hasweights);
		Image<double>* exponentialMultiGradientNorm(bool hasweights);
		
		//	Image<double>* maxApertureGradient(bool hasweights, int nbaperture, double min, double max);
		
		virtual void exponentialblur(double coef, bool hasweights = false);
		virtual void exponentialblur(double coef, unsigned int filter,bool hasweights = false);
		
		void weightrescale();
		void propagablur();
		
		template<class S> void generateMorphedImage(int size, Madstructs::Image<double>* map, Madstructs::Image<S>* &out);
		template<class S> void generateMorphedImage(int size, Madstructs::Image<double>* map, Madstructs::Image<S>* &out, int* rect);
	};
	
	template< > Madstructs::Image<unsigned char>* Image<unsigned char>::makeImageFromFunction(int sizex, int sizey, Madstructs::Evaluatable* T);
	template< > Madstructs::Image<unsigned char>* Image<unsigned char>::makeImageFromFunction(int x, int y,int sizex, int sizey, Madstructs::Evaluatable* T);
	template< >	void Image<unsigned char>::LoadTiffImage(vector<Madstructs::Image<unsigned char>* > &out, const char* path);
	template< >	void Image<unsigned short>::LoadTiffImage(vector<Madstructs::Image<unsigned short>* > &out, const char* path);
	template< >	void Image<unsigned char>::LoadTiffImageBin(vector<Madstructs::Image<unsigned char>* > &out, const char* path);
	template< >	void Image<unsigned char> ::SaveTiffImage(vector<Madstructs::Image<unsigned char>* > &list, const char* path, int ignore, int nbchannels);
	template< >	void Image<unsigned short> ::SaveTiffImage(vector<Madstructs::Image<unsigned short>* > &list, const char* path, int ignore, int nbchannels);
	
	template< >	void Image<float>::LoadTiffImage(vector<Madstructs::Image<float>* > &out, const char* path);
	template< >	void Image<float>::SaveTiffImage(vector<Madstructs::Image<float>* > &list, const char* path, int ignore, int nbchannels);
	template< >	void Image<double>::LoadTiffImage(vector<Madstructs::Image<double>* > &out, const char* path);
	template< >	void Image<double>::SaveTiffImage(vector<Madstructs::Image<double>* > &list, const char* path, int ignore, int nbchannels);
	
	template< >	void Image<unsigned char> ::makeImageFromChannels(vector<Madstructs::Image<unsigned char> *> &srclist, vector<int> &srcchan);
	template< >	void Image<unsigned short> ::makeImageFromChannels(vector<Madstructs::Image<unsigned short> *> &srclist, vector<int> &srcchan);
	template< >	void Image<unsigned int> ::addPixel(int x, int y, double *out, double factor);
	template< >	void Image<double> ::addPixel(int x, int y, double *out, double factor);
	template< >	void Image<unsigned char>::getPixel(int x, int y, double *outp);
	template< >	void Image<unsigned short>::getPixel(int x, int y, double *outp);
	template< >	void Image<unsigned int>::getPixel(int x, int y, double *outp);
	template< >	void Image<float>::getPixel(int x, int y, double *outp);
	template< >	void Image<double>::getPixel(int x, int y, double *outp);
	template< >	void Image<float>::setPixel(int x, int y, double *outp);
	template< >	void Image<double>::setPixel(int x, int y, double *outp);
	template< >	void Image<unsigned char>::setPixel(int x, int y, double *outp);
	template< >	void Image<unsigned short>::setPixel(int x, int y, double *outp);
	template< >	void Image<unsigned int>::setPixel(int x, int y, double *outp);
	template< >	void Image<unsigned char>::GradientFilter(Madstructs::Image<unsigned char>* out);
	template< >	void Image<unsigned short>::GradientFilter(Madstructs::Image<unsigned short>* out);
	template< >	void Image<float>::BFpropagation(double *transit, double *bounds, int dir, Madstructs::Image<float>* out);
	template< >	void Image<float>::BFpropagation(double *transit, double *bounds, int dir, Madstructs::Image<float>* out, int rect[4]);
	template< >	void Image<float>::BFpropagation(double *transit, double *bounds, int dir, Madstructs::Image<float>* out, int rect[4], double *outtransit);
	template< >	void Image<float>::computeDistances(Madstructs::Image<float>* out, Madstructs::Image<float>* rprob);
	template< >	void Image<unsigned char>::groupPixel(Madstructs::Image<unsigned char>* origred,Madstructs::Image<unsigned char>* origreen, char* name);
	template< >	void Image<float>::groupCovPixel(Madstructs::Image<unsigned char>* origred,Madstructs::Image<unsigned char>* origreen, char* name);
	template< >	void Image<float>::groupCovPixel(Madstructs::Image<unsigned char>* origred,Madstructs::Image<unsigned char>* origreen, FILE* name);
	template< >	void Image<float>::groupCovPixel(Madstructs::Image<unsigned short>* origred,Madstructs::Image<unsigned short>* origreen, FILE* name);
	template< >	void Image<float>::groupCovPixel(Madstructs::Image<float>* origreen, FILE* name);
	template< >	void Image<double>::initMorphMap(Madstructs::Evaluatable* source, Madstructs::Evaluatable* sink);
	template< > void Image<double>::morphingStepSym(Madstructs::Evaluatable* source, Madstructs::Evaluatable* sink);
	template< > void Image<double>::morphingStep(Madstructs::Evaluatable* source, Madstructs::Evaluatable* sink, double magnitude);
	template< > void Image<double>::probmorphingStep(Madstructs::Evaluatable* source, Madstructs::Evaluatable* sink, double magnitude);
	template< > double Image<double>::evalMap(Madstructs::Evaluatable* source, Madstructs::Evaluatable* sink, double &var);
	template< > void Image<double>::loadRawImage(char* path);
	template< > void Image<double>::saveRawImage(char* path);
	template< >	void Image<double>::genStdDevs(int start, int nbchannels);
	template< > void Image<double>::writeTransformedPixels(Madstructs::Evaluatable* source, Madstructs::Transformed* xfodrm, Madstructs::Image<double>* target, double weight);
//	template< > void Image<double>::writeTransformedPixels(Madstructs::Evaluatable* source, Madstructs::Transformed* xform, Madstructs::ImageArray<double>* target, double weight);
	template< >	double Image<double>::computeArea(Madstructs::Evaluatable* source);
	template< > void Image<double>::vectorfieldRescale(double factor);
	
	#include "Madstructs_t.h"


}




