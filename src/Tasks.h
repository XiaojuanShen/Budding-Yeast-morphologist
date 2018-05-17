/*
 * Tasks.h
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

using namespace std;
#include "./primitive.h"
using namespace LFHPrimitive;
#include "Madstructs.h"
enum task {
	TASK_NULL=0,
	
	// in TiffManip.cpp
	TASK_REMOVE_FACTORS, 
	TASK_TIFF_FILE_OPERATION, 
	TASK_TIFF_FILE_MANIPULATION,
	// in segmentation.cpp
	TASK_HMM_SEGMENTATION_AND_DISTANCE,
	TASK_GEN_META_PIXELS,
	TASK_FIND_CELL_COVER_UPGRADE, 
	TASK_EXTRACT_HIDDENMAP_DIRECT,  
	// in Extraction.cpp
	TASK_HIDDENMAP_BASE_DATAEXTRACTION,
	TASK_MAKE_DISPLAYABLE,
	TASK_GETLABEL_CONFIDENCE,
	// in Pvalue.cpp
	TASK_CDTGTR_PVALUES, // Pvalue.cpp
	TASK_RANK_PVALUES, // Pvalue.cpp
	TASK_MULTI_PAIR_DISTANCE_PVAL,
	// in Modelling.cpp 
	TASK_PROFILES,
	TASK_HIDDENMAP_BASE_DATAEXTRACTION_V2

};

	template<task TASKTYPE>
	class Taskscope{
	public:
		static int submain(int argc, char * const argv[]);
		static int submain_extra(int argc, char * const argv[]);
	};
	
#define TASK_MEMBER_DEFINITIONS Taskscope(); int defstore(char* const * token, int nbtoken); void store(char* const * token, int nbtoken); void nbaddtoken(char const * const token, int& min, int& max); void help();
	
	template< > class Taskscope<TASK_HMM_SEGMENTATION_AND_DISTANCE>: public LFHPrimitive::ArgumentParser{public:
		const char* file_in;
		const char* file_out;
		bool show;
		const char* file_average;
		const char* file_stddev;
		const char* file_scope;
		const char* file_dist;
		const char* file_nohmm;
		bool bright;
		unsigned int crop_rect[4];
		double min,max;
		double artfrac;
		
		double blur_wind;
		int nbclasses;
		TASK_MEMBER_DEFINITIONS
	};	
	
	template< > class Taskscope<TASK_FIND_CELL_COVER_UPGRADE>: public LFHPrimitive::ArgumentParser{public:
		const char* file_in;
		const char* file_in2;
		const char* file_in3;
		const char* out_mcv;
		const char* out_tif;
		const char* out_txt;
		
		const char* in_mcv_known;
		
		//	const char* out_prvw_tif;
		
		int hint_max_width;
		bool doem;
		int whichheu;
		bool quick;
		const char* scope_file;
		
		bool show;
		int max_dist;
		
		double contour[2];
		
		TASK_MEMBER_DEFINITIONS
	};	
	

	template< > class Taskscope<TASK_GEN_META_PIXELS> : public LFHPrimitive::ArgumentParser{public:
		const char* file_in;
		const char* file_out;
		const char* file_mcv;
		const char* mask_file;
		const char* out_tif;
		const char* out_dist;
		const char* out_partition;
		const char* out_cpartition;
		bool knight;
		bool distpen;
		//	const char* out_prvw_tif;
		bool mark_da_maxima;
		bool quick;
		bool show;
		bool comp;
		float valid_int_range[2];
		int min_nbpixel;
		const char* grametapix;
		double blur_std;
		TASK_MEMBER_DEFINITIONS
	};	
	
	template< > class Taskscope<TASK_REMOVE_FACTORS> : public LFHPrimitive::ArgumentParser{public:
		TASK_MEMBER_DEFINITIONS
	};	
	
	template< > class Taskscope<TASK_TIFF_FILE_MANIPULATION> : public LFHPrimitive::ArgumentParser{public:
		const char* file_in;
		const char* file_out;
		bool show;
		vector<unsigned int> channels;
		vector<unsigned int> frame;
		int mode;
		const char* scale_factor;	
		bool make_normalization_image;
		bool use_fourrier;
		float val;
		float clamprange[2];
		const char* trans;
		char type_out;
		bool wantrgb;
		unsigned int ndim[2];
		bool fastresize;
		TASK_MEMBER_DEFINITIONS
	};
	
	template< > class Taskscope<TASK_TIFF_FILE_OPERATION> : public LFHPrimitive::ArgumentParser{public:
		const char* file_in;
		const char* file_out;
		const char* file_dev;
		const char* concat_file;
		bool show;
		vector<unsigned int> channels;
		vector<unsigned int> frame;
		int mode;
		double clamprange[2];
		bool clamp;
		char type_out;
		char oper_ation;
		TASK_MEMBER_DEFINITIONS
	};
	
	template< > class Taskscope<TASK_EXTRACT_HIDDENMAP_DIRECT> : public LFHPrimitive::ArgumentParser{public:
		const char* file_in;
		const char* file_in2;
		const char* file_mcv;
		
		//	bool file_additionnal;
		
		
		
		const char* out_mcv;
		const char* out_tif;
		const char* out_dist;
		
		const char* out_partition;
		const char* out_cpartition;
		const char* out_ccpartition;
		
		//	const char* out_prvw_tif;
		
		bool doem;
		bool quick;
		
		bool show;
		int min_nbpixel;
		
		bool da_max_is_marked;
		
		const char* grametapix;
		
		const char* ellipsedev;
		
		double contour;
		bool do_refine_area;
		const char* Z_mcv;
		const char* Z_txt;
		const char* Z_tif;
		TASK_MEMBER_DEFINITIONS
	};
	
	template< > class Taskscope<TASK_HIDDENMAP_BASE_DATAEXTRACTION_V2> : public LFHPrimitive::ArgumentParser{public:
		const char* file_in[5];
		
		bool show;
		int cellID;
		const char* out_p;
		const char* seg_im_path;
		const char* residual;
		const char* table_out;
		
		double red_substract_fact;
		
		const char* mcv_in; 
		double distance_thr;
		const char* file_dascp;
		double artifact_class[5];
		bool ignore_red_in_confidence;
		const char* filter_cell_class; 
		bool recursive;
		
		double art_k;
		TASK_MEMBER_DEFINITIONS
	};
	
	
	template< > class Taskscope<TASK_HIDDENMAP_BASE_DATAEXTRACTION> : public LFHPrimitive::ArgumentParser{public:
		const char* file_in[5];
		
		bool show;
		int cellID;
		const char* out_p;
		const char* seg_im_path;
		const char* residual;
		const char* table_out;
		bool ignore_red_in_confidence;
		double red_substract_fact;
		
		const char* mcv_in; 
		double distance_thr;
		const char* file_dascp;
		double artifact_class[5];
		
		const char* filter_cell_class; 
		bool recursive;
		
		double art_k;
		
		TASK_MEMBER_DEFINITIONS
	};
	
	
	template< > class Taskscope<TASK_MAKE_DISPLAYABLE> : public LFHPrimitive::ArgumentParser{public:
		const char* file_in[5];
		bool show;bool inter;
		int torun;
		double k_dens;
		
		TASK_MEMBER_DEFINITIONS
	};
	
template< > class Taskscope<TASK_GETLABEL_CONFIDENCE> : public LFHPrimitive::ArgumentParser{public:
	const char* file_in[2];
	bool show;
	int torun;
	const char* table_list;
	
	TASK_MEMBER_DEFINITIONS
};

	
	template< > class Taskscope<TASK_PROFILES> : public LFHPrimitive::ArgumentParser{public:
		const char* file_in[3];
		const char* hypercl;
		const char* hypercl_out;
		bool hypercl_near;
		bool show,cluster;
		bool inter;
		bool no_daugther;
		double no_confidence;
		
		bool does_bining;
		double daval_sepsep;
		double daval_sepsepi;
		
		const char* pair_wise;
		const char* svm_joachims;
		double factor[6];
		
		
		TASK_MEMBER_DEFINITIONS
	};
	
	
	template< > class Taskscope<TASK_CDTGTR_PVALUES> : public LFHPrimitive::ArgumentParser{public:	
		const char* filter_file;
		const char* permute_file;
		bool show, inter, permute, w_name,gray,scaleR,scaleG,scaleB,doopen;
		double factR,factG,factB;
		const char* cdtcolname;
		
		unsigned int nb_annot_range[2];
		unsigned int nb_permutations;
		int torun;
		unsigned int dims[2];
		bool dodel;
		const char* threshold_path;
		double threshold;
		TASK_MEMBER_DEFINITIONS
	};	
	
	template< > class Taskscope<TASK_RANK_PVALUES> : public LFHPrimitive::ArgumentParser{public:	
		const char* filter_file;
		const char* permute_file;
		bool show, inter, permute, w_name,gray,doopen;
		const char* cdtcolname;
		bool reverse;
		unsigned int nb_annot_range[2];
		unsigned int nb_permutations;
		int torun;
		unsigned int dims[2];
		bool dodel;
		const char* threshold_path;
		double threshold;
		TASK_MEMBER_DEFINITIONS
	};	

	template< > class Taskscope<TASK_MULTI_PAIR_DISTANCE_PVAL> : public LFHPrimitive::ArgumentParser{public:	
		const char* filter_file;
		const char* permute_file;
		bool show, inter, permute, w_name,gray,scaleR,scaleG,scaleB,doopen;
		double factR,factG,factB;
		const char* cdtcolname;
		
		unsigned int nb_annot_range[2];
		unsigned int nb_permutations;
		int torun;
		unsigned int dims[2];
		bool dodel;
		const char* threshold_path;
		double threshold;
		TASK_MEMBER_DEFINITIONS
	};

