/*
 * Segmentation.cpp
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

#include "Tasks.h"
#include "Madstructs.h"

	bool isToLabel_seg(const double & a){return a > 0.1f;}
	
		Taskscope<TASK_HMM_SEGMENTATION_AND_DISTANCE>::Taskscope(): show(false), file_average(NULL), file_stddev(NULL), min(0.0f) , max(0.0f),nbclasses(2),file_dist(NULL),file_out(NULL),file_scope(NULL),artfrac(0.001f),file_nohmm(NULL),blur_wind(2.0f),bright(false){
			crop_rect[0] = 0xFFFFFFFF;
		}
		void Taskscope<TASK_HMM_SEGMENTATION_AND_DISTANCE>::nbaddtoken(char const * const token, int& min, int& max){
			switch(*token){
				case '\0': min =1; max=2; break;
				case 's': min =0; break;
				case 'a': min =1; break;
				case 'v': min =1; break;
				case 'M': min =1; break;
				case 'K': min =4; break;
				case 'd': min =1; break;
				case 'n': min =1; break;
				case 'c': min =2; break;
				case 'b': min =1; break;
				case 'u': min =1; break;
				case 'B': min =1; break;
				case 'D': min =0; break;
			}
		}
		void Taskscope<TASK_HMM_SEGMENTATION_AND_DISTANCE>::store(char* const * token, int nbtoken){ 
			
			switch(token[0][1]){
				case 's': show =true;break;
				case 'a': file_average = token[1];break;
				case 'v': file_stddev = token[1];break;
				case 'M': file_nohmm = token[1];break;
				case 'c': min = atof(token[1]);
					max = atof(token[2]);
					break;
				case 'd': file_dist = token[1];
					break;
				case 'n': nbclasses = atoi(token[1]);break;
				case 'b': file_scope = token[1];break;
					break;
				case 'B': blur_wind = atof(token[1]);break;
					break;
				case 'u': artfrac = atof(token[1]);
					break;
				case 'K': crop_rect[0] = atoi(token[1]);crop_rect[1] = atoi(token[2]);crop_rect[2] = atoi(token[3]);crop_rect[3] = atoi(token[4]); break;
				case 'D': bright = true; break;
			}
			
		}
		int Taskscope<TASK_HMM_SEGMENTATION_AND_DISTANCE>::defstore(char* const * token, int nbtoken){ 
			if (nbtoken) file_in = token[0];
			if (nbtoken>1) file_out = token[1];
			
			unsigned int i,j,k,l;
			char buffer[65536];
			//			 char* fakeargs[] = {"./PureMadness","-a","../../../../HOwt/Average.tif", "-b", "../../../../superscope.scp", "-d", "../../../../eviltwo.tif" , "../../../../HOwt/HOwt_plate01_001015_00.tif", "../../../../evilseg.tif"};
			//			 mmm(9,fakeargs);
			//"-a","../../../../Average_002.tif", 
			//	char* fakeargs[] = {"./PureMadness","-b", "../../../../003003_sg.scp" ,"-d", "../../../../evilthree.tif","-M", "../../../../evilfour.tif","../../../../003003_00.tif", "../../../../eviltwo.tif"};
			//	mmm(9,fakeargs);
			
			SerialStore< unsigned int > scope(file_scope);
			
			
			TiffFile bypassfile(file_nohmm, true);

			unsigned int coor[2];
			
			TiffFile distfile(file_dist, true);
			
			
			DataGrid<double,2> transition;
			if (file_dist){
				coor[0] = 2;
				coor[1] = 2;
				transition.setSizes(coor);
			}
			DataGrid<double,2> profile[2];
			
			strcpy(buffer, file_in);
			
			TiffFile input_file( buffer);
			DataGrid<double,2> raw_image;
			
			TiffFile ave_file( file_average);
			TiffFile std_file( file_stddev);
			
			vector<Madstructs::Image<float> *> output;
			
			
			if (file_scope){
				
				if (scope.has(9)){
					scope.load(9, blur_wind);
					printf("Using stored Blur window size %f\n", blur_wind);
				}else{
					scope.save(9, blur_wind);
				}
				
			}
			
			for(k=0; input_file.fetch_grayscale(raw_image);k++){
				HMM<2> hmm;
				
				LFHPrimitive::Classifier< Tuple<double, 2> ,2 > classif;
				
				if (file_average) ave_file.fetch(profile[0]);
				if (file_stddev) ave_file.fetch(profile[1]);
				
				
				GaussianDistribution<2> foreg;
				
				
				GaussianDistribution<2> backg;
				classif[0] = &backg; 
				classif[1] = &foreg;
				Tuple<double, 2> tmpdp;
				Tuple<double,1> tmptmptmp;
				double pix[32];
				DataGrid<Tuple<double, 1>, 2> image;
				DataGrid< Tuple<double, 2> , 2> classimage;
				
				DataGrid<Tuple<double, 3>, 2> rmap;
				unsigned int dims[3];
				dims[0] = raw_image.dims[0];
				dims[1] = raw_image.dims[1];
                                Tuple< double, 3> prepbound;
                                Tuple< double, 2> bound;

                                image.setSizes(dims);
				
				
				
				classimage.setSizes(dims);
				
				double minmax[2];
				
				minmax[0] = 10000.0f;
				minmax[1] = 0.0f;
				Tuple< double, 1> pixel;
				
				for(coor[1]=0;coor[1]<raw_image.dims[1];coor[1]++){
					for(coor[0]=0;coor[0]<raw_image.dims[0];coor[0]++){
						pixel[0] = raw_image(coor);
						
						if (file_average){
							dims[0] = (int)((coor[0] * profile[0].dims[0]) / raw_image.dims[0]);
							dims[1] = (int)((coor[1] * profile[0].dims[1]) / raw_image.dims[1]);
							pixel[0] -=  profile[0](dims);
						}
						if (file_stddev){
							dims[0] = (int)((coor[0] * profile[1].dims[0]) / raw_image.dims[0]);
							dims[1] = (int)((coor[1] * profile[1].dims[1]) / raw_image.dims[1]);
							pixel[0] /= profile[1](dims);
						}
						
						if (min != max){
							if (min > pixel[0]) pixel[0] = min;
							if (pixel[0] > max) pixel[0] = max;
						}
						image(coor) = pixel; //if (k == 1) printf("%f\n", pixel[0]);
					}
				}
				DataGrid<Tuple<double, 2>, 2> t_image;
				dims[0] = raw_image.dims[0]; 
				dims[1] = raw_image.dims[1];
				t_image.setSizes(dims);
				DataGrid<Tuple<double, 1>, 2> b_image = image.crudeGaussianBlur(blur_wind);
				if (bright){
					DataGrid<double , 2> image_grad = b_image.makeGradientNorm();image_grad *= -1.0f;
					DataGrid<unsigned int , 2> image_index = image_grad.SegementClimbToMax(j);
					
					//Tuple<double, 2>* minmax = new Tuple<double, 2>[j];
					//for(l=0;l<j;l++){ExOp::toMax(minmax[l][0]); ExOp::toMin(minmax[l][1]);}
					
					WeightElem<double, 2>* minmax_st = new WeightElem<double, 2>[j];
					for(l=0;l<j;l++) ExOp::toZero(minmax_st[l]);
					DataGrid<unsigned int , 2>::KeyIterator k_ite_i = image_index.getKeyIterator();
					Tuple<double,2> input_data;
					
					
					if ( k_ite_i.first())do{
						l = image_index(k_ite_i());// printf("%i\n", l);
						//if (minmax_st[l][0] > image(k_ite_i())[0]) minmax_st[l][0] = image(k_ite_i())[0];
						//if (minmax_st[l][1] < image(k_ite_i())[0]) minmax_st[l][1] = image(k_ite_i())[0];
						minmax_st[l] += WeightElem<double, 2>(image(k_ite_i())[0]);
					} while ( k_ite_i.next());
					
					if ( k_ite_i.first())do{
						l = image_index(k_ite_i());
						input_data[0] = image(k_ite_i())[0];
					//	input_data[1] = minmax_st[l][1] - minmax_st[l][0];
						input_data[1] = sqrt(minmax_st[l].getVar());
						if (!ExOp::isValid(input_data[1])) input_data[1] = 0.0f;
						image_grad(k_ite_i()) = input_data[1];
						t_image(k_ite_i()) = input_data;
					} while ( k_ite_i.next());
					

					delete[](minmax_st);
					
			//		TiffFile tf_tf_tf("./tmptmp.tif");
			//		tf_tf_tf.put(image_grad, (float)0.0f, (float)1.0f);
					}else{
					t_image( Tuple<double,2>::Concatenate<1>(), image, b_image);				
				}
				

				
				GaussianDistribution<2> totalDistr;
				totalDistr.EMinit();
				for(coor[1]=0;coor[1] < dims[1];coor[1]++) for(coor[0]=0;coor[0] < dims[0];coor[0]++){
					tmpdp = t_image(coor);
					
					for(j=0;j<2;j++) totalDistr.EMregist(tmpdp,1.0f);
					
				}
				totalDistr.EMfinit();
				printf("Total Distribution:\n");
				totalDistr.show(stdout);
				
				
				if ((file_scope)&&(scope.has(k*10))){
					printf("Reusing saved HMM parameters for frame %i.\n",k);
					scope.load(k*10, backg);
					scope.load(1+k*10, foreg);
					if (artfrac != 0.0f){
						classif.setUnknownProbability(0.001f);
						scope.load(2+k*10,classif.unknown_scope[0]);
					}
					scope.load(3+k*10,hmm);
				}else{
					printf("Finding HMM parameters for frame %i:\n",k);
					
					
					
					if (crop_rect[0] != 0xFFFFFFFF){
					
					
					}
					foreg = totalDistr;
					backg =totalDistr;
					foreg.mean[0] *=1.2f;
					backg.mean[0] *=0.8f;
					
					bound[0] =0.9f;
					bound[1] =0.1f;
					
					hmm.init(0.1f,bound);
					
					if (artfrac != 0.0f) classif.setUnknownProbability(artfrac);
					double LLout;
					
					printf("EM start:\n");
					classif.EMinit();
					for(l=0,i=10;i<0x80000000;i--,l++){
						
						classimage(classif, t_image);
						
						hmm.EMinit();
						
						DataGrid<Tuple<double, 2>, 2> rmap2 =  hmm.runHMM(hmm.runHMM(classimage,0),1,true);
						DataGrid<Tuple<double, 2>, 2> rmap3 =  hmm.runHMM(hmm.runHMM(classimage,1),0,true);
						
						hmm.EMfinit();
						
						classif.EMAlphainit(0.5f);
						
						for(coor[1]=0;coor[1] < dims[1];coor[1]++) for(coor[0]=0;coor[0] < dims[0];coor[0]++){
							bound = rmap2(coor);
							bound += rmap3(coor);
							bound *= 0.5f;
							tmpdp = t_image(coor);
							classif.EMregist(tmpdp, bound);
						} 
						
						
						LLout = classif.EMfinit();
						if (classif.unknown_scope == NULL) printf("Log-Likelyhood: %e\n", LLout);
						else{
						if ((classif.relerr > 0.2f)||(classif.relerr < -0.5f)) {
							i++;
						}
							printf("Log-Likelyhood: %e", LLout);
							if (classif.relerr > 0.2f) printf("(too many identified artifacts)\n");
							else printf("\n");
						}
						if ((!LFHPrimitive::ExCo<double>::isValid(hmm.boundary[0]))||(l > 200)){
							fflush(stdout);
							LFHPrimitive::static_warning_handdle << LFH_WARNING_MAXITE;
							break;
						}
					}
					printf("EM end:\n");
					
					//	saveBMP("./tmptmp.bmp",out);
					
					if (((GaussianDistribution<2> *)classif[0])->mean[0] > ((GaussianDistribution<2> *)classif[1])->mean[0]){
						classif[0] = &foreg; 
						classif[1] = &backg;
						hmm.swapstates(1,0);
					}
					
					if (file_scope){
						scope.save(k*10, *((GaussianDistribution<2> *)classif[0]));
						scope.save(1+k*10, *((GaussianDistribution<2> *)classif[1]));
						scope.save(2+k*10, classif.unknown_scope[0]);
						scope.save(3+k*10, hmm);
						scope.flush();
					}
				}
				
				
				//		LFHPrimitive::Classifier<Tuple<double, 2>, 2> classif_final;
				
				
				hmm.transition.show(stdout);	
				hmm.boundary.show(stdout);
				
				printf("Background Distribution:\n");((GaussianDistribution<2> *)classif[0])->show(stdout);
				printf("Foreground Distribution:\n");((GaussianDistribution<2> *)classif[1])->show(stdout);
				
				classimage(classif, t_image);
				if (file_nohmm != NULL){
					DataGrid<double, 3> da_da_conv(classimage);
					DataGrid<double, 2> da_da_conv2 = da_da_conv.makeSlice(0,0);
					bypassfile.put(da_da_conv2, (float)0.0f, (float)1.0f, false);
				}
				
				DataGrid<Tuple<double, 2>, 2> rmap2 =  hmm.runHMM(hmm.runHMM(classimage,0),1,true);
				DataGrid<Tuple<double, 2>, 2> rmap3 =  hmm.runHMM(hmm.runHMM(classimage,1),0,true);
				
				
				rmap2 += rmap3;
				rmap2 *= 0.5f;
				
				for(coor[1]=0;coor[1] < dims[1];coor[1]++) for(coor[0]=0;coor[0] < dims[0];coor[0]++){
					bound = rmap2(coor);
					tmpdp = t_image(coor);
					pix[1] = classif.UnknownProbability(tmpdp, bound);
					pix[2] = bound[0] * (1.0f - pix[1]);
					pix[0] = bound[1] * (1.0f - pix[1]);
					
					//		test_output[k]->setPixel(coor[0],coor[1],pix);
					
				}
				
				
				b_image(Tuple<double,1>::Deconcatenate<1>(), image,rmap2);
				Madstructs::Image<float>* im_out = new Madstructs::Image<float>(image);
				output.push_back(im_out);
				
				
				
				
				if (file_dist){
					
					// im_out segmented image!
					DataGrid<double,2> distance;
					DataGrid<double,3> classimage_conv(classimage);
					
					coor[0] = 0;
					coor[1] = 0;
					transition(coor) = hmm.transition.data[0];
					coor[0] = 1;
					transition(coor) = hmm.transition.data[1];
					coor[1] = 1;
					transition(coor) = hmm.transition.data[3];
					coor[0] = 0;
					transition(coor) = hmm.transition.data[2];
					distance.ExpectedDistanceToOther(classimage_conv,transition,true);
					
					//		Vector<KeyElem<double, double > > datapts;
					//		typename DataGrid<double,2>::KeyIterator dataite =  distance.getKeyIterator();
					
					//		if (dataite.first()) do{
					//			tmptmptmp = image(dataite());
					//			printf("%f\t%f\n", distance(dataite()), tmptmptmp[0]);
					//		}while(dataite.next());
					distfile.put(distance, (float)0.0f, (float)100.0f);
					
				}
				
			}
			
			if (file_out) {
				strcpy(buffer, file_out);
				
			Madstructs::Image<float>::SaveTiffImage(output, buffer);}
			//	Madstructs::Image<unsigned char>::SaveTiffImage(test_output,"./wannasee.tif");
			if ((show)&&(file_out)){
				strcpy( buffer, "open -a Preview \"");
				strcpy(buffer + 17 , file_out);
				strcpy( buffer + 17 + strlen(file_out), "\"");
				system(buffer);
			}
			
		return 0;}
		void Taskscope<TASK_HMM_SEGMENTATION_AND_DISTANCE>::help(){
			printf("Foreground/Background Segmentation, for background and foreground of distinct intensity. 16bit images assumed\n");
			printf("\n");
			printf("Step 1: Deterniming starting guess for pixel classes.\n");
			printf("Step 2: Running EM-HMM.\n");
			printf("Arguments [1 required]:\n\n");
			printf("\t[in](tiff file) (grayscale)\n");
			printf("\t[out](tiff file) (grayscale, floating point):\n");
			printf("\nPreprocessing Flags\n");
			printf("\t-a (file path):\tMean Intensity (float) Image: to be subtracted to curent image prior to segmentation\n");
			printf("\t-v (file path):\tStdDev Intensity (float) Image: to normalize the image prior to segmentation\n");
			printf("\t-K (int x) (int y) (int width) (int height): Crop rectangle for 2D-HMM parameter inference\n");
			printf("\t-c (float min) (float max):\t clamp range which is applied to input image\n");
			printf("\nOther Flags\n\n");
			printf("\t-s:\tShow output in Preview\n");
			printf("\t-u (float):\t fraction of artifact pixels (default =0.001)\n");
			printf("\t-b (FILE ):\t scope file [store/reuse information for quick recalculations]\n");
			printf("\t-B (float):\t Blur Window size (2 pixel default, additionnal info to capture punctate noise by the 2D-HMM)\n");
			printf("\t-D: DIC (bright field) microscopy attempt [not-working]\n");
			printf("\nAdditionnal Output Flags\n\n");
			printf("\t-d (FILE ):\t compute physical distances under HMM\n");
			printf("\t-M (FILE ):\t marginal, (bypasses the pseudo 2D-HMM)\n");
			printf("Version 1.0\n");
		}

	
	
	
Taskscope<TASK_GEN_META_PIXELS>::Taskscope(): show(false), min_nbpixel(0), mask_file(NULL), out_tif(NULL), comp(false),quick(false),out_partition(NULL),out_cpartition(NULL),grametapix(NULL),blur_std(0.0f),knight(false), distpen(false),mark_da_maxima(false){
	valid_int_range[0] = 1.0f;
	valid_int_range[1] = 0.0f;
}

void Taskscope<TASK_GEN_META_PIXELS>::nbaddtoken(char const * const token, int& min, int& max){
	switch(*token){
		case '\0': min =2; break;
		case 's': min =0; break;
		case 'm': min =1;break;
		case 'H': min =0;break;
		case 'D': min =1;break;
		case 'p': min =0;break;
		case 'k': min =0;break;
		case 'C': min =0;break;
		case 'b': min =1;break;
		case 'M': min =0;break;
		case 'R': min =2;break;
			//			case 'p': min =1;break;
	}
}
 void Taskscope<TASK_GEN_META_PIXELS>::store(char* const * token, int nbtoken){ 
	
	switch(token[0][1]){
		case 'C':
			comp = true;
			break;
		case 'k': knight = true; break;
		case 'p': distpen = true; break;
		case 'D':
			out_dist = token[1];
			break;
		case 'P':
			out_partition = token[1];
			break;
			//			case 'p':
			//				out_prvw_tif = token[1];
			//				break;
		case 'b': blur_std = atof(token[1]); break;
		case 'm': mask_file =  token[1]; break;
		case 'f': min_nbpixel = atoi(token[1]); break;
		case 'M': mark_da_maxima = true;break;
		case 'G': grametapix = token[1]; break;
		case 'R': valid_int_range[0] = atof( token[1]); valid_int_range[1] = atof( token[2]); break;
			
	}
	
}

void Taskscope<TASK_GEN_META_PIXELS>::help(){
	printf("Generates Unseeded Watershed Segmention.\n");
	printf("\n");
	printf("\t(in) Input Image File\n"); 
	printf("\t(in) Output Image File (32-bit TIFF)\n"); 
	printf("Flags\n\n");
	printf("\t-C: cluster metapixel further\n");
	printf("\t-b (double): blur image before clustering (gausian with input STD)\n");
	printf("\t-m (Tiff File): Appy mask before clustering\n");
	printf("\t-k: consider knight moves\n");
	printf("\t-p: penalize by transition distance\n");
	printf("\t-M: Mark Local Maxima\n");
	printf("\t-R (float min) (float max) : Computes the Watershed only on pixel with intensity in the given range: any other pixel has its own catchment bassin.\n");
	
	printf("Version 1.0\n");
	
}

int Taskscope<TASK_GEN_META_PIXELS>::defstore(char* const * token, int nbtoken){ 
	file_in = token[0];
	file_out = token[1];
	
	TiffFile tf_in(file_in);
	TiffFile tf_out(file_out, true);
	TiffFile tf_mask(mask_file);
	//	TiffFile tf_test("./testest.tif");
	
	
	DataGrid<double,3> im_in;
	
	DataGrid<double,2> im_in_inter;
	DataGrid<double,2> im_in_inter_blur;
	
	DataGrid<bool,2> im_filter;
	
	DataGrid<unsigned int,2> im_out;
	unsigned int coor[2];
	while(tf_in.fetch(im_in)){
		im_in_inter = im_in.makeSlice(0,0);
		unsigned int nbmetpix;
		DataGrid<double,2>& curent = im_in_inter;
		
		
		if (mask_file) {
			tf_mask.fetch(im_in);
			im_in_inter_blur = im_in.makeSlice(0,0);
			curent *= im_in_inter_blur;
		}
		
		if (valid_int_range[0] < valid_int_range[1]){
			im_filter.setSizes(curent.dims);
			for(coor[1]=0;coor[1] < curent.dims[1];coor[1]++) for(coor[0]=0;coor[0] < curent.dims[0];coor[0]++) im_filter(coor) = ((curent(coor) > valid_int_range[0])&&(curent(coor) < valid_int_range[1]));
		}
		
		if (blur_std != 0) {im_in_inter_blur = curent.crudeGaussianBlur(blur_std,false);curent = im_in_inter_blur;}
		//		tf_test.put(curent, (float)0, (float)65500);
		
		//printf("dfasdfa %c", mark_da_maxima ? 'Y' : 'N');

		
		if (comp) im_out = curent.SegementClimbToMaxMax(nbmetpix);
		else im_out = curent.SegementClimbToMax(nbmetpix, distpen , knight, false, mark_da_maxima, (valid_int_range[0] < valid_int_range[1]) ? &im_filter : NULL);
		
		
		tf_out.put(im_out, (unsigned int)0, nbmetpix);
	}
	return 0;
}

	
	Taskscope<TASK_FIND_CELL_COVER_UPGRADE>::Taskscope(): show(false),max_dist(100), out_mcv(NULL), out_txt(NULL), hint_max_width(0), out_tif(NULL), doem(true),quick(false),scope_file(NULL),in_mcv_known(NULL),whichheu(0){contour[1] =0.0f;}
	void Taskscope<TASK_FIND_CELL_COVER_UPGRADE>::nbaddtoken(char const * const token, int& min, int& max){
		switch(*token){
			case '\0': min =3; break;
			case 's': min =0; break;
			case 'm': min =1;break;
			case 'o': min =1;break;
			case 'M': min =1;break;
			case 'H': min =0;break;
			case 'B': min =0;break;
			case 'R': min =2;break;
				//			case 'p': min =1;break;
		}
	}
	void Taskscope<TASK_FIND_CELL_COVER_UPGRADE>::store(char* const * token, int nbtoken){ 
		
		switch(token[0][1]){
			case 'o':show = true;
			{
				int l = strlen(token[1])-4;
				if (strcmp(token[1]+ l, ".mcv") == 0) out_mcv = token[1];
				else if (strcmp(token[1] + l, ".tif") == 0) out_tif = token[1];
				else if (strcmp(token[1] + l, ".txt") == 0) out_txt = token[1];
			}
				break;
			case 'M':	hint_max_width =atoi(token[1]); break;
			case 'm':
				in_mcv_known = token[1];
				break;
			case 'H':
				doem = false;
				whichheu = (token[0][2] == '2') ? 2 : 1;
				break;
			case 'B':
				quick = true;
				break;
				//			case 'p':
				//				out_prvw_tif = token[1];
				//				break;
			case 'R': 
				contour[0] = atof(token[1]); 
				contour[1] = atof(token[2]);
				whichheu = (token[0][2] == '2') ? 4 : 3;
				break;
				
		}
		
	}
	int Taskscope<TASK_FIND_CELL_COVER_UPGRADE>::defstore(char* const * token, int nbtoken){ 
		file_in = token[0];
		file_in2 = token[1];
		file_in3 = token[2];
		ProgressBarPrint daprogbar; daprogbar.lenght = 20;
		unsigned int i,j,k,l;
		char buffer[65336];
		
		double pix[32];
		unsigned int dims[2];
		Tuple<unsigned int, 2> coors;
		Tuple<unsigned int, 2> acoors;
		vector<Madstructs::Image<float>* > seg;
		vector<Madstructs::Image<unsigned short>* > ima;
		
		vector<Madstructs::Image<float>* > hiddenmap;
		
		printf("Loading Images!\n"); fflush(stdout);
		Madstructs::Image<unsigned short>::LoadTiffImage(ima, file_in);
		printf("Raw Image Read!\n"); fflush(stdout);
		Madstructs::Image<float>::LoadTiffImage(seg, file_in2);
		printf("Segmented File Read!\n"); fflush(stdout);
		printf("Loading Images Done!\n"); fflush(stdout);
		
		
		
		TiffFile tf_dist(file_in3);
		
		// Z-state from 2D HMM
		// greedy map from dynamic programming
		// 
		
		Tuple<unsigned int, 3> metacoor; metacoor[0] = 0;
		
		FILE * fmcvin = (in_mcv_known) ? fopen(in_mcv_known,"rb+") : NULL; 
		
		Tuple<double, 2> prob;
		Tuple<double, 1> inten;
		
		FILE* mcvf = (out_mcv) ? fopen(out_mcv,"wb+") : NULL;
		
		FILE* mcvtxt = (out_txt) ? fopen(out_txt,"w+") : NULL;
		
		printf("Initialization Done!\n"); fflush(stdout);
		
		DataGrid<double, 2> distanceBeta;
		for(k=0;k<seg.size();k++){
			DataGrid<double, 2> segmented;
			DataGrid<double, 2> distance;
			DataGrid<double, 3> distance_in;
			DataGrid<unsigned int, 2> labels; // = segmented.LabelConnected(&isIn,NULL,NULL);
			dims[0] = seg[k]->sizex;
			dims[1] = seg[k]->sizey;
			
			if (!tf_dist.fetch(distance_in)) {printf("Incohenrent number of frame in tiff files: '%s' '%s' (critical)\n", file_in2, file_in3);exit(1);}
			segmented.setSizes(dims);
			distance.setSizes(dims);
			labels.setSizes(dims);
			
			for(coors[1] = 0;coors[1]<dims[1];coors[1]++) for(coors[0] = 0;coors[0]<dims[0];coors[0]++) {
				ima[k]->getPixel(coors[0],coors[1],pix);
				seg[k]->getPixel(coors[0],coors[1],pix);
				segmented(coors) = pix[0];
				metacoor[1] = coors[0];
				metacoor[2] = coors[1];
				distance(coors) = distance_in(metacoor); 
			}
			
			printf("done!\n"); fflush(stdout);
			

			

				
				
			for(coors[1] = 0;coors[1]<dims[1];coors[1]++) for(coors[0] = 0;coors[0]<dims[0];coors[0]++) distance(coors) *= segmented(coors);
			
			Madstructs::MultiCover mc;
			
			if (in_mcv_known) {
				mc.load(fmcvin);
				printf("Loaded Multicover File!\n"); fflush(stdout);
			}else{
				
				
				ExOp::toZero(labels);
				Tuple<unsigned int,4> bounds;
				vector<Tuple<unsigned int,4> > rects;
				vector<int > areas;
				vector<double> maxdist;
				vector< Tuple<unsigned int,2> > buf;
				
				
				
				
				daprogbar.start("Cutting Foreground into clumps");

				for(coors[1] = 0;coors[1]<dims[1];coors[1]++) {daprogbar.update(((double)coors[1]) / dims[1]); for(coors[0] = 0;coors[0]<dims[0];coors[0]++){
					
					
					
					seg[k]->getPixel(coors[0],coors[1],pix);
					
					if ((pix[0]  > 0.1f)&&(labels(coors) == 0)){
						pix[31] = pix[0];
						buf.push_back(coors);
						bounds[0] = coors[0];
						bounds[1] = coors[1];
						bounds[2] = coors[0];
						bounds[3] = coors[1];
						l=0;
						while(buf.size()!= 0){
							acoors = *(buf.end()-1);buf.pop_back();
							seg[k]->getPixel(acoors[0],acoors[1],pix);
							
							if  ((pix[0]  > 0.1f)&&(labels(acoors)  == 0)){
								labels(acoors) = rects.size()+1;
								if (pix[31] < pix[0]) pix[31] = pix[0];
								l++;
								//		setPixel(x,y,pix);
								if (acoors[0] > 0) {buf.push_back(acoors); (*(buf.end()-1))[0]--;}
								if (acoors[0] < dims[0]-1) {buf.push_back(acoors);(*(buf.end()-1))[0]++;}
								if (acoors[1] > 0) {buf.push_back(acoors);(*(buf.end()-1))[1]--;}
								if (acoors[1] < dims[1]-1) {buf.push_back(acoors);(*(buf.end()-1))[1]++;}
								if (acoors[0]< bounds[0]) bounds[0]= acoors[0];
								if (acoors[1]< bounds[1]) bounds[1] =acoors[1];
								if (acoors[0]> bounds[2]) bounds[2]=acoors[0];
								if (acoors[1]> bounds[3]) bounds[3]=acoors[1];
							}
						}
						areas.push_back(l);
						bounds[2] -=  bounds[0] -1;
						bounds[3] -=  bounds[1] -1;
						rects.push_back(bounds);
						maxdist.push_back(pix[31]);
					}
				} }daprogbar.finish();

				printf("Finding Cells in frame %i! there is %i groups\n",k,rects.size()); fflush(stdout);
				
				// first pass, find cell
				daprogbar.start("Cutting clumps into cells");
				for( l=0;l<rects.size();l++){
					daprogbar.update(((double)l) / rects.size());
					//					printf("%f\n", pix[0]);
					if (areas[l] <25){  // is too small
					//if ((areas[l] <25)||(rects[l][2] == seg[k]->sizex)||(rects[l][3] == seg[k]->sizey) ){ // is too small, or whole/row or collumn
						// printf("%i\t%i\t%i\t%i\t (area %i, maxdist %f) rejected!\n", rects[l][0],rects[l][1],rects[l][2],rects[l][3],areas[l],maxdist[l]);fflush(stdout);
					}else{
						//		printf("rect%i:\t",l);
						//		printf("%i\t%i\t%i\t%i (area %i, maxdist %f) \n", rects[l][0],rects[l][1],rects[l][2],rects[l][3],areas[l],maxdist[l]);
						Madstructs::Image<double> subim;
						subim.sizex = rects[l][2];
						subim.sizey = rects[l][3];
						subim.channels =4;
						subim.allocateBuffer();
						
						for(coors[1] = rects[l][1];coors[1]<rects[l][1] + rects[l][3];coors[1]++) for(coors[0] = rects[l][0];coors[0]<rects[l][0] +rects[l][2];coors[0]++){
							if (labels(coors) == l+1){
								ima[k]->getPixel(coors[0],coors[1],pix+4);
								seg[k]->getPixel(coors[0],coors[1],pix);
								pix[2] = distance(coors);
								pix[1] = pix[4];
								if (!ExCo<double>::isValid(pix[1])) {pix[1] =0; printf("imposibble!!!!!!\n");}
								subim.setPixel(coors[0] - rects[l][0],coors[1] - rects[l][1], pix);
							}else{
								memset(pix,'\0', sizeof(double)*4);
								subim.setPixel(coors[0] - rects[l][0],coors[1] - rects[l][1], pix);
							}
						}
						
						
						Madstructs::CellCover cc;
						switch(whichheu){
							case 1:	cc.findCellCoverHeuristic_Sep2011(&subim,0.25f,NULL);break;
							case 2: cc.findCellCoverHeuristic_agglo2(&subim, 12.0f);break;
							case 3: 
								//	printf("in\n"); fflush(stdout);
								cc.findCellCoverAlgerabric(&subim, contour[0], contour[1]); 
								//	printf("out\n"); fflush(stdout);
								break;
							case 4: 
								cc.findCellCoverAlgerabricHeuristic(&subim, contour[0], contour[1],0.25f,100, hint_max_width); 
								break;
						}
						
						//		if (contour != 0.0f) cc.update_from_contour(&subim, contour);
						
						
						
						
						//	if (!quick) cc.findCellCover_hiddenmap2(&subim,false, (out_tif == NULL) ? NULL : &hiddenmap, false);
						
						
						
						
						if (cc.cover.size() > 0) {
							cc.rect[0]= rects[l][0];
							cc.rect[1]= rects[l][1];
							cc.rect[2]= rects[l][2];
							cc.rect[3]= rects[l][3];
							mc.group.push_back(cc);
						}
						
						
						
					}
					
				}daprogbar.finish();
				
				
				
				
				
				
				
				printf("Finding Cells done!\n"); fflush(stdout);
				
			}
			
			for(l=0;l< mc.group.size();l++){
				
				
			}
			
			
			
			if (out_txt) {
				fprintf(mcvtxt, "Frame No %i:\n", k+1);
				mc.show(mcvtxt);
			}
			if (out_mcv) {mc.save(mcvf);
				printf("%i groups recorded!\n",mc.group.size());
			}
			
			
			
			
			
			
			
		}
		
		//		if (out_prvw_tif){
		//			TiffFile tfpo(out_prvw_tif);
		
		
		//		}
		
		if (out_mcv) fclose(mcvf);
		if (out_txt) fclose(mcvtxt);
		
		
		if (out_tif){
			Madstructs::Image<float>::SaveTiffImage(hiddenmap, out_tif);
		}
		
		
		
	return 0;}
	void Taskscope<TASK_FIND_CELL_COVER_UPGRADE>::help(){
		printf("Finds the circle cover in the segmented image\n modeling the pixel intensity within cells to be different in the raw image.\n");
		printf("\n");
		printf("\t(in) Raw Image\n");
		printf("\t(in) Segmented Image File\n");
		printf("\t(in) Distance Image File\n");
		printf("Flags\n\n");
		//					printf("\t-m (int):\tMaximum distance searched (default = 100)\n");
		printf("\t-o (filepath).mcv: Multicover file output\n");
		printf("\t-o (filepath).txt: Table of found cells output\n");
		printf("\t-m (filepath).mcv: provides the mcv, does recompute\n");
		
		printf("\t-H1: Heuritic centers, max areas for Edge dist\n");
		printf("\t-H2: Heuritic centers, anti-agglometive from Intensity\n");
		printf("\t-R1 (float min) (float max): Robust regression Algeabric centers, random sampling, min-max range of distance needed from distance file to generate contours.\n");
		printf("\t-R2 (float min) (float max): Robust regression Algeabric centers, Heuristic center sampling, min-max range of distance needed from distance file to generate contours.\n");
		printf("\t-B: Pure Heuritic, (center and radii from dist alone, hidden map cannot be recovered if this flag is used)");
		printf("\t-M (int) : Hint for maximal width of a cell (works with -R2)\n");
		printf("Predecated:\n");
		printf("\t-o (filepath).tif: Tiff Hidden Map file output [predecated]\n");
		printf("\t-p (filepath).tif: Hidden Map Preview [predecated]\n");
		
		printf("Version 1.0\n");
	}



		Taskscope<TASK_EXTRACT_HIDDENMAP_DIRECT>::Taskscope(): do_refine_area(false),show(false), min_nbpixel(0), out_mcv(NULL), out_tif(NULL), doem(true),quick(false),out_partition(NULL),out_cpartition(NULL),out_ccpartition(NULL),grametapix(NULL),da_max_is_marked(false),ellipsedev(NULL),Z_mcv(NULL),Z_tif(NULL){}
		void Taskscope<TASK_EXTRACT_HIDDENMAP_DIRECT>::nbaddtoken(char const * const token, int& min, int& max){
			switch(*token){
				case '\0': min =3; break;
				case 's': min =0; break;
				case 'm': min =1;break;
				case 'H': min =0;break;
				case 'D': min =1;break;
				case 'P': min =1;break;
				case 'G': min =1;break;
				case 'C': min =1;break;
				case 'b': min =1;break;
				case 'c': min =1;break;
				case 'M': min =0; break;
				case 'R': min =0; break;
				case 'd': min =1; break;
				case 'Z': min =3; break;
					//			case 'p': min =1;break;
			}
		}
		void Taskscope<TASK_EXTRACT_HIDDENMAP_DIRECT>::store(char* const * token, int nbtoken){ 
			
			switch(token[0][1]){
				case 'Z':
					Z_mcv= token[1];Z_tif= token[2];Z_txt= token[3];
					break;
				case 'H':
					doem = false;
					break;
				case 'D':
					out_dist = token[1];
					break;
				case 'P':
					out_partition = token[1];
					break;
				case 'd':
					ellipsedev =  token[1];
					break;

				case 'C':
					if (token[0][2] == '\0') out_cpartition = token[1];
					else out_ccpartition = token[1];
					break;
					
					//			case 'p':
					//				out_prvw_tif = token[1];
					//				break;
					//		case 'c': contour = atof(token[1]); break;
				case 'm': out_mcv =  token[1]; break;
				case 'f': min_nbpixel = atoi(token[1]); break;
				case 'G': grametapix = token[1]; break;
				case 'M':da_max_is_marked = true; break;
				case 'R':do_refine_area = true; break;
			}
			
		}
		int Taskscope<TASK_EXTRACT_HIDDENMAP_DIRECT>::defstore(char* const * token, int nbtoken){ 
			file_in = token[0];
			file_in2 = token[1];
			file_mcv = token[2];
			unsigned int i,j,k,l;
			char buffer[65536];
			
			//char* fakeargs[] = {"./PMHiddenMapDirect", "-d", "../../../../report/$_jeb_dev.txt", "-C", "../../../../report/$_jeb_hid.tif", "-M", "-G", "../../../../report/$_t.tif","../../../../report/$_tsg.tif","../../../../report/$_dst.tif", "../../../../report/$_jeb.mcv"};
			//mmm(sizeof(fakeargs) / sizeof(char*) ,fakesubstitution(sizeof(fakeargs) / sizeof(char*) , fakeargs, "HOwt_plate01_001012"));
			
			//			char* fakeargs[] = {"./PMHiddenMapDirect", "-G", "../../../../autmptmp.tif", "-C", "../../../../autmptmp2.tif", "../../../../jebdir/HOwt_plate01_002002_sg.tif", "../../../../jebdir/HOwt_plate01_002002_dt.tif", "../../../../jebdir/HOwt_plate01_002002_ellfit.mcv"};
			//			mmm(8,fakeargs);
			//char* fakeargs[] = {"./PMv2_Circle_fit_heuristic","-o","../../../../autmptmp_new.mcv","-o","../../../../autmptmp_hid.tif","-c", "5.0", "-H","../../../../HOwt/HOwt_plate01_001024_tt.tif","../../../../HOwt/HOwt_plate01_001024_tg.tif","../../../../HOwt/HOwt_plate01_001024_td.tif"};
			//mmm(11,fakeargs);
			//		char* fakeargs[] = {"./PMv2_Circle_fit_heuristic","-o","../../../../autmptmp_new.mcv","-o","../../../../autmptmp_hid.tif","-c", "5.0", "-H","../../../../HOwt/HOwt_plate01_001024_tt.tif","../../../../HOwt/HOwt_plate01_001024_tg.tif","../../../../HOwt/HOwt_plate01_001024_td.tif"};
			
			
			//		char* fakeargs[] = {"./PureMadness", "-P", "../../../../partition.tif", "-C", "../../../../cellpartition.tif", "-CC", "../../../../cellconfpartition.tif","../../../../HOwt/HOwt_plate01_001024_tg.tif","../../../../HOwt/HOwt_plate01_001024_td.tif", "../../../../autmptmp_new.mcv"};
			//		mmm(10,fakeargs);
			
			
			//		mmm(12,fakeargs);
			
			vector<Madstructs::Image<float>* > hiddenmap;
			
			
			Tuple<unsigned int, 3> metacoor; metacoor[0] = 0;
			
			FILE* mcvf = fopen(file_mcv,"rb+");
			FILE* mcvf_Z;
			FILE* txtf_Z;
			if (Z_mcv){
				mcvf_Z = fopen(Z_mcv,"rb+");
				txtf_Z = fopen(Z_txt,"wb+");
			}
			
			FILE* out_mcvf = (out_mcv) ? fopen(out_mcv,"wb+") : NULL;
			TiffFile tfZ = TiffFile(Z_tif);
			
			
			FILE* ell_dev = (ellipsedev) ? fopen(ellipsedev,"w+") : NULL;
			double *ellipse_buf = (ellipsedev) ? new double[6* 100] : NULL;
			if (ellipsedev) fprintf(ell_dev,"CELL_ID\tNB_IN_GROUP\tCENTER_DEVIATION\tELLIPSE_VEC_DEVIATION\tWIDTH_LOG_DEVIATION\tJEB_AREA\tJEB_TYPE\n");
			Vector< Vector< Tuple<double, 2> > > datapts;
			Tuple<double, 2> tmpcenter;
			
			Madstructs::MultiCover mcc;
			Madstructs::CellPose tmpcellpose;
			
			
			Madstructs::MultiCover mcc2;
			
			
			double pix[32]; 
			
			TiffFile tf = TiffFile(file_in);
			TiffFile tf_d = TiffFile(file_in2);
			Tuple<unsigned int ,3> coor;
			
			Tuple<unsigned int ,2> coors;
			
			TiffFile tf_p = TiffFile(out_partition, true);
			TiffFile tf_cp = TiffFile(out_cpartition, true);
			TiffFile tf_ccp = TiffFile(out_ccpartition, true);
			TiffFile tf_metix = TiffFile(grametapix);
			
			DataGrid<double, 3> im_in;
			DataGrid<double, 2> im;
			DataGrid<double, 2> im_d;
			
			DataGrid<unsigned int, 2> im_metix;
			
			DataGrid<unsigned int, 2> labels;
			
			DataGrid<unsigned int, 2> out_hidden;
			
			Vector< pair < Tuple< double, 2> , double> > pts;
			
			
			Vector< Tuple<unsigned int, 4 > > labelrect;
			unsigned int nbrect;
			k=0;
			unsigned int curlab;
			unsigned int metalab;
			unsigned int ID_base=1;
			unsigned int OLDID_base=1;
			unsigned int OLDID_base_alt=1;
			map<unsigned int , unsigned int> indexes; 
			map<unsigned int , unsigned int>::iterator ite;
			vector <double*> distances;
			
			Tuple<Tuple<unsigned int, 2>,(TEMPLATE_INT_POWER<3,2>::ans) -1 > ncoor;int nbncoor;
			
			//	Vector < Tuple<double, 2> > point_to_check;
			
			while(tf.fetch(im_in)){	im = im_in.makeSlice(0,0);
				printf("Fetch Segemented image\n");fflush(stdout);
				mcc.load(mcvf);
				printf("Multicover loaded\n");fflush(stdout);
				tf_d.fetch(im_in);im_d = im_in.makeSlice(0,0);
				printf("Fetch Distance image\n");fflush(stdout);
				out_hidden.setSizes(im.dims); 
				ExOp::toZero(out_hidden);
				Vector<unsigned int> lab_to_group;
				
				if (grametapix) {DataGrid<unsigned int, 3> im_metix_in; tf_metix.fetch(im_metix_in); im_metix = im_metix_in.makeSlice(0,0);}
				printf("Watershed loaded\n"); fflush(stdout);
				DataGrid<unsigned int, 2> lab = im.LabelConnected( isToLabel_seg ,&nbrect,&labelrect);
				
				
				for(i = 0 ; i <mcc.group.size(); i++){
					printf("Process group no %i\n",i); fflush(stdout);// ExOp::show(mcc.group[i].rect);ExOp::show(lab.dims);
					for(j=0;j<mcc.group[i].cover.size();j++){
						mcc.group[i].cover[j].show();
						coors[0] = (unsigned int) (mcc.group[i].cover[j].center[0] +  mcc.group[i].rect[0]);
						coors[1] = (unsigned int) (mcc.group[i].cover[j].center[1] +  mcc.group[i].rect[1]);
						if (coors[0] > 0x80000000) coors[0]=0;
						if (coors[1] > 0x80000000) coors[1]=0;
						if (coors[0] >= lab.dims[0]) coors[0]=lab.dims[0]-1;
						if (coors[1] >= lab.dims[1]) coors[1]=lab.dims[1]-1;
						
						
						curlab = lab(coors);
						if (curlab != 0) break;
						nbncoor = lab.get_indirectNeightbor(coors, ncoor);
						for(k=0;k<nbncoor;k++) {ExOp::show(ncoor[k]);curlab = lab(ncoor[k]); if (curlab != 0) break;}
						if (curlab != 0) break;
					}
				//	printf("inner\n"); fflush(stdout);
					if (j<mcc.group[i].cover.size()) {
						
						if (mcc.group[i].cover.size() > 1){
							
							if (grametapix){
						//		printf("inner grame\n"); fflush(stdout);
								
								if  (da_max_is_marked){
									
									if (mcc.group[i].rect[0] < 0) exit(1);
									if (mcc.group[i].rect[1] < 0) exit(1);
									if (mcc.group[i].rect[0] + mcc.group[i].rect[2] > lab.dims[0]) mcc.group[i].rect[2] = lab.dims[0] - mcc.group[i].rect[0];
									if (mcc.group[i].rect[1] + mcc.group[i].rect[3] > lab.dims[1]) mcc.group[i].rect[3] = lab.dims[1] - mcc.group[i].rect[1];
									for( coors[0] = mcc.group[i].rect[0] ;  coors[0] < mcc.group[i].rect[0] + mcc.group[i].rect[2] ; coors[0]++)
										for( coors[1] = mcc.group[i].rect[1] ;  coors[1] < mcc.group[i].rect[1] + mcc.group[i].rect[3] ; coors[1]++){
											
											
											if (lab(coors) == curlab) { // incontiguous area
												metalab = im_metix(coors);
												if (metalab & 0x80000000){
													metalab = distances.size();
													distances.push_back( new double[mcc.group[i].cover.size()]);
													for(j = 0; j< mcc.group[i].cover.size(); j++){
														distances[metalab][j] = hypot( mcc.group[i].cover[j].center[0] + mcc.group[i].rect[0] - coors[0], mcc.group[i].cover[j].center[1] + mcc.group[i].rect[1] - coors[1]);
													}
													indexes[im_metix(coors) & 0x7FFFFFFF ] = metalab;
												}
											}
										}
									}else{
									
									if (mcc.group[i].rect[0] < 0) exit(1);
									if (mcc.group[i].rect[1] < 0) exit(1);
									if (mcc.group[i].rect[0] + mcc.group[i].rect[2] > lab.dims[0]) mcc.group[i].rect[2] = lab.dims[0] - mcc.group[i].rect[0];
									if (mcc.group[i].rect[1] + mcc.group[i].rect[3] > lab.dims[1]) mcc.group[i].rect[3] = lab.dims[1] - mcc.group[i].rect[1];
									
									for( coors[0] = mcc.group[i].rect[0] ;  coors[0] < mcc.group[i].rect[0] + mcc.group[i].rect[2] ; coors[0]++)
										for( coors[1] = mcc.group[i].rect[1] ;  coors[1] < mcc.group[i].rect[1] + mcc.group[i].rect[3] ; coors[1]++){
											if (lab(coors) == curlab) { // incontiguous area
												metalab = im_metix(coors);
												if ( (ite = indexes.find(metalab)) == indexes.end()){
													metalab = distances.size();
													distances.push_back( new double[mcc.group[i].cover.size()]);
													for(j = 0; j< mcc.group[i].cover.size(); j++){
														distances[metalab][j] = hypot( mcc.group[i].cover[j].center[0] + mcc.group[i].rect[0] - coors[0], mcc.group[i].cover[j].center[1] + mcc.group[i].rect[1] - coors[1]);
													}
													indexes[im_metix(coors)] = metalab;
												} else {metalab = ite->second;
													for(j = 0; j< mcc.group[i].cover.size(); j++){
														pix[0] = hypot( mcc.group[i].cover[j].center[0] + mcc.group[i].rect[0] - coors[0], mcc.group[i].cover[j].center[1] + mcc.group[i].rect[1] - coors[1]);
														if (pix[0]  < distances[metalab][j]) distances[metalab][j] = pix[0];
													}
												}
											}
										}}
								
								
								
								for(l=0;l< distances.size();l++) {
									pix[0] = distances[l][0] - 0.5f * mcc.group[i].cover[0].width ;
									coors[0] =0;
									//		printf("%f\t",pix[0] );
									for(j = 1; j< mcc.group[i].cover.size(); j++){
										//			printf("%f\t",distances[l][j] );
										pix[1] = distances[l][j] - 0.5f * mcc.group[i].cover[j].width;
										if (pix[1]< pix[0]){pix[0] =pix[1] ; coors[0] =j;}
									}
									distances[l][0] = (double) coors[0];
									//		printf("\n");
								}
								
								
								if (do_refine_area){
									// starts with largest , frees undesirables
									// force undesirable to next, then frees undesirables
									
									
									
									
									
									
									
								}
								
								// write to image the found hidden state
								if (ellipsedev) memset(ellipse_buf,'\0',sizeof(double)*6 * mcc.group[i].cover.size());
								
								for( coors[0] = mcc.group[i].rect[0] ;  coors[0] < mcc.group[i].rect[0] + mcc.group[i].rect[2] ; coors[0]++)
									for( coors[1] = mcc.group[i].rect[1] ;  coors[1] < mcc.group[i].rect[1] + mcc.group[i].rect[3] ; coors[1]++){
										if (lab(coors) == curlab) {
											if ( (ite = indexes.find(im_metix(coors) & 0x7FFFFFFF)) == indexes.end()){
												
											}else{
												metalab = ite->second;
												
												out_hidden(coors) = ID_base + ((unsigned int)distances[metalab][0]); 
												if (ellipsedev){
													pix[0] = im(coors);
													ellipse_buf[0 + ((unsigned int)distances[metalab][0])*6 ] += pix[0];
													ellipse_buf[1 + ((unsigned int)distances[metalab][0])*6 ] += pix[0] * coors[0];
													ellipse_buf[2 + ((unsigned int)distances[metalab][0])*6 ] += pix[0] * coors[1];
													ellipse_buf[3 + ((unsigned int)distances[metalab][0])*6 ] += pix[0] * coors[0]* coors[0];
													ellipse_buf[4 + ((unsigned int)distances[metalab][0])*6 ] += pix[0] * coors[1]* coors[1];
													ellipse_buf[5 + ((unsigned int)distances[metalab][0])*6 ] += pix[0] * coors[1]* coors[0];
												}
											}
										}
									}
								
								if (ellipsedev){
									
									
									
									for(l=0;l< mcc.group[i].cover.size();l++) {
										tmpcellpose.initFromCummul(ellipse_buf + l*6);
										pix[0] = hypot( mcc.group[i].cover[l].eccentric[0] - tmpcellpose.eccentric[0],mcc.group[i].cover[l].eccentric[1] - tmpcellpose.eccentric[1]);
										pix[1] = hypot( mcc.group[i].cover[l].eccentric[0] + tmpcellpose.eccentric[0],mcc.group[i].cover[l].eccentric[1] + tmpcellpose.eccentric[1]);
										fprintf(ell_dev, "%i\t%i\t%f\t%f\t%f\t%f\t%c\n", ID_base + l, mcc.group[i].cover.size()
												, hypot(mcc.group[i].rect[0] + mcc.group[i].cover[l].center[0] - tmpcellpose.center[0],mcc.group[i].rect[1] + mcc.group[i].cover[l].center[1] - tmpcellpose.center[1]) 
												, pix[0] < pix[1] ? pix[0] : pix[1]
												, log(mcc.group[i].cover[l].width) - log(tmpcellpose.width)
												, mcc.group[i].cover[l].cmpArea()
												, (mcc.group[i].cover[l].error[0] <=2.5f) ? ((mcc.group[i].cover[l].error[0] <= 0.5f) ? ((mcc.group[i].cover[l].error[0] <= -0.5f) ? 'x' : 'a') : ((mcc.group[i].cover[l].error[0] <= 1.5f) ? 'm' : 'b')) : ((mcc.group[i].cover[l].error[0] <=3.5f) ? 'd' : 'l')
												);
										
									}
									
								}
								
							}else{
								
							}
							
							for(l=0;l< distances.size();l++) delete[](distances[l]);
							distances.clear();
							indexes.clear();
						}else{
							
							printf("%i\t%i\n",mcc.group[i].rect[0] + mcc.group[i].rect[2],mcc.group[i].rect[1] + mcc.group[i].rect[3]);
							ExOp::show(out_hidden.dims);
							if (ellipsedev) {
								memset(ellipse_buf,'\0',sizeof(double)*6);
								for( coors[0] = mcc.group[i].rect[0] ;  coors[0] < mcc.group[i].rect[0] + mcc.group[i].rect[2] ; coors[0]++)
									for( coors[1] = mcc.group[i].rect[1] ;  coors[1] < mcc.group[i].rect[1] + mcc.group[i].rect[3] ; coors[1]++){
										if (lab(coors) == curlab){
											out_hidden(coors) = ID_base;
											pix[0] = im(coors);
											ellipse_buf[0] += pix[0];
											ellipse_buf[1] += pix[0] * coors[0];
											ellipse_buf[2] += pix[0] * coors[1];
											ellipse_buf[3] += pix[0] * coors[0]* coors[0];
											ellipse_buf[4] += pix[0] * coors[1]* coors[1];
											ellipse_buf[5] += pix[0] * coors[1]* coors[0];
										}	}
							}else{
								for( coors[0] = mcc.group[i].rect[0] ;  coors[0] < mcc.group[i].rect[0] + mcc.group[i].rect[2] ; coors[0]++)
									for( coors[1] = mcc.group[i].rect[1] ;  coors[1] < mcc.group[i].rect[1] + mcc.group[i].rect[3] ; coors[1]++){
										if (lab(coors) == curlab) out_hidden(coors) = ID_base;
									}	}	
							if (ellipsedev) {
								tmpcellpose.initFromCummul(ellipse_buf);
								pix[0] = hypot( mcc.group[i].cover[0].eccentric[0] - tmpcellpose.eccentric[0],mcc.group[i].cover[0].eccentric[1] - tmpcellpose.eccentric[1]);
								pix[1] = hypot( mcc.group[i].cover[0].eccentric[0] + tmpcellpose.eccentric[0],mcc.group[i].cover[0].eccentric[1] + tmpcellpose.eccentric[1]);
								l=0;
								fprintf(ell_dev, "%i\t%i\t%f\t%f\t%f\t%f\t%c\n", ID_base, mcc.group[i].cover.size()
										, hypot(mcc.group[i].rect[0] + mcc.group[i].cover[0].center[0] - tmpcellpose.center[0],mcc.group[i].rect[1] + mcc.group[i].cover[0].center[1] - tmpcellpose.center[1]) 
										, pix[0] < pix[1] ? pix[0] : pix[1]
										, log(mcc.group[i].cover[0].width) - log(tmpcellpose.width) 
										, mcc.group[i].cover[0].cmpArea()
										, (mcc.group[i].cover[0].error[0] <=2.5f) ? ((mcc.group[i].cover[0].error[0] <= 0.5f) ? ((mcc.group[i].cover[l].error[0] <= -0.5f) ? 'x' : 'a') : ((mcc.group[i].cover[l].error[0] <= 1.5f) ? 'm' : 'b')) : ((mcc.group[i].cover[l].error[0] <=3.5f) ? 'd' : 'l')
										);
							}								
							
						}
						
						
						
						
						
					}
					ID_base += mcc.group[i].cover.size();
				}
				
				
				
				
				
				
				
				if (out_cpartition){
					DataGrid<unsigned int, 3> cp_conv = out_hidden.fromSlice(0);
					tf_cp.put(cp_conv, (unsigned int ) 0, (unsigned int ) nbrect);
				}
				
				if (Z_mcv){
					
					{//:
						mcc2.load(mcvf_Z);
						
						map<unsigned int, unsigned int> layer_map;
						map<unsigned int, unsigned int>::iterator mite;
						
						
						
						
						
						Vector< pair< Tuple<double, 2>, pair<Tuple<unsigned int, 2> , char> > > lista;
						Vector< pair< Tuple<double, 2>, pair<Tuple<unsigned int, 2> , char> > > listb;
						
						Vector< KeyElem<double, Tuple<unsigned int, 2> > > list_or;
						
						pair<Tuple<double, 2> , pair< Tuple<unsigned int, 2> , char> > tmpitem;
						KeyElem<double,  Tuple<unsigned int, 2> > tmporitem;
						
						
						for(i=0;i<mcc.group.size();i++){
							tmpitem.second.first[0] = i;
							for(j=0;j<mcc.group[i].cover.size();j++){
								tmpitem.second.first[1] = j;
								tmpitem.first[0] = mcc.group[i].cover[j].center[0] + mcc.group[i].rect[0];
								tmpitem.first[1] = mcc.group[i].cover[j].center[1] + mcc.group[i].rect[1];
								lista.push_back(tmpitem);
							}
						}
						k=0;
						for(i=0;i<mcc2.group.size();i++){
							tmpitem.second.first[0] = i;
							
							for(j=0;j<mcc2.group[i].cover.size();j++){
								tmpitem.second.first[1] = j;
								tmpitem.first[0] = mcc2.group[i].cover[j].center[0] + mcc2.group[i].rect[0];
								tmpitem.first[1] = mcc2.group[i].cover[j].center[1] + mcc2.group[i].rect[1];
								listb.push_back(tmpitem);
								
								if (mcc2.group[i].cover[j].error[0] <= -0.5f) fprintf(txtf_Z,"%i\td\n", k+OLDID_base_alt);
								else if (mcc2.group[i].cover[j].error[0] <= 0.5f) fprintf(txtf_Z,"%i\ta\n", k+OLDID_base_alt);
								else if (mcc2.group[i].cover[j].error[0] <= 1.5f) fprintf(txtf_Z,"%i\tm\n", k+OLDID_base_alt);
								else if (mcc2.group[i].cover[j].error[0] <= 2.5f) fprintf(txtf_Z,"%i\tb\n", k+OLDID_base_alt);
								else if (mcc2.group[i].cover[j].error[0] <= 3.5f) fprintf(txtf_Z,"%i\td\n", k+OLDID_base_alt);
								else if (mcc2.group[i].cover[j].error[0] <= 4.5f) fprintf(txtf_Z,"%i\tl\n", k+OLDID_base_alt);
								k++;
							}
						}
						
						for(tmporitem.d[0]=0;tmporitem.d[0]<lista.size();tmporitem.d[0]++){
							for(tmporitem.d[1]=0;tmporitem.d[1]<listb.size();tmporitem.d[1]++){
								tmporitem.k = hypot( lista[tmporitem.d[0]].first[0] - listb[tmporitem.d[1]].first[0] , lista[tmporitem.d[0]].first[1] - listb[tmporitem.d[1]].first[1]);
								list_or.push_back(tmporitem);
							}
						}
						list_or.sort();
						
						for(i=0;i<list_or.size();i++){
							if (lista[list_or[i].d[0] ].first[0] <= -1.0f) continue;
							if (listb[list_or[i].d[1] ].first[0] <= -1.0f) continue;
							
							lista[list_or[i].d[0] ].first[0] = -1.0f; // mark matched
							listb[list_or[i].d[1] ].first[0] = -1.0f;
							
							
							layer_map[list_or[i].d[0] +OLDID_base] = list_or[i].d[1] +OLDID_base_alt;
							
						}
						
						for(i=0;i<mcc2.group.size();i++){
							for(j=0;j<mcc2.group[i].cover.size();j++) OLDID_base_alt++;
						}
						
						
						
						
						
						
						
						
						DataGrid<unsigned int, 2>::KeyIterator ite =  out_hidden.getKeyIterator();
						
						if (ite.first()) do{
							
							if ((mite= layer_map.find(out_hidden(ite()))) == layer_map.end()) out_hidden(ite()) = 0;
							else out_hidden(ite()) = mite->second;
							
						}while(ite.next());
						
						DataGrid<unsigned int, 3> cp_conv = out_hidden.fromSlice(0);
						tfZ.put(cp_conv, (unsigned int ) 0, (unsigned int ) nbrect);
					}//:
					
				}
				
				
				
				
				labelrect.clear();
				k++;
				OLDID_base = ID_base;
			}
			
			if (Z_mcv){
				fclose(mcvf_Z);
				fclose(txtf_Z);
			}
			
			
			if (ellipsedev) {fclose(ell_dev); delete[](ellipse_buf);}
			
			
		return 0;}
		void Taskscope<TASK_EXTRACT_HIDDENMAP_DIRECT>::help(){
			printf("Finds the circle cover in the segmented image\nmodeling the pixel intensity within cells to be different in the raw image.\n");
			printf("\n");
			printf("\t(in) Segmented Image File\n"); // foreground probability
			printf("\t(in) Distance  Image File\n"); // foreground probability
			printf("\t(in) .mcv file\n");
			printf("Flags\n\n");
			//					printf("\t-m (int):\tMaximum distance searched (default = 100)\n");
			
			printf("\t-D (filepath).tif: Distance to Othercell image\n");
			
			//		printf("\t-m (filepath).mcv: output multicoverfile\n");
			printf("\t-p (filepath).tif: Hidden Map Preview\n");
			printf("\t-P (filepath).tif: output clump partition\n");
			printf("\t-C (filepath).tif: output cell partition (unsigned int)\n");
			
			printf("\t-M: use Marked point to determine ownership\n");
			printf("\t-R : refine area to maximize density in ellipse\n");
			
			
			//	printf("\t-CC (filepath).tif: output cell partition confidence (double)\n");
			
			printf("\t-G (filepath).tif: use metapixel constrain\n");
			printf("\t-f (int x): filter out cells with nb of pixels <= x\n");
			printf("\t-d (filepath).txt: report deviation from input to output ellipse\n");
			
			printf("\t-Z (input filepath).mcv (output filepath).tif (type).txt: in parallell, do MCV intersection with current MCV, and make a hiddenmap file for input MCV\n");
			printf("Version 1.0\n");
		}

	
	


