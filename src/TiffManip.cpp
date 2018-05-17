/*
 * TiffManip.cpp
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


	Taskscope<TASK_TIFF_FILE_MANIPULATION>::Taskscope(): show(false), mode(0), val(1.0f),wantrgb(false),scale_factor(NULL),trans(NULL),use_fourrier(false), type_out('Z'), make_normalization_image(false) {ndim[0] =0; ndim[1] =0;
		clamprange[1] = 0.0f;clamprange[0] =1.0f;
	}
	
void Taskscope<TASK_TIFF_FILE_MANIPULATION>::nbaddtoken(char const * const token, int& min, int& max){
	switch(*token){
		case '\0': min =2; break;
		case 's': min =0; break;
		case 'r': min =0; break;
		case 'f': min =1; break;
		case 'c': min =1; break;
		case 't': min=1;break;
		case 'T': max=1;min=0;break;
		case 'C': min =2; break;
		case 'R': min =2; break;
		case 'N': min =0; break;
		case 'F': min =0; break;
		case 'S': min =1; break;
	}
}

void Taskscope<TASK_TIFF_FILE_MANIPULATION>::store(char* const * token, int nbtoken){ 
	
	switch(token[0][1]){
		case 's':show = true;
			break;
		case 'f': mode =1;
			readList(token[1], frame);
			break;
		case 'c': mode =1;
			readList(token[1], channels);
			break;
		case 't': 
			trans = token[1];
			break;
		case 'r':wantrgb = true;
			break;
		case 'C': clamprange[0] = atof(token[1]);clamprange[1] = atof(token[2]);break;
		case 'T': 
			type_out = token[0][2] == '\0' ? token[1][0] : token[0][2];
			break;
		case 'R': 
			ndim[0] = atoi(token[1]);
			ndim[1] = atoi(token[2]);
			break;
		case 'N':
			make_normalization_image = true; break;
			//		case 'm':
		case 'S': scale_factor = token[1]; break;
		case 'F': use_fourrier = true; break;
	}
	
}
	
/*
int Taskscope<TASK_TIFF_FILE_MANIPULATION>::defstore(char* const * token, int nbtoken){ 
	file_in = token[0];
	file_out = token[1];
	int i,j,l,k;
	char buffer[65536];
	
	double pix[32];
	double pox[32];
	vector<Madstructs::Image<TEMPLATE_IMAGE_TYPE> *> input;
	strcpy(buffer, file_in);
	Madstructs::Image<TEMPLATE_IMAGE_TYPE>::LoadTiffImage(input, buffer);
	
	Madstructs::Image<TEMPLATE_IMAGE_TYPE>* tmptmp;
	
	switch(mode){
		case 1:{
			
			vector<Madstructs::Image<TEMPLATE_IMAGE_TYPE> *> output;
			
			for(i=0;i< frame.size();i++) {
				if ((frame[i] == 0)||(frame[i] > input.size())) fprintf(stderr, "Warning! frame %i does not exist!\n", frame[i]);
				else {
					if (!wantrgb){
						if (trans){ 
							tmptmp = new Madstructs::Image<TEMPLATE_IMAGE_TYPE>();
							tmptmp->sizex = input[frame[i]-1]->sizex;
							tmptmp->sizey = input[frame[i]-1]->sizey;
							
							if (channels.size() == 0) {
							tmptmp->channels = input[frame[i]-1]->channels;
							tmptmp->allocateBuffer();
							tmptmp->initBlack();
							memset(pix,'\0',sizeof(double)*32);
							for(l=0;l<tmptmp->sizex ;l++){
								for(j=0;j<tmptmp->sizey ;j++){
									input[frame[i]-1]->getPixel(l,j,pix);
									if (pix[0] < 1.0f) pix[0] =0.0f;
									else if (pix[0] > 2.5f) pix[0] =0.0f;
									else pix[0] =1.0f;
									tmptmp->setPixel(l,j,pix);
								}}
							
							}else{
								tmptmp->channels = channels.size();
								tmptmp->allocateBuffer();
								tmptmp->initBlack();
								for(l=0;l<tmptmp->sizex ;l++){
									for(j=0;j<tmptmp->sizey ;j++){
										input[frame[i]-1]->getPixel(l,j,pix);
										for(k=0;k<channels.size();k++) pox[k] = pix[channels[k]];
										tmptmp->setPixel(l,j,pox);
								}}
							}
							output.push_back(tmptmp);
						}else output.push_back(input[frame[i]-1]);
						
					}else{
						tmptmp = new Madstructs::Image<TEMPLATE_IMAGE_TYPE>();
						tmptmp->sizex = input[frame[i]-1]->sizex;
						tmptmp->sizey = input[frame[i]-1]->sizey;
						tmptmp->channels =3 ;
						tmptmp->allocateBuffer();
						tmptmp->initBlack();
						memset(pix,'\0',sizeof(double)*32);
						for(l=0;l<tmptmp->sizex ;l++){
							for(j=0;j<tmptmp->sizey ;j++){
								input[frame[i]-1]->getPixel(l,j,pix);
								if (trans){
									if (pix[0] < 1.0f) pix[0] =0.0f;
									else if (pix[0] > 2.5f) pix[0] =0.0f;
									else pix[0] =1.0f;
								}
								tmptmp->setPixel(l,j,pix);
							}}
						output.push_back(tmptmp);
					}
					
				}
				
			}
			strcpy(buffer, file_out);
			Madstructs::Image<TEMPLATE_IMAGE_TYPE>::SaveTiffImage(output, buffer);
		}break;
			
		default:
			
		strcpy(buffer, file_out);
		Madstructs::Image<TEMPLATE_IMAGE_TYPE>::SaveTiffImage(input, buffer);
			
			
			
	}
	
	
	if (show){
		strcpy( buffer, "open -a Preview \"");
		strcpy(buffer + 17 , file_out);
		strcpy( buffer + 17 + strlen(file_out), "\"");
		system(buffer);
	}
	
	
	
return 0;}*/

	int Taskscope<TASK_TIFF_FILE_MANIPULATION>::defstore(char* const * token, int nbtoken){ 
		file_in = token[0];
		file_out = token[1];
		unsigned int i,j,l,k;
		char buffer[65536];
		double pix[32];
		double pox[32];
		TiffFile tfi(file_in);
		TiffFile tfo(file_out, true);
		DataGrid<double, 3> image_in;
		
		DataGrid<double, 3> image_out;
		Tuple<unsigned int,3> dims, coor,acoor;
		DataGrid<double, 3> tmp_image;
		Vector< DataGrid<double, 3> > image_vector;
		char imageType;
		unsigned int curfr;
		if (frame.size() != 0) image_vector.setSize(frame.size()); 
		
		for(unsigned int frno=1;tfi.fetch(image_in, &imageType);frno++){
			if ((frno == 1)&&(type_out == 'Z')) type_out = imageType;
			if (frame.size() != 0){
				for(curfr=0;curfr<frame.size() ;curfr++) if (frame[curfr] == frno) break;
				if (curfr == frame.size()) continue;
			}else{curfr = image_vector.size(); image_vector.push_back();}
			//printf("process!\n");
			dims[1] = image_in.dims[1];
			dims[2] = image_in.dims[2];
			if (channels.size() != 0){
				for(k=0;k<channels.size();k++) if (channels[k] >= image_in.dims[0]) {printf("Image does not have a %ith channel!\n", channels[k]); exit(1); }
				dims[0] = channels.size();
			}else dims[0] = image_in.dims[0];
			
			image_vector[curfr].setSizes(dims);
			
			if (channels.size() != 0){
			for(coor[2]=0;coor[2] < dims[2];coor[2]++) for(coor[1]=0;coor[1] < dims[1];coor[1]++){
					acoor = coor;
					for(acoor[0]=0; acoor[0]< channels.size();acoor[0]++) {coor[0] = channels[acoor[0]]; 
						if (clamprange[1] < clamprange[0]) image_vector[curfr](acoor) = image_in(coor);
						else if (clamprange[1]< image_in(coor)) image_vector[curfr](acoor) =clamprange[1];
						else if (clamprange[0]> image_in(coor)) image_vector[curfr](acoor) =clamprange[0];
						else image_vector[curfr](acoor) = image_in(coor);
					}
			}
			} else image_vector[curfr] = image_in;
			
			if (ndim[0] != 0){
				coor[0] = image_vector[curfr].dims[0];
				coor[1] = ndim[0];
				coor[2] = ndim[1];
				if (use_fourrier) image_vector[curfr].toresize_crude(coor);
				else image_vector[curfr].toresize(coor);
			}
		}
		
		TiffFile tf_scale(scale_factor);
		
		for(k=0;k<image_vector.size();k++) {
			if (make_normalization_image){
				double* pixbuffer = new double[image_vector[k].dims[0]];
				for(coor[0]=0;coor[0] < image_vector[k].dims[0];coor[0]++) pixbuffer[coor[0]] =0;
				for(coor[2]=0;coor[2] < image_vector[k].dims[2];coor[2]++) for(coor[1]=0;coor[1] < image_vector[k].dims[1];coor[1]++) for(coor[0]=0;coor[0] < image_vector[k].dims[0];coor[0]++){
					pixbuffer[coor[0]] += image_vector[k](coor);
				}
				for(coor[0]=0;coor[0] < image_vector[k].dims[0];coor[0]++) pixbuffer[coor[0]] /= image_vector[k].dims[2] * image_vector[k].dims[1];
				
				for(coor[2]=0;coor[2] < image_vector[k].dims[2];coor[2]++) for(coor[1]=0;coor[1] < image_vector[k].dims[1];coor[1]++) for(coor[0]=0;coor[0] < image_vector[k].dims[0];coor[0]++){
					image_vector[k](coor) = pixbuffer[coor[0]] / image_vector[k](coor);
				}
				delete[](pixbuffer);
			}
			
			if (scale_factor){
				if (!tf_scale.fetch(tmp_image)) {fprintf(stderr, "Missing frame in tiff file %s\n", scale_factor); exit(1);}
				if (use_fourrier) {tmp_image.toresize_crude(image_vector[k].dims);}
				else {tmp_image.toresize(image_vector[k].dims);}
				if ((type_out =='C')||(type_out =='S')||(type_out =='I')) {
					for(coor[2]=0;coor[2] < tmp_image.dims[2];coor[2]++) for(coor[1]=0;coor[1] < tmp_image.dims[1];coor[1]++) for(coor[0]=0;coor[0] < tmp_image.dims[0];coor[0]++) if (tmp_image(coor) < 0.0f) tmp_image(coor) = 0.0f;
				}
				image_vector[k] *= tmp_image;
			}
			
			
			switch(type_out){
				case 'c': tfo.put( image_vector[k], (char) -128,(char) 127); break;
				case 'C': tfo.put( image_vector[k], (unsigned char) 0,(unsigned char) 255); break;
				case 's': tfo.put( image_vector[k], (short) 0x8000,(short) 0x7FFF); break;
				case 'S': tfo.put( image_vector[k], (unsigned short) 0,(unsigned short) 65535); break;
				case 'i': tfo.put( image_vector[k], (int) 0x80000000,(int) 0x7FFFFFFF); break;
				case 'I': tfo.put( image_vector[k], (unsigned int) 0,(unsigned int) 0xFFFFFFFF); break;
				case 'f': tfo.put( image_vector[k], (float) 0.0f,(float) 1.0f); break;
				case 'd': tfo.put( image_vector[k], (double) 0.0f,(double) 1.0f); break;
			}
		}
		
		
		/*
		
		switch(mode){
			case 1:{
				
				vector<Madstructs::Image<TEMPLATE_IMAGE_TYPE> *> output;
				
				for(i=0;i< frame.size();i++) {
					if ((frame[i] == 0)||(frame[i] > input.size())) fprintf(stderr, "Warning! frame %i does not exist!\n", frame[i]);
					else {
						if (!wantrgb){
							if (trans){ 
								tmptmp = new Madstructs::Image<TEMPLATE_IMAGE_TYPE>();
								tmptmp->sizex = input[frame[i]-1]->sizex;
								tmptmp->sizey = input[frame[i]-1]->sizey;
								
								if (channels.size() == 0) {
									tmptmp->channels = input[frame[i]-1]->channels;
									tmptmp->allocateBuffer();
									tmptmp->initBlack();
									memset(pix,'\0',sizeof(double)*32);
									for(l=0;l<tmptmp->sizex ;l++){
										for(j=0;j<tmptmp->sizey ;j++){
											input[frame[i]-1]->getPixel(l,j,pix);
											if (pix[0] < 1.0f) pix[0] =0.0f;
											else if (pix[0] > 2.5f) pix[0] =0.0f;
											else pix[0] =1.0f;
											tmptmp->setPixel(l,j,pix);
										}}
									
								}else{
									tmptmp->channels = channels.size();
									tmptmp->allocateBuffer();
									tmptmp->initBlack();
									for(l=0;l<tmptmp->sizex ;l++){
										for(j=0;j<tmptmp->sizey ;j++){
											input[frame[i]-1]->getPixel(l,j,pix);
											for(k=0;k<channels.size();k++) pox[k] = pix[channels[k]];
											tmptmp->setPixel(l,j,pox);
										}}
								}
								output.push_back(tmptmp);
							}else output.push_back(input[frame[i]-1]);
							
						}else{
							tmptmp = new Madstructs::Image<TEMPLATE_IMAGE_TYPE>();
							tmptmp->sizex = input[frame[i]-1]->sizex;
							tmptmp->sizey = input[frame[i]-1]->sizey;
							tmptmp->channels =3 ;
							tmptmp->allocateBuffer();
							tmptmp->initBlack();
							memset(pix,'\0',sizeof(double)*32);
							for(l=0;l<tmptmp->sizex ;l++){
								for(j=0;j<tmptmp->sizey ;j++){
									input[frame[i]-1]->getPixel(l,j,pix);
									if (trans){
										if (pix[0] < 1.0f) pix[0] =0.0f;
										else if (pix[0] > 2.5f) pix[0] =0.0f;
										else pix[0] =1.0f;
									}
									tmptmp->setPixel(l,j,pix);
								}}
							output.push_back(tmptmp);
						}
						
					}
					
				}
				strcpy(buffer, file_out);
				Madstructs::Image<TEMPLATE_IMAGE_TYPE>::SaveTiffImage(output, buffer);
			}break;
				
			default:
				
				strcpy(buffer, file_out);
				Madstructs::Image<TEMPLATE_IMAGE_TYPE>::SaveTiffImage(input, buffer);
				
				
				
		}
		*/
		
		if (show){
			strcpy( buffer, "open -a Preview \"");
			strcpy(buffer + 17 , file_out);
			strcpy( buffer + 17 + strlen(file_out), "\"");
			system(buffer);
		}
		
		
		
	return 0;}
	
void Taskscope<TASK_TIFF_FILE_MANIPULATION>::help(){
	printf("Makes an Tiff file from a Tiff file, perfoming some basic operation.\n");
	printf("Arguments: [1-2]\n");
	printf("	(file) input tif image (compressed format disallowed)\n");
	printf("	(file) output tif image\n");
	printf("\n");
	printf("Flags\n\n");
	printf("\t-s:\tShow output in Preview (works on MacOS)\n");
	printf("\t-f x,x,...,x :\tFrame selection,\n");
	printf("\t-c x,x,...,x :\tChannel selection,\n");
	printf("Operation Flags\n\n");
	printf("\t-S (tiff File): Multiply by image, input image is rescaled if it does not match input file\n");

	printf("\t-N :\t make Normalization image, that is, if I is the Input, and O is the Output, if K is the average value of pixels in I, then K = O(x,y) * I(x,y) for all pixels\n");
	
	printf("\t-R (int width) (int height): Resize images\n");

	printf("Special Flags\n\n");
	printf("\t-F: Use Fourier Transform for any resizing (slow)\n");
//	printf("\t-t (string): applies transformation, (uses X,R,G,B variables)\n");
//	printf("\t-r : Forces output format to be RGB.\n");
	printf("\t-T (char) : Type for Output image {c,C,s,S,i,I,f,d}.\n");
	printf("Version 1.0\n");
}


	
	Taskscope<TASK_REMOVE_FACTORS>::Taskscope(){}
	
	void Taskscope<TASK_REMOVE_FACTORS>::nbaddtoken(char const * const token, int& min, int& max){
		switch(*token){
			case '\0': min =1; break;
			case 's': min =0; break;
			case 'r': min =0; break;
			case 'f': min =1; break;
			case 't': min=1;break;
		}
	}
	
	void Taskscope<TASK_REMOVE_FACTORS>::store(char* const * token, int nbtoken){ 
		

	}
	int Taskscope<TASK_REMOVE_FACTORS>::defstore(char* const * token, int nbtoken){ 
		char buffer[65536];
		
		Vector<DataGrid<double, 3>* > data;
		Vector<double> dafacts;
		DataGrid<double, 3> tmp;
		{//:
		TiffFile tfin(token[0]);
		while(tfin.fetch(tmp)){
			data.push_back(new DataGrid<double, 3>(tmp));
			}
		}//:
		float tmpf;
		FILE* f = fopen(token[0], "r+");
		
		while ( 1 == fscanf(f,"%c", buffer)){
			if (buffer[0] != 'F') continue;
			if (1 != fscanf(f,"%c", buffer)) break;
			if (buffer[0] != 'a') continue;
			if (1 != fscanf(f,"%c", buffer)) break;
			if (buffer[0] != 'c') continue;
			if (1 != fscanf(f,"%c", buffer)) break;
			if (buffer[0] != 't') continue;
			if (1 != fscanf(f,"%c", buffer)) break;
			if (buffer[0] != 'o') continue;
			if (1 != fscanf(f,"%c", buffer)) break;
			if (buffer[0] != 'r') continue;
			if (1 != fscanf(f,"%c", buffer)) break;
			if (buffer[0] != '=') continue;
			if (1 != fscanf(f,"%c", buffer)) break;
			if (buffer[0] != '"') continue;
			if (1 != fscanf(f,"%f", &tmpf)) break;
			dafacts.push_back(tmpf);
			if (data.size() == dafacts.size()) break;
		}
		fclose(f);
		printf("%i frames\t %i factors\n", data.size(), dafacts.size());
		
		strcpy(buffer, token[0]);
		unsigned int i = strlen(token[0])-1;
		while(buffer[i] != '.') i--;
		memcpy(buffer +i, ".tif", sizeof(char)*5);
		
		TiffFile tfout(buffer);
		for(i=0;i<data.size();i++){
			(*data[i]) *= dafacts[i];
			tfout.put( (*data[i]), (unsigned short) 0, (unsigned short) 65535);
			}
		
		
		
	return 0;}
	
	void Taskscope<TASK_REMOVE_FACTORS>::help(){
		printf("Makes an Tiff file from Flex file, uses 16-bit after factor extraction.\n");
		printf("   Name of the file *.flex");
		printf("\n");
		printf("Flags\n\n");
		printf("Version 1.0\n");
	}
	







Taskscope<TASK_TIFF_FILE_OPERATION>::Taskscope(): show(false), mode(0),file_dev(NULL),clamp(false),type_out('C'),oper_ation('a'),concat_file(NULL) {}
void Taskscope<TASK_TIFF_FILE_OPERATION>::nbaddtoken(char const * const token, int& min, int& max){
	switch(*token){
		case '\0': min =1; max =2; break;
		case 's': min =0; break;
		case 'a': min =0; break;
		case 'd': min =1; break;
		case 'c': min =2; break;
		case 'T': min =0; max =1; break;
		case 'o': min =0; max =1; break;
		case 'C': min =1; break;
	}
}
void Taskscope<TASK_TIFF_FILE_OPERATION>::store(char* const * token, int nbtoken){ 
	
	switch(token[0][1]){
		case 's':show = true;
			break;
		case 'd': file_dev = token[1];break;
		case 'C': concat_file = token[1];break;
		case 'c': 
			clamprange[0] = atof(token[1]);
			clamprange[1] = atof(token[2]);
			clamp = true;
			break;
		case 'T': 
			type_out = token[0][2] == '\0' ? token[1][0] : token[0][2];
			break;
		case 'o':
			oper_ation = token[0][2] == '\0' ? token[1][0] : token[0][2];
			break;
	}
	
}
int Taskscope<TASK_TIFF_FILE_OPERATION>::defstore(char* const * token, int nbtoken){ 
	file_out = token[0];
	file_in = (nbtoken == 2) ? token[1] : NULL;
	int i;
	char buffer[65536];
	
	//		char* fakeargs[] = {"./PureMadness","-s","-d", "../../../../Average.tif", "../../../../tmplist", "../../../../autmptmp11.tif"};
	//		mmm(6,fakeargs);
	
	
	FILE *flist = (file_in) ? fopen(file_in,"r+") : stdin;
	TiffFile tfout(file_out, true);
	int scfout;
	
	TiffFile concat_tf(concat_file);
	Tuple<unsigned int, 3> coor;
	vector< DataGrid< double ,3>* > images_final;
	
	switch(oper_ation){
		case 'a':{//:
			vector< DataGrid< WeightElem<double,1> ,3>* > images;
			
		};break ; //:
		case 'v':{//:
			vector< DataGrid< WeightElem<double,2> ,3>* > images;
			
		};break ; //:
		case 's':
		case 'm':
		case 'M':{//:
			
			
		};break ; //:
	}
	
	unsigned int nbimages =0;
	DataGrid<double,3> dainim;
	
	for(scfout = fscanf(flist, "%[^;\t\n ]", buffer); (scfout!= -1);scfout = fscanf(flist,"%*[;\t\n ]%[^;\t\n ]", buffer)){
		printf("Processing %s\n", buffer);fflush(stdout);
		nbimages++;
		TiffFile tfin(buffer);
		i=0xFFFFFFFF;

		char imageType;
		for(i=0;tfin.fetch(dainim, &imageType);i++){
			
			if (i == images_final.size()){
				images_final.push_back(new DataGrid<double,3>());
				images_final[i]->setSizes(dainim.dims);
				ExOp::toZero(* images_final[i]);
			}else{
				if (((*images_final[i]).dims[2] < dainim.dims[2])||((*images_final[i]).dims[1] < dainim.dims[1])||((*images_final[i]).dims[0] < dainim.dims[0])) {
					coor[2] = ((*images_final[i]).dims[2] < dainim.dims[2]) ? dainim.dims[2] : (*images_final[i]).dims[2];
					coor[1] = ((*images_final[i]).dims[1] < dainim.dims[1]) ? dainim.dims[1] : (*images_final[i]).dims[1];
					coor[0] = ((*images_final[i]).dims[0] < dainim.dims[0]) ? dainim.dims[0] : (*images_final[i]).dims[0];
					(*images_final[i]).toresize_crude(&(coor[0]));
				}
			}
			
			dainim.toresize_crude((*images_final[i]).dims);
			if (clamp){
				for(coor[2] =0;coor[2]< dainim.dims[2]; coor[2]++)
					for(coor[1] =0;coor[1]< dainim.dims[1]; coor[1]++)
						for(coor[0] =0;coor[0]< dainim.dims[0]; coor[0]++){
							if (dainim(coor) < clamprange[0]) dainim(coor) = clamprange[0];
							else if (dainim(coor) > clamprange[1]) dainim(coor) = clamprange[1];
						}
			}
			
			
			(*images_final[i]) += dainim;
			
			if (concat_file){
				switch(imageType){
					case 'c': concat_tf.put( dainim, (char) -128,(char) 127); break;
					case 'C': concat_tf.put( dainim, (unsigned char) 0,(unsigned char) 255); break;
					case 's': concat_tf.put( dainim, (short) 0x8000,(short) 0x7FFF); break;
					case 'S': concat_tf.put( dainim, (unsigned short) 0,(unsigned short) 65535); break;
					case 'i': concat_tf.put( dainim, (int) 0x80000000,(int) 0x7FFFFFFF); break;
					case 'I': concat_tf.put( dainim, (unsigned int) 0,(unsigned int) 0xFFFFFFFF); break;
					case 'f': concat_tf.put( dainim, (float) 0.0f,(float) 1.0f); break;
					case 'd': concat_tf.put( dainim, (double) 0.0f,(double) 1.0f); break;
				}
			}
		}
	}
	
	for(i=0;i<images_final.size();i++) (*images_final[i]) /= nbimages;
	
	for(i=0;i<images_final.size();i++){
		printf("Save frame %i, ", i); ExOp::show((*images_final[i]).dims);
		switch(type_out){
			case 'c': tfout.put( *(images_final[i]), (char) -128,(char) 127); break;
			case 'C': tfout.put( *(images_final[i]), (unsigned char) 0,(unsigned char) 255); break;
			case 's': tfout.put( *(images_final[i]), (short) 0x8000,(short) 0x7FFF); break;
			case 'S': tfout.put( *(images_final[i]), (unsigned short) 0,(unsigned short) 65535); break;
			case 'i': tfout.put( *(images_final[i]), (int) 0x80000000,(int) 0x7FFFFFFF); break;
			case 'I': tfout.put( *(images_final[i]), (unsigned int) 0,(unsigned int) 0xFFFFFFFF); break;
			case 'f': tfout.put( *(images_final[i]), (float) 0.0f,(float) 1.0f); break;
			case 'd': tfout.put( *(images_final[i]), (double) 0.0f,(double) 1.0f); break;
		}
	}
	
	
	
	exit(0);/*
	 for(i= strlen(file_in); buffer[i] != '/' ;i--) ;
	 
	 buffer[i+1] = '\0';
	 path.setChunk(0, buffer);
	 
	 
	 
	 
	 vector<Madstructs::Image<unsigned short> *> input;
	 
	 vector<Madstructs::Image<unsigned short> *> output;
	 vector<Madstructs::Image<unsigned short> *> output2;
	 FILE *f = fopen(file_in,"r+");
	 fscanf(f,"%[^\n]\n",buffer);
	 
	 path.setChunk(1, buffer);
	 Madstructs::Image<unsigned short>::LoadTiffImage(input, path());
	 printf("%s\n",path() );
	 
	 vector<DataGrid<Tuple<WeightElem<double, 2>, 1 > ,2> > compile;
	 
	 for(k=0;k<input.size();k++){
	 compile.push_back(DataGrid<Tuple<double, 1>,2>() );
	 }
	 unsigned int dims[2];
	 unsigned int coor[2];
	 Tuple<WeightElem<double, 2> ,1> in;
	 double pix[32];
	 for(k=0;k<input.size();k++){
	 dims[0] = input[k]->sizex;
	 dims[1] = input[k]->sizey;
	 compile[k].setSizes(dims);
	 
	 for(coor[1]=0;coor[1]<dims[1];coor[1]++) for(coor[0]=0;coor[0]<dims[0];coor[0]++){
	 input[k]->getPixel(coor[0],coor[1],pix);
	 if (clamp){
	 if (pix[0] < clamprange[0]) pix[0] = clamprange[0];
	 if (pix[0] > clamprange[1]) pix[0] = clamprange[1];
	 }
	 in[0] = WeightElem<double, 2 >(pix[0]);
	 compile[k](coor) = in;
	 }
	 
	 delete(input[k]);
	 }
	 printf("%i\n", compile.size());
	 input.clear();
	 
	 while(fscanf(f,"%[^\n]\n",buffer)){
	 path.setChunk(1, buffer);
	 Madstructs::Image<unsigned short>::LoadTiffImage(input, path());
	 printf("%s\n",path() );
	 for(k=0;k<input.size();k++){
	 
	 for(coor[1]=0;coor[1]<dims[1];coor[1]++) for(coor[0]=0;coor[0]<dims[0];coor[0]++){
	 input[k]->getPixel(coor[0],coor[1],pix);
	 in[0] = WeightElem<double, 2 >(pix[0]);
	 compile[k](coor) += in;
	 }
	 
	 delete(input[k]);
	 }
	 input.clear();
	 
	 if (feof(f)) break;
	 }
	 
	 DataGrid<Tuple<double, 1> ,2> tmptmp[2];
	 
	 for(k=0;k<compile.size();k++){
	 tmptmp[0].setSizes(compile[k].dims);				
	 if (file_dev){
	 tmptmp[1].setSizes(compile[k].dims);
	 tmptmp[1](WeightElem<double, 2>::OpStd(), compile[k] );
	 output2.push_back(new Madstructs::Image<unsigned short>(tmptmp[1]));	
	 }
	 tmptmp[0](WeightElem<double, 2>::OpMean(), compile[k] );
	 output.push_back(new Madstructs::Image<unsigned short>(tmptmp[0]));			
	 }
	 
	 strcpy(buffer, file_out);
	 Madstructs::Image<unsigned short>::SaveTiffImage(output, buffer);
	 
	 if (file_dev){
	 strcpy(buffer, file_dev);
	 Madstructs::Image<unsigned short>::SaveTiffImage(output2, buffer);
	 }
	 */
	
	if (show){
		strcpy( buffer, "open -a ImageJ \"");
		strcpy(buffer + 17 , file_out);
		strcpy( buffer + 17 + strlen(file_out), "\"");
		system(buffer);
	}

	
	
	
return 0;}
void Taskscope<TASK_TIFF_FILE_OPERATION>::help(){
	printf("Makes an Tiff file from a list of Tiff files (stored in a file), performs average.\n");
	printf("Arguments: [1-2]\n");
	printf("	(file) output image\n");
	printf("	(file) list of paths (stdin if absent)\n");
	printf("\n");
	printf("Flags\n\n");
	printf("\t-s:\tShow output in Preview (works on MacOS)\n");
	printf("\t-d (filepath): \tOutput std dev Image\n");
	printf("\t-c (float) (float): \t Any read pixel value is clamped in the range\n");
	printf("\t-C (FILE): Concatenate into a single tiff file\n");
	printf("\t-T(char) or -T (char): \t type of image for output (c C s S i I l L f d D t-z), default = 'C'\n");
	printf("\t-o(char) or -o (char): \t operation: ('a' average(default), 's' sum, 'm' minimum, 'M' maximum, 'v' variance)\n");
	
	printf("Version 1.0\n");
	
}



