/*
 * Madstructs_t.h
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


// cell segmentation stuff unsigned char
template<int flag>
double Classifier<flag>::getPriorMean(int inp,int state){
	if (inp == 0){
		if (state == 0) return(4.0f/255.0f);
		if (state == 1) return(64.0f/255.0f);
		if (state == 2) return(68.0/255.0f);
	}else{
		if (state == 0) return(64.0f/255.0f);
		if (state == 1) return(64.0f/255.0f);
		if (state == 2) return(64.0f/255.0f);
	}
	return(0.0f);
}

template<int flag>
double Classifier<flag>::getPriorStd(int inp,int state){
	if (inp == 0){
		if (state == 0) return(3.0 /255.0f);
		if (state == 1) return(6.0f/255.0f);
		if (state == 2) return(6.0f/255.0f);
	}else{
		if (state == 0) return(3.0f/255.0f);
		if (state == 1) return(3.0f/255.0f);
		if (state == 2) return(3.0f/255.0f);
	}
	return(0.0f);
	
}


template<int flag>
void Classifier<flag>::applyConstraints(){
	// state 0 to 2 is forbidden!
	 double tmp;	

	/*	transit[2] = 0.0f;
	 transit[6] = 0.0f;

	 if (means[1] < means[2]) { tmp = means[1]; means[1] = means[2];means[2] = tmp;}
	 if (means[1] < means[0]) { tmp = means[1]; means[1] = means[0];means[0] = tmp;}
	 if (means[2] < means[0]) { tmp = means[2]; means[2] = means[0];means[0] = tmp;}
	 if (boundspr[0] < boundspr[2]) {means[0]*=2; means[2] *=2;}*/
}


template<int flag>
void Classifier<flag>::showparam(){
	int z,w;
	for(w=0;w<inputsize;w++){
		for(z=0;z<nbstates;z++){
			printf("%i,%i: mean(%f) std(%f)\n",w,z,means[z+w*nbstates],stds[z+w*nbstates]);
		}
	}
	for(w=0;w<nbstates;w++){
		printf("%f%c",counts[w], w == nbstates-1 ? '\n' : '\t');
	}
	for(w=0;w<nbstates;w++){
		printf("%f%c",boundspr[w], w == nbstates-1 ? '\n' : '\t');
	}
	printf("\n\n");
	for(z=0;z<nbstates;z++){
		for(w=0;w<nbstates;w++){
			printf("%f\t",transit[w + z *nbstates]);
		}
		printf("\n");
	}
	printf("\n");
}

template<int flag>
Classifier<flag>::Classifier(int n_nbstate,int n_nbinput) : nbstates(n_nbstate), inputsize(n_nbinput), means(new double[n_nbstate*n_nbinput]), stds(new double[n_nbstate*n_nbinput]), transit(new double[n_nbstate*n_nbstate]),boundspr(new double[n_nbstate]), counts(new double[n_nbstate]){
	int w,z;
	for(w=0;w<inputsize;w++)
		for(z=0;z<nbstates;z++) {means[z+w*nbstates] = getPriorMean(w,z); stds[z+w*nbstates] =getPriorStd(w,z);}
	
}

template<int flag>
void Classifier<flag>::permState(int a, int b){
	if (a == b) return;
	
	double tmp;
	int w;
	for(w=0;w<inputsize;w++){
		tmp = means[a+w*nbstates];
		means[a+w*nbstates] = means[b+w*nbstates];
		means[b+w*nbstates] = tmp;
		tmp = stds[a+w*nbstates];
		stds[a+w*nbstates] = stds[b+w*nbstates];
		stds[b+w*nbstates] = tmp;
	}
	tmp = counts[a];
	counts[a] = counts[b];
	counts[b] = tmp;
	tmp = boundspr[a];
	boundspr[a] = boundspr[b];
	boundspr[b] = tmp;
	
	for(w=0;w<nbstates;w++){
		tmp = transit[a+w*nbstates];
		transit[a+w*nbstates] = transit[b+w*nbstates];
		transit[b+w*nbstates] = tmp;
	}
	for(w=0;w<nbstates;w++){
		tmp = transit[w+a*nbstates];
		transit[w+a*nbstates] = transit[w+b*nbstates];
		transit[w+b*nbstates] = tmp;
	}
}

template<int flag>
template<class T>
Image<float>* Classifier<flag>::initializeOutput(vector<Image<T>*> &input){
	Image<float>* outi = new Image<float>();
	outi->sizex = input[0]->sizex;
	outi->sizey = input[0]->sizey;
	outi->channels = nbstates;
	outi->allocateBuffer();
	return(outi);
}

template<int flag>
template<class T>
void Classifier<flag>::applyClassifier(vector<Image<T>* > &input,Image<float>* outi){
	
	double buffer[32];
	double tmpbuf[32];
	double ibuffer[32];
	double tmp;
	int x,y,z,w;
	
//		double prior[2];
//	prior[0] = 0.033f;
//	prior[1] = 1 - prior[0];
	
	for(w=0;w<inputsize;w++){
		for(z=0;z<outi->channels;z++){
			if (stds[z+w*nbstates] <= 0) stds[z+w*nbstates] =0.00000000001f;
		}
	}
	
	for(x=0;x<outi->sizex;x++){
		for(y=0;y<outi->sizey;y++){
			
			for(w=0;w<inputsize;w++){
				input[w]->getPixel(x,y,ibuffer);
				for(z=0;z<outi->channels;z++){
					if (stds[z+w*nbstates] > 0){
						tmpbuf[z] = (ibuffer[0] - means[z+w*nbstates]);
						tmpbuf[z] = exp(- tmpbuf[z]*tmpbuf[z] / (2*stds[z+w*nbstates]*stds[z+w*nbstates] )) / ( stds[z+w*nbstates]  * sqrt(2 * M_PI));
					}else tmpbuf[z] =0;
					if (w == 0) buffer[z] = tmpbuf[z];
					else buffer[z] *= tmpbuf[z];
				}
			}
			tmp = buffer[0];
			for(z=1;z<outi->channels;z++){tmp += buffer[z];}
			if (tmp != 0){
				for(z=0;z<outi->channels;z++){buffer[z] /= tmp;}
			}else buffer[0] =1.0f;
			outi->setPixel(x,y,buffer);
		}
	}
	
}

template<int flag>
template<class T>
void Classifier<flag>::updateClassifier(vector<Image<T>* > &input,Image<float>* outi, Image<float>* filter){
	memset(counts,'\0',sizeof(double)*nbstates);
	memset(boundspr,'\0',sizeof(double)*nbstates);
	memset(means,'\0',sizeof(double)*nbstates* inputsize);
	memset(stds,'\0',sizeof(double)*nbstates*inputsize);
	memset(transit,'\0',sizeof(double)*nbstates*nbstates);
	
	double buffer[32];
	double ibuffer[32];
	double obuffer[32];
	double wbuffer[32];
	int x,y,z,w;
	for(x=0;x<32;x++) wbuffer[x] =1.0f;

	
	

	for(x=0;x<outi->sizex;x++){
		for(y=0;y<outi->sizey;y++){
			outi->getPixel(x,y,buffer);
			
			if (filter != NULL) {filter->getPixel(x,y,wbuffer);	for(z=0;z<nbstates;z++)buffer[z] *= wbuffer[0];	}
			
			for(z=0;z<nbstates;z++){
				if (! isnan( buffer[z]))counts[z] += buffer[z];
				else buffer[z] =0.0f;

			}
			
			for(w=0;w<inputsize;w++){
				input[w]->getPixel(x,y,ibuffer);			
				for(z=0;z<nbstates;z++){
					
					means[z+w*nbstates] += ibuffer[0] * buffer[z];
					stds[z+w*nbstates] += ibuffer[0] * ibuffer[0] * buffer[z];
				}
			}
			if (x>0){
				outi->getPixel(x-1,y,obuffer);
				for(z=0;z<nbstates;z++){
					if (isnan(obuffer[z])) obuffer[z] =0.0f;
				}
				
				if (filter != NULL) {filter->getPixel(x-1,y,wbuffer);for(w=z;w<nbstates;w++) obuffer[w] *= wbuffer[0];}
				for(z=0;z<nbstates;z++){

					for(w=z;w<nbstates;w++){
						transit[w + z *nbstates] += obuffer[z] * buffer[w] + obuffer[w] * buffer[z];
					}
				}
			}
			if (y>0){
				outi->getPixel(x,y-1,obuffer);
				for(z=0;z<nbstates;z++){
					if (isnan(obuffer[z])) obuffer[z] =0.0f;
				}
				if (filter != NULL) {filter->getPixel(x,y-1,wbuffer);for(w=z;w<nbstates;w++) obuffer[w] *= wbuffer[0];}
				for(z=0;z<nbstates;z++){
					
					for(w=z;w<nbstates;w++){
						transit[w + z *nbstates] += obuffer[z] * buffer[w] + obuffer[w] * buffer[z];
					}
				}
			}
		}
	}
	
	
	for(w=0;w<inputsize;w++){
		for(z=0;z<outi->channels;z++){
	//		printf("%f\t%f\t%f\n",means[z+w*nbstates],stds[z+w*nbstates], counts[z]);
			if (counts[z] == 0) {means[z+w*nbstates] = getPriorMean(w,z); stds[z+w*nbstates] = getPriorStd(w,z);}
			else{
				
				means[z+w*nbstates] = means[z+w*nbstates] / counts[z];
				stds[z+w*nbstates] = (stds[z+w*nbstates] / counts[z]) - means[z+w*nbstates]*means[z+w*nbstates];
				stds[z+w*nbstates] = sqrt(stds[z+w*nbstates]);
			}
		}
	}
	for(z=0;z<nbstates;z++){
		for(w=0;w<z;w++){
			transit[w + z *nbstates] = transit[z + w *nbstates];
		}
	}
	
//	applyConstraints();
	double tmp;
	for(z=0;z<nbstates;z++){
		for(w=0;w<nbstates;w++){
			boundspr[z] += transit[z + w *nbstates];
		}
		for(w=0;w<nbstates;w++){
			if (boundspr[z] != 0.0f) transit[z + w *nbstates] /= boundspr[z];
		}
	}
	tmp = 0;
	for(w=0;w<nbstates;w++){
		tmp +=counts[w];
	}	
	for(w=0;w<nbstates;w++){
		if (tmp != 0.0f) boundspr[w] = counts[w] / tmp;
	}	
}

template<int flag>
template<class T>
void Classifier<flag>::applyClassifier(vector<Image<T>* > &input,Image<float>* outi, int rect[4]){
	
	double buffer[32];
	double tmpbuf[32];
	double ibuffer[32];
	double tmp;
	int x,y,z,w;
	
	//		double prior[2];
	//	prior[0] = 0.033f;
	//	prior[1] = 1 - prior[0];
	
	for(w=0;w<inputsize;w++){
		for(z=0;z<outi->channels;z++){
			if (stds[z+w*nbstates] <= 0) stds[z+w*nbstates] =0.00000000001f;
		}
	}
	
	for(x=rect[0];x<=rect[2];x++){
		for(y=rect[1];y<=rect[3];y++){
			
			for(w=0;w<inputsize;w++){
				input[w]->getPixel(x,y,ibuffer);
				if (flag & 1) ibuffer[0] = log(ibuffer[0] +0.00001f);
				for(z=0;z<outi->channels;z++){
					if (stds[z+w*nbstates] > 0){
						tmpbuf[z] = (ibuffer[0] - means[z+w*nbstates]);
						tmpbuf[z] = exp(- tmpbuf[z]*tmpbuf[z] / (2*stds[z+w*nbstates]*stds[z+w*nbstates] )) / ( stds[z+w*nbstates]  * sqrt(2 * M_PI));
					}else tmpbuf[z] =0;
					if (w == 0) buffer[z] = tmpbuf[z];
					else buffer[z] *= tmpbuf[z];
				}
			}
			tmp = buffer[0];
			for(z=1;z<outi->channels;z++){tmp += buffer[z];}
			if (tmp != 0){
				for(z=0;z<outi->channels;z++){buffer[z] /= tmp;}
			}else for(z=0;z<outi->channels;z++) buffer[z] =1.0f / (outi->channels);
			if ((isnan(buffer[0]))||(isnan(buffer[1]))){
				printf("honot");
			}
			outi->setPixel(x,y,buffer);
		}
	}
	
}

template<int flag>
template<class T>
void Classifier<flag>::updateClassifier(vector<Image<T>* > &input,Image<float>* outi, Image<float>* filter, int rect[4], int orifilter){
	memset(counts,'\0',sizeof(double)*nbstates);
	memset(boundspr,'\0',sizeof(double)*nbstates);
	memset(means,'\0',sizeof(double)*nbstates* inputsize);
	memset(stds,'\0',sizeof(double)*nbstates*inputsize);
	memset(transit,'\0',sizeof(double)*nbstates*nbstates);
	
	double buffer[32];
	double ibuffer[32];
	double obuffer[32];
	double wbuffer[32];
	int x,y,z,w;
	for(x=0;x<32;x++) wbuffer[x] =1.0f;
	
	
	
	
	for(x=rect[0];x<=rect[2];x++){
		for(y=rect[1];y<=rect[3];y++){
			outi->getPixel(x,y,buffer);
			
			if (filter != NULL) filter->getPixel(x,y,wbuffer);
			if (wbuffer[0] != 0.0f){
			for(z=0;z<nbstates;z++){
				buffer[z] *= wbuffer[0];
				counts[z] += buffer[z];
			}
			for(w=0;w<inputsize;w++){
				input[w]->getPixel(x,y,ibuffer);			
				for(z=0;z<nbstates;z++){
					
					means[z+w*nbstates] += ibuffer[0] * buffer[z];
					stds[z+w*nbstates] += ibuffer[0] * ibuffer[0] * buffer[z];
				}
			}
				
			if (((orifilter & 1) == 0)&&(x>rect[0])){
				outi->getPixel(x-1,y,obuffer);
				if (filter != NULL) {filter->getPixel(x-1,y,wbuffer);	for(w=0;w<nbstates;w++) obuffer[w] *= wbuffer[0];}
				for(z=0;z<nbstates;z++){
					
					for(w=z;w<nbstates;w++){
						transit[w + z *nbstates] += obuffer[z] * buffer[w] + obuffer[w] * buffer[z];
					}
				}
			}
			if (((orifilter & 2) == 0)&&(y>rect[1])){
				outi->getPixel(x,y-1,obuffer);
				if (filter != NULL) {filter->getPixel(x,y-1,wbuffer);		for(w=0;w<nbstates;w++) obuffer[w] *= wbuffer[0];}
				for(z=0;z<nbstates;z++){
					for(w=z;w<nbstates;w++){
						transit[w + z *nbstates] += obuffer[z] * buffer[w] + obuffer[w] * buffer[z];
					}
				}
			}
			}
		}
	}
	for(w=0;w<inputsize;w++){
		for(z=0;z<outi->channels;z++){
			if (counts[z] == 0) {means[z+w*nbstates] = getPriorMean(w,z); stds[z+w*nbstates] = getPriorStd(w,z);}
			else{
				means[z+w*nbstates] = means[z+w*nbstates] / counts[z];
				stds[z+w*nbstates] = (stds[z+w*nbstates] / counts[z]) - means[z+w*nbstates]*means[z+w*nbstates];
				stds[z+w*nbstates] = sqrt(stds[z+w*nbstates]);
			}
		}
	}
	for(z=0;z<nbstates;z++){
		for(w=0;w<z;w++){
			transit[w + z *nbstates] = transit[z + w *nbstates];
		}
	}
	
	//	applyConstraints();
	double tmp;
	for(z=0;z<nbstates;z++){
		for(w=0;w<nbstates;w++){
			boundspr[z] += transit[w + z *nbstates];
		}
		for(w=0;w<nbstates;w++){
			if (boundspr[z] != 0.0f) transit[w + z *nbstates] /= boundspr[z];
		}
	}
	tmp = 0;
	for(w=0;w<nbstates;w++){
		tmp +=counts[w];
	}	
	for(w=0;w<nbstates;w++){
		if (tmp != 0.0f) boundspr[w] = counts[w] / tmp;
	}	
}

template<int flag>
template<class T>
void Classifier<flag>::updateClassifier(vector<Image<T>* > &input,Image<float>* outi, Image<float>* filter, int rect[4], double* expect_transit){
	memset(counts,'\0',sizeof(double)*nbstates);
	memset(boundspr,'\0',sizeof(double)*nbstates);
	memset(means,'\0',sizeof(double)*nbstates* inputsize);
	memset(stds,'\0',sizeof(double)*nbstates*inputsize);
	memset(transit,'\0',sizeof(double)*nbstates*nbstates);
	
	double buffer[32];
	double ibuffer[32];
	double wbuffer[32];
	int x,y,z,w;
	for(x=0;x<32;x++) wbuffer[x] =1.0f;
	
	
	
	
	for(x=rect[0];x<=rect[2];x++){
		for(y=rect[1];y<=rect[3];y++){
			outi->getPixel(x,y,buffer);
			
			if (filter != NULL) filter->getPixel(x,y,wbuffer);
			if (wbuffer[0] != 0.0f){

				for(z=0;z<nbstates;z++){
				//	printf("%f\t%f\n",buffer[z],wbuffer[0]);
					buffer[z] *= wbuffer[0];
					counts[z] += buffer[z];
				}
				for(w=0;w<inputsize;w++){
					input[w]->getPixel(x,y,ibuffer);	
					if (flag & 1) ibuffer[0] = log(ibuffer[0] +0.00001f);
					for(z=0;z<nbstates;z++){
						
						means[z+w*nbstates] += ibuffer[0] * buffer[z];
						stds[z+w*nbstates] += ibuffer[0] * ibuffer[0] * buffer[z];
					}
				}
				
			}
		}
	}
	for(w=0;w<inputsize;w++){
		for(z=0;z<outi->channels;z++){
			if (counts[z] == 0) {means[z+w*nbstates] = getPriorMean(w,z); stds[z+w*nbstates] = getPriorStd(w,z);}
			else{
				means[z+w*nbstates] = means[z+w*nbstates] / counts[z];
				stds[z+w*nbstates] = sqrt((stds[z+w*nbstates] / counts[z]) - means[z+w*nbstates]*means[z+w*nbstates]);
			}
		}
	}
	
	//	applyConstraints();
	double tmp;
	for(z=0;z<nbstates;z++){
		//	expect_transit[z + z *nbstates] *= 10;
		for(w=0;w<nbstates;w++){
			boundspr[z] += (expect_transit[w + z *nbstates] + 1.0f);
		}
		for(w=0;w<nbstates;w++){
			 transit[w + z *nbstates] = (expect_transit[w + z *nbstates] + 1.0f) / boundspr[z];
		}
	}
	tmp = 0;
	for(w=0;w<nbstates;w++){
		tmp +=counts[w];
	}	
	for(w=0;w<nbstates;w++){
		if (tmp != 0.0f) boundspr[w] = counts[w] / tmp;
	}	
}

template<class C>
void MultiCover::drawDarkCircles(Image<C>* where){
	int i,j,x,y;
	CellPose* wh;
	
	double pix[32];
	double val,tmp, tmp2;
	int rad;
	for(i=0;i<group.size();i++) for(j=0;j<group[i].cover.size();j++){
		wh = &(group[i].cover[j]);
		rad = (int)((wh->width+1) /2);
		
//		tmp2 = (wh->error[0] + 0.12f) * 2;
//		tmp = wh->error[1] * 10.0f;
//		goodness = exp(-tmp2*tmp2 -tmp);
//		printf("that's good %f\n", goodness);
		for(x = 1-rad ; x < rad;x++) for(y = 1-rad ; y < rad;y++){
			tmp = x - wh->eccentric[0]; tmp2 = tmp*tmp;
			tmp = y - wh->eccentric[1]; val = sqrt(tmp2 + tmp*tmp);
			tmp = x + wh->eccentric[0]; tmp2 = tmp*tmp;
			tmp = y + wh->eccentric[1]; val += sqrt(tmp2 + tmp*tmp);
			tmp = 4*(val/2 - rad)/rad;   
			val = exp(-tmp*tmp);
			if (val > 0.01f){
				where->getPixel((unsigned int)(x + wh->center[0] + group[i].rect[0]),(unsigned int)(y + wh->center[1] + group[i].rect[1]),pix);
				pix[0] += val * 0.25f;
				
				if (pix[0] > 1.0f) pix[0] = 1.0f;
				if (pix[2] > 1.0f) pix[2] = 1.0f;
				where->setPixel((unsigned int)(x + wh->center[0] + group[i].rect[0]),(unsigned int)(y + wh->center[1] + group[i].rect[1]),pix);
			}
		}
	}
	
}

template<class T>
void MultiCover::drawACircle(Image<T>* where, int channel, double value, int i, int j, bool ignore_conf){
	int x,y;
	CellPose* wh;
	int rect[4];
	double pix[32];
	double val,tmp, tmp2;
	int rad;
	
//	for(i=0;i<group.size();i++) {
		
		if (group[i].rect[0] < 0){
			group[i].rect[2] += group[i].rect[0];group[i].rect[0] =0;
		}
		if (group[i].rect[1] < 0){
			group[i].rect[3] += group[i].rect[1];group[i].rect[1] =0;
		}
		
//		for(j=0;j<group[i].cover.size();j++){
			wh = &(group[i].cover[j]);
			rad = (int)((wh->width+1) /2);
			
			//		printf("circle %i,%i\n",  group[i].rect[2], group[i].rect[0]);
			
			rect[0] = (int)((wh->center[0]) + 1-rad< 0 ? 0 : (wh->center[0]) + 1-rad);
			rect[1] = (int)(1-rad  +(wh->center[1]) < 0  ? 0 : 1-rad  +(wh->center[1])) ;
			rect[2] = (int)(rad + wh->center[0]>  group[i].rect[2] ? group[i].rect[2]: rad + (wh->center[0]));
			rect[3] = (int)(rad + wh->center[1]>  group[i].rect[3] ? group[i].rect[3]: rad + (wh->center[1]));
			//	printf("%i,%i\n", rect[2] +  group[i].rect[0],rect[3] +  group[i].rect[1] );fflush(stdout);
			if (rect[2] + group[i].rect[0] > where->sizex) rect[2]  = where->sizex - group[i].rect[0];
			if (rect[3] + group[i].rect[1] > where->sizey) rect[3]  = where->sizey - group[i].rect[1];
			for(x = rect[0]; x < rect[2];x++) for(y = rect[1] ; y < rect[3];y++){
				tmp = x - wh->center[0] - wh->eccentric[0]; tmp2 = tmp*tmp;
				tmp = y - wh->center[1] - wh->eccentric[1]; val = sqrt(tmp2 + tmp*tmp);
				tmp = x - wh->center[0] + wh->eccentric[0]; tmp2 = tmp*tmp;
				tmp = y - wh->center[1] + wh->eccentric[1]; val += sqrt(tmp2 + tmp*tmp);
				tmp = 4*(val/2 - rad)/rad;   
				val = exp(-tmp*tmp*50.0f);
				if (val > 0.01f){
					where->getPixel(x + group[i].rect[0],y + group[i].rect[1],pix);
					if (pix[channel] < val){
						pix[channel] = val;	
						where->setPixel(x + group[i].rect[0],y + group[i].rect[1],pix);
					}
				}
//			}}
	}
}

template<class T>
void MultiCover::drawCircles(Image<T>* where, int channel, double value, bool ignore_conf){
	int i,j,x,y;
	CellPose* wh;
	int rect[4];
	double pix[32];
	double val,tmp, tmp2;
	int rad;

	for(i=0;i<group.size();i++) {
		
		if (group[i].rect[0] < 0){
			group[i].rect[2] += group[i].rect[0];group[i].rect[0] =0;
		}
		if (group[i].rect[1] < 0){
			group[i].rect[3] += group[i].rect[1];group[i].rect[1] =0;
		}
		
		for(j=0;j<group[i].cover.size();j++){
		wh = &(group[i].cover[j]);
		rad = (int)((wh->width+1) /2);
		
//		printf("circle %i,%i\n",  group[i].rect[2], group[i].rect[0]);

		rect[0] = (int)((wh->center[0]) + 1-rad< 0 ? 0 : (wh->center[0]) + 1-rad);
		rect[1] = (int)(1-rad  +(wh->center[1]) < 0  ? 0 : 1-rad  +(wh->center[1])) ;
		rect[2] = (int)(rad + wh->center[0]>  group[i].rect[2] ? group[i].rect[2]: rad + (wh->center[0]));
		rect[3] = (int)(rad + wh->center[1]>  group[i].rect[3] ? group[i].rect[3]: rad + (wh->center[1]));
	//	printf("%i,%i\n", rect[2] +  group[i].rect[0],rect[3] +  group[i].rect[1] );fflush(stdout);
		if (rect[2] + group[i].rect[0] > where->sizex) rect[2]  = where->sizex - group[i].rect[0];
		if (rect[3] + group[i].rect[1] > where->sizey) rect[3]  = where->sizey - group[i].rect[1];
		for(x = rect[0]; x < rect[2];x++) for(y = rect[1] ; y < rect[3];y++){
			tmp = x - wh->center[0] - wh->eccentric[0]; tmp2 = tmp*tmp;
			tmp = y - wh->center[1] - wh->eccentric[1]; val = sqrt(tmp2 + tmp*tmp);
			tmp = x - wh->center[0] + wh->eccentric[0]; tmp2 = tmp*tmp;
			tmp = y - wh->center[1] + wh->eccentric[1]; val += sqrt(tmp2 + tmp*tmp);
			tmp = 4*(val/2 - rad)/rad;   
			val = exp(-tmp*tmp*50.0f);
			if (val > 0.01f){
				where->getPixel(x + group[i].rect[0],y + group[i].rect[1],pix);
				if (pix[channel] < val){
					pix[channel] = val;	
					where->setPixel(x + group[i].rect[0],y + group[i].rect[1],pix);
				}
			}
		}}
	}
}


template<class T>
void MultiCover::drawCross(Image<T>* where, int channel, double value){ 
	int i,j,x,y;
	CellPose* wh;
	int rect[4];
	double pix[32];
	double val,tmp, tmp2;
	int rad;
	
	for(i=0;i<group.size();i++) for(j=0;j<group[i].cover.size();j++){
		wh = &(group[i].cover[j]);
		rad = (int)(wh->width);
		
		//		rect[0] = (wh->center[0]) + 1-rad< 0 ? 0 : (wh->center[0]) + 1-rad;
		//		rect[1] = 1-rad  +(wh->center[1]) < 0  ? 0 : 1-rad  +(wh->center[1]) ;
		//		rect[2] = rad + wh->center[0]>  group[i].rect[2] -group[i].rect[0] ? group[i].rect[2]-group[i].rect[0] : rad + (wh->center[0]);
		//		rect[3] = rad + wh->center[1]>  group[i].rect[3] -group[i].rect[1] ? group[i].rect[3]-group[i].rect[1]: rad + (wh->center[1]);
				rect[0] = 0;
				rect[1] = 0;
				rect[2] = group[i].rect[2] - group[i].rect[0];
				rect[3] = group[i].rect[3] - group[i].rect[1];
		//	printf("%i,%i\n", rect[2] +  group[i].rect[0],rect[3] +  group[i].rect[1] );fflush(stdout);
		for(x = rect[0]; x < rect[2];x++) for(y = rect[1] ; y < rect[3];y++){
			tmp = x - wh->center[0] - wh->eccentric[0]; tmp2 = tmp*tmp;
			tmp = y - wh->center[1] - wh->eccentric[1]; val = sqrt(tmp2 + tmp*tmp);
			tmp = x - wh->center[0] + wh->eccentric[0]; tmp2 = tmp*tmp;
			tmp = y - wh->center[1] + wh->eccentric[1]; val += sqrt(tmp2 + tmp*tmp);
			tmp = val/rad - 1;   
			val = exp(-tmp*tmp*100.0f);
			if  ((!(((((int)(x - wh->center[0]+1)) & 0xFFFFFFFC)&&((int)(y - wh->center[1]+1))& 0xFFFFFFFC)))&&(tmp <0.0f)){
			where->getPixel(x + group[i].rect[0],y + group[i].rect[1],pix);
				if (pix[channel] < value){
					pix[channel] = value;	
					where->setPixel(x + group[i].rect[0],y + group[i].rect[1],pix);
				}
			
			}else if (val > 0.01f){
				where->getPixel(x + group[i].rect[0],y + group[i].rect[1],pix);
				if (pix[channel] < value*val){
					pix[channel] = value*val;	
					where->setPixel(x + group[i].rect[0],y + group[i].rect[1],pix);
				}
			}
		}
	}
	
}



template<class C>
void MultiCover::drawLabel(Image<C>* where){
	int i,j;
	CellPose* wh;
	

	
	char labels[256];
	int rad;
	for(i=0;i<group.size();i++) {
		group[i].getStageType(labels);
		for(j=0;j<group[i].cover.size();j++){
		wh = &(group[i].cover[j]);
		rad = (wh->width+1) /2;
		where->drawLetter(labels[j],group[i].rect[0] + group[i].cover[j].center[0] -8,group[i].rect[1] + group[i].cover[j].center[1] -8,0);

	}
}
	
}



#undef LFHTEMP
#define LFHTEMP template <class T>


LFHTEMP Image<T>::Image() : data(NULL), bits(sizeof(T)){}
LFHTEMP Image<T>::~Image(){delete[](data);}

LFHTEMP template<class D, unsigned int NBCHAN> Image<T>::Image(LFHPrimitive::DataGrid<LFHPrimitive::Tuple<D, NBCHAN>,2> & other) : data(NULL), bits(sizeof(T)){
	sizex = other.dims[0];
	sizey = other.dims[1];
	channels = NBCHAN;
	allocateBuffer();
	
	int coor[2];
	double pix[32];
	int z;
	
	LFHPrimitive::Tuple<D, NBCHAN> cur;
	
	for(coor[1] =0;coor[1]<sizey;coor[1]++){
		for(coor[0] =0;coor[0]<sizex;coor[0]++){
			cur = other(coor);
			for(z=0;z<NBCHAN;z++){
				pix[z]= (double)cur[z];
			}
			setPixel(coor[0],coor[1],pix);
		}	
	}
	
	
}

LFHTEMP template<class D> Image<T>::Image(LFHPrimitive::DataGrid<D,2> & other) : data(NULL), bits(sizeof(T)){
	sizex = other.dims[0];
	sizey = other.dims[1];
	channels = 1;
	allocateBuffer();
	
	int coor[2];
	double pix[32];
	int z;
	
	
	for(coor[1] =0;coor[1]<sizey;coor[1]++){
		for(coor[0] =0;coor[0]<sizex;coor[0]++){
			pix[0]= other(coor);
			setPixel(coor[0],coor[1],pix);
		}	
	}
	
	
}

LFHTEMP void Image<T>::LoadBmpImage(char* path){
	FILE* f = fopen(path,"rb+");
	char buffer[65536];
	fread(buffer,sizeof(char),54,f);
	memcpy(&sizex,buffer + 18,sizeof(int));
	memcpy(&sizey,buffer + 22,sizeof(int));
	memcpy(&bits,buffer + 26,sizeof(short));
	channels =3;
	int start;
	memcpy(&start,buffer + 10,sizeof(int));
	
	sizex = abs(sizex);
	sizey = abs(sizey);
	data = new unsigned char[bits*3*sizex*sizey];
	int i;
	for(i=0;i< sizey;i++){
		fread(&data[i*sizex*3],sizeof(char)*bits,sizex*3,f);
		if ((sizex*3) % 4 != 0) fread(buffer,sizeof(char),4-((sizex*3)%4),f);
	}
	//  printf("Opened %s sucessfully\n",path);
	//  printf("size: %i,%i\n", sizex, sizey);
	fclose(f);
}

LFHTEMP void Image<T>::SaveBmpImage(char* path){
	char buffer[65536];
	FILE* f = fopen(path,"rb+");
	int i = strlen(path);
	memcpy(buffer,path,i* sizeof(char));
	buffer[i-4] = '_';
	buffer[i-3] = 'o';
	buffer[i-2] = '.';
	buffer[i-1] = 'b';
	buffer[i] = 'm';
	buffer[i+1] = 'p';
	buffer[i+2] = '\0';
	FILE* g = fopen(buffer,"wb+");
	fread(buffer,sizeof(char),54,f);
	fwrite(buffer,sizeof(char),54,g);
	memset(buffer,'\0',sizeof(char)*4);
	for(i=0;i< sizey;i++){
		fwrite(&data[i*sizex*3],sizeof(char)*bits,sizex*3,g);
		if ((sizex*3) % 4 != 0) fwrite(buffer,sizeof(char),4-((sizex*3)%4),g);
	}
	fclose(f);
	fclose(g);
}


LFHTEMP void Image<T>::LoadTiffImage(vector<Image<T>* > &out, const char* path){
	FILE* in = fopen(path,"rb+");
	if (in == NULL) return;
	char buffer[65536];
	fread(buffer,sizeof(char),8,in);
	stack<int> s;
	int cur = *((int*)(buffer + 4));
	int i;
	Image<T>* curimage;
	void* target;
	int type;
	int x,y,w;
	int stripof; // address
	int stripco;
	int stripro;
	double pix[32];
	//	int nbstrip;
	
	int fieldsize;
	while(cur != 0){
		curimage = new Image<T>();
		fseek(in,cur,SEEK_SET);
		fread(buffer,sizeof(short),1,in);
		i = *((short*)(buffer));
		for(i--;i>=0;i--){
			fread(buffer,sizeof(char),12,in);
			
			switch(*((short*)(buffer))){
					//			case 254:memcpy(buffer+512,"NewSubfileType",sizeof(char)*15);break;
				case 256:target = (void*)&(curimage->sizex); type =4; break;//memcpy(buffer+512,"ImageWidth",sizeof(char)*11);break;
				case 257:target = (void*)&(curimage->sizey); type =4; break;//memcpy(buffer+512,"ImageHeight",sizeof(char)*12);break;
					//			case 258:memcpy(buffer+512,"BitsPerSample",sizeof(char)*14);break;
					//			case 259:memcpy(buffer+512,"Compression",sizeof(char)*12);break;
					//			case 262:memcpy(buffer+512,"PhotometricInterpretation",sizeof(char)*26);break; 
				case 273:target = (void*)&(stripof); type =4; break;//memcpy(buffer+512,"SkipOffsets",sizeof(char)*12);break;
					//			case 274:memcpy(buffer+512,"Orientation",sizeof(char)*12);break;
				case 277:target = (void*)&(curimage->channels); type =4; break;//memcpy(buffer+512,"SamplesPerPixel",sizeof(char)*16);break;
				case 278:target = (void*)&(stripro); type =4; break;//memcpy(buffer+512,"RowsPerStrip",sizeof(char)*13);break;
				case 279:target = (void*)&(stripco); type =4; break;//memcpy(buffer+512,"StripByteCounts",sizeof(char)*16);break;
					//			case 282:memcpy(buffer+512,"XResolution",sizeof(char)*12);break;
					//			case 283:memcpy(buffer+512,"YResolution",sizeof(char)*12);break;
					//			case 296:memcpy(buffer+512,"ResolutionUnit",sizeof(char)*15);break;
				default: target = NULL; //memcpy(buffer+512,"",sizeof(char)*1);
			}
			
			switch(*((short*)(buffer+2))){
				case 1:fieldsize = sizeof(char);break;//memcpy(buffer+256,"BYTE",sizeof(char)*5);break;
				case 2:fieldsize = sizeof(char);break;//memcpy(buffer+256,"ASCII",sizeof(char)*6);break;
				case 3:fieldsize = sizeof(short);break;//memcpy(buffer+256,"SHORT",sizeof(char)*6);break;
				case 4:fieldsize = sizeof(int);break;//memcpy(buffer+256,"LONG",sizeof(char)*5);break;
				case 5:fieldsize = sizeof(int)*2;break;//memcpy(buffer+256,"RATIONAL",sizeof(char)*9);break;
				case 6:fieldsize = sizeof(char);break;//memcpy(buffer+256,"SBYTE",sizeof(char)*6);break;
				case 7:fieldsize = sizeof(char);break;//memcpy(buffer+256,"UNDEFINED",sizeof(char)*10);break;
				case 8:fieldsize = sizeof(short);break;//memcpy(buffer+256,"SSHORT",sizeof(char)*7);break;
				case 9:fieldsize = sizeof(int);break;//memcpy(buffer+256,"SLONG",sizeof(char)*6);break;
				case 10:fieldsize = sizeof(int)*2;break;//memcpy(buffer+256,"SRATIONAL",sizeof(char)*10);break;
				case 11:fieldsize = sizeof(float)*2;break;//memcpy(buffer+256,"FLOAT",sizeof(char)*6);break;
				case 12:fieldsize = sizeof(double)*2;break;//memcpy(buffer+256,"DOUBLE",sizeof(char)*7);break;
				default: memcpy(buffer+256,"UNKNOWN",sizeof(char)*8);
			}
			if (target != NULL){
				if ((*((int*)(buffer+4)))*fieldsize <= 4){
					switch(*((short*)(buffer+2))){
						case 1: if (type == 4) *((int*)target) = (int) *((char*)(buffer + 8)); break;
						case 2: if (type == 4) *((int*)target) = (int) *((char*)(buffer + 8)); break;
						case 3: if (type == 4) *((int*)target) = (int) *((short*)(buffer + 8)); break;
						case 4: if (type == 4) *((int*)target) = (int) *((int*)(buffer + 8)); break;
					}
				}else{
					// for now, that's not used
					
				}
			}
			
		}
		fread(&cur,sizeof(int),1,in);
		fseek(in,stripof,SEEK_SET);
		curimage->allocateBuffer();
		
		for(y=0;y<curimage->sizey;y++) for(x=0;x<curimage->sizex;x++){
			fread(buffer,sizeof(char),curimage->channels,in); // assumes a single strip...	
			for(w=0;w<curimage->channels;w++) pix[w] = ((double)buffer[w]) /255.0f;
			curimage->setPixel(x,y,pix);
		}
		
		
		out.push_back(curimage);
	}
	
	fclose(in);
}

LFHTEMP void Image<T> ::SaveTiffImage(vector<Image<T>* > &list, const char* path, int ignore, int nbchannels){
	if (nbchannels == 0) nbchannels = list[0]->channels - ignore;
	FILE* out = fopen(path,"wb+");
	char buffer[65536];
	buffer[0] = 'I';
	buffer[1] = 'I';
	*(short*)(buffer + 2) = 42;
	char *p;
	*(int*)(buffer + 4) = 8;
	fwrite(buffer,sizeof(char),8,out);
	int imc=0;
	int cur = 8;
	int nbflag;
	int extra;
	int x,y,w;
	double pix[32];
	for(imc=0;imc<list.size();imc++){
		p = buffer;
		extra =0;
		nbflag = 14;
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
		*(int*)(p+8) = list[imc]->sizex; 
		p+=12;
		
		*(short*)(p) = 257; 
		*(short*)(p+2) = 4; 
		*(int*)(p+4) = 1; 
		*(int*)(p+8) = list[imc]->sizey; 
		p+=12;
		
		*(short*)(p) = 258; 
		*(short*)(p+2) = 3; 
		if (nbchannels == 1){
			*(int*)(p+4) = 1;
			*(short*)(p+8) = 8; //bitpersample
		}else{
			*(int*)(p+4) = 3;
			*(int*)(p+8) =cur + 6 + 12 * nbflag + extra;
			extra += sizeof(short)*3;
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
		switch(nbchannels){//photothing
			case 1: *(short*)(p+8) = 1; break;
			case 2: *(short*)(p+8) = 1; break;
			default: *(short*)(p+8) = 2; break;
		}
		p+=12;
		
		*(short*)(p) = 273; 
		*(short*)(p+2) = 4; 
		*(int*)(p+4) = 1; 
		*(int*)(p+8) = cur + 6 + 12 * nbflag + 16 + extra; //data start
		p+=12;
		
		*(short*)(p) = 277; 
		*(short*)(p+2) = 3; 
		*(int*)(p+4) = 1; 
		*(short*)(p+8) = nbchannels; 
		p+=12;
		
		*(short*)(p) = 278; 
		*(short*)(p+2) = 4; 
		*(int*)(p+4) = 1; 
		*(int*)(p+8) = list[imc]->sizey; //Rowperstrip
		p+=12;
		
		*(short*)(p) = 279; 
		*(short*)(p+2) = 4; 
		*(int*)(p+4) = 1; 
		*(int*)(p+8) = list[imc]->sizex*list[imc]->sizey*nbchannels; //NBbyte
		p+=12;
		
		*(short*)(p) = 282; 
		*(short*)(p+2) = 5; 
		*(int*)(p+4) = 1; 
		*(int*)(p+8) = cur + 6 + 12 * nbflag + extra;
		extra += sizeof(int)*2;
		p+=12;
		
		*(short*)(p) = 283; 
		*(short*)(p+2) = 5; 
		*(int*)(p+4) = 1; 
		*(int*)(p+8) = cur + 6 + 12 * nbflag + extra;
		extra += sizeof(int)*2;
		p+=12;
		
		*(short*)(p) = 284; 
		*(short*)(p+2) = 3; 
		*(int*)(p+4) = 1; 
		*(short*)(p+8) = 1; //planar_rep
		p+=12;
		
		*(short*)(p) = 296; 
		*(short*)(p+2) = 3; 
		*(int*)(p+4) = 1; 
		*(short*)(p+8) = 3; //unit
		p+=12;
		
		/*	*(short*)(p) = 338; 
		 *(short*)(p+2) = 3; 
		 *(int*)(p+4) = 1; 
		 *(short*)(p+8) = 3; //unit
		 p+=12;*/
		
		
		cur += 6 + 12*nbflag + extra + list[imc]->sizex* list[imc]->sizey* nbchannels;
		if ((imc + 1) == list.size()) cur = 0;
		*(int*)(p) = cur;
		
		if (nbchannels == 1){
			*(int*)(p+4) = 1;
			*(int*)(p+8) = 1;
			*(int*)(p+12) = 1;
			*(int*)(p+16) = 1;
		}else{
			*(short*)(p+4) = 8;
			*(short*)(p+6) = 8;
			*(short*)(p+8) = 8;			
			*(int*)(p+10) = 1;
			*(int*)(p+14) = 1;
			*(int*)(p+18) = 1;
			*(int*)(p+22) = 1;
		}
		fwrite(buffer,sizeof(char), 6 + 12*nbflag + extra,out);
		
		for(y=0;y<list[imc]->sizey;y++) {for(x=0;x<list[imc]->sizex;x++){
			list[imc]->getPixel(x,y,pix);
			for(w=ignore;w<ignore + nbchannels;w++) {
				if (pix[w] >= 1.0f) buffer[w] = 255;
				else if (pix[w] <=0.0f) buffer[w] =0;
				else buffer[w] = (char)(pix[w] *255.1f);
			}
			fwrite(buffer+ignore,sizeof(char), nbchannels,out);
		}}
		
	}
	fclose(out);
	
}


LFHTEMP void Image<T> ::getMaxPixel(double *out){
	out[0] =out[1] =out[2] = 255;
}

LFHTEMP void Image<T>::shiftColor(int x, int y, double factor){
	if ((x<0)||(y<0)||(x>=sizex)||(y>=sizey)) return;
	double modif;
	int i;
	for(i=0;i<3;i++){
		modif = data[(x + y * sizex)*3];
		if (modif < 0) modif += 128;
		modif = (i==0 ? 0 : 255) * factor + modif * (1- factor);
		if (modif < 128) data[(x + y * sizex)*3+i] = modif; else data[(x + y * sizex)*3+i] = modif-128;
	}
}

LFHTEMP void Image<T> ::computeColorDistribution(double center[2], double radius, double* correlval, double &n, double *filter){
	int x,y,i;
	double tmp;
	double dist;
	double value[3];
	double b2 = (radius+1)*(radius+1);
	int bounds[4];
	bounds[0] = (int)(center[0]-radius-1 <0 ? 0 : center[0]-radius-1);
	bounds[1] = (int)(center[1]-radius-1 <0 ? 0 : center[1]-radius-1);
	bounds[2] = (int)(center[0]+radius+2 >=sizex ? sizex : center[0]+radius+2);
	bounds[3] = (int)(center[1]+radius+2 >=sizey ? sizey : center[1]+radius+2);
	for(x = bounds[0];x< bounds[2];x++){
		for(y = bounds[1];y< bounds[3];y++){
			tmp = x - center[0]; dist = tmp * tmp; tmp = y - center[1]; dist += tmp*tmp;
			if (dist < b2) {
				dist = sqrt(dist)/radius;
				getPixel(x,y, value);
				//		correlval[0] += dist; 
				//		correlval[1] += dist*dist;
				dist = (1-dist*dist)/ (1 + exp(1 * filter[x + y*sizex]));
				for(i=0;i<3;i++){
					correlval[0+2*i] += value[i] * dist; 
					correlval[1+2*i] += value[i]*value[i] * dist; 
					//correlval[4+3*i] += dist*value[i]; 
				}
				n += dist;
			}
		}
	}
}

LFHTEMP double Image<T> ::computeWhiteness(double center[2], double radius){
	int x,y,i;
	int n=0;
	double tmp;
	double dist;
	double out =0.0f;
	double value[3];
	double b2 = (radius+1)*(radius+1);
	int bounds[4];
	bounds[0] = center[0]-radius-1 <0 ? 0 : center[0]-radius-1;
	bounds[1] = center[1]-radius-1 <0 ? 0 : center[1]-radius-1;
	bounds[2] = center[0]+radius+2 >=sizex ? sizex : center[0]+radius+2;
	bounds[3] = center[1]+radius+2 >=sizey ? sizey : center[1]+radius+2;
	for(x = bounds[0];x< bounds[2];x++){
		for(y = bounds[1];y< bounds[3];y++){
			tmp = x - center[0]; dist = tmp * tmp; tmp = y - center[1]; dist += tmp*tmp;
			if (dist < b2) {
				dist = sqrt(dist)/radius;
				getPixel(x,y, value); 
				out += (value[0]*value[1]*value[2])/65536;
				if ((value[0]*value[1]*value[2]) < 0) printf("ouch");
				n++;
			}
		}
	}
	return(out / n);
}

/* This method tries to the best coordonates which maximise the background texture occurence on the interior off-spots regions
 *
 */

LFHTEMP double Image<T> ::findGrid_evalroutine(double *coords, int col, int row,double size){
	double buffer[6];
	double n;
	double center[2];
	memset(buffer,'\0',sizeof(double)*6);
	n=0;
	int x,y;
	
	for(x=1;x<col;x++){
		for(y=1;y<row;y++){
			center[0] = (coords[0]* (2.0f -((double)2*x-1)/col-((double)2*y-1)/row)  +  ((2*x-1)*coords[2])/col + ((2*y-1)*coords[4])/row)/2;
			center[1] = (coords[1]* (2.0f -((double)2*x-1)/ col-((double)2*y-1)/row)  +  ((2*x-1)*coords[3])/col + ((2*y-1)*coords[5])/row)/2;
			if ((center[0] < 0)||(center[1] < 0)||(center[0] > sizex)||(center[1] > sizey)) return(0);
			computeColorDistribution(center, size,buffer,n, NULL);
		}
	}
	return(n * 10.0f /((buffer[1]+buffer[3]+buffer[5] - (buffer[0]*buffer[0]+buffer[2]*buffer[2]+buffer[4]*buffer[4]) / n)/n));
}

LFHTEMP double Image<T> ::findGrid_evalroutine2(double *coords, int col, int row,double size,double size2,double size3, double *filter){
	double buffer[6];
	double inbuffer[6];
	double n;
	double m;
	double center[2];
	memset(buffer,'\0',sizeof(double)*6);
	memset(inbuffer,'\0',sizeof(double)*6);
	n=0;m=0;
	int x,y;
	double out =0;
	for(x=0;x<col*2-1;x++){
		for(y=0;y<row*2-1;y++){
			center[0] = (coords[0]* (2.0f -((double)x)/col-((double)y)/row)  +  ((x)*coords[2])/col + ((y)*coords[4])/row)/2;
			center[1] = (coords[1]* (2.0f -((double)x)/col-((double)y)/row)  +  ((x)*coords[3])/col + ((y)*coords[5])/row)/2;
			if ((center[0] < 0)||(center[1] < 0)||(center[0] > sizex)||(center[1] > sizey)) return(0);
			if ((x & 1)||(y&1)){
				if ((x & 1)&&(y&1))	computeColorDistribution(center, size,buffer,n,filter);
				else computeColorDistribution(center, size2,buffer,n,filter);
			} else computeColorDistribution(center, size3,inbuffer,m,filter);
		}
	}
	out = n *( 5.0f + fabs((inbuffer[0]+inbuffer[2]+inbuffer[4])/m-(buffer[0]+buffer[2]+buffer[4])/n)*0.015f) / ((buffer[1]+buffer[3]+buffer[5] - (buffer[0]*buffer[0]+buffer[2]*buffer[2]+buffer[4]*buffer[4]) / n)/n);
	/*
	 memset(buffer,'\0',sizeof(double)*6);
	 n=0;
	 x=-1;
	 for(y=0;y<row*2-1;y++){
	 center[0] = (coords[0]* (2.0f -((double)x)/col-((double)y)/row)  +  ((x)*coords[2])/col + ((y)*coords[4])/row)/2;
	 center[1] = (coords[1]* (2.0f -((double)x)/col-((double)y)/row)  +  ((x)*coords[3])/col + ((y)*coords[5])/row)/2;
	 if ((center[0] < 0)||(center[1] < 0)||(center[0] > sizex)||(center[1] > sizey)) return(0);
	 if ((x & 1)||(y&1)){
	 if ((x & 1)&&(y&1))	computeColorDistribution(center, size,buffer,n);
	 else computeColorDistribution(center, size2,buffer,n);
	 }
	 }
	 x=col*2-1;
	 for(y=0;y<row*2-1;y++){
	 center[0] = (coords[0]* (2.0f -((double)x)/col-((double)y)/row)  +  ((x)*coords[2])/col + ((y)*coords[4])/row)/2;
	 center[1] = (coords[1]* (2.0f -((double)x)/col-((double)y)/row)  +  ((x)*coords[3])/col + ((y)*coords[5])/row)/2;
	 if ((center[0] < 0)||(center[1] < 0)||(center[0] > sizex)||(center[1] > sizey)) return(0);
	 if ((x & 1)||(y&1)){
	 if ((x & 1)&&(y&1))	computeColorDistribution(center, size,buffer,n);
	 else computeColorDistribution(center, size2,buffer,n);
	 }
	 }
	 y=-1;
	 for(x=0;x<col*2-1;x++){
	 center[0] = (coords[0]* (2.0f -((double)x)/col-((double)y)/row)  +  ((x)*coords[2])/col + ((y)*coords[4])/row)/2;
	 center[1] = (coords[1]* (2.0f -((double)x)/col-((double)y)/row)  +  ((x)*coords[3])/col + ((y)*coords[5])/row)/2;
	 if ((center[0] < 0)||(center[1] < 0)||(center[0] > sizex)||(center[1] > sizey)) return(0);
	 if ((x & 1)||(y&1)){
	 if ((x & 1)&&(y&1))	computeColorDistribution(center, size,buffer,n);
	 else computeColorDistribution(center, size2,buffer,n);
	 }
	 }
	 y=row*2-1;
	 for(x=0;x<col*2-1;x++){
	 center[0] = (coords[0]* (2.0f -((double)x)/col-((double)y)/row)  +  ((x)*coords[2])/col + ((y)*coords[4])/row)/2;
	 center[1] = (coords[1]* (2.0f -((double)x)/col-((double)y)/row)  +  ((x)*coords[3])/col + ((y)*coords[5])/row)/2;
	 if ((center[0] < 0)||(center[1] < 0)||(center[0] > sizex)||(center[1] > sizey)) return(0);
	 if ((x & 1)||(y&1)){
	 if ((x & 1)&&(y&1))	computeColorDistribution(center, size,buffer,n);
	 else computeColorDistribution(center, size2,buffer,n, filter);
	 }
	 }
	 out += 200000.0f / ((buffer[1]+buffer[3]+buffer[5] - (buffer[0]*buffer[0]+buffer[2]*buffer[2]+buffer[4]*buffer[4]) / n)/n);*/
	return(out);
}

LFHTEMP double Image<T> ::findGrid_evalroutine3(double *coords, int col, int row,double size){
	double buffer[6];
	double n;
	
	double center[2];
	double value[3];
	double dist;
	memset(buffer,'\0',sizeof(double)*6);
	n=0;
	int x,y,px,py,i;
	
	int mix,max,miy,may;
	
	mix = (int)floor(coords[0]* (1.0f -(-1)/(col-1) -(-1)/(row-1))  +  ((-1)*coords[2])/(col-1) + ((-1)*coords[4])/(row-1));
	miy = (int)floor(coords[1]* (1.0f -(-1)/(col-1) -(-1)/(row-1))  +  ((-1)*coords[3])/(col-1) + ((-1)*coords[5])/(row-1));
	max = (int)floor(coords[0]* (1.0f -(col)/(col-1) -(row)/(row-1))  +  ((col)*coords[2])/(col-1) + ((row)*coords[4])/(row-1));
	may = (int)floor(coords[1]* (1.0f -(col)/(col-1) -(row)/(row-1))  +  ((col)*coords[3])/(col-1) + ((row)*coords[5])/(row-1));
	if (mix <0) mix =0;
	if (miy <0) miy =0;
	if (max >=sizex) max =sizex-1;
	if (may >=sizey) may =sizey-1;
	for(x=mix;x<max;x++){
		for(y=miy;y<may;y++){
			px = (int)floor(0.5f + (col-1)*((x - coords[0]) * (coords[2] - coords[0]) + (y - coords[1]) * (coords[3] - coords[1])) / ((coords[2] - coords[0])*(coords[2] - coords[0]) + (coords[3] - coords[1])*(coords[3] - coords[1])));
			py = (int)floor(0.5f + (row-1)*((x - coords[0]) * (coords[4] - coords[0]) + (y - coords[1]) * (coords[5] - coords[1])) / ((coords[4] - coords[0])*(coords[4] - coords[0]) + (coords[5] - coords[1])*(coords[5] - coords[1])));
			
			if ((px >= 0)&&(py >= 0)&&(px < col)&&(px < row)){
				
				center[0] = coords[0]* (1.0f -((double)px)/(col-1) -((double)py)/(row-1))  +  ((px)*coords[2])/(col-1) + ((py)*coords[4])/(row-1);
				center[1] = coords[1]* (1.0f -((double)px)/(col-1) -((double)py)/(row-1))  +  ((px)*coords[3])/(col-1) + ((py)*coords[5])/(row-1);
				if ((center[0] < 0)||(center[1] < 0)||(center[0] > sizex)||(center[1] > sizey)) return(0);
				
				dist = sqrt((center[0]-x)*(center[0]-x)+(center[1]-y)*(center[1]-y));
				if (dist > size){
					n += 1;
					getPixel(x,y, value);
					for(i=0;i<3;i++){
						buffer[0+2*i] += value[i]; 
						buffer[1+2*i] += value[i]*value[i]; 
					}
				}
			}
		}
	}
	
	
	return(1000000.0f / ((buffer[1]+buffer[3]+buffer[5] - (buffer[0]*buffer[0]+buffer[2]*buffer[2]+buffer[4]*buffer[4]) / n)/n));
}





template<class T>
void Image<T> ::applyLinearityFilter(double radius){
	int x,y,px,py;
	double value[12];
	int tmp, tmp2;
	double dist;
	int irad = floor(radius)+1;
	int k;
	double xi;
	double yi;
	double xxi;
	double yyi;
	double xyi;
	double score;
	double thisvalue;
	double average;
	int n=0;
	for(x=0;x<sizex;x++){
		for(y=0;y<sizey;y++){
			xi=yi=0;
			for(px = MY_MAX(1,x-irad,tmp,tmp2);(px< x+irad) && (px < sizex-1);px++){
				for(py = MY_MAX(1,y-irad,tmp,tmp2);(py< y+irad) && (py < sizey-1);py++){
					dist = (px - x) *(px - x) +(py - y) *(py - y);
					if (dist< radius*radius){
						dist = (1 - dist / (radius*radius));
						getPixel(px+1,py, value);
						getPixel(px,py+1, value+3);
						getPixel(px-1,py, value+6);
						getPixel(px,py-1, value+9);
						thisvalue =0;
						for(k = 0;k<3;k++) thisvalue += fabs(value[k]-value[k+3])+fabs(value[k]-value[k+6])+fabs(value[k]-value[k+9])+fabs(value[k+3]-value[k+6])+fabs(value[k+3]-value[k+9])+fabs(value[k+6]-value[k+9]);
						xi += (px - x) * dist* thisvalue;
						yi += (py - y) * dist* thisvalue;
						xxi += (px - x) * (px - x) * dist* thisvalue;
						yyi += (py - y) * (py - y) * dist* thisvalue;
						xyi += (px - x) * (py - y) * dist* thisvalue;
					}
				}
			}
			score = sqrt(xi*xi *xxi +2*xi*yi*xyi + yi*yi*yyi)/ (radius*radius*radius);
			average = (score + average * n) / (n+1);
			n++;
			if (fabs(score) < average /10){
				shiftColor(px,py,1.0f);
			}
			
		}
	}
}

template<class T>
double* Image<T> ::makecolorDistanceMap(){
	double* out = new double[sizex*sizey];
	int x,y;
	int n;
	double value[6];
	double v=0;
	double vv=0;
	double dist;
	for(x=0;x<sizex;x++){
		for(y=0;y<sizey;y++){
			dist =0; n =0;
			getPixel(x,y, value);
			if (x > 0)      {getPixel(x-1,y, value+3); value[3] -= value[0];value[4] -= value[1];value[5] -= value[2]; dist += sqrt(value[3]*value[3]+value[4]*value[4]+value[5]*value[5]);n++;}
			if (x < sizex-1){getPixel(x+1,y, value+3); value[3] -= value[0];value[4] -= value[1];value[5] -= value[2]; dist += sqrt(value[3]*value[3]+value[4]*value[4]+value[5]*value[5]);n++;}
			if (y > 0)      {getPixel(x,y-1, value+3); value[3] -= value[0];value[4] -= value[1];value[5] -= value[2]; dist += sqrt(value[3]*value[3]+value[4]*value[4]+value[5]*value[5]);n++;}
			if (y < sizey-1){getPixel(x,y+1, value+3); value[3] -= value[0];value[4] -= value[1];value[5] -= value[2]; dist += sqrt(value[3]*value[3]+value[4]*value[4]+value[5]*value[5]);n++;}
			out[x + y* sizex] = dist / n;
			v += dist / n;
			vv += (dist / n)*(dist / n);
		}
	}
	v /= (sizex*sizey);
	vv /= (sizex*sizey);
	
	vv = sqrt(vv - v*v);
	
	for(x=0;x<sizex*sizey;x++) out[x] = (out[x] - v) / vv;
	return(out);
}

template<class T>
double* Image<T> ::makecolorLinearityMap(){
	double* out = new double[sizex*sizey];
	int x,y;
	//	int n[8];
	double value[15];
	double v=0;
	double vv=0;
	double dist[2];
	double tmp[4];
	int i;
	for(x=1;x<sizex-1;x++){
		for(y=1;y<sizey-1;y++){
			memset(dist,'\0',sizeof(double)*2);
			getPixel(x,y, value);
			getPixel(x-1,y, value+3);
			getPixel(x+1,y, value+6);
			getPixel(x,y-1, value+9);
			getPixel(x,y+1, value+12);
			
			for(i=0;i<3;i++) tmp[i] = (value[0+i] - value[3+i]);
			tmp[3] = sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2])+1;
			for(i=0;i<3;i++) tmp[i] = (value[0+i] - value[6+i]);
			tmp[3] += sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
			for(i=0;i<3;i++) tmp[i] = (value[0+i] - value[9+i]);
			dist[0] = sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2])+1;
			for(i=0;i<3;i++) tmp[i] = (value[0+i] - value[12+i]);
			dist[0] += sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
			dist[0] = log(dist[0] / tmp[3]);
			
			for(i=0;i<3;i++) tmp[i] = (value[3+i] - value[9+i]);
			tmp[3] = sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2])+1;
			for(i=0;i<3;i++) tmp[i] = (value[6+i] - value[12+i]);
			tmp[3] += sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
			for(i=0;i<3;i++) tmp[i] = (value[3+i] - value[12+i]);
			dist[1] = sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2])+1;
			for(i=0;i<3;i++) tmp[i] = (value[6+i] - value[9+i]);
			dist[1] += sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
			dist[1] = log(dist[1] / tmp[3]);
			
			tmp[3] = sqrt(dist[0]*dist[0] + dist[1]*dist[1]);
			
			out[x + y* sizex] = tmp[3];
			v += tmp[3];
			vv += tmp[3]*tmp[3];
		}
	}
	v /= ((sizex-2)*(sizey-2));
	vv /= ((sizex-2)*(sizey-2));
	
	vv = sqrt(vv - v*v);
	
	for(x=0;x<sizex*sizey;x++) out[x] = (out[x] - v) / vv;
	return(out);
}



template<class T>
template<class S>
void Image<T>::addImage(Image<S>* other){
	int x,y;
	double buffer[3];
	for(x=0;x<sizex;x++){
		for(y=0;y<sizey;y++){
			other->getPixel(x,y,buffer);
			addPixel(x,y,buffer);
		}
	}
}

template<class T>
void Image<T>::markbymap(double* map, double threshold){
	int x;
	for(x=0;x<sizex*sizey;x++){
		if (map[x] > threshold) shiftColor(x % sizex,x/sizex,1.0f);
	}
}


template<class T>
Image<double>* Image<T>::exponentialblur_x(double coef,bool hasweights){
	Image<double>* out = new Image<double>();
	out->allocateBuffer(this);
	int x,y,z;
	double* buffer = new double[sizex*channels];
	double* upbuffer = new double[channels];
	double pix[32];
	double den;
	for(y=0;y<sizey;y++){
		x = sizex-1;
		memset(buffer + channels*x, '\0', sizeof(double)*channels);
		for(;x>=1;x--){
			getPixel(x,y,pix);
			for(z=0;z<channels;z++) buffer[z+  channels*(x-1)] = coef *(pix[z] + buffer[z+  channels*x]);
		}
		x=0;
		getPixel(0,y,pix);
		for(z=0;z<channels;z++) upbuffer[z] = pix[z];
		if (hasweights){
			for(z=0;z<channels;z++) pix[z] = (upbuffer[z] + buffer[z+x*channels]);
		}else{
			den = (1-coef) / (1 + coef - pow(coef,x+1) - pow(coef,sizex-x)); 
			
			for(z=0;z<channels;z++) pix[z] = (upbuffer[z] + buffer[z+x*channels]) *den;
		}
		out->setPixel(0,y,pix);
		for(x=1;x<sizex;x++){
			getPixel(x,y,pix);
			for(z=0;z<channels;z++) upbuffer[z] = pix[z] + coef * upbuffer[z];
			if (hasweights){
				for(z=0;z<channels;z++) pix[z] = (upbuffer[z] + buffer[z+x*channels]);
			}else{
				den = (1-coef) / (1 + coef - pow(coef,x+1) - pow(coef,sizex-x)); 
				
				for(z=0;z<channels;z++) pix[z] = (upbuffer[z] + buffer[z+x*channels]) *den;
			}
			out->setPixel(x,y,pix);
		}
	}
	delete[](buffer);
	delete[](upbuffer);
	
	return(out);
}

template<class T>
Image<double>* Image<T>::exponentialblur_y(double coef,bool hasweights){
	Image<double>* out = new Image<double>();
	out->allocateBuffer(this);
	int x,y,z;
	double* buffer = new double[sizey*channels];
	double* upbuffer = new double[channels];
	double pix[32];
	double den;
	
	for(x=0;x<sizex;x++){
		
		y = sizey-1;
		
		memset(buffer + channels*y, '\0', sizeof(double)*channels);
		
		for(;y>=1;y--){
			getPixel(x,y,pix);
			for(z=0;z<channels;z++) buffer[z+  channels*(y-1)] = coef *(pix[z] + buffer[z+  channels*y]);
		}
		
		
		y=0;
		getPixel(x,0,pix);
		for(z=0;z<channels;z++) upbuffer[z] = pix[z];
		pix[0] = buffer[channels];
		if (hasweights){
			for(z=0;z<channels;z++) pix[z] = (upbuffer[z] + buffer[z+y*channels]);
		}else{
			den = (1-coef) / (1 + coef - pow(coef,y+1) - pow(coef,sizey-y)); 
			for(z=0;z<channels;z++) pix[z] = (upbuffer[z] + buffer[z+y*channels]) *den ;
		}
		out->setPixel(x,0,pix);
		for(y=1;y<sizey;y++){
			getPixel(x,y,pix);
			for(z=0;z<channels;z++) upbuffer[z] = pix[z] + coef * upbuffer[z];
			if (hasweights){
				for(z=0;z<channels;z++) pix[z] = (upbuffer[z] + buffer[z+y*channels]);
			}else{
				den = (1-coef) / (1 + coef - pow(coef,y+1) - pow(coef,sizey-y)); 
				
				for(z=0;z<channels;z++) pix[z] = (upbuffer[z] + buffer[z+y*channels]) *den ;
			}
			out->setPixel(x,y,pix);
		}
	}
	delete[](buffer);
	delete[](upbuffer);
	return(out);
}

template<class T>
Image<double>* Image<T>::exponentialGradientNorm(double coef,bool hasweights){
	Image<double>* xb = exponentialblur_x(coef,hasweights);
	Image<double>* yb = xb->exponentialblur_y(coef,hasweights);
	Image<double>* out = new Image<double>();
	out->allocateBuffer(this);
	out->initBlack();
	int dist = (int)(1 / (1 - coef));
	double pix[32];
	double pjx[32];
	double pkx[32];
	double plx[32];
	int x,y,z;
	for(y = 0;y<sizey-1;y++){
		for(x = 0;x<sizex-1;x++){
			if ((x < dist) ||(x> sizex-dist-1)) {memset(pix,'\0',sizeof(double)*32);memset(pjx,'\0',sizeof(double)*32);}
			else { yb->getPixel(x-dist,y,pix); yb->getPixel(x+dist,y,pjx);}
			if ((y < dist) ||(y> sizey-dist-1)) {memset(pkx,'\0',sizeof(double)*32);memset(plx,'\0',sizeof(double)*32);}
			else { yb->getPixel(x,y-dist,pkx); yb->getPixel(x,y+dist,plx);}
			for(z=0;z<channels;z++) {
				pjx[z] = (pix[z] - pjx[z]);
				plx[z] = (plx[z] - pkx[z]);
				pix[z] = 3.0f * sqrt(pjx[z]*pjx[z]+plx[z]*plx[z]) * (1-coef) / ((1+coef)*(1 - pow(coef,dist)));
			}
			out->setPixel(x,y,pix);
		}
	}
	delete(xb);
	delete(yb);
	return(out);
}

template<class T>
Image<double>* Image<T>::exponentialMultiGradientNorm(bool hasweights){
	Image<double>* out = new Image<double>();
	out->allocateBuffer(this);
	out->initBlack();
	Image<double>* xb[4];
	xb[0]->exponentialGradientNorm(0.25,hasweights);
	xb[1]->exponentialGradientNorm(0.5,hasweights);
	xb[2]->exponentialGradientNorm(0.75,hasweights);
	xb[3]->exponentialGradientNorm(0.875,hasweights);
	int x,y;
	for(y = 0;y<sizey-1;y++){
		for(x = 0;x<sizex-1;x++){
			
		}
	}
	return(out);
	
}

template<class T>
void Image<T>::exponentialblur(double coef,bool hasweights){
	int x,y,z;
	x = (sizex> sizey) ? sizex:sizey;
	double* buffer = new double[x*channels];
	double* upbuffer = new double[channels];
	double pix[32];
	double den =1.0f;
	for(x=0;x<sizex;x++){
		y = sizey-1;
		
		memset(buffer + channels*y, '\0', sizeof(double)*channels);
		
		for(;y>=1;y--){
			getPixel(x,y,pix);
			for(z=0;z<channels;z++) buffer[z+  channels*(y-1)] = coef *(pix[z] + buffer[z+  channels*y]);
		}
		
		getPixel(x,0,pix);
		for(z=0;z<channels;z++) upbuffer[z] = pix[z];
		y=0;
		
		if (hasweights){
			for(z=0;z<channels;z++) pix[z] = (upbuffer[z] + buffer[z+y*channels]);
		}else{
			den = (1-coef) / (1 + coef - pow(coef,y+1) - pow(coef,sizey-y)); 
			for(z=0;z<channels;z++) pix[z] = (upbuffer[z] + buffer[z+y*channels]) *den ;
		}
		
		setPixel(x,y,pix);
		
		for(y=1;y<sizey;y++){
			getPixel(x,y,pix);
			for(z=0;z<channels;z++) upbuffer[z] = pix[z] + coef * upbuffer[z];
			
			if (hasweights){
				for(z=0;z<channels;z++) pix[z] = (upbuffer[z] + buffer[z+y*channels]);
			}else{
				den = (1-coef) / (1 + coef - pow(coef,y+1) - pow(coef,sizey-y)); 
				for(z=0;z<channels;z++) pix[z] = (upbuffer[z] + buffer[z+y*channels]) *den ;
			}
			setPixel(x,y,pix);
		}
	}
	for(y=0;y<sizey;y++){
		x = sizex-1;
		memset(buffer + channels*x, '\0', sizeof(double)*channels);
		for(;x>=1;x--){
			getPixel(x,y,pix);
			for(z=0;z<channels;z++) buffer[z+  channels*(x-1)] = coef *(pix[z] + buffer[z+  channels*x]);
		}
		x=0;
		getPixel(x,y,pix);
		
		for(z=0;z<channels;z++) upbuffer[z] = pix[z];
		
		if (hasweights){
			for(z=0;z<channels;z++) pix[z] = (upbuffer[z] + buffer[z+x*channels]);
		}else{
			
			den = (1-coef) / (1 + coef - pow(coef,x+1) - pow(coef,sizex-x)); 
			for(z=0;z<channels;z++) pix[z] = (upbuffer[z] + buffer[z+x*channels]) *den;
		}
		setPixel(x,y,pix);
		for(x=1;x<sizex;x++){
			getPixel(x,y,pix);
			
			for(z=0;z<channels;z++) upbuffer[z] = pix[z] + coef * upbuffer[z];
			if (hasweights){
				for(z=0;z<channels;z++) pix[z] = (upbuffer[z] + buffer[z+x*channels]);
			}else{
				den = (1-coef) / (1 + coef - pow(coef,x+1) - pow(coef,sizex-x)); 
				for(z=0;z<channels;z++) pix[z] = (upbuffer[z] + buffer[z+x*channels]) *den;
			}
			setPixel(x,y,pix);
		}
	}
	delete[](buffer);
	delete[](upbuffer);
}


template<class T>
void Image<T>::exponentialblur(double coef, unsigned int filter,bool hasweights){
	int x,y,z;
	x = (sizex> sizey) ? sizex:sizey;
	double* buffer = new double[x*channels];
	double* upbuffer = new double[channels];
	double pix[32];
	double den =1.0f;
	for(x=0;x<sizex;x++){
		y = sizey-1;
		
		memset(buffer + channels*y, '\0', sizeof(double)*channels);
		
		for(;y>=1;y--){
			getPixel(x,y,pix);
			for(z=0;z<channels;z++) buffer[z+  channels*(y-1)] = coef *(pix[z] + buffer[z+  channels*y]);
		}
		
		getPixel(x,0,pix);
		for(z=0;z<channels;z++) upbuffer[z] = pix[z];
		y=0;
		
		if (hasweights){
			z=0;
			den = (pix[z] > 0.0f) ? (upbuffer[z] + buffer[z+y*channels]) / pix[z] : 0.0f;
			pix[z] = (upbuffer[z] + buffer[z+y*channels]);
			for(z=1;z<channels;z++) pix[z] = ((filter >> z) &1) ? (upbuffer[z] + buffer[z+y*channels]) : pix[z] * den;
		}else{
			den = (1-coef) / (1 + coef - pow(coef,y+1) - pow(coef,sizey-y)); 
			for(z=0;z<channels;z++) if ((filter >> z) &1) pix[z] = (upbuffer[z] + buffer[z+y*channels]) *den ;
		}
		
		setPixel(x,y,pix);
		
		for(y=1;y<sizey;y++){
			getPixel(x,y,pix);
			for(z=0;z<channels;z++) upbuffer[z] = pix[z] + coef * upbuffer[z];
			
			if (hasweights){
				z=0;
				den = (pix[z] > 0.0f) ? (upbuffer[z] + buffer[z+y*channels]) / pix[z] : 0.0f;
				pix[z] = (upbuffer[z] + buffer[z+y*channels]);
				for(z=1;z<channels;z++) pix[z] = ((filter >> z) &1) ? (upbuffer[z] + buffer[z+y*channels]) : pix[z] * den;
			}else{
				den = (1-coef) / (1 + coef - pow(coef,y+1) - pow(coef,sizey-y)); 
				for(z=0;z<channels;z++) if ((filter >> z) &1) pix[z] = (upbuffer[z] + buffer[z+y*channels]) *den ;
			}
			setPixel(x,y,pix);
		}
	}
	for(y=0;y<sizey;y++){
		x = sizex-1;
		memset(buffer + channels*x, '\0', sizeof(double)*channels);
		for(;x>=1;x--){
			getPixel(x,y,pix);
			for(z=0;z<channels;z++) buffer[z+  channels*(x-1)] = coef *(pix[z] + buffer[z+  channels*x]);
		}
		x=0;
		getPixel(x,y,pix);
		
		for(z=0;z<channels;z++) upbuffer[z] = pix[z];
		
		if (hasweights){
			z=0;
			den = (pix[z] > 0.0f) ? (upbuffer[z] + buffer[z+x*channels]) / pix[z] : 0.0f;
			pix[z] = (upbuffer[z] + buffer[z+x*channels]);
			for(z=1;z<channels;z++) pix[z] = ((filter >> z) &1) ? (upbuffer[z] + buffer[z+x*channels]) : pix[z] * den;
		}else{
			
			den = (1-coef) / (1 + coef - pow(coef,x+1) - pow(coef,sizex-x)); 
			for(z=0;z<channels;z++) if ((filter >> z) &1) pix[z] = (upbuffer[z] + buffer[z+x*channels]) *den;
		}
		setPixel(x,y,pix);
		for(x=1;x<sizex;x++){
			getPixel(x,y,pix);
			
			for(z=0;z<channels;z++) upbuffer[z] = pix[z] + coef * upbuffer[z];
			if (hasweights){
				z=0;
				den = (pix[z] > 0.0f) ? (upbuffer[z] + buffer[z+x*channels]) / pix[z] : 0.0f;
				pix[z] = (upbuffer[z] + buffer[z+x*channels]);
				for(z=1;z<channels;z++) pix[z] = ((filter >> z) &1) ? (upbuffer[z] + buffer[z+x*channels]) : pix[z] * den;
			}else{
				den = (1-coef) / (1 + coef - pow(coef,x+1) - pow(coef,sizex-x)); 
				for(z=0;z<channels;z++) if ((filter >> z) &1) pix[z] = (upbuffer[z] + buffer[z+x*channels]) *den;
			}
			setPixel(x,y,pix);
		}
	}
	delete[](buffer);
	delete[](upbuffer);
}

// assumes nbaperture >=6;
/*
 template<class T>
 Image<double>* Image<T>::maxApertureGradient(bool hasweights, int nbaperture, double min, double max){
 int i;
 vector<Image<double>*> grd;
 Image<double>* tmp;
 double app;
 for(i=0;i<nbaperture;i++){
 app = min * pow(max/min, i / (nbaperture -1));
 tmp = exponentialblur_x(app, hasweights);
 grd.push_back(tmp->exponentialblur_y(app, hasweights));
 delete(tmp);
 }
 
 int x,y;
 double pix[1024];
 double pox[32];
 double mag[32];
 
 
 LFHPrimitive::PolyThing< LFHPrimitive::Vector<double> > ddpoly;
 LFHPrimitive::PolyThing< LFHPrimitive::Vector<double> > poly;
 
 
 for(y =0;y<sizey;y++){
 for(x =0;x<sizex;x++){
 for(i=0;i<5;i++){
 grd[i]->getPixel(x,y,pix +i* channels);
 }
 memset(pox,'\0',sizeof(double)*32);
 for(;i<nbaperture;i++){
 grd[i]->getPixel(x,y,pix + (i & 7)* channels);
 
 
 
 }
 
 
 }
 }
 return(tmp);
 
 }
 */

template<class T>
void Image<T>::weightrescale(){
	double pix[32];
	int x,y,z;
	for(y=0;y<sizey;y++){
		for(x=0;x<sizex;x++){
			getPixel(x,y,pix);
			if (pix[0] <= 0.0f){
				pix[0] =0.0f;
				for(z=1;z<channels;z++) pix[z] = NAN;
				//	memset(pix, '\0',sizeof(T)*(channels));
			}else for(z=1;z<channels;z++) pix[z] /= pix[0];
			setPixel(x,y,pix);
		}
	}
	
	
}



template<class T>
void Image<T>::propagablur(vector<Image<T>* > im){
	int i,x,y,z,w;
	double pix[32];
	double pout[32];
	double g;
	for(z=0;z<im.size();z++){
		for(y=0;y<im[z]->sizey;y++){
			for(x=0;x<im[z]->sizex;x++){
				im[z]->getPixel(x,y,pix);
				if (pix[0] > 1.0f){
					for(w=1;w<im[0]->channels;w++) pix[w] /= pix[0];
					pix[0] =1.0f;
					im[z]->setPixel(x,y,pix);
				}
			}
		}
	}
	
	Image<double>* tmpim[3];
	for(w=0;w<2;w++) {
		tmpim[w] = new Image<double>();
		tmpim[w]->sizex = im[0]->sizex;
		tmpim[w]->sizey = im[0]->sizey;
		tmpim[w]->channels = im[0]->channels;
		tmpim[w]->allocateBuffer();
	}
	for(i=0;i<50;i++){
		printf("Propablur step %i\n",i);
		for(y=0;y<im[0]->sizey;y++){
			for(x=0;x<im[0]->sizex;x++){
				memset(pout,'\0',sizeof(double)*32);
				im[1]->getPixel(x,y,pix);  for(w=0;w<im[0]->channels;w++) pout[w] += pix[w];
				if (x > 0) {im[0]->getPixel(x-1,y,pix);  for(w=0;w<im[0]->channels;w++) pout[w] += pix[w];}
				if (y > 0) {im[0]->getPixel(x,y-1,pix);  for(w=0;w<im[0]->channels;w++) pout[w] += pix[w];}
				if (x < im[0]->sizex-1) {im[0]->getPixel(x+1,y,pix);  for(w=0;w<im[0]->channels;w++) pout[w] += pix[w];}
				if (y < im[1]->sizey-1) {im[0]->getPixel(x,y+1,pix);  for(w=0;w<im[0]->channels;w++) pout[w] += pix[w];}
				tmpim[1]->setPixel(x, y, pout);
			}
		}
		
		for(z=0;z<im.size();z++){
			tmpim[2] = tmpim[0];
			tmpim[0] = tmpim[1];
			tmpim[1] = tmpim[2];
			for(y=0;y<im[z]->sizey;y++){
				for(x=0;x<im[z]->sizex;x++){
					memset(pout,'\0',sizeof(double)*32);
					if (z <im.size()-1){
						if (z < im.size()-2) {im[z+2]->getPixel(x,y,pix);  for(w=0;w<im[0]->channels;w++) pout[w] += pix[w];}
						if (x > 0) {im[z+1]->getPixel(x-1,y,pix);  for(w=0;w<im[0]->channels;w++) pout[w] += pix[w];}
						if (y > 0) {im[z+1]->getPixel(x,y-1,pix);  for(w=0;w<im[0]->channels;w++) pout[w] += pix[w];}
						if (x < im[0]->sizex-1) {im[z+1]->getPixel(x+1,y,pix);  for(w=0;w<im[0]->channels;w++) pout[w] += pix[w];}
						if (y < im[1]->sizey-1) {im[z+1]->getPixel(x,y+1,pix);  for(w=0;w<im[0]->channels;w++) pout[w] += pix[w];}
					}
					im[z]->getPixel(x,y,pix);  for(w=0;w<im[0]->channels;w++) pout[w] += pix[w];
					tmpim[1]->setPixel(x, y, pout);
					
					
					if (pix[0] < 1.0f){
						tmpim[0]->getPixel(x, y, pout);
						if (pout[0] > 0.0f){
							g = (pout[0] -pix[0]*5)/6.0f;
							if (g > 0){
								if (g + pix[0] > 1.0f) g = 1 - pix[0];
								for(w=1;w<im[0]->channels;w++) pix[w] = pix[w] + g*pout[w] / pout[0] ;
								pix[0] += g;
								im[z]->setPixel(x,y,pix);
							}
						}
					}
					
					
				}
			}		
		}
		
	}
	
	delete(tmpim[0]);
	delete(tmpim[1]);
}

template<class T>
void Image<T>::propagablurexp(vector<Image<T>* > &im){
	int x,y,z,w;
	double coef =0.1f;
	
	for(z=0;z<im.size();z++){
		im[z]->exponentialblur(coef,true);
	}
	
	int channels = im[0]->channels;
	double* buffer = new double[(im.size())*channels];
	double* upbuffer = new double[channels];
	double pix[32];
	
	for(y=0;y<im[0]->sizey;y++){
		for(x=0;x<im[0]->sizex;x++){
			z = im.size()-1;
			
			memset(buffer + channels*z, '\0', sizeof(double)*channels);
			
			for(;z>=1;z--) {
				im[z]->getPixel(x,y,pix);
				for(w=0;w<channels;w++) buffer[w+  channels*(z-1)] = coef *(pix[w] + buffer[w+  channels*z]);
			}
			
			z=0;
			
			im[z]->getPixel(x,y,pix);
			for(w=0;w<channels;w++) upbuffer[w] = pix[w] + coef * upbuffer[w];
			
			for(w=0;w<channels;w++) pix[w] = (upbuffer[w] +  buffer[w+z*channels]);
			
			im[z]->setPixel(x,y,pix);
			
			
			for(z=1;z<im.size();z++){
				im[z]->getPixel(x,y,pix);
				for(w=0;w<channels;w++) upbuffer[w] = pix[w] + coef * upbuffer[w];
				
				for(w=0;w<channels;w++) pix[w] = (upbuffer[w] + buffer[w+z*channels]);
				
				im[z]->setPixel(x,y,pix);
			}
			
			
			
		}
	}
	
}

template<class T>
void Image<T>::weightrescale(vector<Image<T>* > &im){
	int x,y,z,w;
	double pix[32];
	
	for(z=0;z<im.size();z++) for(y=0;y<im[0]->sizey;y++) for(x=0;x<im[0]->sizex;x++){
		im[z]->getPixel(x,y,pix);
		if (pix[0] == 0) memset(pix,'\0',sizeof(double)*im[z]->channels);
		else for(w=1;w<im[z]->channels;w++) pix[w] /= pix[0];
   	    pix[0] =1.0f;
		im[z]->setPixel(x,y,pix);
	}
}

template<class T>
void Image<T>::makeEmptyImage(int n_sizex, int n_sizey, int n_channels, int n_bits){
	sizex = n_sizex;
	sizey = n_sizey;
	channels = n_channels;
	bits = n_bits;
	allocateBuffer();
};

template<class T>
void Image<T>::initBlack(){memset(data,'\0',channels * sizex * sizey * sizeof(T));}

template <class T>
void Image<T>::allocateBuffer(){
	if (data != NULL) delete[](data);
	data = (T*) new char[sizeof(T)*channels * sizex * sizey];
}

template <class T> template<class S>
void Image<T>::allocateBuffer(Image<S>* Toclone){
	sizex = Toclone->sizex;
	sizey = Toclone->sizey;
	channels = Toclone->channels;
	allocateBuffer();
};

template <class T> template <class S>
void Image<T>::allocateBuffer(LFHPrimitive::DataGrid<S,2> &Toclone, int _nbchannels){
	sizex = Toclone.dims[0];
	sizey = Toclone.dims[1];
	channels = _nbchannels;
	allocateBuffer();	
}


template<class T>
void Image<T>::scaleImage(double factor){
	double tmp[16];
	int x,y,i;
	for(x=0;x<sizex;x++){
		for(y=0;y<sizey;y++){
			getPixel(x,y,tmp);
			for(i=0;i<channels;i++){tmp[i] *=factor;}
			setPixel(x,y,tmp);
		}
	}
}

template<class T> template< class S>
void Image<T>::copyFrom(Image<S>* other){
	double tmp[16];
	int x,y;
	delete[](data);
	sizex = other->sizex;
	sizey = other->sizey;
	channels = other->channels;
	allocateBuffer();
	for(x=0;x<sizex;x++){
		for(y=0;y<sizey;y++){
			other->getPixel(x,y,tmp);
			setPixel(x,y,tmp);
		}
	}
}

template<class T>
void Image<T>::fillHistogram(vector<int*> counts){
	int x,y,i;
	double buffer[3];
	for(x=0;x<sizex;x++){
		for(y=0;y<sizey;y++){
			getPixel(x,y,buffer);
			for(i=0;i<channels;i++) counts[i][(int)buffer[i]]++;
		}
	}
	
	
}

template <class T>
void Image<T>::rescale(int x, int y){
	int ox = sizex;
	int oy = sizey;
	T* odata = data;
	sizex = x;
	sizey = y;
	allocateBuffer();
	delete[](odata);
}

template<class T>
Image<double>* Image<T>::makeMorphingInput(){
	Image<double>* outp = new Image<double>();
	outp->sizex = sizex;
	outp->sizey = sizey;
	outp->channels = 3;
	outp->bits = sizeof(double);
	int x,y;
	outp->allocateBuffer();
	double pix[16];
	double pox[16];
	double nei[4][16];
	double tmp, tmp2;
	for(y=0;y<sizey;y++){
		for(x=0;x<sizex;x++){
			getPixel(x,y,pix);
			pox[0] = pix[0];
			if ((x > 0) &&(x<sizex-1)){
				getPixel(x-1,y,nei[0]);
				getPixel(x+1,y,nei[1]);
				tmp = nei[0][0] - nei[1][0];
				tmp2 = tmp*tmp;
			}else tmp2 = 0;
			if ((y > 0) &&(y<sizey-1)){
				getPixel(x,y-1,nei[2]);
				getPixel(x,y+2,nei[3]);
				tmp = nei[0][0] - nei[1][0];
				tmp2 += tmp*tmp;
			}		
			pox[1] = sqrt(tmp2);
			pox[2] = sqrt(tmp2);
			outp->setPixel(x,y,pox);
		}
	}
	outp->blur(4.99f); printf("blur1 done\n");
	outp->blur(4.99f); printf("blur2 done\n");
	outp->blur(4.99f); printf("blur3 done\n");
	outp->blur(4.99f); printf("blur2 done\n");
	outp->blur(4.99f); printf("blur3 done\n");
	return(outp);
}


template<class T>
void Image<T>::blur(double radius){
	Image<T>* clone = new Image<T>();
	memcpy(clone,this,sizeof(Image<T>));
	data = NULL;
	allocateBuffer();
	initBlack();
	int x,y,px,py;
	int irad = (int)floor(radius);
	int iirad;
	double pix[32];
	double factor;
	double totalfactor = 0;
	for(py=-irad;py<=irad;py++){
		for(px=-irad;px<=irad;px++){
			factor = py*py+px*px - radius*radius;
			if (factor <= 0){
				totalfactor += sqrt(-factor);
			}			
		}
	}
	totalfactor = 1.0f / totalfactor;
	for(y=0;y<sizey;y++){
		for(x=0;x<sizex;x++){
			clone->getPixel(x,y,pix);
			for(py= ((y-irad>=0)?-irad:-y) ;(py<=irad) && (y+py<sizey);py++){
				iirad = (int) floor(sqrt(radius*radius-py*py))+1;
				for(px=((x-iirad>=0)?-iirad:-x);(px<=iirad) && (x+px<sizex);px++){
					factor = py*py+px*px - radius*radius;
					if (factor < 0){
						factor = totalfactor * sqrt(-factor);
						addPixel(x+px,y+py,pix,factor);
					}
				}
			}
		}
	}
	
	delete(clone);
}

template<class T>
void Image<T>::eval(double* coor,double* output){
	getPixel((int)coor[0],(int)coor[1],output);
	/*	double pixbuf[4][32];
	 double px = coor[0]*(sizex-1);
	 double py = coor[1]*(sizey-1);
	 int x = floor(px); px -= x;  
	 int y = floor(py); py -= y;
	 if (x >= sizex-1) {x=sizex-2;px =1.0f;}
	 else if (x < 0) {x=0;px =0.0f;}
	 if (y >= sizex-1) {y=sizey-2;py =1.0f;}
	 else if (y < 0) {y=0;py =0.0f;}
	 getPixel(x,y,pixbuf[0]);
	 getPixel(x+1,y,pixbuf[1]);
	 getPixel(x,y+1,pixbuf[2]);
	 getPixel(x+1,y+1,pixbuf[3]);
	 int z;
	 for(z=0;z<channels;z++) output[z] = (pixbuf[0][z] * (1-px) + pixbuf[1][z]* px) *(1-py) + (pixbuf[2][z] * (1-px) + pixbuf[3][z]* px) *py;*/
}

template<class T>
void Image<T>::pixelDistribution(int rect[4], double* _out){
	int x,y,z;
	double pix[32];
	double cum[64];
	memset(cum,'\0',sizeof(double)*64);
	for(y=rect[1]; y <=rect[3];y++) for(x=rect[0]; x <=rect[2];x++) {
		getPixel(x,y,pix);
		if (pix[0] > 0.0f){
			cum[0]+=pix[0];
			cum[1]+=pix[0]*pix[0];
			for(z=1;z<channels;z++){
				cum[z<<1] += pix[0]*pix[z];
				cum[1 + (z<<1)] += pix[0]*pix[z]*pix[z];
			}
		}
	}
	for(z=1;z<channels;z++){
		cum[z<<1] /= cum[0];
		cum[1 + (z<<1)] = sqrt((cum[1 + (z<<1)] /cum[0]) - cum[z<<1]*cum[z<<1]);
	}
	memcpy(_out,cum,sizeof(double)*channels*2);
	
}

template<class T>
void Image<T>::evalDerivatives(double* coor,double* output){
	double pixbuf[4][32];
	double px = coor[0]*(sizex-1);
	double py = coor[1]*(sizey-1);
	int x = (int) floor(px); px -= (double)x;  
	int y = (int) floor(py); py -= (double)y;
	if (x >= sizex-1) {x=sizex-2;px =1.0f;}
	else if (x < 0) {x=0;px =0.0f;}
	if (y >= sizex-1) {y=sizey-2;py =1.0f;}
	else if (y < 0) {y=0;py =0.0f;}
	getPixel(x,y,pixbuf[0]);
	getPixel(x+1,y,pixbuf[1]);
	getPixel(x,y+1,pixbuf[2]);
	getPixel(x+1,y+1,pixbuf[3]);
	int z;
	for(z=0;z<channels;z++){
		output[0+2*z] = (pixbuf[1][z] * (1-py) + pixbuf[3][z]* py) - (pixbuf[0][z] * (1-py) + pixbuf[2][z]* py);
		output[1+2*z] = (pixbuf[2][z] * (1-px) + pixbuf[3][z]* px) - (pixbuf[0][z] * (1-px) + pixbuf[1][z]* px);
	}
}


template<class T>
void Image<T>::propagablur(){
	double pix[32];
	double pout[32];
	int x,y,z,i;
	double g;
	Image<double> *tmpim = new Image<double>();
	tmpim->sizex = sizex;
	tmpim->sizey = sizey;
	tmpim->channels = channels;
	tmpim->allocateBuffer();
	
	for(y=0;y<sizey;y++){
		for(x=0;x<sizex;x++){
			getPixel(x,y,pix);
			if (pix[0] >= 1.0f){
				for(z=1;z<channels;z++) pix[z] /= pix[0];
				pix[0] = 1.0f;
				setPixel(x,y,pix);
			}
		}
	}
	
	for(i=0;i<10;i++){
		for(y=0;y<sizey;y++){
			for(x=0;x<sizex;x++){
				memset(pout,'\0',sizeof(double)*32);
				if (x > 0) {
					getPixel(x-1,y,pix);
					pout[0] += pix[0];
					for(z=1;z<channels;z++) pout[z] += pix[z];
				}
				if (y > 0) {
					getPixel(x,y-1,pix);
					pout[0] += pix[0];
					for(z=1;z<channels;z++) pout[z] += pix[z];
				}
				if (x < sizex-1) {
					getPixel(x+1,y,pix);
					pout[0] += pix[0];
					for(z=1;z<channels;z++) pout[z] += pix[z];
				}
				if (y < sizey-1) {
					getPixel(x,y+1,pix);
					pout[0] += pix[0];
					for(z=1;z<channels;z++) pout[z] += pix[z];
				}
				tmpim->setPixel(x,y,pout);
			}
		}
		for(y=0;y<sizey;y++){
			for(x=0;x<sizex;x++){
				getPixel(x,y,pix);
				if (pix[0] < 1.0f){
					tmpim->getPixel(x,y,pout);
					if (pout[0] > 0.0f){
						g = (pout[0] -pix[0]*3)/4.0f;
						if (g > 0){
							if (g + pix[0] > 1.0f) g = 1 - pix[0];
							for(z=1;z<channels;z++) pix[z] = pix[z] + g*pout[z] / pout[0] ;
							pix[0] += g;
							setPixel(x,y,pix);
						}
					}
				}
			}
		}		 
	}
	delete(tmpim);
}

template<class T> template<class S>
void Image<T>::generateMorphedImage(int size, Image<double>* map, Image<S>* &outim){
	Image<double> *tmpim = new Image<double>();
	tmpim->sizex = size;
	tmpim->sizey = size;
	tmpim->channels = channels + 1;
	tmpim->allocateBuffer();
	tmpim->initBlack();
	double pix[32];
	double pmx[32];
	double pox[32];
	double w;
	int x,y,z;
	int px,py;
	for(y=0;y<sizey;y++){
		for(x=0;x<sizex;x++){
			getPixel(x,y,pix);
			map->getPixel(x,y,pmx);
			pmx[0] *= size-1;
			if (pmx[0] >= size-1) {px = size-2; pmx[0] =1.0f;}
			else if (pmx[0] <=0) {px = 0; pmx[0] =0.0f;}
			else {px = (int)floor(pmx[0]); pmx[0] -= px;}
			pmx[1] *= size-1;
			if (pmx[1] >= size-1) {py = size-2; pmx[1] =1.0f;}
			else if (pmx[1] <=0) {py = 0; pmx[1] =0.0f;}
			else {py = (int)floor(pmx[1]); pmx[1] -= py;}
			
			tmpim->getPixel(px,py,pox); w = (1.0f-pmx[0])* ( 1.0f- pmx[1]);
			pox[0] += w; for(z=0;z<channels;z++) pox[1+z] += pix[z] * w;
			tmpim->setPixel(px,py,pox);
			
			tmpim->getPixel(px+1,py,pox); w = pmx[0]* ( 1.0f- pmx[1]);
			pox[0] += w; for(z=0;z<channels;z++) pox[1+z] += pix[z] * w;
			tmpim->setPixel(px+1,py,pox);
			
			tmpim->getPixel(px,py+1,pox); w = (1.0f-pmx[0])* pmx[1];
			pox[0] += w; for(z=0;z<channels;z++) pox[1+z] += pix[z] * w;
			tmpim->setPixel(px,py+1,pox);
			
			tmpim->getPixel(px+1,py+1,pox); w = pmx[0]* pmx[1];
			pox[0] += w; for(z=0;z<channels;z++) pox[1+z] += pix[z] * w;
			tmpim->setPixel(px+1,py+1,pox);
			
		}
	}
	
	tmpim->propagablur();
	
	outim = new Image<S>(); 
	outim->sizex = size;
	outim->sizey = size;
	outim->channels = channels;
	outim->allocateBuffer();
	for(y=0;y<size;y++){
		for(x=0;x<size;x++){
			tmpim->getPixel(x,y,pox);
			if (pox[0] > 0){
				for(z=0;z<channels;z++) pix[z] = pox[z+1] / pox[0];
			} else memset(pix,'\0',sizeof(double)*channels);
			outim->setPixel(x,y,pix);
		}
	}
	
	delete(tmpim);
}
template<class T>  template<class S>
void Image<T>::generateMorphedImage(int size, Image<double>* map, Image<S>* &outim, int* rect){
	Image<double> *tmpim = new Image<double>();
	tmpim->sizex = size;
	tmpim->sizey = size;
	tmpim->channels = channels + 1;
	tmpim->allocateBuffer();
	tmpim->initBlack();
	double pix[32];
	double pmx[32];
	double pox[32];
	double w;
	int x,y,z;
	int px,py;
	for(y=0;y<=rect[3]-rect[1];y++){
		for(x=0;x<=rect[2]-rect[0];x++){
			getPixel(x+rect[0],y+rect[1],pix);
			map->getPixel(x,y,pmx);
			pmx[0] =pmx[1];
			pmx[1] = pmx[2];
			pmx[0] *= size-1;
			if (pmx[0] >= size-1) {px = size-2; pmx[0] =1.0f;}
			else if (pmx[0] <=0) {px = 0; pmx[0] =0.0f;}
			else {px = (int)floor(pmx[0]); pmx[0] -= px;}
			pmx[1] *= size-1;
			if (pmx[1] >= size-1) {py = size-2; pmx[1] =1.0f;}
			else if (pmx[1] <=0) {py = 0; pmx[1] =0.0f;}
			else {py = (int)floor(pmx[1]); pmx[1] -= py;}
			
			tmpim->getPixel(px,py,pox); w = (1.0f-pmx[0])* ( 1.0f- pmx[1]);
			pox[0] += w; for(z=0;z<channels;z++) pox[1+z] += pix[z] * w;
			tmpim->setPixel(px,py,pox);
			
			tmpim->getPixel(px+1,py,pox); w = pmx[0]* ( 1.0f- pmx[1]);
			pox[0] += w; for(z=0;z<channels;z++) pox[1+z] += pix[z] * w;
			tmpim->setPixel(px+1,py,pox);
			
			tmpim->getPixel(px,py+1,pox); w = (1.0f-pmx[0])* pmx[1];
			pox[0] += w; for(z=0;z<channels;z++) pox[1+z] += pix[z] * w;
			tmpim->setPixel(px,py+1,pox);
			
			tmpim->getPixel(px+1,py+1,pox); w = pmx[0]* pmx[1];
			pox[0] += w; for(z=0;z<channels;z++) pox[1+z] += pix[z] * w;
			tmpim->setPixel(px+1,py+1,pox);
			
		}
	}
	
	tmpim->propagablur();
	
	outim = new Image<T>(); 
	outim->sizex = size;
	outim->sizey = size;
	outim->channels = channels;
	outim->allocateBuffer();
	for(y=0;y<size;y++){
		for(x=0;x<size;x++){
			tmpim->getPixel(x,y,pox);
			if (pox[0] > 0){
				for(z=0;z<channels;z++) pix[z] = pox[z+1] / pox[0];
			} else memset(pix,'\0',sizeof(double)*channels);
			outim->setPixel(x,y,pix);
		}
	}
	
	delete(tmpim);
}

template<class T>
void Image<T>::makeImageFromChannels(int in_sizex, int in_sizey, vector<Evaluatable*> &srclist, vector<int> &srcchan){
	channels = srclist.size();
	//	Evaluatable* cur = NULL;
	int i;
	
	sizex = in_sizex;
	sizey = in_sizey;
	
	channels = srcchan.size();
	data = new T[sizex*sizey*channels];
	int x,y;
	
	double pix[32];
	double pox[32];
	
	double coor[2];
	for(x=0;x<sizex;x++){
		coor[0] =x;
		for(y=0;y<sizey;y++){
			coor[1] =y;
			for(i=0;i<channels;i++){
				if (srclist[i] == NULL) pix[i] = 0.0f;
				else{
					srclist[i]->eval(coor,pox);
					pix[i] = pox[srcchan[i]];
				}
			}
			setPixel(x,y,pix);
		}
	}
}


template<class T>
void Image<T>::scaleweight(float fact){
	int x,y;
	double pix[32];
	for(y=0;y<sizey;y++) for(x=0;x<sizex;x++){
		getPixel(x,y,pix);
		pix[0] *= fact;
		setPixel(x,y,pix);
	} 
}

template<class T>
int Image<T>::inputsize(){return(2);}
template<class T>
int Image<T>::outputsize(){return(channels);}

template<class T>
template<class S> void Image<T>::groupCovPixel(Image<S>* origred,Image<S>* origreen, int plateId, MultiCover &mcc){
	int x,y,r,c;
	
	CellCover cc;
	//	cc.plate = plateId;
	double pix[3]; 
	double pox[3]; 
	//	Image<double>* vecf;
	//	Image<unsigned char>* outsink;
	int* glabel = new int[sizex*sizey];
	memset(glabel,'\0',sizex*sizey*sizeof(int));
	vector<int> buf;
	int cur = 0;
	int bounds[256][4];
	for(r = 0;r<sizex*sizey;r++){
		getPixel(x = (r %sizex),(y = r/sizex),pix);
		if ((pix[0]  > 0.0f)&&(glabel[x+y*sizex] == 0x00)){
			
			buf.push_back(x);
			buf.push_back(y);
			bounds[cur][0] = x;
			bounds[cur][1] = y;
			bounds[cur][2] = x;
			bounds[cur][3] = y;
			while(buf.size()!= 0){
				y = *(buf.end()-1);buf.pop_back();x = *(buf.end()-1);buf.pop_back();
				getPixel(x,y,pix);
				if  ((pix[0]  > 0.0f)&&(glabel[x+y*sizex] == 0x00)){
					glabel[x+y*sizex] = (cur+1);
					//		setPixel(x,y,pix);
					if (x > 0) {buf.push_back(x-1);buf.push_back(y);}
					if (x < sizex-1) {buf.push_back(x+1);buf.push_back(y);}
					if (y > 0) {buf.push_back(x);buf.push_back(y-1);}
					if (y < sizey-1) {buf.push_back(x);buf.push_back(y+1);}
					if (x< bounds[cur][0]) bounds[cur][0]= x;
					if (y< bounds[cur][1]) bounds[cur][1] =y;
					if (x> bounds[cur][2]) bounds[cur][2]=x;
					if (y> bounds[cur][3]) bounds[cur][3]=y;
				}
			}
			//		printf("did %i\t%i\n",cur,pawa);
			bounds[cur][0]--;
			bounds[cur][1]--;
			bounds[cur][2]++;
			bounds[cur][3]++;
			cur++;
		}
		//if (cur == 5) break;
	}
	
	Image<double>* newimage;
	
	vector<Image<unsigned char> *> vachehe;
	//	Image<double>* vecf;
	Image<unsigned char>* outsink;
	outsink  = new Image<unsigned char>();
	outsink->sizex = 50;
	outsink->sizey = 50;
	outsink->channels =3;
	outsink->allocateBuffer();
	//	Transformed* xform;
	c=0;
	for(r=0;r<cur;r++){
		
		
		if ((bounds[r][0] > 0)&&(bounds[r][1] > 0)&&(bounds[r][2] < sizex-1)&&(bounds[r][3] < sizey-1)&&((bounds[r][2]-bounds[r][0])*(bounds[r][3]-bounds[r][1]) > 250)){
			//		printf("\t(%i,%i)-(%i,%i)\n",bounds[r][0],bounds[r][1],bounds[r][2],bounds[r][3]);
			
			newimage = new Image<double>();
			newimage->sizex = bounds[r][2] - bounds[r][0] +1;
			newimage->sizey = bounds[r][3] - bounds[r][1] +1;
			//	printf("%i,%i\n", newimage->sizex,newimage->sizey);
			newimage->channels = 3;
			newimage->allocateBuffer();
			for(x=bounds[r][0];x<=bounds[r][2];x++){
				for(y=bounds[r][1];y<=bounds[r][3];y++){
					getPixel(x,y,pix);
					if (glabel[x+y*sizex]  == (r+1)) pox[2] = 1.0f - pix[0];
					else {pox[2] = 1.0f;}
					origred->getPixel(x,y,pix);
					pox[0] = pix[0];
					origreen->getPixel(x,y,pix);
					pox[1] = pix[0];
					newimage->setPixel(x-bounds[r][0],y- bounds[r][1], pox);
				}
			}
			
			cc.findCellCover(newimage);
			memcpy(cc.rect, bounds[r],sizeof(int)*4); 
			if (cc.cover.size() > 0){
				
				c++;
				mcc.group.push_back(cc);
				
				/*	vecf = new Image<double>();
				 vecf->sizex = newimage->sizex;
				 vecf->sizey = newimage->sizey;
				 vecf->channels = 3;		
				 
				 vecf->allocateBuffer();*/
				// cc.fillVectField(vecf);
				//	if (cc.cover.size() < 3){
				
				//	xform =  cc.genTransformation(51);
				//	xform->ev = newimage; 
				//	xform->outsize = 3;
				/*	for(y=0;y<newimage->sizey;y++){
				 for(x=0;x<newimage->sizex;x++){
				 newimage->getPixel(x,y,pix);
				 pix[2] =0.0f;
				 newimage->setPixel(x,y,pix);
				 }}
				 vachehe.push_back(Image<unsigned char>::makeImageFromFunction(newimage->sizex,newimage->sizey,newimage));
				 //		vachehe.push_back(newimage);
				 //	 Image<unsigned char>::SaveTiffImage(vachehe, buffer);
				 delete(vachehe[0]);
				 vachehe.clear();*/
				
				//	 delete(outsink);
				//	 delete(vecf);
				//	 delete(xform);
				//	}
			}
			cc.cover.clear();
			delete(newimage);
			
		}
	}
	printf(" has %i out of %i frames\n",mcc.group.size(),cur);
	
	delete[](glabel);
	
}


template<class T>
void Image<T>::drawLetter(char what, int x, int y, int c){
	int i,j;
	x -=8;
	y -=8;
	if (x<0) x=0;
	if (y<0) y=0;
	if (x+15>= sizex) x=sizex-16;
	if (y+15>= sizey) y=sizey-16;
	double pix[32];
	for(j=0;j<16;j++){
		int data = font[((int)what)*16 + 15 - j];
		for(i=15;i>=0;i--){
			if (data & 1){
				getPixel(x+i,y+j,pix);
				pix[c] = 1.0f; 
				setPixel(x+i,y+j,pix);
			}
			data >>= 1;
		}		
	}
}

template<class T>
void Image<T>::drawLine(int x, int y, int x2,int y2, int c){
	int i,k;
	if (y2 < y) {
		i = y; y = y2; y2 =i;
		i = x; x = x2; x2 =i;
	} 
	double pix[32];	
	
	if ((x2 - x <= y2 - y)&&(x - x2 <= y2 - y)){
		if (y == y2){
			// a point
			getPixel(x,y,pix);
			pix[c] = 1.0f; 
			setPixel(x,y,pix);
		}else{
			k = x2 - x;
			for(i=y;i<=y2;i++){
				
				getPixel(x + ((i-y) * (x2-x))/(y2 - y) ,i,pix);
				pix[c] = 1.0f; 
				setPixel(x + ((i-y) * (x2-x))/(y2 - y) ,i,pix);
			}
		}
		
	}else{
		if (x2 < x){
			for(i=x2;i<=x;i++){
				getPixel(i,y2 + ((i-x2) * (y-y2))/(x - x2) ,pix);
				pix[c] = 1.0f; 
				setPixel(i,y2 + ((i-x2) * (y-y2))/(x - x2) ,pix);
			}
		}else{
			for(i=x;i<=x2;i++){
				getPixel(i,y + ((i-x) * (y2-y))/(x2 - x) ,pix);
				pix[c] = 1.0f; 
				setPixel(i,y + ((i-x) * (y2-y))/(x2 - x) ,pix);
			}
		}
		
	}
	
	
}

template<class T>
template<class D, unsigned int NBCHAN>
Image<T>::operator LFHPrimitive::DataGrid<LFHPrimitive::Tuple<D, NBCHAN>,2> (){
	LFHPrimitive::DataGrid<LFHPrimitive::Tuple<D, NBCHAN>,2> _out;
	unsigned int coor[2];
	
	coor[0] = sizex;
	coor[1] = sizey;
	_out.setSizes(coor);
	
	double pox[32];
	
	int z;
	LFHPrimitive::Tuple<D, NBCHAN> cur;
	
	for(coor[1]=0;coor[1]<sizey;coor[1]++){
		for(coor[0]=0;coor[0]<sizex;coor[0]++){
			getPixel(coor[0],coor[1],pox);
			for(z=0;z<NBCHAN;z++){
				cur[z] = (D) pox[z];
			}
			_out(coor) = cur;
		}
	}
	
	return(_out);
}


template<class T>
template<class D, unsigned int NBCHAN>
LFHPrimitive::DataGrid<LFHPrimitive::Tuple<D, NBCHAN>,2> Image<T>::makedatagrid(){
	LFHPrimitive::DataGrid<LFHPrimitive::Tuple<D, NBCHAN>,2> _out;
	_out.dims[0] = sizex;
	_out.dims[1] = sizey;
	_out.allocate();
	
	double pox[32];
	int coor[2];
	int x,y,z;
	LFHPrimitive::Tuple<D, NBCHAN> cur;
	
	for(coor[1]=0;coor[1]<sizey;coor[1]++){
		for(coor[0]=0;coor[0]<sizex;coor[0]++){
			getPixel(coor[0],coor[1],pox);
			for(z=0;z<NBCHAN;z++){
				cur[z] = (D) pox[z];
			}
			_out(coor) = cur;
		}
	}
	
	return(_out);
}


template<class T>
void Image<T>::addWeigthedPixel(vector<Image<T>* > im, double x, double y, double z, double* pix, double weight){
	int px, py, pz, d;
	double fx,fy,fz,fw;
	double pox[256];
	if (x < 0){px =0; fx =0.0f;}
	else if (x>=im[0]->sizex-1) {px = im[0]->sizex-2;fx =1.0f;}
	else {px = (int)floor(x); fx = x-px;}
	
	if (y < 0){py =0; fy =0.0f;}
	else if (y>=im[0]->sizey-1) {py = im[0]->sizey-2;fy =1.0f;}
	else {py = (int)floor(y); fy = y-py;}
	
	if (z < 0){pz =0; fz =0.0f;}
	else if (z>=im.size()-1) {pz = im.size()-2;fz =1.0f;}
	else {pz = (int)floor(z); fz = z-pz;}
	
	im[pz]->getPixel(px,py,pox);
	fw = weight*(1-fx)*(1-fy)*(1-fz);
	for(d=1;d<im[0]->channels;d++) pox[d] += pix[d-1] * fw;
	pox[0] += fw;
	im[pz]->setPixel(px,py,pox);
	
	im[pz+1]->getPixel(px,py,pox);
	fw = weight*(1-fx)*(1-fy)*fz;
	for(d=1;d<im[0]->channels;d++) pox[d] += pix[d-1] * fw;
	pox[0] += fw;
	im[pz+1]->setPixel(px,py,pox);
	
	im[pz]->getPixel(px,py+1,pox);
	fw = weight*(1-fx)*fy*(1-fz);
	for(d=1;d<im[0]->channels;d++) pox[d] += pix[d-1] * fw;
	pox[0] += fw;
	im[pz]->setPixel(px,py+1,pox);
	
	im[pz+1]->getPixel(px,py+1,pox);fw = weight*(1-fx)*fy*fz;
	for(d=1;d<im[0]->channels;d++) pox[d] += pix[d-1] * fw;
	pox[0] += fw;
	im[pz+1]->setPixel(px,py+1,pox);
	
	im[pz]->getPixel(px+1,py,pox);fw = weight*fx*(1-fy)*(1-fz);
	for(d=1;d<im[0]->channels;d++) pox[d] += pix[d-1] * fw;
	pox[0] += fw;
	im[pz]->setPixel(px+1,py,pox);
	
	im[pz+1]->getPixel(px+1,py,pox);fw = weight*fx*(1-fy)*fz;
	for(d=1;d<im[0]->channels;d++) pox[d] += pix[d-1] * fw;
	pox[0] += fw;
	im[pz+1]->setPixel(px+1,py,pox);
	
	im[pz]->getPixel(px+1,py+1,pox);fw = weight*fx*fy*(1-fz);
	for(d=1;d<im[0]->channels;d++) pox[d] += pix[d-1] * fw;
	pox[0] += fw;
	im[pz]->setPixel(px+1,py+1,pox);
	
	im[pz+1]->getPixel(px+1,py+1,pox);fw = weight*fx*fy*fz;
	for(d=1;d<im[0]->channels;d++) pox[d] += pix[d-1] * fw;
	pox[0] += fw;
	im[pz+1]->setPixel(px+1,py+1,pox);
};
