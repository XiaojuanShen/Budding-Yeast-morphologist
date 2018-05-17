/*
 *  Render.cpp
 *  PureMadness
 *
 *  Created by Louis-Francois Handfield on 11/09/12.
 *  Copyright 2012 Louis-Francois Handfield. All rights reserved.
 *
 */

#include "Tasks.h"


	

	Taskscope<TASK_Cummulus>::Taskscope(): contour_plot(0), show(false), mode(0), val(1.0f),clampstd_x(0.0f), ruler(0), clampstd_y(0.0f),ispos(false),scale(1.0f),hasrect(0),histo(false),condd(false),ims(1024),col_x(0),col_y(1),gp_flags(0),gp_scale(0.0f),col_w(0), col_m(0),normal_dens(0), batch_mode(false),quality(false),kernel_smooth_factor(1.0f),drawLabel(false),dotssize(0.0f) {
		penalize[0] = penalize[1] = false;
		exponent[0] = exponent[1]= exponent[2]= exponent[3]=exponent[4]= 1.0f;
		valscale[0] = valscale[1]= valscale[2]= valscale[3]=valscale[4]= 1.0f;
		valshift[0] = valshift[1]= valshift[2]= valshift[3]=valshift[4]= 0.0f;
		ignore = false; // outside data
		nanreplace = getexp =0;
		active_collumns.setSize(2);
		active_collumns[0] = NULL;
		active_collumns[1] = NULL;
		daprogbar.lenght = 20;
	}
	void Taskscope<TASK_Cummulus>::nbaddtoken(char const * const token, int& min, int& max){
		switch(*token){
			case '\0': min =0; max =2; break;
			case 'c': min =1; break;
			case 'n': min =0; break;
			case 'p': min =0; break;
			case 's': min =1; break;
			case 'r': min =2; max =4; break;
			case 'o': min =0; break;
			case 'i': min =0; break;
			case 'd': min =0; break;
			case 'H': min =0; break;
			case 'F': min =1; break;
			case 'f': min =1;max =2; break;
			case 'w': min =0;max =1; break;
			case 'l': min =1;max =256; break;
			case 'g': min =0;max =1; break;
			case 'P': min =1; break;
			case 'S': min =1; break;
			case 'R': min =1; break;
			case 'A': min =1; break;
			case 'N': min =1; break;
			case 'E': min =0; break;
			case 'L': min =0; break;
			case 'C': min =1; break;
			case 'M': min =1; break;
			case 'K': min =1; break;
			case 'B': min =0; break;
			case 'Q': min =0; break;
			case 'O': min =3; break;
			case 'q': min =1; break;
//			case 'X': min =0; max =256; break;
			case 'X': min =0; break;
		}
	}

	void Taskscope<TASK_Cummulus>::store(char* const * token, int nbtoken){ 
		switch(token[0][1]){
			case 'B': batch_mode = true; break;
			case 'C': contour_plot = atoi(token[1]);break;
			case 'n': normal_dens = 1;
				if (token[0][2] == 'm') normal_dens |=2;
				break;
			case 'c': 
				if (token[0][2] == 'x'){
					clampstd_x = atof(token[1]);
					
				}else if (token[0][2] == 'y'){
					clampstd_y = atof(token[1]);
					//	if (clampstd_y == 0) clampstd_y = getEXCOL(token[1]);
				}
				break;
			case 'P': 
				if (token[0][2] == 'x'){
					exponent[0] = atof(token[1]);
				}else if (token[0][2] == 'y'){
					exponent[1] = atof(token[1]);
				}else if (token[0][2] == 'w'){
					exponent[2] = atof(token[1]);
				}else if (token[0][2] == 'l'){
					exponent[3] = atof(token[1]);
				}
				break;
			case 'F': 
				if (token[0][2] == 'x'){
					valscale[0] = atof(token[1]);
				}else if (token[0][2] == 'y'){
					valscale[1] = atof(token[1]);
				}else if (token[0][2] == 'w'){
					valscale[2] = atof(token[1]);
				}else if (token[0][2] == 'l'){
					valscale[3] = atof(token[1]);
				}
				break;
			case 'A': 
				if (token[0][2] == 'x'){
					valshift[0] = atof(token[1]);
				}else if (token[0][2] == 'y'){
					valshift[1] = atof(token[1]);
				}else if (token[0][2] == 'w'){
					valshift[2] = atof(token[1]);
				} else if (token[0][2] == 'l'){
					valshift[3] = atof(token[1]);
				} 
				break;
			case 'E': 
				if (token[0][2] == 'x'){
					getexp |= 1;
				}else if (token[0][2] == 'y'){
					getexp |= 2;
				}else if (token[0][2] == 'w'){
					getexp |= 4;
				}else if (token[0][2] == 'l'){
					getexp |= 8;
				}
				break;
			case 'L': 
				if (token[0][2] == 'x'){
					getexp |= 16;
				}else if (token[0][2] == 'y'){
					getexp |= 32;
				}else if (token[0][2] == 'w'){
					getexp |= 64;
				}else if (token[0][2] == 'l'){
					getexp |= 128;
				}
				break;
			case 'N': 
				if (token[0][2] == 'x'){
					nanreplace |= 1; nanval[0] = atof(token[1]);
				}else if (token[0][2] == 'y'){
					nanreplace |= 2; nanval[1] = atof(token[1]);
				}else if (token[0][2] == 'w'){
					nanreplace |= 4; nanval[2] = atof(token[1]);
				}else if (token[0][2] == 'l'){
					nanreplace |= 8; nanval[3] = atof(token[1]);
				}
				break;
			case 'g': 
				if (token[0][2] == 'y') gp_flags |= 2;
				else gp_flags |= 1;
				gp_nbsample =(nbtoken == 0) ? 0 : atoi(token[1]);
				break;	
			case 'Q': quality = true; break;
			case 'q': quality = true; dotssize = atof(token[1]); break;

			case 'p': ispos = true; break;
			case 's': scale = atof(token[1]); break;
			case 'r': 
				if (token[0][2] == '\0'){
					hasrect =15;
					rect[0] = atof(token[1]);
					rect[1] = atof(token[2]);
					rect[2] = atof(token[3]);
					rect[3] = atof(token[4]);
				}else if (token[0][2] =='x'){
					if (token[0][3] =='i'){
						hasrect |=1;
						rect[0] = atof(token[1]);
					}else if (token[0][3] =='a'){
						hasrect |=4;
						rect[1] = atof(token[1]);
					}else{
						hasrect |=5;
						rect[0] = atof(token[1]);
						rect[1] = atof(token[2]);
					}
				}else if (token[0][2] =='y'){
					if (token[0][3] =='i'){
						hasrect |=2;
						rect[2] = atof(token[1]);
					}else if (token[0][3] =='a'){
						hasrect |=8;
						rect[3] = atof(token[1]);
					}else{
						hasrect |=10;
						rect[2] = atof(token[1]);
						rect[3] = atof(token[2]);
					}
				}
				break;
			case 'o': penalize[0] =  penalize[1] = true;break;
			case 'i': 
				ignore =true;
				break;
			case 'H': histo =true; break;
			case 'd': condd= true; break; 
			case 'f':
				if (nbtoken == 1) {
					active_collumns[1] = token[1];
				}else{
					active_collumns[0] = token[1];
					active_collumns[1] = token[2];
				}
				break;
			case 'w':
				col_w = active_collumns.size();
				active_collumns.push_back(token[1]);
				break;
			case 'l': 
				col_m = active_collumns.size();
				active_collumns.push_back(token[1]);
				drawLabel = (token[0][2] == 'R');
				if (nbtoken > 1){
					binsep.setSize(nbtoken-1);
					for(unsigned int i=2;i<=nbtoken;i++) binsep[i-2] = token[i];
				}
				break;
			case 'S': ims= atoi(token[1]); break;
			case 'K': kernel_smooth_factor = atof(token[1]); break;
			case 'R': ruler = atoi(token[1]); break;
			case 'X':
				for(unsigned int i = 1; i < nbtoken;i++) colormap.push_back(token[i]);
				switch(token[0][2]){
					case 'c':  //:
						
						
						
						break;//:
				}
				
				
				break;
		}
		
	}
	int Taskscope<TASK_Cummulus>::defstore(char* const * token, int nbtoken){ 
		unsigned int i,j,k,l;
		char buffer[65536];
		file_in[0] = (nbtoken > 0) ? token[0] : NULL;
		file_in[1] = (nbtoken > 1) ? token[1] : NULL;
		double pix[4];
		
		CharArea<> dachar;dachar.setDefaultStyle();
		Tuple<unsigned int,2> mmmmmm; mmmmmm[0] = ims+14; mmmmmm[1] = ims+14;dachar.initialize(mmmmmm);
		
		
		Tuple<unsigned int,3> coors;
		Vector<KeyElem< double, Tuple<unsigned int,2> > > pixel_order;
		
		Tuple< Tuple<unsigned int, 2u >,8u> neigh;
		
		FILE* out  = stdout ;
		
		TableReader tr = TableReader(file_in[1], active_collumns);
		
		myHashmap<string, unsigned int> col_names;
		char sep;
		
//		printf("Cummulus log file for job within arguments: ");
//		for(i=1;i<argc;i++) printf("%s%c",argv[i], (i == argc-1)? '\n' : ' ' );

		
		vector<int> channels;
		channels.push_back(0);channels.push_back(1);
		
		Tuple<double, 2> val;
		
		Vector<char* > tiff_opt;
		TiffFile tf_out(file_in[0] ? file_in[0] : "./autmptmp.tif", true);
		
		//	tf_out.addopt_Description(tiff_opt, "Cummulus Graph");
		
		
		vector< Tuple<double, 4> > data; 
		
		Tuple<double, 4> topush;
		
	
		
		vector<char*> names;
		
		
		FILE* bin;
		if (col_w == 0) topush[2] = 1.0f;
		if (col_m == 0) topush[3] = 0.0f;
		
		while(tr.nextRow()){
			j = 0;
			if (tr.current_row[0].length() == 0) continue;
			else topush[0] = Madstructs::cleveratof(tr.current_row[0].c_str());
			if (tr.current_row[1].length() == 0) continue;
			else topush[1] = Madstructs::cleveratof(tr.current_row[1].c_str());
			if (col_w) {if (tr.current_row[col_w].length() == 0) continue;
			else topush[2] = Madstructs::cleveratof(tr.current_row[col_w].c_str());}
			if (col_m) if (binsep.size() ==0) {
				for(k=0;k<names.size();k++) if (strcmp(buffer, names[k])==0) break;
				topush[2] = Madstructs::cleveratof(tr.current_row[col_m].c_str());
				if (k == names.size()) names.push_back(Madstructs::cloneString(buffer));
				topush[3] = k;
			}else if (tr.current_row[col_m].length() == 0) continue;
			else topush[3] = Madstructs::cleveratof(tr.current_row[col_m].c_str());
	
			data.push_back(topush);
		}
		
		// transform values:
		
		
		for(i=0;i< data.size();i++){
			
			if  ((getexp & 16)==16) data[i][0] = log(data[i][0]);
			if (exponent[0] != 1.0f) data[i][0] = pow(data[i][0],exponent[0]);
			data[i][0] = data[i][0] * valscale[0] + valshift[0];
			if ((getexp & 1)==1) data[i][0] = exp(data[i][0]);
			if (!(ExCo<double>::isValid(data[i][0]))){
				if (nanreplace & 1) data[i][0] = nanval[0];
				else{data[i] = data[data.size()-1]; i--; data.pop_back(); continue;}	
			} 
			
			
			
			if  ((getexp & 32)==32) {data[i][1] = log(data[i][1]);}
			if (exponent[1] != 1.0f) data[i][1] = pow(data[i][1],exponent[1]);
			data[i][1] = data[i][1] * valscale[1] + valshift[1];
			if ((getexp & 2)==2) data[i][1] = exp(data[i][1]);
			if (!(ExCo<double>::isValid(data[i][1]))){
				if (nanreplace & 2) data[i][1] = nanval[1];
				else{data[i] = data[data.size()-1]; i--; data.pop_back(); continue;}	
			}
			
			if (col_w != 0) {
				if  ((getexp & 64)==64) data[i][2] = log(data[i][2]);
				if (exponent[2] != 1.0f) data[i][2] = pow(data[i][2],exponent[2]);
				data[i][2] = data[i][2] * valscale[2] + valshift[2];
				if ((getexp & 4)==4) data[i][2] = exp(data[i][2]);
				if (!(ExCo<double>::isValid(data[i][2]))){
					if (nanreplace & 4) data[i][2] = nanval[2];
					else{data[i] = data[data.size()-1]; i--; data.pop_back(); continue;}
				}	
			}
			
			/*		if ((col_m > -1)&&(binsep.size() !=0)) {
			 if (ExOp::isValid(data[i][3])){
			 if  ((getexp & 128)==128) data[i][3] = log(data[i][3]);
			 if (exponent[3] != 1.0f) data[i][3] = pow(data[i][3],exponent[3]);
			 data[i][3] = data[i][2] * valscale[3] + valshift[3];
			 if ((getexp & 8)==8) data[i][3] = exp(data[i][3]);
			 }else if (nanreplace & 1) data[i][3] = nanval[3];
			 else{data[i] = data[data.size()-1]; i--; data.pop_back(); continue;}
			 }*/
			
		}
		
		// label binning!
		
		if (col_m == 0) { names.push_back(NULL);}			// find the number of label, 
		else if (binsep.size() !=0) {
			
			for(i=0;i<binsep.size()-1;i++) {
				j = strlen(binsep[i]); buffer[j] = '-';
				memcpy(buffer,binsep[i],sizeof(char)*j);
				strcpy(buffer + j+1,binsep[i+1]);
				names.push_back(cloneString(buffer));
			}
			names.push_back(cloneString("Outside"));
			float* bin_list = new float[binsep.size()];
			for(j=0;j<binsep.size();j++) bin_list[j] = atof(binsep[j]);
			for(i=0;i< data.size();i++){
				for(j=0;j<binsep.size();j++) if (data[i][3] < bin_list[j]) break;
				data[i][3] = (j + binsep.size() - 1) % binsep.size();
			}
			delete[](bin_list);
			printf("Using label %s\n", tr.col_name[col_m].c_str());
		} else printf("Using label %s, found %i labels types\n", tr.col_name[col_m].c_str(), names.size());
		
		
		Madstructs::MomentsGenerator0D* dist = new Madstructs::MomentsGenerator0D(channels);
		
		
		dist->init();
		val[0] = data[0][0];
		val[1] = data[0][1];
		dist->saferecord(&(val[0]));
		for(i=1;i< data.size();i++){
			val[0] = data[i][0];
			val[1] = data[i][1];
		//	if (val[0] < trect[0]) rect[0] =val[0];
		//	else if (val[0] > trect[1]) rect[1] =val[0];
		//	if (val[1] < trect[2]) rect[2] =val[1];
		//	else if (val[1] > trect[3]) rect[3] =val[1];
			dist->saferecord(&(val[0]));
		}
		dist->finish(); 
		
		if ((hasrect & 1) == 0) rect[0] = dist->x[5];
		if ((hasrect & 4) == 0) rect[1] = dist->x[6];	
		if ((hasrect & 2) == 0) rect[2] = dist->x[12];
		if ((hasrect & 8) == 0) rect[3] = dist->x[13];
		
		if (clampstd_x != 0){
			if (rect[0] < dist->x[1] - dist->x[2] * clampstd_x) rect[0] = dist->x[1] - dist->x[2] * clampstd_x;
			if (rect[1] > dist->x[1] + dist->x[2] * clampstd_x) rect[1] = dist->x[1] + dist->x[2] * clampstd_x;
		}
		if (clampstd_y != 0){
			if (rect[2] < dist->x[8] - dist->x[9] * clampstd_y) rect[2] = dist->x[8] - dist->x[9] * clampstd_y;
			if (rect[3] > dist->x[8] + dist->x[9] * clampstd_y) rect[3] = dist->x[8] + dist->x[9] * clampstd_y;
		}
		

		
		if (ispos){
			rect[0] =0.0f;
			rect[2] =0.0f;
		}
		
		for(i=0;i< data.size();i++){
			if ((data[i][0] < rect[0])||(data[i][0] > rect[1])||(data[i][1] < rect[2])||(data[i][1] > rect[3])){
				data[i] = data[data.size()-1];
				i--;
				data.pop_back();
			}
		}
		printf("\n\nDisplay Rect: [%f, %f] x [%f, %f], %i data points in total\n", rect[0], rect[1], rect[2], rect[3], data.size()); fflush(stdout);
		
		GaussianDistribution<2>* correl = new GaussianDistribution<2>[names.size() + 1];
		
		
		for (k =0; k<names.size()+1;k++) correl[k].EMinit();
		daprogbar.start("Compiling Bin Specific Spreads");
		for(i=0;i< data.size();daprogbar.update(((double)i++)/ data.size())){
			val[0] = data[i][0];
			val[1] = data[i][1];
			if ((ignore)&&((val[0] < rect[0])||(val[0] > rect[1])||(val[1] < rect[2])||(val[1] > rect[3]))) continue ;
			correl[(unsigned int)data[i][3]].EMregist(val, data[i][2]);
			correl[names.size()].EMregist(val, data[i][2]); 
		}daprogbar.finish();
		
		if (names.size() > 1){
			for (k =0; k<names.size();k++){
				printf("For label %s:\n", names[k]);
				correl[k].EMfinit();
				correl[k].show();
			}
			printf("For all label:\n");
			correl[k].EMfinit();
			correl[k].show();
		}else{
			correl[0].EMfinit();
			correl[0].show();
		}
		fflush(stdout);
		
		
		DataGrid<double , 2>* daframes = new DataGrid<double , 2>[names.size() == 1 ? 1 : names.size()+1];
		
		DataGrid<double , 2> daframe_copy;
		DataGrid<double , 3> da_colored_frame;
		
		
		DataGrid<double , 1>* daframes_hist = (histo) ? new DataGrid<double , 1>[2* (names.size() == 1 ? 1 : names.size()+1) ] : NULL;
		
		
		Tuple<unsigned int,1> cor_h;
		Tuple<unsigned int,2> cor;
		KeyElem<Tuple<double,2> , double > valin;
		KeyElem<Tuple<double,1> , double > valin_h;
		
		cor[0] = ims;
		cor[1] = ims;
		cor_h[0] = ims; 
		for (k =0; k<names.size()+1;k++){daframes[k].setSizes(cor);ExOp::toZero(daframes[k]);}
		if (histo) for (k =0; k<(names.size()+1)*2;k++) {daframes_hist[k].setSizes(cor_h);ExOp::toZero(daframes_hist[k]);}
		// remove outside pts
		
		
		if (names.size() == 1){
			for(i=0;i< data.size();i++){
				valin.k[0] = (ims -1)*(data[i][0] - rect[0])/(rect[1] - rect[0]);
				valin.k[1] = (ims -1)*(data[i][1] - rect[2])/(rect[3] - rect[2]);
				valin.d = data[i][2];
				daframes[0] += valin;
				if (histo) {valin_h.d = valin_h.d; valin_h.k[0] = valin.k[0]; daframes_hist[0] += valin_h; valin_h.k[0] = valin.k[1]; daframes_hist[1] += valin_h;}
			}
		}else{
			for(i=0;i< data.size();i++){
				valin.k[0] = (ims -1)*(data[i][0] - rect[0])/(rect[1] - rect[0]);
				valin.k[1] = (ims -1)*(data[i][1] - rect[2])/(rect[3] - rect[2]);
				valin.d = data[i][2];
				daframes[(unsigned int) data[i][3]] += valin;
				daframes[names.size()] += valin;
				if (histo) {
					valin_h.d = valin.d;
					valin_h.k[0] = valin.k[0]; daframes_hist[(((unsigned int) data[i][3])<<1)] += valin_h; daframes_hist[((names.size())<<1)] += valin_h;
				valin_h.k[0] = valin.k[1]; daframes_hist[((((unsigned int) data[i][3])<<1)| 1)] += valin_h; daframes_hist[(((names.size())<<1)|1)] += valin_h;}
			}
		}
		
		
		{//:
			
			for (k =0; k< ((names.size() == 1) ? 1 : names.size()+1);k++) {
				Tuple<double,2> dascale;
				dascale[0] = kernel_smooth_factor * (ims / (rect[1] - rect[0])) * pow(correl[k].weight, -1.0f / 6.0f);
				dascale[1] = kernel_smooth_factor * (ims / (rect[3] - rect[2])) * pow( correl[k].weight, -1.0f / 6.0f); 
				GaussianDistribution<2> dadamat = correl[k];
				dadamat.mean[0] =0.0f;
				dadamat.mean[1] =0.0f;
				
								
				dadamat *= dascale;
				
				if (contour_plot) daframe_copy = daframes[k];
				
				if (quality) {
					if (dotssize == 0.0f) daframes[k].blur_zeroborder(dadamat);
					else daframes[k].convolvecircle_zeroborder(dotssize);
				} else daframes[k].blur_crude(dadamat);
				
				if (normal_dens){ 
					if (normal_dens == 1) daframes[k] *= scale / correl[k].weight;
					else{
						{ //:
							double max = 0;
							Tuple<unsigned int ,2 > coor;
							for(coor[1]=0;coor[1]<ims;coor[1]++){
								for(coor[0]=0;coor[0]<ims;coor[0]++){
									if (max < daframes[k](coor)) max = daframes[k](coor);
								}
							}
							if (max != 0.0f) daframes[k] *= scale / max;
							else daframes[k] *= scale / correl[k].weight;
						} //:
					}
				}else daframes[k] *= scale; 
				
				//	tf_out.put(daframes[k], (unsigned char) 0, (unsigned char) 255);
				
		//		if (ruler){
					
		//			LFHPrimitive::addRuler(daframes[k]);
				
		//		}
				
	
				if (contour_plot){
				//	daframe_copy.exponentialBlur(2.0f, false);
					
					unsigned char tmp;
					DataGrid<unsigned int, 2> daslices = daframes[k].makeSlicesIndexes(contour_plot);
					DataGrid<unsigned char, 3u> dacolored;
					Tuple<unsigned int,3u> coor;
					
					coor[2] = daframes[k].dims[1];
					coor[1] = daframes[k].dims[0];
					coor[0] =3;ExOp::show(coor);fflush(stdout);dacolored.setSizes(coor); 
					DataGrid<unsigned int, 2>::KeyIterator dasl_ite = daslices.getKeyIterator();
					if (dasl_ite.first()) do{
						coor[2] = dasl_ite()[1];
						coor[1] = dasl_ite()[0];
						
						l = daslices.get_indirectNeightbor(dasl_ite(),neigh);
						unsigned int nei_count=0;
						for(l--;l!= 0xFFFFFFFF;l--) if (daslices(dasl_ite()) > daslices(neigh[l])) nei_count++; 
						if (nei_count > 0){
							coor[0] =0; dacolored(coor) = 255;
							coor[0] =1; dacolored(coor) = 0;
							coor[0] =2; dacolored(coor) = 0;
						}else{
							tmp = (daframe_copy(dasl_ite()) > 1.0f) ? 0 : 255 - ((unsigned char) (daframe_copy(dasl_ite()) * 255.0f));
							coor[0] =0; dacolored(coor) = tmp;
							coor[0] =1; dacolored(coor) = tmp;
							coor[0] =2; dacolored(coor) = tmp;
						}
						
					} while(dasl_ite.next());
					dachar.setAxes(rect[0],rect[2],rect[1],rect[3]);
					dachar.overwriteFrame(dacolored);
					tf_out.put(dachar.image, (unsigned char) 0.0f, (unsigned char) 255);
				} else {
					dachar.setAxes(rect[0],rect[2],rect[1],rect[3]);
				//	dachar.drawFrame(daframes[k], 0.0f, 255.0f);
					tf_out.put(dachar.image, (unsigned char) 0, (unsigned char) 255);
				}
				
				if (histo){ // make HISTOGRAMS!!
					k = k << 1;
					GaussianDistribution<1> dadamat_h;
					Tuple<double,1> dascale_h;
					cor_h[0] =0;
					dadamat_h = dadamat.makeSubGaussian(cor_h);
					dascale_h[0] = kernel_smooth_factor * (ims / (rect[1] - rect[0])) * pow(0.75f*correl[(k>>1)].weight, -1.0f / 5.0f); // silverman thumb rule!
					
					dadamat_h *= dascale_h;
					dadamat_h.show();
					if (quality) daframes_hist[k].blur_zeroborder(dadamat_h);
					else daframes_hist[k].blur_crude(dadamat_h);
					
					ExOp::toZero(dascale[0]);
					for(cor_h[0] = 0;cor_h[0]< ims;cor_h[0]++) if (daframes_hist[k](cor_h) > dascale[0]) dascale[0] = daframes_hist[k](cor_h);
					
					
					for(cor_h[0] = 0;cor_h[0]< ims;cor_h[0]++){
						dascale[1] = (daframes_hist[k](cor_h) * (ims -1)) / dascale[0];
						
						for(cor[0] =cor_h[0], cor[1] = 0; ((ims -cor[1]> ((unsigned int )dascale[1]+2)))&&(cor[1] < ims);cor[1]++) daframes[(k>>1)](cor) = 0.0f;
						if ((cor[1] < ims)&&(dascale[1]>=0.0f)) daframes[(k>>1)](cor) = (dascale[1] +2.0f + cor[1])  - ims;
						for(cor[1]++; cor[1] < ims;cor[1]++) daframes[(k>>1)](cor) = 1.0f;
					}
					
					tf_out.put(daframes[(k>>1)], (float) 0.0f, (float) 255.0f);
					
					k = k | 1;
					cor_h[0] =1;
					dadamat_h = dadamat.makeSubGaussian(cor_h);
					dascale_h[0] = kernel_smooth_factor * (ims / (rect[3] - rect[1])) * pow(0.75f*correl[((k>>1))].weight, -1.0f / 5.0f); // silverman thumb rule!
					dadamat_h *= dascale_h;
					if (quality) daframes_hist[k].blur_zeroborder(dadamat_h);
					else daframes_hist[k].blur_crude(dadamat_h);
					
					ExOp::toZero(dascale[0]);
					for(cor_h[0] = 0;cor_h[0]< ims;cor_h[0]++) if (daframes_hist[k](cor_h) > dascale[0]) dascale[0] = daframes_hist[k](cor_h);
					
					for(cor_h[0] = 0;cor_h[0]< ims;cor_h[0]++){
						dascale[1] = (daframes_hist[k](cor_h) * (ims -1)) / dascale[0];
						
						for(cor[0] =ims -1 -cor_h[0], cor[1] = 0; ((ims -cor[1] > ((unsigned int )dascale[1]+2)))&&(cor[1] < ims);cor[1]++) daframes[(k>>1)](cor) = 0.0f;
						if ((cor[1] < ims)&&(dascale[1]>=0.0f))  daframes[(k>>1)](cor) =  (dascale[1] +2.0f + cor[1])  - ims;
						for(cor[1]++; cor[1] < ims;cor[1]++) daframes[(k>>1)](cor) = 1.0f;
					}
					
					tf_out.put(daframes[(k>>1)], (float) 0.0f, (float) 255.0f);
					k = k >> 1;
					

					
				}
			}
		}//:
		
		

		
		exit(0);
		if (colormap.size() >0){
			{
				DataGrid<double , 3> colframe = new DataGrid<double , 3>();
				unsigned int coor[3];
				coor[0] = 3;
				coor[1] = daframes[0].dims[0];
				coor[2] = daframes[0].dims[1];
				colframe.setSizes(coor);
				
				tf_out.put(colframe, (float) 0.0f, (float) 255.0f);
			}//:
		}
		
		
		if (gp_flags != 0){
			{//:
				GaussianProcessProbSimple<1,1> simpleGP;
				
				for(j=0;j<10;j++){
					//simpleGP.EMinit();
					//for(i=0;i< data.size();i++) EMregist(KeyElem<double , double>(data[i][0],valin.k[1]));
					//	printf("step %i: GP likelihood %f\n", j+1, simpleGP.EMfinit());
				}
				
				
				
			}//:
		}
		
		
		int ims = ims;
		Madstructs::WImage* w = new Madstructs::WImage(ims,ims,0);
		Madstructs::WImage* w2;
		w->initBlack();
		Madstructs::WImage* h[2];
		unsigned int coor[1];
		DataGrid< WeightElem< double , 2>, 1 > cdist[3];
		
		double maxh[6];
		if (histo){
			h[0] = new Madstructs::WImage(ims,1,0);
			h[1] = new Madstructs::WImage(ims,1,0);
		}
		
		if (condd){
			coor[0] = ims;
			cdist[0].setSizes(coor);
			cdist[1].setSizes(coor);
			ExOp::toZero(cdist[0]);
			ExOp::toZero(cdist[1]);
		}
		
		dist->init();
		
		for (k =0; k<names.size();k++) {correl[k].EMinit();	}
		
		for(i=0;i< data.size();i++){
			val[0] = data[i][0];
			val[1] = data[i][1];
			//	val[1] /= val[0] + 15.0f;
			
			
			pix[0] = (ims -1)*(val[0] - rect[0])/(rect[1] - rect[0]);
			pix[1] = (ims -1)*(rect[3] - val[1])/(rect[3] - rect[2]);
			if (!(isnan(val[0])||isnan(val[1])||isinf(val[0])||isinf(val[1]))) {
				if ((!ignore)||((pix[0] >= 0.0f)&&(pix[1] >= 0.0f)&&(pix[0] <= ims -1)&&(pix[1] <= ims -1))){
					w->addWeigthedPixel(pix[0] ,pix[1],NULL,scale * data[i][2]);
					dist->saferecord(&(val[0]));
					correl[(unsigned int)data[i][3]].EMregist(val,data[i][2]);
				}
				if (condd){
					coor[0] = (int)pix[1];
					if ((coor[0] >= 0)&&(coor[0]< ims))	cdist[0](coor) += WeightElem<double,2>(pix[0], data[i][2]* (1.0f + coor[0] - pix[1]));
					coor[0]++;
					if ((coor[0] >= 0)&&(coor[0]< ims))	cdist[0](coor) += WeightElem<double,2>(pix[0], data[i][2]*(pix[1] - coor[0]));
					coor[0] = (int)pix[0];
					if ((coor[0] >= 0)&&(coor[0]< ims))	cdist[1](coor) += WeightElem<double,2>(pix[1], data[i][2]*(1.0f + coor[0] - pix[0]));
					coor[0]++;
					if ((coor[0] >= 0)&&(coor[0]< ims))	cdist[1](coor) += WeightElem<double,2>(pix[1], data[i][2]*(pix[0] - coor[0]));
				}
				
				if (histo){
					h[0]->addWeigthedPixel(pix[0],NULL, data[i][2] );
					h[1]->addWeigthedPixel(ims - pix[1] -1,NULL, data[i][2]);
				}
			}
		}
		
		if (ignore){
			printf("==========\nData in clamped region:\n==========\n");
			printf("%s\t%s\n", tr.col_name[col_x].c_str(), tr.col_name[col_y].c_str());
			dist->finish();
			dist->show(out);
			for (k =0; k<names.size();k++){
				if (names.size() > 1) printf("Correlations for label %s:\n", names[k]);
				correl[k].EMfinit();
				correl[k].show();
			}
		}
		//	w->scaleweight(100.0f);
		w2 = w->exponentialblur_x(0.8f);
		delete(w);
		w = w2->exponentialblur_y(0.8f);
		delete(w2);
		
		DataGrid<float, 3> iout;
		
		
		int y,z;
		if (histo){
			maxh[0] = maxh[1] =0.0f;
			for(y=0;y<ims;y++){
				h[0]->getPixel(y,0,pix);
				if (maxh[0] < pix[0]) maxh[0] = pix[0];
				h[1]->getPixel(y,0,pix);
				if (maxh[1] < pix[0]) maxh[1] = pix[0];
			}
		}
		maxh[2] = maxh[3] =0.0f;
		
			coors[0]=1;
			coors[1]=ims;
			coors[2]=ims;
			iout.setSizes(coors);coors[0] = 0;
			for(coors[1]=0;coors[1]<ims;coors[1]++){
				for(coors[2]=0;coors[2]<ims;coors[2]++){
					iout(coors) = w->data[coors[1] + coors[2]*ims]; 
				}
			}
		
		
		tf_out.put(iout,(float) 0, (float) 255.0f);
		//			tf_out.put(iout,tiff_opt);
		if (histo){
			for(coors[1]=0;coors[1]<ims;coors[1]++){
				if (histo){
					h[0]->getPixel(coors[1],0,pix);
					maxh[4] = pix[0] / maxh[0];
					maxh[2] += pix[0] / (dist->x[0]);
					h[1]->getPixel(coors[1],0,pix);
					maxh[5] = pix[0] / maxh[1];
					maxh[3] += pix[0] / (dist->x[0]);
				}
				for(coors[2]=0;coors[2]<ims;coors[2]++){
					pix[2] =pix[0] = (ims - coors[2]-1 < maxh[2] * (ims-1)) ? 1.0f : 0.0f;
					pix[1] = (ims - y< maxh[4] * ims) ? 1.0f : 0.0f;
					//	iout[1]->setPixel(x,y,pix);
					pix[2] =pix[0] = (ims - coors[2]-1< maxh[3] * (ims-1)) ? 1.0f : 0.0f;
					pix[1] = (ims - y< maxh[5] * ims) ? 1.0f : 0.0f;
					//	iout[2]->setPixel(x,y,pix);
				}
				
			}
			
		}
		
		/*
		 
		 if  (condd){
		 cdist[2] = cdist[0].crudeGaussianBlur(ims * 0.1f);
		 iout.push_back(new Madstructs::Image<unsigned char>());
		 l = iout.size()-1;
		 iout[l]->sizex = ims;
		 iout[l]->sizey = ims;
		 iout[l]->channels = 3;
		 iout[l]->allocateBuffer();
		 for(x=0;x<ims;x++){
		 coor[0] =x;
		 maxh[0] = cdist[0](coor).getMean();
		 maxh[1] = -0.5f / cdist[0](coor).getVar();
		 for(y=0;y<ims;y++){
		 pix[0] = y - maxh[0];
		 pix[0] =pix[1] =pix[2] = exp(maxh[1] * pix[0] *pix[0]);
		 iout[l]->setPixel(x,y,pix);
		 }
		 }
		 //	cdist[2] = cdist[1].crudeGaussianBlur(ims * 0.1f);
		 iout.push_back(new Madstructs::Image<unsigned char>());
		 l = iout.size()-1;
		 iout[l]->sizex = ims;
		 iout[l]->sizey = ims;
		 iout[l]->channels = 3;
		 iout[l]->allocateBuffer();
		 for(x=0;x<ims;x++){
		 coor[0] =x;
		 maxh[0] = cdist[1](coor).getMean();
		 maxh[1] = -0.5f / cdist[1](coor).getVar();
		 for(y=0;y<ims;y++){
		 pix[0] = y - maxh[0];
		 pix[0] =pix[1] =pix[2] = exp(maxh[1] * pix[0] *pix[0]);
		 iout[l]->setPixel(y,x,pix);
		 }
		 }
		 
		 }*/
		
		//	system("open -a Preview \"./autmptmp.tif\"");

		
	return 0;
		
		}
	void Taskscope<TASK_Cummulus>::help(){
		printf("Make a Cloud Image.\n");
		printf("\n");
		printf("Optionnal Default arguments(2):\n");
		printf("\t(out) produced tiff image (default = stdout)\n");
		printf("\t(in) 2 collumns table (default = stdin) \n");
		printf("Optionnal Flags\n\n");
		printf("\tChooses the collunms in the tabular delimited file to display data\n First collumn is '1' or 'A', the 28th is 28 or 'AB' and so on:\n\n");

		printf("\t-f (int or string) (int or string) [int or string] : use the collumns x and y (and optionnaly z), if it is a string, a the header colunm name is searched for\n");
		printf("\t-w (int or string): collumn index for weights\n" );
		printf("\t-l[R] (int x=4) (min bin 1) (bin1 bin2 sep) ... (max bin n): collumn index for label, optionnaly followed by bin delimiters\n");
		printf("\t  [R] option: draw label on graph\n");
		printf("\t-z : collumn index for 3rd component) value (z)\n");
		printf("\t-ih : has a header row\n");
		
		printf("\t-L{x,y,z,w,l}: take logarithm (before all other mods)\n");
		printf("\t-P{x,y,z,w,l} (float n) : exponentiate values by 'n'\n");
		printf("\t-F{x,y,z,w,l} (float n) : scale values by\n");
		printf("\t-A{x,y,z,w,l} (float n) : add value 'n'\n");
		printf("\t-E{x,y,z,w,l}: take exponential (after all other mods)\n");
		printf("\t-N{x,y,z,w,l} (float n) : replace NAN by value (rows with any NAN are ignored by default)\n");
		
		
		printf("\nData Related:\n");
		printf("\t-B: Batch mode, input file (or std input) is a list of file names (files which contains data tables that are iterativly opened):\n\n");
		printf("\t-c{x,y} (float n) : Clamp image using 'n' STD (default =3), the input 0 disable the clamping\n");
		printf("\t-p : Clamp to positive quadrant (always)\n");
		printf("\t-s (float n) : scale cumulative values by 'n'\n");
		printf("\t-o : penalize non-gaussian outliers (not yet)\n");
		printf("\t-r (float) (float) (float) (float): clamping rectangle\n");
		printf("\t-r{xy} (float) (float): clamping range\n");
		printf("\t-r{xy}{ia} (float): clamping range (min or max)\n");
		printf("\t-i : ignore data outside clamping for statistics\n");
		printf("\t-d : make conditionnal distributions\n");
		printf("\t-g{x,y} (int x=0): Gaussian Process Fit. the default basis is the datapoint, or 'X' equidistant points\n");
		printf("\t-S (int s) : image size (default 1024)\n");
		printf("\t-Y : use special weight for cells table\n");
		printf("\t-C (int X) : makes and RGB contour plot, with X stripes\n");
		printf("\t-n : normalize density\n");
		printf("\t-nm : normalize density so maximum matches scale parameter (modified with -S flag) \n");
		
		printf("\nRender Related:\n");
		
		printf("\t-O (float x) (float y) (float angle): draw Oblique line\n");
		printf("\t-X : Include axes\n");
		printf("\t-Q : High-Quality flag (slower but better)\n");		
		printf("\t-q (float) : draw circles instead, [turns on -Q flag]\n");		
		printf("\t-H : make histogramms and Cummulative\n");
		printf("\t-K (float scale = 1.0) : density smoothing factor\n");
		printf("\t-De (int color) (int nbstripes) (meanx) (meany) (varx) (vary) (correl): render atop the graph\n");
		printf("\tExtra Graphs:\n\n");
	//	printf("\t-Xc (string '(label name or number)=Colorstr') ... : makes a colored density graph\n");
		
		printf("\t-h : Help \n");
	}

	
	


