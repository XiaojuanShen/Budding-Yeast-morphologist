/*
 * Extraction.cpp
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


Taskscope<TASK_HIDDENMAP_BASE_DATAEXTRACTION>::Taskscope(): red_substract_fact(0.0f), show(false),recursive(false),residual(NULL),out_p(NULL),table_out(false),seg_im_path(NULL),mcv_in(NULL),file_dascp(NULL),distance_thr(0.0f){
	art_k = 4.510662f -0.165461f; // ram...
	art_k *= 1.95248835f - 0.035715438f; // intensity
	art_k *= -0.030004794f - -3.258060922f; // density				
	art_k *= -0.039659578f - -4.744727495f; // edge err/ 
	art_k *= 5900.0f;
	art_k =  0.4404983f / (art_k * 0.5595016f);
	//					printf("%e\t%e\n",art_k);
	//art_k = (0.9f*art_k + 0.4f * 0.00001f) * 0.5f;
	//					printf("%e\t%e\n",art_k, atof("2.474696e-06"));
	art_k *= 10.0f;
	artifact_class[0] = 0.05f;
	artifact_class[1] = 10.0f;
	artifact_class[2] = 0.9f;
	artifact_class[3] = 0.1f;
	artifact_class[4] = 50.0;
	ignore_red_in_confidence =false;
}
void Taskscope<TASK_HIDDENMAP_BASE_DATAEXTRACTION>::nbaddtoken(char const * const token, int& min, int& max){
	switch(*token){
		case '\0': min =3; break;
		case 's': min =1;break;
		case 'c': min =1;break;
		case 'S': min =1;break;
		case 't': min=1;break;
		case 'B': min=1;break;
		case 'm': min=1;break;
		case 'M': min=2;break;
		case 'p': min=1;break;
		case 'R': min=0;break;
			//	case 'A': min=0;break;
		case 'a': min=5;break;
		case 'r': min=0;break;
	}
}
void Taskscope<TASK_HIDDENMAP_BASE_DATAEXTRACTION>::store(char* const * token, int nbtoken){ 
	
	switch(token[0][1]){
		case 'r':ignore_red_in_confidence =true;
		case 'c': file_dascp = token[1]; break;
		case 'a':
			artifact_class[0] = atof(token[1]);
			artifact_class[1] = atof(token[2]);
			artifact_class[2] = atof(token[3]);
			artifact_class[3] = atof(token[4]);
			artifact_class[4] = atof(token[5]);
			break;
		case 's':show = true;
			break;
		case 'o':
			out_p = token[1];
			break;
		case 'p':art_k = atof(token[1]);break;
		case 't':
			table_out = token[1];
			break;
		case 'S': seg_im_path = token[1];break;
		case 'B': red_substract_fact= atof(token[1]);break;
		case 'm': mcv_in = token[1];break;
		case 'M': mcv_in = token[1]; distance_thr = atof(token[2]);break;
		case 'R': recursive =true;
			//	case 'A': file_additionnal= true;
			
	}
	
}
int Taskscope<TASK_HIDDENMAP_BASE_DATAEXTRACTION>::defstore(char* const * token, int nbtoken){ 
	file_in[0] = token[0];
	file_in[1] = token[1];
	file_in[2] = token[2];
	file_in[3] = token[3];
	unsigned int i,j,k,l;
	char buffer[65536];
	
	// ./PMExtract_Features_From_Hidden -t $TMPDIR""$3""_feat_data.txt -S $TMPDIR""$3""_tsg.tif $TMPDIR""$3_tred.tif $TMPDIR""$3_tgre.tif $TMPDIR""$3_hid.tif
	
	
	
	//		char* fakeargs[] = {"./PMExtract_Features_From_Hidden", "-a", "0.9", "0.1", "10", "50", "-t", "../../../../report/$_feat_data7.txt", "-S", "../../../../report/$_tsg.tif", "../../../../report/$_tred.tif", "../../../../report/$_tgre.tif", "../../../../report/$_hid.tif"};
	//				char* fakeargs[] = {"./PMExtract_Features_From_Hidden", "-a", "0.9", "0.1", "10", "50", "-t", "../../report/$_feat_data6.txt", "-S", "../../report/$_tsg.tif", "../../report/$_tred.tif", "../../report/$_tgre.tif", "../../report/$_hid.tif"};
	
	
	//	mmm(sizeof(fakeargs) / sizeof(char*) ,fakesubstitution(sizeof(fakeargs) / sizeof(char*) , fakeargs, "HOwt_plate01_009023"));	
	
	//		char* fakeargs[] = {"./PMHiddenMapDirect", "-t", "../../../../autmptable.txt","../../../../jebdir/HOwt_plate01_002004_rw.tif", "../../../../jebdir/HOwt_plate01_002004_rw.tif", "../../../../autmptmp2.tif"};
	//		mmm(6,fakeargs);
	
	
	Tuple< TMatrix<double, 4,4> , 7> confalt; 
	
	Tuple< Tuple< double , 4> , 7> conf_mean; 
	
	TMatrix<double, 4,4> tmpmat, tmpmat2;
	Tuple<double,4> conf_key, conf_key2;
	GaussianDistribution<1> size_key;
	if (file_dascp){
		{//:
			SerialStore<unsigned int> dastore(file_dascp);
			Tuple< KeyElem<double, GaussianDistribution<4> > , 7 > keyframe;
			dastore.load(100,size_key);
			dastore.load(101,keyframe);
			for(i=0;i<7;i++) {confalt[i] = keyframe[i].d.ihvar.inverse(); conf_mean[i] = keyframe[i].d.mean;}
		}//:
	}
	
	ProgressBarPrint daprogbar; daprogbar.lenght = 20;
	
	
	TiffFile tf  = TiffFile(file_in[0]);
	TiffFile tfg = TiffFile(file_in[1]);
	TiffFile tfh = TiffFile(file_in[2]);
	
	
	TiffFile tf_sgim = TiffFile(seg_im_path);
	
	FILE* datable = table_out ? fopen(table_out,"w+") : NULL;
	
	if (table_out){
		fprintf(datable, "FrameID\tCellID\tGuess_Center_x\tGuess_Center_y\tGuess_Eccentric_x\tGuess_Eccentric_y\tGuess_Width\tGuess_Area\tCenter_x\tCenter_y\tEccentric_x\tEccentric_y\tWidth\tArea\tEdgeDistance_Radius\tEdgeDistance_Error\tDensity_In_Ellipse\tContour_pixels\tRamanujan\tRed_Intensity_LogFold_Error\tRed_Green_Correlation\tCell_type\tMassCenter_x\tMassCenter_y\tMassCenter_center_dist\tbudneck_x\tbudneck_y\tbudneck_center_dist\tRel_Cell_major\tRel_Cell_area\tRel_Cell_dist\tRel_EdgeDistance_Error\tRel_Density_In_Ellipse\tRel_contour_pixels\tRamanujan\tRel_Red_Intensity_LogFold_Error\tCell_prob_log_factor\tCell_error\tCell_prob\tRelation_prob_log_factor\tRelation_error\tRelation_prob\tRel_Cell_ID");
		
		strcpy(buffer,"Red"); fprintf(datable,"\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
		
		strcpy(buffer,"Green"); fprintf(datable,"\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
		strcpy(buffer,"G_SELF_DST"); fprintf(datable,"\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
		strcpy(buffer,"G_MASC_DST"); fprintf(datable,"\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
		strcpy(buffer,"G_EDGE_DST"); fprintf(datable,"\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
		strcpy(buffer,"G_CENT_DST"); fprintf(datable,"\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
		strcpy(buffer,"G_BUDN_DST"); fprintf(datable,"\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
		
		strcpy(buffer,"Intensity"); fprintf(datable,"\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
		strcpy(buffer,"N_SELF_DST"); fprintf(datable,"\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
		strcpy(buffer,"N_MASC_DST"); fprintf(datable,"\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
		strcpy(buffer,"N_EDGE_DST"); fprintf(datable,"\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
		strcpy(buffer,"N_CENT_DST"); fprintf(datable,"\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
		strcpy(buffer,"N_BUDN_DST"); fprintf(datable,"\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
		
		strcpy(buffer,"NR_SELF_DST"); fprintf(datable,"\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
		strcpy(buffer,"NR_MASC_DST"); fprintf(datable,"\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
		
		/*			if (file_additionnal){
		 strcpy(buffer,"G_ADDMASC_DST"); fprintf(datable,"\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
		 strcpy(buffer,"G_ADDOVER_INT"); fprintf(datable,"\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
		 strcpy(buffer,"N_ADDMASC_DST"); fprintf(datable,"\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
		 strcpy(buffer,"N_ADDOVER_INT"); fprintf(datable,"\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
		 }
		 */			
		
		fprintf(datable,"\n");
	}
	DataGrid<double, 2> im, img, ima;
	DataGrid<double, 2> seg_im;
	DataGrid<unsigned int, 2> imh;
	
	
	Tuple<double,5> conf_m_m[4];
	j=0;
	i=0;conf_m_m[j][i++]= -5.782724e+00;conf_m_m[j][i++]= -4.997317e+00;conf_m_m[j][i++]= 9.833092e-01;conf_m_m[j][i++]= 1.132820e-01;conf_m_m[j++][i++]= 1.866002e+03;
	i=0;conf_m_m[j][i++]= -5.084202e+00;conf_m_m[j][i++]= -3.968588e+00;conf_m_m[j][i++]= 9.855421e-01;conf_m_m[j][i++]= 1.055802e-01;conf_m_m[j++][i++]= 7.157430e+02;
	i=0;conf_m_m[j][i++]= -5.780408e+00;conf_m_m[j][i++]= -4.896796e+00;conf_m_m[j][i++]= 9.745668e-01;conf_m_m[j][i++]= 1.211551e-01;conf_m_m[j++][i++]= 1.308310e+03;
	i=0;conf_m_m[j][i++]= -5.575009e+00;conf_m_m[j][i++]= -4.775826e+00;conf_m_m[j][i++]= 9.845781e-01;conf_m_m[j][i++]= 9.943387e-02;conf_m_m[j++][i++]= 1.519127e+03;
	
	TMatrix<double,5,5>	conf_m_hiv[4];
	
	j=0;i=0;
	conf_m_hiv[j].data[i++]=-6.912987e-01;conf_m_hiv[j].data[i++]=4.328291e-01;conf_m_hiv[j].data[i++]=1.821810e-01;conf_m_hiv[j].data[i++]=-3.689752e-02;conf_m_hiv[j].data[i++]=-5.828682e-05;
	conf_m_hiv[j].data[i++]=4.328291e-01;conf_m_hiv[j].data[i++]=-1.045308e+00;conf_m_hiv[j].data[i++]=8.075323e+00;conf_m_hiv[j].data[i++]=2.178571e-01;conf_m_hiv[j].data[i++]=1.338076e-04;
	conf_m_hiv[j].data[i++]=1.821810e-01;conf_m_hiv[j].data[i++]=8.075323e+00;conf_m_hiv[j].data[i++]=-1.638523e+02;conf_m_hiv[j].data[i++]=-2.050080e+00;conf_m_hiv[j].data[i++]=4.550996e-03;
	conf_m_hiv[j].data[i++]=-3.689752e-02;conf_m_hiv[j].data[i++]=2.178571e-01;conf_m_hiv[j].data[i++]=-2.050080e+00;conf_m_hiv[j].data[i++]=-5.233483e+00;conf_m_hiv[j].data[i++]=-2.758743e-05;
	conf_m_hiv[j].data[i++]=-5.828682e-05;conf_m_hiv[j].data[i++]=1.338076e-04;conf_m_hiv[j].data[i++]=4.550996e-03;conf_m_hiv[j].data[i++]=-2.758743e-05;conf_m_hiv[j].data[i++]=-2.424722e-06;
	j++;i=0;
	conf_m_hiv[j].data[i++]=-1.258348e+00;conf_m_hiv[j].data[i++]=1.301703e+00;conf_m_hiv[j].data[i++]=-9.984736e-01;conf_m_hiv[j].data[i++]=-4.032630e-03;conf_m_hiv[j].data[i++]=2.665070e-05;
	conf_m_hiv[j].data[i++]=1.301703e+00;conf_m_hiv[j].data[i++]=-1.810458e+00;conf_m_hiv[j].data[i++]=5.079954e+00;conf_m_hiv[j].data[i++]=2.051689e-02;conf_m_hiv[j].data[i++]=-1.355913e-04;
	conf_m_hiv[j].data[i++]=-9.984736e-01;conf_m_hiv[j].data[i++]=5.079954e+00;conf_m_hiv[j].data[i++]=-1.135338e+02;conf_m_hiv[j].data[i++]=-4.585396e-01;conf_m_hiv[j].data[i++]=3.030380e-03;
	conf_m_hiv[j].data[i++]=-4.032630e-03;conf_m_hiv[j].data[i++]=2.051689e-02;conf_m_hiv[j].data[i++]=-4.585396e-01;conf_m_hiv[j].data[i++]=-5.248240e+00;conf_m_hiv[j].data[i++]=3.734151e-04;
	conf_m_hiv[j].data[i++]=2.665070e-05;conf_m_hiv[j].data[i++]=-1.355913e-04;conf_m_hiv[j].data[i++]=3.030380e-03;conf_m_hiv[j].data[i++]=3.734151e-04;conf_m_hiv[j].data[i++]=-2.467812e-06;
	j++;i=0;
	conf_m_hiv[j].data[i++]=-7.332171e-01;conf_m_hiv[j].data[i++]=5.301192e-01;conf_m_hiv[j].data[i++]=-4.762023e-01;conf_m_hiv[j].data[i++]=-7.479692e-02;conf_m_hiv[j].data[i++]=1.075601e-04;
	conf_m_hiv[j].data[i++]=5.301192e-01;conf_m_hiv[j].data[i++]=-1.042033e+00;conf_m_hiv[j].data[i++]=8.768116e+00;conf_m_hiv[j].data[i++]=1.729375e-01;conf_m_hiv[j].data[i++]=-1.430582e-04;
	conf_m_hiv[j].data[i++]=-4.762023e-01;conf_m_hiv[j].data[i++]=8.768116e+00;conf_m_hiv[j].data[i++]=-2.242465e+02;conf_m_hiv[j].data[i++]=-1.568487e+00;conf_m_hiv[j].data[i++]=7.069388e-03;
	conf_m_hiv[j].data[i++]=-7.479692e-02;conf_m_hiv[j].data[i++]=1.729375e-01;conf_m_hiv[j].data[i++]=-1.568487e+00;conf_m_hiv[j].data[i++]=-4.375910e+00;conf_m_hiv[j].data[i++]=2.275300e-05;
	conf_m_hiv[j].data[i++]=1.075601e-04;conf_m_hiv[j].data[i++]=-1.430582e-04;conf_m_hiv[j].data[i++]=7.069388e-03;conf_m_hiv[j].data[i++]=2.275300e-05;conf_m_hiv[j].data[i++]=-3.301774e-06;
	j++;i=0;
	conf_m_hiv[j].data[i++]=-5.925807e-01;conf_m_hiv[j].data[i++]=4.078784e-01;conf_m_hiv[j].data[i++]=-4.036895e-01;conf_m_hiv[j].data[i++]=-9.704348e-02;conf_m_hiv[j].data[i++]=1.288953e-04;
	conf_m_hiv[j].data[i++]=4.078784e-01;conf_m_hiv[j].data[i++]=-8.194436e-01;conf_m_hiv[j].data[i++]=6.570772e+00;conf_m_hiv[j].data[i++]=1.949642e-01;conf_m_hiv[j].data[i++]=-6.533255e-05;
	conf_m_hiv[j].data[i++]=-4.036895e-01;conf_m_hiv[j].data[i++]=6.570772e+00;conf_m_hiv[j].data[i++]=-1.478156e+02;conf_m_hiv[j].data[i++]=-1.563335e+00;conf_m_hiv[j].data[i++]=2.502390e-03;
	conf_m_hiv[j].data[i++]=-9.704348e-02;conf_m_hiv[j].data[i++]=1.949642e-01;conf_m_hiv[j].data[i++]=-1.563335e+00;conf_m_hiv[j].data[i++]=-4.268293e+00;conf_m_hiv[j].data[i++]=1.554409e-05;
	conf_m_hiv[j].data[i++]=1.288953e-04;conf_m_hiv[j].data[i++]=-6.533255e-05;conf_m_hiv[j].data[i++]=2.502390e-03;conf_m_hiv[j].data[i++]=1.554409e-05;conf_m_hiv[j].data[i++]=-1.013794e-06;
	Tuple<double , 5> conf_diff, conf_xfor;
	
	FILE* f_mcv = mcv_in ? fopen(mcv_in,"rb+") : NULL;
	
	double wei = 1.0f; // should not me changed unless seg image is used
	DataGrid<unsigned int, 2>::KeyIterator ite(imh);
	unsigned int dacell;
	
	Madstructs::CellRecord_Header_Reborn_AGAIN header;
	Madstructs::MultiCover damcv;
	unsigned int mcv_i,mcv_j, mcv_index;
	if (f_mcv == NULL) {header.g_center_x = 0.0f;header.g_center_y = 0.0f;header.g_excentric_x = 0.0f;header.g_excentric_y = 0.0f;header.g_width = 0.0f;header.g_area = 0.0f;}
	else  {
		if (distance_thr == 0.0f){
			mcv_index=1;
			do {
				damcv.load(f_mcv);mcv_i=0;mcv_j=0;
				while (mcv_i < damcv.group.size()){
					if (damcv.group[mcv_i].cover.size() == 0) mcv_i++;
					else break;
				}
			}while (mcv_i >= damcv.group.size());
		}
	}
	
	double core_data[100];
	
	Tuple<unsigned int , 2 > coor;
	
	Vector< Madstructs::CellPose > da_extern_cells;
	HeapTree< KeyElem<double, Tuple<unsigned int, 2> > > da_extern_order;
	Vector< unsigned int > da_extern_backmap; 
	Vector< unsigned int > da_extern_backmap_back; 		
	
	Vector< Madstructs::TmpScope > reccells;
	
	Tuple<Tuple<unsigned int, 2>,(TEMPLATE_INT_POWER<3,2>::ans) -1 > ncoor; int nbncoor;
	//Tuple< Tuple<unsigned int , 2>,  2 * (2 -1) * 4   > ncoork;int nbncoork;
	unsigned int fr_id =0;
	
	map<unsigned int, unsigned int> maping;
	map<unsigned int, unsigned int>::iterator ite_map;
	double pix[32];
	Vector< KeyElem<unsigned int, unsigned int> > backmap; 
	Vector< KeyElem<unsigned int, Tuple<unsigned int, 2> > > pixlist;
	while(tf.fetch(im)){
		if (!tfg.fetch(img)) {printf("Incohenrent number of frame in '%s' and '%s' input images (critical)\n",file_in[0], file_in[1]);exit(1);}
		if (!tfh.fetch(imh)) {printf("Incohenrent number of frame in '%s' and '%s' input images (critical)\n",file_in[0], file_in[2]);exit(1);}
		
		
		if (seg_im_path)	{if (!tf_sgim.fetch(seg_im)) {printf("Incohenrent number of frame in '%s' and '%s' input images (critical)\n",file_in[0], seg_im_path);exit(1);}}
		DataGrid<double, 2> im_odist = imh.ExpectedDistanceToOther(0.01f);
		
		//				tfout.put( im_odist, (float) 0.0f, (float)100.0f);
		maping[0] = 0;
		reccells.push_back(Madstructs::TmpScope());// makes the background be 0 beforehand
		
		if (ite.first()) do{
			coor = ite();
			dacell = imh(coor);
			
			if (red_substract_fact != 0.0f){
				img(coor) -= im(coor) * red_substract_fact;
				if (img(coor) < 0.0f) img(coor) = 0.0f;
			}
			
			
			if ( (ite_map = maping.find(dacell)) == maping.end()){
				maping[dacell] = reccells.size();
				backmap.push_back(KeyElem<unsigned int, unsigned int>(dacell,reccells.size()));
				reccells.push_back(Madstructs::TmpScope());
				reccells[dacell].backID = dacell;
				dacell = reccells.size()-1;
				
				
			}else dacell = ite_map->second;
			
			nbncoor = imh.get_indirectNeightbor(coor, ncoor);
			//	nbncoork = imh.get_indirectNeightbor(coor, ncoork);
			
			for(j=0;j<nbncoor;j++){
				
				k = imh(ncoor[j]);
				if ( (ite_map = maping.find(k)) == maping.end()){
					maping[k] = reccells.size();
					backmap.push_back(KeyElem<unsigned int, unsigned int>(k,reccells.size()));
					reccells.push_back(Madstructs::TmpScope());
					reccells[dacell].backID = k;
					k = reccells.size()-1;
				}else k= ite_map->second;
				
				if ((k != dacell)&&(k !=0)){
					for(l=0;l<reccells[dacell].neighbor_size;l++) if (reccells[dacell].neighbor[l] == k) break;
					if ((l >= reccells[dacell].neighbor_size)&&(reccells[dacell].neighbor_size < 32)) reccells[dacell].neighbor[reccells[dacell].neighbor_size++] = k;
				}
			}
			
			
			
			
			if (seg_im_path)  wei = seg_im(coor);
			reccells[dacell].frameID = fr_id;
			reccells[dacell].ellipic[0] += wei;
			reccells[dacell].ellipic[1] += wei*coor[0];
			reccells[dacell].ellipic[2] += wei*coor[1];
			reccells[dacell].ellipic[3] += wei*coor[0]*coor[0];
			reccells[dacell].ellipic[4] += wei*coor[1]*coor[1];
			reccells[dacell].ellipic[5] += wei*coor[0]*coor[1];
			
			reccells[dacell].centerofmass[0] += wei*coor[0] * img(coor);
			reccells[dacell].centerofmass[1] += wei*coor[1] * img(coor);
			reccells[dacell].centerofmass[2] += wei*img(coor);
			reccells[dacell].centerofmass[3] += wei*coor[0] * im(coor);
			reccells[dacell].centerofmass[4] += wei*coor[1] * im(coor);
			
			
			
			reccells[dacell].red += WeightElem<double,4>(im(coor),wei);
			reccells[dacell].green += WeightElem<double,4>(img(coor),wei);
			reccells[dacell].red_green += WeightElem<double,1>(im(coor)*img(coor),wei);
			
			reccells[dacell].r_ed_dist += WeightElem<double,4>(im_odist(coor) ,im(coor)*wei);
			reccells[dacell].ed_dist += WeightElem<double,4>(im_odist(coor) ,img(coor)*wei);
			
			if (dacell != 0) pixlist.push_back(  KeyElem<unsigned int, Tuple<unsigned int, 2> >(dacell,coor) );
		} while(ite.next());
		
		// find relations!
		
		
		
		
		
		//	printf("%i itmes\n",reccells.size());
		
		
		// compute radius error, and nb of countour pixels
		
		pixlist.sort();
		daprogbar.start("Computing the Confidence Measures");
		
		j =0; // min index for pts;
		for(i=1;i<reccells.size();i++){
			daprogbar.update(((double)i) / reccells.size());
			for(;pixlist[j].k != i;j++);
			pix[0] = reccells[i].ellipic[1] / reccells[i].ellipic[0];
			pix[1] = reccells[i].ellipic[2] / reccells[i].ellipic[0];
			ExOp::toZero(reccells[i].distance_error);
			reccells[i].contour =0;
			for(k=j;pixlist[k].k == i;k++){
				if (seg_im_path)  wei = seg_im(pixlist[k].d);
				reccells[i].distance_error += WeightElem<double,2>(hypot(pix[0] -pixlist[k].d[0], pix[1] -pixlist[k].d[1]) + im_odist(pixlist[k].d),wei);
				nbncoor = imh.get_indirectNeightbor(pixlist[k].d, ncoor);
				int nb_ext=0;
				for(l=0;l<nbncoor ;l++){
					if (maping[imh(ncoor[l])] != i) nb_ext++;
				}
				if (nb_ext + 6  >= nbncoor) reccells[i].contour++; // if has 2 neighbor pixel, its a countour pixel!
			}
			
		}daprogbar.finish();
		
		
		daprogbar.start("Cell Confidence and Labeling Artifacts");
		for(i=1;i<reccells.size();i++){
			Madstructs::CellPose tmpcellpose;
			tmpcellpose.initFromCummul(&(reccells[i].ellipic[0]));
			// reccells[i].ellipic[0] = tmpcellpose.area;  guarrantied to be equal
			reccells[i].ellipic[1] = tmpcellpose.center[0];
			reccells[i].ellipic[2] = tmpcellpose.center[1];
			reccells[i].ellipic[3] = tmpcellpose.eccentric[0];
			reccells[i].ellipic[4] = tmpcellpose.eccentric[1];
			reccells[i].ellipic[5] = tmpcellpose.width;
			
			// temporaty minor axis
			reccells[i].ramanujan = sqrt(tmpcellpose.width*tmpcellpose.width - 4.0f * (tmpcellpose.eccentric[0]  * tmpcellpose.eccentric[0] + tmpcellpose.eccentric[1] * tmpcellpose.eccentric[1]));// minor axis
			reccells[i].density_error = tmpcellpose.area / (M_PI * 0.25f* tmpcellpose.width * reccells[i].ramanujan );
			reccells[i].ramanujan = reccells[i].contour / (M_PI * (1.5f * (tmpcellpose.width + reccells[i].ramanujan) - sqrt(2.5f * tmpcellpose.width * reccells[i].ramanujan+ 0.75f*(tmpcellpose.width*tmpcellpose.width + reccells[i].ramanujan*reccells[i].ramanujan))));
			
			
			pix[0] = reccells[i].distance_error.getMean();
			pix[1] = reccells[i].distance_error.getVar();
			
			reccells[i].dist_error =  pix[1] / (pix[0]*pix[0]);
			
			if (file_dascp){
				pix[0] = (reccells[i].ellipic[0] - size_key.mean[0]);
				pix[0] = pix[0] * pix[0] * size_key.ihvar.data[0];
				
				pix[1] = reccells[i].ellipic[0] / 500.0f;
				j = (unsigned int)floor(pix[1]);
				if (j >5) {j=5; pix[1]=1.0f;}
				else pix[1] -=j;
				
				tmpmat = (confalt[j] * (1.0f - pix[1])) + (confalt[j+1] * pix[1]); 
				
				conf_key2 = (conf_mean[j] * (1.0f - pix[1])) + (conf_mean[j+1] * pix[1]);
				conf_key[0] = log(reccells[i].dist_error);
				conf_key[1] = log(1.0f - reccells[i].density_error);
				conf_key[2] = reccells[i].ramanujan;
				conf_key[3] = (ignore_red_in_confidence) ? conf_key2[3]: log( reccells[i].red.getMean() );
				conf_key -= conf_key2;
				conf_key *= conf_key; 
				//		printf("%e\t%e\t%e\t%e\t%e\n", (conf_key[0]  / tmpmat.data[0]), (conf_key[1]  / tmpmat.data[5]) , (conf_key[2] /tmpmat.data[10]), (conf_key[3] /tmpmat.data[15]), pix[0]);
				//		printf("%e\t%e\t%e\t%e\t%e\n", tmpmat.data[0], tmpmat.data[5] ,tmpmat.data[10], tmpmat.data[15], size_key.ihvar.data[0]);
				if (tmpmat.data[0] > 0.0f) {printf("invalid confidence prior parameter!\n");exit(1);}
				if (tmpmat.data[5] > 0.0f) {printf("invalid confidence prior parameter!\n");exit(1);}
				if (tmpmat.data[10] > 0.0f) {printf("invalid confidence prior parameter!\n");exit(1);}
				if (tmpmat.data[15] > 0.0f) {printf("invalid confidence prior parameter!\n");exit(1);}
				if (size_key.ihvar.data[0] > 0.0f) {printf("invalid confidence prior parameter!\n");exit(1);}
				reccells[i].cell_prob_factor = -2.5f * log(M_PI) + 0.5f * (log(fabs(size_key.ihvar.data[0])) - log(fabs(tmpmat.data[0])) - log(fabs(tmpmat.data[5])) - log(fabs(tmpmat.data[10])) - log(fabs(tmpmat.data[15]))   ); 
				reccells[i].cell_prob = (conf_key[0] / tmpmat.data[0])  +(conf_key[1] / tmpmat.data[5]) + (conf_key[2] /tmpmat.data[10]) + (conf_key[3] /tmpmat.data[15]) + pix[0];
				
			}else{
				/*						
				 switch(reccells[i].type){
				 case 'm':j=0;break;
				 case 'b':j=1;break;
				 case 'd':j=2;break;
				 default: j=3;
				 }
				 conf_diff[0] = log(reccells[i].dist_error) - conf_m_m[j][0];
				 conf_diff[1] = log(1.0f - reccells[i].density_error) - conf_m_m[j][1];
				 conf_diff[2] = reccells[i].ramanujan - conf_m_m[j][2];
				 conf_diff[3] = reccells[i].fold_red_err - conf_m_m[j][3];
				 conf_diff[4] = reccells[i].ellipic[0] - conf_m_m[j][4];
				 
				 conf_xfor = conf_m_hiv[j] * conf_diff;
				 reccells[i].cell_prob_factor = 1.0f;
				 reccells[i].cell_prob = exp(0.5f * conf_diff.dotProduct(conf_xfor));*/
			}
			
			
			header.density_error = tmpcellpose.area / (M_PI * 0.5f* tmpcellpose.width * sqrt(0.25f* reccells[i].ellipic[5]*reccells[i].ellipic[5] - reccells[i].ellipic[3]  * reccells[i].ellipic[3] - reccells[i].ellipic[4] * reccells[i].ellipic[4]));
			if (1.0f / (1.0f + art_k *exp (-reccells[i].cell_prob - reccells[i].cell_prob_factor)) < artifact_class[0]) reccells[i].type = 'a';
			else if (header.density_error < artifact_class[2]) reccells[i].type = 'a';
			else if (reccells[i].ellipic[0] <= artifact_class[1]) reccells[i].type = 'a';
			else {
				pix[0] = reccells[i].distance_error.getMean();
				pix[1] = reccells[i].distance_error.getVar();
				if ((pix[1] / (pix[0]*pix[0])) > artifact_class[3] ) reccells[i].type = 'a';
				else {
					// ultimate density test!
					if (1.0f - header.density_error < (artifact_class[4] / tmpcellpose.area)) reccells[i].type = (reccells[i].neighbor_size == 0) ? 'l' : 'c';
					else reccells[i].type = 'a';
				}
			}
			daprogbar.update(((double)i) / reccells.size());
		}daprogbar.finish();
		
		bool dirty;
		unsigned int *ll = new unsigned int[(reccells.size() << 1)];
		do{
			dirty = false;
			daprogbar.start("Comparing Neighboring object sizes");
			
			static_warning_handdle(ll == NULL, LFH_ERROR_CANT_ALLOCATE_MEMORY);
			for(i=1;i<reccells.size();i++) if (reccells[i].type == 'c') {ll[(i<<1)] =0;ll[(i<<1)|1] =0;} 
			for(i=1;i<reccells.size();i++){
				if (reccells[i].type != 'c') continue;
				for(j=0;j< reccells[i].neighbor_size;j++){
					l = reccells[i].neighbor[j];
					if (reccells[l].type != 'c') continue;
					if (reccells[i].ellipic[0] < reccells[l].ellipic[0]){
						// pot mother
						if ((ll[(i<<1)|1] == 0) || ( reccells[ll[(i<<1)|1]].ellipic[0] < reccells[l].ellipic[0])) ll[(i<<1)|1] = l;
						if ((ll[(l<<1)  ] == 0) || ( reccells[ll[(l<<1)  ]].ellipic[0] > reccells[i].ellipic[0])) ll[(l<<1)  ] = i;
					}else{
						// pot daughter
						if ((ll[(i<<1)  ] == 0) || ( reccells[ll[(i<<1)  ]].ellipic[0] > reccells[l].ellipic[0])) ll[(i<<1)  ] = l;
						if ((ll[(l<<1)|1] == 0) || ( reccells[ll[(l<<1)|1]].ellipic[0] < reccells[i].ellipic[0])) ll[(l<<1)|1] = i;
					}
				}
				daprogbar.update(((double)i) / reccells.size());
			}daprogbar.finish();
			
			daprogbar.start("Find Mother-Bud Pairs");
			for(i=1;i<reccells.size();i++){
				if (reccells[i].type != 'c') continue;
				if ((l = ll[(i<<1)]) != 0){ 
					if ((ll[(l<<1)|1] == i)&&(ll[l<<1] == 0)){
						reccells[i].type = 'm';
						reccells[i].otherID = l;
						reccells[l].type = 'b';
						reccells[l].otherID = i;
						dirty = true;
					}
				}
				daprogbar.update(((double)i) / reccells.size());
			}daprogbar.finish();
		}while((dirty)&&(recursive));
		
		
		daprogbar.start("Find Daughters Cells");
		for(i=1;i<reccells.size();i++){
			reccells[i].centerofmass[0] = reccells[i].centerofmass[0] / reccells[i].centerofmass[2];
			reccells[i].centerofmass[1] = reccells[i].centerofmass[1] / reccells[i].centerofmass[2];
			reccells[i].centerofmass[3] = reccells[i].centerofmass[3] / reccells[i].centerofmass[5];
			reccells[i].centerofmass[4] = reccells[i].centerofmass[4] / reccells[i].centerofmass[5];
			if ((reccells[i].type == 'l')||(reccells[i].type == 'a')) {reccells[i].otherID = 0xFFFFFFFF; continue;}
			if (reccells[i].type != 'c') continue;
			reccells[i].otherID = 0xFFFFFFFF;
			for(j=0;j< reccells[i].neighbor_size;j++) {
				if ((reccells[reccells[i].neighbor[j]].type == 'b')||(reccells[reccells[i].neighbor[j]].type == 'm')){
					for(k=j+1; k < reccells[i].neighbor_size;k++) if (reccells[i].neighbor[k] == reccells[reccells[i].neighbor[j]].otherID) break;
					if (k < reccells[i].neighbor_size){
						if (reccells[i].otherID == 0xFFFFFFFF) {reccells[i].otherID = (reccells[reccells[i].neighbor[j]].type == 'b') ? reccells[i].neighbor[j] : reccells[reccells[i].neighbor[j]].otherID; reccells[i].type = 'd';}
					}else{
						for(k=0; k < reccells[reccells[reccells[i].neighbor[j]].otherID].neighbor_size;k++) if (reccells[reccells[reccells[i].neighbor[j]].otherID].neighbor[k] == i) break;
						if (k < reccells[reccells[reccells[i].neighbor[j]].otherID].neighbor_size){
							if (reccells[i].otherID == 0xFFFFFFFF) {reccells[i].otherID = (reccells[reccells[i].neighbor[j]].type == 'b') ? reccells[i].neighbor[j] : reccells[reccells[i].neighbor[j]].otherID; reccells[i].type = 'd';}
						}
					}
				}
			}
			daprogbar.update(((double)i) / reccells.size());
		}daprogbar.finish(); 
		delete[](ll);
		
		daprogbar.start("Computing Bud Neck Coordinates");
		j =0; // min index for pts;
		for(i=1;i<reccells.size();i++){
			
			for(;j<pixlist.size();j++) if (pixlist[j].k == i) break;
			
			if (reccells[i].type == 'b'){
				// find the set of point between cur objects and 
				memset(pix,'\0',sizeof(double)*32);
				for(;((j<pixlist.size()) && (pixlist[j].k == i));j++){
					nbncoor = imh.get_indirectNeightbor(pixlist[j].d, ncoor);
					for(l=0;l<nbncoor ;l++){
						if (maping[imh(ncoor[l])] == reccells[i].otherID){
							pix[0] += pixlist[j].d[0] + ncoor[l][0];
							pix[1] += pixlist[j].d[1] + ncoor[l][1];
							pix[2] += 2.0f;
						}
					}
					reccells[i].neck[0] = pix[0] / pix[2];
					reccells[i].neck[1] = pix[1] / pix[2];
					reccells[reccells[i].otherID].neck[0] = pix[0] / pix[2];
					reccells[reccells[i].otherID].neck[1] = pix[1] / pix[2];
				}
			}else if (reccells[i].type == 'd'){
				// find the set of point between cur objects and 
				memset(pix,'\0',sizeof(double)*32);
				for(;((j<pixlist.size()) && (pixlist[j].k == i));j++){
					nbncoor = imh.get_indirectNeightbor(pixlist[j].d, ncoor);
					for(l=0;l<nbncoor ;l++){
						if ((maping[imh(ncoor[l])] == reccells[i].otherID)||(maping[imh(ncoor[l])] == reccells[reccells[i].otherID].otherID)){
							pix[0] += pixlist[j].d[0] + ncoor[l][0];
							pix[1] += pixlist[j].d[1] + ncoor[l][1];
							pix[2] += 2.0f;
						}
					}
					reccells[i].neck[0] = pix[0] / pix[2];
					reccells[i].neck[1] = pix[1] / pix[2];					
				}					
			}
			daprogbar.update(((double)i) / reccells.size());
		}daprogbar.finish();
		
		
		//		printf("at point R!\n");fflush(stdout);
		
		daprogbar.start("Computing Morphological Distances");
		j =0; // min index for pts;
		for(j=0,i=1;i<reccells.size();i++){
			
			// espected red (0,15) (1000,27) (1650,29.5) (4000,31.5)
			if (reccells[i].ellipic[0] > 1650){
				pix[0] = 29.5f + (reccells[i].ellipic[0] - 1650) / 1175.0f; 
			}else if (reccells[i].ellipic[0] > 1000){
				pix[0] = 27.0f + (reccells[i].ellipic[0] - 1000) / 260.0f; 
			}else{
				pix[0] = 15.0f + 0.012f* reccells[i].ellipic[0]; 
			}
			core_data[1] = pix[0] / reccells[i].red.getMean();
			reccells[i].fold_red_err = log(reccells[i].red.getMean()) - log(pix[0]) ;
			
			
			for(;j<pixlist.size();j++) if (pixlist[j].k == i) break;
			if ((reccells[i].type == 'c')||(reccells[i].type == 'l')||(reccells[i].type == 'a')){
				reccells[i].bn_dist += WeightElem<double,4>(0 ,0);
				reccells[i].bn_dist += WeightElem<double,4>(0 ,0);
				reccells[i].bn_dist += WeightElem<double,4>(0 ,0);
				reccells[i].r_bn_dist += WeightElem<double,4>(0 ,0);
				reccells[i].r_bn_dist += WeightElem<double,4>(0 ,0);
				reccells[i].r_bn_dist += WeightElem<double,4>(0 ,0);
			}
			l=0;
			for(k=j;k<pixlist.size();k++){
				if (pixlist[k].k > i) break;
				
				if (seg_im_path)  wei = seg_im(pixlist[k].d);
				core_data[2] = img(pixlist[k].d);
				core_data[5] = im(pixlist[k].d);
				reccells[i].red_intensity += WeightElem<double,1>(core_data[2] * core_data[1] * im(pixlist[k].d)  ,wei);
				reccells[i].intensity += WeightElem<double,4>(core_data[2] * core_data[1] ,wei);
				core_data[2] *= wei;
				core_data[5] *= wei;
				if (!((reccells[i].type == 'c')||(reccells[i].type == 'l')||(reccells[i].type == 'a'))) {
					reccells[i].bn_dist += WeightElem<double,4>(hypot(reccells[i].neck[0] - pixlist[k].d[0],reccells[i].neck[1] - pixlist[k].d[1]) ,core_data[2]);
					reccells[i].r_bn_dist += WeightElem<double,4>(hypot(reccells[i].neck[0] - pixlist[k].d[0],reccells[i].neck[1] - pixlist[k].d[1]) ,core_data[5]);
				}
				pix[0] = reccells[i].ellipic[1];
				pix[1] = reccells[i].ellipic[2];
				
				reccells[i].cn_dist += WeightElem<double,4>(hypot(pix[0] - (double)(pixlist[k].d[0]),pix[1] - (double)(pixlist[k].d[1])) ,core_data[2]);
				reccells[i].r_cn_dist += WeightElem<double,4>(hypot(pix[0] - (double)(pixlist[k].d[0]),pix[1] - (double)(pixlist[k].d[1])) ,core_data[5]);
				
				reccells[i].cm_dist += WeightElem<double,4>(hypot(reccells[i].centerofmass[0] - (double)(pixlist[k].d[0]),reccells[i].centerofmass[1] - (double)(pixlist[k].d[1])) ,core_data[2]);
				reccells[i].r_gcm_dist += WeightElem<double,4>(hypot(reccells[i].centerofmass[0] - (double)(pixlist[k].d[0]),reccells[i].centerofmass[1] - (double)(pixlist[k].d[1])) ,core_data[5]);
				reccells[i].r_cm_dist += WeightElem<double,4>(hypot(reccells[i].centerofmass[3] - (double)(pixlist[k].d[0]),reccells[i].centerofmass[4] - (double)(pixlist[k].d[1])) ,core_data[5]);
				reccells[i].sd_dist += WeightElem<double,4>(0.0f ,core_data[2] *core_data[2] *0.5f);
				reccells[i].r_sd_dist += WeightElem<double,4>(0.0f ,core_data[5] *core_data[5] *0.5f);
				if (l != 50000){
					for(l=k+1;l< pixlist.size();l++){
						if ((pixlist[l].k > i)||(l > k +10000))  break;							
						if (seg_im_path)  wei = seg_im(pixlist[l].d);
						reccells[i].sd_dist += WeightElem<double,4>(hypot((double)(pixlist[l].d[0]) - (double)(pixlist[k].d[0]),(double)(pixlist[l].d[1]) - (double)(pixlist[k].d[1])) ,img(pixlist[l].d) * wei * core_data[2]);
						reccells[i].r_sd_dist += WeightElem<double,4>(hypot((double)(pixlist[l].d[0]) - (double)(pixlist[k].d[0]),(double)(pixlist[l].d[1]) - (double)(pixlist[k].d[1])) ,im(pixlist[l].d) * wei * core_data[5]);
						reccells[i].rg_sd_dist += WeightElem<double,4>(hypot((double)(pixlist[l].d[0]) - (double)(pixlist[k].d[0]),(double)(pixlist[l].d[1]) - (double)(pixlist[k].d[1])) ,img(pixlist[l].d) * wei * core_data[5]);
					}
					if (l > k +10000) {l = 50000;
						ExOp::toZero(reccells[i].sd_dist); reccells[i].sd_dist+= WeightElem<double,4>(0 ,0);reccells[i].sd_dist+= WeightElem<double,4>(0 ,0);
						ExOp::toZero(reccells[i].r_sd_dist); reccells[i].r_sd_dist+= WeightElem<double,4>(0 ,0);reccells[i].r_sd_dist+= WeightElem<double,4>(0 ,0);
						ExOp::toZero(reccells[i].rg_sd_dist); reccells[i].rg_sd_dist+= WeightElem<double,4>(0 ,0);reccells[i].rg_sd_dist+= WeightElem<double,4>(0 ,0);
			}	}	}
			daprogbar.update(((double)i) / reccells.size());
		}daprogbar.finish();
		
		if ((f_mcv != NULL)&&(distance_thr != 0.0f)) { // load independent MCV, and maps cells
			da_extern_cells.clear();
			da_extern_order.clear();
			damcv.load(f_mcv);
			daprogbar.start("Mapping the Hidden Map to External MCV file");
			for(i=0;i<damcv.group.size();i++)
				for(j=0;j<damcv.group[i].cover.size();j++){
					da_extern_cells.push_back(damcv.group[i].cover[j]);
					da_extern_cells[da_extern_cells.size()-1].center[0] += damcv.group[i].rect[0];
					da_extern_cells[da_extern_cells.size()-1].center[1] += damcv.group[i].rect[1];
				}
			
			KeyElem< double, Tuple<unsigned int,2 > > inin;
			for(inin.d[1] = 0;inin.d[1]< reccells.size();inin.d[1]++){
				daprogbar.update(0.75f * (((double)inin.d[1]) / reccells.size()));
				pix[0] = reccells[inin.d[1]].ellipic[1];
				pix[1] = reccells[inin.d[1]].ellipic[2];
				for(inin.d[0] = 0;inin.d[0]< da_extern_cells.size();inin.d[0]++){
					inin.k = hypot(da_extern_cells[inin.d[0]].center[0] - pix[0],da_extern_cells[inin.d[0]].center[1] - pix[1]);
					da_extern_order.insert(inin);
				}
			}
			da_extern_backmap.setSize(reccells.size()); ExOp::toMax(da_extern_backmap);
			da_extern_backmap_back.setSize(da_extern_cells.size()); ExOp::toOne(da_extern_backmap_back);
			i = (reccells.size() > da_extern_cells.size()) ? da_extern_cells.size() : reccells.size();
			j=0;
			while(!da_extern_order.isempty()){
				inin = da_extern_order.pop();
				if (inin.k > distance_thr) break;
				if (da_extern_backmap[inin.d[1]] != 0xFFFFFFFF) continue;
				if (da_extern_backmap_back[inin.d[0]] != 1) continue;
				j++;
				daprogbar.update(0.25f * (3.0f+ (((double)j) / i)));
				da_extern_backmap_back[inin.d[0]] = 0;
				da_extern_backmap[inin.d[1]] = inin.d[0];
				if (j ==i) break;
			}
			daprogbar.finish();
			printf("mapped %i circles (tothidden=%i,totmcv=%i)\n", j, reccells.size(), da_extern_cells.size());
		}
		
		
		//		printf("backmap size = %i, cellrecsize = %i\n", backmap.size(),reccells.size());
		backmap.sort();
		daprogbar.start("Writing Cell in output table");	
		for(i=0;i<backmap.size();i++){
			dacell = backmap[i].d;
			if (dacell >= reccells.size()){
				printf("got %i here!\n", dacell); exit(1);
			}
			
			header.area = reccells[dacell].ellipic[0];
			header.center_x = reccells[dacell].ellipic[1];
			header.center_y = reccells[dacell].ellipic[2];
			header.excentric_x = reccells[dacell].ellipic[3];
			header.excentric_y = reccells[dacell].ellipic[4];
			header.width = reccells[dacell].ellipic[5];
			
			header.contour = reccells[dacell].contour;
			
			header.ramanujan = reccells[dacell].ramanujan; 
			header.density_error = reccells[dacell].density_error;
			
			pix[0] = reccells[dacell].distance_error.getMean();
			
			header.edgedist_radius = pix[0];
			header.edgedist_raw_error = reccells[dacell].dist_error;
			header.celltype = reccells[dacell].type;
			header.cellID = backmap[i].k;
			header.fold_red_err = reccells[dacell].fold_red_err;
			
			if ((reccells[dacell].type == 'c')||(reccells[dacell].type == 'l')||(reccells[dacell].type == 'a')) {
				header.other_dist = 0.0f; 
				header.other_major= 0.0f; 
				header.other_area = 0.0f; 
				header.otherID = -1;
				header.budneck_x = 0.0f;
				header.budneck_y = 0.0f;
				header.other_dist_err =0.0f; 
				header.other_density=0.0f; 
				header.other_red_dev=0.0f; 
				header.other_contour=0.0f; 
				header.other_ramanujan=0.0f; 
				header.relation_prob =0.0f;
				header.relation_prob_factor =0.0f;
				header.relation_error = 0.0;
			}else {
				pix[0] = reccells[reccells[dacell].otherID].distance_error.getMean(); 
				
				header.other_contour=reccells[reccells[dacell].otherID].contour;
				header.other_dist = pix[0];
				header.other_dist_err = reccells[reccells[dacell].otherID].dist_error;
				
				header.other_major= reccells[reccells[dacell].otherID].ellipic[5];
				header.other_area = reccells[reccells[dacell].otherID].ellipic[0]; 
				
				
				header.other_ramanujan = reccells[reccells[dacell].otherID].ramanujan;
				header.other_density = reccells[reccells[dacell].otherID].density_error;
				// temporaty minor axis
				
				header.other_red_dev = reccells[reccells[dacell].otherID].fold_red_err;
				
				for(j=0;j<backmap.size();j++) if (backmap[j].d == reccells[dacell].otherID) break;
				header.otherID = (j < backmap.size()) ? backmap[j].k : 0xFFFFFFFF;
				header.budneck_x = reccells[dacell].neck[0];
				header.budneck_y = reccells[dacell].neck[1];
				
				header.relation_error = reccells[reccells[dacell].otherID].cell_prob;
				header.relation_prob_factor = reccells[reccells[dacell].otherID].cell_prob_factor;
				header.relation_prob = 1.0f / (1.0f + art_k *exp (- header.relation_error - header.relation_prob_factor));
				
			}
			
			header.cell_error = reccells[dacell].cell_prob;
			header.cell_prob_factor = reccells[dacell].cell_prob_factor;
			header.cell_prob = 1.0f / (1.0f + art_k *exp (-header.cell_error - header.cell_prob_factor));
			
			header.neckdist = ((reccells[dacell].type == 'c')||(reccells[dacell].type == 'l')||(reccells[dacell].type == 'a')) ? 0.0f : hypot(header.budneck_x  - header.center_x,header.budneck_y - header.center_y);
			
			
			
			header.mc_x = reccells[dacell].centerofmass[0];
			header.mc_y = reccells[dacell].centerofmass[1];
			header.mcdist = hypot(header.mc_x - header.center_x,header.mc_y - header.center_y);
			
			//				
			/*
			 strcpy(buffer,"Red"); printf("\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
			 
			 strcpy(buffer,"Green"); printf("\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
			 strcpy(buffer,"G_SELF_DST"); printf("\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
			 strcpy(buffer,"G_MASC_DST"); printf("\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
			 strcpy(buffer,"G_EDGE_DST"); printf("\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
			 strcpy(buffer,"G_CENT_DST"); printf("\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
			 strcpy(buffer,"G_BUDN_DST"); printf("\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
			 
			 strcpy(buffer,"Intensity"); printf("\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
			 strcpy(buffer,"N_SELF_DST"); printf("\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
			 strcpy(buffer,"N_MASC_DST"); printf("\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
			 strcpy(buffer,"N_EDGE_DST"); printf("\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
			 strcpy(buffer,"N_CENT_DST"); printf("\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
			 strcpy(buffer,"N_BUDN_DST"); printf("\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
			 
			 strcpy(buffer,"R_MASC_DST"); printf("\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);
			 strcpy(buffer,"RG_CROSS_DST"); printf("\t%s_MEAN\t%s_STDV\t%s_SKEW\t%s_KURT", buffer, buffer, buffer, buffer);*/
			
			l=0;
			core_data[l++] = reccells[dacell].red.getMean();
			core_data[l++] = sqrt(reccells[dacell].red.getVar());
			core_data[l++] = reccells[dacell].red.getSkew_scaleinv();
			core_data[l++] = reccells[dacell].red.getKurt_scaleinv();
			
			core_data[l++] = reccells[dacell].green.getMean();
			core_data[l++] = sqrt(reccells[dacell].green.getVar());
			core_data[l++] = reccells[dacell].green.getSkew_scaleinv();
			core_data[l++] = reccells[dacell].green.getKurt_scaleinv();
			
			header.red_green_correlation = (reccells[dacell].red_green.getMean() - core_data[0] * core_data[4]) / (core_data[5] * core_data[1]);
			
			core_data[l++] = reccells[dacell].sd_dist.getMean();
			core_data[l++] = sqrt(reccells[dacell].sd_dist.getVar());
			core_data[l++] = reccells[dacell].sd_dist.getSkew_scaleinv();
			core_data[l++] = reccells[dacell].sd_dist.getKurt_scaleinv();
			
			core_data[l++] = reccells[dacell].cm_dist.getMean();
			core_data[l++] = sqrt(reccells[dacell].cm_dist.getVar());
			core_data[l++] = reccells[dacell].cm_dist.getSkew_scaleinv();
			core_data[l++] = reccells[dacell].cm_dist.getKurt_scaleinv();
			
			core_data[l++] = reccells[dacell].ed_dist.getMean();
			core_data[l++] = sqrt(reccells[dacell].ed_dist.getVar());
			core_data[l++] = reccells[dacell].ed_dist.getSkew_scaleinv();
			core_data[l++] = reccells[dacell].ed_dist.getKurt_scaleinv();
			
			core_data[l++] = reccells[dacell].cn_dist.getMean();
			core_data[l++] = sqrt(reccells[dacell].cn_dist.getVar());
			core_data[l++] = reccells[dacell].cn_dist.getSkew_scaleinv();
			core_data[l++] = reccells[dacell].cn_dist.getKurt_scaleinv();
			
			core_data[l++] = reccells[dacell].bn_dist.getMean();
			core_data[l++] = sqrt(reccells[dacell].bn_dist.getVar());
			core_data[l++] = reccells[dacell].bn_dist.getSkew_scaleinv();
			core_data[l++] = reccells[dacell].bn_dist.getKurt_scaleinv();
			
			// normalized
			
			core_data[l++] = reccells[dacell].intensity.getMean();
			core_data[l++] = sqrt(reccells[dacell].intensity.getVar());
			core_data[l++] = reccells[dacell].intensity.getSkew_scaleinv();
			core_data[l++] = reccells[dacell].intensity.getKurt_scaleinv();
			
			core_data[l] = log(core_data[l-24]) - log(reccells[dacell].rg_sd_dist.getMean());l++;
			core_data[l] = log(core_data[l-24]) - 0.5f * log(reccells[dacell].rg_sd_dist.getVar()) ;l++;
			core_data[l] = core_data[l-24] - reccells[dacell].rg_sd_dist.getSkew_scaleinv();l++;
			core_data[l] = core_data[l-24] - reccells[dacell].rg_sd_dist.getKurt_scaleinv();l++;
			
			core_data[l] = log(core_data[l-24]) - log(reccells[dacell].r_gcm_dist.getMean());l++;
			core_data[l] = log(core_data[l-24]) - 0.5f * log(reccells[dacell].r_gcm_dist.getVar()) ;l++;
			core_data[l] = core_data[l-24] - reccells[dacell].r_gcm_dist.getSkew_scaleinv();l++;
			core_data[l] = core_data[l-24] - reccells[dacell].r_gcm_dist.getKurt_scaleinv();l++;
			
			core_data[l] = log(core_data[l-24]) - log(reccells[dacell].r_ed_dist.getMean());l++;
			core_data[l] = log(core_data[l-24]) - 0.5f * log(reccells[dacell].r_ed_dist.getVar()) ;l++; 
			core_data[l] = core_data[l-24] - reccells[dacell].r_ed_dist.getSkew_scaleinv();l++;
			core_data[l] = core_data[l-24] - reccells[dacell].r_ed_dist.getKurt_scaleinv();l++;
			
			core_data[l] = log(core_data[l-24]) - log(reccells[dacell].r_cn_dist.getMean());l++;
			core_data[l] = log(core_data[l-24]) - 0.5f * log(reccells[dacell].r_cn_dist.getVar()) ;l++;
			core_data[l] = core_data[l-24] - reccells[dacell].r_cn_dist.getSkew_scaleinv();l++;
			core_data[l] = core_data[l-24] - reccells[dacell].r_cn_dist.getKurt_scaleinv();l++;
			
			core_data[l] = log(core_data[l-24]) - log(reccells[dacell].r_bn_dist.getMean());l++;
			core_data[l] = log(core_data[l-24]) - 0.5f * log(reccells[dacell].r_bn_dist.getVar()) ;l++;
			core_data[l] = core_data[l-24] - reccells[dacell].r_bn_dist.getSkew_scaleinv();l++;
			core_data[l] = core_data[l-24] - reccells[dacell].r_bn_dist.getKurt_scaleinv();l++;
			
			
			core_data[l] = log(core_data[l-44]) - log(reccells[dacell].r_sd_dist.getMean());l++;
			core_data[l] = log(core_data[l-44]) - 0.5f * log(reccells[dacell].r_sd_dist.getVar()) ;l++;
			core_data[l] = core_data[l-44] - reccells[dacell].r_sd_dist.getSkew_scaleinv();l++;
			core_data[l] = core_data[l-44] - reccells[dacell].r_sd_dist.getKurt_scaleinv();l++;
			
			core_data[l] = log(core_data[l-44]) - log(reccells[dacell].r_cm_dist.getMean());l++;
			core_data[l] = log(core_data[l-44]) - 0.5f * log(reccells[dacell].r_cm_dist.getVar()) ;l++;
			core_data[l] = core_data[l-44] - reccells[dacell].r_cm_dist.getSkew_scaleinv();l++;
			core_data[l] = core_data[l-44] - reccells[dacell].r_cm_dist.getKurt_scaleinv();l++;
			if (f_mcv != NULL){
				if (distance_thr == 0.0f){
					
					while(mcv_index < header.cellID){
						
						if (mcv_j +1 < damcv.group[mcv_i].cover.size()) mcv_j++;
						else{mcv_j =0; mcv_i++;
							while (mcv_i < damcv.group.size()){
								if (damcv.group[mcv_i].cover.size() == 0) mcv_i++;
								else break;
							}
							while (mcv_i >= damcv.group.size()) {
								damcv.load(f_mcv);mcv_i=0;
								while (mcv_i < damcv.group.size()){
									if (damcv.group[mcv_i].cover.size() == 0) mcv_i++;
									else break;
						}	}	}
						mcv_index++;
					}
					
					header.g_center_x = damcv.group[mcv_i].cover[mcv_j].center[0] + damcv.group[mcv_i].rect[0];
					header.g_center_y = damcv.group[mcv_i].cover[mcv_j].center[1] + damcv.group[mcv_i].rect[1];
					header.g_excentric_x = damcv.group[mcv_i].cover[mcv_j].eccentric[0];
					header.g_excentric_y = damcv.group[mcv_i].cover[mcv_j].eccentric[1];
					header.g_width = damcv.group[mcv_i].cover[mcv_j].width;
					header.g_area = damcv.group[mcv_i].cover[mcv_j].cmpArea();
				}else{
					mcv_i = da_extern_backmap[dacell];
					if (mcv_i != 0xFFFFFFFF){
						header.g_center_x = da_extern_cells[mcv_i].center[0];
						header.g_center_y = da_extern_cells[mcv_i].center[1];
						header.g_excentric_x = da_extern_cells[mcv_i].error[0];
						header.g_excentric_y = da_extern_cells[mcv_i].eccentric[1];
						header.g_width = da_extern_cells[mcv_i].width;
						header.g_area = da_extern_cells[mcv_i].cmpArea();
					}else{
						header.g_center_x = 0.0f;
						header.g_center_y = 0.0f;
						header.g_excentric_x = 0.0f;
						header.g_excentric_y = 0.0f;
						header.g_width = 0.0f;
						header.g_area = 0.0f;
			}	}	}
			
			if (table_out){
				fprintf(datable, "%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",fr_id, header.cellID, header.g_center_x, header.g_center_y, header.g_excentric_x,header.g_excentric_y,header.g_width, header.g_area, header.center_x, header.center_y, header.excentric_x,header.excentric_y,header.width, header.area);
				fprintf(datable, "\t%f\t%f\t%f\t%i\t%f\t%f\t%f\t%c",header.edgedist_radius,header.edgedist_raw_error,header.density_error,(int)header.contour,header.ramanujan,header.fold_red_err,header.red_green_correlation,header.celltype);
				fprintf(datable, "\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%i", header.mc_x,header.mc_y,header.mcdist,header.budneck_x,header.budneck_y,header.neckdist,header.other_major,header.other_area , header.other_dist,header.other_dist_err ,header.other_density,(int)header.other_contour,header.other_ramanujan,header.other_red_dev,header.cell_prob_factor ,header.cell_error ,header.cell_prob , header.relation_prob_factor, header.relation_error, header.relation_prob,header.otherID );
				for(j=0;j<15;j++) fprintf(datable,"\t%f\t%f\t%f\t%f", core_data[(j<<2)], core_data[(j<<2)|1], core_data[(j<<2)|2], core_data[(j<<2)|3]);
				fprintf(datable,"\n");
			}
			daprogbar.update(((double)i) / backmap.size());	
		}daprogbar.finish();
		pixlist.clear();
		maping.clear();
		reccells.clear();
		backmap.clear();
		fr_id++;
	}
	if (table_out) fclose(datable);
return 0;}
void Taskscope<TASK_HIDDENMAP_BASE_DATAEXTRACTION>::help(){
	printf("Identify cell types for cell adjacency. Extract measurements in identified cells. Ratio of distances to point of interests are compared from the 'Green' to the 'Red' channel\n");
	printf("\n");
	printf("\t(in) Image Red\n");
	printf("\t(in) Image Green\n");
	printf("\t(in) Hiddenmap Image\n");
	//		printf("\t(out) distance to other image\n");
	printf("Flags\n\n");
	printf("\t-m (file).mcv: enter MCV file (fills the ellipse guess entries)\n");
	printf("\t-M (file).mcv (float X): filter using MCV filer, maps HIDDEN map to ellipses in MCV, centers must be within X pixels.");
	printf("\t-s : show in tabular format\n");
	printf("\t-c (file.scp): confidence matrix scope file\n");
	printf("\t-S (file): use Segmented image as weights\n");
	//		printf("\t-A: use RED channel as an additionnal channel, and normalize distances from segmented image instead, (weighted by Segmented image if provided)\n")
	printf("\t-t (txt file): output data table\n");
	printf("\t-p (float): artifact class density\n");
	printf("\t-R: Recursively assign Mother-Bud pair on remaining Lone and Daughter Cells");
	printf("\t-C {abmldcMBD}: cell classes identified or outputed, default=\"abmldc\" ('a'rtifact, 'b'ud, 'm'other, 'l'one, 'c'lumped (lone if filtered) , 'd'aughter (lone if filtered), Recursed 'M'other (mother if filtered), Recursed 'B'other (bud if filtered), Recursed 'D'aughter (lone or daughter if filtered)\n");
	printf("\t-B (float fraction): substract a fraction of red intensity to green\n");
	printf("\t-a (float min_cell_confidence = 0.1) (float min surface = 10.0]) (float min desity =0.9) (float max edge error =0.1]) (float max ellipdev pixels =50): define artifact class filter\n");
	printf("\t-r: ignore red total intensity in confidence calculation\n");
	printf("Version 1.0\n");
}


Taskscope<TASK_MAKE_DISPLAYABLE>::Taskscope(): show(false),torun(0),inter(false){ 
	k_dens = 4.510662f -0.165461f; // ram...
	k_dens *= 1.95248835f - 0.035715438f; // intensity
	k_dens *= -0.030004794f - -3.258060922f; // density				
	k_dens *= -0.039659578f - -4.744727495f; // edge err/ 
	k_dens *= 5900.0f;
	k_dens =  0.4404983f / (k_dens * 0.5595016f);					
}
void Taskscope<TASK_MAKE_DISPLAYABLE>::nbaddtoken(char const * const token, int& min, int& max){
	switch(*token){
		case '\0': min =5; break;
		case 'i': min =0; break;
		case 'd': min =1; break;
	}
}
void Taskscope<TASK_MAKE_DISPLAYABLE>::store(char* const * token, int nbtoken){ 
	switch(token[0][1]){
		case 's': show = true; break;
		case 'i': inter = true; break;
		case 'd': k_dens = atof(token[1]); break;
	}
}
int Taskscope<TASK_MAKE_DISPLAYABLE>::defstore(char* const * token, int nbtoken){
	unsigned int i,j,k,l;
	char buffer[65536];
	file_in[0] = token[0];
	file_in[1] = token[1];
	file_in[2] = token[2];
	file_in[3] = token[3];
	file_in[4] = token[4];
	
	
	Madstructs::Table tb;
	tb.wiseload(file_in[0]);
	TiffFile tfh(file_in[1]);
	TiffFile tfr(file_in[2]);
	TiffFile tfg(file_in[3]);
	TiffFile tfo(file_in[4], true);
	printf("%i nbrows \n", tb.nbrows);
	int dataindex[12];	
	i=0;dataindex[i] = tb.findCol("FrameID"); tb.writecoltype(buffer,dataindex[i]);printf("FrameID type: %s\n",buffer);
	i=1;dataindex[i] = tb.findCol("budneck_x"); tb.writecoltype(buffer,dataindex[i]);printf("budneck_x type: %s\n",buffer);
	i=2;dataindex[i] = tb.findCol("budneck_y"); tb.writecoltype(buffer,dataindex[i]);printf("budneck_y type: %s\n",buffer);
	i=3;dataindex[i] = tb.findCol("Cell_type"); tb.writecoltype(buffer,dataindex[i]);printf("Cell_type type: %s\n",buffer);
	i=4;dataindex[i] = tb.findCol("CellID"); tb.writecoltype(buffer,dataindex[i]);printf("CellID type: %s\n",buffer);
	i=5;dataindex[i] = tb.findCol("Cell_prob_log_factor"); tb.writecoltype(buffer,dataindex[i]);printf("Cell_prob_log_factor type: %s\n",buffer);
	i=6;dataindex[i] = tb.findCol("Cell_error"); tb.writecoltype(buffer,dataindex[i]);printf("Cell_error type: %s\n",buffer);
	DataGrid<unsigned int, 2> imh;
	DataGrid<unsigned short, 2> imr;
	DataGrid<unsigned short, 2> img;
	DataGrid<unsigned char, 3> imo;
	unsigned int coor[3];
	double tmp,tmp2;
	int nbind;
	Tuple< Tuple<unsigned int, 2>, 8> ncoor;
	map<unsigned int , KeyElem<double, unsigned char> > type_map;
	for(i=0;i<tb.nbrows;i++){
		tmp2 = exp(tb.getValue(dataindex[5],i).f + tb.getValue(dataindex[6],i).f);
		tmp2 = tmp2 / (tmp2 + k_dens);
		type_map[tb.getValue(dataindex[4],i).i] = KeyElem<double, unsigned char>(tmp2,tb.getValue(dataindex[3],i).sc);
	}
	for(k=0;tfh.fetch(imh);k++){
		if (inter){
			if (!tfr.fetch(img)) {printf("Incoherent number of frames in tiff images '%s' and '%s' (critical)\n", file_in[1], file_in[2]);exit(1);}
			if (!tfr.fetch(imr)) {printf("Incoherent number of frames in tiff images '%s' and '%s' (critical)\n", file_in[1], file_in[2]);exit(1);}
		}else{
			if (!tfr.fetch(imr)) {printf("Incoherent number of frames in tiff images '%s' and '%s' (critical)\n", file_in[1], file_in[2]);exit(1);}
			if (!tfg.fetch(img)) {printf("Incoherent number of frames in tiff images '%s' and '%s' (critical)\n", file_in[1], file_in[3]);exit(1);}
		}
		coor[0] =3;
		coor[1] =imh.dims[0];
		coor[2] =imh.dims[1];
		imo.setSizes(coor);
		DataGrid<unsigned short , 2>::KeyIterator ite = imr.getKeyIterator();
		unsigned short max = ExCo<unsigned int>::minimum();
		if (ite.first()) do{
			if (imr(ite()) > max) max = imr(ite());
		}while (ite.next());
		coor[0] = 0;
		if (ite.first()) do{
			coor[2] = ite()[1];
			coor[1] = ite()[0];
			imo(coor) = (unsigned char)((255.0f*(double)imr(ite())) / max);
		}while (ite.next());
		max = ExCo<unsigned int>::minimum();
		if (ite.first()) do{
			if (img(ite()) > max) max = img(ite());
		}while (ite.next());
		coor[0] = 1;
		if (ite.first()) do{
			coor[2] = ite()[1];
			coor[1] = ite()[0];
			imo(coor) = (unsigned char)((255.0f*(double)img(ite())) / max);
		}while (ite.next());
		coor[0] = 2;
		if (ite.first()) do{
			coor[2] = ite()[1];
			coor[1] = ite()[0];
			nbind = imh.get_indirectNeightbor(ite(), ncoor);
			j = imh(ite());
			if (j != 0){
				l=0;
				for(i=0;i<nbind;i++) if (j != imh(ncoor[i])) l++;
				tmp = ((680.0f*l) / nbind);
				if (type_map[j].d == 'a'){
					if ((((coor[1] - coor[2]+5000) % 6) == 0)||(((coor[1] + coor[2]) % 6) == 0)) tmp = 0xFF;
					if (tmp > 255.0f) imo(coor) = 0xFF;
					else imo(coor) = (unsigned char)tmp;
					coor[0] = 0;
					if (255.0f - tmp < imo(coor)) imo(coor) = 0xFF;
					else imo(coor) += (unsigned char)tmp;
					coor[0] = 1;tmp *0.5f;
					if (255.0f - tmp < imo(coor)) imo(coor) = 0xFF;
					else imo(coor) += (unsigned char)tmp;
					coor[0] = 2;
				}else{
					if (tmp > 255.0f) imo(coor) = 0xFF;
					else imo(coor) = (unsigned char)tmp;
					coor[0] = 0; 
					tmp *= (1.0f - type_map[j].k);
					if (255.0f - tmp < imo(coor)) imo(coor) = 0xFF;
					else imo(coor) += (unsigned char)tmp;
					coor[0] = 2;
				}
			} else {coor[0] = 2; imo(coor) = (unsigned char)0;}
		}while (ite.next());
		for(i=0;i<tb.nbrows;i++){
			if (tb.getValue(dataindex[0],i).i == k) {
				if (tb.getValue(dataindex[3],i).sc == 'b'){
					//	printf("%u\t%u\t%u\n",i, (unsigned int) tb.getValue(dataindex[1],i).f, (unsigned int) tb.getValue(dataindex[2],i).f);
					coor[2] = (unsigned int) tb.getValue(dataindex[2],i).f;
					coor[1] = (unsigned int) tb.getValue(dataindex[1],i).f;
					nbind = imh.get_knightNeightbor(coor+1, ncoor);
					if (nbind == 8){
						coor[1] += 2;
						coor[0] = 0;imo(coor) = 0xFF;coor[0] = 1;imo(coor) = 0xFF;coor[0] = 2;imo(coor) = 0xFF;
						coor[1] -= 4;
						coor[0] = 0;imo(coor) = 0xFF;coor[0] = 1;imo(coor) = 0xFF;coor[0] = 2;imo(coor) = 0xFF;
						coor[1] += 2;
						coor[2] += 2;
						coor[0] = 0;imo(coor) = 0xFF;coor[0] = 1;imo(coor) = 0xFF;coor[0] = 2;imo(coor) = 0xFF;
						coor[2] -= 4;
						coor[0] = 0;imo(coor) = 0xFF;coor[0] = 1;imo(coor) = 0xFF;coor[0] = 2;imo(coor) = 0xFF;
					}
					for(j=0;j<nbind;j++) {
						coor[2] = ncoor[j][1];
						coor[1] = ncoor[j][0];
						coor[0] = 0;imo(coor) = 0xFF;coor[0] = 1;imo(coor) = 0xFF;coor[0] = 2;imo(coor) = 0xFF;
					}
				}else if (tb.getValue(dataindex[3],i).sc == 'd'){
					//		printf("%u\t%u\t%u\n",i, (unsigned int) tb.getValue(dataindex[1],i).f, (unsigned int) tb.getValue(dataindex[2],i).f);
					coor[2] = (unsigned int) tb.getValue(dataindex[2],i).f;
					coor[1] = (unsigned int) tb.getValue(dataindex[1],i).f;
					nbind = imh.get_knightNeightbor(coor+1, ncoor);
					if (nbind == 8){
						coor[1] += 2;
						coor[0] = 0;imo(coor) = 0x40;coor[0] = 1;imo(coor) = 0xFF;coor[0] = 2;imo(coor) = 0x40;
						coor[1] -= 4;
						coor[0] = 0;imo(coor) = 0x40;coor[0] = 1;imo(coor) = 0xFF;coor[0] = 2;imo(coor) = 0x40;
						coor[1] += 2;
						coor[2] += 2;
						coor[0] = 0;imo(coor) = 0x40;coor[0] = 1;imo(coor) = 0xFF;coor[0] = 2;imo(coor) = 0x40;
						coor[2] -= 4;
						coor[0] = 0;imo(coor) = 0x40;coor[0] = 1;imo(coor) = 0xFF;coor[0] = 2;imo(coor) = 0x40;
					}
					for(j=0;j<nbind;j++) {
						coor[2] = ncoor[j][1];
						coor[1] = ncoor[j][0];
						coor[0] = 0;imo(coor) = 0x40;coor[0] = 1;imo(coor) = 0xFF;coor[0] = 2;imo(coor) = 0x40;
					}
				}
			}
		}
		tfo.put(imo,(unsigned char)0,(unsigned char)255);
	}
	
	
return 0;}
void Taskscope<TASK_MAKE_DISPLAYABLE>::help(){
	printf("Make a Condifence Profile.\n");
	printf("Default Arguments 1-2:\n");
	printf("(in) data table\n");
	printf("(in) hidden_map image\n");
	printf("(in) raw Red image\n");
	printf("(in) raw Green image\n");
	printf("(out) rendered image\n");
	printf("\n");
	printf("Optionnal Default arguments(3):\n");
	printf("-i: assumes that RED/GREEN are in the same file!\n");
	printf("-d: set artifact time piror ratio constant\n");
	printf("\t-h : Help \n");
}


Taskscope<TASK_GETLABEL_CONFIDENCE>::Taskscope(): show(false),torun(0){
}
void Taskscope<TASK_GETLABEL_CONFIDENCE>::nbaddtoken(char const * const token, int& min, int& max){
	switch(*token){
		case '\0': min =1; max =2; break;
		case 's': min =0; break;
		case 'c': min =0; break;
		case 'l': min =0; break;
	}
}
void Taskscope<TASK_GETLABEL_CONFIDENCE>::store(char* const * token, int nbtoken){ 
	switch(token[0][1]){
		case 's': show = true; break;
		case 'c': torun |=1; break;
		case 'l': torun |=2; break;
	}
}
int Taskscope<TASK_GETLABEL_CONFIDENCE>::defstore(char* const * token, int nbtoken){
	unsigned int i,j,k,l;
	char buffer[65536];
	file_in[0] = (nbtoken > 0) ? token[0] : NULL;
	file_in[1] = (nbtoken > 1) ? token[1] : NULL;
	if ((nbtoken < 2)&&(!show)) {show = true; printf("Missing arguments, showing the current scope (no modifications)\n");}
	SerialStore<unsigned int> dastore(file_in[0]);
	Tuple< KeyElem<double, GaussianDistribution<4> > , 7 > keyframe;
	GaussianDistribution<1> size_key;
	if (show){
		dastore.load(100,size_key);
		size_key.show();
		dastore.load(101,keyframe);
		keyframe.show();
		return(0);
	}
	/*	
	 FILE* flist  = fopen("../../../../feat_data6_list.txt","r+");
	 Tuple<int, 5> tcbuffer;ExOp::toZero(tcbuffer);
	 Tuple<int, 5> ucbuffer;ExOp::toZero(ucbuffer);
	 LFHPrimitive::Classifier<Tuple<double, 5>, 4> mega_classif;
	 LFHPrimitive::Classifier<Tuple<double, 1>, 4> mega_classif_s;
	 GaussianDistribution<5> typecl[4];
	 GaussianDistribution<1> typecl_s[4];
	 Tuple<double, 5> in_data;
	 FILE* g;
	 typecl[0].EMinit();
	 typecl[1].EMinit();
	 typecl[2].EMinit();
	 typecl[3].EMinit();
	 j=0;
	 while (fscanf(flist,"%s\n",buffer)== 1){
	 Madstructs::Table tb(buffer,"i i f f f f f f f f f f f f f c f f f f f f f f f f f f f f i f ");
	 Tuple<int, 5> cbuffer;
	 ExOp::toZero(cbuffer);
	 for(i=0;i<tb.nbrows;i++){
	 in_data[0] = log(tb.getValue(7,i).f);
	 in_data[1] = log(tb.getValue(9,i).f);
	 in_data[2] = log(1.0f - tb.getValue(10,i).f); if (!ExCo<double>::isValid(in_data[2])) in_data[2] = -1000.0f;
	 in_data[3] = tb.getValue(12,i).f;
	 in_data[4] = log(tb.getValue(31,i).f); // alt13
	 if (in_data[2] > -1000.0f){
	 switch(tb.getValue(15,i).sc){
	 case 'm': cbuffer[0]++; typecl[0].EMregist(in_data); break;
	 case 'b': cbuffer[1]++; typecl[1].EMregist(in_data); break;
	 case 'd': cbuffer[2]++; typecl[2].EMregist(in_data); break;
	 case 'l': cbuffer[3]++; typecl[3].EMregist(in_data); break;
	 case 'a': cbuffer[4]++; break;
	 }}
	 }
	 }
	 typecl[0].EMfinit();
	 typecl[1].EMfinit();
	 typecl[2].EMfinit();
	 typecl[3].EMfinit();
	 printf("Mother!\n");
	 typecl[0].show();
	 printf("Bud!\n");
	 typecl[1].show();
	 printf("Daughter!\n");
	 typecl[2].show();
	 printf("Lone!\n");
	 typecl[3].show();
	 //	g = fopen("./4_gaussian_5.bin", "w+");			typecl[0].save(g);			typecl[1].save(g);			typecl[2].save(g);			typecl[3].save(g);			fclose(g);
	 //	g = fopen("./4_gaussian_5.bin", "r+");			typecl[0].load(g);			typecl[1].load(g);			typecl[2].load(g);			typecl[3].load(g);			fclose(g);
	 Tuple<unsigned int, 1 > dachan; dachan[0] = 0;
	 printf("Mother!\n");
	 typecl[0].show();
	 printf("Bud!\n");
	 typecl[1].show();
	 printf("Daughter!\n");
	 typecl[2].show();
	 printf("Lone!\n");
	 typecl[3].show();
	 //	typecl_s[0] = typecl[0].makeSubGaussian(dachan);
	 //	typecl_s[1] = typecl[1].makeSubGaussian(dachan);
	 //	typecl_s[2] = typecl[2].makeSubGaussian(dachan);
	 //	typecl_s[3] = typecl[3].makeSubGaussian(dachan);
	 //	Tuple<unsigned int ,2> hh;
	 //	TMatrix<double, 4 ,5> haha; haha.getDims(hh); ExOp::show(hh);
	 //	printf("Mother!\n");
	 //	typecl_s[0].show();
	 //	printf("Bud!\n");
	 //	typecl_s[1].show();
	 //	printf("Daughter!\n");
	 //	typecl_s[2].show();
	 //	printf("Lone!\n");
	 //	typecl_s[3].show();
	 mega_classif[0] = &(typecl[0]);
	 mega_classif[1] = &(typecl[1]);
	 mega_classif[2] = &(typecl[2]);
	 mega_classif[3] = &(typecl[3]);
	 mega_classif_s[0] = &(typecl_s[0]);
	 mega_classif_s[1] = &(typecl_s[1]);
	 mega_classif_s[2] = &(typecl_s[2]);
	 mega_classif_s[3] = &(typecl_s[3]);
	 Tuple<double, 4> likeli;
	 Tuple<double, 25> likbuf; ExOp::toZero(likbuf);
	 GradientSearchScope grscp; grscp.init(0.000000001f, 4);
	 double scope[4];
	 Tuple<double, 4> deriv;
	 double total;
	 double total_ll;
	 ExOp::toZero(scope);
	 scope[0] = 1.937852e+01;
	 scope[1] = 2.103694e+01;
	 scope[2] = 2.183531e+01;
	 scope[3] = 1.445520e+01;
	 DataGrid<double, 2> size_buf;
	 unsigned int coor[2];
	 coor[0] = 250;
	 coor[1] = 5;
	 size_buf.setSizes(coor); ExOp::toZero( size_buf);
	 double best;
	 int best_cl;
	 for(k=0;k<1;k++){
	 ExOp::toZero(deriv);total_ll =0.0f;
	 fseek( flist , 0 , SEEK_SET );
	 while (fscanf(flist,"%s\n",buffer)== 1){
	 Madstructs::Table tb(buffer,"i i f f f f f f f f f f f f f c f f f f f f f f f f f f f f i f ");
	 for(i=0;i<tb.nbrows;i++){
	 in_data[0] = log(tb.getValue(7,i).f);
	 in_data[1] = tb.getValue(9,i).f;
	 in_data[2] = log(1.0f - tb.getValue(10,i).f); if (!ExCo<double>::isValid(in_data[2])) in_data[2] = -1000.0f;
	 in_data[3] = tb.getValue(12,i).f;
	 in_data[4] = log(tb.getValue(31,i).f); // alt13
	 if (in_data[2] > -1000.0f){
	 (*mega_classif[0])(likeli[0], in_data);
	 (*mega_classif[1])(likeli[1], in_data);
	 (*mega_classif[2])(likeli[2], in_data);
	 (*mega_classif[3])(likeli[3], in_data);
	 total = likeli[0] * exp(scope[0]);
	 total += likeli[1] * exp(scope[1]);
	 total += likeli[2] * exp(scope[2]);
	 total += likeli[3] * exp(scope[3]);
	 //		printf("%e\n",likeli[2] * exp(scope[2]));
	 if (ExCo<double>::isValid(total)){
	 //	if (!ExCo<double>::isValid( log(likeli[0]))) printf("%e\n",typecl[0].LL(in_data));
	 //	if (!ExCo<double>::isValid( log(likeli[1]))) printf("%e\n",typecl[1].LL(in_data));
	 //	if (!ExCo<double>::isValid( log(likeli[2]))) printf("%e\n",typecl[2].LL(in_data));
	 //	if (!ExCo<double>::isValid( log(likeli[3]))) printf("%e\n",typecl[3].LL(in_data));
	 //		printf("%f\t%f\t%f\t%f\t%f\t",  likeli[0] * exp(scope[0]) / (1.0f + total),  likeli[1] * exp(scope[1]) / (1.0f + total),  likeli[2] * exp(scope[2]) / (1.0f + total),  likeli[3] * exp(scope[3]) / (1.0f + total),  1.0f / (1.0f + total));
	 switch(tb.getValue(15,i).sc){
	 case 'm': 
	 //	if( k ==0)	   printf("%c: %f\n", tb.getValue(15,i).sc, likeli[0] * exp(scope[0]) / (1.0f + total));
	 for(l=0;l<4;l++) likbuf[0+l] += likeli[l] * exp(scope[l]);
	 best = 1.0f; best_cl =4;
	 for(l=0;l<4;l++) if (best < likeli[l] * exp(scope[l])) {best_cl =l; best = likeli[l] * exp(scope[l]);}
	 likbuf[0+best_cl] += 1.0f;
	 deriv[0] +=  exp(-scope[0]);
	 total_ll += scope[0] + log(likeli[0]);
	 coor[1] = 0;
	 break;
	 case 'b': for(l=0;l<4;l++) likbuf[4+l] += likeli[l] * exp(scope[l]);
	 //	if( k ==0)	   printf("%c: %f\n", tb.getValue(15,i).sc, likeli[1] * exp(scope[1]) / (1.0f + total));
	 best = 1.0f; best_cl =4;
	 for(l=0;l<4;l++) if (best < likeli[l] * exp(scope[l])) {best_cl =l; best = likeli[l] * exp(scope[l]);}
	 likbuf[5+best_cl] += 1.0f;
	 deriv[1] +=  exp(-scope[1]);
	 total_ll += scope[1] + log(likeli[1]);
	 coor[1] = 1;
	 break;
	 case 'd': for(l=0;l<4;l++) likbuf[8+l] += likeli[l] * exp(scope[l]);
	 //if( k ==0)	   printf("%c: %f\n", tb.getValue(15,i).sc, likeli[2] * exp(scope[2]) / (1.0f + total));
	 best = 1.0f; best_cl =4;
	 for(l=0;l<4;l++) if (best < likeli[l] * exp(scope[l])) {best_cl =l; best = likeli[l] * exp(scope[l]);}
	 likbuf[10+best_cl] += 1.0f;
	 deriv[2] +=  exp(-scope[2]);
	 total_ll += scope[2] + log(likeli[2]);
	 coor[1] = 2;
	 break;
	 case 'l': for(l=0;l<4;l++) likbuf[12+l] += likeli[l] * exp(scope[l]);
	 //		   if( k ==0)	   printf("%c: %f\n", tb.getValue(15,i).sc, likeli[3] * exp(scope[3]) / (1.0f + total));
	 best = 1.0f; best_cl =4;
	 for(l=0;l<4;l++) if (best < likeli[l] * exp(scope[l])) {best_cl =l; best = likeli[l] * exp(scope[l]);}
	 likbuf[15+best_cl] += 1.0f;
	 deriv[3] +=  exp(-scope[3]);
	 total_ll += scope[3]+ log(likeli[3]);
	 coor[1] = 3;
	 break;
	 case 'a': for(l=0;l<4;l++) likbuf[16+l] += likeli[l] * exp(scope[l]);
	 //if( k ==0)	   printf("%c: %f\n", tb.getValue(15,i).sc, 1.0f / (1.0f + total));
	 best = 1.0f; best_cl =4;
	 for(l=0;l<4;l++) if (best < likeli[l] * exp(scope[l])) {best_cl =l; best = likeli[l] * exp(scope[l]);}
	 likbuf[20+best_cl] += 1.0f;
	 coor[1] = 4;
	 }
	 //						coor[0] = tb.getValue(7,i).f / 10.0f;
	 // coor[0] = tb.getValue(9,i).f * 30000.0f;
	 //						coor[0] =  tb.getValue(10,i).f * 230.0f;
	 coor[0] = ((tb.getValue(12,i).f * 200.0f) - 75.0f);
	 //						coor[0] = ((tb.getValue(12,i).f * 1000.0f) + 125.0f);
	 if (coor[0] > 1000249) coor[0] = 0;
	 if (coor[0] > 249) coor[0] = 249;
	 size_buf(coor) += 1.0f;
	 deriv[0] -= likeli[0]  / (total + 1.0f);
	 deriv[1] -= likeli[1]  / (total + 1.0f);
	 deriv[2] -= likeli[2]  / (total + 1.0f);
	 deriv[3] -= likeli[3]  / (total + 1.0f);
	 total_ll  -= log(total+ 1.0f);
	 }
	 }
	 }
	 }
	 ExOp::show(deriv);
	 deriv[0] *= exp(scope[0]);
	 deriv[1] *= exp(scope[1]);
	 deriv[2] *= exp(scope[2]);
	 deriv[3] *= exp(scope[3]);
	 ExOp::show(deriv);
	 ExOp::show(scope);
	 printf("LL = %e step = %e\n",total_ll, grscp.updateAscent(total_ll, scope,  deriv.data));
	 ExOp::show(scope);
	 }
	 ExOp::show(likbuf);
	 //	ExOp::show(size_buf);
	 */
	Madstructs::Table tb;
	tb.wiseload(file_in[1]);
	printf("%i nbrows \n", tb.nbrows);
	//		Madstructs::Table tb(file_in[1] ,"i s c f f f f f f f f f f f f f f f i f f f c s f f f f f f f f f f f f f i f f f f i f " );
	//Madstructs::Table tb("../../../../report/complete_jeb4.txt","i s i f f f f c f f f f f f f f f f f f f c f f f f f f f f f f f f f f f f i f " );
	int dataindex[12];
	if (torun & 1){
		i=0;dataindex[i] = tb.findCol("Area"); tb.writecoltype(buffer,dataindex[i]);printf("Area type: %s\n",buffer);
		i=1;dataindex[i] = tb.findCol("EdgeDistance_Error"); tb.writecoltype(buffer,dataindex[i]);printf("EdgeDistance_Error type: %s\n",buffer);
		i=2;dataindex[i] = tb.findCol("Density_In_Ellipse"); tb.writecoltype(buffer,dataindex[i]);printf("Density_In_Ellipse type: %s\n",buffer);
		i=3;dataindex[i] = tb.findCol("Ramanujan"); tb.writecoltype(buffer,dataindex[i]);printf("Ramanujan type: %s\n",buffer);
		i=4;dataindex[i] = tb.findCol("Red_MEAN"); tb.writecoltype(buffer,dataindex[i]);printf("Red_MEAN type: %s\n",buffer);
		GaussianDistribution<4> typecl[7];
		double mix;
		Tuple<double, 7> dar; ExOp::toZero(dar);
		int first;
		Tuple<double, 4> point;
		Tuple<double, 1> point_size;
		for(i=0;i<7;i++) {typecl[i].EMinit();} 
		size_key.EMinit();
		for(i=0;i<tb.nbrows;i++){
			point_size[0] = tb.getValue(dataindex[0],i).f; size_key.EMregist(point_size);
			mix = point_size[0] / 500.0f;
			first = floor(mix);
			if (first >= 6) {first = 5; mix = 1.0f;}
			else mix -= first;
			point[0] = log(tb.getValue(dataindex[1],i).f);
			point[1] = log(1.0f - tb.getValue(dataindex[2],i).f); 
			point[2] = tb.getValue(dataindex[3],i).f;
			point[3] = log(tb.getValue(dataindex[4],i).f);
			if (ExCo<double>::isValid(point[0])&&ExCo<double>::isValid(point[1])&&ExCo<double>::isValid(point[2])&&ExCo<double>::isValid(point[3])) {
				dar[first] += 1.0f - mix;
				dar[first+1] += mix;
				typecl[first].EMregist(point, 1.0f - mix);
				typecl[first+1].EMregist(point,  mix);
			} else printf("Fail: ");
			ExOp::show(point);
		}
		size_key.EMfinit();
		for(i=0;i<7;i++) {typecl[i].EMfinit();  keyframe[i].k = 500.0f * i;keyframe[i].d = typecl[i]; } 
		dastore.save(100,size_key);
		dastore.save(101,keyframe);
		size_key.show();
		keyframe.show();
		TMatrix<double,4,4> inv = keyframe[0].d.ihvar.inverse(); inv.show(); 
		TMatrix<double,4,4> inv2 = keyframe[1].d.ihvar.inverse(); inv2.show(); 
		TMatrix<double,4,4> inv3 = (inv + inv2) * 0.5f; inv3.show(); 
		ExOp::show(-2.5f*log(M_PI) + 0.5f* (log(fabs(size_key.ihvar.data[0])) -log(fabs(inv3.determinant()))) );
		inv = keyframe[3].d.ihvar.inverse(); inv.show(); 
		inv2 = keyframe[4].d.ihvar.inverse(); inv2.show(); 
		inv3 = (inv*0.1f) + (inv2 * 0.9f); inv3.show(); 
		ExOp::show(-2.5f*log(M_PI) + 0.5f* (log(fabs(size_key.ihvar.data[0])) -log(fabs(inv3.determinant()))));			
	}
	// LEARNING THE MIXING PARAMETER!
	if (torun & 2){
		i=0;dataindex[i] = tb.findCol("Area"); tb.writecoltype(buffer,dataindex[i]);printf("Area type: %s\n",buffer);
		i=1;dataindex[i] = tb.findCol("EdgeDistance_Error"); tb.writecoltype(buffer,dataindex[i]);printf("EdgeDistance_Error type: %s\n",buffer);
		i=2;dataindex[i] = tb.findCol("Density_In_Ellipse"); tb.writecoltype(buffer,dataindex[i]);printf("Density_In_Ellipse type: %s\n",buffer);
		i=3;dataindex[i] = tb.findCol("Ramanujan"); tb.writecoltype(buffer,dataindex[i]);printf("Ramanujan type: %s\n",buffer);
		i=4;dataindex[i] = tb.findCol("Red_MEAN"); tb.writecoltype(buffer,dataindex[i]);printf("Red_MEAN type: %s\n",buffer);
		double mix = 0.0005f;
		double inp_dens;
		double dafr;
		Vector<double> densities;
		printf("list of densities (start)\n");
		double pix[32];
		i=0;
		pix[0] = pix[1] = tb.getValue(dataindex[1],0).f;
		pix[2] = pix[3] = tb.getValue(dataindex[2],0).f;
		pix[4] = pix[5] = tb.getValue(dataindex[3],0).f;
		pix[6] = pix[7] = tb.getValue(dataindex[4],0).f;
		//				while(2 == fscanf(f,"%[^\n\t\r]%c",buffer+prefix_length, &sep)){				}
		for(i=1;i<tb.nbrows;i++){
			for(j=0;j<4;j++){
				if (pix[j*2] >  tb.getValue(dataindex[1],0).f) pix[j*2] =  tb.getValue(dataindex[1],0).f;
				if (pix[1|(j*2)] <  tb.getValue(dataindex[1],0).f) pix[1|(j*2)] =  tb.getValue(dataindex[1],0).f;
			}
			inp_dens = exp(tb.getValue(dataindex[1],i).f + tb.getValue(dataindex[0],i).f);
		}
		double density = 5900.0f;
		for(j=0;j<4;j++) density *= pix[1|(j*2)] - pix[j*2]; 
		density = 1.0f / density;
		printf("%e", log(density));
		printf("list of densities (end)\n");
		for(l=0;l<25;l++){
			double in=0.0f; int count=0;
			double ll=0.0f;
			for(i=0;i<tb.nbrows;i++){
				inp_dens = tb.getValue(dataindex[1],i).f + tb.getValue(dataindex[0],i).f;
				if (ExCo<double>::isValid(inp_dens)){
					dafr = mix *density / ( (1.0f-mix)* exp(inp_dens) + mix *density);
					ll += dafr * log(density) + (1.0f - dafr) * inp_dens;
					if (ExCo<double>::isValid(inp_dens)) {count++; in +=dafr;}
				}
			}
			printf("LL[%e] = %e\n", mix, ll);
			mix = in / count;
		}
		for(i=0;i<tb.nbrows;i++){
			inp_dens = tb.getValue(dataindex[1],i).f + tb.getValue(dataindex[0],i).f;
			if (ExCo<double>::isValid(inp_dens)){
				printf("%e\n", mix *density / ( (1.0f-mix)* exp(inp_dens) + mix *density));
			}
		}
	}
	/*
	 if (torun & 2){
	 i=0;dataindex[i] = tb.findCol("Cell_prob_log_factor");tb.writecoltype(buffer,dataindex[i]);printf("Cell_prob_log_factor (%s)\n",buffer);// 43
	 i=1;dataindex[i] = tb.findCol("Cell_error");tb.writecoltype(buffer,dataindex[i]);printf("Cell_error (%s)\n",buffer);// 43
	 i=2;dataindex[i] = tb.findCol("Cell_type"); tb.writecoltype(buffer,dataindex[i]);printf("Cell_type type: %s\n",buffer);
	 double density = 4.510662f -0.165461f; // ram...
	 density *= 1.95248835f - 0.035715438f; // intensity
	 density *= -0.030004794f - -3.258060922f; // density				
	 density *= -0.039659578f - -4.744727495f; // edge err/ 
	 density *= 5900.0f;
	 density = 1.0f / density;
	 printf("%e", log(density));
	 double mix = 0.0005f;
	 double inp_dens;
	 double dafr;
	 density *= 1.0f;
	 printf("list of densities (start)\n");
	 for(i=0;i<tb.nbrows;i++){
	 inp_dens = exp(tb.getValue(dataindex[1],i).f + tb.getValue(dataindex[0],i).f);
	 printf("%e vs %e\n", inp_dens,density);
	 }
	 printf("list of densities (end)\n");
	 for(l=0;l<25;l++){
	 double in=0.0f; int count=0;
	 double ll=0.0f;
	 for(i=0;i<tb.nbrows;i++){
	 inp_dens = tb.getValue(dataindex[1],i).f + tb.getValue(dataindex[0],i).f;
	 if (ExCo<double>::isValid(inp_dens)){
	 dafr = mix *density / ( (1.0f-mix)* exp(inp_dens) + mix *density);
	 ll += dafr * log(density) + (1.0f - dafr) * inp_dens;
	 if (ExCo<double>::isValid(inp_dens)) {count++; in +=dafr;}
	 }
	 }
	 printf("LL[%e] = %e\n", mix, ll);
	 mix = in / count;
	 }
	 for(i=0;i<tb.nbrows;i++){
	 inp_dens = tb.getValue(dataindex[1],i).f + tb.getValue(dataindex[0],i).f;
	 if (ExCo<double>::isValid(inp_dens)){
	 printf("%e\n", mix *density / ( (1.0f-mix)* exp(inp_dens) + mix *density));
	 }
	 }
	 }
	 */
	/*
	 ExOp::show(size_key);
	 Tuple<double,25> confusion; ExOp::toZero(confusion);
	 for(i=0;i<tb.nbrows;i++){
	 mix = 1.05f / ((1.0e-12 / tb.getValue(dataindex[5],i).f) + 1);
	 switch(tb.getValue(dataindex[7],i).sc){
	 case 'm': k=0; break;
	 case 'b': k=5; break;
	 case 'd': k=10; break;
	 case 'l': k=15; break;
	 case 'a': k=20; break;
	 }
	 switch(tb.getValue(dataindex[8],i).sc){
	 case 'm': confusion[k] += mix; confusion[k+4] += 1.0f-mix; break;
	 case 'b': confusion[k+1] += mix; confusion[k+4] += 1.0f-mix; break;				
	 case 'd': confusion[k+2] += mix; confusion[k+4] += 1.0f-mix; break;				
	 case 'l': confusion[k+3] += mix; confusion[k+4] += 1.0f-mix; break;				
	 case 'a': confusion[k+4] += 1.0f; break;				
	 }
	 }
	 printf("%i\n", tb.nbrows);
	 ExOp::show(confusion);ExOp::toZero(confusion);
	 for(i=0;i<tb.nbrows;i++){
	 //	point[0] = tb.getValue(dataindex[5],i).f / 0.925f;
	 //	point[1] = tb.getValue(dataindex[6],i).f / 0.925f;
	 point[0] =1.04f / ((1.0e-10 / exp(tb.getValue(dataindex[5],i).f+tb.getValue(dataindex[9],i).f)) + 1);
	 point[1] =1.04f / ((1.0e-10 / exp(tb.getValue(dataindex[6],i).f+tb.getValue(dataindex[10],i).f)) + 1);
	 switch(tb.getValue(dataindex[7],i).sc){
	 case 'm': k=0; break;
	 case 'b': k=5; break;
	 case 'd': k=10; break;
	 case 'l': k=15; break;
	 case 'a': k=20; break;
	 }
	 switch(tb.getValue(dataindex[8],i).sc){
	 case 'm': confusion[k] += point[0]*point[1]; confusion[k+3] += point[0] * (1.0f - point[1]); confusion[k+4] += 1.0f-point[0]; break;
	 case 'b': confusion[k+1] += point[0]*point[1]; confusion[k+3] += point[0] * (1.0f - point[1]);confusion[k+4] += 1.0f-point[0]; break;				
	 case 'd': confusion[k+2] += point[0]*point[1]; confusion[k+3] += point[0] * (1.0f - point[1]);confusion[k+4] += 1.0f-point[0]; break;				
	 case 'l': confusion[k+3] += point[0]; confusion[k+4] += 1.0f-point[0]; break;				
	 case 'a': confusion[k+4] += 1.0f; break;				
	 }
	 }
	 ExOp::show(confusion);
	 */
	return 0;}
void Taskscope<TASK_GETLABEL_CONFIDENCE>::help(){
	printf("Make a Condifence Profile.\n");
	printf("Default Arguments 1-2:\n");
	printf("(out) Scope file\n");
	printf("(in) tabular file\n");
	printf("\n");
	printf("Optionnal Default arguments(3):\n");
	printf("-s: Only show the content of the scopefile (no modifications), automatically on no tabular file provided.\n");
	printf("-c: Compute confidence profile, required collumns: (JebType,Cell_type,Area,EdgeDistance_Error,Density_In_Ellipse,Ramanujan,Red_MEAN,Relation_prob_log_factor,Relation_error,Cell_prob_log_factor,Cell_error)\n");	
	printf("-a: Make a Displayable annotated image\n");	
	printf("-l (table file list): learn K, and fill in the cell probabilities: (Cell_type,Area,EdgeDistance_Error,Density_In_Ellipse,Ramanujan,Red_MEAN,Relation_prob_log_factor,Relation_error,Cell_prob_log_factor,Cell_error)\n");	
	printf("\t-h : Help \n");
}








