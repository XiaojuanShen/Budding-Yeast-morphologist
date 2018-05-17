/*
 * Pvalue.cpp
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

	
	Taskscope<TASK_CDTGTR_PVALUES>::Taskscope():threshold(0.0f),show(false),filter_file(NULL),permute_file(NULL),torun(0),inter(false),w_name(false),gray(false),scaleR(false),scaleG(false),scaleB(false),factR(0.0f),factG(0.0f),factB(0.0f),doopen(false),dodel(false),cdtcolname(NULL),nb_permutations(0) {
		nb_annot_range[0] =2; nb_annot_range[1] = 0xFFFFFFFF;
	}
	void Taskscope<TASK_CDTGTR_PVALUES>::nbaddtoken(char const * const token, int& min, int& max){
		switch(*token){
			case '\0': min =0; max =1; break;
			case 'p': min =1; max =1; break;
			case 't': min =2; break;
			case 'A': min =0; break;
			case 'g': min =0; break;
			case 'S': min =0; break;
			case 'F': min =1; break;
			case 'f': min =2; break;
			case 'N': min =1; break;
			case 'P': min =1; max=2; break;
		}
	}
	
	void Taskscope<TASK_CDTGTR_PVALUES>::store(char* const * token, int nbtoken){ 
		switch(token[0][1]){
			case 's': show = true; break;
			case 't': threshold_path = token[1] ; threshold = atof(token[2]); break;
			case 'p': inter = true; break;
			case 'g': gray = true; break;
			case 'A': w_name = true; break;
			case 'S': switch(token[0][2]){
				case 'R':scaleR = true; break;
				case 'G':scaleG = true; break;
				case 'B':scaleB = true; break;
				default: scaleR = true;scaleG = true;scaleB = true;
			}break;
			case 'F': filter_file = token[1];break;
			case 'f': nb_annot_range[0] = atoi(token[1]); nb_annot_range[1] = atoi(token[2]); break;
			case 'P': nb_permutations =atoi(token[1]); if (nbtoken > 1) permute_file = token[2]; break;
			case 'N': cdtcolname = token[1];break;
		}
	}
	int Taskscope<TASK_CDTGTR_PVALUES>::defstore(char* const * token, int nbtoken){ 
		unsigned int i,j,k,l;
		char buffer[65536];
		double pix[32];
		
		SerialStore<unsigned int> rstd("PMGTRpvalues_scope.scp");
		Vector<char> str;
		
		if (inter){
			printf("enter the annotation table path:");
			if (rstd.has(0)){
				rstd.load(0, str);
				printf(" currently =%s\n", &(str[0]));
			} else printf("\n");
			scanf(" %[^\n]", buffer);
			l = strlen(buffer);
			if ((l)&&(strcmp(buffer," ") != 0)){
				str.setSize(l+1);
				memcpy(&(str[0]),buffer, l+1);
				rstd.save(0, str);
			}
			printf("enter the annotation description file:");
			if (rstd.has(1)){
				rstd.load(1, str);
				printf(" currently =%s\n", &(str[0]));
			} else printf("\n");
			scanf(" %[^\n]", buffer);
			l = strlen(buffer);
			if ((l)&&(strcmp(buffer," ") != 0)){
				str.setSize(l+1);
				memcpy(&(str[0]),buffer, l+1);
				rstd.save(1, str);
			}
		}else{
			
			
			i = strlen(token[0]); memcpy(buffer,token[0],i); memcpy(buffer+i,".cdt",5);
			
			FILE* f_names = fopen(buffer,"r+"); if (f_names == NULL) {fprintf(stdout,"cannot open %s! (critical)\n", buffer);exit(1);}
			memcpy(buffer+i,".gtr",5);
			FILE* f_trees = fopen(buffer,"r+"); if (f_trees == NULL) {fprintf(stdout,"cannot open %s! (critical)\n", buffer);exit(1);}
			
			
			
			Vector< pair<double, pair<unsigned int,unsigned int > > > best_pvalue; 
			Vector<unsigned int> annotID_counts;
			myHashmap<string, unsigned int> annot_to_annotID;
			
			double dval;
			double da_super_sum;
			Vector< Vector<unsigned int> > annotIDs;
			myHashmap<string, unsigned int > name_to_internal;
			myHashmap<string, unsigned int > cdtname_to_internal;
			
			char* name = NULL;
			Tuple<char,32> ID;
			char sep;
			Vector<unsigned int> annot_input;
			j=0; // annotID count
			
			myHashmap<string,char> filter;
			
			FILE* g_annot;
			if (filter_file){
				printf("loading filter file..."); fflush(stdout);
				g_annot = fopen(filter_file,"r+"); if (g_annot == NULL) {fprintf(stdout,"cannot open %s! (critical)\n", filter_file);exit(1);}
				while(1 == fscanf(g_annot,"%[^\t\n]%*c", buffer)) filter[cloneString(buffer)] =0;
				fclose(g_annot);
				printf("(DONE)\n"); fflush(stdout);
			}
			
			rstd.load(0, str);
			j = strlen(&(str[0]));
			memcpy(buffer,&(str[0]), j);
			buffer[j] = '\0';
			g_annot = fopen(buffer,"r+"); if (g_annot == NULL) {fprintf(stdout,"cannot open %s! (critical)\n", buffer);exit(1);}
			
			// LOADS ALL ANNOTATIONS! 
			printf("loading annotations...");fflush(stdout);
			k =0;
			while(2 == fscanf(g_annot,"%[^\t\n]%c", buffer,&sep)){
				if (name == NULL){
					name = cloneString(buffer);
				}else{
					if ((filter_file == NULL)||(filter.find(buffer) != 0xFFFFFFFF)){
						i = annot_to_annotID.find(buffer);
						if (i == 0xFFFFFFFF){
							
							annot_to_annotID[cloneString(buffer)] = j++;
							i = annot_to_annotID.find(buffer);
							if (i == 0xFFFFFFFF) {printf("NOSE!\n");exit(0);}
						}
						annot_input.push_back(annot_to_annotID.deref(i));
					}
				}
				if (sep == '\n'){
					name_to_internal[string(name)] = k;
					annot_input.sort();
					annotIDs.push_back(annot_input);
					k++;
					delete[](name);
					name = NULL;
					annot_input.clear();
				}
			}
			annotIDs.push_back(annot_input); // empty annotations
			
		//	for(l=0;l< annotIDs.size();l++) printf("size=%i\n", annotIDs[l].size() );
			
			fclose(g_annot);
			printf("(DONE) (%i annotation total)\n", j); fflush(stdout);
			Vector<char*> annot_descr; 
			
			Vector<char*> node_names; node_names.setSize(k);ExOp::toZero(node_names);
			
			annotID_counts.setSize(j);ExOp::toZero(annotID_counts);
			best_pvalue.setSize(j);ExOp::toZero(best_pvalue);
			annot_descr.setSize(j);ExOp::toZero(annot_descr);
			
			
			rstd.load(1, str);
			j = strlen(&(str[0]));
			memcpy(buffer,&(str[0]), j);
			buffer[j] = '\0';
			g_annot = fopen(buffer,"r+"); if (g_annot == NULL) {fprintf(stdout,"cannot open %s! (critical)\n", buffer);exit(1);}
			while(2 == fscanf(g_annot,"%[^\t\n]%c", buffer,&sep)) {
				if (sep == '\t'){
					l = strlen(buffer);
					j = annot_to_annotID.find(buffer);
					buffer[l] = '\t';
					fscanf(g_annot,"%[^\n]\n",buffer+l+1);
					if (j != 0xFFFFFFFF) annot_descr[annot_to_annotID.deref(j)] = cloneString(buffer);
				}
			}
			
			for(i=0;i<annot_to_annotID.heap.size();i++){
				if (annot_descr[annot_to_annotID.heap[i].first.d] == NULL){
					j = annot_to_annotID.heap[i].first.k.length();
					memcpy(buffer,annot_to_annotID.heap[i].first.k.c_str(),j);
					buffer[j] = '\t';
					buffer[j+1] = '\0'; annot_descr[annot_to_annotID.heap[i].first.d] = cloneString(buffer);
				}
			}
			
		
			// finding collumn
			unsigned int col_name=0;
			col_name--;
			do{ col_name++;if (2 != fscanf(f_names,"%[^\t\n]%c",buffer, &sep)) {printf("could not find collunm %s in cdt file!\n", cdtcolname ? cdtcolname : "NAME"); exit(1);}; 
			}while( strcmp(buffer, cdtcolname ? cdtcolname : "NAME") != 0);
			if (sep != '\n') fscanf(f_names,"%*[^\n]\n%*[^\n]\n"); 
			else fscanf(f_names,"%*[^\n]\n");
			
			printf("Loading genelist (%ith col for name)...\n", col_name);fflush(stdout);
			l=0;
			while(1 == fscanf(f_names,"%[^\t\n]\t", buffer)){
				
				for(i=1;i<col_name;i++) fscanf(f_names,"%*[^\t\n]%*c");
				fscanf(f_names,"%[^\t\n]%c", buffer+1024,&sep); if (sep != '\n') fscanf(f_names,"%*[^\n]\n");
				i = name_to_internal.find(buffer+1024);
				l++;
				if (i != 0xFFFFFFFF){
					k = name_to_internal.deref(i);
					cdtname_to_internal[string(buffer)] = k;
					node_names[k] = cloneString(buffer+1024);
					for(j=0;j < annotIDs[k].size();j++) annotID_counts[annotIDs[k][j]]++; // add count!
				} else {printf("no annotation for %s\n", buffer+1024); cdtname_to_internal[string(buffer)] = annotIDs.size()-1;}
			}
			fclose(f_names);
			
			printf("\t...(DONE)! %i protein in cdt file\nRemoving Annotations with nb occurence outside input range...", cdtname_to_internal.getSize());fflush(stdout);
			
			Vector<unsigned int> annotID_counts_dest;
			annotID_counts_dest.setSize(annotID_counts.size()); ExOp::toZero(annotID_counts_dest);
			
			for(j=0;j<annotIDs.size();j++){
				for(k=0;k<annotIDs[j].size();k++){
					
					if ((annotID_counts[annotIDs[j][k]] < nb_annot_range[0])||(annotID_counts[annotIDs[j][k]] > nb_annot_range[1])) annotIDs[j].pop_swap(k--);
				}
				annotIDs[j].sort();
			}
			
			unsigned int nb_genes = l;
			
		//	for(j=0;j<annot_descr.size();j++) printf("%i\t%i\t%s\n" ,annotID_counts[j], annotID_counts_dest[j], annot_descr[j]);
			
			
			
			
			unsigned int jjj;
		
			
			printf("(DONE)\n");
			
			myHashmap<string, unsigned int> gtrname_to_posis;
			Forest< Vector<  KeyElem< unsigned int,pair<unsigned int, double> > >, 3> super_tree;
			
			
			
			Vector<  KeyElem< unsigned int,pair<unsigned int, double> > > *targ_tree[2];
			
			super_tree.setSize(nb_genes*2-1);
			k = nb_genes*2-2;
			l = nb_genes-2;
			
			Vector< pair<unsigned int, unsigned int> > merge_order; merge_order.setSize(nb_genes-1);
			Vector< unsigned int> internal_to_treeorder; internal_to_treeorder.setSize(annotIDs.size());

			
			while(3 == fscanf(f_trees,"%[^\t]\t%[^\t]\t%[^\t]\t%*[^\n]\n", buffer,buffer+64,buffer+128)){
				
				
				i = gtrname_to_posis.find(buffer+64);
				super_tree.clear(l);
				if (i == 0xFFFFFFFF){
					// must be a gene
					j = cdtname_to_internal.find(buffer+64);
					internal_to_treeorder[j] = k;
					
					super_tree.clear(k);
					if (j != 0xFFFFFFFF) {
						jjj = cdtname_to_internal.deref(j);
						super_tree[k].second.setSize(1+annotIDs[jjj].size());
						super_tree[k].second[0].k = 1;
						super_tree[k].second[0].d.first = jjj; 
						for(i=0;i<annotIDs[jjj].size();i++) {super_tree[k].second[i+1].k = annotIDs[jjj][i];super_tree[k].second[i+1].d.first =1; annotID_counts_dest[annotIDs[jjj][i]]++; }
					}else {printf("gene %s missing in annotation table!\n", buffer+64);	super_tree[k].second.setSize(1); super_tree[k].second[0].k = 1; super_tree[k].second[0].d.first = annotIDs.size()-1; }
					super_tree.makeLeftof(k, l);k--;
				}else super_tree.makeLeftof(gtrname_to_posis.deref(i), l);
				
				i = gtrname_to_posis.find(buffer+128);
				if (i == 0xFFFFFFFF){
					// must be a gene
					j = cdtname_to_internal.find(buffer+128);
					internal_to_treeorder[j] = k;
					
					super_tree.clear(k);
					if (j != 0xFFFFFFFF) {
						jjj = cdtname_to_internal.deref(j);
						super_tree[k].second.setSize(1+annotIDs[jjj].size());
						super_tree[k].second[0].k = 1;
						super_tree[k].second[0].d.first = jjj; 
						for(i=0;i<annotIDs[jjj].size();i++) {super_tree[k].second[i+1].k = annotIDs[jjj][i];super_tree[k].second[i+1].d.first =1; annotID_counts_dest[annotIDs[jjj][i]]++;}
					}else {printf("gene %s missing in annotation table!\n", buffer+128);	super_tree[k].second.setSize(1); super_tree[k].second[0].k = 1; super_tree[k].second[0].d.first = annotIDs.size()-1;}
					super_tree.makeRightof(k, l);k--;
				} else super_tree.makeRightof(gtrname_to_posis.deref(i), l);
				
			//	printf("%i<-%i,%i\n",l,super_tree.getLeft(l),super_tree.getRight(l));
				
				targ_tree[0] = &(super_tree[super_tree.getLeft(l)].second);
				targ_tree[1] = &(super_tree[super_tree.getRight(l)].second);

				
				merge_order[l] = pair<unsigned int, unsigned int>(super_tree.getLeft(l),super_tree.getRight(l) );
				
				// merging time!
				super_tree[l].second.push_back(KeyElem< unsigned int, pair<unsigned int, double> >( (*targ_tree[0])[0].k + (*targ_tree[1])[0].k ,pair<unsigned int, double>(0 , 0.0f)));
				
				i=j=1;
				
				
				while((targ_tree[0]->size() > i)&&(targ_tree[1]->size() > j)){
					
					
					if ((*targ_tree[0])[i].k == (*targ_tree[1])[j].k){ //printf("ismerge\n");
						dval = Madstructs::lncumhypergeo((*targ_tree[0])[i].d.first + (*targ_tree[1])[j].d.first,annotID_counts[(*targ_tree[0])[i].k],super_tree[l].second[0].k,nb_genes) / -log(10);
						super_tree[l].second.push_back(KeyElem< unsigned int, pair<unsigned int, double> >( (*targ_tree[0])[i].k ,pair<unsigned int, double>((*targ_tree[0])[i].d.first + (*targ_tree[1])[j].d.first , dval ) ) );
						
						if (best_pvalue[(*targ_tree[0])[i].k].first < dval) best_pvalue[(*targ_tree[0])[i].k] = pair<double,pair<unsigned int,unsigned int > >(dval,pair<unsigned int,unsigned int >(l,super_tree[l].second.size()-1));
						i++;j++;
					}else if ((*targ_tree[0])[i].k < (*targ_tree[1])[j].k){//printf("ispush\n");
						super_tree[l].second.push_back(KeyElem< unsigned int, pair<unsigned int, double> >( (*targ_tree[0])[i].k ,pair<unsigned int, double>((*targ_tree[0])[i].d.first , 0.0f ) ) );
						i++;
					}else{// printf("ispush\n");
						super_tree[l].second.push_back(KeyElem< unsigned int, pair< unsigned int, double> >( (*targ_tree[1])[j].k ,pair<unsigned int, double>((*targ_tree[1])[j].d.first, 0.0f ) ) );
						j++;
					}
					
				}
				
				
				if (targ_tree[0]->size() > i){
					for(;i<targ_tree[0]->size();i++){
						super_tree[l].second.push_back(KeyElem< unsigned int, pair<unsigned int, double> >( (*targ_tree[0])[i].k ,pair<unsigned int, double>((*targ_tree[0])[i].d.first , 0.0f ) ) );
					}
				}else if (targ_tree[1]->size() > j){
					for(;j<targ_tree[1]->size();j++){
						super_tree[l].second.push_back(KeyElem< unsigned int, pair< unsigned int, double> >( (*targ_tree[1])[j].k ,pair<unsigned int, double>((*targ_tree[1])[j].d.first , 0.0f ) ) );
					}
				}
				
				//		printf("%i\t%i\t%i\n", super_tree[l].second.size(), targ_tree[0]->size(), targ_tree[1]->size());
				
				
				gtrname_to_posis[cloneString(buffer)] = l--;
				
				
			}
			
			if (l != (0-1)) exit(1);
			
			
			//				for(i=0;i<super_tree[0].second.size();i++) printf("%i\t%i\n",super_tree[0].second[i].d.first,annotID_counts[super_tree[0].second[i].k]);

			

			for(dval = 0.0f,i=0;i<best_pvalue.size();i++) if ((annotID_counts[i] >= nb_annot_range[0])&&(annotID_counts[i] <= nb_annot_range[1])) dval += best_pvalue[i].first;
			
			for(i=0,j=0;i<best_pvalue.size();i++) if ((annotID_counts[i] >= nb_annot_range[0])&&(annotID_counts[i] <= nb_annot_range[1])) j++;
		
			printf("sum of log-pvalue: %f (for %i annotations)\n",dval, j);
			da_super_sum =dval; 
			
			stack<unsigned int> stck;
			if (threshold != 0.0f){
				FILE* ggg = fopen(threshold_path,"w+");
				for(i=0;i<best_pvalue.getSize();i++) if (best_pvalue[i].first > threshold) {
					dval = (((double)(super_tree[best_pvalue[i].second.first].second[best_pvalue[i].second.second].d.first * nb_genes)) / annotID_counts[super_tree[best_pvalue[i].second.first].second[best_pvalue[i].second.second].k]) / super_tree[best_pvalue[i].second.first].second[0].k;
					fprintf(ggg,"%f\t%f\t%i\t%i\t%i\t%s\t",best_pvalue[i].first, dval, super_tree[best_pvalue[i].second.first].second[best_pvalue[i].second.second].d.first,annotID_counts[super_tree[best_pvalue[i].second.first].second[best_pvalue[i].second.second].k],super_tree[best_pvalue[i].second.first].second[0].k,
						   annot_descr[super_tree[best_pvalue[i].second.first].second[best_pvalue[i].second.second].k]);
					
					
					stck.push(best_pvalue[i].second.first);
					while(!stck.empty()){
						k = stck.top(); stck.pop(); 
						if (super_tree.hasLeft(k)){
							stck.push(super_tree.getLeft(k));
							stck.push(super_tree.getRight(k));
						}else{
							for(l=1; l<super_tree[k].second.size() ;l++) if (super_tree[k].second[l].k == super_tree[best_pvalue[i].second.first].second[best_pvalue[i].second.second].k) break;
							if (l<super_tree[k].second.size()) fprintf(ggg,"%s ", node_names[super_tree[k].second[0].d.first]);
						}
					}
					
					fprintf(ggg,"\n");
				}
				fclose(ggg);
			}
			
			Vector< Vector<unsigned int> > allowed;
			Vector<unsigned int> group;
			Vector< unsigned int> permute; 
			
			myHashmap<string, unsigned int> daHclasses_names_back;
			vector< string > daHclasses_names;
			
			WeightElem<double, 4> back_dist; ExOp::toZero(back_dist);
			WeightElem<double, 4> back_edist; ExOp::toZero(back_edist);
			if (nb_permutations > 0) {
				permute.setSize(nb_genes); for(i=0;i<nb_genes;i++) permute[i] = i + nb_genes-1;
				if (permute_file){
					printf("Loading Permutation file...\n"); fflush(stdout);
					FILE* ggg = fopen(permute_file,"r+"); if (ggg == NULL) {fprintf(stdout,"cannot open permutation file %s! (critical)\n", permute_file);exit(1);} 
					group.setSize(name_to_internal.heap.size());
					
					j=0;
					while(2 == fscanf(ggg,"%[^\t\n]\t%[^\t\n]\n", buffer,buffer+1024)){
						l = daHclasses_names_back.find(string(buffer+1024));
						if (l == 0xFFFFFFFF) {l = daHclasses_names.size(); daHclasses_names_back[string(buffer+1024)] = l; daHclasses_names.push_back(string(buffer+1024));}
						else l = daHclasses_names_back.deref(l);
						k = name_to_internal.find(buffer);
						if (k != 0xFFFFFFFF) group[name_to_internal.deref(k)] = l;
						if (j <= l) j = l+1;
					}
					fclose(ggg);
					allowed.setSize(j); 
					for(i=0;i<nb_genes;i++){
						l = cdtname_to_internal.heap[i].first.d;
						allowed[group[l]].push_back(internal_to_treeorder[i]-nb_genes +1);
					}
					printf("...(DONE)\n"); fflush(stdout);
					
				}
				unsigned int empirical;
			for(empirical=0 ;nb_permutations > 0; nb_permutations--){

				if (permute_file){
					for(i=0;i<allowed.size();i++){
						for(j=allowed[i].size()-1;j > 0 ;j--){
							k = rand() % (j +1);
							if (k != j) {
								l = permute[allowed[i][j]];
								permute[allowed[i][j]] = permute[allowed[i][k]];
								permute[allowed[i][k]] = l;
								
								l = allowed[i][j];
								allowed[i][j] = allowed[i][k];
								allowed[i][k] = l;
							}
						}
					}
				}else{
					for(i=0;i<nb_genes-1;i++){
						k = rand() % (nb_genes - i);
						if (k != 0){
							j = permute[i];
							permute[i] = permute[i+k];
							permute[i+k] =j;
						}
					}
				}
				
				ExOp::toZero(best_pvalue);
				
				for(l = nb_genes-2; l<nb_genes;l--){
						
					if (merge_order[l].first > nb_genes-2){ targ_tree[0] = &(super_tree[ permute[merge_order[l].first -nb_genes +1] ].second);
					} else targ_tree[0] =&(super_tree[merge_order[l].first].second);
					
					if (merge_order[l].second > nb_genes-2){ targ_tree[1] = &(super_tree[ permute[merge_order[l].second -nb_genes +1  ] ].second);
					} else targ_tree[1] = &(super_tree[merge_order[l].second].second);
					
					super_tree[l].second.clear();
					super_tree[l].second.push_back(KeyElem< unsigned int, pair<unsigned int, double> >( (*targ_tree[0])[0].k + (*targ_tree[1])[0].k ,pair<unsigned int, double>(0 , 0.0f)));

			
					i=j=1;
					
					
					while((targ_tree[0]->size() > i)&&(targ_tree[1]->size() > j)){
						if ((*targ_tree[0])[i].k == (*targ_tree[1])[j].k){ //printf("ismerge\n");
							dval = Madstructs::lncumhypergeo((*targ_tree[0])[i].d.first + (*targ_tree[1])[j].d.first,annotID_counts[(*targ_tree[0])[i].k],super_tree[l].second[0].k,nb_genes) / -log(10);
							super_tree[l].second.push_back(KeyElem< unsigned int, pair<unsigned int, double> >( (*targ_tree[0])[i].k ,pair<unsigned int, double>((*targ_tree[0])[i].d.first + (*targ_tree[1])[j].d.first , dval ) ) );
							
							if (best_pvalue[(*targ_tree[0])[i].k].first < dval) best_pvalue[(*targ_tree[0])[i].k] = pair<double,pair<unsigned int,unsigned int > >(dval,pair<unsigned int,unsigned int >(l,super_tree[l].second.size()-1));
							i++;j++;
						}else if ((*targ_tree[0])[i].k < (*targ_tree[1])[j].k){//printf("ispush\n");
							super_tree[l].second.push_back(KeyElem< unsigned int, pair<unsigned int, double> >( (*targ_tree[0])[i].k ,pair<unsigned int, double>((*targ_tree[0])[i].d.first , 0.0f ) ) );
							i++;
						}else{// printf("ispush\n");
							super_tree[l].second.push_back(KeyElem< unsigned int, pair< unsigned int, double> >( (*targ_tree[1])[j].k ,pair<unsigned int, double>((*targ_tree[1])[j].d.first, 0.0f ) ) );
							j++;
						}
						
					}
					
					
					if (targ_tree[0]->size() > i){
						for(;i<targ_tree[0]->size();i++){
							super_tree[l].second.push_back(KeyElem< unsigned int, pair<unsigned int, double> >( (*targ_tree[0])[i].k ,pair<unsigned int, double>((*targ_tree[0])[i].d.first , 0.0f ) ) );
						}
					}else if (targ_tree[1]->size() > j){
						for(;j<targ_tree[1]->size();j++){
							super_tree[l].second.push_back(KeyElem< unsigned int, pair< unsigned int, double> >( (*targ_tree[1])[j].k ,pair<unsigned int, double>((*targ_tree[1])[j].d.first , 0.0f ) ) );
						}
					}
					
					
					
					
					
				}
				
				for(dval = 0.0f,i=0;i<best_pvalue.size();i++) if ((annotID_counts[i] >= nb_annot_range[0])&&(annotID_counts[i] <= nb_annot_range[1])) dval += best_pvalue[i].first;
				if (dval > da_super_sum) empirical++;
				back_dist += WeightElem<double,4>(dval);
				back_edist += WeightElem<double,4>(exp(dval / da_super_sum));
				
			}
				printf("mean=%f\tvar=%f\t4thmom=%f\tnb_samples=%i\n", back_dist.getMean(), back_dist.getSecondMomment(), back_dist.getFouthMomment(), (unsigned int)back_dist.w[0] );
				//pix[0] = (log(back_dist.getFouthMomment()) - 4.0f * log(da_super_sum - back_dist.getMean())) / log(10);
				
				pix[0] = (log(back_edist.getMean()) + da_super_sum) / log(10);
				pix[1] = (log(back_dist.getVar()) - 2.0f * log(da_super_sum - back_dist.getMean()))/ log(10) ;
				if (pix[1]>= 0.0f) pix[1] =0.0f;
				if (pix[0]>= 0.0f) pix[0] =0.0f;
				printf("4th momment Chebyshev Bound(Chebyshev Bound) = %f\t%f\t%i\n", pix[0], pix[1], empirical);

				
			//	for(j=0;j<annot_descr.size();j++) if (annotID_counts[j] > 1) printf("%i\t%i\t%s\n" ,annotID_counts[j], annotID_counts_dest[j], annot_descr[j]);
			}
			
			
			
			
			
			
		}				
		
		return 0;
		
	}
	void Taskscope<TASK_CDTGTR_PVALUES>::help(){
		printf("Analyses P-values in a hierachical cluster.\n");
		printf("Default Arguments 0-1:\n");
		printf("(in) filename (needs both a .cdt and .gtr file)\n");
		printf("\n");
		printf("Optionnal Default arguments(3):\n");
		printf("-p: change path variables (auto-on if nbargs < 1)! If corrupted, del PMGTRpvalues_scope.scp\n");
		printf("-F (file): filter annotations to the ones in the file (1 collunm)\n");
		printf("-f (int min)-(int max): filter annotations based on the number of instances (default: 2 and more)\n");
		
		printf("-P (int nb_perm) [FILE* permute constraint]: permute the genes for significance testing\n");
		printf("-t (file) (float): display all pvalu abovelog10 pvalue threshold\n");
		printf("-N (string): collunm name in cdt file identifying protein for annotation matching (default = \"NAME\")\n");
		//					printf("-D(int size1, int size2): forces the output image size to match\n");
		//	printf("-O: Orient for tall image, or puts mother-bud verticaly\n");
		//					printf("-g: image has 2 frames per plates, green then red.\n");
		//					printf("-S[RGB]: scale colors independantly to maximize contrast.\n");
		//					printf("-F[RGB]: scale factor, or maximum intensity allowed for automatic scaling (with flag -S).\n");
		//					printf("-A: use Area instead of ID for output file name.\n");
		printf("\t-h : Help \n");
	}
	
	
	Taskscope<TASK_RANK_PVALUES>::Taskscope():threshold(0.0f),reverse(false),show(false),filter_file(NULL),permute_file(NULL),torun(0),inter(false),w_name(false),gray(false),doopen(false),dodel(false),cdtcolname(NULL),nb_permutations(0) {
		nb_annot_range[0] =2; nb_annot_range[1] = 0xFFFFFFFF;
	}
	void Taskscope<TASK_RANK_PVALUES>::nbaddtoken(char const * const token, int& min, int& max){
		switch(*token){
			case '\0': min =0; max =1; break;
			case 'p': min =1; max=1; break;
			case 't': min =2; break;
			case 'A': min =0; break;
			case 'g': min =0; break;
			case 'r': min =0; break;
			case 'd': min =0; break;
			case 'F': min =1; break;
			case 'f': min =2; break;
			case 'N': min =1; break;
			case 'P': min =1; max=2; break;
		}
	}
	
	void Taskscope<TASK_RANK_PVALUES>::store(char* const * token, int nbtoken){ 
		switch(token[0][1]){
			case 's': show = true; break;
			case 't': threshold_path = token[1] ; threshold = atof(token[2]); break;
			case 'p': inter = true; break;
			case 'g': gray = true; break;
			case 'A': w_name = true; break;
			case 'r': reverse = true; break;
			case 'F': filter_file = token[1];break;
			case 'f': nb_annot_range[0] = atoi(token[1]); nb_annot_range[1] = atoi(token[2]); break;
			case 'P': nb_permutations =atoi(token[1]); if (nbtoken > 1) permute_file = token[2]; break;
			case 'N': cdtcolname = token[1];break;
		}
	}
	int Taskscope<TASK_RANK_PVALUES>::defstore(char* const * token, int nbtoken){ 
		unsigned int i,j,k,l;
		char buffer[65536];
		double pix[32];
		
		SerialStore<unsigned int> rstd("PMGTRpvalues_scope.scp");
		Vector<char> str;
		
		if (inter){
			printf("enter the annotation table path:");
			if (rstd.has(0)){
				rstd.load(0, str);
				printf(" currently =%s\n", &(str[0]));
			} else printf("\n");
			scanf(" %[^\n]", buffer);
			l = strlen(buffer);
			if ((l)&&(strcmp(buffer," ") != 0)){
				str.setSize(l+1);
				memcpy(&(str[0]),buffer, l+1);
				rstd.save(0, str);
			}
			printf("enter the annotation description file:");
			if (rstd.has(1)){
				rstd.load(1, str);
				printf(" currently =%s\n", &(str[0]));
			} else printf("\n");
			scanf(" %[^\n]", buffer);
			l = strlen(buffer);
			if ((l)&&(strcmp(buffer," ") != 0)){
				str.setSize(l+1);
				memcpy(&(str[0]),buffer, l+1);
				rstd.save(1, str);
			}
		}else{
			
			

			
			FILE* f_names = fopen(token[0],"r+"); if (f_names == NULL) {fprintf(stdout,"cannot open %s! (critical)\n", token[0]);exit(1);}
			
			Vector< pair<double, pair<unsigned int,unsigned int > > > best_pvalue; 
			Vector<unsigned int> annotID_counts;
			myHashmap<string, unsigned int> annot_to_annotID;
			
			double dval;
			double da_super_sum;
			Vector< Vector<unsigned int> > annotIDs;
			myHashmap<string, unsigned int > name_to_internal;
			
			char* name = NULL;
			Tuple<char,32> ID;
			char sep;
			Vector<unsigned int> annot_input;
			j=0; // annotID count
			
			myHashmap<string,char> filter;
			
			FILE* g_annot;
			if (filter_file){
				printf("loading filter file..."); fflush(stdout);
				g_annot = fopen(filter_file,"r+"); if (g_annot == NULL) {fprintf(stdout,"cannot open %s! (critical)\n", filter_file);exit(1);}
				while(1 == fscanf(g_annot,"%[^\t\n]%*c", buffer)) filter[cloneString(buffer)] =0;
				fclose(g_annot);
				printf("(DONE)\n"); fflush(stdout);
			}
			
			rstd.load(0, str);
			j = strlen(&(str[0]));
			memcpy(buffer,&(str[0]), j);
			buffer[j] = '\0';
			g_annot = fopen(buffer,"r+"); if (g_annot == NULL) {fprintf(stdout,"cannot open %s! (critical)\n", buffer);exit(1);}
			
			// LOADS ALL ANNOTATIONS! 
			printf("loading annotations...");fflush(stdout);
			k =0;
			while(2 == fscanf(g_annot,"%[^\t\n]%c", buffer,&sep)){
				if (name == NULL){
					name = cloneString(buffer);
				}else{
					if ((filter_file == NULL)||(filter.find(buffer) != 0xFFFFFFFF)){
						i = annot_to_annotID.find(buffer);
						if (i == 0xFFFFFFFF){
							
							annot_to_annotID[cloneString(buffer)] = j++;
							i = annot_to_annotID.find(buffer);
							if (i == 0xFFFFFFFF) {printf("NOSE!\n");exit(0);}
						}
						annot_input.push_back(annot_to_annotID.deref(i));
					}
				}
				if (sep == '\n'){
					name_to_internal[string(name)] = k;
					annot_input.sort();
					annotIDs.push_back(annot_input);
					k++;
					delete[](name);
					name = NULL;
					annot_input.clear();
				}
			}
			annotIDs.push_back(annot_input); // empty annotations
			
			//	for(l=0;l< annotIDs.size();l++) printf("size=%i\n", annotIDs[l].size() );
			
			fclose(g_annot);
			printf("(DONE) (%i annotation total)\n", j); fflush(stdout);
			Vector<char*> annot_descr; 
			
			Vector<char*> node_names; node_names.setSize(k);ExOp::toZero(node_names);
			
			annotID_counts.setSize(j);ExOp::toZero(annotID_counts);
			best_pvalue.setSize(j);ExOp::toZero(best_pvalue);
			annot_descr.setSize(j);ExOp::toZero(annot_descr);
			
			
			rstd.load(1, str);
			j = strlen(&(str[0]));
			memcpy(buffer,&(str[0]), j);
			buffer[j] = '\0';
			g_annot = fopen(buffer,"r+"); if (g_annot == NULL) {fprintf(stdout,"cannot open %s! (critical)\n", buffer);exit(1);}
			while(2 == fscanf(g_annot,"%[^\t\n]%c", buffer,&sep)) {
				if (sep == '\t'){
					l = strlen(buffer);
					j = annot_to_annotID.find(buffer);
					buffer[l] = '\t';
					fscanf(g_annot,"%[^\n]\n",buffer+l+1);
					if (j != 0xFFFFFFFF) annot_descr[annot_to_annotID.deref(j)] = cloneString(buffer);
				}
			}
			
			for(i=0;i<annot_to_annotID.heap.size();i++){
				if (annot_descr[annot_to_annotID.heap[i].first.d] == NULL){
					j = annot_to_annotID.heap[i].first.k.length();
					memcpy(buffer,annot_to_annotID.heap[i].first.k.c_str(),j);
					buffer[j] = '\t';
					buffer[j+1] = '\0'; annot_descr[annot_to_annotID.heap[i].first.d] = cloneString(buffer);
				}
			}
			
			
			// finding collumn
			Vector<unsigned int> internal_order;

			l=0;
			while(1 == fscanf(f_names,"%[^\t\n]%*c", buffer)){
				i = name_to_internal.find(buffer);
				internal_order.push_back(i);
				l++;
				if (i != 0xFFFFFFFF){
					k = name_to_internal.deref(i);
					node_names[k] = cloneString(buffer);
					for(j=0;j < annotIDs[k].size();j++) annotID_counts[annotIDs[k][j]]++; // add count!
					
				} else {printf("no annotation for %s\n", buffer);}
			}
			fclose(f_names);
			
			if (reverse) internal_order.reverse();
			
			printf("\t...(DONE)! %i protein in rank file\nRemoving Annotations with nb occurence outside input range...", internal_order.getSize());fflush(stdout);
			

		
			for(j=0;j<annotIDs.size();j++){
				for(k=0;k<annotIDs[j].size();k++){
					if ((annotID_counts[annotIDs[j][k]] < nb_annot_range[0])||(annotID_counts[annotIDs[j][k]] > nb_annot_range[1])) annotIDs[j].pop_swap(k--);
				}
				annotIDs[j].sort();
			}
			
			unsigned int nb_genes = l;
			printf("(DONE)\n");
			
			Vector< KeyElem< unsigned int,pair<unsigned int, double> > >* super_tree;
			
			super_tree = new Vector< KeyElem< unsigned int,pair<unsigned int, double> > >[nb_genes+1];
			
			Vector< unsigned int> internal_to_treeorder; internal_to_treeorder.setSize(annotIDs.size());
			
			k=0;

			if (internal_order[k] == 0xFFFFFFFF)	super_tree[0].setSize(1);
			else{
			super_tree[0].setSize(1+annotIDs[internal_order[k]].size());
			for(i=0;i<annotIDs[internal_order[k]].size();i++) {super_tree[0][i+1].k = annotIDs[internal_order[k]][i];super_tree[0][i+1].d.first =1;}
			}
			super_tree[0][0].k = 1;
			super_tree[0][0].d.first = internal_order[k];
			
			for(k++;k < nb_genes;k++){
				if (internal_order[k] == 0xFFFFFFFF) {
					super_tree[nb_genes].setSize(1);
				}else{
				super_tree[nb_genes].setSize(1+annotIDs[internal_order[k]].size());
				for(i=0;i<annotIDs[internal_order[k]].size();i++) {super_tree[nb_genes][i+1].k = annotIDs[internal_order[k]][i];super_tree[nb_genes][i+1].d.first =1; }
				}
				super_tree[nb_genes][0].k = 1;
				super_tree[nb_genes][0].d.first = internal_order[k];
				
				i=j=1;
				super_tree[k].push_back(KeyElem< unsigned int, pair<unsigned int, double> >( super_tree[k-1][0].k + super_tree[nb_genes][0].k ,pair<unsigned int, double>(0 , 0.0f)));
				while((super_tree[k-1].size() > i)&&(super_tree[nb_genes].size() > j)){
					if (super_tree[k-1][i].k == super_tree[nb_genes][j].k){ //printf("ismerge\n");
						dval = Madstructs::lncumhypergeo(super_tree[k-1][i].d.first + super_tree[nb_genes][j].d.first,annotID_counts[super_tree[k-1][i].k],super_tree[k][0].k,nb_genes) / -log(10);
						
						if (best_pvalue[super_tree[k-1][i].k].first < dval){
							best_pvalue[super_tree[k-1][i].k] = pair<double,pair<unsigned int,unsigned int > >(dval,pair<unsigned int,unsigned int >(k,super_tree[k].size()));
							super_tree[k].push_back(KeyElem< unsigned int, pair<unsigned int, double> >( super_tree[k-1][i].k ,pair<unsigned int, double>(super_tree[k-1][i].d.first + super_tree[nb_genes][j].d.first , dval ) ) );
						}else{
						//	if (best_pvalue[super_tree[k-1][i].k].first < (Madstructs::lncumhypergeo(annotID_counts[super_tree[k-1][i].k],annotID_counts[super_tree[k-1][i].k],super_tree[k][0].k + annotID_counts[super_tree[k-1][i].k] - super_tree[k-1][i].d.first - super_tree[nb_genes][j].d.first,nb_genes) / -log(10))) // there is still hope!
								super_tree[k].push_back(KeyElem< unsigned int, pair<unsigned int, double> >( super_tree[k-1][i].k ,pair<unsigned int, double>(super_tree[k-1][i].d.first + super_tree[nb_genes][j].d.first , dval ) ) );
						}
						
						i++;j++;
					}else if (super_tree[k-1][i].k < super_tree[nb_genes][j].k){//printf("ispush\n");
						super_tree[k].push_back(KeyElem< unsigned int, pair<unsigned int, double> >( super_tree[k-1][i].k ,pair<unsigned int, double>(super_tree[k-1][i].d.first , 0.0f ) ) );
						i++;
					}else{
						super_tree[k].push_back(KeyElem< unsigned int, pair< unsigned int, double> >( super_tree[nb_genes][j].k ,pair<unsigned int, double>(super_tree[nb_genes][j].d.first, 0.0f ) ) );
						j++;
					}
				}
				if (super_tree[k-1].size() > i) for(;i<super_tree[k-1].size();i++) super_tree[k].push_back(KeyElem< unsigned int, pair<unsigned int, double> >( super_tree[k-1][i].k ,pair<unsigned int, double>(super_tree[k-1][i].d.first , 0.0f ) ) );
				else if (super_tree[nb_genes].size() > j) for(;j<super_tree[nb_genes].size();j++) super_tree[k].push_back(KeyElem< unsigned int, pair< unsigned int, double> >( super_tree[nb_genes][j].k ,pair<unsigned int, double>(super_tree[nb_genes][j].d.first , 0.0f ) ) );
			}
			
			
			for(dval = 0.0f,i=0;i<best_pvalue.size();i++) if ((annotID_counts[i] >= nb_annot_range[0])&&(annotID_counts[i] <= nb_annot_range[1])) dval += best_pvalue[i].first;
			
			for(i=0,j=0;i<best_pvalue.size();i++) if ((annotID_counts[i] >= nb_annot_range[0])&&(annotID_counts[i] <= nb_annot_range[1])) j++;
			
			printf("sum of log-pvalue: %f (for %i annotations)\n",dval, j);
			da_super_sum =dval; 
			
			stack<unsigned int> stck;
			
			if (threshold != 0.0f){
				FILE* ggg = fopen(threshold_path,"w+");
				fprintf(ggg,"Log10Pvalue\tFoldEnrichment\tnbHits\tnbAnnotations\tClusterSize\tAnnotation\tDescription\n");
				for(i=0;i<best_pvalue.getSize();i++) if (best_pvalue[i].first > threshold) {
					dval = (((double)(super_tree[best_pvalue[i].second.first][best_pvalue[i].second.second].d.first * nb_genes)) / annotID_counts[super_tree[best_pvalue[i].second.first][best_pvalue[i].second.second].k]) / super_tree[best_pvalue[i].second.first][0].k;
					if (super_tree[best_pvalue[i].second.first][best_pvalue[i].second.second].k != i) printf("warning! %i became %i\n", super_tree[best_pvalue[i].second.first][best_pvalue[i].second.second].k, i);
					fprintf(ggg,"%f\t%f\t%i\t%i\t%i\t%s",best_pvalue[i].first, dval,
							super_tree[best_pvalue[i].second.first][best_pvalue[i].second.second].d.first,
							annotID_counts[i],
							super_tree[best_pvalue[i].second.first][0].k, // group size
							annot_descr[super_tree[best_pvalue[i].second.first][best_pvalue[i].second.second].k]);
					fprintf(ggg,"\n");
				}
				fclose(ggg);
			}
			
			Vector< Vector<unsigned int> > allowed;
			Vector<unsigned int> group;
			Vector< unsigned int> permute; 
			
			myHashmap<string, unsigned int> daHclasses_names_back;
			vector< string > daHclasses_names;
			
			WeightElem<double, 4> back_dist; ExOp::toZero(back_dist);
			WeightElem<double, 4> back_edist; ExOp::toZero(back_edist);
			if (nb_permutations > 0) {
				permute.setSize(nb_genes); for(i=0;i<nb_genes;i++) permute[i] = i;
				if (permute_file){
					printf("Loading Permutation file...\n"); fflush(stdout);
					FILE* ggg = fopen(permute_file,"r+"); if (ggg == NULL) {fprintf(stdout,"cannot open permutation file %s! (critical)\n", permute_file);exit(1);} 
					group.setSize(name_to_internal.heap.size());
					
					j=0;
					while(2 == fscanf(ggg,"%[^\t\n]\t%[^\t\n]\n", buffer,buffer+1024)){
						l = daHclasses_names_back.find(string(buffer+1024));
						if (l == 0xFFFFFFFF) {l = daHclasses_names.size(); daHclasses_names_back[string(buffer+1024)] = l; daHclasses_names.push_back(string(buffer+1024));}
						else l = daHclasses_names_back.deref(l);
						k = name_to_internal.find(buffer);
						if (k != 0xFFFFFFFF) group[name_to_internal.deref(k)] = l;
						if (j <= l) j = l+1;
					}
					fclose(ggg);
					allowed.setSize(j); 
					for(i=0;i<nb_genes;i++){
						l = internal_order[i];
						if (l != 0xFFFFFFFF) allowed[group[l]].push_back(i);
					}
					printf("...(DONE)\n"); fflush(stdout);
					
				}
				unsigned int empirical;
				for(empirical=0 ;nb_permutations > 0; nb_permutations--){
					
					if (permute_file){
						for(i=0;i<allowed.size();i++){
							for(j=allowed[i].size()-1;j > 0 ;j--){
								k = rand() % (j +1);
								if (k != j) {
									l = permute[allowed[i][j]];
									permute[allowed[i][j]] = permute[allowed[i][k]];
									permute[allowed[i][k]] = l;
									
									l = allowed[i][j];
									allowed[i][j] = allowed[i][k];
									allowed[i][k] = l;
								}
							}
						}
					}else{
						for(i=0;i<nb_genes-1;i++){
							k = rand() % (nb_genes - i);
							if (k != 0){
								j = permute[i];
								permute[i] = permute[i+k];
								permute[i+k] =j;
							}
						}
					}
					
					ExOp::toZero(best_pvalue);
					
					k=0;
					if (internal_order[permute[k]] == 0xFFFFFFFF)	super_tree[0].setSize(1);
					else{
						super_tree[0].setSize(1+annotIDs[internal_order[permute[k]]].size());
						for(i=0;i<annotIDs[internal_order[permute[k]]].size();i++) {super_tree[0][i+1].k = annotIDs[internal_order[permute[k]]][i];super_tree[0][i+1].d.first =1;}
					}
					super_tree[0][0].k = 1;
					super_tree[0][0].d.first = internal_order[permute[k]];
					
					for(k++;k < nb_genes;k++){
						if (internal_order[permute[k]] == 0xFFFFFFFF) {
							super_tree[nb_genes].setSize(1);
						}else{
							super_tree[nb_genes].setSize(1+annotIDs[internal_order[permute[k]]].size());
							for(i=0;i<annotIDs[internal_order[permute[k]]].size();i++) {super_tree[nb_genes][i+1].k = annotIDs[internal_order[permute[k]]][i];super_tree[nb_genes][i+1].d.first =1; }
						}
						super_tree[nb_genes][0].k = 1;
						super_tree[nb_genes][0].d.first = internal_order[permute[k]];
						
						i=j=1;
						super_tree[k].clear();
						super_tree[k].push_back(KeyElem< unsigned int, pair<unsigned int, double> >( super_tree[k-1][0].k + super_tree[nb_genes][0].k ,pair<unsigned int, double>(0 , 0.0f)));
						while((super_tree[k-1].size() > i)&&(super_tree[nb_genes].size() > j)){
							if (super_tree[k-1][i].k == super_tree[nb_genes][j].k){ //printf("ismerge\n");
								dval = Madstructs::lncumhypergeo(super_tree[k-1][i].d.first + super_tree[nb_genes][j].d.first,annotID_counts[super_tree[k-1][i].k],super_tree[k][0].k,nb_genes) / -log(10);
								if (best_pvalue[super_tree[k-1][i].k].first < dval){ best_pvalue[super_tree[k-1][i].k] = pair<double,pair<unsigned int,unsigned int > >(dval,pair<unsigned int,unsigned int >(k,super_tree[k].size()));
									super_tree[k].push_back(KeyElem< unsigned int, pair<unsigned int, double> >( super_tree[k-1][i].k ,pair<unsigned int, double>(super_tree[k-1][i].d.first + super_tree[nb_genes][j].d.first , dval ) ) );
								}else{
									dval = Madstructs::lncumhypergeo(annotID_counts[super_tree[k-1][i].k],annotID_counts[super_tree[k-1][i].k],super_tree[k][0].k + annotID_counts[super_tree[k-1][i].k] - super_tree[k-1][i].d.first - super_tree[nb_genes][j].d.first,nb_genes) / -log(10);					
									if (best_pvalue[super_tree[k-1][i].k].first < dval) // there is still hope!
										super_tree[k].push_back(KeyElem< unsigned int, pair<unsigned int, double> >( super_tree[k-1][i].k ,pair<unsigned int, double>(super_tree[k-1][i].d.first + super_tree[nb_genes][j].d.first , dval ) ) );
								}
								i++;j++;
							}else if (super_tree[k-1][i].k < super_tree[nb_genes][j].k){//printf("ispush\n");
								super_tree[k].push_back(KeyElem< unsigned int, pair<unsigned int, double> >( super_tree[k-1][i].k ,pair<unsigned int, double>(super_tree[k-1][i].d.first , 0.0f ) ) );
								i++;
							}else{
								super_tree[k].push_back(KeyElem< unsigned int, pair< unsigned int, double> >( super_tree[nb_genes][j].k ,pair<unsigned int, double>(super_tree[nb_genes][j].d.first, 0.0f ) ) );
								j++;
							}
						}
						if (super_tree[k-1].size() > i) for(;i<super_tree[k-1].size();i++) super_tree[k].push_back(KeyElem< unsigned int, pair<unsigned int, double> >( super_tree[k-1][i].k ,pair<unsigned int, double>(super_tree[k-1][i].d.first , 0.0f ) ) );
						else if (super_tree[nb_genes].size() > j) for(;j<super_tree[nb_genes].size();j++) super_tree[k].push_back(KeyElem< unsigned int, pair< unsigned int, double> >( super_tree[nb_genes][j].k ,pair<unsigned int, double>(super_tree[nb_genes][j].d.first , 0.0f ) ) );
						
					}
					
					
					for(dval = 0.0f,i=0;i<best_pvalue.size();i++) if ((annotID_counts[i] >= nb_annot_range[0])&&(annotID_counts[i] <= nb_annot_range[1])) dval += best_pvalue[i].first;
					if (dval > da_super_sum) empirical++;
					back_dist += WeightElem<double,4>(dval);
					back_edist += WeightElem<double,4>(exp(dval / da_super_sum));
					
				}
				printf("mean=%f\tvar=%f\t4thmom=%f\tnb_samples=%i\n", back_dist.getMean(), back_dist.getSecondMomment(), back_dist.getFouthMomment(), (unsigned int)back_dist.w[0] );
				//pix[0] = (log(back_dist.getFouthMomment()) - 4.0f * log(da_super_sum - back_dist.getMean())) / log(10);
				
				pix[0] = (log(back_edist.getMean()) + da_super_sum) / log(10);
				pix[1] = (log(back_dist.getVar()) - 2.0f * log(da_super_sum - back_dist.getMean()))/ log(10) ;
				if (pix[1]>= 0.0f) pix[1] =0.0f;
				if (pix[0]>= 0.0f) pix[0] =0.0f;
				printf("4th momment Chebyshev Bound(Chebyshev Bound) = %f\t%f\t%i\n", pix[0], pix[1], empirical);
				
				
			}
			
		}				
		
		return 0;
		
	}
	void Taskscope<TASK_RANK_PVALUES>::help(){
		printf("Analyses P-values in a ordered protein list.\n");
		printf("Default Arguments 0-1:\n");
		printf("(in) filename\n");
		printf("\n");
		printf("Optionnal Default arguments:\n");
		printf("-p: change path variables (auto-on if nbargs < 1)! If corrupted, del PMGTRpvalues_scope.scp\n");
		printf("-F (file): filter annotations to the ones in the file (1 collunm)\n");
		printf("-f (int min)-(int max): filter annotations based on the number of instances (default: 2 and more)\n");
		printf("-r : reverses the ordering (find enrichment in the tail)\n");
		printf("-P (int nb_perm) [FILE* permute constraint]: permute the genes for significance testing\n");
		printf("-t (file) (float): display all pvalu abovelog10 pvalue threshold\n");
		printf("-N (string): collunm name in cdt file identifying protein for annotation matching (default = \"NAME\")\n");
		printf("\t-h : Help \n");
	}
	
	Taskscope<TASK_MULTI_PAIR_DISTANCE_PVAL>::Taskscope():threshold(0.0f),show(false),filter_file(NULL),permute_file(NULL),torun(0),inter(false),w_name(false),gray(false),scaleR(false),scaleG(false),scaleB(false),factR(0.0f),factG(0.0f),factB(0.0f),doopen(false),dodel(false),cdtcolname(NULL),nb_permutations(0) {
		nb_annot_range[0] =2; nb_annot_range[1] = 0xFFFFFFFF;
	}
	void Taskscope<TASK_MULTI_PAIR_DISTANCE_PVAL>::nbaddtoken(char const * const token, int& min, int& max){
		switch(*token){
			case '\0': min =0; max =1; break;
			case 'p': min =1; max =1; break;
			case 't': min =2; break;
			case 'A': min =0; break;
			case 'g': min =0; break;
			case 'S': min =0; break;
			case 'o': min =0; break;
			case 'd': min =0; break;
			case 'F': min =1; break;
			case 'f': min =2; break;
			case 'N': min =1; break;
			case 'P': min =1; max=2; break;
		}
	}
	
	void Taskscope<TASK_MULTI_PAIR_DISTANCE_PVAL>::store(char* const * token, int nbtoken){ 
		switch(token[0][1]){
			case 's': show = true; break;
			case 't': threshold_path = token[1] ; threshold = atof(token[2]); break;
			case 'p': inter = true; break;
			case 'g': gray = true; break;
			case 'A': w_name = true; break;
			case 'o': doopen = true; break;case 'd': dodel = true; break;
			case 'S': switch(token[0][2]){
				case 'R':scaleR = true; break;
				case 'G':scaleG = true; break;
				case 'B':scaleB = true; break;
				default: scaleR = true;scaleG = true;scaleB = true;
			}break;
			case 'F': filter_file = token[1];break;
			case 'f': nb_annot_range[0] = atoi(token[1]); nb_annot_range[1] = atoi(token[2]); break;
			case 'P': nb_permutations =atoi(token[1]); if (nbtoken > 1) permute_file = token[2]; break;
			case 'N': cdtcolname = token[1];break;
		}
	}
	int Taskscope<TASK_MULTI_PAIR_DISTANCE_PVAL>::defstore(char* const * token, int nbtoken){ 
		unsigned int i,j,k,l;
		char buffer[65536];
		double pix[32];
	
		SerialStore<unsigned int> rstd("PMMPDpvalues_scope.scp");
		Vector<char> str;
		
		if (inter){
			printf("enter the annotation table path:");
			if (rstd.has(0)){
				rstd.load(0, str);
				printf(" currently =%s\n", &(str[0]));
			} else printf("\n");
			scanf(" %[^\n]", buffer);
			l = strlen(buffer);
			if ((l)&&(strcmp(buffer," ") != 0)){
				str.setSize(l+1);
				memcpy(&(str[0]),buffer, l+1);
				rstd.save(0, str);
			}
			printf("enter the annotation description file:");
			if (rstd.has(1)){
				rstd.load(1, str);
				printf(" currently =%s\n", &(str[0]));
			} else printf("\n");
			scanf(" %[^\n]", buffer);
			l = strlen(buffer);
			if ((l)&&(strcmp(buffer," ") != 0)){
				str.setSize(l+1);
				memcpy(&(str[0]),buffer, l+1);
				rstd.save(1, str);
			}
		}else{
		
		
		}
		return 0;
	}
	
	
	void Taskscope<TASK_MULTI_PAIR_DISTANCE_PVAL>::help(){
		printf("Analyses P-values in a hierachical cluster.\n");
		printf("Default Arguments 0-1:\n");
		printf("(in) filename (needs both a .cdt and .gtr file)\n");
		printf("\n");
		printf("Optionnal Default arguments(3):\n");
		printf("-p: change path variables (auto-on if nbargs < 1)! If corrupted, del PMMakeDisplay_scope.scp\n");
		printf("-F (file): filter annotations to the ones in the file (1 collunm)\n");
		printf("-f (int min)-(int max): filter annotations based on the number of instances (default: 2 and more)\n");
		
		printf("-P (int nb_perm) [FILE* permute constraint]: permute the genes for significance testing\n");
		printf("-t (file) (float): display all pvalu abovelog10 pvalue threshold\n");
		printf("-N (string): collunm name in cdt file identifying protein for annotation matching (default = \"NAME\")\n");
		//					printf("-D(int size1, int size2): forces the output image size to match\n");
		//	printf("-O: Orient for tall image, or puts mother-bud verticaly\n");
		//					printf("-g: image has 2 frames per plates, green then red.\n");
		//					printf("-S[RGB]: scale colors independantly to maximize contrast.\n");
		//					printf("-F[RGB]: scale factor, or maximum intensity allowed for automatic scaling (with flag -S).\n");
		//					printf("-A: use Area instead of ID for output file name.\n");
		printf("OTHER\n");
		printf("-d: delete the last produced image.\n");
		printf("-o: open the image with ImageJ once done.\n");
		printf("\t-h : Help \n");
	}
	




