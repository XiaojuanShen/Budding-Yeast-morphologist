/*
 * Modeling.cpp
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

template<unsigned int size> double batta_cluster_metrix(const GaussElem<Tuple<double,size> >&l, const GaussElem<Tuple<double,size> >&r ){return l.bhattacharryya_dist(r);}
void eight_gaussreport(Vector<double> &fout,  Tuple< GaussElem< Tuple<double, 5> >, 8> &mb_intensity){
	GaussElem< Tuple<double, 3> > mb_residual[8];
	GaussElem< Tuple<double, 4> > mb_residual_par[8];
	double tmpdouble[6];
	Tuple<unsigned int, 3> get_res_coors;get_res_coors[0] =3;get_res_coors[1] =4;
	Tuple<unsigned int, 4> get_res_coors_par;get_res_coors_par[0] =3;get_res_coors_par[1] =4;
	unsigned int i,j,k;
	for(j=0;j<8;j++) fout.push_back(mb_intensity[j].w);
	for(j=0;j<8;j++) {tmpdouble[0] = (mb_intensity[j].w * mb_intensity[j].w) / mb_intensity[j].w2;  fout.push_back((ExOp::isValid(tmpdouble[0])) ? tmpdouble[0] : 0.0f);}
	k=0;for(j=0;j<8;j++) fout.push_back( mb_intensity[j].mean[get_res_coors[k]] / mb_intensity[j].w);
	//	fout.push_back( (mb_intensity[5].mean[get_res_coors[k]] + mb_intensity[6].mean[get_res_coors[k]]) / mb_intensity[7].w );
	//	fout.push_back( (+mb_intensity[0].mean[get_res_coors[k]] + mb_intensity[1].mean[get_res_coors[k]]+mb_intensity[2].mean[get_res_coors[k]] + mb_intensity[3].mean[get_res_coors[k]]+mb_intensity[4].mean[get_res_coors[k]] ) / mb_intensity[7].w );
	k=1;for(j=0;j<8;j++) fout.push_back( mb_intensity[j].mean[get_res_coors[k]] / mb_intensity[j].w);
	//	fout.push_back( (mb_intensity[5].mean[get_res_coors[k]] + mb_intensity[6].mean[get_res_coors[k]]) / mb_intensity[7].w );
	//	fout.push_back( (+mb_intensity[0].mean[get_res_coors[k]] + mb_intensity[1].mean[get_res_coors[k]]+mb_intensity[2].mean[get_res_coors[k]] + mb_intensity[3].mean[get_res_coors[k]]+mb_intensity[4].mean[get_res_coors[k]] ) / mb_intensity[7].w );
	k=0;for(j=0;j<8;j++) fout.push_back( mb_intensity[j].getVar()[get_res_coors[k]]);
	tmpdouble[0] = mb_intensity[7].getVar()[get_res_coors[k]];
	tmpdouble[1] = (mb_intensity[6].getVar()[get_res_coors[k]] * mb_intensity[6].w + mb_intensity[5].getVar()[get_res_coors[k]] * mb_intensity[5].w);
	tmpdouble[2] = (mb_intensity[0].getVar()[get_res_coors[k]] * mb_intensity[0].w + mb_intensity[1].getVar()[get_res_coors[k]] * mb_intensity[1].w+ mb_intensity[2].getVar()[get_res_coors[k]] * mb_intensity[2].w + mb_intensity[3].getVar()[get_res_coors[k]] * mb_intensity[3].w+mb_intensity[4].getVar()[get_res_coors[k]] * mb_intensity[4].w)  ;
	fout.push_back( tmpdouble[1]/ (mb_intensity[7].w * tmpdouble[0])); fout.push_back(tmpdouble[2]/ (mb_intensity[7].w * tmpdouble[0]));
	k=1;for(j=0;j<8;j++) fout.push_back( mb_intensity[j].getVar()[get_res_coors[k]]);
	tmpdouble[3] = mb_intensity[7].getVar()[get_res_coors[k]];
	tmpdouble[4] = (mb_intensity[6].getVar()[get_res_coors[k]] * mb_intensity[6].w + mb_intensity[5].getVar()[get_res_coors[k]] * mb_intensity[5].w);
	tmpdouble[5] = (mb_intensity[0].getVar()[get_res_coors[k]] * mb_intensity[0].w + mb_intensity[1].getVar()[get_res_coors[k]] * mb_intensity[1].w+ mb_intensity[2].getVar()[get_res_coors[k]] * mb_intensity[2].w + mb_intensity[3].getVar()[get_res_coors[k]] * mb_intensity[3].w+mb_intensity[4].getVar()[get_res_coors[k]] * mb_intensity[4].w)  ;
	fout.push_back(tmpdouble[4]/ (mb_intensity[7].w * tmpdouble[3]));  fout.push_back(tmpdouble[5]/ (mb_intensity[7].w * tmpdouble[3]));
	tmpdouble[1] = 	tmpdouble[0] * mb_intensity[7].w;
	tmpdouble[4] =tmpdouble[3] * mb_intensity[7].w;
	get_res_coors_par[2] =1; get_res_coors_par[3] =2;
	for(j=0;j<8;j++) mb_residual_par[j] =  mb_intensity[j].makeResidualGauss(get_res_coors_par);
	k=0;
	for(j=0;j<8;j++) fout.push_back( mb_residual_par[j].getVar()[k]);
	fout.push_back(mb_residual_par[7].getVar()[k] / tmpdouble[0]);fout.push_back((mb_residual_par[6].getVar()[k] * mb_intensity[6].w + mb_residual_par[5].getVar()[k] * mb_intensity[5].w)  / tmpdouble[1]  );
	fout.push_back( (mb_residual_par[0].getVar()[k] * mb_intensity[0].w + mb_residual_par[1].getVar()[k] * mb_intensity[1].w+ mb_residual_par[2].getVar()[k] * mb_intensity[2].w + mb_residual_par[3].getVar()[k] * mb_intensity[3].w+mb_residual_par[4].getVar()[k] * mb_intensity[4].w) / tmpdouble[1] );
	k=1;
	for(j=0;j<8;j++) fout.push_back( mb_residual_par[j].getVar()[k]);
	fout.push_back(mb_residual_par[7].getVar()[k] / tmpdouble[3]);fout.push_back((mb_residual_par[6].getVar()[k] * mb_intensity[6].w + mb_residual_par[5].getVar()[k] * mb_intensity[5].w)  / tmpdouble[4]  );
	fout.push_back( (mb_residual_par[0].getVar()[k] * mb_intensity[0].w + mb_residual_par[1].getVar()[k] * mb_intensity[1].w+ mb_residual_par[2].getVar()[k] * mb_intensity[2].w + mb_residual_par[3].getVar()[k] * mb_intensity[3].w+mb_residual_par[4].getVar()[k] * mb_intensity[4].w) / tmpdouble[4] );
	get_res_coors_par[2] =0; get_res_coors_par[3] =2;
	for(j=0;j<8;j++) mb_residual_par[j] =  mb_intensity[j].makeResidualGauss(get_res_coors_par);
	k=0;for(j=0;j<8;j++) fout.push_back( mb_residual_par[j].getVar()[k]);
	fout.push_back( mb_residual_par[7].getVar()[k] / tmpdouble[0]);fout.push_back((mb_residual_par[6].getVar()[k] * mb_intensity[6].w + mb_residual_par[5].getVar()[k] * mb_intensity[5].w)  / tmpdouble[1]  );
	fout.push_back( (mb_residual_par[0].getVar()[k] * mb_intensity[0].w + mb_residual_par[1].getVar()[k] * mb_intensity[1].w+ mb_residual_par[2].getVar()[k] * mb_intensity[2].w + mb_residual_par[3].getVar()[k] * mb_intensity[3].w+mb_residual_par[4].getVar()[k] * mb_intensity[4].w) / tmpdouble[1] );
	k=1;for(j=0;j<8;j++) fout.push_back( mb_residual_par[j].getVar()[k]);
	fout.push_back( mb_residual_par[7].getVar()[k] / tmpdouble[3]);fout.push_back((mb_residual_par[6].getVar()[k] * mb_intensity[6].w + mb_residual_par[5].getVar()[k] * mb_intensity[5].w)  / tmpdouble[4]  );
	fout.push_back( (mb_residual_par[0].getVar()[k] * mb_intensity[0].w + mb_residual_par[1].getVar()[k] * mb_intensity[1].w+ mb_residual_par[2].getVar()[k] * mb_intensity[2].w + mb_residual_par[3].getVar()[k] * mb_intensity[3].w+mb_residual_par[4].getVar()[k] * mb_intensity[4].w) / tmpdouble[4] );
	get_res_coors_par[2] =0; get_res_coors_par[3] =1;
	for(j=0;j<8;j++) mb_residual_par[j] =  mb_intensity[j].makeResidualGauss(get_res_coors_par);
	k=0;for(j=0;j<8;j++) fout.push_back( mb_residual_par[j].getVar()[k]);
	fout.push_back( mb_residual_par[7].getVar()[k] / tmpdouble[0]);fout.push_back((mb_residual_par[6].getVar()[k] * mb_intensity[6].w + mb_residual_par[5].getVar()[k] * mb_intensity[5].w)  / tmpdouble[1] );
	fout.push_back( (mb_residual_par[0].getVar()[k] * mb_intensity[0].w + mb_residual_par[1].getVar()[k] * mb_intensity[1].w+ mb_residual_par[2].getVar()[k] * mb_intensity[2].w + mb_residual_par[3].getVar()[k] * mb_intensity[3].w+mb_residual_par[4].getVar()[k] * mb_intensity[4].w) / tmpdouble[1] );
	k=1;for(j=0;j<8;j++) fout.push_back( mb_residual_par[j].getVar()[k]);
	fout.push_back( mb_residual_par[7].getVar()[k] / tmpdouble[3]);fout.push_back((mb_residual_par[6].getVar()[k] * mb_intensity[6].w + mb_residual_par[5].getVar()[k] * mb_intensity[5].w)  / tmpdouble[4]  );
	fout.push_back( (mb_residual_par[0].getVar()[k] * mb_intensity[0].w + mb_residual_par[1].getVar()[k] * mb_intensity[1].w+ mb_residual_par[2].getVar()[k] * mb_intensity[2].w + mb_residual_par[3].getVar()[k] * mb_intensity[3].w+mb_residual_par[4].getVar()[k] * mb_intensity[4].w) / tmpdouble[4] );
	get_res_coors[2] = 1;
	for(j=0;j<8;j++) mb_residual[j] =  mb_intensity[j].makeResidualGauss(get_res_coors);
	k=0;for(j=0;j<8;j++) fout.push_back( mb_residual[j].getVar()[k]);
	fout.push_back(mb_residual[7].getVar()[k] / tmpdouble[0]);fout.push_back((mb_residual[6].getVar()[k] * mb_intensity[6].w + mb_residual[5].getVar()[k] * mb_intensity[5].w)  / tmpdouble[1]  );
	fout.push_back( (mb_residual[0].getVar()[k] * mb_intensity[0].w + mb_residual[1].getVar()[k] * mb_intensity[1].w+ mb_residual[2].getVar()[k] * mb_intensity[2].w + mb_residual[3].getVar()[k] * mb_intensity[3].w+mb_residual[4].getVar()[k] * mb_intensity[4].w) / tmpdouble[1] );
	k=1;for(j=0;j<8;j++) fout.push_back( mb_residual[j].getVar()[k]);
	fout.push_back(mb_residual[7].getVar()[k] / tmpdouble[3]);fout.push_back((mb_residual[6].getVar()[k] * mb_intensity[6].w + mb_residual[5].getVar()[k] * mb_intensity[5].w)  / tmpdouble[4]  );
	fout.push_back( (mb_residual[0].getVar()[k] * mb_intensity[0].w + mb_residual[1].getVar()[k] * mb_intensity[1].w+ mb_residual[2].getVar()[k] * mb_intensity[2].w + mb_residual[3].getVar()[k] * mb_intensity[3].w+mb_residual[4].getVar()[k] * mb_intensity[4].w) / tmpdouble[4] );
	get_res_coors[2] = 0;
	for(j=0;j<8;j++) mb_residual[j] =  mb_intensity[j].makeResidualGauss(get_res_coors);
	k=0;for(j=0;j<8;j++) fout.push_back( mb_residual[j].getVar()[k]);
	fout.push_back(mb_residual[7].getVar()[k] / tmpdouble[0]);fout.push_back((mb_residual[6].getVar()[k] * mb_intensity[6].w + mb_residual[5].getVar()[k] * mb_intensity[5].w)  / tmpdouble[1]  );
	fout.push_back( (mb_residual[0].getVar()[k] * mb_intensity[0].w + mb_residual[1].getVar()[k] * mb_intensity[1].w+ mb_residual[2].getVar()[k] * mb_intensity[2].w + mb_residual[3].getVar()[k] * mb_intensity[3].w+mb_residual[4].getVar()[k] * mb_intensity[4].w) / tmpdouble[1] );
	k=1;for(j=0;j<8;j++) fout.push_back( mb_residual[j].getVar()[k]);
	fout.push_back(mb_residual[7].getVar()[k] / tmpdouble[3]);fout.push_back((mb_residual[6].getVar()[k] * mb_intensity[6].w + mb_residual[5].getVar()[k] * mb_intensity[5].w)  / tmpdouble[4] );
	fout.push_back( (mb_residual[0].getVar()[k] * mb_intensity[0].w + mb_residual[1].getVar()[k] * mb_intensity[1].w+ mb_residual[2].getVar()[k] * mb_intensity[2].w + mb_residual[3].getVar()[k] * mb_intensity[3].w+mb_residual[4].getVar()[k] * mb_intensity[4].w) / tmpdouble[4] );
}
// 0: budsize 1:mother size 2: DistanceRFP 3: DistanceGFP
void eight_gaussreport_half(Vector<double> &fout,  Tuple< GaussElem< Tuple<double, 4> >, 8> &b_intensity,  Tuple< GaussElem< Tuple<double, 4> >, 8> &m_intensity){
	GaussElem< Tuple<double, 2> > mb_residual[8];
	GaussElem< Tuple<double, 3> > mb_residual_par[8];
	double tmpdouble[6];
	Tuple<unsigned int, 2> get_res_coors;get_res_coors[0] =3;
	Tuple<unsigned int, 3> get_res_coors_par;get_res_coors_par[0] =3;
	unsigned int i,j,k;
	for(j=0;j<8;j++) fout.push_back(b_intensity[j].w);
	for(j=0;j<8;j++) {tmpdouble[0] = (b_intensity[j].w * b_intensity[j].w) / b_intensity[j].w2;  fout.push_back((ExOp::isValid(tmpdouble[0])) ? tmpdouble[0] : 0.0f);}
	k=0;for(j=0;j<8;j++) fout.push_back( b_intensity[j].mean[get_res_coors[k]] / b_intensity[j].w);
	//	fout.push_back( (b_intensity[5].mean[get_res_coors[k]] + b_intensity[6].mean[get_res_coors[k]]) / b_intensity[7].w );
	//	fout.push_back( (+b_intensity[0].mean[get_res_coors[k]] + b_intensity[1].mean[get_res_coors[k]]+b_intensity[2].mean[get_res_coors[k]] + b_intensity[3].mean[get_res_coors[k]]+b_intensity[4].mean[get_res_coors[k]] ) / b_intensity[7].w );
	k=0;for(j=0;j<8;j++) fout.push_back( m_intensity[j].mean[get_res_coors[k]] / b_intensity[j].w);
	//	fout.push_back( (b_intensity[5].mean[get_res_coors[k]] + b_intensity[6].mean[get_res_coors[k]]) / b_intensity[7].w );
	//	fout.push_back( (+b_intensity[0].mean[get_res_coors[k]] + b_intensity[1].mean[get_res_coors[k]]+b_intensity[2].mean[get_res_coors[k]] + b_intensity[3].mean[get_res_coors[k]]+b_intensity[4].mean[get_res_coors[k]] ) / b_intensity[7].w );
	k=0;for(j=0;j<8;j++) fout.push_back( b_intensity[j].getVar()[get_res_coors[k]]);
	tmpdouble[0] = b_intensity[7].getVar()[get_res_coors[k]];
	tmpdouble[1] = (b_intensity[6].getVar()[get_res_coors[k]] * b_intensity[6].w + b_intensity[5].getVar()[get_res_coors[k]] * b_intensity[5].w);
	tmpdouble[2] = (b_intensity[0].getVar()[get_res_coors[k]] * b_intensity[0].w + b_intensity[1].getVar()[get_res_coors[k]] * b_intensity[1].w+ b_intensity[2].getVar()[get_res_coors[k]] * b_intensity[2].w + b_intensity[3].getVar()[get_res_coors[k]] * b_intensity[3].w+b_intensity[4].getVar()[get_res_coors[k]] * b_intensity[4].w)  ;
	fout.push_back( tmpdouble[1]/ (b_intensity[7].w * tmpdouble[0])); fout.push_back(tmpdouble[2]/ (b_intensity[7].w * tmpdouble[0]));
	k=0;for(j=0;j<8;j++) fout.push_back( m_intensity[j].getVar()[get_res_coors[k]]);
	tmpdouble[3] = m_intensity[7].getVar()[get_res_coors[k]];
	tmpdouble[4] = (m_intensity[6].getVar()[get_res_coors[k]] * m_intensity[6].w + m_intensity[5].getVar()[get_res_coors[k]] * m_intensity[5].w);
	tmpdouble[5] = (m_intensity[0].getVar()[get_res_coors[k]] * m_intensity[0].w + m_intensity[1].getVar()[get_res_coors[k]] * m_intensity[1].w+ m_intensity[2].getVar()[get_res_coors[k]] * m_intensity[2].w + m_intensity[3].getVar()[get_res_coors[k]] * m_intensity[3].w+m_intensity[4].getVar()[get_res_coors[k]] * m_intensity[4].w)  ;
	fout.push_back(tmpdouble[4]/ (m_intensity[7].w * tmpdouble[3]));  fout.push_back(tmpdouble[5]/ (m_intensity[7].w * tmpdouble[3]));
	tmpdouble[1] = 	tmpdouble[0] * b_intensity[7].w;
	tmpdouble[4] =tmpdouble[3] * m_intensity[7].w;
	get_res_coors_par[1] =1; get_res_coors_par[2] =2;
	k=0;for(j=0;j<8;j++) mb_residual_par[j] =  b_intensity[j].makeResidualGauss(get_res_coors_par);
	for(j=0;j<8;j++) fout.push_back( mb_residual_par[j].getVar()[k]);
	fout.push_back(mb_residual_par[7].getVar()[k] / tmpdouble[0]);fout.push_back((mb_residual_par[6].getVar()[k] * b_intensity[6].w + mb_residual_par[5].getVar()[k] * b_intensity[5].w)  / tmpdouble[1]  );
	fout.push_back( (mb_residual_par[0].getVar()[k] * b_intensity[0].w + mb_residual_par[1].getVar()[k] * b_intensity[1].w+ mb_residual_par[2].getVar()[k] * b_intensity[2].w + mb_residual_par[3].getVar()[k] * b_intensity[3].w+mb_residual_par[4].getVar()[k] * b_intensity[4].w) / tmpdouble[1] );
	k=0;for(j=0;j<8;j++) mb_residual_par[j] =  m_intensity[j].makeResidualGauss(get_res_coors_par);
	for(j=0;j<8;j++) fout.push_back( mb_residual_par[j].getVar()[k]);
	fout.push_back(mb_residual_par[7].getVar()[k] / tmpdouble[3]);fout.push_back((mb_residual_par[6].getVar()[k] * b_intensity[6].w + mb_residual_par[5].getVar()[k] * b_intensity[5].w)  / tmpdouble[4]  );
	fout.push_back( (mb_residual_par[0].getVar()[k] * b_intensity[0].w + mb_residual_par[1].getVar()[k] * b_intensity[1].w+ mb_residual_par[2].getVar()[k] * b_intensity[2].w + mb_residual_par[3].getVar()[k] * b_intensity[3].w+mb_residual_par[4].getVar()[k] * b_intensity[4].w) / tmpdouble[4] );
	get_res_coors_par[1] =0; get_res_coors_par[2] =2;
	k=0;for(j=0;j<8;j++) mb_residual_par[j] =  b_intensity[j].makeResidualGauss(get_res_coors_par);
	for(j=0;j<8;j++) fout.push_back( mb_residual_par[j].getVar()[k]);
	fout.push_back( mb_residual_par[7].getVar()[k] / tmpdouble[0]);fout.push_back((mb_residual_par[6].getVar()[k] * b_intensity[6].w + mb_residual_par[5].getVar()[k] * b_intensity[5].w)  / tmpdouble[1]  );
	fout.push_back( (mb_residual_par[0].getVar()[k] * b_intensity[0].w + mb_residual_par[1].getVar()[k] * b_intensity[1].w+ mb_residual_par[2].getVar()[k] * b_intensity[2].w + mb_residual_par[3].getVar()[k] * b_intensity[3].w+mb_residual_par[4].getVar()[k] * b_intensity[4].w) / tmpdouble[1] );
	k=0;for(j=0;j<8;j++) mb_residual_par[j] =  m_intensity[j].makeResidualGauss(get_res_coors_par);
	for(j=0;j<8;j++) fout.push_back( mb_residual_par[j].getVar()[k]);
	fout.push_back( mb_residual_par[7].getVar()[k] / tmpdouble[3]);fout.push_back((mb_residual_par[6].getVar()[k] * b_intensity[6].w + mb_residual_par[5].getVar()[k] * b_intensity[5].w)  / tmpdouble[4]  );
	fout.push_back( (mb_residual_par[0].getVar()[k] * b_intensity[0].w + mb_residual_par[1].getVar()[k] * b_intensity[1].w+ mb_residual_par[2].getVar()[k] * b_intensity[2].w + mb_residual_par[3].getVar()[k] * b_intensity[3].w+mb_residual_par[4].getVar()[k] * b_intensity[4].w) / tmpdouble[4] );
	get_res_coors_par[1] =0; get_res_coors_par[2] =1;
	k=0;for(j=0;j<8;j++) mb_residual_par[j] =  b_intensity[j].makeResidualGauss(get_res_coors_par);
	for(j=0;j<8;j++) fout.push_back( mb_residual_par[j].getVar()[k]);
	fout.push_back( mb_residual_par[7].getVar()[k] / tmpdouble[0]);fout.push_back((mb_residual_par[6].getVar()[k] * b_intensity[6].w + mb_residual_par[5].getVar()[k] * b_intensity[5].w)  / tmpdouble[1] );
	fout.push_back( (mb_residual_par[0].getVar()[k] * b_intensity[0].w + mb_residual_par[1].getVar()[k] * b_intensity[1].w+ mb_residual_par[2].getVar()[k] * b_intensity[2].w + mb_residual_par[3].getVar()[k] * b_intensity[3].w+mb_residual_par[4].getVar()[k] * b_intensity[4].w) / tmpdouble[1] );
	k=0;for(j=0;j<8;j++) mb_residual_par[j] =  m_intensity[j].makeResidualGauss(get_res_coors_par);
	for(j=0;j<8;j++) fout.push_back( mb_residual_par[j].getVar()[k]);
	fout.push_back( mb_residual_par[7].getVar()[k] / tmpdouble[3]);fout.push_back((mb_residual_par[6].getVar()[k] * b_intensity[6].w + mb_residual_par[5].getVar()[k] * b_intensity[5].w)  / tmpdouble[4]  );
	fout.push_back( (mb_residual_par[0].getVar()[k] * b_intensity[0].w + mb_residual_par[1].getVar()[k] * b_intensity[1].w+ mb_residual_par[2].getVar()[k] * b_intensity[2].w + mb_residual_par[3].getVar()[k] * b_intensity[3].w+mb_residual_par[4].getVar()[k] * b_intensity[4].w) / tmpdouble[4] );
	get_res_coors[1] = 1;
	k=0;for(j=0;j<8;j++) mb_residual[j] =  b_intensity[j].makeResidualGauss(get_res_coors);
	for(j=0;j<8;j++) fout.push_back( mb_residual[j].getVar()[k]);
	fout.push_back(mb_residual[7].getVar()[k] / tmpdouble[0]);fout.push_back((mb_residual[6].getVar()[k] * b_intensity[6].w + mb_residual[5].getVar()[k] * b_intensity[5].w)  / tmpdouble[1]  );
	fout.push_back( (mb_residual[0].getVar()[k] * b_intensity[0].w + mb_residual[1].getVar()[k] * b_intensity[1].w+ mb_residual[2].getVar()[k] * b_intensity[2].w + mb_residual[3].getVar()[k] * b_intensity[3].w+mb_residual[4].getVar()[k] * b_intensity[4].w) / tmpdouble[1] );
	k=0;for(j=0;j<8;j++) mb_residual[j] =  m_intensity[j].makeResidualGauss(get_res_coors);
	for(j=0;j<8;j++) fout.push_back( mb_residual[j].getVar()[k]);
	fout.push_back(mb_residual[7].getVar()[k] / tmpdouble[3]);fout.push_back((mb_residual[6].getVar()[k] * b_intensity[6].w + mb_residual[5].getVar()[k] * b_intensity[5].w)  / tmpdouble[4]  );
	fout.push_back( (mb_residual[0].getVar()[k] * b_intensity[0].w + mb_residual[1].getVar()[k] * b_intensity[1].w+ mb_residual[2].getVar()[k] * b_intensity[2].w + mb_residual[3].getVar()[k] * b_intensity[3].w+mb_residual[4].getVar()[k] * b_intensity[4].w) / tmpdouble[4] );
	get_res_coors[1] = 0;
	k=0;for(j=0;j<8;j++) mb_residual[j] =  b_intensity[j].makeResidualGauss(get_res_coors);
	for(j=0;j<8;j++) fout.push_back( mb_residual[j].getVar()[k]);
	fout.push_back(mb_residual[7].getVar()[k] / tmpdouble[0]);fout.push_back((mb_residual[6].getVar()[k] * b_intensity[6].w + mb_residual[5].getVar()[k] * b_intensity[5].w)  / tmpdouble[1]  );
	fout.push_back( (mb_residual[0].getVar()[k] * b_intensity[0].w + mb_residual[1].getVar()[k] * b_intensity[1].w+ mb_residual[2].getVar()[k] * b_intensity[2].w + mb_residual[3].getVar()[k] * b_intensity[3].w+mb_residual[4].getVar()[k] * b_intensity[4].w) / tmpdouble[1] );
	k=0;for(j=0;j<8;j++) mb_residual[j] =  m_intensity[j].makeResidualGauss(get_res_coors);
	for(j=0;j<8;j++) fout.push_back( mb_residual[j].getVar()[k]);
	fout.push_back(mb_residual[7].getVar()[k] / tmpdouble[3]);fout.push_back((mb_residual[6].getVar()[k] * b_intensity[6].w + mb_residual[5].getVar()[k] * b_intensity[5].w)  / tmpdouble[4] );
	fout.push_back( (mb_residual[0].getVar()[k] * b_intensity[0].w + mb_residual[1].getVar()[k] * b_intensity[1].w+ mb_residual[2].getVar()[k] * b_intensity[2].w + mb_residual[3].getVar()[k] * b_intensity[3].w+mb_residual[4].getVar()[k] * b_intensity[4].w) / tmpdouble[4] );
}

	Taskscope<TASK_PROFILES>::Taskscope(): show(false),inter(false),cluster(false),no_daugther(false),hypercl(NULL),no_confidence(2.0f),does_bining(false){
		file_in[2]=NULL;
		for(unsigned int i=0; i< 6;i++) factor[i] = 1.0f;
	}
	void Taskscope<TASK_PROFILES>::nbaddtoken(char const * const token, int& min, int& max){
		switch(*token){
			case '\0': min =1; max =2; break;
			case 's': min =0; break;
			case 'C': min =0; break;
			case 'i': min =1; break;
			case 'p': min =1; max =1; break;
			case 'D': min =1; break;
			case 'P': min =1; break;
			case 'L': min =0; min =1; break;
			case 'd': min =0; break;
			case 'B': min =2; break;
			case 'J': min =1; break;
			case 'S': min =6; break;
		}
	}
	void Taskscope<TASK_PROFILES>::store(char* const * token, int nbtoken){ 
		switch(token[0][1]){
			case 's': show = true; break;
			case 'p': inter = true; break;
			case 'P': pair_wise = token[1]; break;
			case 'i': no_confidence = atof(token[1]); break;
			case 'C': cluster = true; break;
			case 'D': file_in[2] = token[1]; break;
			case 'd': no_daugther = true; break;
			case 'L': hypercl_out = (nbtoken < 1) ? NULL : token[1]; hypercl_near = (token[0][2] == 'N');  break;
			case 'B': does_bining =true; daval_sepsep = atof(token[1]); daval_sepsepi = atof(token[2]); break;
			case 'J':  svm_joachims = token[1]; break;
			case 'S': factor[0] =  1.0f / atof(token[1]); factor[1] =  1.0f / atof(token[2]) ;factor[2] =  1.0f / atof(token[3]); factor[3] =  1.0f / atof(token[4]) ;factor[4] =  1.0f / atof(token[5]); factor[5] =  1.0f / atof(token[6]) ; break;
		}
	}
	int Taskscope<TASK_PROFILES>::defstore(char* const * token, int nbtoken){ 
		file_in[0] = (nbtoken > 0) ? token[0] : NULL;
		file_in[1] = (nbtoken > 1) ? token[1] : NULL;
		if ((nbtoken < 1)&&(!inter)) {inter = true; printf("Missing arguments, modifying scope: (try \"./PMProfiles -h\" for help)\n");}
		ProgressBarPrint daprogbar; daprogbar.lenght = 20;
		unsigned int i,j,k,l;
		char buffer[65536];
		//	 char* fakeargs[] = {"./PMProfiles", "../../../../dataA_list", "../../../../dataA_table.txt"}; // 
		//	(sizeof(fakeargs) / sizeof(char*) ,fakesubstitution(sizeof(fakeargs) / sizeof(char*) , fakeargs, "HOwt_plate01_009023"));	
		double tmpdouble;
		SerialStore<unsigned int> rstd("./PMProfiles_scope.scp");
		Vector<char> str;
		FILE* clust_cdt = (cluster) ? fopen("auto_cluster.cdt", "w+") : NULL;
		Forest<double,2> da_clust;
		Vector< Tuple< WeightElem<double,2> , 120> > to_clust;
		Tuple<double, 120> clust_simple_in;
		Tuple< WeightElem<double,2> , 120> clust_in;
		myHashmap<unsigned int, pair< Tuple<double, 6> ,double > > dual_info;
		Vector< KeyElem< pair< unsigned int, Tuple<double, 4 > > , Tuple<double, 12 > > > mother_buds;
		Vector< KeyElem< pair< unsigned int, Tuple<double, 3 > > , Tuple<double, 5 > > > loners;
		KeyElem< pair< unsigned int, Tuple<double, 4 > > , Tuple<double, 12 > > mother_buds_input;
		KeyElem< pair< unsigned int, Tuple<double, 3 > > , Tuple<double, 5 > > loners_input;
		Tuple<double,4> supercount; ExOp::toZero(supercount);
		Vector< char* > to_clust_header;
		unsigned int resol = 10;
		Weight gaus_pix[32];
		gaus_pix[resol]= 1.0f;
		for(i =1;i<=resol ;i++) {gaus_pix[resol-i] = Weight(exp(-3.0f*i*i));gaus_pix[resol+i] = Weight(exp(-3.0f*i*i));}
		myHashmap<string, unsigned int> daHclasses;
		Vector< GaussElem<Tuple<double, 120 > > > daHGauss;
		WeightElem<double,2> da_fix;
		if (inter){
			if (rstd.has(0)){
				rstd.load(0, str);
				printf("enter the data table path and prefix: currently =%s\n", &(str[0]));
			} else printf("enter the data table path and prefix:\n");
			scanf(" %[^\n]", buffer);
			l = strlen(buffer);
			if (l){
				str.setSize(l+1);
				memcpy(&(str[0]),buffer, l+1);
				rstd.save(0, str);
			}
			if (rstd.has(1)){
				rstd.load(1, str);
				printf("enter the data table suffix: currently =%s\n", &(str[0]));
			} else printf("enter the data table suffix:\n");
			scanf(" %[^\n]", buffer);
			l = strlen(buffer);
			if (l){
				str.setSize(l+1);
				memcpy(&(str[0]),buffer, l+1);
				rstd.save(1, str);
			}
			if (rstd.has(2)){
				rstd.load(2, str);
				printf("enter the path to classification table: currently =%s\n", &(str[0]));
			} else printf("enter the data table suffix:\n");
			scanf(" %[^\n]", buffer);
			l = strlen(buffer);
			if (l){
				str.setSize(l+1);
				memcpy(&(str[0]),buffer, l+1);
				rstd.save(2, str);
			}
		}else{
			myHashmap<string, unsigned int> daHclasses;
			myHashmap<string, unsigned int> daHclasses_names_back;
			vector< string > daHclasses_names;
			rstd.load(2, str);
			if (str.size() > 0){
				FILE *hyper = fopen(&(str[0]),"r+");
				while(fscanf(hyper,"%[^\t\n]\t%[^\t\n]\n", buffer, buffer + 1024)  == 2){
					i = daHclasses_names_back.find(string(buffer+1024));
					if (i == 0xFFFFFFFF) {i = daHclasses_names.size(); daHclasses_names_back[string(buffer+1024)] = i; daHclasses_names.push_back(string(buffer+1024));}
					else i = daHclasses_names_back.deref(i);
					daHclasses[string(buffer)] = i;
				}
				daHclasses_names.push_back(string("Undefined Class"));
				daHclasses_names.push_back(string("Total Class"));
				daHGauss.setSize(daHclasses_names.size());
				ExOp::toZero(daHGauss);
			}
			FILE* f = fopen(file_in[0], "r+");
			rstd.load(0, str);
			unsigned int prefix_length = strlen(&(str[0]));
			FILE* g = file_in[1] ?  fopen(file_in[1], "w+") : stdout;
			int dataindex[64];	
			double pix[32];
			l=0;
			Madstructs::Table tb;
			double lone_prob;
			double class_prob;
			double size_stage;
			double cell_stage;
			fprintf(g,"ID\tNAME");
			memcpy(buffer,"BM_CNT_", sizeof(char)*8);for(i=0;i<4*(resol+1);i++) if ((i % (resol +1)) !=0) fprintf(g,"\t%s%c%i", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L') ,(i-1)%(resol+1)); else fprintf(g,"\t%s%c", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L'));
			memcpy(buffer,"BM_INT_", sizeof(char)*8);for(i=0;i<4*(resol+1);i++) if ((i % (resol +1)) !=0) fprintf(g,"\t%s%c%i", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L') ,(i-1)%(resol+1)); else fprintf(g,"\t%s%c", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L'));
			memcpy(buffer,"BM_SEF_", sizeof(char)*8);for(i=0;i<4*(resol+1);i++) if ((i % (resol +1)) !=0) fprintf(g,"\t%s%c%i", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L') ,(i-1)%(resol+1)); else fprintf(g,"\t%s%c", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L'));
			memcpy(buffer,"BM_MCT_", sizeof(char)*8);for(i=0;i<4*(resol+1);i++) if ((i % (resol +1)) !=0) fprintf(g,"\t%s%c%i", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L') ,(i-1)%(resol+1)); else fprintf(g,"\t%s%c", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L'));
			memcpy(buffer,"BM_EDG_", sizeof(char)*8);for(i=0;i<4*(resol+1);i++) if ((i % (resol +1)) !=0) fprintf(g,"\t%s%c%i", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L') ,(i-1)%(resol+1)); else fprintf(g,"\t%s%c", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L'));
			memcpy(buffer,"BM_CEN_", sizeof(char)*8);for(i=0;i<4*(resol+1);i++) if ((i % (resol +1)) !=0) fprintf(g,"\t%s%c%i", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L') ,(i-1)%(resol+1)); else fprintf(g,"\t%s%c", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L'));
			memcpy(buffer,"BM_NEC_", sizeof(char)*8);for(i=0;i<3*(resol+1);i++) if ((i % (resol +1)) !=0) fprintf(g,"\t%s%c%i", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L') ,(i-1)%(resol+1)); else fprintf(g,"\t%s%c", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L'));
			memcpy(buffer,"BV_INT_", sizeof(char)*8);for(i=0;i<4*(resol+1);i++) if ((i % (resol +1)) !=0) fprintf(g,"\t%s%c%i", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L') ,(i-1)%(resol+1)); else fprintf(g,"\t%s%c", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L'));
			memcpy(buffer,"BV_SEF_", sizeof(char)*8);for(i=0;i<4*(resol+1);i++) if ((i % (resol +1)) !=0) fprintf(g,"\t%s%c%i", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L') ,(i-1)%(resol+1)); else fprintf(g,"\t%s%c", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L'));
			memcpy(buffer,"BV_MCT_", sizeof(char)*8);for(i=0;i<4*(resol+1);i++) if ((i % (resol +1)) !=0) fprintf(g,"\t%s%c%i", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L') ,(i-1)%(resol+1)); else fprintf(g,"\t%s%c", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L'));
			memcpy(buffer,"BV_EDG_", sizeof(char)*8);for(i=0;i<4*(resol+1);i++) if ((i % (resol +1)) !=0) fprintf(g,"\t%s%c%i", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L') ,(i-1)%(resol+1)); else fprintf(g,"\t%s%c", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L'));
			memcpy(buffer,"BV_CEN_", sizeof(char)*8);for(i=0;i<4*(resol+1);i++) if ((i % (resol +1)) !=0) fprintf(g,"\t%s%c%i", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L') ,(i-1)%(resol+1)); else fprintf(g,"\t%s%c", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L'));
			memcpy(buffer,"BV_NEC_", sizeof(char)*8);for(i=0;i<3*(resol+1);i++) if ((i % (resol +1)) !=0) fprintf(g,"\t%s%c%i", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L') ,(i-1)%(resol+1)); else fprintf(g,"\t%s%c", buffer, (i <(2*resol+2)) ? ((i <(resol+1)) ? 'B' : 'M') : ((i <(3*resol+3)) ? 'D' : 'L'));
			if (cluster){
				fprintf(clust_cdt,"GID\tID\tNAME\tGWEIGHT");
				/*      memcpy(buffer,"BM_INT_", sizeof(char)*8);for(i=0;i<4*resol;i+=resol) fprintf(clust_cdt,"\t%s%c", buffer, (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L'));
				 memcpy(buffer,"BM_SEF_", sizeof(char)*8);for(i=0;i<4*resol;i+=resol) fprintf(clust_cdt,"\t%s%c", buffer, (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L'));
				 memcpy(buffer,"BM_MCT_", sizeof(char)*8);for(i=0;i<4*resol;i+=resol) fprintf(clust_cdt,"\t%s%c", buffer, (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L'));
				 memcpy(buffer,"BM_EDG_", sizeof(char)*8);for(i=0;i<4*resol;i+=resol) fprintf(clust_cdt,"\t%s%c", buffer, (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L'));
				 memcpy(buffer,"BM_CEN_", sizeof(char)*8);for(i=0;i<4*resol;i+=resol) fprintf(clust_cdt,"\t%s%c", buffer, (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L'));
				 memcpy(buffer,"BM_NEC_", sizeof(char)*8);for(i=0;i<2*resol;i+=resol) fprintf(clust_cdt,"\t%s%c", buffer, (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L'));
				 
				 memcpy(buffer,"BV_INT_", sizeof(char)*8);for(i=0;i<4*resol;i+=resol) fprintf(clust_cdt,"\t%s%c", buffer, (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L'));
				 memcpy(buffer,"BV_SEF_", sizeof(char)*8);for(i=0;i<4*resol;i+=resol) fprintf(clust_cdt,"\t%s%c", buffer, (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L'));
				 memcpy(buffer,"BV_MCT_", sizeof(char)*8);for(i=0;i<4*resol;i+=resol) fprintf(clust_cdt,"\t%s%c", buffer, (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L'));
				 memcpy(buffer,"BV_EDG_", sizeof(char)*8);for(i=0;i<4*resol;i+=resol) fprintf(clust_cdt,"\t%s%c", buffer, (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L'));
				 memcpy(buffer,"BV_CEN_", sizeof(char)*8);for(i=0;i<4*resol;i+=resol) fprintf(clust_cdt,"\t%s%c", buffer, (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L'));
				 memcpy(buffer,"BV_NEC_", sizeof(char)*8);for(i=0;i<2*resol;i+=resol) fprintf(clust_cdt,"\t%s%c", buffer, (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L'));
				 */	
				/*
				 //time profiles
				 fprintf(clust_cdt,"GID\tID\tNAME\tGWEIGHT");
				 
				 memcpy(buffer,"BM_INT_", sizeof(char)*8);for(i=0;i<3;i++) fprintf(clust_cdt,"\t%s%c", buffer, (i <2) ? ((i <1) ? 'B' : 'M') : ((i <3) ? 'L' : 'L') );
				 memcpy(buffer,"BM_SEF_", sizeof(char)*8);for(i=0;i<3;i++) fprintf(clust_cdt,"\t%s%c", buffer, (i <2) ? ((i <1) ? 'B' : 'M') : ((i <3) ? 'L' : 'L') );
				 memcpy(buffer,"BM_MCT_", sizeof(char)*8);for(i=0;i<3;i++) fprintf(clust_cdt,"\t%s%c", buffer, (i <2) ? ((i <1) ? 'B' : 'M') : ((i <3) ? 'L' : 'L') );
				 memcpy(buffer,"BM_EDG_", sizeof(char)*8);for(i=0;i<3;i++) fprintf(clust_cdt,"\t%s%c", buffer, (i <2) ? ((i <1) ? 'B' : 'M') : ((i <3) ? 'L' : 'L') );
				 memcpy(buffer,"BM_CEN_", sizeof(char)*8);for(i=0;i<3;i++) fprintf(clust_cdt,"\t%s%c", buffer, (i <2) ? ((i <1) ? 'B' : 'M') : ((i <3) ? 'L' : 'L') );
				 memcpy(buffer,"BM_NEC_", sizeof(char)*8);for(i=0;i<2;i++) fprintf(clust_cdt,"\t%s%c", buffer, (i <2) ? ((i <1) ? 'B' : 'M') : ((i <3) ? 'L' : 'L') );
				 */
				memcpy(buffer,"BM_INT_", sizeof(char)*8);for(i=0;i<2*resol;i++) fprintf(clust_cdt,"\t%s%c%i", buffer, (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%(resol));
				memcpy(buffer,"BM_SEF_", sizeof(char)*8);for(i=0;i<2*resol;i++) fprintf(clust_cdt,"\t%s%c%i", buffer, (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%(resol));
				memcpy(buffer,"BM_MCT_", sizeof(char)*8);for(i=0;i<2*resol;i++) fprintf(clust_cdt,"\t%s%c%i", buffer, (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%(resol));
				memcpy(buffer,"BM_EDG_", sizeof(char)*8);for(i=0;i<2*resol;i++) fprintf(clust_cdt,"\t%s%c%i", buffer, (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%(resol));
				memcpy(buffer,"BM_CEN_", sizeof(char)*8);for(i=0;i<2*resol;i++) fprintf(clust_cdt,"\t%s%c%i", buffer, (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%(resol));
				memcpy(buffer,"BM_NEC_", sizeof(char)*8);for(i=0;i<2*resol;i++) fprintf(clust_cdt,"\t%s%c%i", buffer, (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%(resol));
				memcpy(buffer,"BV_INT_", sizeof(char)*8);for(i=0;i<2*resol;i++) fprintf(clust_cdt,"\t%s%c%i", buffer, (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%(resol));
				memcpy(buffer,"BV_SEF_", sizeof(char)*8);for(i=0;i<2*resol;i++) fprintf(clust_cdt,"\t%s%c%i", buffer, (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%(resol));
				memcpy(buffer,"BV_MCT_", sizeof(char)*8);for(i=0;i<2*resol;i++) fprintf(clust_cdt,"\t%s%c%i", buffer, (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%(resol));
				memcpy(buffer,"BV_EDG_", sizeof(char)*8);for(i=0;i<2*resol;i++) fprintf(clust_cdt,"\t%s%c%i", buffer, (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%(resol));
				memcpy(buffer,"BV_CEN_", sizeof(char)*8);for(i=0;i<2*resol;i++) fprintf(clust_cdt,"\t%s%c%i", buffer, (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%(resol));
				memcpy(buffer,"BV_NEC_", sizeof(char)*8);for(i=0;i<2*resol;i++) fprintf(clust_cdt,"\t%s%c%i", buffer, (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%(resol));
				fprintf(clust_cdt,"\nEWEIGHT\t\t\t");
				for(i=0;i<clust_in.getSize()*2;i++) fprintf(clust_cdt,"\t%f",1.0f);
				fprintf(clust_cdt,"\n");
			}
			/*
			 for(i=0;i<4*resol;i++) fprintf(g,"\tLM_CNT_%c%i", (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%resol);
			 for(i=0;i<4*resol;i++) fprintf(g,"\tLM_INT_%c%i", (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%resol);
			 for(i=0;i<4*resol;i++) fprintf(g,"\tLM_SEF_%c%i", (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%resol);
			 for(i=0;i<4*resol;i++) fprintf(g,"\tLM_MCT_%c%i", (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%resol);
			 for(i=0;i<4*resol;i++) fprintf(g,"\tLM_EDG_%c%i", (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%resol);
			 for(i=0;i<4*resol;i++) fprintf(g,"\tLM_CEN_%c%i", (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%resol);
			 for(i=0;i<3*resol;i++) fprintf(g,"\tLM_NEC_%c%i", (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%resol);
			 
			 for(i=0;i<4*resol;i++) fprintf(g,"\tLV_INT_%c%i", (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%resol);
			 for(i=0;i<4*resol;i++) fprintf(g,"\tLV_SEF_%c%i", (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%resol);
			 for(i=0;i<4*resol;i++) fprintf(g,"\tLV_MCT_%c%i", (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%resol);
			 for(i=0;i<4*resol;i++) fprintf(g,"\tLV_EDG_%c%i", (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%resol);
			 for(i=0;i<4*resol;i++) fprintf(g,"\tLV_CEN_%c%i", (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%resol);
			 for(i=0;i<3*resol;i++) fprintf(g,"\tLV_NEC_%c%i", (i <(2*resol)) ? ((i <(resol)) ? 'B' : 'M') : ((i <(3*resol)) ? 'D' : 'L') ,i%resol);*/
			fprintf(g,"\n");
			char sep;
			//		DataGrid< WeightElem<Tuple<double, 6>, 2> , 1> profiles[2];
			//		DataGrid< WeightElem<Tuple<double, 5>, 2> , 1> profiles_l;
			//		DataGrid< WeightElem<double, 2> , 1> profiles_d;
			WeightElem<Tuple<double, 6>, 2> prior[2];
			WeightElem<Tuple<double, 5>, 2> prior_l;
			WeightElem<double, 2> prior_d;
			Tuple<double, 6> prior_for_var;
			prior_for_var[0] = exp(-6.0f);prior_for_var[1] = exp(-6.0f);prior_for_var[2] = exp(-6.0f);prior_for_var[3] = exp(-6.0f);prior_for_var[4] = exp(-6.0f);prior_for_var[5] = exp(-6.0f);
			Tuple<unsigned int,1u> gauss_dims;unsigned int coor[1];
			Tuple<double, 1> tmptmp;
			GaussianDistribution<1> fun_blur;  
			fun_blur.EMinit();
			tmptmp[0]  = -64.0f; fun_blur.EMregist(tmptmp);tmptmp[0]  = 64.0f; fun_blur.EMregist(tmptmp);
			fun_blur.EMfinit();
			gauss_dims[0] = 1024;
			//		profiles_m.setSizes(coor);profiles_d.setSizes(coor);profiles_b.setSizes(coor);profiles_l.setSizes(gauss_dims);
			//	GaussianDistribution<3> covar; 
			//	GaussianDistribution<3> covcov; 
			//	covcov.EMinit();
			Tuple<double, 6> data;
			Tuple<double, 5> data_l;
			Tuple<double, 12 > data_dual;
			double data_d;
			l=0;
			Tuple<double, 3> data_cov;
			Vector<unsigned int> locind;
			vector<string> prot_orf;
			vector<string> prot_name;
			while(2 == fscanf(f,"%[^\n\t\r]%c",buffer+prefix_length, &sep)){	
				rstd.load(0, str);
				memcpy(buffer,&(str[0]), prefix_length);
				i = strlen(buffer);
				rstd.load(1, str);
				j = strlen(&(str[0]));
				memcpy(buffer+i,&(str[0]), j);
				buffer[i+j] = '\0';				
				printf("path= %s",buffer);
				tb.wiseload(buffer);
				fscanf(f,"%[^\n]\n",buffer);
				for(i=0;buffer[i] !='\t';i++);
				buffer[i] ='\0';
				prot_name.push_back(string(buffer));
				buffer[i] ='\t';
				for(i++;buffer[i] !=';';i++);
				for(j=i+1;buffer[j] !=';';j++);
				buffer[j] = '\0';
				unsigned int hashite = daHclasses.find(string(buffer+i+1));
				prot_orf.push_back(string(buffer+i+1));
				buffer[j] = ';';
				locind.push_back((hashite <  daHclasses.getSize()) ? ((daHclasses.deref(hashite) != 0xFFFFFFFF) ? daHclasses.deref(hashite) : daHclasses_names.size() - 2 )   : daHclasses_names.size() - 2);
				if (daHGauss.size() > 0) printf("(hyperlabel = %s)\n",daHclasses_names[locind[locind.size()-1]].c_str());
				else printf("\n");
				fprintf(g,"%s", buffer,tb.nbrows);
				if ((cluster) || (file_in[2] != NULL)) {
					sprintf(buffer+32000,"GENE%iX\t%s\t%f", to_clust.size()+1u, buffer,1.0f);
					to_clust_header.push_back(cloneString(buffer+32000));
				}
				i=0;dataindex[i] = tb.findCol("Cell_type");
				i=1;dataindex[i] = tb.findCol("Cell_prob");
				i=2;dataindex[i] = tb.findCol("Relation_prob");
				i=3;dataindex[i] = tb.findCol("Area");
				i=4;dataindex[i] = tb.findCol("Rel_Cell_area"); 				
				i=5;dataindex[i] = tb.findCol("Intensity_MEAN"); if (tb.rowformat[2*dataindex[i] ] != 'f') {printf("missing Intensity_MEAN!\n");exit(1);}
				i=6;dataindex[i] = tb.findCol("N_SELF_DST_MEAN"); if (tb.rowformat[2*dataindex[i] ] != 'f') {printf("missing N_SELF_DST_MEAN!\n");exit(1);}
				//i=7;dataindex[i] = tb.findCol("G_SELF_DST_MEAN"); if (tb.rowformat[2*dataindex[i] ] != 'f') {printf("missing N_MASC_DST_MEAN!\n");exit(1);}
				i=7;dataindex[i] = tb.findCol("N_MASC_DST_MEAN"); if (tb.rowformat[2*dataindex[i] ] != 'f') {printf("missing N_MASC_DST_MEAN!\n");exit(1);}
				i=8;dataindex[i] = tb.findCol("N_EDGE_DST_MEAN"); if (tb.rowformat[2*dataindex[i] ] != 'f') {printf("missing N_EDGE_DST_MEAN!\n");exit(1);}
				i=9;dataindex[i] = tb.findCol("N_CENT_DST_MEAN"); if (tb.rowformat[2*dataindex[i] ] != 'f') {printf("missing N_CENT_DST_MEAN!\n");exit(1);}
				i=10;dataindex[i] = tb.findCol("N_BUDN_DST_MEAN"); if (tb.rowformat[2*dataindex[i] ] != 'f') {printf("missing N_BUDN_DST_MEAN!\n");exit(1);	}
				i=11;dataindex[i] = tb.findCol("Green_MEAN"); if (tb.rowformat[2*dataindex[i] ] != 'f')  {printf("missing Green_MEAN!\n");exit(1);	}	
				i=12;dataindex[i] = tb.findCol("CellID");  if (tb.rowformat[2*dataindex[i] ] != 'i')  {printf("missing CellID!\n");exit(1);	}	
				i=13;dataindex[i] = tb.findCol("Rel_Cell_ID");  if (tb.rowformat[2*dataindex[i] ] != 'i')  {printf("missing Rel_Cell_ID!\n");exit(1);	}	
				i=14;dataindex[i] = tb.findCol("Red_MEAN"); if (tb.rowformat[2*dataindex[i] ] != 'f')  {printf("missing Red_MEAN!\n");exit(1);	}	
				Tuple< WeightElem<Tuple< double, 6>,2> , 33> profile, profile_c; ExOp::toZero(profile_c); 
				Tuple< WeightElem<Tuple< double, 5>,2> , 11> profile_l, profile_lc; ExOp::toZero(profile_lc); 
				Tuple< WeightElem<double,2> , 11> profile_d, profile_dc; ExOp::toZero(profile_dc);
				//		ExOp::toZero(profiles_b);ExOp::toZero(profiles_m);ExOp::toZero(profiles_d);	ExOp::toZero(profiles_l);
				ExOp::toZero(prior); ExOp::toZero(prior_l); ExOp::toZero(prior_d); 
				//	covar.EMinit();
				mother_buds_input.k.first = to_clust.size();
				loners_input.k.first = to_clust.size();
				daprogbar.start("Reading Cell Attributes: ");
				for(i =0;i < tb.nbrows;i++){
					daprogbar.update(((double)i) / tb.nbrows);
					pix[0] = no_confidence == 2.0f ? tb.getValue(dataindex[1],i).f : (((no_confidence == 0.0f)||(no_confidence <= tb.getValue(dataindex[1],i).f)) ? 1.0f : 0.0f);
					pix[1] = no_confidence == 2.0f ? tb.getValue(dataindex[2],i).f : (((no_confidence == 0.0f)||(no_confidence <= tb.getValue(dataindex[2],i).f)) ? 1.0f : 0.0f);
					if (tb.getValue(dataindex[0],i).sc == 'a') continue;
					if (pix[0] == 0.0f) continue;
					if ((pix[1] == 0.0f)&&(tb.getValue(dataindex[0],i).sc != 'l')) continue;
					size_stage = pow(tb.getValue(dataindex[3],i).f / 1200.0f, 1.5f)  / 1.7f;
					if (tb.getValue(dataindex[0],i).sc == 'l'){
						lone_prob = pix[0]; class_prob = 0.0f;
					}else {
						lone_prob = pix[0] * (1 - pix[1]);class_prob =  pix[0] *  pix[1];
						if (tb.getValue(dataindex[0],i).sc == 'm') cell_stage = pow(tb.getValue(dataindex[4],i).f / 1200.0f, 1.5f);
						else if (tb.getValue(dataindex[0],i).sc == 'b')  cell_stage = pow(tb.getValue(dataindex[3],i).f / 1200.0f, 1.5f);
						else cell_stage = size_stage;
					}
					if ((no_daugther)&&(tb.getValue(dataindex[0],i).sc == 'd')){	lone_prob += class_prob;class_prob =0.0f;	} // DISCARDING DAUGHTER LABEL!
					if (lone_prob > 0.0f){
						supercount[0] += lone_prob;
						k = (unsigned int)(size_stage * resol);
						if (size_stage < 0) {k = 0;coor[0]=0; pix[0] = 0.0f;}
						else if (k >= resol-1) {k = resol-2;pix[0] = 1.0f;coor[0] = gauss_dims[0]-1;}
						else { pix[0] = size_stage * resol - k;coor[0] = (unsigned int)(size_stage * gauss_dims[0]);}
						for(j=0;j<5;j++) data_l[j] = tb.getValue(dataindex[5+j],i).f * factor[j];
						//	printf("l_entrie: ");ExOp::show(data_l);
						profile_lc[k+1] += WeightElem<Tuple<double,5>,2>(data_l, lone_prob*(1.0f-pix[0]));
						profile_lc[k+2] += WeightElem<Tuple<double,5>,2>(data_l, lone_prob*pix[0]);
						profile_lc[0] += WeightElem<Tuple<double,5>,2>(data_l, lone_prob);
						//	profiles_l(coor) += WeightElem<Tuple<double,5>,2>(data_l,lone_prob); 
						prior_l += WeightElem<Tuple<double,5>,2>(data_l,lone_prob);
						if (tb.getValue(dataindex[0],i).sc == 'd') lone_prob += class_prob;
						loners_input.d = data_l;
						loners_input.k.second[0] = lone_prob;
						loners_input.k.second[1] = size_stage;
						loners_input.k.second[2] = tb.getValue(dataindex[14],i).f ;
						loners.push_back(loners_input);
					}
					if (class_prob > 0.0f){
						if (tb.getValue(dataindex[0],i).sc == 'b') supercount[1] += class_prob;
						if (tb.getValue(dataindex[0],i).sc == 'm') supercount[2] += class_prob;
						if (tb.getValue(dataindex[0],i).sc == 'd'){
							supercount[0] += class_prob;
							k = (unsigned int)(size_stage * resol);
							if (size_stage < 0) {k = 0;coor[0]=0; pix[0] = 0.0f;}
							else if (k >= resol-1) {k = resol-2;pix[0] = 1.0f;coor[0] = gauss_dims[0]-1;}
							else { pix[0] = size_stage * resol - k;coor[0] = (unsigned int)(size_stage * gauss_dims[0]);}
							for(j=0;j<5;j++) data_l[j] = tb.getValue(dataindex[5+j],i).f * factor[j];
							profile_lc[k+1] += WeightElem<Tuple<double,5>,2>(data_l,  class_prob*(1.0f-pix[0]));
							profile_lc[k+2] += WeightElem<Tuple<double,5>,2>(data_l,  class_prob*pix[0]);
							profile_lc[0] += WeightElem<Tuple<double,5>,2>(data_l,  class_prob);
							//	profiles_l(coor) += WeightElem<Tuple<double,5>,2>(data_l, class_prob); 
							prior_l += WeightElem<Tuple<double,5>,2>(data_l, class_prob);
							data_d = tb.getValue(dataindex[10],i).f;
							profile_dc[k+1] += WeightElem<double,2>(data_d,  class_prob*(1.0f-pix[0]));
							profile_dc[k+2] += WeightElem<double,2>(data_d,  class_prob*pix[0]);
							profile_dc[0] += WeightElem<double,2>(data_d,  class_prob);
							//	profiles_d(coor) += WeightElem<double>,2>(data_d, class_prob);
							prior_d += WeightElem<double,2>(data_d, class_prob);
							k = (unsigned int)(cell_stage * resol);
							if (cell_stage < 0) {k = 0;coor[0]=0; pix[0] =0.0f;}
							else if (k >= resol-1) {k = resol-2;coor[0] = gauss_dims[0]-1;pix[0] =1.0f;}
							else {coor[0] = (unsigned int)(cell_stage * gauss_dims[0]);pix[0] =(cell_stage * resol) -k;}
							for(j=0;j<6;j++) data[j] = tb.getValue(dataindex[5+j],i).f* factor[j];
							profile_c[resol*2+2] += WeightElem<Tuple<double,6>,2>(data, class_prob);
							profile_c[resol*2+3+k] += WeightElem<Tuple<double,6>,2>(data, class_prob*(1.0f-pix[0]));
							profile_c[resol*2+4+k] += WeightElem<Tuple<double,6>,2>(data, class_prob*pix[0]);//profiles_d(coor) += WeightElem<Tuple<double,6>,2>(data,class_prob); prior[2] += WeightElem<Tuple<double,6>,2>(data,class_prob);
						}else{
							k = (unsigned int)(cell_stage * resol);
							if (cell_stage < 0) {k = 0;coor[0]=0; pix[0] =0.0f;}
							else if (k >= resol-1) {k = resol-2;coor[0] = gauss_dims[0]-1;pix[0] =1.0f;}
							else {coor[0] = (unsigned int)(cell_stage * gauss_dims[0]);pix[0] =(cell_stage * resol) -k;}
							for(j=0;j<6;j++) data[j] = tb.getValue(dataindex[5+j],i).f* factor[j];
							//	printf("entrie: ");ExOp::show(data);
							if (tb.getValue(dataindex[0],i).sc == 'm') {
								profile_c[resol+1] += WeightElem<Tuple<double,6>,2>(data, class_prob);
								profile_c[resol+2+k] += WeightElem<Tuple<double,6>,2>(data, class_prob*(1.0f-pix[0]));
								profile_c[resol+3+k] += WeightElem<Tuple<double,6>,2>(data, class_prob*pix[0]);// profiles_m(coor) += WeightElem<Tuple<double,6>,2>(data,class_prob); prior[1] += WeightElem<Tuple<double,6>,2>(data,class_prob);
								l = dual_info.find(tb.getValue(dataindex[12],i).i);
								if (l == 0xFFFFFFFF) dual_info[tb.getValue(dataindex[13],i).i] = pair<Tuple<double, 6>, double>(data, tb.getValue(dataindex[14],i).f);
								else{
									for(j=0;j<6;j++) {mother_buds_input.d[1 | (j<<1)] = data[j]; mother_buds_input.d[ (j<<1)] = dual_info.deref(l).first[j];}
									mother_buds_input.k.second[0] = class_prob;
									mother_buds_input.k.second[1] = tb.getValue(dataindex[4],i).f;
									mother_buds_input.k.second[2] = tb.getValue(dataindex[3],i).f;
									mother_buds_input.k.second[3] = (tb.getValue(dataindex[14],i).f * tb.getValue(dataindex[3],i).f + dual_info.deref(l).second * tb.getValue(dataindex[4],i).f)  / ( tb.getValue(dataindex[3],i).f + tb.getValue(dataindex[4],i).f);
									mother_buds.push_back(mother_buds_input);
								}
								//		if ((tb.getValue(dataindex[3],i).f < 2000.0f)&& (tb.getValue(dataindex[3],i).f > 1500.0f)){
								//		data_cov[0] = tb.getValue(dataindex[3],i).f;  data_cov[1] = tb.getValue(dataindex[4],i).f; data_cov[2] = tb.getValue(dataindex[11],i).f;  covar.EMregist(data_cov, class_prob);
								//		}
							}else {
								 profile_c[0] += WeightElem<Tuple<double,6>,2>(data, class_prob);
								 profile_c[k+1] += WeightElem<Tuple<double,6>,2>(data, class_prob*(1.0f-pix[0]));
								 profile_c[k+2] += WeightElem<Tuple<double,6>,2>(data, class_prob*pix[0]);//profiles_b(coor) += WeightElem<Tuple<double,6>,2>(data,class_prob); prior[0] += WeightElem<Tuple<double,6>,2>(data,class_prob);
							     
								 
								l = dual_info.find(tb.getValue(dataindex[12],i).i);
								if (l == 0xFFFFFFFF) dual_info[tb.getValue(dataindex[13],i).i] = pair<Tuple<double, 6>, double>(data, tb.getValue(dataindex[14],i).f);
								else{
									for(j=0;j<6;j++) {mother_buds_input.d[ (j<<1)] = data[j]; mother_buds_input.d[1 | (j<<1)] = dual_info.deref(l).first[j];}
									mother_buds_input.k.second[0] = class_prob;
									mother_buds_input.k.second[1] = tb.getValue(dataindex[3],i).f;
									mother_buds_input.k.second[2] = tb.getValue(dataindex[4],i).f;
									mother_buds_input.k.second[3] = (tb.getValue(dataindex[14],i).f * tb.getValue(dataindex[3],i).f + dual_info.deref(l).second * tb.getValue(dataindex[4],i).f)  / ( tb.getValue(dataindex[3],i).f + tb.getValue(dataindex[4],i).f);
									mother_buds.push_back(mother_buds_input);
								}
							 }
						}
					}
				}daprogbar.finish();
				// 475, 800
				dual_info.clear();
				for(k=1;k<22;k+= resol+1) {
					profile[k-1] = profile_c[k-1];
					for(i=0;i<resol;i++) {
						profile[k+i] = profile_c[k+i];
						for(j=0;j<resol;j++) {if (j==i) continue;profile[k+i] +=  profile_c[k+j]  * gaus_pix[10+i-j];}
					}
				}
				k=1;profile_l[0] = profile_lc[0];
				for(i=0;i<resol;i++) {
					profile_l[k+i] = profile_lc[k+i];
					for(j=0;j<resol;j++) {if (j==i) continue; profile_l[k+i] += profile_lc[k+j] * gaus_pix[10+i-j];}
				}
				k=1;profile_d[0] = profile_dc[0];
				for(i=0;i<resol;i++) {
					profile_d[k+i] = profile_dc[k+i];
					for(j=0;j<resol;j++) {if (j==i) continue; profile_d[k+i] += profile_dc[k+j] * gaus_pix[10+i-j];}
				}
				//			if (covar.weight > 10.0f) printf("%e\t%e\n",covar.getCorrelationCoef(0,2),covar.getCorrelationCoef(1,2));
				for(i=0;i<2*(resol+1);i++) {fprintf(g,"\t%f", profile[i].w[0]);}
				for(i=0;i<resol+1;i++) {fprintf(g,"\t%f", profile_l[i].w[0]);}
				for(i=0;i<resol+1;i++) {fprintf(g,"\t%f", profile_d[i].w[0]);}
				for(j=0;;j++){
					for(i=0;i<2*(resol+1);i++) {	data = profile[i].getMean();	fprintf(g,"\t%f", data[j]); 	}
					if (j ==5) break;
					for(i=0;i<resol+1;i++) {	data_l = profile_l[i].getMean();	fprintf(g,"\t%f", data_l[j]);	}
				}
				for(i=0;i<resol+1;i++) {	data_d = profile_d[i].getMean();	fprintf(g,"\t%f", data_d);	}
				for(j=0;;j++){
					for(i=0;i<2*(resol+1);i++) {	data = profile[i].getVar();	fprintf(g,"\t%f", sqrt(data[j]) );}
					if (j ==5) break;
					for(i=0;i<resol+1;i++) {	data_l = profile_l[i].getVar();	fprintf(g,"\t%f", sqrt(data_l[j]));}
				}
				for(i=0;i<resol+1;i++) {data_d = profile_d[i].getVar();	fprintf(g,"\t%f", sqrt(data_d));}
				if ((cluster) || (file_in[2] != NULL) || (hypercl) || (pair_wise)){
					l=0;
					for(j=0;j<6;j++){
						// time
						for(i=0;i<resol;i++,l++) { // buds
							data = profile[i+1].getMean(); clust_in[l].e[0] = data[j];
							//data = profile[i+1].getVar(); clust_in[l].e[1] = data[j] + clust_in[l].e[0]*clust_in[l].e[0];
							data = profile[i+1].getVar_underPrior(prior_for_var); clust_in[l].e[1] = data[j] + clust_in[l].e[0]*clust_in[l].e[0];
						} // 	clust_in[i+resol*2*j] = data[j];	}
						for(i=0;i<resol;i++,l++) { // mothers
							data = profile[i+resol+2].getMean(); clust_in[l].e[0] = data[j];
							//data = profile[i+resol+2].getVar(); clust_in[l].e[1] = data[j] + clust_in[l].e[0]*clust_in[l].e[0];
							data = profile[i+resol+2].getVar_underPrior(prior_for_var); clust_in[l].e[1] = data[j] + clust_in[l].e[0]*clust_in[l].e[0];
						} //	clust_in[i+resol*(1+2*j)] = data[j];	}
						/*
						if (j < 5){
							for(i=0;i<resol;i++,l++) {
								data_l = profile_l[i+1].getMean(); clust_in[l].e[0] = data_l[j];
								//	data = profile[i+resol*2+3].getVar() ; clust_in[l].e[1] = data[j] + clust_in[l].e[0]*clust_in[l].e[0];
								data_l = profile_l[i+1].getVar_underPrior(prior_for_var ); clust_in[l].e[1] = data_l[j] + clust_in[l].e[0]*clust_in[l].e[0];
							}
						}else{
							for(i=0;i<resol;i++,l++) {
								data = profile[i+resol*2+3].getMean(); clust_in[l].e[0] = data[j];
								//	data = profile[i+resol*2+3].getVar() ; clust_in[l].e[1] = data[j] + clust_in[l].e[0]*clust_in[l].e[0];
								data = profile[i+resol*2+3].getVar_underPrior(prior_for_var); clust_in[l].e[1] = data[j] + clust_in[l].e[0]*clust_in[l].e[0];
							}
						} */
					}
					for(i=0;i<clust_in.getSize();i++) {
						clust_in[i].w[0] = 1.0f;clust_in[i].w[1] = 0.0f;
						if (!(ExCo<double>::isValid(clust_in[i].getMean()))) clust_in[i].e[0] =0.0f;
						if (!(ExCo<double>::isValid(clust_in[i].getVar()))) clust_in[i].e[1] =exp(-6.0f);
					}
					to_clust.push_back(clust_in);
					if (locind[locind.size()-1] != 0xFFFFFFFF) {daHGauss[locind[locind.size()-1]] += GaussElem<Tuple<double, 120 > >(clust_in);}
				}
				fprintf(g,"\n");
				//		l++; printf("nbrow: %i\n", tb.nbrows);
				//		if (l == 2) break;
			}
			//	covcov.EMfinit();
		//	if (pair_wise){
		//		FILE* f_pairwise = fopen(pair_wise, "w+");
	/*			WeightElem<double,2> dist[3]; ExOp::toZero(dist);
				for(i=0;i<to_clust.size();i++) {
					ExOp::toZero(dist[1]);
					ExOp::toZero(dist[2]);
//					ExtremumScope<unsigned int, double, false> e_scope; ExOp::toZero(e_scope);
					for(j=0;j<to_clust.size();j++) {
						double sum = to_clust[i][0].getBattDist(to_clust[j][0]);
						for(k=1;k<clust_in.size();k++) sum += to_clust[i][k].getBattDist(to_clust[j][k]);
						//e_scope.regist(locind[j], sum);
						dist[ (locind[j] == locind[i] ) ? 1 : 2] += WeightElem<double,2>(sum);
					}
//					printf("%s\t%s\n", daHclasses_names[locind[i]].c_str(),daHclasses_names[locind[e_scope.best_key]].c_str() );
					if (dist[1].w[0] > 0.0f) dist[0] += WeightElem<double,2>(dist[1].getMean() / dist[2].getMean());
				}*/
			// superm
			DataGrid<double , 2> daBhattMap;
			Tuple<unsigned int ,2 > coors;
			WeightElem<double,2> st_inn[4];
			WeightElem<double,2> st_outt[4];
			WeightElem<double,2> st_inno[4]; ExOp::toZero(st_inno);
			WeightElem<double,2> st_outto[4]; ExOp::toZero(st_outto);
			Vector<unsigned int> damap;
			damap.setSize(daHGauss.size()-2);
			coors[0] = damap.size();
			coors[1] = damap.size();
			daBhattMap.setSizes(coors);
			for(coors[0]=0; coors[0]< daBhattMap.dims[0];coors[0]++)
				for(coors[1]=0; coors[1]< daBhattMap.dims[1];coors[1]++)
					daBhattMap(coors) = daHGauss[coors[0]].bhattacharryya_dist(daHGauss[coors[1]],false);
			for(l=0;l<damap.size();l++) damap[l] = l;
			unsigned int fail = 0;
			double da_dif[3];
			for(l=0;l<1000000;l++){
				ExOp::toZero(st_inn);
				ExOp::toZero(st_outt);
				for(k=0;k<10;k++){
					coors[0] = damap[k];
					for(j=0;j<damap.size();j++){
						if (j==k) continue;
						coors[1] = damap[j];
						if (k < 5){
							if (k < 3){
								if (j < 3) st_inn[0] += WeightElem<double,2>(daBhattMap(coors));
								else  st_outt[0] += WeightElem<double,2>(daBhattMap(coors));
							}else{
								if ((j == 3)||(j == 4)) st_inn[1] += WeightElem<double,2>(daBhattMap(coors));
								else st_outt[1] += WeightElem<double,2>(daBhattMap(coors));
							}
						}else if (k > 6){
							if ((j > 5)&&(j < 10)) st_inn[3] += WeightElem<double,2>(daBhattMap(coors));
							else st_outt[3] += WeightElem<double,2>(daBhattMap(coors));
						}else if (k == 6){
							if ((j == 5)||(j == 6)) st_inn[2] += WeightElem<double,2>(daBhattMap(coors));
							else st_outt[2] += WeightElem<double,2>(daBhattMap(coors));
							if ((j > 5)&&(j < 10)) st_inn[3] += WeightElem<double,2>(daBhattMap(coors));
							else st_outt[3] += WeightElem<double,2>(daBhattMap(coors));
						}else{
							if ((j == 5)||(j == 6)) st_inn[2] += WeightElem<double,2>(daBhattMap(coors));
							else st_outt[2] += WeightElem<double,2>(daBhattMap(coors));							
						}
				}	}
				if (l ==0){
					st_inn[0].show();
					st_inn[1].show();
					st_inn[2].show();
					st_inn[3].show();
					st_outt[0].show();
					st_outt[1].show();
					st_outt[2].show();
					st_outt[3].show();
					da_dif[0] =0;
					for(k=0;k<4;k++) da_dif[0] += st_outt[k].getMean() - st_inn[k].getMean();
				}else{
					da_dif[1] = 0;
					for(k=0;k<4;k++) da_dif[1] += st_outt[k].getMean() - st_inn[k].getMean();
					if (da_dif[1] < da_dif[0]) fail++;
					for(k=0;k<4;k++){
					st_inno[k] +=  WeightElem<double,2>(st_inn[k].getMean());
					st_outto[k] +=  WeightElem<double,2>(st_outt[k].getMean());
					}
				}
				damap.random_permute();
			}
			printf("empirical %i\n", fail);
			st_inno[0].show();
			st_inno[1].show();
			st_inno[2].show();
			st_inno[3].show();
			st_outto[0].show();
			st_outto[1].show();
			st_outto[2].show();
			st_outto[3].show();
			// superm */
		//		fclose(f_pairwise);
		//	}
			if (file_in[2]){
				FILE* d2f = fopen(file_in[2],"w+");
				for(j=0;j<to_clust.size();j++)	{
					fprintf(d2f,"%s\t", to_clust_header[j]);
					for(i=0;i<to_clust.size();i++){
						double sum = to_clust[i][0].getBattDist(to_clust[j][0]);
						for(k=1;k<120;k++) sum += to_clust[i][k].getBattDist(to_clust[j][k]);
						fprintf(d2f,"%f%c", sum, i== to_clust.size() -1 ? '\n':'\t');
					}
				}
				fclose(d2f);
			}
			printf("Bhattacharrya Distance");
			for(i=0;i<daHGauss.size();i++) printf("\t%s", daHclasses_names[i].c_str() );
			for(i=0;i<daHGauss.size();i++){
				printf("\n%s", daHclasses_names[i].c_str() );
				for(j=0;j<daHGauss.size();j++){
					printf("\t%e", daHGauss[i].bhattacharryya_dist(daHGauss[j],false,true) );
				}
			} printf("\n");
			for(i=0;i<daHGauss.size();i++){
				printf("%s\t", daHclasses_names[i].c_str() );
				ExOp::show(daHGauss[i].getMean());
			} printf("\n");
			if (cluster){
				da_clust.cluster_likelihood_ratio(to_clust,batta_cluster_metrix);
				pix[0] = ExCo<double>::maximum(); pix[1] = ExCo<double>::minimum();
				for(i=0;i< da_clust.size();i++){ if (ExCo<double>::isValid(da_clust[i].second)) {
					if (da_clust[i].second <=0.0f) continue;
					if (da_clust[i].second < pix[0]) pix[0] = da_clust[i].second;
					if (da_clust[i].second > pix[1]) pix[1] = da_clust[i].second;
				}}
				for(i=0;i< da_clust.size();i++) da_clust[i].second = (ExCo<double>::isValid(da_clust[i].second)) ? ((da_clust[i].second <=0.0f) ? 1.0f: (log(pix[1]) - log(da_clust[i].second)) / (log(pix[1]) - log(pix[0]) ) ) : 0.0f;
				da_clust.saveGTRfile("auto_cluster.gtr", "GENE");
				da_clust.saveCDTfile_W(clust_cdt, to_clust_header, to_clust);
				fclose(clust_cdt);
			}
		/*	 if (daHGauss.size() != 0){
				Tuple<unsigned int, 0u> confusion; confusion.setSize(daHGauss.size()*daHGauss.size());ExOp::toZero(confusion);
				if (hypercl_near){
					for(i=0;i< to_clust.size();i++){
						if ((locind[i] == 0xFFFFFFFF)||(locind[i] == 0)) continue;
						GaussElem<Tuple<double, 120 > > tmp(to_clust[i]);
						ExtremumScope<unsigned int, double, false> e_scope; ExOp::toZero(e_scope);
						for(j=0;j< to_clust.size();j++){
							if ((locind[j] == 0xFFFFFFFF)||(i==j)||(locind[j] == 0)) continue;
							GaussElem<Tuple<double, 120 > > tmp2(to_clust[j]);
							e_scope.regist(locind[j], tmp.likelihoodratio_dist( tmp2,false,false) );
						}
						printf("%i ->%i\n",locind[i],e_scope.best_key);
						confusion[locind[i]*daHGauss.size() + e_scope.best_key]++;
					}
				}else{
				for(i=0;i<daHGauss.size();i++){
					printf("%i", i );
					for(j=0;j<daHGauss.size();j++){
						printf("\t%e", daHGauss[i].bhattacharryya_dist(daHGauss[j],false) / (daHGauss[i].w + daHGauss[j].w) );
					}
					printf("\n");
				}
				for(i=0;i<daHGauss.size();i++){
					printf("%i\t%f", i, daHGauss[i].w);
					for(j=0;j<120;j++){
						clust_simple_in[j] = daHGauss[i].mean[j] / daHGauss[i].w;
						clust_simple_in[j] = daHGauss[i].cov.getDiagonal()[j] - clust_simple_in[j] * clust_simple_in[j];
						printf("\t%e", daHGauss[i].mean[j] / daHGauss[i].w);
					}
					printf("\t");
					ExOp::show(clust_simple_in);
				}
				double best;
				for(i=0;i< to_clust.size();i++){
					if (locind[i] == 0xFFFFFFFF) continue;
					GaussElem<Tuple<double, 120 > > tmp(to_clust[i]);
					GaussElem<Tuple<double, 120 > > back = daHGauss[locind[i]];
					daHGauss[locind[i]] -= tmp;
					ExtremumScope<unsigned int, double, false> e_scope;
					tmpdouble = tmp.likelihoodratio_dist(daHGauss[1],false,true); printf("%i\t%e", i,tmpdouble);
					e_scope.init(1u,tmpdouble);
					for(j=2;j<daHGauss.size();j++) {tmpdouble = tmp.likelihoodratio_dist(daHGauss[j],false,true);e_scope.regist(j,tmpdouble);	printf("\t%e", tmpdouble);}			
					confusion[locind[i]*daHGauss.size() + e_scope.best_key]++;
					daHGauss[locind[i]] = back;
					printf("%i is best!\n",e_scope.best_key);
				}
				}
				FILE *f_out_classif = (hypercl_out) == NULL ? stdout : fopen(hypercl_out, "w+");
				for(i=0,k=0;i<daHGauss.size();i++){
					for(j=0;j<daHGauss.size()-1;j++){
						fprintf(f_out_classif,"%i\t",confusion[k++]);
					}
					fprintf(f_out_classif,"%i\n",confusion[k++]);
				}
				fclose(f_out_classif);
			}*/
			if (file_in[1]) fclose(g);
		}
	return 0;}
	void Taskscope<TASK_PROFILES>::help(){
		printf("Make a Condifence Profile. This version explicitly ignore cell-to-cell covariances in feature measures\n");
		printf("Default Arguments 1-2:\n");
		printf("(in) file of table paths\n");
		printf("(out) Tabular file (stdout otherwise)\n");
		printf("\n");
		printf("Optionnal Default arguments(3):\n");
		printf("-p: change path variables (auto-on if nbargs < 2)! If corrupted, del PMProfiles_scope.scp\n");
		printf("-s: Only show the content of the scopefile (no modifications), automatically on no tabular file provided.\n");
		printf("-a: Make a Displayable annotated image\n");	
		printf("-D (file) : 2D distance map\n");	
		printf("-C: cluster profiles\n");
		printf("-i (float t): ignore confidence measure, use threshold t on confidence to accept!\n");
		printf("-d: discard daughter label, they become lone cells\n");	
		printf("-L(N) (file in) {file out = stdout}: classification! file of (plate_ID label) pairs, (add 'N' for nearest neighbor classification)\n");
		printf("-B (float sep) (float plusminus): bin approach! \n");
		printf("-J (file) : Produces a input file for SVM classification by Joachims's SVMlite program\n");
		printf("-S (float) (float) (float) (float) (float) (float) : scale factors for output distances\n");
		printf("\t-h : Help \n");
	}
