/************************************************************************
Copyright [2021] [Christian B. Huebschle chuebsch@moliso.de]

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
************************************************************************/
#ifndef AUTOHFIX_H
#define AUTOHFIX_H
#define LM 2000000
#include "kissfft/kiss_fftnd.h"
#include "molecule.h"
#if (QT_VERSION >= 0x050e00)
#define skipEmptyParts Qt::SkipEmptyParts
#else
#define skipEmptyParts QString::SkipEmptyParts
#endif

//! Rec is a reflection type of a fcf 6 file.
typedef struct {
  int  ih,//!< h
       ik,//!< k
       il;//!< l
  float fo,//!< F observed
        so,//!< \f$\sigma(observed)\f$
        fc,//!< F calculated
        phi; //!< \f$\varphi\f$
} Rec;
int HKLMX=200;
float sigma[3];//!<sigma values
int n1;//!< dimension of the map in a diRection
int n2;//!< dimension of the map in b diRection
int n3;//!< dimension of the map in c diRection
int n4;//!< \f$ n4 = n1 \times n2\f$
int n5;//!< \f$ n5 = n1 \times n2 \times n3\f$

V3 dx,//!< vector in a direction for each map voxel
   dy,//!< vector in b direction for each map voxel
   dz;//!< vector in c direction for each map voxel
double C[15],D[9],sy[12][192],wave;
Rec lr[LM];

char cen,git;
int nr,nc,ns;
char titl[80];/*fcmax=0,f000=0,resmax=99999.0,*/
float *datfo_fc;//!<data pointer for the fobs-fcalc map in real space
float f000;
/*V3 dx,//!< vector in a direction for each map voxel
   dy,//!< vector in b direction for each map voxel
   dz;//!< vector in c direction for each map voxel
   */

#ifdef FFTW3_H
    fftwf_plan  fwd_plan;//!!!
    fftwf_complex *B;//!!!
#else
    kiss_fftnd_cfg fwd_plan;//!!!
    kiss_fft_cpx *B;//!!!
#endif
    Molecule mole;
    QList<int> sfac;//!<List of Scattering factors.
    QStringList resLines;
    QList<double> fvar;//!<List of Free Variables.
    QMap<int,int> fvarCntr;//!<Free Variable counter QMap.
#endif // AUTOHFIX_H
