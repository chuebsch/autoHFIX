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
#include <QCoreApplication>

#include <QtCore>
#include "autohfix.h"
#include "scatt.h"
#include <math.h>
void trimm(char s[]){
  /*! a trimm function for c-strings.
  */
  char sc[409];
  size_t j=0;
  size_t len=strlen(s);
  strncpy(sc,s,400);
  for (size_t i=0; i<len; i++) if ((sc[i]!='\'')&&(!isspace(sc[i]))) s[j++]=static_cast<char>(toupper(sc[i]));
  s[j]='\0';
}

int readHeader(const char *filename){
  /*! reads the header of an fcf file
  */
  FILE *f=NULL;
  char line[122],*dum=NULL;
  //size_t zlen=120;
  int ok=0;
  int i;double T,V;
  f=fopen(filename,"r");
  if (f==NULL) return 3;
  ns=0;
  sy[0][ns]=1.0;
  sy[1][ns]=0.0;
  sy[2][ns]=0.0;

  sy[3][ns]=0.0;
  sy[4][ns]=1.0;
  sy[5][ns]=0.0;

  sy[6][ns]=0.0;
  sy[7][ns]=0.0;
  sy[8][ns]=1.0;

  sy[9][ns]=0.0;
  sy[10][ns]=0.0;
  sy[11][ns]=0.0;
  ns=1;
  int listcode=0;
  do{
    dum=fgets(line,120,f);
    //printf("%d\n",__LINE__);
    if (dum==NULL){fclose(f);return 2;}
    if (feof(f)){fclose(f);return 2;}
    while (dum[0]==' ') dum++;
    if (!strncmp(dum,"_shelx_title",12)) {
      sscanf(line,"_shelx_title %[^\r\n]",titl);
      trimm(titl);
    }
    if (!strncmp(dum,"_exptl_crystal_F_000",20)){
      sscanf(line,"_exptl_crystal_F_000 %f",&f000);
    }
    if (!strncmp(dum,"_shelx_refln_list_code",22)) {
      sscanf(line,"_shelx_refln_list_code %d",&listcode);
      //qDebug()<<listcode;
      if (listcode!=6) {fclose(f);return 1;}
    }
    if (!strncmp(dum,"_cell_length_a",14)) {
      sscanf(line,"_cell_length_a %lf",&C[0]);
    }
    if (!strncmp(dum,"_cell_length_b",14)) {
      sscanf(line,"_cell_length_b %lf",&C[1]);
    }
    if (!strncmp(dum,"_cell_length_c",14)) {
      sscanf(line,"_cell_length_c %lf",&C[2]);
    }
    if (!strncmp(dum,"_cell_angle_alpha",17)) {
      sscanf(line,"_cell_angle_alpha %lf",&C[3]);
    }
    if (!strncmp(dum,"_cell_angle_beta",16)) {
      sscanf(line,"_cell_angle_beta %lf",&C[4]);
    }
    if (!strncmp(dum,"_cell_angle_gamma",17)) {
      sscanf(line,"_cell_angle_gamma %lf",&C[5]);
      for (i=0;i<3;i++){
        if (C[i]<0.1) return 2;
        T=.0174533*C[i+3];
        if (T<0.001) return 2;
        D[i]=sin(T);
        D[i+3]=cos(T);
        C[i+6]=(D[i]/(C[i]*C[i]));
      }
      V=1.-D[3]*D[3]-D[4]*D[4]-D[5]*D[5]+2.*D[3]*D[4]*D[5];
      C[6]/=V;
      C[7]/=V;
      C[8]/=V;
      C[9]= 2.*sqrt(C[7]*C[8])*(D[4]*D[5]-D[3])/(D[2]*D[2]);
      C[10]=2.*sqrt(C[6]*C[8])*(D[3]*D[5]-D[4])/(D[0]*D[2]);
      C[11]=2.*sqrt(C[6]*C[7])*(D[3]*D[4]-D[5])/(D[0]*D[1]);
      C[14]=C[1]*C[2]*C[0]*sqrt(V);
      D[6]=C[1]*C[2]*D[0]/C[14];
      D[7]=C[0]*C[2]*D[1]/C[14];
      D[8]=C[0]*C[1]*D[2]/C[14];
    }
    if ((!strncmp(dum,"_symmetry_equiv_pos_as_xyz",26))||(!strncmp(dum,"_space_group_symop_operation_xyz",32))){
      dum=fgets(line,120,f);
      trimm(line);
      while (strchr(line,'Y')) {
        QString sc=QString(line).toUpper().remove("SYMM").trimmed();
        sc=sc.remove("'");
        sc=sc.remove(" ");
        QStringList axe=sc.split(",");
        QStringList bruch;
        if (axe.size()!=3) return false;
        double _sx[3],_sy[3],_sz[3],t[3];
        for (int i=0; i<3; i++){
          _sx[i]=0;_sy[i]=0;_sz[i]=0;t[i]=0;
          if (axe.at(i).contains("-X")) {_sx[i]=-1.0;axe[i].remove("-X");}
          else if (axe.at(i).contains("X")) {_sx[i]=1.0;axe[i].remove("X");}
          if (axe.at(i).contains("-Y")) {_sy[i]=-1.0;axe[i].remove("-Y");}
          else if (axe.at(i).contains("Y")) {_sy[i]=1.0;axe[i].remove("Y");}
          if (axe.at(i).contains("-Z")) {_sz[i]=-1.0;axe[i].remove("-Z");}
          else if (axe.at(i).contains("Z")) {_sz[i]=1.0;axe[i].remove("Z");}
          if (axe.at(i).endsWith("+")) axe[i].remove("+");
          if (axe.at(i).contains("/")) {
            bruch=axe.at(i).split("/");
            if (bruch.size()==2) t[i]=bruch.at(0).toDouble() / bruch.at(1).toDouble();
          }
          else if (!axe.at(i).isEmpty()) t[i]=axe.at(i).toDouble();
        }
        sy[0][ns]=_sx[0];
        sy[1][ns]=_sy[0];
        sy[2][ns]=_sz[0];

        sy[3][ns]=_sx[1];
        sy[4][ns]=_sy[1];
        sy[5][ns]=_sz[1];

        sy[6][ns]=_sx[2];
        sy[7][ns]=_sy[2];
        sy[8][ns]=_sz[2];


        sy[9][ns]=t[0];
        sy[10][ns]=t[1];
        sy[11][ns]=t[2];


        strcpy(line,"");
        dum=fgets(line,120,f);
        trimm(line);
        ns++;

      }
    }
    if (!strncmp(dum,"_refln_phase_calc",17)) ok=1;
  }while((!ok)&&(!feof(f)));

  if (listcode!=6) return 1;
  for (int i=0; i<ns;i++){
    for (int n=i+1; n<ns;n++){
      int u=0,v=0;
      for (int j=0; j<9; j++){
        u+=abs(sy[j][n]-sy[j][i]);
        v+=abs(sy[j][n]+sy[j][i]);
      }
      if (fmin(u,v)>0.01) continue;
      for (int j=0; j<12; j++){
        sy[j][n]=sy[j][ns-1];
      }
      ns--;
    }
  }
  fclose(f);
  return 0;
}

void sorthkl(int nr, Rec r[]){
  /*! sorts the reflection list
  */
  Rec *hilf= (Rec*) malloc(sizeof(Rec)*nr);
  if (hilf==NULL)return ;
  int i,j,k,nj,ni,spalte;int index[4096];
  for (spalte=0; spalte<3; spalte++){
    j=-999999;
    k=999999;
    switch (spalte) {
      case 0: for (i=0; i<nr; i++){ j=(j<r[i].ih)?r[i].ih:j; k=(k>r[i].ih)?r[i].ih:k; } break;
      case 1: for (i=0; i<nr; i++){ j=(j<r[i].ik)?r[i].ik:j; k=(k>r[i].ik)?r[i].ik:k; } break;
      case 2: for (i=0; i<nr; i++){ j=(j<r[i].il)?r[i].il:j; k=(k>r[i].il)?r[i].il:k; } break;
    }
    nj=-k;
    ni=(nj+j+1);
    for (i=0; i<=ni; i++) index[i]=0;
    for (i=0; i<nr; i++){
      switch (spalte){
        case 0: j=r[i].ih+nj; break;
        case 1: j=r[i].ik+nj; break;
        case 2: j=r[i].il+nj; break;
      }
      index[j]++;/*brauch ich das? -->JA!*/
      hilf[i].ih=r[i].ih;
      hilf[i].ik=r[i].ik;
      hilf[i].il=r[i].il;
      hilf[i].fo=r[i].fo;
      hilf[i].so=r[i].so;
      hilf[i].fc=r[i].fc;
      hilf[i].phi=r[i].phi;
    }/*/4*/
    j=0;
    for (i=0; i<ni; i++){
      k=j;
      j+=index[i];
      index[i]=k;
    }/*/5*/
    for (i=0; i<nr;i++){
      switch (spalte) {
        case 0: j=hilf[i].ih +nj;break;
        case 1: j=hilf[i].ik +nj;break;
        case 2: j=hilf[i].il +nj;break;
      }
      index[j]++;
      j=index[j]-1;
      r[j].ih=hilf[i].ih;
      r[j].ik=hilf[i].ik;
      r[j].il=hilf[i].il;
      r[j].fo=hilf[i].fo;
      r[j].so=hilf[i].so;
      r[j].fc=hilf[i].fc;
      r[j].phi=hilf[i].phi;
    }/*/6*/
  }/*/spalten*/
  free(hilf);
}

void loadFouAndPerform(const char filename[]){
  /*! loads a fcf file prepares the reciprocal data for the fourier transpormation and performs it.
  */
  const int it[143]= {2,3,4,5,6,8,9,10,12,15,16,18,20,24,25,27,30,32,36,40,45,48,50,54,60,64,72,75,80,81,90,96,100,
    108,120,125,128,135,144,150,160,162,180,192,200,216,225,240,243,250,256,270,288,300,320,324,360,375,384,400,405,
    432,450,480,486,500,512,540,576,600,625,640,648,675,720,729,750,768,800,810,864,900,960,972,1000,1024,1080,1125,
    1152,1200,1215,1250,1280,1296,1350,1440,1458,1500,1536,1600,1620,1728,1800,1875,1920,1944,2000,2025,2048,2160,
    2187,2250,2304,2400,2430,2500,2560,2592,2700,2880,2916,3000,3072,3125,3200,3240,3375,3456,3600,3645,3750,3840,
    3888,4000,4050,4096,4320,4374,4500,4608,4800,4860,5000};//!multiples of 2 3 and 5
  int ok;
  if (strstr(filename,".fcf")==NULL) return;
  FILE *f;
  double rr=6.0;
  ok= readHeader(filename);
  if (ok) {
    switch (ok){
      case 1: fprintf(stderr, "Map generation failed. SHELXL LIST code was not 6.\n");break;
      case 2: fprintf(stderr, "Map generation failed. File %s corrupted.\n",filename);break;
      case 3: fprintf(stderr, "Map generation failed. Cannot open file %s.\n",filename);break;
      case 4: fprintf(stderr, "Map generation failed. No reflection data in file %s.\n",filename);break;
    }
    return;
  }
  f=fopen(filename,"rb");
  if (f==NULL)return;
  char line[122]="";
  while (strstr(line,"_refln_phase_calc")==NULL) {
    if (fgets(line,120,f)){}//just a stupid supression for annoying warnig
  }
  nr=0;
  lr[nr].ih=0;
  lr[nr].ik=0;
  lr[nr].il=0;
  lr[nr].fo=f000*f000;
  lr[nr].so=0.5f;
  lr[nr].fc=f000;
  lr[nr].phi=0.0f;
  nr=1;
  int hmin=0, hmax=0;
  int kmin=0, kmax=0;
  int lmin=0, lmax=0;
  fprintf(stderr,"%s\n",filename);

  while (!feof(f)){
    if (fgets(line,120,f)){}//just a stupid supression for annoying warnig
    int rdi=
      sscanf(line,"%d %d %d %f %f %f %f",&lr[nr].ih,&lr[nr].ik, &lr[nr].il ,&lr[nr].fo, &lr[nr].so, &lr[nr].fc, &lr[nr].phi);
    if (rdi==7) {        
      if ((abs(lr[nr].ih)<HKLMX)&&
          (abs(lr[nr].ik)<HKLMX)&&
          (abs(lr[nr].il)<HKLMX)&&
          ((lr[nr].ih|lr[nr].ik|lr[nr].il)!=0)){
          hmin=lr[nr].ih<hmin?lr[nr].ih:hmin;
          hmax=lr[nr].ih>hmax?lr[nr].ih:hmax;

          kmin=lr[nr].ik<kmin?lr[nr].ik:kmin;
          kmax=lr[nr].ik>kmax?lr[nr].ik:kmax;


          lmin=lr[nr].il<lmin?lr[nr].il:lmin;
          lmax=lr[nr].il>lmax?lr[nr].il:lmax;
        nr++;
      }
    }
    if (nr>=LM) {
      fprintf(stderr,"to many reflections in fcf file\n");
      return;
    }
  }
  fclose(f);
  if (nr<2) {
    fprintf(stderr,"Map generation failed. No reflection data in file %s.",filename);
    return;
  }
  printf("%d<h<%d  %d<h<%d  %d<h<%d\n",hmin,hmax,kmin,kmax,lmin,lmax);
  for (int i=0;i<nr;i++){
    double u=lr[i].ih,v=lr[i].ik,w=lr[i].il;
    int mh=lr[i].ih,mk=lr[i].ik,ml=lr[i].il;
    double p,q=lr[i].phi/180.0*M_PI;
    lr[i].phi=fmod(4*M_PI+q,2*M_PI);
    for (int k=0; k<ns; k++){
      int nh,nk,nl;
      double t=1.0;
      nh=(int) (u*sy[0][k]+ v*sy[3][k] + w*sy[6][k]);
      nk=(int) (u*sy[1][k]+ v*sy[4][k] + w*sy[7][k]);
      nl=(int) (u*sy[2][k]+ v*sy[5][k] + w*sy[8][k]);
      if((nl<0)||((nl==0)&&(nk<0))||((nl==0)&&(nk==0)&&(nh<0)))
      {nh*=-1;nk*=-1;nl*=-1;t=-1.0;}
      if ((nl<ml)||((nl==ml)&&(nk<mk))||((nl==ml)&&(nk==mk)&&(nh<=mh))) continue;
      mh=nh;mk=nk;ml=nl;
      p=u*sy[9][k]+v*sy[10][k]+w*sy[11][k];
      lr[i].phi=fmod(4*M_PI+t*fmod(q-2*M_PI*p,2*M_PI)-0.01,2*M_PI)+0.01;

    }
    lr[i].ih=mh;
    lr[i].ik=mk;
    lr[i].il=ml;
  }
  sorthkl(nr,lr);
  int n=-1;
  {int i=0;
    while(i<nr){
      double t=0.;
      double u=0.;
      double v=0.;
      double z=0.;
      double p=0.;
      int m;
      int k=i;
      while ((i<nr)&&(lr[i].ih==lr[k].ih)&&(lr[i].ik==lr[k].ik)&&(lr[i].il==lr[k].il)) {
        t=t+1.;
        u+=lr[i].fo;
        v+=1./(lr[i].so*lr[i].so);
        z+=lr[i].fc;
        p=lr[i].phi;
        i++;
      }
      m=n+1;
      lr[m].fo=sqrt(fmax(0.,u/t));
      lr[m].so=sqrt(lr[m].fo*lr[m].fo+sqrt(1./v))-lr[m].fo;
      lr[m].fc=z/t;
      lr[m].phi=p;
      n=m;
      lr[n].ih=lr[k].ih;
      lr[n].ik=lr[k].ik;
      lr[n].il=lr[k].il;
    }
  }
  n++;
  nr=n;
  //  printf("%4d%4d%4d %g %g %g %g\n",lr[0].ih,lr[0].ik,lr[0].il,lr[0].fo,lr[0].so,lr[0].fc,lr[0].phi);
  {
    float DX;
    float DY;
    float DZ;
    {
      int mh=0, mk=0, ml=0,j;
      for (int n=0; n<nr; n++){
        double u=lr[n].ih,v=lr[n].ik,w=lr[n].il;
        int a,b,c;
        for (int k=0; k<ns;k++){
          a=abs((int)(u*sy[0][k]+v*sy[3][k]+w*sy[6][k]));
          b=abs((int)(u*sy[1][k]+v*sy[4][k]+w*sy[7][k]));
          c=abs((int)(u*sy[2][k]+v*sy[5][k]+w*sy[8][k]));
          mh=(mh<a)?a:mh;
          mk=(mk<b)?b:mk;
          ml=(ml<c)?c:ml;
        }
      }
      j=(int)(rr*mh+.5);
      for (int i=0; it[i]< j; i++)n1=it[i+1];
      j=(int)(rr*mk+.5);
      for (int i=0; it[i]< j; i++)n2=it[i+1];
      j=(int)(rr*ml+.5);
      for (int i=0; (it[i]< j)||((nc)&&(it[i]%2)); i++) n3=it[i+1];
      n4=n2*n1;
      n5=n3*n4;
      datfo_fc=(float*) malloc(sizeof(float)*n5);
      if (datfo_fc==NULL){
          fprintf(stderr,"Too less memory for F_observed - F_calculated!\n");
          datfo_fc=NULL;
          fprintf(stderr,"Could not allocate sufficient memory for density maps. Defaults have been applied. \n");
          return;
      }
      DX=1.0f/n1;
      DY=1.0f/n2;
      DZ=1.0f/n3;
    }
    //for (int typ=0; typ<2;typ++)
    {
      double miZ=99999.99,maZ=-99999.99;
      int nbytes,dims[3];
      dims[0]=n3;
      dims[1]=n2;
      dims[2]=n1;
#ifdef FFTW3_H
      B=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*n5);
      if (B==NULL)return false;
      for (int i=0; i<n5; i++){B[i][0]=0;B[i][1]=0;}
#else
      B=(kiss_fft_cpx*)KISS_FFT_MALLOC(nbytes = (sizeof(kiss_fft_cpx)*n5));
      if (B==NULL){
          fprintf(stderr,"Too less memory for FFT map!\n");
          free(datfo_fc);
          datfo_fc=NULL;
          fprintf(stderr,"Could not allocate sufficient memory for density maps. Defaults have been applied. \n");
          return;
      }
      for (int i=0; i<n5; i++){B[i].r=0;B[i].i=0;}
#endif
      for (int i=0; i<nr;i++){
        float  u,v,w;
        u=lr[i].ih;
        v=lr[i].ik;
        w=lr[i].il;
        float  ss,s=0,t=0,q,p;
        for (int n=0; n<ns;n++){
          int j,k,l;
          j=(int) (u*sy[0][n]+ v*sy[3][n] + w*sy[6][n]);
          k=(int) (u*sy[1][n]+ v*sy[4][n] + w*sy[7][n]);
          l=(int) (u*sy[2][n]+ v*sy[5][n] + w*sy[8][n]);
          if((abs(j-lr[i].ih)+abs(k-lr[i].ik)+abs(l-lr[i].il))==0)s+=1.0f;
          if(abs(j+lr[i].ih)+abs(k+lr[i].ik)+abs(l+lr[i].il)==0)t+=1.0f;
        }
        if (i==0) {//printf("v%f s%f t%f\n",C[14],s,t);
          s=1;t=0;//f000
        }
        ss=(lr[i].fo-lr[i].fc)/(C[14]*(s+t));
        for (int n=0; n<ns;n++){
          int j,k,l,m;
          j=(int) (u*sy[0][n]+ v*sy[3][n] + w*sy[6][n]);
          k=(int) (u*sy[1][n]+ v*sy[4][n] + w*sy[7][n]);
          l=(int) (u*sy[2][n]+ v*sy[5][n] + w*sy[8][n]);
          //          q=(-2*M_PI*(u*sy[9][n]+v*sy[10][n]+w*sy[11][n]))-M_PI*(j*DX+k*DY+l*DZ);
          q=(lr[i].phi-2*M_PIf*(u*sy[9][n]+v*sy[10][n]+w*sy[11][n]))-M_PI*(j*DX+k*DY+l*DZ);
          j=(999*n1+j)%n1;
          k=(999*n2+k)%n2;
          l=(999*n3+l)%n3;
          m=j+n1*(k+n2*l);
          p=ss*cosf(q);
          q=ss*sinf(q);
#ifdef FFTW3_H
          B[m][0]=p;
          B[m][1]=q;
#else
          B[m].r=p;
          B[m].i=q;
#endif
          j*=-1;
          if(j<0)j=n1+j;
          k*=-1;
          if(k<0)k=n2+k;
          l*=-1;
          if(l<0)l=n3+l;
          m=j+n1*(k+n2*l);
#ifdef FFTW3_H
          B[m][0]=p;
          B[m][1]=-q;
#else
          B[m].r=p;
          B[m].i=-q;
#endif
        }
      }
#ifdef FFTW3_H
      fwd_plan = fftwf_plan_dft_3d(n3,n2,n1,B,B,FFTW_FORWARD,FFTW_ESTIMATE);
      fftwf_execute(fwd_plan);
      fftwf_destroy_plan(fwd_plan);
#else
      fwd_plan = kiss_fftnd_alloc(dims,3,0,0,0);
      kiss_fftnd( fwd_plan,B,B);
      free(fwd_plan);
#endif
      float t=0;
      double DM=0.,  DS=0., DD  ;
      for (int i=0; i<n5;i++){
#ifdef FFTW3_H
        DD=B[i][0];
#else
        DD=B[i].r;
#endif
        miZ=fmin(miZ,DD);
        maZ=fmax(maZ,DD);
        DM+=DD;
        DS+=DD*DD;
#ifdef FFTW3_H
        datfo_fc[i]=B[i][0];
#else
        datfo_fc[i]=B[i].r;
#endif
      }
      sigma[0]=t=sqrt((DS/n5)-((DM/n5)*(DM/n5)));
      printf("sigma %f max %f min %f\n",t,maZ,miZ);
      free(B);
    }//1
  }//2
  dx=V3(1.0/(n1),0,0);
  dy=V3(0,1.0/(n2),0);
  dz=V3(0,0,1.0/(n3));
  mole.frac2kart(dx,dx);
  mole.frac2kart(dy,dy);
  mole.frac2kart(dz,dz);
  fprintf(stderr,"Uniq Reflections: %d\nFourier grid dimensions: %d X %d X %d = %d grid points.\n",nr,n1,n2,n3,n5);
}

void addNewScatteringFactor(int oz){
  if (sfac.contains(oz)) {
    return;
  }
  sfac.append(oz);
  int sfacline = resLines.lastIndexOf(QRegExp("^sfac ",Qt::CaseInsensitive));
  QString sf = resLines.at(sfacline).trimmed();
  if (sf.size()<77) sf.append(QString(" %1").arg(mole.pse(oz)));
  else sf.append(QString("\nSFAC %1").arg(mole.pse(oz)));
  resLines[sfacline]=sf;
  int unitline = resLines.lastIndexOf(QRegExp("^unit ",Qt::CaseInsensitive));
  QString un = resLines.at(unitline).trimmed();
  if (un.size()<77) un.append(" 1");
  else un.append(QString("= \n 1"));
  resLines[unitline]=un;
}

MyAtom findOH(V3 donor, V3 acceptor,int dindex,QStringList alab){
  CEnvironment peaks;
  QString nn;
  peaks.clear();
  double minGoofy=100000.0,g=0;
  int minG=0;
  V3 mdz=(0.5*dx)+(0.5*dy)+(0.5*dz);
  V3 df=donor;
  V3 af=acceptor;
  if (mole.fl(af.x-df.x, af.y-df.y,af.z-df.z)>3.3){
    V3 floors;
    floors = af-df+ V3(0.5,0.5,0.5);
    floors = V3(floor(floors.x),floor(floors.y),floor(floors.z));
    af+=-1.0*floors;
  }
  double ad=mole.fl(af.x-df.x, af.y-df.y,af.z-df.z);
  MyAtom htest;
  htest.an=0;
  htest.part=0;
  htest.resiNr=0;
  htest.hidden=0;
  htest.Label=mole.asymm.at(dindex).Label.section('_',0,0);
  htest.part=mole.asymm.at(dindex).part;
  htest.resiNr=mole.asymm.at(dindex).resiNr;
  htest.Label.replace(0,1,'H');
  int hindex=0;
  QString nnn=htest.Label;
  do {
    if (hindex) nnn=htest.Label.left(3).append(QString::number(hindex+9,36));
    if (nnn.size()>4) nnn=htest.Label.left(2).append(QString::number(hindex+9,36));
    if (htest.resiNr)
      nn=QString("%1_%2")//for shelxl same names in different parts are no problem but for platon so this should fix it.
        .arg(nnn)
        .arg((htest.resiNr)?QString::number(htest.resiNr):"");
    else nn=nnn;
    hindex++;
  }while(alab.contains(nn,Qt::CaseInsensitive));
  htest.Label=nnn;
  alab.append(nn);
  htest.symmGroup=0;
  htest.sg=0;
  htest.scod=555;//the identity
  htest.auidx=-1;
  htest.frac=V3(0,0,0);
  htest.pos=V3(0,0,0);
  htest.isIso=true;
  htest.ufiso_org=QString("-1.5");
  double uIso = mole.ueq(mole.asymm.at(dindex).uf)*1.5;
  htest.uf.m11 = htest.uf.m22 = htest.uf.m33 = uIso;
  htest.uf.m23 = uIso * mole.cell.cosra;
  htest.uf.m13 = uIso * mole.cell.cosrb;
  htest.uf.m12 = uIso * mole.cell.cosrg;
  mole.Uf2Uo(htest.uf,htest.uc);
  if (n5==0) return htest;
  double x,y,z,u,v,w,p,q,r,t;
  for (int k=0;k<n3;k++){
    for (int j=0;j<n2;j++){
      for (int i=0;i<n1;i++){
        V3 o= (i*dx+j*dy+k*dz) + mdz;
        mole.kart2frac(o,o);
        V3 floors;
        floors = o-df+ V3(0.5,0.5,0.5);
        floors = V3(floor(floors.x),floor(floors.y),floor(floors.z));
        o+=-1.0*floors;
        double ha,hd;
        hd=mole.fl(df.x-o.x, df.y-o.y,df.z-o.z);
        ha=mole.fl(af.x-o.x, af.y-o.y,af.z-o.z);
        if (1.2<hd) continue;//donor H must be smaller 1.2 A
        if (2.6<ha) continue;//donor H must be smaller 2.6 A
        if (hd<0.6) continue;//donor H must be greater 0.6 A
        if (((ha*ha+hd*hd-ad*ad)/(2*hd*ha))>-0.57)continue; //angle greater 125deg ok!
        t = datfo_fc[i+ n1*(j+k*n2)];
        x = t -  datfo_fc[((i+1)%n1) + n1*(j+k*n2)];
        y = t -  datfo_fc[i+ n1*(((j+1)%n2)+k*n2)];
        z = t -  datfo_fc[i+ n1*(j+((k+1)%n3)*n2)];
        u = t -  datfo_fc[((i+n1-1)%n1)+ n1*(j+k*n2)];
        v = t -  datfo_fc[i+ n1*(((j+n2-1)%n2)+k*n2)];
        w = t -  datfo_fc[i+ n1*(j+((k+n2-1)%n3)*n2)];
        if (fmin(x,fmin(y,fmin(z,fmin(u,fmin(v,w)))))<0) continue;//if this is not a maximum then skip
        p=0.5*(u-x)/(u+x);
        q=0.5*(v-y)/(v+y);
        r=0.5*(w-z)/(w+z);
        t=t+0.25*(p*(u-x)+q*(v-y)+r*(w-z));
        if (t<fabs(sigma[0])*3.3) continue;// if this is below the green iso surface then skip
        htest.peakHeight=t;
        htest.frac.x = ((double) i + p +0.5)/n1;
        htest.frac.y = ((double) j + q +0.5)/n2;
        htest.frac.z = ((double) k + r +0.5)/n3;
        hd=mole.fl(df.x-htest.frac.x, df.y-htest.frac.y,df.z-htest.frac.z);
        ha=mole.fl(af.x-htest.frac.x, af.y-htest.frac.y,af.z-htest.frac.z);
        if (hd<0.6) continue;//to close to donor
        g=(ha+hd)/htest.peakHeight;
        minG=(g<minGoofy)?peaks.size():minG;
        minGoofy=(g<minGoofy)?g:minGoofy;
        mole.frac2kart( htest.frac, htest.pos);
        peaks.append(htest);
      }
    }
  }
  if (peaks.isEmpty()) {
    htest.frac=htest.pos=V3(0,0,0);
    return htest;
  }
  V3 D,floorD,prime,nf,dp;
  double min=1000,dk=0;
  for (int nn=0;nn<mole.cell.symmops.size();  nn++){
    prime=mole.cell.symmops.at(nn) * peaks.at(minG).frac + mole.cell.trans.at(nn);
    D=prime - mole.asymm.at(dindex).frac + V3(0.5,0.5,0.5) ;
    floorD=V3(floor(D.x),floor(D.y),floor(D.z));
    dp=D - floorD - V3(0.5,0.5,0.5);
    dk=mole.fl(dp.x,dp.y,dp.z);
    if (dk<min){
      min=dk;
      nf=prime-floorD;
    }
  }
  peaks[minG].frac=nf;
  mole.frac2kart(peaks[minG].frac,peaks[minG].pos);
  return peaks[minG];
}

void autoHFix(){
  if (mole.asymm.isEmpty()) {printf("no atoms?\n");return;}
  QStringList alab;
  for (int i =0; i<mole.asymm.size(); i++){
    alab.append(mole.asymm.at(i).Label);
  }
  QString nn;

  mole.grow();

  MyAtom htest,qtest;

  double xh[8]={
    0.98,//1
    0.97,//2
    0.96,//3
    0.93,//4
    1.10,//15 B
    0.82,//8 O
    0.93,//9
    0.93// 16
  };

  if (sfac.indexOf(0)==-1)addNewScatteringFactor(0);
  htest.an=0;
  htest.part=0;
  htest.resiNr=0;
  htest.hidden=0;
  htest.symmGroup=0;
  htest.sg=0;
  htest.scod=555;//the identity
  htest.auidx=-1;
  htest.uc.m12=htest.uc.m23=htest.uc.m13=0.0;
  htest.uc.m21=htest.uc.m32=htest.uc.m31=0.0;
  htest.uf.m12=htest.uf.m23=htest.uf.m13=0.0;
  htest.uf.m21=htest.uf.m32=htest.uf.m31=0.0;
  qtest=htest;
  qtest.an=-1;
  QMap<int,int> keinDonor;
  QMap<int,int> istDonor;
  for (int j=0;j<mole.contact.size();j++){
    double ch1,ch2;
    if (mole.knoepfe.at(mole.contact.at(j).a1).neighbors.size()==0) keinDonor[mole.contact.at(j).a1]++;
    if (mole.knoepfe.at(mole.contact.at(j).a2).neighbors.size()==0) keinDonor[mole.contact.at(j).a2]++;
    if ((mole.asymm.at(mole.contact.at(j).a1).an==7)&&(mole.knoepfe.at(mole.contact.at(j).a1).neighbors.size()>1)) keinDonor[mole.contact.at(j).a1]++;
    if ((mole.asymm.at(mole.contact.at(j).a2).an==7)&&(mole.knoepfe.at(mole.contact.at(j).a2).neighbors.size()>1)) keinDonor[mole.contact.at(j).a2]++;
    if ((mole.contact.at(j).covalent)){//
      ch1=ch2=0;
      for (int k=0; k<mole.knoepfe.at(mole.contact.at(j).a1).neighbors.size(); k++){
        ch1 += sqrt(Distance(mole.asymm.at(mole.contact.at(j).a1).pos,
              mole.showatoms.at(mole.knoepfe.at(mole.contact.at(j).a1).neighbors.at(k)).pos))-
          ((mole.Kovalenz_Radien[mole.asymm.at(mole.contact.at(j).a1).an]+
            mole.Kovalenz_Radien[mole.showatoms.at(mole.knoepfe.at(mole.contact.at(j).a1).neighbors.at(k)).an])*0.01);
      }
      if (mole.knoepfe.at(mole.contact.at(j).a1).neighbors.size())
        ch1*=1.0/mole.knoepfe.at(mole.contact.at(j).a1).neighbors.size();
      for (int k=0; k<mole.knoepfe.at(mole.contact.at(j).a2).neighbors.size(); k++){
        ch2 += sqrt(Distance(mole.asymm.at(mole.contact.at(j).a2).pos,
              mole.showatoms.at(mole.knoepfe.at(mole.contact.at(j).a2).neighbors.at(k)).pos))-
          ((mole.Kovalenz_Radien[mole.asymm.at(mole.contact.at(j).a2).an]+
            mole.Kovalenz_Radien[mole.showatoms.at(mole.knoepfe.at(mole.contact.at(j).a2).neighbors.at(k)).an])*0.01);
      }
      if (mole.knoepfe.at(mole.contact.at(j).a2).neighbors.size())
        ch2*=1.0/mole.knoepfe.at(mole.contact.at(j).a2).neighbors.size();
      if (ch1<-0.2 ) keinDonor[mole.contact.at(j).a1]++;
      if (ch2<-0.2 ) keinDonor[mole.contact.at(j).a2]++;
      if ((ch1-0.01)>ch2) istDonor[mole.contact.at(j).a1]++;
      if ((ch2-0.01)>ch1) istDonor[mole.contact.at(j).a2]++;
    }
  }
  for (int i =0; i<mole.asymm.size(); i++){
    if ((mole.asymm.at(i).an==5)||(mole.asymm.at(i).an==6)){
      double nh= (mole.asymm.at(i).an==6) ? 0.07 : 0;
      double plan=1.0;
      if (mole.knoepfe.at(i).neighbors.size()==1) {
        if ((mole.asymm.at(i).an==6)&&(!istDonor.contains(i))) continue;
        V3 v=V3(0,0,0),w=V3(0,0,0),u=V3(0,0,0) ,a=V3(0,0,0);//,b=V3(0,0,0),c=V3(0,0,0);
        a = mole.showatoms.at(mole.knoepfe.at(i).neighbors.at(0)).pos - mole.asymm.at(i).pos;
        double chi = sqrt(Norm(a))-((mole.Kovalenz_Radien[mole.showatoms.at(mole.knoepfe.at(i).neighbors.at(0)).an]+mole.Kovalenz_Radien[mole.asymm.at(i).an])/100.0);
        if (chi>-0.09) {//afix 137 127
          int j,k,l;
          a=Normalize(a);
          double jx,kx,lx,mx,my,mz;
          double val[24],merg[8],wink=0,t=0,mt=0;
          for (int z=0; z<8;z++) merg[z]=0;
          if (n5){
            for (int z=0; z<24; z++){
              if ((a==V3(0,0,-1))||(a==V3(0,0,1))) {
                w=Normalize(V3(1,1,1)%a);
              } else {
                w=Normalize(V3(0,0,1)%a);
              }
              u=Normalize(w%a);
              v=(-(xh[2]-nh)*0.33333*a)+ mole.asymm.at(i).pos+(xh[2]-nh)*0.942809041582*(sin(wink)*u+cos(wink)*w);
              qtest.pos=v;
              mole.kart2frac(v,v);
              int zplu=0;
              do {
                jx=((zplu+v.x)*n1+0.5);
                zplu++;
              }while (jx<0);
              zplu=0;
              do {
                kx=((zplu+v.y)*n2+0.5);
                zplu++;
              }while (kx<0);
              zplu=0;
              do {
                lx=((zplu+v.z)*n3+0.5);
                zplu++;
              }while (lx<0);
              mx=jx-floor(jx);
              my=kx-floor(kx);
              mz=lx-floor(lx);
              j=((int)jx)%n1;//
              k=((int)kx)%n2;//
              l=((int)lx)%n3;//
              merg[z%8]+=val[z]=(
                  (1-mx)*datfo_fc[(j+n1*(k+l*n2))%n5]+
                  mx*datfo_fc[((j+1)+n1*(k+l*n2))%n5]+
                  (1-my)*datfo_fc[(j+n1*(k+l*n2))%n5]+
                  my*datfo_fc[(j+n1*((k+1)+l*n2))%n5]+
                  (1-mz)*datfo_fc[(j+n1*(k+l*n2))%n5]+
                  mz*datfo_fc[(j+n1*(k+(l+1)*n2))%n5])/3.0;
              qtest.peakHeight=val[z];
              qtest.frac=v;
              qtest.Label=QString("Q%1#").arg(z,-2);
              mole.pmin=qMin(mole.pmin, qtest.peakHeight);
              mole.pmax=qMax(mole.pmax, qtest.peakHeight);
              mole.asymm.append(qtest);
              wink+=0.261799387799;
            }
            for (int z=0; z<8;z++) {
              t =merg[(z+8)%8] ;
              double x = t - merg[(z+9)%8];
              double u = t - merg[(z+7)%8];
              if (fmin(x,u)<0) continue;
              double p=0.5*(u-x)/(u+x);
              t=t+0.25*(p*(u-x));
              mt=(mt>t)?mt:t;
              if (t<mt)continue;
              wink=((double) z + p )*0.261799387799;
            }
          } else {
            int zzz=0;
            int z=mole.knoepfe.at(mole.knoepfe.at(i).neighbors.at(0)).neighbors.at(zzz++);
            while ((z==i)&&(zzz<mole.knoepfe.at(mole.knoepfe.at(i).neighbors.at(0)).neighbors.size())) {
              z=mole.knoepfe.at(mole.knoepfe.at(i).neighbors.at(0)).neighbors.at(zzz++);
            }
            if ((a==V3(0,0,-1))||(a==V3(0,0,1))) {
              w=Normalize(V3(1,1,1)%a);
            } else {
              w=Normalize(V3(0,0,1)%a);
            }
            u=Normalize(w%a);
            v=(-(xh[2]-nh)*0.33333*a)+ mole.asymm.at(i).pos+(xh[2]-nh)*0.942809041582*(sin(wink)*u+cos(wink)*w);
            V3 increm = V3(0.0, 0.0, 0.0);
            if ( z == i ) {
              increm = V3(0.00001, 0.00001, 0.00001);
            }
            wink=mole.dieder(mole.showatoms.at(mole.knoepfe.at(i).neighbors.at(0)).pos+increm-mole.asymm.at(z).pos,
                mole.asymm.at(i).pos-mole.showatoms.at(mole.knoepfe.at(i).neighbors.at(0)).pos,
                v-mole.asymm.at(i).pos)/180.0*M_PI;
          }
          v=(-(xh[2]-nh)*0.33333*a)+ mole.asymm.at(i).pos+(xh[2]-nh)*0.942809041582*(sin(wink)*u+cos(wink)*w);
          htest.Label=mole.asymm.at(i).Label.section('_',0,0);
          htest.part=mole.asymm.at(i).part;
          htest.resiNr=mole.asymm.at(i).resiNr;
          htest.Label.replace(0,1,'H');
          QString baslab=htest.Label;
          int hindex=1;
          QString nnn=htest.Label;
          do {
            htest.Label=baslab;
            if (hindex) nnn=htest.Label.left(3).append(QString::number(hindex+9,36));
            if (nnn.size()>4) nnn=htest.Label.left(2).append(QString::number(hindex+9,36));
            if (htest.resiNr)
              nn=QString("%1_%2")//for shelxl same names in different parts are no problem but for platon so this should fix it.
                .arg(nnn)
                .arg((htest.resiNr)?QString::number(htest.resiNr):"");
            else nn=nnn;
            hindex++;
          }while(alab.contains(nn,Qt::CaseInsensitive));
          htest.Label=nnn;
          alab.append(nn);
          double uIso = mole.ueq(mole.asymm.at(i).uf)*1.5;
          htest.uf.m11 = htest.uf.m22 = htest.uf.m33 = uIso;
          htest.uf.m23 = uIso * mole.cell.cosra;
          htest.uf.m13 = uIso * mole.cell.cosrb;
          htest.uf.m12 = uIso * mole.cell.cosrg;
          htest.isIso=true;
          htest.ufiso_org=QString("-1.5");
          mole.Uf2Uo(htest.uf,htest.uc);
          htest.pos=v;
          mole.kart2frac(htest.pos,htest.frac);
          htest.orginalLine=QString("%1  %2  %3 %4 %5 %6 %7")
            .arg(htest.Label)
            .arg(qMax(1,sfac.indexOf(0)+1))
            .arg(htest.frac.x,9,'f',5)
            .arg(htest.frac.y,9,'f',5)
            .arg(htest.frac.z,9,'f',5)
            .arg(mole.asymm.at(i).sof_org)
            .arg("-1.5");
          mole.asymm.append(htest);
          wink+=2.09439510239;
          v=(-(xh[2]-nh)*0.33333*a)+ mole.asymm.at(i).pos+(xh[2]-nh)*0.942809041582*(sin(wink)*u+cos(wink)*w);
          do {
            htest.Label=baslab;
            if (hindex) nnn=htest.Label.left(3).append(QString::number(hindex+9,36));
            if (nnn.size()>4) nnn=htest.Label.left(2).append(QString::number(hindex+9,36));
            if (htest.resiNr)
              nn=QString("%1_%2")//for shelxl same names in different parts are no problem but for platon so this should fix it.
                .arg(nnn)
                .arg((htest.resiNr)?QString::number(htest.resiNr):"");
            else nn=nnn;
            hindex++;
          }while(alab.contains(nn,Qt::CaseInsensitive));
          htest.Label=nnn;
          alab.append(nn);
          htest.pos=v;
          mole.kart2frac(htest.pos,htest.frac);
          htest.orginalLine=QString("%1  %2  %3 %4 %5 %6 %7")
            .arg(htest.Label)
            .arg(qMax(1,sfac.indexOf(0)+1))
            .arg(htest.frac.x,9,'f',5)
            .arg(htest.frac.y,9,'f',5)
            .arg(htest.frac.z,9,'f',5)
            .arg(mole.asymm.at(i).sof_org)
            .arg("-1.5");
          mole.asymm.append(htest);
          wink+=2.09439510239;
          v=(-(xh[2]-nh)*0.33333*a)+ mole.asymm.at(i).pos+(xh[2]-nh)*0.942809041582*(sin(wink)*u+cos(wink)*w);
          do {
            htest.Label=baslab;
            if (hindex) nnn=htest.Label.left(3).append(QString::number(hindex+9,36));
            if (nnn.size()>4) nnn=htest.Label.left(2).append(QString::number(hindex+9,36));
            if (htest.resiNr)
              nn=QString("%1_%2")//for shelxl same names in different parts are no problem but for platon so this should fix it.
                .arg(nnn)
                .arg((htest.resiNr)?QString::number(htest.resiNr):"");
            else nn=nnn;
            hindex++;
          }while(alab.contains(nn,Qt::CaseInsensitive));
          htest.Label=nnn;
          alab.append(nn);
          htest.pos=v;
          mole.kart2frac(htest.pos,htest.frac);
          htest.orginalLine=QString("%1  %2  %3 %4 %5 %6 %7")
            .arg(htest.Label)
            .arg(qMax(1,sfac.indexOf(0)+1))
            .arg(htest.frac.x,9,'f',5)
            .arg(htest.frac.y,9,'f',5)
            .arg(htest.frac.z,9,'f',5)
            .arg(mole.asymm.at(i).sof_org)
            .arg("-1.5");
          mole.asymm.append(htest);
          QString label = mole.showatoms.at(i).orginalLine;
          int cpos = resLines.indexOf(label);
          while (resLines.at(cpos).endsWith("=")) cpos++;
          resLines.insert(cpos+1,QString("AFIX 137\n%1\n%2\n%3\nAFIX %4")
              .arg(mole.asymm.at(mole.asymm.size()-3).orginalLine)
              .arg(mole.asymm.at(mole.asymm.size()-2).orginalLine)
              .arg(mole.asymm.at(mole.asymm.size()-1).orginalLine)
              .arg(mole.afix5(i))
              );
        }else if (chi> -0.27){//afix 93
          int j,k,l;
          a=Normalize(a);
          double jx,kx,lx,mx,my,mz;
          double val[24],merg[12],wink=0,t=0,mt=0;
          if (n5){
            for (int z=0; z<12;z++) merg[z]=0;
            for (int z=0; z<24; z++){
              wink=z*0.261799387799;
              w=Normalize(V3(0,0,1)%a);
              u=Normalize(w%a);
              v=(-(xh[6]-nh)*0.33333*a)+ mole.asymm.at(i).pos+(xh[6]-nh)*0.942809041582*(sin(wink)*u+cos(wink)*w);
              qtest.pos=v;
              mole.kart2frac(v,v);
              jx=(v.x>0)?(v.x*n1+0.5):((1.0+v.x)*n1+0.5);
              kx=(v.y>0)?(v.y*n2+0.5):((1.0+v.y)*n2+0.5);
              lx=(v.z>0)?(v.z*n3+0.5):((1.0+v.z)*n3+0.5);
              mx=jx-floor(jx);
              my=kx-floor(kx);
              mz=lx-floor(lx);
              j=((int)jx)%n1;//
              k=((int)kx)%n2;//
              l=((int)lx)%n3;//
              merg[z%12]+=val[z]=(
                  (1-mx)*datfo_fc[(j+n1*(k+l*n2))%n5]+
                  mx*datfo_fc[((j+1)+n1*(k+l*n2))%n5]+
                  (1-my)*datfo_fc[(j+n1*(k+l*n2))%n5]+
                  my*datfo_fc[(j+n1*((k+1)+l*n2))%n5]+
                  (1-mz)*datfo_fc[(j+n1*(k+l*n2))%n5]+
                  mz*datfo_fc[(j+n1*(k+(l+1)*n2))%n5])/3.0;
              qtest.peakHeight=val[z];
              qtest.frac=v;
              qtest.Label=QString("Q%1").arg(z+100);
              mole.pmin=qMin(mole.pmin, qtest.peakHeight);
              mole.pmax=qMax(mole.pmax, qtest.peakHeight);
              mole.asymm.append(qtest);
            }
            for (int z=0; z<12;z++) {
              t =merg[(z+12)%12] ;
              double x = t - merg[(z+13)%12];
              double u = t - merg[(z+11)%12];
              if (fmin(x,u)<0) continue;
              double p=0.5*(u-x)/(u+x);
              t=t+0.25*(p*(u-x));
              mt=(mt>t)?mt:t;
              if (t<mt)continue;
              wink=((double) z + p )*0.261799387799;
              /////////////////
            }
          }else{
            int zzz=0;
            int z=mole.knoepfe.at(mole.knoepfe.at(i).neighbors.at(0)).neighbors.at(zzz++);
            while ((z==i)&&(zzz<mole.knoepfe.at(mole.knoepfe.at(i).neighbors.at(0)).neighbors.size()))
              z=mole.knoepfe.at(mole.knoepfe.at(i).neighbors.at(0)).neighbors.at(zzz++);
            w=Normalize(V3(0,0,1)%a);
            u=Normalize(w%a);
            v=(-(xh[2]-nh)*0.33333*a)+ mole.asymm.at(i).pos+(xh[2]-nh)*0.942809041582*(sin(wink)*u+cos(wink)*w);
            V3 increm = V3(0.0, 0.0, 0.0);
            if ( z == i ) {
              increm = V3(0.00001, 0.00001, 0.00001);
            }
            wink=mole.dieder(mole.showatoms.at(mole.knoepfe.at(i).neighbors.at(0)).pos+increm-mole.asymm.at(z).pos,
                mole.asymm.at(i).pos-mole.showatoms.at(mole.knoepfe.at(i).neighbors.at(0)).pos,
                v-mole.asymm.at(i).pos)/180.0*M_PI;
          }
          v=(-(xh[6]-nh)*0.5*a)+ mole.asymm.at(i).pos+(xh[6]-nh)*0.866025403784*(sin(wink)*u+cos(wink)*w);
          htest.Label=mole.asymm.at(i).Label.section('_',0,0);
          htest.part=mole.asymm.at(i).part;
          htest.resiNr=mole.asymm.at(i).resiNr;
          htest.Label.replace(0,1,'H');
          int hindex=1;
          QString nnn=htest.Label;
          QString baslab=htest.Label;
          do {
            htest.Label=baslab;
            if (hindex) nnn=htest.Label.left(3).append(QString::number(hindex+9,36));
            if (nnn.size()>4) nnn=htest.Label.left(2).append(QString::number(hindex+9,36));
            if (htest.resiNr)
              nn=QString("%1_%2")//for shelxl same names in different parts are no problem but for platon so this should fix it.
                .arg(nnn)
                .arg((htest.resiNr)?QString::number(htest.resiNr):"");
            else nn=nnn;
            hindex++;
          }while(alab.contains(nn,Qt::CaseInsensitive));
          htest.Label=nnn;
          alab.append(nn);
          htest.isIso=true;
          htest.ufiso_org=QString("-1.2");
          double uIso = mole.ueq(mole.asymm.at(i).uf)*1.2;
          htest.uf.m11 = htest.uf.m22 = htest.uf.m33 = uIso;
          htest.uf.m23 = uIso * mole.cell.cosra;
          htest.uf.m13 = uIso * mole.cell.cosrb;
          htest.uf.m12 = uIso * mole.cell.cosrg;
          mole.Uf2Uo(htest.uf,htest.uc);
          htest.pos=v;
          mole.kart2frac(htest.pos,htest.frac);
          htest.orginalLine=QString("%1  %2  %3 %4 %5 %6 %7")
            .arg(htest.Label)
            .arg(qMax(1,sfac.indexOf(0)+1))
            .arg(htest.frac.x,9,'f',5)
            .arg(htest.frac.y,9,'f',5)
            .arg(htest.frac.z,9,'f',5)
            .arg(mole.asymm.at(i).sof_org)
            .arg("-1.2");
          mole.asymm.append(htest);
          wink+=3.14159265359;
          v=(-(xh[6]-nh)*0.5*a)+ mole.asymm.at(i).pos+(xh[6]-nh)*0.866025403784*(sin(wink)*u+cos(wink)*w);
          do {
            htest.Label=baslab;
            if (hindex) nnn=htest.Label.left(3).append(QString::number(hindex+9,36));
            if (nnn.size()>4) nnn=htest.Label.left(2).append(QString::number(hindex+9,36));
            if (htest.resiNr)
              nn=QString("%1_%2")//for shelxl same names in different parts are no problem but for platon so this should fix it.
                .arg(nnn)
                .arg((htest.resiNr)?QString::number(htest.resiNr):"");
            else nn=nnn;
            hindex++;
          }while(alab.contains(nn,Qt::CaseInsensitive));
          htest.Label=nnn;
          alab.append(nn);
          htest.pos=v;
          mole.kart2frac(htest.pos,htest.frac);
          htest.orginalLine=QString("%1  %2  %3 %4 %5 %6 %7")
            .arg(htest.Label)
            .arg(qMax(1,sfac.indexOf(0)+1))
            .arg(htest.frac.x,9,'f',5)
            .arg(htest.frac.y,9,'f',5)
            .arg(htest.frac.z,9,'f',5)
            .arg(mole.asymm.at(i).sof_org)
            .arg("-1.2");
          mole.asymm.append(htest);
          QString label = mole.showatoms.at(i).orginalLine;
          int cpos = resLines.indexOf(label);
          while (resLines.at(cpos).endsWith("=")) cpos++;
          resLines.insert(cpos+1,QString("AFIX 93\n%1\n%2\nAFIX %3")
              .arg(mole.asymm.at(mole.asymm.size()-2).orginalLine)
              .arg(mole.asymm.at(mole.asymm.size()-1).orginalLine)
              .arg(mole.afix5(i))
                          );
        }else{//afix 163
          if (mole.asymm.at(i).an==5){
            a = Normalize(mole.showatoms.at(mole.knoepfe.at(i).neighbors.at(0)).pos - mole.asymm.at(i).pos);
            v=(-(xh[7]-nh)*a)+ mole.asymm.at(i).pos;
            htest.Label=mole.asymm.at(i).Label.section('_',0,0);
            htest.part=mole.asymm.at(i).part;
            htest.resiNr=mole.asymm.at(i).resiNr;
            htest.Label.replace(0,1,'H');
            int hindex=0;
            QString nnn=htest.Label;
            do {
              if (hindex) nnn=htest.Label.left(3).append(QString::number(hindex+9,36));
              if (nnn.size()>4) nnn=htest.Label.left(2).append(QString::number(hindex+9,36));
              if (htest.resiNr)
                nn=QString("%1_%2")//for shelxl same names in different parts are no problem but for platon so this should fix it.
                  .arg(nnn)
                  .arg((htest.resiNr)?QString::number(htest.resiNr):"");
              else nn=nnn;
              hindex++;
            }while(alab.contains(nn,Qt::CaseInsensitive));
            htest.Label=nnn;
            alab.append(nn);
            htest.isIso=true;
            htest.ufiso_org=QString("-1.2");
            double uIso = mole.ueq(mole.asymm.at(i).uf)*1.2;
            htest.uf.m11 = htest.uf.m22 = htest.uf.m33 = uIso;
            htest.uf.m23 = uIso * mole.cell.cosra;
            htest.uf.m13 = uIso * mole.cell.cosrb;
            htest.uf.m12 = uIso * mole.cell.cosrg;
            mole.Uf2Uo(htest.uf,htest.uc);
            htest.pos=v;
            mole.kart2frac(htest.pos,htest.frac);
            htest.orginalLine=QString("%1  %2  %3 %4 %5 %6 %7")
              .arg(htest.Label)
              .arg(qMax(1,sfac.indexOf(0)+1))
              .arg(htest.frac.x,9,'f',5)
              .arg(htest.frac.y,9,'f',5)
              .arg(htest.frac.z,9,'f',5)
              .arg(mole.asymm.at(i).sof_org)
              .arg("-1.2");
            mole.asymm.append(htest);
            QString label = mole.showatoms.at(i).orginalLine;
            int cpos = resLines.indexOf(label);
            while (resLines.at(cpos).endsWith("=")) cpos++;
            resLines.insert(cpos+1,QString("AFIX 163\n%1\nAFIX %2")
                .arg(mole.asymm.at(mole.asymm.size()-1).orginalLine)
                .arg(mole.afix5(i))
                );
          }
        }
      }
      if (mole.knoepfe.at(i).neighbors.size()==2){//afix 43 oder 23
        if ((mole.asymm.at(i).an==6)&&(!istDonor.contains(i))) continue;
        double angle=0;
        V3 v=V3(0,0,0),w=V3(0,0,0),u=V3(0,0,0) ,a=V3(0,0,0),b=V3(0,0,0),c=V3(0,0,0);
        a = mole.showatoms.at(mole.knoepfe.at(i).neighbors.at(0)      ).pos - mole.asymm.at(i).pos;
        b = mole.showatoms.at(mole.knoepfe.at(i).neighbors.at(1)      ).pos - mole.asymm.at(i).pos;
        double chi1 = sqrt(Norm(a))-
          ((mole.Kovalenz_Radien[mole.showatoms.at(mole.knoepfe.at(i).neighbors.at(0)).an]+
            mole.Kovalenz_Radien[mole.asymm.at(i).an])/100.0);
        double chi2 = sqrt(Norm(b))-
          ((mole.Kovalenz_Radien[mole.showatoms.at(mole.knoepfe.at(i).neighbors.at(1)).an]+
            mole.Kovalenz_Radien[mole.asymm.at(i).an])/100.0);
        a = c = Normalize(a) ;
        c = c+ (b= Normalize(b));
        angle=mole.winkel(a,b);
        c = Normalize(c);
        double andi=(mole.showatoms.at(mole.knoepfe.at(i).neighbors.at(1)).an+mole.showatoms.at(mole.knoepfe.at(i).neighbors.at(0)).an)/2.0-5.0;
        andi= (andi!=0.0)?0.8:1.0;
        if (angle<160) {
          if ((andi*(chi1+chi2))<-0.13) {//43
            c = (-(xh[3]-nh) * c)+mole.asymm.at(i).pos;
            htest.Label=mole.asymm.at(i).Label.section('_',0,0);
            htest.part=mole.asymm.at(i).part;
            htest.resiNr=mole.asymm.at(i).resiNr;
            htest.Label.replace(0,1,'H');
            int hindex=0;
            QString nnn=htest.Label;
            do {
              if (hindex) nnn=htest.Label.left(3).append(QString::number(hindex+9,36));
              if (nnn.size()>4) nnn=htest.Label.left(2).append(QString::number(hindex+9,36));
              if (htest.resiNr)
                nn=QString("%1_%2")//for shelxl same names in different parts are no problem but for platon so this should fix it.
                  .arg(nnn)
                  .arg((htest.resiNr)?QString::number(htest.resiNr):"");
              else nn=nnn;
              hindex++;
            }while(alab.contains(nn,Qt::CaseInsensitive));
            htest.Label=nnn;
            alab.append(nn);
            htest.isIso=true;
            htest.ufiso_org=QString("-1.2");
            double uIso = mole.ueq(mole.asymm.at(i).uf)*1.2;
            htest.uf.m11 = htest.uf.m22 = htest.uf.m33 = uIso;
            htest.uf.m23 = uIso * mole.cell.cosra;
            htest.uf.m13 = uIso * mole.cell.cosrb;
            htest.uf.m12 = uIso * mole.cell.cosrg;
            mole.Uf2Uo(htest.uf,htest.uc);
            htest.pos=c;
            mole.kart2frac(htest.pos,htest.frac);
            htest.orginalLine=QString("%1  %2  %3 %4 %5 %6 %7")
              .arg(htest.Label)
              .arg(qMax(1,sfac.indexOf(0)+1))
              .arg(htest.frac.x,9,'f',5)
              .arg(htest.frac.y,9,'f',5)
              .arg(htest.frac.z,9,'f',5)
              .arg(mole.asymm.at(i).sof_org)
              .arg("-1.2");
            mole.asymm.append(htest);
            QString label = mole.showatoms.at(i).orginalLine;
            int cpos = resLines.indexOf(label);
            while (resLines.at(cpos).endsWith("=")) cpos++;
            resLines.insert(cpos+1,QString("AFIX 43\n%1\nAFIX %2")
                .arg(mole.asymm.at(mole.asymm.size()-1).orginalLine)
                .arg(mole.afix5(i)));
          }
          else{//23
            w = Normalize(a);
            u = Normalize(b);
            double si=1.0376-0.0346*Norm(w-u);    //=1.91063323625-((acos(u*w)-1.91063323625)/5.0);//-1.91063323625 = 109.47 grad ; 3.7 empirischer Wert
            c = (-cos(si)*(xh[1]-nh) * c);
            v = sin(si)  *(xh[1]-nh) * Normalize(c%a);
            a = mole.asymm.at(i).pos + c + v;
            b = mole.asymm.at(i).pos + c - v;
            c = c + mole.asymm.at(i).pos;
            htest.Label=mole.asymm.at(i).Label.section('_',0,0);
            htest.part=mole.asymm.at(i).part;
            htest.resiNr=mole.asymm.at(i).resiNr;
            htest.Label.replace(0,1,'H');
            int hindex=1;
            QString nnn=htest.Label;
            QString ann=htest.Label;
            do {
              if (hindex) nnn=htest.Label.left(3).append(QString::number(hindex+9,36));
              if (nnn.size()>4) nnn=htest.Label.left(2).append(QString::number(hindex+9,36));
              if (htest.resiNr)
                nn=QString("%1_%2")//for shelxl same names in different parts are no problem but for platon so this should fix it.
                  .arg(nnn)
                  .arg((htest.resiNr)?QString::number(htest.resiNr):"");
              else nn=nnn;
              hindex++;
            }while(alab.contains(nn,Qt::CaseInsensitive));
            htest.Label=nnn;
            alab.append(nn);
            htest.isIso=true;
            htest.ufiso_org=QString("-1.2");
            double uIso = mole.ueq(mole.asymm.at(i).uf)*1.2;
            htest.uf.m11 = htest.uf.m22 = htest.uf.m33 = uIso;
            htest.uf.m23 = uIso * mole.cell.cosra;
            htest.uf.m13 = uIso * mole.cell.cosrb;
            htest.uf.m12 = uIso * mole.cell.cosrg;
            mole.Uf2Uo(htest.uf,htest.uc);
            htest.pos=a;
            mole.kart2frac(htest.pos,htest.frac);
            htest.orginalLine=QString("%1  %2  %3 %4 %5 %6 %7")
              .arg(htest.Label)
              .arg(qMax(1,sfac.indexOf(0)+1))
              .arg(htest.frac.x,9,'f',5)
              .arg(htest.frac.y,9,'f',5)
              .arg(htest.frac.z,9,'f',5)
              .arg(mole.asymm.at(i).sof_org)
              .arg("-1.2");
            mole.asymm.append(htest);;
            do {
              if (hindex) nnn=ann.left(3).append(QString::number(hindex+9,36));
              if (nnn.size()>4) nnn=ann.left(2).append(QString::number(hindex+9,36));
              if (htest.resiNr)
                nn=QString("%1_%2")//for shelxl same names in different parts are no problem but for platon so this should fix it.
                  .arg(nnn)
                  .arg((htest.resiNr)?QString::number(htest.resiNr):"");
              else nn=nnn;
              hindex++;
            }while(alab.contains(nn,Qt::CaseInsensitive));
            htest.Label=nnn;
            alab.append(nn);
            htest.pos=b;
            mole.kart2frac(htest.pos,htest.frac);
            htest.orginalLine=QString("%1  %2  %3 %4 %5 %6 %7")
              .arg(htest.Label)
              .arg(qMax(1,sfac.indexOf(0)+1))
              .arg(htest.frac.x,9,'f',5)
              .arg(htest.frac.y,9,'f',5)
              .arg(htest.frac.z,9,'f',5)
              .arg(mole.asymm.at(i).sof_org)
              .arg("-1.2");
            mole.asymm.append(htest);
            QString label = mole.showatoms.at(i).orginalLine;
            int cpos = resLines.indexOf(label);
            while (resLines.at(cpos).endsWith("=")) cpos++;
            resLines.insert(cpos+1,QString("AFIX 23\n%1\n%2\nAFIX %3")
                .arg(mole.asymm.at(mole.asymm.size()-2).orginalLine)
                .arg(mole.asymm.at(mole.asymm.size()-1).orginalLine)
                .arg(mole.afix5(i)));
          }
        }
      }
      if (mole.knoepfe.at(i).neighbors.size()==3) {
        if ((mole.asymm.at(i).an==6)&&(!istDonor.contains(i))) continue;
        V3 v,w=V3(0,0,0),u ,a,b,c=V3(0,0,0);
        double ch1=0;
        for (int k=0; k<3; k++){
          a = Normalize(mole.showatoms.at(mole.knoepfe.at(i).neighbors.at(k)      ).pos - mole.asymm.at(i).pos);
          b = Normalize(mole.showatoms.at(mole.knoepfe.at(i).neighbors.at((k+1)%3)).pos - mole.asymm.at(i).pos);
          c = c + a;
          v = Normalize(a%b);
          if (k) {plan*=fabs(w*v);//1 for perfect flat
          }
          w = v;
          ch1 += sqrt(Distance(mole.asymm.at(i).pos,
                mole.showatoms.at(mole.knoepfe.at(i).neighbors.at(k)).pos))-
            ((mole.Kovalenz_Radien[mole.asymm.at(i).an]+
              mole.Kovalenz_Radien[mole.showatoms.at(mole.knoepfe.at(i).neighbors.at(k)).an])*0.01);
        }
        c = Normalize(c);
        ch1*=1.0/3;
        if ((plan<0.8)&&(ch1>-0.09)){// not flat and not to short bonds
          c *=-(xh[0]-nh);
          htest.Label=mole.asymm.at(i).Label.section('_',0,0);
          htest.part=mole.asymm.at(i).part;
          htest.resiNr=mole.asymm.at(i).resiNr;
          htest.Label.replace(0,1,'H');
          int hindex=0;
          QString nnn=htest.Label;
          do {
            if (hindex) nnn=htest.Label.left(3).append(QString::number(hindex+9,36));
            if (nnn.size()>4) nnn=htest.Label.left(2).append(QString::number(hindex+9,36));
            if (htest.resiNr)
              nn=QString("%1_%2")//for shelxl same names in different parts are no problem but for platon so this should fix it.
                .arg(nnn)
                .arg((htest.resiNr)?QString::number(htest.resiNr):"");
            else nn=nnn;
            hindex++;
          }while(alab.contains(nn,Qt::CaseInsensitive));
          htest.Label=nnn;
          alab.append(nn);
          htest.isIso=true;
          htest.ufiso_org=QString("-1.2");
          double uIso = mole.ueq(mole.asymm.at(i).uf)*1.2;
          htest.uf.m11 = htest.uf.m22 = htest.uf.m33 = uIso;
          htest.uf.m23 = uIso * mole.cell.cosra;
          htest.uf.m13 = uIso * mole.cell.cosrb;
          htest.uf.m12 = uIso * mole.cell.cosrg;
          mole.Uf2Uo(htest.uf,htest.uc);
          htest.pos = mole.asymm.at(i).pos + c;
          mole.kart2frac(htest.pos,htest.frac);
          QString label = mole.showatoms.at(i).orginalLine;

          htest.orginalLine=QString("%1  %2  %3 %4 %5 %6 %7")
            .arg(htest.Label)
            .arg(qMax(1,sfac.indexOf(0)+1))
            .arg(htest.frac.x,9,'f',5)
            .arg(htest.frac.y,9,'f',5)
            .arg(htest.frac.z,9,'f',5)
            .arg(mole.asymm.at(i).sof_org)
            .arg("-1.2");

          int cpos = resLines.indexOf(label);
          while (resLines.at(cpos).endsWith("=")) cpos++;
          resLines.insert(cpos+1,QString("AFIX 13\n%1\nAFIX %2").arg(htest.orginalLine).arg(mole.afix5(i)));
          mole.asymm.append(htest);
        }
      }
    }
    if (((mole.asymm.at(i).an==4)||(mole.asymm.at(i).an==5))&&((mole.knoepfe.at(i).neighbors.size()==4)||(mole.knoepfe.at(i).neighbors.size()==5))){//boron with 4 or 5 neighbors AFIX 153
        V3 vsum = V3(0.,0.,0.);
        int nbs=mole.knoepfe.at(i).neighbors.size();
        int nbr=0;
        for (int k=0; k<mole.knoepfe.at(i).neighbors.size();k++){
            if (abs(mole.showatoms.at(mole.knoepfe.at(i).neighbors.at(k)).an-5)<2){
                nbr++;
            vsum+=(mole.showatoms.at(mole.knoepfe.at(i).neighbors.at(k)      ).pos - mole.asymm.at(i).pos);
            }
        }
        //printf("153 test length = %f %d==%d %s\n",Norm(vsum), nbs, nbr, mole.asymm.at(i).Label.toStdString().c_str());
        if ((nbs==nbr)&&(Norm(vsum)>10)){
        vsum=Normalize(vsum);
        vsum*=-1.0*xh[4];
        htest.Label=mole.asymm.at(i).Label.section('_',0,0);
        htest.part=mole.asymm.at(i).part;
        htest.resiNr=mole.asymm.at(i).resiNr;
        htest.Label.replace(0,1,'H');
        int hindex=0;
        QString nnn=htest.Label;
        do {
          if (hindex) nnn=htest.Label.left(3).append(QString::number(hindex+9,36));
          if (nnn.size()>4) nnn=htest.Label.left(2).append(QString::number(hindex+9,36));
          if (htest.resiNr)
            nn=QString("%1_%2")//for shelxl same names in different parts are no problem but for platon so this should fix it.
              .arg(nnn)
              .arg((htest.resiNr)?QString::number(htest.resiNr):"");
          else nn=nnn;
          hindex++;
        }while(alab.contains(nn,Qt::CaseInsensitive));
        htest.Label=nnn;
        alab.append(nn);
        htest.isIso=true;
        htest.ufiso_org=QString("-1.2");
        double uIso = mole.ueq(mole.asymm.at(i).uf)*1.2;
        htest.uf.m11 = htest.uf.m22 = htest.uf.m33 = uIso;
        htest.uf.m23 = uIso * mole.cell.cosra;
        htest.uf.m13 = uIso * mole.cell.cosrb;
        htest.uf.m12 = uIso * mole.cell.cosrg;
        mole.Uf2Uo(htest.uf,htest.uc);
        htest.pos = mole.asymm.at(i).pos + vsum;
        mole.kart2frac(htest.pos,htest.frac);
        QString label = mole.showatoms.at(i).orginalLine;
        htest.orginalLine=QString("%1  %2  %3 %4 %5 %6 %7")
          .arg(htest.Label)
          .arg(qMax(1,sfac.indexOf(0)+1))
          .arg(htest.frac.x,9,'f',5)
          .arg(htest.frac.y,9,'f',5)
          .arg(htest.frac.z,9,'f',5)
          .arg(mole.asymm.at(i).sof_org)
          .arg("-1.2");
        int cpos = resLines.indexOf(label);
        while (resLines.at(cpos).endsWith("=")) cpos++;
        resLines.insert(cpos+1,QString("AFIX 153\n%1\nAFIX %2").arg(htest.orginalLine).arg(mole.afix5(i)));
        mole.asymm.append(htest);
        }
    }
  }
  for (int j=0;j<mole.contact.size();j++){
    if ((istDonor.contains(mole.contact.at(j).a1)||(!keinDonor.contains(mole.contact.at(j).a1)))
        &&(mole.asymm.at(mole.contact.at(j).a1).an==7)) {
      if (keinDonor.contains(mole.contact.at(j).a1)) continue;
      if (n5) {
        htest = findOH(mole.cell.symmops.at(mole.contact.at(j).sn) * mole.asymm.at(mole.contact.at(j).a1).frac
            + mole.cell.trans.at(mole.contact.at(j).sn) - mole.contact.at(j).floorD,
            mole.asymm.at(mole.contact.at(j).a2).frac,
            mole.contact.at(j).a1,alab);
        if (Distance(htest.pos,V3(0,0,0))<0.001){
          htest = findOH(mole.asymm.at(mole.contact.at(j).a1).frac,
              mole.cell.symmops.at(mole.contact.at(j).sn) * mole.asymm.at(mole.contact.at(j).a2).frac
              + mole.cell.trans.at(mole.contact.at(j).sn) ,
              mole.contact.at(j).a1,alab);
        } else {
          keinDonor[mole.contact.at(j).a1]++;
        }
        if (Distance(htest.pos,V3(0,0,0))){
          htest.orginalLine=QString("%1  %2  %3 %4 %5 %6 %7")
            .arg(htest.Label)
            .arg(qMax(1,sfac.indexOf(0)+1))
            .arg(htest.frac.x,9,'f',5)
            .arg(htest.frac.y,9,'f',5)
            .arg(htest.frac.z,9,'f',5)
            .arg(mole.asymm.at(mole.contact.at(j).a1).sof_org)
            .arg("-1.5");
          mole.asymm.append(htest);
          QString label = mole.showatoms.at(mole.contact.at(j).a1).orginalLine;
          int cpos = resLines.indexOf(label);
          while (resLines.at(cpos).endsWith("=")) cpos++;
          resLines.insert(cpos+1,QString("AFIX 148\n%1\nAFIX %2")
              .arg(mole.asymm.last().orginalLine)
              .arg(mole.afix5(mole.contact.at(j).a1))
              );
        }//Dist
      }//n5
    }// is O is in list of donors and not in list of no donors

  }// all contacts
}



double getNumber(double v, const QList<double> fv,int idx,int &fixFlag){
  double av=qAbs(v),res=0.0,var=1.0;
  int m=0;
  while ((-10*m+av) > 5){m++;}
  if (m>1) m=qMin(m,fv.size());
  if (m>1) var=fv.at(m-1);
  if (!m) res=v;
  else if (v>0) res=(av-(10*m))*var;
  else res=(av-(10*m))*(1.0-var);
  fvarCntr[m]++;
  if (m>0) fixFlag|=1<<idx; else if (fixFlag&1<<idx) fixFlag-=1<<idx;
  // printf("idx: %d shift %d  flag %d  m%d\n", idx,1<<idx,fixFlag,m);
  return res;
}

double getNumber(double v, const QList<double> fv, double uiso){
  double av=qAbs(v),res=0.0,var=1.0;
  if ((v<-0.5)&& (v>-5.0)) return res = av * uiso ;
  int m=0;
  while ((-10*m+av) > 5){m++;}
  if (m>1) m=qMin(m,fv.size());
  if (m>1) var=fv.at(m-1);
  if (!m) return v;
  else if (v>=0) res=(av-(10*m))*var;
  else res=(av-(10*m))*(1.0-var);
  fvarCntr[m]++;
  return res;
}

int isacommand(QString command){
  QStringList keywords;
  keywords <<
    "ACTA"<<//0
    "AFIX"<<//1
    "MPLA"<<//2
    "ANIS"<<//3
    "BASF"<<//4
    "BIND"<<//5
    "BLOC"<<//6
    "BOND"<<//7 AGS4 to MPLA
    "BUMP"<<//8
    "CELL"<<//9
    "CGLS"<<//10
    "CHIV"<<//11
    "CONF"<<//12
    "CONN"<<//13
    "DAMP"<<//14
    "DANG"<<//15
    "DEFS"<<//16
    "DELU"<<//17
    "DFIX"<<//18
    "DISP"<<//19
    "EADP"<<//20
    "EGEN"<<//21
    "END" <<//22
    "EQIV"<<//23
    "ESEL"<<//24
    "EXTI"<<//25
    "EXYZ"<<//26
    "FEND"<<//27
    "FLAT"<<//28
    "FMAP"<<//29
    "FRAG"<<//30
    "FREE"<<//31
    "FVAR"<<//32
    "GRID"<<//33
    "HFIX"<<//34
    "HKLF"<<//35
    "HOPE"<<//36
    "HTAB"<<//37
    "INIT"<<//38
    "ISOR"<<//39
    "LAST"<<//40
    "LATT"<<//41
    "LAUE"<<//42
    "LIST"<<//43
    "L.S."<<//44
    "MERG"<<//45
    "MOLE"<<//46
    "MORE"<<//47
    "MOVE"<<//48
    "NCSY"<<//49
    "OMIT"<<//50
    "PART"<<//51
    "PATT"<<//52
    "PHAN"<<//53
    "PHAS"<<//54
    "PLAN"<<//55
    "PSEE"<<//56
    "REM"<< //57
    "RESI"<<//58
    "RTAB"<<//59
    "SADI"<<//60
    "SAME"<<//61
    "SFAC"<<//62
    "SHEL"<<//63
    "SIMU"<<//64
    "SIZE"<<//65
    "SPEC"<<//66
    "SPIN"<<//67
    "STIR"<<//68
    "SUMP"<<//69
    "SWAT"<<//70
    "SYMM"<<//71
    "TEMP"<<//72
    "TEXP"<<//73
    "TIME"<<//74
    "TITL"<<//75
    "TREF"<<//76
    "TWIN"<<//77
    "UNIT"<<//78
    "VECT"<<//79
    "WPDB"<<//80
    "WGHT"<<//81
    "ZERR"<<//82
    "XNPD"<<//83
    "REST"<<//84
    "CHAN"<<//85
    "RIGU"<<//86
    "FLAP"<<//87
    "RNUM"<<//88
    "SOCC"<<//89
    "PRIG"<<//90
    "WIGL"<<//91
    "RANG"<<//92
    "TANG"<<//93
    "ADDA"<<//94
    "STAG"<<//95
    "ATOM"<<//96PDB dummy commands ...
    "HETA"<<//97
    "SCAL"<<//98
    "ABIN"<<//99
    "ANSC"<<//100
    "ANSR"<<//101
    "NOTR"<<//102
    "NEUT"<<//103
    "TWST"<<//104
    "BEDE"<<//105
    "LONE"//106<<
    ;//
  QString c=command;
  c=c.toUpper();
  c.remove(QRegExp("_\\S*"));
  if (c.startsWith("+")) return 666;//file inclusion
  return keywords.indexOf(c);
}

CEnvironment lonesome(QString s, int stnum){
  const QRegExp sep=QRegExp("\\s+");
  const QRegExp num=QRegExp("^\\d+");
  const QRegExp alp=QRegExp("^[\\$A-Za-z]+");
  QStringList tok=s.split(sep,skipEmptyParts);
  int m=0,j=0,k=0,ll=stnum;
  MyAtom horse;
  CEnvironment l;
  l.clear();
  MyAtom at;
  at.an=-42;
  at.auidx=-1;
  at.scod=555;
  int ignore=0;
  while (!tok.at(j).contains(num)) {j++;}
  m=tok.at(j).toInt();
  k=j;
  QString such;
  j=1;
  while ((!tok.at(j).contains(alp))) {j++;if (j>=tok.size()) break;}
  if (j==tok.size()) {
    qDebug()<<"oh-oh!";
    return l;//big problem
  }
  such=tok.at(j);
  j = k + 1;
  while (tok.at(j).contains(alp)) {j++;}
  double a   = getNumber(tok.at(j).toDouble(),fvar,0,ignore);
  j++;
  while (tok.at(j).contains(alp)) {j++;}
  double b1  = getNumber(tok.at(j).toDouble(),fvar,0,ignore);
  j++;
  while (tok.at(j).contains(alp)) {j++;}
  double b2  = getNumber(tok.at(j).toDouble(),fvar,0,ignore);
  bool general=true;
  for (int i=0; i<mole.asymm.size(); i++){
    if (such==mole.asymm.at(i).Label) {
      general=false;
      horse=mole.asymm[i];
      break;}
  }
  if (general){
    QString such2=such;
    such2.remove('$');
    int an1= mole.getOZ(such2);
    int jjj=stnum;
    for (int recurs=0; recurs<mole.asymm.size();recurs++){
      if (mole.asymm.at(recurs).an==an1){
        QString s2=s;
        s2=s2.replace(such,mole.asymm.at(recurs).Label);
        CEnvironment lrec;
        lrec=lonesome(s2, jjj);//recursion!
        jjj+=lrec.size();
        for (int dii=0; dii<lrec.size();dii++) {
          l.append(lrec[dii]);
        }
      }
    }
    return l;

  }else{
    switch (m) {
      case 1:{
               MyAtom n1,n2,n3;
               j=0;
               double d,soll=0;
               for (;j<mole.asymm.size();j++) {
                 if (mole.asymm.at(j).an<0) continue;
                 soll=(mole.Kovalenz_Radien[horse.an]+ mole.Kovalenz_Radien[mole.asymm.at(j).an])*0.012;
                 d= sqrt(Distance(horse.pos,mole.asymm.at(j).pos));
                 if ((d>0.1)&&(d<soll)){
                   n1=mole.asymm[j];
                   j++;
                   break;
                 }
               }
               for (;j<mole.asymm.size();j++) {
                 if (mole.asymm.at(j).an<0) continue;
                 soll=(mole.Kovalenz_Radien[horse.an]+ mole.Kovalenz_Radien[mole.asymm.at(j).an])*0.012;
                 d= sqrt(Distance(horse.pos,mole.asymm.at(j).pos));
                 if ((d>0.1)&&(d<soll)){
                   n2=mole.asymm[j];
                   j++;
                   break;
                 }
               }
               for (;j<mole.asymm.size();j++) {
                 if (mole.asymm.at(j).an<0) continue;
                 soll=(mole.Kovalenz_Radien[horse.an]+ mole.Kovalenz_Radien[mole.asymm.at(j).an])*0.012;
                 d= sqrt(Distance(horse.pos,mole.asymm.at(j).pos));
                 if ((d>0.1)&&(d<soll)){
                   n3=mole.asymm[j];
                   j++;
                   break;
                 }
               }
               V3 va= Normalize(n1.pos-horse.pos);
               V3 vb= Normalize(n2.pos-horse.pos);
               V3 vc= Normalize(n3.pos-horse.pos);
               V3 trisect=Normalize(va+vb+vc);
               double rad=tok.last().toDouble();
               if (rad<0) {
                 rad=-rad;
               }
               at.pos=horse.pos-(trisect*rad);
               at.uf=mole.u2b(horse.uf);//possible problem if U is not anisotrop
               mole.kart2frac(at.pos,at.frac);
               at.peakHeight=b1;
               at.sof=a;
               at.auidx=-1;
               at.Label=QString("L%1").arg(ll++);
               l.append(at);
               at.pos=horse.pos;
               at.frac=horse.frac;
               at.peakHeight=b2;
               at.sof=-a*l.size();
               at.Label=QString("L%1").arg(ll++);
               l.append(at);
               if (tok.size()>7) {
                 QString s2=s;
                 s2.remove(horse.Label);
                 l.append(lonesome(s2,ll));
               }
             }
             break;
      case 2:{
               MyAtom n1,n2;
               j=0;
               double d,soll=0;
               for (;j<mole.asymm.size();j++) {
                 if (mole.asymm.at(j).an<0) continue;
                 soll=(mole.Kovalenz_Radien[horse.an]+ mole.Kovalenz_Radien[mole.asymm.at(j).an])*0.012;
                 d= sqrt(Distance(horse.pos,mole.asymm.at(j).pos));
                 if ((d>0.1)&&(d<soll)){
                   n1=mole.asymm[j];
                   j++;
                   break;
                 }
               }
               for (;j<mole.asymm.size();j++) {
                 if (mole.asymm.at(j).an<0) continue;
                 soll=(mole.Kovalenz_Radien[horse.an]+ mole.Kovalenz_Radien[mole.asymm.at(j).an])*0.012;
                 d= sqrt(Distance(horse.pos,mole.asymm.at(j).pos));
                 if ((d>0.1)&&(d<soll)){
                   n2=mole.asymm[j];
                   j++;
                   break;
                 }
               }
               V3 va= Normalize(n1.pos-horse.pos);
               V3 vb= Normalize(n2.pos-horse.pos);
               V3 bisect=Normalize(va+vb);
               V3 ax=Normalize(vb%va);//LLLL
               double phi=tok.last().toDouble();
               double rad=tok.at(tok.size()-2).toDouble();
               phi*=M_PI/360.0;//we want /2 !
               V3 c=-cos(phi)*bisect;
               if (rad<0) {
                 c=-1.0*c;
                 rad=-rad;
               }
               V3 v=sin(phi)*ax;
               at.pos=horse.pos+(c+v)*rad;
               at.uf=mole.u2b(horse.uf);//possible problem if U is not anisotrop
               mole.kart2frac(at.pos,at.frac);
               at.peakHeight=b1;
               at.sof=a;
               at.auidx=-1;
               at.Label=QString("L%1").arg(ll++);
               l.append(at);
               at.pos=horse.pos+(c-v)*rad;
               mole.kart2frac(at.pos,at.frac);
               at.Label=QString("L%1").arg(ll++);
               l.append(at);
               at.pos=horse.pos;
               at.frac=horse.frac;
               at.peakHeight=b2;
               at.sof=-a*l.size();
               at.Label=QString("L%1").arg(ll++);
               l.append(at);
               if (tok.size()>8) {
                 QString s2=s;
                 s2.remove(horse.Label);
                 l.append(lonesome(s2,ll));
               }
             }
             break;

      case 4:{
               MyAtom n1,n2;
               j=0;
               double d,soll=0;
               for (;j<mole.asymm.size();j++) {
                 if (mole.asymm.at(j).an<0) continue;
                 soll=(mole.Kovalenz_Radien[horse.an]+ mole.Kovalenz_Radien[mole.asymm.at(j).an])*0.012;
                 d= sqrt(Distance(horse.pos,mole.asymm.at(j).pos));
                 if ((d>0.1)&&(d<soll)){
                   n1=mole.asymm[j];
                   j++;
                   break;
                 }
               }
               for (;j<mole.asymm.size();j++) {
                 if (mole.asymm.at(j).an<0) continue;
                 soll=(mole.Kovalenz_Radien[horse.an]+ mole.Kovalenz_Radien[mole.asymm.at(j).an])*0.012;
                 d= sqrt(Distance(horse.pos,mole.asymm.at(j).pos));
                 if ((d>0.1)&&(d<soll)){
                   n2=mole.asymm[j];
                   j++;
                   break;
                 }
               }
               V3 va= Normalize(n1.pos-horse.pos);
               V3 vb= Normalize(n2.pos-horse.pos);
               V3 bisect=Normalize(va+vb);
               //V3 ax=Normalize(vb%va);//LLLL
               double rad=tok.last().toDouble();
               at.pos=horse.pos-bisect*rad;
               at.uf=mole.u2b(horse.uf);//possible problem if U is not anisotrop
               mole.kart2frac(at.pos,at.frac);
               at.peakHeight=b1;
               at.sof=a;
               at.Label=QString("L%1").arg(ll++);
               l.append(at);
               at.pos=horse.pos;
               at.frac=horse.frac;
               at.peakHeight=b2;
               at.sof=-a;
               at.Label=QString("L%1").arg(ll++);
               l.append(at);
               if (tok.size()>7) {
                 QString s2=s;
                 s2.remove(horse.Label);
                 l.append(lonesome(s2,ll));
               }
             }
             break;
      case 6:{
               MyAtom n1,n2;
               j=0;
               double d,soll=0;
               for (;j<mole.asymm.size();j++) {
                 if (mole.asymm.at(j).an<0) continue;
                 soll=(mole.Kovalenz_Radien[horse.an]+ mole.Kovalenz_Radien[mole.asymm.at(j).an])*0.012;
                 d= sqrt(Distance(horse.pos,mole.asymm.at(j).pos));
                 if ((d>0.1)&&(d<soll)){
                   n1=mole.asymm[j];
                   j++;
                   break;
                 }
               }
               for (;j<mole.asymm.size();j++) {
                 if (mole.asymm.at(j).an<0) continue;
                 soll=(mole.Kovalenz_Radien[horse.an]+ mole.Kovalenz_Radien[mole.asymm.at(j).an])*0.012;
                 d= sqrt(Distance(horse.pos,mole.asymm.at(j).pos));
                 if ((d>0.1)&&(d<soll)){
                   n2=mole.asymm[j];
                   j++;
                   break;
                 }
               }
               V3 va= Normalize(n1.pos-horse.pos);
               V3 vb= Normalize(n2.pos-horse.pos);
               V3 ax=Normalize(vb%va);//LLLL
               double rad=tok.last().toDouble();
               if (rad<0) {
                 rad=-rad;
               }
               at.pos=horse.pos+(ax)*rad;
               at.uf=mole.u2b(horse.uf);//possible problem if U is not anisotrop
               mole.kart2frac(at.pos,at.frac);
               at.peakHeight=b1;
               at.sof=a;
               at.Label=QString("L%1").arg(ll++);
               l.append(at);
               at.pos=horse.pos-(ax)*rad;
               mole.kart2frac(at.pos,at.frac);
               at.Label=QString("L%1").arg(ll++);
               l.append(at);
               at.pos=horse.pos;
               at.frac=horse.frac;
               at.peakHeight=b2;
               at.sof=-a*l.size();
               at.Label=QString("L%1").arg(ll++);
               l.append(at);
               if (tok.size()>7) {
                 QString s2=s;
                 s2.remove(horse.Label);
                 l.append(lonesome(s2,ll));
               }
             }
             break;
      case 7:{
               MyAtom n1,n2;
               j=0;
               double d,soll=0,d2;
               for (;j<mole.asymm.size();j++) {
                 if (mole.asymm.at(j).an<0) continue;
                 soll=(mole.Kovalenz_Radien[horse.an]+ mole.Kovalenz_Radien[mole.asymm.at(j).an])*0.012;
                 d= sqrt(Distance(horse.pos,mole.asymm.at(j).pos));
                 if ((d>0.1)&&(d<soll)){
                   n1=mole.asymm[j];
                   j++;
                   break;
                 }
               }
               for (j=0;j<mole.asymm.size();j++) {
                 if (mole.asymm.at(j).an<0) continue;
                 soll=(mole.Kovalenz_Radien[n1.an]+ mole.Kovalenz_Radien[mole.asymm.at(j).an])*0.012;
                 d= sqrt(Distance(n1.pos,mole.asymm.at(j).pos));
                 d2= sqrt(Distance(horse.pos,mole.asymm.at(j).pos));
                 if ((d2>1.1)&&(d>0.1)&&(d<soll)){
                   n2=mole.asymm[j];
                   j++;
                   break;
                 }
               }
               V3 va= Normalize(n1.pos-horse.pos);
               V3 vb= Normalize(n2.pos-n1.pos);
               V3 ax=Normalize(vb%va);//LLLL
               double phi=tok.last().toDouble();
               double rad=tok.at(tok.size()-2).toDouble();
               phi*=M_PI/360.0;//we want /2 !
               V3 c=-cos(phi)*va;
               if (rad<0) {
                 c=-1.0*c;
                 rad=-rad;
               }
               V3 v=sin(phi)*ax;
               at.pos=horse.pos+(c+v)*rad;
               at.uf=mole.u2b(horse.uf);//possible problem if U is not anisotrop
               mole.kart2frac(at.pos,at.frac);
               at.peakHeight=b1;
               at.sof=a;
               at.Label=QString("L%1").arg(ll++);
               l.append(at);
               at.pos=horse.pos+(c-v)*rad;
               mole.kart2frac(at.pos,at.frac);
               at.Label=QString("L%1").arg(ll++);
               l.append(at);
               at.pos=horse.pos;
               at.frac=horse.frac;
               at.peakHeight=b2;
               at.sof=-a*l.size();
               at.Label=QString("L%1").arg(ll++);
               l.append(at);
               if (tok.size()>8) {
                 QString s2=s;
                 s2.remove(horse.Label);
                 l.append(lonesome(s2,ll));
               }
             }
             break;
      case 9:{
               MyAtom n1,n2;
               j=0;
               double d,soll=0,d2;
               for (;j<mole.asymm.size();j++) {
                 if (mole.asymm.at(j).an<0) continue;
                 soll=(mole.Kovalenz_Radien[horse.an]+ mole.Kovalenz_Radien[mole.asymm.at(j).an])*0.012;
                 d= sqrt(Distance(horse.pos,mole.asymm.at(j).pos));
                 if ((d>0.1)&&(d<soll)){
                   n1=mole.asymm[j];
                   j++;
                   break;
                 }
               }
               for (j=0;j<mole.asymm.size();j++) {
                 if (mole.asymm.at(j).an<0) continue;
                 soll=(mole.Kovalenz_Radien[n1.an]+ mole.Kovalenz_Radien[mole.asymm.at(j).an])*0.012;
                 d= sqrt(Distance(n1.pos,mole.asymm.at(j).pos));
                 d2= sqrt(Distance(horse.pos,mole.asymm.at(j).pos));
                 if ((d2>1.1)&&(d>0.1)&&(d<soll)){
                   n2=mole.asymm[j];
                   j++;
                   break;
                 }
               }
               V3 va= Normalize(n1.pos-horse.pos);
               V3 vb= Normalize(n2.pos-n1.pos);
               V3 ax1=Normalize(vb%va);//LLLL
               V3 ax=Normalize(va%ax1);//LLLL
               double phi=tok.last().toDouble();
               double rad=tok.at(tok.size()-2).toDouble();
               phi*=M_PI/360.0;//we want /2 !
               V3 c=-cos(phi)*va;
               if (rad<0) {
                 c=-1.0*c;
                 rad=-rad;
               }
               V3 v=sin(phi)*ax;
               at.pos=horse.pos+(c+v)*rad;
               at.uf=mole.u2b(horse.uf);//possible problem if U is not anisotrop
               mole.kart2frac(at.pos,at.frac);
               at.peakHeight=b1;
               at.sof=a;
               at.Label=QString("L%1").arg(ll++);
               l.append(at);
               at.pos=horse.pos+(c-v)*rad;
               mole.kart2frac(at.pos,at.frac);
               at.Label=QString("L%1").arg(ll++);
               l.append(at);
               at.pos=horse.pos;
               at.frac=horse.frac;
               at.peakHeight=b2;
               at.sof=-a*l.size();
               at.Label=QString("L%1").arg(ll++);
               l.append(at);
               if (tok.size()>8) {
                 QString s2=s;
                 s2.remove(horse.Label);
                 l.append(lonesome(s2,ll));
               }
             }
             break;
      default:
             fprintf(stderr,"%d das kann ich noch nicht! %s\n",m,such.toStdString().c_str());
    }
  }
  return l;
}


void load_sheldrick(QString fileName, QString &inhalt){
  //! Parses the SHELX res/ins file into data structures.
  //printf("%d\n",__LINE__);
  bool neutrons=false;
  QString problem;
  int minp=1,maxp=0;
  bool qbeforehkl=false,hklf=false;
  mole.cell.symmops.clear();
  mole.cell.trans.clear();
  mole.asymm.clear();
  mole.sdm.clear();
  mole.showbonds.clear();
  mole.showatoms.clear();
  mole.newstructure();
  mole.lbonds.clear();
  mole.envi_sdm.clear();
  mole.contact.clear();
  mole.envi_sdm.clear();
  mole.contact.clear();
  mole.symmopsEQIV.clear();
  mole.labelEQIV.clear();
  mole.transEQIV.clear();
  mole.freeatoms.clear();
  mole.bindatoms.clear();
  for (int i=0; i<mole.knoepfe.size();i++) mole.knoepfe[i].neighbors.clear();
  mole.knoepfe.clear();
  mole.showbonds.clear();
  mole.showatoms.clear();
  V3 nl(0,0,0);
  QStringList bedeRloneLY;
  mole.cell.trans.append(nl);
  mole.cell.symmops.append(Matrix(1,0,0, 0,1,0, 0,0,1));
  mole.pmin=10000;
  mole.pmax=-10000;
  mole.lmin=10000;
  mole.lmax=-10000;
  int firstAtomLine = 0;
  int ls_cycls = 0;
  MyAtom newAtom;
  newAtom.part=0;
  newAtom.resiNr=0;
  newAtom.hidden=0;
  newAtom.symmGroup=0;
  newAtom.sg=0;
  newAtom.scod=555;//the identity
  newAtom.auidx=-1;
  newAtom.sof=0;
  newAtom.sof_org=0;
  newAtom.afix=0;
  newAtom.isIso =false;
  //bool virg=true;
  Scatt s(&mole);
  int part=0,afix=0;
  double part_fvar=11.0;
  double uiso=0.05;
  int lattice=1;
  int isoat=0,adpat=0,qpeaks=0;
  int afixparent=0;
  mole.exti=-666.0;
  mole.swat=-666.0;
  sfac.clear();
  fvar.clear();
  mole.osf=0.0;
  mole.hklScale=mole.hklSigmaScale=1.0;
  mole.hklMat=Matrix(1,0,0,0,1,0,0,0,1);
  mole.hklOmitSig=-2;
  mole.hklOmit2th=180.0;
  mole.hklShellLow=99999.0;
  mole.hklShellHig=0.0;
  fvarCntr.clear();
  mole.transEQIV.clear();
  mole.symmopsEQIV.clear();
  mole.labelEQIV.clear();
  mole.freeatoms.clear();
  mole.bindatoms.clear();
  mole.envi_sdm.clear();
  mole.bindPart.clear();

  QRegExp sep=QRegExp("\\s+");
  // Disabled SkipEmptyParts because otherwise i is not in sync with the line number!
  QStringList lines = inhalt.split(QRegExp("[\\r\\n]")); //,skipEmptyParts);
  bool fragIgnore=false;
  for (int i=0; i<lines.size(); i++){
    if (lines.at(i).startsWith(" ")) continue;
    if (!lines.at(i).isEmpty()){
      newAtom.orginalLine=lines.at(i);//.section("=",0,0);
      lines[i].remove("=");
      int cmd=isacommand(lines.at(i).section(sep,0,0));
      QString resiSpez=lines.at(i).section(sep,0,0).section('_',1,1);
      if (cmd==30) fragIgnore=true;
      if (cmd==27) fragIgnore=false;
      if (!fragIgnore)
        switch (cmd) {
          case 1: //AFIX
            afix=lines.at(i).section(sep,1,1).toInt();
            newAtom.afix=afix;
            break;
          case 5: {//bind
                    if (!lines.at(i).section(sep,1,2).contains(QRegExp("[A-Za-z]"))){
                      int pn01=0,pn02=0;
                      bool ok1,ok2;
                      pn01=lines.at(i).section(sep,1,1).toInt(&ok1);
                      pn02=lines.at(i).section(sep,2,2).toInt(&ok2);
                      if ((ok1)&&(ok2)) mole.bindPart[pn01]=pn02;
                      if ((ok1)&&(ok2)) mole.bindPart[pn02]=pn01;
                    }
                    MyBind aa;
                    aa.Lab1=lines.at(i).section(sep,1,1).toUpper();
                    aa.Lab2=lines.at(i).section(sep,2,2).toUpper();
                    mole.bindatoms.append(aa);
                  }
                  break;
          case 9: { //CELL
                    mole.cell.wave = lines.at(i).section(sep,1,1).toDouble();
                    mole.cell.a = lines.at(i).section(sep,2,2).toDouble();
                    mole.cell.b = lines.at(i).section(sep,3,3).toDouble();
                    mole.cell.c = lines.at(i).section(sep,4,4).toDouble();
                    mole.cell.al= lines.at(i).section(sep,5,5).toDouble();
                    mole.cell.be = lines.at(i).section(sep,6,6).toDouble();
                    mole.cell.ga = lines.at(i).section(sep,7,7).toDouble();
                    mole.cellSetup();
                  }
                  break;
          case 10: { // L.S.
                     ls_cycls = lines.at(i).section(sep, 1, 1).toInt();
                   }
                   break;
          case 19:{//DISP
                    QStringList tok = lines.at(i).split(sep,skipEmptyParts);
                    s.setDISP(tok);
                  }
                  break;

          case 23:{ //EQIV
                    mole.decodeSymmCardEQIV(lines.at(i));
                  }
                  break;
          case 25:{//EXTI
                    mole.exti=lines.at(i).section(sep,1,1).toDouble();
                  }
                  break;
          case 31:{//free
                    MyBind aa;
                    aa.Lab1=lines.at(i).section(sep,1,1).toUpper();
                    aa.Lab2=lines.at(i).section(sep,2,2).toUpper();
                    mole.freeatoms.append(aa);
                  }
                  break;
          case 32:{ //FVAR
                    QStringList tok = lines.at(i).split(sep,skipEmptyParts);
                    for (int ifv=1; ifv<tok.size();ifv++){
                      fvar.append(tok.at(ifv).toDouble());
                    }
                    if (!fvar.isEmpty()) {
                      mole.osf=fvar.at(0);
                    }
                  }
                  break;
          case 35: //HKLF
                  {
                    hklf=true;
                    QStringList tok = lines.at(i).split(sep,skipEmptyParts);
                    if (tok.size()>1)mole.hklf=tok.at(1).toInt();
                    if (tok.size()>2)mole.hklScale=tok.at(2).toDouble();
                    if (tok.size()>11){
                      mole.hklMat.m11=tok.at(3).toDouble();
                      mole.hklMat.m12=tok.at(4).toDouble();
                      mole.hklMat.m13=tok.at(5).toDouble();
                      mole.hklMat.m21=tok.at(6).toDouble();
                      mole.hklMat.m22=tok.at(7).toDouble();
                      mole.hklMat.m23=tok.at(8).toDouble();
                      mole.hklMat.m31=tok.at(9).toDouble();
                      mole.hklMat.m32=tok.at(10).toDouble();
                      mole.hklMat.m33=tok.at(11).toDouble();

                    }
                    if (tok.size()>12)mole.hklSigmaScale=tok.at(12).toDouble();
                  }
                  break;
          case 41: { //LATT
                     lattice=lines.at(i).section(sep,1,1).toInt();
                     break;
                   }
          case 44: { // CGLS
                     ls_cycls = lines.at(i).section(sep, 1, 1).toInt();
                   }
                   break;
          case 50://OMIT
                   {
                     QStringList tok = lines.at(i).split(sep,skipEmptyParts);
                     QRegExp num=QRegExp("$[A-Za-z]");
                     if ((tok.size()==3)&&(!tok.at(1).contains(num))&&(!tok.at(2).contains(num))){
                       //hklOmitSig,hklOmit2th,hklShellLow,hklShellHig
                       mole.hklOmitSig=tok.at(1).toDouble();
                       mole.hklOmit2th=tok.at(2).toDouble();
                       mole.hklShellHig=0.5*mole.cell.wave/sin(mole.hklOmit2th/360.0*M_PI);
                     }
                   }
                   break;
          case 51: { //PART n sof
                     part = lines.at(i).section(sep,1,1).toInt();
                     // part_fvar stores sof (and fvar) to decide during atom parsing
                     part_fvar = lines.at(i).section(sep, 2, 2).toDouble();
                     minp=(minp<part)?minp:part;
                     maxp=(maxp>part)?maxp:part;
                   }
                   break;
          case 58: {//RESI
                     if (lines.at(i).section(sep,1,1).indexOf(QRegExp("^[0-9]+"))!=-1) {
                       newAtom.resiNr=lines.at(i).section(sep,1,1).toInt();
                       newAtom.ResiClass=lines.at(i).section(sep,2,2);
                     } else {
                       newAtom.resiNr=lines.at(i).section(sep,2,2).toInt();
                       newAtom.ResiClass=lines.at(i).section(sep,1,1);
                     }
                     break;
                   }
          case 62: { //SFAC
                     QStringList tok = lines.at(i).split(sep,skipEmptyParts);
                     if ((tok.size()>4)&&(!tok.at(2).contains(QRegExp("[A-Za-z]")))) {
                       s.setSCAT(tok);
                       if ((tok.at(2).toDouble()==0.0)&&((tok.at(3).toDouble()==1.0))) {
                         neutrons=true;
                       }
                       sfac.append(mole.getOZ(tok.at(1)));
                     } else {
                       for (int isf=1; isf<tok.size();isf++) {
                         sfac.append(mole.getOZ(tok.at(isf)));
                       }
                     }
                   }
                   break;
          case 63:{//shel
                    QStringList tok = lines.at(i).split(sep,skipEmptyParts);
                    if (tok.size()>2){
                      mole.hklShellLow=tok.at(1).toDouble();
                      mole.hklShellHig=tok.at(2).toDouble();
                    }
                  }break;
          case 69:{ //SUMP
                    QStringList forml= lines.at(i).split(sep);
                    int fvix=4;
                    while (forml.size()>fvix){
                      int m=forml.at(fvix).toInt();
                      fvarCntr[m]++;
                      fvix+=2;
                    }
                  }
                  break;
          case 70:{
                    mole.swat=0.0;
                    mole.exti=2.0;
                    QStringList tok = lines.at(i).split(sep,skipEmptyParts);
                    if (tok.size()>1) mole.swat=tok.at(1).toDouble();
                    if (tok.size()>2) mole.exti=tok.at(2).toDouble();
                  }
                  break;

          case 71: {// SYMM
                     mole.decodeSymmCard(lines.at(i).section(" ",1,130));
                     break;
                   }
          case 75: {// TITL
                     mole.titl=lines.at(i).section(" ",1,-1).simplified();
                     break;
                   }
          case 78: {// UNIT
                     mole.applyLatticeCentro(lattice);

                     break;
                   }
          case 81: {// WGHT
                     break;
                   }
          case 82: {// ZERR
                     mole.cell.Z= lines.at(i);
                     break;
                   }
          case 103: {// NEUT
                      neutrons=true;
                      break;
                    }
          case 105:// BEDE
          case 106:{// LONE
                     bedeRloneLY.append(lines.at(i));
                     break;
                   }
          case -1:{//an atom or an error!
                    newAtom.isIso=false;
                    newAtom.ufiso_org="";
                    if (firstAtomLine == 0) {
                      firstAtomLine = i;
                    }
                    if ((newAtom.afix>170)&&((newAtom.afix<180))) break;
                    QStringList tok = lines.at(i).split(sep,skipEmptyParts);
                    newAtom.part=part;
                    if (newAtom.resiNr<0){
                      newAtom.resiNr=0;
                      newAtom.ResiClass="";
                    }
                    if (tok.size()==7){
                      // We have an isotropic atom!
                      if ((newAtom.part!=0)||(newAtom.resiNr!=0)){
                        newAtom.Label=QString("%1_%3%4")
                          .arg(tok.at(0))
                          .arg((newAtom.resiNr)?QString::number(newAtom.resiNr):"")
                          .arg((newAtom.part)?QString::number((newAtom.part<0)?
                                36+newAtom.part:
                                newAtom.part+9,36):"");
                      } else {
                        newAtom.Label=tok.at(0);
                      }
                      int fac=tok.at(1).toInt()-1;
                      newAtom.fixFlag=0;
                      newAtom.an=((fac<0)||(fac>=sfac.size()))?-2:sfac.at(fac);
                      newAtom.frac.x = getNumber(tok.at(2).toDouble(),fvar,0,newAtom.fixFlag);
                      newAtom.frac.y = getNumber(tok.at(3).toDouble(),fvar,1,newAtom.fixFlag);
                      newAtom.frac.z = getNumber(tok.at(4).toDouble(),fvar,2,newAtom.fixFlag);
                      if (qAbs(part_fvar) > 11.0 ){
                        newAtom.sof = getNumber(part_fvar, fvar, 3, newAtom.fixFlag);
                        newAtom.sof_org = part_fvar;
                      } else {
                        newAtom.sof = getNumber(tok.at(5).toDouble(), fvar, 3, newAtom.fixFlag);
                        newAtom.sof_org = tok.at(5).toDouble();
                      }
                      newAtom.isIso=true;
                      double uIso = getNumber(tok.at(6).toDouble(),fvar,uiso);
                      //printf("UISO? %f %f\n",uIso,uiso);
                      newAtom.uf.m11 = newAtom.uf.m22 = newAtom.uf.m33 = uIso;
                      newAtom.uf.m32 = newAtom.uf.m23 = uIso * mole.cell.cosra;
                      newAtom.uf.m31 = newAtom.uf.m13 = uIso * mole.cell.cosrb;
                      newAtom.uf.m21 = newAtom.uf.m12 = uIso * mole.cell.cosrg;
                      mole.Uf2Uo(newAtom.uf, newAtom.uc);
                      newAtom.ufiso_org=tok.at(6);
                      newAtom.peakHeight=-666;
                      if (newAtom.an==-66) newAtom.peakHeight=9999.99;
                      newAtom.afixParent=-1;
                      if (newAtom.an>0) afixparent=mole.asymm.size();
                      else newAtom.afixParent=afixparent;
                      mole.asymm.append(newAtom);
                      isoat++;
                    }
                    if (tok.size()==8) {
                      // We have a q-peak!
                      newAtom.Label=tok.at(0);
                      if (newAtom.Label.startsWith('Q',Qt::CaseInsensitive)) {
                        newAtom.an=-1;
                        newAtom.resiNr=-1;
                        newAtom.ResiClass="Q-Peak";
                        qpeaks++;
                      } else {
                        int fac=tok.at(1).toInt()-1;
                        newAtom.an=((fac<0)||(fac>=sfac.size()))?-2:sfac.at(fac);
                        isoat++;
                      }
                      newAtom.fixFlag=0;
                      newAtom.frac.x = tok.at(2).toDouble();
                      newAtom.frac.y = tok.at(3).toDouble();
                      newAtom.frac.z = tok.at(4).toDouble();
                      if (qAbs(part_fvar) > 11.0 ){
                        // fvar is defined for entire part:
                        newAtom.sof = getNumber(part_fvar, fvar, 3, newAtom.fixFlag);
                        newAtom.sof_org = part_fvar;
                      } else {
                        // fvar is defined at each atom:
                        newAtom.sof = getNumber(tok.at(5).toDouble(), fvar, 3, newAtom.fixFlag);
                        newAtom.sof_org = tok.at(5).toDouble();
                      }
                      newAtom.isIso=true;
                      double uIso = getNumber(tok.at(6).toDouble(),fvar,uiso);
                      newAtom.uf.m11 = newAtom.uf.m22 = newAtom.uf.m33 = uIso;
                      newAtom.uf.m32 = newAtom.uf.m23 = uIso * mole.cell.cosra;
                      newAtom.uf.m31 = newAtom.uf.m13 = uIso * mole.cell.cosrb;
                      newAtom.uf.m21 = newAtom.uf.m12 = uIso * mole.cell.cosrg;
                      newAtom.ufiso_org=tok.at(6);
                      newAtom.peakHeight=tok.at(7).toDouble();
                      if (!hklf) qbeforehkl=true;
                      if(newAtom.an==-1){
                        mole.pmin=qMin(mole.pmin, newAtom.peakHeight);
                        mole.pmax=qMax(mole.pmax, newAtom.peakHeight);
                      }
                      newAtom.afixParent=-1;
                      mole.asymm.append(newAtom);
                    }
                    if (tok.size()==12){
                      // We have an ellipsoid atom!
                      if ((newAtom.part!=0)||(newAtom.resiNr!=0)){
                        newAtom.Label=QString("%1_%3%4")
                          .arg(tok.at(0))
                          .arg((newAtom.resiNr)?QString::number(newAtom.resiNr):"")
                          .arg((newAtom.part)?QString::number((newAtom.part<0)?
                                36+newAtom.part:
                                newAtom.part+9,36):"");
                      }else newAtom.Label=tok.at(0);
                      int fac=tok.at(1).toInt()-1;
                      newAtom.an=((fac<0)||(fac>=sfac.size()))?-2:sfac.at(fac);
                      newAtom.frac.x = getNumber(tok.at(2).toDouble(),fvar,0,newAtom.fixFlag);
                      newAtom.frac.y = getNumber(tok.at(3).toDouble(),fvar,1,newAtom.fixFlag);
                      newAtom.frac.z = getNumber(tok.at(4).toDouble(),fvar,2,newAtom.fixFlag);
                      if (qAbs(part_fvar) > 11.0 ){
                        newAtom.sof = getNumber(part_fvar, fvar, 3, newAtom.fixFlag);
                        newAtom.sof_org = part_fvar;
                      } else {
                        newAtom.sof = getNumber(tok.at(5).toDouble(), fvar, 3, newAtom.fixFlag);
                        newAtom.sof_org = tok.at(5).toDouble();
                      }
                      newAtom.uf.m11 = getNumber(tok.at(6).toDouble(),fvar,4,newAtom.fixFlag);
                      newAtom.uf.m22 = getNumber(tok.at(7).toDouble(),fvar,5,newAtom.fixFlag);
                      newAtom.uf.m33 = getNumber(tok.at(8).toDouble(),fvar,6,newAtom.fixFlag);
                      newAtom.uf.m32 = newAtom.uf.m23 = getNumber(tok.at(9).toDouble(),fvar,7,newAtom.fixFlag);
                      newAtom.uf.m31 = newAtom.uf.m13 = getNumber(tok.at(10).toDouble(),fvar,8,newAtom.fixFlag);
                      newAtom.uf.m21 = newAtom.uf.m12 = getNumber(tok.at(11).toDouble(),fvar,9,newAtom.fixFlag);
                      newAtom.peakHeight=-666;
                      uiso=mole.ueq(newAtom.uf);
                      if (newAtom.an>0) afixparent=mole.asymm.size();
                      newAtom.afixParent=-1;
                      mole.asymm.append(newAtom);
                      adpat++;
                    }
                  }
                  break;
        }
    }
  }
  mole.qboMax=1.67;
  mole.qbeforehkl=qbeforehkl;
  for (int i=0;i<mole.asymm.size();i++){
    mole.Uf2Uo(mole.asymm.at(i).uf,mole.asymm[i].uc);
    mole.frac2kart(mole.asymm.at(i).frac,mole.asymm[i].pos);
    mole.showatoms.append(mole.asymm.at(i));
  }
  if (mole.asymm.size()){
    mole.packer(mole.sdmcompleter());
    mole.showbonds = mole.connecting(mole.showatoms);
    int lnrr=0;
    CEnvironment belo;
    for (int bl=0; bl<bedeRloneLY.size(); bl++){
      int ignore=0;
      QStringList tok = bedeRloneLY.at(bl).split(sep,skipEmptyParts);
      if ((tok.at(0).endsWith("BEDE",Qt::CaseInsensitive))&&(tok.size()>6)){
        lnrr=belo.size();
        QString aton1=tok.at(1),aton2=tok.at(2);
        double rad = getNumber(tok.at(3).toDouble(),fvar,0,ignore);
        double a   = getNumber(tok.at(4).toDouble(),fvar,0,ignore);
        double b1  = getNumber(tok.at(5).toDouble(),fvar,0,ignore);
        double b2  = getNumber(tok.at(6).toDouble(),fvar,0,ignore);
        newAtom.an=-42;
        newAtom.symmGroup=0;
        newAtom.sg=0;
        newAtom.scod=555;//the identity
        newAtom.auidx=-1;
        int an1=-1,an2=-2;
        for (int bl1=0; bl1<mole.asymm.size(); bl1++)
          for (int bl2=0; bl2<mole.asymm.size(); bl2++){
            if ((mole.asymm.at(bl1).an<0)||(mole.asymm.at(bl2).an<0)) continue;
            double soll=-1,d= sqrt(Distance(mole.asymm.at(bl1).pos,mole.asymm.at(bl2).pos));
            if ((mole.asymm.at(bl1).an==0)&&(mole.asymm.at(bl2).an==0)) continue;
            if ((mole.asymm.at(bl1).an>-1)&&(mole.asymm.at(bl2).an>-1)&&
                ((!mole.asymm.at(bl1).part)||(!mole.asymm.at(bl2).part)||(mole.asymm.at(bl1).part==mole.asymm.at(bl2).part)))
              soll=(mole.Kovalenz_Radien[mole.asymm.at(bl1).an]+ mole.Kovalenz_Radien[mole.asymm.at(bl2).an])*0.012;
            else
              soll=-1;
            if ((d>0.1)&&(d<soll)) {
              an1=-1;
              an2=-2;
              if (aton1[0]=='$') an1=mole.getOZ(aton1);
              if (aton2[0]=='$') an2=mole.getOZ(aton2);
              if (((an1>=0)&&(mole.asymm.at(bl1).an==an1))||(aton1==mole.asymm.at(bl1).Label))
                if (((an2>=0)&&(mole.asymm.at(bl2).an==an2))||(aton2==mole.asymm.at(bl2).Label)){
                  V3 pos = mole.asymm.at(bl1).pos+Normalize(mole.asymm.at(bl2).pos-mole.asymm.at(bl1).pos)*rad;
                  newAtom.pos=pos;
                  mole.kart2frac(pos,pos);
                  newAtom.Label=QString("L%1").arg(++lnrr);
                  newAtom.frac=pos;
                  newAtom.auidx=-1;
                  newAtom.peakHeight=b1;
                  newAtom.sof=a;
                  newAtom.uf=mole.u2b(mole.asymm.at(bl1).uf);//possible problem if U is not anisotrop
                  belo.append(newAtom);
                  pos=mole.asymm.at(bl1).frac;
                  newAtom.Label=QString("L%1").arg(++lnrr);
                  newAtom.frac=mole.asymm.at(bl1).frac;
                  mole.frac2kart(newAtom.frac,newAtom.pos);
                  newAtom.sof=-a;
                  newAtom.peakHeight=b2;
                  belo.append(newAtom);

                  mole.lmin=qMin(mole.lmin,-a);
                  mole.lmax=qMax(mole.lmax, a);

                }
            }
          }

      }else
        if (tok.at(0).endsWith("LONE",Qt::CaseInsensitive)){
          belo.append(lonesome(bedeRloneLY.at(bl),belo.size()+1));
        }
    }
    if (!belo.isEmpty()) {
      for (int bl=0; bl<belo.size(); bl++){
        belo[bl].orginalLine=QString("%1    2 %2 %3 %4").arg(belo.at(bl).Label.section('_',0,0),-4)
          .arg(belo.at(bl).frac.x,10,'f',5)
          .arg(belo.at(bl).frac.y,10,'f',5)
          .arg(belo.at(bl).frac.z,10,'f',5);
      }
      mole.asymm.append(belo);
    }
    QString info;
    for (int ifv=0; ifv<fvar.size();ifv++){
      if (ifv==0) {
        info.append("OSF: ");
        info.append(QString("<font color=red>%1</font>").arg(fvar.at(ifv)));
      }
      else {
        info.append(QString("FVAR<font color=blue>%1</font>: ").arg(ifv+1));
        if (fvarCntr.value(ifv+1)<1)info.append(QString("<font color=red>%1</font>(Never Used!)").arg(fvar.at(ifv)));
        else info.append(QString("<font color=red>%1</font>(used %2 time%3.)").arg(fvar.at(ifv))
            .arg(fvarCntr.value(ifv+1)).arg((fvarCntr.value(ifv+1)>1)?"s":""));
      }
      if ((ifv+1)<fvar.size()) info.append(", ");
      info.append("<br>");
    }

  }  
  QString fcfname = s.myfcf(fileName);//compute fcf6 from hkl and structure
  loadFouAndPerform(fcfname.toLocal8Bit());
  QMap<int,double> unit;
  for (int i=0; i<mole.asymm.size(); i++){
    if (mole.asymm.at(i).an>=0) unit[mole.asymm.at(i).an]+= mole.asymm.at(i).sof * mole.cell.symmops.size();
  }
}//shelx

int main(int argc, char *argv[]){
    HKLMX=200;
    bool recheck=false;
    QCoreApplication a(argc, argv);
    QString fileName="";
    if (argc>1){
      for (int i=1;i<argc;i++){
        if ((QCoreApplication::arguments().at(i).contains(".ins",Qt::CaseInsensitive)) ||
            (QCoreApplication::arguments().at(i).contains(".res",Qt::CaseInsensitive)))          
            fileName=QCoreApplication::arguments().at(i);
        if ((QCoreApplication::arguments().at(i)=="-HKLMAX")) HKLMX = QCoreApplication::arguments().at(i+1).toInt();
        if ((QCoreApplication::arguments().at(i)=="-recheck")) recheck=true;
        }
      }
    if (fileName.isEmpty()) return 1;
    bool gut=false;
    QFile test;
    QString alltest;
    test.setFileName(fileName);
    if (test.exists()&&test.size()){
      test.open(QIODevice::ReadOnly|QIODevice::Text);
      alltest=test.readAll().replace('\0','~');// title can be empty. Fortran then has a string with many \0 this would leed to end of string in Rfactor REMARK section as the title is there

      gut = ((alltest.contains("UNIT",Qt::CaseInsensitive))&&(alltest.contains("HKLF",Qt::CaseInsensitive)||alltest.contains("END",Qt::CaseInsensitive)));
      test.close();
    }
    if (!gut){
      fileName.replace(QRegExp("res$",Qt::CaseInsensitive),"ins");
      test.setFileName(fileName);
      if (test.exists()&&test.size()){
        test.open(QIODevice::ReadOnly|QIODevice::Text);
        alltest=test.readAll();
        gut = ((alltest.contains("UNIT",Qt::CaseInsensitive))&&(alltest.contains("HKLF",Qt::CaseInsensitive)||alltest.contains("END",Qt::CaseInsensitive)));
        test.close();
      }
    }
    if (!gut) {
      test.open(QIODevice::ReadOnly|QIODevice::Text);
      alltest=test.readAll();
      test.close();
      return 2;
    }
    alltest.replace(QRegExp("=\\s*[\\r\\n]+\\s{1,}"),"==");
    resLines = alltest.split('\n');
    alltest.replace(QRegExp("REM[^\\n]*\n"),"\n");
    alltest.replace(QRegExp("![^\\n]*\n"),"\n");
    while (alltest.contains(QRegExp("[\\n^]\\++[^\\n]*\n"))) {
      QString incl=alltest.section(QRegExp("[\\n^]\\++"),1,1,QString::SectionSkipEmpty);
      incl=incl.section('\n',0,0);
      QString pre=fileName.section('/',0,-2);
      if (QFileInfo(pre+"/"+incl).exists()) {
        QFile include(pre+"/"+incl);
        QString inst;
        if (include.open(QIODevice::ReadOnly|QIODevice::Text)) inst=include.readAll();
        alltest.replace("++" + incl, inst);
        alltest.replace("+" + incl, inst);
      } else {
        alltest.remove("++" + incl);
        alltest.remove("+" + incl);
      }
    }

    load_sheldrick(fileName, alltest);
    printf("\nautoHfix ...%d %d\n",resLines.size(),mole.asymm.size());
    autoHFix();
    printf("\n?autoHfix %d %d\n",resLines.size(),mole.asymm.size());
    alltest.clear();
    alltest=resLines.join("\n");
    QString insname = fileName;
    insname.replace(QRegExp(".res$|.ins$",Qt::CaseInsensitive),".ins");
    if (recheck) load_sheldrick(insname, alltest);
    alltest=alltest.replace("==","=\n   ");
    QFile insfi(insname);
    insfi.open(QIODevice::WriteOnly|QIODevice::Text);
    insfi.write(alltest.toAscii());
    insfi.close();
}
