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
#include "molecule.h"
const double g2r=180.0/M_PI;

Molecule::Molecule(){
  qboMin=0.97;
  qboMax=1.67;
  HAWink=135;
  HAMax=2.5;
  theseAreAcceptors.append(6);//NOF
  theseAreAcceptors.append(7);
  theseAreAcceptors.append(8);

  theseAreDonors.append(5);//CNOF
  theseAreDonors.append(6);
  theseAreDonors.append(7);
  theseAreDonors.append(8);

  Kovalenz_Radien[0  ]=32  ;Kovalenz_Radien[1  ]=1   ;Kovalenz_Radien[2  ]=123 ;Kovalenz_Radien[3  ]=90  ;
  Kovalenz_Radien[4  ]=80  ;Kovalenz_Radien[5  ]=77  ;Kovalenz_Radien[6  ]=74  ;Kovalenz_Radien[7  ]=66  ;
  Kovalenz_Radien[8  ]=72  ;Kovalenz_Radien[9  ]=1   ;Kovalenz_Radien[10 ]=154 ;Kovalenz_Radien[11 ]=149 ;
  Kovalenz_Radien[12 ]=118 ;Kovalenz_Radien[13 ]=111 ;Kovalenz_Radien[14 ]=106 ;Kovalenz_Radien[15 ]=102 ;
  Kovalenz_Radien[16 ]=99  ;Kovalenz_Radien[17 ]=1   ;Kovalenz_Radien[18 ]=203 ;Kovalenz_Radien[19 ]=174 ;
  Kovalenz_Radien[20 ]=144 ;Kovalenz_Radien[21 ]=132 ;Kovalenz_Radien[22 ]=122 ;Kovalenz_Radien[23 ]=118 ;
  Kovalenz_Radien[24 ]=117 ;Kovalenz_Radien[25 ]=117 ;Kovalenz_Radien[26 ]=116 ;Kovalenz_Radien[27 ]=124 ;
  Kovalenz_Radien[28 ]=117 ;Kovalenz_Radien[29 ]=125 ;Kovalenz_Radien[30 ]=126 ;Kovalenz_Radien[31 ]=122 ;
  Kovalenz_Radien[32 ]=120 ;Kovalenz_Radien[33 ]=116 ;Kovalenz_Radien[34 ]=114 ;Kovalenz_Radien[35 ]=1   ;
  Kovalenz_Radien[36 ]=218 ;Kovalenz_Radien[37 ]=191 ;Kovalenz_Radien[38 ]=162 ;Kovalenz_Radien[39 ]=145 ;
  Kovalenz_Radien[40 ]=134 ;Kovalenz_Radien[41 ]=130 ;Kovalenz_Radien[42 ]=127 ;Kovalenz_Radien[43 ]=125 ;
  Kovalenz_Radien[44 ]=125 ;Kovalenz_Radien[45 ]=128 ;Kovalenz_Radien[46 ]=134 ;Kovalenz_Radien[47 ]=148 ;
  Kovalenz_Radien[48 ]=144 ;Kovalenz_Radien[49 ]=141 ;Kovalenz_Radien[50 ]=140 ;Kovalenz_Radien[51 ]=136 ;
  Kovalenz_Radien[52 ]=133 ;Kovalenz_Radien[53 ]=1   ;Kovalenz_Radien[54 ]=235 ;Kovalenz_Radien[55 ]=198 ;
  Kovalenz_Radien[56 ]=169 ;Kovalenz_Radien[57 ]=165 ;Kovalenz_Radien[58 ]=165 ;Kovalenz_Radien[59 ]=164 ;
  Kovalenz_Radien[60 ]=164 ;Kovalenz_Radien[61 ]=162 ;Kovalenz_Radien[62 ]=185 ;Kovalenz_Radien[63 ]=161 ;
  Kovalenz_Radien[64 ]=159 ;Kovalenz_Radien[65 ]=159 ;Kovalenz_Radien[66 ]=157 ;Kovalenz_Radien[67 ]=157 ;
  Kovalenz_Radien[68 ]=156 ;Kovalenz_Radien[69 ]=170 ;Kovalenz_Radien[70 ]=156 ;Kovalenz_Radien[71 ]=144 ;
  Kovalenz_Radien[72 ]=134 ;Kovalenz_Radien[73 ]=130 ;Kovalenz_Radien[74 ]=128 ;Kovalenz_Radien[75 ]=126 ;
  Kovalenz_Radien[76 ]=127 ;Kovalenz_Radien[77 ]=130 ;Kovalenz_Radien[78 ]=134 ;Kovalenz_Radien[79 ]=149 ;
  Kovalenz_Radien[80 ]=148 ;Kovalenz_Radien[81 ]=147 ;Kovalenz_Radien[82 ]=146 ;Kovalenz_Radien[83 ]=146 ;
  Kovalenz_Radien[84 ]=145 ;Kovalenz_Radien[85 ]=1   ;Kovalenz_Radien[86 ]=0   ;Kovalenz_Radien[87 ]=1   ;
  Kovalenz_Radien[88 ]=188 ;Kovalenz_Radien[89 ]=165 ;Kovalenz_Radien[90 ]=161 ;Kovalenz_Radien[91 ]=142 ;
  Kovalenz_Radien[92 ]=130 ;Kovalenz_Radien[93 ]=151 ;Kovalenz_Radien[94 ]=182 ;
  PSE<<"H"<<"He"<<"Li"<<"Be"<<"B"<<"C"<<"N"<<"O"<<"F"<<"Ne"<<"Na"<<"Mg"<<"Al"<<"Si"<<"P"<<"S"<<"Cl"<<"Ar"<<
    "K"<<"Ca"<<"Sc"<<"Ti"<<"V"<<"Cr"<<"Mn"<<"Fe"<<"Co"<<"Ni"<<"Cu"<<"Zn"<<"Ga"<<"Ge"<<"As"<<"Se"<<"Br"<<"Kr"<<
    "Rb"<<"Sr"<<"Y"<<"Zr"<<"Nb"<<"Mo"<<"Tc"<<"Ru"<<"Rh"<<"Pd"<<"Ag"<<"Cd"<<"In"<<"Sn"<<"Sb"<<"Te"<<"I"<<"Xe"<<
    "Cs"<<"Ba"<< "La"<<"Ce"<<"Pr"<<"Nd"<<"Pm"<<"Sm"<<"Eu"<<"Gd"<<"Tb"<<"Dy"<<"Ho"<<"Er"<<"Tm"<<"Yb"<<"Lu"<<
    "Hf"<<"Ta"<<"W"<<"Re"<<"Os"<<"Ir"<<"Pt"<<"Au"<<"Hg"<<"Tl"<<"Pb"<<"Bi"<<"Po"<<"At"<<"Rn"<<"Fr"<<"Ra"<<
    "Ac"<<"Th"<<"Pa"<<"U"<<"Np"<<"Pu"<<"Am"<<"Cm"<<"Bk"<<"Cf"<<"Es"<<"Fm"<<"Md"<<"No"<<"Lr";//no need of these I guess--> <<"Ku"<<"Ha"<<"Rf"<<"Ns"<<"Hs"<<"Mt";
}

int Molecule::afix5(int idx){
  int f5 = 0,m,cnt = 0,aafx=-1;
  if (showatoms.at(idx).afix==0) return 0;
  printf("AFIXm5 %d (%s)\n",showatoms.at(idx).afix,showatoms.at(idx).Label.toStdString().c_str());
  for (int i = 0; (i <= idx) && (i < showatoms.size()); i++){
    if (showatoms.at(i).an > 0){
      if (aafx != showatoms.at(i).afix) cnt = 0;// afix is different than before => reset the counter
      aafx = showatoms.at(i).afix;
      m = showatoms.at(i).afix / 10;
      switch (m) {
        case 0://afix 6 afix 9
          printf("ehm m is zero?\n");
          f5=5;
          break;
        case 6://66 for phenyl rings 
          cnt++;
          f5 = (cnt % 6)? m * 10 + 5: 0;
          //        printf("cnt%d %s\n",cnt,showatoms.at(i).Label.toStdString().c_str());
          break;
        case 5://56 planar regular 5 membered ring (cp)
          cnt++;
          f5 = (cnt % 5)? m * 10 + 5: 0;
          break;
        case 10://cp*
          cnt++;
          f5 = (cnt % 10)? m * 10 + 5: 0;
          break;
        case 11://naphtyl
          cnt++;
          f5 = (cnt % 10)? m * 10 + 5: 0;
          break;
      }
    }
  }
  //  printf("==> %d aafx%d m%d cnt%d idx%d %s\n",f5,aafx,m,cnt,idx,showatoms.at(idx).Label.toStdString().c_str());
  return f5;
}

double Molecule::dieder(V3 a,V3 b, V3 c){

  static double erg;
  double A[3],B[3],sig;
  sig=a.x*b.y*c.z - a.z*b.y*c.x +
    a.z*b.x*c.y - a.x*b.z*c.y +
    a.y*b.z*c.x - a.y*b.x*c.z;

  A[0]= a.y*b.z - a.z*b.y;
  A[1]=-a.x*b.z + a.z*b.x;
  A[2]= a.x*b.y - a.y*b.x;

  B[0]=-b.y*c.z + b.z*c.y;
  B[1]= b.x*c.z - b.z*c.x;
  B[2]=-b.x*c.y + b.y*c.x;
  //  printf("A%f B%f\n",sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]),sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]));
  erg=(A[0]*B[0]+A[1]*B[1]+A[2]*B[2])/(sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2])*sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]));
  erg=acos(erg)/M_PI*180.0;
  return (sig>0)?(erg):(-erg);
}

double Molecule::winkel(V3 a,V3 b){

  static double erg;
  if ((Norm(a)<0.001)||(Norm(b)<0.001)) return 0;
  erg= a*b /(sqrt(Norm(a))*sqrt(Norm(b)));
  //  erg=(a.x*b.x+a.y*b.y+a.z*b.z)/(sqrt(a.x*a.x+a.y*a.y+a.z*a.z)*sqrt(b.x*b.x+b.y*b.y+b.z*b.z));
  erg=acos(erg)/M_PI*180.0;
  return(erg);
}

V3 Molecule::kreuzX(double x1,double y1,double z1,double x2,double y2,double z2) {

  static V3 erg;
  erg.x=-y1*z2+z1*y2;
  erg.y=-z1*x2+x1*z2;
  erg.z=-x1*y2+y1*x2;
  return(erg);
}


bool Molecule::decodeSymmCard(QString symmCard){ 
  /*! decodes a symmetry card like 'SYMM -X, 1/2+Y, -Z' and feeds cell.symmops and cell.trans lists.
   * @param symmCard like 'SYMM -X, 1/2+Y, -Z'.
   * \returns true on sucess.
   * */ 
  QString sc=symmCard.toUpper().remove("SYMM").trimmed();
  sc=sc.remove("'");
  sc=sc.remove(" ");
  QStringList axe=sc.split(",");
  QStringList bruch;
  if (axe.size()!=3) return false;
  double sx[3],sy[3],sz[3],t[3];
  for (int i=0; i<3; i++){
    sx[i]=0;sy[i]=0;sz[i]=0;t[i]=0;
    if (axe.at(i).contains("-X")) {sx[i]=-1.0;axe[i].remove("-X");}
    else if (axe.at(i).contains("X")) {sx[i]=1.0;axe[i].remove("X");}
    if (axe.at(i).contains("-Y")) {sy[i]=-1.0;axe[i].remove("-Y");}
    else if (axe.at(i).contains("Y")) {sy[i]=1.0;axe[i].remove("Y");}
    if (axe.at(i).contains("-Z")) {sz[i]=-1.0;axe[i].remove("-Z");}
    else if (axe.at(i).contains("Z")) {sz[i]=1.0;axe[i].remove("Z");}
    if (axe.at(i).endsWith("+")) axe[i].remove("+");
    if (axe.at(i).contains("/")) {
      bruch=axe.at(i).split("/");
      if (bruch.size()==2) t[i]=bruch.at(0).toDouble() / bruch.at(1).toDouble();
    }
    else if (!axe.at(i).isEmpty()) t[i]=axe.at(i).toDouble();
  }
  Matrix sm = Matrix(sx[0],sy[0],sz[0],	  sx[1],sy[1],sz[1],  sx[2],sy[2],sz[2]);

  cell.symmops.append(sm);
  cell.trans.append(V3(t[0],t[1],t[2]));
  return true;
}

bool Molecule::decodeSymmCardEQIV(QString symmCard){
  /*! decodes a symmetry card like 'SYMM -X, 1/2+Y, -Z' and feeds Molecule.symmopsEQIV, Molecule.transEQIV and Molecule.labelEQIV.
   * \returns true on success.
   * @param symmCard like 'SYMM -X, 1/2+Y, -Z'.
   */
  QString sc=symmCard.toUpper().remove("EQIV").trimmed();
  labelEQIV.append(sc.section(' ',0,0));
  sc.remove(labelEQIV.last());
  sc=sc.remove("'");
  sc=sc.remove(" ");
  QStringList axe=sc.split(",");
  QStringList bruch;
  if (axe.size()!=3) return false;
  double sx[3],sy[3],sz[3],t[3];
  for (int i=0; i<3; i++){
    sx[i]=0;sy[i]=0;sz[i]=0;t[i]=0;
    if (axe.at(i).contains("-X")) {sx[i]=-1.0;axe[i].remove("-X");}
    else if (axe.at(i).contains("X")) {sx[i]=1.0;axe[i].remove("X");}
    if (axe.at(i).contains("-Y")) {sy[i]=-1.0;axe[i].remove("-Y");}
    else if (axe.at(i).contains("Y")) {sy[i]=1.0;axe[i].remove("Y");}
    if (axe.at(i).contains("-Z")) {sz[i]=-1.0;axe[i].remove("-Z");}
    else if (axe.at(i).contains("Z")) {sz[i]=1.0;axe[i].remove("Z");}
    if (axe.at(i).endsWith("+")) axe[i].remove("+");
    if (axe.at(i).contains("/")) {
      bruch=axe.at(i).split("/");
      if (bruch.size()==2) t[i]=bruch.at(0).toDouble() / bruch.at(1).toDouble();
    }
    else if (!axe.at(i).isEmpty()) t[i]=axe.at(i).toDouble();
  }
  Matrix sm = Matrix(sx[0],sy[0],sz[0],	  sx[1],sy[1],sz[1],  sx[2],sy[2],sz[2]);
  symmopsEQIV.append(sm);
  transEQIV.append(V3(t[0],t[1],t[2]));
  //  qDebug()<<sx[0]<<sy[0]<<sz[0]<<   sx[1]<<sy[1]<<sz[1]<<  sx[2]<<sy[2]<<sz[2]<<t[0]<<t[1]<<t[2];
  return true;
}

QString Molecule::symmCard2Code(QString symmCard){
  /*! Creates a internal symmetry code from a human readable symmetry card.
   * @param symmCard human readable symmetry card like "-X, 1/2+Y, -Z"
   * \returns internal symmetry code like 'n_555:1'
   */
  QString sc=symmCard.toUpper().trimmed();
  sc=sc.remove("'");
  sc=sc.remove(" ");
  QStringList axe=sc.split(",");
  QStringList bruch;
  if (axe.size()!=3) return "";
  double sx[3],sy[3],sz[3],t[3];
  for (int i=0; i<3; i++){
    sx[i]=0;sy[i]=0;sz[i]=0;t[i]=0;
    if (axe.at(i).contains("-X")) {sx[i]=-1.0;axe[i].remove("-X");}
    else if (axe.at(i).contains("X")) {sx[i]=1.0;axe[i].remove("X");}
    if (axe.at(i).contains("-Y")) {sy[i]=-1.0;axe[i].remove("-Y");}
    else if (axe.at(i).contains("Y")) {sy[i]=1.0;axe[i].remove("Y");}
    if (axe.at(i).contains("-Z")) {sz[i]=-1.0;axe[i].remove("-Z");}
    else if (axe.at(i).contains("Z")) {sz[i]=1.0;axe[i].remove("Z");}
    if (axe.at(i).endsWith("+")) axe[i].remove("+");
    if (axe.at(i).contains("/")) {
      bruch=axe.at(i).split("/");
      if (bruch.size()==2) t[i]=bruch.at(0).toDouble() / bruch.at(1).toDouble();
    }
    else if (!axe.at(i).isEmpty()) t[i]=axe.at(i).toDouble();
  }
  Matrix sm = Matrix(sx[0],sy[0],sz[0],	  sx[1],sy[1],sz[1],  sx[2],sy[2],sz[2]);
  if (!cell.symmops.contains(sm)) return "";
  V3 r;
  r.x=fmod(t[0]+10,1.0);
  r.y=fmod(t[1]+10,1.0);
  r.z=fmod(t[2]+10,1.0);
  for (int i=0; i<cell.symmops.size();i++){
    if ((cell.symmops.at(i)==sm)&&(cell.trans.at(i)==r)){
      QString ss=QString("%1_%2%3%4:1,")
        .arg(i+1)
        .arg(static_cast<int>(t[0]-r.x) +5)
        .arg(static_cast<int>(t[1]-r.y) +5)
        .arg(static_cast<int>(t[2]-r.z) +5);
      return ss;
    }
  }
  return "";
}

double Molecule::shortestDistance(QString sc){
  /*! computes the shortest distance for the given internal symmetry code.
  */
  double erg=100000;
  V3 frac;
  int s,h,k,l,symmgroup;
  sscanf(sc.toLatin1(),"%d_%1d%1d%1d:%d",&s,&h,&k,&l,&symmgroup);
  //printf("BS:!%s! %d h%d k%d l%d sg%d\n",brauchSymm.at(j).toStdString().c_str(),s,h,k,l,symmgroup);
  h-=5;
  k-=5;
  l-=5;
  s--;
  for (int i=0;i<asymm.size();i++){
    if ((asymm[i].molindex==symmgroup)&&(asymm[i].an>-1)){
      frac=cell.symmops.at(s)*asymm[i].frac+cell.trans.at(s)+V3(h,k,l);
      for (int j=0;j<asymm.size();j++){
        if (asymm.at(j).an<0)continue;
        erg=qMin(erg,fl(frac.x-asymm.at(j).frac.x,frac.y-asymm.at(j).frac.y,frac.z-asymm.at(j).frac.z));
      }
    }
  }
  return erg;
}

int commonFaces(Vert a,Vert b, int f[2]){
  int matches=0; 
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      if (a.faces[i]==b.faces[j])f[matches++]=a.faces[i];
  return matches;
}

void Molecule::voronoij(CEnvironment au, int intat){
  double dk,range=5.0;
  voroMsg.clear();
  voroMsg.append("<b>Voronoi:</b><br>");
  QTime speedTest;
  speedTest.start();
  vtriangles.clear();
  VPoly triangle;
  V3 prime,dp,D,floorD;
  QList<SdmItem> sdm;
  SdmItem sdmItem;
  sdmItem.a1=0;
  sdmItem.a2=1;
  sdmItem.sn=0;
  //sdmItem.p1=sdmItem.p2=V3(0,0,0);
  for (int i=0; i<au.size(); i++){ 
    for (int j=0; j<au.size(); j++ ){
      if ((au.at(i).sg)||(au.at(j).sg)) continue;
      //  bool hma=false;
      for (int n=0;n<cell.symmops.size();  n++){
        prime=cell.symmops.at(n) * au.at(i).frac + cell.trans.at(n);
        D=prime - au.at(j).frac+ V3(0.5,0.5,0.5) ;
        floorD=V3(floor(D.x),floor(D.y),floor(D.z));
        for (int h=-1; h<2; h++){
          for (int k=-1; k<2; k++){
            for (int l=-1; l<2; l++){
              V3 fD=floorD+V3(h,k,l);  
              dp=D - fD - V3(0.5,0.5,0.5);
              dk=fl(dp.x,dp.y,dp.z);
              if ((dk>0.01)&&((range)>=dk)){
                sdmItem.d=dk;
                sdmItem.floorD=fD;
                sdmItem.a1=i;
                sdmItem.a2=j;
                sdmItem.sn=n;
                sdm.append(sdmItem);
              }
            }
          }
        }
      }
    }
  }
  qSort(sdm.begin(),sdm.end());
  //SDMprint(sdm,au);
  QList<V3> n,m,doneLine;
  QList<int> intra, altcol;
  QList<Vert> v;
  Vert vert;
  V3 pf,pc,mx,nx,of,oc;
  Matrix mat;
  double vol=0.0,avol,tvol;
  for (int i=(intat>-1)?intat:0; i<au.size(); i++){
    if (au.at(i).an<0) continue;
    //find faces
    n.clear();
    m.clear();
    v.clear();
    doneLine.clear();
    intra.clear();
    avol=0.0;
    int faci=0;
    //int vcnt=0;
    frac2kart(au.at(i).frac,oc);
    int tris=0;
    for (int j=0; j<sdm.size();j++){
      if(sdm.at(j).d<0.1) continue;
      if ((sdm.at(j).a2==i)&&(au[sdm.at(j).a1].an>-1)) {
        pf = cell.symmops.at(sdm.at(j).sn) * au[sdm.at(j).a1].frac + cell.trans.at(sdm.at(j).sn) - sdm.at(j).floorD;
        frac2kart(pf,pc);
        /*
           printf("%9.5f %9.5f %9.5f  %9.5f %9.5f %9.5f\n"
           ,oc.x
           ,oc.y
           ,oc.z
           ,pc.x
           ,pc.y
           ,pc.z
           );// */
        mx =0.5*(oc+pc);
        nx = Normalize(pc-oc);
        bool exists = false; 
        double dpl;
        for (int tp=0; tp<m.size(); tp++){
          dpl=n.at(tp)*(mx-m.at(tp));
          if (fabs(dpl)<0.001){exists=true;break;}
        }
        if (!exists){
          m.append(mx);
          n.append(nx);
          intra.append((sdm.at(j).sn==0)&&(Norm(sdm.at(j).floorD)==0.0));
          if (intat>-1)altcol.append(au[sdm.at(j).a1].an);
        }
      }//is neighbour of i
    }// j sdm
    for (int h=0; h<n.size()-2; h++)
      for (int k=h+1; k<n.size()-1; k++)
        for (int l=k+1; l<n.size(); l++){
          mat=Matrix(n.at(h),n.at(k),n.at(l));
          double d=determinant(mat),dpl;
          //printf("%d,%d,%d %9.5f\n",h,k,l,d);
          if (fabs(d)<0.0002) continue;

          //p = (dot(p1,n1)*cross(n2,n3)-dot(p2,n2)*cross(n1,n3)+dot(p3,n3)*cross(n1,n2))/d 
          bool out=false;
          vert.pos=(
              (m.at(h)*n.at(h))*(n.at(k)%n.at(l))+
              (m.at(k)*n.at(k))*(n.at(l)%n.at(h))+
              (m.at(l)*n.at(l))*(n.at(h)%n.at(k)))*(1.0/d);
          vert.faces[0]=h;
          vert.faces[1]=k;
          vert.faces[2]=l;
          for (int pli=0; pli<n.size(); pli++){//teste lage zu allen ebenen
            dpl=n.at(pli)*(vert.pos-m.at(pli));
            if (dpl>0.00001) {
              out=true;
              break;
            }
          }
          if (!out) {
            v.append(vert);
            //  printf("[%d,%d,%d] %9.5f%9.5f%9.5f  \n",h,k,l,vert.pos.x,vert.pos.y,vert.pos.z);
          }
        }
    for (int mi=0; mi<m.size();mi++){// this is to get the midpoint of a face in its actual center
      V3 mnew=V3(0,0,0);
      int poli=0;
      for (int vi=0; vi<v.size(); vi++){
        if ((v.at(vi).faces[0]==mi)||(v.at(vi).faces[1]==mi)||(v.at(vi).faces[2]==mi)){
          mnew+=v.at(vi).pos;
          poli++;
        }
      }
      mnew*=1.0/poli;
      //      printf("%dsoll 0=%f\n",poli,n.at(mi)*(mnew-m.at(mi)));
      if (poli) {m[mi]=mnew;faci++;}
    }
    /*    glBegin(GL_LINES);
          for(int vi=0; vi<v.size()-1; vi++)
          for(int vj=vi+1; vj<v.size(); vj++){
          int fc[2];
          if (commonFaces(v.at(vi),v.at(vj),fc)==2){
          glColor4fv(Acol[au[i].an]);
          glVertex3d(v.at(vi).pos.x,v.at(vi).pos.y,v.at(vi).pos.z);
          glVertex3d(v.at(vj).pos.x,v.at(vj).pos.y,v.at(vj).pos.z);
          }  
          }
          glEnd();*/
    //      glBegin(GL_TRIANGLES);
    //glBegin(GL_LINES);

    for(int vi=0; vi<v.size()-1; vi++)
      for(int vj=vi+1; vj<v.size(); vj++){
        int fc[2];
        if (commonFaces(v.at(vi),v.at(vj),fc)==2){
          V3 lin=v.at(vi).pos+v.at(vj).pos;
          double dis;
          bool done=false;
          for (int di=0; di<doneLine.size(); di++){
            dis=Distance(lin,doneLine.at(di));
            if (dis<0.000001) {done=true; break;}
          }
          if (done) continue;
          doneLine.append(lin);
          tvol=determinant(Matrix(m.at(fc[0])-oc, v.at(vi).pos-oc,v.at(vj).pos-oc));
          int v1=(tvol<0)?vj:vi,
              v2=(tvol<0)?vi:vj;
          triangle.mid=m.at(fc[0]);
          triangle.nor=n.at(fc[0]);
          triangle.acol=(intat>-1)?altcol.at(fc[0]):au[i].an;
          triangle.intra=(intat>-1)?0:intra.at(fc[0]);
          triangle.verts0=v[v1].pos;
          triangle.verts1=v[v2].pos;

          vtriangles.append(triangle);
          /*
             GLfloat colo[4]={ Acol[au[i].an][0], Acol[au[i].an][1], Acol[au[i].an][2], 0.4};
          //GLfloat colo[4]={ Acol[tris%100][0], Acol[tris%100][1], Acol[tris%100][2], 0.4};
          colo[3]=intra.at(fc[0])?1.0:0.4;
          glColor4fv(colo);
          glNormal3d(n.at(fc[0]).x,n.at(fc[0]).y,n.at(fc[0]).z);
          glVertex3d(m.at(fc[0]).x,m.at(fc[0]).y,m.at(fc[0]).z);
          glVertex3d(v.at(v1).pos.x,v.at(v1).pos.y,v.at(v1).pos.z);
          glVertex3d(v.at(v2).pos.x,v.at(v2).pos.y,v.at(v2).pos.z);
          */
          tris++;
          double tv=fabs(tvol)/6.0;
          avol+=tv;
          vol+=tv*au[i].sof;
          tvol=determinant(Matrix(m.at(fc[1])-oc, v.at(vi).pos-oc,v.at(vj).pos-oc));
          v1=(tvol<0)?vj:vi;
          v2=(tvol<0)?vi:vj;
          triangle.mid=m.at(fc[1]);
          triangle.nor=n.at(fc[1]);

          triangle.acol=(intat>-1)?altcol.at(fc[1]):au[i].an;
          triangle.intra=(intat>-1)?0:intra.at(fc[1]);
          //         triangle.col[3]=((intat>-1)&&(intra.at(fc[1])))?0.8:0.4;
          triangle.verts0=v[v1].pos;
          triangle.verts1=v[v2].pos;
          vtriangles.append(triangle);
          /*
             colo[3]=intra.at(fc[1])?1.0:0.4;
             glColor4fv(colo);

             glNormal3d(n.at(fc[1]).x,n.at(fc[1]).y,n.at(fc[1]).z);
             glVertex3d(m.at(fc[1]).x,m.at(fc[1]).y,m.at(fc[1]).z);
             glVertex3d(v.at(v1).pos.x,v.at(v1).pos.y,v.at(v1).pos.z);
             glVertex3d(v.at(v2).pos.x,v.at(v2).pos.y,v.at(v2).pos.z);
             */
          tris++;
          tv=fabs(tvol)/6.0;
          avol+=tv;
          vol+=tv*au[i].sof;
        }  
      } 
    //    glEnd();
    //    printf("%d neighbours found for %s. Voronoij polyeder has %d verices, %d edges and %d faces. Triangles %d Volume %f Total Volume= %f cell volume %f\n",
    //        m.size(),au.at(i).Label,v.size(),v.size()+m.size()-2,m.size(),
    printf("%-12s: Triangles %6d Volume %18.5f %d\n",au.at(i).Label.toStdString().c_str(),tris,avol,-tris/2+faci+v.size());
    voroMsg.append(QString("%1 : EulerTest: %2 Volume %3 &Aring;<sup>3</sup>.<br>").arg(au.at(i).Label).arg(-tris/2+faci+v.size()).arg(avol));

    if (intat>-1) break;
  }//i atoms au
  if (intat==-1) {
    printf("Total Volume= %18.5f Total Volume * z = %18.5f cell volume %18.5f delta =%f%%\n",vol,vol*cell.symmops.size(),cell.V,(vol*cell.symmops.size()-cell.V)/cell.V*100.0);
    voroMsg.append(QString("Total Volume= %1 &Aring;<sup>3</sup> cell volume %2 &Aring;<sup>3</sup> &Delta;V %3%<br>").arg(vol*cell.symmops.size()).arg(cell.V).arg((vol*cell.symmops.size()-cell.V)/cell.V*100.0));
  }
  //  printf("drawing %dms\n",speedTest.restart());
}
V3 eye;
bool vbyeye(VPoly a,VPoly b){
  return (Distance(a.mid,eye)>Distance(b.mid,eye));
}
void Molecule::setHBondMaxDist(double d){if ((d>0.0)&&(d<4.0)) HAMax=d;}
void Molecule::setHBondMaxAngl(double w){if ((w>0.0)&&(w<180.0)) HAWink=w;}

bool Molecule::canbeDonor(int an){
  return theseAreDonors.contains(an);
}

bool Molecule::canbeAcceptor(int an){
  return theseAreAcceptors.contains(an);
}
//QList<int> theseAreAcceptors;//
//QList<int> theseAreDonors;
double Molecule::hbdist(){return HAMax;}
double Molecule::hbangl(){return HAWink;}
#include <stdio.h>
QString Molecule::symmcode2human(int s){
  /*! @param s the Nth symmetry operation of the space group.
   *  \returns human readable string. 
   */
  QString erg;
  Matrix m;
  V3 t;
  int h=0,k=0,l=0;
  m=cell.symmops.at(s);
  QString symstrX,symstrY,symstrZ;
  if (m.m11==1) symstrX.append("+x");
  if (m.m11==-1) symstrX.append("-x");
  if (m.m21==1) symstrX.append("+y");
  if (m.m21==-1) symstrX.append("-y");
  if (m.m31==1) symstrX.append("+z");
  if (m.m31==-1) symstrX.append("-z");

  if (m.m12==1) symstrY.append("+x");
  if (m.m12==-1) symstrY.append("-x");
  if (m.m22==1) symstrY.append("+y");
  if (m.m22==-1) symstrY.append("-y");
  if (m.m32==1) symstrY.append("+z");
  if (m.m32==-1) symstrY.append("-z");

  if (m.m13==1) symstrZ.append("+x");
  if (m.m13==-1) symstrZ.append("-x");
  if (m.m23==1) symstrZ.append("+y");
  if (m.m23==-1) symstrZ.append("-y");
  if (m.m33==1) symstrZ.append("+z");
  if (m.m33==-1) symstrZ.append("-z");
  t=cell.trans.at(s);
  V3 zaehler,nenner;
  double egal;
  for (int g=1;g<7; g++){
    nenner.x=(t.x<0)?-g:g;
    zaehler.x=static_cast<int>(round((t.x+h)*g));
    if (fabs(modf(t.x*g,&egal))<0.05) break;
  }
  for (int g=1;g<7; g++){
    nenner.y=(t.y<0)?-g:g;
    zaehler.y=static_cast<int>(round((t.y+k)*g));
    if (fabs(modf(t.y*g,&egal))<0.05) break;
  }
  for (int g=1;g<7; g++){
    nenner.z=(t.z<0)?-g:g;
    zaehler.z=static_cast<int>(round((t.z+l)*g));
    if (fabs(modf(t.z*g,&egal))<0.05) break;
  }
  erg=(QString(" %1/%2%3, %4/%5%6, %7/%8%9")
      .arg(zaehler.x)
      .arg(nenner.x)
      .arg(symstrX)
      .arg(zaehler.y)
      .arg(nenner.y)
      .arg(symstrY)
      .arg(zaehler.z)
      .arg(nenner.z)
      .arg(symstrZ));
  erg.remove(QRegExp("0/\\d"));
  erg.replace("1/1","1");
  erg.replace("2/1","2");
  erg.replace("3/1","3");
  erg.replace("4/1","4");
  erg.replace("5/1","5");
  erg.replace("6/1","6");
  return erg;
}

QString Molecule::symmcode2human(QString brauchSymm){
  /*! @param brauchSymm internal symmetry code n_555
   *  \returns human readable string 
   */
  QString erg;
  Matrix m;
  V3 t;
  int h,k,l,s;
  sscanf(brauchSymm.toLatin1(),"%d_%1d%1d%1d",&s,&h,&k,&l);
  h-=5;
  k-=5;
  l-=5;
  s--;
  m=cell.symmops.at(s);
  QString symstrX,symstrY,symstrZ;
  if (m.m11==1) symstrX.append("+x");
  if (m.m11==-1) symstrX.append("-x");
  if (m.m21==1) symstrX.append("+y");
  if (m.m21==-1) symstrX.append("-y");
  if (m.m31==1) symstrX.append("+z");
  if (m.m31==-1) symstrX.append("-z");

  if (m.m12==1) symstrY.append("+x");
  if (m.m12==-1) symstrY.append("-x");
  if (m.m22==1) symstrY.append("+y");
  if (m.m22==-1) symstrY.append("-y");
  if (m.m32==1) symstrY.append("+z");
  if (m.m32==-1) symstrY.append("-z");

  if (m.m13==1) symstrZ.append("+x");
  if (m.m13==-1) symstrZ.append("-x");
  if (m.m23==1) symstrZ.append("+y");
  if (m.m23==-1) symstrZ.append("-y");
  if (m.m33==1) symstrZ.append("+z");
  if (m.m33==-1) symstrZ.append("-z");
  t=cell.trans.at(s);
  V3 zaehler,nenner;
  double egal;
  for (int g=1;g<7; g++){
    nenner.x=(t.x<0)?-g:g;
    zaehler.x=static_cast<int>(round((t.x+h)*g));
    if (fabs(modf(t.x*g,&egal))<0.05) break;
  }
  for (int g=1;g<7; g++){
    nenner.y=(t.y<0)?-g:g;
    zaehler.y=static_cast<int>(round((t.y+k)*g));
    if (fabs(modf(t.y*g,&egal))<0.05) break;
  }
  for (int g=1;g<7; g++){
    nenner.z=(t.z<0)?-g:g;
    zaehler.z=static_cast<int>(round((t.z+l)*g));
    if (fabs(modf(t.z*g,&egal))<0.05) break;
  }
  erg=(QString("%1/%2%3, %4/%5%6, %7/%8%9")
      .arg(zaehler.x)
      .arg(nenner.x)
      .arg(symstrX)
      .arg(zaehler.y)
      .arg(nenner.y)
      .arg(symstrY)
      .arg(zaehler.z)
      .arg(nenner.z)
      .arg(symstrZ));


  erg.remove(QRegExp("0/\\d"));
  erg.replace("1/1","1");
  erg.replace("2/1","2");
  erg.replace("3/1","3");
  erg.replace("4/1","4");
  erg.replace("5/1","5");
  erg.replace("6/1","6");
  return erg;
}

QString Molecule::symmcode2human(QString brauchSymm, int j){
  /*! @param brauchSymm internal symmetry code n_555
   *  @param j symmetry number to prepend the human readable string returned.
   *  \returns human readable string 
   */
  QString erg;
  Matrix m;
  V3 t;
  int h,k,l,s;
  sscanf(brauchSymm.toLatin1(),"%d_%1d%1d%1d",&s,&h,&k,&l);
  h-=5;
  k-=5;
  l-=5;
  s--;
  m=cell.symmops.at(s);
  QString symstrX,symstrY,symstrZ;
  if (m.m11==1) symstrX.append("+x");
  if (m.m11==-1) symstrX.append("-x");
  if (m.m21==1) symstrX.append("+y");
  if (m.m21==-1) symstrX.append("-y");
  if (m.m31==1) symstrX.append("+z");
  if (m.m31==-1) symstrX.append("-z");

  if (m.m12==1) symstrY.append("+x");
  if (m.m12==-1) symstrY.append("-x");
  if (m.m22==1) symstrY.append("+y");
  if (m.m22==-1) symstrY.append("-y");
  if (m.m32==1) symstrY.append("+z");
  if (m.m32==-1) symstrY.append("-z");

  if (m.m13==1) symstrZ.append("+x");
  if (m.m13==-1) symstrZ.append("-x");
  if (m.m23==1) symstrZ.append("+y");
  if (m.m23==-1) symstrZ.append("-y");
  if (m.m33==1) symstrZ.append("+z");
  if (m.m33==-1) symstrZ.append("-z");
  t=cell.trans.at(s);
  V3 zaehler,nenner;
  double egal;
  for (int g=1;g<7; g++){
    nenner.x=(t.x<0)?-g:g;
    zaehler.x=static_cast<int>(round((t.x+h)*g));//+0.0001
    if (fabs(modf(t.x*g,&egal))<0.05) break;
  }
  for (int g=1;g<7; g++){
    nenner.y=(t.y<0)?-g:g;
    zaehler.y=static_cast<int>(round((t.y+k)*g));//+0.001
    if (fabs(modf(t.y*g,&egal))<0.05) break;
  }
  for (int g=1;g<7; g++){
    nenner.z=(t.z<0)?-g:g;
    zaehler.z=static_cast<int>(round((t.z+l)*g));//+0.001
    if (fabs(modf(t.z*g,&egal))<0.05) break;
  }

  erg=(QString("&raquo;%1 :<b> %2/%3%4, %5/%6%7, %8/%9%10</b><br>")
      .arg(j)
      .arg(zaehler.x)
      .arg(nenner.x)
      .arg(symstrX)
      .arg(zaehler.y)
      .arg(nenner.y)
      .arg(symstrY)
      .arg(zaehler.z)
      .arg(nenner.z)
      .arg(symstrZ));

  erg.remove(QRegExp("0/\\d"));
  erg.replace("1/1","1");
  erg.replace("2/1","2");
  erg.replace("3/1","3");
  erg.replace("4/1","4");
  erg.replace("5/1","5");
  erg.replace("6/1","6");
  return erg;
}

QString Molecule::symmcode2human(QStringList brauchSymm){
  /*! @param brauchSymm list of internal symmetry codes 
   *  \returns human readable string sepearated by <br>
   */
  QString erg;
  Matrix m;
  V3 t;
  int h,k,l,s;
  for (int j=0; j<brauchSymm.size();j++){
    sscanf(brauchSymm[j].toLatin1(),"%d_%1d%1d%1d",&s,&h,&k,&l);
    h-=5;
    k-=5;
    l-=5;
    s--;
    m=cell.symmops.at(s);
    QString symstrX,symstrY,symstrZ;
    if (m.m11==1) symstrX.append("+x");
    if (m.m11==-1) symstrX.append("-x");
    if (m.m21==1) symstrX.append("+y");
    if (m.m21==-1) symstrX.append("-y");
    if (m.m31==1) symstrX.append("+z");
    if (m.m31==-1) symstrX.append("-z");

    if (m.m12==1) symstrY.append("+x");
    if (m.m12==-1) symstrY.append("-x");
    if (m.m22==1) symstrY.append("+y");
    if (m.m22==-1) symstrY.append("-y");
    if (m.m32==1) symstrY.append("+z");
    if (m.m32==-1) symstrY.append("-z");

    if (m.m13==1) symstrZ.append("+x");
    if (m.m13==-1) symstrZ.append("-x");
    if (m.m23==1) symstrZ.append("+y");
    if (m.m23==-1) symstrZ.append("-y");
    if (m.m33==1) symstrZ.append("+z");
    if (m.m33==-1) symstrZ.append("-z");
    t=cell.trans.at(s);
    V3 zaehler,nenner;
    double egal;
    for (int g=1;g<7; g++){
      nenner.x=(t.x<0)?-g:g;
      zaehler.x=static_cast<int>(round((t.x+h)*g));
      if (fabs(modf(t.x*g,&egal))<0.05) break;
    }
    for (int g=1;g<7; g++){
      nenner.y=(t.y<0)?-g:g;
      zaehler.y=static_cast<int>(round((t.y+k)*g));
      if (fabs(modf(t.y*g,&egal))<0.05) break;
    }
    for (int g=1;g<7; g++){
      nenner.z=(t.z<0)?-g:g;
      zaehler.z=static_cast<int>(round((t.z+l)*g));
      if (fabs(modf(t.z*g,&egal))<0.05) break;
    }
    erg.append(QString("&raquo;%1 :<b> %2/%3%4, %5/%6%7, %8/%9%10</b><br>")
        .arg(j+1)
        .arg(zaehler.x)
        .arg(nenner.x)
        .arg(symstrX)
        .arg(zaehler.y)
        .arg(nenner.y)
        .arg(symstrY)
        .arg(zaehler.z)
        .arg(nenner.z)
        .arg(symstrZ));

  }
  erg.remove(QRegExp("0/\\d"));
  erg.replace("1/1","1");
  erg.replace("2/1","2");
  erg.replace("3/1","3");
  erg.replace("4/1","4");
  erg.replace("5/1","5");
  erg.replace("6/1","6");
  return erg;
}

void Molecule::frac2kart(V3 x, V3 &y){
  /*!
   * converts a fractional coordinate into a cartesian using cell information.
   * @param[in] x fractional coordinate.
   * @param[out] y cartesian coordinate.
   */
  Matrix u;
  u.m11 = cell.a;
  u.m21 = 0.0;
  u.m31 = 0.0;
  u.m12 = cell.b * cell.cosga;
  u.m22 = cell.b * cell.singa;
  u.m32 = 0.0;
  u.m13 = cell.c * cell.cosbe;
  u.m23 = cell.tau;
  u.m33 = cell.c * cell.phi / cell.singa;

  y.x = x.x * u.m11 + x.y * u.m12 + x.z * u.m13;
  y.y = x.x * u.m21 + x.y * u.m22 + x.z * u.m23;
  y.z = x.x * u.m31 + x.y * u.m32 + x.z * u.m33;

}

void Molecule::kart2frac(V3 x, V3 &y){
  /*!
   * converts a cartesian coordinate into a fractional using cell information.
   * @param[in] x cartesian coordinate.
   * @param[out] y fractional coordinate.
   */
  Matrix u;
  y.x =0.0 ;
  y.y =0.0 ;
  y.z =0.0;
  u.m11 = 1.0/cell.a;
  u.m21 = 0.0;
  u.m31 = 0.0;
  u.m12 = -1.0/(cell.a * cell.tanga);
  u.m22 = 1.0/(cell.b * cell.singa);
  u.m32 = 0.0;
  u.m13 = (cell.cosal*cell.cosga-cell.cosbe)/(cell.a*cell.phi*cell.singa);
  u.m23 = (cell.cosbe*cell.cosga-cell.cosal)/(cell.b*cell.phi*cell.singa);
  u.m33 = cell.singa/(cell.c*cell.phi);
  // Wird jetzt hier genauso wie in int Tab B S.345 gerechnet (M^-1^).
  y.x = x.x * u.m11 + x.y * u.m12 + x.z * u.m13;
  y.y = x.x * u.m21 + x.y * u.m22 + x.z * u.m23;
  y.z = x.x * u.m31 + x.y * u.m32 + x.z * u.m33;
}

QString Molecule::pse(int oz){


  if ((oz>-1)&&(oz<PSE.size())) return PSE.at(oz);
  switch (oz){
  case -1: return "Q";
  case -42: return "BELO";
  case -66: return "Cnt";
  default: return "";
  }
}

int Molecule::getOZ(QString S1){
  QString s=S1;
  s=s.remove('$');
  s=s.section(QRegExp("[^A-Za-z]"),0,0);
  s=s.toUpper();
  if (s=="CNT") return -66;
  if (s=="D") return 0;
  QRegExp r = QRegExp(s,Qt::CaseInsensitive);
  return PSE.indexOf(r);
}

void Molecule::uniqueInCell(){
  /*! Moves centers of gravity of each fragment into the unit cell box and close to the origin.
   * */
  QList<int> flags;
  int toflipp=0;
  for (int i=0; i<asymm.size(); i++) if (asymm.at(i).an>-1){
    flags.append(-1);
    toflipp++;
  }else flags.append(0);
  flags[0]=1;
  //  printf("MaxMol = %d\n",maxmol);
  double dk;
  V3 prime,dp,D,floorD;
  //  V3 D,floorD;
  QList<V3> minMol;
  for (int i=0; i<maxmol; i++) minMol.append(V3(0,0,0));
  SdmItem sdmItem;
  sdmItem.a1=0;
  sdmItem.a2=1;
  sdmItem.sn=0;
  sdmItem.covalent=true;
  QList<SdmItem> starter;
  QList<int> mii;
  for (int j=0; j<maxmol; j++){
    int mi=0;
    for (int i=0; i<asymm.size(); i++){
      if (asymm.at(i).an<0) continue;
      if ((asymm.at(i).molindex-1)!=j) continue;
      minMol[j]+=asymm.at(i).frac;
      mi++;
    }
    mii.append(mi);
    minMol[j]*=1.0/mi;
    // printf("mole#%-2d has %3d atoms. Center of mass is at %9.4f %9.4f %9.4f \n",j,mi, minMol[j].x, minMol[j].y, minMol[j].z);
  }
  for (int j=0; j<maxmol; j++){
    if (mii.at(j)==0){
      sdmItem.sn=0;
      sdmItem.floorD=V3(floor(minMol.at(j).x),floor(minMol.at(j).y),floor(minMol.at(j).z));
      sdmItem.a2=sdmItem.a1=j;
      starter.append(sdmItem);
    }else{
      double min=99999.0;
      for (int n=0;n<cell.symmops.size();  n++){
        prime=cell.symmops.at(n) * minMol.at(j) + cell.trans.at(n);
        D=prime ;//+ V3(0.5,0.5,0.5) ;
        floorD=V3(floor(D.x),floor(D.y),floor(D.z));
        dp=D - floorD - V3(0.5,0.5,0.5);
        dk=fl(dp.x,dp.y,dp.z);
        if (dk<min){
          sdmItem.d=dk;
          min=dk;
          sdmItem.floorD=floorD;
          sdmItem.a1=j;
          sdmItem.a2=j;
          sdmItem.sn=n;
        }
      }
      //      printf("Mole#%-2d floorD %3.0f %3.0f %3.0f %d\n",j,sdmItem.floorD.x,sdmItem.floorD.y,sdmItem.floorD.z,sdmItem.sn);
      starter.append(sdmItem);
    }
  }
  for (int i=0; i<asymm.size(); i++){
    if (asymm.at(i).an<0) continue;
    if (asymm.at(i).molindex<1) continue;
    int n = starter.at(asymm.at(i).molindex-1).sn;
    V3 f = starter.at(asymm.at(i).molindex-1).floorD;
    asymm[i].frac= cell.symmops.at(n) * asymm.at(i).frac + cell.trans.at(n) - f;
    frac2kart(asymm[i].frac,asymm[i].pos);

  }
  for (int j=0; j<maxmol; j++){
    int mi=0;
    minMol[j]=V3(0,0,0);
    for (int i=0; i<asymm.size(); i++){
      if ((asymm.at(i).molindex-1)!=j) continue;
      minMol[j]+=asymm.at(i).frac;
      mi++;
    }
    minMol[j]*=1.0/mi;
    // printf("mole#%-2d has %2d atoms. Center of mass is at %9.4f %9.4f %9.4f \n",j,mi, minMol[j].x, minMol[j].y, minMol[j].z);
  }
  packer(sdmcompleter());
  showatoms.clear();
  for (int i=0; i<asymm.size();i++){
    showatoms.append(asymm[i]);
    showatoms[i].molindex=asymm[i].molindex;
  }
  showbonds.clear();
  showbonds=connecting(showatoms);
}

void Molecule::enviSDM(double range){
  /*! Calculates the shortes distance matrix for the ENVIronment functionality in the given range.
   * @param range maximal distances around each atom to generate a matrix entry. 
   */

  // George Sheldrick Seminar ideas
  double dk,dddd;
  V3 prime,dp,D,floorD;
  envi_sdm.clear();
  SdmItem sdmItem;
  sdmItem.a1=0;
  sdmItem.a2=1;
  sdmItem.sn=0;
  //sdmItem.p1=sdmItem.p2=V3(0,0,0);
  for (int i=0; i<asymm.size(); i++){ 
    for (int j=0; j<asymm.size(); j++ ){
      //  bool hma=false;
      for (int n=0;n<cell.symmops.size();  n++){
        prime=cell.symmops.at(n) * asymm.at(i).frac + cell.trans.at(n);
        D=prime - asymm.at(j).frac+ V3(0.5,0.5,0.5) ;
        floorD=V3(floor(D.x),floor(D.y),floor(D.z));
        for (int h=-1; h<2; h++){
          for (int k=-1; k<2; k++){
            for (int l=-1; l<2; l++){
              V3 fD=floorD+V3(h,k,l);  
              dp=D - fD - V3(0.5,0.5,0.5);
              dk=fl(dp.x,dp.y,dp.z);
              if ((dk>0.01)&&((range)>=dk)){
                sdmItem.d=dk;
                sdmItem.floorD=fD;
                sdmItem.a1=i;
                sdmItem.a2=j;
                sdmItem.sn=n;
                if ((asymm[i].an>-1)&&(asymm[j].an>-1)&&
                    ((asymm[i].part*asymm[j].part==0)||
                     (asymm[i].part==asymm[j].part)))
                  dddd=(Kovalenz_Radien[asymm[i].an]+ Kovalenz_Radien[asymm[j].an])*0.012;
                else dddd=0;
                if (sdmItem.d<dddd){
                  sdmItem.covalent=true;
                }else sdmItem.covalent=false;
                envi_sdm.append(sdmItem);
              }
            }
          }
        }
      }
    }
  }
#if (QT_VERSION >= 0x050000)
  std::sort(envi_sdm.begin(),envi_sdm.end());
#else
  qSort(envi_sdm.begin(),envi_sdm.end());
#endif
  //histogram
  /*
     double hmi,hma;
     hmi=envi_sdm.first().d;
     hma=envi_sdm.last().d;
     double hstep=(hma-hmi)/20.0;
     printf("%g %g %g %d\n",hma,hmi,hstep,envi_sdm.size());
     if (hstep<0.01) return;
     int hbox[20]={
     0,0,0,0,0,
     0,0,0,0,0,
     0,0,0,0,0,
     0,0,0,0,0};
     for (int i=0; i<envi_sdm.size();i++){
     int hind=(int) (envi_sdm.at(i).d/hstep);
     hbox[hind]++;
     }
     for (int i=0;i<20;i++){
     printf("%g-%g %d  (%g)\n",hmi+hstep*i,hmi+i*hstep+hstep,hbox[i],hma);
     }
  // */
}

void Molecule::neighborSort(Knopf &kn){//sort neighbors by atomic number
  // I think bublle sort is fine for such a sort list
  for (int i = 0; i < kn.neighbors.size(); i++){ 
    for (int j = i+1; j < kn.neighbors.size(); j++){
      if ((kn.neighbors.at(i)>=asymm.size())||(kn.neighbors.at(j)>=asymm.size())) continue;
      if (asymm.at(kn.neighbors.at(i)).an < asymm.at(kn.neighbors.at(j)).an){
        int swap=kn.neighbors.at(i);
        kn.neighbors[i]=kn.neighbors.at(j);
        kn.neighbors[j]=swap; 
      }
    }
  }
}

QList<SdmItem> Molecule::computeSDM(CEnvironment au){
  double dk,dddd;
  V3 prime,dp,D,floorD;
  SdmItem sdmItem;
  sdmItem.a1=0;
  sdmItem.a2=0;
  sdmItem.sn=0;
  QList<SdmItem> _sdm;
  for (int i=0; i<au.size(); i++){
    for (int j=i+1; j<au.size(); j++ ){
      double min=1000000;
      bool hma=false;
      for (int n=0;n<cell.symmops.size();  n++){
        //for (int n=0;n<1;  n++){
        prime=cell.symmops.at(n) * au.at(i).frac + cell.trans.at(n);
        D=prime - au.at(j).frac+ V3(0.5,0.5,0.5) ;
        floorD=V3(floor(D.x),floor(D.y),floor(D.z));
        dp=D - floorD - V3(0.5,0.5,0.5);
        dk=fl(dp.x,dp.y,dp.z);
        if (n) dk+=0.0001;
        if ((dk>0.01)&&((min+0.00)>=dk)){
          min=fmin(dk,min);
          sdmItem.d=min;
          sdmItem.floorD=floorD;
          sdmItem.a1=i;
          sdmItem.a2=j;
          sdmItem.sn=n;
          hma=true;
        }
      }
      if ((au[sdmItem.a1].an>-1)&&(au[sdmItem.a2].an>-1)&& ((au[sdmItem.a1].part*au[sdmItem.a2].part==0)||
           (au[sdmItem.a1].part==au[sdmItem.a2].part))){
        dddd=(Kovalenz_Radien[au[sdmItem.a1].an]+ Kovalenz_Radien[au[sdmItem.a2].an])*0.012;
      }
      else {
        dddd=0;
      }
      sdmItem.covalent=(sdmItem.d<dddd);
      if (hma) {
        _sdm.append(sdmItem);
      }
    }

  }
  return _sdm;
}

QStringList Molecule::sdmcompleter(){
  /*! Calculates the shortest distance matrix. 
   * Counts the number of fragmets in the asymmetric unit. 
   * Finds possible hydrogen bond contacts.
   * \returns internal symmtry code list for grown structures.
   */
  // George Sheldrick Seminar ideas
  //  printf("sdm1 %d\n",__LINE__);
  //  QTime zeit; zeit.start();
  //printf("sdmcompleter\n");
  double dk,dddd;
  V3 prime,dp,D,floorD;
  contact.clear();
  sdm.clear();
  SdmItem sdmItem;
  sdmItem.a1=0;
  sdmItem.a2=0;
  sdmItem.sn=0;
  SdmItem sdmItem2;
  sdmItem2.a1=0;
  sdmItem2.a2=0;
  sdmItem2.sn=0;

  knoepfe.clear();
  Knopf kn;
  kn.neighbors.clear();
  for (int i=0; i<asymm.size(); i++){
    kn.neighbors.clear();
    for (int j=0; j<asymm.size(); j++ ){
      double min=1000000;
      bool hma=false;
      for (int n=0;n<cell.symmops.size();  n++){
        //for (int n=0;n<1;  n++){
        //if ((n>0)&&((asymm[i].part<0)||(asymm[j].part<0))) continue;
        prime=cell.symmops.at(n) * asymm.at(i).frac + cell.trans.at(n);
        D=prime - asymm.at(j).frac+ V3(0.5,0.5,0.5) ;
        floorD=V3(floor(D.x),floor(D.y),floor(D.z));
        dp=D - floorD - V3(0.5,0.5,0.5);
        dk=fl(dp.x,dp.y,dp.z);
        if (n) dk+=0.0001;
        if ((dk>0.01)&&((min+0.00)>=dk)){
          min=fmin(dk,min);
          sdmItem.d=min;
          sdmItem.floorD=floorD;
          sdmItem.a1=i;
          sdmItem.a2=j;
          sdmItem.sn=n;
          hma=true;
        }
      }
      if ((asymm[sdmItem.a1].an>-1)&&(asymm[sdmItem.a2].an>-1)&&
          ((asymm[sdmItem.a1].part*asymm[sdmItem.a2].part==0)||
           (asymm[sdmItem.a1].part==asymm[sdmItem.a2].part)))
        dddd=(Kovalenz_Radien[asymm[sdmItem.a1].an]+ Kovalenz_Radien[asymm[sdmItem.a2].an])*0.0115;
      else dddd=0;
      if (sdmItem.d<dddd){
        if (hma) kn.neighbors.append(j);
        sdmItem.covalent=true;
      } else {
        sdmItem.covalent=false;
      }
      //if ((sdmItem.sn!=0)&&((asymm[sdmItem.a1].part<0)||(asymm[sdmItem.a2].part<0))) sdmItem.covalent=false;
      if (hma) sdm.append(sdmItem);
      }
      neighborSort(kn);
      knoepfe.append(kn);
    }

    for (int i=0; i<asymm.size(); i++){ //wasserstoffbrueckensuche
      for (int j=0; j<asymm.size(); j++ ){
        for (int nn=0;nn<cell.symmops.size();  nn++){
          prime=cell.symmops.at(nn) * asymm.at(i).frac + cell.trans.at(nn);
          D=prime - asymm.at(j).frac+ V3(0.5,0.5,0.5) ;
          floorD=V3(floor(D.x),floor(D.y),floor(D.z));
          dp=D - floorD - V3(0.5,0.5,0.5);
          dk=fl(dp.x,dp.y,dp.z);
          if((1)&&(abs(asymm[i].an-7)<2)&&(abs(asymm[j].an-7)<2)&&(fabs(dk-2.725)<0.275)){//N O F mit abstanden zwischen 2.45 bis 3 A koennten H=Bruecken sein
            sdmItem2.d=dk;
            sdmItem2.floorD=floorD;
            sdmItem2.a1=i;
            sdmItem2.a2=j;
            sdmItem2.sn=nn;
            sdmItem2.covalent=true;
            double ang1=361,ang2=361,an;
            for (int oo=0; oo<knoepfe.at(i).neighbors.size(); oo++){
              an=wl(asymm.at(j).frac,
                  cell.symmops.at(nn) * asymm.at(i).frac + cell.trans.at(nn)-floorD,
                  cell.symmops.at(nn) * asymm.at(knoepfe.at(i).neighbors.at(oo)).frac + cell.trans.at(nn)-floorD);
              ang1=(an<ang1)?an:ang1;
            }
            for (int oo=0; oo<knoepfe.at(j).neighbors.size(); oo++){
              an=wl(
                  cell.symmops.at(nn) * asymm.at(i).frac  + cell.trans.at(nn)-floorD,
                  asymm.at(j).frac,
                  asymm.at(knoepfe.at(j).neighbors.at(oo)).frac);
              ang2=(an<ang2)?an:ang2;
            }
            /*printf("Paarschippen: %s %s %f %f %f\n",
              asymm.at(i).Label.toStdString().c_str(),
              asymm.at(j).Label.toStdString().c_str(),
              dk, ang1, ang2);*/
            if ((qMax(ang1,ang2)!=361)&&((qMin(ang1,ang2)<85)||(qMax(ang1,ang2)>180))) sdmItem2.covalent=false;
            else {
              contact.append(sdmItem2);
            }
          }
        }
      }
    }
    //  printf("sdm3 %d\n",zeit.elapsed());//\n"); zeit.restart();
    QString bs;
    QStringList brauchSymm;
    int someleft=0,nextmol=0;
    maxmol=1;
    for (int i=0; i<asymm.size();i++) asymm[i].molindex=-1;
    asymm[0].molindex=1;//starter;
    do {
      nextmol=0;
      do {
        someleft=0;
        for  (int i=0; i<sdm.size();i++){
          if ((sdm.at(i).covalent)&&(asymm[sdm.at(i).a1].molindex*asymm[sdm.at(i).a2].molindex<0)) {
            asymm[sdm.at(i).a1].molindex=maxmol;
            asymm[sdm.at(i).a2].molindex=maxmol;
            someleft++;
          }
        }
      }while (someleft);
      for ( int i=0; i<asymm.size();i++) {
        if ((asymm[i].an>-1)&&(asymm[i].molindex<0)) {nextmol=i;break;}
      }
      if (nextmol) {
        asymm[nextmol].molindex=(++maxmol);
      }
    }while (nextmol);
    Fragments=QString("<b>The asymmetric unit contains %1 fragments.</b><br>").arg(maxmol);
    for (int k =0; k<sdm.size();k++){
      if (sdm.at(k).covalent){
        if (asymm[sdm.at(k).a1].molindex<1) continue;
        if (asymm[sdm.at(k).a1].molindex>6) continue;
        for (int n=0;n<cell.symmops.size();  n++){
          if (((asymm[sdm.at(k).a1].part!=0)&&(asymm[sdm.at(k).a2].part!=0)&&(asymm[sdm.at(k).a1].part!=asymm[sdm.at(k).a2].part)))continue;
          if ((asymm[sdm.at(k).a1].an==asymm[sdm.at(k).a2].an)&&(asymm[sdm.at(k).a1].an==0)) continue;
          prime=cell.symmops.at(n) * asymm[sdm.at(k).a1].frac + cell.trans.at(n);
          D=prime - asymm[sdm.at(k).a2].frac+ V3(0.5,0.5,0.5) ;
          floorD=V3(floor(D.x),floor(D.y),floor(D.z));
          dp=D - floorD - V3(0.5,0.5,0.5);
          if ((n==0)&&(V3(0,0,0)==floorD)) continue;
          dk=fl(dp.x,dp.y,dp.z);
          dddd=(sdm.at(k).d+0.2);
          if ((dk>0.001)&&(dddd>=dk)) {
            bs=QString("%1_%2%3%4:%5,").arg(n+1).arg(5-static_cast<int>(floorD.x)).arg(5-static_cast<int>(floorD.y)).arg(5-static_cast<int>(floorD.z)).arg(asymm[sdm.at(k).a1].molindex);
            if  ((!brauchSymm.contains(bs))) {
              brauchSymm.append(bs);
              //            printf("%d %d %d %d [%d]\n",n+1,(int)floorD.x,(int)floorD.y,(int)floorD.z ,brauchSymm.size());
            }
          }
        }
      }
    }
    //  printf("sdm4 %d\n",zeit.elapsed()); zeit.restart();
#if (QT_VERSION >= 0x050000)
    std::sort(sdm.begin(),sdm.end());
#else
    qSort(sdm.begin(),sdm.end());
#endif
    //  printf("sdm5 %d\n",zeit.elapsed());//\n"); zeit.restart();
    enviSDM(2.5);
    //  printf("sdm5 %d\n",zeit.elapsed());//\n"); zeit.restart();

    QList<int> flags;
    for (int i=0; i<asymm.size(); i++) flags.append((asymm.at(i).an<0)?-1:1);
    for (int i=0; i<envi_sdm.size(); i++){
      if ((flags.at(envi_sdm.at(i).a1)*flags.at(envi_sdm.at(i).a2))==-1) {
        //if ((qbeforehkl)||((flags.at(sdm.at(i).a1)==-1)&&((asymm.at(sdm.at(i).a2).an>-1)))){
        //   if (((asymm[envi_sdm.at(i).a1].Label=="Q11")||(asymm[envi_sdm.at(i).a2].Label=="Q11"))&&(envi_sdm.at(i).d<2.6))
        //       qDebug()<<asymm[envi_sdm.at(i).a1].Label<<asymm[envi_sdm.at(i).a2].Label<<envi_sdm.at(i).d<<flags[envi_sdm.at(i).a1]<<flags[envi_sdm.at(i).a2];
        if (asymm[envi_sdm.at(i).a1].an>-1) continue;
        if (asymm[envi_sdm.at(i).a1].an==-66) continue;
        if ((envi_sdm.at(i).sn==0)&&(envi_sdm.at(i).floorD==V3(0,0,0))) {flags[envi_sdm.at(i).a1]=1;continue;}
        if (envi_sdm.at(i).d>2.4) continue;
        asymm[envi_sdm.at(i).a1].frac = cell.symmops.at(envi_sdm.at(i).sn) * asymm[envi_sdm.at(i).a1].frac + cell.trans.at(envi_sdm.at(i).sn) - envi_sdm.at(i).floorD;
        flags[envi_sdm.at(i).a1]=1;
        frac2kart(asymm[envi_sdm.at(i).a1].frac,asymm[envi_sdm.at(i).a1].pos);
        //}
      }
    }
    //  printf("sdm6 %d\n",zeit.elapsed()); zeit.restart();
    return brauchSymm;
  }

  /*
     qDebug()<<neededSymm.size()<<ss;
  // */

  Connection Molecule::connecting(const CEnvironment &atom,bool eac){
    /*! Creates bonds and lbonds depending on covalent radii.
     * If edit atom style flag (eac) is not set FREE and BIND instructions are considered.
     * Feeds also Molecule.knoepfe list.
     * \returns the Connection bond list.
     * */
    Connection bonds;
    bonds.clear();
    lbonds.clear();
    MyBond bond;
    double soll=0;
    for (int i=0; i<atom.size(); i++){
      for (int j=0; j<=i; j++ ){
        double d = sqrt(Distance(atom.at(i).pos,atom.at(j).pos));
        if ((atom.at(i).an>-1)&&(atom.at(j).an>-1)&&(atom.at(i).part!=0)&&(atom.at(j).part!=0)&&((bindPart.key(atom.at(i).part)==atom.at(j).part)||(bindPart.value(atom.at(i).part)==atom.at(j).part)))
          soll=(Kovalenz_Radien[atom.at(i).an]+ Kovalenz_Radien[atom.at(j).an])*0.012;
        else {
          if ((atom.at(i).part<0)&&(atom.at(j).part!=0)&&(atom.at(i).symmGroup!=atom.at(j).symmGroup)) continue;
          if ((atom.at(j).part<0)&&(atom.at(i).part!=0)&&(atom.at(i).symmGroup!=atom.at(j).symmGroup)) continue;
          if ((atom.at(i).an==0)&&(atom.at(j).an==0)) continue;
          if ((atom.at(i).an>-1)&&(atom.at(j).an>-1)&&
              ((!atom.at(i).part)||(!atom.at(j).part)||(atom.at(i).part==atom.at(j).part)))
            soll=(Kovalenz_Radien[atom.at(i).an]+ Kovalenz_Radien[atom.at(j).an])*0.012;
          else
            soll=-1;
        }
        if ((d>0.1)&&(d<soll)) {
          bond.ato1=&atom[i];
          bond.ato2=&atom[j];
          bond.a1=i;
          bond.a2=j;
          bond.length=d;
          bonds.append(bond);
        }
        else if ((d<qboMax)&&(d>qboMin)){
          if ((atom.at(i).an==-42)||(atom.at(j).an==-42)) continue;//belo: wuff!
          if ((atom.at(i).an==0)||(atom.at(j).an==0)) continue;
          if  ((!atom.at(i).part)||(!atom.at(j).part)||(atom.at(i).part==atom.at(j).part)){ 
            bond.ato1=&atom[i];
            bond.ato2=&atom[j];
            bond.a1=i;
            bond.a2=j;
            bond.length=d;
            lbonds.append(bond);
          }
        }
      }
    }

    if (!eac) for (int i=0; i<freeatoms.size(); i++){
      MyAtom a1, a2;
      bool bo=false;
      if (freeatoms.at(i).Lab1.contains('_')) a1.resiNr=freeatoms.at(i).Lab1.section('_',1,1).toInt(&bo);
      if (!bo) a1.resiNr=-1;
      if (freeatoms.at(i).Lab2.contains('_'))a2.resiNr=freeatoms.at(i).Lab2.section('_',1,1).toInt(&bo);
      if (!bo)a2.resiNr=-1;
      a1.Label = freeatoms.at(i).Lab1.section('_',0,0);
      a2.Label = freeatoms.at(i).Lab2.section('_',0,0);
      //printf("#%d %d\n",a2.resiNr,a1.resiNr);

      for (int j=0;j<bonds.size();j++){
        if (((a1.resiNr>-1)&& (bonds.at(j).ato1->resiNr!=a1.resiNr))|| 
            ((a2.resiNr>-1)&& (bonds.at(j).ato2->resiNr!=a2.resiNr))) continue;
        if (((a1.resiNr>-1)&& (bonds.at(j).ato2->resiNr!=a1.resiNr))|| 
            ((a2.resiNr>-1)&& (bonds.at(j).ato1->resiNr!=a2.resiNr))) continue;
        if (((bonds.at(j).ato1->Label.section('_',0,0).toUpper()== a1.Label)&&
              (bonds.at(j).ato2->Label.section('_',0,0).toUpper()== a2.Label))||
            ((bonds.at(j).ato2->Label.section('_',0,0).toUpper()== a1.Label)&&
             (bonds.at(j).ato1->Label.section('_',0,0).toUpper()== a2.Label))){
          bonds.removeAt(j);
          j=qMax(j-1,0);
        }
      }
    }
    if (!eac) for (int i=0; i<bindatoms.size(); i++){
      MyAtom a1, a2;
      bool bo=false;
      int sy1=-1,sy2=-1;
      if (bindatoms.at(i).Lab1.contains("_$"))sy1=bindatoms.at(i).Lab1.section("_$",1,1).toInt(&bo);
      if (bindatoms.at(i).Lab2.contains("_$"))sy2=bindatoms.at(i).Lab2.section("_$",1,1).toInt(&bo);
      bo=false;
      if (bindatoms.at(i).Lab1.contains('_'))a1.resiNr=bindatoms.at(i).Lab1.section('_',1,1).toInt(&bo);
      if (!bo) a1.resiNr=-1;
      if (bindatoms.at(i).Lab2.contains('_'))a2.resiNr=bindatoms.at(i).Lab2.section('_',1,1).toInt(&bo);
      if (!bo)a2.resiNr=-1;
      a1.Label = bindatoms.at(i).Lab1.section('_',0,0);
      a2.Label = bindatoms.at(i).Lab2.section('_',0,0);
      //printf("#>%d %d %d %d [%d]\n",a2.resiNr,a1.resiNr,sy1,sy2,i);
      double mindd=50.0;
      //printf("-Q %d\n",eqivMap.size());
      if ((!eqivMap.isEmpty())&&(qMax(sy1,sy2)>0)) {
        //    printf("sg %d\n",eqivMap.value(QString("$%1").arg(sy2)));
        //    printf("%s %s\n",a1.Label.toUpper().toStdString().c_str(),a2.Label.toUpper().toStdString().c_str());
        for (int iii=0; iii<atom.size(); iii++){
          for (int jjj=0; jjj<atom.size(); jjj++){
            if ((a1.resiNr>-1)&& (atom.at(iii).resiNr!=a1.resiNr)) continue;
            if ((a2.resiNr>-1)&& (atom.at(jjj).resiNr!=a2.resiNr)) continue;
            if ((sy1>-0)&&(atom.at(iii).symmGroup==0)) continue;
            if ((sy2>-0)&&(atom.at(jjj).symmGroup==0)) continue;
            //printf("sg %d\n",eqivMap.value(QString("$%1").arg(sy1)));
            if ((sy1>-0)&&(eqivMap.value(QString("$%1").arg(sy1))!=atom.at(iii).symmGroup)) continue;
            if ((sy2>-0)&&(eqivMap.value(QString("$%1").arg(sy2))!=atom.at(jjj).symmGroup)) continue;
            if ((atom.at(iii).Label.section('_',0,0).section(QString::fromUtf8("»"),0,0).toUpper() == a1.Label.toUpper())&&
                (atom.at(jjj).Label.section('_',0,0).section(QString::fromUtf8("»"),0,0).toUpper() == a2.Label.toUpper())){
              bond.ato1=&atom[iii];
              bond.ato2=&atom[jjj];
              bond.a1=iii;
              bond.a2=jjj;
              bond.length=sqrt(Distance(atom[iii].pos,atom[jjj].pos));
              //printf("%g %f\n",bond.length,mindd);
              if ((mindd+0.01)>=bond.length) bonds.append(bond);
              if (bond.length>0.7) mindd=qMin(bond.length, mindd);

            }
          }
        }
      }
      for (int j=0; j < bonds.size(); j++){
        if (((bonds.at(j).ato1->Label.section('_',0,0).section(QString::fromUtf8("»"),0,0).toUpper()==a1.Label.toUpper())&&
              (bonds.at(j).ato2->Label.section('_',0,0).section(QString::fromUtf8("»"),0,0).toUpper()==a2.Label.toUpper()))||
            ((bonds.at(j).ato1->Label.section('_',0,0).section(QString::fromUtf8("»"),0,0).toUpper()==a2.Label.toUpper())&&
             (bonds.at(j).ato2->Label.section('_',0,0).section(QString::fromUtf8("»"),0,0).toUpper()==a1.Label.toUpper()))) 
          if (bonds.at(j).length>(mindd+0.1)) {
            double d=bonds.at(j).length;
            printf("-- %f %f\n",d,mindd);
            bonds.removeAt(j);
            if (j) j--;
          }
      }
    }
    //-->Verknopfung   !!!
    if (!eac){
      for (int i=0; i<knoepfe.size();i++) knoepfe[i].neighbors.clear();
      knoepfe.clear();
      Knopf kn;
      kn.neighbors.clear();
      //for (int i=0; i<asymm.size(); i++){
      for (int i=0; i<atom.size(); i++){
        kn.neighbors.clear();
        for (int j=0; j<bonds.size(); j++){
          if (bonds.at(j).a1==i){
            kn.neighbors.append(bonds.at(j).a2);
          }
          if (bonds.at(j).a2==i){
            kn.neighbors.append(bonds.at(j).a1);
          }
        }
        neighborSort(kn);
        knoepfe.append(kn);
      }
    }
    return bonds;
    }

    double Molecule::fl(double x,double y,double z){
      double a,b,c;
      a = (cell.ga==90.0)?0.0:2.0*x*y*cell.a*cell.b*cell.cosga;
      b = (cell.be==90.0)?0.0:2.0*x*z*cell.a*cell.c*cell.cosbe;
      c = (cell.al==90.0)?0.0:2.0*y*z*cell.b*cell.c*cell.cosal;
      double erg=sqrt(x*x*cell.a*cell.a+
          y*y*cell.b*cell.b+
          z*z*cell.c*cell.c+
          a+b+c);
      return erg;
    }

    double Molecule::wl(V3 p2, V3 p1, V3 p3){
      V3 f1=p2-p1;//12
      V3 f2=p3-p1;//13
      V3 f3=p3-p2;//23
      double a=0,b=0,c=0;
      a =  (cell.ga==90.0)?0.0:2.0*f1.x*f1.y*cell.a*cell.b*cell.cosga;
      a += (cell.be==90.0)?0.0:2.0*f1.x*f1.z*cell.a*cell.c*cell.cosbe;
      a += (cell.al==90.0)?0.0:2.0*f1.y*f1.z*cell.b*cell.c*cell.cosal;

      b =  (cell.ga==90.0)?0.0:2.0*f2.x*f2.y*cell.a*cell.b*cell.cosga;
      b += (cell.be==90.0)?0.0:2.0*f2.x*f2.z*cell.a*cell.c*cell.cosbe;
      b += (cell.al==90.0)?0.0:2.0*f2.y*f2.z*cell.b*cell.c*cell.cosal;

      c = ( cell.ga==90.0)?0.0:2.0*f3.x*f3.y*cell.a*cell.b*cell.cosga;
      c += (cell.be==90.0)?0.0:2.0*f3.x*f3.z*cell.a*cell.c*cell.cosbe;
      c += (cell.al==90.0)?0.0:2.0*f3.y*f3.z*cell.b*cell.c*cell.cosal;
      double d12,d13,d23;
      d12 = f1.x*f1.x*cell.a*cell.a+
        f1.y*f1.y*cell.b*cell.b+
        f1.z*f1.z*cell.c*cell.c+a;	
      d13 = f2.x*f2.x*cell.a*cell.a+
        f2.y*f2.y*cell.b*cell.b+
        f2.z*f2.z*cell.c*cell.c+b;	
      d23 = f3.x*f3.x*cell.a*cell.a+
        f3.y*f3.y*cell.b*cell.b+
        f3.z*f3.z*cell.c*cell.c+c;
      double erg=acos((d12+d13-d23)/(2.0*sqrt(d12)*sqrt(d13)))*g2r;
      return erg;
    }

    void Molecule::cellSetup(){
      cell.cosal=cos(cell.al/g2r);
      cell.cosbe=cos(cell.be/g2r);
      cell.cosga=cos(cell.ga/g2r);
      cell.sinal=sin(cell.al/g2r);
      cell.sinbe=sin(cell.be/g2r);
      cell.singa=sin(cell.ga/g2r);
      cell.tanga=tan(cell.ga/g2r);
      cell.phi=  sqrt(1-(cell.cosal*cell.cosal)-
          (cell.cosbe*cell.cosbe)-(cell.cosga*cell.cosga)
          +2*cell.cosal*cell.cosbe*cell.cosga);
      cell.tau=cell.c *((cell.cosal- cell.cosbe* cell.cosga)/ cell.singa);
      cell.V = cell.a*cell.b*cell.c*cell.phi;
      cell.as=cell.c*cell.b*cell.sinal/cell.V;
      cell.bs=cell.c*cell.a*cell.sinbe/cell.V;
      cell.cs=cell.a*cell.b*cell.singa/cell.V;
      cell.cosra = (cell.cosbe*cell.cosga-cell.cosal)/(cell.sinbe*cell.singa);
      cell.cosrb = (cell.cosga*cell.cosal-cell.cosbe)/(cell.singa*cell.sinal);
      cell.cosrg = (cell.cosal*cell.cosbe-cell.cosga)/(cell.sinal*cell.sinbe);
      cell.G.m11=cell.a*cell.a;
      cell.G.m22=cell.b*cell.b;
      cell.G.m33=cell.c*cell.c;
      cell.G.m12=cell.G.m21=cell.a*cell.b*cell.cosga;
      cell.G.m13=cell.G.m31=cell.a*cell.c*cell.cosbe;
      cell.G.m23=cell.G.m32=cell.b*cell.c*cell.cosal;
      cell.Gi=inverse(cell.G);
    }

    void Molecule::applyLatticeCentro(int gitter){
      /*! Adds centrings and inversion symmetry to the symmops and trans list. 
       * @param gitter the parameter of the shelx LATT instruction. ==> see Shelxl manual.
       */
      int z=cell.symmops.size();
      cell.ns0=z;
      Matrix inv(-1.0,0.0,0.0, 0.0,-1.0,0.0, 0.0,0.0,-1.0);  
      cell.centeric=false;
      if (gitter>0){ 
        for (int i=0; i<z;i++){
          Matrix m=cell.symmops.at(i)*inv;
          cell.symmops.append(m);
          cell.trans.append(cell.trans.at(i));
        }
        cell.centeric=true;
      }
      latt=gitter;
      gitter=(gitter>0)?gitter:-gitter;
      z=cell.symmops.size();
      cell.centered=false;
      switch (gitter){
        case 5 :
          for (int i=0; i<z;i++){
            V3 tt = cell.trans.at(i)+V3(0.0, 0.5, 0.5);
            tt.x=(tt.x>1)?tt.x-1:tt.x;
            tt.y=(tt.y>1)?tt.y-1:tt.y;
            tt.z=(tt.z>1)?tt.z-1:tt.z;
            cell.symmops.append(cell.symmops.at(i));
            cell.trans.append(tt);
            cell.centered=true;
          }
          break;
        case 6 :
          for (int i=0; i<z;i++){
            V3 tt = cell.trans.at(i)+V3(0.5, 0.0, 0.5);
            tt.x=(tt.x>1)?tt.x-1:tt.x;
            tt.y=(tt.y>1)?tt.y-1:tt.y;
            tt.z=(tt.z>1)?tt.z-1:tt.z;
            cell.symmops.append(cell.symmops.at(i));
            cell.trans.append(tt);
            cell.centered=true;
          }
          break;
        case 7 :
          for (int i=0; i<z;i++){
            V3 tt = cell.trans.at(i)+V3(0.5, 0.5, 0.0);
            tt.x=(tt.x>1)?tt.x-1:tt.x;
            tt.y=(tt.y>1)?tt.y-1:tt.y;
            tt.z=(tt.z>1)?tt.z-1:tt.z;
            cell.symmops.append(cell.symmops.at(i));
            cell.trans.append(tt);
            cell.centered=true;
          }
          break;
        case 4 :
          for (int i=0; i<z;i++){
            V3 tt = cell.trans.at(i)+V3(0.0, 0.5, 0.5);
            tt.x=(tt.x>1)?tt.x-1:tt.x;
            tt.y=(tt.y>1)?tt.y-1:tt.y;
            tt.z=(tt.z>1)?tt.z-1:tt.z;
            cell.symmops.append(cell.symmops.at(i));
            cell.trans.append(tt);
            tt = cell.trans.at(i)+V3(0.5, 0.0, 0.5);
            tt.x=(tt.x>1)?tt.x-1:tt.x;
            tt.y=(tt.y>1)?tt.y-1:tt.y;
            tt.z=(tt.z>1)?tt.z-1:tt.z;
            cell.symmops.append(cell.symmops.at(i));
            cell.trans.append(tt);
            tt = cell.trans.at(i)+V3(0.5, 0.5, 0.0);
            tt.x=(tt.x>1)?tt.x-1:tt.x;
            tt.y=(tt.y>1)?tt.y-1:tt.y;
            tt.z=(tt.z>1)?tt.z-1:tt.z;
            cell.symmops.append(cell.symmops.at(i));
            cell.trans.append(tt);
            cell.centered=true;
          }
          break;
        case 2 :
          for (int i=0; i<z;i++){
            V3 tt = cell.trans.at(i)+V3(0.5, 0.5, 0.5);
            tt.x=(tt.x>1)?tt.x-1:tt.x;
            tt.y=(tt.y>1)?tt.y-1:tt.y;
            tt.z=(tt.z>1)?tt.z-1:tt.z;
            cell.symmops.append(cell.symmops.at(i));
            cell.trans.append(tt);
            cell.centered=true;
          }
          break;
        case 3 :
          for (int i=0; i<z;i++){
            V3 tt = cell.trans.at(i)+V3(2.0/3.0, 1.0/3.0, 1.0/3.0);
            tt.x=(tt.x>1)?tt.x-1:tt.x;
            tt.y=(tt.y>1)?tt.y-1:tt.y;
            tt.z=(tt.z>1)?tt.z-1:tt.z;
            cell.symmops.append(cell.symmops.at(i));
            cell.trans.append(tt);
            tt = cell.trans.at(i)+V3(1.0/3.0, 2.0/3.0, 2.0/3.0);
            tt.x=(tt.x>1)?tt.x-1:tt.x;
            tt.y=(tt.y>1)?tt.y-1:tt.y;
            tt.z=(tt.z>1)?tt.z-1:tt.z;
            cell.symmops.append(cell.symmops.at(i));
            cell.trans.append(tt);
            cell.centered=true;

          }
          break;
        case 0 :break;  
      }

    }

    double Molecule::dimension(){
      /*! calculates the maximal interatomic distance (Ang) of showatoms environment if this is 0 then 10.0 is returned instead.
      */
      double max=0,gg=0;
      for (int i=0;i<showatoms.size();i++)
        for (int j=i+1;j<showatoms.size();j++){
          if ((showatoms[j].an>=0)&&(showatoms[j].an>=0)) 
            gg=sqrt( Distance(showatoms[i].pos,showatoms[j].pos));
          max=(max<gg)?gg:max;
        }
      if (max==0) return 10.0;
      return max;
    }

    double Molecule::dimension(CEnvironment ce){
      /*! calculates the maximal interatomic distance (Ang) of showatoms environment if this is 0 then 10.0 is returned instead.
      */
      double max=0,gg=0;
      for (int i=0;i<ce.size();i++)
        for (int j=i+1;j<ce.size();j++){
          if ((ce[j].an>=0)&&(ce[j].an>=0)) 
            gg=sqrt( Distance(ce[i].pos,ce[j].pos));
          max=(max<gg)?gg:max;
        }
      if (max==0) return 10.0;
      return max;
    }

#define Ato4d(arr)       arr[0], arr[1], arr[2], arr[3]

    void Molecule::Uf2Uo(const Matrix x, Matrix & y) {
      /*! Turns fractional Uij's into cartesinan based Uij's.
       * @param[in] x fractional Uij matrix.
       * @param[out] y cartesinan Uij matrix.
       */
      Matrix o,a,w;   /*Cholesky decomposition of the real space Metric tensor
                        Wird fuer die Umrechnung von fraktionellen in kartesischen Korrdinaten benoetigt.*/
      a.m11 =cell.as;
      a.m12 = 0;
      a.m13 = 0;
      a.m21 = 0;
      a.m22 = cell.bs;
      a.m23 = 0;
      a.m31 = 0;
      a.m32 = 0;
      a.m33 = cell.cs;
      w=(a*x)*a;
      o.m11 =cell.a;
      o.m12 = 0.0;
      o.m13 = 0.0;
      o.m21 = cell.b*cell.cosga;
      o.m22 = cell.b*cell.singa;
      o.m23 = 0.0;
      o.m31 = cell.c* cell.cosbe;
      o.m32 = cell.tau;
      o.m33 = cell.c* cell.phi / cell.singa ;
      //  Matrix to=;
      y=(o*w)*transponse(o);
      /*  printf("%f %f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f %f\n%f %f %f %f %f %f %f %f %f\n",
          o.m11,o.m12,o.m13,o.m21,o.m22,o.m23,o.m31,o.m32,o.m33,
          to.m11,to.m21,to.m31,to.m12,to.m22,to.m32,to.m13,to.m23,to.m33,  
          x.m11,x.m12,x.m13,x.m21,x.m22,x.m23,x.m31,x.m32,x.m33,
          w.m11,w.m12,w.m13,w.m21,w.m22,w.m23,w.m31,w.m32,w.m33
          );*/
      //  printf("%f %f %f %f %f %f\n",y.m11*10000,y.m22*10000,y.m33*10000,y.m12*10000,y.m13*10000,y.m23*10000);//		  */
    }

    void Molecule::Usym (Matrix x,Matrix sym, Matrix & y){
      /*! Applies the symmetry matrix sym to Uij's 
       * @param[in] x Uij matrix.
       * @param[in] sym symmtry matrix.
       * @param[out] y resulting Uij matrix. 
       */
      //y=(transponse(sym)*x)*sym;
      y=(sym*x)*transponse(sym);
    }

    void Molecule::USymCode(Matrix x,int scod, Matrix & y){
      int s = scod / 1000;
      y=(cell.symmops.at(s)*x)*transponse(cell.symmops.at(s));


    }

//void Molecule::

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau); a[k][l]=h+s*(g-h*tau);
//NumericalRecepies....
double * Molecule::jacobi(const Matrix &uij, V3 &ev) { 
      int j,iq,ip,i,n=3,nrot; 
      double tresh=0,theta,tau,t,sm,s,h,g,c; 
      double a[3][3],b[3],z[3],v[3][3],d[3];
      a[0][0]=uij.m11;
      a[0][1]=uij.m12;
      a[0][2]=uij.m13;
      a[1][0]=uij.m21;
      a[1][1]=uij.m22;
      a[1][2]=uij.m23;
      a[2][0]=uij.m31;
      a[2][1]=uij.m32;
      a[2][2]=uij.m33;
      /*    printf("%12.6f %12.6f %12.6f\n%12.6f %12.6f %12.6f\n%12.6f %12.6f %12.6f\n"
            ,a[0][0],a[1][0],a[2][0]
            ,a[0][1],a[1][1],a[2][1]
            ,a[0][2],a[1][2],a[2][2]);*/
      static double erg[4]={0.0,1.0,0.0,0.0};
      for (ip=1;ip<=n;ip++) {   
        for (iq=1;iq<=n;iq++) v[ip-1][iq-1]=0.0; 
        v[ip-1][ip-1]=1.0; 
      } 
      for (ip=1;ip<=n;ip++) {  
        b[ip-1]=d[ip-1]=a[ip-1][ip-1]; 
        z[ip-1]=0.0; 
      } 
      nrot=0; 
      for (i=1;i<=150;i++) {
        sm=0.0; 
        for (ip=1;ip<=n-1;ip++) { 
          for (iq=ip+1;iq<=n;iq++) 
            sm += fabs(a[ip-1][iq-1]);       
        } 

        //printf("sm =%20.19f\n",sm);

        if (float(sm) < tresh) { 
          if ((v[0][0]+v[1][1]+v[2][2])!=3.0) {
            erg[0]=acos((v[0][0]+v[1][1]+v[2][2]-1.0)/2.0);
            erg[1]=(v[2][1]-v[1][2])/(2.0*sin(erg[0]));
            erg[2]=(v[0][2]-v[2][0])/(2.0*sin(erg[0]));
            erg[3]=(v[1][0]-v[0][1])/(2.0*sin(erg[0]));
            erg[0]*=180.0/M_PI;}
          else {erg[0]=0.0;erg[1]=1.0;erg[2]=0.0;erg[3]=0.0; }
          //printf("%d??ERG:%f %f %f %f\n",i,Ato4d(erg));
          /*
             printf("=%d======================================\n%8.5f %8.5f %8.5f \n%8.5f %8.5f %8.5f \n%8.5f %8.5f %8.5f \n%8.5f %8.5f %8.5f\n========================================\n",i,
             d[0],d[1],d[2],v[0][0],v[0][1],v[0][2]
             ,v[1][0],v[1][1],v[1][2]
             ,v[2][0],v[2][1],v[2][2]
             );*/
          ev=V3(d[0],d[1],d[2]);
          /*    printf("\n%12.6f %12.6f %12.6f\n%12.6f %12.6f %12.6f\n%12.6f %12.6f %12.6f\n"
                ,v[0][0],v[1][0],v[2][0]
                ,v[0][1],v[1][1],v[2][1]
                ,v[0][2],v[1][2],v[2][2]);
                if ((ev.x>0)&&(ev.y>0)&&(ev.z>0)){;}else{
                printf("NPD: %g %g %g \n",d[0],d[1],d[2]);
                }*/
          return static_cast<double*>(erg);
        }
        if (i < 4) tresh=0.00001;  
        else tresh=0.0001;
        for (ip=1;ip<=n-1;ip++) { 
          for (iq=ip+1;iq<=n;iq++) { 
            //printf("\np:%i q:%i i:%i nrot:%i\n",ip,iq,i,nrot);
            g=100.0*fabs(a[ip-1][iq-1]);  
            if ((i > 4) && ((fabs(d[ip-1])+g) == fabs(d[ip-1])) && ((fabs(d[iq-1])+g) == fabs(d[iq-1]))) {a[ip-1][iq-1]=0.0;}
            else if (fabs(a[ip-1][iq-1]) >= tresh) { 
              h=d[iq-1]-d[ip-1]; 
              if ((fabs(h)+g) == fabs(h)) {t=(a[ip-1][iq-1])/h; } 
              else { theta=0.5*h/(a[ip-1][iq-1]);  
                t=1.0/(fabs(theta)+sqrt(1.0+theta*theta)); 
                if (theta < 0.0) {t = -1.0*t;}
              } 
              c=1.0/sqrt(1+t*t); 
              s=t*c; 
              tau=s/(1.0+c); 
              h=t*a[ip-1][iq-1];
              z[ip-1] -= h;
              z[iq-1] += h;
              d[ip-1] -= h;
              d[iq-1] += h;
              a[ip-1][iq-1]=0.0;
              for (j=1;j<=ip-1;j++) { 
                ROTATE(a,j-1,ip-1,j-1,iq-1)
                  //printf("%i %i %i %i",j,ip,j,iq);
              } 
              for (j=ip+1;j<=iq-1;j++) { 
                ROTATE(a,ip-1,j-1,j-1,iq-1)
                  //printf("%i %i %i %i ",ip,j,j,iq);
              } 
              for (j=iq+1;j<=n;j++) {  
                ROTATE(a,ip-1,j-1,iq-1,j-1)
                  //printf("%i %i %i %i",ip,j,iq,j);
              } 
              for (j=1;j<=n;j++) { 
                ROTATE(v,j-1,ip-1,j-1,iq-1)
              } 
              ++(nrot);	  
              //    printf("U|\n%f %f %f  \n%f %f %f\n%f %f %f\nV|\n%f %f %f  \n%f %f %f\n%f %f %f\n\n",a[0][0],a[1][0],a[2][0],a[0][1],a[1][1],a[2][1],a[0][2],a[1][2],a[2][2],v[0][0],v[1][0],v[2][0],v[0][1],v[1][1],v[2][1],v[0][2],v[1][2],v[2][2]);
            } //else ;//printf("nix:%f p%i q%i",fabs(a[ip-1][iq-1]),ip,iq);
          } 
        } 
        for (ip=1;ip<=n;ip++) { 
          b[ip-1] += z[ip-1];
          d[ip-1] =b[ip-1];
          z[ip-1] =0.0;
        } 
      } 
      erg[0]=acos((v[0][0]+v[1][1]+v[2][2]-1.0)/2.0);
      if (erg[0]==0) {
        erg[1]=1.0;
        erg[2]=0.0;
        erg[3]=0.0;
      }else{
        erg[1]=(v[2][1]-v[1][2])/(2.0*sin(erg[0]));
        erg[2]=(v[0][2]-v[2][0])/(2.0*sin(erg[0]));
        erg[3]=(v[1][0]-v[0][1])/(2.0*sin(erg[0]));
        erg[0]*=180.0/M_PI;
      }
      /*printf("=%d=======================================\n%8.5f %8.5f %8.5f \n%8.5f %8.5f %8.5f \n%8.5f %8.5f %8.5f \n%8.5f %8.5f %8.5f\n========================================\n",i,
        d[0],d[1],d[2],v[0][0],v[0][1],v[0][2]
        ,v[1][0],v[1][1],v[1][2]
        ,v[2][0],v[2][1],v[2][2]

        );*/
      ev=V3(d[0],d[1],d[2]);
      /* if ((ev.x>0)&&(ev.y>0)&&(ev.z>0)){;}else{
         printf("NPD: %g %g %g \n",d[0],d[1],d[2]);
         }*/
      return static_cast<double*>(erg);
    }

void Molecule::fuse(){
  /*! reduces the showatoms list to the asymmetric unit.
   * */
  HumanSymmetry=QString("Used Symmetry:<br>%1").arg(symmcode2human(QStringList()));
  usedSymmetry.clear();
  //printf("FUSE\n");
  showatoms.clear();
  for (int i=0; i<asymm.size();i++){
    showatoms.append(asymm[i]);
    showatoms[i].molindex=asymm[i].molindex;
  }
  showbonds.clear();
  showbonds=connecting(showatoms);
}

void Molecule::grow(){
  /*! completes covalent bound molecule fragments. 
   *
   * */
  showatoms.clear();
  for (int i=0; i<asymm.size();i++){
    showatoms.append(asymm[i]);
    showatoms[i].molindex=asymm[i].molindex;
  }
  showbonds.clear();
  QStringList brauchSymm;
  QString bs;
  brauchSymm.clear();
  V3 prime,dp,D,floorD;
  double dk,dddd;
  for (int k =0; k<sdm.size();k++){
    if ((sdm.at(k).covalent)){
      if ((asymm[sdm.at(k).a1].molindex<1)) continue;
      for (int n=0;n<cell.symmops.size();  n++){
        if (((asymm[sdm.at(k).a1].part!=0)&&(asymm[sdm.at(k).a2].part!=0)&&(asymm[sdm.at(k).a1].part!=asymm[sdm.at(k).a2].part)))continue;
        if ((asymm[sdm.at(k).a1].an==asymm[sdm.at(k).a2].an)&&(asymm[sdm.at(k).a1].an==0)) continue;
        prime=cell.symmops.at(n) * asymm[sdm.at(k).a1].frac + cell.trans.at(n);
        D=prime - asymm[sdm.at(k).a2].frac+ V3(0.5,0.5,0.5) ;
        floorD=V3(floor(D.x),floor(D.y),floor(D.z));
        dp=D - floorD - V3(0.5,0.5,0.5);
        if ((n==0)&&(V3(0,0,0)==floorD)) continue;
        dk=fl(dp.x,dp.y,dp.z);
        //printf ("%f n%d\n",dk,n);
        dddd=(sdm.at(k).covalent)?(sdm.at(k).d+0.2):0;
        // an<0 means hydrogen atom (this is a hydrogen contact):
        if ((asymm[sdm.at(k).a1].an<0)&&(asymm[sdm.at(k).a2].an<0)) dddd=1.8;//||
        if ( (dk>0.001)&&(dddd>=dk)) {
          bs=QString("%1_%2%3%4:%5,").arg(n+1).arg(5-static_cast<int>(floorD.x)).arg(5-static_cast<int>(floorD.y)).arg(5-static_cast<int>(floorD.z)).arg(asymm[sdm.at(k).a1].molindex);
          //	     printf("%s %s %g %g\n",asymm[sdm.at(k).a1].Label.toStdString().c_str(),asymm[sdm.at(k).a2].Label.toStdString().c_str(),sdm.at(k).d,dddd);
          //             printf("(grow)%s\n",bs.toStdString().c_str());
          if  ((!brauchSymm.contains(bs))) {
            brauchSymm.append(bs);
          }
        }
      }
    }
  }
  for (int i=0; i<bindatoms.size(); i++){
    MyAtom a1, a2;
    bool bo=false;
    int sy1=-1,sy2=-1;
    if (bindatoms.at(i).Lab1.contains("_$"))sy1=bindatoms.at(i).Lab1.section("_$",1,1).toInt(&bo);
    if (bindatoms.at(i).Lab2.contains("_$"))sy2=bindatoms.at(i).Lab2.section("_$",1,1).toInt(&bo);
    bo=false;
    if (bindatoms.at(i).Lab1.contains('_'))a1.resiNr=bindatoms.at(i).Lab1.section('_',1,1).toInt(&bo);
    if (!bo) a1.resiNr=-1;
    if (bindatoms.at(i).Lab2.contains('_'))a2.resiNr=bindatoms.at(i).Lab2.section('_',1,1).toInt(&bo);
    if (!bo)a2.resiNr=-1;
    a1.Label = bindatoms.at(i).Lab1.section('_',0,0);
    a2.Label = bindatoms.at(i).Lab2.section('_',0,0);
    int frag1=0,frag2=0;
    for (int ii=0; ii<asymm.size();ii++){
      if ((a1.resiNr>-1)&& (asymm[ii].resiNr!=a1.resiNr)) continue;
      if ((a2.resiNr>-1)&& (asymm[ii].resiNr!=a2.resiNr)) continue;
      if (asymm[ii].Label.section('_',0,0).toUpper() == a1.Label.toUpper()){
        frag1=asymm[ii].molindex;
      }
      if (asymm[ii].Label.section('_',0,0).toUpper() == a2.Label.toUpper()){
        frag2=asymm[ii].molindex;
      }
    }
    QString ss;
    QList<Matrix> smm;
    QList<V3> tmm;
    int syy1=-1;
    if (sy1>-1) syy1=labelEQIV.indexOf(QString("$%1").arg(sy1));
    if (syy1>-1) {
      smm.append(symmopsEQIV.at(syy1));
      tmm.append(transEQIV.at(syy1));
    }
    int syy2=-1;
    if (sy2>-1) syy2=labelEQIV.indexOf(QString("$%1").arg(sy2));
    if (syy2>-1) {
      smm.append(symmopsEQIV.at(syy2));
      tmm.append(transEQIV.at(syy2));
    }
    for (int ui=0; ui<smm.size(); ui++){
      if (cell.symmops.contains(smm.at(ui))) {
        V3 r;
        r.x=fmod(tmm.at(ui).x+10,1.0);
        r.y=fmod(tmm.at(ui).y+10,1.0);
        r.z=fmod(tmm.at(ui).z+10,1.0);

        for (int is=0; is<cell.symmops.size();is++){
          if ((cell.symmops.at(is)==smm.at(ui))&&(cell.trans.at(is)==r))
            for (int yi=1; yi<=5; yi++){
              ss=QString("%1_%2%3%4:%5,")
                .arg(is+1)
                .arg(static_cast<int>(tmm.at(ui).x-r.x) +5)
                .arg(static_cast<int>(tmm.at(ui).y-r.y) +5)
                .arg(static_cast<int>(tmm.at(ui).z-r.z) +5)
                .arg(((sy1>-1)&&(!ui))?frag1:frag2);
              if  (!brauchSymm.contains(ss)) {
                brauchSymm.append(ss);
                eqivMap[labelEQIV.at(((sy1>-1)&&(!ui))?sy1-1:sy2-1)]=brauchSymm.size();
                printf("%d ==> %d %s\n",sy2-1, brauchSymm.size(),labelEQIV.at(((sy1>-1)&&(!ui))?sy1-1:sy2-1).toStdString().c_str());
              }
            }
        }
      }
    }
  }

  /*  for (int i=0; i<sdm.size();i++){
      if ((sdm.at(i).sn==0)&&(V3(0,0,0)==sdm.at(i).floorD)) continue;
      double dddd=((Kovalenz_Radien[asymm[sdm.at(i).a1].an]+ Kovalenz_Radien[asymm[sdm.at(i).a2].an])*0.01+0.5);
      if (dddd>sdm.at(i).d){
      bs=QString("%1_%2%3%4:%5,").arg(sdm.at(i).sn+1).arg(5-(int)sdm.at(i).floorD.x).arg(5-(int)sdm.at(i).floorD.y).arg(5-(int)sdm.at(i).floorD.z).arg(asymm[sdm.at(i).a1].molindex);
      / *        printf("%s %f %s %s\n",bs.toStdString().c_str(),sdm.at(i).d,
      asymm[sdm.at(i).a1].Label.toStdString().c_str(),
      asymm[sdm.at(i).a2].Label.toStdString().c_str());
   * /
   if  (!brauchSymm.contains(bs)) {
   brauchSymm.append(bs);
   }
   }
   } */
  //packer(brauchSymm);
  usedSymmetry+=brauchSymm;    
  //complete();
  // printf("GROW\n");
  packer(brauchSymm);
  showbonds.clear();
  showbonds=connecting(showatoms);
}

    void Molecule::grow_plus(){
      /*! completes covalent bound molecule fragments. 
       *
       * */
      showatoms.clear();
      for (int i=0; i<asymm.size();i++){
        showatoms.append(asymm[i]);
        showatoms[i].molindex=asymm[i].molindex;
      }
      showbonds.clear();
      QStringList brauchSymm;
      QString bs;
      brauchSymm.clear();
      V3 prime,dp,D,floorD;
      double dk,dddd;
      for (int k =0; k<sdm.size();k++){
        if ((sdm.at(k).covalent)){
          if ((asymm[sdm.at(k).a1].molindex<1)) continue;
          for (int n=0;n<cell.symmops.size();  n++){
            if (((asymm[sdm.at(k).a1].part!=0)&&(asymm[sdm.at(k).a2].part!=0)&&(asymm[sdm.at(k).a1].part!=asymm[sdm.at(k).a2].part)))continue;
            if ((asymm[sdm.at(k).a1].an==asymm[sdm.at(k).a2].an)&&(asymm[sdm.at(k).a1].an==0)) continue;
            prime=cell.symmops.at(n) * asymm[sdm.at(k).a1].frac + cell.trans.at(n);
            D=prime - asymm[sdm.at(k).a2].frac+ V3(0.5,0.5,0.5) ;
            floorD=V3(floor(D.x),floor(D.y),floor(D.z));
            dp=D - floorD - V3(0.5,0.5,0.5);
            if ((n==0)&&(V3(0,0,0)==floorD)) continue;
            dk=fl(dp.x,dp.y,dp.z);
            //printf ("%f n%d\n",dk,n);
            dddd=(sdm.at(k).covalent)?(sdm.at(k).d+0.2):0;
            if ((asymm[sdm.at(k).a1].an<0)&&(asymm[sdm.at(k).a2].an<0)) dddd=1.8;//||
            if ( (dk>0.001)&&(dddd>=dk)) {
              bs=QString("%1_%2%3%4:%5,").arg(n+1).arg(5-static_cast<int>(floorD.x)).arg(5-static_cast<int>(floorD.y)).arg(5-static_cast<int>(floorD.z)).arg(asymm[sdm.at(k).a1].molindex);
              //	     printf("%s %s %g %g\n",asymm[sdm.at(k).a1].Label.toStdString().c_str(),asymm[sdm.at(k).a2].Label.toStdString().c_str(),sdm.at(k).d,dddd);
              //             printf("(grow)%s\n",bs.toStdString().c_str());
              if  ((!brauchSymm.contains(bs))) {
                brauchSymm.append(bs);
              }
            }
          }
        }
      }
      for (int i=0; i<bindatoms.size(); i++){
        MyAtom a1, a2;
        bool bo=false;
        int sy1=-1,sy2=-1;
        if (bindatoms.at(i).Lab1.contains("_$"))sy1=bindatoms.at(i).Lab1.section("_$",1,1).toInt(&bo);
        if (bindatoms.at(i).Lab2.contains("_$"))sy2=bindatoms.at(i).Lab2.section("_$",1,1).toInt(&bo);
        bo=false;
        if (bindatoms.at(i).Lab1.contains('_'))a1.resiNr=bindatoms.at(i).Lab1.section('_',1,1).toInt(&bo);
        if (!bo) a1.resiNr=-1;
        if (bindatoms.at(i).Lab2.contains('_'))a2.resiNr=bindatoms.at(i).Lab2.section('_',1,1).toInt(&bo);
        if (!bo)a2.resiNr=-1;
        a1.Label = bindatoms.at(i).Lab1.section('_',0,0);
        a2.Label = bindatoms.at(i).Lab2.section('_',0,0);
        int frag1=0,frag2=0;
        for (int ii=0; ii<asymm.size();ii++){
          if ((a1.resiNr>-1)&& (asymm[ii].resiNr!=a1.resiNr)) continue;
          if ((a2.resiNr>-1)&& (asymm[ii].resiNr!=a2.resiNr)) continue;
          if (asymm[ii].Label.section('_',0,0).toUpper() == a1.Label.toUpper()){
            frag1=asymm[ii].molindex;
          }
          if (asymm[ii].Label.section('_',0,0).toUpper() == a2.Label.toUpper()){
            frag2=asymm[ii].molindex;
          }
        }
        QString ss;
        QList<Matrix> smm;
        QList<V3> tmm;
        int syy1=-1;
        if (sy1>-1) syy1=labelEQIV.indexOf(QString("$%1").arg(sy1));
        if (syy1>-1) {
          smm.append(symmopsEQIV.at(syy1));
          tmm.append(transEQIV.at(syy1));
        }
        int syy2=-1;
        if (sy2>-1) syy2=labelEQIV.indexOf(QString("$%1").arg(sy2));
        if (syy2>-1) {
          smm.append(symmopsEQIV.at(syy2));
          tmm.append(transEQIV.at(syy2));
        }
        for (int ui=0; ui<smm.size(); ui++)
          if (cell.symmops.contains(smm.at(ui))) {
            V3 r;
            r.x=fmod(tmm.at(ui).x+10,1.0);
            r.y=fmod(tmm.at(ui).y+10,1.0);
            r.z=fmod(tmm.at(ui).z+10,1.0);

            for (int is=0; is<cell.symmops.size();is++){
              if ((cell.symmops.at(is)==smm.at(ui))&&(cell.trans.at(is)==r))
                for (int yi=1; yi<=5; yi++){
                  ss=QString("%1_%2%3%4:%5,")
                    .arg(is+1)
                    .arg(static_cast<int>(tmm.at(ui).x-r.x) +5)
                    .arg(static_cast<int>(tmm.at(ui).y-r.y) +5)
                    .arg(static_cast<int>(tmm.at(ui).z-r.z) +5)
                    .arg(((sy1>-1)&&(!ui))?frag1:frag2);
                  if  (!brauchSymm.contains(ss)) {
                    brauchSymm.append(ss);
                    //		printf("### %d  = %d\n",((sy1>-1)&&(!ui))?sy1:sy2,brauchSymm.size());
                    eqivMap[labelEQIV.at(((sy1>-1)&&(!ui))?sy1-1:sy2-1)]=brauchSymm.size();
                    //		qDebug()<<ss;
                  }
                }
            }
          }
      }

      /*  for (int i=0; i<sdm.size();i++){
          if ((sdm.at(i).sn==0)&&(V3(0,0,0)==sdm.at(i).floorD)) continue;
          double dddd=((Kovalenz_Radien[asymm[sdm.at(i).a1].an]+ Kovalenz_Radien[asymm[sdm.at(i).a2].an])*0.01+0.5);
          if (dddd>sdm.at(i).d){
          bs=QString("%1_%2%3%4:%5,").arg(sdm.at(i).sn+1).arg(5-(int)sdm.at(i).floorD.x).arg(5-(int)sdm.at(i).floorD.y).arg(5-(int)sdm.at(i).floorD.z).arg(asymm[sdm.at(i).a1].molindex);
          / *        printf("%s %f %s %s\n",bs.toStdString().c_str(),sdm.at(i).d,
          asymm[sdm.at(i).a1].Label.toStdString().c_str(),
          asymm[sdm.at(i).a2].Label.toStdString().c_str());
       * /
       if  (!brauchSymm.contains(bs)) {
       brauchSymm.append(bs);
       }
       }
       } */
      //packer(brauchSymm);
      usedSymmetry+=brauchSymm;    
      complete();
      showbonds.clear();
      showbonds=connecting(showatoms);
    }

    void Molecule::packInLimits(double ami, double ama, double bmi, double bma, double cmi, double cma){//!< Pack inside given limits (eg multiple unit cells)

      MyAtom  newAtom;
      newAtom.hidden=0;
      V3 prime,floorD;
      double dawars=1000.0,dl;
      for (int n=0; n<cell.symmops.size();n++){
        for (int i=0; i<asymm.size();i++){
          //if ((asymm.at(i).an>-1)&&(asymm.at(i).molindex>0))
          {
            prime=cell.symmops.at(n) * asymm.at(i).frac + cell.trans.at(n);
            floorD=V3(floor(prime.x),floor(prime.y),floor(prime.z));
            //prime=prime -floorD;
            dawars=1000.0;
            for (int g=0; g<showatoms.size();g++){
              if ((asymm[i].an*showatoms[g].an<0)) continue;
              dl=fl(prime.x-showatoms[g].frac.x, prime.y-showatoms[g].frac.y, prime.z-showatoms[g].frac.z);
              dawars=(dl<dawars)?dl:dawars;
            }
            if (dawars<0.01) continue; 
            for (int hh=-11; hh<12; hh++)
              for (int kk=-11; kk<12; kk++)
                for (int ll=-11; ll<12; ll++){
                  if ((prime.x+hh)<ami) continue;
                  if ((prime.y+kk)<bmi) continue;
                  if ((prime.z+ll)<cmi) continue;
                  if ((prime.x+hh)>ama) continue;
                  if ((prime.y+kk)>bma) continue;
                  if ((prime.z+ll)>cma) continue;
                  newAtom.part=asymm[i].part;
                  newAtom.frac=prime+V3(static_cast<double>(hh),static_cast<double>(kk),static_cast<double>(ll));
                  frac2kart(newAtom.frac,newAtom.pos);
                  newAtom.an=asymm[i].an;
                  newAtom.molindex=asymm[i].molindex;
                  newAtom.peakHeight=asymm[i].peakHeight;
                  newAtom.resiNr=asymm[i].resiNr;

                  newAtom.sof=asymm[i].sof;
                  newAtom.Label=asymm[i].Label;
                  newAtom.isIso=asymm[i].isIso;
                  newAtom.ufiso_org=asymm[i].ufiso_org;

                  newAtom.sof_org=asymm[i].sof_org;
                  newAtom.symmGroup=1;
                  newAtom.sg=n;
                  //printf("A%d %d %d %d %g %g %g\n", n, hh, kk, ll, floorD.x, floorD.y, floorD.z);
                  newAtom.scod = genSymCode(n,hh,kk,ll);//genSymCode(n,(int)(hh+floorD.x),(int)(kk+floorD.x),(int)(ll+floorD.x));//genSymCode(n,hh,kk,ll);
                  newAtom.auidx = i;
                  if ((asymm[i].uc.m12==0.0)&&(asymm[i].uc.m23==0.0)&&(asymm[i].uc.m13==0.0)){
                    newAtom.uf.m11=newAtom.uc.m11=newAtom.uc.m22=newAtom.uc.m33=asymm[i].uf.m11;
                    newAtom.uc.m12=newAtom.uc.m13=newAtom.uc.m23=newAtom.uc.m21=newAtom.uc.m31=newAtom.uc.m32=0.0;
                  }
                  else {
                    Usym(asymm[i].uf,cell.symmops[n],newAtom.uf);
                    Uf2Uo(newAtom.uf,newAtom.uc);
                  }
                  showatoms.append(newAtom);
                }//ll
          }//packable 
          //    else if ((asymm.at(i).an<0)&&(n==0)) showatoms.append(asymm[i]);//non atoms BCPs, Dummys, etc.
        }//i++
      }//n++
      showbonds.clear();
      showbonds=connecting(showatoms);

    }

    CEnvironment Molecule::packInLimits(CEnvironment &au, double ami, double ama, double bmi, double bma, double cmi, double cma){//!< Pack inside given limits (eg multiple unit cells)
      MyAtom  newAtom;
      CEnvironment sa;
      newAtom.hidden=0;
      V3 prime,floorD;
      //double dawars=1000.0,dl;
      for (int n=0; n<cell.symmops.size();n++){
        for (int i=0; i<au.size();i++){
          prime=cell.symmops.at(n) * au.at(i).frac + cell.trans.at(n);
          floorD=V3(floor(prime.x),floor(prime.y),floor(prime.z));
          prime=prime -floorD;
          /*dawars=1000.0;
            for (int g=0; g<sa.size();g++){
            dl=fl(prime.x-sa[g].frac.x, prime.y-sa[g].frac.y, prime.z-sa[g].frac.z);
            dawars=(dl<dawars)?dl:dawars;
            }
            if (dawars<0.01) continue;*/
          for (int hh=-11; hh<12; hh++)
            for (int kk=-11; kk<12; kk++)
              for (int ll=-11; ll<12; ll++){
                if ((prime.x+hh)<ami) continue;
                if ((prime.y+kk)<bmi) continue;
                if ((prime.z+ll)<cmi) continue;
                if ((prime.x+hh)>ama) continue;
                if ((prime.y+kk)>bma) continue;
                if ((prime.z+ll)>cma) continue;
                newAtom.frac=prime+V3(static_cast<double>(hh), static_cast<double>(kk), static_cast<double>(ll));
                frac2kart(newAtom.frac,newAtom.pos);
                newAtom.an=au[i].an;
                newAtom.peakHeight=au[i].peakHeight;
                newAtom.sof=au[i].sof;
                newAtom.Label=au[i].Label;
                newAtom.isIso=au[i].isIso;
                newAtom.ufiso_org=au[i].ufiso_org;
                newAtom.sof_org=au[i].sof_org;
                newAtom.symmGroup=1;
                newAtom.sg=n;
                newAtom.scod = genSymCode(n,hh,kk,ll);
                printf("B%d %d %d %d %g %g %g\n", n, hh, kk, ll, floorD.x, floorD.y, floorD.z);
                //newAtom.scod = genSymCode(n,(int)(hh-floorD.x),(int)(kk-floorD.x),(int)(ll-floorD.x));
                newAtom.auidx = i;
                Usym(au[i].uf,cell.symmops[n],newAtom.uf);
                sa.append(newAtom);
              }//ll
        }//i++
      }//n++
      return sa;
    }




    void Molecule::fillCell(){
      /*! Fills the unit cell with molecules.
       *
       * */
      packInLimits(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
      return;

      showatoms.clear();
      for (int i=0; i<asymm.size();i++){
        showatoms.append(asymm[i]);
        showatoms[i].molindex=asymm[i].molindex;
      }
      showbonds.clear();
      QStringList brauchSymm;
      QString bs;
      brauchSymm.clear();
      V3 prime,dp,D,floorD;
      for (int k =0; k<sdm.size();k++){
        if (sdm.at(k).covalent){
          if (asymm[sdm.at(k).a1].molindex<1) continue;
          for (int n=0;n<cell.symmops.size();  n++){
            D=prime=cell.symmops.at(n) * asymm[sdm.at(k).a1].frac + cell.trans.at(n);
            floorD=V3(floor(D.x),floor(D.y),floor(D.z));
            if ((n==0)&&(V3(0,0,0)==floorD)) continue;
            bs=QString("%1_%2%3%4:%5,").arg(n+1).arg(5-static_cast<int>(floorD.x)).arg(5-static_cast<int>(floorD.y)).arg(5-static_cast<int>(floorD.z)).arg(asymm[sdm.at(k).a1].molindex);
            if  (!brauchSymm.contains(bs)) {
              brauchSymm.append(bs);
            }
          }
        }
      }
      packer(brauchSymm);
      showbonds.clear();
      showbonds=connecting(showatoms);
    }

    void Molecule::expandAll(){
      /*! like Molecule.expandAt but this searches vor neighboring molecules arround every atom of the asymmetric unit and adds them to the showatoms list.
      */
      //if (index>=showatoms.size()) return;
      //  V3 expander= showatoms.at(index).frac;
      V3 prime,floorD,D;
      double SUCHRAD=3.2;
      for (int i=0; i<asymm.size();i++){//Rapunzel lass dein H herunter
        if (asymm.at(i).an==0) {
          SUCHRAD-=0.7;
          break;
        }
      }
      //  printf(" Rapunzels H ist %f lang\n",SUCHRAD);
      QStringList brauchSymm;
      QString bs;
      brauchSymm.clear();
      //  brauchSymm+=(usedSymmetry);
      double dk;
      for (int j=0; j<asymm.size();j++){
        if (asymm.at(j).an<0) continue;
        for (int n=0; n<cell.symmops.size();n++){
          for (int i=0; i<asymm.size();i++){
            if (asymm.at(i).an<0) continue;
            prime=cell.symmops.at(n) * asymm.at(i).frac + cell.trans.at(n);
            D=prime - asymm.at(j).frac + V3(0.5,0.5,0.5) ;
            floorD=V3(floor(D.x),floor(D.y),floor(D.z));
            if ((n==0)&&(floorD==V3(0,0,0))) {continue;}
            D=D - floorD - V3(0.5,0.5,0.5);
            dk=fl(D.x,D.y,D.z);
            if (dk<SUCHRAD){
              bs=QString("%1_%2%3%4:%5,").arg(n+1).arg(5-static_cast<int>(floorD.x)).arg(5-static_cast<int>(floorD.y)).arg(5-static_cast<int>(floorD.z)).arg(asymm.at(i).molindex);
              if  ((!brauchSymm.contains(bs))&&(asymm.at(i).molindex>0)) {
                brauchSymm.append(bs);
              }
            }
          }
        }
      }
      showatoms.clear();
      for (int i=0; i<asymm.size();i++){
        showatoms.append(asymm[i]);
        showatoms[i].molindex=asymm[i].molindex;
      }
      packer(brauchSymm);
      showbonds.clear();
      showbonds=connecting(showatoms);
    }


    void Molecule::complete(){
      /*! If the current visible atoms are covalently bound to  symmetry equivalents the latter are added to the showatoms list.
       * */
      //if (index>=showatoms.size()) return;
      //  V3 expander= showatoms.at(index).frac;
      V3 prime,floorD,D;
      QStringList brauchSymm;
      QString bs;
      //  brauchSymm.clear();
      brauchSymm+=(usedSymmetry);
      double dk;
      int s,h,k,l,symmgroup;
      V3 pos0;
      for (int sy=0; sy<usedSymmetry.size(); sy++){ 
        sscanf(usedSymmetry.at(sy).toLatin1(),"%d_%1d%1d%1d:%d",&s,&h,&k,&l,&symmgroup);
        h-=5;
        k-=5;
        l-=5;
        s--;
        for (int j=0; j<asymm.size();j++){
          if (asymm.at(j).an<0) continue;
          for (int n=0; n<cell.symmops.size();n++){
            for (int i=0; i<asymm.size();i++){
              if (asymm.at(i).an<1) continue;
              if (asymm[i].molindex!=symmgroup) continue;
              pos0=cell.symmops.at(s)*asymm.at(i).frac+cell.trans.at(s)+V3(h,k,l);
              prime=cell.symmops.at(n) * asymm.at(j).frac + cell.trans.at(n);
              D=prime - pos0 + V3(0.5,0.5,0.5) ;
              floorD=V3(floor(D.x),floor(D.y),floor(D.z));
              if ((n==0)&&(floorD==V3(0,0,0))) {continue;}
              D=D - floorD - V3(0.5,0.5,0.5);
              dk=fl(D.x,D.y,D.z);
              double dddd=qMin(1.7,(Kovalenz_Radien[asymm.at(i).an]+ Kovalenz_Radien[asymm.at(j).an])*0.012);
              if (dk<dddd){
                bs=QString("%1_%2%3%4:%5,").arg(n+1).arg(5-static_cast<int>(floorD.x)).arg(5-static_cast<int>(floorD.y)).arg(5-static_cast<int>(floorD.z)).arg(asymm.at(i).molindex);
                if  ((!brauchSymm.contains(bs))&&(asymm.at(i).molindex>0)) {
                  brauchSymm.append(bs);
                }
              }
            }
          }
        }
      }
      showatoms.clear();
      for (int i=0; i<asymm.size();i++){
        showatoms.append(asymm[i]);
        showatoms[i].molindex=asymm[i].molindex;
      }
      packer(brauchSymm);
      showbonds.clear();
      showbonds=connecting(showatoms);
    }

    void Molecule::expandAt(int index){
      /*! searches SUCHRAD=3.2 A around the given atom for symmetry equivalent molecules and adds them to the showatoms list.
       *  if Hydrogens are inside the structure the SUCHRAD is reduced about 0.7 A.
       *  @param index index of the specified atom from the current showatoms list.
       */
      if (index>=showatoms.size()) return;
      V3 expander= showatoms.at(index).frac;
      V3 prime,floorD,D;
      double SUCHRAD=3.2;
      for (int i=0; i<asymm.size();i++){//Rapunzel lass dein H herunter
        if (asymm.at(i).an==0) {//haha if there are H atoms th search radius is reduced by 0.7 A
          SUCHRAD-=0.7;
          break;
        }
      }
      //printf(" Rapunzels H ist %f lang\n",SUCHRAD);
      QStringList brauchSymm;
      QString bs;
      brauchSymm.clear();
      brauchSymm+=(usedSymmetry);
      double dk;
      for (int n=0; n<cell.symmops.size();n++){
        for (int i=0; i<asymm.size();i++){
          prime=cell.symmops.at(n) * asymm.at(i).frac + cell.trans.at(n);
          D=prime - expander + V3(0.5,0.5,0.5) ;
          floorD=V3(floor(D.x),floor(D.y),floor(D.z));
          if ((n==0)&&(floorD==V3(0,0,0))) {continue;}
          D=D - floorD - V3(0.5,0.5,0.5);
          dk=fl(D.x,D.y,D.z);
          if (dk<SUCHRAD){
            bs=QString("%1_%2%3%4:%5,").arg(n+1).arg(5-static_cast<int>(floorD.x)).arg(5-static_cast<int>(floorD.y)).arg(5-static_cast<int>(floorD.z)).arg(asymm.at(i).molindex);
            if  ((!brauchSymm.contains(bs))&&(asymm.at(i).molindex>0)) {
              brauchSymm.append(bs);
            }
          }
        }
      }
      showatoms.clear();
      for (int i=0; i<asymm.size();i++){
        showatoms.append(asymm[i]);
        showatoms[i].molindex=asymm[i].molindex;
      }
      packer(brauchSymm);
      showbonds.clear();
      showbonds=connecting(showatoms);
    }

    void Molecule::packer(QStringList brauchSymm){
      /*! Packs symmetry equivalent atoms according to the given list of internal symmetry codes and adds them to the showatoms list
       * HumanSymmetry is feeeded with a human readalble list of used symmetry.
       * @param brauchSymm list of internal symmetry codes.
       */
      //  printf("MOEBEL? packer\n");
      usedSymmetry.clear();
      usedSymmetry+=(brauchSymm);
      MyAtom  newAtom;
      newAtom.hidden=0;
      int s,h,k,l,gibscho=0,symmgroup;
      //      balken->setMinimum(0);
      //      balken->setMaximum(brauchSymm.size());
      //      balken->show();
      HumanSymmetry=QString("Used Symmetry:<br>%1").arg(symmcode2human(brauchSymm));
      QString pre,suff;
      for (int j=0;j<brauchSymm.size();j++){
        //	balken->setValue(j);

        sscanf(brauchSymm.at(j).toLatin1(),"%d_%1d%1d%1d:%d",&s,&h,&k,&l,&symmgroup);
        //printf("BS:!%s! %d h%d k%d l%d sg%d\n",brauchSymm.at(j).toStdString().c_str(),s,h,k,l,symmgroup);
        h-=5;
        k-=5;
        l-=5;
        s--;
        for (int i=0;i<asymm.size();i++){
          // if (asymm[i].molindex) printf ("doch %d %d %s\n",asymm[i].an,asymm[i].molindex,asymm.at(i).Label.toStdString().c_str());
          if ((asymm[i].molindex==symmgroup)&&(asymm[i].an>-1)){
            newAtom.frac=cell.symmops.at(s)*asymm[i].frac+cell.trans.at(s)+V3(h,k,l);
            newAtom.part=asymm[i].part;
            frac2kart(newAtom.frac,newAtom.pos);
            newAtom.Label=QString("%1%2%3")
              .arg(asymm.at(i).Label)
              .arg(QString::fromUtf8("»"))
              .arg(j+1);// */
            newAtom.an=asymm[i].an;
            newAtom.symmGroup=j+1;
            newAtom.sg=s;
            newAtom.scod = genSymCode(s,h,k,l);
            newAtom.auidx = i;
            //printf("scod=%d %d %d %d %d\n",newAtom.scod, s,h+5,k+5,l+5);
            newAtom.sof_org=asymm[i].sof_org;
            newAtom.isIso=asymm[i].isIso;
            newAtom.ufiso_org=asymm[i].ufiso_org;
            newAtom.molindex=asymm[i].molindex;
            newAtom.ResiClass=asymm[i].ResiClass;
            newAtom.resiNr=asymm[i].resiNr;
            if ((asymm[i].uc.m12==0.0 )&&(asymm[i].uc.m23==0.0)&&(asymm[i].uc.m13==0.0)){
              newAtom.uc.m11=newAtom.uc.m22=newAtom.uc.m33=asymm[i].uf.m11;
              newAtom.uc.m12=newAtom.uc.m13=newAtom.uc.m23=newAtom.uc.m21=newAtom.uc.m31=newAtom.uc.m32=0.0;
              newAtom.uf=asymm[i].uf;
            }
            else {
              Usym(asymm[i].uf,cell.symmops[s],newAtom.uf);
              Uf2Uo(newAtom.uf,newAtom.uc);
            }
            gibscho=0;
            if (newAtom.part>=0){
              for(int gbt=0;gbt<showatoms.size();gbt++){
                if (showatoms.at(gbt).an<0) continue;
                if (showatoms.at(gbt).part!=newAtom.part)continue;
                if (fl(newAtom.frac.x-showatoms[gbt].frac.x,
                      newAtom.frac.y-showatoms[gbt].frac.y,
                      newAtom.frac.z-showatoms[gbt].frac.z)<0.2) gibscho=1;
              }
            }
            if (!gibscho) {
              showatoms.append(newAtom);
            }
          }// /*
          else if ((asymm[i].an==-1)&&(asymm[i].molindex==symmgroup)){
            newAtom.frac=cell.symmops.at(s)*asymm[i].frac+cell.trans.at(s)+V3(h,k,l);
            newAtom.part=asymm[i].part;        
            newAtom.isIso=asymm[i].isIso;
            frac2kart(newAtom.frac,newAtom.pos);
            newAtom.Label=asymm.at(i).Label;
            newAtom.an=asymm[i].an;
            newAtom.symmGroup=j+1;
            newAtom.sg=s;
            newAtom.scod = genSymCode(s,h,k,l);
            newAtom.auidx = i;
            newAtom.sof_org=asymm[i].sof_org;
            newAtom.sof=asymm[i].sof;
            newAtom.ufiso_org=asymm[i].ufiso_org;
            newAtom.molindex=asymm[i].molindex;
            newAtom.ResiClass=asymm[i].ResiClass;
            newAtom.resiNr=asymm[i].resiNr;
            newAtom.peakHeight=asymm[i].peakHeight;
            newAtom.orginalLine= asymm[i].orginalLine;
            gibscho=0;
            if (newAtom.part>=0){
              for(int gbt=0;gbt<showatoms.size();gbt++){
                if (fl(newAtom.frac.x-showatoms[gbt].frac.x,
                      newAtom.frac.y-showatoms[gbt].frac.y,
                      newAtom.frac.z-showatoms[gbt].frac.z)<0.1) gibscho=1;
              }
            }
            if (!gibscho) {
              showatoms.append(newAtom);
            }
          }// */
        }
      }
      //statusBar()->showMessage(tr("Neighbor search is finished"));
    }

    Matrix Molecule::b2u(const Matrix u){
      const double twopi2=2.0*M_PI*M_PI;
      Matrix erg;
      /*
       *         99      Atoms(ii)%Bij(1)=2*Pisqr*Atoms(ii)%Uij(1)*(rcell(1)*rcell(1))
       *        100      Atoms(ii)%Bij(2)=2*Pisqr*Atoms(ii)%Uij(2)*(rcell(2)*rcell(2))
       *        101      Atoms(ii)%Bij(3)=2*Pisqr*Atoms(ii)%Uij(3)*(rcell(3)*rcell(3))
       *        102      Atoms(ii)%Bij(4)=2*Pisqr*Atoms(ii)%Uij(4)*(rcell(1)*rcell(2))
       *        103      Atoms(ii)%Bij(5)=2*Pisqr*Atoms(ii)%Uij(5)*(rcell(1)*rcell(3))
       *        104      Atoms(ii)%Bij(6)=2*Pisqr*Atoms(ii)%Uij(6)*(rcell(2)*rcell(3))
       * */
      erg.m11 =          u.m11 /( twopi2 * cell.as*cell.as);
      erg.m22 =          u.m22 /( twopi2 * cell.bs*cell.bs);
      erg.m33 =          u.m33 /( twopi2 * cell.cs*cell.cs);
      erg.m12 = erg.m21 =u.m12 /( twopi2 * cell.as*cell.bs);
      erg.m13 = erg.m31 =u.m13 /( twopi2 * cell.as*cell.cs);
      erg.m23 = erg.m32 =u.m23 /( twopi2 * cell.bs*cell.cs);

      //  printf("Uij: %9.6f%9.6f%9.6f%9.6f%9.6f%9.6f\n",u.m11,u.m22,u.m33,u.m12,u.m13,u.m23);
      //  printf("Bij: %9.6f%9.6f%9.6f%9.6f%9.6f%9.6f\n",erg.m11,erg.m22,erg.m33,erg.m12,erg.m13,erg.m23);
      //printf("erg:\n%9.5f %9.5f %9.5f\n%9.5f %9.5f %9.5f\n%9.5f %9.5f %9.5f\n",erg.m11,erg.m12,erg.m13,erg.m12,erg.m22,erg.m23,erg.m13,erg.m23,erg.m33);
      //printf("erg2:\n%9.5f %9.5f %9.5f\n%9.5f %9.5f %9.5f\n%9.5f %9.5f %9.5f\n",erg2.m11,erg2.m12,erg2.m13,erg2.m12,erg2.m22,erg2.m23,erg2.m13,erg2.m23,erg2.m33);
      return erg;

    }

    Matrix Molecule::u2b(const Matrix u){
      const double twopi2=2.0*M_PI*M_PI;
      Matrix erg;
      /*
       *         99      Atoms(ii)%Bij(1)=2*Pisqr*Atoms(ii)%Uij(1)*(rcell(1)*rcell(1))
       *        100      Atoms(ii)%Bij(2)=2*Pisqr*Atoms(ii)%Uij(2)*(rcell(2)*rcell(2))
       *        101      Atoms(ii)%Bij(3)=2*Pisqr*Atoms(ii)%Uij(3)*(rcell(3)*rcell(3))
       *        102      Atoms(ii)%Bij(4)=2*Pisqr*Atoms(ii)%Uij(4)*(rcell(1)*rcell(2))
       *        103      Atoms(ii)%Bij(5)=2*Pisqr*Atoms(ii)%Uij(5)*(rcell(1)*rcell(3))
       *        104      Atoms(ii)%Bij(6)=2*Pisqr*Atoms(ii)%Uij(6)*(rcell(2)*rcell(3))
       * */
      erg.m11 =          u.m11 * twopi2 * cell.as*cell.as;
      erg.m22 =          u.m22 * twopi2 * cell.bs*cell.bs;
      erg.m33 =          u.m33 * twopi2 * cell.cs*cell.cs;
      erg.m12 = erg.m21 =u.m12 * twopi2 * cell.as*cell.bs;
      erg.m13 = erg.m31 =u.m13 * twopi2 * cell.as*cell.cs;
      erg.m23 = erg.m32 =u.m23 * twopi2 * cell.bs*cell.cs;

      //  printf("Uij: %9.6f%9.6f%9.6f%9.6f%9.6f%9.6f\n",u.m11,u.m22,u.m33,u.m12,u.m13,u.m23);
      //  printf("Bij: %9.6f%9.6f%9.6f%9.6f%9.6f%9.6f\n",erg.m11,erg.m22,erg.m33,erg.m12,erg.m13,erg.m23);
      //printf("erg:\n%9.5f %9.5f %9.5f\n%9.5f %9.5f %9.5f\n%9.5f %9.5f %9.5f\n",erg.m11,erg.m12,erg.m13,erg.m12,erg.m22,erg.m23,erg.m13,erg.m23,erg.m33);
      //printf("erg2:\n%9.5f %9.5f %9.5f\n%9.5f %9.5f %9.5f\n%9.5f %9.5f %9.5f\n",erg2.m11,erg2.m12,erg2.m13,erg2.m12,erg2.m22,erg2.m23,erg2.m13,erg2.m23,erg2.m33);
      return erg;

    }

    int Molecule::genSymCode(int s, int h, int k, int l){//!< creates a platon style symmetry code as an integer (s555)
      return s*1000+(h+5)*100+(k+5)*10+(l+5);
    }

    V3 Molecule::applySymCode(const V3 &frac, int scode){//!< applies a platon style symmetry code as an integer (s555)
      int s, h, k, l;
      s = scode / 1000;
      if (s >= cell.symmops.size())return frac;// silent error handling
      l = (scode % 10) - 5 ;
      k = ((scode % 100) / 10) - 5;   
      h = ((scode % 1000) / 100) - 5;   
      return cell.symmops.at(s) * frac + cell.trans.at(s) + V3(h, k, l);
    }

    int eccall=0;

    void Molecule::extendChain(QList<Ring> &rings,QList<int> &currentChain, QList<QList<int> > &neighbors,QSet<int> &nogo){
      eccall++;

      /*for (int k=0; k<currentChain.size(); k++){
        printf("%s-",asymm.at(currentChain.at(k)).Label.toStdString().c_str());
        }// */
      //printf("EC[%d]nogosize%d#current chain%d\n",eccall,nogo.size(),currentChain.size());

      int last=currentChain.last();
      for (int ri=0; ri<rings.size();ri++){
        if (rings.at(ri).done) continue;
        if (rings[ri].members.contains(last)){
          int bbv=currentChain.size();
          int strt=rings[ri].members.indexOf(last);
          for (int rri=0,rrj=0; rri<rings[ri].members.size(); rri++){
            rrj=(strt>=0)?(strt+rri)%rings[ri].members.size():rri;
            if (!currentChain.contains(rings[ri].members.at(rrj))) currentChain.append(rings[ri].members.at(rrj));
          }
          bbv-=currentChain.size();
          if (bbv) {
            rings[ri].done=true;
            //    printf("extRING%d %d %d %d\n",ri,bbv,currentChain.size(),rings[ri].members.size());
          }
        }
      }


      if (neighbors.at(last).size()==1) return;
      QList<int> c;
      int ic=-1;
      int manb=0;
      for (int i=0; i<neighbors.at(last).size(); i++){
        if (currentChain.contains(neighbors.at(last).at(i))) continue;
        if (nogo.contains(neighbors.at(last).at(i))) continue;
        int z=neighbors.at(neighbors.at(last).at(i)).size();
        if(z>manb){
          manb=z;
          ic=i;
        }
      }
      if (ic!=-1){
        c=currentChain;
        currentChain.append(neighbors.at(last).at(ic));
        extendChain(rings, currentChain, neighbors, nogo);
      }/*else{ // forking
         chains.append(c);
         int lc = chains.size() - 1;
         chains[lc].append(neighbors.at(last).at(i));
         extendChain(chains, chains[lc], neighbors, ng);
         ic++;
         ng.unite(chains[lc].toSet());      
         }// */
      c.clear();
    }

    void Molecule::inventNewLabels(QList<int> &result){

      /*
       * Idea / draft:
       * for atoms in a molecule
       * find chain ends  (E) neighbors 1
       * find chain atoms (C) neighbors 2   
       * find nodes       (N) neighbors >2
       * start at an E and append indexes to a QList chain
       * when a node N appears copy the list neighbors - 1 times (forking)
       * when no atoms left then compare number members in each chain
       * the longest chain wins
       * the longest part with carbons wins
       * the E with the greatest distance from center of mass wins
       * E-C-C-C-N-C-C-N-N-N-C-C-E
       *         |     | | |
       *         C-C-C-C E E
       *
       * 1-2-3-4-5-g-h-a-b-c-d-e-f
       *         |     | | |
       *         6-7-8-9 i j
       */  
      //return;
      QList<QList<int> > neighbors; 
      QList<QList<int> > chains;
      QList<int> emptyChain;
      QSet<int> nogo;
      QList<int> molsizes;
      QList<bool> moleIsDone;
      QList<double> avan, di2com;// avarage atomic number, distance to center of mass
      QList<V3> centerOfMass;
      QList<int>cmax;//maximal chain length
      for (int i=0; i<asymm.size(); i++){
        if (asymm.at(i).an<0 )continue;
        neighbors.append(emptyChain);
        if (asymm.at(i).molindex> molsizes.size()) {
          molsizes.append(0);
          cmax.append(0);
          moleIsDone.append(false);
          centerOfMass.append(V3(0.0,0.0,0.0));
          avan.append(999990.0);
          di2com.append(0.0);
        }
        molsizes[asymm.at(i).molindex-1]++;
        centerOfMass[asymm.at(i).molindex-1]+=asymm.at(i).pos;
      }
      for (int i=0; i<molsizes.size(); i++){
        centerOfMass[i]*=1.0/molsizes.at(i);
      }
      for (int i=0; i<sdm.size(); i++){
        if ((sdm.at(i).a1>=neighbors.size())||(sdm.at(i).a2>=neighbors.size())) continue;
        if (sdm.at(i).covalent) neighbors[sdm.at(i).a1].append(sdm.at(i).a2); 
      }
      // find rings first
      QList<QList<int> > rings; 
      QList<Tripel> tripels;
      Tripel t;//midcenter
      QList<Ring> ringe;
      QSet<int>nonMetals;
      //B, C, N, O, F,  Si, P, S, Cl, Ge, As, Se, Br,  Sb, Te, I
      nonMetals<<4<<5<<6<<7<<8<<13<<14<<15<<16<<31<<32<<33<<34<<50<<51<<52;
      for (int i=0; i<neighbors.size();i++){
        if (!nonMetals.contains(asymm.at(i).an)) continue; //no rings with metals
        if (neighbors.at(i).size()>1) 
          for (int j=0; j<neighbors.at(i).size();j++){
            for (int k=j+1; k<neighbors.at(i).size();k++){
              if (!nonMetals.contains(asymm.at(neighbors.at(i).at(j)).an)) continue; //no rings with metals
              if (!nonMetals.contains(asymm.at(neighbors.at(i).at(k)).an)) continue; //no rings with metals
              t.n[0]=i;
              t.n[1]=neighbors.at(i).at(j);
              t.n[2]=neighbors.at(i).at(k);
              double w = winkel(asymm.at(t.n[1]).pos - asymm.at(t.n[0]).pos, asymm.at(t.n[2]).pos - asymm.at(t.n[0]).pos) / 180.0 * M_PI;
              double scal = sin(w * 0.5) * sin(w * 0.5) / (sin(w) * sin(w));
              t.midcenter = asymm.at(t.n[0]).pos + scal * (asymm.at(t.n[1]).pos + asymm.at(t.n[2]).pos -(asymm.at(t.n[0]).pos * 2.0));
              tripels.append(t);
              //printf("tripel#%d %s %s %s\n",tripels.size(),asymm.at(t.n[1]).Label.toStdString().c_str(), asymm.at(t.n[0]).Label.toStdString().c_str(), asymm.at(t.n[2]).Label.toStdString().c_str());
              //      if (neighbors.at(i).size()==2) break;
            }
          }
      }
      printf("atoms%d tripel%d\n",neighbors.size(),tripels.size());
      const double RINGTEST=0.7;
      for (int i=0; i<tripels.size(); i++){
        /*  print the midcenter as q peaks
            V3 f;
            kart2frac(tripels.at(i).midcenter,f);
            printf("Q%03d 1 %9.5f %9.5f %9.5f 11.0 0.05 %5.2f\n",i,f.x,f.y,f.z,i*0.01);
            */
        for (int j=0; j<i; j++){
          double dist=Distance(tripels.at(i).midcenter,tripels.at(j).midcenter);
          if (dist<RINGTEST){
            bool alreadyHaveIt=false;
            for (int k = 0; k < ringe.size(); k++){
              //printf("j%d i%d k%d\n",j,i,k);
              //          if (ringe[k].isComplete()) continue;
              V3 r1=ringe[k].theCenter();
              double 
                d1=Distance(tripels.at(i).midcenter,r1),
                d2=Distance(tripels.at(j).midcenter,r1);
              if ((d1<RINGTEST)||(d2<RINGTEST)){
                ringe[k].addTripel(tripels.at(i));
                ringe[k].addTripel(tripels.at(j));
                //printf(":: %d %f %f ::\n", k, d1, d2);
                alreadyHaveIt=true;
                //  continue;
              }
              //else printf(":: %f %f ::\n",d1,d2);
            }
            if (!alreadyHaveIt){
              //printf("tic\n");
              ringe.append(Ring(tripels.at(i),tripels.at(j)));
              //ringe.last().printRing(this);
            }
          }//else {printf("DIST %g\n",dist);}
        //printf("j%d %d\n",j,i);
        }
        //printf("i%d \n",i);
      }
      //return;
      printf("Rings %d\n",ringe.size());
      //  int fix=0;
      for (int ri=0; ri<ringe.size(); ri++){
        bool fertig = ringe[ri].isComplete();
        if ((!fertig)&&(ringe[ri].members.size()<7)) {
          //ringe[ri].printRing(this);
          int u=ringe[ri].members.size();
          for (int w=0; w<u; w++) {
            for (int v=0; v<tripels.size();v++){
              if (
                  (ringe[ri].members.contains(tripels.at(v).n[0]))&&
                  (ringe[ri].members.contains(tripels.at(v).n[1]))&&
                  (ringe[ri].members.contains(tripels.at(v).n[2]))) {
                ringe[ri].addTripel(tripels.at(v));
                continue;
              }
            }
          }
        }
        //ringe[ri].printRing(this);
        fertig = ringe[ri].isComplete();
        ringe[ri].isItARealRing();
        fertig = ringe[ri].isComplete();
        //printf("%d %d %d\n",ringe[ri].t.size(),ringe[ri].members.size(),fertig);
        if (!fertig) {ringe.removeAt(ri);ri--;}
      }
      /*
         for (int ri=0; ri<ringe.size(); ri++){
         for (int rj=1; rj<ri; rj++){
         if (rj==ri)continue;
         if (Distance(ringe[ri].theCenter(),ringe[rj].theCenter())<0.2){ringe.removeAt(ri);ri--;break;}
         }
         }
         */
      printf("Rings %d are real and complete.\n",ringe.size());
      //  return;
      for (int mi=0; mi<centerOfMass.size();mi++){
        printf("Molecule#%d has %d atoms\n", mi, molsizes.at(mi));
        double away=-10.0;
        int itakethisring=0;
        for (int ri=0; ri<ringe.size(); ri++){
          if ((asymm.at(ringe.at(ri).members.at(0)).molindex-1)!=mi) continue;
          bool fertig = ringe[ri].isComplete();
          //    if (fertig) fix++;
          ringe[ri].printRing(this);
          printf("Ring #%d size = %d/%d is complete %d\n",ri,ringe.at(ri).t.size(),ringe.at(ri).members.size(),fertig);
          if (ringe.at(ri).members.size()==7){
            ringe[ri].debugRing(this);
            ringe[ri].isItARealRing();

          }
          double d=Distance(centerOfMass[asymm.at(ringe.at(ri).members.at(0)).molindex-1],ringe[ri].theCenter());
          V3 dir = (asymm.at(ringe.at(ri).members.at(0)).pos - ringe[ri].theCenter())%(asymm.at(ringe.at(ri).members.at(1)).pos - ringe[ri].theCenter());
          V3 com2r = ringe[ri].theCenter()-centerOfMass[asymm.at(ringe.at(ri).members.at(0)).molindex-1];
          double test= dir*com2r;
          QList<int> rcpy=ringe[ri].members;
          if (test>0.0){//flip ring order
            ringe[ri].members.clear();
            for (int kk=rcpy.size()-1; kk>=0; kk--){
              ringe[ri].members.append(rcpy.at(kk));
            }
            dir = (asymm.at(ringe.at(ri).members.at(0)).pos - ringe[ri].theCenter())%(asymm.at(ringe.at(ri).members.at(1)).pos - ringe[ri].theCenter());
            test= dir*com2r;
            //printf("fliped ");
          }
          ringe[ri].printRing(this);
          if (away<d){
            away=d;
            itakethisring=ri;
          }
          printf("AWAYOFCOM %f %d  test%g   \n",d,ri,test);
        }
        //printf("Complete Rings %d\n",fix);
        printf("Rings %d\n",ringe.size());
        if ((!ringe.isEmpty())&&((asymm.at(ringe.at(itakethisring).members.at(0)).molindex-1)==mi)){
          ringe[itakethisring].printRing(this);
          printf("most far AWAY from COM %f %d\n",away,itakethisring);
          int starthere=0,mn=0;
          double away2=1000;
          for (int i=0; i<ringe[itakethisring].members.size(); i++){
            if ((neighbors.at(ringe[itakethisring].members.at(i)).size()>=mn)&&
                (Distance(asymm.at(ringe[itakethisring].members.at(i)).pos,centerOfMass[asymm.at(ringe.at(itakethisring).members.at(i)).molindex-1])<away2)){
              mn=neighbors.at(ringe[itakethisring].members.at(i)).size();
              away2=Distance(asymm.at(ringe[itakethisring].members.at(i)).pos,centerOfMass[asymm.at(ringe.at(itakethisring).members.at(i)).molindex-1]);
              starthere=i;
            }
          }      
          chains.append(emptyChain);
          for (int i=0; i<ringe[itakethisring].members.size(); i++){
            chains.last().append(ringe[itakethisring].members.at(starthere%ringe[itakethisring].members.size()));
            starthere++;
          }
          ringe[itakethisring].done=true;
          //aneal ring
          for (int ci=0; ci<chains.last().size(); ci++){
            for (int ri=0; ri<ringe.size();ri++){
              if (ringe.at(ri).done) continue;
              if ((ringe[ri].members.contains(chains.last().at(ci)))&&(ringe[ri].members.contains(chains.last().at((ci+1)%chains.last().size())))){
                int bbv=chains.last().size();
                int ll=chains.last().last();
                int strt=ringe[ri].members.indexOf(ll);
                if (strt<0)strt = ringe[ri].members.indexOf(chains.last().at(ci));
                for (int rri=0,rrj=0; rri<ringe[ri].members.size(); rri++){
                  rrj=(strt>=0)?(strt+rri)%ringe[ri].members.size():rri;
                  if (!chains.last().contains(ringe[ri].members.at(rrj))) chains.last().append(ringe[ri].members.at(rrj));
                }          
                bbv-=chains.last().size();
                if (bbv) ringe[ri].done=true;
                //printf("RING%d %d %d %d %d\n",ri,bbv,ci,chains.last().size(),ringe[ri].members.size());
              }
            }
          }
          do{
            int ce=chains.last().size();
            for (int cc=0; cc<ce; cc++){
              for (int nn=0; nn<neighbors.at(chains.last().at(cc)).size();nn++){
                if (!chains.last().contains(neighbors.at(chains.last().at(cc)).at(nn))) {
                  chains.last().append(neighbors.at(chains.last().at(cc)).at(nn));
                  //printf("ext %s\n",asymm.at(neighbors.at(chains.last().at(cc)).at(nn)).Label.toStdString().c_str());
                  extendChain(ringe,chains.last(), neighbors,nogo);
                }
              }
            }
          }while(chains.last().size()<molsizes.at(mi));
          if (chains.last().size()==molsizes.at(mi))moleIsDone[mi]=true;
        }
        if ((!chains.isEmpty())&&(!result.contains(chains.last().last()))) {
          result.append(chains.last());
        }
      }//mi

      /*
         for (int cs=0; cs<chains.size(); cs++){
         chains[cs].clear();
         }
         chains.clear();
         */
      //for(int i=0; i<result.size(); i++){printf("%s[%d]=",asymm.at(result.at(i)).Label.toStdString().c_str(),i+1);if(!(i%10))printf("\n");}printf("\n");
      //return;
      //
      //printf("Molecules with rings are done\n");
      for (int mi=0; mi<molsizes.size(); mi++){
        if (moleIsDone.at(mi)) continue;
        //printf("mol%d %d %d %d resultsize%d chains%d\n",mi,molsizes.at(mi),molsizes.size(),centerOfMass.size(),
        //     result.size(),chains.size());
        int starthere=-1;
        double away=-10;
        for (int i=0; i<neighbors.size(); i++){
          if (asymm.at(i).an<0)continue;
          if ((asymm.at(i).molindex-1 == mi)&&(neighbors.at(i).size()==0)){
            result.append(i);//isolated atoms
            continue;
          }
          if ((asymm.at(i).molindex-1 == mi)&&(neighbors.at(i).size()==1)){//chain end
            //printf("%d %d %d %d %d\n",__LINE__,i,asymm.size(),neighbors.size());
            double d=Distance(asymm.at(i).pos,centerOfMass.at(mi));
            //printf("%d\n",__LINE__);
            if (d>away){
              away=d;
              starthere=i;
            }
          }
          // printf("%d\n",__LINE__);
        }
        //printf("%d\n",__LINE__);
        if (starthere<0) continue;
        //printf("start at %d if not done%d\n",starthere,!result.contains(starthere));

        if (!result.contains(starthere)){
          chains.append(emptyChain);
          chains.last().append(starthere);
          extendChain(ringe,chains.last(),neighbors,nogo);
        }
        do{
          int ce=chains.last().size();
          for (int cc=0; cc<ce; cc++){
            for (int nn=0; nn<neighbors.at(chains.last().at(cc)).size();nn++){
              if (!chains.last().contains(neighbors.at(chains.last().at(cc)).at(nn))) {
                chains.last().append(neighbors.at(chains.last().at(cc)).at(nn));
                //printf("ext %s\n",asymm.at(neighbors.at(chains.last().at(cc)).at(nn)).Label.toStdString().c_str());
                extendChain(ringe,chains.last(), neighbors,nogo);
              }
            }
          }
        }while(chains.last().size()<molsizes.at(mi));
        if ((!chains.isEmpty())&&(!result.contains(chains.last().last()))) {
          result.append(chains.last());
        }

      }
      for (int cs=0; cs<chains.size(); cs++){
        chains[cs].clear();
      }
      chains.clear();
    }

    /*

       if (neighbors.at(i).size()==1){//start a chain
       printf("start a chain @ %s\n",asymm.at(i).Label.toStdString().c_str() );
       chains.append(emptyChain);
       chains.last().append(i);
       chains.last().append(neighbors.at(i).at(0));
       extendChain(ringe,chains.last(), neighbors,nogo);
       }
       }
       for (int l = 0; l < chains.size(); l++){
       if (!chains.at(l).isEmpty()){
       int idx = asymm.at(chains.at(l).at(0)).molindex-1;
       cmax[idx]=(cmax.at(idx)>chains.at(l).size())?cmax.at(idx):chains.at(l).size();
       }
       }
       for (int l = 0; l < chains.size(); l++){
       if (!chains.at(l).isEmpty()){
       int idx = asymm.at(chains.at(l).at(0)).molindex-1;
       if (chains.at(l).size()<cmax.at(idx)) chains[l].clear();
       }
       if (!chains.at(l).isEmpty()){
       printf("[%d]: ",chains.at(l).size());
       for (int k = 0; k < chains.at(l).size(); k++){
       printf("%s-",asymm.at(chains.at(l).at(k)).Label.toStdString().c_str());
       }
       printf("\n");
       int ansum = 0;
       for (int k = 0; k < chains.at(l).size(); k++){
       ansum+= asymm.at(chains.at(l).at(k)).an + 1;
       printf("-%d-",neighbors.at(chains.at(l).at(k)).size());
       }
       int idx=asymm.at(chains.at(l).at(0)).molindex-1;
       double avan0=fabs((double)ansum / chains.at(l).size()-6.0);
       avan[idx] = fmin(avan[idx],avan0);
       double di2com0 = Distance(centerOfMass.at(asymm.at(chains.at(l).at(0)).molindex-1),asymm.at(chains.at(l).at(0)).pos);
       printf(" %f %d %f, %f [%d]\n",avan0, 
       idx,
       di2com0,
       avan[idx],
       idx
       );
       }
       }
       for (int l = 0; l < chains.size(); l++){
       if (!chains.at(l).isEmpty()){
       int idx = asymm.at(chains.at(l).at(0)).molindex-1;

       int ansum = 0;
       for (int k = 0; k < chains.at(l).size(); k++){
       ansum+= asymm.at(chains.at(l).at(k)).an + 1;
       }
       double avan0=fabs((double)ansum / chains.at(l).size()-6.0);
       if (avan0>avan.at(idx)) chains[l].clear();
       if (!chains.at(l).isEmpty()){
       double di2com0 = Distance(centerOfMass.at(asymm.at(chains.at(l).at(0)).molindex-1),asymm.at(chains.at(l).at(0)).pos);
       di2com[idx] = fmax(di2com0,di2com[idx]);
       }
       }
       }
       for (int l = 0; l < chains.size(); l++){
       if (!chains.at(l).isEmpty()){
       int idx = asymm.at(chains.at(l).at(0)).molindex-1;
       double di2com0 = Distance(centerOfMass.at(asymm.at(chains.at(l).at(0)).molindex-1),asymm.at(chains.at(l).at(0)).pos);
       if (di2com0<di2com[idx]) chains[l].clear();
       if (!chains.at(l).isEmpty()){

       printf("Winner is[%d]: ",chains.at(l).size());
       for (int k = 0; k < chains.at(l).size(); k++){
       printf("%s-",asymm.at(chains.at(l).at(k)).Label.toStdString().c_str());
       }
    printf("!== %f %d\n",avan.at(idx),idx);
}
}
}
*/




    double Molecule::ueq(const Matrix m){
      double erg=0;
      erg+=m.m11*cell.as*cell.a*cell.a*cell.as;
      erg+=m.m12*cell.as*cell.a*cell.b*cell.bs;
      erg+=m.m13*cell.as*cell.a*cell.c*cell.cs;
      erg+=m.m21*cell.bs*cell.b*cell.a*cell.as;
      erg+=m.m22*cell.bs*cell.b*cell.b*cell.bs;
      erg+=m.m23*cell.bs*cell.b*cell.c*cell.cs;
      erg+=m.m31*cell.cs*cell.c*cell.a*cell.as;
      erg+=m.m32*cell.cs*cell.c*cell.b*cell.bs;
      erg+=m.m33*cell.cs*cell.c*cell.c*cell.cs;
      erg*=1/3.0;
      return erg;
    }
