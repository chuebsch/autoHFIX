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
#ifndef MOLECULE_H
#define MOLECULE_H 1
#include <math.h>
#ifndef fmin
#define fmin(x, y) (((x) < (y)) ? (x) : (y))
#endif
#ifndef fmax
#define fmax(x, y) (((x) > (y)) ? (x) : (y))
#endif
#if defined _MSC_VER &&  _MSC_VER == 1500
#ifndef round 
#define round(x) (x<0?ceil((x)-0.5):floor((x)+0.5))
#endif
#endif
#include <QString>
#include <QStringList>
#include <QList>
#include <QtCore>
#ifndef M_PI
#define	M_PI		3.14159265358979323846
#endif
#ifndef M_PIf
#define	M_PIf		3.14159265358979323846f
#endif
//! V3 is a three dimensional vector in cartesian space
struct V3 {
  double x//! x is the X coordinate 
    ,y//! y is the Y coordinate 
    ,z//! z is the Z coordinate
    ;
  //  int rc;
  inline V3( void ){}
  inline V3( const double& _x, const double& _y, const double& _z ) : 
    x(_x), y(_y), z(_z)//!< initializer
    //,rc(0) 
  {

  }
  inline V3& operator *= ( const double& d ){
    x *= d;
    y *= d;
    z *= d;
    return *this;
  }//!< The *= operator to scale by a scalar
  inline V3& operator += ( const V3& v ){
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
  }//!< The += operator to add a V3  
  inline V3& operator += ( const double& v ){
    x += v;
    y += v;
    z += v;
    return *this;
  }//!< The += operator to add a scalar
};
inline V3 operator + ( const V3& v1, const V3& v2 ) {
  V3 t;
  t.x = v1.x + v2.x;
  t.y = v1.y + v2.y;
  t.z = v1.z + v2.z;
  return t;
}//!< The + operator to add two V3
inline V3 operator - ( const V3& v1, const V3& v2 ) {
  V3 t;
  t.x = v1.x - v2.x;
  t.y = v1.y - v2.y;
  t.z = v1.z - v2.z;
  return t;
}//!< The + operator to subtract two V3
inline V3 operator * ( const V3& v, const double& d ) {
  V3 t;
  t.x = v.x*d;
  t.y = v.y*d;
  t.z = v.z*d;
  return t;
}//!< The * to scale a V3
inline V3 operator * ( const double& d, const V3& v ) {
  V3 t;
  t.x = v.x*d;
  t.y = v.y*d;
  t.z = v.z*d;
  return t;
}//!< The * to scale a V3
inline V3 operator % ( const V3& v1, const V3& v2 ) {
  V3 t;
  t.x = v1.y*v2.z - v2.y*v1.z;
  t.y = v1.z*v2.x - v2.z*v1.x;
  t.z = v1.x*v2.y - v2.x*v1.y;
  return t;
}//!< The % operator the cross product of two V3
inline double operator * ( const V3& v1, const V3& v2 ) {
  return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}//!< The * operator the scalar product of two V3
inline double Norm( const V3& v ) {
  return v.x*v.x + v.y*v.y + v.z*v.z;
}//!< The squared lenght of a V3
inline double Distance( const V3& v1, const V3& v2 ) {
  return Norm(v1 - v2);
}//!< The squared distance between two V3
inline bool operator == (const V3& v1, const V3& v2 ) {
  //  return ((v1.x==v2.x)&&(v1.y==v2.y)&&(v1.z==v2.z));
  return (Distance(v1,v2)<0.005);
}
inline V3& Normalize( V3 v ) {
  static V3 erg=V3(1,0,0);
  if (Norm(v))  erg= (v * (1.0/sqrt(Norm(v)))); 
  return erg; 
}
//! Matrix is a 3 x 3 Matrix with all needed operators
struct Matrix{
  double m11, m21, m31, m12, m22, m32, m13, m23, m33;
  inline Matrix(void){}
  inline Matrix( const V3 &a, const V3 &b, const V3 &c):
    m11(a.x), m21(b.x), m31(c.x),
    m12(a.y), m22(b.y), m32(c.y),
    m13(a.z), m23(b.z), m33(c.z){;}
  inline Matrix( const double& x11, const double& x21, const double& x31,
      const double& x12, const double& x22, const double& x32,
      const double& x13, const double& x23, const double& x33):
    m11(x11), m21(x21), m31(x31),
    m12(x12), m22(x22), m32(x32),
    m13(x13), m23(x23), m33(x33){;}

};
inline Matrix transponse (Matrix a){//transponse
  return Matrix(
      a.m11, a.m12, a.m13,
      a.m21, a.m22, a.m23,
      a.m31, a.m32, a.m33);
}
inline double determinant (Matrix a){
  return a.m11*a.m22*a.m33 - a.m11*a.m23*a.m32 - a.m12*a.m21*a.m33 + a.m12*a.m23*a.m31 +a.m13*a.m21*a.m32 -a.m13*a.m22*a.m31;
}
inline Matrix inverse (Matrix A){
  double D = determinant(A);
  if (D==0) return A;
  D=1.0/D;
  return Matrix(
      D*(A.m22*A.m33-A.m23*A.m32),//x11
      D*(A.m13*A.m32-A.m12*A.m33),//x21
      D*(A.m21*A.m23-A.m13*A.m22),//x31
      D*(A.m23*A.m31-A.m21*A.m33),//x12
      D*(A.m11*A.m33-A.m13*A.m31),//x22
      D*(A.m13*A.m21-A.m11*A.m23),//x32
      D*(A.m21*A.m32-A.m22*A.m31),//x13
      D*(A.m12*A.m31-A.m11*A.m32),//x23
      D*(A.m11*A.m22-A.m12*A.m21)//x33
      );
}
inline bool operator == (const Matrix &a,const Matrix &b){
  return ((a.m11 == b.m11)&&(a.m21 == b.m21)&&(a.m31 == b.m31)&&
      (a.m12 == b.m12)&&(a.m22 == b.m22)&&(a.m23 == b.m23)&&
      (a.m13 == b.m13)&&(a.m32 == b.m32)&&(a.m33 == b.m33));
}
inline Matrix operator * (const Matrix &a,const Matrix &b){
  Matrix erg;
  erg.m11 = a.m11 * b.m11 + a.m21 * b.m12 + a.m31 * b.m13;
  erg.m21 = a.m11 * b.m21 + a.m21 * b.m22 + a.m31 * b.m23;
  erg.m31 = a.m11 * b.m31 + a.m21 * b.m32 + a.m31 * b.m33;

  erg.m12 = a.m12 * b.m11 + a.m22 * b.m12 + a.m32 * b.m13;
  erg.m22 = a.m12 * b.m21 + a.m22 * b.m22 + a.m32 * b.m23;
  erg.m32 = a.m12 * b.m31 + a.m22 * b.m32 + a.m32 * b.m33;

  erg.m13 = a.m13 * b.m11 + a.m23 * b.m12 + a.m33 * b.m13;
  erg.m23 = a.m13 * b.m21 + a.m23 * b.m22 + a.m33 * b.m23;
  erg.m33 = a.m13 * b.m31 + a.m23 * b.m32 + a.m33 * b.m33;
  return erg;
}
inline V3 operator * (const Matrix &a, const V3 &b){
  V3 erg;
  erg.x = a.m11*b.x + a.m21*b.y + a.m31*b.z;
  erg.y = a.m12*b.x + a.m22*b.y + a.m32*b.z;
  erg.z = a.m13*b.x + a.m23*b.y + a.m33*b.z;
  return erg;
}
inline V3 operator * (const V3 &a, const Matrix &b){
  V3 erg;
  erg.x = b.m11*a.x + b.m12*a.y + b.m13*a.z;
  erg.y = b.m21*a.x + b.m22*a.y + b.m23*a.z;
  erg.z = b.m31*a.x + b.m32*a.y + b.m33*a.z;
  return erg;
}
inline Matrix operator * (const Matrix &a,const double &b){
  Matrix erg;
  erg.m11 = a.m11*b;
  erg.m21 = a.m21*b;
  erg.m31 = a.m31*b;

  erg.m12 = a.m12*b;
  erg.m22 = a.m22*b;
  erg.m32 = a.m32*b;

  erg.m13 = a.m13*b;
  erg.m23 = a.m23*b;
  erg.m33 = a.m33*b;
  return erg;
}
inline Matrix operator + (const Matrix &a, const Matrix &b){
  Matrix erg;
  erg.m11 = a.m11+b.m11;
  erg.m21 = a.m21+b.m21;
  erg.m31 = a.m31+b.m31;

  erg.m12 = a.m12+b.m12;
  erg.m22 = a.m22+b.m22;
  erg.m32 = a.m32+b.m32;

  erg.m13 = a.m13+b.m13;
  erg.m23 = a.m23+b.m23;
  erg.m33 = a.m33+b.m33;
  return erg;
}
//! MyAtom an atom object
typedef struct MyAtom{
  //! Label is the atom label with residue nuber seperated by '_' and the symmetry operation number seperated by french quotes.
  QString Label;
  //! ResiClass is the four or three letters residue class
  QString ResiClass;
  //! The original line in the res file.
  QString orginalLine;
  //! fragment number of the structure part
  int molindex;
  //! Binary flag for fixed (tied to first FVAR) atom parameters.
  int fixFlag;
  //! the symmetry number  
  int symmGroup ;
  //!symmetry operation number
  int sg;
  //! platon style symetry code (s555)
  int scod;
  //! index of the source atom in the asymmetric unit
  int auidx;
  //! a is hidden flag
  int hidden;
  //! the site occupancy factor as it is the file e.g. 11.000
  double sof_org;
  //! the site occupancy factor e.g. 1.0 or 0.5
  double sof;
  //! isotropric flag
  bool isIso;
  //! the Uiso string as it appears in the res file eg -1.5 for 150% of the ueq of the horse atom 
  QString ufiso_org;
  //! the fractional Uij
  Matrix uf;
  //! the cartesian Uij
  Matrix uc;
  //! the cartesian coordinates.
  V3 pos;
  //! the fractional coordinates.
  V3 frac;
  //! the electrondensity value of a Q-Peak.
  double peakHeight;
  //! the residue number
  int resiNr;
  //! the atomic number
  int an;
  //! the (dissordered) part number
  int part;
  //! the style type of the atom
  int style;
  //! the afix type of the parent atom
  int afixParent;
  //! the afix type
  int afix;
}MyAtom;
inline bool  operator == (const MyAtom &a1,const MyAtom &a2){
  return ((a1.Label == a2.Label)&&(a1.resiNr == a2.resiNr)&&(a1.part == a2.part)&&(a1.symmGroup == a2.symmGroup));
} 
inline bool operator < (const MyAtom &a1, const MyAtom &a2){
  return (a1.Label < a2.Label);
} 
inline bool operator < (MyAtom &a1, MyAtom &a2){
  return (a1.Label < a2.Label);
} 
//! MyBond a bond object
typedef struct MyBond{
  //! a pointer to the atom where the bond starts
  MyAtom const *ato1;  
  //! a pointer to the atom where the bond ends
  MyAtom const *ato2;
  //! the bond length
  double length;
  //! the index of the atom where the bond starts  
  int a1;
  //! the index of the atom where the bond ends
  int a2;
}MyBond;
inline bool operator ==(const MyBond &a1,const MyBond &a2){
    return ((a1.a1==a2.a1)&&(a1.a2==a2.a2))||((a1.a2==a2.a1)&&(a1.a1==a2.a2));
}
//! for the BIND instruction 
typedef struct MyBind{
  QString Lab1,Lab2;
}MyBind;

struct Vert {
  int faces[3];
  V3 pos;
};

typedef struct VPoly{
 V3 mid,nor,verts0,verts1;
 int acol,intra;
}VPoly;
typedef QList<MyBond> Connection;
typedef QList<MyAtom> CEnvironment;
//! Unit cell parameters
struct Cell {
  //! the dimension a in Angstrom
  double a;
  //! the dimension b in Angstrom
  double b;
  //! the dimension c in Angstrom
  double c;
  //! the angle alpha in degrees
  double al;
  //! the angle beta in degrees
  double be;
  //! the angle gamma in degrees
  double ga;
  //! \f$\varphi =  \sqrt(1 - (\cos(\alpha)^2) - (\cos(\beta)^2) - (\cos(\gamma)^2) + 2\cos(\alpha)\cos(\beta)\cos(\gamma))\f$
  double phi;
  //! the cell volume in Angstrom^3
  double V;
  //! the reciprocal dimension a 
  double as;
  //! the reciprocal dimension b 
  double bs;
  //! the reciprocal dimension c 
  double cs;

  //! \f$ \tau = c ((\cos(\alpha) - \cos(\beta)  \cos(\gamma)) / \sin(\gamma))\f$
  double tau;
  //! \f$ \cos(\gamma)\f$ 
  double cosga;
  //! \f$ \cos(\alpha) \f$ 
  double cosal;
  //! \f$ \cos(\beta) \f$ 
  double cosbe;
  //! \f$ \sin(\alpha) \f$ 
  double sinal;
  //! \f$ \sin(\beta) \f$ 
  double sinbe;
  //! \f$ \sin(\gamma) \f$ 
  double singa;
  //! \f$ \tan(\gamma) \f$ 
  double tanga;
  double cosra;//! cos reciprocal alpha
  double cosrb;//! cos reciprocal beta
  double cosrg;//! cos reciprocal gamma
  //! List of symmetry operators
  QList<Matrix> symmops;
  //! List of translations
  QList<V3> trans;
  //!number of symmops without centring and inversion applyed
  int ns0;
  //! space group is centred (A,B,C,I,F)
  bool centered;
  //! space group is centro symmetric
  bool centeric;
  //! ZERR string as it is from the res file
  QString Z;
  //! the wavelenth  \f$ \lambda  \f$ in Angstroms
  double wave;
  Matrix G,Gi;//metric tensor and its inverse
};
//! shortest distance matrix item type
struct SdmItem{
  //! shortest contact distance
  double d;
  //! atom index of the starting atom
  int a1;
  //! atom index of the ending atom
  int a2;
  //! symmetry number of the contact
  int sn;
  V3 floorD;
  //! contact is covalent
  bool covalent;
};

inline bool operator < (const SdmItem &a1, const SdmItem &a2){
  return (a1.d<a2.d);
}

struct Tripel{
  int n[3];
  V3 midcenter;
};
inline bool operator == (const Tripel& v1, const Tripel& v2 ) {
  return ((v1.n[0]==v2.n[0])&&(v1.n[1]==v2.n[1])&&(v1.n[2]==v2.n[2]));
}
//!Knopf is a list of neighbors of an atom
struct Knopf{
  //! neighbors is list of atom indices of neigbours of that atom
  QList<int> neighbors;
};
class Ring;

/*! \brief Molecule is a crystallography and molecular structure object of shelXle. 
 *
 * It has routines to convert coordinates and 
 * ADPs from fractional to cartesian space, symmetry operations etc. It provides also routintes to draw atoms bonds and q-peaks. 
 * It is thought that there is only one instance of this class in ShelXle during runtime. Molecule accesses the settings file of ShelXle.
 */
class Molecule{
  public:    
    bool qbeforehkl;//!< There are Q-Peaks before HKLF instruction 
    double  qboMin,qboMax;
    double exti,swat,osf;
    int hklf;
    double hklScale,hklSigmaScale,hklOmitSig,hklOmit2th,hklShellLow,hklShellHig;
    Matrix hklMat;
    QString titl;
    QString mynewnameis;    
    Molecule();//!<The constructor
    int icocnt;
    int latt;
    int afix5(int idx);//!< retuns fo example 65 if less then 6 non H atoms (in AFIX 66) are before this line.
    void setHBondMaxDist(double d);
    void setHBondMaxAngl(double w);
    int genSymCode(int s, int h, int k, int l);//!< creates a platon style symmetry code as an integer (s555)
    V3 applySymCode(const V3 &frac, int scode);//!< applies a platon style symmetry code as an integer (s555)
    double hbdist();
    double hbangl();
    Cell cell;//!<The unit cell parameters
    QString Fragments;//!<A html string of the kind "<b>The asymmetric unit contains N fragments.</b><br>" were N is the number of fragments.
    QList<SdmItem> sdm;//!< shortest distance matrix.
    QList<SdmItem> envi_sdm;//!< shortest distance matrix for the ENVIron ment listing functionality.
    QList<SdmItem> contact;//!< shortest distance matrix with possible hydrogen bonds for auto HFIX.
    QList<Matrix> symmopsEQIV;//!< list of symmetry operations from  EQIV instructions.
    QStringList labelEQIV;//!< list of atom lable from  EQIV instructions.
    QList<V3> transEQIV;//!< list of translation from  EQIV instructions.
    QMap<QString,int> eqivMap;//!< maps labelEQIV and symmetry operation indices.
    QList<MyBind> freeatoms,//!< bond list for FREE instructions
      bindatoms;//!< bond list for BIND instructions 
    QList<Knopf> knoepfe;//!< neighbors list of each atom
    QMap<int,int> bindPart;
    bool canbeDonor(int an);
    bool canbeAcceptor(int an);
    QList<int> theseAreAcceptors;
    QList<int> theseAreDonors;
    unsigned short Kovalenz_Radien[109];
    double  pmin,//!< minimum peak height of the Q-Peaks.
            pmax;//!< maximum peak height of the Q-Peaks.
    double  lmin,//!< minimum a for belo.
            lmax;//!< maximum a for belo.
    CEnvironment asymm;//!< atom list containing the asymertic unit
    CEnvironment showatoms;//!< atom list containing all atoms to be drawn.
    Connection showbonds;//!< list of visible bonds.
    Connection lbonds;//!< list of line bonds (for Q-Peaks)
    Connection theH_Bonds;
    QStringList usedSymmetry;//!< List of internal symmetry codes of symmetry operations curently applied.
    QString HumanSymmetry;//!< A human readable list of symmetry operations in use.
    QSettings *einstellung;//!< the QSetttings for ShelXle.
    QString voroMsg;
    QList<VPoly> vtriangles;
    void voronoij(CEnvironment au, int intat=-1);
    double ueq(Matrix m);
    //unsigned short ElNeg[83];
    //
    //
    QStringList sdmcompleter();
    QList<SdmItem> computeSDM(CEnvironment au);
    Connection connecting(const CEnvironment &atom, bool eac=false);
    bool decodeSymmCard(const QString symmCard);
    bool decodeSymmCardEQIV(const QString symmCard);
    //void countMols(QList<INP> & xdinp);
    //bool applyLatticeCentro(const QChar latt,const bool centro);
    QString symmcode2human(QStringList brauchSymm);//!< converts a list of internal symmetry codes in human readable ones.
    QString symmcode2human(QString brauchSymm);//!< converts an internal symmetry code into human readable symmetry operation description
    QString symmcode2human(QString brauchSymm, int j);//!< converts an internal symmetry code into human readable symmetry operation description prepends french quotes j: 
    QString symmcode2human(int s);
    QString symmCard2Code(QString symmCard);
    void frac2kart (V3 x, V3 & y);
    void kart2frac (V3 x, V3 & y);
    int getOZ(QString S1);//!< returns the atomic number (starting from 0) for a given chmical elememt symbol S1. @param S1 element symbol.
    double shortestDistance(QString sc);
    static double winkel(V3 a,V3 b);//!< calculates the angle in degrees between two vectors 
    double dieder(V3 a,V3 b, V3 c);//!< calculates a torsion angle from 3 given vectors
    static V3 kreuzX(double x1,double y1,double z1,double x2,double y2,double z2) ;//!< a cross product 
    double fl(double x,double y, double z);//!< calculates the length of a vector given by x,y,z in fractional coordinates in Angstroms @param x,y,z fractional coordinates of a vector. 
    double wl(V3 f2, V3 f1, V3 f3);//!< calculates angle in degrees of the 3 positional vectors in fractional space. 
    void cellSetup();//!< calculates unit cell parameters from a,b,c, alpha, beta, gammma
    void fuse();
    void grow();
    void grow_plus();
    void fillCell();
    void uniqueInCell();
    void packer(QStringList brauchSymm);
    void packInLimits(double ami, double ama, double bmi, double bma, double cmi, double cma);//!< Pack inside given limits (eg multiple unit cells)
    CEnvironment packInLimits(CEnvironment &au, double ami, double ama, double bmi, double bma, double cmi, double cma);//!< Pack inside given limits (eg multiple unit cells)
    void applyLatticeCentro(int gitter);
    void Uf2Uo(const Matrix x, Matrix & y) ;
    void Usym (Matrix x,Matrix sym, Matrix & y);
    void USymCode(Matrix x,int scod, Matrix & y);
    double dimension();
    double dimension(CEnvironment ce);
    void expandAt(int index);
    void expandAll();
    void complete();
    void enviSDM(double range);
    QString pse(int oz);//!< returns the chemical symbol for the given atomic number (starting from 0)
    int maxmols(){return maxmol;}//!< returns the number of fragments (molecules)
    void newstructure(){maxmol=0;}//!< new structure initialize maxmol to 0
    const QStringList thepse(){return PSE;}
    int pseSize(){return PSE.size();}
    Matrix u2b(Matrix u);//!< compute Bij values from Uij 
    Matrix b2u(Matrix u);//!< compute Uij values from Bij 
    void extendChain(QList<Ring > &rings, QList<int> &currentChain, QList<QList<int> > &neighbors,QSet<int> &nogo);

    void inventNewLabels(QList<int> &result);//!< lets see where it goes
    double * jacobi(const Matrix &uij, V3 &ev);
  private:
    void neighborSort(Knopf &kn);
    double HAMax,HAWink;
    int maxmol;
    QStringList PSE;
};

class Ring{
  private:
    bool complete;

  public:
    bool done;
    QList<Tripel> t;
    QList<int> members;
    V3 center;
    Ring(Tripel t1,Tripel t2){
      complete=false;
      done=false;
      t.append(t1);
      t.append(t2);
      center = t1.midcenter + t2.midcenter;
      members.append(t1.n[1]);
      members.append(t1.n[0]);
      members.append(t1.n[2]);
      if (!members.contains(t2.n[1])) members.append(t2.n[1]);
      if (!members.contains(t2.n[0])) members.append(t2.n[0]);
      if (!members.contains(t2.n[2])) members.append(t2.n[2]);

    }
    ~Ring(){
      t.clear();
      members.clear();
    }

    bool isComplete(){
      if (complete) return complete;
      complete = (members.size()==t.size());
      if (complete) bringInOrder();
      //printf("CPLTE%d %d %d\n",complete,members.size(),t.size());
      return complete;
    }

    V3 theCenter(){
      return center*(1.0/t.size());
    }

    void bringInOrder(){
      //printf("bio\n");
      //for (int aai=0; aai<members.size(); aai++)printf("Members %d\n",members.at(aai));
      members.clear();
      members.append(t.at(0).n[0]);
      int loop=0;
      while (members.size()!=t.size()){
        loop++;
        for (int i=0; i<t.size(); i++){
          if ((members.last()==t.at(i).n[1])&&(!members.contains(t.at(i).n[0]))&&((members.size()+1==t.size())||(members.first()!=t.at(i).n[2])))
            members.append(t.at(i).n[0]); else
              if ((members.last()==t.at(i).n[2])&&(!members.contains(t.at(i).n[0]))&&((members.size()+1==t.size())||(members.first()!=t.at(i).n[1])))
                members.append(t.at(i).n[0]);
        }
        if (loop>50){
          /*printf("ALARM? %d %d %d\n",members.size(),t.size(),loop);
            for (int aai=0; aai<members.size(); aai++)printf("members %d\n",members.at(aai));
            for (int tai=0; tai<t.size(); tai++)printf("t %d %d %d\n",t.at(tai).n[0],t.at(tai).n[1],t.at(tai).n[2]);*/
          complete=false;
          done=true;
          members.clear();
          return;
        }
      }  
      isItARealRing();
    }

    void isItARealRing(){
      //printf("complete %d length %d %d\n",complete,t.size(),members.size());
      if (complete==false) return;
      bool killme=false;
      killme = (t.size()!=members.size());
      for (int i=0; i<t.size(); i++){
        int u,v,w;
        u=members.indexOf(t.at(i).n[0]);
        v=members.indexOf(t.at(i).n[1]);
        w=members.indexOf(t.at(i).n[2]);
        killme|=!((((u+1)%t.size()==v)&&((u-1+t.size())%t.size()==w))||(((u+1)%t.size()==w)&&((u-1+t.size())%t.size()==v)));
        killme|=(u*v*w<0);
        /*printf("%d?%d %d %d %d ?%d\n",i,u,v,w,
          (((u+1)%t.size()==v)&&((u-1+t.size())%t.size()==w))||
          (((u+1)%t.size()==w)&&((u-1+t.size())%t.size()==v)),
          u*v*w<0);// */

      }
      if (killme) {
        //printf("NOT A REAL RING!\n");
        complete=false;
        members.clear();
        done=true;
      } 
    }

    void printRing(Molecule *mol){
      QString str="+";
      for (int i = 0; i < members.size(); i++){
        str.append(mol->asymm.at(members.at(i)).Label);
        str.append((i==members.size()-1)?"+\n+":"=");
      }
      int l=str.size()-4;
      for (int i=0; i<l; i++) str.append("=");
      str.append("+\n");
      printf("%s",str.toStdString().c_str());
      str.clear();
    }

    void debugRing(Molecule *mol){
      printRing(mol);
      QString str="+";
      for (int i = 0; i < t.size(); i++){
        str.append("|");
        str.append(mol->asymm.at(t.at(i).n[1]).Label);
        str.append("~");
        str.append(mol->asymm.at(t.at(i).n[0]).Label);
        str.append("~");
        str.append(mol->asymm.at(t.at(i).n[2]).Label);
        str.append((i==t.size()-1)?"+\n+":"=");
      }
      int l=str.size()-4;
      for (int i=0; i<l; i++) str.append("=");
      str.append("+\n");
      printf("%s",str.toStdString().c_str());
      str.clear();

    }

    void addTripel(Tripel nt){
      if (t.contains(nt)) return;
      t.append(nt);
      center += nt.midcenter;
      if (!members.contains(nt.n[1])) members.append(nt.n[1]);
      if (!members.contains(nt.n[0])) members.append(nt.n[0]);
      if (!members.contains(nt.n[2])) members.append(nt.n[2]);
      isComplete();
    }

};
#endif
