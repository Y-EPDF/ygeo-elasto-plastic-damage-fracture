/* Copyright (C) 2000, Dr. Antonio Munjiza
 *
 * This code is provided as part of the book entitled "The Combined
 * Finite Discrete Element Method". It is distributed WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. Inclusion of a part or whole of this code into any other
 * commercial or research or other purpose code is not granted without author's
 * written explicit permission.
 * When results using whole or any part of this code
 * are published, Y code must be mentioned and acknowledgement to Dr Munjiza must be made.
 * Should you modify this source code, the Copyright (C) on the modified code
 * as a whole belongs to Dr. A. Munjiza regardless of the extent or nature
 * of modifications.
 * Copyright (C) to whole of any code containing any part of this code
 * also belongs to Dr. A.Munjiza.
 * Any code comprising any part of this source code
 * must be called Y program.
 * If you do not agree with this, you are not allowed to do
 * any modifications to any part of this source code or included
 * any part of it in any other program.
 */
/* file Yd.h  Y data base description */
#include "Ytypes.h"
#ifndef FRAMEINCL
#include "frame.h"
#define FRAMEINCL
#endif
#ifndef YTYPESINCL
#include "Ytypes.h"
#define YTYPESINCL
#endif

typedef struct YDC_struct *YDC;
struct YDC_struct
{ YINT  mcstep, ncstep;       /* maximum/current number of time steps                  */
  FILE *finp, *fcheck;

  DBL  dcgray;               /* gravity y                                              */
  DBL  dcsizc;               /* size coord.                                            */
  DBL  dcsizf;               /* size force                                             */
  DBL  dcsizs;               /* size stress                                            */
  DBL  dcsizv;               /* size velocity                                          */
  DBL  dcstec;               /* current time step size                                 */
  DBL  dctime;               /* current time                                           */
  YINT  icoutf;               /* write output frequency                                 */
  YINT  icoutrf;              /* write output reduced frequency                         */
  YINT  icoutnf;              /* threshold number of fractures beyond which use icoutrf */
  YINT  icouti;               /* current write output No                                */
  YINT  icoutp;               /* output precision - digits per number                   */
  YINT  icresf;               /* output frequency of restart data                       */
  YINT iuptrimesh;
  YINT iuseeczm;


};




typedef struct YDE_struct *YDE;
struct YDE_struct
{ YINT melem, nelem;          /* maximum (actual) number of elements                   */
  YINT melst, nelst;          /* maximum (actual) number of elemen. states var.        */
  YINT melno, nelno;          /* maximum (actual) number of elemen. nodes              */
  YINT nelemst;               /* actual number of elements @ the beginning             */
  YINT nebrk;                 /* actual number of joint elements broken                */
  YINT netbrk;                /* number of joint elements broken in one time step      */
  YINT nesft;                 /* number of joint elements softened                     */
  YINT netsft;                /* number of joint elements softened in one time step    */

  YINT   *i1elcf;             /*[melem]    contacting couple first                     */
  YINT   *i1elpr;             /*[melem]    element property                            */
  YINT   *i1elprtmp;          /*[melem]    element property  used for meshing          */
  DBL   *d1elfs;             /*[melem]    shear strength at joint (Mohr-Coulomb)      */
  DBL   **d2elst;            /*[melst][melem]    - element state                      */
  YINT   **i2elto;            /*[melno][melem]    - element topology                   */
  YINT   **i2eltost;          /*[melno][melem]    - element topology @ the beginning   */
  YINT   *i1ebrk;             /*[melem]         broken joint elements                  */
  DBL   **d2ecbrk;           /*[mnodim][melem] joint element breakage coordinates     */
  DBL   **d2ecbrk_NEW;           /*[mnodim][melem] joint element breakage coordinates     */
  DBL   *d1etbrk;            /*[melem]         joint element breakage time            */
  DBL   *d1elbrk;            /*[melem]         length of broken joint elements        */
  DBL   *d1efe;              /*[melem]         joint element fracture energy          */
  YINT   *i1esft;             /*[melem]         softened joint elements                */
  DBL   **d2ecsft;           /*[mnodim][melem] joint element yielding coordinates     */
  DBL   *d1etsft;            /*[melem]         joint element yielding time            */
  YINT   *i1esftf;            /*[melem]         softened element flag                  */
  DBL   *d1ebrkf;            /*[melem]         broken element flag                    */
  DBL   *d1eike;             /*[melem]    joint element initial kinetic energy        */
  DBL   *d1edke;             /*[melem]    joint element differential kinetic energy   */
  DBL   *d1etmke;            /*[melem]    joint element maximum kinetic energy time   */
  DBL   *d1elfr;             /*[melem] friction angle of elements next to DFN cracks  */
  DBL   *d1elpe;             /*[melem] normal penalty of elements next to DFN cracks  */
  DBL   *d1elpt;             /*[melem] tangent penalty of elements next to DFN cracks */
  YINT   **i2elnext;          /*[4][melem] element & edge indeces next to joints       */
  YINT   **i2eledge;          /*[nelno][melem] flags for intact/broken element edges   */
  YINT   *i1edfnf;            /*[melem] flag for joint element belonging to a DFN      */
  YINT   *i1edft;             /*[melem] flag for triangles next to DFN crack type 3    */
  DBL   *d1etike;            /*[melem] initial time of KE monitoring window           */
  
  DBL   **d2eldmg;           /*[3][melem] maximum damage coefficient of joint elements*/
  
  DBL   **d2elstr;           /*[4][melem] stress tensor of elements                   */



  /* --- 新增：无旋转构型转换 --- */

DBL **d2elR;        // [9][melem]
DBL **d2elV;        // [9][melem]
DBL **d2elAlpha;    // [9][melem]
DBL **d2elState;    // [nstate][melem]







  /* --- 新增：乘法分解状态量 --- */
  DBL **Fe00, **Fe01, **Fe10, **Fe11, **Fe33;         /* 塑性形变梯度 Fp */
  DBL **Fp00, **Fp01, **Fp10, **Fp11, **Fp33, **Fp31, **Fp32;         /* 塑性形变梯度 Fp */
  DBL **Fpkd00, **Fpkd01, **Fpkd10, **Fpkd11, **Fpkd33; /* 动硬化耗散分量 */
  DBL **Fpcd00, **Fpcd01, **Fpcd10, **Fpcd11, **Fpcd33; /* 弹性核耗散分量 */
  DBL **ss00, **ss01, **ss10, **ss11, **ss33 , **ss31, **ss32;  
  DBL **sigma00, **sigma10, **sigma11, **sigma33,**sigma31,**sigma32 ;                  /* 第二 Piola-Kirchhoff 应力张量 S */
  //DBL **Fpks00, **Fpks01, **Fpks10, **Fpks11; /* 动硬化存储分量 */
  //DBL **Fpcs00, **Fpcs01, **Fpcs10, **Fpcs11; /* 弹性核存储分量 */

  /* --- 新增：Mandel-like 张量（对称存三个分量） --- */
  DBL **M00, **M01,**M10, **M11, **M33, **M31, **M32;
  DBL **Mk00, **Mk01, **Mk11,**Mk33;
  DBL **Mc00, **Mc01, **Mc11,**Mc33;

  /* --- 新增：硬化标量 --- */
  DBL **H;    /* 等向硬化 */
  DBL **R;    /* 次载面比 */
  DBL **Rc;   /* 弹性核比 */
  DBL **EPS; /* 塑性变形 */
  /* --- 可选诊断 --- */
  DBL **detFe, **detFp, **eqp;


  /* --- 新增：ECZM --- */

  DBL   **d2centroid;    // 每个原始三角形的质心坐标 [2][n_elem]
  DBL   *d1area;            /*[melem] area          */

  DBL   *d1inst;          /*[melem]    节理单元激活计算时使用的拉应力       */
  DBL   *d1inss;          /*[melem]    节理单元激活计算时使用的剪应力         */
  DBL   **d2dmg;          /*[melem]    节理单元激活计算时使用的拉应力       */
DBL   **d2D;   
  YINT   *i1elfr;          /*[melem]    element property  used for meshing          */


 DBL   *d1J; 
  DBL   *d1J2; 
    DBL   *d1volce; 
 YINT** i2elel0;          /*[4][melem]单元边对应节理       */


};

typedef struct YDI_struct *YDI;
struct YDI_struct
{ YINT micoup, nicoup;        /* maximum possible number of contacting couples         */
  YINT    iiecff;             /* interaction element contact. couple free first        */

  DBL    diedi;              /* travel since last detection                           */
  DBL    diezon;             /* buffer zone size                                      */
  DBL   *d1iesl;             /*[mcoup] contact sliding                                */
  YINT   *i1iecn;             /*[mcoup] couple next                                    */
  YINT   *i1iect;             /*[mcoup] couple target                                  */
  YINT   mistate;             /* number of states for d2sldis                          */
  DBL   **d2sldis;           /*[mistate][mcoup] sliding distance btw  couples         */
};

typedef struct YDN_struct *YDN;
struct YDN_struct
{ YINT mnodim, nnodim;        /* max(actual) nodal dimensions number                   */
  YINT mnopo, nnopo;          /* maximum (actual) number of nodal points               */
  YINT nnopst;                /* actual number of nodal points at the beginning        */

  DBL   *d1nmct;             /* [mnopo] nodal mass current translation                */
  DBL  **d2ncc;              /* [mnodim][mnopo] nodal coordinate current              */
  DBL  **d2nci;              /* [mnodim][mnopo] nodal coordinate initial              */
  DBL  **d2nfc;              /* [mnodim][mnopo] nodal force current                   */
  DBL  **d2nvc;              /* [mnodim][mnopo] nodal velocity current                */
  YINT  *i1nobf;              /* [mnopo] nodal boundary >0 is boundary                 */
  YINT  *i1nopr;              /* [mnopo] nodal boundary condition                      */
  YINT  *i1nowe;              /* [mnopo] hydrofrac boundary, >0 node is wet            */
  YINT  **i2noid;             /* [2][mnopo] hydrofrac, IDs of cw and ccw nodes         */
  DBL   *d1nfp;              /* [mnopo] hydrofrac, nodal fluid pressure               */
  DBL  **d2nc0;              /* [mnodim][mnopo] nodal coordinate at time step 0       */


  /* --- 新增：ECZM --- */
  YINT  *i1nobf0;              /* [mnopo]  >0哦判断是不是边界节点   */

  YINT *i1getmaster;                         /* [mnopo]  主节点   */
  YINT *i1getmasterfem;                      /* [mnopo]  节点对应有限元节点   */
  YINT  *i1notoel;                      /* [mnopo]  节点对应三角   */
  YINT  *i1remeshf;                      /* [mnopo]  节理破坏后节点标志   */



};

typedef struct YDMN_struct *YDMN;
struct YDMN_struct
{          
  YINT mnopo, nnopo;          /* maximum (actual) number of nodal points               */
  DBL  **d2ncc;              /* [mnodim][mnopo] nodal coordinate current              */
  DBL  **d2nci; 
  YINT **i2elsort;    /* [主节点][50] 每个主节点连接的所有三角/排序后结果             */
  YINT *i1elcount;          /* [主节点]主节点连接三角数量             */
  YINT **i2elno;  /* [主节点][50] 每个主节点连接的所有三角对应的从节点            */
  YINT **i2elnext;     // 循环链表的下一个索引 (nnopo × 50)
  YINT *i1head;        // 链表的起始索引 (nnopo)
  YINT *i1boundaryflag;   // 记录每个主节点是否为边界节点（0-否，1-是，2-两个端点）

  //体积锁定
  DBL  *d1noa0; //主节点的初始体积
  DBL  *d1noac; //主节点的当前体积

};







typedef struct YDB_struct *YDB;
struct YDB_struct
{ YINT mborh, nborh;          /* maximum (actual) number of boreholes                  */
  YINT mbdim, nbdim;          /* maximum (actual) number of dimensions of boreholes    */
  YINT nbpaf;                 /* number of amplitude factors                           */

  DBL   **d2bca;             /* [mbdim][mborh] coordinates of point A in borehole     */
  DBL   **d2bcb;             /* [mbdim][mborh] coordinates of point B in borehole     */
  DBL   *d1brad;             /* [mborh] radii of boreholes                            */
  DBL   *d1bpaf;             /* [nbpf] amplitude factor of pressure amplitude         */
  DBL   *d1bpts;             /* [mborh] start time of pressure load on boreholes      */
  DBL   *d1bpte;             /* [mborh] end time of pressure load on boreholes        */
  DBL   *d1bvdt;             /* [mborh] velocity of detonation on boreholes           */
  DBL   *d1bprs;             /* [mborh] amplitudes of pressure for each borehole      */
  DBL   dblmax;              /* max length of all boreholes                           */
  DBL   dbbuf;               /* buffer (max. size of element)                         */
};

typedef struct YDS_struct *YDS;
struct YDS_struct
{ YINT msour, nsour;          /* maximum (actual) number of sources                    */
  YINT msdim, nsdim;          /* maximum (actual) number of dimensions of sources      */
  YINT nspaf;                 /* number of pressure amplitude factors (p=p(t) )        */
  YINT nssaf;                 /* number of pressure amplitude factors (s=s(r) )        */

  DBL   **d2scs;             /* [msdim][msour] coordinates of sources                 */
  DBL   *d1spaf;             /* [nspaf] pressure amplitude factors (p=p(t) )          */
  DBL   *d1ssaf;             /* [nssaf] pressure amplitude factors (s=s(r) )          */
  DBL   *d1spts;             /* [msour] start time of pressure load                   */
  DBL   *d1spte;             /* [msour] end time of pressure load                     */
  DBL   *d1svpr;             /* [msour] velocity of pressure propagation              */
  DBL   *d1sprs;             /* [msour] amplitudes of pressures (p=p(t) )             */
  DBL   *d1ssir;             /* [msour] source initial radius                         */
  DBL   dsbuf;               /* buffer (max. size of element)                         */
};

typedef struct YDO_struct *YDO;
struct YDO_struct
{ YINT mohys, nohys;          /* maximum (actual) number of history variables          */

  DBL     dohyp;             /* output history accuracy                               */

  DBL   *d1ohyf;             /*[mohys] output history factor to scale state           */
  DBL   *d1ohyc;             /*[mohys] output history factor to scale time            */
  DBL   *d1ohys;             /*[mohys] output history state                           */
  DBL   *d1ohyt;             /*[mohys] output history time                            */
  DBL   *d1ohyx;             /*[mohys] output history x coordinate of the point       */
  DBL   *d1ohyy;             /*[mohys] output history y coordinate of the point       */
  DBL   *d1ohyz;             /*[mohys] output history z coordinate of the point       */
  FILE  **f2ohyf;            /*[mohys] output history files                           */
  YINT   *i1ohyt;             /*[mohys] output history type, i.e. which variable       */
  FILE  **f2orsf;            /*[...] restart data files                               */
  FILE  *foebrk;             /* output time, position and energy of broken elements   */
  FILE  *foesft;             /* output time, position and energy of softened elements */
  FILE  *fohyfr;             /* time, flow rate and fluid pressure for hydrofrac      */
  FILE  *fofrac;             /* total fracture area in the model                      */
};

/* Y Database Properties - Elements */
typedef struct YDPE_struct *YDPE;
struct YDPE_struct
{ YINT mprop, nprop;          /* maximum (actual) number of properties                 */

  DBL   *d1peca;             /*[mprop] child age - procreation                        */
  DBL   *d1pecl;             /*[mprop] child life - interval for procreation          */
  DBL   *d1peks;             /*[mprop] dpeks=2hbeta*sqrt(E*ro) in 2D or 3D,0<beta<1   */
  DBL   *d1pela;             /*[mprop] property lamda - Lame elastic constant         */
  DBL   *d1pemu;             /*[mprop] property mu    - Lame elastic constant         */
  DBL   *d1pepe;             /*[mprop] contact penalty parameter for Y-RC             */
  //DBL   *d1pept;             /*[mprop] tangential penalty param., (0.1*ydpj->d1pjpe)  */
  //DBL   *d1pefr;             /*[mprop] Coloumb friction                               */
  DBL   *d1pera;             /*[mprop] property radius of sphere                      */
  DBL   *d1pero;             /*[mprop] property ro    - density                       */
  DBL   *d1pevi;             /*[mprop] viscosity for  granular flow                   */
  
  YINT    mperow;             /* Possible combinations of Pr sets for phi, p_n, p_t    */
  DBL   **d2peint;           /*[mprop][5] Coulomb friction, normal penalty, tangential penalty between property sets */
  
  DBL   *d1peem;             /*[mprop] Young's modulus                                */
  DBL   *d1penu;             /*[mprop] Poisson's ratio                                */

  DBL   *d1psem;             /*[mprop] maximum tensile stretch                        */

  YINT   *i1pecn;             /*[mprop] property No to be assigned to child's nodes    */
  YINT   *i1pecp;             /*[mprop] permanent property to be assigned to child     */
  YINT   *i1pect;             /*[mprop] temporary property to be assigned to child     */
  YINT   *i1pefr;             /*[mprop] if >, fracture                                 */
  YINT   *i1pejp;             /*[mprop] joint property; if<0, no joints                */
  YINT   *i1pemb;             /*[mprop] mark boundary nodes 1 yes 0 no                 */
  YINT   *i1pemn;             /*[mprop] number of mesh refinements                     */

  YINT   *i1pnib;             /*[mprop] 1 if borhole regardles boundary, else 0        */
  YINT   *i1ptyp;             /*[mprop] property type                                  */ 
  YINT   *i1psde;             /*[mprop] state damage elastic id                        */
  
  YINT   *i1pexc;             /*[mprop] excavation flag                                */
  
  YINT   *i1usan;             /*[mprop] flag for using anisotropic elasticity (=1:use) */
  DBL   *d1peex;             /*[mprop] Young's modulus E_x                            */
  DBL   *d1peey;             /*[mprop] Young's modulus E_y                            */
  DBL   *d1pemx;             /*[mprop] Poisson's ratio mu_xy                          */
  DBL   *d1pemy;             /*[mprop] Poisson's ratio mu_yx                          */
  DBL   *d1peg;              /*[mprop] Shear modulus G                                */
  
  YINT   *i1psup;             /*[mprop] excavation support flag                        */ 

/* ---------------- 显式塑性 + 扩展次载面模型材料参数（DP版） ---------------- */



/* Drucker–Prager 屈服参数（有限变形） */
DBL   *d1phi;      /* 摩擦角 φ（rad 或 deg，按实现约定） */
DBL   *d1DPc;        /* 内聚力 c（同单位制） */
DBL   *d1sig_t;    /* 抗拉截断 σ_t（可选，<=0 表示不启用） */

/* 等向硬化参数（F(H)） */
DBL   *d1F0;       /* 初始屈服函数值（等效初始强度） */
DBL   *d1h1;       /* 非线性硬化参数 h1（控制饱和幅度） */
DBL   *d1h2;       /* 非线性硬化参数 h2（控制硬化速率） */
DBL   *d1Hc;       /* 线性硬化模量（可选，0 表示不用） */

/* 动硬化能量参数（可选，Mk） */
DBL   *d1Ck;       /* 动硬化能量系数（ψk） */

/* 弹性核能量参数（Mc） */
DBL   *d1Cc;       /* 弹性核能量系数（ψc） */

/* 动硬化演化与极限（可选） */
DBL   *d1bk;       /* 动硬化饱和值 */
//DBL   *d1ck;       /* 动硬化速率参数 */
//DBL   *d1nk;       /* 动硬化自旋系数 */

/* 次载面参数（R, Rc 演化） */
DBL   *d1URu;        /* U(R) 放大系数（式(68)） */
DBL   *d1uc;       /* 重装载放大系数（式(111)） */
DBL   *d1Re;       /* 准弹性阈值 R_e */
DBL   *d1Chi;      /* 弹性核极限比 χ */

/* 塑性自旋参数（式(113)） */
DBL   *d1nn;       /* 塑性自旋系数 n */
DBL   *d1nP;       /* 塑性自旋系数 n^p */
DBL   *d1n2;       /* 弹性核自旋系数 n^c */

/* 算法选项 */
YINT  *i1use_ext_subloading; /* 1=使用扩展次载面（Mc），0=初始次载面 */


YINT  *i1use_kinhard;       /* 1=启用 Mk/ψk 与 Fpkd 更新；0=关闭，直接跳过 Mk、Fpkd 的更新 */





};

/* Y Database Properties for Joints */
typedef struct YDPJ_struct *YDPJ;
struct YDPJ_struct
{ YINT mpjset, npjset;        /* maximum (actual) number of Joint Property sets       */

  DBL *d1pjfs;               /* [mpjset] ultimate shear strength at joint            */
  DBL *d1pjft;               /* [mpjset] ultimate tensile strength at joint          */
  DBL *d1pjgf;               /* [mpjset] ultimate fracture energy at joint in tension*/
  DBL *d1pjgs;               /* [mpjset] ultimate fracture energy at joint in shear  */
  DBL *d1pjco;               /* [mpjset] cohesion at joint                           */
  DBL *d1pjfr;               /* [mpjset] friction at joint                           */

  DBL *d1pjpe;               /* [mpjset] fracture penalty parameter                  */

  YINT *i1psde;               /* [mpjset] state damage elastic id                     */
  YINT *i1ptyp;               /* [mpjset] Property Joint Type                         */
  
  DBL *d1usaf;               /* [mpjset] flag for using anisotropic fracture (=1or2) */
  DBL *d1pjcr;               /* [mpjset] reduced cohesion at joint                   */
  DBL *d1pjfd;               /* [mpjset] reduced friction at joint                   */
  DBL *d1pjtr;               /* [mpjset] reduced ultimate tensile strength at joint  */
  DBL *d1pjgr;               /* [mpjset] reduced ultimate fracture energy at joint in tension */
  DBL *d1pjsr;               /* [mpjset] reduced ultimate fracture energy at joint in shear */
  DBL *d1pjal;               /* [mpjset] layering orientation (0-180)                */
 
  YINT iusehy;                /* flag for using hysteretic fracture model             */
};

/* Y Database Properties - Nodes, Nodal Properties or Boundary Conditions for nodes  */
typedef struct YDPN_struct *YDPN;
struct YDPN_struct
{ YINT mpnset, npnset;        /* maximum (actual) number of node property sets        */
  YINT mpnfact, npnfact;      /* maximum (actual) number of factors                   */
  DBL ***d3pnfac;            /*[2][mbc][mfact] time and amplitude factor             */

  YINT   *i1pnfx;             /*[mpnset] fixity x direction 1 force; 2 acc. 3 vel.    */
  YINT   *i1pnfy;             /*[mpnset] fixity y direction 1 force; 2 acc. 3 vel.    */
  YINT   *i1pnfz;             /*[mpnset] fixity z direction 1 force; 2 acc. 3 vel.    */

  DBL   *d1pnaf;             /*[mpnset] amplitude factor all ampltd multp by it      */
  DBL   *d1pnap;             /*[mpnset] amplitude of element surface pressure        */
  DBL   *d1pnat;             /*[mpnset] amplitude of element surface traction        */
  DBL   *d1pnax;             /*[mpnset] amplitude of force/velocity x                */
  DBL   *d1pnay;             /*[mpnset] amplitude of force/velocity y                */
  DBL   *d1pnaz;             /*[mpnset] amplitude of force/velocity z                */

  DBL   *d1pnxx;             /*[mpnset] direction of local x                         */
  DBL   *d1pnxy;             /*[mpnset] direction of local x                         */
  DBL   *d1pnxz;             /*[mpnset] direction of local x                         */
  DBL   *d1pnyx;             /*[mpnset] direction of local y                         */
  DBL   *d1pnyy;             /*[mpnset] direction of local y                         */
  DBL   *d1pnyz;             /*[mpnset] direction of local y                         */
  DBL   *d1pnzx;             /*[mpnset] direction of local z                         */
  DBL   *d1pnzy;             /*[mpnset] direction of local z                         */
  DBL   *d1pnzz;             /*[mpnset] direction of local z                         */
};

/* Y Database Properties - Meshing */
typedef struct YDPM_struct *YDPM;
struct YDPM_struct
{ YINT mpmcom;                /* combination of pr. sets (rows in I2PMSET)            */
  YINT mpmcol;                /* Number of columns  in I2PMSET                        */
  YINT **i2pmset;             /*[mpmcol][mpmcom] mesh pr.sets i & j together          */
  YINT mpmrow;                /* Number of rows in I2PMIJ                             */
  YINT **i2pmij;              /* Meshing combinations defining joints btw 2 pr. sets  */
};

/* Y Database Properties - Meshing */
typedef struct YDFN_struct *YDFN;
struct YDFN_struct
{ YINT iusefn;                /* flag for using discrete fracture network (=1: use)        */
  YINT mdfnfr;                /* maximum number of fractures in the DFN                    */
  YINT mdfnno;                /* maximum number of nodes per fracture                      */
  YINT **i2dfnn;              /* [mdfnfr][mdfnno] node IDs of each fracture                */
  DBL *d1dffr;               /* [mdfnfr] friction coefficient of each fracture            */
  DBL *d1dfpe;               /* [mdfnfr] normal penalty of each fracture                  */
  DBL *d1dfpt;               /* [mdfnfr] tangential penalty of each fracture              */
  DBL ddfnft;                /* Tensile strength of DFN fractures (when iusefn==2)        */
  DBL ddfnco;                /* Cohesion of DFN fractures (when iusefn==2)                */ 
  DBL ddfngf;                /* Mode I fracture energy of DFN fractures (when iusefn==2)  */ 
  DBL ddfngs;                /* Mode II fracture energy of DFN fractures (when iusefn==2) */
  YINT *i1dfft;               /* [mdfnfr] type of DFN fracture (when iusefn==3, 1=broken, 2=cohesive */
};

typedef struct YDIS_struct *YDIS;
struct YDIS_struct
{ 
  YINT  iuseis;               /* flag for using in-situ stress (=1: use)               */
  DBL  dcstxx;               /* in-situ stress tensor xx component                    */
  DBL  dcstxy;               /* in-situ stress tensor xy component                    */ 
  DBL  dcstyy;               /* in-situ stress tensor yy component                    */ 
  DBL  dcsyxx;               /* in-situ stress tensor xx component y gradient         */	
  DBL  dcsyxy;               /* in-situ stress tensor xy component y gradient         */	
  DBL  dcsyyy;               /* in-situ stress tensor yy component y gradient         */
  DBL  dcsrfy;               /* y coordinate of the (flat) topographic surface        */
};

typedef struct YDHF_struct *YDHF;
struct YDHF_struct
{ 
  YINT  iusehf;               /* flag for using in-situ stress (=1: use)               */
  YINT  ihftyp;               /* hydro-frac input type, 1 = pressure, 2 = flow rate    */
  DBL  dhfflp;               /* input fluid pressure                                  */
  DBL  dhfflq;               /* input flow rate                                       */
  YINT  hfarow;               /* number of amplitude factors                           */
  DBL **d2hfaf;              /* pressure or flow rate amplitude factor vs time        */
  DBL  fluvol;               /* Fluid volume                                          */
  DBL  flupres;              /* Fluid pressure                                        */
  DBL  flumass;              /* Fluid mass                                            */
  DBL  flurho0;              /* Fluid density at reference pressure p_0               */
  DBL  flupres0;             /* Reference fluid pressure p_0                          */
  DBL  flubulk;              /* Bulk modulus                                          */
  YINT  fradim;               /* Dimensionality of fractures (2 = 2D, 3 = 3D)          */
  YINT  ihfmsin;              /* Flag for fluid mass initialization                    */
  DBL  gravacc;              /* Gravitational acceleration for hydrostatic pressure   */
  DBL **d2wtlev;             /* Upstream and downstream x-coord and water level       */
};

typedef struct YDSM_struct *YDSM;
struct YDSM_struct
{ 
  YINT  iusesm;               /* flag for using alternative monitoring (=1: use)       */
  DBL  dctwle;               /* maximum duration of monitoring window                 */
};

/* reference points */
typedef struct YDR_struct *YDR;
struct YDR_struct
{ YINT mnodim, nnodim; /* max(actual) nodal dimensions number                   */
  YINT mrdim, nrdim;   /* max(actual) reference point dimensions number         */
  YINT mrldm, nrldm;   /* max(actual) local reference point dim. number         */
  YINT mrepo, nrepo;   /* maximum (actual) number of ref. nodal points          */
  YINT nbrjointrb;
  DBL  **d2rcig;  /* [mnodim][mrepo] coordinate initial global                 */
  DBL  **d2rccg;  /* [mnodim][mrepo] coordinate current global                 */
  DBL  **d2rccl;  /* [mrldm][mrepo]  coordinate current local                  */
  DBL  **d2rvcg;  /* [mnodim][mrepo] velocity current global                   */
  DBL  **d2riLc;  /* [mrldm][mrepo] referent coordinate initial local          */ 
  DBL   **d2rsctr;/* [2][mrepo]carent vector ref.points                        */ 
  YINT  **i2relto; /* [2][mrepo]  1D joint element topology                     */
  YINT   *i1rmyel; /* [mrepo]         my element                                */
  YINT   *i1rrpn;  /* [mrepo]  ref. point next on sing list                     */
  YINT   *i1rprop;  /*[mrepo]  properties of ref. points                        */
  YINT   *i1refbar; /*[mrepo]  ref. points - bar                                */
  YINT   *i1myjoint;   /* [mrepo]  steel joint                                  */
  DBL   *d1rbsig;   /* [mrepo]  stress in rebar                                */
  DBL   *d1rbfrc;   /* [mrepo]  force in rebar                                 */
  DBL   *d1rbstr;   /* [mrepo]  strain in rebar                                */
  DBL   *d1rjsig;   /* [mrepo]  normal (axial) stress in 1D joint              */
  DBL   *d1rjtau;   /* [mrepo]  tangential (shear) stress in 1D joint          */
  DBL   *d1rjslpnor;/* [mrepo]  normal (axial) strain in 1D joint              */
  DBL   *d1rjdelnor;/* [mrepo]  tangential (shear) strain in 1D joint          */
  YINT   *i1rjstnor; /* [mrepo]  axial state of 1D joint  (1 = ela, 2 = yielded */
  DBL   **d2rjfrv;  /* [2][mrepo] force vector in 1D joint                     */
  YINT   **i2rbedn;  /* [2][mrepo] edge nodes associated with reference point   */
};

/* steel reinforcement elements*/
typedef struct YDSB_struct *YDSB;
struct YDSB_struct
{ YINT msdim, nsdim; /* max(actual) reference point dimensions number    */
  YINT msbar, nsbar; /* maximum (actual) number of bar elements          */
  YINT isfirst;      /* iffirst step for steel bar    yes/no             */   

  DBL    **d2sic;    /*[msdim][msbar] initial steel bar coordinate        */
  YINT    *i1srpf;    /*[nsbar]        reference point first               */ 
  YINT    *i1sbpr;    /*[msbar]        bar element property                */
  DBL    *d1spea;    /*[msbar]        bar area                            */
  DBL    *d1sdiam;   /*[msbar]        bar diametar                        */ 
  DBL    *d1smmdiam; /*[msbar]        bar diametar  in (mm)               */
  DBL    *d1crlcr;   /*[msbar]        interval of discrete cracks in (mm) */
  YINT    *i1sbty;    /*[msbar]        type of bar, 0 = elastic, 1 = plastic (Y-RC formulation), 5 = elastic (Andrea's formulation) */
  YINT    *i1sbac;    /*[msbar]        0 = inactivated, 1 = activated      */
};

/* steel reinforcement property*/
typedef struct YDPS_struct *YDPS;
struct YDPS_struct
{ YINT mprop, nprop; /* max(actual) number of steel property                                        */
  DBL    *d1young;    /*[mprop] property young modulus of elasticity                               */
  DBL    *d1sfc;      /*[mprop] concrete compressive strength                                      */
  DBL    *d1mpsfc;    /*[mprop] concrete compressive strength in MPa                               */
  DBL    *d1sfy;      /*[mprop] yielding strength of steel                                         */
  DBL    *d1epssh;    /*[mprop] strain - hardening point of steel                                  */
  DBL    *d1sfu;      /*[mprop] ultimate strength  of steel                                        */
  DBL    *d1epsu;     /*[mprop] ultimate strain  of steel                                          */
  DBL    *d1sfbr;     /*[mprop] break strength  of steel                                           */
  DBL    *d1epsbr;    /*[mprop] break strain  of steel                                             */
  DBL    *d1stkn;     /*[mprop] normal stiffness of 1D joint element (used if i1sbty = 5 or 6)     */
  DBL    *d1stkt;     /*[mprop] tangential stiffness of 1D joint element (used if i1sbty = 5 or 6) */
  DBL    *d1styns;    /*[mprop] yield stress in the axial direction  (used if i1sbty = 6)          */
  DBL    *d1strns;    /*[mprop] rupture strain in the axial direction  (used if i1sbty = 6)        */
  DBL    *d1stcoh;    /*[mprop] cohesion of rebar-rock interface                                   */
  DBL    *d1stfri;    /*[mprop] friction coefficient of rebar-rock interface                       */
};



typedef struct YD_struct *YD;
struct YD_struct
{ struct YDC_struct ydc;     /* control structure                                    */
  struct YDE_struct yde;     /* element description structure                        */
  struct YDI_struct ydi;     /* interaction structure                                */
  struct YDN_struct ydn;     /* node description structure                           */
  struct YDMN_struct ydmn;     /*mapping node description structure      for     eczm  hydro-fracturing              */
  struct YDB_struct ydb;     /* borehole description structure                       */
  struct YDS_struct yds;     /* inter. fluid (source) description structure          */
  struct YDO_struct ydo;     /* output description structure                         */
  struct YDPE_struct ydpe;   /* property description structure for elements          */
  struct YDPN_struct ydpn;   /* property description structure for nodes             */
  struct YDPJ_struct ydpj;   /* property description structure for joints            */
  struct YDPM_struct ydpm;   /* property description structure for meshing           */
  struct YDFN_struct ydfn;   /* discrete fracture network structure                  */
  struct YDIS_struct ydis;   /* in-situ stress parameters                            */
  struct YDHF_struct ydhf;   /* hydro-fracturing parameters                          */
  struct YDSM_struct ydsm;   /* seismic monitoring parameters                        */
  struct YDR_struct ydr;     /* reference points description structure Y-RC          */
  struct YDSB_struct ydsb;   /* reinforcement description structure Y-RC             */
  struct YDPS_struct ydps;   /* property description structure for steel Y-RC        */
};

