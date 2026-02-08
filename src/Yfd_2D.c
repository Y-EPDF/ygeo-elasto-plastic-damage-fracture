
  /* File   Yfd.c */
#include "Yproto.h"
#include <assert.h>

#include <math.h>
#include <string.h> // 用于 memset

void calc_exsub(double F0, double S_r, double CH, double Ck, double Bk, double UA, double EL, double POI,
                double* SIG,     // [in/out] 应力数组
                double* AL,      // [in/out] 背应力数组
                double* De_m,    // [in]     应变增量数组
                int*    C_m,     // [in]     控制模式数组
                double* Ep_M,    // [in/out] 塑性应变数组
                double* H_ptr,   // [in/out] 累积塑性应变 H
                double* Hdh_ptr) // [out]    本次增量 Hdh
{
    // ==========================================
    // 1. 局部变量定义 (替代原全局变量)
    // ==========================================
    int i, j, k;
    
    // --- 定义原代码中的常量 ---
    // I_M: 单位张量 (原全局变量)
    double I_M[6] = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0}; 
    // Ds_m: 应力增量目标 (原全局变量，显式算法中默认为0)
    double Ds_m[6] = {0.0}; 

    // --- 获取当前状态 H ---
    double H = *H_ptr; 

    // --- 前置计算 F 和 Fdh (原代码是在上一步末尾算的，这里改为开头算) ---
    double F, Fdh, fhn;
    // 屈服面大小 F
    F = F0 * (1.0 + S_r * (1.0 - exp(-CH * H)));
    // 硬化模量辅助参数 Fdh
    Fdh = F0 * S_r * CH * exp(-CH * H);
    // 常数
    fhn = sqrt(2.0 / 3.0);

    /* 临时量定义 (保持原样) */
    double Pivot, Factor, SumVal;
    double RHS_Vec[6];
    double Mat_K[6][6];

    double SIGM, SD_M[6], SD;
    double ALM, ALD_M[6];
    double SA[6], SAM, SAD_M[6], SAD;
    double SG, KV, EML1, EML2, EML3;
    double EL_M[6][6], IEL_M[6][6];
    double N[6], EN[6], trNEN;
    double MP1[6], MP2[6], MP3[6], MP_N[6], MP;
    double EPM_M[6][6];
    double DX_M[6], DW_M[6];
    double DWM, DWD_M[6];
    double trNDS, trNED, LLM, LM;
    double Da_M[6];
    double DP_M[6];
    
    // 新增局部变量 (原全局变量或未定义变量)
    double CLU = 0.0; 
    double DPX = 0.0;
    double Hdh = 0.0;
    double R = 0.0;
    double UR = 0.0;
    double R_eff = 0.0;

    // ==========================================
    // 2. 核心逻辑 (与原代码完全一致)
    // ==========================================

    SIGM = (SIG[0] + SIG[1] + SIG[2]) / 3.0;
    for (i = 0; i < 6; ++i) {
        SD_M[i] = SIG[i] - SIGM * I_M[i];
    }

    SD = 0.0;
    for (i = 0; i < 3; ++i) {
        SD += SD_M[i] * SD_M[i];
    }
    for (i = 3; i < 6; ++i) {
        SD += 2.0 * SD_M[i] * SD_M[i];
    }
    SD = sqrt(SD);

    /* --- 背应力计算 --- */
    ALM = (AL[0] + AL[1] + AL[2]) / 3.0;
    for (i = 0; i < 6; ++i) {
        ALD_M[i] = AL[i] - ALM * I_M[i];
    }

    for (i = 0; i < 6; ++i) {
        SA[i] = SIG[i] - AL[i];
    }
    SAM = (SA[0] + SA[1] + SA[2]) / 3.0;

    for (i = 0; i < 6; ++i) {
        SAD_M[i] = SA[i] - SAM * I_M[i];
    }

    SAD = 0.0;
    for (i = 0; i < 3; ++i) {
        SAD += SAD_M[i] * SAD_M[i];
    }
    for (i = 3; i < 6; ++i) {
        SAD += 2.0 * SAD_M[i] * SAD_M[i];
    }
    SAD = sqrt(SAD);

    if (SAD < 1.0e-12) SAD = 1.0e-12;

    R = sqrt(3.0 / 2.0) * SAD / F;
    if (R < 1.0) CLU = 0.0;

    /*  弹性矩阵 */
    SG = EL / (2.0 * (1.0 + POI));
    KV = EL * SG / (3.0 * (3.0 * SG - EL));

    EML1 = KV + (4.0 / 3.0) * SG;
    EML2 = KV - (2.0 / 3.0) * SG;
    EML3 = SG;

    for (i = 0; i < 6; ++i) {
        for (j = 0; j < 6; ++j) {
            EL_M[i][j]  = 0.0;
            IEL_M[i][j] = 0.0;
        }
    }

    EL_M[0][0] = EML1; EL_M[1][1] = EML1; EL_M[2][2] = EML1;
    EL_M[3][3] = EML3; EL_M[4][4] = EML3; EL_M[5][5] = EML3;
    EL_M[0][1] = EML2; EL_M[1][2] = EML2; EL_M[2][0] = EML2;
    EL_M[1][0] = EML2; EL_M[2][1] = EML2; EL_M[0][2] = EML2;

    {
        double IEML1 = 1.0 / (9.0 * KV) + 1.0 / (3.0 * SG);
        double IEML2 = 1.0 / (9.0 * KV) - 1.0 / (6.0 * SG);
        double IEML3 = 1.0 / SG;

        IEL_M[0][0] = IEML1; IEL_M[1][1] = IEML1; IEL_M[2][2] = IEML1;
        IEL_M[3][3] = IEML3; IEL_M[4][4] = IEML3; IEL_M[5][5] = IEML3;
        IEL_M[0][1] = IEML2; IEL_M[1][2] = IEML2; IEL_M[2][0] = IEML2;
        IEL_M[1][0] = IEML2; IEL_M[2][1] = IEML2; IEL_M[0][2] = IEML2;
    }

    /*  N */
    for (i = 0; i < 3; ++i) {
        N[i] = SAD_M[i] / SAD;
    }
    for (i = 3; i < 6; ++i) {
        N[i] = SAD_M[i] / SAD * 2.0;
    }

    /*  UR */
    R_eff = (R > 1.0e-2) ? R : 1.0e-2;
    UR = UA * (1.0 - R_eff) * exp(1.0 / R_eff - 1.0);

    /* E:n */
    for (i = 0; i < 6; ++i) {
        EN[i] = 0.0;
        for (j = 0; j < 6; ++j) {
            EN[i] += EL_M[i][j] * N[j];
        }
    }

    /* n:E:n */
    trNEN = 0.0;
    for (i = 0; i < 6; ++i) {
        trNEN += N[i] * EN[i];
    }

    /*  MP */
    for (i = 0; i < 6; ++i) {
        MP1[i] = (Fdh / F) * fhn * SAD_M[i];
        MP2[i] = 0.0;
        if (Bk > 1.0e-12) {
            MP2[i] = R * Ck * (SAD_M[i] / SAD - (1.0 / (Bk / F)) * AL[i]);
        }
        MP3[i] = (UR / R) * SIG[i];
        MP_N[i] = MP1[i] + MP2[i] + MP3[i];
    }

    MP = 0.0;
    for (i = 0; i < 6; ++i) {
        MP += N[i] * MP_N[i];
    }

    /* ----------------- 循环起点 ----------------- */
L20:

    /* [3] 构建刚度矩阵 EPM_M */
    for (i = 0; i < 6; ++i) {
        for (j = 0; j < 6; ++j) {
            EPM_M[i][j] = EL_M[i][j] - (CLU * EN[i]) * (EN[j] / (MP + trNEN));
        }
    }

    /* [4] 混合控制求解器：高斯消元求 DX_M */
    /* 初始化 */
    for (i = 0; i < 6; ++i) {
        RHS_Vec[i] = Ds_m[i];
        for (j = 0; j < 6; ++j) {
            Mat_K[i][j] = EPM_M[i][j];
        }
    }

    /* 把应变控制方向移项到 RHS */
    for (j = 0; j < 6; ++j) {
        if (C_m[j] == 0) {
            for (i = 0; i < 6; ++i) {
                RHS_Vec[i] -= Mat_K[i][j] * De_m[j];
                Mat_K[i][j] = 0.0;
            }
        }
    }

    /* 只求解未知分量 */
    for (i = 0; i < 6; ++i) {
        if (C_m[i] == 0) {
            for (j = 0; j < 6; ++j) {
                Mat_K[i][j] = 0.0;
            }
            Mat_K[i][i] = 1.0;
            RHS_Vec[i]  = 0.0;
        }
    }

    /* 高斯消元 */
    for (k = 0; k < 5; ++k) {
        for (i = k + 1; i < 6; ++i) {
            Factor = Mat_K[i][k] / Mat_K[k][k];
            for (j = k; j < 6; ++j) {
                Mat_K[i][j] -= Factor * Mat_K[k][j];
            }
            RHS_Vec[i] -= Factor * RHS_Vec[k];
        }
    }

    /* 回代 */
    for (i = 5; i >= 0; --i) {
        SumVal = 0.0;
        for (j = i + 1; j < 6; ++j) {
            SumVal += Mat_K[i][j] * DX_M[j];
        }
        DX_M[i] = (RHS_Vec[i] - SumVal) / Mat_K[i][i];
    }

    /* 组合最终应变增量 */
    for (i = 0; i < 6; ++i) {
        if (C_m[i] == 0) DX_M[i] = De_m[i];
    }

    /* 计算应力增量 */
    for (i = 0; i < 6; ++i) {
        DW_M[i] = 0.0;
        for (j = 0; j < 6; ++j) {
            DW_M[i] += EPM_M[i][j] * DX_M[j];
        }
    }

    /* ------------ Loading index 计算与卸载判定 ------------ */
    DWM = (DW_M[0] + DW_M[1] + DW_M[2]) / 3.0;
    for (i = 0; i < 6; ++i) {
        DWD_M[i] = DW_M[i] - DWM * I_M[i];
    }

    trNDS = 0.0;
    trNED = 0.0;
    for (i = 0; i < 6; ++i) {
        trNDS += N[i] * DWD_M[i];
        trNED += EN[i] * DX_M[i];
    }

    LLM = trNED / (MP + trNEN);
    LM  = trNDS / MP;

    if (CLU == 0.0) {
        LM  = 0.0;
        LLM = 0.0;
        goto L10;
    }
    if (LLM < 0.0) {
        LM  = 0.0;
        CLU = 0.0;
        LLM = 0.0;
        goto L20;
    }

L10:
    /* 塑性应变增量 */
    for (i = 0; i < 6; ++i) {
        DP_M[i] = CLU * LLM * N[i];
    }

    /* 等效塑性应变 */
    DPX = 0.0;
    for (i = 0; i < 3; ++i) {
        DPX += DP_M[i] * DP_M[i];
    }
    for (i = 3; i < 6; ++i) {
        DPX += (DP_M[i] * DP_M[i]) / 2.0;
    }
    DPX = sqrt(DPX);

    /* (已删除) 弹性应变增量 Dex_M 的计算，因为没有传 Ee_M */

    /* 更新应力和应变 */
    for (i = 0; i < 6; ++i) {
        SIG[i]  += DW_M[i];
        // E_M[i]  += DX_M[i];  // (已删除) 你没传总应变数组
        // Ee_M[i] += Dex_M[i]; // (已删除) 你没传弹性应变数组
        Ep_M[i] += DP_M[i];
    }

    /* 更新背应力 */
    for (i = 0; i < 6; ++i) {
        if (Bk > 1.0e-12) {
            Da_M[i] = Ck * F * (SAD_M[i] / SAD - (1.0 / (Bk * F)) * AL[i]) * DPX;
            AL[i]  += Da_M[i];
        }
    }

    /* 计算 Hdh (用于返回) */
    Hdh = sqrt(2.0 / 3.0) * DPX;

    // ==========================================
    // 3. 结果写回 (Post-Processing)
    // ==========================================
    
    // 更新累积塑性应变 H
    H += Hdh;
    *H_ptr = H; // 写回指针指向的内存
    
    // 返回增量 (如果需要)
    if (Hdh_ptr != 0) {
        *Hdh_ptr = Hdh;
    }
    
    // 注意：原代码末尾更新 F 和 Fdh 的代码已删除
    // 因为它们会在下一次进入函数时，在开头被重新计算
}





static void Yfd2TRIELS_EP(  /* small strain elastic triangle  */
            nelem,
            iprop,
            npnfact,mprop,nprop,
            d3pnfac,
            i1ptyp,
            dpeks,dpela, dpemu, dpero ,
            dpsem,
            dpeem,dpenu,
            d1nccx,d1nccy,d1ncix,d1nciy,d1nfcx,
            d1nfcy,d1nmct,d1nvcx,d1nvcy,d1pnaf,
            d1pnap,d1pnat,
            i1elpr,i1nopr,i2elto,
            nohys, dohyp, dctime,
            d1ohys, d1ohyt, d1ohyx, d1ohyy,
            i1ohyt, npnset, d1elfr,
            i1usan, d1peex, d1peey, d1pemx, d1pemy, d1peg,
            iuseis, dcstxx, dcstxy, dcstyy,
            dcsyxx, dcsyxy, dcsyyy, dcsrfy,
            i1pnfx, i1pnfy,
            i1pexc, i1nowe, iusehf,
            nsbar, d2elstr, ncstep,
            d2elR, d2elV, d2elAlpha,d1area,dcstec,
            d1F0,d1h2,d1HC,d1Ck,d1bk,d1URU,
            sigma00, sigma10, sigma11,sigma33, sigma31, sigma32,
ss00,   ss10,    ss11,   ss33,  ss31,  ss32,   
M00,     M10,    M11,    M33,   M31,   M32,
H, R, Rc
            )
  YINT    nelem;
  YINT    iprop;
  YINT   npnfact; YINT    mprop; YINT    nprop;
  YINT npnset;
  DBL ***d3pnfac;
  YINT i1ptyp;
  DBL    dpeks; DBL   dpela; DBL    dpemu; DBL   dpero;
  DBL   dpsem; DBL   dpeem; DBL    dpenu;
  DBL *d1nccx; DBL  *d1nccy; DBL *d1ncix; DBL  *d1nciy; DBL *d1nfcx;
  DBL *d1nfcy; DBL  *d1nmct; DBL *d1nvcx; DBL  *d1nvcy; DBL *d1pnaf;
  DBL *d1pnap; DBL  *d1pnat;
  YINT *i1elpr; YINT  *i1nopr; YINT **i2elto;
  YINT   nohys; DBL    dohyp; DBL   dctime;
  DBL *d1ohys; DBL  *d1ohyt; DBL  *d1ohyx; DBL *d1ohyy;
  YINT *i1ohyt;
  DBL *d1elfr;

  YINT i1usan;
  DBL d1peex; DBL d1peey; DBL d1pemx; DBL d1pemy; DBL d1peg;
  
  YINT iuseis; DBL  dcstxx; DBL  dcstxy;  DBL  dcstyy;
  DBL  dcsyxx; DBL  dcsyxy;  DBL  dcsyyy; DBL dcsrfy; 
  YINT *i1pnfx; YINT *i1pnfy;
  YINT *i1pexc; YINT *i1nowe; YINT iusehf;
  YINT nsbar; DBL **d2elstr; YINT ncstep; 
DBL **d2elR;
DBL **d2elV;
DBL **d2elAlpha;
DBL *d1area,*sigma00, *sigma10, *sigma11,*sigma33, *sigma31, *sigma32,
*ss00, *  ss10,   * ss11,  * ss33, * ss31, * ss32,   
*M00,    * M10,    *M11,  *  M33,  * M31,   *M32;
DBL dcstec,d1F0,d1h2,d1HC,d1Ck,d1bk,d1URU;
DBL *H, * R, * Rc;

{ DBL nx,ny,voli,volc;

DBL  V[3];

/* ---------------- F&T + plane strain 临时变量 ---------------- */

/* 2D TRI3 的基础几何变量（保留） */
DBL F2D[2][2];      /* 2D deformation gradient */
DBL L2D[2][2];      /* 2D velocity gradient */
DBL F0[2][2];       /* initial local base */
DBL FX[2][2];       /* current local base */
DBL F0inv[2][2];    /* inverse of initial base */
DBL FXinv[2][2];    /* inverse of current base */
  DBL LX[2][2]; 
/* ---------------- 3D F&T 变量 ---------------- */
//DBL F3D[3][3];      /* 3D deformation gradient */
DBL L3D[3][3];      /* 3D velocity gradient */
DBL D3D[3][3];      /* rate of deformation */
DBL W3D[3][3];      /* spin tensor */

DBL Rloc[3][3], Vloc[3][3];  /* 当前单元的 R、V（从 d2elR/d2elV 读入） */
DBL wvec[3];                 /* w 向量 */



DBL z[3];           /* F&T: z vector */
DBL w[3];           /* F&T: w vector */
DBL omega[3];       /* F&T: ω vector */
DBL Omega[3][3];    /* F&T: Ω tensor */
DBL Q[3][3];        /* F&T: rotation increment */

DBL dNR[3][3];      /* non-rotational d = R^T D R */
DBL DepsNR[3][3];   /* Δε_nr = dNR * dt */

DBL T3D[3][3];      /* 3D Cauchy stress */

/* ---------------- 其他变量（保留） ---------------- */

DBL v0, v1, v2, rpx, rpy, r0x, r0y, r1x, r1y, r2x, r2y, stprev;

YINT ielem;
YINT i, j, k, in, jn, kn, jnopr, knopr;
YINT ihys;

DBL **d2fact;
DBL **d2time;

DBL Tinsitu[2][2];
DBL cp, cs, mu, lambda;

/* plane strain TI constants */
DBL d1peex_pstrain;
DBL d1peey_pstrain;
DBL d1pemx_pstrain;
DBL d1pemy_pstrain;

DBL Delta;



  if(d1pnaf == DBL1NULL)
  { d1pnaf = TalDBL1(mprop);
  }
  for(i=0; i<nprop; i++)
  { d1pnaf[i]=R1;
  }

  if(d3pnfac != DBL3NULL)
  { if(d3pnfac[0][0][0] != -R1)
    { d2time = d3pnfac[0];
      d2fact = d3pnfac[1];
      for(j=0; j<npnset; j++)
      { for(i=1; i<npnfact; i++)
        { if((dctime>=d2time[j][i-1])&&(dctime<=d2time[j][i]))
          { d1pnaf[j]=d2fact[j][i-1]-((d2fact[j][i-1]-d2fact[j][i])*
                      ((dctime-d2time[j][i-1])/(d2time[j][i]-d2time[j][i-1])));
  } } } } }
  
  for(ielem=0;ielem<nelem;ielem++)
  { if(i1elpr[ielem]==iprop)
    { for(i=1;i<3;i++)
      {F0[0][i-1]=d1ncix[(i2elto[i][ielem])]-d1ncix[(i2elto[0][ielem])];
        F0[1][i-1]=d1nciy[(i2elto[i][ielem])]-d1nciy[(i2elto[0][ielem])];
        FX[0][i-1]=d1nccx[(i2elto[i][ielem])]-d1nccx[(i2elto[0][ielem])];
        FX[1][i-1]=d1nccy[(i2elto[i][ielem])]-d1nccy[(i2elto[0][ielem])];
        LX[0][i-1]=d1nvcx[(i2elto[i][ielem])]-d1nvcx[(i2elto[0][ielem])];
        LX[1][i-1]=d1nvcy[(i2elto[i][ielem])]-d1nvcy[(i2elto[0][ielem])];
      }
      YMATINV2(F0,F0inv,voli);
      YMATINV2(FX,FXinv,volc);
       d1area[ielem]=volc/2;






      for(i=0;i<2;i++)
      { for(j=0;j<2;j++)
        { //F2D[i][j]=R0;
          L2D[i][j]=R0;
          for(k=0;k<2;k++)
          {// F2D[i][j]=F2D[i][j]+FX[i][k]*F0inv[k][j];


            L2D[i][j]=L2D[i][j]+LX[i][k]*FXinv[k][j];
      } } }

/* 
F3D[0][0] = F2D[0][0];
F3D[0][1] = F2D[0][1];
F3D[0][2] = 0.0;

F3D[1][0] = F2D[1][0];
F3D[1][1] = F2D[1][1];
F3D[1][2] = 0.0;

F3D[2][0] = 0.0;
F3D[2][1] = 0.0;
F3D[2][2] = 1.0;
 */
/* 扩展到 3D plane strain 的 L3D */
L3D[0][0] = L2D[0][0];
L3D[0][1] = L2D[0][1];
L3D[0][2] = 0.0;

L3D[1][0] = L2D[1][0];
L3D[1][1] = L2D[1][1];
L3D[1][2] = 0.0;

L3D[2][0] = 0.0;
L3D[2][1] = 0.0;
L3D[2][2] = 0.0;


for (i=0; i<3; ++i) {
  for (j=0; j<3; ++j) {
    D3D[i][j] = 0.5 * (L3D[i][j] + L3D[j][i]);
    W3D[i][j] = 0.5 * (L3D[i][j] - L3D[j][i]);
  }
}


/* ---------- 读入该单元上一时刻的 R、V ---------- */
/* d2elR[0..8][ielem] 存的是 3×3 R 依次按 (0,0),(0,1),(0,2),...,(2,2) 展开 */

for (i = 0; i < 3; ++i) {
  for (j = 0; j < 3; ++j) {
    Rloc[i][j] = d2elR[i*3 + j][ielem];
    Vloc[i][j] = d2elV[i*3 + j][ielem];
  }
}


/* ---------- 从 W3D 构造自旋向量 w ---------- */
/* w_x = 0.5 * (W23 - W32)
   w_y = 0.5 * (W31 - W13)
   w_z = 0.5 * (W12 - W21)
*/

wvec[0] = 0.5 * (W3D[1][2] - W3D[2][1]);  /* wx */
wvec[1] = 0.5 * (W3D[2][0] - W3D[0][2]);  /* wy */
wvec[2] = 0.5 * (W3D[0][1] - W3D[1][0]);  /* wz */


/* -------------------------------------------------------------------------
   Compute the F&T "z‑vector" (stretch‑related spin)
   原文公式:  z_i = ε_{ikj} D_{jm} V_{mk}
   ------------------------------------------------------------------------- */

DBL A[3][3], B[3][3];
DBL zvec[3];

/* A = D * V */
for (i = 0; i < 3; ++i) {
  for (j = 0; j < 3; ++j) {
    A[i][j] = 0.0;
    for (k = 0; k < 3; ++k) {
      A[i][j] += D3D[i][k] * Vloc[k][j];
    }
  }
}

/* B = A - A^T （反对称部分）*/
for (i = 0; i < 3; ++i) {
  for (j = 0; j < 3; ++j) {
    B[i][j] = A[i][j] - A[j][i];
  }
}

/* z = 1/2 * axial(B) */
zvec[0] = 0.5 * (B[2][1] - B[1][2]);  /* zx */
zvec[1] = 0.5 * (B[0][2] - B[2][0]);  /* zy */
zvec[2] = 0.5 * (B[1][0] - B[0][1]);  /* zz */

// 计算Vloc的迹（trace(V)）
DBL trace_V = Vloc[0][0] + Vloc[1][1] + Vloc[2][2];

DBL M[3][3];  // 定义矩阵M = trace(V)*I - V
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        if (i == j) {
            M[i][j] = trace_V * 1.0 - Vloc[i][j];  // I的对角线元素为1
        } else {
            M[i][j] = 0.0 - Vloc[i][j];  // I的非对角线元素为0
        }
    }
}


// ---------------- 步骤3：调用现有YMATINV3宏求M的逆（核心修改）----------------
DBL M_inv[3][3];  // 逆矩阵
DBL detM_inv;          // 矩阵行列式（宏的输出参数）

// 调用已定义的求逆宏：参数依次为「原矩阵M、逆矩阵M_inv、行列式det」
YMATINV3(M, M_inv, detM_inv);



// ---------------- 步骤4：计算M_inv · zvec（矩阵×向量，不变）----------------
DBL M_inv_z[3] = {0.0};
for (int i = 0; i < 3; i++) {
    for (int k = 0; k < 3; k++) {
        M_inv_z[i] += M_inv[i][k] * zvec[k];
    }
}

// ---------------- 步骤5：叠加得到ω向量（不变）----------------
DBL omega[3];
for (int i = 0; i < 3; i++) {
    omega[i] = wvec[i] + M_inv_z[i];
    // 数值降噪：极小值置0
    if (fabs(omega[i]) < 1e-12) omega[i] = 0.0;
}

// 定义Ω张量（3×3反对称，对应原文式(3-9)）
DBL Omega[3][3] = {0.0};  // 初始化所有元素为0，简化对角线赋值

// 按置换张量规则赋值非对角线元素（核心逻辑）
// Ω[0][1] = -ω_z, Ω[1][0] = ω_z
Omega[0][1] = -omega[2];
Omega[1][0] = omega[2];//这里赋值务必重思考

// Ω[0][2] = ω_y, Ω[2][0] = -ω_y
Omega[0][2] = omega[1];
Omega[2][0] = -omega[1];

// Ω[1][2] = -ω_x, Ω[2][1] = ω_x
Omega[1][2] = -omega[0];
Omega[2][1] = omega[0];

// 数值降噪（可选，与omega向量处理逻辑一致）
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        if (fabs(Omega[i][j]) < 1e-12) {
            Omega[i][j] = 0.0;
        }
    }
}

// 步骤11：根据ω向量计算标量Ω（原文式(3-15)）
DBL Omega_scalar;  // 注意：区分之前的Ω张量（Omega[3][3]），命名加scalar避免混淆
// 1. 计算ω向量各分量的平方和
DBL omega_sq_sum = omega[0]*omega[0] + omega[1]*omega[1] + omega[2]*omega[2];
// 2. 开方得到标量Ω（补充数值稳定性：极小值置0，避免后续除以0）
if (omega_sq_sum < 1e-24) {  // 比1e-12更严格，避免开方后出现微小噪声
    Omega_scalar = 0.0;
} else {
    Omega_scalar = sqrt(omega_sq_sum);
}

// 步骤12：计算正交张量Q（原文式(3-14)）
DBL Q[3][3];          // 目标Q张量
DBL Omega_sq[3][3];   // 存储Ω张量的平方（Ω·Ω）
DBL theta;            // Δt·Ω（旋转角度）
DBL sin_theta_over_O; // sin(theta)/Ω
DBL one_minus_cos_over_O2; // (1 - cos(theta))/Ω²

// ---------------- 步骤1：初始化Q为单位矩阵（基础项I）----------------
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        Q[i][j] = (i == j) ? 1.0 : 0.0; // 单位矩阵I
    }
}

// ---------------- 步骤2：处理Ω=0的特殊情况（直接返回单位矩阵）----------------
if (Omega_scalar < 1e-12) {
    // 无旋转，Q=I，无需后续计算
    goto Q_COMPLETE; // 跳转到结束标记，避免无效计算
}

// ---------------- 步骤3：计算theta和两个核心系数 ----------------
theta = dcstec * Omega_scalar; // Δt·Ω（dctime是时间步长）
sin_theta_over_O = sin(theta) / Omega_scalar;
one_minus_cos_over_O2 = (1.0 - cos(theta)) / (Omega_scalar * Omega_scalar);

// ---------------- 步骤4：计算Ω张量的平方（Ω·Ω）----------------
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        Omega_sq[i][j] = 0.0;
        for (int k = 0; k < 3; k++) {
            Omega_sq[i][j] += Omega[i][k] * Omega[k][j]; // 矩阵乘法
        }
    }
}

// ---------------- 步骤5：按公式叠加项到Q张量 ----------------
// 项1：[sin(theta)/Ω]·Ω
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        Q[i][j] += sin_theta_over_O * Omega[i][j];
    }
}

// 项2：-[(1 - cos(theta))/Ω²]·Ω²
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        Q[i][j] -= one_minus_cos_over_O2 * Omega_sq[i][j];
    }
}

// ---------------- 标记：Q计算完成 ----------------
Q_COMPLETE:
// 可选：数值降噪（极小值置0，保持张量正交性）
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        if (fabs(Q[i][j]) < 1e-12) {
            Q[i][j] = 0.0;
        }
    }
}



// 步骤14：更新旋转张量R_t = Q · R_{t-Δt}（原文式(3-13)）
DBL R_new[3][3] = {0.0};  // 存储当前时刻的R_t

// ---------------- 步骤1：执行矩阵乘法 Q · Rloc ----------------
// 矩阵乘法规则：R_new[i][j] = Σ(k=0→2) Q[i][k] * Rloc[k][j]
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
            R_new[i][j] += Q[i][k] * Rloc[k][j];
        }
    }
}

// ---------------- 步骤2：数值稳定性处理（可选）----------------
// 1. 极小值置0，避免浮点噪声
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        if (fabs(R_new[i][j]) < 1e-12) {
            R_new[i][j] = 0.0;
        }
    }
}
// 【新增】平面应变强制约束：
R_new[0][2] = 0.0; R_new[1][2] = 0.0; // 第3列 (0, 1)
R_new[2][0] = 0.0; R_new[2][1] = 0.0; // 第3行 (0, 1)
R_new[2][2] = 1.0;                    // z轴方向保持为1 (无旋转分量混入)
// 2. 正交性修正（可选，确保R_t是正交张量：R·R^T=I）
// 若需严格正交，可添加Gram-Schmidt正交化，新手阶段可先省略

// ---------------- 步骤3：将R_new写回d2elR数组 ----------------
// d2elR的存储规则：d2elR[i*3 + j][ielem] = R_new[i][j]（行优先展开）
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        d2elR[i*3 + j][ielem] = R_new[i][j];
    }
}


// 步骤15：更新左拉伸张量V_t（原文式(3-8)+显式积分）
DBL dV[3][3] = {0.0};    // 存储V的时间导数dot(V)
DBL V_new[3][3] = {0.0}; // 存储当前时刻的V_t
DBL L_V[3][3] = {0.0};   // 临时存储L·V
DBL V_Omega[3][3] = {0.0}; // 临时存储V·Ω

// ---------------- 步骤1：计算第一项 L·V ----------------
// 矩阵乘法：L_V[i][j] = Σ(k=0→2) L3D[i][k] * Vloc[k][j]
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
            L_V[i][j] += L3D[i][k] * Vloc[k][j];
        }
    }
}

// ---------------- 步骤2：计算第二项 V·Ω ----------------
// 矩阵乘法：V_Omega[i][j] = Σ(k=0→2) Vloc[i][k] * Omega[k][j]
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
            V_Omega[i][j] += Vloc[i][k] * Omega[k][j];
        }
    }
}

// ---------------- 步骤3：计算dot(V) = L·V - V·Ω ----------------
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        dV[i][j] = L_V[i][j] - V_Omega[i][j];
    }
}

// ---------------- 步骤4：显式积分更新V_t = V_{t-Δt} + dot(V)*Δt ----------------
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        V_new[i][j] = Vloc[i][j] + dV[i][j] * dcstec;
    }
}

// ---------------- 步骤5：数值稳定性处理 ----------------
// 1. 极小值置0，避免浮点噪声（平面应变下V_new[2][2]应≈1）
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        if (fabs(V_new[i][j]) < 1e-12) {
            V_new[i][j] = 0.0;
        }
        // 平面应变特殊处理：强制z方向拉伸为1（避免数值漂移）
        if (i == 2 && j == 2) {
            V_new[i][j] = 1.0;
        }
    }
}
// 【新增】平面应变强制约束：
V_new[0][2] = 0.0; V_new[1][2] = 0.0; // 剪切 Vxz, Vyz = 0
V_new[2][0] = 0.0; V_new[2][1] = 0.0; // 对称部分
V_new[2][2] = 1.0;                    // z方向拉伸比恒为 1.0
// ---------------- 步骤6：将V_new写回d2elV数组 ----------------
// d2elV存储规则与d2elR一致：d2elV[i*3 + j][ielem] = V_new[i][j]
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        d2elV[i*3 + j][ielem] = V_new[i][j];
    }
}


// 步骤17：计算非旋转构型下的应变率张量d = R_t^T · D3D · R_t（原文式(3-6)）
DBL R_T[3][3] = {0.0};       // 存储R_t的转置（R_new^T）
DBL temp1[3][3] = {0.0};     // 临时存储R_T · D3D的结果

// ---------------- 步骤1：计算R_t的转置（R_new^T）----------------
// 转置规则：R_T[i][j] = R_new[j][i]
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        R_T[i][j] = R_new[j][i];
    }
}

// ---------------- 步骤2：计算第一步乘法：temp1 = R_T · D3D ----------------
// 矩阵乘法规则：temp1[i][j] = Σ(k=0→2) R_T[i][k] * D3D[k][j]
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
            temp1[i][j] += R_T[i][k] * D3D[k][j];
        }
    }
}

// ---------------- 步骤3：计算第二步乘法：dNR = temp1 · R_new ----------------
// dNR即为非旋转构型下的d张量，初始化清零
memset(dNR, 0, sizeof(dNR)); // 确保初始值为0（需包含<string.h>头文件，若未包含可手动赋值）
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
            dNR[i][j] += temp1[i][k] * R_new[k][j];
        }
    }
}

// ---------------- 步骤4：数值稳定性处理 ----------------
// 1. 极小值置0，避免浮点噪声
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        if (fabs(dNR[i][j]) < 1e-12) {
            dNR[i][j] = 0.0;
        }
    }
}
// 2. 平面应变特殊处理：z方向应变率强制为0（无out-of-plane变形）
dNR[2][0] = dNR[2][1] = dNR[2][2] = 0.0;
dNR[0][2] = dNR[1][2] = 0.0;



// 步骤18：计算非旋转构型下的应变增量Δε_nr = d · Δt
// 初始化DepsNR为0（避免残留值干扰）
memset(DepsNR, 0, sizeof(DepsNR)); // 需包含<string.h>，若未包含可手动循环赋值

// 核心计算：逐分量乘以时间步长
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        DepsNR[i][j] = dNR[i][j] * dcstec;
        
        // 数值稳定性：极小值置0，避免浮点噪声
        if (fabs(DepsNR[i][j]) < 1e-12) {
            DepsNR[i][j] = 0.0;
        }
        
        // 平面应变特殊处理：强制z方向应变增量为0
        if ((i == 2) || (j == 2)) {
            DepsNR[i][j] = 0.0;
        }
    }
}


// （之前已完成）步骤1-18：计算得到非旋转构型应变增量DepsNR[3][3]

// ---------------- 适配步骤：将3×3的DepsNR转换为你算法的6维应变向量 ----------------
// 你的算法用6维向量存储应变（对应工程应变：[εxx, εyy, εzz, γxy, γyz, γzx]）
// 非旋转构型下DepsNR是3×3对称张量，转换规则：
double De_m[6];  // 对应你算法中的应变增量输入
De_m[0] = DepsNR[0][0];          // εxx
De_m[1] = DepsNR[1][1];          // εyy
De_m[2] = DepsNR[2][2];          // εzz（平面应变下为0）
De_m[3] = 2.0 * DepsNR[0][1];    // γxy = 2εxy（工程剪应变定义）
De_m[4] = 2.0 * DepsNR[1][2];    // γyz = 2εyz（平面应变下为0）
De_m[5] = 2.0 * DepsNR[0][2];    // γzx = 2εzx（平面应变下为0）
double SIG[6];  // 应力输入
SIG[0] = sigma00[ielem];         
SIG[1] = sigma11[ielem];       
SIG[2] = sigma33[ielem];          
SIG[3] = sigma10[ielem];   
SIG[4] = sigma31[ielem];    
SIG[5] = sigma32[ielem];    
double Ss[6];  // 背应力
Ss[0] = ss00[ielem];          
Ss[1] = ss11[ielem];         
Ss[2] = ss33[ielem];        
Ss[3] = ss10[ielem];    
Ss[4] = ss31[ielem];    
Ss[5] = ss32[ielem];    
double EP_M[6];  //塑性应变增量
EP_M[0] = M00[ielem];       
EP_M[1] = M11[ielem];        
EP_M[2] = M33[ielem];        
EP_M[3] = M10[ielem];    
EP_M[4] = M31[ielem];  
EP_M[5] = M32[ielem];   
 
int C_m[6] = {0, 0, 0, 0, 0, 0};
// ---------------- 调用你的弹塑性本构算法 ----------------
// 注意：此时你的算法计算的是「非旋转构型下的应力σ」（对应文档中的σ）
calc_exsub(
  
  d1F0,d1h2,d1HC,d1Ck,d1bk,d1URU,dpeem,dpenu,

 SIG,Ss, De_m,C_m,EP_M,
&H[ielem], &Rc[ielem]); 


// 用计算出的真实应力 T3D 更新全局存储
sigma00[ielem] = SIG[0];
sigma11[ielem] = SIG[1];
sigma33[ielem] = SIG[2];
sigma10[ielem] = SIG[3]; // 注意对称性 xy
sigma31[ielem] = SIG[4]; // yz (注意你的索引对应关系)
sigma32[ielem] = SIG[5]; // zx

ss00[ielem] = Ss[0];
ss11[ielem] = Ss[1];
ss33[ielem] = Ss[2];
ss10[ielem] = Ss[3];
ss31[ielem] = Ss[4];
ss32[ielem] = Ss[5];

M00[ielem] = EP_M[0];
M11[ielem] = EP_M[1];
M33[ielem] = EP_M[2];
M10[ielem] = EP_M[3];
M31[ielem] = EP_M[4];
M32[ielem] = EP_M[5];

// ---------------- 适配步骤：将6维σ转换为3×3非旋转应力张量 ----------------
double sigma_nr[3][3] = {0.0};  // 非旋转构型下的σ张量
sigma_nr[0][0] = SIG[0];
sigma_nr[1][1] = SIG[1];
sigma_nr[2][2] = SIG[2];
sigma_nr[0][1] = sigma_nr[1][0] = SIG[3] ; 
sigma_nr[1][2] = sigma_nr[2][1] = SIG[4] ;
sigma_nr[0][2] = sigma_nr[2][0] = SIG[5] ;


// ---------------- 步骤24：计算真实柯西应力T = R·σ·R^T + η·D ----------------
DBL T3D[3][3] = {0.0};  // 真实柯西应力（最终用于节点力计算）
DBL sigma_R[3][3] = {0.0};  // 临时存储σ·R^T
DBL R_sigma_R[3][3] = {0.0};  // 临时存储R·(σ·R^T)
DBL eta_D[3][3] = {0.0};  // 粘性阻尼项（η·D）


// 1. 计算σ·R^T（R^T是R_new的转置，已在步骤17中计算为R_T）
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
            sigma_R[i][j] += sigma_nr[i][k] * R_T[k][j];
        }
    }
}

// 2. 计算R·(σ·R^T) → R·σ·R^T
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
            R_sigma_R[i][j] += R_new[i][k] * sigma_R[k][j];
        }
    }
}

// 3. 计算粘性阻尼项η·D（D=D3D，变形率张量）
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        eta_D[i][j] = dpeks * D3D[i][j];
    }
}

// 4. 叠加得到真实应力T
for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
        T3D[i][j] = R_sigma_R[i][j] + eta_D[i][j];
    }
}

// ---------------- 步骤25：转换为节点力（与原FDEM逻辑一致） ----------------
// 后续按原FDEM流程，将T3D转换为等效节点力，叠加到d1elfr中，无需修改

      if((ncstep%1==0)&& (ielem==89))
      { 

      { printf("T:%.10f,%.10f,%.10f,%.10f,%.10f,%1d\n", T3D[0][0], T3D[0][1], T3D[1][1], T3D[2][2],  ncstep);
  } 
  }
      /* output history states */
      r0x = d1ncix[(i2elto[0][ielem])];
      r0y = d1nciy[(i2elto[0][ielem])];
      r1x = d1ncix[(i2elto[1][ielem])];
      r1y = d1nciy[(i2elto[1][ielem])];
      r2x = d1ncix[(i2elto[2][ielem])];
      r2y = d1nciy[(i2elto[2][ielem])];
      for(ihys=0; ihys<nohys; ihys++)
      { 
	if(i1ohyt[ihys]==(YFLEE)) /* Strain energy */
	{
	  d1ohyt[ihys] = dctime;    /* output history time  */
          if(i1usan==1) /* Transversely isotropic elasticity */
	  {   //d1ohys[ihys] += 0.5 * ((E[0][0] * ((d1peex / (1 - d1pemx * d1pemy)) * (E[0][0] + d1pemy * E[1][1]))) +      /* output history state */
        //                           (E[0][1] * ((1 / (2 * d1peg)) * E[0][1])) +
        //                           (E[1][0] * ((1 / (2 * d1peg)) * E[0][1])) +
        //                           (E[1][1] * ((d1peey / (1 - d1pemx * d1pemy)) * (E[1][1] + d1pemx * E[0][0])))) * (volc/2);    
	  }
	  else if(i1usan==2) /* Transversely isotropic elasticity (plane strain) */
	  { //d1peex_pstrain = d1peex / (1 - d1pemx * d1pemy);
       //     d1peey_pstrain = d1peey / (1 - d1pemy * d1pemx);
	    //d1pemx_pstrain = (d1pemx + d1pemx * d1pemy) / (1 - d1pemx * d1pemy);
	    //d1pemy_pstrain = (d1pemy + d1pemx * d1pemy) / (1 - d1pemx * d1pemy);
	    //d1ohys[ihys] += 0.5 * ((E[0][0] * ((d1peex_pstrain / (1 - d1pemx_pstrain * d1pemy_pstrain)) * (E[0][0] + d1pemy * E[1][1]))) +      /* output history state */
      //                             (E[0][1] * ((1 / (2 * d1peg)) * E[0][1])) +
      //                             (E[1][0] * ((1 / (2 * d1peg)) * E[0][1])) +
       //                            (E[1][1] * ((d1peey_pstrain / (1 - d1pemx_pstrain * d1pemy_pstrain)) * (E[1][1] + d1pemx_pstrain * E[0][0])))) * (volc/2);    
	  }
	  else /* Isotropic elasticity */
	  { //d1ohys[ihys] += 0.5 * ((E[0][0] * (2.0*dpemu*E[0][0]*(voli/volc) + dpela*(volc/voli-voli/volc))) +      /* output history state */
      //                           (E[0][1] * (2.0*dpemu*E[0][1]*(voli/volc))) +
       //                          (E[1][0] * (2.0*dpemu*E[1][0]*(voli/volc))) +
        //                         (E[1][1] * (2.0*dpemu*E[1][1]*(voli/volc) + dpela*(volc/voli-voli/volc)))) * (volc/2);    
	  }
	}
	rpx=d1ohyx[ihys];  /* x coordinate of point P */
        rpy=d1ohyy[ihys];  /* y coordinate of point P */
        V2DCro(v0,(r1x-r0x),(r1y-r0y),(rpx-r0x),(rpy-r0y));
        V2DCro(v1,(r2x-r1x),(r2y-r1y),(rpx-r1x),(rpy-r1y));
        V2DCro(v2,(r0x-r2x),(r0y-r2y),(rpx-r2x),(rpy-r2y));

        if((v0>R0)&&(v1>R0)&&(v2>R0))    /* if point is inside the triangle */
        { if(i1ohyt[ihys]==(YFLDSXX))
          { stprev=MAXIM((EPSILON),(ABS(d1ohys[ihys])));
            if((ABS(T3D[0][0]-stprev))>=dohyp)
            { d1ohyt[ihys] = dctime;    /* output history time  */
              d1ohys[ihys] = T3D[0][0];    /* output history state */
          } }
          else if(i1ohyt[ihys]==(YFLDSXY))
          { stprev=MAXIM((EPSILON),(ABS(d1ohys[ihys])));
            if((ABS(T3D[0][1]-stprev))>=dohyp)  /* if((ABS(R1-T[0][1]/stprev))>dohyp) */
            { d1ohyt[ihys] = dctime;
              d1ohys[ihys] = T3D[0][1];
          } }
          else if(i1ohyt[ihys]==(YFLDSYY))
          { stprev=MAXIM((EPSILON),(ABS(d1ohys[ihys])));
            if((ABS(T3D[1][1]-stprev))>=dohyp)
            { d1ohyt[ihys] = dctime;
              d1ohys[ihys] = T3D[1][1];
          } }
          else if(i1ohyt[ihys]==(YFLDSZZ))
          { d1ohyt[ihys] = dctime;
            d1ohys[ihys] = R0;
          }
          else if(i1ohyt[ihys]==(YFLDSZX))
          { d1ohyt[ihys] = dctime;
            d1ohys[ihys] = R0;
          }
          else if(i1ohyt[ihys]==(YFLDSZY))
          { d1ohyt[ihys] = dctime;
            d1ohys[ihys] = R0;
          }
          else if(i1ohyt[ihys]==(YFLDVEL))
          { stprev=MAXIM((EPSILON),(ABS(d1ohys[ihys])));
            for(i=0; i<3; i++)
            { V[i] = SQRT((d1nvcx[(i2elto[i][ielem])]*d1nvcx[(i2elto[i][ielem])])
                         +(d1nvcy[(i2elto[i][ielem])]*d1nvcy[(i2elto[i][ielem])]));
            }
            if((ABS(((V[0]+V[1]+V[2])/R3)-stprev))>=dohyp)
            { d1ohyt[ihys] = dctime;
              d1ohys[ihys] = (V[0]+V[1]+V[2])/R3;  /* average velocity of element */
          } }
          else if(i1ohyt[ihys]==(YFLDVEX))
          { stprev=MAXIM((EPSILON),(ABS(d1ohys[ihys])));
            for(i=0; i<3; i++)
            { V[i] = d1nvcx[(i2elto[i][ielem])];
            }
            if((ABS(((V[0]+V[1]+V[2])/R3)-stprev))>=dohyp)
            { d1ohyt[ihys] = dctime;
              d1ohys[ihys] = (V[0]+V[1]+V[2])/R3;  /* average velocity x of element */
          } }
          else if(i1ohyt[ihys]==(YFLDVEY))
          { stprev=MAXIM((EPSILON),(ABS(d1ohys[ihys])));
            for(i=0; i<3; i++)
            { V[i] = d1nvcy[(i2elto[i][ielem])];
            }
            if((ABS(((V[0]+V[1]+V[2])/R3)-stprev))>=dohyp)
            { d1ohyt[ihys] = dctime;
              d1ohys[ihys] = (V[0]+V[1]+V[2])/R3;  /* average velocity y of element */
            } 
          } 
      } }
       
      /* Store element stress tensor only if rebars are used */ 
      if(nsbar>0)
      { d2elstr[0][ielem]=T3D[0][0];
        d2elstr[1][ielem]=T3D[0][1];
        d2elstr[2][ielem]=T3D[1][0]; 
        d2elstr[3][ielem]=T3D[1][1];
      }

      /* Nodal Forces */
      for(i=0;i<3;i++)
      { j=i+1; if(j>2)j=0;
        k=j+1; if(k>2)k=0;
        in=i2elto[i][ielem];
        jn=i2elto[j][ielem];
        kn=i2elto[k][ielem];

        nx=d1nccy[kn]-d1nccy[jn];
        ny=d1nccx[jn]-d1nccx[kn];
        d1nmct[in]=d1nmct[in]+dpero*voli/YR6;
        d1nfcx[in]=d1nfcx[in]+(T3D[0][0]*nx+T3D[0][1]*ny)/R2;
        d1nfcy[in]=d1nfcy[in]+(T3D[1][0]*nx+T3D[1][1]*ny)/R2;

	if (iuseis==1) /* Apply element in-situ stress */
        { Tinsitu[0][0]=-dcstxx-dcsyxx*((d1nciy[i2elto[0][ielem]]+d1nciy[i2elto[1][ielem]]+d1nciy[i2elto[2][ielem]])/3.0 - dcsrfy);
          Tinsitu[0][1]=-dcstxy-dcsyxy*((d1nciy[i2elto[0][ielem]]+d1nciy[i2elto[1][ielem]]+d1nciy[i2elto[2][ielem]])/3.0 - dcsrfy);
          Tinsitu[1][0]=-dcstxy-dcsyxy*((d1nciy[i2elto[0][ielem]]+d1nciy[i2elto[1][ielem]]+d1nciy[i2elto[2][ielem]])/3.0 - dcsrfy);
          Tinsitu[1][1]=-dcstyy-dcsyyy*((d1nciy[i2elto[0][ielem]]+d1nciy[i2elto[1][ielem]]+d1nciy[i2elto[2][ielem]])/3.0 - dcsrfy);
	  /* Nodal Forces due to in-situ stress */
	  d1nfcx[in]=d1nfcx[in]+(Tinsitu[0][0]*nx+Tinsitu[0][1]*ny)/R2;
	  d1nfcy[in]=d1nfcy[in]+(Tinsitu[1][0]*nx+Tinsitu[1][1]*ny)/R2;
	}
	
       /* Nodal Forces due to edge force*/
        jnopr=i1nopr[jn];
        knopr=i1nopr[kn];
        if( ((DABS(d1pnap[jnopr]))>EPSILON)&&
            ((DABS(d1pnap[knopr]))>EPSILON) )
        { d1nfcx[jn]=d1nfcx[jn]-
          d1pnap[jnopr]*d1pnaf[jnopr]*nx/R3-
          d1pnap[knopr]*d1pnaf[knopr]*nx/YR6;
          d1nfcy[jn]=d1nfcy[jn]-
          d1pnap[jnopr]*d1pnaf[jnopr]*ny/R3-
          d1pnap[knopr]*d1pnaf[knopr]*ny/YR6;
          d1nfcx[kn]=d1nfcx[kn]-
          d1pnap[jnopr]*d1pnaf[jnopr]*nx/YR6-
          d1pnap[knopr]*d1pnaf[knopr]*nx/R3;
          d1nfcy[kn]=d1nfcy[kn]-
          d1pnap[jnopr]*d1pnaf[jnopr]*ny/YR6-
          d1pnap[knopr]*d1pnaf[knopr]*ny/R3;
        }
        if( ((DABS(d1pnat[jnopr]))>EPSILON)&&
            ((DABS(d1pnat[knopr]))>EPSILON) )
        { d1nfcx[jn]=d1nfcx[jn]-
          d1pnat[jnopr]*d1pnaf[jnopr]*ny/R3-
          d1pnat[knopr]*d1pnaf[knopr]*ny/YR6;
          d1nfcy[jn]=d1nfcy[jn]+
          d1pnat[jnopr]*d1pnaf[jnopr]*nx/R3+
          d1pnat[knopr]*d1pnaf[knopr]*nx/YR6;
          d1nfcx[kn]=d1nfcx[kn]-
          d1pnat[jnopr]*d1pnaf[jnopr]*ny/YR6-
          d1pnat[knopr]*d1pnaf[knopr]*ny/R3;
          d1nfcy[kn]=d1nfcy[kn]+
          d1pnat[jnopr]*d1pnaf[jnopr]*nx/YR6+
          d1pnat[knopr]*d1pnaf[knopr]*nx/R3;
        }
        /* Nodal forces due to absorbing boundary condition */
        if (i1ptyp == YTE2TRIELS)
        { cp=SQRT((2*dpemu+dpela)/dpero);
          cs=SQRT(dpemu/dpero);
        }
        else if (i1ptyp == YTE2PLANESTRESS)
        { lambda=(dpeem*dpenu)/((1+dpenu)*(1-2*dpenu));
          mu=(dpeem)/(2*(1+dpenu));
          cp=SQRT((2*mu+lambda)/dpero);
          cs=SQRT(mu/dpero);
        }  
        else if (i1ptyp == YTE2PLANESTRAIN)
        { lambda=(dpeem*dpenu)/((1+dpenu)*(1-2*dpenu));
          mu=(dpeem)/(2*(1+dpenu));
          cp=SQRT((2*mu+lambda)/dpero);
          cs=SQRT(mu/dpero);
        }
        if ((i1pnfx[jnopr]==4)&&(i1pnfx[knopr]==4))
        { cp=SQRT((2*dpemu+dpela)/dpero);
          cs=SQRT(dpemu/dpero);
          d1nfcx[jn]=d1nfcx[jn]-dpero*cp*d1nvcx[jn]*ABS(nx)/R3-
                     dpero*cp*d1nvcx[kn]*ABS(nx)/YR6-
                     dpero*cs*d1nvcx[jn]*ABS(ny)/R3-
                     dpero*cs*d1nvcx[kn]*ABS(ny)/YR6;
          d1nfcx[kn]=d1nfcx[kn]-dpero*cp*d1nvcx[jn]*ABS(nx)/YR6-
                     dpero*cp*d1nvcx[kn]*ABS(nx)/R3-
                     dpero*cs*d1nvcx[jn]*ABS(ny)/YR6-
                     dpero*cs*d1nvcx[kn]*ABS(ny)/R3;
        }
        if ((i1pnfy[jnopr]==4)&&(i1pnfy[knopr]==4))
        { cp=SQRT((2*dpemu+dpela)/dpero);
          cs=SQRT(dpemu/dpero);
          d1nfcy[jn]=d1nfcy[jn]-dpero*cp*d1nvcy[jn]*ABS(ny)/R3-
                     dpero*cp*d1nvcy[kn]*ABS(ny)/YR6-
                     dpero*cs*d1nvcy[jn]*ABS(nx)/R3-
                     dpero*cs*d1nvcy[kn]*ABS(nx)/YR6;
          d1nfcy[kn]=d1nfcy[kn]-dpero*cp*d1nvcy[jn]*ABS(ny)/YR6-
          dpero*cp*d1nvcy[kn]*ABS(ny)/R3-
          dpero*cs*d1nvcy[jn]*ABS(nx)/YR6-
          dpero*cs*d1nvcy[kn]*ABS(nx)/R3;
        }
      }
      //! If element is "excavated" set nodal boundary condition to v_x = 0 and v_y = 0 (hardcoded bc ID)
      //! and translate by a constant vector (hardcoded, 200)
      if(i1pexc[i1elpr[ielem]]==1)
      { i1nopr[i2elto[0][ielem]] = 1; /* 1 is a hard-coded value */
        i1nopr[i2elto[1][ielem]] = 1; /* 1 is a hard-coded value */
	i1nopr[i2elto[2][ielem]] = 1; /* 1 is a hard-coded value */
	Delta = 200;
	d1nccx[(i2elto[0][ielem])]=d1ncix[(i2elto[0][ielem])]+Delta;
	d1nccx[(i2elto[1][ielem])]=d1ncix[(i2elto[1][ielem])]+Delta;
	d1nccx[(i2elto[2][ielem])]=d1ncix[(i2elto[2][ielem])]+Delta;
	d1nccy[(i2elto[0][ielem])]=d1nciy[(i2elto[0][ielem])]+Delta;
	d1nccy[(i2elto[1][ielem])]=d1nciy[(i2elto[1][ielem])]+Delta;
	d1nccy[(i2elto[2][ielem])]=d1nciy[(i2elto[2][ielem])]+Delta;
	if(iusehf==1)
	{ //! For hydrofrac: set nodes to be non-wettable
	  i1nowe[i2elto[0][ielem]]=3;
	  i1nowe[i2elto[1][ielem]]=3;
	  i1nowe[i2elto[2][ielem]]=3;
	}
      }
    }
  }
  FREE(d1pnaf);
}



static void Yfd2TRIELS(  /* small strain elastic triangle  */
            nelem,
            iprop,
            npnfact,mprop,nprop,
            d3pnfac,
            i1ptyp,
            dpeks,dpela, dpemu, dpero ,
            dpsem,
            dpeem,dpenu,
            d1nccx,d1nccy,d1ncix,d1nciy,d1nfcx,
            d1nfcy,d1nmct,d1nvcx,d1nvcy,d1pnaf,
            d1pnap,d1pnat,
            i1elpr,i1nopr,i2elto,
            nohys, dohyp, dctime,
            d1ohys, d1ohyt, d1ohyx, d1ohyy,
            i1ohyt, npnset, d1elfr,
            i1usan, d1peex, d1peey, d1pemx, d1pemy, d1peg,
            iuseis, dcstxx, dcstxy, dcstyy,
            dcsyxx, dcsyxy, dcsyyy, dcsrfy,
            i1pnfx, i1pnfy,
            i1pexc, i1nowe, iusehf,
            nsbar, d2elstr, ncstep,
            sigma00, sigma10, sigma11, sigma33,d1area,dcstec,
            Fe00, Fe01, Fe10, Fe11,Fe33,d1noa0,d1noac,i1getmaster,mcstep,i1nobf0,d1J
            )
  YINT    nelem;
  YINT    iprop;
  YINT   npnfact; YINT    mprop; YINT    nprop;
  YINT npnset;
  DBL ***d3pnfac;
  YINT i1ptyp;
  DBL    dpeks; DBL   dpela; DBL    dpemu; DBL   dpero;
  DBL   dpsem; DBL   dpeem; DBL    dpenu;
  DBL *d1nccx; DBL  *d1nccy; DBL *d1ncix; DBL  *d1nciy; DBL *d1nfcx;
  DBL *d1nfcy; DBL  *d1nmct; DBL *d1nvcx; DBL  *d1nvcy; DBL *d1pnaf;
  DBL *d1pnap; DBL  *d1pnat;
  YINT *i1elpr; YINT  *i1nopr; YINT **i2elto;
  YINT   nohys; DBL    dohyp; DBL   dctime;
  DBL *d1ohys; DBL  *d1ohyt; DBL  *d1ohyx; DBL *d1ohyy;
  YINT *i1ohyt;
  DBL *d1elfr;

  YINT i1usan;
  DBL d1peex; DBL d1peey; DBL d1pemx; DBL d1pemy; DBL d1peg;
  
  YINT iuseis; DBL  dcstxx; DBL  dcstxy;  DBL  dcstyy;
  DBL  dcsyxx; DBL  dcsyxy;  DBL  dcsyyy; DBL dcsrfy; 
  YINT *i1pnfx; YINT *i1pnfy;
  YINT *i1pexc; YINT *i1nowe; YINT iusehf;
  YINT nsbar; DBL **d2elstr; YINT ncstep; 
  DBL   **sigma00; DBL **sigma10; DBL **sigma11; DBL **sigma33; DBL *d1area;DBL dcstec;
   DBL   **Fe00; DBL **Fe01; DBL **Fe10; DBL **Fe11;DBL **Fe33;
   
   DBL *d1noa0;DBL *d1noac;YINT* i1getmaster;YINT mcstep;YINT *i1nobf0;DBL *d1J;
{ DBL nx,ny,voli,volc;
  DBL v0, v1, v2, rpx, rpy, r0x, r0y, r1x, r1y, r2x, r2y, stprev;

  DBL  jbarno[3];//tijisuoding
  DBL  Jbar;//tijisuoding
  DBL  V[3];
  DBL  B[2][2];     /* left Cauchy-Green strain tensor                        */
  DBL  D[2][2];     /* rate of deformation (stretching) tensor                */
  DBL  E[2][2];     /* strain tensor (small strains)                          */
  DBL  F[2][2];     /* deformation gradient in global base                    */
  DBL F0[2][2];     /* initial local base                                     */
  DBL FX[2][2];     /* current local base                                     */
  DBL F0inv[2][2];  /* global base in initial local base                      */
  DBL FXinv[2][2];  /* global base in current local base                      */
  DBL  L[2][2];     /* velocity gradient in global base                       */
  DBL LX[2][2];     /* vel. gradient in current local base = delta x/delta X  */
  DBL  T[2][2];     /* Cauchy stress                                          */
  YINT ielem;
  YINT i,j,k,in,jn,kn,jnopr,knopr;
  YINT ihys;
  DBL **d2fact;
  DBL **d2time;
  
  DBL Tinsitu[2][2]; /* element in-situ stress */
  DBL cp; /* element p-wave velocity */
  DBL cs; /* element s-wave velocity */
  DBL mu;
  DBL lambda;
  
  /* Elastic constants for plane strain transverse isotropy */
  DBL d1peex_pstrain; 
  DBL d1peey_pstrain;
  DBL d1pemx_pstrain;
  DBL d1pemy_pstrain;
  
  DBL Delta;
  
  if(d1pnaf == DBL1NULL)
  { d1pnaf = TalDBL1(mprop);
  }
  for(i=0; i<nprop; i++)
  { d1pnaf[i]=R1;
  }

  if(d3pnfac != DBL3NULL)
  { if(d3pnfac[0][0][0] != -R1)
    { d2time = d3pnfac[0];
      d2fact = d3pnfac[1];
      for(j=0; j<npnset; j++)
      { for(i=1; i<npnfact; i++)
        { if((dctime>=d2time[j][i-1])&&(dctime<=d2time[j][i]))
          { d1pnaf[j]=d2fact[j][i-1]-((d2fact[j][i-1]-d2fact[j][i])*
                      ((dctime-d2time[j][i-1])/(d2time[j][i]-d2time[j][i-1])));
  } } } } }
  
  for(ielem=0;ielem<nelem;ielem++)
  { if(i1elpr[ielem]==iprop)
    { for(i=1;i<3;i++)
      { F0[0][i-1]=d1ncix[(i2elto[i][ielem])]-d1ncix[(i2elto[0][ielem])];
        F0[1][i-1]=d1nciy[(i2elto[i][ielem])]-d1nciy[(i2elto[0][ielem])];
        FX[0][i-1]=d1nccx[(i2elto[i][ielem])]-d1nccx[(i2elto[0][ielem])];
        FX[1][i-1]=d1nccy[(i2elto[i][ielem])]-d1nccy[(i2elto[0][ielem])];
        LX[0][i-1]=d1nvcx[(i2elto[i][ielem])]-d1nvcx[(i2elto[0][ielem])];
        LX[1][i-1]=d1nvcy[(i2elto[i][ielem])]-d1nvcy[(i2elto[0][ielem])];
      }

     YMATINV2(F0,F0inv,voli);
      YMATINV2(FX,FXinv,volc);

      DBL JJ=volc/voli;

assert(volc>1e-9);



if(d1J[ielem]==-1000000000){
  
  d1J[ielem]=JJ;
}
double vol_corr = sqrt(d1J[ielem] / JJ);
       d1area[ielem]=volc/2;
double F_old[2][2], F_new[2][2];

F_old[0][0] = Fe00[0][ielem]; F_old[0][1] = Fe01[0][ielem];
F_old[1][0] = Fe10[0][ielem]; F_old[1][1] = Fe11[0][ielem];
DBL F33_old = Fe33[0][ielem];




      for(i=0;i<2;i++)
      { for(j=0;j<2;j++)
        { F[i][j]=R0;
          L[i][j]=R0;
          for(k=0;k<2;k++)
          { F[i][j]=F[i][j]+FX[i][k]*F0inv[k][j];


            L[i][j]=L[i][j]+LX[i][k]*FXinv[k][j];
      } } }



          for(i=0;i<2;i++) {
        for(j=0;j<2;j++) {
  //   vol_corr=1.0;
            F[i][j] = vol_corr * F[i][j];

        }
    }
//push_Fp_expmap(L, F_old, dcstec, F_new);

//Fe00[0][ielem] =F_new[0][0] ; Fe01[0][ielem] = F_new[0][1] ;
//Fe10[0][ielem] = F_new[1][0] ; Fe11[0][ielem] =F_new[1][1]  ;
//Fe33[0][ielem] = 1.0;


      for(i=0;i<2;i++)
      { for(j=0;j<2;j++)
        { B[i][j]=R0;
          for(k=0;k<2;k++)
          { B[i][j]=B[i][j]+F[i][k]*F[j][k]; /* left Cauchy-Green strain */
          }
          D[i][j]=RP5*(L[i][j]+L[j][i]);     /* rate of deformation      */
          E[i][j]=RP5*B[i][j];               /* small strain             */
          if(i==j)E[i][j]=E[i][j]-RP5;
      } }
      if(i1usan==1) /* Apply transversely isotropic elastic constitutive law (plane stress)*/
      { T[0][0] = (d1peex / (1 - d1pemx * d1pemy)) * (E[0][0] + d1pemy * E[1][1]) + dpeks * D[0][0];
        T[1][1] = (d1peey / (1 - d1pemx * d1pemy)) * (E[1][1] + d1pemx * E[0][0]) + dpeks * D[1][1];
        T[0][1] =  (2 * d1peg) * E[0][1] + dpeks * D[0][1];
        T[1][0] =  (2 * d1peg) * E[1][0] + dpeks * D[1][0];
      }
      else if(i1usan==2) /* Apply transversely isotropic elastic constitutive law (plane strain)*/
      { d1peex_pstrain = d1peex / (1 - d1pemx * d1pemy);
        d1peey_pstrain = d1peey / (1 - d1pemy * d1pemx);
	d1pemx_pstrain = (d1pemx + d1pemx * d1pemy) / (1 - d1pemx * d1pemy);
	d1pemy_pstrain = (d1pemy + d1pemx * d1pemy) / (1 - d1pemx * d1pemy);
	
	T[0][0] = (d1peex_pstrain / (1 - d1pemx_pstrain * d1pemy_pstrain)) * (E[0][0] + d1pemy_pstrain * E[1][1]) + dpeks * D[0][0];
        T[1][1] = (d1peey_pstrain / (1 - d1pemx_pstrain * d1pemy_pstrain)) * (E[1][1] + d1pemx_pstrain * E[0][0]) + dpeks * D[1][1];
        T[0][1] =  (2 * d1peg) * E[0][1] + dpeks * D[0][1];
        T[1][0] =  (2 * d1peg) * E[1][0] + dpeks * D[1][0];
      }
      else
      { if (i1ptyp == YTE2TRIELS)
        {
          for(i=0;i<2;i++)    
          {  
              for(j=0;j<2;j++)
              { T[i][j]=R2*dpemu*E[i][j]*(voli/volc)+dpeks*D[i][j];
              }
              T[i][i]=T[i][i]+dpela*(volc/voli-voli/volc);
          }
        }
        /* Plane Stress formulation based on E,nu */
        else if (i1ptyp == YTE2PLANESTRESS)
        {        
          T[0][0] = (dpeem / (1 - dpenu*dpenu)) * (E[0][0] + dpenu * E[1][1]) + dpeks * D[0][0];
          T[1][1] = (dpeem / (1 - dpenu*dpenu)) * (E[1][1] + dpenu * E[0][0]) + dpeks * D[1][1];
          T[0][1] = (dpeem / (1 + dpenu)) * E[0][1] + dpeks * D[0][1];
          T[1][0] = (dpeem / (1 + dpenu)) * E[1][0] + dpeks * D[1][0];
        }
        /* Plane Strain formulation based on E,nu */
        else if (i1ptyp == YTE2PLANESTRAIN)
        {

           T[0][0] = dpeem/((1+dpenu)*(1-2*dpenu)) * ( (1-dpenu)*E[0][0] +  dpenu*E[1][1]) + dpeks * D[0][0];
           T[1][1] = dpeem/((1+dpenu)*(1-2*dpenu)) * ( (1-dpenu)*E[1][1] +  dpenu*E[0][0]) + dpeks * D[1][1];
           T[0][1] = (dpeem / (1+dpenu)) * E[0][1] + dpeks * D[0][1];
           T[1][0] = (dpeem / (1+dpenu)) * E[1][0] + dpeks * D[1][0];

        }
      }
      if((ncstep%1000==0)&& (ielem==0))
    //if((ncstep%100==0))
      { 

      { printf("T:%.10f,%.10f,%.10f,%.10f,%.10f,%1d\n", T[0][0], T[0][1], T[1][0], T[1][1],  ncstep);
  } 
  }
      /* output history states */
      r0x = d1ncix[(i2elto[0][ielem])];
      r0y = d1nciy[(i2elto[0][ielem])];
      r1x = d1ncix[(i2elto[1][ielem])];
      r1y = d1nciy[(i2elto[1][ielem])];
      r2x = d1ncix[(i2elto[2][ielem])];
      r2y = d1nciy[(i2elto[2][ielem])];
      for(ihys=0; ihys<nohys; ihys++)
      { 
	if(i1ohyt[ihys]==(YFLEE)) /* Strain energy */
	{
	  d1ohyt[ihys] = dctime;    /* output history time  */
          if(i1usan==1) /* Transversely isotropic elasticity */
	  { d1ohys[ihys] += 0.5 * ((E[0][0] * ((d1peex / (1 - d1pemx * d1pemy)) * (E[0][0] + d1pemy * E[1][1]))) +      /* output history state */
                                   (E[0][1] * ((1 / (2 * d1peg)) * E[0][1])) +
                                   (E[1][0] * ((1 / (2 * d1peg)) * E[0][1])) +
                                   (E[1][1] * ((d1peey / (1 - d1pemx * d1pemy)) * (E[1][1] + d1pemx * E[0][0])))) * (volc/2);    
	  }
	  else if(i1usan==2) /* Transversely isotropic elasticity (plane strain) */
	  { d1peex_pstrain = d1peex / (1 - d1pemx * d1pemy);
            d1peey_pstrain = d1peey / (1 - d1pemy * d1pemx);
	    d1pemx_pstrain = (d1pemx + d1pemx * d1pemy) / (1 - d1pemx * d1pemy);
	    d1pemy_pstrain = (d1pemy + d1pemx * d1pemy) / (1 - d1pemx * d1pemy);
	    d1ohys[ihys] += 0.5 * ((E[0][0] * ((d1peex_pstrain / (1 - d1pemx_pstrain * d1pemy_pstrain)) * (E[0][0] + d1pemy * E[1][1]))) +      /* output history state */
                                   (E[0][1] * ((1 / (2 * d1peg)) * E[0][1])) +
                                   (E[1][0] * ((1 / (2 * d1peg)) * E[0][1])) +
                                   (E[1][1] * ((d1peey_pstrain / (1 - d1pemx_pstrain * d1pemy_pstrain)) * (E[1][1] + d1pemx_pstrain * E[0][0])))) * (volc/2);    
	  }
	  else /* Isotropic elasticity */
	  { d1ohys[ihys] += 0.5 * ((E[0][0] * (2.0*dpemu*E[0][0]*(voli/volc) + dpela*(volc/voli-voli/volc))) +      /* output history state */
                                 (E[0][1] * (2.0*dpemu*E[0][1]*(voli/volc))) +
                                 (E[1][0] * (2.0*dpemu*E[1][0]*(voli/volc))) +
                                 (E[1][1] * (2.0*dpemu*E[1][1]*(voli/volc) + dpela*(volc/voli-voli/volc)))) * (volc/2);    
	  }
	}
	rpx=d1ohyx[ihys];  /* x coordinate of point P */
        rpy=d1ohyy[ihys];  /* y coordinate of point P */
        V2DCro(v0,(r1x-r0x),(r1y-r0y),(rpx-r0x),(rpy-r0y));
        V2DCro(v1,(r2x-r1x),(r2y-r1y),(rpx-r1x),(rpy-r1y));
        V2DCro(v2,(r0x-r2x),(r0y-r2y),(rpx-r2x),(rpy-r2y));

        if((v0>R0)&&(v1>R0)&&(v2>R0))    /* if point is inside the triangle */
        { if(i1ohyt[ihys]==(YFLDSXX))
          { stprev=MAXIM((EPSILON),(ABS(d1ohys[ihys])));
            if((ABS(T[0][0]-stprev))>=dohyp)
            { d1ohyt[ihys] = dctime;    /* output history time  */
              d1ohys[ihys] = T[0][0];    /* output history state */
          } }
          else if(i1ohyt[ihys]==(YFLDSXY))
          { stprev=MAXIM((EPSILON),(ABS(d1ohys[ihys])));
            if((ABS(T[0][1]-stprev))>=dohyp)  /* if((ABS(R1-T[0][1]/stprev))>dohyp) */
            { d1ohyt[ihys] = dctime;
              d1ohys[ihys] = T[0][1];
          } }
          else if(i1ohyt[ihys]==(YFLDSYY))
          { stprev=MAXIM((EPSILON),(ABS(d1ohys[ihys])));
            if((ABS(T[1][1]-stprev))>=dohyp)
            { d1ohyt[ihys] = dctime;
              d1ohys[ihys] = T[1][1];
          } }
          else if(i1ohyt[ihys]==(YFLDSZZ))
          { d1ohyt[ihys] = dctime;
            d1ohys[ihys] = R0;
          }
          else if(i1ohyt[ihys]==(YFLDSZX))
          { d1ohyt[ihys] = dctime;
            d1ohys[ihys] = R0;
          }
          else if(i1ohyt[ihys]==(YFLDSZY))
          { d1ohyt[ihys] = dctime;
            d1ohys[ihys] = R0;
          }
          else if(i1ohyt[ihys]==(YFLDVEL))
          { stprev=MAXIM((EPSILON),(ABS(d1ohys[ihys])));
            for(i=0; i<3; i++)
            { V[i] = SQRT((d1nvcx[(i2elto[i][ielem])]*d1nvcx[(i2elto[i][ielem])])
                         +(d1nvcy[(i2elto[i][ielem])]*d1nvcy[(i2elto[i][ielem])]));
            }
            if((ABS(((V[0]+V[1]+V[2])/R3)-stprev))>=dohyp)
            { d1ohyt[ihys] = dctime;
              d1ohys[ihys] = (V[0]+V[1]+V[2])/R3;  /* average velocity of element */
          } }
          else if(i1ohyt[ihys]==(YFLDVEX))
          { stprev=MAXIM((EPSILON),(ABS(d1ohys[ihys])));
            for(i=0; i<3; i++)
            { V[i] = d1nvcx[(i2elto[i][ielem])];
            }
            if((ABS(((V[0]+V[1]+V[2])/R3)-stprev))>=dohyp)
            { d1ohyt[ihys] = dctime;
              d1ohys[ihys] = (V[0]+V[1]+V[2])/R3;  /* average velocity x of element */
          } }
          else if(i1ohyt[ihys]==(YFLDVEY))
          { stprev=MAXIM((EPSILON),(ABS(d1ohys[ihys])));
            for(i=0; i<3; i++)
            { V[i] = d1nvcy[(i2elto[i][ielem])];
            }
            if((ABS(((V[0]+V[1]+V[2])/R3)-stprev))>=dohyp)
            { d1ohyt[ihys] = dctime;
              d1ohys[ihys] = (V[0]+V[1]+V[2])/R3;  /* average velocity y of element */
            } 
          } 
      } }
       
      /* Store element stress tensor only if rebars are used */ 
      if(nsbar>0)
      { d2elstr[0][ielem]=T[0][0];
        d2elstr[1][ielem]=T[0][1];
        d2elstr[2][ielem]=T[1][0]; 
        d2elstr[3][ielem]=T[1][1];
      }

      /* Nodal Forces */
      for(i=0;i<3;i++)
      { j=i+1; if(j>2)j=0;
        k=j+1; if(k>2)k=0;
        in=i2elto[i][ielem];
        jn=i2elto[j][ielem];
        kn=i2elto[k][ielem];

        nx=d1nccy[kn]-d1nccy[jn];
        ny=d1nccx[jn]-d1nccx[kn];
        d1nmct[in]=d1nmct[in]+dpero*voli/YR6;
        d1nfcx[in]=d1nfcx[in]+(T[0][0]*nx+T[0][1]*ny)/R2;
        d1nfcy[in]=d1nfcy[in]+(T[1][0]*nx+T[1][1]*ny)/R2;
	sigma00[0][ielem] = T[0][0];
  sigma10[0][ielem] = T[0][1];
  sigma11[0][ielem] = T[1][1];
  sigma33[0][ielem] = dpenu*(T[0][0]+T[1][1]);
	if (iuseis==1) /* Apply element in-situ stress */
        { Tinsitu[0][0]=-dcstxx-dcsyxx*((d1nciy[i2elto[0][ielem]]+d1nciy[i2elto[1][ielem]]+d1nciy[i2elto[2][ielem]])/3.0 - dcsrfy);
          Tinsitu[0][1]=-dcstxy-dcsyxy*((d1nciy[i2elto[0][ielem]]+d1nciy[i2elto[1][ielem]]+d1nciy[i2elto[2][ielem]])/3.0 - dcsrfy);
          Tinsitu[1][0]=-dcstxy-dcsyxy*((d1nciy[i2elto[0][ielem]]+d1nciy[i2elto[1][ielem]]+d1nciy[i2elto[2][ielem]])/3.0 - dcsrfy);
          Tinsitu[1][1]=-dcstyy-dcsyyy*((d1nciy[i2elto[0][ielem]]+d1nciy[i2elto[1][ielem]]+d1nciy[i2elto[2][ielem]])/3.0 - dcsrfy);
	  /* Nodal Forces due to in-situ stress */
	  d1nfcx[in]=d1nfcx[in]+(Tinsitu[0][0]*nx+Tinsitu[0][1]*ny)/R2;
	  d1nfcy[in]=d1nfcy[in]+(Tinsitu[1][0]*nx+Tinsitu[1][1]*ny)/R2;
	}
	
       /* Nodal Forces due to edge force*/
        jnopr=i1nopr[jn];
        knopr=i1nopr[kn];
        if( ((DABS(d1pnap[jnopr]))>EPSILON)&&
            ((DABS(d1pnap[knopr]))>EPSILON) )
        { 
                    DBL d1pnap_jnopr= d1pnap[jnopr];
          DBL d1pnap_knopr= d1pnap[knopr];

YINT numm=0.01*mcstep;
if (ncstep <= (numm)) {
    // 前10%步数：切向载荷从0 线性递增到 输入D1PNAT满值
    d1pnap_jnopr = d1pnat[jnopr] * ncstep / (numm);
    d1pnap_knopr = d1pnat[knopr] * ncstep / (numm); 

} 
          
          
          
          
          
          d1nfcx[jn]=d1nfcx[jn]-
          d1pnap_jnopr*d1pnaf[jnopr]*nx/R3-
          d1pnap_knopr*d1pnaf[knopr]*nx/YR6;
          d1nfcy[jn]=d1nfcy[jn]-
          d1pnap_jnopr*d1pnaf[jnopr]*ny/R3-
          d1pnap_knopr*d1pnaf[knopr]*ny/YR6;
          d1nfcx[kn]=d1nfcx[kn]-
          d1pnap_jnopr*d1pnaf[jnopr]*nx/YR6-
          d1pnap_knopr*d1pnaf[knopr]*nx/R3;
          d1nfcy[kn]=d1nfcy[kn]-
          d1pnap_jnopr*d1pnaf[jnopr]*ny/YR6-
          d1pnap_knopr*d1pnaf[knopr]*ny/R3;
        }
        if( ((DABS(d1pnat[jnopr]))>EPSILON)&&
            ((DABS(d1pnat[knopr]))>EPSILON) )
        { 
          DBL d1pnat_jnopr= d1pnat[jnopr];
          DBL d1pnat_knopr= d1pnat[knopr];

YINT numm=0.1*mcstep;
if (ncstep <= (numm)) {
    // 前10%步数：切向载荷从0 线性递增到 输入D1PNAT满值
    d1pnat_jnopr = d1pnat[jnopr] * ncstep / (numm);
    d1pnat_knopr = d1pnat[knopr] * ncstep / (numm); 

} 



          d1nfcx[jn]=d1nfcx[jn]-
          d1pnat_jnopr*d1pnaf[jnopr]*ny/R3-
          d1pnat_knopr*d1pnaf[knopr]*ny/YR6;
          d1nfcy[jn]=d1nfcy[jn]+
          d1pnat_jnopr*d1pnaf[jnopr]*nx/R3+
          d1pnat_knopr*d1pnaf[knopr]*nx/YR6;
          d1nfcx[kn]=d1nfcx[kn]-
          d1pnat_jnopr*d1pnaf[jnopr]*ny/YR6-
          d1pnat_knopr*d1pnaf[knopr]*ny/R3;
          d1nfcy[kn]=d1nfcy[kn]+
          d1pnat_jnopr*d1pnaf[jnopr]*nx/YR6+
          d1pnat_knopr*d1pnaf[knopr]*nx/R3;
        }
        /* Nodal forces due to absorbing boundary condition */
        if (i1ptyp == YTE2TRIELS)
        { cp=SQRT((2*dpemu+dpela)/dpero);
          cs=SQRT(dpemu/dpero);
        }
        else if (i1ptyp == YTE2PLANESTRESS)
        { lambda=(dpeem*dpenu)/((1+dpenu)*(1-2*dpenu));
          mu=(dpeem)/(2*(1+dpenu));
          cp=SQRT((2*mu+lambda)/dpero);
          cs=SQRT(mu/dpero);
        }  
        else if (i1ptyp == YTE2PLANESTRAIN)
        { lambda=(dpeem*dpenu)/((1+dpenu)*(1-2*dpenu));
          mu=(dpeem)/(2*(1+dpenu));
          cp=SQRT((2*mu+lambda)/dpero);
          cs=SQRT(mu/dpero);
        }
        if ((i1pnfx[jnopr]==4)&&(i1pnfx[knopr]==4))
        { cp=SQRT((2*dpemu+dpela)/dpero);
          cs=SQRT(dpemu/dpero);
          d1nfcx[jn]=d1nfcx[jn]-dpero*cp*d1nvcx[jn]*ABS(nx)/R3-
                     dpero*cp*d1nvcx[kn]*ABS(nx)/YR6-
                     dpero*cs*d1nvcx[jn]*ABS(ny)/R3-
                     dpero*cs*d1nvcx[kn]*ABS(ny)/YR6;
          d1nfcx[kn]=d1nfcx[kn]-dpero*cp*d1nvcx[jn]*ABS(nx)/YR6-
                     dpero*cp*d1nvcx[kn]*ABS(nx)/R3-
                     dpero*cs*d1nvcx[jn]*ABS(ny)/YR6-
                     dpero*cs*d1nvcx[kn]*ABS(ny)/R3;
        }
        if ((i1pnfy[jnopr]==4)&&(i1pnfy[knopr]==4))
        { cp=SQRT((2*dpemu+dpela)/dpero);
          cs=SQRT(dpemu/dpero);
          d1nfcy[jn]=d1nfcy[jn]-dpero*cp*d1nvcy[jn]*ABS(ny)/R3-
                     dpero*cp*d1nvcy[kn]*ABS(ny)/YR6-
                     dpero*cs*d1nvcy[jn]*ABS(nx)/R3-
                     dpero*cs*d1nvcy[kn]*ABS(nx)/YR6;
          d1nfcy[kn]=d1nfcy[kn]-dpero*cp*d1nvcy[jn]*ABS(ny)/YR6-
          dpero*cp*d1nvcy[kn]*ABS(ny)/R3-
          dpero*cs*d1nvcy[jn]*ABS(nx)/YR6-
          dpero*cs*d1nvcy[kn]*ABS(nx)/R3;
        }
      }
      //! If element is "excavated" set nodal boundary condition to v_x = 0 and v_y = 0 (hardcoded bc ID)
      //! and translate by a constant vector (hardcoded, 200)
      if(i1pexc[i1elpr[ielem]]==1)
      { i1nopr[i2elto[0][ielem]] = 1; /* 1 is a hard-coded value */
        i1nopr[i2elto[1][ielem]] = 1; /* 1 is a hard-coded value */
	i1nopr[i2elto[2][ielem]] = 1; /* 1 is a hard-coded value */
	Delta = 200;
	d1nccx[(i2elto[0][ielem])]=d1ncix[(i2elto[0][ielem])]+Delta;
	d1nccx[(i2elto[1][ielem])]=d1ncix[(i2elto[1][ielem])]+Delta;
	d1nccx[(i2elto[2][ielem])]=d1ncix[(i2elto[2][ielem])]+Delta;
	d1nccy[(i2elto[0][ielem])]=d1nciy[(i2elto[0][ielem])]+Delta;
	d1nccy[(i2elto[1][ielem])]=d1nciy[(i2elto[1][ielem])]+Delta;
	d1nccy[(i2elto[2][ielem])]=d1nciy[(i2elto[2][ielem])]+Delta;
	if(iusehf==1)
	{ //! For hydrofrac: set nodes to be non-wettable
	  i1nowe[i2elto[0][ielem]]=3;
	  i1nowe[i2elto[1][ielem]]=3;
	  i1nowe[i2elto[2][ielem]]=3;
	}
      }
    }
  }
  FREE(d1pnaf);
}



static void Yfd2JOINTS(  /* joint element */
            nelem, iprop,
            dpeft, dpegf, dpegs,
            dpeco, dpefr, dpepe,
            d1nccx,d1nccy,d1nfcx,d1nfcy,d1nvcx,
            d1nvcy,
            i1elpr,i2elto,d1elfs,
            dctime, dcstec, nebrk,netbrk,i1ebrk,d2ecbrk,d2ecbrk_NEW,
            d1etbrk,d1elbrk,d1efe,
            nesft,netsft,i1esft,d2ecsft,d1etsft,i1esftf,d1ebrkf,
            d1eike,d1edke,d1nmct,d1etmke,
            i2elnext,i2eledge,
            i1pexc,
            dusaf,dpealp,
            dpecor,dpefrrd,dpeftr,dpegfr,dpegsr,
            i1nowe,i2noid, 
            iusefn,i1edfnf,ddfnft,ddfnco,ddfngf,ddfngs,iusehf,
            i1edft,d1etike,iusesm,dctwle,
            iusehy,d2eldmg,EPS,d2dmg,d2D
            )
  YINT   nelem; YINT   iprop;
  DBL   dpeft; DBL   dpegf; DBL   dpegs;
  DBL   dpeco; DBL   dpefr; DBL   dpepe;
  DBL *d1nccx; DBL *d1nccy; DBL *d1nfcx; DBL *d1nfcy;
  DBL *d1nvcx; DBL *d1nvcy;
  YINT *i1elpr; YINT **i2elto; DBL *d1elfs;
  DBL dctime; DBL dcstec;
  YINT  *nebrk; YINT  *netbrk; YINT *i1ebrk;
  DBL  **d2ecbrk; DBL  **d2ecbrk_NEW; DBL *d1etbrk; DBL *d1elbrk; DBL *d1efe;
  YINT  *nesft; YINT  *netsft; YINT *i1esft;
  DBL  **d2ecsft; DBL *d1etsft; YINT *i1esftf; DBL *d1ebrkf;
  DBL *d1eike; DBL *d1edke; DBL *d1nmct; DBL *d1etmke;
  YINT **i2elnext; YINT **i2eledge;
  YINT *i1pexc;
  DBL dusaf; DBL dpealp; DBL dpecor; DBL dpefrrd; DBL dpeftr; DBL dpegfr; DBL dpegsr;
  YINT *i1nowe; YINT **i2noid; 
  YINT iusefn; YINT *i1edfnf; DBL ddfnft; DBL ddfnco; DBL ddfngf; DBL ddfngs; YINT iusehf;
  YINT *i1edft; DBL *d1etike; YINT iusesm; DBL dctwle;
  YINT iusehy; DBL **d2eldmg; DBL **EPS; DBL **d2dmg; DBL **d2D;
  
{ DBL dpefa=0.63;
  DBL dpefb=1.8;
  DBL dpefc=6.0;
  DBL dpefm=0.0;
  DBL small,sabs,o,s,o1,o2,s1,s2,op,sp,ot,st,dmg,z,sigma,tau;
  DBL e1x,e1y,h,area;
  YINT ielem,integ,i0,i1,i2,i3,nfail,el1,el2,el1edge,el2edge,flag;
  YINT nsoft;
  DBL dpefs;
  DBL joint_ke; /* joint kinetic energy (i.e. kinetic energy of the four joint nodes) */
  DBL delta_ke; /* difference between joint_ke and the joint kinetic energy calculated when it first yields */
  
  DBL beta;  /* joint orientation angle (0-180) */
  DBL gamma; /* relative angle between layering and joint element (0-90) */
  
  //! Assigning strength input values to "peak" values
  DBL dpeftp = dpeft;
  DBL dpegfp = dpegf;
  DBL dpegsp = dpegs;
  DBL C; /* exponent for power-law function */
  
  //DBL Fxi0;
  //DBL Fyi0;
  
  /*static FILE *out1=FILENULL;
  if(out1 == FILENULL)
  { out1=fopen("Yfd2JOINTS.txt", "a");
  }*/
  









  small=EPSILON; 
  for(ielem=0;ielem<nelem;ielem++)
  { if(i1elpr[ielem]==iprop)
    {// if((ielem%75)==0){i1elpr[ielem]=iprop-YIPROPMAX;}
      dpefs=d1elfs[ielem];
      i0=i2elto[0][ielem];
      i1=i2elto[1][ielem];
      i2=i2elto[2][ielem];
      i3=i2elto[3][ielem];
    
      
      e1x=RP5*(d1nccx[i1]+d1nccx[i2]-d1nccx[i0]-d1nccx[i3]);
      e1y=RP5*(d1nccy[i1]+d1nccy[i2]-d1nccy[i0]-d1nccy[i3]);
      h=SQRT(e1x*e1x+e1y*e1y);
      e1x=e1x/(h+small);
      e1y=e1y/(h+small);
      s1=(d1nccy[i0]-d1nccy[i3])*e1y+(d1nccx[i0]-d1nccx[i3])*e1x;
      s2=(d1nccy[i1]-d1nccy[i2])*e1y+(d1nccx[i1]-d1nccx[i2])*e1x;
      o1=(d1nccy[i0]-d1nccy[i3])*e1x-(d1nccx[i0]-d1nccx[i3])*e1y;
      o2=(d1nccy[i1]-d1nccy[i2])*e1x-(d1nccx[i1]-d1nccx[i2])*e1y;
      //! Anisotropic fracture model
      if(dusaf>0.0)
      {
        //! Calculation of joint orientation (i.e. angle beta)
        //beta = atan( (((d1nccy[i1]+d1nccy[i2])/2)-((d1nccy[i0]+d1nccy[i3])/2))/(((d1nccx[i1]+d1nccx[i2])/2)-((d1nccx[i0]+d1nccx[i3])/2))) * 180/MYPI;
        beta = atan((d1nccy[i1]-d1nccy[i0])/(d1nccx[i1]-d1nccx[i0])) * 180/MYPI; // simplified formula
        if (beta < 0.0)
        { beta = beta + 180.0; }
      
        //! Calculation of relative angle between layering and joint element
        gamma = ABS(dpealp - beta);
        if(gamma > 90.0)
        { 
          gamma = 180.0 - gamma; 
        }
        
        //! Power-law variation with exponent C = IUSAF
        C = dusaf;
        dpeft = dpeftr + (dpeftp-dpeftr) * pow((gamma/90.0),C);
        dpegf = dpegfr + (dpegfp-dpegfr) * pow((gamma/90.0),C);
        dpegs = dpegsr + (dpegsp-dpegsr) * pow((gamma/90.0),C);
        // dpefs is updated below according to Mohr-Coulomb */        

      }
      
      op=R2*h*dpeft/dpepe;
      sp=R2*h*dpefs/dpepe;
      ot=MAXIM(EPSILON,(R3*dpegf/dpeft));
      st=MAXIM(EPSILON,(R3*dpegs/dpefs));
      //ot=MAXIM((R2*op),(R3*dpegf/dpeft));
      //st=MAXIM((R2*sp),(R3*dpegs/dpefs));
      
      //! Use "cohesive" DFN properties to calculate op, sp, ot, and st
      if((iusefn == 2) && (i1edfnf[ielem]==1))
      { //! If joint element belongs to DFN
        op=R2*h*ddfnft/dpepe;
        sp=R2*h*dpefs/dpepe;
        ot=MAXIM(EPSILON,(R3*ddfngf/ddfnft));
        st=MAXIM(EPSILON,(R3*ddfngs/dpefs));
        // dpefs is updated below according to Mohr-Coulomb 
      }
      
      nfail=0;
      nsoft=0;
      
      el1=i2elnext[0][ielem];     /* 1st element next to the joint */
      el2=i2elnext[1][ielem];     /* 2nd element next to the joint */
      el1edge=i2elnext[2][ielem]; /* edge number (0,1,2) of the 1st element (el1) */
      el2edge=i2elnext[3][ielem]; /* edge number (0,1,2) of the 2nd element (el2) */
      if((el1==-1)||(el2==-1)) continue; /* no need to do further computations if this is an external edge */
      //! Applying mixed DFN (i.e., DFN type 3) using the flag assigned to the joint element
      if(iusefn == 3)
      { if(i1edft[ielem] == 1) //! Broken-type DFN crack
      	{ i1elpr[ielem]=iprop-YIPROPMAX;
          d1ebrkf[ielem]=5.0;
          //! For hydrofrac
          i2noid[0][i0]=i1;
          i2noid[1][i1]=i0;
          i2noid[0][i2]=i3;
          i2noid[1][i3]=i2;
        }
        if(i1edft[ielem] == 2) //! Cohesive-type DFN crack
        { op=R2*h*ddfnft/dpepe;
          sp=R2*h*dpefs/dpepe;
          ot=MAXIM(EPSILON,(R3*ddfngf/ddfnft));
          st=MAXIM(EPSILON,(R3*ddfngs/dpefs));
          // dpefs is updated below according to Mohr-Coulomb
        }
      }
         
      // 体单元塑性状态（外部每步更新）
     // 每个体单元的等效塑性应变（或塑性耗散能）

// 门控控制参数
DBL epeq_crit=0.001;         // 开始触发界面软化的塑性阈值   取材料屈服后、但还没明显局部化时的等效塑性应变。常见范围 0.0005–0.002
DBL epeq_band=0.2*epeq_crit;         // 平滑过渡带宽 Δεp 控制 α 从 0 增长到 1 的“带宽”，越大过渡越缓，初值可取 epeq_band ≈ (0.5–1.0) × epeq_crit，这样 α 在 epeq_crit 附近有一个平滑的 S 型增长。如果想让界面延迟更久，取大一些
DBL epeqL = EPS[0][el1];
DBL epeqR = EPS[0][el2];
DBL epeq_edge = MAXIM(epeqL, epeqR); // 或平均

// --- smooth plasticity weight: w ∈ [0,1], S-shaped unlock around epeq_crit ---
DBL w_arg = (epeq_edge - epeq_crit) / (epeq_band + EPSILON);
DBL w = 0.5*(1.0 + tanh(w_arg));       // soft gating (α)



// --- damage evolution parameters (new) ---
DBL k_rate   = 1.0;                  // global scale; calibrate to match Gc   初值1.0
DBL p_exp    = 1.0;                    // displacement sensitivity (1–2)
DBL q_exp    = 3.0;                    // plasticity sensitivity (2–4)
DBL k_floor  = 1e-3*k_rate;            // tiny floor once beyond peak; 1e-4–1e-2*k_rate

      //! Performing excavation: set joint element state to broken if between at least one "excavated" element
      if((i1pexc[i1elpr[el1]]==1)||(i1pexc[i1elpr[el2]]==1))
      { i1elpr[ielem]=iprop-YIPROPMAX; 
        d1ebrkf[ielem]=4.0;
        if(iusehf==1)
	{ //! For hydrofrac
          i2noid[0][i0]=i1;
          i2noid[1][i1]=i0;
          i2noid[0][i2]=i3;
          i2noid[1][i3]=i2;
        }
      }
      else 
      {
      //Fxi0=0.0;
      //Fyi0=0.0;
      for(integ=0;integ<3;integ++)
      {
        
        flag=0;
        if(integ==0)
        { o=o1; s=s1;
        }
        else if(integ==2)
        { o=o2; s=s2;
        }
        else
        { o=RP5*(o1+o2); s=RP5*(s1+s2);
        }
        sabs=ABS(s);


if ((o > op) && (sabs > sp)) {
  dmg=SQRT(((o-op)/ot)*((o-op)/ot)+((sabs-sp)/st)*((sabs-sp)/st));
} else if (o > op) {
  dmg=(o-op)/ot;
} else if (sabs > sp) {
  dmg=(sabs-sp)/st;
}
else
        { dmg=R0;
        }

          if((integ==0) && (dmg>=d2dmg[0][ielem])) 
          {  flag=1;
            d2dmg[0][ielem]=dmg;
            }
          else if((integ==1) && (dmg>=d2dmg[1][ielem])) 
          { flag=1;
            d2dmg[1][ielem]=dmg;
            }
          else if((integ==2) && (dmg>=d2dmg[2][ielem])) 
          { flag=1;
            d2dmg[2][ielem]=dmg;
            }




// --- equivalent separation for geometric drive (choose form consistent with your law) ---
    // Option A: use projections (o for normal opening; sabs for shear), weighted:
    DBL alpha_n = 1.0, alpha_s = 1.0;   // tune if anisotropic mix needed
    DBL del_bar = SQRT(alpha_n*o*o + alpha_s*sabs*sabs);
    ////DBL del_c   = MAXIM(EPSILON, MAXIM(op, sp));   // scale for non-dimensionalization

    //// --- geometric and plastic drive terms for rate ---
    ////DBL geom = pow(MAXIM(EPSILON, del_bar/del_c), p_exp);
    DBL del_bar_n;
if (sabs > o) {del_bar_n = del_bar / (sp + EPSILON);} else {del_bar_n = del_bar / (op + EPSILON);}
DBL geom = pow(del_bar_n, p_exp);

    DBL plast = pow(w, q_exp);

    // --- previous damage state: reuse d2eldmg as "peak" storage; add D store if you prefer clarity ---
    DBL D_prev = d2D[integ][ielem];  // previously used for dmg peak; reinterpret as D peak

    // --- damage rate with soft gating ---
    DBL D_dot = k_rate * geom * plast;

    // --- tiny floor beyond peak to avoid infinite platform (still geometry-controlled) ---
    if (((o >= op) || (sabs >= sp))) {
        D_dot = MAXIM(D_dot, k_floor * geom);
    }
     if (flag==0) { /* unloading branch */
        D_dot = k_rate * geom * plast * (-1.0);
        //D_dot = 0.0;


    }
    // --- integrate and clamp ---
    DBL D = D_prev + D_dot * dcstec;
    if (fabs(1.0 - D) < 1.0e-2)  D = 1.0;
    if (D < 0.0) D = 0.0;
    if (D > 1.0){
D = 1.0;


    } 
// --- hysteresis peak tracking (loading branch) ---
    if (D >= d2eldmg[integ][ielem]) {
        d2eldmg[integ][ielem] = D;
    }

     
       d2D[integ][ielem] = D;



    // --- compute z from D (replace dmg→D everywhere) ---
    if(iusehy==0) {
        z = (R1 - ((dpefa+dpefb-R1)/(dpefa+dpefb)) *
            exp(D*(dpefa+dpefc*dpefb)/((dpefa+dpefb)*(R1-dpefa-dpefb))))
            * (dpefa*(R1-D) + dpefb*pow((R1-D),dpefc));
    } else {
        // Loading branch (peak update already done above)
        DBL Dpeak = d2eldmg[integ][ielem];
        z = (R1 - ((dpefa+dpefb-R1)/(dpefa+dpefb)) *
            exp(Dpeak*(dpefa+dpefc*dpefb)/((dpefa+dpefb)*(R1-dpefa-dpefb))))
            * (dpefa*(R1-Dpeak) + dpefb*pow((R1-Dpeak), dpefc));

        // Unloading branch (keep your linear scaling idea, swap dmg→D)
        if (D < Dpeak) {
            //DBL scale = (op + D*ot) / (op + Dpeak*ot); // mirrors your ((op + dmg*ot)/(op + dmg_peak*ot))
           DBL scale =  D /  Dpeak;
           
            z *= scale;
        }
    }




    /* 8) 统计“软化检测计数”（保留你的 nsoft 用法：超过一定次数登记软化事件）
      原来是 dmg>0 的三种分支里递增；现在用“几何已经跨峰值”判断递增，
      并保留你登记 i1esft[], d1etsft[], KE等的逻辑。 */
/* 软化检测：几何已跨峰值（统一使用 op, sp） */
if ((o > op) && (sabs > sp)) {

    nsoft = nsoft + 1;
    if ((nsoft > 2) && (i1esftf[ielem] == 0)) {
        i1esft[*nesft] = ielem;
        d1etsft[*nesft] = dctime;
        d1etmke[ielem] = dctime;
        d2ecsft[0][*nesft] = (d1nccx[i0]+d1nccx[i1]+d1nccx[i2]+d1nccx[i3])/R4;
        d2ecsft[1][*nesft] = (d1nccy[i0]+d1nccy[i1]+d1nccy[i2]+d1nccy[i3])/R4;
        d1eike[ielem] = 0.5*(d1nmct[i0]*(d1nvcx[i0]*d1nvcx[i0]+d1nvcy[i0]*d1nvcy[i0])+
                             d1nmct[i1]*(d1nvcx[i1]*d1nvcx[i1]+d1nvcy[i1]*d1nvcy[i1])+
                             d1nmct[i2]*(d1nvcx[i2]*d1nvcx[i2]+d1nvcy[i2]*d1nvcy[i2])+
                             d1nmct[i3]*(d1nvcx[i3]*d1nvcx[i3]+d1nvcy[i3]*d1nvcy[i3]));
        d1etike[ielem] = dctime;
        i1esftf[ielem] = 1;
        (*nesft)++;
        (*netsft)++;
    }
} else if (o > op) {
  
    nsoft = nsoft + 1;
    if ((nsoft > 2) && (i1esftf[ielem] == 0)) {
        i1esft[*nesft] = ielem;
        d1etsft[*nesft] = dctime;
        d1etmke[ielem] = dctime;
        d2ecsft[0][*nesft] = (d1nccx[i0]+d1nccx[i1]+d1nccx[i2]+d1nccx[i3])/R4;
        d2ecsft[1][*nesft] = (d1nccy[i0]+d1nccy[i1]+d1nccy[i2]+d1nccy[i3])/R4;
        d1eike[ielem] = 0.5*(d1nmct[i0]*(d1nvcx[i0]*d1nvcx[i0]+d1nvcy[i0]*d1nvcy[i0])+
                             d1nmct[i1]*(d1nvcx[i1]*d1nvcx[i1]+d1nvcy[i1]*d1nvcy[i1])+
                             d1nmct[i2]*(d1nvcx[i2]*d1nvcx[i2]+d1nvcy[i2]*d1nvcy[i2])+
                             d1nmct[i3]*(d1nvcx[i3]*d1nvcx[i3]+d1nvcy[i3]*d1nvcy[i3]));
        d1etike[ielem] = dctime;
        i1esftf[ielem] = 2;
        (*nesft)++;
        (*netsft)++;
    }
} else if (sabs > sp) {
  
    nsoft = nsoft + 1;
    if ((nsoft > 2) && (i1esftf[ielem] == 0)) {
        i1esft[*nesft] = ielem;
        d1etsft[*nesft] = dctime;
        d1etmke[ielem] = dctime;
        d2ecsft[0][*nesft] = (d1nccx[i0]+d1nccx[i1]+d1nccx[i2]+d1nccx[i3])/R4;
        d2ecsft[1][*nesft] = (d1nccy[i0]+d1nccy[i1]+d1nccy[i2]+d1nccy[i3])/R4;
        d1eike[ielem] = 0.5*(d1nmct[i0]*(d1nvcx[i0]*d1nvcx[i0]+d1nvcy[i0]*d1nvcy[i0])+
                             d1nmct[i1]*(d1nvcx[i1]*d1nvcx[i1]+d1nvcy[i1]*d1nvcy[i1])+
                             d1nmct[i2]*(d1nvcx[i2]*d1nvcx[i2]+d1nvcy[i2]*d1nvcy[i2])+
                             d1nmct[i3]*(d1nvcx[i3]*d1nvcx[i3]+d1nvcy[i3]*d1nvcy[i3]));
        d1etike[ielem] = dctime;
        i1esftf[ielem] = 3;
        (*nesft)++;
        (*netsft)++;
    }
}



/* 9) 断裂判据：用 D≥1 代替原 dmg≥1；其余登记逻辑保留 */
if (D >= R1) /* joint element broken */ 
{
    nfail = nfail + 1;
    if ((nfail > 1) && (i1elpr[ielem] >= 0)) {
        i1elpr[ielem] = iprop - YIPROPMAX;
        i1ebrk[*nebrk] = ielem;
        d1etbrk[*nebrk] = dctime;

        /* Mode 分类仍用你的 o/sabs 与残余门槛判据 */
        if ((o >= (op + ot)) && (sabs >= (sp + st))) {       // Mode I + II
            d1ebrkf[ielem] = 3.0;
        } else if (o >= (op + ot)) {                         // Mode I
            d1ebrkf[ielem] = 1.0;
        } else if (sabs >= (sp + st)) {                      // Mode II
            d1ebrkf[ielem] = 2.0;
        } else {                                             // 混合模式（你的向量式）
            d1ebrkf[ielem] = 1.0 + (sabs - sp) / (st);
        }

        /* 断裂位置与长度登记（保留） */
        d2ecbrk[0][*nebrk]=(d1nccx[i0]+d1nccx[i1]+d1nccx[i2]+d1nccx[i3])/R4;
        d2ecbrk[1][*nebrk]=(d1nccy[i0]+d1nccy[i1]+d1nccy[i2]+d1nccy[i3])/R4;
        d2ecbrk_NEW[0][ielem]=(d1nccx[i0]+d1nccx[i1]+d1nccx[i2]+d1nccx[i3])/R4;
        d2ecbrk_NEW[1][ielem]=(d1nccy[i0]+d1nccy[i1]+d1nccy[i2]+d1nccy[i3])/R4;
        d1elbrk[*nebrk]=SQRT((d1nccx[i0]-d1nccx[i1])*(d1nccx[i0]-d1nccx[i1])+(d1nccy[i0]-d1nccy[i1])*(d1nccy[i0]-d1nccy[i1]));
        (*nebrk)++;
        (*netbrk)++;

        //! For hydrofrac
        i2noid[0][i0]=i1;
        i2noid[1][i1]=i0;
        i2noid[0][i2]=i3;
        i2noid[1][i3]=i2;
    }
    /* 将当前 D 钳到 1（与原 dmg=R1 一致） */
    D = R1;
}

    // --- intact/broken edges logic unchanged, just use D instead of dmg for the condition ---
    if (i0!=i3 && i1!=i2) {
        if (D >= R1) {
            i2eledge[el1edge][el1] = -1;
            i2eledge[el2edge][el2] = -1;
        } else {
            i2eledge[el1edge][el1] = 1;
            i2eledge[el2edge][el2] = 1;
        }
    }






  /* --- 法向应力 σ 更新 --- */
if (o < R0) {  
    // 压缩状态：线性接触刚度
    sigma = R2 * o * dpeft / op;   // op_eff = op*(1-κw)，若不用微前移则直接 op
} 
else if (o > op) {  
    // 超过峰值：由退化因子 z(D) 控制
    sigma = dpeft * z;  
} 
else {  
    // 峰值前：二次多项式上升段
    DBL ratio = o / op;
    sigma = (R2*ratio - ratio*ratio) * z * dpeft;
}

/* --- Mohr-Coulomb 剪切强度更新 --- */
if (dpeco > R0) {
    if (sigma > R0) {
        // 拉开状态
        if (dusaf > 0) { // 各向异性修正
            C = dusaf;
            dpefs = dpecor + (dpeco - dpecor) * pow((gamma/90.0), C);
        } else {
            dpefs = dpeco;
        }
        // DFN 裂缝覆盖
        if (((iusefn == 2) && (i1edfnf[ielem] == 1)) || ((iusefn == 3) && (i1edft[ielem] == 2))) {
            dpefs = ddfnco;
        }
    } else {
        // 压缩状态
        if (dusaf > 0) {
            C = dusaf;
            dpefs = dpecor + (dpeco - dpecor) * pow((gamma/90.0), C)
                  - sigma * (dpefrrd + (dpefr - dpefrrd) * pow((gamma/90.0), C));
        } else {
            dpefs = dpeco - sigma * dpefr;
        }
        // DFN 裂缝覆盖
        if (((iusefn == 2) && (i1edfnf[ielem] == 1)) || ((iusefn == 3) && (i1edft[ielem] == 2))) {
            dpefs = ddfnco - sigma * dpefr;
        }
    }
}

/* --- 剪切应力 τ 更新 --- */
if ((sigma > R0) && (sabs > sp)) {
    tau = z * dpefs;
} 
else if (sigma > R0) {
    DBL ratio = sabs / sp;
    tau = (R2*ratio - ratio*ratio) * z * dpefs;
} 
else if (sabs > sp) {
    tau = z * dpefs - dpefm * sigma;
} 
else {
    DBL ratio = sabs / sp;
    tau = (R2*ratio - ratio*ratio) * (z * dpefs - dpefm * sigma);
}

/* 更新 Mohr-Coulomb 剪切强度 */
d1elfs[ielem] = dpefs;

/* 剪切方向修正 */
if (s < R0) tau = -tau;


    // --- energy accounting (optional but recommended) ---
    // Ed += (sigma*delta_n_inc + tau*delta_s_inc) - elastic_storage_change;
    // Ep += solid_plastic_work_increment_from_el1_el2();


        

  
        if(integ==0)  /* nodal forces */
        { area=h/YR6; /* area=h/6.0; */
          d1nfcx[i0]=d1nfcx[i0]-area*(tau*e1x-sigma*e1y);
          d1nfcx[i3]=d1nfcx[i3]+area*(tau*e1x-sigma*e1y);
          d1nfcy[i0]=d1nfcy[i0]-area*(tau*e1y+sigma*e1x);
          d1nfcy[i3]=d1nfcy[i3]+area*(tau*e1y+sigma*e1x);
          /* Element fracture energy */
          d1efe[ielem]=d1efe[ielem]+dcstec*((-area*(tau*e1x-sigma*e1y))*d1nvcx[i0]+(-area*(tau*e1y+sigma*e1x))*d1nvcy[i0]+
                        (area*(tau*e1x-sigma*e1y))*d1nvcx[i3]+(area*(tau*e1y+sigma*e1x))*d1nvcy[i3]);
          
          /*if(ielem==3)
          { Fxi0=Fxi0-area*(tau*e1x-sigma*e1y);
            Fyi0=Fyi0-area*(tau*e1y+sigma*e1x); }*/
        }
        else if(integ==1)
        { area=h/R3;  /* area=h/3.0; */
          d1nfcx[i0]=d1nfcx[i0]-area*(tau*e1x-sigma*e1y);
          d1nfcx[i3]=d1nfcx[i3]+area*(tau*e1x-sigma*e1y);
          d1nfcy[i0]=d1nfcy[i0]-area*(tau*e1y+sigma*e1x);
          d1nfcy[i3]=d1nfcy[i3]+area*(tau*e1y+sigma*e1x);
          d1nfcx[i1]=d1nfcx[i1]-area*(tau*e1x-sigma*e1y);
          d1nfcx[i2]=d1nfcx[i2]+area*(tau*e1x-sigma*e1y);
          d1nfcy[i1]=d1nfcy[i1]-area*(tau*e1y+sigma*e1x);
          d1nfcy[i2]=d1nfcy[i2]+area*(tau*e1y+sigma*e1x);
          /* Element fracture energy */
          d1efe[ielem]=d1efe[ielem]+dcstec*((-area*(tau*e1x-sigma*e1y))*d1nvcx[i0]+(-area*(tau*e1y+sigma*e1x))*d1nvcy[i0]+
                       (-area*(tau*e1x-sigma*e1y))*d1nvcx[i1]+(-area*(tau*e1y+sigma*e1x))*d1nvcy[i1]+
                       (+area*(tau*e1x-sigma*e1y))*d1nvcx[i2]+(+area*(tau*e1y+sigma*e1x))*d1nvcy[i2]+
                       (+area*(tau*e1x-sigma*e1y))*d1nvcx[i3]+(+area*(tau*e1y+sigma*e1x))*d1nvcy[i3]);

          /*if(ielem==3)
          { Fxi0=Fxi0-area*(tau*e1x-sigma*e1y);  
            Fyi0=Fyi0-area*(tau*e1y+sigma*e1x); }*/
          //if(ielem==3)
          //{ fprintf(out1,"%.6f \t %.6f \t %.6f \n",dctime,sigma,tau); }
        }
        else
        { area=h/YR6; /* area=h/6.0; */
          d1nfcx[i1]=d1nfcx[i1]-area*(tau*e1x-sigma*e1y);
          d1nfcx[i2]=d1nfcx[i2]+area*(tau*e1x-sigma*e1y);
          d1nfcy[i1]=d1nfcy[i1]-area*(tau*e1y+sigma*e1x);
          d1nfcy[i2]=d1nfcy[i2]+area*(tau*e1y+sigma*e1x);
          /* Element fracture energy */
          d1efe[ielem]=d1efe[ielem]+dcstec*((-area*(tau*e1x-sigma*e1y))*d1nvcx[i1]+(-area*(tau*e1y+sigma*e1x))*d1nvcy[i1]+
                       (+area*(tau*e1x-sigma*e1y))*d1nvcx[i2]+(+area*(tau*e1y+sigma*e1x))*d1nvcy[i2]);
        }
      }
      }
      /* If joint is yielded compute differential kinetic energy */
      if(i1esftf[ielem]>0)
      { joint_ke = 0.5*(d1nmct[i0]*(d1nvcx[i0]*d1nvcx[i0]+d1nvcy[i0]*d1nvcy[i0])+
                        d1nmct[i1]*(d1nvcx[i1]*d1nvcx[i1]+d1nvcy[i1]*d1nvcy[i1])+
                        d1nmct[i2]*(d1nvcx[i2]*d1nvcx[i2]+d1nvcy[i2]*d1nvcy[i2])+
                        d1nmct[i3]*(d1nvcx[i3]*d1nvcx[i3]+d1nvcy[i3]*d1nvcy[i3]));
        delta_ke = joint_ke - d1eike[ielem];
        /* If the differential kinetic energy is greater than that from previous timestep update the joint differential kin energy */
        if(delta_ke > d1edke[ielem])
        { d1edke[ielem] = delta_ke; 
          d1etmke[ielem] = dctime;
        }
        if(iusesm == 1) //! Use maximum time window duration
        { if (dctime >= d1etike[ielem] + dctwle)
          { d1eike[ielem] = joint_ke;
            d1etike[ielem] = dctime; 
        } }
      }
      if(iusesm == 2) //! Record energy at the time of failure
      { if(d1ebrkf[ielem]>0)
        { joint_ke = 0.5*(d1nmct[i0]*(d1nvcx[i0]*d1nvcx[i0]+d1nvcy[i0]*d1nvcy[i0])+
                          d1nmct[i1]*(d1nvcx[i1]*d1nvcx[i1]+d1nvcy[i1]*d1nvcy[i1])+
                          d1nmct[i2]*(d1nvcx[i2]*d1nvcx[i2]+d1nvcy[i2]*d1nvcy[i2])+
                          d1nmct[i3]*(d1nvcx[i3]*d1nvcx[i3]+d1nvcy[i3]*d1nvcy[i3]));
          d1edke[ielem] = joint_ke;
          d1etmke[ielem] = dctime;
      } }
    }
  }
}

static void Yfd2JOINTS_ECZM(  /* joint element */
  YINT   nelem, YINT   iprop,
  DBL   dpeft, DBL   dpegf, DBL   dpegs,
  DBL   dpeco, DBL   dpefr, DBL   dpepe,
  DBL *d1nccx, DBL *d1nccy, DBL *d1nfcx, DBL *d1nfcy,
  DBL *d1nvcx, DBL *d1nvcy,
  YINT *i1elpr, YINT **i2elto, DBL *d1elfs,
  DBL dctime, DBL dcstec,
  YINT  *nebrk, YINT  *netbrk, YINT *i1ebrk,
  DBL  **d2ecbrk, DBL **d2ecbrk_NEW, DBL *d1etbrk, DBL *d1elbrk, DBL *d1efe,
  YINT  *nesft, YINT  *netsft, YINT *i1esft,
  DBL **d2ecsft, DBL *d1etsft, YINT *i1esftf, DBL *d1ebrkf,
  DBL *d1eike, DBL *d1edke, DBL *d1nmct, DBL *d1etmke,
  YINT **i2elnext, YINT **i2eledge,
  YINT *i1pexc,
  DBL dusaf, DBL dpealp, DBL dpecor, DBL dpefrrd, DBL dpeftr, DBL dpegfr, DBL dpegsr,
  YINT *i1nowe, YINT **i2noid,
  YINT iusefn, YINT *i1edfnf, DBL ddfnft, DBL ddfnco, DBL ddfngf, DBL ddfngs, YINT iusehf,
  YINT *i1edft, DBL *d1etike, YINT iusesm, DBL dctwle,
  YINT iusehy, DBL **d2eldmg,

   DBL **sigma00, DBL **sigma10, DBL **sigma11, DBL *d1area,DBL  ncstep,DBL *d1inst,DBL *d1inss,YINT *i1elprtmp,YINT *i1elfr ,YINT *iuptrimesh ,YINT *i1remeshf
  ){
   DBL dpefa=0.63;
  DBL dpefb=1.8;
  DBL dpefc=6.0;
  DBL dpefm=0.0;
  DBL sabs;
  DBL o;
  DBL s;
  DBL ot;
  DBL st;
  DBL dmg;
  DBL z;
  DBL sigma=0.0;
  DBL tau=0.0;
  DBL sp;
  DBL op;
  DBL s2;
  DBL s1;
  DBL o1;
  DBL o2;
  DBL small=EPSILON;
   DBL e1x,e1y,h,area;
  YINT ielem,integ,i0,i1,i2,i3,nfail,el1,el2,el1edge,el2edge;
  YINT nsoft,i,j;
  DBL dpefs;
  DBL joint_ke; /* joint kinetic energy (i.e. kinetic energy of the four joint nodes) */
  DBL delta_ke; /* difference between joint_ke and the joint kinetic energy calculated when it first yields */
  
  DBL beta;  /* joint orientation angle (0-180) */
  DBL gamma; /* relative angle between layering and joint element (0-90) */
  DBL Tedge[2][2],TT0[2][2],TT1[2][2];
  //! Assigning strength input values to "peak" values
  DBL dpeftp = dpeft;
  DBL dpegfp = dpegf;
  DBL dpegsp = dpegs;
  DBL C; /* exponent for power-law function */

  small=EPSILON; 
  for(ielem=0;ielem<nelem;ielem++)
  { 
    
    
    if((i1elpr[ielem]==iprop))
    {// if((ielem%75)==0){i1elpr[ielem]=iprop-YIPROPMAX;}
      dpefs=d1elfs[ielem];
      i0=i2elto[0][ielem];
      i1=i2elto[1][ielem];
      i2=i2elto[2][ielem];
      i3=i2elto[3][ielem];     

      el1=i2elnext[0][ielem];     /* 1st element next to the joint */
      el2=i2elnext[1][ielem];     /* 2nd element next to the joint */
      el1edge=i2elnext[2][ielem]; /* edge number (0,1,2) of the 1st element (el1) */
      el2edge=i2elnext[3][ielem]; /* edge number (0,1,2) of the 2nd element (el2) */
      
      if((el1==-1)||(el2==-1)) continue; /* no need to do further computations if this is an external edge */
             e1x=RP5*(d1nccx[i1]+d1nccx[i2]-d1nccx[i0]-d1nccx[i3]);
      e1y=RP5*(d1nccy[i1]+d1nccy[i2]-d1nccy[i0]-d1nccy[i3]);
      h=SQRT(e1x*e1x+e1y*e1y);
      e1x=e1x/(h+small);//qiexiang
      e1y=e1y/(h+small);
      s1=(d1nccy[i0]-d1nccy[i3])*e1y+(d1nccx[i0]-d1nccx[i3])*e1x;
      s2=(d1nccy[i1]-d1nccy[i2])*e1y+(d1nccx[i1]-d1nccx[i2])*e1x;
      o1=(d1nccy[i0]-d1nccy[i3])*e1x-(d1nccx[i0]-d1nccx[i3])*e1y;
      o2=(d1nccy[i1]-d1nccy[i2])*e1x-(d1nccx[i1]-d1nccx[i2])*e1y;
      //! Anisotropic fracture model
      if(dusaf>0.0)
      {
        //! Calculation of joint orientation (i.e. angle beta)
        //beta = atan( (((d1nccy[i1]+d1nccy[i2])/2)-((d1nccy[i0]+d1nccy[i3])/2))/(((d1nccx[i1]+d1nccx[i2])/2)-((d1nccx[i0]+d1nccx[i3])/2))) * 180/MYPI;
        beta = atan((d1nccy[i1]-d1nccy[i0])/(d1nccx[i1]-d1nccx[i0])) * 180/MYPI; // simplified formula
        if (beta < 0.0)
        { beta = beta + 180.0; }
      
        //! Calculation of relative angle between layering and joint element
        gamma = ABS(dpealp - beta);
        if(gamma > 90.0)
        { 
          gamma = 180.0 - gamma; 
        }
        
        //! Power-law variation with exponent C = IUSAF
        C = dusaf;
        dpeft = dpeftr + (dpeftp-dpeftr) * pow((gamma/90.0),C);
        dpegf = dpegfr + (dpegfp-dpegfr) * pow((gamma/90.0),C);
        dpegs = dpegsr + (dpegsp-dpegsr) * pow((gamma/90.0),C);
      }
      op=R2*h*dpeft/dpepe;
      sp=R2*h*dpefs/dpepe;
      ot=MAXIM(EPSILON,(R3*dpegf/dpeft));
      st=MAXIM(EPSILON,(R3*dpegs/dpefs));
      //ot=MAXIM((R2*op),(R3*dpegf/dpeft));
      //st=MAXIM((R2*sp),(R3*dpegs/dpefs));
      
      //! Use "cohesive" DFN properties to calculate op, sp, ot, and st
      if((iusefn == 2) && (i1edfnf[ielem]==1))
      { //! If joint element belongs to DFN
        op=R2*h*ddfnft/dpepe;
        sp=R2*h*dpefs/dpepe;
        ot=MAXIM(EPSILON,(R3*ddfngf/ddfnft));
        st=MAXIM(EPSILON,(R3*ddfngs/dpefs));
        // dpefs is updated below according to Mohr-Coulomb 
      }
      
      nfail=0;
      nsoft=0;
                TT0[0][0] = sigma00[0][el1];
                 TT0[1][1] =sigma11[0][el1];
                 TT0[0][1] = sigma10[0][el1];
                 TT0[1][0] = sigma10[0][el1];
                TT1[0][0] = sigma00[0][el2];
                 TT1[1][1] =sigma11[0][el2];
                 TT1[0][1] = sigma10[0][el2];
                 TT1[1][0] = sigma10[0][el2];

      for(i=0;i<2;i++)   
      { for(j=0;j<2;j++)   
        { Tedge[i][j]=0.0;
          Tedge[i][j]=((d1area[el1]*TT0[i][j])+(d1area[el2]*TT1[i][j]))/(d1area[el1]+d1area[el2]);
        }
      }
          
          // 法向应力 σ_n = n · σ · n
          DBL   sigma_n = (-e1y) * (Tedge[0][0] * (-e1y)+ Tedge[0][1] * e1x) 
                        + e1x * (Tedge[1][0] * (-e1y) + Tedge[1][1] * e1x);

          // 切向应力 τ_s = t · σ · n
          DBL  tau_s = e1x * (Tedge[0][0] * (-e1y) + Tedge[0][1] * e1x) 
                     + e1y * (Tedge[1][0] * (-e1y) + Tedge[1][1] * e1x);

if(((sigma_n >= dpeft)||(fabs(tau_s) >= dpefs))||(i1elprtmp[ielem]==1111)||(i1elfr[ielem]==1))
    {
                        if(i1elprtmp[ielem]<0){
                        i1remeshf[i2elto[0][ielem]]=1;
                        i1remeshf[i2elto[1][ielem]]=1;
                        i1remeshf[i2elto[2][ielem]]=1;
                        i1remeshf[i2elto[3][ielem]]=1;
                        

                        d1inst[ielem]=sigma_n;
                        d1inss[ielem]=tau_s;
                         i1elfr[ielem]=1;
                        i1elprtmp[ielem]=1111;
      if(sigma_n>0){d1inst[ielem]=MINIM(dpeft,d1inst[ielem]);}
      if(tau_s>0){
        d1inss[ielem]=MINIM(dpefs,d1inss[ielem]);}
      if(tau_s<0){
        d1inss[ielem]=MAXIM(-dpefs,d1inss[ielem]);}
                       (*iuptrimesh)=1;
                      
                        }


      //! Applying mixed DFN (i.e., DFN type 3) using the flag assigned to the joint element
      if(iusefn == 3)
      { if(i1edft[ielem] == 1) //! Broken-type DFN crack
      	{ i1elpr[ielem]=iprop-YIPROPMAX;
          d1ebrkf[ielem]=5.0;
          //! For hydrofrac
          i2noid[0][i0]=i1;
          i2noid[1][i1]=i0;
          i2noid[0][i2]=i3;
          i2noid[1][i3]=i2;
        }
        if(i1edft[ielem] == 2) //! Cohesive-type DFN crack
        { op=R2*h*ddfnft/dpepe;
          sp=R2*h*dpefs/dpepe;
          ot=MAXIM(EPSILON,(R3*ddfngf/ddfnft));
          st=MAXIM(EPSILON,(R3*ddfngs/dpefs));
          // dpefs is updated below according to Mohr-Coulomb
        }
      }
                  
      //! Performing excavation: set joint element state to broken if between at least one "excavated" element
      if((i1pexc[i1elpr[el1]]==1)||(i1pexc[i1elpr[el2]]==1))
      { i1elpr[ielem]=iprop-YIPROPMAX; 
        d1ebrkf[ielem]=4.0;
        if(iusehf==1)
	{ //! For hydrofrac
          i2noid[0][i0]=i1;
          i2noid[1][i1]=i0;
          i2noid[0][i2]=i3;
          i2noid[1][i3]=i2;
        }
      }
      else 
      {
      //Fxi0=0.0;
      //Fyi0=0.0;
      for(integ=0;integ<3;integ++)
      { if(integ==0)
        { o=o1+(h*d1inst[ielem]/dpepe); s=s1+(h*d1inss[ielem]/dpepe);
        }
        else if(integ==2)
        { o=o2+(h*d1inst[ielem]/dpepe); s=s2+(h*d1inss[ielem]/dpepe);
        }
        else
        { o=RP5*(o1+o2)+(h*d1inst[ielem]/dpepe); s=RP5*(s1+s2)+(h*d1inss[ielem]/dpepe);
        }
        sabs=ABS(s);
        if((o>op)&&(sabs>sp))
        { dmg=SQRT(((o-op)/ot)*((o-op)/ot)+((sabs-sp)/st)*((sabs-sp)/st));
          nsoft=nsoft+1;
          if((nsoft>2)&&(i1esftf[ielem]==0))  
          { i1esft[*nesft]=ielem;
            d1etsft[*nesft]=dctime;
            d1etmke[ielem]=dctime;
            d2ecsft[0][*nesft]=(d1nccx[i0]+d1nccx[i1]+d1nccx[i2]+d1nccx[i3])/R4;
            d2ecsft[1][*nesft]=(d1nccy[i0]+d1nccy[i1]+d1nccy[i2]+d1nccy[i3])/R4;
            /* Kinetic energy of the joint as soon as it yields */ 
            d1eike[ielem]= 0.5*(d1nmct[i0]*(d1nvcx[i0]*d1nvcx[i0]+d1nvcy[i0]*d1nvcy[i0])+
                                d1nmct[i1]*(d1nvcx[i1]*d1nvcx[i1]+d1nvcy[i1]*d1nvcy[i1])+
                                d1nmct[i2]*(d1nvcx[i2]*d1nvcx[i2]+d1nvcy[i2]*d1nvcy[i2])+
                                d1nmct[i3]*(d1nvcx[i3]*d1nvcx[i3]+d1nvcy[i3]*d1nvcy[i3]));
            d1etike[ielem] = dctime; //! initial time of KE monitoring window
            i1esftf[ielem]=1;
            (*nesft)++;
            (*netsft)++;
          }
        }
        else if(o>op)
        { dmg=(o-op)/ot;
          nsoft=nsoft+1;
          if((nsoft>2)&&(i1esftf[ielem]==0))  
          { 
            i1esft[*nesft]=ielem;
            d1etsft[*nesft]=dctime;
            d1etmke[ielem]=dctime;
            d2ecsft[0][*nesft]=(d1nccx[i0]+d1nccx[i1]+d1nccx[i2]+d1nccx[i3])/R4;
            d2ecsft[1][*nesft]=(d1nccy[i0]+d1nccy[i1]+d1nccy[i2]+d1nccy[i3])/R4;
            /* Kinetic energy of the joint as soon as it yields */ 
            d1eike[ielem]= 0.5*(d1nmct[i0]*(d1nvcx[i0]*d1nvcx[i0]+d1nvcy[i0]*d1nvcy[i0])+
                                d1nmct[i1]*(d1nvcx[i1]*d1nvcx[i1]+d1nvcy[i1]*d1nvcy[i1])+
                                d1nmct[i2]*(d1nvcx[i2]*d1nvcx[i2]+d1nvcy[i2]*d1nvcy[i2])+
                                d1nmct[i3]*(d1nvcx[i3]*d1nvcx[i3]+d1nvcy[i3]*d1nvcy[i3]));
            d1etike[ielem] = dctime; //! initial time of KE monitoring window
            i1esftf[ielem]=2;
            (*nesft)++;
            (*netsft)++;
          }
        }
        else if(sabs>sp)
        { dmg=(sabs-sp)/st;
          nsoft=nsoft+1;
          if((nsoft>2)&&(i1esftf[ielem]==0))  
          { 
            i1esft[*nesft]=ielem;
            d1etsft[*nesft]=dctime;
            d1etmke[ielem]=dctime;
            d2ecsft[0][*nesft]=(d1nccx[i0]+d1nccx[i1]+d1nccx[i2]+d1nccx[i3])/R4;
            d2ecsft[1][*nesft]=(d1nccy[i0]+d1nccy[i1]+d1nccy[i2]+d1nccy[i3])/R4;
            /* Kinetic energy of the joint as soon as it yields */ 
            d1eike[ielem]= 0.5*(d1nmct[i0]*(d1nvcx[i0]*d1nvcx[i0]+d1nvcy[i0]*d1nvcy[i0])+
                                d1nmct[i1]*(d1nvcx[i1]*d1nvcx[i1]+d1nvcy[i1]*d1nvcy[i1])+
                                d1nmct[i2]*(d1nvcx[i2]*d1nvcx[i2]+d1nvcy[i2]*d1nvcy[i2])+
                                d1nmct[i3]*(d1nvcx[i3]*d1nvcx[i3]+d1nvcy[i3]*d1nvcy[i3]));
            d1etike[ielem] = dctime; //! initial time of KE monitoring window
            i1esftf[ielem]=3;
            (*nesft)++;
            (*netsft)++;
          }
        }
        else
        { dmg=R0;
        }
        if(dmg>=R1) /* joint element broken */
        { nfail=nfail+1;
          if((nfail>1)&&(i1elpr[ielem]>=0))
          {
            i1elpr[ielem]=iprop-YIPROPMAX;
            i1ebrk[*nebrk]=ielem;
            d1etbrk[*nebrk]=dctime;
            //i1ebrkf[ielem]=i1esftf[ielem];
            if((o>=(op+ot))&&(sabs>=(sp+st))) // Mode 1 + mode 2 failure
            { d1ebrkf[ielem]=3.0; 
            }
            else if(o>=(op+ot)) // Mode 1 failure
            { d1ebrkf[ielem]=1.0; 
            }
            else if (sabs>=(sp+st)) // Mode 2 failure
            { d1ebrkf[ielem]=2.0; 
            }
            else 
            //{ d1ebrkf[ielem]=1.0+(sabs-sp)/(st-sp); } // Mode 1 + mode 2 failure (vectorial sum of s and o overcomes residual value)
            { d1ebrkf[ielem]=1.0+(sabs-sp)/(st); 
            } // Mode 1 + mode 2 failure (vectorial sum of s and o overcomes residual value)
            d2ecbrk[0][*nebrk]=(d1nccx[i0]+d1nccx[i1]+d1nccx[i2]+d1nccx[i3])/R4;
            d2ecbrk[1][*nebrk]=(d1nccy[i0]+d1nccy[i1]+d1nccy[i2]+d1nccy[i3])/R4;
            d2ecbrk_NEW[0][ielem]=(d1nccx[i0]+d1nccx[i1]+d1nccx[i2]+d1nccx[i3])/R4;
            d2ecbrk_NEW[1][ielem]=(d1nccy[i0]+d1nccy[i1]+d1nccy[i2]+d1nccy[i3])/R4;
            d1elbrk[*nebrk]=SQRT((d1nccx[i0]-d1nccx[i1])*(d1nccx[i0]-d1nccx[i1])+(d1nccy[i0]-d1nccy[i1])*(d1nccy[i0]-d1nccy[i1]));
            (*nebrk)++;
            (*netbrk)++;
            //! For hydrofrac
             i2noid[0][i0]=i1;
             i2noid[1][i1]=i0;
             i2noid[0][i2]=i3;
             i2noid[1][i3]=i2;
          }
          dmg=R1;
        }
        
        /* Specify intact/broken element edges */
        if(i0!=i3 && i1!=i2)
        { if(dmg==R1)
          { i2eledge[el1edge][el1]=-1;  /* -1: broken */
            i2eledge[el2edge][el2]=-1;  
          }
          else
          { i2eledge[el1edge][el1]=1;   /* 1: intact  */ 
            i2eledge[el2edge][el2]=1;
        } }
        
        /* Calculation of stress multiplier (z) from damage coefficient (dmg) */
        
        if(iusehy==0) /* Loading curve = unloading curve (classic formulation of Y-code) */
        { z=(R1-((dpefa+dpefb-R1)/(dpefa+dpefb))*exp(dmg*(dpefa+dpefc*dpefb)/((dpefa+dpefb)*(R1-dpefa-dpefb))))*(dpefa*(R1-dmg)+dpefb*pow((R1-dmg),dpefc));
        }
        else /* Hysteretic model with linear unloading */
        { if((integ==0) && (dmg>=d2eldmg[0][ielem])) /* Loading for integration point 0 */
          { d2eldmg[0][ielem]=dmg;
            z=(R1-((dpefa+dpefb-R1)/(dpefa+dpefb))*exp(d2eldmg[0][ielem]*(dpefa+dpefc*dpefb)/((dpefa+dpefb)*(R1-dpefa-dpefb))))*(dpefa*(R1-d2eldmg[0][ielem])+dpefb*pow((R1-d2eldmg[0][ielem]),dpefc));
          }
          else if((integ==1) && (dmg>=d2eldmg[1][ielem])) /* Loading for integration point 1 */
          { d2eldmg[1][ielem]=dmg;
            z=(R1-((dpefa+dpefb-R1)/(dpefa+dpefb))*exp(d2eldmg[1][ielem]*(dpefa+dpefc*dpefb)/((dpefa+dpefb)*(R1-dpefa-dpefb))))*(dpefa*(R1-d2eldmg[1][ielem])+dpefb*pow((R1-d2eldmg[1][ielem]),dpefc));
          }
          else if((integ==2) && (dmg>=d2eldmg[2][ielem]))  /* Loading for integration point 2 */
          { d2eldmg[2][ielem]=dmg;
            z=(R1-((dpefa+dpefb-R1)/(dpefa+dpefb))*exp(d2eldmg[2][ielem]*(dpefa+dpefc*dpefb)/((dpefa+dpefb)*(R1-dpefa-dpefb))))*(dpefa*(R1-d2eldmg[2][ielem])+dpefb*pow((R1-d2eldmg[2][ielem]),dpefc));
          }
          if((integ==0) && (dmg<d2eldmg[0][ielem])) /* Unloading for integration point 0 */
          {  z=(R1-((dpefa+dpefb-R1)/(dpefa+dpefb))*exp(d2eldmg[0][ielem]*(dpefa+dpefc*dpefb)/((dpefa+dpefb)*(R1-dpefa-dpefb))))*(dpefa*(R1-d2eldmg[0][ielem])+dpefb*pow((R1-d2eldmg[0][ielem]),dpefc));
             z=((op+dmg*ot)/(op+d2eldmg[0][ielem]*ot))*z;
          }
          else if((integ==1) && (dmg<d2eldmg[1][ielem])) /* Unloading for integration point 1 */
          {  z=(R1-((dpefa+dpefb-R1)/(dpefa+dpefb))*exp(d2eldmg[1][ielem]*(dpefa+dpefc*dpefb)/((dpefa+dpefb)*(R1-dpefa-dpefb))))*(dpefa*(R1-d2eldmg[1][ielem])+dpefb*pow((R1-d2eldmg[1][ielem]),dpefc));
             z=((op+dmg*ot)/(op+d2eldmg[1][ielem]*ot))*z;
          }
          else if((integ==2) && (dmg<d2eldmg[2][ielem])) /* Unloading for integration point 2 */
          {  z=(R1-((dpefa+dpefb-R1)/(dpefa+dpefb))*exp(d2eldmg[2][ielem]*(dpefa+dpefc*dpefb)/((dpefa+dpefb)*(R1-dpefa-dpefb))))*(dpefa*(R1-d2eldmg[2][ielem])+dpefb*pow((R1-d2eldmg[2][ielem]),dpefc));
             z=((op+dmg*ot)/(op+d2eldmg[2][ielem]*ot))*z;
          }
        }  
  
        if(o<R0)           /* normal stress */
        { sigma=R2*o*dpeft/op; /* sigma=R0; */
        }
        else if(o>op)
        { sigma=dpeft*z; nsoft=nsoft+1;
        }
        else
        { sigma=(R2*o/op-(o/op)*(o/op))*z*dpeft;
        }
        /* take into account Mohr-Coulomb   */
        if(dpeco>R0)
        { if(sigma>R0)
          { 
            if(dusaf>0) //! Anisotropic fracture model
            { 
              C = dusaf;
              dpefs = dpecor + (dpeco-dpecor) * pow((gamma/90.0),C); 
            }
            else
            { 
              dpefs = dpeco; 
            }
            //! Use "cohesive" DFN properties if joint element belongs to DFN
            //if((iusefn == 2)&&(i1edfnf[ielem]==1))
            if(((iusefn == 2)&&(i1edfnf[ielem]==1))||((iusefn == 3) && (i1edft[ielem] == 2)))
            { dpefs = ddfnco; }
          }
          else
          { 
            if(dusaf>0) //! Anisotropic fracture model
            { 
              C = dusaf; 
              dpefs = dpecor + (dpeco-dpecor) * pow((gamma/90.0),C) - sigma * (dpefrrd + (dpefr - dpefrrd) * pow((gamma/90.0),C));               
            }
            else
            { 
              dpefs = dpeco-sigma*dpefr; 
            }
            //! Use "cohesive" DFN properties if joint element belongs to DFN
            //if((iusefn == 2)&&(i1edfnf[ielem]==1))
            if(((iusefn == 2)&&(i1edfnf[ielem]==1))||((iusefn == 3) && (i1edft[ielem] == 2)))
            { dpefs = ddfnco-sigma*dpefr;}
        } }
        if((sigma>R0)&&(sabs>sp))           /* shear stress */
        { tau=z*dpefs;
        }
        else if(sigma>R0)
        { 
           tau=R2*sabs*dpefs/sp;
          //tau=(R2*(sabs/sp)-(sabs/sp)*(sabs/sp))*z*dpefs;
        }
        else if(sabs>sp)
        { tau=z*dpefs-dpefm*sigma;
        }
        else
        { tau=R2*sabs*dpefs/sp;
          //tau=(R2*(sabs/sp)-(sabs/sp)*(sabs/sp))*(z*dpefs-dpefm*sigma);
        }
        d1elfs[ielem]=dpefs;    /* update fs to Mohr-Coulomb    */
        if(s<R0)tau=-tau;
        if(integ==0)  /* nodal forces */
        { area=h/YR6; /* area=h/6.0; */
          d1nfcx[i0]=d1nfcx[i0]-area*(tau*e1x-sigma*e1y);
          d1nfcx[i3]=d1nfcx[i3]+area*(tau*e1x-sigma*e1y);
          d1nfcy[i0]=d1nfcy[i0]-area*(tau*e1y+sigma*e1x);
          d1nfcy[i3]=d1nfcy[i3]+area*(tau*e1y+sigma*e1x);
          /* Element fracture energy */
          d1efe[ielem]=d1efe[ielem]+dcstec*((-area*(tau*e1x-sigma*e1y))*d1nvcx[i0]+(-area*(tau*e1y+sigma*e1x))*d1nvcy[i0]+
                        (area*(tau*e1x-sigma*e1y))*d1nvcx[i3]+(area*(tau*e1y+sigma*e1x))*d1nvcy[i3]);
          
          /*if(ielem==3)
          { Fxi0=Fxi0-area*(tau*e1x-sigma*e1y);
            Fyi0=Fyi0-area*(tau*e1y+sigma*e1x); }*/
        }
        else if(integ==1)
        { area=h/R3;  /* area=h/3.0; */
          d1nfcx[i0]=d1nfcx[i0]-area*(tau*e1x-sigma*e1y);
          d1nfcx[i3]=d1nfcx[i3]+area*(tau*e1x-sigma*e1y);
          d1nfcy[i0]=d1nfcy[i0]-area*(tau*e1y+sigma*e1x);
          d1nfcy[i3]=d1nfcy[i3]+area*(tau*e1y+sigma*e1x);
          d1nfcx[i1]=d1nfcx[i1]-area*(tau*e1x-sigma*e1y);
          d1nfcx[i2]=d1nfcx[i2]+area*(tau*e1x-sigma*e1y);
          d1nfcy[i1]=d1nfcy[i1]-area*(tau*e1y+sigma*e1x);
          d1nfcy[i2]=d1nfcy[i2]+area*(tau*e1y+sigma*e1x);
          /* Element fracture energy */
          d1efe[ielem]=d1efe[ielem]+dcstec*((-area*(tau*e1x-sigma*e1y))*d1nvcx[i0]+(-area*(tau*e1y+sigma*e1x))*d1nvcy[i0]+
                       (-area*(tau*e1x-sigma*e1y))*d1nvcx[i1]+(-area*(tau*e1y+sigma*e1x))*d1nvcy[i1]+
                       (+area*(tau*e1x-sigma*e1y))*d1nvcx[i2]+(+area*(tau*e1y+sigma*e1x))*d1nvcy[i2]+
                       (+area*(tau*e1x-sigma*e1y))*d1nvcx[i3]+(+area*(tau*e1y+sigma*e1x))*d1nvcy[i3]);

          /*if(ielem==3)
          { Fxi0=Fxi0-area*(tau*e1x-sigma*e1y);  
            Fyi0=Fyi0-area*(tau*e1y+sigma*e1x); }*/
          //if(ielem==3)
          //{ fprintf(out1,"%.6f \t %.6f \t %.6f \n",dctime,sigma,tau); }
        }
        else
        { area=h/YR6; /* area=h/6.0; */
          d1nfcx[i1]=d1nfcx[i1]-area*(tau*e1x-sigma*e1y);
          d1nfcx[i2]=d1nfcx[i2]+area*(tau*e1x-sigma*e1y);
          d1nfcy[i1]=d1nfcy[i1]-area*(tau*e1y+sigma*e1x);
          d1nfcy[i2]=d1nfcy[i2]+area*(tau*e1y+sigma*e1x);
          /* Element fracture energy */
          d1efe[ielem]=d1efe[ielem]+dcstec*((-area*(tau*e1x-sigma*e1y))*d1nvcx[i1]+(-area*(tau*e1y+sigma*e1x))*d1nvcy[i1]+
                       (+area*(tau*e1x-sigma*e1y))*d1nvcx[i2]+(+area*(tau*e1y+sigma*e1x))*d1nvcy[i2]);
        }
      }
      }
      /* If joint is yielded compute differential kinetic energy */
      if(i1esftf[ielem]>0)
      { joint_ke = 0.5*(d1nmct[i0]*(d1nvcx[i0]*d1nvcx[i0]+d1nvcy[i0]*d1nvcy[i0])+
                        d1nmct[i1]*(d1nvcx[i1]*d1nvcx[i1]+d1nvcy[i1]*d1nvcy[i1])+
                        d1nmct[i2]*(d1nvcx[i2]*d1nvcx[i2]+d1nvcy[i2]*d1nvcy[i2])+
                        d1nmct[i3]*(d1nvcx[i3]*d1nvcx[i3]+d1nvcy[i3]*d1nvcy[i3]));
        delta_ke = joint_ke - d1eike[ielem];
        /* If the differential kinetic energy is greater than that from previous timestep update the joint differential kin energy */
        if(delta_ke > d1edke[ielem])
        { d1edke[ielem] = delta_ke; 
          d1etmke[ielem] = dctime;
        }
        if(iusesm == 1) //! Use maximum time window duration
        { if (dctime >= d1etike[ielem] + dctwle)
          { d1eike[ielem] = joint_ke;
            d1etike[ielem] = dctime; 
        } }
      }
      if(iusesm == 2) //! Record energy at the time of failure
      { if(d1ebrkf[ielem]>0)
        { joint_ke = 0.5*(d1nmct[i0]*(d1nvcx[i0]*d1nvcx[i0]+d1nvcy[i0]*d1nvcy[i0])+
                          d1nmct[i1]*(d1nvcx[i1]*d1nvcx[i1]+d1nvcy[i1]*d1nvcy[i1])+
                          d1nmct[i2]*(d1nvcx[i2]*d1nvcx[i2]+d1nvcy[i2]*d1nvcy[i2])+
                          d1nmct[i3]*(d1nvcx[i3]*d1nvcx[i3]+d1nvcy[i3]*d1nvcy[i3]));
          d1edke[ielem] = joint_ke;
          d1etmke[ielem] = dctime;
      } }
    


} 
    }
  }
}







/*********************PUBLIC*************************************/
void Yfd(   ydc,  yde,  ydn, ydmn, ydo,  ydpe, ydpn, ydpj, ydis, ydfn, ydhf, ydsm, ydsb    /***  nodal forces  ***/
        )
  YDC ydc; YDE yde; YDN ydn;YDMN ydmn; YDO ydo; YDPE ydpe; YDPN ydpn; YDPJ ydpj; YDIS ydis; YDFN ydfn; YDHF ydhf; YDSM ydsm; YDSB ydsb;
{ YINT iprop,inopo,jprop,i,j;
  YINT ielem;
  static YINT pmcstep=0; /* previous maximum number of time steps */
  YINT ihys;
  
  YINT s,k,r;
  DBL xi_0,yi_0,xi_1,yi_1;
  
  /* zero model strain energy */
  for(ihys=0;ihys<ydo->nohys;ihys++)
  { if(ydo->i1ohyt[ihys]==(YFLEE)) 
    { ydo->d1ohys[ihys] = 0.0;
    }
  }
    if (ydn->i1remeshf == INT1NULL)
  { ydn->i1remeshf = TalINT1(ydn->mnopo);
    for(i=0;i<ydn->mnopo;i++)
    ydn->i1remeshf[i]=-1;
  }

 /* Initializing d1area */
  if(yde->d1area==DBL1NULL)
  { yde->d1area=TalDBL1(yde->melem);
    for(i=0;i<yde->melem;i++)
       yde->d1area[i]=R0;
  }
 /* Initializing d1area */
  if(yde->d1inst==DBL1NULL)
  { yde->d1inst=TalDBL1(yde->melem);
    yde->d1inss=TalDBL1(yde->melem);

    for(i=0;i<yde->melem;i++)
       yde->d1inst[i]=1000000000000000000.0;
       yde->d1inss[i]=1000000000000000000.0;
  }

  if(yde->d2dmg==DBL2NULL)
  { yde->d2dmg=TalDBL2(3,yde->melem);
    yde->d2D=TalDBL2(3,yde->melem);
    for(i=0;i<3;i++)
    { for(j=0;j<yde->melem;j++)
        yde->d2dmg[i][j]=R0;
        yde->d2D[i][j]=R0;
    }
  }


    /* init. shear strength from joint database */
  if(ydc->ncstep==0 || ydc->ncstep==pmcstep)
  { pmcstep=ydc->mcstep;
    if(ydpj->npjset>0)
    { yde->d1elfs=TalDBL1(yde->melem);
      for(ielem=0;ielem<yde->melem;ielem++)
      { if(yde->i1elpr[ielem]>=ydpe->nprop)  /* joints   */
        { jprop=yde->i1elpr[ielem]-ydpe->nprop;
          if((ydpj->i1ptyp[jprop])==(YTE2JOINTS))
          { if(ydpj->d1pjfs[jprop]>R0)
            { yde->d1elfs[ielem]=ydpj->d1pjfs[jprop];
            }
            else
            { yde->d1elfs[ielem]=ydpj->d1pjco[jprop];  
        } } }
        else                                /* trians   */
        { iprop=yde->i1elpr[ielem];
	  // This part has to be commented out (still to figure out why :-)
	  /*if((ydpe->i1ptyp[iprop])==(YTE2TRISOF))
          { yde->d1elfs[ielem]=ydpj->d1pjfs[ydpe->i1pejp[iprop]];
          } */
  } } } }
  /* zero nodal forces and masses */
  for(inopo=0;inopo<ydn->nnopo;inopo++)
  { ydn->d1nmct[inopo]=R0;
    if(ydn->nnodim>0)ydn->d2nfc[0][inopo]=R0;
    if(ydn->nnodim>1)ydn->d2nfc[1][inopo]=R0;
    if(ydn->nnodim>2)ydn->d2nfc[2][inopo]=R0;
  }
  /* zero number of joint elements broken in the timestep */
  yde->netbrk=0;
  /* zero number of joint elements softened in the timestep */
  yde->netsft=0;
  /* Initializing i1ebrk */
  if(yde->i1ebrk == INT1NULL)
  { yde->i1ebrk=TalINT1(yde->melem);
    for (i=0;i<yde->melem;i++)
	{ yde->i1ebrk[i]=-1;
    }
  }
  /* Initializing i1esft */
  if(yde->i1esft == INT1NULL)
  { yde->i1esft=TalINT1(yde->melem);
    for (i=0;i<yde->melem;i++)
	{ yde->i1esft[i]=-1;
    }
  }
  /* Initializing d2ecbrk */
  if(yde->d2ecbrk==DBL2NULL)
  { yde->d2ecbrk=TalDBL2(ydn->mnodim,yde->melem);
    for(i=0;i<ydn->mnodim;i++)
    { for(j=0;j<yde->melem;j++)
        yde->d2ecbrk[i][j]=R0;
    }
  }
  /* Initializing d2ecbrk_NEW */
  if(yde->d2ecbrk_NEW==DBL2NULL)
  { yde->d2ecbrk_NEW=TalDBL2(ydn->mnodim,yde->melem);
    for(i=0;i<ydn->mnodim;i++)
    { for(j=0;j<yde->melem;j++)
        yde->d2ecbrk_NEW[i][j]=R0;
    }
  }
  /* Initializing d2ecsft */
  if(yde->d2ecsft==DBL2NULL)
  { yde->d2ecsft=TalDBL2(ydn->mnodim,yde->melem);
    for(i=0;i<ydn->mnodim;i++)
    { for(j=0;j<yde->melem;j++)
        yde->d2ecsft[i][j]=R0;
    }
  }
  /* initializing d1etbrk */
  if(yde->d1etbrk==DBL1NULL)
  { yde->d1etbrk=TalDBL1(yde->melem);
    for(i=0;i<yde->melem;i++)
      yde->d1etbrk[i]=R0;
  }
  /* initializing d1etsft */
  if(yde->d1etsft==DBL1NULL)
  { yde->d1etsft=TalDBL1(yde->melem);
    for(i=0;i<yde->melem;i++)
      yde->d1etsft[i]=R0;
  }
  /* initializing d1etmke */
  if(yde->d1etmke==DBL1NULL)
  { yde->d1etmke=TalDBL1(yde->melem);
    for(i=0;i<yde->melem;i++)
      yde->d1etmke[i]=R0;
  }
  /* initializing d1elbrk */
  if(yde->d1elbrk==DBL1NULL)
  { yde->d1elbrk=TalDBL1(yde->melem);
    for(i=0;i<yde->melem;i++)
      yde->d1elbrk[i]=R0;
  }
  /* Initializing d1efe */
  if(yde->d1efe==DBL1NULL)
  { yde->d1efe=TalDBL1(yde->melem);
    for(i=0;i<yde->melem;i++)
	  yde->d1efe[i]=R0;
  }
  /* Initializing d1esftf */
  if(yde->i1esftf==INT1NULL)
  { yde->i1esftf=TalINT1(yde->melem);
    for(i=0;i<yde->melem;i++)
	  yde->i1esftf[i]=0;
  }
  /* Initializing d1ebrkf */
  if(yde->d1ebrkf==DBL1NULL)
  { yde->d1ebrkf=TalDBL1(yde->melem);
    for(i=0;i<yde->melem;i++)
	  yde->d1ebrkf[i]=R0;
  }
  /* Initializing d1eike */
  if(yde->d1eike==DBL1NULL)
  { yde->d1eike=TalDBL1(yde->melem);
    for(i=0;i<yde->melem;i++)
       yde->d1eike[i]=R0;
  }
  /* Initializing d1edke */
  if(yde->d1edke==DBL1NULL)
  { yde->d1edke=TalDBL1(yde->melem);
    for(i=0;i<yde->melem;i++)
       yde->d1edke[i]=R0;
  }
  /* Initializing d1etike */
  if(yde->d1etike==DBL1NULL)
  { yde->d1etike=TalDBL1(yde->melem);
    for(i=0;i<yde->melem;i++)
       yde->d1etike[i]=R0;
  }

  /* Initializing i1psup */
  if(ydpe->i1psup==INT1NULL)
  { ydpe->i1psup=TalINT1(ydpe->mprop);
    for (i=0;i<ydpe->mprop;i++)
    { yde->i1ebrk[i]=0; }
  }
  
//   /* Initializing i1pjhy with 0 if not specified in the input file */
//   if(ydpj->i1pjhy==INT1NULL)
//   { ydpj->i1pjhy=TalINT1(ydpj->mpjset);
//     for (i=0;i<ydpj->mpjset;i++)
//     { ydpj->i1pjhy[i]=0; }
//   }
  

  /* Initializing d2eldmg only if hysteretic joint model is used */
  if(ydpj->iusehy>0)
  { if(yde->d2eldmg==DBL2NULL)
    { yde->d2eldmg=TalDBL2(3,yde->melem);
      for(i=0;i<3;i++)
      { for(j=0;j<yde->melem;j++)
        yde->d2eldmg[i][j]=R0;
      }
    }
  }
  
  /* Initializing d2elstr only if rebars are used */
  if(ydsb->nsbar>0)
  { if(yde->d2elstr==DBL2NULL)
    { yde->d2elstr=TalDBL2(4,yde->melem);
      for(i=0;i<4;i++)
      { for(j=0;j<yde->melem;j++)
        yde->d2elstr[i][j]=R0;
      }
    }
  }
  
  /* Initializing d2nc0 (current coordinates at timestep 0, used to output displacements */
  if(ydn->d2nc0==DBL2NULL)
  { ydn->d2nc0=TalDBL2(4,ydn->mnopo);
    for(i=0;i<2;i++)
    { for(j=0;j<ydn->mnopo;j++)
      ydn->d2nc0[i][j]=ydn->d2ncc[i][j];
    }
  }
    


/* --- 非旋转构型 --- */
if (ydc->ncstep == 0) {

    /* 旋转张量 R */
    if (yde->d2elR == DBL2NULL) {
        yde->d2elR = TalDBL2(9, yde->melem);
        for (i = 0; i < 9; i++)
            for (ielem = 0; ielem < yde->melem; ielem++)
                yde->d2elR[i][ielem] = (i==0 || i==4 || i==8) ? 1.0 : 0.0;
    }

    /* 左伸长张量 V */
    if (yde->d2elV == DBL2NULL) {
        yde->d2elV = TalDBL2(9, yde->melem);
        for (i = 0; i < 9; i++)
            for (ielem = 0; ielem < yde->melem; ielem++)
                yde->d2elV[i][ielem] = (i==0 || i==4 || i==8) ? 1.0 : 0.0;
    }

    /* 非旋转构型应力 alpha */
    if (yde->d2elAlpha == DBL2NULL) {
        yde->d2elAlpha = TalDBL2(9, yde->melem);
        for (i = 0; i < 9; i++)
            for (ielem = 0; ielem < yde->melem; ielem++)
                yde->d2elAlpha[i][ielem] = 0.0;
    }

    /* 内部变量 state（如果需要） */
    if (yde->d2elState == DBL2NULL) {
        yde->d2elState = TalDBL2(0, yde->melem);
        for (i = 0; i < 0; i++)
            for (ielem = 0; ielem < yde->melem; ielem++)
                yde->d2elState[i][ielem] = 0.0;
    }
}




/* --- 初始化显式塑性 + 扩展次载面（DP版）状态量，只在 ncstep==0 时执行 --- */
if (ydc->ncstep == 0) {

  /* 塑性形变梯度 Fp */
  if (yde->Fp00 == DBL2NULL) {
    yde->Fe00 = TalDBL2(NINT, yde->melem);
    yde->Fe01 = TalDBL2(NINT, yde->melem);
    yde->Fe10 = TalDBL2(NINT, yde->melem);
    yde->Fe11 = TalDBL2(NINT, yde->melem);
    yde->Fe33 = TalDBL2(NINT, yde->melem);

    
    
    yde->Fp00 = TalDBL2(NINT, yde->melem);
    yde->Fp01 = TalDBL2(NINT, yde->melem);
    yde->Fp10 = TalDBL2(NINT, yde->melem);
    yde->Fp11 = TalDBL2(NINT, yde->melem);
    yde->Fp33 = TalDBL2(NINT, yde->melem);
    yde->Fp32 = TalDBL2(NINT, yde->melem);
    yde->Fp31 = TalDBL2(NINT, yde->melem);

    yde->sigma00 = TalDBL2(NINT, yde->melem);
    yde->sigma10 = TalDBL2(NINT, yde->melem);
    yde->sigma11 = TalDBL2(NINT, yde->melem);
    yde->sigma33 = TalDBL2(NINT, yde->melem);
    yde->sigma31 = TalDBL2(NINT, yde->melem);
    yde->sigma32 = TalDBL2(NINT, yde->melem);
    yde->ss00 = TalDBL2(NINT, yde->melem);
    yde->ss01 = TalDBL2(NINT, yde->melem);
    yde->ss10 = TalDBL2(NINT, yde->melem);
    yde->ss11 = TalDBL2(NINT, yde->melem);
    yde->ss33 = TalDBL2(NINT, yde->melem);
    yde->ss31 = TalDBL2(NINT, yde->melem);
    yde->ss32 = TalDBL2(NINT, yde->melem);
    for (i = 0; i < NINT; i++)
      for (ielem = 0; ielem < yde->melem; ielem++) {
        yde->Fp00[i][ielem] = 1.0;
        yde->Fp01[i][ielem] = 0.0;
        yde->Fp10[i][ielem] = 0.0;
        yde->Fp11[i][ielem] = 1.0;
        yde->Fp33[i][ielem] = 1.0; 
        yde->Fp31[i][ielem] = 0.0;
        yde->Fp32[i][ielem] = 0.0;
        
        yde->Fe00[i][ielem] = 1.0;
        yde->Fe01[i][ielem] = 0.0;
        yde->Fe10[i][ielem] = 0.0;
        yde->Fe11[i][ielem] = 1.0;
        yde->Fe33[i][ielem] = 1.0;

        yde->sigma00[i][ielem] = 1.0e-9;
        yde->sigma10[i][ielem] = 0.0;
        yde->sigma11[i][ielem] = 0.0;
        yde->sigma33[i][ielem] = 0.0;
        yde->sigma31[i][ielem] = 0.0;
        yde->sigma32[i][ielem] = 0.0;
        yde->ss00[i][ielem] = 1.0e-12;
        yde->ss01[i][ielem] = -yde->ss00[i][ielem]/2.0;
        yde->ss10[i][ielem] = -yde->ss00[i][ielem]/2.0;
        yde->ss11[i][ielem] = 0.0;
        yde->ss33[i][ielem] = 0.0;
        yde->ss31[i][ielem] = 0.0;
        yde->ss32[i][ielem] = 0.0;
      }
  }




  /* 动硬化耗散分量 Fpkd */
  if (yde->Fpkd00 == DBL2NULL) {
    yde->Fpkd00 = TalDBL2(NINT, yde->melem);
    yde->Fpkd01 = TalDBL2(NINT, yde->melem);
    yde->Fpkd10 = TalDBL2(NINT, yde->melem);
    yde->Fpkd11 = TalDBL2(NINT, yde->melem);
    yde->Fpkd33 = TalDBL2(NINT, yde->melem);
    for (i = 0; i < NINT; i++)
      for (ielem = 0; ielem < yde->melem; ielem++) {
        yde->Fpkd00[i][ielem] = 1.0;
        yde->Fpkd01[i][ielem] = 0.0;
        yde->Fpkd10[i][ielem] = 0.0;
        yde->Fpkd11[i][ielem] = 1.0;
        yde->Fpkd33[i][ielem] = 1.0;
      }
  }

  /* 弹性核耗散分量 Fpcd */
  if (yde->Fpcd00 == DBL2NULL) {
    yde->Fpcd00 = TalDBL2(NINT, yde->melem);
    yde->Fpcd01 = TalDBL2(NINT, yde->melem);
    yde->Fpcd10 = TalDBL2(NINT, yde->melem);
    yde->Fpcd11 = TalDBL2(NINT, yde->melem);
    yde->Fpcd33 = TalDBL2(NINT, yde->melem);
    for (i = 0; i < NINT; i++)
      for (ielem = 0; ielem < yde->melem; ielem++) {
        yde->Fpcd00[i][ielem] = 1.0;
        yde->Fpcd01[i][ielem] = 0.0;
        yde->Fpcd10[i][ielem] = 0.0;
        yde->Fpcd11[i][ielem] = 1.0;
        yde->Fpcd33[i][ielem] = 1.0;
      }
  }


  /* Mk / Mc （Mandel-like 背应力与弹性核中心） */
  if (yde->Mk00 == DBL2NULL) {
    yde->Mk00 = TalDBL2(NINT, yde->melem);
    yde->Mk01 = TalDBL2(NINT, yde->melem);
    yde->Mk11 = TalDBL2(NINT, yde->melem);
    yde->Mk33 = TalDBL2(NINT, yde->melem);


    yde->Mc00 = TalDBL2(NINT, yde->melem);
    yde->Mc01 = TalDBL2(NINT, yde->melem);
    yde->Mc11 = TalDBL2(NINT, yde->melem);
    yde->Mc33 = TalDBL2(NINT, yde->melem);
    
    
    
    yde->M00 = TalDBL2(NINT, yde->melem);
    yde->M01 = TalDBL2(NINT, yde->melem);
    yde->M10 = TalDBL2(NINT, yde->melem);
    yde->M11 = TalDBL2(NINT, yde->melem);
    yde->M33 = TalDBL2(NINT, yde->melem);
    yde->M31 = TalDBL2(NINT, yde->melem);
    yde->M32 = TalDBL2(NINT, yde->melem);




    for (i = 0; i < NINT; i++)
      for (ielem = 0; ielem < yde->melem; ielem++) {
        yde->Mk00[i][ielem] = 0.0;
        yde->Mk01[i][ielem] = 0.0;
        yde->Mk11[i][ielem] = 0.0;
        yde->Mk33[i][ielem] = 0.0;


        yde->Mc00[i][ielem] = 0.0;
        yde->Mc01[i][ielem] = 0.0;
        yde->Mc11[i][ielem] = 0.0;
        yde->Mc33[i][ielem] = 0.0;




        yde->M00[i][ielem] = 0.0;
        yde->M01[i][ielem] = 0.0;
        yde->M10[i][ielem] = 0.0;
        yde->M11[i][ielem] = 0.0; 
        yde->M33[i][ielem] = 0.0;
        yde->M31[i][ielem] = 0.0;
        yde->M32[i][ielem] = 0.0;
      }
  }

  /* 硬化标量 H / R / Rc */

  if (yde->H == DBL2NULL) {
    yde->H  = TalDBL2(NINT, yde->melem);
    yde->R  = TalDBL2(NINT, yde->melem);
    yde->Rc = TalDBL2(NINT, yde->melem);
    yde->EPS= TalDBL2(NINT, yde->melem);
    for (i=0; i<NINT; i++)
      for (ielem=0; ielem<yde->melem; ielem++) {
        yde->H[i][ielem]  = 0.0;
        yde->R[i][ielem]  = 0.0; /* 或材料参数 Re 初值 */
        yde->Rc[i][ielem] = 0.0;
        yde->EPS[i][ielem]= 0.0;
      }
  }

  /* 可选诊断 detFp / eqp */

  if (yde->detFe == DBL2NULL) {
    yde->detFe = TalDBL2(NINT, yde->melem);
    yde->detFp = TalDBL2(NINT, yde->melem);
    yde->eqp   = TalDBL2(NINT, yde->melem);
    for (i = 0; i < NINT; i++) {
        for (ielem = 0; ielem < yde->melem; ielem++) {
            yde->detFe[i][ielem] = 1.0;  /* 初始 Fe 行列式 */
            yde->detFp[i][ielem] = 1.0;  /* 初始 Fp 行列式 */
            yde->eqp[i][ielem]   = 0.0;  /* 初始等效塑性应变 */
        }
    }
}
}














  /* At time step zero, if a mixed DFN type is used, find and assign crack type to all joint elements belonging to the DFN */
  if(ydc->ncstep == 0) 
  { if(ydfn->iusefn == 3)
    { for(jprop=0;jprop<ydpj->npjset;jprop++) // Loop over joint elements
      { if((ydpj->i1ptyp[jprop])==(YTE2JOINTS))
        { for(ielem=0;ielem<yde->nelem;ielem++)
          { if(yde->i1elpr[ielem]==(jprop+(ydpe->nprop)))
            // Check if two edge nodes (of the joint element) belong to a DFN crack
       	    { xi_0 = ydn->d2nci[0][yde->i2elto[0][ielem]];  
      	      yi_0 = ydn->d2nci[1][yde->i2elto[0][ielem]];                
	      xi_1 = ydn->d2nci[0][yde->i2elto[1][ielem]];
	      yi_1 = ydn->d2nci[1][yde->i2elto[1][ielem]];
	      for(s=0; s<ydfn->mdfnfr; s++)
	      { for(k=0; k<ydfn->mdfnno; k++)
	        { if(ydfn->i2dfnn[k][s] >= 0)  
	          { if(xi_0 == ydn->d2nci[0][ydfn->i2dfnn[k][s]])
	            { if(yi_0 == ydn->d2nci[1][ydfn->i2dfnn[k][s]])
	              { for(r=0; r<ydfn->mdfnno; r++)
	                { if(ydfn->i2dfnn[r][s] >= 0)
		          { if(xi_1 == ydn->d2nci[0][ydfn->i2dfnn[r][s]])
	                    { if(yi_1 == ydn->d2nci[1][ydfn->i2dfnn[r][s]]) 
	                      { yde->i1edft[ielem] = ydfn->i1dfft[s]; //! Assign crack type to joint element (1 = broken, 2 = cohesive)
  } } } } } } } } } } } } } } }
  
  
   
  /* Apply support */
  for(iprop=0;iprop<ydpe->nprop;iprop++)
  { if(ydpe->i1psup[iprop]==1)
    { for(ielem=0;ielem<yde->nelem;ielem++)
      { if(yde->i1elpr[ielem]==iprop) //! Set initial coordinates equal to current coordinates (i.e., reset elastic deformation)
        { ydn->d2nci[0][(yde->i2elto[0][ielem])]=ydn->d2ncc[0][(yde->i2elto[0][ielem])];
          ydn->d2nci[0][(yde->i2elto[1][ielem])]=ydn->d2ncc[0][(yde->i2elto[1][ielem])];
          ydn->d2nci[0][(yde->i2elto[2][ielem])]=ydn->d2ncc[0][(yde->i2elto[2][ielem])];
          ydn->d2nci[1][(yde->i2elto[0][ielem])]=ydn->d2ncc[1][(yde->i2elto[0][ielem])];
          ydn->d2nci[1][(yde->i2elto[1][ielem])]=ydn->d2ncc[1][(yde->i2elto[1][ielem])];
          ydn->d2nci[1][(yde->i2elto[2][ielem])]=ydn->d2ncc[1][(yde->i2elto[2][ielem])];
      } } 
      ydpe->i1psup[iprop]=0;
   } }
   
  for(iprop=0;iprop<ydpe->nprop;iprop++)
  { if( 
        (ydpe->i1ptyp[iprop])==(YTE2PLANESTRAIN) )
    { 

//  /*
Yfd2TRIELS(   
      yde->nelem,
      iprop,
      ydpn->npnfact,ydpn->mpnset,ydpn->npnset,
      ydpn->d3pnfac,
      ydpe->i1ptyp[iprop],
      ydpe->d1peks[iprop],ydpe->d1pela[iprop],
      ydpe->d1pemu[iprop],ydpe->d1pero[iprop],ydpe->d1psem[iprop],
      ydpe->d1peem[iprop],ydpe->d1penu[iprop],
      ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nci[0],ydn->d2nci[1],ydn->d2nfc[0],
      ydn->d2nfc[1],ydn->d1nmct  ,ydn->d2nvc[0],ydn->d2nvc[1],
      ydpn->d1pnaf , ydpn->d1pnap  ,ydpn->d1pnat  ,
      yde->i1elpr,ydn->i1nopr,yde->i2elto,
      ydo->nohys, ydo->dohyp , ydc->dctime,
      ydo->d1ohys, ydo->d1ohyt, ydo->d1ohyx, ydo->d1ohyy,
      ydo->i1ohyt, ydpn->npnset,yde->d1elfr,
      ydpe->i1usan[iprop], ydpe->d1peex[iprop], ydpe->d1peey[iprop],
      ydpe->d1pemx[iprop], ydpe->d1pemy[iprop], ydpe->d1peg[iprop],
      ydis->iuseis, ydis->dcstxx, ydis->dcstxy, ydis->dcstyy,
      ydis->dcsyxx, ydis->dcsyxy, ydis->dcsyyy, ydis->dcsrfy,
      ydpn->i1pnfx, ydpn->i1pnfy,
      ydpe->i1pexc, ydn->i1nowe, ydhf->iusehf,
      ydsb->nsbar,yde->d2elstr,ydc->ncstep,
yde->sigma00, yde->sigma10, yde->sigma11, yde->sigma33,yde->d1area,ydc->dcstec,
yde->Fe00, yde->Fe01, yde->Fe10, yde->Fe11,yde->Fe33,ydmn->d1noa0,ydmn->d1noac,ydn->i1getmaster,ydc->mcstep, ydn->i1nobf0,yde->d1J
      );

  //   */

  
     /*
      Yfd2TRIELS_EP(
    yde->nelem,
      iprop,
      ydpn->npnfact,ydpn->mpnset,ydpn->npnset,
      ydpn->d3pnfac,
      ydpe->i1ptyp[iprop],
      ydpe->d1peks[iprop],ydpe->d1pela[iprop],
      ydpe->d1pemu[iprop],ydpe->d1pero[iprop],ydpe->d1psem[iprop],
      ydpe->d1peem[iprop],ydpe->d1penu[iprop],
      ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nci[0],ydn->d2nci[1],ydn->d2nfc[0],
      ydn->d2nfc[1],ydn->d1nmct  ,ydn->d2nvc[0],ydn->d2nvc[1],
      ydpn->d1pnaf , ydpn->d1pnap  ,ydpn->d1pnat  ,
      yde->i1elpr,ydn->i1nopr,yde->i2elto,
      ydo->nohys, ydo->dohyp , ydc->dctime,
      ydo->d1ohys, ydo->d1ohyt, ydo->d1ohyx, ydo->d1ohyy,
      ydo->i1ohyt, ydpn->npnset,yde->d1elfr,
      ydpe->i1usan[iprop], ydpe->d1peex[iprop], ydpe->d1peey[iprop],
      ydpe->d1pemx[iprop], ydpe->d1pemy[iprop], ydpe->d1peg[iprop],
      ydis->iuseis, ydis->dcstxx, ydis->dcstxy, ydis->dcstyy,
      ydis->dcsyxx, ydis->dcsyxy, ydis->dcsyyy, ydis->dcsrfy,
      ydpn->i1pnfx, ydpn->i1pnfy,
      ydpe->i1pexc, ydn->i1nowe, ydhf->iusehf,
      ydsb->nsbar,yde->d2elstr,ydc->ncstep,
      yde->d2elR , yde->d2elV , yde->d2elAlpha,yde->d1area,ydc->dcstec,

      ydpe->d1F0[iprop], ydpe->d1h2[iprop], ydpe->d1Hc[iprop],
      ydpe->d1Ck[iprop], 
      ydpe->d1bk[iprop], 
      ydpe->d1URu[iprop], 
yde->sigma00[0], yde->sigma10[0], yde->sigma11[0], yde->sigma33[0], yde->sigma31[0], yde->sigma32[0],
yde->ss00[0],    yde->ss10[0],    yde->ss11[0],    yde->ss33[0],    yde->ss31[0],    yde->ss32[0],   
yde->M00[0],     yde->M10[0],      yde->M11[0],     yde->M33[0] ,    yde->M31[0],     yde->M32[0],

    yde->H[0], yde->R[0], yde->Rc[0]


);

    */


    }
   
  }
  for(jprop=0;jprop<ydpj->npjset;jprop++)
  { if((ydpj->i1ptyp[jprop])==(YTE2JOINTS))
    { 
      
       if((ydc->iuseeczm)==1)
    { 
   /*
Yfd2JOINTS_ECZM( 
      yde->nelem,(jprop+ydpe->nprop),
      ydpj->d1pjft[jprop], ydpj->d1pjgf[jprop], ydpj->d1pjgs[jprop],
      ydpj->d1pjco[jprop], ydpj->d1pjfr[jprop], ydpj->d1pjpe[jprop],
      ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nfc[0],ydn->d2nfc[1],
      ydn->d2nvc[0],ydn->d2nvc[1],
      yde->i1elpr,yde->i2elto,yde->d1elfs,ydc->dctime,ydc->dcstec,&(yde->nebrk),&(yde->netbrk),yde->i1ebrk,yde->d2ecbrk,yde->d2ecbrk_NEW,
      yde->d1etbrk,yde->d1elbrk,yde->d1efe,&(yde->nesft),&(yde->netsft),yde->i1esft,yde->d2ecsft,
      yde->d1etsft,yde->i1esftf,yde->d1ebrkf,yde->d1eike,yde->d1edke,ydn->d1nmct,yde->d1etmke,
      yde->i2elnext,yde->i2eledge,
      ydpe->i1pexc,
      ydpj->d1usaf[jprop],ydpj->d1pjal[jprop],
      ydpj->d1pjcr[jprop],ydpj->d1pjfd[jprop],ydpj->d1pjtr[jprop],ydpj->d1pjgr[jprop],ydpj->d1pjsr[jprop],ydn->i1nowe,ydn->i2noid,
      ydfn->iusefn,yde->i1edfnf,ydfn->ddfnft,ydfn->ddfnco,ydfn->ddfngf,ydfn->ddfngs,ydhf->iusehf,
      yde->i1edft,yde->d1etike,ydsm->iusesm,ydsm->dctwle,
      ydpj->iusehy,yde->d2eldmg,   
      yde->sigma00, yde->sigma10, yde->sigma11,yde->d1area,ydc->ncstep, yde->d1inst,yde->d1inss,     yde->i1elprtmp,yde->i1elfr,&(ydc->iuptrimesh),ydn->i1remeshf
    );

       */

    }
      else
    {
      
      Yfd2JOINTS(  /* joint element  */
      yde->nelem,(jprop+ydpe->nprop),
      ydpj->d1pjft[jprop], ydpj->d1pjgf[jprop], ydpj->d1pjgs[jprop],
      ydpj->d1pjco[jprop], ydpj->d1pjfr[jprop], ydpj->d1pjpe[jprop],
      ydn->d2ncc[0],ydn->d2ncc[1],ydn->d2nfc[0],ydn->d2nfc[1],
      ydn->d2nvc[0],ydn->d2nvc[1],
      yde->i1elpr,yde->i2elto,yde->d1elfs,ydc->dctime,ydc->dcstec,&(yde->nebrk),&(yde->netbrk),yde->i1ebrk,yde->d2ecbrk,yde->d2ecbrk_NEW,
      yde->d1etbrk,yde->d1elbrk,yde->d1efe,&(yde->nesft),&(yde->netsft),yde->i1esft,yde->d2ecsft,
      yde->d1etsft,yde->i1esftf,yde->d1ebrkf,yde->d1eike,yde->d1edke,ydn->d1nmct,yde->d1etmke,
      yde->i2elnext,yde->i2eledge,
      ydpe->i1pexc,
      ydpj->d1usaf[jprop],ydpj->d1pjal[jprop],
      ydpj->d1pjcr[jprop],ydpj->d1pjfd[jprop],ydpj->d1pjtr[jprop],ydpj->d1pjgr[jprop],ydpj->d1pjsr[jprop],ydn->i1nowe,ydn->i2noid,
      ydfn->iusefn,yde->i1edfnf,ydfn->ddfnft,ydfn->ddfnco,ydfn->ddfngf,ydfn->ddfngs,ydhf->iusehf,
      yde->i1edft,yde->d1etike,ydsm->iusesm,ydsm->dctwle,
      ydpj->iusehy,yde->d2eldmg,yde->EPS,yde->d2dmg,yde->d2D);




  }

    }
  }
}

