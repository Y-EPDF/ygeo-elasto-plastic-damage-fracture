
/* File   Ymd.c */
#include "Yproto.h"
/**************GRAINS G2 R=0.5**********/
/* File   Ymd.c */
#include <assert.h>
#include "Yproto.h"

/**************MESH GRAINS***********/
// 结构体用于存储三角形索引和对应的预计算角度
typedef struct {
  YINT ielem;
  double clockwise_angle;
} TriAngle;

// 归并排序函数
void merge(TriAngle *arr, int l, int m, int r) {
  int n1 = m - l + 1;
  int n2 = r - m;

  TriAngle *L = malloc(n1 * sizeof(TriAngle));
  TriAngle *R = malloc(n2 * sizeof(TriAngle));

  for (int i = 0; i < n1; i++)
      L[i] = arr[l + i];
  for (int j = 0; j < n2; j++)
      R[j] = arr[m + 1 + j];

  int i = 0, j = 0, k = l;
  while (i < n1 && j < n2) {
      if (L[i].clockwise_angle <= R[j].clockwise_angle) {
          arr[k++] = L[i++];
      } else {
          arr[k++] = R[j++];
      }
  }

  while (i < n1) arr[k++] = L[i++];
  while (j < n2) arr[k++] = R[j++];

  free(L);
  free(R);
}

void merge_sort(TriAngle *arr, int l, int r) {
  if (l < r) {
      int m = l + (r - l) / 2;
      merge_sort(arr, l, m);
      merge_sort(arr, m + 1, r);
      merge(arr, l, m, r);
  }
}


/*********************PUBLIC********************************************************/
void Ymdupdate(   ydc, yde, ydi, ydn,ydmn, ydpe, ydpn, ydpm, ydfn , ydpj   /***  mesh elements  ***/
        )
  YDC ydc; YDE yde; YDI ydi;  YDN ydn;YDMN ydmn; YDPE ydpe; YDPN ydpn; YDPM ydpm; YDFN ydfn;YDPJ ydpj;
{ YINT nelest, nnopst;
  YINT iprop;
  YINT elprflag;
  YINT icom, i,j, ielem, irow;
  YINT inopo,jprop;
 
  YINT i0, i1, i2,u, i3;

if (ydc->iuptrimesh == 1||ydc->ncstep == 0)
  {

    if(ydpj->npjset>0&&ydc->iuptrimesh == 1)
    { 
      for(ielem=0;ielem<yde->melem;ielem++)
      { if(yde->i1elpr[ielem]>=ydpe->nprop&&yde->d1elfs[ielem]==0)  /* joints   */
        { jprop=yde->i1elpr[ielem]-ydpe->nprop;
          if((ydpj->i1ptyp[jprop])==(YTE2JOINTS))
          { if(ydpj->d1pjfs[jprop]>R0)
            { yde->d1elfs[ielem]=ydpj->d1pjfs[jprop];
            }
            else
            { yde->d1elfs[ielem]=ydpj->d1pjco[jprop];  
        } } }
 } } 


if(ydmn->i2elno==INT2NULL)
  { ydmn->i2elno=TalINT2(ydmn->mnopo,50);
    ydmn->i2elsort=TalINT2(ydmn->mnopo,50);
    ydmn->i2elnext = TalINT2(ydmn->mnopo, 50); // 新增链表指针数组
    ydmn->i1head = TalINT1(ydmn->mnopo);       // 新增链表头指针
    ydmn->i1elcount=TalINT1(ydmn->mnopo);


    }





    for (i=0;i<ydmn->mnopo;i++)

    {
      ydmn->i1head[i] = -1; // 初始链表头为 -1 (空链表)
      ydmn->i1elcount[i]=0;
      for(j=0; j<50; j++)
      { ydmn->i2elno[i][j]=-1;
        ydmn->i2elsort[i][j]=-1;
        ydmn->i2elnext[i][j] = -1; // 初始化为 -1 (无后继)

    } }
  
  for(int ielem=yde->nelemst;ielem<yde->nelem;ielem++)
  { 

  if((yde->i2elto[3][ielem]==yde->i2elto[0][ielem]&&  yde->i2elto[2][ielem]==yde->i2elto[1][ielem])||(yde->i1elpr[ielem]<-100))
   { 
         yde->i1elfr[ielem]=1;//标记边界重合节理单元
  } }


    for(int ielem=0;ielem<yde->nelem;ielem++)
    { if (yde->i2elto[3][ielem]<0 &&yde->i1elpr[ielem]<1e8)
      {
     // yde->d2centroid[0][ielem] = (ydn->d2ncc[0][yde->i2elto[0][ielem]] + ydn->d2ncc[0][yde->i2elto[1][ielem]] + ydn->d2ncc[0][yde->i2elto[2][ielem]]) / 3.0;
     // yde->d2centroid[1][ielem] = (ydn->d2ncc[1][yde->i2elto[0][ielem]] + ydn->d2ncc[1][yde->i2elto[1][ielem]] + ydn->d2ncc[1][yde->i2elto[2][ielem]]) / 3.0;
      yde->d2centroid[0][ielem] = (ydn->d2nci[0][yde->i2elto[0][ielem]] + ydn->d2nci[0][yde->i2elto[1][ielem]] + ydn->d2nci[0][yde->i2elto[2][ielem]]) / 3.0;
      yde->d2centroid[1][ielem] = (ydn->d2nci[1][yde->i2elto[0][ielem]] + ydn->d2nci[1][yde->i2elto[1][ielem]] + ydn->d2nci[1][yde->i2elto[2][ielem]]) / 3.0;       
     
      
      for(int i=0;i<3;i++)
         {
           ydmn->i2elsort[ydn->i1getmasterfem[yde->i2elto[i][ielem]]][ydmn->i1elcount[ydn->i1getmasterfem[yde->i2elto[i][ielem]]]] = ielem;
           ydmn->i1elcount[ydn->i1getmasterfem[yde->i2elto[i][ielem]]]++;
         }
    }
if (yde->i2elto[3][ielem]>=0 )
      {
      //yde->d2centroid[0][ielem] = (ydn->d2ncc[0][yde->i2elto[0][ielem]] + ydn->d2ncc[0][yde->i2elto[1][ielem]] + ydn->d2ncc[0][yde->i2elto[2][ielem]]+ ydn->d2ncc[0][yde->i2elto[3][ielem]]) / 4.0;
     // yde->d2centroid[1][ielem] = (ydn->d2ncc[1][yde->i2elto[0][ielem]] + ydn->d2ncc[1][yde->i2elto[1][ielem]] + ydn->d2ncc[1][yde->i2elto[2][ielem]]+ ydn->d2ncc[1][yde->i2elto[3][ielem]]) / 4.0;
     yde->d2centroid[0][ielem] = (ydn->d2nci[0][yde->i2elto[0][ielem]] + ydn->d2nci[0][yde->i2elto[1][ielem]] + ydn->d2nci[0][yde->i2elto[2][ielem]]+ ydn->d2nci[0][yde->i2elto[3][ielem]]) / 4.0;
      yde->d2centroid[1][ielem] = (ydn->d2nci[1][yde->i2elto[0][ielem]] + ydn->d2nci[1][yde->i2elto[1][ielem]] + ydn->d2nci[1][yde->i2elto[2][ielem]]+ ydn->d2nci[1][yde->i2elto[3][ielem]]) / 4.0;

           ydmn->i2elsort[ydn->i1getmasterfem[yde->i2elto[0][ielem]]][ydmn->i1elcount[ydn->i1getmasterfem[yde->i2elto[0][ielem]]]] = ielem;
           ydmn->i1elcount[ydn->i1getmasterfem[yde->i2elto[0][ielem]]]++;
           ydmn->i2elsort[ydn->i1getmasterfem[yde->i2elto[1][ielem]]][ydmn->i1elcount[ydn->i1getmasterfem[yde->i2elto[1][ielem]]]] = ielem;
           ydmn->i1elcount[ydn->i1getmasterfem[yde->i2elto[1][ielem]]]++;
         
    }

  }
    
  for (int i = 0; i < ydmn->nnopo; i++) 
  {if(ydmn->i1elcount[i] > 0) {
      YINT count = ydmn->i1elcount[i]; // 实际相邻单元数
      TriAngle *tri_angles = malloc(count * sizeof(TriAngle));
      
      // 1. 计算每个单元的相对角度
      for(int t = 0; t < count; t++) {
        YINT ielem = ydmn->i2elsort[i][t];
        //DBL x0 = ydmn->d2ncc[0][i];
        //DBL y0 = ydmn->d2ncc[1][i];
        DBL x0 = ydmn->d2nci[0][i];
        DBL y0 = ydmn->d2nci[1][i];

        DBL dx = yde->d2centroid[0][ielem] - x0;
        DBL dy = yde->d2centroid[1][ielem] - y0;
        double angle = atan2(dy, dx);
        if(angle < 0) angle += 2 * MYPI;
        tri_angles[t].clockwise_angle = fmod(2 * MYPI - (angle - (MYPI*0.5)), 2 * MYPI);
        tri_angles[t].ielem = ielem;
      }
 // 2. 归并排序
      merge_sort(tri_angles, 0, count - 1);
      
      // 3. 更新排序后的单元列表
      for(int t = 0; t < count; t++) {
        ydmn->i2elsort[i][t] = tri_angles[t].ielem;
      }
      
      // 4. 构建循环链表
      for(int t = 0; t < count; t++) {
        if(t == count - 1) {
          // 最后一个指向第一个
          ydmn->i2elnext[i][t] = 0; 
        } else {
          // 其他指向下一个
          ydmn->i2elnext[i][t] = t + 1; 
        }
      }
      
      // 5. 设置链表头
      ydmn->i1head[i] = 0; // 以第一个单元为头节点
      
      free(tri_angles);
    }
 int master=-1;
        for (int t = 0; t < ydmn->i1elcount[i]; t++) {
         if (ydmn->i2elno[i][t]<0) {
            ydn->i1getmaster[ydmn->i2elno[i][t]] = -1;
         }
              if (yde->i2elto[3][ydmn->i2elsort[i][t]]<0) {
            for (int e = 0; e < 3; e++) {

              if (yde->i2elto[e][ydmn->i2elsort[i][t]] != -1 && ydn->i1getmasterfem[yde->i2elto[e][ydmn->i2elsort[i][t]]] == i) {
                  ydmn->i2elno[i][t] = yde->i2elto[e][ydmn->i2elsort[i][t]];
                                if (t==0||t==1||(t==2&&ydmn->i2elno[i][0]<0&&ydmn->i2elno[i][1]<0)) {master=yde->i2elto[e][ydmn->i2elsort[i][t]];}
                  ydn->i1getmaster[ydmn->i2elno[i][t]] = master;
              }
          }
          }
      }

//     printf("Master fem/add Node %d: Slave Order\n", i);
//       for (int j = 0; j < ydmn->i1elcount[i]; j++)
//      {
//          printf("  Slave element%d: (count=%d, elno=%d,master=%d)\n", ydmn->i2elsort[i][j], ydmn->i1elcount[i] ,ydmn->i2elno[i][j], ydn->i1getmaster[ydmn->i2elno[i][j]] );
//        }
        

   }






// -------------------------- 新增边界节点批量更新逻辑 --------------------------
// 遍历所有主节点（ydmn->nnopo为总主节点数，与原代码遍历一致）
for (int i = 0; i < ydmn->nnopo; i++) 
{
    // 仅处理有有效关联单元/从节点的主节点，跳过空主节点
    if (ydmn->i1elcount[i] <= 0) continue;

    YINT has_boundary_slave = 0; // 标记：该主节点下是否存在边界从节点（0=无，1=有）
    YINT slave_node;            // 临时变量：存储当前遍历的从节点索引

    // 第一步：检测当前主节点下是否存在任意一个边界从节点
    for (int t = 0; t < ydmn->i1elcount[i]; t++) 
    {
        slave_node = ydmn->i2elno[i][t];
        // 过滤无效从节点（原代码中未分配的从节点记为-1，需跳过）
        if (slave_node == -1) continue;
        // 找到任意一个边界从节点，立即标记并退出检测循环（无需继续遍历）
        if (ydn->i1nobf0[slave_node] == 1) 
        {
            has_boundary_slave = 1;
            break;
        }
    }

    // 第二步：若存在边界从节点，将该主节点下所有有效从节点设为边界
    if (has_boundary_slave == 1) 
    {
        for (int t = 0; t < ydmn->i1elcount[i]; t++) 
        {
            slave_node = ydmn->i2elno[i][t];
            // 仅处理有效从节点，避免对-1索引赋值导致内存越界
            if (slave_node != -1) 
            {
                ydn->i1nobf0[slave_node] = 1; // 统一置为边界节点标记
            }
        }
    }
}
// -------------------------- 新增逻辑结束 --------------------------
/**/









 for (int i = 0; i < ydmn->nnopo; i++)  // 遍历每个主节点
  {
    int count = ydmn->i1elcount[i];
    if (count <= 0) continue;
    
    int start = ydmn->i1head[i];
    int current = start;
    int processed = 0;
    
    // 第一步：处理连续的破坏节理单元（边界情况）
    do {
      int elem_index = ydmn->i2elsort[i][current];
      
      // 检查当前单元是否为破坏节理单元
      if (yde->i2elto[3][elem_index] >= 0 && yde->i1elfr[elem_index] == 1)
      {
        // 找到破坏节理后的下一个单元
        int next_index = ydmn->i2elnext[i][current];
        int next_elem = ydmn->i2elsort[i][next_index];
        
        // 检查下一个单元是否也是破坏节理单元
        if (yde->i2elto[3][next_elem] >= 0 ) {
          assert(yde->i1elfr[next_elem] == 1);
          current = next_index;
          processed++;
          continue;
        }
        
        // 检查下一个单元是否为三角形单元
        if (yde->i2elto[3][next_elem] < 0) // 三角形单元
        {
          // 找到三角形单元对应的从节点作为新的主节点
          int new_master = -1;
          for (int e = 0; e < 3; e++) {
            if (yde->i2elto[e][next_elem] != -1 && 
                ydn->i1getmasterfem[yde->i2elto[e][next_elem]] == i) {
              new_master = yde->i2elto[e][next_elem];
              assert(new_master>=0);
              break;
            }
          }
          
          if (new_master >= 0) {
            // 创建临时链表索引以处理环形结构
            int process_index = next_index;
            int first_processed = process_index;
            
            // 处理后续三角形直到遇到破坏节理
            while (1) {
              int tri_elem = ydmn->i2elsort[i][process_index];
              
              // 确保处理的是三角形单元
              if (yde->i2elto[3][tri_elem] < 0) {
                // 更新该三角形单元的所有从节点
                for (int e = 0; e < 3; e++) {
                  int node = yde->i2elto[e][tri_elem];
                  if (node != -1 && ydn->i1getmasterfem[node] == i) {
                    ydn->i1getmaster[node] = new_master;
                  }
                }
              }
              
              // 移动到下一个单元
              int next = ydmn->i2elnext[i][process_index];
              
              // 检查终止条件：回到起点或遇到破坏节理
              if (yde->i2elto[3][ydmn->i2elsort[i][next]] >= 0 && 
                  yde->i1elfr[ydmn->i2elsort[i][next]] == 1) {
                break;
              }
              
              // 检查是否完成整个环形
              if (next == first_processed) {
                break;
              }
              
              process_index = next;
            }
          }
        }
      }
      
      // 移动到链表中的下一个单元
      current = ydmn->i2elnext[i][current];
      processed++;
    } while (processed < count);
  
  //     printf("Master fem/add Node %d: Slave Order\n", i);
  //     for (int j = 0; j < ydmn->i1elcount[i]; j++)
  //   {
 //          printf("  Slave element%d: (count=%d, elno=%d,master=%d)\n", ydmn->i2elsort[i][j], ydmn->i1elcount[i] ,ydmn->i2elno[i][j], ydn->i1getmaster[ydmn->i2elno[i][j]] );
  //     }
  
  
  
  }



  /* At time step zero, if a mixed DFN type is used, find and assign crack type to all joint elements belonging to the DFN */
  if(ydc->ncstep == 0) 
  { if(ydfn->iusefn == 3)
    { for(int jprop=0;jprop< ydpj->npjset;jprop++) // Loop over joint elements
      { if((ydpj->i1ptyp[jprop])==(YTE2JOINTS))
        { for(ielem=0;ielem<yde->nelem;ielem++)
          { if(yde->i1elpr[ielem]==(jprop+(ydpe->nprop)))
            // Check if two edge nodes (of the joint element) belong to a DFN crack
       	    {DBL xi_0 = ydn->d2nci[0][yde->i2elto[0][ielem]];  
              DBL yi_0 = ydn->d2nci[1][yde->i2elto[0][ielem]];                
              DBL xi_1 = ydn->d2nci[0][yde->i2elto[1][ielem]];
              DBL yi_1 = ydn->d2nci[1][yde->i2elto[1][ielem]];
	      for(int  s=0; s<ydfn->mdfnfr; s++)
	      { for(int k=0; k<ydfn->mdfnno; k++)
	        { if(ydfn->i2dfnn[k][s] >= 0)  
	          { if(xi_0 == ydn->d2nci[0][ydfn->i2dfnn[k][s]])
	            { if(yi_0 == ydn->d2nci[1][ydfn->i2dfnn[k][s]])
	              { for(int r=0; r<ydfn->mdfnno; r++)
	                { if(ydfn->i2dfnn[r][s] >= 0)
		          { if(xi_1 == ydn->d2nci[0][ydfn->i2dfnn[r][s]])
	                    { if(yi_1 == ydn->d2nci[1][ydfn->i2dfnn[r][s]]) 
	                      { 
                          
                          
                          yde->i1edft[ielem] = ydfn->i1dfft[s]; //! Assign crack type to joint element (1 = broken, 2 = cohesive)
  
  
                        } } } } } } } } } } } } } } 
  

}




} 
ydc->iuptrimesh=0;









if(ydmn->d1noa0==DBL1NULL)
  {  
ydmn->d1noa0 = TalDBL1(ydmn->mnopo); 
ydmn->d1noac = TalDBL1(ydmn->mnopo);


    }
for (i=0;i<ydmn->mnopo;i++)

    {ydmn->d1noa0[i] = 0.0;
    ydmn->d1noac[i] = 0.0;
    }


  /* Initializing d1edke */
  if(yde->d1J==DBL1NULL)
  { yde->d1J=TalDBL1(yde->melem);
yde->d1J2=TalDBL1(yde->melem);
yde->d1volce=TalDBL1(yde->melem);
}
    for(i=0;i<yde->melem;i++) {
      yde->d1J[i]=R0;
    yde->d1J2[i]=R0;
      yde->d1volce[i]=R0;
    }
//体积锁定1



 //  /*
 DBL F0[2][2];     
  DBL FX[2][2];     
  DBL F0inv[2][2];  
  DBL FXinv[2][2];  
  DBL voli, volc;
    for(int ielem=0;ielem<yde->nelem;ielem++)
    { if (yde->i2elto[3][ielem]<0 &&yde->i1elpr[ielem]<1e8)
      {
     
      for(i=1;i<3;i++)
      { F0[0][i-1]=ydn->d2nci[0][(yde->i2elto[i][ielem])]-ydn->d2nci[0][(yde->i2elto[0][ielem])];
        F0[1][i-1]=ydn->d2nci[1][(yde->i2elto[i][ielem])]-ydn->d2nci[1][(yde->i2elto[0][ielem])];
        FX[0][i-1]=ydn->d2ncc[0][(yde->i2elto[i][ielem])]-ydn->d2ncc[0][(yde->i2elto[0][ielem])];
        FX[1][i-1]=ydn->d2ncc[1][(yde->i2elto[i][ielem])]-ydn->d2ncc[1][(yde->i2elto[0][ielem])];

      }
      YMATINV2(F0,F0inv,voli);
      YMATINV2(FX,FXinv,volc);
     
    ydmn->d1noa0[ ydn->i1getmaster[yde->i2elto[0][ielem]]] += voli / 6.0;
    ydmn->d1noa0[ydn->i1getmaster[yde->i2elto[1][ielem]]] += voli / 6.0;
    ydmn->d1noa0[ ydn->i1getmaster[yde->i2elto[2][ielem]]] += voli / 6.0;
     ydmn->d1noac[ ydn->i1getmaster[yde->i2elto[0][ielem]]] += volc / 6.0;
    ydmn->d1noac[ydn->i1getmaster[yde->i2elto[1][ielem]]] += volc / 6.0;
    ydmn->d1noac[ ydn->i1getmaster[yde->i2elto[2][ielem]]] += volc / 6.0;      

    }
  }



  


    for(int ielem=0;ielem<yde->nelem;ielem++)
    { if (yde->i2elto[3][ielem]<0 &&yde->i1elpr[ielem]<1e8)
      {
     
  DBL  jbarno[3];//tijisuoding
int nobf=0;
//Jbar=0.0;





for(i=0;i<3;i++)
      { 
        assert(ydmn->d1noac[ydn->i1getmaster[(yde->i2elto[i][ielem])]]>0&&ydmn->d1noa0[ydn->i1getmaster[(yde->i2elto[i][ielem])]]>1e-9);
        jbarno[i]=ydmn->d1noac[ydn->i1getmaster[(yde->i2elto[i][ielem])]]/ydmn->d1noa0[ydn->i1getmaster[(yde->i2elto[i][ielem])]];
      if(ydn->i1nobf0[yde->i2elto[i][ielem]]==1){ nobf=nobf+1; }
    }


            if(nobf==0)
               // if(nobf==0||nobf==3)        
            {
                                            for(i=0;i<3;i++)
                                                  { 

                                                yde->d1J[ielem]=yde->d1J[ielem]+(jbarno[i])/3.0;

       }
                                                  }
   else if (nobf==3)
    {

 yde->d1J[ielem]=-1000000000;
    }
    else 
    {assert(nobf!=3);
     for(i=0;i<3;i++)
       { if(ydn->i1nobf0[yde->i2elto[i][ielem]]!=1){

 yde->d1J[ielem]=yde->d1J[ielem]+(jbarno[i]/(3-nobf));
  
   }}}   




      for(i=1;i<3;i++)
      { F0[0][i-1]=ydn->d2nci[0][(yde->i2elto[i][ielem])]-ydn->d2nci[0][(yde->i2elto[0][ielem])];
        F0[1][i-1]=ydn->d2nci[1][(yde->i2elto[i][ielem])]-ydn->d2nci[1][(yde->i2elto[0][ielem])];
        FX[0][i-1]=ydn->d2ncc[0][(yde->i2elto[i][ielem])]-ydn->d2ncc[0][(yde->i2elto[0][ielem])];
        FX[1][i-1]=ydn->d2ncc[1][(yde->i2elto[i][ielem])]-ydn->d2ncc[1][(yde->i2elto[0][ielem])];

      }
      YMATINV2(F0,F0inv,voli);
   yde->d1volce[ielem]=yde->d1J[ielem]*voli;




    }
  }






















   /*
for (i=0;i<ydmn->mnopo;i++)

    {
    ydmn->d1noac[i] = 0.0;
    }




    for(int ielem=0;ielem<yde->nelem;ielem++)
    { if (yde->i2elto[3][ielem]<0 &&yde->i1elpr[ielem]<1e8)
      {
     
     ydmn->d1noac[ ydn->i1getmaster[yde->i2elto[0][ielem]]] += yde->d1volce[ielem] / 6.0;
    ydmn->d1noac[ydn->i1getmaster[yde->i2elto[1][ielem]]] += yde->d1volce[ielem] / 6.0;
    ydmn->d1noac[ ydn->i1getmaster[yde->i2elto[2][ielem]]] += yde->d1volce[ielem] / 6.0;      

    }
  }


    for(int ielem=0;ielem<yde->nelem;ielem++)
    { if (yde->i2elto[3][ielem]<0 &&yde->i1elpr[ielem]<1e8)
      {
     
  DBL  jbarno[3];//tijisuoding
int nobf=0;
//Jbar=0.0;





for(i=0;i<3;i++)
      { 
        assert(ydmn->d1noac[ydn->i1getmaster[(yde->i2elto[i][ielem])]]>0&&ydmn->d1noa0[ydn->i1getmaster[(yde->i2elto[i][ielem])]]>1e-9);
        jbarno[i]=ydmn->d1noac[ydn->i1getmaster[(yde->i2elto[i][ielem])]]/ydmn->d1noa0[ydn->i1getmaster[(yde->i2elto[i][ielem])]];
      if(ydn->i1nobf0[yde->i2elto[i][ielem]]==1){ nobf=nobf+1; }
    }


// if(nobf==0){
      for(i=0;i<3;i++)
            { 

          yde->d1J2[ielem]=yde->d1J[ielem]+(jbarno[i])/3.0;

//       }
}
//    else if (nobf==3)
//    {

// Jbar=JJ;
//    }
//    else 
//    {assert(nobf==2);
//     for(i=0;i<3;i++)
//       { if(i1nobf0[i2elto[i][ielem]]!=1){

//   Jbar=Jbar+(jbarno[i]/(3-nobf));
//    }}}   

    }
  }


 //   */




















  //  */






//体积锁定2

/*
 DBL F0[2][2];     
  DBL FX[2][2];     
  DBL F0inv[2][2];  
  DBL FXinv[2][2];  
  DBL voli, volc, voli_self, volc_self;  


for(int ielem=0; ielem<yde->nelem; ielem++)
{ 
    if (yde->i2elto[3][ielem]<0 && yde->i1elpr[ielem]<1e8)
    {




      for(int m=1;m<3;m++)
      { F0[0][m-1]=ydn->d2nci[0][(yde->i2elto[m][ielem ])]-ydn->d2nci[0][(yde->i2elto[0][ielem])];
        F0[1][m-1]=ydn->d2nci[1][(yde->i2elto[m][ielem ])]-ydn->d2nci[1][(yde->i2elto[0][ielem ])];
        FX[0][m-1]=ydn->d2ncc[0][(yde->i2elto[m][ielem ])]-ydn->d2ncc[0][(yde->i2elto[0][ielem])];
        FX[1][m-1]=ydn->d2ncc[1][(yde->i2elto[m][ielem ])]-ydn->d2ncc[1][(yde->i2elto[0][ielem ])];

      }
      YMATINV2(F0,F0inv,voli_self);
      YMATINV2(FX,FXinv,volc_self);







        int countedge = 0; // 邻居单元数量
        double sum_neighbor_voli = 0.0; 
         double sum_neighbor_volc = 0.0; 
        int isbound=0;
        // 第一步：收集所有邻居单元的d1J
        for(int i=0; i<3; i++)
        { 
            if(yde->i2elnext[0][yde->i2elel0[i][ielem]] == ielem)
            {
                if(yde->i2elnext[1][yde->i2elel0[i][ielem]] != -1)
                {
                    int neighbor_elem = yde->i2elnext[1][yde->i2elel0[i][ielem]];


      for(int m=1;m<3;m++)
      { F0[0][m-1]=ydn->d2nci[0][(yde->i2elto[m][neighbor_elem ])]-ydn->d2nci[0][(yde->i2elto[0][neighbor_elem ])];
        F0[1][m-1]=ydn->d2nci[1][(yde->i2elto[m][neighbor_elem ])]-ydn->d2nci[1][(yde->i2elto[0][neighbor_elem ])];
        FX[0][m-1]=ydn->d2ncc[0][(yde->i2elto[m][neighbor_elem ])]-ydn->d2ncc[0][(yde->i2elto[0][neighbor_elem ])];
        FX[1][m-1]=ydn->d2ncc[1][(yde->i2elto[m][neighbor_elem ])]-ydn->d2ncc[1][(yde->i2elto[0][neighbor_elem ])];

      }
      YMATINV2(F0,F0inv,voli);
      YMATINV2(FX,FXinv,volc);





sum_neighbor_volc += volc; 
                    sum_neighbor_voli += voli; // 累加邻居
                    countedge++;
                }
                else{
 isbound=1;

                }
            }
            if(yde->i2elnext[1][yde->i2elel0[i][ielem]] == ielem)
            {
                if(yde->i2elnext[0][yde->i2elel0[i][ielem]] != -1)
                {
                                       int neighbor_elem = yde->i2elnext[0][yde->i2elel0[i][ielem]];


      for(int m=1;m<3;m++)
      { F0[0][m-1]=ydn->d2nci[0][(yde->i2elto[m][neighbor_elem ])]-ydn->d2nci[0][(yde->i2elto[0][neighbor_elem ])];
        F0[1][m-1]=ydn->d2nci[1][(yde->i2elto[m][neighbor_elem ])]-ydn->d2nci[1][(yde->i2elto[0][neighbor_elem ])];
        FX[0][m-1]=ydn->d2ncc[0][(yde->i2elto[m][neighbor_elem ])]-ydn->d2ncc[0][(yde->i2elto[0][neighbor_elem ])];
        FX[1][m-1]=ydn->d2ncc[1][(yde->i2elto[m][neighbor_elem ])]-ydn->d2ncc[1][(yde->i2elto[0][neighbor_elem ])];

      }
      YMATINV2(F0,F0inv,voli);
      YMATINV2(FX,FXinv,volc);





sum_neighbor_volc += volc; 
                    sum_neighbor_voli += voli; // 累加邻居
                    countedge++;
                }
                                else{
 isbound=1;

                }
            }
        }

        // 第二步：混合平均（自身+邻居），区分边界/内部单元

        double d1J2_val = 0.0;
sum_neighbor_volc += volc_self; 
                    sum_neighbor_voli += voli_self; 
        // 判断是否为边界单元（假设i1elpr[ielem]标记边界，可根据你的实际标记调整）
        if(isbound==1) // 比如加载边单元标记为1，内部为0
        {            d1J2_val = sum_neighbor_volc / sum_neighbor_voli;
            // 边界单元：
            d1J2_val =(1.0-(countedge/3.0))  * (volc_self/voli_self )+ (countedge/3.0) * d1J2_val;
        }
        else
        {
            assert(countedge==3);
            d1J2_val = sum_neighbor_volc / sum_neighbor_voli;
        }

        // 赋值给最终的d1J2
        yde->d1J2[ielem] = d1J2_val;
    }
}

   */














}
