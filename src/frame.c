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
 * Should you modify this source code, the Copyright (C) on the mdified code
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
#include "frame.h"
/* Arrays input output and allocation routines */ 
 
/***************PUBLIC***************************************/
void CHRcpynoext(c1,c2)  /* copy no extension */
  CHR *c1; CHR *c2;
{ YINT i;
  i=0;
  while((c2[i]!='\0')&&(c2[i]!='.')&&(i<300))
  { c1[i]=c2[i];
    i=i+1;
  }
  c1[i]='\0';
}

DBL *TalDBL1(m1)
  YINT m1; 
{ YINT isize=m1*sizeof(DBL);
  if(isize==0)return DBL1NULL;
  return (DBL*)MALLOC(isize);
}

DBL **TalDBL2(m2,m1)
  YINT m2; YINT m1; 
{ YINT isize,i2;
  DBL     *p1;
  DBL    **p2;
  void    *v1;

  isize=sizeof(DBL*)*(m2+3)+
        sizeof(DBL )*(m2*m1+3);
  if(isize==0)return DBL2NULL;
  v1=MALLOC(isize);
  p2=(DBL**)v1;
  p1=(DBL*)v1;
  p1=p1+((m2+1)*sizeof(DBL**))/sizeof(DBL)+2;    
  for(i2=0;i2<m2;i2++)
  { p2[i2]=p1+i2*m1; 
  }
  return p2; 
}
DBL ***TalDBL3(m3,m2,m1)
  YINT m3; YINT m2; YINT m1; 
{ YINT isize,i2,i3;
  DBL     *p1b;
  DBL    **p2b,  **p2e, **p2;
  DBL   ***p3b, ***p3e;
  void    *v1;

  isize=sizeof(DBL**)*(m3+3)+
        sizeof(DBL* )*(m3*m2+3)+
        sizeof(DBL )*(m3*m2*m1+3);
  if(isize==0)return DBL3NULL;
  v1=MALLOC(isize);
  p3b=(DBL***)v1; p3e=p3b+m3+1;
  p2b=(DBL**)p3e; p2e=p2b+m3*m2+1;
  p1b=(DBL*)p2e;
  p2=p2b;    
  for(i3=0;i3<m3;i3++)
  { p2=p2b+i3*m2;
    p3b[i3]=p2;
    for(i2=0;i2<m2;i2++)
    { p2[i2]=p1b+i3*m2*m1+i2*m1;
  } }
  return p3b; 
}
YINT *TalINT1(m1)
  YINT m1; 
{ YINT isize=m1*sizeof(YINT);
  if(isize==0)return INT1NULL;
  return (YINT*)MALLOC(isize);   
} 

YINT **TalINT2(m2,m1)
  YINT m2; YINT m1; 
{ YINT isize,i2;
  YINT     *p1;
  YINT    **p2;
  void    *v1;

  isize=sizeof(YINT*)*(m2+3)+
        sizeof(YINT )*(m2*m1+3);
  if(isize==0)return INT2NULL;
  v1=MALLOC(isize);
  p2=(YINT**)v1;
  p1=(YINT*)v1;
  p1=p1+((m2+1)*sizeof(YINT**))/sizeof(YINT)+2;
  for(i2=0;i2<m2;i2++)
  { p2[i2]=p1+i2*m1; 
  }
  return p2; 
}
YINT ***TalINT3(m3,m2,m1)
  YINT m3; YINT m2; YINT m1; 
{ YINT isize,i2,i3;
  YINT     *p1b;
  YINT    **p2b,  **p2e, **p2;
  YINT   ***p3b, ***p3e;
  void    *v1;

  isize=sizeof(YINT**)*(m3+3)+
        sizeof(YINT* )*(m3*m2+3)+
        sizeof(YINT )*(m3*m2*m1+3);
  if(isize==0)return INT3NULL;
  v1=MALLOC(isize);
  p3b=(YINT***)v1; p3e=p3b+m3+1;
  p2b=(YINT**)p3e; p2e=p2b+m3*m2+1;
  p1b=(YINT*)p2e;
  p2=p2b;    
  for(i3=0;i3<m3;i3++)
  { p2=p2b+i3*m2;
    p3b[i3]=p2;
    for(i2=0;i2<m2;i2++)
    { p2[i2]=p1b+i3*m2*m1+i2*m1;
  } }
  return p3b; 
}
/*--------------------------------------------*/
YINT  Getname(argc, argv, name)
  YINT argc; char **argv; CHR *name;
{ YINT i;
  for(i=0;i<argc;i++)
  { if(CHRcmp(argv[i],name,2)==0)return (i+1);
  };
  return i;
}

void TreadDBL1(fptr,n1,d1aray) 
  FILE  *fptr;  YINT n1; DBL *d1aray;
{ YINT i;
  DBL dnum;
  for(i=0;i<n1;i++)
  { DBLr(fptr,&dnum);
    d1aray[i]=dnum;
} } 
 
void TreadDBL2(fptr,n1,n2,d2aray) 
  FILE  *fptr;  YINT n1; YINT n2; DBL **d2aray;
{ YINT i1,i2;
  DBL dnum;
  for(i1=0;i1<n1;i1++)
  { for(i2=0;i2<n2;i2++)
    { DBLr(fptr,&dnum);
      d2aray[i2][i1]=dnum;
} } }
 
void TreadINT1(fptr,n1,i1aray) 
  FILE  *fptr;  YINT n1; YINT *i1aray;
{ YINT i; YINT inum;
  for(i=0;i<n1;i++)
  { INTr(fptr,&inum);
    i1aray[i]=inum;
} } 
 
void TreadINT2(fptr,n1,n2,i2aray) 
  FILE  *fptr;  YINT n1; YINT n2; YINT **i2aray;
{ YINT i1,i2;
  YINT inum;
  for(i1=0;i1<n1;i1++)
  { for(i2=0;i2<n2;i2++)
    { INTr(fptr,&inum);
      i2aray[i2][i1]=inum;
} } }

void TformDBL1(fptr,dinit,m1,d1aray) 
  FILE  *fptr;  YINT m1;  DBL dinit; DBL **d1aray; 
{ YINT i1,n1;
  DBL dnum;
  DBL *d1=*d1aray;
  n1=0;
  if(d1==DBL1NULL)d1=TalDBL1(m1);
  *d1aray=d1;
  for(i1=0;i1<m1;i1++)
  { d1[i1]=dinit;
  }
  if(fptr!=FILENULL)
  { INTr(fptr,&n1);
    for(i1=0;i1<n1;i1++)
    { DBLr(fptr,&dnum);
      d1[i1]=dnum;
  } }
}
void TformDBL2(fptr,dinit,m1,m2,d2aray) 
  FILE  *fptr;  YINT m1; YINT m2; DBL dinit; DBL ***d2aray; 
{ YINT i1,i2,n1,n2,how;
  DBL dnum;
  DBL **d2=*d2aray;
  n1=0; n2=0; how=0;
  if(d2==DBL2NULL)d2=TalDBL2(m1,m2);
  *d2aray=d2;
  for(i1=0;i1<m1;i1++)
  { for(i2=0;i2<m2;i2++)
    { d2[i1][i2]=dinit;
  } }
  if(fptr!=FILENULL)
  { INTr(fptr,&how);
    INTr(fptr,&n1);
    INTr(fptr,&n2);
    if(how==12)
    { for(i1=0;i1<n1;i1++)
      { for(i2=0;i2<n2;i2++)
        { DBLr(fptr,&dnum);
          d2[i1][i2]=dnum;
    } } }
    else
    { for(i2=0;i2<n2;i2++)
      { for(i1=0;i1<n1;i1++)
        { DBLr(fptr,&dnum);
          d2[i1][i2]=dnum;
} } } } }
void TformDBL3(fptr,dinit,m1,m2,m3,d3aray) 
  FILE  *fptr;  YINT m1; YINT m2; YINT m3; DBL dinit; DBL ****d3aray; 
{ YINT i1,i2,i3,n1,n2,n3,how;
  DBL dnum;
  DBL ***d3=*d3aray;
  n1=0; n2=0; n3=0; how=0;
  if(d3==DBL3NULL)d3=TalDBL3(m1,m2,m3);
  *d3aray=d3;
  for(i1=0;i1<m1;i1++)
  { for(i2=0;i2<m2;i2++)
    { for(i3=0;i3<m3;i3++)
      {  d3[i1][i2][i3]=dinit;     
  } } }
  if(fptr!=FILENULL)
  { INTr(fptr,&how);
    INTr(fptr,&n1);
    INTr(fptr,&n2);
    INTr(fptr,&n3);
    if(how==123)
    { for(i1=0;i1<n1;i1++)
      { for(i2=0;i2<n2;i2++)
        { for(i3=0;i3<n3;i3++)
          { DBLr(fptr,&dnum);
            d3[i1][i2][i3]=dnum;
    } } } }
    else if(how==132)
    { for(i1=0;i1<n1;i1++)
      { for(i3=0;i3<n3;i3++)
        { for(i2=0;i2<n2;i2++)
          { DBLr(fptr,&dnum);
            d3[i1][i2][i3]=dnum;
    } } } }
    else if(how==213)
    { for(i2=0;i2<n2;i2++)
      { for(i1=0;i1<n1;i1++)
        { for(i3=0;i3<n3;i3++)
          { DBLr(fptr,&dnum);
            d3[i1][i2][i3]=dnum;
    } } } }
    else if(how==231)
    { for(i2=0;i2<n2;i2++)
      { for(i3=0;i3<n3;i3++)
        { for(i1=0;i1<n1;i1++)
          { DBLr(fptr,&dnum);
            d3[i1][i2][i3]=dnum;
    } } } }
    else if(how==312) 
    { for(i3=0;i3<n3;i3++)
      { for(i1=0;i1<n1;i1++)
        { for(i2=0;i2<n2;i2++)
          { DBLr(fptr,&dnum);
            d3[i1][i2][i3]=dnum;
    } } } }
    else if(how==321) 
    { for(i3=0;i3<n3;i3++)
      { for(i2=0;i2<n2;i2++)
        { for(i1=0;i1<n1;i1++)
          { DBLr(fptr,&dnum);
            d3[i1][i2][i3]=dnum;
} } } } } }

void TformINT1(fptr,iinit,m1,i1aray) 
  FILE  *fptr;  YINT m1;  YINT iinit; YINT **i1aray; 
{ YINT i1,n1;
  YINT inum;
  YINT *i1a=*i1aray;
  n1=0;
  if(i1a==INT1NULL)i1a=TalINT1(m1);
  *i1aray=i1a;
  for(i1=0;i1<m1;i1++)
  { i1a[i1]=iinit;
  }
  if(fptr!=FILENULL)
  { INTr(fptr,&n1);
    for(i1=0;i1<n1;i1++)
    { INTr(fptr,&inum);
      i1a[i1]=inum;
  } }
}
void TformINT2(fptr,iinit,m1,m2,i2aray) 
  FILE  *fptr;  YINT m1; YINT m2; YINT iinit; YINT ***i2aray; 
{ YINT i1,i2,n1,n2,how;
  YINT inum;
  YINT **i2a=*i2aray;
  n1=0; n2=0; how=0;
  if(i2a==INT2NULL)i2a=TalINT2(m1,m2);
  *i2aray=i2a;
  for(i1=0;i1<m1;i1++)
  { for(i2=0;i2<m2;i2++)
    { i2a[i1][i2]=iinit;
  } }
  if(fptr!=FILENULL)
  { INTr(fptr,&how);
    INTr(fptr,&n1);
    INTr(fptr,&n2);
    if(how==12)
    { for(i1=0;i1<n1;i1++)
      { for(i2=0;i2<n2;i2++)
        { INTr(fptr,&inum);
          i2a[i1][i2]=inum;
    } } }
    else
    { for(i2=0;i2<n2;i2++)
      { for(i1=0;i1<n1;i1++)
        { INTr(fptr,&inum);
          i2a[i1][i2]=inum;
} } } } }
void TformINT3(fptr,iinit,m1,m2,m3,i3aray) 
  FILE  *fptr;  YINT m1; YINT m2; YINT m3; YINT iinit; YINT ****i3aray; 
{ YINT i1,i2,i3,n1,n2,n3,how;
  YINT inum;
  YINT ***i3a=*i3aray;
  n1=0; n2=0; n3=0; how=0;
  if(i3a==INT3NULL)i3a=TalINT3(m1,m2,m3);
  *i3aray=i3a;
  for(i1=0;i1<m1;i1++)
  { for(i2=0;i2<m2;i2++)
    { for(i3=0;i3<m3;i3++)
      { i3a[i1][i2][i3]=iinit;
  } } }
  if(fptr!=FILENULL)
  { INTr(fptr,&how);
    INTr(fptr,&n1);
    INTr(fptr,&n2);
    INTr(fptr,&n3);
    if(how==123)
    { for(i1=0;i1<n1;i1++)
      { for(i2=0;i2<n2;i2++)
        { for(i3=0;i3<n3;i3++)
          { INTr(fptr,&inum);
            i3a[i1][i2][i3]=inum;
    } } } }
    else if(how==132)
    { for(i1=0;i1<n1;i1++)
      { for(i3=0;i3<n3;i3++)
        { for(i2=0;i2<n2;i2++)
          { INTr(fptr,&inum);
            i3a[i1][i2][i3]=inum;
    } } } }
    else if(how==213)
    { for(i2=0;i2<n2;i2++)
      { for(i1=0;i1<n1;i1++)
        { for(i3=0;i3<n3;i3++)
          { INTr(fptr,&inum);
            i3a[i1][i2][i3]=inum;
    } } } }
    else if(how==231)
    { for(i2=0;i2<n2;i2++)
      { for(i3=0;i3<n3;i3++)
        { for(i1=0;i1<n1;i1++)
          { INTr(fptr,&inum);
            i3a[i1][i2][i3]=inum;
    } } } }
    else if(how==312)
    { for(i3=0;i3<n3;i3++)
      { for(i1=0;i1<n1;i1++)
        { for(i1=0;i1<n1;i1++)
          { INTr(fptr,&inum);
            i3a[i1][i2][i3]=inum;
    } } } }
    { for(i3=0;i3<n3;i3++)
      { for(i2=0;i2<n2;i2++)
        { for(i3=0;i3<n3;i3++)
          { INTr(fptr,&inum);
            i3a[i1][i2][i3]=inum;
} } } } } }
/*--------------------------------------------*/
FILE *TwriteCHR2(namep,fpt,ndigit,nperli,n1,c2aray,name)
  CHR   *namep; FILE *fpt; YINT ndigit; YINT nperli; YINT  n1;
  CHR **c2aray; CHR *name;
{ YINT i1;
  YINT counter;
  FILE *fptr=fpt;   

  if(fptr==FILENULL)fptr=fopen(namep,"w");
  if((c2aray!=CHR2NULL)&&(fptr!=FILENULL))
  { CHRw(fptr,name) ; CHRwcr(fptr);
    counter=0;   
    for(i1=0;i1<n1;i1++)
    { CHRwsp(fptr); CHRwsp(fptr);  CHRw(fptr,c2aray[i1]);  counter++;
      if(counter>=nperli)
      {  counter=0; CHRwcr(fptr);
    } }
    if(counter>0)CHRwcr(fptr);
  }
  return fptr;
}

FILE *TwriteDBL1(namep,fpt,ndigit,nperli,n1,d1aray,name)
  CHR *namep; FILE *fpt; YINT ndigit; YINT nperli; YINT  n1;
  DBL  *d1aray; CHR *name;
{ YINT i1;
  YINT counter;
  DBL dnum;
  FILE *fptr=fpt;   

  if(fptr==FILENULL)fptr=fopen(namep,"w");
  if((d1aray!=DBL1NULL)&&(fptr!=FILENULL))
  { CHRw(fptr,name) ; CHRwsp(fptr); INTw(fptr,n1,8); CHRwcr(fptr);
    counter=0;   
    for(i1=0;i1<n1;i1++)
    { dnum=d1aray[i1];
      CHRwsp(fptr); DBLw(fptr,dnum,ndigit);  counter++;
      if(counter>=nperli)
      {  counter=0; CHRwcr(fptr);
    } }
    if(counter>0)CHRwcr(fptr);
  }
  return fptr;
}

FILE *TwriteDBL2(namep,fpt,ndigit,nperli,ihow,n1,n2,d2aray,name)
  CHR *namep; FILE *fpt; YINT ndigit; YINT nperli; YINT ihow; YINT  n1; YINT n2; 
  DBL **d2aray; CHR *name;
{ YINT i1,i2;
  YINT counter;
  DBL dnum;   
  FILE *fptr=fpt;  

  if(fptr==FILENULL)fptr=fopen(namep,"w");
  if((d2aray!=DBL2NULL)&&(fptr!=FILENULL))
  { CHRw(fptr,name); INTw(fptr,ihow,8); 
    CHRwsp(fptr); INTw(fptr,n1,8);
    CHRwsp(fptr); INTw(fptr,n2,8);
    CHRwcr(fptr);
    counter=0; 
    if(ihow==12)  
    { for(i1=0;i1<n1;i1++)
      { for(i2=0;i2<n2;i2++)
        { dnum=d2aray[i1][i2];
          CHRwsp(fptr); DBLw(fptr,dnum,ndigit);  counter++;
          if(counter>=nperli)
          { counter=0; CHRwcr(fptr);
    } } } }
    else  
    { for(i2=0;i2<n2;i2++)
      { for(i1=0;i1<n1;i1++)
        { dnum=d2aray[i1][i2];
          CHRwsp(fptr); DBLw(fptr,dnum,ndigit);  counter++;
          if(counter>=nperli)
          { counter=0; CHRwcr(fptr);
    } } } }
    if(counter>0)CHRwcr(fptr);
  }
  return fptr;
}
FILE *TwriteDBL3(namep,fpt,ndigit,nperli,ihow,n1,n2,n3,d3aray,name)
  CHR *namep; FILE *fpt; YINT ndigit; YINT nperli; YINT ihow;
  YINT  n1; YINT n2; YINT n3;
  DBL ***d3aray; CHR *name;
{ YINT i1,i2,i3;
  YINT counter;
  DBL dnum;   
  FILE *fptr=fpt;  

  if(fptr==FILENULL)fptr=fopen(namep,"w");
  if((d3aray!=DBL3NULL)&&(fptr!=FILENULL))
  { CHRw(fptr,name); INTw(fptr,ihow,8); 
    CHRwsp(fptr); INTw(fptr,n1,8);
    CHRwsp(fptr); INTw(fptr,n2,8);
    CHRwsp(fptr); INTw(fptr,n3,8);
    CHRwcr(fptr);
    counter=0; 
    if(ihow==123)  
    { for(i1=0;i1<n1;i1++)
      { for(i2=0;i2<n2;i2++)
        { for(i3=0;i3<n3;i3++)
          { dnum=d3aray[i1][i2][i3];
            CHRwsp(fptr); DBLw(fptr,dnum,ndigit);  counter++;
            if(counter>=nperli)
            { counter=0; CHRwcr(fptr);
    } } } } }
    if(ihow==132)  
    { for(i1=0;i1<n1;i1++)
      { for(i3=0;i3<n3;i3++)
        { for(i2=0;i2<n2;i2++)
          { dnum=d3aray[i1][i2][i3];
            CHRwsp(fptr); DBLw(fptr,dnum,ndigit);  counter++;
            if(counter>=nperli)
            { counter=0; CHRwcr(fptr);
    } } } } }    
    if(ihow==213)  
    { for(i2=0;i2<n2;i2++)
      { for(i1=0;i1<n1;i1++)
        { for(i3=0;i3<n3;i3++)
          { dnum=d3aray[i1][i2][i3];
            CHRwsp(fptr); DBLw(fptr,dnum,ndigit);  counter++;
            if(counter>=nperli)
            { counter=0; CHRwcr(fptr);
    } } } } }    
    if(ihow==231)  
    { for(i2=0;i2<n2;i2++)
      { for(i3=0;i3<n3;i3++)
        { for(i1=0;i1<n1;i1++)
          { dnum=d3aray[i1][i2][i3];
            CHRwsp(fptr); DBLw(fptr,dnum,ndigit);  counter++;
            if(counter>=nperli)
            { counter=0; CHRwcr(fptr);
    } } } } }    
    if(ihow==312)  
    { for(i3=0;i3<n3;i3++)
      { for(i1=0;i1<n1;i1++)
        { for(i2=0;i2<n2;i2++)
          { dnum=d3aray[i1][i2][i3];
            CHRwsp(fptr); DBLw(fptr,dnum,ndigit);  counter++;
            if(counter>=nperli)
            { counter=0; CHRwcr(fptr);
    } } } } }    
    if(ihow==321)  
    { for(i3=0;i3<n3;i3++)
      { for(i2=0;i2<n2;i2++)
        { for(i1=0;i1<n1;i1++)
          { dnum=d3aray[i1][i2][i3];
            CHRwsp(fptr); DBLw(fptr,dnum,ndigit);  counter++;
            if(counter>=nperli)
            { counter=0; CHRwcr(fptr);
    } } } } }    
    if(counter>0)CHRwcr(fptr);
  }
  return fptr;
}
FILE *TwriteINT1(namep,fpt,ndigit,nperli,n1,i1aray,name)
  CHR *namep; FILE *fpt; YINT ndigit; YINT nperli; YINT n1;
  YINT  *i1aray; CHR *name;
{ YINT i1;
  YINT counter,istring;
  YINT inum;   
  FILE *fptr=fpt;  

  if(fptr==FILENULL)fptr=fopen(namep,"w");
  if((i1aray!=INT1NULL)&&(fptr!=FILENULL))
  { CHRw(fptr,name) ; CHRwcr(fptr);
    counter=0;
    for(i1=0;i1<n1;i1++)
    { counter=MAXIM(counter,ABS(i1aray[i1]));
    }
    istring=0;
    while(counter>0)
    { counter=counter/10;
      istring++;
    }
    istring=MAXIM(istring+1,ndigit);

    counter=0;   
    for(i1=0;i1<n1;i1++)
    { inum=i1aray[i1];
      CHRwsp(fptr); INTw(fptr,inum,istring);  counter++;
      if(counter>=nperli)
      {  counter=0; CHRwcr(fptr);
    } }
    if(counter>0)CHRwcr(fptr);
  }
  return fptr;
}

FILE *TwriteINT2(namep,fpt,ndigit,nperli,ihow,n1,n2,i2aray,name)
  CHR *namep; FILE *fpt; YINT ndigit; YINT nperli; YINT ihow; 
  YINT n1; YINT  n2;
  YINT **i2aray; CHR *name;
{ YINT i1,i2;
  YINT counter,istring;
  YINT inum; 
  FILE *fptr=fpt;  

  if(fptr==FILENULL)fptr=fopen(namep,"w");
  if((i2aray!=INT2NULL)&&(fptr!=FILENULL))
  { CHRw(fptr,name); INTw(fptr,ihow,8); 
    CHRwsp(fptr); INTw(fptr,n1,8);
    CHRwsp(fptr); INTw(fptr,n2,8);
    CHRwcr(fptr);

    counter=0;
    for(i1=0;i1<n1;i1++)
    { for(i2=0;i2<n2;i2++)
      { counter=MAXIM(counter,ABS(i2aray[i1][i2]));
    } }
    istring=0;
    while(counter>0)
    { counter=counter/10;
      istring++;
    }
    istring=MAXIM(istring+1,ndigit); 

    counter=0; 
    if(ihow==12)
    { counter=0;   
      for(i1=0;i1<n1;i1++)
      { for(i2=0;i2<n2;i2++)
        { inum=i2aray[i1][i2];
          CHRwsp(fptr); INTw(fptr,inum,istring);  counter++;
          if(counter>=nperli)
          {  counter=0; CHRwcr(fptr);
    } } } }
    else
    { counter=0;   
      for(i2=0;i2<n2;i2++)
      { for(i1=0;i1<n1;i1++)
        { inum=i2aray[i1][i2];
          CHRwsp(fptr); INTw(fptr,inum,istring);  counter++;
          if(counter>=nperli)
          {  counter=0; CHRwcr(fptr);
    } } } }
    if(counter>0)CHRwcr(fptr);
  }
  return fptr;
} 

FILE *TwriteINT3(namep,fpt,ndigit,nperli,ihow,n1,n2,n3,i3aray,name)
  CHR *namep; FILE *fpt; YINT ndigit; YINT nperli; YINT ihow; 
  YINT n1; YINT  n2; YINT n3;
  YINT ***i3aray; CHR *name;
{ YINT i1,i2,i3;
  YINT counter,istring;
  YINT inum; 
  FILE *fptr=fpt;  

  if(fptr==FILENULL)fptr=fopen(namep,"w");
  if((i3aray!=INT3NULL)&&(fptr!=FILENULL))
  { CHRw(fptr,name); INTw(fptr,ihow,8); 
    CHRwsp(fptr); INTw(fptr,n1,8);
    CHRwsp(fptr); INTw(fptr,n2,8);
    CHRwsp(fptr); INTw(fptr,n3,8);
    CHRwcr(fptr);

    counter=0;
    for(i1=0;i1<n1;i1++)
    { for(i2=0;i2<n2;i2++)
      { for(i3=0;i3<n3;i3++)
        { counter=MAXIM(counter,ABS(i3aray[i1][i2][i3]));
    } } }
    istring=0;
    while(counter>0)
    { counter=counter/10;
      istring++;
    }
    istring=MAXIM(istring+1,ndigit); 

    counter=0; 
    if(ihow==123)
    { counter=0;   
      for(i1=0;i1<n1;i1++)
      { for(i2=0;i2<n2;i2++)
        { for(i3=0;i3<n3;i3++)
          { inum=i3aray[i1][i2][i3];
            CHRwsp(fptr); INTw(fptr,inum,istring);  counter++;
            if(counter>=nperli)
            {  counter=0; CHRwcr(fptr);
    } } } } }
    else if(ihow==132)
    { counter=0;   
      for(i1=0;i1<n1;i1++)
      { for(i3=0;i3<n3;i3++)
        { for(i2=0;i2<n2;i2++)
          { inum=i3aray[i1][i2][i3];
            CHRwsp(fptr); INTw(fptr,inum,istring);  counter++;
            if(counter>=nperli)
            {  counter=0; CHRwcr(fptr);
    } } } } }
    else if(ihow==213)
    { counter=0;   
      for(i2=0;i2<n2;i2++)
      { for(i1=0;i1<n1;i1++)
        { for(i3=0;i3<n3;i3++)
          { inum=i3aray[i1][i2][i3];
            CHRwsp(fptr); INTw(fptr,inum,istring);  counter++;
            if(counter>=nperli)
            {  counter=0; CHRwcr(fptr);
    } } } } }
    else if(ihow==231)
    { counter=0;   
      for(i2=0;i2<n2;i2++)
      { for(i3=0;i3<n3;i3++)
        { for(i1=0;i1<n1;i1++)
          { inum=i3aray[i1][i2][i3];
            CHRwsp(fptr); INTw(fptr,inum,istring);  counter++;
            if(counter>=nperli)
            {  counter=0; CHRwcr(fptr);
    } } } } }
    else if(ihow==312)
    { counter=0;   
      for(i3=0;i3<n3;i3++)
      { for(i1=0;i1<n1;i1++)
        { for(i2=0;i2<n2;i2++)
          { inum=i3aray[i1][i2][i3];
            CHRwsp(fptr); INTw(fptr,inum,istring);  counter++;
            if(counter>=nperli)
            {  counter=0; CHRwcr(fptr);
    } } } } }
    else
    { counter=0;   
      for(i3=0;i3<n3;i3++)
      { for(i2=0;i2<n2;i2++)
        { for(i1=0;i1<n1;i1++)
          { inum=i3aray[i1][i2][i3];
            CHRwsp(fptr); INTw(fptr,inum,istring);  counter++;
            if(counter>=nperli)
            {  counter=0; CHRwcr(fptr);
    } } } } }    
    if(counter>0)CHRwcr(fptr);
  }
  return fptr;
}
/****************SORTING ARRAYS BY SIZE - smallest ... largest number ****/
void TsortDBL1(n,d1x,i1x)      /* Sorts double array smallest...largest */
  YINT n;  DBL *d1x;  YINT *i1x; /*rearanges i1x in the same way          */
{ YINT iblock,i,j,k,itmp;
  DBL tmp,m;
  YINT ibig[50];
  YINT iend[50];
  ibig[0]=0;
  iend[0]=n-1;
  iblock=0;
  while(iblock>=0)
  { i=ibig[iblock];
    j=iend[iblock];   
    iblock=iblock-1; 
    if(j>i)
    { tmp=d1x[j];
      m=tmp;
      for(k=i;k<j;k++)
      { tmp=MINIM(tmp,d1x[k]);
        m=MAXIM(m,d1x[k]);
      }
      if(tmp!=m)
      { m=(tmp+m)/2;
        while(i<=j)
        { while((i<=j)&&(d1x[i]<=m)) {i=i+1;}
          while((j>=i)&&(d1x[j]>m))  {j=j-1;}
          if(j>i)
          { tmp=d1x[j];
            d1x[j]=d1x[i];
            d1x[i]=tmp;
            if(i1x!=INT1NULL)
            { itmp=i1x[j];
              i1x[j]=i1x[i];
              i1x[i]=itmp;
        } } }
        ibig[iblock+2]=ibig[iblock+1];
        iend[iblock+2]=j;       
        ibig[iblock+1]=i;
        iblock=iblock+2;
} } } } 
 
void TsortINT(n,nsort,nrear,i2)      /* Sorts nsort arrays     - smallest...largest */
                                     /* and rearange nrear-nsort arrays  same way   */
  YINT n; YINT nsort; YINT nrear;  YINT **i2; 
{ YINT iblock,i,j,k,isort,tmp,m;
  YINT *i1,*i1r;
  YINT ibig[50];
  YINT iend[50];
  ibig[0]=0;
  iend[0]=n-1;
  iblock=0;
  while(iblock>=0)
  { i=ibig[iblock];
    j=iend[iblock];   
    iblock=iblock-1; 
    if(j>i)
    { tmp=0; m=0; isort=0;
      while((isort<nsort)&&(m==tmp))
      { i1=i2[isort];
        isort=isort+1;
        tmp=i1[j]; 
        m=tmp;
        for(k=i;k<j;k++)
        { tmp=MINIM(tmp,i1[k]);
          m=MAXIM(m,i1[k]);
      } }
      if(tmp!=m)
      { m=(tmp+m)/2;
        while(i<=j)
        { while((i<=j)&&(i1[i]<=m)) {i=i+1;}
          while((j>=i)&&(i1[j]>m))  {j=j-1;}
          if(j>i)
          { for(isort=0;isort<nrear;isort++)
            { i1r=i2[isort];
              tmp=i1r[j];
              i1r[j]=i1r[i];
              i1r[i]=tmp;
        } } }
        ibig[iblock+2]=ibig[iblock+1];
        iend[iblock+2]=j;       
        ibig[iblock+1]=i;
        iblock=iblock+2;
} } } }
 
/**************SAVING SPACE BY CODED ARRAYS**************/
#define ICODEBASE 90
static CHR *c1diga="0123456789-=qwertyuiop[]^asdfg";
static CHR *c1digb="hjkl;zxcvbnm,./~!@#$%&*()_+QWE";
static CHR *c1digc="RTYUIOP{}|ASDFGHJKL:ZXCVBNM<>?";
static CHR  c1dig[ICODEBASE+1];
static YINT  i1dig[1000];
static YINT  i1max[10];
static YINT iffirst=0;

static void initcode()
{ YINT i,j;

  iffirst=1;
  i1max[0]=1;
  for(i=1;i<10;i++)
  { i1max[i]=i1max[i-1]*ICODEBASE;
  }
  for(j=0;j<1000;j++)
  { i1dig[j]=0;
  }    
  CHRcpy(c1dig,c1diga);
  CHRcat(c1dig,c1digb);
  CHRcat(c1dig,c1digc);
  for(i=0;i<ICODEBASE;i++)
  { j=(YINT)c1dig[i];
    if(j<0)j=j+500;
    if(j>1000)
    { CHRw(stderr,"Unexpected digit value");
      exit(1);
    }
    i1dig[j]=i;
} }
  
void codeCHRtoINT(c1code,i1num)
   CHR *c1code; YINT *i1num;
{ YINT inum,idig,ndigit,nnum,icar;

  if(iffirst==0)initcode();
  icar=(YINT)c1code[0];
  if(icar<0)icar=icar+500;
  ndigit=i1dig[icar];
  icar=(YINT)c1code[1];
  if(icar<0)icar=icar+500;
  nnum=i1dig[icar];
  i1num[0]=ndigit;
  i1num[1]=nnum; 
  for(inum=2;inum<nnum;inum++)
  { i1num[inum]=0;
    for(idig=ndigit-1;idig>=0;idig--)
    { icar=(YINT)c1code[(inum-2)*ndigit+2+idig];
      if(icar<0)icar=icar+500;
      i1num[inum]=i1num[inum]+i1dig[icar]*i1max[idig];
} } }
 
void codeINTtoCHR(c1code,i1num)
   CHR *c1code; YINT *i1num;
{ YINT inum,idig,ival,ndigit,nnum;

  if(iffirst==0)initcode();
  ndigit=i1num[0];
  nnum=i1num[1];
  c1code[0]=c1dig[ndigit];
  c1code[1]=c1dig[nnum];
  c1code[(nnum-2)*ndigit+2]=CHRTERMINATE; 
  for(inum=2;inum<nnum;inum++)
  { ival=i1num[inum];
    for(idig=ndigit-1;idig>=0;idig--)
    { c1code[(inum-2)*ndigit+2+idig]=c1dig[(ival/i1max[idig])];
      ival=ival%i1max[idig];
} } }

void codeINTtoDBL(d1num,i1num)
  DBL *d1num; YINT *i1num;
{ YINT inum,nnum,ndigit;
  DBL dmax,dval;

  if(iffirst==0)initcode();
  ndigit=i1num[0];
  nnum=i1num[1];
  dmax=(DBL)(i1max[ndigit]-1);
  for(inum=0;inum<nnum;inum++)
  { dval=(DBL)i1num[inum];
    d1num[inum]=R2*dval/dmax-R1;
} }
 
void codeDBLtoINT(d1num,i1num)
  DBL *d1num; YINT *i1num;
{ YINT inum,nnum,ndigit;
  DBL dmax,dval;

  if(iffirst==0)initcode();
  ndigit=i1num[0];
  nnum=i1num[1];
  dmax=(DBL)(i1max[ndigit]-1);
  for(inum=2;inum<nnum;inum++)
  { dval=RP5*(d1num[inum]+R1)*dmax;
    i1num[inum]=(YINT)dval;
} }

void UnVec(DBL *e1x,DBL *e1y,DBL rx,DBL ry)
{ DBL mag;
  mag=SQRT((rx*rx)+(ry*ry));
  *e1x=rx/mag;
  *e1y=ry/mag;
}

void Rot90(DBL *e1x,DBL *e1y,DBL e2x,DBL e2y)
{ *e1x=e2y;
  *e1y=-e2x;
}

void V2DTranToLoc(DBL *u,DBL *v,DBL rx,DBL ry,DBL e1x,DBL e1y,DBL e2x,DBL e2y)
{ *u=(rx)*(e1x)+(ry)*(e1y);
  *v=(rx)*(e2x)+(ry)*(e2y);
}

