
/* File   Y.c */
#include "Yproto.h"
#include <time.h>
#include "paraview.h" //! PARAVIEW

int main(argc, argv)
  YINT argc; char **argv;
{ CHR c1name[300];         /* name of the problem i.e. input file */
  struct YD_struct yd;     /* Y database                          */
  YDC ydc=&(yd.ydc);       /* Y control database                  */
  YDE yde=&(yd.yde);       /* Y element database                  */
  YDI ydi=&(yd.ydi);       /* Y interaction database              */
  YDN ydn=&(yd.ydn);       /* Y node database                     */
YDMN ydmn=&(yd.ydmn);         /* = nodal database                      */
  YDB ydb=&(yd.ydb);       /* Y borehole database                 */
  YDS yds=&(yd.yds);       /* Y source (inter. fluid) database    */
  YDO ydo=&(yd.ydo);       /* Y output database                   */
  YDPE ydpe=&(yd.ydpe);    /* Y property database  for elements   */
  YDPN ydpn=&(yd.ydpn);    /* Y property database  for nodes (BC) */
  YDPJ ydpj=&(yd.ydpj);    /* Y property database  for joints     */
  YDPM ydpm=&(yd.ydpm);    /* Y property database  for meshing    */
  YDFN ydfn=&(yd.ydfn);    /* Y property database for DFN         */
  YDIS ydis=&(yd.ydis);    /* Y in-situ stress parameters         */
  YDHF ydhf=&(yd.ydhf);    /* Y hydro-frac parameters             */
  YDSM ydsm=&(yd.ydsm);    /* Y seismic monitoring parameters     */
  YDR ydr=&(yd.ydr);       /* Y reference points database Y-RC    */
  YDSB ydsb=&(yd.ydsb);    /* Y steel-bar element database Y-RC   */
  YDPS ydps=&(yd.ydps);    /* Y property database  for steel Y-RC */
  
  time_t beginTime, endTime;
  time(&beginTime);
  
  /* get name of the problem */
  if(argv[1]!=NULL)
  { CHRcpy(c1name,argv[1]); 

/* process data */
    CHRw(stdout,"doi: http://dx.doi.org/10.1061/(ASCE)GM.1943-5622.0000216 \n\n");
//Mahabadi, O., Lisjak, A., Munjiza, A., and Grasselli, G. (2012). 
//Y-Geo: New Combined Finite-Discrete Element Numerical Code for Geomechanical Applications.
//International Journal of Geomechanics. 12, SPECIAL ISSUE: Advances in Modeling Rock Engineering Problems, 
//(doi: http://dx.doi.org/10.1061/(ASCE)GM.1943-5622.0000216)


    ydc->finp=FILENULL; ydc->fcheck=FILENULL;
    while(Yrd(c1name,&yd)>0)                 /* Process while any input */
    { CHRw(stdout,"NEW INPUT: "); CHRw(stdout, c1name); CHRwcr(stdout);
      for(ydc->ncstep=ydc->ncstep;ydc->ncstep<ydc->mcstep;ydc->ncstep++)
      { if(ydc->ncstep==0)Ymd(ydc,yde,ydi,ydn,ydmn,ydpe,ydpn,ydpm,ydfn,ydpj);                      /* mesh elements            */
      Ymdupdate(ydc,yde,ydi,ydn,ydmn,ydpe,ydpn,ydpm,ydfn,ydpj);
        Yfd(ydc,yde,ydn,ydmn,ydo,ydpe,ydpn, ydpj,ydis,ydfn,ydhf,ydsm,ydsb); /* nodal forces             */
        
        Ycd(ydc,yde,ydi,ydn,ydpe,ydpn);                                /* contact detection        */
        //Yrb(ydc,yde,ydn,ydpe,ydpj,ydps,ydr,ydsb,ydo);                  /* RBAR steel elements  Y-RC*/
        Yid(ydc,yde,ydi,ydn,ydo,ydpe,ydpn, ydpj,ydpm,ydfn);            /* interaction              */
         Ysd(ydc,yde,ydn,ydmn,ydo,ydpe,ydpn);                                /* solve equations          */
	    Yod(c1name,&yd);                                               /* output results           */
        Yfrd(ydc,yde,ydi,ydn,ydpe,ydpn,ydpj,ydpm);                     /* fracture                 */

        ydc->dctime=ydc->dctime+ydc->dcstec;                           /* update time              */
      }
    }
    pv_free();
    time(&endTime);
    printf("\nStart time: ");
    printf(ctime(&beginTime));
    printf("End time:   ");
    printf(ctime(&endTime));
    printf("Total run time (seconds) = %.1lf \n\n", difftime(endTime, beginTime));
    CHRw(stderr,"   ***** Y HAS ORDERLY FINISHED *****");
    CHRwcr(stderr);
    CHRw(stderr,"Press a key to continue");
    CHRwcr(stderr);
    getchar();
  }
  else
  { CHRw(stdout,"Double click the data file to run the program.\n");
  } 
}
