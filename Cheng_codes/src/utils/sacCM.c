/************************************************************
*	Cut and Merge  SAC files
*	Usage
*		sacCM o yyyy:mm:dd:hh:MM:ss.s b 5 e 30 f sac_files ...
*
*       History
*             2008  Weijun Wang, 
*sacCM 2008:05:27:17:40:00.0 2 20 f YUB.CQ.2008148173357.BHZ > a.sac
*  gcc -m32  sacCM.c time2double.c -g -o ~/seisbin/sacCM -lm
*************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "SacHeaderFile.h"


FILE * read_sachead2(const char *, struct SACheader *);
double time2double(char *);
int main(int argc, char **argv)
{
	int		n1;
	struct SACheader	hd,hdpre=NullSacHeader;
//	putenv("TZ=GMT+0");
	if(argc < 6) {
		fprintf(stderr,"Usage: sacCM yyyy:mm:dd:hh:MM:ss.s b e f file_lists\n");
		fprintf(stderr,"   ex. sacCM 2008:05:12:6:28.10 -300 18000 f YUB.BHZ.*.sac\n");
		fprintf(stderr,"   ex. sacCM 2008:133:6:28.10.j.m -300 18000 f YUB.BHZ.*.sac; j means julian day,.m means msec\n");
		fprintf(stderr,"       Weijun Wang\n"); 
		return -1;
	}
	double odt=0.0;	
	double b=-300, e=18000;
	n1=0; argv++; argc--;
	while ( *argv[0] != 'f' ) {
		if(n1==0) {
			odt=time2double(argv[0]);
//			fprintf(stderr,"odt=%f\n",odt);
		}
		else if(n1==1)
			b=atof(argv[0]);
		else if(n1==2)
			e=atof(argv[0]);
		n1++; argv++; argc--;
	}
	double cutbeg,cutend,sacbeg,sacend;
	
//	fprintf(stderr,"cutbeg=%f cutend=%f\n",cutbeg,cutend);
	int i;
	float *data;
	FILE *sacf;
	for(i = 1; i < argc; i++) {
	 	if( (sacf = read_sachead2(argv[i], &hd)) != NULL ){
			if(i==1){
				
				struct tm *tgmt;
				time_t todt = (time_t)odt;
				tgmt = gmtime(&todt);
				//hdpre =  NullSacHeader;
				hdpre.nzyear = tgmt->tm_year + 1900; 
				hdpre.nzjday = tgmt->tm_yday+1;
				hdpre.nzhour = tgmt->tm_hour;
				hdpre.nzmin = tgmt->tm_min ;
				hdpre.nzsec = tgmt->tm_sec;
				double tmp= odt*1000 - floor(odt)*1000.0;
				hdpre.nzmsec = hd.nzmsec; // set the nzmsec same as  will 
				odt = odt  - tmp/1000.0 + hd.nzmsec/1000.0;
				cutbeg = odt + b;
				cutend = odt + e;
				data=calloc(((cutend-cutbeg)/hd.delta+1), sizeof(float))	;
				if(!data) {fprintf(stderr,"Can't calloc memory for data\n");return(-1);}

//fprintf(stderr,"odt=%f,nzmsec=%d\n",odt,hdpre.nzmsec);
				hdpre.b = b;
				hdpre.e = e;
				hdpre.delta = hd.delta;
				hdpre.nvhdr = hd.nvhdr;
				hdpre.iftype = hd.iftype;
				hdpre.leven = hd.leven;
				hdpre.npts = (e-b)/hdpre.delta + 1;
//fprintf(stderr,"b=%f,e=%f,npts=%ld,%ld\n",hdpre.b,hdpre.e,hdpre.npts,length(data));
				strcpy(hdpre.kcmpnm,hd.kcmpnm);
				strcpy(hdpre.kstnm,hd.kstnm);
				strcpy(hdpre.knetwk,hd.knetwk);
				hdpre.stla = hd.stla;
				hdpre.stlo = hd.stlo;
				hdpre.stel = hd.stel;
				hdpre.evla = hd.evla;
				hdpre.evlo = hd.evlo;
				hdpre.evdp = hd.evdp;
			}
			if( strcmp(hd.kcmpnm,hdpre.kcmpnm)==0 && strcmp(hd.kstnm,hdpre.kstnm)==0 ){				
				char ts[40];
				sprintf(ts, "%d:%d:%d:%d:%d.%d.j",hd.nzyear,hd.nzjday,hd.nzhour,hd.nzmin,hd.nzsec,hd.nzmsec);
//				fprintf(stderr,"%s\n",ts);
				double reft = time2double(ts);
				sacbeg = reft + hd.b;
				sacend = reft + hd.e;
//fprintf(stderr,"sacbeg=%f, sacend=%f,npts=%f\n",sacbeg,sacend,(sacend-sacbeg)/hd.delta+1);
				double from, to;
				long fromp,locatep,top;
				if(cutbeg>=sacbeg && sacend >= cutbeg){
					// have data to cut
						from = cutbeg   ;
						to = sacend > cutend ? cutend : sacend ;
						
				}else if( sacbeg >= cutbeg && sacbeg < sacend){
						from = sacbeg;
						to = sacend > cutend ? cutend:sacend;
				}else 
					continue;
				fromp = ((from - sacbeg)/hd.delta +0.5);
				locatep = (from - cutbeg)/hd.delta +0.5;

				top  =  (to - from )/hd.delta+ 0.5+1;
//fprintf(stderr,"fromp-f:%f,locatep-f:%f,top-f:%f\n",(from - sacbeg)/hd.delta,(from - cutbeg)/hd.delta,(to - from )/hd.delta);
//fprintf(stderr,"sacf position:%ld\n",ftell(sacf));			
				fseek(sacf,fromp*sizeof(float),SEEK_CUR);
//fprintf(stderr,"sacf position:%ld\n",ftell(sacf));
				fread(data+locatep,sizeof(float),top,sacf);
//fprintf(stderr,"sacf position:%ld\n",ftell(sacf));
//fprintf(stderr,"from=%f,to=%f,fromp=%ld,locatep=%ld,top=%ld\n",from,to,fromp,locatep,top);
				}
		  } //end of if i==1
	
	}//end for all sac files
//fprintf(stderr,"stdout position:%ld\n",ftell(stdout));
	fwrite(&hdpre,sizeof(struct SACheader),1,stdout);
//fprintf(stderr,"stdout position:%ld\n",ftell(stdout));
//	fprintf(stderr,"npts:%ld\n",hdpre.npts);
	fwrite(data,sizeof(float),hdpre.npts,stdout);
//fprintf(stderr,"stdout position:%ld\n",ftell(stdout));

}  

FILE* 
read_sachead2(const char *name, struct SACheader *hd ) {
  FILE		*strm;
  
  if ((strm = fopen(name, "rb")) == NULL) {
    fprintf(stderr, "Unable to open %s\n",name);
    return NULL;
  }
  
  if (fread(hd, sizeof(struct SACheader), 1, strm) != 1) {
    fprintf(stderr, " Error in reading SAC header %s\n",name);
    fclose(strm);
    return NULL;
  }
  if(  hd->nvhdr != 6  ) {
      fprintf(stderr, " Error in SAC byteOrder or version %s\n",name);
      fclose(strm);
      hd = NULL;
      return(NULL);
    }
  
  // fclose(strm);
  return strm;
}


