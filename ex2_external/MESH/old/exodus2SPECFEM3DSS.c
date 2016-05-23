/* This program converts the Binary (provided that ncdump command exists)
or ASCII exodus file exported from the CUBIT to several mesh files required
by the SPECFEM3D package. Basic steps starting from the CUBIT:
in CUBIT:
- blocks define only the material regions
  -> actual material properties should not be defined within the CUBIT.
  actual material properties can be listed later corresponding to each block 
  (i.e., material region).
- surface boundary conditions must be defined using side set
  -> each side set name must contain the corresponding BC names as defined in 
  char *ss_bcname[] below
  e.g., side set name can be front_ssbcux or front_ssbcux_ssbcuy etc.
- nodal boundary conditions must be defined using node set
  -> each node set name must contain the corresponding BC names as defined in 
  char *ns_bcname[] below
  e.g., node set name can be front_nsbcux or front_nsbcux_nsbcuy etc.

step1: export mesh file as exodus file say "tunnel.e" (use 3D option)
step2: convert "tunnel.e" to ASCII file using
  >>ncdump tunnel.e > tunnel.txt
step3: produce mesh and BC files
  >>exodus2SPECFEM3DSS tunnel.txt
  OR
  >>exodus2SPECFEM3DSS tunnel.txt 1000.0

There will be several output files:
*_nodes_coords_file : total number of nodes followed by nodal coordinate ? (? -> x, y, z)
*_mesh_file : total number of elements followed by connectivity list
*_materials_file : total number of elements followed by material IDs
*_??bcu? : node IDs which have u? = 0 as the boundary conditions (?? -> ns or ss, ? -> x, y, z) 
------------------------------------------------------
DEVELOPER:
  Hom Nath Gharti
  Princeton University
  homnath_AT_princeton_DOT_edu
  Hom Nath Gharti
  NORSAR
  homnath_AT_norsar_DOT_no
DEPENDENCY:
  stringmanip.c: string manipulation routines
COMPILE:
  gcc exodus2SPECFEM3DSS.c -o exodus2SPECFEM3DSS
USAGE: 
  exodus2SPECFEM3DSS <inputfile> <OPTIONS>
  Example: exodus2SPECFEM3DSS tunnel.txt
  or
  exodus2SPECFEM3DSS tunnel.e -fac=0.001 -bin=1
OPTIONS:
  -fac: use this option to multiply coordinates. this is importantn for unit 
        conversion, e.g., to convert m to km use -fac=0.001
  -bin: use this option if you want to convert exodus binary directly, provided
        that the command ncdump is in the path. ncdump is a part of netCDF library
        that can be downloaded freely from 
        http://www.unidata.ucar.edu/downloads/netcdf/index.jsp. use -bin=1 for binary 
        or -bin=0 for ascii file.
  -norm: use this option to check the normal of the faces. use -norm=1 for checking or -norm=0 (default) for no checking
HISTORY:
  HNG,Jan 12,2015: fixed hex-face ordering and tested for latest version of SPECFEM3D
  HNG,Jan 18,2011: modified to convert to SPECFEM3D
  HNG,Apr 23,2010;HNG,Apr 17,2010;HNG,Feb 08,2009
-------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "stringmanip.c"

#define OFF 0
#define ON 1

/* auxiliary routines */
void removeExtension(char *, char *);
int get_int(int *, char *, char *);
int look_int(int *, char *, char *);
int getfirstquote(char *, char *);


/* main routine */
int main(int argc,char **argv){
int i,ielmt,iface,inode,itmp,j,k;
int ndim;	/* geometry dimension */ 
int nnode,nelmt; /* number of nodes, number of elements */
int nblk,nss; /* number of blocks, number of side sets */
int elmt_count,node_count; /* element, node count */
int blk_count,ss_count; /* block, node set count */
int dim_stat,ss_stat,con_stat,coord_stat,mat_stat; /* status */
int *blk_nelmt,*blk_nenod; /* number of elements, number of nodes per element in each node */
int *ss_nside; /* number of nodes in each node set */

double fac,ftmp; /* multiplication factor for coordinates, temporary float */
double **coord;
int **elmt_node;
char *bulk,line[100],token[62],dumc[62],stag[62];
char **coord_name; /* coordinates name */
char fonly[62],infname[62],outfname[62];
char **ss_name; /* node set names */

/* change this line to look for other BCs in side sets */
double normal[3],p[3][4];
/* int fnod[6][4]={{1,2,6,5},{2,3,7,6},{3,4,8,7},{4,1,5,8},{2,1,4,3},{5,6,7,8}}; // face nodes */
int fnod[6][4]={{1,2,6,5},{2,3,7,6},{3,4,8,7},{1,5,8,4},{1,4,3,2},{5,6,7,8}}; // face nodes
int igfnod,gfnod[4]; // global face nodes
int nface_change;

int count_sselmt,count_ssside,sumss_nelmt,*ss_elmt,*ss_side; /* total number of element in absorbing BC and free surface */

int isbin; /* test if binary */
int ndir,isnorm; /* normal direction, test if normal has to checked */

FILE *inf,*outf_mat,*outf_con,*outf_coord,*outf_side;

/* default factor and binary switch*/    
fac=1.0; isbin=OFF; isnorm=OFF; ndir=1;

if(argc<2){
  fprintf(stderr,"ERROR: input file not entered!\n");
  exit(-1);
}

/* scan command */
if(argc>2){
  for(i=2;i<argc;i++){
    if(look_double(&ftmp,"-fac=",argv[i])==0){
	  fac=ftmp;		  
	  continue;
	}
	else if(look_int(&itmp,"-bin=",argv[i])==0){
	  isbin=itmp;
	  continue;
	}else if(look_int(&itmp,"-norm=",argv[i])==0){
	  isnorm=itmp;
	  continue;
	}else{
	  printf("ERROR: unrecognized option \"%s\"",argv[i]);
	  exit(-1);
	}
  }
}

printf("input file: %s\n",argv[1]);
printf("fac: %f\n",fac);
printf("isbin: %d\n",isbin);
printf("isnorm: %d\n",isnorm);

printf("--------------------------------\n");
/* default input file name is argv[1]*/
strcpy(infname,argv[1]);
removeExtension(argv[1],fonly);

if (isbin){
  printf("converting binary to ascii...");
  /* set input file name */
  sprintf(infname,"%s.txt",fonly);

  /* convert binary netCDF file to ascii file */
  sprintf(dumc,"ncdump %s > %s.txt",argv[1],fonly);  
  if (system(dumc)!=0){
	printf("ERROR: command \"%s\" cannot be executed! use -bin=0 or no option for ascii input file! \n",dumc);
	exit(-1);
  }
  printf("complete!\n");
}

/* open input file */
inf=fopen(infname,"r");
if(inf==NULL){
  fprintf(stderr,"ERROR: file \"%s\" not found!",argv[1]);
  exit(-1);
}
/*printf("--------------------------------\n");*/

bulk=malloc(1000); /* bulk string */
      
/* initialize some variables to 0 */
ndim=0; nblk=0; nnode=0; nelmt=0; nss=0;

/* intialize count to 0 */
blk_count=0; ss_count=0; node_count=0; elmt_count=0;	

/* set default status to OFF */
dim_stat=OFF; ss_stat=OFF; con_stat=OFF; coord_stat=OFF;

fscanf(inf,"%s",token);
if(strcmp(token,"netcdf")!=0){
  printf("ERROR: invalid exodus file or wrong -bin value!\n");
  printf("HINT: try correct value for -bin option or use valid exodus file!\n");
  exit(-1);
}

count_sselmt=0; count_ssside=0;
while(!feof(inf)){
  fscanf(inf,"%s",token);
  /* read dimensions */
  if(dim_stat!=ON && strcmp(token,"dimensions:")==0){
    printf("reading dimensions...");
    while(strstr(fgets(line,100,inf),"variables:") == NULL){		
      strncat(bulk,line,strcspn(line,";")+1);              
    }    
    get_int(&ndim,"num_dim =",bulk);
    if(ndim>0){	  
      /* allocate memory */
      coord_name=malloc(ndim*sizeof(char *));
      for(i=0;i<ndim;i++){	    
      coord_name[i]=malloc(62*sizeof(char)); /* each name has maximum of 62 characters */		
      }	  
    }else{
      printf("ERROR: illegal value of dimension!\n");
      exit(-1);
    }
    get_int(&nnode,"num_nodes =",bulk);
    get_int(&nelmt,"num_elem =",bulk);
    get_int(&nblk,"num_el_blk =",bulk);
    
    /* allocate memory */
    blk_nelmt=malloc(nblk*sizeof(int));
    blk_nenod=malloc(nblk*sizeof(int));
	elmt_node=malloc(8*sizeof(int *));
	if(elmt_node==NULL){
		printf("ERROR: not enough memory!");
		exit(-1);
	}
	for(i=0;i<8;i++){
		/*printf("%d %d\n",nelmt,i);*/
        elmt_node[i]=malloc(nelmt*sizeof(int)); /* need to make general */
		if(elmt_node[i]==NULL){
			printf("ERROR: not enough memory!");
			exit(-1);
		}
		/*printf("%d %d\n",nelmt,i);*/
    }
	coord=malloc(3*sizeof(double *));
	for(i=0;i<3;i++){
        coord[i]=malloc(nnode*sizeof(double));
		if(coord[i]==NULL){
			printf("ERROR: not enough memory!");
			exit(-1);
		}
      }

    /* sideset information */
    if (look_int(&nss,"num_side_sets =",bulk)!=0){
      nss=0;
    }else{
      /* allocate memory */
      ss_name=malloc(nss*sizeof(char *));
      for(i=0;i<nss;i++){
        ss_name[i]=malloc(62*sizeof(char)); /* each name has maximum of 62 characters */
      }
      ss_nside=malloc(nss*sizeof(int));
    }
    
    if(nss>0){ /* This segment has a significance only if nss has legitimate value */
      for(i=0;i<nss;i++){
        sprintf(stag,"num_side_ss%d =",i+1);
        get_int(&ss_nside[i],stag,bulk);					
      }
    }

    /* block information */        
    for(i=0;i<nblk;i++){
      sprintf(stag,"num_el_in_blk%d =",i+1);
      get_int(&blk_nelmt[i],stag,bulk); 
      
      sprintf(stag,"num_nod_per_el%d =",i+1);
      get_int(&blk_nenod[i],stag,bulk);
    }    
    	
    dim_stat=ON;
    free(bulk);
    printf("complete!\n");
    printf(" geometry dimension: %d\n",ndim);    
    printf(" number of nodes: %d\n",nnode);
    printf(" number of elements: %d\n",nelmt);
    printf(" number of blocks: %d\n",nblk);
    continue;
  }  
  
  /* read coordinate names */
  if(strcmp(token,"coor_names")==0){    
    fscanf(inf,"%s",dumc); /* = */			
    for (i=0; i<ndim; i++){
      fscanf(inf,"%s",dumc);	  
	  getfirstquote(dumc,coord_name[i]);	  
    }	
    continue;
  }  
  
  /* read and store side boundary conditions */
  if(strcmp(token,"ss_names")==0){
    //printf("saving side BCs...");
    fscanf(inf,"%s",dumc); /* = */			
    for (i=0; i<nss; i++){  
	  fscanf(inf,"%s",dumc);	  
	  getfirstquote(dumc,ss_name[i]);
	}
	/* total number of sideset elements */
	sumss_nelmt=0;
	for(i=0;i<nss;i++)sumss_nelmt+=ss_nside[i];	

	ss_elmt=malloc(sumss_nelmt*sizeof(int));
	ss_side=malloc(sumss_nelmt*sizeof(int)); 
    continue;
  }
 
  if(ss_stat!=ON){
    for(i=0;i<nss;i++){
      sprintf(stag,"elem_ss%d",i+1);
      if(strcmp(token,stag)==0){

        fscanf(inf,"%s",dumc); /* = */
        for(j=0;j<ss_nside[i]; j++){
          fscanf(inf,"%d,",&ss_elmt[count_sselmt]); /* read comma separated data */
		  count_sselmt++; 
        }
        
		fscanf(inf,"%s",dumc); /* ; */
		fscanf(inf,"%s",token);
	    sprintf(stag,"side_ss%d",i+1);
	    if(strcmp(token,stag)==0){
	      fscanf(inf,"%s",dumc); /* = */
	      for(j=0;j<ss_nside[i]; j++){
	        fscanf(inf,"%d,",&ss_side[count_ssside]); /* read comma separated data */
			count_ssside++;
	      }
	    } 
        continue;
      }          
    }
  }
 
  /* Connectivity */
  if(nblk>0 && con_stat!=ON){   
   
    /* write connectivity and material id */  
    for(i=0;i<nblk;i++){
      sprintf(stag,"connect%d",i+1);
      if(strcmp(token,stag)==0){
        blk_count++;
        
        /* open connectivity and material files */
        if(blk_count==1){          
          printf("writing connectivity and materials..."); 
          sprintf(outfname,"%s_mesh_file",fonly);
          outf_con=fopen(outfname,"w");
          fprintf(outf_con,"%d\n",nelmt);
          
          sprintf(outfname,"%s_materials_file",fonly);
          outf_mat=fopen(outfname,"w");
          //fprintf(outf_mat,"%d\n",nelmt);
        }
        
        fscanf(inf,"%s",dumc); /* = */
        for(j=0;j<blk_nelmt[i];j++){
	      fprintf(outf_con,"%d ",elmt_count+1);
          for(k=0;k<blk_nenod[i];k++){
            fscanf(inf,"%d,",&itmp);
            fprintf(outf_con,"%d ",itmp);
			elmt_node[k][elmt_count]=itmp;
          }
		  fprintf(outf_con,"\n");/* new line */
          fprintf(outf_mat,"%d %d\n",elmt_count+1,i+1);
          elmt_count++;
          
        }        
        
        if(blk_count==nblk){
          con_stat=ON;
          mat_stat=ON;
          free(blk_nelmt);
          free(blk_nenod);
          printf("complete!\n");
          fclose(outf_con);
          fclose(outf_mat);
        }
        continue;
      }
    }
  }
  
  /* Coordinates */
  if(strcmp(token,"coord")==0){
    printf("reading coordinates...");
    fscanf(inf,"%s",dumc);
    for(i=0;i<ndim;i++){ 
      for(j=0;j<nnode;j++){
        fscanf(inf,"%lf,",&ftmp); /* read comma separated data */ 
		/*printf("%d %d\n",nnode,node_count);*/
		//printf("%d %d %d\n",i,j,nnode);
		coord[i][j]=fac*ftmp;
        if(i==0)node_count++;
      }
    }
	coord_stat=ON;
	printf("complete!\n");
	break;
  }  
}
fclose(inf);

/* write coordinates file */
printf("writing coordinates...");
sprintf(outfname,"%s_nodes_coords_file",fonly);
outf_coord=fopen(outfname,"w");
fprintf(outf_coord,"%d\n",nnode);
for(i=0;i<nnode;i++)fprintf(outf_coord,"%d %.15lf %.15lf %.15lf\n",i+1,coord[0][i],coord[1][i],coord[2][i]);
fclose(outf_coord);
printf("complete!\n");


/* Write sidesets */
if(ss_stat=ON){
	printf("writing side sets...");
	count_sselmt=0;
	nface_change=0;
	for(i=0;i<nss;i++){
		ndir=1;
		sprintf(outfname,"%s_%s",fonly,ss_name[i]);
		outf_side=fopen(outfname,"w");
		fprintf(outf_side,"%d\n",ss_nside[i]);
		/* assign normal */
		if (strstr(ss_name[i],"xmin")!=NULL){
			normal[0]=-1.0; normal[1]=0.0; normal[2]=0.0;
		}else if(strstr(ss_name[i],"xmax")!=NULL){
			normal[0]=1.0; normal[1]=0.0; normal[2]=0.0;
		}else if(strstr(ss_name[i],"ymin")!=NULL){
			normal[0]=0.0; normal[1]=-1.0; normal[2]=0.0;
		}else if(strstr(ss_name[i],"ymax")!=NULL){
			normal[0]=0.0; normal[1]=1.0; normal[2]=0.0;
		}else if(strstr(ss_name[i],"zmin")!=NULL || strstr(ss_name[i],"bottom")!=NULL){
			normal[0]=0.0; normal[1]=0.0; normal[2]=-1.0;
		}else if(strstr(ss_name[i],"zmax")!=NULL || strstr(ss_name[i],"top")!=NULL){
			normal[0]=0.0; normal[1]=0.0; normal[2]=1.0;
		}else{
			/* do not check normal for unknown-normal faces */
			isnorm=0;
			ndir=1;
		}
		for(j=0;j<ss_nside[i];j++){
			fprintf(outf_side,"%d ",ss_elmt[count_sselmt]);
			ielmt=ss_elmt[count_sselmt]-1; /* index starts from 0!*/
			iface=ss_side[count_sselmt]-1; /* index starts from 0!*/			
			for(k=0;k<4;k++){
			  gfnod[k]=elmt_node[fnod[iface][k]-1][ielmt]; /* index starts from 0!*/
			  igfnod=gfnod[k]-1; /* index starts from 0!*/
			  p[0][k]=coord[0][igfnod]; p[1][k]=coord[1][igfnod]; p[2][k]=coord[2][igfnod];
			}
			if(isnorm){
				ndir=check_normal(p,normal);
			}
			//printf("Hello %d %d %d %d\n",gfnod[0],gfnod[1],gfnod[2],gfnod[3]);
			if(ndir==1){
				fprintf(outf_side,"%d %d %d %d\n",gfnod[0],gfnod[1],gfnod[2],gfnod[3]);
			}else if(ndir==-1){
				fprintf(outf_side,"%d %d %d %d\n",gfnod[0],gfnod[3],gfnod[2],gfnod[1]);
				nface_change++;
			}else{
				printf("ERROR: absurd value of ndir:%d for normal orientation!\n",ndir);
				exit(-1);
			}
			count_sselmt++;			
		}
		fclose(outf_side);
		ss_stat=ON;
	}
	printf("complete!\n");
	printf(" number of side sets written: %d\n",nss);
	printf(" number of reversed normals: %d\n",nface_change);
}
 
for(i=0;i<3;i++)free(coord[i]);
free(coord);
for(i=0;i<8;i++)free(elmt_node[i]);
free(elmt_node);
free(ss_elmt);
free(ss_side);
/* check status */
if(nnode!=node_count){
  printf("ERROR: number of nodes inconsistent!\n");
  exit(-1);
}
if(nelmt!=elmt_count){
  printf("ERROR: number of elements inconsistent!\n");
  exit(-1);
}
if(dim_stat!=ON){
  printf("ERROR: dimensions cannot be read!\n");
  exit(-1);
}

if(ss_stat!=ON){
  printf("WARNING: side boundary conditions cannot be read!\n");
  exit(-1);
}
if(con_stat!=ON){
  printf("ERROR: connectivity cannot be read!\n");
  exit(-1);
}
if(mat_stat!=ON){
  printf("ERROR: material IDs cannot be read!\n");
  exit(-1);
}
if(coord_stat!=ON){
  printf("ERROR: coordinates cannot be read!\n");
  exit(-1);
}
printf("--------------------------------\n");
return(0);
}
/*======================================*/

/* this function check the direction of normal to the given normal
and returns 1 if the direction matches, returns -1 if does not. */
int check_normal(double p[3][4],double normal[3])
{	
	int i;
	double a[3],b[3],comp_normal[3],dot;

	for(i=0;i<3;i++){
		/* compute 2 vectors */
		a[i]=p[i][1]-p[i][0];
		b[i]=p[i][2]-p[i][1];
	}
	/* compute normal */
	comp_normal[0]=(a[1]*b[2]-a[2]*b[1]);
	comp_normal[1]=(a[2]*b[0]-a[0]*b[2]);
	comp_normal[2]=(a[0]*b[1]-a[1]*b[0]);
	dot=0.0;
	for(i=0;i<3;i++)dot=dot+comp_normal[i]*normal[i];

	if(dot>0.0){
		/* directions match */
		return(1);
	}else if(dot<0.0){
		/* directions do not match */
		return(-1);
	}else{
		printf("ERROR: perpendicular normal orientation !\n");
		exit(-1);
	}
}
