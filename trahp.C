/*Author LUKE KRISTOPHER DAVIS

This code was originally written for MPHYS project at Swansea University. 

PARALLEL REPLICA EXCHANGE VERSION
*/

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <mpi.h>

using namespace std;

int ranseq(int N);
int buddycheck(int latlength, int mono1, int mono2);
double rndnum(void);
void  rini(int iran);
void placelattice(int N,int latlength,int size);
void oldposition(int N);
void pullmove(int N, int seqmono,int latlength, int size);
void pivot(int N, int axismono,int latlength, int size);
void kinkflip(int N,int latlength, int size);
void chainterminal(int N, int latlength);
void frw(int N, int latlength);
void typeii(int N,int latlength); 
void printposarr(int N);
void internalenergy(int hsize,int N);
void heatcapacity(int hsize,int N);
int compenergy(int N, int latlength);
int uppermove(int posmono2,int latlength);
int lowermove(int posmono2, int latlength);
int rightmove(int posmono2, int latlength);
int leftmove(int posmono2, int latlength);
int maxi(int a,int b);
int mini(int a,int b);
int checkthis(int N,int checkthis);
int row(int pos,int latlength);
int exclvol(int N);
int exclvol2(int N, int checkme, int upperlim, int lowerlim);
int col(int pos,int latlength);
int histoflat(int hsize, float p);
int SAW(int N,int latlength);

FILE *ifp, *ofp, *ofp2;
int numprocs,myid,visited[1],H[1],old[100],POS[100], HP[100],lowE=0, i;
long double logG[1],expo,T,Tlim=1, Tinc=0.001,Tstart=0.000001,lambda=0;
clock_t start,finish;
float p=0.7,elapsed_time;
long double k_b= 1,e=2.71828;


int main(int argc,char *argv[])
{
 int *sendbuf,*recvbuf;
 long int numsecs=86400,enesize,possize,seedid,hsize=1,saw, mv,k,y,x,mon,i,N=100,hp,latlength, size, seed=357913;
  long double *sendbufD,*recvbufD,ratio,f=1,test,minilog=0;
  enesize=hsize+1;possize=N+1;

  std::cout.precision(10);

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  seedid= (seed+myid)*2*(myid+1);
  rini(seedid);//initialize the seed
  latlength = 1000; size = latlength*latlength;

  if(myid==0)
    {
  ifp= fopen("hpseq2D100b","r");
  for(i=1; i<=N; i++)//writes seqeunce to an array which stores the sequence
    {fscanf(ifp,"%d\n", &HP[i]);/*printf("\nMonomer %d has HP value = %d\n",i,HP[i]);*/}
  fclose(ifp);

  for(i=0;i<=hsize;i++)
    {visited[i]=0;}
    }
  else;
  
  MPI_Bcast(HP,possize,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(visited,enesize, MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&numsecs,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  
  //Places initial position of the linear chain on the lattice
  placelattice(N,latlength,size);
  printf("Hello from %d my seed is %d\n", myid,seedid);
  
  int source,dest,mlog,brcnt=0,sizecnt,oldenergy=0, energy=0, minenergy=0,nosaw=0, count=0,F=0,inter=0;
  float ran,ki,pu;

  for(i=0; i<=hsize; i++)//initialise the DOS and histogram
    {logG[i]=0; H[i]=0;}  
  start=clock();
  numsecs=(long int)3120*2*numsecs;//1 days
  for(mv=1;mv<numsecs;mv++)
    {
      brcnt=0;
      
      for(i=0;i<=hsize;i++)//
	{if(i==0)
	    {mlog=logG[i];}
	  else
	    {if(mlog>logG[i])
		mlog=logG[i];
	      else;}
	}
      // MPI_Barrier(MPI_COMM_WORLD);
      //======================= REPLICA EXCHANGING ===========
      if(mv%10000==0)
	{
	 
	  for(i=0;i<=(numprocs-1);i++)
	    {if(myid==0)//if I am master thread
		{source=(int)(rndnum()*(numprocs-1));//printf("source= %d\n",source);
		  dest=i;//printf("dest= %d\n",dest);
		}
	      else;
	      MPI_Barrier(MPI_COMM_WORLD);
	      MPI_Bcast(&source,1,MPI_INT,0,MPI_COMM_WORLD);
	      MPI_Bcast(&dest,1,MPI_INT,0,MPI_COMM_WORLD);
	      if(source != dest)
		{
		  if(myid==source)
		    {
		      MPI_Send(POS,possize,MPI_INT,dest,1,MPI_COMM_WORLD);//printf("sending to %d\n",dest);
		    }
		  else if(myid==dest)
		    {
		      MPI_Recv(POS,possize,MPI_INT,source,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);//printf("receiving from %d\n",source);
		    }
		  else;
		}
	      else;
	    }
	}
      else;
      //======================================================
      
      oldenergy= compenergy(N,latlength);
      //anneal=-(int)(1*hsize);
      
      {  ran=rndnum();
	if(ran>=0.30)
	  {oldposition(N);
	    pullmove(N,ranseq(N),latlength,size);
	    // frw(N,latlength);
	      //      printf("PULL\n");
	      //  printposarr(N);
	    }
	  else if(ran<0.30 && ran>=0.25)
	    {oldposition(N);
	      // kinkflip(N,latlength,size);
	      frw(N,latlength);
	      // printf("FRW \n");
	      // printposarr(N);
	    }
	  else if(ran<0.25 && ran>=0.23)
	    {oldposition(N);
	      kinkflip(N,latlength,size);
	    }
	  else if(ran<0.23 && ran>0.04)
	    {
	      oldposition(N);
	      ran=rndnum();
	      if(ran>0.70)
		{chainterminal(N,latlength);}
	      else
		{typeii(N,latlength);} 
	    }
	  else
	    {oldposition(N);
	      pivot(N,ranseq(N),latlength,size);
	    }
	}
      
	  energy = compenergy(N,latlength);
      
	  nosaw=exclvol(N);
	  
	  saw=SAW(N,latlength)+nosaw;//check connections between neighbours to ensure they are okay
	  
	  if(saw==0)//if it is self avoiding
	    {
	      oldenergy=oldenergy*-1; energy=energy*-1;
	      
	      ratio=expl(logG[oldenergy]-logG[energy]);
	      //std::cout << "ratio=" << ratio <<'\n';
	      ran=rndnum();
	      if(ratio > (long double)ran)
		{;}//accept the new position array values
	      else 
		{//go back to old positions
		  for(i=1;i<=N;i++)
		    {POS[i]= old[i];}
		}
	      energy=(-1)*compenergy(N,latlength);//the resulting energy
	      if(visited[energy]==0)//i.e. new
		{
		  logG[energy]=mlog;
		  for(i=0;i<=hsize;i++)//reset the histogram
		    {H[i]=0;}
		  H[energy]=H[energy]+1;
		  logG[energy]=logG[energy]+f;
		  // continue;
		}
	      else;
	      visited[energy]=1;
	      H[energy]=H[energy]+1;//update the histogram
	      // g[energy]=g[energy]*f;//update DOS
	      logG[energy]=logG[energy]+f;
	      energy=(-1)*energy;
	      if(energy < minenergy)//is it a new minimum?
		{
		  minenergy =energy; 
		  
		}
	      else;
	      if(count==0)
		{
		  if(mv%100000==0)//checking flatness every 10^5
		    {
		      if(histoflat(hsize,p)==0)//it isn't flat enough yet
			{;}
		      else
			{f=f/2;
			  // std::cout << "modi f=" << f << '\n' << std::scientific;
			  for(i=0;i<=hsize;i++)//reset the histogram
			    {H[i]=0;}
			  //F=log(f);
			  inter=((long double)(hsize*N))/(long double)mv;//monte carlo time
			  if(f <= inter)
			    {f=inter;count++;}
			    else;
			}
		    }
		  else;
		}
	      else
		{inter=((long double)hsize)/(long double)mv;
		  f= inter;
		  std::cout << "modi f=" << f << '\n' << std::scientific;
		}
	      
	    }
	  else
	    {
	      for(i=1;i<=N;i++)
		{POS[i]=old[i];}
	      oldenergy=oldenergy*-1;
	      H[oldenergy]++;
	      logG[oldenergy]=logG[oldenergy]+f;
	    }
	  energy=0, oldenergy=0;
	    
	  if(myid==0 && mv%100000==0)
	    {std::cout << "Approximate time passed in MC iters = " << mv << " out of a total:"<< numsecs << '\n' << std::scientific;
	      std::cout << "Minenergy = " << minenergy << '\n' << std::scientific;
	    }
	  else;
	  
	  finish=clock()-start;
	  elapsed_time=((float)finish/((float)CLOCKS_PER_SEC));
    }
  printf("I am process %d and my minene = %d\n",myid,minenergy);
  sendbuf = new int; recvbuf = new int;
  sendbuf=&minenergy;
    MPI_Reduce(sendbuf,recvbuf,1,MPI_INT,MPI_MIN,0,MPI_COMM_WORLD);
    //printf("Minenergy= %d from process: %d\n",minenergy,myid);
    std::cout << "final modi f=" << f << "\t from process:" << myid <<  '\n' << std::scientific;
    if(myid==0)
      {
    std::cout << "From all processes min energy =" << *recvbuf << '\n' <<std::scientific;
    printf("Time taken: %f\n",elapsed_time);
      }
    else;
    
    //      heatcapacity(hsize,N);
    //  internalenergy(hsize,N);
    
    //  MPI_Barrier(MPI_COMM_WORLD);  
  MPI_Finalize();
  return 0;
}

int exclvol(int N)
{int y,n=0;
  for(i=1;i<=N;i++)//check self avoidance (not >1 bead per site
    {for(y=1;y<=N;y++)
	{if(i==y)
	    {;}
	  else
	    {
	      if(POS[i]==POS[y] || POS[i]-POS[y]==0)
		{n++;}
	      else{;}
	    }
	}
    }
  if(n==0)
    {return 0;}
  else
    {return 1;}
}

int exclvol2(int N,int checkme,int upperlim, int lowerlim)
{int j,n=0;
  if(lowerlim==1)
    {
      for(j=lowerlim;j<=upperlim;j++)
	{
	  if(POS[j]==checkme || POS[j]-checkme==0)
	    {n++;}
	  else{;}
	}
      
      if(n==0)
	{return 0;}
      else
	{return 1;}
    }
  else
    {
      for(j=lowerlim;j>=upperlim;j--)
	{
	  if(POS[j]==checkme || POS[j]-checkme==0)
	    {n++;}
	  else{;}
	}
      
      if(n==0)
	{return 0;}
      else
	{return 1;}
    }
}
void placelattice(int N, int latlength,int size)
{int i, spos;
  spos=(size/2)+(latlength/2);
  for(i=1;i<=N;i++)
    {POS[i]=spos;spos=spos+latlength;
    }
}

void heatcapacity(int hsize,int N)
{
  long double lar=0,max=0,expo2=0,temp=0, C_v=0,num1=0,num2=0,de1=0,de2=0,tem=0,total=0;
  char filename[sizeof "dataheatcapacityPARATEST1.txt"];
  sprintf(filename, "heat64PLE%03d.txt",myid); 
  ofstream myfile; expo=0;
  T=0.000001;
  myfile.open(filename);
  
  while(T<=4)
    {
      
      for(i=lowE;i<=hsize;i++)
	{
	  if(logG[i] > 0)
	    {if(i==0)
		{max=logG[i]-(i/(k_b*T));}
	      else
		{
		  lambda=logG[i]-(i/(k_b*T));
		}
	      if(lambda > max)
		{max=lambda;}
	      else;
	    }
	  else;
	  
	}
      for(i=lowE;i<=hsize;i++)
	{
	  if(logG[i] > 0)
	    {
	      temp=logG[i]-(i/(k_b*T))-max;
	      // tem=logG[i];
	      // expo=powl(e,temp);
	      expo=expl(temp);
	      // expo2=powl(e,temp);
	      num1=num1+ ((i*i))*expo;
	      de1=de1+expo;
	      num2=num2+(i*expo);
	      
	      //de2=de2+(expo);
	      //de2=de1*de1;
	    }
	  else {continue;}
	  
	}
      num2=num2*num2;de2=de1*de1;
      total=((num1/de1) - (num2/de2));
      C_v = (1/(T*T))*total;
      /*if(C_v>lar)
	{lar=C_v;}
	else;*/
      myfile <<T <<'\t'<< C_v/N << '\n';
      // fprintf(ofp,"%.8Lf\t %.8Lf\n",T,C_v);
      total=0; C_v=0;num1=0;num2=0;de1=0;de2=0;tem=0;temp=0;
      T=T+0.01;
    }
  
        myfile.close();
}

void internalenergy(int hsize, int N)
{
  long double expo=0,total=0,S=0,F=0,C_v=0,lar=0,max,pow,U=0,T=0.000001,heanum1=0,dehea1=0, num1=0,num2=0,de2=0,de1=0,tem=0;
  ofstream energyfile,freefile,entropyfile;
  char internal[sizeof "datainternalenergyPARATEST1.txt"];
  char free[sizeof "datafreePARATEST1.txt"];
  char entropy[sizeof "dataentropyPARATEST1.txt"];
  sprintf(internal, "internal64PLE%03d.txt",myid);
  sprintf(free, "free64PLE%03d.txt",myid);
  sprintf(entropy, "entropy64PLE%03d.txt",myid);
  energyfile.open(internal);
  freefile.open(free);
  entropyfile.open(entropy);
  //  heat.open("dataheatcapacity85test14");
  
  while(T<=4)
    {
      for(i=lowE;i<=hsize;i++)
	{
	  if(logG[i] > 0)
	    {if(i==0)
		{max=logG[i]-(i/T);}
	      else
		{
		  lambda=logG[i]-(i/T);
		}
	      if(lambda > max)
		{max=lambda;}
	      else;
	    }
	  else;

	}
      for(i=lowE;i<=hsize;i++)
	{
	  if(logG[i] > 0)
	    { 
	      tem =logG[i]-(i/T)-max;
	      pow=expl(tem);
	      // expo=pow;
	      // expo2=powl(e,temp);
	      heanum1=heanum1+(i*i)*pow;
	      // dehea1=dehea1+pow;
	      num1=num1+(i*pow);
	      de1=de1+pow;
	      de2=de1*de1;
	      
	    }
	  else {continue;}
	}
      //      num2=num1*num1;
      total=((heanum1/dehea1) - (num2/de2));
      // C_v = (1/(T*T))*total;
      U=num1/de1;
      F=-1*T*logl(de1);
      S=(U-F)/T;
      energyfile << T <<'\t'<< U/N<< '\n';
      freefile << T <<'\t'<< F/N<< '\n';
      entropyfile << T <<'\t'<<S/N<< '\n';
      // heat << T << '\t' << C_v/N<< '\n';
      num2=0;de2=0;U=0;F=0;S=0;C_v=0;heanum1=0;dehea1=0;total=0;expo=0;num1=0;num2=0;de2=0;de1=0;tem=0;
      T=T+0.01;
    }
  energyfile.close();
  freefile.close();
  entropyfile.close();
  // heat.close();
}

int histoflat(int hsize, float p)
{int flatcount=0,Htotal=0,Hcount=0;
  long double q,Haverage=0;
  
  for(i=lowE;i<=hsize;i++)
    {if(visited[i]>0)
	{
	  Htotal = Htotal + H[i];
	  Hcount++;
	}
      else;
    }

  Haverage = ((long double)Htotal/(long double)Hcount);

  for(i=lowE;i<=hsize; i++)
    {if(H[i]>0)
	{
	  q=p*Haverage;
	  if((long double)H[i]>=q)
	    {flatcount++;}
	  else;
	}
      else;
    }

  if(Hcount > flatcount)//not all visited energy historgrams fall within the flatness region
    {return 0;}
  else//it is flat so return 1
    {return 1;}
  //else {return 0;}
  
  /* for(i=lowE;i<=hsize;i++)
    {
      if(H[i]>0)
	{flatcount++;}
      else;
    }
   if(flatcount==(hsize-lowE)+1)
    {return 1;}
  else
  {return 0;}*/
    
}

int SAW(int N, int latlength)
{int sw, bc, error=0,upper,lower,right,left;

  for(sw=1;sw<=N;sw++)
    {
      upper=POS[sw]-latlength;
      lower=POS[sw]+latlength;
      right=POS[sw]+1;
      left=POS[sw]-1;
      if(sw==1)
	{
	  if(POS[2]==upper)
	    {;}
	  else if(POS[2]==lower)
	    {;}
	  else if(POS[2]==right)
	    {;}
	  else if(POS[2]==left)
	    {;}
	  else
	    {error=1;}

	}
      else if(sw==N)
	{
	  if(POS[N-1]==upper)
	    {;}
	  else if(POS[N-1]==lower)
	    {;}
	  else if(POS[N-1]==right)
	    {;}
	  else if(POS[N-1]==left)
	    {;}
	  else
	    {error=1;}
	}
      else
	{
	  if(POS[sw-1]==upper)
	    {;}
	  else if(POS[sw-1]==lower)
	    {;}
	  else if(POS[sw-1]==right)
	    {;}
	  else if(POS[sw-1]==left)
	    {;}
	  else
	    {error=1;}

	  if(POS[sw+1]==upper)
	    {;}
	  else if(POS[sw+1]==lower)
	    {;}
	  else if(POS[sw+1]==right)
	    {;}
	  else if(POS[sw+1]==left)
	    {;}
	  else
	    {error=1;}
	}
    }
  if(error>0)
    {return 1;}
  else {return 0;}
  
}

int checkpos(int N,int checkthis)//ensures future positions are actually free
{int ct,counter=0;

    for(ct=1;ct<=N;ct++)
      {if(old[ct]==checkthis)
	  {counter++;}
	else;
      }

    if(counter > 0)
      {return 1;}//not free
    else return 0; //all free


}

void printposarr(int N)
{
  int i;
  for(i=1;i<=N;i++)
    {printf("POS[%d]=%d\n",i,POS[i]);}
  return;
}

void frw(int N, int latlength)
{int ranint,um=0,nosaw=0,go,mono,upper,lower,right,left,above=0,below=0,rc=0,lc=0,uc=0,dc=0,tc=0;
  float ran;
  mono=ranseq(N);
  if(mono==1)
    {above=1;}
  else if(mono==N)
    {below=1;}
  else
    {ran=rndnum();
      if(ran >=0.5)
	{above=1;}
      else
	{below=1;}
    }

  if(above==1 && below==0)
    {
      for(i=mono+1;i<=N;i++)
	{if(um>0)
	    {i=i-1;um=0;}
	  else;
	  upper=POS[i-1]-latlength; lower=POS[i-1]+latlength; right=POS[i-1]+1; left=POS[i-1]-1;

	  ranint=(int)(rndnum()*4)+1;
	  //	  printf("ranint = %d\n",ranint);
	  if(ranint==1)
	    {POS[i]=rightmove(POS[i-1],latlength);nosaw=exclvol2(N,POS[i],i-1,1);rc=1;}
	  else if(ranint==2)
	    {POS[i]=leftmove(POS[i-1],latlength);nosaw=exclvol2(N,POS[i],i-1,1);lc=1;}
	  else if(ranint==3)
	    {POS[i]=lowermove(POS[i-1],latlength);nosaw=exclvol2(N,POS[i],i-1,1);dc=1;}
	  else if(ranint==4)
	    {POS[i]=uppermove(POS[i-1],latlength);nosaw=exclvol2(N,POS[i],i-1,1);uc=1;}
	  tc=rc+lc+dc+uc;
	  if(nosaw==1 && tc < 4)
	    {POS[i]=old[i];nosaw=0;um++;}
	  else if(nosaw==1 && tc==4)//give up?
	    {POS[i]=old[i];return;}
	  else if(nosaw==0)//go to next monomer
	    {rc=lc=uc=dc=tc=nosaw=0;}
	}
      return;
    }
  else if(above==0 && below==1)
    {
      for(i=mono-1;i>=1;i--)
	{if(um>0)
	      {i=i+1;um=0;}
	    else;
	  upper=POS[i+1]-latlength; lower=POS[i+1]+latlength; right=POS[i+1]+1; left=POS[i+1]-1;

	  ranint=(int)(rndnum()*4)+1;
	  //  printf("ranint = %d\n",ranint);
	  if(ranint==1)
	    {POS[i]=rightmove(POS[i+1],latlength);nosaw=exclvol2(N,POS[i],i+1,N);rc=1;}
	  else if(ranint==2)
	    {POS[i]=leftmove(POS[i+1],latlength);nosaw=exclvol2(N,POS[i],i+1,N);lc=1;}
	  else if(ranint==3)
	    {POS[i]=lowermove(POS[i+1],latlength);nosaw=exclvol2(N,POS[i],i+1,N);dc=1;}
	  else if(ranint==4)
	    {POS[i]=uppermove(POS[i+1],latlength);nosaw=exclvol2(N,POS[i],i+1,N);uc=1;}
	  tc=rc+lc+dc+uc;
	  if(nosaw==1 && tc < 4)
	    {POS[i]=old[i];nosaw=0;um++;}
	  else if(nosaw==1 && tc==4)
	    {POS[i]=old[i];return;}
	  else if(nosaw==0)
	    {rc=lc=uc=dc=tc=nosaw=0;}
	}
      return;

    }
  else return;
  
}

void typeii(int N, int latlength)
{int done=0,coun,m=0,coldiff=0,rowdiff=0,bcseq=0,bcconbud=0,conbud=0,seqcount=0,seq=0,toposeq=0,topoconbud=0,smallest=0,largest=0,exch=0,i,j=0;
  float ran;
  for(coun=0;coun<=N;coun++)
  { 
  seq=ranseq(N);//random monomer
  ran=rndnum();

  if(seq==1)
    {conbud=2;}
  else if(seq==N)
    {conbud=N-1;}
  else
    {
      if(ran>0.5)
	{conbud=seq-1;}
      else
	{conbud=seq+1;}
    }

  coldiff=col(POS[seq],latlength)-col(POS[conbud],latlength);
  rowdiff=row(POS[seq],latlength)-row(POS[conbud],latlength);
  coldiff=abs(coldiff);
  rowdiff=abs(rowdiff);

  if(coldiff==0)//in same column so only look for topo neighboours (r  or l)
    {
      for(i=1;i<=N;i++)
	{
	  /* if(i==seq || i==conbud || i==seq+1 || i==seq-1 || i==conbud+1 || i==conbud-1)
	    {continue;}
	    else;*/

	  bcseq=buddycheck(latlength,POS[seq],POS[i]);
	  bcconbud=buddycheck(latlength,POS[conbud],POS[i]);

	  if(bcseq==1)//looking for right neightbours
	    {toposeq=i;}
	  else;
	  if(bcconbud==1)
	    {topoconbud=i;}
	  else;
	  if(toposeq>0 && topoconbud > 0 && toposeq!= topoconbud)//found a square pair!
	    {break;}
	  else;

	  if(bcseq==2)//looking for left neighbours
	    {toposeq=i;}
	  else;
	  if(bcconbud==2)
	    {topoconbud=i;}
	  else;
	  if(toposeq>0 && topoconbud >0 && toposeq != topoconbud)//found a square pair!
	    {break;}
	  else;

	}

      if(toposeq==0 || topoconbud==0)//couldn't form the square
	{topoconbud=0;toposeq=0;continue;}
      else;

      if(topoconbud-toposeq < 0 && conbud-seq < 0)//impose the ascension rule
	{;}
      else if(topoconbud-toposeq > 0 && conbud-seq > 0)
	{;}
      else {topoconbud=0;toposeq=0;continue;}//exit move function (nothing has been change)

      smallest=mini(mini(mini(toposeq,conbud),topoconbud),seq);
      largest=maxi(maxi(maxi(toposeq,conbud),topoconbud),seq);
      POS[smallest]=old[smallest];
      POS[largest]=old[largest];
      POS[1]=old[1];
      POS[N]=old[N];
      for(j=2;j<smallest;j++)//deal with smallest positions
	{if(smallest==1)
 	    {break;}
	  else
	    {POS[j]=old[j];}}
      for(j=largest+1;j<N;j++)//deal with largest positions
	{POS[j]=old[j];}
      m=largest-1;
      for(j=smallest+1;j<largest;j++)
	{POS[j]=old[m];m--;}
      done++;
      if(done>0)
	{return;}
      else;

    }

  else if(rowdiff==0)
    {

      for(i=1;i<=N;i++)
	{
	  /*if(i==seq || i==conbud || i==seq+1 || i==seq-1 || i==conbud+1 || i==conbud-1)
	    {continue;}
	    else;*/

	  bcseq=buddycheck(latlength,POS[seq],POS[i]);
	  bcconbud=buddycheck(latlength,POS[conbud],POS[i]);

	  if(bcseq==4)//looking for upper neightbours
	    {toposeq=4;}
	  else;
	  if(bcconbud==4)
	    {topoconbud=i;}
	  else;
	  if(toposeq>0 && topoconbud > 0 && toposeq!= topoconbud)//found a square pair!
	    {break;}
	  else;

	  if(bcseq==3)//looking for lower neighbours
	    {toposeq=i;}
	  else;
	  if(bcconbud==3)
	    {topoconbud=i;}
	  else;
	  if(toposeq>0 && topoconbud >0 && toposeq != topoconbud)//found a square pair!
	    {break;}
	  else;

	}
      
      if(toposeq==0 || topoconbud==0)//couldn't form the square
	{topoconbud=0;toposeq=0;continue;}
      else;

      if(topoconbud-toposeq < 0 && conbud-seq < 0)//impose the ascension rule
	{;}
      else if(topoconbud-toposeq > 0 && conbud-seq > 0)
	{;}
      else {topoconbud=0;toposeq=0;continue;}//exit move function (nothing has been change)

      smallest=mini(mini(mini(toposeq,conbud),topoconbud),seq);
      largest=maxi(maxi(maxi(toposeq,conbud),topoconbud),seq);
      POS[smallest]=old[smallest];
      POS[largest]=old[largest];
      POS[1]=old[1];
      POS[N]=old[N];
      for(j=2;j<smallest;j++)//deal with smallest positions
	{if(smallest==1)
	    {break;}
	  else
	    {POS[j]=old[j];}}
      for(j=largest+1;j<N;j++)//deal with largest positions
	{POS[j]=old[j];}
      m=largest-1;
      for(j=smallest+1;j<largest;j++)
	{POS[j]=old[m];m--;}
      done++;
      if(done>0)
	{return;}
      else;
      
    }
  else {topoconbud=0;toposeq=0;continue;}

  done=0;topoconbud=0;toposeq=0;
  }
}

int mini(int a, int b)
{
  if(a-b > 0)
    {return b;}
  else
    {return a;}

}

int maxi(int a, int b)
{
  if(a-b > 0)
    {return a;}
  else
    {return b;}

}
void chainterminal(int N, int latlength)
{int topo=0,exch,i,j=0,bead,right,left,upper,lower;
  float ran;


  ran=rndnum();
  if(ran>0.5)
    {bead=1;}
  else
    {bead=N;}
  upper=POS[bead]-latlength;
  lower=POS[bead]+latlength;
  right=POS[bead]+1;
  left=POS[bead]-1;

  if(bead==1)
    {
      for(i=3;i<=N;i++)//scan through positions of monomers looking for topo neigbhbours
	{
	  if(POS[i]==upper)
	    {topo=i;ran=rndnum();
	      if(ran>0.5)
		break;
	      else;
	    }
	  else if(POS[i]==lower)
	    {topo=i;ran=rndnum();
	      if(ran>0.5)
		break;
	      else;
	    }
	  else if(POS[i]==left)
	    {topo=i;ran=rndnum();
	      if(ran>0.5)
		break;
	      else;
	    }
	  else if(POS[i]==right)
	    {topo=i;ran=rndnum();
	      if(ran>0.5)
		break;
	      else;
	    }
	  else;
	  
	}
      if(topo>0)
	{;}
      else {return;}
      
      exch= topo-1;
      j=exch;
      POS[1]=old[exch];
      j=j-1;
      for(i=2;i<=exch;i++)
	{
	  POS[i]=old[j];
	  j=j-1;
	}
    }
  else if(bead==N)
    { for(i=N-2;i>=1;i--)//scan through positions of monomers looking for topo neigbhbours
	{
	  if(POS[i]==upper)
	    {topo=i;ran=rndnum();
	      if(ran>0.5)
		break;
	      else;
	    }
	  else if(POS[i]==lower)
	    {topo=i;ran=rndnum();
	      if(ran>0.5)
		break;
	      else;
	    }
	  else if(POS[i]==left)
	    {topo=i;ran=rndnum();
	      if(ran>0.5)
		break;
	      else;
	    }
	  else if(POS[i]==right)
	    {topo=i;ran=rndnum();
	      if(ran>0.5)
		break;
	      else;
	    }
	  else;
	  
	}
      if(topo>0)
	{;}
      else {return;}
      exch= topo+1;
      j=exch;
      POS[N]=old[exch];
      for(i=N-1;i>=exch;i--)
	{
	  POS[i]=old[j];
	  j=j+1;
	}

    }
  else
    {std::cout << "ERROR bead is not 1 or N" << '\n';}
    
}

int compenergy(int N, int latlength)//computes the total energy of chain and prints(for now)
{int posi,i,j,ene=0,enetot=0,upper,lower,right,left;

  for(i=1;i<=N;i++)//cycle through monomers
    {
      if(HP[i]==1)//only hydrophobic monomers contribute to the energy
	{ upper=POS[i]-latlength;
	  lower=POS[i]+latlength;
	  right=POS[i]+1;
	  left=POS[i]-1;
	  for(j=1;j<=N;j++)//cycle through monomers
	    {if(j==i+1 || j==i-1 || j==i)//monomers obeying this if won't be topo connected
		{continue;}
	      else if(HP[j]==0)//hydrophilic monomers don't contribute to the energy
		{continue;}
	      else
		{posi=POS[j];
		  if(posi==upper && HP[j]==1)//if topo hydrophobic guy above me local contributions= -1
		      {ene=-1;}
		  else if(posi==lower && HP[j]==1)
		      {ene=-1;}
		  else if(posi==right && HP[j]==1)
		      {ene=-1;}
		  else if(posi==left && HP[j]==1)
		      {ene=-1;}
		    else
		      {continue;}
		  enetot=enetot+ene;ene=0;//add the contribution to the total energy
		}
	    }
	}
      else continue;
    }	

  return enetot/2;//total energy =enetot/2 (due to double counting)

}

void kinkflip(int N,int latlength,int size)//performs a kinkflip move
{
  int check,seqmono, upperguy, lowerguy,rowdiff,coldiff,bcu,bcl;
  // cleanlattice(size,N);
  for(seqmono=2;seqmono<N;seqmono++)//goes through chain looking for kinks
    {
      upperguy=seqmono+1;//bead # ahead of me (on the sequence)
      lowerguy=seqmono-1;//bead # behind me " "
      rowdiff=abs(row(POS[upperguy],latlength)-row(POS[lowerguy],latlength));//get absolute value of row difference
      coldiff=abs(col(POS[upperguy],latlength)-col(POS[lowerguy],latlength));//get absolute value of col difference
      
      if(rowdiff==1 && coldiff==1)
	{bcu=buddycheck(latlength,POS[seqmono],POS[upperguy]);//check orientation of seqmono with ahead bead
	  bcl=buddycheck(latlength,POS[seqmono],POS[lowerguy]);//check oritentation of seqmono with behind bead
	  if(abs(bcu-bcl)==3 || abs(bcu-bcl)==1)//technically a move top right quadrant or bottom left quadrant
	    {if(bcu==4 ||bcl==4)//would be a top right quadrant move
		{check=checkpos(N,POS[seqmono]-latlength+1);
		  if(check==0)//if the position is avaliable
		    {
		      POS[seqmono]=POS[seqmono]-(latlength+1);return;
		    }//job done monomer moved
		  else continue;
		}
	      else if(bcu==3 || bcl==3)//bottom left quadrant move
		{check=checkpos(N,POS[seqmono]+(latlength-1));
		  if(check==0)//is the position available?
		    {
		      POS[seqmono]=POS[seqmono]+(latlength-1);return;
		    }//job done monomer moved
		  else continue;//cant do anything so give up
		}
	      else continue;
	    }
	  
	  else if(abs(bcu-bcl)==2)//technically a move top left quadrant or bottom right quadrant
	    {if(bcu==4 || bcl==4)//would be a top left quadrant move
		{check=checkpos(N,POS[seqmono]-(latlength-1));
		  if(check==0)//if the space is free
		    {
		      POS[seqmono]=POS[seqmono]-(latlength-1);return;
		    }//job done
		  else continue;//cant do anything
		}
	      else if(bcu==3 || bcl==3)//bottom right quadrant move
		{check=checkpos(N,POS[seqmono]+(latlength+1));
		  if(check==0)//is the position free?
		    {
		      POS[seqmono]=POS[seqmono]+(latlength+1);return;
		    }
		  else continue;
		}
	    }
	  else continue;
	}
      else continue;
    }
  return;
  
}

void pivot(int N,int axismono,int latlength, int size)
{//cleanlattice(size,N);
  int cnt=0,g,arrow,check;
  int *direct;
  direct =(int *) malloc((N+1) * sizeof(int));
  if(direct==NULL)
    {printf("OUT OF MEMORY FOR PIVOT MOVE\n");exit(1);}
  
  else//memory is free
    {
      if(axismono==1 || axismono==N)
	{free(direct);return;}//rotating around the ends is pointless
      
      else//monomer is between 1 and N and hence a pivot move is significant
	{
	  if(rndnum()<=0.5)//focus on monomers with bead # < axismono
	    {
	      for(g=axismono-1;g>=1;g--)
		{
		  direct[g]=buddycheck(latlength,POS[g+1],POS[g]);
		  // printf("Direction value of %d wrt %d = %d\n",g,g+1,direct[g]);
		}
	      //**********************************************************************
	      if(rndnum()<=0.5)//perform an anticlockwise move
		{
		  for(g=axismono-1;g>=1;g--)
		    {arrow=direct[g];
		      switch (arrow)
			{
			case 1://originally was right to g-1
			  {direct[g]=4;break;}
			case 2:
			  {direct[g]=3;break;}
			case 3:
			  {direct[g]=1;break;}
			case 4:
			  {direct[g]=2;break;}
			default:
			  {continue;}
			}
		    }
		  for(g=axismono-1;g>=1;g--)//check future positions are actually free
		    {arrow=direct[g];
		      switch (arrow)
			{
			case 1:
			  {check=checkpos(N,POS[g+1]+1);
			    if(check==0)
			      {break;}
			    else
			      {cnt++;break;}
			  }
			case 2:
			  {check=checkpos(N,POS[g+1]-1);
			    if(check==0)
			      {break;}
			    else
				{cnt++;break;}
			  }
			case 3:
			  {check=checkpos(N,POS[g+1]+latlength);
			    if(check==0)
			      {break;}
			    else
			      {cnt++;break;}
			  }
			case 4:
			  {check=checkpos(N,POS[g+1]-latlength);
			    if(check==0)
			      {break;}
			    else
			      {cnt++;break;}
			  }
			default:
			  {cnt++;break;}
			}
		    }
		  if(cnt==0)//positions are free so can actually move them
		    {for(g=axismono-1;g>=1;g--)
			{arrow=direct[g];
			  switch (arrow)
			    {
			    case 1:
			      {POS[g]=POS[g+1]+1;break;}
			    case 2:
			      {POS[g]=POS[g+1]-1;break;}
			    case 3:
			      {POS[g]=POS[g+1]+latlength;break;}
			    case 4:
			      {POS[g]=POS[g+1]-latlength;break;}
			    default:
			      {break;}
			    }
			}
		    }
		  else {;}//do nothing since no positions are free
		}
	      //************************************************************************
	      else//perform a clockwise move
		{
		  for(g=axismono-1;g>=1;g--)//change directions due to clockwise move
		    {arrow=direct[g];
		      switch (arrow)
			{
			case 1:
			  {direct[g]=3;break;}
			case 2:
			  {direct[g]=4;break;}
			case 3:
			  {direct[g]=2;break;}
			case 4:
			  {direct[g]=1;break;}
			default:
			  {continue;}
			}
		    }
		  for(g=axismono-1;g>=1;g--)//check positions are free
		    {arrow=direct[g];
		      switch (arrow)
			{
			case 1:
			  {check=checkpos(N,POS[g+1]+1);
			    if(check==0)
			      {break;}
			    else
			      {cnt++;break;}
			  }
			case 2:
			  {check=checkpos(N,POS[g+1]-1);
			    if(check==0)
			      {break;}
			    else
			      {cnt++;break;}
			  }
			case 3:
			  {check=checkpos(N,POS[g+1]+latlength);
			    if(check==0)
			      {break;}
			    else
			      {cnt++;break;}
			  }
			case 4:
			  {check=checkpos(N,POS[g+1]-latlength);
			    if(check==0)
			      {break;}
			    else
			      {cnt++;break;}
			  }
			default:
			  {cnt++;break;}
			}
		    }
		  if(cnt==0)
		    {for(g=axismono-1;g>=1;g--)
			{arrow=direct[g];
			  switch (arrow)
			    {
			    case 1:
			      {POS[g]=POS[g+1]+1;break;}
			    case 2:
			      {POS[g]=POS[g+1]-1;break;}
			    case 3:
			      {POS[g]=POS[g+1]+latlength;break;}
			    case 4:
			      {POS[g]=POS[g+1]-latlength;break;}
			    default:
			      {break;}
			    }
			}
		    }
		  else {;}//do nothing since no positions are free
		}
	    }



       //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	  else//rndnum()>0.5 => focus on monomers with bead # > axismono
	    {
	      for(g=axismono+1;g<=N;g++)
		{
		  direct[g]=buddycheck(latlength,POS[g-1],POS[g]);
		  // printf("Direction value of %d wrt %d = %d\n",g,g-1,direct[g]); 
		}
	      if(rndnum()<=0.5)//perform an anticlockwise move
		{
		  for(g=axismono+1;g<=N;g++)
		    {arrow=direct[g];
		      switch (arrow)
			{
			case 1://originally was right
			  {direct[g]=4;break;}
			case 2://originally was left
			  {direct[g]=3;break;}
			case 3://originally was below
			  {direct[g]=1;break;}
			case 4://originally was above
			  {direct[g]=2;break;}
			default:
			  {continue;}
			}
		    }
		  for(g=axismono+1;g<=N;g++)//check future positions are actually free
		    {arrow=direct[g];
		      switch (arrow)
			{
			case 1://would be to the right of g-1 guy
			  {check=checkpos(N,POS[g-1]+1);
			    if(check==0)//if it's free do nothing
			      {break;}
			    else
			      {cnt++;break;}//if it's not free add it to the cnt
			  }
			case 2://would be to the left of g-1 guy
			  {check=checkpos(N,POS[g-1]-1);
			    if(check==0)
			      {break;}
			    else
			      {cnt++;break;}
			  }
			case 3://would be below the g-1 guy
			  {check=checkpos(N,POS[g-1]+latlength);
			    if(check==0)
			      {break;}
			    else
			      {cnt++;break;}
			  }
			case 4://would be above the g-1 guy
			  {check=checkpos(N,POS[g-1]-latlength);
			    if(check==0)
			      {break;}
			    else
			      {cnt++;break;}
			  }
			default:
			  {cnt++;}
			}
		    }
		  if(cnt==0)//positions are free
		    {for(g=axismono+1;g<=N;g++)//implement new positions
			{arrow=direct[g];
			  switch (arrow)
			    {
			    case 1://move it to the right of g-1
			      {POS[g]=POS[g-1]+1;break;}
			    case 2://move it to the left of g-1
			      {POS[g]=POS[g-1]-1;break;}
			    case 3://move it below g-1
			      {POS[g]=POS[g-1]+latlength;break;}
			    case 4://move it above g-1
			      {POS[g]=POS[g-1]-latlength;break;}
			    default:
			      {break;}
			    }
			}
		    }
		  else {;}//do nothing since not all the positions are free
		}
	      
	      else//perform a clockwise move
		{
		  for(g=axismono+1;g<=N;g++)
		    {arrow=direct[g];
		      switch (arrow)
			{
			case 1:
			  {direct[g]=3;break;}
			case 2:
			  {direct[g]=4;break;}
			case 3:
			  {direct[g]=2;break;}
			case 4:
			  {direct[g]=1;break;}
			default:
			  {continue;}
			}
		    }
		  for(g=axismono+1;g<=N;g++)//check positions are acutally free
		    {arrow=direct[g];
		      switch (arrow)
			{
			case 1: 
			  {check=checkpos(N,POS[g-1]+1);
			    if(check==0)
			      {break;}
			    else
			      {cnt++;break;}
			  }
			case 2:
			  {check=checkpos(N,POS[g-1]-1);
			    if(check==0)
			      {break;}
			    else
			      {cnt++;break;}
			  }
			case 3:
			  {check=checkpos(N,POS[g-1]+latlength);
			    if(check==0)
			      {break;}
			    else
			      {cnt++;break;}
			  }
			case 4:
			  {check=checkpos(N,POS[g-1]-latlength);
			    if(check==0)
			      {break;}
			    else
			      {cnt++;break;}
			  }
			default:
			  {cnt++;break;}
			}
		    }
		  if(cnt==0)//positions are free so actually move them
		    {for(g=axismono+1;g<=N;g++)
			{arrow=direct[g];
			  switch (arrow)
			    {
			    case 1:
			      {POS[g]=POS[g-1]+1;break;}
			    case 2:
			      {POS[g]=POS[g-1]-1;break;}
			    case 3:
			      {POS[g]=POS[g-1]+latlength;break;}
			    case 4:
			      {POS[g]=POS[g-1]-latlength;break;}
			    default:
			      {break;}
			    }
			}
		    }
		  else {;}//do nothing since no positions are free
		}
	    }
	}
    }
  
  free(direct);
  

}

int buddycheck(int latlength,int mono1, int mono2)/*(mono1 and mono2 are position values checks connection type)*/
{int  d,rl,ud;

  d=mono1-mono2;
  rl=latlength-1;ud=latlength*rl;

  if(d==-1)
      {return 1;}
  else if(d==1)
      {return 2;}
  else if(d==-latlength)
      {return 3;}
  else if(d==latlength)
      {return 4;}
  else if(d==-ud)
      {return 4;}
  else if(d==ud)
      {return 3;}
  else if(d==rl)
      {return 1;}
  else if (d==-rl)
      {return 2;}
  else
      {return 5;}


}

void pullmove(int N, int seqmono, int latlength,int size)
{//cleanlattice(size,N);
  int sw,check1,check2,g,bc,choice1=0,choice2=0,proposed,seqmono2;
    float k;
    /* printf("seqmono = %d\n", seqmono);*/
    oldposition(N);//stores old positions in an array

    //========== Chooses the 'axis' of the move which means this bead is fixed and everything moves around it in a pull fashion === 
    if(seqmono==1)// if I am 1 my only sequence neighbour is bead 2
      {seqmono2=2;/*printf("seqmono2 = %d\n",seqmono2);*/}

    else if (seqmono==N)// if I am N my only sequence neighbour is bead N-1
    {seqmono2=N-1;/*printf("seqmono2 = %d\n",seqmono2);*/}

    else// else pick with 50/50 chance a bead to be the fixed bead as I have 2 neighbours 
    {k=rndnum();
      if(k<=0.5)
	{seqmono2=seqmono-1;/*printf("seqmono2 = %d\n", seqmono2);*/}
      else if(k>0.5)
	  {seqmono2=seqmono+1;/*printf("seqmono2=%d\n",seqmono2);*/}
    }
  //===============================================================================
  bc=buddycheck(latlength, POS[seqmono], POS[seqmono2]);//check relative neighbor position
  /* printf("buddycheck gives %d\n",bc);*/
  switch (bc)//depending on neighbour position and orientation dictates how the pull move will be carried out
    {
      //&&&&&&&&&&&&&&&&  HORIZANTAL CASE WITH AXIS BEAD ON THE RIGHT &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    case 1://horizontal case where mono2 is on the right
      {check1=checkpos(N,POS[seqmono2]-latlength);check2=checkpos(N,POS[seqmono2]-latlength-1); 
	if(check1==0 && check2==0)//if the position above me is free do..
	  {choice1++;}
	else; 
	check1=checkpos(N,POS[seqmono2]+latlength);check2=checkpos(N,POS[seqmono2]+latlength-1);
	if(check1==0 && check2==0)//if the position below me is free do..
	  {choice2++;}
	else;

	if(choice1>0 && choice2>0)//if pos above and below is free make the choice to go either way random
	  {
	    if(rndnum()<=0.5)//choice1 moveup
	      {if (seqmono==1 || seqmono==N)
		  {
		    POS[seqmono]=uppermove(POS[seqmono2],latlength);break;//just need to move it above easy peasy
		  } 
		else if(seqmono2>seqmono)//if new axis bead # > current bead #
		  {
		POS[seqmono]=uppermove(POS[seqmono2],latlength);
		POS[seqmono-1]=POS[seqmono]-1; //move penultimate bead to the left of seqmono
		for(g=1;g<seqmono-1;g++)//drag chain along
		  {POS[g]=old[g+2];
		    sw=SAW(N,latlength);
		    if(sw==0)
		      {break;}
		    else
		      {continue;}
		  }
		break;
		  }
		else//if new axis bead # < current bead #
		  {
		    POS[seqmono]=uppermove(POS[seqmono2],latlength);
		    POS[seqmono+1]=POS[seqmono]-1;//move my other neighbor to the correct spot
		    for(g=(seqmono+2);g<=N;g++)//drag chain along
		      {POS[g]=old[g-2];
			sw=SAW(N,latlength);
			if(sw==0)
			  {break;}
			else
			  {continue;}
		      }
		    break;
		  }
	      }

	    else//choice2 movedown
	      {if(seqmono==1 || seqmono==N)
		  {
		    POS[seqmono]=lowermove(POS[seqmono2],latlength);break;//just need to move it down easy peasy
		  }
		else if(seqmono2>seqmono)//if new axis bead # > current bead #
		  {
		    POS[seqmono]=lowermove(POS[seqmono2],latlength);
		    POS[seqmono-1]=POS[seqmono]-1;//move penultimate bead to the left of seqmono
		    for(g=1;g<seqmono-1;g++)//drag chain along
		      {POS[g]=old[g+2];
			sw=SAW(N,latlength);
			if(sw==0)
			  {break;}
			else
			  {continue;}
		      } break;//e.g 1 goes to where 3 was
		  }
		else//if new axis bead # < current bead #
		  {
		    POS[seqmono]=lowermove(POS[seqmono2],latlength);
		    POS[seqmono+1]=POS[seqmono]-1;//move my other neighbor to the correct spot
		    for(g=(seqmono+2);g<=N;g++)
		      {POS[g]=old[g-2];
			sw=SAW(N,latlength);
			if(sw==0)
			  {break;}
			else
			  {continue;}
		      } break;// e.g N goes to where N-2 was
		  }
	      }
	  }

	else if(choice1>0 && choice2==0)//choice 1 moveup unconditionally since it is the only one free
	  {if(seqmono==1 || seqmono==N)
	      {
		POS[seqmono]=uppermove(POS[seqmono2],latlength);break;//just need to move it up easy peasy
	      }
	    else if(seqmono2>seqmono)//if new axis bead # > current bead #
	      {
		POS[seqmono]=uppermove(POS[seqmono2],latlength);
		POS[seqmono-1]=POS[seqmono]-1; //move penultimate bead to the left of seqmono
		for(g=1;g<seqmono-1;g++)//drag chain along
		  {POS[g]=old[g+2];
		    sw=SAW(N,latlength);
		    if(sw==0)
		      {break;}
		    else
		      {continue;}
		  }break;
	      }
	    else//if new axis bead # < current bead #
	      {
		POS[seqmono]=uppermove(POS[seqmono2],latlength);
		POS[seqmono+1]=POS[seqmono]-1;//move my other neightbor to the correct spot
		for(g=(seqmono+2);g<=N;g++)//drag chain along
		  {POS[g]=old[g-2];
		    sw=SAW(N,latlength);
		    if(sw==0)
		      {break;}
		    else
		      {continue;}
		  }break;
	      }
	  }

	else if(choice1==0 && choice2>0)//choice 2 movedown unconditionally since it is the only one free
	  {if(seqmono==1 || seqmono==N)
	      {
		POS[seqmono]=lowermove(POS[seqmono2],latlength);break;//just need to move it down easy peasy
	      }
	    else if(seqmono2>seqmono)//if new axis bead # > current bea #
	      {
		POS[seqmono]=lowermove(POS[seqmono2],latlength);
		POS[seqmono-1]=POS[seqmono]-1; //move penultimate bead to the left of seqmono
		for(g=1;g<seqmono-1;g++)//drag chain along
		  {POS[g]=old[g+2];
		    sw=SAW(N,latlength);
		    if(sw==0)
		      {break;}
		    else
		      {continue;}
		  }break;
	      }
	    else//if new axis bead # < current bead #
	      {
		POS[seqmono]=lowermove(POS[seqmono2],latlength);
		POS[seqmono+1]=POS[seqmono]-1;//move my other neighbor to the correct spot
		for(g=(seqmono+2);g<=N;g++)//drag chain along
		  {POS[g]=old[g-2];
		    sw=SAW(N,latlength);
		    if(sw==0)
		      {break;}
		    else
		      {continue;}
		  }break;

	      }
	  }

	else break; //do nothing since below and above are not available

      }

      //&&&&&&&&&&&&&&&&&&&&&&&&&& HORIZONTAL CASE WHERE AXIS BEAD IS ON THE LEFT &&&&&&&&&&&&&&&
    case 2: //horizontal partner to my left
      {check1=checkpos(N,POS[seqmono2]-latlength);check2=checkpos(N,POS[seqmono2]-latlength+1);
	if(check1==0 && check2==0)//to check position above and to the right
	  {choice1++;}
	else;
	check1=checkpos(N,POS[seqmono2]+latlength);check2=checkpos(N,POS[seqmono2]+latlength+1);
	if(check1==0 && check2==0)//to check position below and to the right
	  {choice2++;}
	else;

	if(choice1>0 && choice2>0)
	  {if(rndnum()<=0.5)//randomly move up
	      {if(seqmono==1 || seqmono==N)
		  {
		    POS[seqmono]=uppermove(POS[seqmono2],latlength);break;//just move it above no messing with neighbours
		  }
		else if(seqmono2>seqmono)//if new axis bead # > current bead #
		  {
		    POS[seqmono]=uppermove(POS[seqmono2],latlength);
		    POS[seqmono-1]=POS[seqmono]+1;
		    for(g=1;g<seqmono-1;g++)
		      {POS[g]=old[g+2];
			sw=SAW(N,latlength);
			if(sw==0)
			  {break;}
			else
			  {continue;}
		      }break;
		  }
		else//if new axis bead # < current bead #
		  {
		    POS[seqmono]=uppermove(POS[seqmono2],latlength);
		    POS[seqmono+1]=POS[seqmono]+1;
		    for(g=(seqmono+2);g<=N;g++)
		      {POS[g]=old[g-2];
			sw=SAW(N,latlength);
			if(sw==0)
			  {break;}
			else
			  {continue;}
		      }break;
		  }
	      }

	    else//random number >0.5 then move down
	      {if(seqmono==1 || seqmono==N)
		  {
		    POS[seqmono]=lowermove(POS[seqmono2],latlength);break;//just move it down no messing with neighbours
		  }
		else if(seqmono2>seqmono)//if new axis bead # > current bead #
		  {
		    POS[seqmono]=lowermove(POS[seqmono2],latlength);
		    POS[seqmono-1]=POS[seqmono]+1;
		    for(g=1;g<seqmono-1;g++)
		      {POS[g]=old[g+2];
			sw=SAW(N,latlength);
			if(sw==0)
			  {break;}
			else
			  {continue;}
		      }break;
		  }
		else//if new axis bead # < current bead #
		  {
		    POS[seqmono]=lowermove(POS[seqmono2],latlength);
		    POS[seqmono+1]=POS[seqmono]+1;
		    for(g=(seqmono+2);g<=N;g++)
		      {POS[g]=old[g-2];
			sw=SAW(N,latlength);
			if(sw==0)
			  {break;}
			else
			  {continue;}
		      }break;
		  }
	      }
	  }

	else if(choice1>0 && choice2==0)// move up unconditionally
	  {if(seqmono==1 || seqmono==N)
	      {
		POS[seqmono]=uppermove(POS[seqmono2],latlength);break;//just move it above no sweat
	      }
	    else if(seqmono2>seqmono)//if new axis bead # > current bead #
	      {
		POS[seqmono]=uppermove(POS[seqmono2],latlength);
		POS[seqmono-1]=POS[seqmono]+1;
		for(g=1;g<seqmono-1;g++)
		  {POS[g]=old[g+2];
		    sw=SAW(N,latlength);
		    if(sw==0)
		      {break;}
		    else
		      {continue;}
		  }break;
	      }
	    else//if new axis bead # < current bead #
	      {
		POS[seqmono]=uppermove(POS[seqmono2],latlength);
		POS[seqmono+1]=POS[seqmono]+1;
		for(g=(seqmono+2);g<=N;g++)
		  {POS[g]=old[g-2];
		    sw=SAW(N,latlength);
		    if(sw==0)
		      {break;}
		    else
		      {continue;}
		  }break;
	      }
	  }

	else if(choice1==0 && choice2>0)//move down unconditionally
	  {if(seqmono==1 || seqmono==N)
	      {
		POS[seqmono]=lowermove(POS[seqmono2],latlength);break;//just move it down no messing with neighbours
	      }
	    else if(seqmono2>seqmono)//if new axis bead # > current bead #
	      {
		POS[seqmono]=lowermove(POS[seqmono2],latlength);
		POS[seqmono-1]=POS[seqmono]+1;
		for(g=1;g<seqmono-1;g++)
		  {POS[g]=old[g+2];
		    sw=SAW(N,latlength);
		    if(sw==0)
		      {break;}
		    else
		      {continue;}
		  }break;
	      }
	    else//if new axis bead # < current bead #
	      {
		POS[seqmono]=lowermove(POS[seqmono2],latlength);
		POS[seqmono+1]=POS[seqmono]+1;
		for(g=(seqmono+2);g<=N;g++)
		  {POS[g]=old[g-2];
		    sw=SAW(N,latlength);
		    if(sw==0)
		      {break;}
		    else
		      {continue;}
		  }break;
	      }
	  }

	else break;//choice1==0 and choice2==0 which means SAW conditions will be violated

      }
   
 //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& VERTICAL CASE WITH AXIS BEAD BELOW &&&&&&&&&&&&&&&&&&&&&&

    case 3://vertical partner below me
      {check1=checkpos(N,POS[seqmono2]+1);check2=checkpos(N,POS[seqmono2]+1-latlength);
	if(check1==0 && check2==0)
	  {choice1++;}//check position to right and above is free for neighbour
	else;
	check1=checkpos(N,POS[seqmono2]-1); check2=checkpos(N,POS[seqmono2]-1-latlength);
	if(check1==0 && check2==0)
	  {choice2++;}
	else;

	if(choice1>0 && choice2>0)//it is valid for me to move to the left or right so choose randomly
	  {if(rndnum()<=0.5)//randomly move to the right
	      {if(seqmono==1 || seqmono==N)
		  {
		    POS[seqmono]=rightmove(POS[seqmono2],latlength);break;//just move to the right no problem
		  }
		else if(seqmono2>seqmono)//if new axis bead # > current bead #
		  {
		    POS[seqmono]=rightmove(POS[seqmono2],latlength);
		    POS[seqmono-1]=POS[seqmono]-latlength;
		    for(g=1;g<seqmono-1;g++)
		      {POS[g]=old[g+2];
			sw=SAW(N,latlength);
			if(sw==0)
			  {break;}
			else
			  {continue;}
		      }break;
		  }
		else//if new axis bead # < current bead number
		  {
		    POS[seqmono]=rightmove(POS[seqmono2],latlength);
		    POS[seqmono+1]=POS[seqmono]-latlength;
		    for(g=seqmono+2;g<=N;g++)
		      {POS[g]=old[g-2];
			sw=SAW(N,latlength);
			if(sw==0)
			  {break;}
			else
			  {continue;}
		      }break;
		  }
	      }
	   
	    else//random number > 0.5 then move left
	      {if(seqmono==1 || seqmono==N)
		  {
		    POS[seqmono]=leftmove(POS[seqmono2],latlength);break;
		  }

		else if(seqmono2>seqmono)//if new axis bead # > current bead #
		  {
		    POS[seqmono]=leftmove(POS[seqmono2],latlength);
		    POS[seqmono-1]=POS[seqmono]-latlength;
		    for(g=1;g<(seqmono-1);g++)
		      {POS[g]=old[g-2];
			sw=SAW(N,latlength);
			if(sw==0)
			  {break;}
			else
			  {continue;}
		      }break;
		  }
		else//if new axis bead # < current bead number
		  {
		    POS[seqmono]=leftmove(POS[seqmono2],latlength);
		    POS[seqmono+1]=POS[seqmono]-latlength;
		    for(g=seqmono+2;g<=N;g++)
		      {POS[g]=old[g-2];
			sw=SAW(N,latlength);
			if(sw==0)
			  {break;}
			else
			  {continue;}
		      }break;
		  }
	      }
	  }

	else if(choice1>0 && choice2==0)//make a rightmove unconditionally 
	  {if(seqmono==1 || seqmono==N)
	      {
		POS[seqmono]=rightmove(POS[seqmono2],latlength);break;
	      }
	    else if(seqmono2>seqmono)//if new axis bead > current bead #
	      {
		POS[seqmono]=rightmove(POS[seqmono2],latlength);
		POS[seqmono-1]=POS[seqmono]-latlength;
		for(g=1;g<seqmono-1;g++)
		  {POS[g]=old[g+2];
		    sw=SAW(N,latlength);
		    if(sw==0)
		      {break;}
		    else
		      {continue;}
		  }break;
	      }
	    else//if new axis bead # < current bead # 
	      {
		POS[seqmono]=rightmove(POS[seqmono2],latlength);
		POS[seqmono+1]=POS[seqmono]-latlength;
		for(g=seqmono+2;g<=N;g++)
		  {POS[g]=old[g-2];
		    sw=SAW(N,latlength);
		    if(sw==0)
		      {break;}
		    else
		      {continue;}
		  }break;
	      }
	  }

	else if(choice1==0 && choice2>0)//can only move to the left
	  {if(seqmono==1 || seqmono==N)
	      {
		POS[seqmono]=leftmove(POS[seqmono2],latlength);break;// just move to left no sweat
	      }
	    else if(seqmono2>seqmono)// if new axis bead #> current bead # (seqmono)
	      {
		POS[seqmono]=leftmove(POS[seqmono2],latlength);
		POS[seqmono-1]=POS[seqmono]-latlength;
		for(g=1;g<seqmono-1;g++)
		  {POS[g]=old[g+2];
		    sw=SAW(N,latlength);
		    if(sw==0)
		      {break;}
		    else
		      {continue;}
		  }break;
	      }
	    else//if new axis bead # < current bead #
	      {
		POS[seqmono]=leftmove(POS[seqmono2],latlength);
		POS[seqmono+1]=POS[seqmono]-latlength;
		for(g=seqmono+2;g<=N;g++)
		  {POS[g]=old[g-2];
		    sw=SAW(N,latlength);
		    if(sw==0)
		      {break;}
		    else
		      {continue;}
		  }break;
	      }
	  }
	else break;// do nothing since no position is legally available
	  

      }//closes case 3

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& VERTICAL CASE WHERE SEQMONO2 IS ABOVE ME &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    case 4:
      {check1=checkpos(N,POS[seqmono2]+1);check2=checkpos(N,POS[seqmono2]+1+latlength);
	if(check1==0 && check2==0)
	{choice1++;}//check position to right and below is free
	else;
	check1=checkpos(N,POS[seqmono2]-1);check2=checkpos(N,POS[seqmono2]-1+latlength);
	if(check1==0 && check2==0)
	  {choice2++;}//check position to left and below is free
	else;

	if(choice1>0 && choice2>0)//choose randomly between right and left
	  {if(rndnum()<=0.5)//move right
	      {if(seqmono==1 || seqmono==N)
		  {
		    POS[seqmono]=rightmove(POS[seqmono2],latlength);break;
		  }//move right no sweating about the neighbours
		else if(seqmono2>seqmono)
		  {
		    POS[seqmono]=rightmove(POS[seqmono2],latlength);
		    POS[seqmono-1]=POS[seqmono]+latlength;
		    for(g=1;g<seqmono-1;g++)
		      {POS[g]=old[g+2];
			sw=SAW(N,latlength);
			if(sw==0)
			  {break;}
			else
			  {continue;}
		      }break;
		  }
		else//if semono2 < seqmono
		  {
		    POS[seqmono]=rightmove(POS[seqmono2],latlength);
		    POS[seqmono+1]=POS[seqmono]+latlength;
		    for(g=seqmono-2;g<=N;g++)
		      {POS[g]=old[g-2];
			sw=SAW(N,latlength);
			if(sw==0)
			  {break;}
			else
			  {continue;}
		      }break;
		  }
	      }

	    else//random number > 0.5 then move to the left
	      {if(seqmono==1 || seqmono==N)
		  {
		    POS[seqmono]=leftmove(POS[seqmono2],latlength);break;
		  }
		else if(seqmono2>seqmono)
		  {
		    POS[seqmono]=leftmove(POS[seqmono2],latlength);
		    POS[seqmono-1]=POS[seqmono]+latlength;
		    for(g=1;g<seqmono-1;g++)
		      {POS[g]=old[g+2];
			sw=SAW(N,latlength);
			if(sw==0)
			  {break;}
			else
			  {continue;}
		      }break;
		  }
		else
		  {
		    POS[seqmono]=leftmove(POS[seqmono2],latlength);
		    POS[seqmono+1]=POS[seqmono]+latlength;
		    for(g=seqmono+2;g<=N;g++)
		      {POS[g]=old[g-2];
			sw=SAW(N,latlength);
			if(sw==0)
			  {break;}
			else
			  {continue;}
		      }break;
		  }
	      }
	  }

	else if(choice1>0 && choice2==0)//make a right turn unconditionally
	  {if(seqmono==1 || seqmono==N)
	      {
		POS[seqmono]=rightmove(POS[seqmono2],latlength);break;
	      }

	    else if(seqmono2> seqmono)
	      {
		POS[seqmono]=rightmove(POS[seqmono2],latlength);
		POS[seqmono-1]=POS[seqmono]+latlength;
		for(g=1;g<seqmono-1;g++)
		  {POS[g]=old[g+2];
		    sw=SAW(N,latlength);
		    if(sw==0)
		      {break;}
		    else
		      {continue;}
		  }break;
	      }
	    else//seqmono2 < seqmono
	      {
		POS[seqmono]=rightmove(POS[seqmono2],latlength);
		POS[seqmono+1]=POS[seqmono]+latlength;
		for(g=seqmono+2;g<=N;g++)
		  {POS[g]=old[g-2];
		    sw=SAW(N,latlength);
		    if(sw==0)
		      {break;}
		    else
		      {continue;}
		  }break;
	      }
	  }

	else if(choice1==0 && choice2>0)//move left unconditionally
	  {if(seqmono==1 || seqmono==N)
	      {
		POS[seqmono]=leftmove(POS[seqmono2],latlength);break;
	      }
	    else if(seqmono2>seqmono)
	      {
		POS[seqmono]=leftmove(POS[seqmono2],latlength);
		POS[seqmono-1]=POS[seqmono]+latlength;
		for(g=1;g<seqmono-1;g++)
		  {POS[g]=old[g+2];
		    sw=SAW(N,latlength);
		    if(sw==0)
		      {break;}
		    else
		      {continue;}
		  }break;
	      }
	    else//seqmono2 < seqmono
	      {
		POS[seqmono]=leftmove(POS[seqmono2],latlength);
		POS[seqmono+1]=POS[seqmono]+latlength;
		for(g=seqmono+2;g<=N;g++)
		  {POS[g]=old[g-2];
		    sw=SAW(N,latlength);
		    if(sw==0)
		      {break;}
		    else
		      {continue;}
		  }break;
	      }
	  }

	else break;
      }
	  
    }//closes switch
}//function ends
		  

int ranseq(int N)//produces a random sequence number to pick a monomer for moves etc.
{
  int k;

  k= (int)(rndnum()*N)+1;
  return k;
}

int uppermove(int posmono2, int latlength)//returns new upper position if you want to move above + checks pbc
{int newpos;

   if(row(posmono2,latlength)==0)
    {
  newpos=posmono2+(latlength*(latlength-1));
    }
  else newpos=posmono2-latlength;

  return newpos;
  //newpos=posmono2-latlength;
  return newpos;
}

int lowermove(int posmono2, int latlength)//returns new lower position if you want to move below + checks pbc
{int newpos;

  if(row(posmono2,latlength)==(latlength-1))
    {
      newpos=posmono2-(latlength*(latlength-1));
    }
  else newpos=posmono2+latlength;
  
  //newpos=posmono2+latlength;
  return newpos;
}

int leftmove(int posmono2, int latlength)
{int newpos;
  
  if(col(posmono2,latlength)==0)
    {newpos=posmono2+(latlength-1);
    }
  else newpos=posmono2-1;
  
  //newpos=posmono2-1;
  return newpos;
}

int rightmove(int posmono2, int latlength)
{int newpos;
  
  if(col(posmono2,latlength)==(latlength-1))
    {
      newpos=posmono2-(latlength-1);
    }
  else newpos=posmono2+1;
  
  // newpos=posmono2+1;
  return newpos;
}

void oldposition(int N)
{
  for(i=1;i<=N;i++)
    {
      old[i]=POS[i];
    }
}

int row(int pos, int latlength) 
{
  return (int)(pos/latlength);

}

int col(int pos, int latlength)
{
  return pos%latlength;

}

