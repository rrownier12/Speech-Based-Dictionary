#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <iostream>
#include <fstream>
#include <conio.h>
#include <string>
#include "util.h"
#include "test.h"
int choice=0;

#define N 5
#define M 32
#define T 60 //number of time frame
long double threshold =1e-30;
using namespace std;
long double delta[T][N];
long double p_star;
long double a_matrix[6][6];   //array to store the a matrix values
long double b_matrix[6][33];  //store the values of the b matrix
int pi[6];

long double avg_a[6][6];
long double avg_b[6][33];
int avg_pi[6];


long double gamma[T][N];
long double zeta[T][N][N];
int pi_bar[6];
long double a_bar[N][N];
long double b_bar[N][M];
int obs_seq[61];
int xsi[T][N];
int Qstar[T];
double alpha[T][N];
double beta[T][N];
double p=0;


void feed_forward_model()      
{
	for(int i=0;i<N;i++)       
		if(i==0)             
			pi[i]=1.0;
		else
			pi[i]=0;
    for(int i=0;i<N;i++)        
        for(int j=0;j<N;j++)
			if(i==j&&i!=N-1)
				a_matrix[i][j]=0.8; 
			else if(i==j&&i==N-1)
				a_matrix[i][j]=1; 
			else if(j==i+1)
				a_matrix[i][j]=0.2;
			else
				a_matrix[i][j]=0; 
    for(int i=0;i<N;i++) 
        for(int j=0;j<M;j++)
            b_matrix[i][j]=1.0/M;
}


void forward_procedure()
{
	//for (int o=0;o<1;o++)//for each observation sequence
	{
	//initialize
	for(int i=0;i<N;i++)
	{
	  alpha[0][i]=pi[i]*b_matrix[i][obs_seq[0]];
	 // cout<<"alp:"<<alpha[0][i]<<" pi:"<<pi[i]<<"bmat:"<<b_matrix[i][obs_seq[o][0]-1]<<endl;
	}

	for(int t=0;t<T-1;t++)
	{
	  for(int j=0;j<N;j++)
	  {
	   double sum=0;
	   for(int i=0;i<N;i++)
	   {
	     sum+=alpha[t][i]*a_matrix[i][j];
	   }
	   alpha[t+1][j]=sum*b_matrix[j][obs_seq[t+1]];
	   
	  }
	 
	}
	p=0;
	for(int i=0;i<N;i++) //estimate the probability
	{
	p+=alpha[T-1][i];
	 //cout<<alpha[T-1][i]<<" ";
	}
//print alpha values
/*for(int j=0;j<N;j++)	
	{
	  for(int i=0;i<T;i++)
		 {
			cout<<alpha[j][i]<<" ";			 
	  }
	 cout<<endl;;
	}
*/
	
//cout<<endl;
	cout<<"Probability:"<<p<<endl;    //probability

	}
}

void backward_procedure()
{
	//for(int o=0;o<1;o++)
	{
		{
		for (int i=0;i<N;i++)  //intialization
		{
		 beta[T-1][i]=1.0;
		}
		//induction
		for(int t=T-2;t>=0;t--)
		{
		for(int i=0;i<N;i++)
		{
		   double sum=0;
		   for(int j=0;j<N;j++)
		   {
			   sum+=a_matrix[i][j]*b_matrix[j][obs_seq[t+1]]*beta[t+1][j];		   
		   }
		   beta[t][i]=sum;
		}
		}
		//print beta values
 /*	for(int t=0;t<T;t++)
		  {
		    for(int i=0;i<N;i++)
		     {
		        cout<<beta[t][i]<<" ";
		     }
		   cout<<endl;
		   }	*/	
		}
	}


}

void call_viterbi_Algorithm()
{
 int argmax =0;
// for (int o=0;o<1;o++)     //iterate for every observation sequence
 {
	for(int i=0;i<T;i++)
	{
	   for(int j=0;j<N;j++)
	   {
		delta[i][j]=0;
	    xsi[i][j]=0;
	   }
	}
	
   for(int i=0;i<N;i++)
     {
       delta[0][i]=pi[i]*b_matrix[i][obs_seq[0]];
       xsi[0][i]=0;         //initial no state 
      }

   for (int t=1;t<T;t++)
   {
       for(int j=0;j<N;j++)
	    {
	       argmax= 0;
		   for (int i=0;i<N;i++)    //////////////
		   {
		      if((delta[t-1][i]*a_matrix[i][j])>(delta[t-1][argmax]*a_matrix[argmax][j]))
		              argmax=i;
		   }

		   delta[t][j]=delta[t-1][argmax]*a_matrix[argmax][j]*b_matrix[j][obs_seq[t]];
	       xsi[t][j]=argmax; 
	   }
   }


   argmax=0;
   for(int i=0;i<N;i++)
   {
     if (delta[T-1][i]>delta[T-1][argmax])
		 argmax=i;  
   }
   p_star = delta[T-1][argmax];

   /*cout<<"Observation sequence"<<endl;;
   for(j=0;j<T;j++)
		   cout<<obs_seq[j]<<" ";
	   cout<<endl;*/
   Qstar[T-1]=argmax;
   printf("\nP* value: %e \n",p_star);
   printf("Qstar value: %d \n",Qstar[T-1]+1);
   Qstar[T-1]=argmax;
	//cout << T-1+1 << " ---> " << Qstar[T-1]+1 << endl;
	for(int t=T-2;t>=0;t--) //back tracking the path
	{
		Qstar[t]=xsi[t+1][Qstar[t+1]];
	}
	for(int i=0;i<T;i++)
		{
			printf("%d  ",Qstar[i]+1);
		}
 }
}

void gamma_calc()
{
	int *q,argmax=0; //q--> store the state which has maximum probability of occurence at time t.
	q=new int[T];
	double devider=0; //used as devider in baye's theorem for computation of gamma
	for(int t=0;t<T;t++)
	{
		for(int i=0;i<N;i++) //compute it once for t
		{
			devider+=alpha[t][i]*beta[t][i];
		}
		argmax=0;
		for(int i=0;i<N;i++)
		{
			gamma[t][i]=alpha[t][i]*beta[t][i]/devider;
			if(gamma[t][argmax]<gamma[t][i])
				argmax=i;
		}
		q[t]=argmax;
		devider=0;
	}
	
	
/*	for(int t=0;t<T;t++)
	{
		for(int i=0;i<N;i++) 
		{
			cout<<gamma[t][i]<<" ";
		}
		cout<<endl;
	}*/

}

void baum_welch()
{
	
	double devider=0;      
	for(int t=0;t<T-1;t++)    
	{
		devider=0;
		for(int i=0;i<N;i++)
		{
			for(int j=0;j<N;j++)
				devider+=alpha[t][i]*a_matrix[i][j]*b_matrix[j][obs_seq[t+1]]*beta[t+1][j];
		}
		for(int i=0;i<N;i++)
		{
			for(int j=0;j<N;j++)
				{
				   zeta[t][i][j]=alpha[t][i]*a_matrix[i][j]*b_matrix[j][obs_seq[t+1]]*beta[t+1][j]/devider;
			     //  cout<<zeta[t][i][j]<<" ";
			    }
			//cout<<endl;
		}
		//cout<<endl;
	}

}

void re_estimation()
{
	long double count=0;
	long double numerator=0, denominator=0; //for re-estimation of transition probabilities
	for(int i=0;i<N;i++) //re-estimation of pi as pi_bar
	{	pi_bar[i]=gamma[0][i];
	  // cout<<pi_bar[i]<<endl;
	}
	for(int i=0;i<N;i++) //re-estimation of a as a_bar
	{
		for(int j=0;j<N;j++)
		{
			numerator=0;
			denominator=0;
			for(int t=0;t<T-2;t++)
			{
				numerator+=zeta[t][i][j];
				denominator+=gamma[t][i];
			}
			a_bar[i][j]=(numerator/denominator);
		}
		cout<<endl;
	}


	for(int j=0;j<N;j++) //re-estimation of b as b_bar
	{
		for(int k=0;k<M;k++)
		{
			numerator=0;
			denominator=0;
			for(int t=0;t<T;t++)
			{
				if(obs_seq[t]==k)
					numerator+=gamma[t][j];
			}
			for(int t=0;t<T;t++)
			{
				denominator+=gamma[t][j];
			}
			b_bar[j][k]=numerator/denominator;
			if ((b_bar[j][k])==0)
				b_bar[j][k]=threshold;

		}

	}

}
void replace()
{
    for(int i=0;i < N;i++) //assign the given values
        pi[i]=pi_bar[i];
    for(int i=0;i < N;i++) //assign the given values for transition probability distribution
        for(int j=0;j < N;j++)
            a_matrix[i][j]=a_bar[i][j];
    for(int i=0;i < N;i++) //assign the given values for observation symbol probability distribution
        for(int j=0;j < M;j++)
		{
		//	if(j==N-1)
			//	b_matrix[i][j]=threshold;
		//	else
            b_matrix[i][j]=b_bar[i][j];
		}
}

void displayab()    //test for checking a and b matrix
{
	long double c=0;
     for(int i=0;i<N;i++) //re-estimation of a as a_bar
	     {
			
		  for(int j=0;j<N;j++)
	     	{
		      cout<<a_matrix[i][j]<<" ";	
	    	}
		    cout<<endl;
	    }

		cout<<"b:"<<endl;
		for(int i=0;i<N;i++) 
	     {
		  for(int j=0;j<M;j++)
	     	{
		      cout<<b_matrix[i][j]<<" ";
	    	}
           }
}



void dump_modelInformation(int i,int j)
{
 char store_model[140];


	sprintf(store_model, "./Models/Digit_%d.txt",i);
	FILE *stModel;
	if(j==1)
	stModel=fopen(store_model,"w");
	else
    stModel=fopen(store_model,"a");

	fprintf(stModel,"Sample %d.......... \n",j);

	fprintf(stModel,"P star:%e \n",p_star);
	//fprintf(stModel,"P :%Lf \n",p);
	for(int i=0;i<T;i++)
		{
			fprintf(stModel," %d  ",Qstar[i]+1);
		}
	//print a matrix
	fprintf(stModel,"\nConverged A matrix ......... \n");
	 for(int i=0;i<N;i++) 
        {
			for(int j=0;j<N;j++)
             {
				 fprintf(stModel,"%lf ",a_bar[i][j]);
		     }
			fprintf(stModel,"\n");
	    }
	 //print b matrix
	 fprintf(stModel,"Converged B matrix ......... \n");
	  for(int i=0;i < N;i++) 
	{
        for(int j=0;j < M;j++)
		{
           fprintf(stModel,"%e ",b_bar[i][j]);
		}
		fprintf(stModel,"\n");
	}

	  fprintf(stModel,"                             **************************************************************\n");
	  fclose(stModel);
}



void call_functions()
{
    forward_procedure();
	backward_procedure();
	call_viterbi_Algorithm();
	gamma_calc();
	baum_welch();
	re_estimation();
}


void store_average_model(int digit)
{
  char store_a[140];
  char store_b[140];
  char store_pi[140];

	sprintf(store_a, "./average model/A_%d.txt",digit);
	sprintf(store_b, "./average model/B_%d.txt",digit);
	sprintf(store_pi, "./average model/pi_%d.txt",digit);

	FILE *a,*b,*p;
	a=fopen(store_a,"w");
	b=fopen(store_b,"w");
	p=fopen(store_pi,"w");

	for(int i=0;i < N;i++) 
       fprintf(p,"%d ",avg_pi[i]/20);
    for(int i=0;i<N;i++) 
        {
			for(int j=0;j<N;j++)
             {
				 fprintf(a,"%lf ",avg_a[i][j]/20);
		     }
			fprintf(a,"\n");
	    }
    for(int i=0;i < N;i++) 
	{
        for(int j=0;j < M;j++)
		{
		   avg_b[i][j]=avg_b[i][j]/20;
           fprintf(b,"%e ",avg_b[i][j]);
		}
		fprintf(b,"\n");
	}

	fclose(a);
	fclose(b);
	fclose(p);
	//reset
	for(int i=0;i < N;i++) 
       avg_pi[i]=0;
    for(int i=0;i < N;i++) 
        for(int j=0;j < N;j++)
            avg_a[i][j]=0;
    for(int i=0;i < N;i++) 
        for(int j=0;j < M;j++)
		{
           avg_b[i][j]=0;
		}
	
}

void addto_average()
{
	for(int i=0;i < N;i++) 
       avg_pi[i]+= pi_bar[i];
    for(int i=0;i < N;i++) 
        for(int j=0;j < N;j++)
            avg_a[i][j]+=a_bar[i][j];
    for(int i=0;i < N;i++) 
        for(int j=0;j < M;j++)
		{
           avg_b[i][j]+=b_bar[i][j];
		}
}


void train_model()
{
	cout<<"\nTraining.....\n";
	char store_obs[140];
	for(int ind=1;ind<=11;ind++)
	{
	for (int index=1;index<21;index++)
	{

	if(ind==1)
	{
	sprintf(store_obs, "./obs/ambigious/214101024_ambigious_%d.txt",ind);
	}
	else if(ind==2)
	{
	sprintf(store_obs, "./obs/artificial/214101024_artificial_%d.txt",ind);
	}
	else if(ind==3)
	{
	sprintf(store_obs, "./obs/broadband/214101024_broadband_%d.txt",ind);
	}
	else if(ind==4)
	{
	sprintf(store_obs, "./obs/coding/214101024_coding_%d.txt",ind);
	}
	else if(ind==5)
	{
	sprintf(store_obs, "./obs/covid/214101024_covid_%d.txt",ind);
	}
	else if(ind==6)
	{
	sprintf(store_obs, "./obs/deepfake/214101024_deepfake_%d.txt",ind);
	}
	else if(ind==7)
	{
	sprintf(store_obs, "./obs/internet/214101024_internet_%d.txt",ind);
	}
	else if(ind==8)
	{
	sprintf(store_obs, "./obs/mandatory/214101024_mandatory_%d.txt",ind);
	}
	else if(ind==9)
	{
	sprintf(store_obs, "./obs/podcast/214101024_podcast_%d.txt",ind);
	}
	else if(ind==10)
	{
	sprintf(store_obs, "./obs/ransomware/214101024_ransomware_%d.txt",ind);
	}
	else if(ind==11)
	{
	sprintf(store_obs, "./obs/ripple/214101024_ripple_%d.txt",ind);
	}


	FILE *inp=fopen(store_obs,"r");
	cout<<store_obs<<endl;
	feed_forward_model();
	if(!inp)
		cout<<"Error\n";
	int i=0;
	int num=0,count=0;
	while(fscanf(inp,"%d",&num)!=EOF)
	{
	obs_seq[i]=num;
	i++;
	}
	fclose(inp);
	while (count<100)
	{
	count++;
	call_functions();
	replace();               //replace matrix/update
	}
	//dump_modelInformation(ind,index);
	addto_average();
	}
//store average model for this digit
	store_average_model(ind);	
}
}



//take input from mic
void take_input_from_mic()
{
string s="Recording_Module.exe 3 ";
	string s1="live_test_sample";
	for(int i=1;i<=1;i++)
	{
		for(int j=1;j<=1;j++)
		{
			{
			char buffer[10];
			//sprintf(buffer,"%d",j);
			string s2=s1 +".txt";
			string s3=s1 +".wav";
			cout<<"\n" << s2<<endl;
			string sfin= s+s3+" "+s2;
			system(sfin.c_str());
			}
		}
	}
}

void train_data_input()
{

	read_data();
	calculate_Ci();
    generate_obs();

}


int main()
{
	
	cout<<"1.Train pre-recorded data\n";
	cout<<"2.Test on pre-recorded data\n";
	cout<<"3.Live testing\n";
cin>>choice;

	switch(choice)
	{
	  
	case 1:
		train_data_input();
	    train_model();
		getchar();
		break;
	case 2:
	  cout<<"Testing program.............\n";
      call_recorded(); 
	  getchar();
	  break;
	case 3:
	   //live testing
	  take_input_from_mic();
	  //test for the given input
	  live_testing();
	  getchar();
	  break;
	}
	cout<<"END";
	getchar();
	return 0;

}