#include <fstream>
int obs[61];
double hamming[322];   //store values from the pre-calculated file
//double universe[18000][13];
using namespace std;
int check=0;
double codeBook[32][13]={0};
//read the training data and normalize it
void read_data()
{
	cout<<"Reading training data......\n";
	char pathR[150];
	char pathW[150];
	char pathW2[150];
	FILE *readFile;
	int num=0; 
	//std::string dictionary[11]={"ambigious","artificial","broadband","coding","covid","deepfake","internet","mandatory","podcast","ransomware","ripple"};

	for(int i=1;i<=11;i++)
	{
	for (int j=1;j<=20;j++)
	{
	int max=-999999;
	if(i==1)
	{
	sprintf(pathR, "./words/ambigious/214101024_ambigious_%d.txt",j);
	sprintf(pathW, "./Normalised/ambigious/214101024_ambigious_%d.txt",j);
	}
	else if(i==2)
	{
	sprintf(pathR, "./words/artificial/214101024_artificial_%d.txt",j);
	sprintf(pathW, "./Normalised/artificial/214101024_artificial_%d.txt",j);
	}
	else if(i==3)
	{
	sprintf(pathR, "./words/broadband/214101024_broadband_%d.txt",j);
	sprintf(pathW, "./Normalised/broadband/214101024_broadband_%d.txt",j);
	}
	else if(i==4)
	{
	sprintf(pathR, "./words/coding/214101024_coding_%d.txt",j);
	sprintf(pathW, "./Normalised/coding/214101024_coding_%d.txt",j);
	}
	else if(i==5)
	{
	sprintf(pathR, "./words/covid/214101024_covid_%d.txt",j);
	sprintf(pathW, "./Normalised/covid/214101024_covid_%d.txt",j);
	}
	else if(i==6)
	{
	sprintf(pathR, "./words/deepfake/214101024_deepfake_%d.txt",j);
	sprintf(pathW, "./Normalised/deepfake/214101024_deepfake_%d.txt",j);
	}
	else if(i==7)
	{
	sprintf(pathR, "./words/internet/214101024_internet_%d.txt",j);
	sprintf(pathW, "./Normalised/internet/214101024_internet_%d.txt",j);
	}
	else if(i==8)
	{
	sprintf(pathR, "./words/mandatory/214101024_mandatory_%d.txt",j);
	sprintf(pathW, "./Normalised/mandatory/214101024_mandatory_%d.txt",j);
	}
	else if(i==9)
	{
	sprintf(pathR, "./words/podcast/214101024_podcast_%d.txt",j);
	sprintf(pathW, "./Normalised/podcast/214101024_podcast_%d.txt",j);
	}
	else if(i==10)
	{
	sprintf(pathR, "./words/ransomware/214101024_ransomware_%d.txt",j);
	sprintf(pathW, "./Normalised/ransomware/214101024_ransomware_%d.txt",j);
	}
	else if(i==11)
	{
	sprintf(pathR, "./words/ripple/214101024_ripple_%d.txt",j);
	sprintf(pathW, "./Normalised/ripple/214101024_ripple_%d.txt",j);
	}
	//cout<<pathR<<endl;
	
	readFile=fopen(pathR,"r");
	if(!readFile)
		printf("Cannot open file\n");
	//find the maximum
	double dc=0;
	int index=0;
	while (fscanf(readFile,"%d",&num)!=EOF)
	{
		index++;
	   if(index<320*6)
			dc+=num;

	    if(abs(num)>max)
		   max=abs(num);
	}

	//cout<<dc/6<<" is the dc\n";
	dc/=6;
	rewind(readFile);
	FILE *normalizedValues;
	normalizedValues = fopen(pathW,"w");
	if (!normalizedValues)
		cout<<"error\n";
	int temp1=0;
	//ofstream out(pathW);
	while (fscanf(readFile,"%d",&num)!=EOF)      //normalization of values
	{
	    	temp1=((num-dc)*5000)/max ;
		 
		//	cout<<temp1<<endl;
			fprintf(normalizedValues,"%d\n",temp1);
			//out<<temp1<<endl;
			
	}
	fclose(normalizedValues);
    }

	}

}
//load codebook
void load_codebook()
{
	cout<<"Loading CodeBook\n";
	FILE *in = fopen("codebook.txt","r");
	long double num=0;
	int i=0;
	int j=0;
	while(fscanf(in,"%Lf",&num)!=EOF)
	{
		codeBook[i][j]=num;
		j++;
		if(j==12)
		{
		i++;
		j=0;
		}
	
	}
}

//laod Hamming
void load_Hamming()
{
FILE *inp=fopen("Hamming_window.txt","r");
int i=0;
double num=0;
while(fscanf(inp,"%lf",&num)!=EOF)
		{                                    
			hamming[i]=num;      //store in array whole 85 frame
           i++;
		 //  cout<<hamming[i-1]<<" ";
	}
}



//calculate Ci's
void calculate_Ci()
{
	cout<<"Calculating Ci values.........\n";
	int start=20*320;
	int end=80*320;
	load_Hamming();

	for(int ind=1;ind<=11;ind++)
	{
    for (int file_count=1;file_count<=20;file_count++)    //for 20 samples of a particular digit
    {
    int count=0,num;
	int totalframesize[27400]={0};
	double Ai[14];
	double Ci[86][14];
	char pathW[150];
	char strCi[150];
	double alpha[14][14]={0.0};
	double E[14]={0.0};

   if(ind==1)
	{
	sprintf(strCi, "./Ci/ambigious/214101024_ambigious_%d.txt",file_count);
	sprintf(pathW, "./Normalised/ambigious/214101024_ambigious_%d.txt",file_count);
	}
	else if(ind==2)
	{
	sprintf(strCi, "./Ci/artificial/214101024_artificial_%d.txt",file_count);
	sprintf(pathW, "./Normalised/artificial/214101024_artificial_%d.txt",file_count);
	}
	else if(ind==3)
	{
	sprintf(strCi, "./Ci/broadband/214101024_broadband_%d.txt",file_count);
	sprintf(pathW, "./Normalised/broadband/214101024_broadband_%d.txt",file_count);
	}
	else if(ind==4)
	{
	sprintf(strCi, "./Ci/coding/214101024_coding_%d.txt",file_count);
	sprintf(pathW, "./Normalised/coding/214101024_coding_%d.txt",file_count);
	}
	else if(ind==5)
	{
	sprintf(strCi, "./Ci/covid/214101024_covid_%d.txt",file_count);
	sprintf(pathW, "./Normalised/covid/214101024_covid_%d.txt",file_count);
	}
	else if(ind==6)
	{
	sprintf(strCi, "./Ci/deepfake/214101024_deepfake_%d.txt",file_count);
	sprintf(pathW, "./Normalised/deepfake/214101024_deepfake_%d.txt",file_count);
	}
	else if(ind==7)
	{
	sprintf(strCi, "./Ci/internet/214101024_internet_%d.txt",file_count);
	sprintf(pathW, "./Normalised/internet/214101024_internet_%d.txt",file_count);
	}
	else if(ind==8)
	{
	sprintf(strCi, "./Ci/mandatory/214101024_mandatory_%d.txt",file_count);
	sprintf(pathW, "./Normalised/mandatory/214101024_mandatory_%d.txt",file_count);
	}
	else if(ind==9)
	{
	sprintf(strCi, "./Ci/podcast/214101024_podcast_%d.txt",file_count);
	sprintf(pathW, "./Normalised/podcast/214101024_podcast_%d.txt",file_count);
	}
	else if(ind==10)
	{
	sprintf(strCi, "./Ci/ransomware/214101024_ransomware_%d.txt",file_count);
	sprintf(pathW, "./Normalised/ransomware/214101024_ransomware_%d.txt",file_count);
	}
	else if(ind==11)
	{
	sprintf(strCi, "./Ci/ripple/214101024_ripple_%d.txt",file_count);
	sprintf(pathW, "./Normalised/ripple/214101024_ripple_%d.txt",file_count);
	}



   int index_ci=0;
   FILE *norm =fopen(pathW,"r");
   int j=0;
	while(fscanf(norm,"%d",&num)!=EOF)
		{                                  
        count++;
		if(count>=start && count<end)
		{
			totalframesize[j]=num;      //store in array whole 85 frame
             j++;
			//cout<<totalframesize[j-1]<<endl;
		}
	}
///////////////////////////////
	//for(int i=0;i<27200;i++)
		//cout<<totalframesize[i]<<" ";
	//calculate Ri
	
	int tValues[322];
	int x=0;

		for (int frame=1;frame<61;frame++)     
		{
		
		double Ri[14]={0};
		//load 1 frame
		for(int i=0;i<320;i++)
		{
			tValues[i]=totalframesize[x]*hamming[i];  //apply hamming at the same time
	        x++;
		}

		for(int i=0;i<=12;i++) //p=12
		{
		 for(int j=0;j<320-i;j++)    //less then N-i
		 {
			  Ri[i]+=tValues[j]*tValues[j+i];

		 }
		}
		//calculate Ai

	double temp=0;
	double k[14]={0};
	E[0]=Ri[0];
	for(int i=1;i<=12;i++)
	{
	temp=0;
	for(int j=1;j<=i-1;j++)
	{
	   temp+=(alpha[i-1][j]*Ri[i-j]);
	}
	
	k[i]=(Ri[i]-temp)/E[i-1];
	
	alpha[i][i]=k[i];

	for(int j=1;j<=i-1;j++)
	{
	alpha[i][j] = alpha[i-1][j] - (k[i]*alpha[i-1][i-j]);
	}
	E[i]=(1-k[i]*k[i])*E[i-1];
	}
     for(int j=1;j<=12;j++)
     Ai[j]=alpha[12][j]; 

	 //calculate Ci
	Ci[index_ci][0]=log(Ri[0]*Ri[0]); //C0 = log(sigma)^2        sigma = R0     gain term of LPC             
	int m=0;
	for (m=1;m<=12;m++)
	{
		temp=0.0;
		for(int i=1;i<=m-1;i++)
		{
		    temp+=(i*Ci[index_ci][i]*Ai[m-i])/m;
		}
		Ci[index_ci][m]=Ai[m]+temp;	
	}
	index_ci++;

	}
/////////////////
	int sinee[17]={0};
	for(int m=1;m<=12;m++)
	{
	   sinee[m]=1+ (6*sin((3.14*m)/12));
	}
	
	for(int i=0;i<85;i++){
	for(int m=1;m<=12;m++)
	{
	Ci[i][m]=Ci[i][m]*sinee[m];                                               //apply sine window
	}
	}

	//  FILE *outp = fopen(strCi,"w");
	 // cout<<strCi<<endl;
	  ofstream out;
	  out.open(strCi);
      int check=0;
	  for(int i=0;i<60;i++){
		for(int j=1;j<=12;j++)
	     {
		//  check++;
		 // printf("%lf ",Ci[i][j]);                              
		 // fprintf(outp,"%lf ",Ci[i][j]);
		  out<<Ci[i][j]<<" ";
	      }
		// cout<<endl;
		//fprintf(outp,"\n");
		out<<"\n";
	  }
	  out.close();
	 // cout<<check;
  }
  }
  cout<<"Ci values calculated\n";
}





//generate observation sequence


void generate_obs()
{
	cout<<"\nGenerating observation sequence\n";
	double wt[13]={1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
	double codebook[33][12];
	double num=0;
	FILE *rdcb = fopen("codeBook.txt","r");

	int i=0,j=0;
	while(fscanf(rdcb,"%lf",&num)!=EOF)
	{
	   codebook[i][j]=num;
	   //cout<<num;
	   j++;
	   if(j==12)
	   {
		   //cout<<endl;
	     j=0;
		 i++;
	   }
	}

	char obs_file[140];

	for(int ind=1;ind<=11;ind++)
	{
	for (int index=1;index<21;index++)
	{
	//cout<<endl<<endl;
	//find the minimum distance
	double Ci[86][13];
	char store_obs[140];
    i=0,j=0;                           //////


	if(ind==1)
	{
	sprintf(obs_file, "./Ci/ambigious/214101024_ambigious_%d.txt",index);
	sprintf(store_obs, "./obs/ambigious/214101024_ambigious_%d.txt",index);
	}
	else if(ind==2)
	{
	sprintf(obs_file, "./Ci/artificial/214101024_artificial_%d.txt",index);
	sprintf(store_obs, "./obs/artificial/214101024_artificial_%d.txt",index);
	}
	else if(ind==3)
	{
	sprintf(obs_file, "./Ci/broadband/214101024_broadband_%d.txt",index);
	sprintf(store_obs, "./obs/broadband/214101024_broadband_%d.txt",index);
	}
	else if(ind==4)
	{
	sprintf(obs_file, "./Ci/coding/214101024_coding_%d.txt",index);
	sprintf(store_obs, "./obs/coding/214101024_coding_%d.txt",index);
	}
	else if(ind==5)
	{
	sprintf(obs_file, "./Ci/covid/214101024_covid_%d.txt",index);
	sprintf(store_obs, "./obs/covid/214101024_covid_%d.txt",index);
	}
	else if(ind==6)
	{
	sprintf(obs_file, "./Ci/deepfake/214101024_deepfake_%d.txt",index);
	sprintf(store_obs, "./obs/deepfake/214101024_deepfake_%d.txt",index);
	}
	else if(ind==7)
	{
	sprintf(obs_file, "./Ci/internet/214101024_internet_%d.txt",index);
	sprintf(store_obs, "./obs/internet/214101024_internet_%d.txt",index);
	}
	else if(ind==8)
	{
	sprintf(obs_file, "./Ci/mandatory/214101024_mandatory_%d.txt",index);
	sprintf(store_obs, "./obs/mandatory/214101024_mandatory_%d.txt",index);
	}
	else if(ind==9)
	{
	sprintf(obs_file, "./Ci/podcast/214101024_podcast_%d.txt",index);
	sprintf(store_obs, "./obs/podcast/214101024_podcast_%d.txt",index);
	}
	else if(ind==10)
	{
	sprintf(obs_file, "./Ci/ransomware/214101024_ransomware_%d.txt",index);
	sprintf(store_obs, "./obs/ransomware/214101024_ransomware_%d.txt",index);
	}
	else if(ind==11)
	{
	sprintf(obs_file, "./Ci/ripple/214101024_ripple_%d.txt",index);
	sprintf(store_obs, "./obs/ripple/214101024_ripple_%d.txt",index);
	}



	FILE *inp ;
	inp=fopen(obs_file,"r");
	while(fscanf(inp,"%lf",&num)!=EOF)
	{
		Ci[i][j]=num;
		j++;
		if(j==12)
		{
		i++;
		j=0;
		}
	}
	fclose(inp);
	//check cb and frame
/*
	for(int i=0;i<85;i++)
	{
	  for(int j=0;j<12;j++)
		  cout<<Ci[i][j]<<" ";
	  cout<<endl;
	}
	cout<<"codebook\n";

	for(int i=0;i<32;i++)
	{
	  for(int j=0;j<12;j++)
		  cout<<codebook[i][j]<<" ";
	  cout<<endl;
	}

	*/


	int cnt=0;
	int minIndex=0;
	for (int i=0;i<60;i++)
	{
	//	cout<<i<<"frame\n";
	  double min=999999;
	  double diff=0,d=0;
	  for(int j=0;j<32;j++)
	  {
		  d=0;
		  diff=0;
	    for(int k=0;k<12;k++)
		{
		  diff=Ci[i][k]-codebook[j][k];
		  d+=wt[k]*diff*diff;
		}
		//cout<<"d value:"<<d<<" jvalue"<<j<<endl;
		  if(d<min)         /////////
		{
			min=d;
			minIndex=j;
			//cout<<"obs:"<<obs[j]<<" min:"<<min<<" j value:"<<j<<endl; ////////////////			
	     }
		}
	  obs[i]=minIndex+1;
	  //cout<<"min:"<<min<<endl;
	}
	
	//cout<<endl<<endl;
	ofstream out;
	out.open(store_obs);
	for(int i=0;i<60;i++)
		{
			//cout<<obs[i]<<" ";  //print observation sequence
			out<<obs[i]<<" ";
	    }
	out.close();

}
}
}