//this file handles the testing part of the progam

#define N 5
#define M 32
#define T 60 
std::string dictionary[12]={"ambigious","artificial","broadband","coding","covid","deepfake","internet","mandatory","podcast","ransomware","ripple"};

using namespace std;
int count_check=0;    //check total correct
long double max_prob=-999;
int output_digit;
int obs_seq_test[61];
int totalDigit=0;
//int choice=0;
int digit=0;
long double prob=0;
//int obs[61];
long double avg_atest[6][6]={0};
long double avg_btest[6][33]={0};
int avg_pitest[6]={1,0,0,0,0};
double alpha_test[T][N];




//calculate Ci's for test file
void calculate_Ci_test()
{
	int start=20*320;
	int end=80*320;
	load_Hamming();
	for(int ind=1;ind<=11;ind++)
	{
     for (int file_count=1;file_count<=10;file_count++)    //for 10 samples of a particular digit
    {
    int count=0,num;
	int totalframesize[27400]={0};
	double Ai[14];
	double Ci[86][14];
	char pathW[150];
	char strCi[150];
	double alpha[14][14]={0.0};
	double E[14]={0.0};

   //sprintf(pathW, "./Normalised_test/214101024_%d.txt",ind,file_count);
   //sprintf(strCi, "./Ci_test/214101024_%d_%d.txt",ind,file_count);
   ////

   if(ind==1)
	{
	sprintf(strCi, "./Ci_test/214101024_ambigious_%d.txt",file_count);
	sprintf(pathW, "./Normalised_test/ambigious/214101024_ambigious_%d.txt",file_count);
	}
	else if(ind==2)
	{
	sprintf(strCi, "./Ci_test/214101024_artificial_%d.txt",file_count);
	sprintf(pathW, "./Normalised_test/artificial/214101024_artificial_%d.txt",file_count);
	}
	else if(ind==3)
	{
	sprintf(strCi, "./Ci_test/214101024_broadband_%d.txt",file_count);
	sprintf(pathW, "./Normalised_test/broadband/214101024_broadband_%d.txt",file_count);
	}
	else if(ind==4)
	{
	sprintf(strCi, "./Ci_test/214101024_coding_%d.txt",file_count);
	sprintf(pathW, "./Normalised_test/coding/214101024_coding_%d.txt",file_count);
	}
	else if(ind==5)
	{
	sprintf(strCi, "./Ci_test/214101024_covid_%d.txt",file_count);
	sprintf(pathW, "./Normalised_test/covid/214101024_covid_%d.txt",file_count);
	}
	else if(ind==6)
	{
	sprintf(strCi, "./Ci_test/214101024_deepfake_%d.txt",file_count);
	sprintf(pathW, "./Normalised_test/deepfake/214101024_deepfake_%d.txt",file_count);
	}
	else if(ind==7)
	{
	sprintf(strCi, "./Ci_test/214101024_internet_%d.txt",file_count);
	sprintf(pathW, "./Normalised_test/internet/214101024_internet_%d.txt",file_count);
	}
	else if(ind==8)
	{
	sprintf(strCi, "./Ci_test/214101024_mandatory_%d.txt",file_count);
	sprintf(pathW, "./Normalised_test/mandatory/214101024_mandatory_%d.txt",file_count);
	}
	else if(ind==9)
	{
	sprintf(strCi, "./Ci_test/214101024_podcast_%d.txt",file_count);
	sprintf(pathW, "./Normalised_test/podcast/214101024_podcast_%d.txt",file_count);
	}
	else if(ind==10)
	{
	sprintf(strCi, "./Ci_test/214101024_ransomware_%d.txt",file_count);
	sprintf(pathW, "./Normalised_test/ransomware/214101024_ransomware_%d.txt",file_count);
	}
	else if(ind==11)
	{
	sprintf(strCi, "./Ci_test/214101024_ripple_%d.txt",file_count);
	sprintf(pathW, "./Normalised_test/ripple/214101024_ripple_%d.txt",file_count);
	}

   ////


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
	
	for(int i=0;i<60;i++){
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
}


//read the test data and normalize it
void read_testfiles()
{
	cout<<"\nReading Digits\n";
	char pathR[150];
	char pathW[150];
	FILE *readFile;
	int num=0; 

	for(int ind=1;ind<=11;ind++)
	{
	for (int i=1;i<=10;i++)
	{
	int max=-999999;
	//////////
	if(ind==1)
	{
	sprintf(pathR, "./test/ambigious/214101024_ambigious_%d.txt",i);
	sprintf(pathW, "./Normalised_test/ambigious/214101024_ambigious_%d.txt",i);
	}
	else if(ind==2)
	{
	sprintf(pathR, "./test/artificial/214101024_artificial_%d.txt",i);
	sprintf(pathW, "./Normalised_test/artificial/214101024_artificial_%d.txt",i);
	}
	else if(ind==3)
	{
	sprintf(pathR, "./test/broadband/214101024_broadband_%d.txt",i);
	sprintf(pathW, "./Normalised_test/broadband/214101024_broadband_%d.txt",i);
	}
	else if(ind==4)
	{
	sprintf(pathR, "./test/coding/214101024_coding_%d.txt",i);
	sprintf(pathW, "./Normalised_test/coding/214101024_coding_%d.txt",i);
	}
	else if(ind==5)
	{
	sprintf(pathR, "./test/covid/214101024_covid_%d.txt",i);
	sprintf(pathW, "./Normalised_test/covid/214101024_covid_%d.txt",i);
	}
	else if(ind==6)
	{
	sprintf(pathR, "./test/deepfake/214101024_deepfake_%d.txt",i);
	sprintf(pathW, "./Normalised_test/deepfake/214101024_deepfake_%d.txt",i);
	}
	else if(ind==7)
	{
	sprintf(pathR, "./test/internet/214101024_internet_%d.txt",i);
	sprintf(pathW, "./Normalised_test/internet/214101024_internet_%d.txt",i);
	}
	else if(ind==8)
	{
	sprintf(pathR, "./test/mandatory/214101024_mandatory_%d.txt",i);
	sprintf(pathW, "./Normalised_test/mandatory/214101024_mandatory_%d.txt",i);
	}
	else if(ind==9)
	{
	sprintf(pathR, "./test/podcast/214101024_podcast_%d.txt",i);
	sprintf(pathW, "./Normalised_test/podcast/214101024_podcast_%d.txt",i);
	}
	else if(ind==10)
	{
	sprintf(pathR, "./test/ransomware/214101024_ransomware_%d.txt",i);
	sprintf(pathW, "./Normalised_test/ransomware/214101024_ransomware_%d.txt",i);
	}
	else if(ind==11)
	{
	sprintf(pathR, "./test/ripple/214101024_ripple_%d.txt",i);
	sprintf(pathW, "./Normalised_test/ripple/214101024_ripple_%d.txt",i);
	}
	///////////

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
	//load_Hamming();
	calculate_Ci_test();

}

//generate observation sequence
void generate_obs_test()
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
	for (int index=1;index<11;index++)
	{
	//cout<<endl<<endl;
	//find the minimum distance
	double Ci[86][13];
	char store_obs[140];
    i=0,j=0;

	//sprintf(store_obs, "./obs/214101024_%d_%d.txt",ind,index);
	//sprintf(obs_file, "./Ci/214101024_%d_%d.txt",ind,index);
	//cout<<store_obs<<endl;

	 if(ind==1)
	{
	sprintf(obs_file, "./Ci_test/214101024_ambigious_%d.txt",index);
	sprintf(store_obs, "./obs_test/214101024_ambigious_%d.txt",index);
	}
	else if(ind==2)
	{
	sprintf(obs_file, "./Ci_test/214101024_artificial_%d.txt",index);
	sprintf(store_obs, "./obs_test/214101024_artificial_%d.txt",index);
	}
	else if(ind==3)
	{
	sprintf(obs_file, "./Ci_test/214101024_broadband_%d.txt",index);
	sprintf(store_obs, "./obs_test/214101024_broadband_%d.txt",index);
	}
	else if(ind==4)
	{
	sprintf(obs_file, "./Ci_test/214101024_coding_%d.txt",index);
	sprintf(store_obs, "./obs_test/214101024_coding_%d.txt",index);
	}
	else if(ind==5)
	{
	sprintf(obs_file, "./Ci_test/214101024_covid_%d.txt",index);
	sprintf(store_obs, "./obs_test/214101024_covid_%d.txt",index);
	}
	else if(ind==6)
	{
	sprintf(obs_file, "./Ci_test/214101024_deepfake_%d.txt",index);
	sprintf(store_obs, "./obs_test/214101024_deepfake_%d.txt",index);
	}
	else if(ind==7)
	{
	sprintf(obs_file, "./Ci_test/214101024_internet_%d.txt",index);
	sprintf(store_obs, "./obs_test/214101024_internet_%d.txt",index);
	}
	else if(ind==8)
	{
	sprintf(obs_file, "./Ci_test/214101024_mandatory_%d.txt",index);
	sprintf(store_obs, "./obs_test/214101024_mandatory_%d.txt",index);
	}
	else if(ind==9)
	{
	sprintf(obs_file, "./Ci_test/214101024_podcast_%d.txt",index);
	sprintf(store_obs, "./obs_test/214101024_podcast_%d.txt",index);
	}
	else if(ind==10)
	{
	sprintf(obs_file, "./Ci_test/214101024_ransomware_%d.txt",index);
	sprintf(store_obs, "./obs_test/214101024_ransomware_%d.txt",index);
	}
	else if(ind==11)
	{
	sprintf(obs_file, "./Ci_test/214101024_ripple_%d.txt",index);
	sprintf(store_obs, "./obs_test/214101024_ripple_%d.txt",index);
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

void read_average_model(int iteration)
{
	char readA[150];
	char readB[150];
	//char readPi[150];


	sprintf(readA, "./average model/A_%d.txt", iteration);
	sprintf(readB, "./average model/B_%d.txt", iteration);
	//sprintf(store_pi, "./average model/pi_%d.txt", iteration);
	double num=0;

	FILE  *rA,*rB,*rPi;
	rA=fopen(readA,"r");
	rB=fopen(readB,"r");
	//rPi=fopen(readPi,"r");
	//if (!rPi )
		//cout<<"erooor\n";


	//read a matrix
	int i=0,j=0;
	while (fscanf(rA,"%lf",&num)!=EOF)
	{
	//	cout<<num<<" ";
	 avg_atest[i][j]=num;
	 j++;
	 if(j==5)
	 {
	 j=0;
	 i++;
	// cout<<endl;
	 }
	}
	fclose(rA);
/*	//read pi
	int n=0;
	i=0;
	while (fscanf(rPi,"%d",&n)!=EOF)
	{
	   avg_pitest[i]=n;
	   i++;
	}*/
	//fclose(rPi);
	ifstream in;

	in.open(readB); ////
		string temp;
	for (int i = 0; i < N; i++){
		for (int j = 0; j < M; j++){
			in >> temp;
			avg_btest[i][j] = stold(temp);
		}
	}
}


void displayab_test()    
{
	for(int i=0;i<=4;i++)
		cout<<avg_pitest[i]<<" ";
	cout<<endl;
for(int i=0;i<N;i++) 
	     {
			
		  for(int j=0;j<N;j++)
	     	{
		      cout<<avg_atest[i][j]<<" ";	
	    	}
		    cout<<endl;
	    }

		cout<<"b:"<<endl;
		for(int i=0;i<N;i++) 
	     {
		  for(int j=0;j<M;j++)
	     	{
		      cout<<avg_btest[i][j]<<" ";
	    	}
		  cout<<endl;
           }
}




void forward_procedure(int iter)
{
	//for (int o=0;o<1;o++)//for each observation sequence
	{
	//initialize
	for(int i=0;i<N;i++)
	{
	  alpha_test[0][i]=avg_pitest[i]*avg_btest[i][obs_seq_test[0]];
	 // cout<<"alp:"<<alpha[0][i]<<" pi:"<<avg_pitest[i]<<"bmat:"<<avg_b[i][obs_seq_test[0]-1]<<endl;
	}

	for(int t=0;t<T-1;t++)
	{
	  for(int j=0;j<N;j++)
	  {
	   double sum=0;
	   for(int i=0;i<N;i++)
	   {
	     sum+=alpha_test[t][i]*avg_atest[i][j];
	   }
	   alpha_test[t+1][j]=sum*avg_btest[j][obs_seq_test[t+1]];  
	  }
	}
	prob=0;
	for(int i=0;i<N;i++) //estimate the probability
	{
	prob+=alpha_test[T-1][i];
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
	if(prob>max_prob)
	{
	max_prob=prob;
	output_digit=iter;	  /////////////
	}
	//print alpha values
/*for(int j=0;j<N;j++)	
	{
	  for(int i=0;i<T;i++)
		 {
			cout<<alpha[j][i]<<" ";			 
	  }
	 cout<<endl;
	}
*/

//cout<<endl;
	cout<<"Probability for model "<<iter<<" is:"<<prob<<endl;    //probability

	}
}

void load_obs(int i,int j)
{
	char obs[140];
	//sprintf(obs, "./obs/214101024_%d_%d.txt",i,j);


	 if(i==1)
	{
	sprintf(obs, "./obs_test/214101024_ambigious_%d.txt",j);
	}
	else if(i==2)
	{
	sprintf(obs, "./obs_test/214101024_artificial_%d.txt",j);
	}
	else if(i==3)
	{
	sprintf(obs, "./obs_test/214101024_broadband_%d.txt",j);
	}
	else if(i==4)
	{
	sprintf(obs, "./obs_test/214101024_coding_%d.txt",j);
	}
	else if(i==5)
	{
	sprintf(obs, "./obs_test/214101024_covid_%d.txt",j);
	}
	else if(i==6)
	{
	sprintf(obs, "./obs_test/214101024_deepfake_%d.txt",j);
	}
	else if(i==7)
	{
	sprintf(obs, "./obs_test/214101024_internet_%d.txt",j);
	}
	else if(i==8)
	{
	sprintf(obs, "./obs_test/214101024_mandatory_%d.txt",j);
	}
	else if(i==9)
	{
	sprintf(obs, "./obs_test/214101024_podcast_%d.txt",j);
	}
	else if(i==10)
	{
	sprintf(obs, "./obs_test/214101024_ransomware_%d.txt",j);
	}
	else if(i==11)
	{
	sprintf(obs, "./obs_test/214101024_ripple_%d.txt",j);
	}


	////////// 
	//cout<<obs;
	int ind=0,num=0;
	FILE *inp ;
	inp=fopen(obs,"r");
	if(!inp)
		cout<<"\nError opening obs file\n";
	while(fscanf(inp,"%d",&num)!=EOF)
	{
		//cout<<num<<" ";
		obs_seq_test[ind]=num;
		ind++;
		//cout<<obs_seq_test[ind-1];
	}

	fclose(inp);

}

void check_detection(int digit){	
	if (digit == output_digit){
		count_check++;
	}
}


void call_recorded()
{
        // read_testfiles();
         // generate_obs_test();
	
	//test
	     int index=0;
	     for(int i=1;i<=11;i++)
	        {
		     cout<<"Input word:"<<dictionary[i-1]<<endl;
		     for(int j=1;j<=10;j++)
		       {  
		         totalDigit++;
		         cout<<"Input word:"<<dictionary[i-1]<<" sample "<<j<<endl;
		         max_prob=-9;
		         output_digit=-1;
		         load_obs(i,j);
		      for(int k=1;k<12;k++)
		       {
         	     read_average_model(k);
			   // displayab();
			   // cout<<"/////////////\n";
			      forward_procedure(k);
			   //cout<<"MaxP:"<<max_prob<<endl;
		        }
		       cout<<"Input word:"<<dictionary[i-1]<<endl;
		       cout<<"Output:"<<dictionary[output_digit-1]<<endl<<endl<<endl;
		       check_detection(i);
		        }
			 cout<<endl<<endl;
	         }
         cout<<"Total test cases:"<<totalDigit<<endl;
         cout<<"Correctly predicted total:"<<count_check<<endl;
}


void live_testing()
{
	//read the input digit for the live testing part
	cout<<"\nReading Digits\n";
	char pathR[150];
	char pathW[150];
	FILE *readFile;
	int num=0; 

	int max=-999999;
	sprintf(pathR, "live_test_sample.txt");
	sprintf(pathW, "Normalised_test_sample.txt");

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
			fprintf(normalizedValues,"%d\n",temp1);		
	}
	fclose(normalizedValues);
	//calculate ci


	int start=20*320;
	int end=80*320;
	load_Hamming();

    int count=0;
	int totalframesize[27400]={0};
	double Ai[14];
	double Ci[86][14];
//	char pathW[150];
	char strCi[150];
	double alpha[14][14]={0.0};
	double E[14]={0.0};

 	sprintf(strCi, "Ci_test.txt");
	sprintf(pathW, "Normalised_test_sample.txt");

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
	
	for(int i=0;i<60;i++){
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
  

//generate observation sequence 

   cout<<"\nGenerating observation sequence\n";
	double wt[13]={1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
	double codebook[33][12];
	double n=0;
	FILE *rdcb = fopen("codeBook.txt","r");

	 int i=0;
	 j=0;
	while(fscanf(rdcb,"%lf",&n)!=EOF)
	{
	   codebook[i][j]=n;
	 //  cout<<n;
	   j++;
	   if(j==12)
	   {
	     j=0;
		 i++;
	   }
	}
	char obs_file[140];
	//find the minimum distance
//	double Ci[86][13];              //////
	char store_obs[140];
    i=0,j=0;
	sprintf(obs_file, "Ci_test.txt");
	sprintf(store_obs, "obs_test.txt");

	FILE *inp ;
	inp=fopen(obs_file,"r");
	///////
	double why=0;
	while(fscanf(inp,"%lf",&why)!=EOF)
	{
		Ci[i][j]=why;
	//	cout<<why<<" ";
		j++;
		if(j==12)
		{
		i++;
		j=0;
		}
	}
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
	ofstream out7;
	out7.open(store_obs);
	for(int i=0;i<60;i++)
		{
			obs_seq_test[i]=obs[i];    //
			//cout<<obs[i]<<" ";  //print observation sequence
			out7<<obs[i]<<" ";
	    }
	out7.close();
///////////

//prediction  of word 
//compare with all the models
	for(int k=1;k<12;k++)
       {
          read_average_model(k);
	      forward_procedure(k);
	   //cout<<"MaxP:"<<max_prob<<endl;
       }
       cout<<"Output:"<<dictionary[output_digit-1]<<endl<<endl<<endl;
     //  check_detection(i);



}
