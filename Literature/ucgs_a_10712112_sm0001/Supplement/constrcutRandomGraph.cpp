

#include<iostream.h>
#include<iomanip.h>
#include<fstream.h>			
#include<stdlib.h>          
#include<ctime>             
#include<math.h> 
#include<string.h>
#include <stdio.h> 

#include"graph.h"



void main()


{
	my_rand_seeds();
	clock_t start, finish;
	start=clock();
	

	//The number of graphcial models for each case
	int num=5;
	//The number of vertices of graph
	int num_vert[4]={50, 100, 200, 300};
	//the probability of adding edge between two vertices
	double prob[6]={0.2, 0.15, 0.1,  0.05,	0.01, 0.005};
	char * prob_grate[]={"0.2", "0.15", "0.1", "0.05", "0.01", "0.005"};


	for(int i=0; i<4; i++)
	{//int i=1;
		for(int j=2; j<6; j++)
		{//int j=4;
			
			//for(int k=0; k<1; k++){	
			for(int k=0; k<num; k++){			
				//Generate randomly a connect graph G by "g3(num_vert[i], prob[j])".
				ALGraph<int> g3(num_vert[i], prob[j]);	
				//ALGraph<int> g3(num_vert[i], "cycle");	
				char graphName[20];
				sprintf(graphName, "%d", num_vert[i]);
				strcat(graphName, "_");

				char GraphName[60]="graphs/";
				strcat(GraphName, graphName);
			
				strcat(GraphName, prob_grate[j]);
				strcat(GraphName, "_");
			
				char num[10];
				sprintf(num, "%d", k+1);
				strcat(GraphName, num);			
				strcat(GraphName, ".txt");

				g3.Output(GraphName);
			}


		}
	}

	finish=clock();
	double totaltime=double(finish-start)/CLOCKS_PER_SEC;	
	cout<<"total time is "<<totaltime<<" seconds\n";

}




