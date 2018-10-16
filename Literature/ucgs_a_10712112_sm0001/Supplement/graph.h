//Mcs_M's madj set is adjusted to reduce the memory, 
//But Lex_M is not, we will do that later
//#include "matlib.h"
//#include "ips_load_cliques.h"

#include<iostream.h>
#include<iomanip.h>
#include<fstream.h>			
#include<stdlib.h>          
#include<ctime>             
#include<math.h> 
#include<string>
#include <stdio.h>


int n_d=0;

const int M=3000;//maximun number of vertex

struct ArcNode
{
	int adjvex;
	ArcNode * next;
};

template <class T>
struct VertexNode 
{
	T vertex;
	ArcNode * firstedge;
	int num;
};



int comp[M];//connected component





/******************************************************************************/
//global variables for Mcs_M 
struct Madj
{

	ArcNode * first;
	ArcNode * tail;

};

	
int McsM_order[M];
int weight[M], weight_update[M];

struct chain1 //save element of the sets madj[v_i] and reach[j]
{
	ArcNode * first;
};

Madj madj[M]; // monotony adjcent neighbors, global variables
chain1 reach[M];//reach(j) set, the notation in the orginal paper

int Lex_order[M];
int labeled[M]; //if vertex i is numbered, labeled[i-1]=its number, otherwise =0



int num_blocks;

time_t t_lex, t_Mcs; 

time_t lf1, lf2, mf1, mf2;	int lfnum=0, mfnum=0;
time_t lc1, lc2, mc1, mc2;	double lc=0.0, mc=0.0;
time_t ld1, ld2, md1, md2;  double ld=0.0, md=0.0;
time_t lco1, lco2, mco1, mco2; double lco=0.0, mco=0.0;
time_t ldg1, ldg2, mdg1, mdg2; double ldg=0.0, mdg=0.0;

time_t lg1, lg2, mg1, mg2;	double lg=0.0, mg=0.0;

/******************************************************************************/
// the global variables used to find all cliques 
int compsub[M];
int ccc;
int number_cliques;
/*
struct chain
{
	int num;
	chain *next;
};*/
//chain * pc1, * pc2, * pc3, * pc4;

ArcNode * pc1, * pc2, * pc3, * pc4;

/******************************************************************************/




template <class T>
class ALGraph
{
public:
	ALGraph(double a[]);   //construct a graph from shuzu, i.e., a[]  
	ALGraph(char * p);//construt a graph from disk         
	ALGraph(int num); ///construct a tree, where num is the number of vertices, but must smaller than M
	ALGraph(int num, double prob);//construct a connected graph from a tree, where num is the number of vertices, ，prob is the probability of adding edges
    ALGraph(int num, char *cycle); //a tree////construct an n_cycle
	int In_edge(int k, int i);//called by ALGraph(int num, double prob), judge wether k is in the adjacent of i
	ALGraph(const ALGraph<T>&);//copy function
	
	~ALGraph();

	void Num(){cout<<"vertexNum="<<vertexNum<<"  "<<"arcNum="<<arcNum<<endl;}
	int  Num_of_vertex(){ return vertexNum;}//return the number of vertices
	int  Num_of_edge(){return arcNum;}//return the number of edge
	void Print(int i); // print all parents of i to screem
	void Print();//print all parents of all vertices
	void Output();// output the graph to disk
	void Output(char * name);//output the graph to disk, name.txt
	void Output_to_matlab(int number);//this function is simllar to Output(char * name)
							//the input "number" is the same 
							//as the parameter of input_graph(cases, number)
	int * Output_to_matlab();//
	int * GraphToArray( );
	void OutputVertex();//output vertices
	void Compchain(int i);

	int  Chain(int m, int n);

/***************************************************************************************************************************************/

	//Mcs_M copy the paper of Lex_M
	void Mcs_M2(ALGraph<T> & g0, int largest);
	void Decom_McsM_4(ALGraph<T> & g0,   ArcNode * f_t, int n);////don't copy g1 and g2
	void Print_block(int number);///called by Decom_McsM_4(ALGraph<T> & g) to print all prime blocks to screen
									//number is number_th prime block
	ArcNode * Mcs_M5(double subset[], int num);//changed from Mcs_M5(int largest);
	                                        //it is called to construct a junction tree whose first node containing subset[]
    ArcNode * Mcs_M52(double subset[], int num);
	ArcNode * Mcs_M5(int largest);//return a point of f_ts//reduce the memory of madj[i];  //largest is largest numbered vertex
	void Decom_McsM_5(ALGraph<T> & g0,   ArcNode * f_t, int n);
/***************************************************************************************************************************************/
	//copy the orginal paper of Lex_M and compile the procedure
	void Lex_M(ALGraph<T> & g0, int largest);//largest is the numbered vertex

	
	void Decom_LexM_4(ALGraph<T> & g0, ArcNode * e, ArcNode * f, int num);//don't copy g1 and g2



/***************************************************************************************************************************************/
	int  adjacent(int i, int j);////for a graph g, judge whether i and j are adjacent, 
										//if i=j or i,j are adjacent, then return 1, otherwise return 0
	void All_Cliques(ALGraph<T> & g);//find all cliques of undirected graph
	void All_Cliques();//find all cliques of undirected graph
	void All_Cliques_To_Matlab();//find all cliques_To_Matlab() of undirected graph
	void extend_version_2(int old[], int ne, int ce, ALGraph<T> & g);//called by ALL cliques
	void extend_version_2(int old[], int ne, int ce);//called by ALL cliques
	void extend_version_2_To_Matlab(int old[], int ne, int ce);//called by ALL cliques_To_Matlab()
	void extend_version_2(int old[], int ne, int ce, ALGraph<T> & g, ArcNode *num_vertex_cliques);

//for finding all cliques of prime block
	void All_Cliques_1(ALGraph<T> & g, int number);//find all cliques of prime block
	void All_Cliques_1(int number);//find all cliques of prime block

	void junction_clique(int n);
	/***************************************************************************************************************************************/
	void  jun_tree();//construct junction tree of global graph after order the vertices by ArcNode * Mcs_M5(int largest);
	
private:
	VertexNode <T> adjlist[M]; 
	int vertexNum, arcNum;// the numbers of vertices and edges
};

template<class T>                
ALGraph<T>::ALGraph(double a[])// a[] is a shuzu,// dim is shuzu's dimension
{
	//	VertexNode <T> adjlist[M];	
	
	int  j,  e=0, i=0;
//	infile>>vertexNum;
//	infile>>arcNum;
	vertexNum = int(a[i]);	i++;
	arcNum = int(a[i]);	 i++;

//	cout<<"vertexNum="<<vertexNum<<endl;
	struct ArcNode * s, *s1;



	for(int k=0; k<vertexNum; k++)  
	{
	//	infile>>adjlist[k].vertex;
		adjlist[k].vertex = int(a[i]); i++;
		adjlist[k].firstedge=NULL;
		adjlist[k].num=0;

//		infile>>j;  
		j = int(a[i]);	i++;  
		int m=0;
		while(j!=0)//(j>k)
		{
			e++;  
			s= new ArcNode; n_d++;
			s->adjvex=j;
			adjlist[k].num++;
			if(m==0)
			{
				adjlist[k].firstedge=s;
				m++;
			}
			else 
				s1->next=s;
			s1=s;
	//		infile>>j;
			j = int(a[i]);	i++;  
		}
		if(m==0)
			adjlist[k].firstedge=NULL;
		else
			s1->next=NULL;

	}//	arcNum=e;	  cout<<"arcNum="<<arcNum<<endl; cout<<"vertexNum="<<vertexNum<<endl;


}

template<class T>                
ALGraph<T>::ALGraph(char *p)// p is a pointer pointing to filename
{
	//	VertexNode <T> adjlist[M];	
	ifstream infile(p, ios::in);
	if(infile.fail())
	{
		cout<<"file no exist!"<<endl;
		return;
	}
	int j,  e=0;
	infile>>vertexNum;
	infile>>arcNum;
	cout<<"vertexNum="<<vertexNum<<endl;
	struct ArcNode * s, *s1;



	for(int k=0; k<vertexNum; k++)  
	{
		infile>>adjlist[k].vertex;
		adjlist[k].firstedge=NULL;
		adjlist[k].num=0;

		infile>>j;  
		int m=0;
		while(j!=0)//(j>k)
		{
			e++;  
			s= new ArcNode; n_d++;
			s->adjvex=j;
			adjlist[k].num++;
			if(m==0)
			{
				adjlist[k].firstedge=s;
				m++;
			}
			else 
				s1->next=s;
			s1=s;
			infile>>j;
		}
		if(m==0)
			adjlist[k].firstedge=NULL;
		else
			s1->next=NULL;

	}//	arcNum=e;	  cout<<"arcNum="<<arcNum<<endl; cout<<"vertexNum="<<vertexNum<<endl;
	infile.close();
}



template<class T>
ALGraph<T>::ALGraph(const ALGraph<T>&g)//copy function
{
	vertexNum=g.vertexNum;// why?
	arcNum=g.arcNum;
	for(int i=0; i< vertexNum; i++)
	{
		adjlist[i].vertex=g.adjlist[i].vertex; 
		adjlist[i].num=g.adjlist[i].num;
		ArcNode *bs, *bs1, *as;
		as=g.adjlist[i].firstedge;
		if(as!=NULL)
		{
			bs=new ArcNode; n_d++;
			bs->adjvex=as->adjvex;
			adjlist[i].firstedge=bs;
			bs1=bs;
			as=as->next;		
			while(as!=NULL)
			{
				bs=new ArcNode; n_d++;
				bs->adjvex=as->adjvex;
				bs1->next=bs;
				bs1=bs;
				as=as->next;

			}
			bs->next=NULL;

		}
		else
		{
			adjlist[i].firstedge=NULL;
		}
	}
//	cout<<"n_d="<<n_d<<endl;
}

template<class T>
ALGraph<T>::~ALGraph()   
{
//	int e=arcNum;
//	cout<<"e="<<e<<endl;
	ArcNode *p, *q;
	for(int k=0; k< vertexNum; k++)
	{
		p=adjlist[k].firstedge;
		while(p)
		{
			q=p;
			p=p->next;
			delete q; n_d--;
//			e--;
		}
		adjlist[k].firstedge=NULL;
	}
//	cout<<"e="<<e<<endl;
//	cout<<"jieshu"<<endl;
//	cout<<"n_d="<<n_d<<endl;
}


template <class T>
void ALGraph<T>::Print(int i)//print all parents of vertex i on the screen 
{
	struct ArcNode * p;
	//	if(adjlist[i-1].vertex != 0)
	{
		cout<<adjlist[i-1].vertex<<" ";         
		p= adjlist[i-1].firstedge;
		while(p!=NULL)
		{
		
			cout<<p->adjvex<<" ";
			p=p->next;
		}
		cout<<endl;
	}
}


void Print(ArcNode * head)
{
	cout<<endl;
	ArcNode * p;
	p=head;
	while(p !=NULL)
	{
		cout<<setw(2)<<p->adjvex<<" ";
		p=p->next;
	}
	cout<<endl;
}

template <class T>
void ALGraph<T>::Print()//print all parents of all vertices on the screen
{
	int n=vertexNum;
	for(int i=1; i<=vertexNum/*vertexNum*/; i++)
	{
		if(adjlist[i-1].vertex != 0)
		{
			Print(i);
			cout<<i<<" has "<<adjlist[i-1].num<<" edges.\n";
		}
		else
			n--;
	}//	cout<<endl;
	cout<<"vertexNum="<<n<<"  "<<"arcNum="<<arcNum<<endl;
}
template <class T>
void ALGraph<T>::Output()
{

	fstream outfile;      
	outfile.open("a5.txt", ios::app);
	if(!outfile)
	{
		cout<<"can,t open.\n";
		abort();
	}
	outfile<<"vertexNum="<<vertexNum<<"  "<<"arcNum="<<arcNum<<endl;
	cout<<"vertexNum="<<vertexNum<<"  "<<"arcNum="<<arcNum<<endl;                    
	struct ArcNode *p;

	for(int i=0; i<vertexNum; i++)
	{
		if(adjlist[i].vertex != 0)
		{                  
			outfile<<adjlist[i].vertex<<" ";
			p=adjlist[i].firstedge;
			while(p!=NULL)
			{
				outfile<<p->adjvex<<" ";
				p=p->next;
			}
			outfile<<endl;//there isn't 0 at the end
		}
		
	}
	outfile<<endl;
	
	outfile.close();
}	

template <class T>
void ALGraph<T>::Output(char * name)//name is file name
{
	fstream outfile;      
	outfile.open(name, ios::app);//outfile.open(name, ios::out);
	if(!outfile)
	{
		cout<<"can,t open.\n";
		abort();
	}
	outfile<<vertexNum<<"  "<<arcNum<<endl;//outfile<<"vertexNum="<<vertexNum<<"  "<<"arcNum="<<arcNum<<endl;              
	struct ArcNode *p;

	for(int i=0; i<vertexNum; i++)
	{
	//	if(adjlist[i].vertex != 0)
		{                  
			outfile<<adjlist[i].vertex<<" ";
			p=adjlist[i].firstedge;
			while(p!=NULL)
			{
				outfile<<p->adjvex<<" ";
				p=p->next;
			}
			outfile<<"0"<<endl;//add 0 to the end
		}
		
	}
	outfile<<endl;
	
	outfile.close();
}

template <class T>
void ALGraph<T>::Output_to_matlab(int number)//this function is simllar to Output(char * name)//name is file name
{
	ofstream outfile;
	char filename[50]="graph_";
	char name[20];
	sprintf(name, "%d", number);
	strcat(name, ".txt");
	strcat(filename, name);

	outfile.open(filename);
	//outfile.open(name, ios::out);
	if(!outfile)
	{
		cout<<"can,t open.\n";
		abort();
	}


	outfile<<vertexNum<<" ";//<<arcNum<<endl;//outfile<<"vertexNum="<<vertexNum<<"  "<<"arcNum="<<arcNum<<endl;              
	struct ArcNode *p;

	for(int i=0; i<vertexNum; i++)
	{
		if(adjlist[i].vertex != 0)
		{                  
		//	outfile<<adjlist[i].vertex<<" ";
			outfile<<adjlist[i].num<<" ";
			p=adjlist[i].firstedge;
			while(p!=NULL)
			{
				outfile<<p->adjvex<<" ";
				p=p->next;
			}
		//	outfile<<"0"<<endl;//add 0 to the end
		}
		
	}
//	outfile<<endl;
	
	outfile.close();
}

template <class T>
int * ALGraph<T>::Output_to_matlab( )//for mex function
{									//this function is simllar to Output(char * name)//name is file name
	struct ArcNode *p;
	int num = vertexNum;
	for(int i=0; i<vertexNum; i++)         
		num = num + adjlist[i].num;
	int * adj = new int [num +1];	
	int k = 0;
	for(i=0; i<vertexNum; i++){
		adj[k] = adjlist[i].num;
		k++;
		p=adjlist[i].firstedge;
		while(p!=NULL){
			adj[k] = p->adjvex;
			k++;
			p=p->next;
		}
	}
	return adj;
}

template <class T>
int * ALGraph<T>::GraphToArray( )//for mex function
{									//this function is the same as Output(char * name)
	struct ArcNode *p;
	int num = 3 + 2*vertexNum;
	int i;
	for(i=0; i<vertexNum; i++)         
		num = num + adjlist[i].num;
	int * adj = new int [num];	
	adj[0] = num;
	adj[1] = vertexNum;
	adj[2] = arcNum;
	int k = 3;
	for(i=0; i<vertexNum; i++){
		adj[k] = i+1;
		k++;
		p=adjlist[i].firstedge;
		while(p!=NULL){
			adj[k] = p->adjvex;
			k++;
			p=p->next;
		}
		adj[k] = 0;
		k++;
	}
	return adj;
}

template <class T>
void ALGraph<T>::OutputVertex()//output vertices
{

	fstream outfile;      
	outfile.open("vertex.txt", ios::app);
	if(!outfile)
	{
		cout<<"can,t open.\n";
		abort();
	}	           
	for(int i=0; i<vertexNum; i++)
	{
		if(adjlist[i].vertex != 0)
		{                  
			outfile<<adjlist[i].vertex<<" ";
		}
		
	}
	outfile<<endl;
	outfile<<"*************************"<<endl<<endl;
	outfile.close();
}





template<class T>
void ALGraph<T>::Compchain(int i)//Depth first search
{
	comp[i-1]=1;
	ArcNode * p=adjlist[i-1].firstedge;
	int j;
	while(p)
	{
		j=p->adjvex;
		if(comp[j-1]==-1)Compchain(j);
		p=p->next;

	}

}


template<class T>
int ALGraph<T>::Chain(int m, int n)
{
	ArcNode * p0;
	p0=adjlist[m-1].firstedge;
	while(p0 !=NULL)
	{
		if(p0->adjvex==n) return 0;//is not the chain
		p0=p0->next;
	}
	return 1; //not in the chain
}



/**************************************************************************************************************************/
//Begin MCS-M algorithm for computing minimal elimination ordering and minimal triangulation


template<class T>
int complete(ArcNode * first, ALGraph<T> & g)//called by 
{
	ArcNode * p, *p1;
	p=first;
	while(p != NULL)
	{
		p1=p->next;
		while(p1 !=NULL)
		{
			if( g.Chain(p->adjvex, p1->adjvex)) return 0;//not complete
			p1=p1->next;
		}
		p=p->next;
	}
	return 1;//complete
}


/**************************************************************************************************************************/
//Begin to construct an n_cycle

template<class T>                
ALGraph<T>::ALGraph(int num, char *cycle) //a tree////construct an n_cycle 
{
	struct ArcNode * s, * p, * p1;
	int   e=0;// the number of edges
	int a[M], b[M];
	if(num>M)
		cout<<"Please input a larger M\n"<<endl;
	vertexNum=num;  
//	cout<<"vertexNum="<<vertexNum<<endl;
	int k;
	for(k=2; k<num; k++)
	{
		adjlist[k-1].num=2;
		adjlist[k-1].vertex=k;
		p=new ArcNode;
		p->adjvex=k-1;
		adjlist[k-1].firstedge=p;
		p1=new ArcNode;
		p1->adjvex=k+1;
		p1->next=NULL;
		p->next=p1;		
	}
	
	k=1;
	adjlist[k-1].num=2;
		adjlist[k-1].vertex=k;
		p=new ArcNode;
		p->adjvex=k+1;
		adjlist[k-1].firstedge=p;
		p1=new ArcNode;
		p1->adjvex=num;
		p1->next=NULL;
		p->next=p1;

	k=num;
	adjlist[k-1].num=2;
		adjlist[k-1].vertex=k;
		p=new ArcNode;
		p->adjvex=1;
		adjlist[k-1].firstedge=p;
		p1=new ArcNode;
		p1->adjvex=num-1;
		p1->next=NULL;
		p->next=p1;

	arcNum=num;	

}


//construct a connect graph, num is the number of vertex，
double myRand()//Uniform(0, 1)
{
	//	return   (double)((double)rand() / ((double)RAND_MAX   +   1.0)   *   1); 
	return   (double)rand()/(RAND_MAX+1.0); 
}


int random(int k)
{
	return rand() % k; //draw a random number from 0 to k-1, not k
}



template<class T>                
ALGraph<T>::ALGraph(int num) //a tree////construct a tree 
{
	struct ArcNode * s, * p, * p1;
	int   e=0;// the number of edges
	int a[M], b[M];
	if(num>M)
		cout<<"Please input a larger M"<<endl;
	vertexNum=num;  
//	cout<<"vertexNum="<<vertexNum<<endl;
	int k;
	for(k=1; k<=num; k++)
	{
		adjlist[k-1].vertex=k;
		adjlist[k-1].num=0;
		adjlist[k-1].firstedge=NULL;
		a[k-1]=k;
	}
	k--;
//	cout<<endl;
	int k1, j2, j1=random(k);


		b[0]=a[j1];
		a[j1]=a[k-1];

	int num2;
	num2=num-1;
	for(k=num2; k>0; k--)
	{
		k1=num-k;
		j1=random(k);
		j2=random(k1);
		e++;
		
		p=adjlist[a[j1]-1].firstedge;
		s= new ArcNode;
		s->adjvex=b[j2];
		s->next=NULL;
		adjlist[a[j1]-1].num++;
		if(p==NULL)
			adjlist[a[j1]-1].firstedge=s;
		else
		{
			p1=p->next;
			while(p1)
			{
				p=p1;
				p1=p->next;
			}
			p->next=s;
		}
		
		p=adjlist[b[j2]-1].firstedge;
		s= new ArcNode;
		s->adjvex=a[j1];
		s->next=NULL;
		adjlist[b[j2]-1].num++;
		if(p==NULL)
			adjlist[b[j2]-1].firstedge=s;
		else
		{
			p1=p->next;
			while(p1)
			{
				p=p1;
				p1=p->next;
			}
			p->next=s;
		}
		
		
		b[k1]=a[j1];
		a[j1]=a[k-1];
		
	}
	for(k=0; k<num; k++)
	{
	}
	arcNum=e;	

}


template<class T>                
int ALGraph<T>::In_edge(int k, int i)//whether or not k is in adjacent set of i, called by ALGraph(int num, double prob)
{
	ArcNode * p, * p1=adjlist[i-1].firstedge;
	while(p1)
	{
		if(p1->adjvex==k)
			return 1;//k is in adjacent set of i
		else
		{
			p=p1;
			p1=p->next;
		}
	}
	return 0;//k is not in adjacent set of i

}



template<class T>                
ALGraph<T>::ALGraph(int num, double prob) //based on the tree, construct a connected graph, where num is the number of vertices, prob is the probality
{
	struct ArcNode * s, * p, * p1;
	ALGraph<T> g(num);// costruct a tree named tree
	//copy a tree
	vertexNum=num; 
	int	e=g.arcNum;//copy a tree
	for(int i=0; i< vertexNum; i++)
	{
		adjlist[i].vertex=g.adjlist[i].vertex; 
		adjlist[i].num=g.adjlist[i].num;
		ArcNode *bs, *bs1, *as;
		as=g.adjlist[i].firstedge;
		if(as!=NULL)
		{
			bs=new ArcNode;
			bs->adjvex=as->adjvex;
			adjlist[i].firstedge=bs;
			bs1=bs;
			as=as->next;		
			while(as!=NULL)
			{
				bs=new ArcNode;
				bs->adjvex=as->adjvex;
				bs1->next=bs;
				bs1=bs;
				as=as->next;

			}
			bs->next=NULL;

		}
		else
		{
			adjlist[i].firstedge=NULL;
		}
	}	
//_______________________________________________________________________________________________-/
//based on the tree, construct a connected graph, where num is the number of vertices, prob is the probality
	double r;
	for(int k=1; k<num; k++)
		for(int i=k+1; i<=num; i++)
		{
			r=myRand();//cout<<k<<", "<<i<<"r="<<r<<endl;
			if(r<=prob)
			{
				if(g.In_edge(k, i)==0)
				{
	//___________________________________________________________________________/
					e++;
					//add i to k's adjacent set
					p=adjlist[k-1].firstedge;
					s= new ArcNode;
					s->adjvex=i;
					s->next=NULL;
					adjlist[k-1].num++;
					if(p==NULL)
						adjlist[k-1].firstedge=s;
					else
					{
						p1=p->next;
						while(p1)//if p1==NULL, then p points to the last neighbor
						{
							p=p1;
							p1=p->next;
						}
						p->next=s;
					}
					//add k to i's adjacent set
					p=adjlist[i-1].firstedge;
					s= new ArcNode;
					s->adjvex=k;
					s->next=NULL;
					adjlist[i-1].num++;
				 	if(p==NULL)
						adjlist[i-1].firstedge=s;
					else
					{
						p1=p->next;
						while(p1)//if p1==NULL, then p points to the last neighbor
						{
							p=p1;
							p1=p->next;
						}
						p->next=s;
					}
	//_________________________________________________________________/

				}//End if(g.In_edge(k, i)==0)//
			}
		}
	

//_________________________________________________________/
	arcNum=e;	//  cout<<"arcNum="<<arcNum<<endl; cout<<"vertexNum="<<vertexNum<<endl;

}



//End costruct a connented graph
/**************************************************************************************************************************/



/**************************************************************************************************************************/
//copy the orginal paper of Lex_M and compile the procedure

 


void add_reach(int w, int i)// add w to reach[i-1];
{
	ArcNode * p, * p1;
	p=reach[i-1].first;
	p1=new ArcNode; 
	p1->adjvex=w;
	p1->next=p;
	reach[i-1].first=p1;
}

void add_madj(int v_max, int w)// add v_max into w's monotony adjcent neighbors
{
	ArcNode * p1, *p2;
	p1=new ArcNode; 
	p1->adjvex=v_max;
	p1->next=NULL;
	if(madj[w-1].first==NULL)
	{		
		madj[w-1].first=p1;
		madj[w-1].tail=p1;
	}
	else
	{
		p2=madj[w-1].tail;
		p2->next=p1;
		madj[w-1].tail=p1;
	}
	

}

void delete_reach( int j)//delete reach[j-1] the first vertex
{
	ArcNode * p;
	p=reach[j-1].first;
//	p1=p->next;
	reach[j-1].first=p->next;
	delete p; 

}



template<class T>
void ALGraph<T>::Lex_M(ALGraph<T> & g, int largest)//largest is the largest numbered vertex
{
	ArcNode *p, *p1;

	float la[M];//save label numbers
	int reached[M];//	ALGraph<T> g(g0);
	int n=g.Num_of_vertex();
	int i;
	for(i=0; i<n; i++)
	{
		labeled[i]=0;
		la[i]=1;
		reach[i].first=NULL;
		reached[i]=0;
		madj[i].first=NULL;
	}
	int k=1;
	int v_max=largest;
	Lex_order[n-1]=v_max;
	labeled[v_max-1]=n;//v_max is n_th vertex
	reached[v_max-1]=1;
	for(i=n; i>1; i--)
	{
		//select://mark all unnumbered vertices unreached
		
		for(int j=0; j<n; j++)
		{
			if(labeled[j]==0)
				reached[j]=0;
		}//updage v_max's neighbor la
		p=g.adjlist[v_max-1].firstedge;
		int w, nei;
		while(p)
		{
			w=p->adjvex;
			if(labeled[w-1]==0)
			{
				add_reach(w, int(la[w-1]));
				reached[w-1]=1;
				la[w-1]=la[w-1]+float(0.5);
				add_madj(v_max, w);
			}
			p=p->next;	
		}
		//search: fill-path 
		{for(int j=1; j<=k; j++)
		{
			p=reach[j-1].first;
			while(p)
			{
				//delete w from reach[j-1]
				w=p->adjvex;
				delete_reach(j);//delete reach[j-1]'s first vertex
				
				p1=g.adjlist[w-1].firstedge;
				while(p1)
				{
					nei=p1->adjvex;
					if(reached[nei-1]==0)
					{
						reached[nei-1]=1;
						if(la[nei-1] > j)
						{
							add_reach(nei, int(la[nei-1]));
							la[nei-1]=la[nei-1]+float(0.5);
							add_madj(v_max, nei);
						}
						else
							add_reach(nei, j);
					}
					p1=p1->next;
				}
				p=reach[j-1].first;
			}
		}}
		//sort: sort unnumbered vertices by la(w) value
		k=0;
		{for(int j=0; j<n; j++)
		{
			if(labeled[j]==0)
			{
				la[j]=float(ceil(la[j]));
				if(la[j] > k)
				{
					k=int(la[j]);
					v_max=j+1;
				}
			}
		}}
		Lex_order[i-2]=v_max;
		labeled[v_max-1]=i-1;
		reached[v_max-1]=1;
	}
	int mm=0;
	for(i=0; i<n; i++)
	{
		p=madj[i].first;
		cout<<"madj["<<i+1<<"]=";
		while(p)
		{
			p1=p;
			cout<<p->adjvex<<" ";
			p=p->next;
			
			mm++;
		}
		cout<<endl<<i+1<<" mm="<<mm<<endl;

	}
	cout<<"mm="<<mm<<endl;


}


int compare1(ArcNode *p1, ArcNode *p2) //return 1 时， p1 >= p2 ; return 0 时 p1< p2

{
	do
	{
		if(p2==NULL)return 1;
		else if (p1==NULL) return 0;
		else if (labeled[p1->adjvex -1] > labeled[p2->adjvex -1])	return 1;
		else if (labeled[p1->adjvex -1] < labeled[p1->adjvex -1])   return 0;
		else
		{
			p1=p1->next; p2=p2->next;
		}
	}while(1);
}


ArcNode * F1(int verticesNum)//F saves the numbers, not vertices
{
	int ff=0;
	ArcNode * p, * p0,  * head;
	int n=0;	
	for(int i=0; i<verticesNum-1; i++)
	{ 
		if(compare1(madj[Lex_order[i+1]-1].first,  madj[Lex_order[i]-1].first))
		{
			if(n==0)
			{
				p=new ArcNode; n_d++;
				ff++;
				p->adjvex=i+1;//is number, not vertex
				p->next=NULL;
				p0=p;
				n++;
			}
			else 
			{
				p=new ArcNode; n_d++;
				ff++;
				p->adjvex=i+1;   
				p->next=NULL;
				p0->next=p;
				if(n==1) 
				{
					head=p0; n++;
				}
				p0=p;
			}
		}
	}
	p=new ArcNode; n_d++;
	ff++;
	p->adjvex=verticesNum;//is number, not vertex
	p->next=NULL;
	if(n==0)//for the case: complete graph
		head=p;
	else if(n==1)
	{
		p0->next=p;
		head=p0;		
	}
	else
		p0->next=p;
	cout<<"ff="<<ff<<endl;

	return head;
}


ArcNode * E1(ArcNode *head, int verticesNum)
{
	int ee=0;
	ArcNode * p, * e, * e1, *h;
	int n=0;
	e=new ArcNode; n_d++;
	ee++;
	//e->adjvex=order[0];
	e->adjvex=Lex_order[0];
	e->next=NULL;
	e1=e;
	p=head;
	while(p->adjvex != verticesNum)
	{
		e=new ArcNode; n_d++;
		ee++;
	//	e->adjvex=order[p->adjvex];
		e->adjvex=Lex_order[p->adjvex];
		e->next=NULL;
		e1->next=e;
		if(n==0)
		{
			h=e1; n++;
		}
		e1=e;
		p=p->next;
	}
	if(n==0) h=e1;//for the case: complete graph
	cout<<"ee="<<ee<<endl;
	return h;

}

int relabel_1(ArcNode * e, ArcNode * f)
{
	ArcNode * p1, * p2, * p, *p0;
	int m, n, k, num=0;
	p1=e; p2=f;
	while(p1 != NULL && p2 != NULL)
	{
		num++;
		m=p1->adjvex;  n=p2->adjvex;//m is vertex, n is a number
	//	p=L[m-1].first;
		p=madj[m-1].first;
		if(p !=NULL)
		{
		//	k=p->adjvex;//where p->adjvex is a number, not vertex
			k=labeled[p->adjvex -1];
			p0=p;
			while(k > n)
			{
				p0=p;
				p=p->next;
				if(p==NULL) break;
			//	k=p->adjvex;
				k=labeled[p->adjvex -1];
			}
			if(p !=NULL) p0->next=NULL; 
		}
		p1=p1->next;
		p2=p2->next;
	}
	return num;
	
}



/**************************************************************************************************************************/


//Mcs_M copy the paper of lex_M
template<class T>
void ALGraph<T>::Mcs_M2(ALGraph<T> & g, int largest)//largest is largest numbered vertex  
{
	ArcNode *p, *p1;
	int reached[M]; 

	int n=g.Num_of_vertex();
	int i;
	for(i=0; i<n; i++)
	{
		labeled[i]=0;
	//	la[i]=1;
		weight[i]=0;
		reach[i].first=NULL;
		reached[i]=0;
		madj[i].first=NULL;
	}
	
	int k=1;
	int v_max=largest;
	McsM_order[n-1]=v_max;

	labeled[v_max-1]=n;//v_max is the n_th vertex；
	reached[v_max-1]=1;
	for(i=n; i>1; i--)
	{
		//select:
		//mark all unnumbered vertices unreached
	
	
		{for(int j=0; j<n; j++)
		{
	
			if(labeled[j]==0)
				reached[j]=0;
		}}
		//update v_max's neighbor la
		p=g.adjlist[v_max-1].firstedge;
		int w, nei;
		while(p)
		{
			w=p->adjvex;
			if(labeled[w-1]==0)
			{
			
				add_reach(w, weight[w-1]+1);
				reached[w-1]=1;
				weight[w-1]++;			
				add_madj(v_max, w);
			}
			p=p->next;	
		}
		//search:fill-in path
	
		for(int j=1; j<=k; j++)
		{
			p=reach[j-1].first;
			while(p)
			{
				//delete w from reach[j-1]
				w=p->adjvex;
				delete_reach(j);//delete reach[j-1]'s first edge
				
				p1=g.adjlist[w-1].firstedge;
				while(p1)
				{
					nei=p1->adjvex;
					if(reached[nei-1]==0)
					{
						reached[nei-1]=1;
					//	if(la[nei-1] > j)
						if(weight[nei-1]+1 > j)
						{
						
							add_reach(nei, weight[nei-1]+1);
						
							weight[nei-1]++;			
							add_madj(v_max, nei);
						}
						else
							add_reach(nei, j);
					}
					p1=p1->next;
				}
				p=reach[j-1].first;
			}
		}
	
		//sort: sort unnumbered vertices by la(w) value
		k=0;
		for(j=0; j<n; j++)
		{
			if(labeled[j]==0)
			{
			//	la[j]=float(ceil(la[j]));
			//	if(la[j] > k)
				if(weight[j] > k)
				{
				//	k=int(la[j]);
					k=weight[j];
					v_max=j+1;
				}
			}
		}
		McsM_order[i-2]=v_max;
		labeled[v_max-1]=i-1;
		reached[v_max-1]=1;

	}


}


/**************************************************************************************************************************/



ArcNode *F_t(int n)//this F_t doesn't contain the largest numbered vertex
{
	ArcNode *s, *s1, *head;
	int h=0;//head's flag
	for(int i=0; i<n-1; i++ )
	{
		if(weight[/*v_1*/McsM_order[i]-1]<=weight[/*v_2*/McsM_order[i+1]-1])
		{
			if(h==0)
			{
				h++;
				s=new ArcNode;n_d++;
				s->adjvex=/*v_1*/McsM_order[i];
			//	s->next=NULL;
				s1=s;
				head=s;

			}
			else
			{
				s=new ArcNode; n_d++;
				s->adjvex=/*v_1*/McsM_order[i];
		//		s->next=NULL;
				s1->next=s;
				s1=s;
			}
	 		
		}

	}//	s1->next=McsM_order[n-1];
	if(h==0)
		head=NULL;
	else
		s1->next=NULL;

	return head;

}

template<class T>
void ALGraph<T>::Print_block(int number)// called by Decom_McsM_4(ALGraph<T> & g), print prime blocks to the screen, i.e., only print j such that comp[j]=0 or comp[j]=1
{
	ArcNode * p;
	int j;
	cout<<"  g "<<number<<" ="<<endl;
	for(int i=0; i<vertexNum; i++)
	{
		if(comp[i]==1 || comp[i]==0)
		{
			cout<<i+1<<" ";
			p=adjlist[i].firstedge;
			while(p)
			{
				j=p->adjvex;
				if(comp[j-1]==1 || comp[j-1]==0)
					cout<<j<<" ";
				p=p->next;
			}
			cout<<endl;
		}
	}
	cout<<endl;
}

/**************************************************************************************************************************/


//Begin Decom_McsM_4 version 


template<class T>
void ALGraph<T>::Decom_McsM_4(ALGraph<T> & g0,   ArcNode * f_t, int n)//called by Decom_McsM_3(ALGraph<T> & g), don't copy g1，g2
{
	ArcNode *pt, * p;//*pt point to f_t
	int number=1;//number of cliques


	pt=f_t;
	while(pt)
	{
		mfnum++;
	//	v_1=pt->adjvex;
		p=madj[/*v_1*/pt->adjvex -1].first;
		mc1=clock();
		int c=complete(p, g0);
		mc2=clock();
		mc=mc+double(mc2-mc1)/CLOCKS_PER_SEC;
		if(/*complete(p, g0)*/c)
			{
				
				mg1=clock();
				mg2=clock();
				mg=mg+double(mg2-mg1)/CLOCKS_PER_SEC;
				
				md1=clock();
				while(p !=NULL)
				{
					comp[p->adjvex-1]=0;
					p=p->next;
				}
				md2=clock();
				md=md+double(md2-md1)/CLOCKS_PER_SEC;

				mco1=clock();//	g2.Compchain(/*v_1*/pt->adjvex);
				g0.Compchain(/*v_1*/pt->adjvex);
				mco2=clock();
				mco=mco+double(mco2-mco1)/CLOCKS_PER_SEC;


								
				g0.Print_block(number);//only print j such that comp[j]=0 or comp[j]=1
				

				mdg1=clock();
				for(int j=0; j<n; j++)
				{	                       //cout<<"comp["<<j<<"]="<<comp[j-1]<<" ";
					if(comp[j]==0)	comp[j]=-1;					    	
					else if (comp[j]==1)  comp[j]++;//comp[j]=2, the vertex j-1 is deleted 
										
				}
				mdg2=clock();
				mdg=mdg+double(mdg2-mdg1)/CLOCKS_PER_SEC;
				number++;
			}		

		pt=pt->next;
	}
	mfnum++;

	for(int i=0; i<n; i++)
		if(comp[i]==-1)
			comp[i]=1;
	 
	g0.Print_block(number);cout<<endl;

	cout<<"the number of the prime blocks is	"<<number<<"	by Decom_McsM_4 "<<endl;
	//realse the memory of madj[];
	ArcNode *p1;
	for( i=0; i<n; i++)
	{
		p=madj[i].first;
		while(p)
		{
			p1=p;
			p=p->next;
			delete p1; n_d--;
		}

	}
	p=f_t;//realse the memory of f_t
	while(p)
	{
		p1=p;
		p=p->next;
		delete p1;  n_d--;
	}//realse the memory 
}

template<class T>
void Decom_McsM_4(ALGraph<T> & g)
										
{
	int n=g.Num_of_vertex();
	for(int i=0; i<n; i++)
		comp[i]=-1;
			
	g.Mcs_M2(g, 1);
	//_________________________
	t_Mcs=clock();
	//Mcs order time 
	
	mf1=clock();
	ArcNode *f_t=F_t(n);
	mf2=clock();
	ArcNode * p=f_t;


	g.Decom_McsM_4(g, f_t, n );
}





//End Decom_McsM_4 version
/**************************************************************************************************************************/

// Begin Lex_M_4 version

template<class T>
void ALGraph<T>::Decom_LexM_4(ALGraph<T> & g0, ArcNode * e, ArcNode * f, int num)
{		
	ArcNode * p1, * p2, * p;
	p1=e; p2=f; int number=0;
	int i;
	for( i=1; i< num; i++ )
//	while(p1)
	{
		p=madj[p1->adjvex -1].first;
		lc1=clock();
		int c=complete(p, g0);
		lc2=clock();
		lc=lc+double(lc2-lc1)/CLOCKS_PER_SEC;

		if(/*complete(p, g0)*/c)
		{
			lg1=clock();
		//	ALGraph<T> g1(g0), g2(g0);
			lg2=clock();
			lg=lg+double(lg2-lg1)/CLOCKS_PER_SEC;

			number++;

			ld1=clock();
			while(p !=NULL)
			{
				comp[p->adjvex-1]=0;
				p=p->next;
			}
			ld2=clock();
			ld=ld+double(ld2-ld1)/CLOCKS_PER_SEC;

			lco1=clock();
			g0.Compchain(Lex_order[p2->adjvex-1]);
			lco2=clock();
			lco=lco+double(lco2-lco1)/CLOCKS_PER_SEC;

			//print prime blocks
							
			g0.Print_block(number);
			


			ldg1=clock();

			for(int j=0; j<vertexNum; j++)
			{	                      
				if(comp[j]==0)	comp[j]=-1;					    	
				else if (comp[j]==1)  comp[j]++;
									
			}
			ldg2=clock();
			ldg=ldg+double(ldg2-ldg1)/CLOCKS_PER_SEC;
	

		}
		p1=p1->next;  p2=p2->next;
	}
	number++;
	
	num_blocks=number;//global variable, return the number of prime blocks

	for(i=0; i<vertexNum; i++)
		if(comp[i]==-1)
			comp[i]=1;
	 
	g0.Print_block(number);


	cout<<"the number of the mp-subgraphs is	"<<number<<"	by Decom_LexM_4 version"<<endl;
}



ArcNode * E2(ArcNode *head, int verticesNum)
{
	int ee=0;
	ArcNode * p, * e, * e1, *h=NULL;
	int n=0;
	e=new ArcNode; n_d++;
	ee++;
	//e->adjvex=order[0];
	e->adjvex=Lex_order[0];
	e->next=NULL;
	e1=e;
	p=head;
	while(p->adjvex != verticesNum)
	{
		e=new ArcNode; n_d++;
		ee++;
	//	e->adjvex=order[p->adjvex];
		e->adjvex=Lex_order[p->adjvex];
		e->next=NULL;
		e1->next=e;
		if(n==0)
		{
			h=e1; n++;
		}
		e1=e;
		p=p->next;
	}
	if(n==0) h=e1;//for the case: complete graph
	cout<<"ee="<<ee<<endl;
	return h;

}


template<class T>
void Decom_LexM_4(/*ALGraph<T> & g2,*/ ALGraph<T> & g)

{
	ArcNode *p, *p1;
	cout<<"Begin LexM, n_d="<<n_d<<endl;
//	ALGraph<T> g2(g), g3(g);//copy g2, g3 for ordering and decompostion's delete edges
	int n=g.Num_of_vertex();
	int i=0;
	for( i=0; i<n; i++)
		comp[i]=-1;

//	g2.Lex_M(g2, 1);
	g.Lex_M(g, 1);
	for(i=0; i<n; i++)
	{
		p=madj[i].first;
		cout<<"madj["<<i+1<<"]=";
		while(p)
		{
			p1=p;
			cout<<p->adjvex<<" ";
			p=p->next;
		}
		cout<<endl;
	}

	//Lex ordering time

	t_lex=clock();
	
	//Lex ordering time
	
	lf1=clock();
	cout<<"n_d="<<n_d<<endl;
	ArcNode *f=F1(n);
	cout<<"n_d="<<n_d<<endl;
	for(i=0; i<n; i++)
	{
		p=madj[i].first;
		cout<<"madj["<<i+1<<"]=";
		while(p)
		{
			p1=p;
			cout<<p->adjvex<<" ";
			p=p->next;
		}
		cout<<endl;
	}
	
	ArcNode *e=E1(f, n);
	cout<<"n_d="<<n_d<<endl;
	for(i=0; i<n; i++)
	{
		p=madj[i].first;
		cout<<"madj["<<i+1<<"]=";
		while(p)
		{
			p1=p;
			cout<<p->adjvex<<" ";
			p=p->next;
		}
		cout<<endl;
	}

	int num=relabel_1(e, f); //number of cliques of fill-in graph

	lf2=clock();
	lfnum=num;//number of cliques of fill-in graph
	cout<<"relable 之后\n";
	for(i=0; i<n; i++)
	{
		p=madj[i].first;
		cout<<"madj["<<i+1<<"]=";
		while(p)
		{
			p1=p;
			cout<<p->adjvex<<" ";
			p=p->next;
		}
		cout<<endl;
	}
	
//	g3.Decom_LexM_4(g3, e, f, num); 
	cout<<"n_d="<<n_d<<endl;
    g.Decom_LexM_4(g, e, f, num); 
	cout<<"n_d="<<n_d<<endl;
	for(i=0; i<n; i++)
	{
		p=madj[i].first;
		cout<<"madj["<<i+1<<"]=";
		while(p)
		{
			p1=p;
			cout<<p->adjvex<<" ";
			p=p->next;
		}
		cout<<endl;
	}

	//releas memory 
	
	//releas memory of f
	p=f;

	int ff=0;
	while(p)
	{
		p1=p;
		p=p->next;
		delete p1; n_d--;
		ff++;
	}
	cout<<"ff="<<ff<<endl;
	cout<<"n_d="<<n_d<<endl;
	//releas memory of e
	p=e;
	int ee=0;
	while(p)
	{
		p1=p;
		p=p->next;
		delete p1; n_d--;
		ee++;
	}
	cout<<"ee="<<ee<<endl;
	cout<<"n_d="<<n_d<<endl;
	//releas memory of madj[];
	int mm=0;
	for(i=0; i<n; i++)
	{
		p=madj[i].first;
		cout<<"madj["<<i+1<<"]=";
		while(p)
		{
			p1=p;
			cout<<p->adjvex<<" ";
			p=p->next;
			delete p1; n_d--;
			mm++;
		}
		cout<<endl<<i+1<<"  mm="<<mm<<endl;

	}
	cout<<"mm="<<mm<<endl;
	cout<<"End LexM, n_d="<<n_d<<endl;
//	cout<<"releas memory..."<<endl;
	
}

// End Lex_M_4 version 
/*************************************************************************/
//Begin  construct a junction tree of fill_in graph of prime block


void junction_tree(int n, int number)//n is the number of vertices of glabel graph
{						//print the junction tree of fill_in graph of the number_th prime block, 
						//the first number, i.e., 0, is the parent  of the first node, 
						//then the number of vertices of the first node, and then vertices...
						//Then begin print the second clique's parent, its verticies number, and its vertices
						//the following is the third, forth, ...
						//And the last number is the number of nodes of the junction tree
	int j, v1, v2, num_clique=0, label=n;//v1, v2 are two vertices, and label=n is used to number the vertices in the first clique 
	int clique[M];// is the funciton defined in orginal paper
	ArcNode *p;
	for(int i=n-1; i>=0; i--)
	{
		if(comp[McsM_order[i]-1]==1)
		{
			v2=McsM_order[i];
			j=i-1;
			i=-1;
		}
	}
/*	cout<<"\n ___________________________________________"<<endl;	
	for(i=0; i<n; i++)
	{
		cout<<i+1<<"  dan diao lin ji：";
		p=madj[i].first;
		if(p)
		{
		while(p)
		{
			cout<<p->adjvex<<", ";cout<<"p="<<p<<" ";
			p=p->next;
		}
		cout<<endl;
		p=madj[i].tail;
		cout<<" p="<<p;
		cout<<" "<<p->adjvex<<", ";
		}cout<<endl;
	}

cout<<"\n ___________________________________________"<<endl;//*/


	ofstream outfile; 
	char filename[50]="junction_tree_";
	char name[20];
	sprintf(name, "%d", number);
	strcat(name, ".txt");
	strcat(filename, name);

	outfile.open(filename);
	//outfile.open(name, ios::out);
	if(!outfile)
	{
		cout<<"can,t open.\n";
		abort();
	}
	outfile<<0<<" ";//print the number of parent clique of the first clique,
	int r_num;//r is the number of vertices in R=C\S, where C is the clique and S is the separator
	while(j>=0)
	{
		if(comp[McsM_order[j]-1]==1)
		{
			v1=McsM_order[j];
			if(weight[v1-1]<=weight[v2-1])
			{
		/*		num_clique++;
				//print the number of vertices of the clique, and all vertices
				outfile<<weight[v2-1]+1<<" "<<v2<<" ";
				clique[v2-1]=num_clique;
				p=madj[v2-1].first;
				while(p)
				{
					outfile<<p->adjvex<<" ";
					if(labeled[p->adjvex -1]<=label)
						clique[p->adjvex -1]=num_clique;
					p=p->next;
				}
				label=labeled[v1-1];
				//print the number of parent clique of next clique
				p=madj[v1-1].tail;
				outfile<<clique[p->adjvex -1]<<" ";//*/

				num_clique++;
				//print the number of vertices of the clique, and all vertices
				outfile<<weight[v2-1]+1<<" ";
				clique[v2-1]=num_clique;
				r_num=1;
				p=madj[v2-1].first;
				while(p)
				{
					outfile<<p->adjvex<<" ";
					if(labeled[p->adjvex -1]<=label)
					{
						clique[p->adjvex -1]=num_clique;
						r_num++;
					}
					p=p->next;
				}
				outfile<<v2<<" "<<r_num<<" ";
				label=labeled[v1-1];
				//print the number of parent clique of next clique
				p=madj[v1-1].tail;
				outfile<<clique[p->adjvex -1]<<" ";
			
			}
			v2=v1;
		}
		j--;
	}

/*	num_clique++;
	outfile<<weight[v2-1]+1<<" "<<v2<<" ";
	p=madj[v2-1].first;
	while(p)
	{
		outfile<<p->adjvex<<" ";
	//	if(labeled[p->adjvex -1]<=label)
	//		clique[p->adjvex -1]=num_clique;
		p=p->next;
	}
	outfile<<num_clique<<endl;//*/

	num_clique++;
	outfile<<weight[v2-1]+1<<" ";
	r_num=1;
	p=madj[v2-1].first;
	while(p)
	{
		outfile<<p->adjvex<<" ";
		if(labeled[p->adjvex -1]<=label)
		{
			clique[p->adjvex -1]=num_clique;
			r_num++;
		}
		p=p->next;
	}
	outfile<<v2<<" "<<r_num<<" "<<num_clique<<endl;

	
	outfile.close();//*/

	
}





//End  construct a junction tree 
/**************************************************************************************************************************/
/*************************************************************************/
//Begin  construct a junction tree of the global fill_in graph 

template<class T>
void  ALGraph<T>::jun_tree()//n is the number of vertices of glabel graph
						//print the junction tree of fill_in graph, 
						//the first number, i.e., 0, is the parent  of the first node, 
						//then the number of vertices of the first node, and then vertices...
						//Then begin print the second clique's parent, its verticies number, and its vertices
						//the following is the third, forth, ...
						//And the last number is the number of nodes of the junction tree
{
	Mcs_M5(1);//firstly we order the vertices, then construct junction tree
	int n=vertexNum;
	int j, v1, v2, num_clique=0, label=n;//v1, v2 are two vertices, and label=n is used to number the vertices in the first clique 
	int clique[M];// is the funciton defined in orginal paper
	ArcNode *p;
	
	ofstream outfile; 
	char filename[50]="jun_tree_";
	char name[30];
	sprintf(name, "%d", n);
	strcat(name, ".txt");
	strcat(filename, name);

	outfile.open(filename);
	//outfile.open(name, ios::out);
	if(!outfile)
	{
		cout<<"can,t open.\n";
		abort();
	}
	outfile<<0<<" ";//print the number of parent clique of the first clique,

	v2=McsM_order[n-1];
	j=n-2;
	int r_num;//r is the number of vertices in R=C\S, where C is the clique and S is the separator
	while(j>=0)
	{
		{
			v1=McsM_order[j];
			if(weight[v1-1]<=weight[v2-1])
			{
				num_clique++;
				//print the number of vertices of the clique, and all vertices
				outfile<<weight[v2-1]+1<<" ";
				clique[v2-1]=num_clique;
				r_num=1;
				p=madj[v2-1].first;
				while(p)
				{
					outfile<<p->adjvex<<" ";
					if(labeled[p->adjvex -1]<=label)
					{
						clique[p->adjvex -1]=num_clique;
						r_num++;
					}
					p=p->next;
				}
				outfile<<v2<<" "<<r_num<<" ";
				label=labeled[v1-1];
				//print the number of parent clique of next clique
				p=madj[v1-1].tail;
				outfile<<clique[p->adjvex -1]<<" ";
			
			}
			v2=v1;
		}
		j--;
	}
	num_clique++;
	outfile<<weight[v2-1]+1<<" ";
	r_num=1;
	p=madj[v2-1].first;
	while(p)
	{
		outfile<<p->adjvex<<" ";
		if(labeled[p->adjvex -1]<=label)
		{
			clique[p->adjvex -1]=num_clique;
			r_num++;
		}
		p=p->next;
	}
	outfile<<v2<<" "<<r_num<<" "<<num_clique<<endl;


	
	outfile.close();//*/
	//realse the memory of madj[];
	
	for(int i=0; i<n; i++)
	{
		p=madj[i].first;
		while(p)
		{
			p1=p;
			p=p->next;			
	//		cout<<p1->adjvex<<"	";
			delete p1; n_d--;
		}
//		cout<<endl;

	}//*/

	
}


//End  construct a junction tree of global graph
/**************************************************************************************************************************/
/*************************************************************************/
//Begin to reduce the memory of Mcs_M algorithm, reduce the memory of madj[i];  //Mcs_M copy the paper of lex_M


template<class T>
ArcNode * ALGraph<T>::Mcs_M5(int largest)//reduce the memory of madj[i];  //largest is largest numbered vertex  
{
	ArcNode *s, * head=NULL;
	int h=0;//head's flag

	ArcNode *p, *p1;//, *p3, *p4;
	int reached[M]; 

	int n=vertexNum;
	int i;
	for(i=0; i<n; i++)
	{
		labeled[i]=0;
	//	la[i]=1;
		weight[i]=0;
		reach[i].first=NULL;
		reached[i]=0;
		madj[i].first=NULL;
	}
	
	int k=1;
	int v_max=largest, v_max_old;
	McsM_order[n-1]=v_max;

	labeled[v_max-1]=n;//v_max is the n_th vertex；
	reached[v_max-1]=1;
	for(i=n; i>1; i--)
	{
		//select:
        int j;
        //for j:=1 until k do reach(j): = \emptyset
        for(j=1; j<=vertexNum; j++){
            p=reach[j-1].first;
            ArcNode *q;
            while (p){
                q = p;
                p = p->next;
                delete q;
            }
            reach[j-1].first=NULL;
        }//*/
		//mark all unnumbered vertices unreached
		for(j=0; j<n; j++)
		{
	
			if(labeled[j]==0)
				reached[j]=0;
		}
		//update v_max's neighbor la
		//p=g.adjlist[v_max-1].firstedge;
		p=adjlist[v_max-1].firstedge;
		int w, nei;
		while(p)
		{
			w=p->adjvex;
			if(labeled[w-1]==0)
			{
			
				add_reach(w, weight[w-1]+1);
				reached[w-1]=1;
				weight[w-1]++;			
				add_madj(v_max, w);
			}
			p=p->next;	
		}
		//search:fill-in path
	
		for(j=1; j<=k; j++)
		{
			p=reach[j-1].first;
			while(p)
			{
				//delete w from reach[j-1]
				w=p->adjvex;
				delete_reach(j);//delete reach[j-1]'s first edge
				
			//	p1=g.adjlist[w-1].firstedge;
				p1=adjlist[w-1].firstedge;
				while(p1)
				{
					nei=p1->adjvex;
					if(reached[nei-1]==0)
					{
						reached[nei-1]=1;
					//	if(la[nei-1] > j)
						if(weight[nei-1]+1 > j)
						{
						
							add_reach(nei, weight[nei-1]+1);
						
							weight[nei-1]++;			
							add_madj(v_max, nei);
						}
						else
						{
							add_reach(nei, j);
						}
					}
					p1=p1->next;
				}
				p=reach[j-1].first;
			}
		}
		//in order to determine f_t
		v_max_old=v_max;
		//sort: sort unnumbered vertices by la(w) value
		k=0;
		for(j=0; j<n; j++)
		{
			if(labeled[j]==0)
			{
				if(weight[j] > k)
				{
				//	k=int(la[j]);
					k=weight[j];
					v_max=j+1;
				}
			}
		}
				
		McsM_order[i-2]=v_max;
		labeled[v_max-1]=i-1;
		reached[v_max-1]=1;

		//determine f_t		
		if(weight[v_max-1] <= weight[v_max_old-1])//If the conditoion is satisfied, then v_max a f_t
		{
			if(h==0)
			{
				s=new ArcNode; n_d++;
				s->adjvex=v_max;
				s->next=NULL;
				head=s;
				h++;
			}
			else
			{
				s=new ArcNode; n_d++;
				s->adjvex=v_max;
				s->next=head;
				head=s;

			}	 		
		}
		//in order to construct junction tree of prime block, we save the madj
	/*	else//otherwise, v_max's madj is not useful, in order to reduce memory, we delete its madj 
		{
			p3=madj[v_max-1].first;
			while(p3)
			{
				p4=p3;
				p3=p3->next;
				delete p4; n_d--;
			}
			madj[v_max-1].first=NULL;
		}*/
		

	}

/*
	cout<<endl;
	for(i=n-1; i>=0; i--)
	{
		cout<<"McsM_order["<<i+1<<"]= "<<McsM_order[i]<<endl;
	}
	cout<<"___________________"<<endl;
	for(i=0; i<n; i++)
	{
		cout<<i+1<<"  dan diao lin ji：";
		p1=madj[i].first;
		if(p1)
		{
		while(p1)
		{
			cout<<p1->adjvex<<", ";cout<<"p1="<<p1<<" ";
			p1=p1->next;
		}
		cout<<endl;
		p1=madj[i].tail;
		cout<<"	"<<p1->adjvex<<", ";
		cout<<"p1="<<p1<<endl;
		}
	}
	
	cout<<"____________________"<<endl;
/*	for(i=0; i<n; i++)
	{
		cout<<"labeled["<<i+1<<"]= "<<labeled[i]<<endl;
	}
	for(i=0; i<n; i++)
	{
		cout<<"weight["<<i+1<<"]="<<weight[i];
		cout<<"	 "<<i+1<<"  labels为：";
		p1=madj[i].first;
		while(p1)
		{
			cout<<labeled[p1->adjvex-1]<<", ";
			p1=p1->next;
		}
		cout<<endl;
	}//*/
/*	
	p=head;
	while(p)
	{
		cout<<p->adjvex<<"	";
		p=p->next;
	}//*/
    int j;
        //for j:=1 until k do reach(j): = \emptyset
        for(j=1; j<=vertexNum; j++){
            p=reach[j-1].first;
            ArcNode *q;
            while (p){
                q = p;
                p = p->next;
                delete q;
            }
            reach[j-1].first=NULL;
        }//*/

	return head;



}





template<class T>
ArcNode * ALGraph<T>::Mcs_M5(double subset[], int num)//reduce the memory of madj[i];  //largest is largest numbered vertex  
{
	ArcNode *s, * head=NULL;
	int h=0;//head's flag

	ArcNode *p, *p1;//, *p3, *p4;
	int reached[M]; 

	int n=vertexNum;
	int i;
	for(i=0; i<n; i++)
	{
		labeled[i]=0;
	//	la[i]=1;
		weight[i]=0;
		reach[i].first=NULL;
		reached[i]=0;
		madj[i].first=NULL;
	}
/*	//==========outfile======
		ofstream outfile; 
	char filename[50]="Mcs_M5_";
	char name[30];
	sprintf(name, "%d", n);
	strcat(name, ".txt");
	strcat(filename, name);

	outfile.open(filename);
	//outfile.open(name, ios::out);
	if(!outfile)
	{
		cout<<"can,t open.\n";
		abort();
	}

	for(i=0; i<num; i++)
		outfile<<subset[i]<<"  ";
	outfile<<endl<<"---------------------\n";
	//============outfile===========*/

	int k=1;
	int v_max = int(subset[0]), v_max_old;
	McsM_order[n-1]=v_max;

	labeled[v_max-1]=n;//v_max is the n_th vertex；
	reached[v_max-1]=1;

	int index = 1;

	for(i=n; i>1; i--)
	{
		//select:
        int j;
        //for j:=1 until k do reach(j): = \emptyset
        for(j=1; j<=k; j++){
            p=reach[j-1].first;
            ArcNode *q;
            while (p){
                q = p;
                p = p->next;
                delete q;
            }
            reach[j-1].first=NULL;
        }//*/
		//mark all unnumbered vertices unreached
		for(j=0; j<n; j++)
		{	
			if(labeled[j]==0)
				reached[j]=0;
		}
		//update v_max's neighbor la
		//p=g.adjlist[v_max-1].firstedge;
		p=adjlist[v_max-1].firstedge;
		int w, nei;
		while(p)
		{
			w=p->adjvex;
			if(labeled[w-1]==0)
			{
			
				add_reach(w, weight[w-1]+1);
				reached[w-1]=1;
				weight[w-1]++;			
				add_madj(v_max, w);
			}
			p=p->next;	
		}
		//search:fill-in path
	
		for(j=1; j<=k; j++)
		{
			p=reach[j-1].first;
            /*    outfile<<"reachreach("<<j<<") = ";
                ArcNode *q;
                q = p;
                while (q){
                    outfile<<q->adjvex<<"  ";
                    q = q->next;
                }outfile<<endl;*/
			while(p)
			{
				//delete w from reach[j-1]
				w=p->adjvex;
                
				delete_reach(j);//delete reach[j-1]'s first edge
				
			//	p1=g.adjlist[w-1].firstedge;
				p1=adjlist[w-1].firstedge;
				while(p1)
				{
					nei=p1->adjvex;
					if(reached[nei-1]==0)
					{
						reached[nei-1]=1;
					//	if(la[nei-1] > j)
						if(weight[nei-1]+1 > j)
						{
						
							add_reach(nei, weight[nei-1]+1);
						
							weight[nei-1]++;			
							add_madj(v_max, nei);
						}
						else
						{
							add_reach(nei, j);
						}
					}
					p1=p1->next;
				}
				p=reach[j-1].first;
			}
		}
     /*   outfile<<endl;
        for(j=1; j<=7; j++)
		{
			p=reach[j-1].first;
                outfile<<"reach("<<j<<") = ";
                ArcNode *q;
                q = p;
                while (q){
                    outfile<<q->adjvex<<"  ";
                    q = q->next;
                }outfile<<endl;
        }
        outfile<<"__________________________"<<endl;//*/
		//in order to determine f_t
		v_max_old=v_max;
		//sort: sort unnumbered vertices by la(w) value
		v_max = int(subset[ index ]);
		k = weight[v_max -1];
		index++;
	//	outfile<<"McsM_order["<<i-1<<"]= "<<v_max<<endl;
	/*	k=0;
		for(j=0; j<n; j++)
		{
			if(labeled[j]==0)
			{
				if(weight[j] > k)
				{
				//	k=int(la[j]);
					k=weight[j];
					v_max=j+1;
				}
			}
		}*/
		
				
		McsM_order[i-2]=v_max;
		labeled[v_max-1]=i-1;
		reached[v_max-1]=1;

		//determine f_t		
		if(weight[v_max-1] <= weight[v_max_old-1])//If the conditoion is satisfied, then v_max a f_t
		{
			if(h==0)
			{
				s=new ArcNode; 
				s->adjvex=v_max;
				s->next=NULL;
				head=s;
				h++;
			}
			else
			{
				s=new ArcNode; 
				s->adjvex=v_max;
				s->next=head;
				head=s;

			}	 		
		}

		if(index == num){
			i--;
			break;
		}
		

	}

	for( ; i>1; i--)
	{
		//select:
        int j;
        //for j:=1 until k do reach(j): = \emptyset
        for(j=1; j<=k; j++){
            p=reach[j-1].first;
            ArcNode *q;
            while (p){
                q = p;
                p = p->next;
                delete q;
            }
            reach[j-1].first=NULL;
        }//*/
		//mark all unnumbered vertices unreached
		for(j=0; j<n; j++)
		{
	
			if(labeled[j]==0)
				reached[j]=0;
		}
		//update v_max's neighbor la
		//p=g.adjlist[v_max-1].firstedge;
		p=adjlist[v_max-1].firstedge;
		int w, nei;
		while(p)
		{
			w=p->adjvex;
			if(labeled[w-1]==0)
			{
			
				add_reach(w, weight[w-1]+1);
				reached[w-1]=1;
				weight[w-1]++;			
				add_madj(v_max, w);
			}
			p=p->next;	
		}
		//search:fill-in path
	
		for(j=1; j<=k; j++)
		{
			p=reach[j-1].first;
            /*    outfile<<"reachreach("<<j<<") = ";
                ArcNode *q;
                q = p;
                while (q){
                    outfile<<q->adjvex<<"  ";
                    q = q->next;
                }outfile<<endl;*/
			while(p)
			{
				//delete w from reach[j-1]
				w=p->adjvex;
				delete_reach(j);//delete reach[j-1]'s first edge
				
			//	p1=g.adjlist[w-1].firstedge;
				p1=adjlist[w-1].firstedge;
				while(p1)
				{
					nei=p1->adjvex;
					if(reached[nei-1]==0)
					{
						reached[nei-1]=1;
					//	if(la[nei-1] > j)
						if(weight[nei-1]+1 > j)
						{
						
							add_reach(nei, weight[nei-1]+1);
						
							weight[nei-1]++;			
							add_madj(v_max, nei);
						}
						else
						{
							add_reach(nei, j);
						}
					}
					p1=p1->next;
				}
				p=reach[j-1].first;
			}
		}
      /*  outfile<<endl;
        for(j=1; j<=7; j++)
		{
			p=reach[j-1].first;
                outfile<<"reach("<<j<<") = ";
                ArcNode *q;
                q = p;
                while (q){
                    outfile<<q->adjvex<<"  ";
                    q = q->next;
                }outfile<<endl;
        }
        outfile<<"__________________________"<<endl;//*/
		//in order to determine f_t
		v_max_old=v_max;
		//sort: sort unnumbered vertices by la(w) value
		k=0;
		for(j=0; j<n; j++)
		{
			if(labeled[j]==0)
			{
				if(weight[j] > k)
				{
				//	k=int(la[j]);
					k=weight[j];
					v_max=j+1;
				}
			}
		}
		
	//	outfile<<"McsM_order["<<i-1<<"]= "<<v_max<<endl;
				
		McsM_order[i-2]=v_max;
		labeled[v_max-1]=i-1;
		reached[v_max-1]=1;

		//determine f_t		
		if(weight[v_max-1] <= weight[v_max_old-1])//If the conditoion is satisfied, then v_max a f_t
		{
			if(h==0)
			{
				s=new ArcNode; n_d++;
				s->adjvex=v_max;
				s->next=NULL;
				head=s;
				h++;
			}
			else
			{
				s=new ArcNode; n_d++;
				s->adjvex=v_max;
				s->next=head;
				head=s;

			}	 		
		}
		//in order to construct junction tree of prime block, we save the madj
	/*	else//otherwise, v_max's madj is not useful, in order to reduce memory, we delete its madj 
		{
			p3=madj[v_max-1].first;
			while(p3)
			{
				p4=p3;
				p3=p3->next;
				delete p4; n_d--;
			}
			madj[v_max-1].first=NULL;
		}*/
		

	}
	

/*
	outfile<<endl;
	for(i=n-1; i>=0; i--)
	{
		outfile<<"McsM_order["<<i+1<<"]= "<<McsM_order[i]<<endl;
	}
	outfile<<"___________________"<<endl;
	for(i=0; i<n; i++)
	{
		outfile<<i+1<<"  dan diao lin ji：";
		p1=madj[i].first;
		if(p1)
		{
		while(p1)
		{
			outfile<<p1->adjvex<<", ";outfile<<"p1="<<p1<<" ";
			p1=p1->next;
		}
		outfile<<endl;
		p1=madj[i].tail;
		outfile<<"	"<<p1->adjvex<<", ";
		outfile<<"p1="<<p1<<endl;
		}
	}
	
	outfile<<"____________________"<<endl;
	for(i=0; i<n; i++)
	{
		outfile<<"labeled["<<i+1<<"]= "<<labeled[i]<<endl;
	}
	for(i=0; i<n; i++)
	{
		outfile<<"weight["<<i+1<<"]="<<weight[i];
		outfile<<"	 "<<i+1<<"  labels为：";
		p1=madj[i].first;
		while(p1)
		{
			outfile<<labeled[p1->adjvex-1]<<", ";
			p1=p1->next;
		}
		outfile<<endl;
	}
    outfile.close();//*/
/*	
	p=head;
	while(p)
	{
		outfile<<p->adjvex<<"	";
		p=p->next;
	}//*/
    
    int j;
        //for j:=1 until k do reach(j): = \emptyset
        for(j=1; j<=vertexNum; j++){
            p=reach[j-1].first;
            ArcNode *q;
            while (p){
                q = p;
                p = p->next;
                delete q;
            }
            reach[j-1].first=NULL;
        }//*/
	return head;



}



//Begin Decom_McsM_5 version 


template<class T>
void ALGraph<T>::Decom_McsM_5(ALGraph<T> & g0,   ArcNode * f_t, int n)//called by Decom_McsM_5(ALGraph<T> & g), don't copy g1，g2
{
	ArcNode *pt, * p, *p1;;//*pt point to f_t
	int number=1;//number of cliques


	pt=f_t;
	while(pt)
	{
		mfnum++;
	//	v_1=pt->adjvex;
		p=madj[/*v_1*/pt->adjvex -1].first;
		mc1=clock();
		int c=complete(p, g0);
		
		mc2=clock();
		mc=mc+double(mc2-mc1)/CLOCKS_PER_SEC;
		if(/*complete(p, g0)*/c)
		{
				
				mg1=clock();
				mg2=clock();
				mg=mg+double(mg2-mg1)/CLOCKS_PER_SEC;
				
				md1=clock();
				while(p !=NULL)
				{
					comp[p->adjvex-1]=0;
					p=p->next;
				}
				md2=clock();
				md=md+double(md2-md1)/CLOCKS_PER_SEC;

				mco1=clock();//	g2.Compchain(/*v_1*/pt->adjvex);
				g0.Compchain(/*v_1*/pt->adjvex);
				mco2=clock();
				mco=mco+double(mco2-mco1)/CLOCKS_PER_SEC;

	for(int i=0; i<vertexNum; i++)
	{
		int j;
		if(comp[i]==1 )
		{
		//	cout<<i+1<<" ";
			p=madj[i].first;
			while(p)
			{
				j=p->adjvex;
				if(comp[j-1]!=1 && comp[j-1]!=0)
					cout<<j<<" error! ";
				p=p->next;
			}
		//	cout<<endl;
		}
	}
				junction_tree(n, number);
								
			//	g0.Print_block(number);//only print j such that comp[j]=0 or comp[j]=1
				

				//find all cliques of prime block
				g0.All_Cliques_1(number);

				mdg1=clock();
				for(int j=0; j<n; j++)
				{	                       //cout<<"comp["<<j<<"]="<<comp[j-1]<<" ";
					if(comp[j]==0)	comp[j]=-1;					    	
					else if (comp[j]==1)  comp[j]++;//comp[j]=2, the vertex j-1 is deleted 
										
				}
				mdg2=clock();
				mdg=mdg+double(mdg2-mdg1)/CLOCKS_PER_SEC;
				number++;
		}
		//release the memory of madj[];
	/*	p=madj[pt->adjvex -1].first;
		while(p)
		{
			p1=p;
			p=p->next;
			delete p1; n_d--;
		}
		madj[pt->adjvex -1].first=NULL;//*/
		
		p=pt;//realse the memory of f_t
		pt=pt->next;//go to the next f_t
		delete p; n_d--;//realse the memory of f_t
	}
	mfnum++;

	for(int i=0; i<n; i++)
		if(comp[i]==-1)
			comp[i]=1;

	junction_tree(n, number);
//	g0.Print_block(number);
	//find all cliques of prime block
	g0.All_Cliques_1(number);

	cout<<"the number of the prime blocks is	"<<number<<"	by Decom_McsM_5 "<<endl;
	//realse the memory of madj[];
	
	for( i=0; i<n; i++)
	{
		p=madj[i].first;
		while(p)
		{
			p1=p;
			p=p->next;			
	//		cout<<p1->adjvex<<"	";
			delete p1; n_d--;
		}
//		cout<<endl;

	}//*/
/*	p=f_t;//realse the memory of f_t
	while(p)
	{
		p1=p;
		p=p->next;
		delete p1; n_d--;
	}//*///realse the memory 
}

template<class T>
void Decom_McsM_5(ALGraph<T> & g)
										
{
	int n=g.Num_of_vertex();
	for(int i=0; i<n; i++)
		comp[i]=-1;
			
	ArcNode *f_t = g.Mcs_M5(1);
	//_________________________
	t_Mcs=clock();
	//Mcs order time 


//	g.junction_clique(n);

	mf1=clock();
//	ArcNode *f_t=F_t(n);
	mf2=clock();

	g.Decom_McsM_5(g, f_t, n );
}


//End to reduce the memory of Mcs_M algorithm, reduce the memory of madj[i]; 
/**************************************************************************************************************************/

//Begin
//find all cliques of an undirected graph
template <class T>
int ALGraph<T>::adjacent(int i, int j)//for a graph g, judge whether i and j are adjacent, 
										//if i=j or i,j are adjacent, then return 1, otherwise return 0
{
	if(i==j)
		return 1;
	struct ArcNode * p;
	if(adjlist[i-1].vertex != 0)
	{
//		cout<<adjlist[i-1].vertex<<" ";         
		p= adjlist[i-1].firstedge;
		while(p!=NULL)
		{
			if(p->adjvex==j)		
				return 1;
			p=p->next;
		}
		
	}
	return 0;

}


template<class T>
void ALGraph<T>::All_Cliques()//Find all cliques of undirected graph g
{
	// output maximal complete subgraphs 2(cnnected, N)
	//value N; integer N;
	//Boolean array connected;
	//Commnet: the input graph is expected in the form of a symmetrical Boolean matrix connected. 
	//N is the number of nodes in the graph. 
	//The values of the diagonal elements should be true;
	
	//the inputted graph must be undirected, and
	
//	cout<<"Start to find all cliques"<<endl;
	ofstream outfile; //print the vertices of prime blocks, 
						//the first number is the number of vertices of prime blocks,  then the vertices... and then
						//print all the cliques of the prime block in disk, 
						//the first number is the number of cliques, i.e., number_cliques
						//from the 2nd to number_cliques+1, the number are the numbers of vertices of all clques
						//the following number is the vertices of all cliques
						//but the last number is not a vertex of clique, because the pointer (needle) is not given any value
	
	char filename[50]="cliques_";
	char name[30];
	sprintf(name, "%d", vertexNum);
	strcat(name, ".txt");
	strcat(filename, name);

	outfile.open(filename);
	//outfile.open(name, ios::out);
	if(!outfile)
	{
		cout<<"can,t open.\n";
		abort();
	}

	int ALL[M];
	outfile<<vertexNum<<"	";//print the number of vertices of prime blocks,

	for(int c1=1; c1<=vertexNum; c1++)
	{
		ALL[c1-1]=c1;
		outfile<<c1<<"	";
	}

	
	//	int compsub[4];
	
//	cout<<"Num_vertex="<<vertexNum<<endl;


		
	
	ccc=0;
	number_cliques=0;
	

//	number_of_recursives=0;

//	pc1=new chain;
	pc1=new ArcNode; 
	pc1->next=NULL;
//	chain *p;
	ArcNode *p;
	p=pc1;
//	pc3=new chain;
	pc3=new ArcNode; 
	pc3->next=NULL;
//	chain *pp;
	ArcNode *pp;
	pp=pc3;
	extend_version_2(ALL, 0, vertexNum);
//	p->num=number_cliques;cout<<endl;
	p->adjvex=number_cliques;//cout<<endl;
/*	while(p!=NULL)
	{
		cout<<p->num<<"	";
		p=p->next;
	}*/
	/******************************************************************/

	ArcNode * q;
	while(p!=NULL)
	{
	//	outfile<<p->num<<"	";
		outfile<<p->adjvex<<"	";
		q=p;
		p=p->next;
		delete q; 
	}
	while(pp!=NULL)
	{
//		outfile<<pp->num<<"	";
		outfile<<pp->adjvex<<"	";
		q=pp;
		pp=pp->next;
		delete q; 
	}
	                  
	

	outfile<<endl;
	
	outfile.close();
	/******************************************************************/

}




template<class T>
void ALGraph<T>::extend_version_2(int old[], int ne, int ce/*, ALGraph<T> & g, ArcNode *number_vertex*/)//called by All cliques
{
	
//	number_of_recursives++;
//	cout<<"number_of_recursives = "<<number_of_recursives<<endl;
	int new1[M];
	int nod, fixp;
	int newne, newce, i, j, count, pos, p, s, sel, minnod;
	//Commnent: The latter set of integers is local in scope but nee not be declared recursively;
	minnod=ce;  	
	nod=0;
	//determine each counter value and look for minimum:
	i=1;
	while(i<=ce && minnod!=0)
	{
		{
			p=old[i-1]; //cout<<"p=old[i-1]; =  "<<p<<endl;
			count=0;
			j=ne+1; //in orignal paper: j=ne;
			//count disconnections:
			while(j<=ce && count<minnod)//for(; j<=ce; j++)//in orignal paper: j=ne; for j:=j+1 while j<=ce 且 minnod!=0 do
			{
				if(adjacent(p, old[j-1])==0)//(connected[p-1][old[j-1]-1]==0)
				{
					count++;   
					//Save position of potential candidate:
					pos=j;   
				}
				j++;
			}

			if(count<minnod)
			{
				fixp=p;		
				minnod=count;
				if(i<=ne)
				{
					s=pos; 
				}
				else
				{
					s=i; 
					//Preincr
					nod=1;
				}//End new minimum;
			}
		}

		i++;
	}//End i;
	//Comment: If fixed point initially chosen from candidates the number of disconnecctions will be preincreased by one;
	
	//Backtrackcycle:
	
	nod=minnod+nod;		
	for( ; nod>=1; nod--) //for nod:=minnod+nod step -1 until 1 do
	{
		//Intechaange:
		p=old[s-1];	old[s-1]=old[ne];	sel=p;	old[ne]=p;//it is as following in the orignal paper: p=old[s];	old[s]=old[ne+1];	sel=p;	old[ne+1]=p;
		
		//Fill new set not:
		newne=0;	
		i=0;i++;
		for(; i<=ne; i++)  //in orignal paper: for i:=i+1 while i<=ne do
		{
			 
			if(adjacent(sel, old[i-1])==1)//(connected[sel-1][old[i-1]-1]==1)
			{
				newne++;
				new1[newne-1]=old[i-1];
				if(old[i-1]==0)
				{
					int x;
					cout<<"old[i-1]=0;   Input a number x=";
					cin>>x;
					cout<<endl;
				}
			}
		}
		//Fillnew set cand:
		newce=newne;
		i=ne+1;i++;
		 
		for( ; i<=ce; i++)
		{
			
			if(adjacent(sel, old[i-1])==1)//(connected[sel-1][old[i-1]-1]==1)
			{
				newce++;
				new1[newce-1]=old[i-1];
				if(old[i-1]==0)
				{
					int x;
					cout<<"old[i-1]=0; 222222 Input a number x=";
					cin>>x;
					cout<<endl;
				}
			}
		}
		//Add to compsub:
		ccc++;  
		compsub[ccc-1]=sel;
		
		if(newce==0)
		{
			int loc;
			number_cliques++;
		//	pc2=new chain;
			pc2=new ArcNode; n_d++;
		//	pc2->num=ccc;//record the number of vertices of cliques
			pc2->adjvex=ccc;//record the number of vertices of cliques
			pc2->next=NULL;
			pc1->next=pc2;
			pc1=pc2;
			
	/*		b=new AroNode;
			b->adjvex=ccc;
			b->next=NULL;
			number_vertex->next=b;*/


			//outstring(1, "clique = ");
	//		cout<<endl<<"clique = "<<endl;
	//		cout<<"compsub["<<number_cliques<<"]= {";
			for(loc=1; loc<=ccc; loc++)
			{
				//outinterge(1, compsub[loc]);
	//			cout<<compsub[loc-1]<<", ";//print vertices of clique in the screen
			//	pc4=new chain;
				pc4=new ArcNode; n_d++;
				pc4->next=NULL;
		//		pc3->num=compsub[loc-1];
				pc3->adjvex=compsub[loc-1];
				pc3->next=pc4;
				pc3=pc4;
			}//end output of clique
	//		cout<<"}"<<endl<<"***************************************************"<<endl;
	
		}
		else if(newne<newce)
		{
			extend_version_2(new1, newne, newce/*, b*/);
		}
		//Remove from compsub:
		ccc--; 
		//Add to not:
		ne++;
		
		if(nod>1)
		{
			//Select a candidate disconnected to the fixed point:
			s=ne;
			//Look: for candidate:
			do
			{
				s++;
			}while(adjacent(fixp, old[s-1])==1 /*&& s<ce*/);
	//		cout<<"nod="<<nod<<"  s_4 = "<<s<<"	ce_4="<<ce<<endl;
			if(s>ce)
			{
				int x;
				cout<<"s > ce "<<"  Input x=";
				cin>>x;
				cout<<endl;
			}
		}//end selection
	
	}//end Backtrackcycle
	
	
}


//End find all cliques of an undiretced graph
/*****************************************************************************************************/




template<class T>
void ALGraph<T>::extend_version_2_To_Matlab(int old[], int ne, int ce)//called by All cliques_To_Matlab()
{
	
//	number_of_recursives++;
//	cout<<"number_of_recursives = "<<number_of_recursives<<endl;
	int new1[M];
	int nod, fixp;
	int newne, newce, i, j, count, pos, p, s, sel, minnod;
	//Commnent: The latter set of integers is local in scope but nee not be declared recursively;
	minnod=ce;  	
	nod=0;
	//determine each counter value and look for minimum:
	i=1;
	while(i<=ce && minnod!=0)
	{
		{
			p=old[i-1]; //cout<<"p=old[i-1]; =  "<<p<<endl;
			count=0;
			j=ne+1; //in orignal paper: j=ne;
			//count disconnections:
			while(j<=ce && count<minnod)//for(; j<=ce; j++)//in orignal paper: j=ne; for j:=j+1 while j<=ce 且 minnod!=0 do
			{
				if(adjacent(p, old[j-1])==0)//(connected[p-1][old[j-1]-1]==0)
				{
					count++;   
					//Save position of potential candidate:
					pos=j;   
				}
				j++;
			}

			if(count<minnod)
			{
				fixp=p;		
				minnod=count;
				if(i<=ne)
				{
					s=pos; 
				}
				else
				{
					s=i; 
					//Preincr
					nod=1;
				}//End new minimum;
			}
		}

		i++;
	}//End i;
	//Comment: If fixed point initially chosen from candidates the number of disconnecctions will be preincreased by one;
	
	//Backtrackcycle:
	
	nod=minnod+nod;		
	for( ; nod>=1; nod--) //for nod:=minnod+nod step -1 until 1 do
	{
		//Intechaange:
		p=old[s-1];	old[s-1]=old[ne];	sel=p;	old[ne]=p;//it is as following in the orignal paper: p=old[s];	old[s]=old[ne+1];	sel=p;	old[ne+1]=p;
		
		//Fill new set not:
		newne=0;	
		i=0;i++;
		for(; i<=ne; i++)  //in orignal paper: for i:=i+1 while i<=ne do
		{
			 
			if(adjacent(sel, old[i-1])==1)//(connected[sel-1][old[i-1]-1]==1)
			{
				newne++;
				new1[newne-1]=old[i-1];
				if(old[i-1]==0)
				{
					int x;
					cout<<"old[i-1]=0;   Input a number x=";
					cin>>x;
					cout<<endl;
				}
			}
		}
		//Fillnew set cand:
		newce=newne;
		i=ne+1;i++;
		 
		for( ; i<=ce; i++)
		{
			
			if(adjacent(sel, old[i-1])==1)//(connected[sel-1][old[i-1]-1]==1)
			{
				newce++;
				new1[newce-1]=old[i-1];
				if(old[i-1]==0)
				{
					int x;
					cout<<"old[i-1]=0; 222222 Input a number x=";
					cin>>x;
					cout<<endl;
				}
			}
		}
		//Add to compsub:
		ccc++;  
		compsub[ccc-1]=sel;
		
		if(newce==0)
		{
			int loc;
			number_cliques++;
		//	pc2=new chain;
			pc2=new ArcNode; n_d++;
		//	pc2->num=ccc;//record the number of vertices of cliques
			pc2->adjvex=ccc;//record the number of vertices of cliques
			pc2->next=NULL;
			pc1->next=pc2;
			pc1=pc2;
			
	/*		b=new AroNode;
			b->adjvex=ccc;
			b->next=NULL;
			number_vertex->next=b;*/


			//outstring(1, "clique = ");
	//		cout<<endl<<"clique = "<<endl;
	//		cout<<"compsub["<<number_cliques<<"]= {";
			for(loc=1; loc<=ccc; loc++)
			{
				//outinterge(1, compsub[loc]);
	//			cout<<compsub[loc-1]<<", ";//print vertices of clique in the screen
			//	pc4=new chain;
				pc4=new ArcNode; n_d++;
				pc4->next=NULL;
		//		pc3->num=compsub[loc-1];
				pc3->adjvex=compsub[loc-1];
				pc3->next=pc4;
				pc3=pc4;
			}//end output of clique
	//		cout<<"}"<<endl<<"***************************************************"<<endl;
	
		}
		else if(newne<newce)
		{
			extend_version_2(new1, newne, newce/*, b*/);
		}
		//Remove from compsub:
		ccc--; 
		//Add to not:
		ne++;
		
		if(nod>1)
		{
			//Select a candidate disconnected to the fixed point:
			s=ne;
			//Look: for candidate:
			do
			{
				s++;
			}while(adjacent(fixp, old[s-1])==1 /*&& s<ce*/);
	//		cout<<"nod="<<nod<<"  s_4 = "<<s<<"	ce_4="<<ce<<endl;
			if(s>ce)
			{
				int x;
				cout<<"s > ce "<<"  Input x=";
				cin>>x;
				cout<<endl;
			}
		}//end selection
	
	}//end Backtrackcycle
	
	
}


//End find all cliques of an undiretced graph
/*****************************************************************************************************/

//Begin to find all cliques of prime blocks
//We will apply the golable variable comp[M] in the decomposition step
 
template<class T>
void ALGraph<T>::All_Cliques_1( int number)//Find all cliques of prime blocks
{
	// output maximal complete subgraphs 2(cnnected, N)
	//value N; integer N;
	//Boolean array connected;
	//Commnet: the input graph is expected in the form of a symmetrical Boolean matrix connected. 
	//N is the number of nodes in the graph. 
	//The values of the diagonal elements should be true;
	
	//the inputted graph must be undirected, and

//	cout<<"Start to find all cliques"<<endl;

	ArcNode *p, *q;
	ofstream outfile; //print the vertices of prime blocks, 
						//the first number is the number of vertices of prime blocks,  then the vertices... and then
						//print all the cliques of the prime block in disk, 
						//the first number is the number of cliques, i.e., number_cliques
						//from the 2nd to number_cliques+1, the number are the numbers of vertices of all clques
						//the following number is the vertices of all cliques
						//but the last number is not a vertex of clique, because the pointer (needle) is not given any value
	char filename[50]="clique_";
	char name[20];
	sprintf(name, "%d", number);
	strcat(name, ".txt");
	strcat(filename, name);

	outfile.open(filename);
	//outfile.open(name, ios::out);
	if(!outfile)
	{
		cout<<"can,t open.\n";
		abort();
	}

	int ALL[M];

	
//	for(int c1=1; c1<=vertexNum; c1++)
//		ALL[c1-1]=c1;
	int n=0;
	for(int i=0; i<vertexNum; i++)
	{
		if(comp[i]==1 || comp[i]==0)
		{
			ALL[n]=i+1;
			n++;
		}
	}


	outfile<<n<<"	";//print the number of vertices of prime blocks,

	for(i=0; i<n; i++)
		outfile<<ALL[i]<<"	";//print the vertices of prime blocks,
		
	
/*	for(i=0; i<n; i++)
		cout<<ALL[i]<<"	";
	cout<<endl;*/

	ccc=0;
	number_cliques=0;
	


	pc1=new ArcNode;n_d++;
	pc1->next=NULL;
	
	p=pc1;
	pc3=new ArcNode; n_d++;
	pc3->next=NULL;
	ArcNode *pp;
	pp=pc3;
	extend_version_2(ALL, 0, n);
//	p->num=number_cliques;cout<<endl;
	p->adjvex=number_cliques;//cout<<endl;

	

	while(p!=NULL)
	{
		outfile<<p->adjvex<<"	";
		q=p;
		p=p->next;
		delete q; n_d--;
	}
	while(pp!=NULL)
	{
		outfile<<pp->adjvex<<"	";
		q=pp;
		pp=pp->next;
		delete q; n_d--;
	}
	                  
	outfile<<endl;
	
	outfile.close();
	/******************************************************************/

}

//End to find all cliques of prime blocks
/*****************************************************************************************************/
/*************************************************************************/
//Begin to construct a junction tree of the global fill_in graph
//and begin to find cliques in each node of junction tree 

template<class T>
void ALGraph<T>::junction_clique(int n)//n is the number of vertices of glabel graph
						//print the junction tree of fill_in graph, 
						//the first number, i.e., 0, is the parent  of the first node, 
						//then the number of vertices of the first node, and then vertices...
						//Then begin print the second clique's parent, its verticies number, and its vertices
						//the following is the third, forth, ...
						//And the last number is the number of nodes of the junction tree 
						
						//and print cliques in each node of junction tree
{
	int j, v1, v2, num_clique=0, label=n;//v1, v2 are two vertices, and label=n is used to number the vertices in the first clique 
	int clique[M];// is the funciton defined in orginal paper
	ArcNode *p, *q;
//_____________________________________________________________________________________________//*/
	//find cliques of each node of junction tree
	ofstream print; //print cliques of node of junction tree
					//before which, we print all vertices of the graph
	char filename1[50]="jun_cliques_";
	char name1[20];
	sprintf(name1, "%d", n);
	strcat(name1, ".txt");
	strcat(filename1, name1);

	print.open(filename1);
	//outfile.open(name, ios::out);
	if(!print)
	{
		cout<<"can,t open.\n";
		abort();
	}
	print<<n<<" ";
	for(j=1; j<=n; j++)
		print<<j<<" ";
	int ALL[M];
	int nn;//nn is index of ALL[]
//	int num_vertex_each;//number of shuzi in each printing cliques
//_____________________________________________________________________________________________//*/
	
	ofstream outfile; 
	char filename[50]="jun_cliques_tree_";
	char name[20];
	sprintf(name, "%d", n);
	strcat(name, ".txt");
	strcat(filename, name);

	outfile.open(filename);
	//outfile.open(name, ios::out);
	if(!outfile)
	{
		cout<<"can,t open.\n";
		abort();
	}
	outfile<<0<<" ";//print the number of parent clique of the first clique,

	v2=McsM_order[n-1];
	j=n-2;
	int r_num;//r is the number of vertices in R=C\S, where C is the clique and S is the separator
	while(j>=0)
	{
		{
			v1=McsM_order[j];
			if(weight[v1-1]<=weight[v2-1])
			{
				num_clique++;
				//print the number of vertices of the clique, and all vertices
				nn=0;//nn is index of ALL[]

				outfile<<weight[v2-1]+1<<" ";
				clique[v2-1]=num_clique;
				r_num=1;
				p=madj[v2-1].first;
				while(p)
				{
					ALL[nn]=p->adjvex;
					nn++;

					outfile<<p->adjvex<<" ";					
					if(labeled[p->adjvex -1]<=label)
					{
						clique[p->adjvex -1]=num_clique;
						r_num++;
					}
					p=p->next;
				}
				outfile<<v2<<" "<<r_num<<" ";
				label=labeled[v1-1];
				//print the number of parent clique of next clique
				p=madj[v1-1].tail;
				outfile<<clique[p->adjvex -1]<<" ";
				
				ALL[nn]=v2;
				nn++;

//_____________________________________________________________________________________________//*/
/*	int nn=0;
	for(int i=0; i<vertexNum; i++)
	{
		if(comp[i]==1 || comp[i]==0)
		{
			ALL[nn]=i+1;
			n++;
		}
	}//*/

//	num_vertex_each=0;

	ccc=0;
	number_cliques=0;
	


	pc1=new ArcNode;n_d++;
	pc1->next=NULL;
	
	p=pc1;
	pc3=new ArcNode; n_d++;
	pc3->next=NULL;
	ArcNode *pp;
	pp=pc3;
	extend_version_2(ALL, 0, nn);
//	p->num=number_cliques;cout<<endl;
	p->adjvex=number_cliques;//cout<<endl;

	

	while(p!=NULL)
	{
	//	num_vertex_each++;
		print<<p->adjvex<<" ";
		q=p;
		p=p->next;
		delete q; n_d--;
	}
	while(pp!=NULL)
	{
	//	num_vertex_each++;
		print<<pp->adjvex<<" ";
		q=pp;
		pp=pp->next;
		delete q; n_d--;
	}
//	print<<num_vertex_each<<" ";
//_____________________________________________________________________________________________//*/
			
			}
			v2=v1;
		}
		j--;
	}

	nn=0;
	num_clique++;
	outfile<<weight[v2-1]+1<<" ";
	r_num=1;
	p=madj[v2-1].first;
	while(p)
	{
		ALL[nn]=p->adjvex;
		nn++;

		outfile<<p->adjvex<<" ";
		if(labeled[p->adjvex -1]<=label)
		{
			clique[p->adjvex -1]=num_clique;
			r_num++;
		}
		p=p->next;
	}
	outfile<<v2<<" "<<r_num<<" "<<num_clique<<endl;	
	outfile.close();//*/

//_____________________________________________________________________________________________//*/
	ALL[nn]=v2;
	nn++;

//	num_vertex_each=0;
	ccc=0;
	number_cliques=0;
	pc1=new ArcNode;n_d++;
	pc1->next=NULL;
	
	p=pc1;
	pc3=new ArcNode; n_d++;
	pc3->next=NULL;
	ArcNode *pp;
	pp=pc3;
	extend_version_2(ALL, 0, nn);
	p->adjvex=number_cliques;//cout<<endl;
	

	while(p!=NULL)
	{
	//	num_vertex_each++;
		print<<p->adjvex<<" ";
		q=p;
		p=p->next;
		delete q; n_d--;
	}
	while(pp!=NULL)
	{
	//	num_vertex_each++;
		print<<pp->adjvex<<" ";
		q=pp;
		pp=pp->next;
		delete q; n_d--;
	}
//	print<<num_vertex_each<<" ";
//_____________________________________________________________________________________________//*/
//	print<<num_clique<<endl;
	
	print.close();//*/

	
}

double my_rand_seeds()
{
	time_t ttt; 
	srand((unsigned) time(&ttt)); 
	return rand();
}


