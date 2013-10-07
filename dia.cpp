#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <random>
#include <time.h>
#include <vector>
#include <string>
#include <omp.h>
using namespace std;
long double dx[8], x[8]={1e10, 0,0,1e3,0,0,0,0};
//  /*
long long jca=0;
bool DIAnoSSA=true;
long double t0=1e10, tmax=1e12, pw=750/10.0, pm=500/10.0;
long double cor=0.6, cod=3.6,cod2=3.6, cain=6, cain2=5, cor2=0.5; // 0.6 3.6 3.6 6 5 0.5
long double kUni = 0.005/1.0, kSpec1wi = 0.085/1.0, kSpec1mu = 0.001/1.0, kSpec2wi = 0.005/1.0, kSpec2mu = 0.085/1.0;
long double ncoeff = 2*5*1e-12*t0, cmutwtm = 0.3, cmutmtw = 0.3;
long double d=2.5, delta=0.02, cv=5*0.001*1e4;
long double ATW[3];
long double c1w = kSpec1wi/(kSpec1wi + kSpec1mu), c1m = kSpec1mu/(kSpec1wi + kSpec1mu), c2w = kSpec2wi/(kSpec2wi + kSpec2mu), c2m = kSpec2mu/(kSpec2wi + kSpec2mu);
long double speed = 0.3, MAX_POP_VIR = 1e10;
long double bph1=1, bph2=1;
//  */
vector<long double> data[9];
void ofdata(long double T)
{
	for(int i=0; i<8; i++)data[i].push_back(1+x[i]);
	data[8].push_back(T);
}
long double uplim(long double x, long double y){x/=y; if(!x>1.0)return x; return 1.0;}
void o2fdata()
{
	ofstream out;
	out.open("datascript.m", ofstream::out);
	out << "AX = [";
	for(long long i=0; i<8; i++)
	{
		for(long long j=0; j<data[i].size(); j++)
		{
			out << data[i][j]<<" ";
		}
		out << ";"<<endl;
	}
	out<<"];"<<endl<<endl;
	out << "T = [";
	for(long long i=0; i<data[8].size(); i++)out<<data[8][i]<<" ";
	out<<"];"<<endl;
	out << "AX = transpose(AX);"<<endl;
	out <<"figure"<<endl;
	out <<"semilogy(T,AX);";
	out.close();
}

inline void x2dx(){for(int i=0; i<8; i++){x[i] = 0; if(dx[i]>0)x[i] = dx[i];}}

inline double geo(long double x){if(x>0)return 1.0; return 0.0;}
inline long double sign(long double x){if(x>0)return 1; if(x<0)return -1; return 0;}
inline long double mdl(long double x){x*=sign(x); return x;}
inline long double min(long double x,long double y){if(x<y)return x; return y;}
long double rnd(){return ((long double)rand())/((long double)(RAND_MAX));}
void mdf_dh(){bph1=bph2=1;}
// nonuse
void h_make_d(long double SA)
{
	//bph1 = ((cyt[0]+cyt[1]+cyt[2])/1e8);
	//bph2 = cyt[3]/1e8;
	// уберём влияние цитокинов
	//bph1 = 1;
	//bph2 = 1;
	dx[0] += SA*((bph1*d*x[0]*(1-(x[0]/t0))-bph2*ncoeff*x[0]*(x[3]+x[4])));
	dx[1] += SA*(bph2*ncoeff*x[0]*x[3] - delta*x[1] - kUni*x[5]*(x[1]) - kSpec1wi*x[6]*(x[1]) - kSpec2wi*x[7]*(x[1]));
	dx[2] += SA*(bph2*ncoeff*x[0]*x[4] - delta*x[2] - kUni*x[5]*(x[2]) - kSpec1mu*x[6]*(x[2]) - kSpec2mu*x[7]*(x[2]));
	dx[3] += SA*(cmutmtw*pm*x[2] + pw*x[1]*(1-cmutwtm)-cv*x[3]-ncoeff*x[0]*x[3]);
	dx[4] += SA*(cmutwtm*pw*x[1] + pm*x[2]*(1-cmutmtw)-cv*x[4]-ncoeff*x[0]*x[4]);
	dx[5] += SA*(cain+bph1*cor*min(1,(0.5*(x[1]+x[2]))/1e4)*x[5]*(1-x[5]/tmax)-cod*x[5]);
	dx[6] += SA*(cain2+bph1*cor2*min(1,(c1w*x[1] + c1m*x[2])/1e4)*x[6]*(1-x[6]/tmax)-cod2*x[6]);
	dx[7] += SA*(cain2+bph1*cor2*min(1,(c2m*x[2] + c2w*x[1])/1e4)*x[7]*(1-x[7]/tmax)-cod2*x[7]);
	return;
}
inline void make_d(long double SA)
{
	dx[0] += SA*((d*x[0]*(1-(x[0]/t0))-ncoeff*uplim(x[0],(x[3]+x[4]))*(x[3]+x[4])));
	dx[1] += SA*(((ncoeff*uplim(x[0],(x[3]+1))*x[3]-delta*x[1]-kUni*x[1]*x[5]/(x[1]+1) - kSpec1wi*x[1]*x[6]/(x[1]+1) - kSpec2wi*x[1]*x[7]/(x[1]+1))));
	dx[2] += SA*(((ncoeff*uplim(x[0],(x[4]+1))*x[4]-delta*x[2]-kUni*x[2]*x[5]/(x[2]+1) - kSpec1mu*x[2]*x[6]/(x[2]+1) - kSpec2mu*x[2]*x[7]/(x[2]+1))));
	dx[3] += SA*(((cmutmtw*pm*x[2] + pw*x[1]*(1-cmutwtm)-cv*x[3]-ncoeff*uplim(x[0],(x[3]+1))*x[3])));
	dx[4] += SA*(((cmutwtm*pw*x[1] + pm*x[2]*(1-cmutmtw)-cv*x[4]-ncoeff*uplim(x[0],(x[4]+1))*x[4])));
	dx[5] += SA*(cain+cor*(((x[1]+x[2])/2.0)/1e2)*x[5]*(1-x[5]/tmax)-cod*x[5]);
	dx[6] += SA*(cain2+cor2*((c1w*x[1] + c1m*x[2])/1e2)*x[6]*(1-x[6]/tmax)-cod2*x[6]);
	dx[7] += SA*(cain2+cor2*((c2m*x[2] + c2w*x[1])/1e2)*x[7]*(1-x[7]/tmax)-cod2*x[7]);
	
}
// выделим память под алгоритм Гилеспи
long double ret, sum, RND, RND2;
long double SREC[100];
inline long double make_ssa()
{
	long double ret, sum, RND = rnd(), RND2;
	int select=-1;
	RND2 = rnd();
	while(RND2==0){RND2 = rnd();}
	SREC[0] = d*x[0]*(1.0-(x[0]/t0)); // replic

	SREC[1] = ncoeff*uplim(x[0],x[3]+1)*x[3];
	SREC[2] = delta*x[1];
	SREC[3] = kUni*x[1]*x[5]/(x[1]+1);
	SREC[4] = kSpec1wi*x[1]*x[6]/(x[1]+1);
	SREC[5] = kSpec2wi*x[1]*x[7]/(x[1]+1);
	
	SREC[6] = ncoeff*uplim(x[0],x[4]+1)*x[4];
	SREC[7] = delta*x[1];
	SREC[8] = kUni*x[2]*x[5]/(x[2]+1);
	SREC[9] = kSpec1mu*x[2]*x[6]/(x[2]+1);
	SREC[10] = kSpec2mu*x[2]*x[7]/(x[2]+1);
	
	SREC[11] = cmutmtw*pm*x[2];
	SREC[12] = pw*x[1]-pw*x[1]*cmutwtm;
	SREC[13] = cv*x[3];
	
	SREC[14] = cmutwtm*pw*x[1];
	SREC[15] = pm*x[2]-pm*x[2]*cmutmtw;
	SREC[16] = cv*x[4];
	
	SREC[17] = cain;
	SREC[18] = cor*(((x[1]+x[2])/2.0)/1e2)*x[5]*(1-x[5]/tmax);
	SREC[19] = cod*x[5];
	
	SREC[20] = cain2;
	SREC[21] = cor2*((c1w*x[1] + c1m*x[2])/1e2)*x[6]*(1-x[6]/tmax);
	SREC[22] = cod2*x[6];
	
	SREC[23] = cain2;
	SREC[24] = cor2*((c2m*x[2] + c2w*x[1])/1e2)*x[7]*(1-x[7]/tmax);
	SREC[25] = cod2*x[7];
	sum = 0;//modul(SREC[0])+modul(SREC[1])+modul(SREC[2]);
	for(int i=0; i<26; i++){sum+=SREC[i];}
	//if(sum<1e-10)return 0;
	ret = (1.0/sum)*log(1.0/RND2);
	sum*=RND;
	while(sum>0)
	{
		select++;
		sum-=SREC[select];
	}
	// #1
	if(select==0){dx[0]++;}
	// #2
	if(select==1){dx[1]++; dx[0]--; dx[3]--;}
	if(select==2){dx[1]--;}
	if(select==3){dx[1]--;}
	if(select==4){dx[1]--;}
	if(select==5){dx[1]--;}
	// #3
	if(select==6){dx[2]++; dx[0]--; dx[4]--;}
	if(select==7){dx[2]--;}
	if(select==8){dx[2]--;}
	if(select==9){dx[2]--;}
	if(select==10){dx[2]--;}
	// #4
	if(select==11){dx[3]++;}
	if(select==12){dx[3]++;}
	if(select==13){dx[3]--;}
	// #5
	if(select==14){dx[4]++;}
	if(select==15){dx[4]++;}
	if(select==16){dx[4]--;}
	// #6
	if(select==17){dx[5]++;}
	if(select==18){dx[5]++;}
	if(select==19){dx[5]--;}
	// #7
	if(select==20){dx[6]++;}
	if(select==21){dx[6]++;}
	if(select==22){dx[6]--;}
	// #8
	if(select==23){dx[7]++;}
	if(select==24){dx[7]++;}
	if(select==25){dx[7]--;}
	return ret;
}

long double o_make_ssa()
{
	RND = rnd();
	RND2 = rnd();
	while(RND2==0){RND2 = rnd();}
	// 1
	SREC[0] = bph1*d*x[0]*(1-x[0]/t0); // replic
	// 2
	SREC[1] = bph2*ncoeff*x[0]*x[3]; // заражение здоровых клеток
	SREC[2] = delta*x[1]; // гибель отработавших ресурс
	SREC[3] = kUni*x[5]*(x[1]); // уничтожение неспецифическими лимфоцитами
	SREC[4] = kSpec1wi*(x[1])*x[6]; // уничтожение W-специфичными лимфоцитами
	SREC[5] = kSpec2wi*x[7]*(x[1]); // уничтожение M-специфичными лимфоцитами
	// 3
	SREC[6] = bph2*ncoeff*x[0]*x[4];
	SREC[7] = delta*x[2];
	SREC[8] = kUni*x[5]*(x[2]);
	SREC[9] = kSpec1mu*x[6]*(x[2]);
	SREC[10] = kSpec2mu*x[7]*(x[2]);
	// 4
	SREC[11] = cmutmtw*pm*x[2] + pw*x[1]-pw*x[1]*cmutwtm;
	SREC[12] = cv*x[3];
	// 5
	SREC[13] = cmutwtm*pw*x[1] + (pm*x[2]-pm*x[2]*cmutmtw);
	SREC[14] = cv*x[4];
	// 6
	SREC[15] = cain;
	SREC[16] = cor*min(1,((x[1]+x[2])/2.0)/1e4)*x[5]*(1-x[5]/tmax);
	SREC[17] = cod*x[5];
	// 7
	SREC[18] = cain2;
	SREC[19] = cor2*min(1,(c1w*x[1] + c1m*x[2])/1e4)*x[6]*(1-x[6]/tmax);
	SREC[20] = cod2*x[6];
	// 8
	SREC[21] = cain2;
	SREC[22] = cor2*min(1,(c2m*x[2] + c2w*x[1])/1e4)*x[7]*(1-x[7]/tmax);
	SREC[23] = cod2*x[7];
	sum = 0;
	for(int i=0; i<24; i++){sum+=SREC[i];}
	ret = (1.0/sum)*log(1.0/RND2);
	sum*=RND;
	int select=-1;
	while(sum>0){select++;sum-=SREC[select];}
	// #1
	if(select==0){dx[0]++;}
	// #2
	if(select==1){dx[1]++; dx[0]--;}
	if(select==2){dx[1]--;}
	if(select==3){dx[1]--;}
	if(select==4){dx[1]--;}
	if(select==5){dx[1]--;}
	// #3
	if(select==6){dx[2]++; dx[0]--;}
	if(select==7){dx[2]--;}
	if(select==8){dx[2]--;}
	if(select==9){dx[2]--;}
	if(select==10){dx[2]--;}
	// #4
	if(select==11){dx[3]++;}
	if(select==12){dx[3]--;}
	// #5
	if(select==13){dx[4]++;}
	if(select==14){dx[4]--;}
	// #6
	if(select==15){dx[5]++;}
	if(select==16){dx[5]++;}
	if(select==17){dx[5]--;}
	// #7
	if(select==18){dx[6]++;}
	if(select==19){dx[6]++;}
	if(select==20){dx[6]--;}
	// #8
	if(select==21){dx[7]++;}
	if(select==22){dx[7]++;}
	if(select==23){dx[7]--;}
	return ret;
}

long double get_ssa_fs()
{
	long double ret, sum, RND = random(), RND2;
	int select=-1;
	RND2 = random();
	while(RND2==0){RND2 = random();}
	SREC[0] = d*x[0]*(1.0-(x[0]/t0)); // replic
	
	SREC[1] = ncoeff*uplim(x[0],x[3]+1)*x[3];
	SREC[2] = delta*x[1];
	SREC[4-1] = kUni*x[1]*x[5]/(x[1]+1);
	SREC[5-1] = kSpec1wi*x[1]*x[6]/(x[1]+1);
	SREC[6-1] = kSpec2wi*x[1]*x[7]/(x[1]+1);
	
	SREC[7-1] = ncoeff*uplim(x[0],x[4]+1)*x[4];
	SREC[7] = delta*x[1];
	SREC[8] = kUni*x[2]*x[5]/(x[2]+1);
	SREC[9] = kSpec1mu*x[2]*x[6]/(x[2]+1);
	SREC[10] = kSpec2mu*x[2]*x[7]/(x[2]+1);
	
	SREC[11] = cmutmtw*pm*x[2];
	SREC[12] = pw*x[1]-pw*x[1]*cmutwtm;
	SREC[13] = cv*x[3];
	
	SREC[14] = cmutwtm*pw*x[1];
	SREC[15] = pm*x[2]-pm*x[2]*cmutmtw;
	SREC[16] = cv*x[4];
	
	SREC[17] = cain;
	SREC[18] = cor*(((x[1]+x[2])/2.0)/1e2)*x[5]*(1-x[5]/tmax);
	SREC[19] = cod*x[5];
	
	SREC[20] = cain2;
	SREC[21] = cor2*((c1w*x[1] + c1m*x[2])/1e2)*x[6]*(1-x[6]/tmax);
	SREC[22] = cod2*x[6];
	
	SREC[23] = cain2;
	SREC[24] = cor2*((c2m*x[2] + c2w*x[1])/1e2)*x[7]*(1-x[7]/tmax);
	SREC[25] = cod2*x[7];
	sum = 0;
	for(int i=0; i<26; i++){sum+=SREC[i];}
	return sum;
}

long double o_get_ssa_fs()
{
	RND = rnd();
	RND2 = rnd();
	while(RND2==0){RND2 = rnd();}
	// 1
	SREC[0] = bph1*d*x[0]*(1-(x[0]/t0)); // replic
	// 2
	SREC[1] = bph2*ncoeff*x[0]*x[3]; // заражение здоровых клеток
	SREC[2] = delta*x[1]; // гибель отработавших ресурс
	SREC[3] = kUni*x[5]*x[1]; // уничтожение неспецифическими лимфоцитами
	SREC[4] = kSpec1wi*x[1]*x[6]; // уничтожение W-специфичными лимфоцитами
	SREC[5] = kSpec2wi*x[7]*(x[1]); // уничтожение M-специфичными лимфоцитами
	// 3
	SREC[6] = bph2*ncoeff*x[0]*x[4];
	SREC[7] = delta*x[2];
	SREC[8] = kUni*x[5]*(x[2]);
	SREC[9] = kSpec1mu*x[6]*(x[2]);
	SREC[10] = kSpec2mu*x[7]*(x[2]);
	// 4
	SREC[11] = cmutmtw*pm*x[2] + pw*x[1]-pw*x[1]*cmutwtm;
	SREC[12] = cv*x[3];
	// 5
	SREC[13] = cmutwtm*pw*x[1] + (pm*x[2]-pm*x[2]*cmutmtw);
	SREC[14] = cv*x[4];
	// 6
	SREC[15] = cain;
	SREC[16] = cor*min(1,((x[1]+x[2])/2.0)/1e2)*x[5]*(1-x[5]/tmax);
	SREC[17] = cod*x[5];
	// 7
	SREC[18] = cain2;
	SREC[19] = cor2*min(1,(c1w*x[1] + c1m*x[2])/1e2)*x[6]*(1-x[6]/tmax);
	SREC[20] = cod2*x[6];
	// 8
	SREC[21] = cain2;
	SREC[22] = cor2*min(1,(c2m*x[2] + c2w*x[1])/1e2)*x[7]*(1-x[7]/tmax);
	SREC[23] = cod2*x[7];
	sum = 0;
	for(int i=0; i<24; i++){sum+=SREC[i];}
	return sum;
}

bool autocheck(long double LSS)
{
	if(x[0]<1e3)return true;
	if(get_ssa_fs()<LSS)return true;
	if(x[1]+x[2]+x[3]+x[4]<1e3)return true;
	return false;
}

int main(int argc, char** argv)
{
	for(int i=0; i<100; i++)SREC[i]=0;
	ifstream fin;
	fin.open("iData.txt", ifstream::in);
	for(int i=0; i<8; i++){fin >> x[i]; cout << "x["<<i<<"] = "<<x[i]<<endl;}
	fin.close();
	for(int i=0; i<8; i++)dx[i]=x[i];
	double FP, SA, LSS,TOI, O2F;
	//cout << "FinalPoint, StepAccuracy, LevelStochastickSpeed"<<endl;
	//cin >> FP >> SA >> LSS;
	string nrand;
	char c_nrand[10];
	if(argc!=7)return -1;
	sscanf(argv[1], "%lf", &FP);
	sscanf(argv[2], "%lf", &SA);
	sscanf(argv[3], "%lf", &LSS);
	sscanf(argv[4], "%lf", &O2F);
	sscanf(argv[5], "%lf", &TOI);
	sscanf(argv[6], "%s", &c_nrand[0]);
	bool b_nrand = true;
	nrand = string(c_nrand);
	if(nrand=="nrand"){b_nrand = false;}
	cout << FP << ":" << SA << ":" << O2F << ":" << LSS << endl;
	long double TM[4]={0,0,0,0};
	long double GTO[2]={omp_get_wtime(),GTO[0]};
	long double SP=0, O2FP=0;
	long double hdi;
	cout <<"lqs_#1"<< get_ssa_fs() << endl;
	make_ssa();
	cout <<"lqs_#2"<< get_ssa_fs() << endl;
	while(SP<FP)
	{
		if(SP>=O2FP)
		{
			ofdata(SP);
			O2FP+=O2F;
		}
		GTO[1] = omp_get_wtime();
		if(GTO[1]-GTO[0]>TOI)
		{
			system("clear");
			cout << "Time: "<<SP<<" / "<<FP<<" sec"<<endl;
			get_ssa_fs();
			cout << "MODE: ";
			if((hdi=get_ssa_fs())<LSS){cout << "SSA";}else{cout << "DIA";}
			cout << "("<<hdi<<" ev/day)"<<"["<<jca<<" jumps]"<<endl;
			GTO[0] = omp_get_wtime();
			for(int i=0; i<8; i++)
			{
				cout <<"["<<i<<"]: "<< x[i]<<";"<<endl;
			}
			for(int i=0; i<26; i++)
			{
				cout <<"{"<<i<<"}: "<< SREC[i]<<";"<<endl;
			}
		}
		if(b_nrand && autocheck(LSS))
		{
			// SSA algorythm
			SP+=make_ssa();
			if(DIAnoSSA){jca++; DIAnoSSA=false;}
		}
		else
		{
			//TM[0] = omp_get_wtime();
			// Determin algorythm
			make_d(SA);
			if(!DIAnoSSA){jca++; DIAnoSSA=true;}
			SP+=SA;
			//TM[1] = omp_get_wtime();
			//TM[2]++;
			//TM[3]*=TM[2]-1;
			//TM[3]+=TM[1]-TM[0];
			//TM[3]/=TM[2];
			//cout << "Midtime ~" << TM[3] << " sec."<<endl;
		}
		x2dx();
	}
	o2fdata();
	return 0;
}
