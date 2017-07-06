#include <iostream>
#include <stdlib.h>         	
#include <math.h>
#include <time.h>  
#include <iomanip>
using namespace std;



const int MAX_REQUESTS		= 200;
const int MAX_PROCESSORS	= 5;
const int MUT_NUM			= 3;
const int CROSS_POINT		= 3;
const int pop_size = 60;
const int CC	   = 100;
const int MAX_NUM			= 2000;
const double FOG_EDGE_DELTA = 0;


typedef struct Chrom
{
	int bit[MAX_REQUESTS *MAX_PROCESSORS];
	int fsb;
	double avg_ltncy;
	double fit;
	double norm_fit;
	int theta[MAX_REQUESTS*MAX_PROCESSORS];
	double sTime[MAX_REQUESTS];
	double pTime[MAX_REQUESTS];
	double tTime[MAX_REQUESTS];
	double ltncy[MAX_REQUESTS];
} chrom;

typedef struct Request
{
	double init_time;		//referenced to 0
	double size;
	double deadline;			//as a duration
	double priority;
	double sTime;			
	double pTime;		
	double tTime;	
	double ltncy;
} request;

typedef struct Processor
{
	double speed;
	double delta;			
} processor;

//----------------------------------------------------------------------------

void define_prob(int num_requests, int num_processors);
void initialize_pop(chrom * popcurrent, int pop_size, int num_requests, int num_processors);
void calc_fittness(chrom * pop, int pop_size, int num_requests, int num_processors);
int check_fsb(chrom * pop, int ind, int num_requests, int num_processors);
double latency(chrom * pop, int ind, int num_requests, int num_processors);
void display_pop(chrom * pop, int pop_size, int num_requests, int num_processors);
double fsb_percentage(chrom * pop, int pop_size);
void sort(chrom * popcurrent, int pop_size);
int select_chrom(chrom * popnext, int pop_size);
void crossover(chrom * popnext, int pop_size, int num_requests, int num_processors);
void mutation(chrom * popnext, int pop_size, int num_requests, int num_processors);
request * req_set;
processor * proc_set;
int * Beta;

//----------------------------------------------------------------------------

double Init_0[MAX_REQUESTS] = { 2, 4, 8, 10, 13, 18, 21, 23, 25, 27, 30, 32, 34, 36, 36, 37, 38, 39, 40, 43, 45, 49, 52, 53, 56, 60, 61, 62, 63, 66, 70, 70, 72, 73, 75, 78, 80, 80, 83, 88, 88, 89, 90, 90, 91, 92, 93, 96, 96, 99, 101, 104, 106, 109, 112, 118, 122, 124, 128, 128, 133, 134, 137, 137, 140, 142, 143, 144, 147, 149, 151, 151, 154, 156, 162, 167, 169, 172, 173, 175, 178, 180, 181, 181, 183, 186, 189, 191, 193, 196, 198, 200, 200, 205, 210, 211, 214, 217, 219, 220, 221, 222, 225, 225, 226, 228, 233, 237, 240, 241, 243, 246, 248, 250, 252, 252, 255, 260, 261, 263, 266, 267, 269, 271, 271, 273, 275, 276, 278, 283, 284, 287, 287, 290, 290, 292, 295, 298, 298, 299, 299, 301, 302, 305, 305, 306, 308, 310, 313, 316, 316, 317, 318, 319, 321, 322, 322, 326, 331, 331, 331, 331, 333, 334, 335, 339, 343, 345, 346, 350, 354, 357, 358, 359, 360, 361, 366, 368, 372, 376, 379, 383, 384, 386, 388, 389, 391, 392, 394, 396, 399, 403, 406, 407, 412, 413, 414, 418, 419, 420 };
double Size_0[MAX_REQUESTS] = { 5392, 4185, 3761, 4575, 6146, 5203, 5059, 5936, 4989, 3726, 4350, 7209, 4649, 4028, 5331, 5195, 4652, 5672, 5370, 3766, 5088, 6945, 5424, 6482, 4676, 3481, 5182, 6752, 5416, 4023, 4222, 6310, 4395, 5024, 5095, 4803, 5009, 5131, 5394, 6111, 3901, 4606, 4833, 5178, 5320, 5002, 3439, 5171, 4727, 4885, 6590, 4534, 5596, 3707, 5970, 4761, 5791, 4063, 3479, 5341, 4838, 5630, 4237, 5932, 4180, 2965, 4452, 3744, 4366, 5515, 5201, 5484, 4825, 4793, 5481, 5429, 4167, 5420, 3913, 5758, 4278, 5497, 6333, 6835, 5414, 4699, 4327, 5089, 5263, 4027, 4024, 5228, 6236, 6622, 2520, 3805, 6313, 4588, 5255, 5287, 3349, 4366, 5706, 4899, 5296, 3743, 5535, 1467, 3001, 4792, 4608, 5929, 4172, 4642, 2716, 5530, 4537, 4618, 4478, 5664, 4805, 5093, 3722, 6792, 4849, 4025, 4978, 4600, 4084, 6397, 2357, 4505, 2701, 4384, 3254, 3737, 6092, 5571, 5085, 6963, 4023, 5900, 4665, 4414, 5472, 4044, 4881, 3742, 5658, 4773, 5604, 5571, 5708, 4510, 5737, 4691, 5222, 3848, 5881, 5085, 5148, 6864, 5496, 4530, 3931, 6245, 6189, 5349, 5346, 4035, 5801, 6162, 5136, 6382, 4685, 6043, 3738, 4698, 4716, 4896, 3038, 5949, 5098, 5921, 4134, 3988, 4890, 4140, 5437, 6462, 5014, 3827, 4376, 6255, 7339, 7072, 5219, 6756, 4469, 6112 };
double Dead_0[MAX_REQUESTS] = { 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000, 1000000 };
double Pri_0[MAX_REQUESTS] = { 2, 8, 15, 1, 10, 5, 19, 19, 3, 5, 6, 6, 2, 8, 2, 12, 16, 3, 8, 17, 12, 5, 3, 14, 13, 3, 2, 17, 19, 16, 8, 7, 12, 19, 10, 13, 8, 20, 16, 15, 4, 12, 3, 14, 14, 5, 2, 12, 14, 9, 8, 5, 3, 18, 18, 20, 4, 2, 10, 19, 17, 16, 11, 3, 9, 7, 1, 3, 5, 9, 7, 6, 11, 10, 11, 11, 7, 2, 14, 9, 10, 4, 5, 15, 17, 1, 7, 17, 12, 9, 5, 20, 7, 4, 18, 19, 19, 3, 10, 2, 14, 16, 20, 19, 5, 11, 18, 7, 14, 7, 2, 6, 5, 13, 11, 10, 18, 14, 18, 13, 7, 11, 2, 17, 16, 8, 16, 15, 12, 13, 11, 11, 2, 5, 7, 11, 8, 12, 8, 18, 18, 8, 14, 4, 6, 10, 10, 19, 2, 9, 3, 7, 7, 11, 14, 9, 1, 12, 3, 16, 11, 20, 5, 18, 9, 4, 16, 2, 3, 11, 12, 17, 15, 1, 17, 2, 9, 20, 9, 5, 2, 15, 14, 20, 19, 19, 1, 9, 8, 8, 9, 14, 9, 4, 8, 2, 11, 18, 14, 15 };

//----------------------------------------------------------------------------

double Spd_0[MAX_PROCESSORS] = { 1000 , 1000, 1000, 1000, 1000 };						//			//packets/s
double Delta_0[MAX_PROCESSORS] = { 0.00005 , 0.00005, 0.00005, 0.00005, 0.00005 };		//250 ms				//s/packet

//----------------------------------------------------------------------------

void define_prob(int num_requests, int num_processors)
{
	for (int j = 0; j < num_requests; j++)
	{
		req_set[j].init_time	= Init_0[j];
		req_set[j].size			= Size_0[j];
		req_set[j].deadline		= Dead_0[j];
		req_set[j].priority		= Pri_0[j];
		req_set[j].sTime = 0;
		req_set[j].pTime = 0;
		req_set[j].tTime = 0;
		req_set[j].ltncy = 0;
	}

	for (int i = 0; i < num_processors; i++)
	{
		proc_set[i].speed = Spd_0[i];
		proc_set[i].delta = Delta_0[i];
	}



	double pri_sum = 0;
	for (int j = 0; j < num_requests; j++)
		pri_sum += req_set[j].priority;


	for (int j = 0; j < num_requests; j++)
		req_set[j].priority = (double)req_set[j].priority / pri_sum;

}

void initialize_pop(chrom * pop, int pop_size, int num_requests, int num_processors)
{

	//Filling zeros
	for (int i = 0; i<pop_size; i++)
		for (int j = 0; j < (num_requests*num_processors); j++)
		{
			pop[i].bit[j] = 0;
			pop[i].theta[j] = -1;
		}


	//Allocate, check feasibility and calculate fitness
	for (int i = 0; i < pop_size; i++)
	{
		for (int j = 0; j < num_requests; j++)
		{
			pop[i].bit[(j*num_processors) + (rand() % num_processors)] = 1;
		}
	}


	calc_fittness(pop, pop_size, num_requests, num_processors);
}

void calc_fittness(chrom * pop, int pop_size, int num_requests, int num_processors)
{
	double fit_sum = 0;
	for (int i = 0; i<pop_size; i++)
	{
		//cout << "latency Good of " << i << "\n";
		pop[i].avg_ltncy = latency(pop, i, num_requests, num_processors);
		pop[i].fit = (1.0 / pop[i].avg_ltncy);
		pop[i].fsb = check_fsb(pop, i, num_requests, num_processors);
		//cout << pop[i].fsb << "  ";
		fit_sum += pop[i].fit;
	}
	//exit(0);
	//Normalized fittness
	for (int i = 0; i < pop_size; i++)
		pop[i].norm_fit = (pop[i].fit / fit_sum);

	fit_sum = 0;
	for (int i = 0; i < pop_size; i++)
	{
	fit_sum+=pop[i].fit;
	}cout << fit_sum << endl; exit(0);
}

int check_fsb(chrom * pop, int ind, int num_requests, int num_processors)
{
	int sum;

	//A request must be served and served only one time
	for (int i = 0; i < num_requests; i++)
	{
		sum = 0;
		for (int j = 0; j < num_processors; j++)
			sum += pop[ind].bit[(i*num_processors) + j];
		if (sum != 1)
		{
			//cout << "\nInfeasible-1\n";
			return 0;
		}
	}

	//Satisfying deadline requirements
	for (int j = 0; j < num_requests; j++)
	{
		//cout << req_set[j].ltncy << " <?> " << req_set[j].deadline << endl;
		if (req_set[j].ltncy > req_set[j].deadline   &&     req_set[j].deadline!=0)
		{
			//cout << "\nInfeasible-2\n";
			return 0;
		}
	}


	return 1;
}

double latency(chrom * pop, int ind, int num_requests, int num_processors)
{
	//for (int i = 0; i < num_processors; i++)
	//{
	//	for (int j = 0; j < num_requests; j++)
	//	{
	//		cout << (pop[ind].bit[(j*num_processors) + i]) << "  ";
	//	}
	//	cout << endl;
	//}
	double last_camein = req_set[0].init_time;
	for (int j = 0; j < num_requests; j++)
	{
		if (req_set[j].init_time > last_camein)
		{
			last_camein = req_set[j].init_time + req_set[j].size*FOG_EDGE_DELTA;;
		}
	}
	//cout << "last_camein = " << last_camein << endl;
	//cout <<"\nBeginning\n";
	//for (int j = 0; j < num_requests; j++)
	//{
	//cout << "sTime[ "<<j<<" ] = " << req_set[j].sTime << "\tpTime = " << req_set[j].pTime<<endl;
	//}
	//cout << "\n======================================================================\n";

	//A request cannot start processing before initiation time and transmission time 
	for (int j = 0; j < num_requests; j++)
	{
		for (int i = 0; i < num_processors; i++)
		{
			if (pop[ind].bit[(j*num_processors) + i] == 1)
			{

				req_set[j].tTime = proc_set[i].delta * req_set[j].size;
				req_set[j].sTime = req_set[j].init_time + req_set[j].tTime;
				//cout <<"sssss   "<< req_set[j].sTime << endl;
				req_set[j].pTime =  req_set[j].size / proc_set[i].speed;

				break;
			}
		}
	}
	//cout << endl << endl;
	//cout << "\n*******************************************";
	//cout << endl << "One request at a time\n";
	//for (int j = 0; j < num_requests; j++)
	//{
	//	cout << "sTime[ " << j << " ] = " << req_set[j].sTime << "\tpTime = " << req_set[j].pTime << "\ttTime = " << req_set[j].tTime << "\tltncy = " << req_set[j].ltncy << "\tinit_time = " << req_set[j].init_time << endl;
	//}
	//cout << "\n======================================================================\n";

	//A processor can process only one request at a time
	int * indx_arr = new int[num_requests];
	int * indx_arr2 = new int[num_requests];
	int tm;
	int * arrive_dur = new int[num_requests];

	for (int i = 0; i < num_processors*num_requests; i++)
		pop[ind].theta[i] = -1;

	for (int i = 0; i < num_processors; i++)
	{
		int tht_ind = i * num_requests;
		int req_cnt = 0;
		int Indx = 0;
		for (int j = 0; j < num_requests; j++)
			indx_arr[j] = -1;


		//cout << "\n*******************************************";
		//cout << "\nP[" << i << "]: ";
		for (int j = 0; j < num_requests; j++)
		{
			//cout << x.bit[i + (j*num_processors)];
			if (pop[ind].bit[i + (j*num_processors)])
			{
				req_cnt++;
				indx_arr[Indx] = j;
				//indx_arr2[Indx] = j;
				Indx++;
			}
		}

		float tF_next_at = 0 ;
		float p_next_at = 0;
		float tb_next_at = 0, sTime_next=0;
		int cnt_back = req_cnt;

		for (int j = 0; j < req_cnt; j++)
		{
			int r = rand() % cnt_back;
			int J = indx_arr[r];
			pop[ind].theta[tht_ind++] = J;
			//cout << "J= " << J<<endl;
			//Execute
			if (pop[ind].bit[i + (J*num_processors)])
			{
				/*
				if (req_set[J].sTime < tF_next_at)
					req_set[J].sTime = tF_next_at;
				req_set[J].sTime = req_set[J].sTime + req_set[J].tTime;
				tF_next_at = req_set[J].sTime;
				//------------------------------------------------------
				if (req_set[J].sTime < p_next_at)
					req_set[J].sTime = p_next_at;
				req_set[J].sTime = req_set[J].sTime + req_set[J].pTime;
				p_next_at = req_set[J].sTime;
				//------------------------------------------------------
				if (req_set[J].sTime < tb_next_at)
					req_set[J].sTime = tb_next_at;
				req_set[J].sTime = req_set[J].sTime + req_set[J].tTime;
				tb_next_at = req_set[J].sTime;*/
				

				if (req_set[J].sTime < sTime_next)
					req_set[J].sTime = sTime_next;

				sTime_next = req_set[J].sTime + req_set[J].pTime;

				req_set[J].ltncy = req_set[J].sTime + req_set[J].pTime + req_set[J].tTime - req_set[J].init_time;

				/*cout << "req_set[J].init_time = " << req_set[J].init_time << endl;
				cout << "req_set[J].sTime = " << req_set[J].sTime << endl;
				cout << "req_set[J].tTime = " << req_set[J].tTime << endl;
				cout << "req_set[J].size = " << req_set[J].size << endl;
				cout << "req_set[J].pTime = " << req_set[J].pTime << endl;
				cout << "req_set[J].ltncy = " << req_set[J].ltncy << endl;
				cout << "sTime_next = " << sTime_next << endl;
				cout << endl<<endl;*/
				//exit(0);


			}
			
			//pack up
			for (int z = 0; z < num_requests; z++)
				indx_arr2[z] = -1;

			int Indx = 0;
			for (int z = 0; z < num_requests; z++)
			{
				if (indx_arr[z] != J)
					indx_arr2[Indx++] = indx_arr[z];
			}

			//cout << endl << "+++++++++" << endl;
			for (int z = 0; z < num_requests; z++)
			{
				indx_arr[z] = indx_arr2[z];
				//cout << indx_arr[z] << "  ";
			}
			cnt_back--;

			//------------------------------------------------------

			//cout << "\n&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
		}
		//exit(0);
		//cout << "\n&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
	}

	for (int j = 0; j < num_requests; j++)
	{
		pop[ind].tTime[j] = req_set[j].tTime;
		pop[ind].pTime[j] = req_set[j].pTime;
		pop[ind].ltncy[j] = req_set[j].ltncy;
	}

	//Calculate the latency cost
	double latency = 0;

	for (int j = 0; j < num_requests; j++)
	{
		latency += (req_set[j].ltncy*req_set[j].priority);
		//cout << "req_set[j].ltncy: " << req_set[j].ltncy << endl;
		//latency += req_set[j].sTime;
	}
	//cout << "latency: " << latency << endl;
	//exit(0);
	//latency = latency / 1000.0;
	//cout << "\nlatency  " << latency<<endl;
	//cout << "-------------------------";

	delete[] indx_arr;
	delete[] indx_arr2;
	delete[] arrive_dur;

	return latency;
}

void display_pop(chrom * pop, int pop_size, int num_requests, int num_processors)
{
	double sum_fit = 0;
	for (int i = 0; i < pop_size; i++)
	{
		cout << "\nchrom[" << i << "]  ";
		if (i < 10)cout << " ";
		cout<<"=  ";
		for (int j = 0; j < (num_requests*num_processors); j++)
		{
			cout << pop[i].bit[j];
			if (((j + 1) % num_processors) == 0) cout << " ";
		}
		cout << endl;
		for (int j = 0; j < (num_requests*num_processors); j++)
		{
		if(pop[i].theta[j] !=-1)	cout << pop[i].theta[j]+1<<", ";
		else cout << "-, ";
		if (((j+1) % num_requests) == 0) cout << endl;
		}
		cout << "fsb = " << pop[i].fsb;
		cout << "\t\tfit = " << pop[i].fit;
		cout << "\t\tnorm_fit = " << pop[i].norm_fit;
		cout << "\t\tlatency = " << pop[i].avg_ltncy;
		//cout << setprecision(10) << "\t\tfit = " << pop[i].fit;
		//cout << setprecision(4) << "\t\tlatency = " << pop[i].avg_ltncy;
		//cout << endl;
		sum_fit += pop[i].norm_fit;
	}
	//cout << "\nSum of norm_fit = " << sum_fit << endl;
	cout << "\nFeasibility = " << fsb_percentage(pop, pop_size) * 100 << " %" << endl;
}

void sort(chrom * pop, int pop_size)
{
	chrom * temp = new chrom[1];

	for (int i = 0; i < pop_size; i++)
	{
		for (int j = 0; j < pop_size - 1; j++)
		{
			if (pop[j + 1].fit > pop[j].fit)
			{
				temp[0] = pop[j + 1];
				pop[j + 1] = pop[j];
				pop[j] = temp[0];
			}
		}
	}

	delete[] temp;
}

int select_chrom(chrom * pop, int pop_size)
{
	//cout << "select_chrom\n";
	double R = (rand() % 100000) / 100000.0;


	double fit_sum = 0;
	int i = 0;
	for (i = 0; i < pop_size; i++)
	{
		fit_sum += pop[i].norm_fit;
		//cout << "fit_sum = "<<fit_sum << endl;
		if ((fit_sum >= R) && (pop[i].fsb))
		{
			//cout << "Elitist\n";
			return i;
		}
	}
	//cout << "Random\n";
	return (rand() % pop_size);
}

void display_chrom(chrom * chrom, int num_requests, int num_processors)
{
	for (int j = 0; j < (num_requests*num_processors); j++)
	{
		cout << chrom[0].bit[j];
		if (((j + 1) % num_processors) == 0) cout << "  ";
	}

	cout << endl;
	for (int j = 0; j < (num_requests*num_processors); j++)
	{
		if (chrom[0].theta[j] != -1)
			cout << chrom[0].theta[j] + 1 << " ";
		else cout << "0 ";
		if (((j + 1) % num_requests) == 0 && num_requests>1) cout << endl;
	}
	if (num_requests==1) cout << endl;

	cout << "fsb = " << chrom[0].fsb;
	cout << "\t\tfit = " << chrom[0].fit;
	cout << "\t\tnorm_fit = " << chrom[0].norm_fit;
	cout << "\t\tlatency = " << chrom[0].avg_ltncy;
}

void cross(chrom * chrom1, chrom * chrom2, int num_requests, int num_processors)
{
	int cross_point = (rand() % CROSS_POINT) + 1;    //-------------------------------->>>> Tuning

	//cout << "\nCrossOver:";
	
	chrom * temp = new chrom[1];
	temp[0] = chrom1[0];

	for (int j = 0; j < cross_point; j++)
	{
		for (int i = 0; i < num_processors; i++)
		{
			chrom1[0].bit[(j*num_processors) + i] = chrom2[0].bit[(j*num_processors) + i];
			chrom2[0].bit[(j*num_processors) + i] = temp[0].bit[(j*num_processors) + i];
		}
	}
	
	delete[] temp;
}

void crossover(chrom * pop, int pop_size, int num_requests, int num_processors)
{
	//cout << "\ncrossover1\n" << pop_size;

	chrom * temp_pop = new chrom[pop_size];
	chrom * chrom1 = new chrom[1];
	chrom * chrom2 = new chrom[1];
	

	//cout << "crossover2\n";

	int index = 0;
	//cout << endl;
	for (int i = 0; i < pop_size / 2; i++)
	{
		int x = select_chrom(pop, pop_size);
		int y = select_chrom(pop, pop_size);

		//cout << "\n============================================================\n";
		//cout << "\ncrosssing number " << i << endl;
		//cout << x << " with " << y << endl << endl;

		chrom1[0] = pop[x];
		chrom2[0] = pop[y];

		//display_chrom(chrom1, num_requests, num_processors);
		//display_chrom(chrom2, num_requests, num_processors);
		
		cross(chrom1, chrom2, num_requests, num_processors);
		
		//display_chrom(chrom1, num_requests, num_processors);
		//display_chrom(chrom2, num_requests, num_processors);

		temp_pop[index++] = chrom1[0];
		temp_pop[index++] = chrom2[0];
	

	}
	//display_pop(temp_pop, pop_size, num_requests, num_processors);
	
	for (int i = 0; i < pop_size; i++)
		pop[i] = temp_pop[i];

	calc_fittness(pop, pop_size, num_requests, num_processors);

	delete[] temp_pop;
	delete[] chrom1;
	delete[] chrom2;
}

void mutation(chrom * pop, int pop_size, int num_requests, int num_processors)
{
	for (int i = 0; i < pop_size; i++)
	{
		if ((rand() % 4) == 0)//25				//Probability of mutation for each chromosome is 10%
		{
			for (int w = 0; w < MUT_NUM; w++)
			{
				int req_indx = rand() % num_requests;
				int prc_indx = rand() % num_processors;

				for (int c = 0; c < num_processors; c++)
					pop[i].bit[(req_indx*num_processors) + c] = 0;

				pop[i].bit[(req_indx*num_processors) + prc_indx] = 1;
				//cout << "\nMutate ["<<i<<"]:  i = " << prc_indx << " j = " << req_indx ;
			}

		}
	}

	calc_fittness(pop, pop_size, num_requests, num_processors);
}

double fsb_percentage(chrom * pop, int pop_size)
{
	int cnt = 0;

	for (int i = 0; i < pop_size; i++)
		if (pop[i].fsb == 1) cnt++;

	//return ((double)cnt / (double)pop_size);
	return cnt;

}

void main()
{
	srand((unsigned)time(NULL));
	clock_t t;
	t = clock();
	int mostucky = 0;
	double best_sofar;
	double prev_best;


	int num = MAX_NUM, num_requests = MAX_REQUESTS, num_processors = MAX_PROCESSORS; 


	req_set = new request[num_requests];
	proc_set = new processor[num_processors];
	Beta = new int[num_requests*num_requests];


	define_prob(num_requests, num_processors);


	chrom * popcurrent = new chrom[pop_size];
	chrom * popnext = new chrom[pop_size];
	chrom * best_sln = new chrom[1];

	//cout << "Initial population:";
	initialize_pop(popcurrent, pop_size, num_requests, num_processors);
	//display_pop(popcurrent, pop_size, num_requests, num_processors); 

	prev_best   = popcurrent[0].avg_ltncy;
	best_sofar  = popcurrent[0].avg_ltncy;
	best_sln[0] = popcurrent[0];

	//if (fsb_percentage(popcurrent, pop_size) == 0)
	//{
	//	cout << "\nTermination\n";
	//	exit(0);
	//}
	cout << popcurrent[0].avg_ltncy << endl;
	int i, updated;
	for (i = 1; i <= num; i++)
	{
		//cout << "\ni = " << i << "\n";

		//cout << "best_sofar = " << best_sofar << endl;
		
		for (int j = 0; j < pop_size; j++)
			popnext[j] = popcurrent[j];            	//copy popcurrent to popnext in order to adjust it

		//cout << "\nSorting:";
		sort(popnext, pop_size);
		cout << popnext[0].avg_ltncy << " \t " << popnext[0].fsb << endl;
		//display_pop(popnext, pop_size, num_requests, num_processors);
		//cout << endl;

		//-----------------------
		updated = 0 ;
		for (int j = 0; j < pop_size; j++)
		{
			if ((popnext[j].avg_ltncy < best_sofar) && (popnext[j].fsb))
			{
				best_sln[0] = popnext[j];
				best_sofar = popnext[j].avg_ltncy;
				mostucky = 0;
				//cout << "Better :)\n";
				//cout << "Latency: " << best_sofar << " ,FSB:" << best_sln[0].fsb<<endl;
				updated = 1;
			}
		}
		if (!updated)		
		{
			mostucky++;
			//cout << "No better (" << mostucky << ") :(\n";
			if (mostucky == CC)
				break;
		}
		/*cout << "avg_ltncy = " << popnext[0].avg_ltncy << endl;
		if (popnext[0].avg_ltncy < best_sofar && (popnext[0].fsb))
		{
			best_sln[0] = popnext[0];
			best_sofar = popnext[0].avg_ltncy;
			mostucky = 0;
			cout << "Better :)\n";
		}
		else
		{
			mostucky++;
			cout << "Bad (" << mostucky << ") :(\n";
			if (mostucky == CC)
				break;
		}*/
		//prev_best = popnext[0].avg_ltncy;
		//-----------------------
		
		if (num_requests>1)
			crossover(popnext, pop_size, num_requests, num_processors);  
		
		//display_pop(popnext, pop_size, num_requests, num_processors);
		//cout << endl;


		//cout << "\nMutation:";
		mutation(popnext, pop_size, num_requests, num_processors);
		//display_pop(popnext, pop_size, num_requests, num_processors);
		//cout << endl;

		for (int j = 0; j < pop_size; j++)
			popcurrent[j] = popnext[j];

		//cout << "\n----------------------------------------------------------------------------------------------------\n\n";
	}


	cout << "\n\n==============================================================================================================\n\n";

	cout << "Best solution: " << endl;
	//display_chrom(best_sln, num_requests, num_processors);
	cout << endl;

	t = clock() - t;
	//cout << "\n**************************************************" << endl;
	cout << best_sofar << "\n";// << ((double)t / CLOCKS_PER_SEC) << endl;
	cout << endl << endl;
	delete[] popcurrent;
	delete[] popnext;
	delete[] best_sln;

	delete[] req_set;
	delete[] proc_set;
	delete[] Beta;

}

