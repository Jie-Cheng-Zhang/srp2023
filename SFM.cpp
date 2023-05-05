// SFM.cpp : Defines the entry point for the console application.
//
#include "iostream"
#include "fstream"
#include "sstream"
#include "stdafx.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include <algorithm>
#include <vector>
#include <random>
#include <list>
using namespace std;
#define MAXOBSLINES	1000
#define	MAXAGENTS	1000
double	tao = 0.5;
double	v0 = 1.02;
double	A = 1000;
double	B = 0.08;
double	k1 = 120000; //120000
double	k2 = 240000;
double	maxv = 3;
double	c_mass = 320.; // m/cmass = r-->0.25 ==> cmass = m/0.25 = 320;
double	sense_range = 10;
double  col_range=382;
double  row_range=679;
int  obstical_line_num=0;
double  map_factor=10;
int num_ticks=100000;
double tick = 0.02;
int plan_range=5;
int A_step=10000;
double arrive_range=0.3;
vector<vector<bool>> map_matrix;

struct node{
	bool flag;
	double h;
	double g;
	int x;
	int y;
	node* parent;

	node()=default;
	node(int a, int b, double c=0){
		flag=0;
		x=a;
		y=b;
		h=c;
	}
	static bool compare(node& a, node& b){
		if ((a.g+a.h)>(b.g+b.h)){
			return true;
		}
		return false;
	}
};

vector<vector<node>> map_matrix_A;

struct cordinate{
	double x;
	double y;
	cordinate()=default;
	cordinate(double a, double b){
		x=a;
		y=b;
	}
};

struct AGENT
{
	double	m;
	double	x;
	double	y;
	double	vx;
	double	vy;
	double	gx;
	double	gy;
	double	lgx;
	double	lgy;
	double	v0;
	int		id;
	list <cordinate> path; 
};
struct AGENT agents[MAXAGENTS];
//==========================================================================================================
static int phase = 0;
double gaussian()                                                       
{
    static double V1, V2, S; 
    double X;    
    if ( phase == 0 ) {
        do {       
            double U1 = (double)rand() /(double)RAND_MAX;
            double U2 = (double)rand() /(double)RAND_MAX;
            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);        
        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);
        
    phase = 1 - phase;
    return X;
}
double Point_to_line_distance(double x0, double y0, 
	double sx, double sy, double ex, double ey, double d, double* crx, double* cry)
{

	double t0 = ((ex - sx) * (x0 - sx) + (ey - sy) * (y0 - sy)) / (d * d);
	if(t0 < 0){
		d = sqrt((x0 - sx) * (x0 - sx) + (y0 - sy) * (y0 - sy));
	}else if(t0 > 1){
		d = sqrt((x0 - ex) * (x0 - ex) + (y0 - ey) * ( y0 - ey));
	}else{
		d = sqrt(
			(x0 - (sx + t0 * ( ex  - sx))) * (x0 - (sx + t0 * ( ex  - sx))) +
			(y0 - (sy + t0 * ( ey  - sy))) * (y0 - (sy + t0 * ( ey  - sy)))
			);
	}
	*crx = sx + t0 * (ex - sx);
	*cry = sy + t0 * (ey - sy);
	return d;
}

double randval(double a, double b)
{
	return a + (b - a) * rand() /(double)RAND_MAX;
}
struct OBSLINE
{
	double sx;
	double sy;
	double ey;
	double ex;
	double d;
};
struct OBSLINE obstical_lines[MAXOBSLINES];

void initialize() 
{   
	fstream infile("./map/obstacles3.txt", ios::in);
    string buf;
    int num = 0;
    getline(infile, buf);
    while (!infile.eof()) {
        stringstream ss(buf);
        ss >> obstical_lines[num].sy >> obstical_lines[num].sx >> obstical_lines[num].ey >> obstical_lines[num].ex;
		obstical_lines[num].sx/=map_factor;
		obstical_lines[num].sy/=map_factor;
		obstical_lines[num].ex/=map_factor;
		obstical_lines[num].ey/=map_factor;
		obstical_lines[num].d = sqrt((obstical_lines[num].sx - obstical_lines[num].ex) * (obstical_lines[num].sx - obstical_lines[num].ex) +
		(obstical_lines[num].sy - obstical_lines[num].ey) * (obstical_lines[num].sy - obstical_lines[num].ey));
		num++;
		getline(infile,buf);
	}
	obstical_line_num=num;
	infile.close();

	fstream infile2("./map/matrix3.txt", ios::in);
	buf="";
	getline(infile2, buf);
	while(!infile2.eof()){
		stringstream ss(buf);
		map_matrix.push_back(vector<bool>());
		int temp;
		while (!ss.eof())
		{
			ss>>temp;
			map_matrix[map_matrix.size()-1].push_back(temp);
		}
		getline(infile2, buf);
	}
	infile2.close();

	for (int i = 0; i < MAXAGENTS; i++) {
		agents[i].id = i;
		do{
			agents[i].x = randval(0, 1) * col_range/map_factor;
			agents[i].y = randval(0,1) * row_range/map_factor;
		}while(!map_matrix[int(agents[i].y*map_factor)][int(agents[i].x*map_factor)]);
		do{
			agents[i].gx = 15;
			agents[i].gy = 5;
		}while(!map_matrix[int(agents[i].gy*map_factor)][int(agents[i].gx*map_factor)]);
		agents[i].v0 = v0;
		agents[i].vx = randval(-2, 2);
		agents[i].vy = randval(-2, 2);
		agents[i].m = 80;        			
	}

	for (int i=0; i<row_range; i++){
		map_matrix_A.push_back(vector<node>());
		for (int j=0; j<col_range; j++){
			map_matrix_A[i].push_back(node(j,i));
		}
	}

}


void compute_agent_force(AGENT*temp, AGENT*temp2, double *fx, double *fy)
{
	double d = sqrt((temp->x - temp2->x) * (temp->x - temp2->x) + (temp->y - temp2->y) * (temp->y - temp2->y));
	if(d == 0) {printf("here d == 0 error but fixed..."); d = 1e-10;}
	double delta_d = temp->m / c_mass + temp2->m / c_mass - d;
	double fexp = A * exp(delta_d / B);
	double fkg = delta_d < 0? 0: k1 * delta_d;
	double nijx = (temp->x - temp2->x) / d;
	double nijy = (temp->y - temp2->y) / d;
	double fnijx = (fexp + fkg) * nijx;
	double fnijy = (fexp + fkg) * nijy;
	double fkgx = 0;
	double fkgy = 0;
	if(delta_d>0){
		double tix = -nijy;
		double tiy = nijx;
		fkgx = k2 * delta_d;
		fkgy = k2 * delta_d;
		double delta_vij = (temp2->vx - temp->vx) * tix + (temp2->vy - temp->vy) * tiy;
		fkgx = fkgx * delta_vij * tix;
		fkgy = fkgy * delta_vij * tiy; 
	}			
	*fx += fnijx + fkgx;
	*fy += fnijy + fkgy;
}

int Intersection_2LineSegmen(double p0_x, double p0_y, double p1_x, double p1_y, double p2_x, double p2_y, double p3_x, double p3_y, double* i_x, double* i_y)
{
    double s1_x, s1_y, s2_x, s2_y;
    s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
    s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;

    double s, t;
    s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
    t = ( s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);

    if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
    {
        // Collision detected
        if (i_x != NULL)
            *i_x = p0_x + (t * s1_x);
        if (i_y != NULL)
            *i_y = p0_y + (t * s1_y);
        return 1;
    }
    return 0; // No collision
}

struct cross_info
{
	double diw;
	double crx;
	double cry;
	static bool compare(cross_info& a, cross_info& b){
		return a.diw<b.diw;
	}
};

vector<cordinate> direction={{1,0}, {-1,0}, {0,1}, {0,-1}};

bool check_in_map(int x, int y){
	if (x<0 || x>=col_range){
		return false;
	}
	if (y<0 || y>=row_range){
		return false;
	}
	return true;
}

void A_algorithm(){

	for (int i=0; i<MAXAGENTS; i++){
		list<node*> temp_list;
		vector<node*> delete_list;
		temp_list.push_back(&map_matrix_A[int(agents[i].y*map_factor)][int(agents[i].x*map_factor)]);
		delete_list.push_back(&map_matrix_A[int(agents[i].y*map_factor)][int(agents[i].x*map_factor)]);
		map_matrix_A[int(agents[i].y*map_factor)][int(agents[i].x*map_factor)].flag=1;
		map_matrix_A[int(agents[i].y*map_factor)][int(agents[i].x*map_factor)].h=0;
		map_matrix_A[int(agents[i].y*map_factor)][int(agents[i].x*map_factor)].parent=nullptr;
		while (temp_list.size())
		{
			node* temp_node= temp_list.front();
			temp_list.pop_front();
			if(temp_node->x==int(agents[i].gx*map_factor) && temp_node->y==int(agents[i].gy*map_factor)){
				agents[i].path.clear();
				int count=0;
				while(temp_node->parent!=nullptr){
					if(count%plan_range==0){
						agents[i].path.push_front(cordinate(temp_node->x/map_factor,temp_node->y/map_factor));
					}
					temp_node=temp_node->parent;
					count++;
				}
				while (delete_list.size())
				{
					delete_list.back()->flag=0;
					delete_list.pop_back();
				}
				AGENT* temp_agent=&agents[i];
				break;
			}
			for (auto j: direction){
				if (check_in_map(temp_node->x+j.x, temp_node->y+j.y) && !map_matrix_A[temp_node->y+j.y][temp_node->x+j.x].flag){
					if(!map_matrix[temp_node->y+j.y][temp_node->x+j.x]){
						continue;
					}
					bool temp_flag=false;
					for(auto k: direction){
						if(check_in_map(temp_node->x+j.x+k.x, temp_node->y+j.y+k.y)&&!map_matrix[temp_node->y+j.y+k.y][temp_node->x+j.x+k.x]){
							temp_flag=true;
							break;
						}
						if(check_in_map(temp_node->x+j.x+2*k.x, temp_node->y+j.y+2*k.y)&&!map_matrix[temp_node->y+j.y+2*k.y][temp_node->x+j.x+2*k.x]){
							temp_flag=true;
							break;
						}
					}
					if(temp_flag){
						continue;
					}
					map_matrix_A[temp_node->y+j.y][temp_node->x+j.x].h=temp_node->h+1;
					map_matrix_A[temp_node->y+j.y][temp_node->x+j.x].g=sqrt(pow(map_matrix_A[temp_node->y+j.y][temp_node->x+j.x].x-agents[i].gx*map_factor,2)+pow(map_matrix_A[temp_node->y+j.y][temp_node->x+j.x].y-agents[i].gy*map_factor,2));
					map_matrix_A[temp_node->y+j.y][temp_node->x+j.x].flag=1;
					map_matrix_A[temp_node->y+j.y][temp_node->x+j.x].parent=temp_node;
					if(!temp_list.size()){
						temp_list.push_front(&map_matrix_A[temp_node->y+j.y][temp_node->x+j.x]);
					}
					else{
						auto k=temp_list.begin();
						for(k; k!=temp_list.end() && node::compare(map_matrix_A[temp_node->y+j.y][temp_node->x+j.x],**k); k++);
						temp_list.insert(k, &map_matrix_A[temp_node->y+j.y][temp_node->x+j.x]);
						delete_list.push_back(&map_matrix_A[temp_node->y+j.y][temp_node->x+j.x]);
					}
				}
			}
		}
	}
}

void check_next_tg(){
	for (int i=0; i<MAXAGENTS; i++){
		if((agents[i].path.front().x-agents[i].x<arrive_range && agents[i].path.front().x-agents[i].x>-arrive_range) && (agents[i].path.front().y-agents[i].y<arrive_range && agents[i].path.front().y-agents[i].y>-arrive_range)){
			if (agents[i].path.size()>1){
				agents[i].path.pop_front();
			}
		}
	}
}

void step() {

	int i, j, k;
	double dx, dy, d;
	double wfx, wfy;
	double pfx, pfy;
	

	for (i = 0; i < MAXAGENTS; i++){	
		AGENT *temp = &agents[i];
		double dvtx = 0, dvty = 0;
		//first compute the desired direction;
		double x0 = temp->x;
		double y0 = temp->y;
		double d0 = sqrt((x0 -temp->path.front().x) * (x0 - temp->path.front().x) + (y0 - temp->path.front().y) * (y0 - temp->path.front().y));
		dx = temp->v0 *(temp->path.front().x - x0) / d0;
		dy = temp->v0 *(temp->path.front().y - y0) / d0;
		dvtx = (dx - temp->vx) / tao;
		dvty = (dy - temp->vy) / tao;

		//compute the force of other agents;
		double total_fx = 0;
		double total_fy = 0;

		for (j = 0; j < MAXAGENTS; j++){
			if(j == i) continue;			
			AGENT *temp2 = &agents[j];
			double dis = (temp->x - temp2->x) *  (temp->x - temp2->x) + (temp->y - temp2->y) *  (temp->y - temp2->y);
			if(dis > sense_range * sense_range) continue;
			compute_agent_force(temp, temp2, &total_fx, &total_fy);
		}

		// compute the wall....;

		vector <cross_info> cross_info_list;
		for(j = 0; j < obstical_line_num; j++){
			cross_info temp;
			temp.diw = Point_to_line_distance(x0, y0, obstical_lines[j].sx, obstical_lines[j].sy, obstical_lines[j].ex, obstical_lines[j].ey, obstical_lines[j].d, &temp.crx, &temp.cry);
			cross_info_list.push_back(temp);
		}
		
		sort(cross_info_list.begin(),cross_info_list.end(),cross_info::compare);
		for (auto i=cross_info_list.begin(); i != cross_info_list.end(); i++){
			if(i->diw>sense_range){
				cross_info_list.erase(i, cross_info_list.end());
				break;
			}
		}

		for (j=0; j < cross_info_list.size();j++){
			double vir_diw = sqrt((x0 - cross_info_list[j].crx) * (x0 - cross_info_list[j].crx) + (y0 - cross_info_list[j].cry) * (y0 - cross_info_list[j].cry));
			if(vir_diw == 0) {printf("here vir_diw == 0 error but fixed..."); vir_diw = 1e-10;}
			double niwx = (x0 - cross_info_list[j].crx)/vir_diw;
			double niwy = (y0 - cross_info_list[j].cry)/vir_diw;
			double drw = temp->m / c_mass - cross_info_list[j].diw;
			double fiw_1 = A * exp(drw/B);
				
			if(drw > 0){
				fiw_1 += k1 * drw;
			}
			

			double fniwx = fiw_1 * niwx;
			double fniwy = fiw_1 * niwy;

			//---the second_part.
			double fiw_kgx = 0;
			double fiw_kgy = 0;
			
			if(drw > 0){
				double fiw_kg = k2 * drw * (temp->vx * (-niwy) + temp->vy * niwx);
				fiw_kgx = fiw_kg * (-niwy);
				fiw_kgy = fiw_kg * (niwx); 
			}
			total_fx += fniwx - fiw_kgx;
			total_fy += fniwy - fiw_kgy;
		}

		//===================================
		dvtx = dvtx + total_fx / temp->m;
		dvty = dvty + total_fy / temp->m; 

		temp->vx = temp->vx + dvtx * tick;
		temp->vy = temp->vy + dvty * tick;

		double dv = sqrt(temp->vx * temp->vx + temp->vy * temp->vy);
		if(dv > maxv){
			temp->vx = temp->vx * maxv / dv;
			temp->vy = temp->vy * maxv / dv;
		}

        /* double mint = 1;
		for(j = 0; j < obstical_line_num; j++){
			double crx, cry, tt;
                        int ret = Intersection_2LineSegmen(x0, y0, x0+0.5 * temp->vx * tick, y0 + 0.5 * temp->vy * tick, obstical_lines[j].sx, obstical_lines[j].sy, obstical_lines[j].ex, obstical_lines[j].ey, &crx, &cry);
                        if(ret == 1){
                             if(fabs(crx - x0) > 0)
				tt = (crx - x0) / (temp->vx * tick);
                             else tt = (cry - y0) / (temp->vy * tick + 1e-20);
                             if(tt < mint) mint = tt;
			}
		}
        
		temp->vx = temp->vx * mint;
        temp->vy = temp->vy * mint;
		*/

		temp->x = temp->x + temp->vx * tick;
		temp->y = temp->y + temp->vy * tick;        
		if(!check_in_map(int(temp->x*map_factor),int(temp->y*map_factor))){
			temp->vx=-temp->vx;
			temp->vy=-temp->vy;
			temp->x=temp->x+2*temp->vx*tick;
			temp->y=temp->y+2*temp->vy*tick;
		}       
	}	
}

int _tmain(int argc, _TCHAR* argv[])
{
	time_t now = time(NULL);
	tm* tm_t = localtime(&now);
	std::stringstream ss;
	ss  <<"_"<< tm_t->tm_year + 1900 << "-" << tm_t->tm_mon + 1 << "-" << tm_t->tm_mday
		<< "_" << tm_t->tm_hour << "-"<< tm_t->tm_min << "-"<< tm_t->tm_sec;
	string fileName = "./result/demo"+ss.str()+".txt";
	const char * fileNameChar = fileName.c_str();
	FILE *f = fopen(fileNameChar,"w");
	//fprintf(f,"%g,%g,%g,%g,%d\n", 1/tick, map_factor, col_range, row_range, MAXAGENTS);
	initialize();
	for(int i = 0; i <num_ticks; i++){
		if(i%A_step==0){
			A_algorithm();
		}
		check_next_tg();
		step();
		fprintf(f,"%d\n", MAXAGENTS);
		for(int j = 0; j < MAXAGENTS; j++){
			fprintf(f,"%g,%g,%g,%g,%g,%g\n", agents[j].x, agents[j].y, agents[j].vx, agents[j].vy,agents[j].path.front().x,agents[j].path.front().y);
			//fprintf(f,"%g,%g,%g,%g\n", agents[j].x, agents[j].y, agents[j].vx, agents[j].vy);
		}
	}
	fclose(f);

	return 0;
}

