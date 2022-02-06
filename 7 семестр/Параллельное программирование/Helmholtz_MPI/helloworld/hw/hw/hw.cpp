#include<iostream>
#include <fstream>
#include <cmath>
#include<mpi.h>

//решаем уравнение Гельмгольца вида -d^2(u)/dx^2-d^2(u)/dy^2+k^2*u=f с нулевыми НУ

using namespace std;

typedef double mytipe;  //тип данных, использующийся во всей программе

const mytipe EPS = 10e-8;//константа сравнения с 0
const mytipe Pi = 3.14159265358979323846264338328;

mytipe f(mytipe x, mytipe y, mytipe k)
{
	return 2.*sin(Pi*y) + k*k*(1. - x)*x*sin(Pi*y) + Pi*Pi*(1. - x)*x*sin(Pi*y);
}

mytipe u(mytipe x, mytipe y)
{
	return (1. - x)*x*sin(Pi*y);
}


void print(int N, mytipe** A); 
void print_file(int N, mytipe h, mytipe** A);

int main(int argc, char **argv)
{
	const int n_glob = 1024;//кол-во узлов
	const int N_glob = n_glob - 1;//кол-во отрезков
	mytipe k = 150.;
	mytipe h = 1. / mytipe(N_glob); 
	mytipe c = 1. / (k*k*h*h + 4.);
	mytipe d = h*h*c;
	const mytipe acc = 10e-4;//точность вычисления

	MPI_Init(&argc, &argv); 

	MPI_Status stat;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);//номер исполняемового процесса
	int nproc;
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);//кол-во исполняемых процессов

	int* nLoc = nullptr;//кол-во разбиений матрицы U для каждого процесса
	int n_loc;//Кол-во строк марицы U_loc для текущего процесса
	mytipe** U = nullptr;//Численное решение, найденное с заданной точностью

	int* sdvig = nullptr;
	int sdvig_proc;

	mytipe t_init = MPI_Wtime();

	if (rank == 0)
	{
		U = new mytipe*[n_glob];//решение
		for (int i = 0; i < n_glob; i++) U[i] = new mytipe[n_glob];

		nLoc = new int[nproc];
		for (int i = 0; i < nproc; i++) nLoc[i] = n_glob / nproc;//кол-во строк под каждый процесс
		for (int i = 0; i < n_glob  % nproc; i++) ++nLoc[i];

		sdvig = new int[nproc];
		sdvig[0] = 0;
		for (int i = 1; i < nproc; i++) sdvig[i] = sdvig[i - 1] + nLoc[i - 1];
	}


	MPI_Scatter(nLoc, 1, MPI_INT, &n_loc, 1, MPI_INT, 0, MPI_COMM_WORLD);//отправляем каждому процессу свое значение nLoc, записывая в переменную n_loc
	MPI_Scatter(sdvig, 1, MPI_INT, &sdvig_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int N_loc = n_loc - 1;

	mytipe** yk = new mytipe*[n_loc]; 
	mytipe** yk1 = new mytipe*[n_loc]; 
	mytipe** right = new mytipe*[n_loc];

	for (int i = 0; i < n_loc; i++) {
		yk[i] = new mytipe[n_glob];
		yk1[i] = new mytipe[n_glob];
		right[i] = new mytipe[n_glob];
		for (int j = 0; j < n_glob; j++) { right[i][j] = d*f((sdvig_proc+i)*h, j*h, k); yk[i][j] = 0.; yk1[i][j] = 0.; }//нач приближ
	}

	mytipe* y_up = new mytipe[n_glob];
	mytipe* y_down = new mytipe[n_glob];

	

	MPI_Request * reqOut;
	MPI_Request * reqIn;
	if (rank == 0 || rank == nproc - 1)
	{
		reqOut = new MPI_Request;
		reqIn = new MPI_Request;
	}
	else
	{
		reqOut = new MPI_Request[2];
		reqIn = new MPI_Request[2];
	}
	int nreqIn = 0;
	int nreqOut = 0;

	if (rank < nproc - 1)
	{
		MPI_Send_init(yk[N_loc], n_glob, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, reqOut + nreqOut);
		MPI_Recv_init(y_up, n_glob, MPI_DOUBLE, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, reqIn + nreqIn);
		nreqIn++;
		nreqOut++;
	}
	if (rank > 0)
	{
		MPI_Send_init(yk[0], n_glob, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, reqOut + nreqOut);
		MPI_Recv_init(y_down, n_glob, MPI_DOUBLE, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, reqIn + nreqIn);
		nreqIn++;
		nreqOut++;
	}

	MPI_Status * statIn = new  MPI_Status[2];
	MPI_Status * statOut = new  MPI_Status[2];

	int iter = 0;

	//Якоби

	do {

		MPI_Startall(nreqOut, reqOut);
		MPI_Startall(nreqIn, reqIn);
		
		for (int i = 1; i < N_loc; i++)
			for (int j = 1; j < N_glob; j++)
				yk1[i][j] = c*(yk[i + 1][j] + yk[i - 1][j] + yk[i][j + 1] + yk[i][j - 1]) + right[i][j];

		MPI_Waitall(nreqIn, reqIn, statIn);//ждём, когда получим граничные значения 

		if (rank < nproc - 1) {
			for (int j = 1; j < N_glob; j++)
				yk1[N_loc][j] = c*(y_up[j] + yk[N_loc - 1][j] + yk[N_loc][j + 1] + yk[N_loc][j - 1]) + right[N_loc][j];
		}
		if (rank > 0) {
			for (int j = 1; j < N_glob; j++)
				yk1[0][j] = c*(yk[1][j] + y_down[j] + yk[0][j + 1] + yk[0][j - 1]) + right[0][j];
		}

		MPI_Waitall(nreqOut, reqOut, statOut);
		
		swap(yk, yk1);

		iter++;

	} while (iter < 1500);
	


	//Зейдель

	/*
	do {

		MPI_Startall(nreqOut, reqOut);
		MPI_Startall(nreqIn, reqIn);

		for (int i = 1; i < N_loc; i++)
			for (int j = (i % 2 + 1); j < N_glob; j += 2)
				yk[i][j] = c*(yk[i + 1][j] + yk[i - 1][j] + yk[i][j + 1] + yk[i][j - 1]) + right[i][j];

		for (int i = 1; i < N_loc; i++)
			for (int j = ((i + 1) % 2 + 1); j < N_glob; j += 2)
				yk[i][j] = c*(yk[i + 1][j] + yk[i - 1][j] + yk[i][j + 1] + yk[i][j - 1]) + right[i][j];

		MPI_Waitall(nreqIn, reqIn, statIn);//ждём, когда получим граничные значения 

		if (rank < nproc - 1) {
			for (int j = (N_loc % 2 + 1); j < N_glob; j += 2)
				yk[N_loc][j] = c*(y_up[j] + yk[N_loc - 1][j] + yk[N_loc][j + 1] + yk[N_loc][j - 1]) + right[N_loc][j];

			for (int j = ((N_loc + 1) % 2 + 1); j < N_glob; j += 2)
				yk[N_loc][j] = c*(y_up[j] + yk[N_loc - 1][j] + yk[N_loc][j + 1] + yk[N_loc][j - 1]) + right[N_loc][j];
		}
		if (rank > 0) {
			for (int j = 1; j < N_glob; j += 2)
				yk[0][j] = c*(yk[1][j] + y_down[j] + yk[0][j + 1] + yk[0][j - 1]) + right[0][j];

			for (int j = 2; j < N_glob; j += 2)
				yk[0][j] = c*(yk[1][j] + y_down[j] + yk[0][j + 1] + yk[0][j - 1]) + right[0][j];
		}

		MPI_Waitall(nreqOut, reqOut, statOut);

		iter++;

	} while (iter < 1500);
	*/

	//Сбор решения на 0-й процесс
	if (rank == 0)
	{
		for (int i = 0; i < n_loc; ++i)
			for (int j = 0; j < n_glob; ++j) U[i][j] = yk[i][j];

		for (int p = 1; p < nproc; ++p) 
			for (int i = 0; i < nLoc[p]; ++i)
				MPI_Recv(U[i + sdvig[p]], n_glob, MPI_DOUBLE, p, 42, MPI_COMM_WORLD, &stat);
	}
	else
	{
		for (int i = 0; i < n_loc; ++i)
			MPI_Send(yk[i], n_glob, MPI_DOUBLE, 0, 42, MPI_COMM_WORLD);
	}

	mytipe t_end = MPI_Wtime();

	if (rank == 0) {
		cout << "dt = " << t_end - t_init << endl; 
		/*print(n_glob, U);*/ 
		mytipe MAX = 0.0;

		for (int i = 0; i < n_glob; i++)
			for (int j = 0; j < n_glob; j++) if (fabs(U[i][j] - u(i*h,j*h)) > MAX) MAX = fabs(U[i][j] - u(i*h, j*h));
			cout << MAX;

		print_file(n_glob, h, U);

		for (int i = 0; i < n_glob; i++)  delete[] U[i];
		delete[] U; delete[] nLoc; delete[] sdvig;}
			
	for (int i = 0; i < n_loc; i++) {
		delete[] yk1[i]; delete[] yk[i]; delete[] right[i];
	}

	delete[] yk1; delete[] yk; delete[] y_up; delete[] y_down; delete[] right;

	for (int i = 0; i < nreqOut; i++)
		MPI_Request_free(reqOut + i);
	for (int i = 0; i < nreqIn; i++)
		MPI_Request_free(reqIn + i);

	

	MPI_Finalize();


	return 0;
}


void print(int N, mytipe** A)
{
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) cout << A[i][j] << " ";
		cout << endl;
	}
}



void print_file(int N, mytipe h, mytipe** A)
{
	ofstream fout;    // создали переменную для записи в файл
	fout.open("Solve.dat", ios_base::out | ios_base::trunc);
	for (int i = 0; i < N; i += 32) 
		for (int j = 0; j < N; j+=32) fout <<i*h<<" "<<j*h<<" "<<A[i][j]<<endl;
		
	
	fout.close();
	fout.clear();
}