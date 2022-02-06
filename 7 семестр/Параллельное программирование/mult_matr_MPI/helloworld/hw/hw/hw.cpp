#include<iostream>
#include <fstream>
#include <cmath>
#include<mpi.h>

using namespace std;

typedef double mytipe;  //тип данных, использующийся во всей программе

const mytipe EPS = 10e-8;//константа сравнения с 0

void print(int N, mytipe** A);

int main(int argc, char **argv)
{
	const int N_glob = 1024; 

	mytipe** A = new mytipe*[N_glob];
	mytipe** B = new mytipe*[N_glob]; 

	for (int i = 0; i < N_glob; i++) {
		B[i] = new mytipe[N_glob]; A[i] = new mytipe[N_glob];
		for (int j = 0; j < N_glob; j++) { A[i][j] = i + j; B[i][j] = i + j; }
	}

	MPI_Init(&argc, &argv); 

	MPI_Status stat;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);//номер исполняемового процесса
	int nproc;
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);//кол-во исполняемых процессов

	int* nLoc = nullptr;//кол-во разбиений матрицы A для каждого процесса
	int n_loc;//Кол-во строк марицы A_loc для текущего процесса
	int* sdvig = nullptr;//начало считывания для каждого процесса из A в A_loc
	int sdvig_proc;//начало считывания для текущего процесса из A в A_loc

	
	mytipe** C = nullptr;

	mytipe t_init = MPI_Wtime();

	if (rank == 0)
	{
		C = new mytipe*[N_glob];

		for (int i = 0; i < N_glob; i++) C[i] = new mytipe[N_glob];

		nLoc = new int[nproc];
		for (int i = 0; i < nproc; i++) nLoc[i] = N_glob / nproc;//локальная размерность матрицы для каждого процесса
		for (int i = 0; i < N_glob % nproc; i++) nLoc[i] += 1;
		
		sdvig = new int[nproc];
		sdvig[0] = 0; 
		for (int i = 1; i < nproc; i++) sdvig[i] = sdvig[i - 1] + nLoc[i - 1];

	}

	//for (int i = 0; i < N_glob; i++) MPI_Bcast(B[i], N_glob, MPI_DOUBLE, 0, MPI_COMM_WORLD);//рассылвем B на все процессы

	MPI_Scatter(nLoc, 1, MPI_INT, &n_loc, 1, MPI_INT, 0, MPI_COMM_WORLD);//отправляем каждому процессу свое значение nLoc, записывая в переменную n_loc
	MPI_Scatter(sdvig, 1, MPI_INT, &sdvig_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);

	mytipe** A_loc = new mytipe*[n_loc]; mytipe** C_loc = new mytipe*[n_loc];

	for (int i = 0; i < n_loc; i++) {
		A_loc[i] = new mytipe[N_glob];
		C_loc[i] = new mytipe[N_glob];
		for (int j = 0; j < N_glob; j++) C_loc[i][j] = 0.;
	}

	/*if (rank == 0)
	{
		for (int i = 0; i < n_loc; ++i)
			for (int j = 0; j < N_glob; ++j) A_loc[i][j] = A[i][j];

		int vspom = 0;
		for (int p = 1; p < nproc; ++p) {
			vspom += nLoc[p - 1];
			for (int i = vspom; i < vspom + nLoc[p]; ++i)
				MPI_Send(A[i], N_glob, MPI_DOUBLE, p, 42, MPI_COMM_WORLD);
		}
	}
	else
	{
		for (int i = 0; i < n_loc; ++i)
			MPI_Recv(A_loc[i], N_glob, MPI_DOUBLE, 0, 42, MPI_COMM_WORLD, &stat);
	}*/

	for (int i = 0; i < n_loc; ++i)
		for (int j = 0; j < N_glob; ++j) A_loc[i][j] = A[i + sdvig_proc][j];


	for (int i = 0; i < n_loc; i++)
		for (int k = 0; k < N_glob; k++)
			for (int j = 0; j < N_glob; j++)
				C_loc[i][j] += A_loc[i][k]*B[k][j];


	//Сбор решения на 0-й процесс
	if (rank == 0)
	{
		for (int i = 0; i < n_loc; ++i)
			for (int j = 0; j < N_glob; ++j) C[i][j] = C_loc[i][j];

		for (int p = 1; p < nproc; ++p) 
			for (int i = 0; i < nLoc[p]; ++i)
				MPI_Recv(C[i + sdvig[p]], N_glob, MPI_DOUBLE, p, 42, MPI_COMM_WORLD, &stat);
	}
	else
	{
		for (int i = 0; i < n_loc; ++i)
			MPI_Send(C_loc[i], N_glob, MPI_DOUBLE, 0, 42, MPI_COMM_WORLD);
	}

	mytipe t_end = MPI_Wtime();

	if (rank == 0) { cout << "dt = " << t_end - t_init << endl; /*print(N_glob, C); */
	for (int i = 0; i < N_glob; i++)  delete[] C[i];
	delete[] C; delete[] nLoc; delete[] sdvig;}
			
	for (int i = 0; i < n_loc; i++) {
		delete[] A_loc[i];
		delete[] C_loc[i];
	}

	delete[] A_loc; delete[] C_loc;

	MPI_Finalize();

	for (int i = 0; i < N_glob; i++) { delete[] A[i]; delete[] B[i];  }
	
	delete[] A; delete[] B; 

	return 0;
}


void print(int N, mytipe** A)
{
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) cout << A[i][j] << " ";
		cout << endl;
	}
}

