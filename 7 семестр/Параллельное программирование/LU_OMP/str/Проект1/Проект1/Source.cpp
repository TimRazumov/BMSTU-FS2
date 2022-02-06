#include "iostream"
#include <fstream>
#include <cmath>
#include <omp.h>
#include "Matr.h"

using namespace std;

typedef double mytipe;  //тип данных, использующийся во всей программе
const mytipe EPSILON = 10e-8;//константа сравнения с 0

void LU(int N, int b, mytipe* A, int l);
void LU_Paral(int N, int b, mytipe* A, int l);
void block_LU_Paral(int N, mytipe* A, int b);


int main() {
	setlocale(LC_ALL, "Russian"); //подключаем русский язык
	int n = 2*2048;//размерность матрицы

	mytipe* A = new mytipe[n*n];
	mytipe* L = new mytipe[n*n];
	mytipe* A1 = new mytipe[n*n];


	for (int i = 0; i <n*n; i++) A[i] = rand() / 10000.;



	int b = 64;//размер блока
	mytipe t_par, t, t1, t2;

	for (int i = 1; i <= 4; i++)
	{
		omp_set_num_threads(i);
		copy(n*n, A1, A);
		cout << "\nПотоков = " << i << endl;
		cout << "Classic LU\n";
		t1 = omp_get_wtime();
		LU_Paral(n, n, A1, 0);
		t2 = omp_get_wtime();
		t_par = t2 - t1;
		cout << "t = " << t_par << endl;
	}

	for (int i = 4; i <= 4; i++)
	{
		omp_set_num_threads(i);
		cout << "\nПотоков = " << i << endl;
		cout << "Block LU\nb = " << b << endl;
		copy(n*n, A1, A);
		t1 = omp_get_wtime();
		block_LU_Paral(n, A1, b);
		t2 = omp_get_wtime();
		t_par = t2 - t1;
		cout << "t = " << t_par << endl;
	}



	//cout << "\n\nA = ";
	//print(n*n, A);
	//cout << "\n\nL|U = ";
	//print(n*n, A1);


	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			if (i>j)
			{
				L[i*n + j] = A1[i*n + j];
				A1[i*n + j] = 0.;
			}
			if (i == j)
				L[i*n + j] = 1.;
			if (i < j)
				L[i*n + j] = 0.0;
		}

	A_is_LU(n, L, A1, A);


	delete[] A;
	delete[] L;
	delete[] A1;

	cout << "\n end";
	cin.get();
	return 0;
}


void LU_Paral(int N, int b, mytipe* A, int l)
{
	for (int i = 0; i < b; i++) {
#pragma omp parallel for 
		for (int j = i + 1; j < b; j++)
		{
			A[(j + l)*N + i + l] /= A[(i + l)*N + i + l];

			for (int k = i + 1; k < b; k++) A[(j + l)*N + k + l] -= A[(i + l)*N + k + l] * A[(j + l)*N + i + l];
		}
	}

}

void LU(int N, int b, mytipe* A, int l)
{
	for (int i = 0; i < b; i++) {

		for (int j = i + 1; j < b; j++)
		{
			A[(j + l)*N + i + l] /= A[(i + l)*N + i + l];

			for (int k = i + 1; k < b; k++) A[(j + l)*N + k + l] -= A[(i + l)*N + k + l] * A[(j + l)*N + i + l];
		}
	}

}


void block_LU_Paral(int N, mytipe* A, int b)
{

	int b1 = b;

	mytipe* L_down = new mytipe[(N - b1)*b1];
	mytipe* U_right = new mytipe[(N - b1)*b1];

	for (int i = 0; i < N; i += b)
	{
		if (N - i < b) b1 = N - i;

		LU(N, b1, A, i);


		Ldown_Uright(b1, N - i - b1, A, i);


		//-------Копирование матриц-----------
		for (int k = 0; k < N - i - b1; k++) {
			for (int j = 0; j < b1; j++) {
				L_down[k*b1 + j] = A[(k + i + b1)*N + i + j];
			}
		}
		for (int k = 0; k < b1; k++) {
			for (int j = 0; j < N - i - b1; j++) {
				U_right[k*(N - i - b1) + j] = A[(k + i)*N + j + i + b1];
			}
		}

		multmatr_minus(b1, N - i - b1, A, L_down, U_right, i);


	}

	delete[] L_down;
	delete[] U_right;

}