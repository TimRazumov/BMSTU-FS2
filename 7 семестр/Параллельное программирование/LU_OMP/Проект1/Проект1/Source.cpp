#include "iostream"
#include <fstream>
#include <cmath>
#include <omp.h>
#include "Matr.h"

using namespace std;

typedef double mytipe;  //тип данных, использующийся во всей программе
const mytipe EPSILON = 10e-8;//константа сравнения с 0

void LU(int N, mytipe** A, int l);
void block_LU(int N, mytipe** A, int b);
void LU_Paral(int N, mytipe** A, int l);
void block_LU_Paral(int N, mytipe** A, int b);

void printfile(int N, mytipe** A) {
	FILE* f;
	fopen_s(&f, "LU_matrix.txt", "w");
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			fprintf(f,"%lf ",A[i][j]);
		}
		fprintf(f,"\n");
	}
	fclose(f);
}

int main() {
	setlocale(LC_ALL, "Russian"); //подключаем русский язык
	int n = 2*2048;//размерность матрицы
	mytipe** A = new mytipe*[n];
	mytipe** L = new mytipe*[n];
	mytipe** A1 = new mytipe*[n];

	for (int i = 0; i < n; i++) {
		A[i] = new mytipe[n];
		A1[i] = new mytipe[n];
	}


	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			A[i][j] = rand() / 10000.;
		}

	cout << "A = ";
	//print(n, A);
	
	int threadsNum = 4;
	omp_set_num_threads(threadsNum);

	int b = 64;//размер блока
	mytipe t_par, t, t1, t2;

	for (int i = 1; i <= 4; i++)
	{
		omp_set_num_threads(i);
		copy(n, A1, A);
		cout << "\nПотоков = " << i << endl;
		cout << "Classic LU\n";
		t1 = omp_get_wtime();
		LU_Paral(n, A1, 0);
		t2 = omp_get_wtime();
		t_par = t2 - t1;
		cout <<  "t = " << t_par  << endl;
	}

	for (int i = 1; i <= 4; i++)
	{
		omp_set_num_threads(i);
		cout << "\nПотоков = " << i << endl;
		cout << "Block LU\nb = " << b << endl;
		copy(n, A1, A);
		t1 = omp_get_wtime();
		block_LU_Paral(n, A1, b);
		t2 = omp_get_wtime();
		t_par = t2 - t1;
		cout << "t = " << t_par << endl;
	}


	//while (b < n / 4) {
	//	copy(n, A1, A);
	//	b *= 2;
	//	cout << "\nBlock LU\nb = " << b << endl;
	//	t1 = omp_get_wtime();
	//	block_LU(n, A1, b);
	//	t2 = omp_get_wtime();
	//	t = t2 - t1;
	//	t1 = omp_get_wtime();
	//	copy(n, A1, A);
	//	block_LU_Paral(n, A1, b);
	//	t2 = omp_get_wtime();
	//	t_par = t2 - t1;
	//	cout << "t = " << t << " , t_parallel = " << t_par << " , t/t_parallel = " << t / t_par << endl;
	//}




	/*cout << "\n\nL|U = ";
	print(n, A);*/



	for (int i = 0; i < n; i++) {
		L[i] = new mytipe[n];
	}

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			if (i>j)
			{
				L[i][j] = A1[i][j];
				A1[i][j] = 0.;
			}
			if (i == j)
				L[i][j] = 1.;
			if (i < j)
				L[i][j] = 0.0;
		}

	A_is_LU(n, L, A1, A);

	for (int i = 0; i < n; i++) {
		delete[] A[i];
		delete[] L[i];
		delete[] A1[i];
	}

	delete[] A;
	delete[] L;
	delete[] A1;

	cout << "\n end";
	cin.get();
	return 0;
}

void LU_Paral(int N, mytipe** A, int l)
{
	for (int i = 0; i < N; i++)
#pragma omp parallel for 
		for (int j = i + 1; j < N; j++)
		{
			A[j + l][i + l] /= A[i + l][i + l];

			for (int k = i + 1; k < N; k++) A[j + l][k + l] -= A[i + l][k + l] * A[j + l][i + l];

		}
}

void LU(int N, mytipe** A, int l)
{
	for (int i = 0; i < N; i++)
		for (int j = i + 1; j < N; j++)
		{
			A[j + l][i + l] /= A[i + l][i + l];

			for (int k = i + 1; k < N; k++) A[j + l][k + l] -= A[i + l][k + l] * A[j + l][i + l];

		}
}


void block_LU_Paral(int N, mytipe** A, int b)
{

	int b1 = b;


	for (int i = 0; i < N; i += b)
	{
		if (N - i < b) b1 = N - i;

		mytipe** L_down = new mytipe*[N - i - b1];
		mytipe** U_right = new mytipe*[b1];

		LU(b1, A, i);

		Ldown_Uright_Paral(b1, N - i - b1, A, i);

		//-------Копирование указателей------------
		for (int k = 0; k < N - i - b1; k++) {
			L_down[k] = &(A[k + i + b1][i]);
		}
		for (int k = 0; k < b1; k++) {
			U_right[k] = &(A[k + i][i + b1]);
		}

		//-------Копирование матриц-----------
		//for (int k = 0; k < N - i - b1; k++) {
		//	L_down[k] = new mytipe[b1];
		//	for (int j = 0; j < b1; j++) {
		//		L_down[k][j] = A[k + i + b1][i + j];
		//	}
		//}
		//for (int k = 0; k < b1; k++) {
		//	U_right[k] = new mytipe[N - i - b1];
		//	for (int j = 0; j < N - i - b1; j++) {
		//		U_right[k][j] = A[k + i][j + i + b1];
		//		}
		//}

		multmatr_minus_Paral(b1, N - i - b1, A, L_down, U_right, i);

		delete[] L_down;
		delete[] U_right;
	}

}

void block_LU(int N, mytipe** A, int b)
{

	int b1 = b;

	mytipe** L_down = new mytipe*[N - b1];
	mytipe** U_right = new mytipe*[b1];

	for (int i = 0; i < N; i += b)
	{
		if (N - i < b) b1 = N - i;
		

		LU(b1, A, i);

		Ldown_Uright(b1, N - i - b1, A, i);

		//-------Копирование указателей------------
		for (int k = 0; k < N - i - b1; k++) {
			L_down[k] = &(A[k + i + b1][i]);
		}
		for (int k = 0; k < b1; k++) {
			U_right[k] = &(A[k + i][i + b1]);
		}

		//-------Копирование матриц-----------
		//for (int k = 0; k < N - i - b1; k++) {
		//	L_down[k] = new mytipe[b1];
		//	for (int j = 0; j < b1; j++) {
		//		L_down[k][j] = A[k + i + b1][i + j];
		//	}
		//}
		//for (int k = 0; k < b1; k++) {
		//	U_right[k] = new mytipe[N - i - b1];
		//	for (int j = 0; j < N - i - b1; j++) {
		//		U_right[k][j] = A[k + i][j + i + b1];
		//		}
		//}

		multmatr_minus(b1, N - i - b1, A, L_down, U_right, i);

		
	}

	delete[] L_down;
	delete[] U_right;

}






