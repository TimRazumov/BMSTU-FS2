#pragma once//директива дл€ того, чтобы файл подкл строго 1 раз

#include <iostream>

typedef double mytipe;


using namespace std;


mytipe norm(const unsigned int DIM, const mytipe* const b, char flag) //3 нормы вектора с запросом варианта нормы
{
	mytipe norm = 0.;  //получаема€ норма вектора
	if (flag == 'k')  //кубическа€ норма (max)
	{
		mytipe* thread_max = new mytipe[4];
		for (int i = 0; i < 4; i++)
			thread_max[i] = 0.;

#pragma omp parallel for num_threads(4)
		for (int i = 0; i < DIM; i++)
		{
			if (thread_max[omp_get_thread_num()] < abs(b[i]))  thread_max[omp_get_thread_num()] = abs(b[i]);   //находим максимум
		}
		for (int i = 0; i < 4; i++)
		{
			if (norm < thread_max[omp_get_thread_num()])  norm = thread_max[omp_get_thread_num()];   //находим максимум
		}
		return norm;
	}

	if (flag == '1') //октаэдрическа€ норма (сумма элементов)
	{
		for (int j = 0; j < DIM; j++)
		{
			norm += abs(b[j]);  //производим сложение элементов
		}

		return norm;
	}

	if (flag == '2')  //шарова€ норма  (≈вклидова)
	{
		for (int i = 0; i<DIM; i++)
		{
			norm += b[i] * b[i];  //суммируем квадраты элементов
		}

		return sqrt(norm); //возвращаем корень квадратный из суммы квадратов
	}
}


void print(const unsigned int DIM, mytipe** const a) //процедура вывода матрицы
{
	cout << endl;
	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			cout << a[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}


void print(const unsigned int DIM, const mytipe* const b) //процедура вывода столбца
{
	cout << '(';
	for (int i = 0; i < DIM; i++)
	{
		cout << b[i];
		if (i < DIM - 1) { cout << ','; } //красива€ запись
	}
	cout << ")^т\n";
}



void copy(const unsigned int DIM, mytipe** a, mytipe** b) //копирование матриц a=b
{
	for (int i = 0; i < DIM; i++)
		for (int j = 0; j < DIM; j++) a[i][j] = b[i][j];

}

void copy(const unsigned int DIM, mytipe* a, mytipe* b) //копирование векторов
{
	for (int j = 0; j < DIM; j++)	a[j] = b[j];

}


void multmatr_minus(int n, int m, mytipe* A, mytipe* L, mytipe* U, int l)//умножение матрицы A размерностью (n_a)х(m_ab)
{
#pragma omp parallel for//на матрицу B размерностью (m_ab)х(n_b) и вычитаетс€ из матрицы ј начина€ с индекса l + n
	for (int i = 0; i < m; i++)
		for (int k = 0; k < n; k++)
			for (int j = 0; j < m; j++)
				A[(i + l + n)*(m + n + l) + j + l + n] -= L[i*n + k] * U[k*m + j];

}


void Ldown_Uright(int N, int M, mytipe* A, int l) //N <= M; N, M - размерности пр€моугольных матриц Ldown и Uright
{
#pragma omp parallel for
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < j; k++)
			{
				A[(i + l + N)*(M + N + l) + j + l] -= A[(k + l)*(M + N + l) + j + l] * A[(i + l + N)*(M + N + l) + k + l]; //Ldown
				A[(j + l)*(M + N + l) + i + l + N] -= A[(j + l)*(M + N + l) + k + l] * A[(k + l)*(M + N + l) + i + l + N]; //Uright
			}
			A[(i + l + N)*(M + N + l) + j + l] /= A[(j + l)*(M + N + l) + j + l]; //делим на элемент на главной диагонали U
		}
	}
}



void A_is_LU(int N, mytipe* L, mytipe* U, mytipe* A)  //проверка по максимальной норме на равенство исходной матрицы и матрицы, получаемой путем перемножени€ L и U
{
#pragma omp parallel for
	for (int i = 0; i < N; i++)
		for (int k = 0; k < N; k++)
			for (int j = 0; j < N; j++)
				A[i*N + j] -= L[i*N + k] * U[k*N + j];

	//print(N*N, A);

	cout << "\n||A-LU|| = " << N*norm(N*N, A, 'k') << endl;
}
