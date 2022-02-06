#pragma once//директива для того, чтобы файл подкл строго 1 раз

#include <iostream>

typedef double mytipe;


using namespace std;

mytipe norm(const unsigned int DIM, mytipe** const A, char flag) //4 нормы матрицы с запросом варианта нормы
{
	mytipe norm = 0;  //получаемая норма матрицы
	mytipe vspom = 0;  //вспомогательная переменная
	if (flag == '1')  //октаэдрическая норма матрицы (max суммируя по столбцам)
	{
		for (int i = 0; i < DIM; i++)
		{
			vspom = 0;
			for (int j = 0; j<DIM; j++)
			{
				vspom += fabs(A[j][i]); //суммируем элементы i-ой строки
			}
			if (norm < vspom) { norm = vspom; };   //проверяем на максимум
		}
		return norm;
	}

	if (flag == 'k')  //кубическая норма матрицы (max суммирую по строкам)
	{
		for (int j = 0; j < DIM; j++)
		{
			vspom = 0;
			for (int i = 0; i<DIM; i++)
			{
				vspom += fabs(A[j][i]);  //суммируем элементы j-ой строки
			}
			if (norm < vspom) { norm = vspom; };   //проверяем на максимум
		}

		return norm;
	}

	if (flag == '2')  //шаровая норма матрицы (корень из суммы квадратов всех элементов)
	{
		for (int j = 0; j < DIM; j++)
		{
			for (int i = 0; i<DIM; i++)
			{
				norm += A[j][i] * A[j][i];  //складываем квадраты
			}
		}

		return sqrt(norm);  //возвращаем корень из суммы квадратов
	}

	if (flag == 'm')  //максимальная норма матрицы (n*max)
	{
		mytipe* Tred_max = new mytipe[4];
		for (int i = 0; i < 4; i++) Tred_max[i] = 0;
#pragma omp parallel for 
		for (int j = 0; j < DIM; j++)
		{
			for (int i = 0; i<DIM; i++)
			{
				if (fabs(A[j][i]) > Tred_max[omp_get_thread_num()]) { Tred_max[omp_get_thread_num()] = fabs(A[j][i]); } //находим максимальный элемент
			}
		}
		for (int i = 0; i < 4; i++)
			if (norm < Tred_max[i]) { norm = Tred_max[omp_get_thread_num()]; } //находим максимальный элемент;
		return DIM*norm;
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
		if (i < DIM - 1) { cout << ','; } //красивая запись
	}
	cout << ")^т\n";
}



void copy(const unsigned int DIM, mytipe** a, mytipe** b) //копирование матриц a=b
{

	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			a[i][j] = b[i][j];
		}
	}
}

void copy(const unsigned int DIM, mytipe* a, mytipe* b) //копирование векторов
{
	for (int j = 0; j < DIM; j++)
	{
		a[j] = b[j];
	}
}


void trans(const unsigned int DIM, mytipe** a) //транспонирование матрицы
{
	mytipe vspom;
	for (int i = 0; i < DIM; i++)
	{

		for (int j = i + 1; j < DIM; j++)
		{
			vspom = a[i][j];
			a[i][j] = a[j][i];
			a[j][i] = vspom;
		}
	}
}


void singlematr(const unsigned int DIM, mytipe** const E) //заполнение матрицы единичной
{
	for (int i = 0; i < DIM; i++)      //|делаем
	{									//|матрицу Т
		for (int j = 0; j < DIM; j++)  //|единичной
		{							    //|для
			E[i][j] = 0.;				//|дальнейших
		}								//|преобразований
		E[i][i] = 1.;					//|
	}
}
void zeromatr(const int n, const int m, mytipe** const Z) //заполнение матрицы нулевой
{
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			Z[i][j] = 0.;
}

void multmatr(int n_a, int m_ab, int n_b, mytipe** A, mytipe** B, mytipe** C)//умножение матрицы A размерностью (n_a)х(m_ab)
{                                                                            //на матрицу B размерностью (m_ab)х(n_b)
	zeromatr(n_a, n_b, C);
#pragma omp parallel for
	for (int i = 0; i < n_a; i++)
		for (int k = 0; k < m_ab; k++)
			for (int j = 0; j < n_b; j++)
				C[i][j] += A[i][k] * B[k][j];

}

void multmatr_minus(int n, int m, mytipe** A, mytipe** L, mytipe** U, int l)//умножение матрицы A размерностью (n_a)х(m_ab)
{                                                   //на матрицу B размерностью (m_ab)х(n_b) и вычитается из матрицы А начиная с индекса l + n
	for (int i = 0; i < m; i++)
		for (int k = 0; k < n; k++)
			for (int j = 0; j < m; j++)
				A[i + l + n][j + l + n] -= L[i][k] * U[k][j];

}


void Ldown_Uright(int N, int M, mytipe** A, int l) //N <= M; N, M - размерности прямоугольных матриц Ldown и Uright
{
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < j; k++)
			{
				A[i + l + N][j + l] -= A[k + l][j + l] * A[i + l + N][k + l]; //Ldown
				A[j + l][i + l + N] -= A[j + l][k + l] * A[k + l][i + l + N]; //Uright
			}
			A[i + l + N][j + l] /= A[j + l][j + l]; //делим на элемент на главной диагонали U
		}
	}
}

void multmatr_minus_Paral(int n, int m, mytipe** A, mytipe** L, mytipe** U, int l)//умножение матрицы A размерностью (n_a)х(m_ab)
{
#pragma omp parallel for                            //на матрицу B размерностью (m_ab)х(n_b) и вычитается из матрицы А начиная с индекса l + n
	for (int i = 0; i < m; i++)
		for (int k = 0; k < n; k++)
			for (int j = 0; j < m; j++)
				A[i + l + n][j + l + n] -= L[i][k] * U[k][j];

}


void Ldown_Uright_Paral(int N, int M, mytipe** A, int l) //N <= M; N, M - размерности прямоугольных матриц Ldown и Uright
{
#pragma omp parallel for
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < j; k++)
			{
				A[i + l + N][j + l] -= A[k + l][j + l] * A[i + l + N][k + l]; //Ldown
				A[j + l][i + l + N] -= A[j + l][k + l] * A[k + l][i + l + N]; //Uright
			}
			A[i + l + N][j + l] /= A[j + l][j + l]; //делим на элемент на главной диагонали U
		}
	}
}


void A_is_LU(int N, mytipe** L, mytipe** U, mytipe** A)  //проверка по максимальной норме на равенство исходной матрицы и матрицы, получаемой путем перемножения L и U
{
#pragma omp parallel for
	for (int i = 0; i < N; i++)
		for (int k = 0; k < N; k++)
			for (int j = 0; j < N; j++)
				A[i][j] -= L[i][k] * U[k][j];

	//print(N*N, A);

	cout << "\n||A-LU|| = " << norm(N, A, 'm') << endl;
}