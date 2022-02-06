#pragma once//��������� ��� ����, ����� ���� ����� ������ 1 ���

#include <iostream>

typedef double mytipe;


using namespace std;

mytipe norm(const unsigned int DIM, mytipe** const A, char flag) //4 ����� ������� � �������� �������� �����
{
	mytipe norm = 0;  //���������� ����� �������
	mytipe vspom = 0;  //��������������� ����������
	if (flag == '1')  //�������������� ����� ������� (max �������� �� ��������)
	{
		for (int i = 0; i < DIM; i++)
		{
			vspom = 0;
			for (int j = 0; j<DIM; j++)
			{
				vspom += fabs(A[j][i]); //��������� �������� i-�� ������
			}
			if (norm < vspom) { norm = vspom; };   //��������� �� ��������
		}
		return norm;
	}

	if (flag == 'k')  //���������� ����� ������� (max �������� �� �������)
	{
		for (int j = 0; j < DIM; j++)
		{
			vspom = 0;
			for (int i = 0; i<DIM; i++)
			{
				vspom += fabs(A[j][i]);  //��������� �������� j-�� ������
			}
			if (norm < vspom) { norm = vspom; };   //��������� �� ��������
		}

		return norm;
	}

	if (flag == '2')  //������� ����� ������� (������ �� ����� ��������� ���� ���������)
	{
		for (int j = 0; j < DIM; j++)
		{
			for (int i = 0; i<DIM; i++)
			{
				norm += A[j][i] * A[j][i];  //���������� ��������
			}
		}

		return sqrt(norm);  //���������� ������ �� ����� ���������
	}

	if (flag == 'm')  //������������ ����� ������� (n*max)
	{
		mytipe* Tred_max = new mytipe[4];
		for (int i = 0; i < 4; i++) Tred_max[i] = 0;
#pragma omp parallel for 
		for (int j = 0; j < DIM; j++)
		{
			for (int i = 0; i<DIM; i++)
			{
				if (fabs(A[j][i]) > Tred_max[omp_get_thread_num()]) { Tred_max[omp_get_thread_num()] = fabs(A[j][i]); } //������� ������������ �������
			}
		}
		for (int i = 0; i < 4; i++)
			if (norm < Tred_max[i]) { norm = Tred_max[omp_get_thread_num()]; } //������� ������������ �������;
		return DIM*norm;
	}
}




void print(const unsigned int DIM, mytipe** const a) //��������� ������ �������
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


void print(const unsigned int DIM, const mytipe* const b) //��������� ������ �������
{
	cout << '(';
	for (int i = 0; i < DIM; i++)
	{
		cout << b[i];
		if (i < DIM - 1) { cout << ','; } //�������� ������
	}
	cout << ")^�\n";
}



void copy(const unsigned int DIM, mytipe** a, mytipe** b) //����������� ������ a=b
{

	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			a[i][j] = b[i][j];
		}
	}
}

void copy(const unsigned int DIM, mytipe* a, mytipe* b) //����������� ��������
{
	for (int j = 0; j < DIM; j++)
	{
		a[j] = b[j];
	}
}


void trans(const unsigned int DIM, mytipe** a) //���������������� �������
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


void singlematr(const unsigned int DIM, mytipe** const E) //���������� ������� ���������
{
	for (int i = 0; i < DIM; i++)      //|������
	{									//|������� �
		for (int j = 0; j < DIM; j++)  //|���������
		{							    //|���
			E[i][j] = 0.;				//|����������
		}								//|��������������
		E[i][i] = 1.;					//|
	}
}
void zeromatr(const int n, const int m, mytipe** const Z) //���������� ������� �������
{
#pragma omp parallel for
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			Z[i][j] = 0.;
}

void multmatr(int n_a, int m_ab, int n_b, mytipe** A, mytipe** B, mytipe** C)//��������� ������� A ������������ (n_a)�(m_ab)
{                                                                            //�� ������� B ������������ (m_ab)�(n_b)
	zeromatr(n_a, n_b, C);
#pragma omp parallel for
	for (int i = 0; i < n_a; i++)
		for (int k = 0; k < m_ab; k++)
			for (int j = 0; j < n_b; j++)
				C[i][j] += A[i][k] * B[k][j];

}

void multmatr_minus(int n, int m, mytipe** A, mytipe** L, mytipe** U, int l)//��������� ������� A ������������ (n_a)�(m_ab)
{                                                   //�� ������� B ������������ (m_ab)�(n_b) � ���������� �� ������� � ������� � ������� l + n
	for (int i = 0; i < m; i++)
		for (int k = 0; k < n; k++)
			for (int j = 0; j < m; j++)
				A[i + l + n][j + l + n] -= L[i][k] * U[k][j];

}


void Ldown_Uright(int N, int M, mytipe** A, int l) //N <= M; N, M - ����������� ������������� ������ Ldown � Uright
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
			A[i + l + N][j + l] /= A[j + l][j + l]; //����� �� ������� �� ������� ��������� U
		}
	}
}

void multmatr_minus_Paral(int n, int m, mytipe** A, mytipe** L, mytipe** U, int l)//��������� ������� A ������������ (n_a)�(m_ab)
{
#pragma omp parallel for                            //�� ������� B ������������ (m_ab)�(n_b) � ���������� �� ������� � ������� � ������� l + n
	for (int i = 0; i < m; i++)
		for (int k = 0; k < n; k++)
			for (int j = 0; j < m; j++)
				A[i + l + n][j + l + n] -= L[i][k] * U[k][j];

}


void Ldown_Uright_Paral(int N, int M, mytipe** A, int l) //N <= M; N, M - ����������� ������������� ������ Ldown � Uright
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
			A[i + l + N][j + l] /= A[j + l][j + l]; //����� �� ������� �� ������� ��������� U
		}
	}
}


void A_is_LU(int N, mytipe** L, mytipe** U, mytipe** A)  //�������� �� ������������ ����� �� ��������� �������� ������� � �������, ���������� ����� ������������ L � U
{
#pragma omp parallel for
	for (int i = 0; i < N; i++)
		for (int k = 0; k < N; k++)
			for (int j = 0; j < N; j++)
				A[i][j] -= L[i][k] * U[k][j];

	//print(N*N, A);

	cout << "\n||A-LU|| = " << norm(N, A, 'm') << endl;
}