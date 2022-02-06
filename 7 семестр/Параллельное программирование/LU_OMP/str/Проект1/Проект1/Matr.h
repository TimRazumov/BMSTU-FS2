#pragma once//��������� ��� ����, ����� ���� ����� ������ 1 ���

#include <iostream>

typedef double mytipe;


using namespace std;


mytipe norm(const unsigned int DIM, const mytipe* const b, char flag) //3 ����� ������� � �������� �������� �����
{
	mytipe norm = 0.;  //���������� ����� �������
	if (flag == 'k')  //���������� ����� (max)
	{
		mytipe* thread_max = new mytipe[4];
		for (int i = 0; i < 4; i++)
			thread_max[i] = 0.;

#pragma omp parallel for num_threads(4)
		for (int i = 0; i < DIM; i++)
		{
			if (thread_max[omp_get_thread_num()] < abs(b[i]))  thread_max[omp_get_thread_num()] = abs(b[i]);   //������� ��������
		}
		for (int i = 0; i < 4; i++)
		{
			if (norm < thread_max[omp_get_thread_num()])  norm = thread_max[omp_get_thread_num()];   //������� ��������
		}
		return norm;
	}

	if (flag == '1') //�������������� ����� (����� ���������)
	{
		for (int j = 0; j < DIM; j++)
		{
			norm += abs(b[j]);  //���������� �������� ���������
		}

		return norm;
	}

	if (flag == '2')  //������� �����  (���������)
	{
		for (int i = 0; i<DIM; i++)
		{
			norm += b[i] * b[i];  //��������� �������� ���������
		}

		return sqrt(norm); //���������� ������ ���������� �� ����� ���������
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
		for (int j = 0; j < DIM; j++) a[i][j] = b[i][j];

}

void copy(const unsigned int DIM, mytipe* a, mytipe* b) //����������� ��������
{
	for (int j = 0; j < DIM; j++)	a[j] = b[j];

}


void multmatr_minus(int n, int m, mytipe* A, mytipe* L, mytipe* U, int l)//��������� ������� A ������������ (n_a)�(m_ab)
{
#pragma omp parallel for//�� ������� B ������������ (m_ab)�(n_b) � ���������� �� ������� � ������� � ������� l + n
	for (int i = 0; i < m; i++)
		for (int k = 0; k < n; k++)
			for (int j = 0; j < m; j++)
				A[(i + l + n)*(m + n + l) + j + l + n] -= L[i*n + k] * U[k*m + j];

}


void Ldown_Uright(int N, int M, mytipe* A, int l) //N <= M; N, M - ����������� ������������� ������ Ldown � Uright
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
			A[(i + l + N)*(M + N + l) + j + l] /= A[(j + l)*(M + N + l) + j + l]; //����� �� ������� �� ������� ��������� U
		}
	}
}



void A_is_LU(int N, mytipe* L, mytipe* U, mytipe* A)  //�������� �� ������������ ����� �� ��������� �������� ������� � �������, ���������� ����� ������������ L � U
{
#pragma omp parallel for
	for (int i = 0; i < N; i++)
		for (int k = 0; k < N; k++)
			for (int j = 0; j < N; j++)
				A[i*N + j] -= L[i*N + k] * U[k*N + j];

	//print(N*N, A);

	cout << "\n||A-LU|| = " << N*norm(N*N, A, 'k') << endl;
}
