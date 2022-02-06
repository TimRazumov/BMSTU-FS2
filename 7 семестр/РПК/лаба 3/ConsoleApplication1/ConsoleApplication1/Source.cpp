//d/dx(k(x)du/dx)=f(x)

#include<iostream>
#include <fstream>

using namespace std;

const double Eps = 1e-8;//����������� ����

double k(double x)
{
	return exp(x);
}


double f(double x)
{
	return -100*sin(10*x);
}

bool findx(const unsigned int DIM, double* b, double** a);
void hod(const unsigned int DIM, double* b, double** a);

int main()
{
	double a = 0.; 
	double b = 1.; 
	int N = 100;//���-�� ��������� 
	int n = N + 1;//���-�� �����
	double h = (b - a) / double(N);
	double Ua = 10.; double Ub = 20.;//��������� �������
	double K_loc[2][2] = { {1.,-1.}, {-1.,1.} };

	double** K = new double*[n];
	for (int i = 0; i < n; i++)
	{
		K[i] = new double[n]; 
		for (int j = 0; j < n; j++) K[i][j] = 0.;
	}

	double* F = new double [n];
	for (int j = 0; j < n; j++) F[j] = 0.;

	double* x = new double[n];
	for (int j = 0; j < n; j++) x[j] = j*h;
	
	int** Elems = new int* [2];//������� ���������
	for (int i = 0; i < 2; i++) Elems[i] = new int [N];

	for (int j = 0; j < N; j++) { Elems[0][j] = j; Elems[1][j] = j + 1; };

	double M = 10e+10;//�������� ������
	double k_e;
	double h_e;
	int I; int J;
	double c;

	for (int i = 0; i < N; i++) {
		I = Elems[0][i]; J = Elems[1][i];
		h_e = x[J] - x[I];
		k_e = k(0.5*(x[J] + x[I]));
		c = k_e / h_e;
		K[I][I] += c*K_loc[0][0];
		K[I][J] += c*K_loc[0][1];
		K[J][I] += c*K_loc[1][0];
		K[J][J] += c*K_loc[1][1];
		F[I] += 0.5*h_e*f(x[I]);
		F[J] += 0.5*h_e*f(x[J]);
	}

	K[0][0] += M; K[N][N] += M;
	F[0] += M*Ua; F[N] += M*Ub;

	findx(n, F, K);
	hod(n, F, K);

	ofstream fout;    // ������� ���������� ��� ������ � ����
	fout.open("sol.txt", ios_base::out | ios_base::trunc);

	for (int i = 0; i < n; i++) { fout << x[i] << " " << F[i] << endl; }

	fout.close();
	fout.clear();

	for (int i = 0; i < n; i++) delete[] K[i];
	for (int i = 0; i < 2; i++) delete[] Elems[i];

	delete[] K;
	delete[] F;
	delete[] x;
	delete[] Elems;

	cout << "end";
	cin.get();
	return 0;
}


bool findx(const unsigned int DIM, double* b, double** a) //������ ��� ������ � ��������� �� �������������
{
	double max;     //��������������� ���������� ��� ��������� ��������� (���������� max)
	double* vspom;  //��������������� ������ ��� ������������ ����� �������
	double vspom2;   //��������������� ���������� ��� ������������ ��������� �������
	bool flag = true;  //���������� ����������(�������� ������������� �������)


	for (int k = 0, h; k < (DIM); k++)
	{
		max = fabs(a[k][k]);  //������ ������������� �������� � ���������� ���������
		vspom = a[k];   //������ ������ �� ������ ������������� ��������(���� ������� � ����� �� ���� �� ����������)
		h = k;   //������ ������ ������ (�� �������, ��������� ����)

		if (fabs(a[k][k]) < Eps) { flag = true; } //���� ������ ������� �������, �� ������������ ����
		else { flag = false; } //���� �� �������, �� �������������� ����

		for (int i = (k + 1); i < DIM; i++)
		{
			if (fabs(a[i][k]) > max && fabs(a[i][k])>Eps) //���� ������� ������ ����������� max � ������� � �� ����
			{
				max = fabs(a[i][k]); vspom = a[i];  h = i; flag = false;
			}; //��������� ������, �������� max � �������������� ����
		}

		if (!flag) //���� ������� �� �������
		{
			a[h] = a[k];  //����� ��������
			a[k] = vspom;
			vspom2 = b[h]; //����� ���������� �������
			b[h] = b[k];
			b[k] = vspom2;
		}
		else { cout << "������� ���������\n";   return false; }

		for (int i = k + 1; i < DIM; i++)
		{
			b[i] = b[i] - b[k] / a[k][k] * a[i][k]; //��������� �� ��������� �������
			for (int j = (DIM - 1); j >= 0; j--)
			{
				a[i][j] = a[i][j] - a[k][j] / a[k][k] * a[i][k]; //��������� ���� ��������� ��� ������������
			}
		}

	}
	return true;
}

void hod(const unsigned int DIM, double* b, double** a) //�������� ��� ������
{

	for (int i = (DIM - 1); i >= 0; i--)
	{
		for (int j = i + 1; j < DIM; j++)
		{
			b[i] = b[i] - a[i][j] * b[j]; //��������� �� ������� ������ ����� ���� ��������� ������� ���� �� ������ ����� �������� �� �������� ���������
		}
		b[i] = b[i] / a[i][i]; //������� �������� ������� ������ ����� �� �������� ������� ���������
	}
}
