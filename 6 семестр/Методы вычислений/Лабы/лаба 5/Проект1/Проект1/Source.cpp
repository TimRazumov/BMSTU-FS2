#include "iostream"
#include <fstream>
#include <typeinfo>
#include <cmath>

using namespace std;

const double EPSILON = 10e-18;   //��������� ��������� � 0
const double Pi = 3.1415926535897932385;
typedef double mytipe;  //��� ������, �������������� �� ���� ���������


//part 1
void Quadrature_method(mytipe a, mytipe b, int N);
void Simple_iter_method(mytipe a, mytipe b, int N);
void Degenerate_kernel(mytipe a, mytipe b, int N, int m);
void Degenerate_kernel_Cheb(mytipe a, mytipe b, int N, int m);

//part 2
void Singular_kernel(int N);

void hod(const unsigned int DIM, mytipe* b, mytipe** a); //�������� ��� ������
bool findx(const unsigned int DIM, mytipe* b, mytipe** a);  //������ ��� ������ � ������� ������������� �������

void copy(const unsigned int DIM, mytipe** a, mytipe** b); //����������� ������
void copy(const unsigned int DIM, mytipe* a, mytipe* b); //����������� ��������
void sol_to_file(int N, mytipe* x, mytipe* Y);
void sol_to_file(int N, mytipe* x, mytipe* y, mytipe* g);

void vivod(const unsigned int DIM, mytipe** const a);  //����� �������
void vivod(const unsigned int DIM, const mytipe* const b);  //����� �������
void error(int N, mytipe* y);
mytipe phiC(mytipe x, int i);
mytipe psiC(mytipe s, int i);

mytipe f(mytipe x)
{
	//part1
	return 0.5*(1. + sin(x));
	//return  x*x + sqrt(x);

	//part2
	//return cos(7*x);
}

mytipe K(mytipe x, mytipe s)
{
	return 0.5*(1. - x*cos(x*s));
}

mytipe phi(mytipe x, int i);
mytipe psi(mytipe x, int i);

int main() {
	setlocale(LC_ALL, "Russian"); //���������� ������� ����

								  //part 1
	mytipe a = 0.;
	mytipe b = 1.;
	mytipe h = 0.1;
	int N = round((b - a) / h) + 1;

	//Quadrature_method(a, b, N);
	//Simple_iter_method(a, b, N);
	//Degenerate_kernel(a, b, N, 4);
	Degenerate_kernel_Cheb(a, b, N, 4);


	//part 2
	//N = 50;
	//Singular_kernel(N);

	cout << "end:)";
	cin.get();
	return 0;
}


void Quadrature_method(mytipe a, mytipe b, int N)
{
	mytipe h = (b - a) / (N - 1.);

	mytipe** A = new mytipe*[N];
	for (int i = 0; i < N; i++) { A[i] = new mytipe[N]; }

	mytipe* x = new mytipe[N];
	mytipe* d = new mytipe[N];//������ �����

	for (int i = 0; i < N; i++) { x[i] = a + i*h; d[i] = f(x[i]); }


	for (int i = 0; i < N; i++) {

		A[i][0] = -h / 2.*K(x[i], x[0]);
		for (int k = 1; k < N - 1; k++) { A[i][k] = -h*K(x[i], x[k]); }
		A[i][N - 1] = -h / 2.*K(x[i], x[N - 1]);

		A[i][i] += 1.;
	}

	findx(N, d, A);
	hod(N, d, A);

	error(N, d);

	sol_to_file(N, x, d);

	for (int i = 0; i < N; i++) { delete[] A[i]; }
	delete[] A;
	delete[] d;
	delete[] x;

}

void Simple_iter_method(mytipe a, mytipe b, int N)
{

	mytipe h = (b - a) / (N - 1.);
	mytipe I;//�������� ���������
	mytipe max, vspom;
	mytipe eps = 1e-8;//��������
	int k = 0;//���-�� ����

	mytipe* x = new mytipe[N];
	mytipe* uk = new mytipe[N];
	mytipe* uk1 = new mytipe[N];
	for (int i = 0; i < N; i++) { x[i] = a + i*h; uk[i] = f(x[i]); }

	do {
		max = 0.; vspom = 0.;

		for (int i = 0; i < N; i++) {

			I = 0.;
			for (int j = 0; j < N - 1; j++) { I += 0.5*h*(K(x[i], x[j])*uk[j] + K(x[i], x[j + 1])*uk[j + 1]); }
			uk1[i] = f(x[i]) + I;

			vspom = fabs(uk1[i] - uk[i]);
			if (max < vspom) max = vspom;
		}

		copy(N, uk, uk1);
		k++;

	} while (max > eps);

	cout << k << " ��������\n";
	error(N, uk1);
	sol_to_file(N, x, uk1);

	delete[] uk1;
	delete[] uk;
	delete[] x;

}


void Degenerate_kernel(mytipe a, mytipe b, int N, int m)//m-����� ������� ��, ���
{
	mytipe h = (b - a) / (N - 1.);

	mytipe** alfa = new mytipe*[m];
	for (int i = 0; i < m; i++) alfa[i] = new mytipe[m];

	mytipe* betta = new mytipe[m];//������ �����

	mytipe* x = new mytipe[N];//���-�� ����� ��� ���������� ���������
	mytipe* u = new mytipe[N];//�������

	for (int i = 0; i < N; i++)  x[i] = a + i*h;

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			alfa[i][j] = 0.;
			for (int k = 0; k < N - 1; k++) alfa[i][j] -= 0.5*h*(psi(x[k], i)*phi(x[k], j) + psi(x[k + 1], i)*phi(x[k + 1], j));
		}

		alfa[i][i] += 1.;
		betta[i] = 0.;
		for (int k = 0; k < N - 1; k++) betta[i] += 0.5*h*(psi(x[k], i)*f(x[k]) + psi(x[k + 1], i)*f(x[k + 1]));
	}

	findx(m, betta, alfa);
	hod(m, betta, alfa);

	for (int i = 0; i < N; i++)
	{
		u[i] = f(x[i]);
		for (int j = 0; j < m; j++) u[i] += betta[j] * phi(x[i], j);
	}
	error(N, u);
	sol_to_file(N, x, u);

	for (int i = 0; i < m; i++)  delete[] alfa[i];
	delete[] alfa;
	delete[] betta;
	delete[] x;
	delete[] u;

}


void Degenerate_kernel_Cheb(mytipe a, mytipe b, int N, int m)//m-����� ������� ��, ���
{
	mytipe h = (b - a) / (N - 1.);

	mytipe** alfa = new mytipe*[m];
	for (int i = 0; i < m; i++) alfa[i] = new mytipe[m];

	mytipe* betta = new mytipe[m];//������ �����

	mytipe* x = new mytipe[N];//���-�� ����� ��� ���������� ���������
	mytipe* u = new mytipe[N];//�������

	for (int i = 0; i < N; i++)  x[i] = a + i*h;

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			alfa[i][j] = 0.;
			for (int k = 0; k < N - 1; k++) alfa[i][j] -= 0.5*h*(psiC(x[k], i)*phiC(x[k], j) + psiC(x[k + 1], i)*phiC(x[k + 1], j));
		}

		alfa[i][i] += 1.;
		betta[i] = 0.;
		for (int k = 0; k < N - 1; k++) betta[i] += 0.5*h*(psiC(x[k], i)*f(x[k]) + psiC(x[k + 1], i)*f(x[k + 1]));
	}

	findx(m, betta, alfa);
	hod(m, betta, alfa);

	for (int i = 0; i < N; i++)
	{
		u[i] = f(x[i]);
		for (int j = 0; j < m; j++) u[i] += betta[j] * phi(x[i], j);
	}
	error(N, u);
	sol_to_file(N, x, u);

	for (int i = 0; i < m; i++)  delete[] alfa[i];
	delete[] alfa;
	delete[] betta;
	delete[] x;
	delete[] u;

}



void Singular_kernel(int N)
{
	mytipe** A = new mytipe*[N + 1];
	for (int i = 0; i < N + 1; i++) { A[i] = new mytipe[N + 1]; }
	mytipe* d = new mytipe[N + 1];//������ �����

	mytipe h = 2.*Pi / N;

	mytipe* cx = new mytipe[N];
	mytipe* cy = new mytipe[N];
	mytipe* kx = new mytipe[N];//=nx
	mytipe* ky = new mytipe[N];//=ny



	for (int i = 0; i < N; i++) {
		cx[i] = cos(h*i); cy[i] = sin(h*i);
		kx[i] = cos(h*(i + 0.5)); ky[i] = sin(h*(i + 0.5)); d[i] = f(h*(i + 0.5));
	}

	d[N] = 0.; A[N][N] = 0.;

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {//nxQx+nyQy
			A[i][j] = 1. / N*(ky[i] * (kx[i] - cx[j]) - kx[i] * (ky[i] - cy[j])) / ((kx[i] - cx[j])*(kx[i] - cx[j]) + (ky[i] - cy[j])*(ky[i] - cy[j]));
		}
		A[i][N] = 1.; //��� �������
		A[N][i] = 1.;//��� ������
	}

	findx(N + 1, d, A);
	hod(N + 1, d, A);

	mytipe vspom = 0.;
	mytipe err = 0.;

	for (int i = 0; i < N; i++) {
		vspom = fabs(d[i] - 2 * sin(7 * i*h));
		if (err < vspom) err = vspom;
	}

	cout << "������ = " << err << endl;

	sol_to_file(N, cx, cy, d);

	for (int i = 0; i < N + 1; i++) { delete[] A[i]; }

	delete[] A;
	delete[] d;
	delete[] cx; delete[] cy;
	delete[] kx; delete[] ky;
}


void error(int N, mytipe* y)
{
	mytipe vspom = 0.;
	mytipe err = 0.;

	for (int i = 0; i < N; i++) {
		vspom = fabs(y[i] - 1.);
		if (err < vspom) {
			err = vspom;
			cout << "���� = " << i << "������ = " << err << endl;
		}
	}

	cout << "������ = " << err << endl;
}

bool findx(const unsigned int DIM, mytipe* b, mytipe** a) //������ ��� ������ � ��������� �� �������������
{
	mytipe max;     //��������������� ���������� ��� ��������� ��������� (���������� max)
	mytipe* vspom;  //��������������� ������ ��� ������������ ����� �������
	mytipe vspom2;   //��������������� ���������� ��� ������������ ��������� �������
	bool flag = true;  //���������� ����������(�������� ������������� �������)


	for (int k = 0, h; k < (DIM); k++)
	{
		max = fabs(a[k][k]);  //������ ������������� �������� � ���������� ���������
		vspom = a[k];   //������ ������ �� ������ ������������� ��������(���� ������� � ����� �� ���� �� ����������)
		h = k;   //������ ������ ������ (�� �������, ��������� ����)

		if (fabs(a[k][k]) < EPSILON) { flag = true; } //���� ������ ������� �������, �� ������������ ����
		else { flag = false; } //���� �� �������, �� �������������� ����

		for (int i = (k + 1); i < DIM; i++)
		{
			if (fabs(a[i][k]) > max && fabs(a[i][k])>EPSILON) //���� ������� ������ ����������� max � ������� � �� ����
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

void hod(const unsigned int DIM, mytipe* b, mytipe** a) //�������� ��� ������
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

mytipe phiC(mytipe x, int i) {

	if (i == 0) return 0.5*(1. - x);

	if (i == 1) return 0.25* pow(x, 3);

	if (i == 2) return -0.0208332* pow(x, 5);

	if (i == 3) return 0.000694148*pow(x, 7);

}



mytipe psiC(mytipe s, int i) {

	if (i == 0) return 1.;

	if (i == 1) return pow(s, 2);

	if (i == 2) return  pow(s, 4);

	if (i == 3) return pow(s, 6);
}

mytipe phi(mytipe x, int i) {

	if (i == 0) return 0.5*(1. - x);

	if (i == 1) return 0.25* pow(x, 3);

	if (i == 2) return -1. / 48.* pow(x, 5);

	if (i == 3) return  1. / 1440.*pow(x, 7);

}



mytipe psi(mytipe s, int i) {

	if (i == 0) return 1.;

	if (i == 1) return pow(s, 2);

	if (i == 2) return  pow(s, 4);

	if (i == 3) return pow(s, 6);
}



void copy(const unsigned int DIM, mytipe** a, mytipe** b) //����������� ������
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

void sol_to_file(int N, mytipe* x, mytipe* Y)
{
	ofstream fout;    // ������� ���������� ��� ������ � ����
	fout.open("sol.txt", ios_base::out | ios_base::trunc);

	for (unsigned int i = 0; i < N; i++) { fout << x[i] << " " << Y[i] << endl; }
	fout.close();
	fout.clear();
};

void sol_to_file(int N, mytipe* x, mytipe* y, mytipe* g)
{
	ofstream fout;    // ������� ���������� ��� ������ � ����
	fout.open("sol_SINGL.txt", ios_base::out | ios_base::trunc);

	for (unsigned int i = 0; i < N; i++) { fout << x[i] << " " << y[i] << " " << g[i] << endl; }
	fout.close();
	fout.clear();
};

void vivod(const unsigned int DIM, mytipe** const a) //��������� ������ �������
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

void vivod(const unsigned int DIM, const mytipe* const b) //��������� ������ �������
{
	cout << '(';
	for (int i = 0; i < DIM; i++)
	{
		cout << b[i];
		if (i < DIM - 1) { cout << ','; } //�������� ������
	}
	cout << ")^�\n";
}
