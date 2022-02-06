#include "iostream"
#include <fstream>
#include <typeinfo>
#include <cmath>

using namespace std;

const double EPSILON = 1e-8;   //��������� ��������� � 0

typedef double mytipe;  //��� ������, �������������� �� ���� ���������

mytipe *x_real;


void vivod(const unsigned int DIM, const mytipe* const b, mytipe** const a); //����� ������� � �������
void vivod(const unsigned int DIM, mytipe** const a);  //����� �������
void vivod(const unsigned int DIM, const mytipe* const b);  //����� �������
mytipe neviaz(const unsigned int DIM, const mytipe* const x, mytipe** const A, const mytipe* const b, const char k); //�������
void copy(const unsigned int DIM, mytipe** B, mytipe** const A); //����������� ������
void copy(const unsigned int DIM, mytipe* b, const mytipe* a); //����������� ��������
void trans(const unsigned int DIM, mytipe** a);  // ���������������� �������
void multimatrix(const unsigned int DIM, mytipe** const a, mytipe** const b); // ������������ ������ ��� �������� ����������
void mult_matr(const unsigned int DIM, mytipe** const A1, mytipe** const A2, mytipe** X); // ������������ ������ ��� �������� ����������
void mult_matrvect(const unsigned int DIM, mytipe** const A, const mytipe* const b, mytipe* x); // ������������ ������
void sum_vect(const unsigned int DIM, const mytipe* const b1, const mytipe* const b2, mytipe* x); //�������� ��������
void sum_matr(const unsigned int DIM, const mytipe** A1, const mytipe** A2, mytipe** X); //�������� ������
void mult_matrnum(const unsigned int DIM, mytipe** const A, const mytipe a, mytipe** X); //�������� ������� �� �����
void mult_vectnum(const unsigned int DIM, const mytipe* const b, const mytipe a, mytipe* x); //��������� ������� �� �����
mytipe norm(const unsigned int DIM, mytipe** const A, const char flag); // 1 �� 4 ���� �������
mytipe norm(const unsigned int DIM, const mytipe* const b, const char flag);  // 1 �� 3 ���� �������
void singlematr(const unsigned int DIM, mytipe** const E);  //���������� ������� ���������
void matrvec(const unsigned int DIM, mytipe** A, mytipe* b);  //��������� ������� �� �������

mytipe SimpleIter_C(const unsigned int DIM, mytipe** const A, const mytipe* const b, const mytipe TAU, const char); //���������� ������� � � ������ ������� ��������
mytipe Yakobi_C(const unsigned int DIM, mytipe** const A, const mytipe* const b, const char); //���������� ������� � � ������ �����
mytipe Seidel_C(const unsigned int DIM, mytipe** const A, const mytipe* const b, const char); //���������� ������� � � ������ �������
mytipe Relax_C(const unsigned int DIM, mytipe** const A, const mytipe* const b, const mytipe W, const char); //���������� ������� � � ������ ����������

void SimpleIter(const unsigned int DIM, mytipe** const A, const mytipe* const b, const mytipe TAU, const char NORM, mytipe EPS, const mytipe normC, mytipe* xk); //����� ������� ��������
void Yakobi(const unsigned int DIM, mytipe** const A, const mytipe* const b, const char NORM, mytipe EPS, const mytipe normC, mytipe* xk); //����� �����
void Seidel(const unsigned int DIM, mytipe** const A, const mytipe* const b, const char NORM, mytipe EPS, const mytipe normC, mytipe* xk); //����� �������
void Relax(const unsigned int DIM, mytipe** const A, const mytipe* const b, const mytipe W, const char NORM, mytipe EPS, const mytipe normC, mytipe* xk); //����� ����������

void C_Iter(const unsigned int DIM, mytipe** const A, const mytipe* const b, mytipe** const C, const mytipe* const y, const char NORM); //������� ������������ �������, ��������� � � y

void relax3d(const int var, const char NORM, const mytipe w, const mytipe EPS);

void C_l_d_u(mytipe** C, const int DIM, const char NORM);

int main() {
	setlocale(LC_ALL, "Russian"); //���������� ������� ����

	ifstream fin;    // ������� ���������� ��� ���������� �� �����
	fin.open("1.txt", ios_base::in | ios_base::app | ios_base::binary);  //������� ����
	unsigned int a;  // ������ ������������
	fin >> a;  //������� �� ����� ������ ������������
	const unsigned int DIM1 = a;  //������ ������������� ������� �����������


	mytipe **A;  //�������, ��������� � �����
	A = new mytipe *[DIM1];            // |��������
	for (int i = 0; i < DIM1; i++) {   // |����� ���
		A[i] = new mytipe[DIM1];       // |������� �
	}

	mytipe *b;  //�������, ���������� �� �����
	b = new mytipe[DIM1];     //|�������� ����� ��� ������� b


	for (int i = 0; i < DIM1; i++)  //��������� �� ����� ������� � �������
	{
		for (int j = 0; j < DIM1; j++)
		{
			fin >> A[i][j];  //�������
		}
		fin >> b[i];  //�������
	}

	fin.close();  //��������� ����
	fin.clear();  //���������� ���� �� �������� ����������

	cout << "Dim=" << a << endl;
	vivod(DIM1, b, A);    //������� �������

	const char NORM = 'k';
	const mytipe EPS = 1e-3;

	mytipe *x;  //������� ��������� �����������
	x = new mytipe[DIM1];


	x_real = new mytipe[DIM1];
	x_real[0] = 5;
	x_real[1] = -7;
	x_real[2] = 12;
	x_real[3] = 4;

	//����� ������� ��������
	for (int i = 0; i < DIM1; i++) { x[i] = 0.0; }
	const mytipe TAU = 0.06;
	mytipe normC_SimplIter = SimpleIter_C(DIM1, A, b, TAU, NORM);
	SimpleIter(DIM1, A, b, TAU, NORM, EPS, normC_SimplIter, x);
	cout << "������ x = \n"; vivod(DIM1, x);
	cout << "����� ������� = " << neviaz(DIM1, x, A, b, NORM) << endl;



	// ����� �����	
     for (int i = 0; i < DIM1; i++) { x[i] = 0.0; }
	mytipe normC_Yakobi = Yakobi_C(DIM1, A, b, NORM);
	Yakobi(DIM1, A, b, NORM, EPS, normC_Yakobi, x);
	cout << "������ x = \n"; vivod(DIM1, x);
	cout << "����� ������� = " << neviaz(DIM1, x, A, b, NORM) << endl;


	// ����� �������
	for (int i = 0; i < DIM1; i++) { x[i] = 0.0; }
	mytipe normC_Seidel = Seidel_C(DIM1, A, b, NORM);
	Seidel(DIM1, A, b, NORM, EPS, normC_Seidel, x);
	cout << "������ x = \n"; vivod(DIM1, x);
	cout << "����� ������� = " << neviaz(DIM1, x, A, b, NORM) << endl;


	//����� ����������
	for (int i = 0; i < DIM1; i++) { x[i] = 0.0; }
	const mytipe W = 0.9;
	mytipe normC_Relax = Relax_C(DIM1, A, b, W, NORM);
	Relax(DIM1, A, b, W, NORM, EPS, normC_Relax, x);
	cout << "������ x = \n"; vivod(DIM1, x);
	cout << "����� ������� = " << neviaz(DIM1, x, A, b, NORM) << endl;


	//3-x ������������ �������
	relax3d(14, NORM, 1.9, EPS);


	delete[] b; //������� ������������ �������

	for (int i = 0; i < DIM1; i++) //������� ������������ ������ �������
	{
		delete[] A[i];
	}
	delete[] A;//������� ������������ �������
	delete[]x; delete[]x_real;

	cin.get();
	return 0;
}

void vivod(const unsigned int DIM, const mytipe* const b, mytipe** const a) //��������� ������ ������� � �������
{
	cout << endl;
	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			cout << a[i][j] << " ";
		}
		cout << '|' << b[i] << endl;
	}
	cout << endl;
}

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




mytipe neviaz(const unsigned int DIM, const mytipe* const x, mytipe** const A, const mytipe* const b, const char k) //�������
{
	mytipe* b1;          //��������� ������ ���
	b1 = new mytipe[DIM];  //������, ������� ��������� ��� ����������� ������� � �������
	mytipe nor;  //����� ������� ��������

	for (int i = 0; i < DIM; i++)
	{
		b1[i] = 0;
		for (int j = 0; j < DIM; j++)
		{
			b1[i] += A[i][j] * x[j]; //�������, ������������ �������, ������� ������ �����

		}
		b1[i] -= b[i]; //�������� i-�� ���������� ������������� ������� � ���������
	}

	nor = norm(DIM, b1, k);
	delete[] b1; //������������ ������ ������������� �������
	return nor;
}


void copy(const unsigned int DIM, mytipe** B, mytipe** const A) //����������� ������
{

	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			B[i][j] = A[i][j];
		}
	}
}

void copy(const unsigned int DIM, mytipe* b, const mytipe* a) //����������� ��������
{
	for (int j = 0; j < DIM; j++)
	{
		b[j] = a[j];
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

void multimatrix(const unsigned int DIM, mytipe** const a, mytipe** const b) //������������ ������
{
	double t;
	mytipe** E;
	E = new mytipe *[DIM];              // �������� �����
	for (int i = 0; i < DIM; i++)       // ��� �������
	{                                   //  ���������� �����������
		E[i] = new mytipe[DIM];         //  ������������ ������ a � b
	}

	for (int i = 0; i<DIM; i++) {
		for (int l = 0; l<DIM; l++) {
			t = 0;
			for (int j = 0; j<DIM; j++) {
				t += a[i][j] * b[j][l];  //������������ ������ ������� a �� ������� ������� b
			}
			E[i][l] = t;   //������ ���������� ������������ � �������������� �������
		}
	}

	vivod(DIM, E);  //������� ���������

	for (int i = 0; i < DIM; i++) {   //������� ������������ ������ �������
		delete[] E[i];
	}
	delete[] E;             // ������� ������������ �������
}

void mult_matr(const unsigned int DIM, mytipe** const A1, mytipe** const A2, mytipe** X) //������������ ������
{
	for (int i = 0; i<DIM; i++) {
		for (int l = 0; l<DIM; l++) {
			X[i][l] = 0;
			for (int j = 0; j<DIM; j++) {
				X[i][l] += A1[i][j] * A2[j][l];  //������������ ������ ������� a �� ������� ������� b
			}
		}
	}
}

void mult_matrvect(const unsigned int DIM, mytipe** const A, const mytipe* const b, mytipe* x) //��������� ������� �� ������
{
	for (int i = 0; i<DIM; i++) {
		x[i] = 0;
		for (int j = 0; j < DIM; j++) {
			x[i] += A[i][j] * b[j];
		}
	}
}

void mult_matrnum(const unsigned int DIM, mytipe** const A, const mytipe a, mytipe** X) //�������� ������� �� �����
{
	for (int i = 0; i<DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			X[i][j] = a * A[i][j];
		}
	}
}

void mult_vectnum(const unsigned int DIM, const mytipe* const b, const mytipe a, mytipe* x) //��������� ������� �� �����
{
	for (int i = 0; i<DIM; i++) {
		x[i] = a * b[i];
	}
}

void sum_vect(const unsigned int DIM, const mytipe* const b1, const mytipe* const b2, mytipe* x) //�������� ��������
{
	for (int i = 0; i < DIM; i++) {
		x[i] = b1[i] + b2[i];
	}
}

void sum_matr(const unsigned int DIM, const mytipe** A1, const mytipe** A2, mytipe** X) //�������� ������
{
	for (int i = 0; i<DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			X[i][j] = A1[i][j] + A2[i][j];
		}
	}
}

mytipe norm(const unsigned int DIM, mytipe** const A, const char flag) //4 ����� ������� � �������� �������� �����
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
				vspom += abs(A[j][i]); //��������� �������� i-�� ������
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
				vspom += abs(A[j][i]);  //��������� �������� j-�� ������
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

		for (int j = 0; j < DIM; j++)
		{
			for (int i = 0; i<DIM; i++)
			{
				if (abs(A[j][i]) > norm) { norm = abs(A[j][i]); } //������� ������������ �������
			}
		}
		return DIM*norm;
	}
}

mytipe norm(const unsigned int DIM, const mytipe* const b, const char flag) //3 ����� ������� � �������� �������� �����
{
	mytipe norm = 0;  //���������� ����� �������
	if (flag == 'k')  //���������� ����� (max)
	{

		for (int i = 0; i < DIM; i++)
		{
			if (norm < abs(b[i])) { norm = abs(b[i]); }  //������� ��������
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


void singlematr(const unsigned int DIM, mytipe** const E) //���������� ������� ���������
{
	for (int i = 0; i < DIM; i++)      //|������
	{									//|������� �
		for (int j = 0; j < DIM; j++)  //|���������
		{							    //|���
			E[i][j] = 0;				//|����������
		}								//|��������������
		E[i][i] = 1;					//|
	}
}


void matrvec(const unsigned int DIM, mytipe** A, mytipe* b) //��������� ������� �� �������
{
	mytipe* c;
	c = new mytipe[DIM];
	for (int i = 0; i < DIM; i++)
	{
		c[i] = 0;
		for (int j = 0; j < DIM; j++)
		{
			c[i] += A[i][j] * b[j];
		}
	}
	copy(DIM, b, c);
	delete[] c;
}

mytipe SimpleIter_C(const unsigned int DIM, mytipe** const A, const mytipe* const b, const mytipe tau, const char NORM) //���������� ������� � � ������ ������� ��������
{
	mytipe **C;  //������� �� ������ ���� ��������
	C = new mytipe *[DIM];
	for (int i = 0; i < DIM; i++) {
		C[i] = new mytipe[DIM];
	}

	mytipe *y = new mytipe[DIM];           //������ ����������� �� (�+1)-�� ����

	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			C[i][j] = -tau*A[i][j];
		}
		C[i][i] += 1;
	}
	copy(DIM, y, b);
	mult_vectnum(DIM, y, tau, y);

	cout << "\n������� C=";
	vivod(DIM, C);  //������� ������� C
	cout << "������ y=\n";
	vivod(DIM, y);
	mytipe norma = norm(DIM, C, NORM);
	cout << "\n����� ������� C=" << norma << endl;

	//C_Iter(DIM, A, b, C, y, NORM);

	for (int i = 0; i < DIM; i++) //������� ������������ ������ �������
	{
		delete[] C[i];
	}

	delete[] C;//������� ������������ �������
	delete[] y;
	return norma;
}

mytipe Yakobi_C(const unsigned int DIM, mytipe** const A, const mytipe* const b, const char NORM) //���������� ������� � � ������ �����
{
	mytipe **C;  //������� �� ������ ���� ��������
	C = new mytipe *[DIM];
	for (int i = 0; i < DIM; i++) {
		C[i] = new mytipe[DIM];
	}

	mytipe *y = new mytipe[DIM];           //������ ����������� �� (�+1)-�� ����


	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			C[i][j] = -A[i][j] / A[i][i];
		}
		C[i][i] += 1;
	}
	for (int i = 0; i < DIM; i++)
	{
		y[i] = b[i] / A[i][i];
	}

	cout << "\n������� C=";
	vivod(DIM, C);  //������� ������� C
	cout << "������ y=\n";
	vivod(DIM, y);
	mytipe norma = norm(DIM, C, NORM);
	cout << "\n����� ������� C=" << norma << endl;

	//C_Iter(DIM, A, b, C, y, NORM);

	for (int i = 0; i < DIM; i++) //������� ������������ ������ �������
	{
		delete[] C[i];
	}

	delete[] C;//������� ������������ �������
	delete[] y;
	return norma;
}

mytipe Seidel_C(const unsigned int DIM, mytipe** const A, const mytipe* const b, const char NORM) //���������� ������� � � ������ �������
{
	mytipe** vspom = new mytipe *[DIM];;                  //�������� ������ ��� �������� ������� (A1+D)
	mytipe **C = new mytipe *[DIM];  //������� �� ������ ���� ��������
	for (int i = 0; i < DIM; i++) { vspom[i] = new mytipe[DIM]; C[i] = new mytipe[DIM]; }

	mytipe *y = new mytipe[DIM];           //������ ����������� �� (�+1)-�� ����

	mytipe vspom1;
	for (int j = 0; j < DIM; j++) {
		vspom[j][j] = 1 / A[j][j];
		for (int i = j + 1; i < DIM; i++) {
			vspom1 = 0;
			for (int k = j; k < i; k++) {
				vspom1 -= vspom[k][j] * A[i][k];
			}
			vspom[i][j] = vspom1 / A[i][i];
		}

		for (int i = 0; i < j; i++) {
			vspom[i][j] = 0;
		}
	}

	for (int i = 0; i < DIM; i++) {  // -(A1+D)^(-1)*A2
		for (int j = 1; j < DIM; j++) {
			C[i][j] = 0;
			for (int k = 0; k < j; k++) {
				C[i][j] -= vspom[i][k] * A[k][j];
			}
		}
		C[i][0] = 0;
	}
	for (int i = 0; i < DIM; i++) {
		y[i] = 0;
		for (int j = 0; j <= i; j++) {
			y[i] += vspom[i][j] * b[j];
		}
	}
	cout << "\n������� C=";
	vivod(DIM, C);  //������� ������� C
	cout << "������ y=\n";
	vivod(DIM, y);
	mytipe norma = norm(DIM, C, NORM);
	cout << "\n����� ������� C=" << norma << endl;

	//C_l_d_u(C, DIM, NORM);
	//C_Iter(DIM, A, b, C, y, NORM);

	for (int i = 0; i < DIM; i++) //������� ������������ ������ �������
	{
		delete[] C[i]; delete[] vspom[i];
	}

	delete[] C;//������� ������������ �������
	delete[] y;
	return norma;
}

mytipe Relax_C(const unsigned int DIM, mytipe** const A, const mytipe* const b, const mytipe W, const char NORM) //���������� ������� � � ������ ����������
{
	mytipe** vspom = new mytipe *[DIM];;                  //�������� ������ ��� �������� ������� (A1+D)
	mytipe **C = new mytipe *[DIM];  //������� �� ������ ���� ��������
	for (int i = 0; i < DIM; i++) { vspom[i] = new mytipe[DIM]; C[i] = new mytipe[DIM]; }

	mytipe *y = new mytipe[DIM];           //������ ����������� �� (�+1)-�� ����
	mytipe vspom1;
	for (int j = 0; j < DIM; j++) {
		vspom[j][j] = 1 / A[j][j];
		for (int i = j + 1; i < DIM; i++) {
			vspom1 = 0;
			for (int k = j; k < i; k++) {
				vspom1 -= vspom[k][j] * A[i][k];
			}
			vspom[i][j] = vspom1 / A[i][i] * W;
		}

		for (int i = 0; i < j; i++) {
			vspom[i][j] = 0;
		}
	}

	for (int i = 0; i < DIM; i++) {  // -(A1+D)^(-1)*A2
		for (int j = 0; j < DIM; j++) {
			C[i][j] = 0;
			for (int k = 0; k < j; k++) {
				C[i][j] -= vspom[i][k] * A[k][j];
			}
			C[i][j] = W*C[i][j] + (1 - W)*A[j][j] * vspom[i][j];
		}
	}
	for (int i = 0; i < DIM; i++) {
		y[i] = 0;
		for (int j = 0; j <= i; j++) {
			y[i] += vspom[i][j] * b[j];
		}
		y[i] = W*y[i];
	}

	cout << "\n������� C=";
	vivod(DIM, C);  //������� ������� C
	cout << "������ y=\n";
	vivod(DIM, y);
	mytipe norma = norm(DIM, C, NORM);
	cout << "\n����� ������� C=" << norma << endl;

	//C_l_d_u(C,DIM,NORM);
	//C_Iter(DIM, A, b, C, y, NORM);

	for (int i = 0; i < DIM; i++) //������� ������������ ������ �������
	{
		delete[] C[i]; delete[] vspom[i];
	}

	delete[] C;//������� ������������ �������
	delete[] y;
	return norma;
}

void SimpleIter(const unsigned int DIM, mytipe** const A, const mytipe* const b, const mytipe TAU, const char NORM, mytipe EPS, mytipe normC, mytipe* xk) //����� ������� ��������
{
	cout << "\n\t\t\t----------------------\n\t\t\t����� ������� ��������\n\t\t\t----------------------\n ";
	mytipe* xk1 = new mytipe[DIM];
	unsigned int iter = 0;
	mytipe* V = new mytipe[DIM];
	EPS *= fabs(1 - normC) / normC;
	mytipe norm_new;

	do {
	
		iter++;
		for (int i = 0; i < DIM; i++)
		{
			xk1[i] = 0.0;
			for (int j = 0; j < DIM; j++)
			{
				xk1[i] += A[i][j] * xk[j];
			};
			xk1[i] = xk[i] + TAU*(-xk1[i] + b[i]);
		}

		for (int i = 0; i < DIM; i++) { V[i] = x_real[i] - xk[i]; }
		cout << "����� ����������� �� " << iter << "-�� ��������:  " << norm(DIM, V, NORM) << endl;
		for (int i = 0; i < DIM; i++) { V[i] = xk1[i] - xk[i]; }
		norm_new = norm(DIM, V, NORM);
		copy(DIM, xk, xk1);
	} while (norm_new > EPS);
	
	
	cout << "����� ������� �������� �������� \n ���������� �������� = " << iter << endl; 
	delete[] xk1;
	delete[] V;
}

void Yakobi(const unsigned int DIM, mytipe** const A, const mytipe* const b, const char NORM, mytipe EPS, const mytipe normC, mytipe* xk) //����� �����
{
	cout << "\n\t\t\t-----------\n\t\t\t����� �����\n\t\t\t-----------\n ";
	mytipe* xk1 = new mytipe[DIM];
	unsigned int iter = 0;
	//while (neviaz(DIM, xk, A, b, NORM) > EPS) {
	mytipe* V = new mytipe[DIM];
	EPS *= (1 - normC) / normC;
	mytipe norm_new = norm(DIM, V, NORM);

	do {
		iter++;
		for (unsigned int j = 0; j < DIM; j++) {
			xk1[j] = b[j];
			for (unsigned int i = 0; i < j; i++) {
				xk1[j] -= A[j][i] * xk[i];
			}
			for (unsigned int i = j + 1; i < DIM; i++) {
				xk1[j] -= A[j][i] * xk[i];
			}
			xk1[j] = xk1[j] / A[j][j];
		}

		for (int i = 0; i < DIM; i++) { V[i] = x_real[i] - xk[i]; }
		cout << "����� ����������� �� " << iter << "-�� ��������:  " << norm(DIM, V, NORM) << endl;

		for (int i = 0; i < DIM; i++) { V[i] = xk1[i] - xk[i]; }
		norm_new = norm(DIM, V, NORM);
		copy(DIM, xk, xk1);
	} while (norm_new > EPS);
	cout << "����� ����� �������� \n ���������� �������� = " << iter << endl;
	delete[] xk1;
	delete[] V;
}

void Seidel(const unsigned int DIM, mytipe** const A, const mytipe* const b, const char NORM, mytipe EPS, const mytipe normC, mytipe* xk) //����� �������
{
	cout << "\n\t\t\t-------------\n\t\t\t����� �������\n\t\t\t-------------\n ";
	mytipe* xk1 = new mytipe[DIM];
	unsigned int iter = 0;
	mytipe* V = new mytipe[DIM];
	EPS *= (1 - normC) / normC;
	mytipe norm_new;

	do {
		iter++;
		for (unsigned int j = 0; j < DIM; j++) {
			xk1[j] = b[j];
			for (unsigned int i = 0; i < j; i++) {
				xk1[j] -= A[j][i] * xk1[i];
			}
			for (unsigned int i = j + 1; i < DIM; i++) {
				xk1[j] -= A[j][i] * xk[i];
			}
			xk1[j] = xk1[j] / A[j][j];
		}
		for (int i = 0; i < DIM; i++) { V[i] = x_real[i] - xk[i]; }
		cout << "����� ����������� �� " << iter << "-�� ��������:  " << norm(DIM, V, NORM) << endl;

		for (int i = 0; i < DIM; i++) { V[i] = xk1[i] - xk[i]; }
		norm_new = norm(DIM, V, NORM);
		copy(DIM, xk, xk1);
	} while (norm_new > EPS);

	cout << "����� ������� �������� \n ���������� �������� = " << iter << endl;
	delete[] xk1;
	delete[] V;
}

void Relax(const unsigned int DIM, mytipe** const A, const mytipe* const b, const mytipe W, const char NORM, mytipe EPS, const mytipe normC, mytipe* xk) //����� ����������
{
	cout << "\n\t\t\t----------------\n\t\t\t����� ����������\n\t\t\t----------------\n ";
	mytipe* xk1 = new mytipe[DIM];
	unsigned int iter = 0;
	mytipe vspom;
	mytipe* V = new mytipe[DIM];
	EPS *= (1 - normC) / normC;
	mytipe norm_new;


	do {
		//while (neviaz(DIM, xk, A, b, NORM) > EPS) {
		iter++;
		for (unsigned int j = 0; j < DIM; j++) {
			xk1[j] = (b[j] / A[j][j] - xk[j])*W + xk[j];
			vspom = 0;
			for (unsigned int i = 0; i < j; i++) {
				vspom -= A[j][i] * xk1[i];
			}
			xk1[j] += vspom*W / A[j][j];
			vspom = 0;
			for (unsigned int i = j + 1; i < DIM; i++) {
				vspom -= A[j][i] * xk[i];
			}
			xk1[j] += vspom*W / A[j][j];
		}
		for (int i = 0; i < DIM; i++) { V[i] = x_real[i] - xk[i]; }
		cout << "����� ����������� �� " << iter << "-�� ��������:  " << norm(DIM, V, NORM) << endl;
		for (int i = 0; i < DIM; i++) { V[i] = xk1[i] - xk[i]; }
		norm_new = norm(DIM, V, NORM);
		copy(DIM, xk, xk1);
	} while (norm_new > EPS);

	cout << "����� ���������� �������� \n ���������� �������� = " << iter << endl;
	delete[] xk1;
	delete[] V;
}

void C_Iter(const unsigned int DIM, mytipe** const A, const mytipe* const b, mytipe** const C, const mytipe* const y, const char NORM) //������� ������������ �������, ��������� � � y
{
	cout << "\n\t\t\t------------\n\t\t\tx_k+1=Cx_k+y\n\t\t\t------------\n ";
	mytipe EPS = 1e-3;
	unsigned int n = 0;
	mytipe* xk1 = new mytipe[DIM];
	mytipe* xk = new mytipe[DIM]; for (int i = 0; i < DIM; i++) { xk1[i] = 0.0; xk[i] = 0.0; }
	while (neviaz(DIM, xk, A, b, NORM) > EPS) {
		n++;
		mult_matrvect(DIM, C, xk, xk1);
		sum_vect(DIM, xk1, y, xk1);
		n++;
		mult_matrvect(DIM, C, xk1, xk);
		sum_vect(DIM, xk, y, xk);
	};
	cout << "\nx = ";
	vivod(DIM, xk1);
	cout << n << " ��������\n";
	delete[] xk1;
}

void relax3d(const int var, const char NORM, const mytipe w, const mytipe EPS)
{
	const unsigned int DIM = 200 + var;  //������ ������������� ������� �����������
	int iter = 0;
	mytipe *a;  //�������, ���������� �� �����
	a = new mytipe[DIM];     //|�������� ����� ��� ������� b
	mytipe *b;  //�������, ���������� �� �����
	b = new mytipe[DIM];     //|�������� ����� ��� ������� b
	mytipe *c;  //�������, ���������� �� �����
	c = new mytipe[DIM];     //|�������� ����� ��� ������� b
	mytipe *d;  //�������, ���������� �� �����
	d = new mytipe[DIM];     //|�������� ����� ��� ������� b
	mytipe *xk1;                //������ ����������� �� (�+1)-�� ����
	xk1 = new mytipe[DIM];
	mytipe *xk;                //������ ����������� �� (�)-�� ����
	xk = new mytipe[DIM];//��������� �����������

	d[0] = 6.0; a[0] = 0.0;

	for (int i = 0; i < DIM - 1; i++) {
		xk[i] = 0.0; xk1[i] = 0.0;
		a[i + 1] = 1.0; b[i] = 4.0; c[i] = 1.0; d[i + 1] = 10.0 - 2.0*((i + 2) % 2);
	}

	b[DIM - 1] = 4.0;
	d[DIM - 1] = 9.0 - 3.0*((DIM) % 2);
	c[DIM - 1] = 0.0;
	xk[DIM - 1] = 0.0;
	xk1[DIM - 1] = 0.0;

	mytipe *V = new mytipe[DIM];//������ ����������� �� (�)-�� ����

	do
	{
		copy(DIM, xk, xk1);
		xk1[0] = w*(d[0] - c[0] * xk[1]) / b[0] + (1.0 - w)*xk[0];
		for (int j = 1; j < DIM - 1; j++) { xk1[j] = w*(d[j] - a[j] * xk1[j - 1] - c[j] * xk[j + 1]) / b[j] + (1 - w)*xk[j]; }
		xk1[DIM - 1] = w*(d[DIM - 1] - a[DIM - 1] * xk1[DIM - 2]) / b[DIM - 1] + (1.0 - w)*xk[DIM - 1];

		for (int i = 0; i < DIM; i++) { V[i] = xk1[i] - xk[i]; }
		iter++;

	} while (norm(DIM, V, NORM) > EPS*(2.0 - w) / w);//(norm(DIM1, V, NORM) / (norm(DIM1, xk, NORM) + EPSILON) > EPS);

	cout << "\n\t\t\t------------------------\n\t\t\t3-� ������������ �������\n\t\t\t------------------------\n\n ";
	cout << "\n\t\t\t\tDim = " << DIM;
	cout << "\n\na = ";
	vivod(DIM, a);
	cout << "\nb = ";
	vivod(DIM, b);
	cout << "\nc = ";
	vivod(DIM, c);
	cout << "\nd = ";
	vivod(DIM, d);
	cout << "\nx = ";
	vivod(DIM, xk1);
	cout << "\n\t\t\t\t" << iter << " ��������.\n";

	delete[] a; //������� ������������ �������
	delete[] b; //������� ������������ �������
	delete[] c; //������� ������������ �������
	delete[] d; //������� ������������ �������
	delete[] xk; //������� ������������ �������
	delete[] xk1; //������� ������������ �������
	delete[] V; //������� ������������ �������
}

void C_l_d_u(mytipe** C, const int DIM, const char NORM)
{
	mytipe **Cl = new mytipe *[DIM];   //�������, ��������� � �����
	mytipe **Cd = new mytipe *[DIM];   //�������, ��������� � �����
	mytipe **Cu = new mytipe *[DIM];   //�������, ��������� � �����
	for (int i = 0; i < DIM; i++) {   // |����� ���
		Cl[i] = new mytipe[DIM];       // |������� �
		Cd[i] = new mytipe[DIM];
		Cu[i] = new mytipe[DIM];
	}
	for (int i = 0; i < DIM; i++)
	{
		for (int j = i + 1; j < DIM; j++)
		{
			Cl[j][i] = C[j][i]; Cu[i][j] = C[i][j];
			Cl[i][j] = 0.0; Cu[j][i] = 0.0; Cd[i][j] = 0.0; Cd[j][i] = 0.0;
		}
		Cl[i][i] = 0.0; Cu[i][i] = 0.0; Cd[i][i] = C[i][i];
	}

	mytipe cl, cu;
	cout << "\n������� Cl=";
	vivod(DIM, Cl);  //������� ������� C

	cl = norm(DIM, Cl, NORM);
	cout << "����� ������� Cl=" << cl << endl;

	cout << "\n������� Cd=";
	vivod(DIM, Cd);  //������� ������� C

	cout << "����� ������� Cd=" << norm(DIM, Cd, NORM) << endl;

	cout << "\n������� Cu=";
	vivod(DIM, Cu);  //������� ������� C

	cu = norm(DIM, Cu, NORM);
	cout << "����� ������� Cu=" << cu << "\n\n||Cl||+||Cu|| = " << cu + cl;

	for (int i = 0; i < DIM; i++) //������� ������������ ������ �������
	{
		delete[] Cl[i]; delete[] Cd[i]; delete[] Cu[i];
	}
	delete[] Cl; delete[] Cd; delete[] Cu;
}

