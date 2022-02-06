#include "iostream"
#include <fstream>
#include <typeinfo>
#include <cmath>

using namespace std;

const double EPSILON = 10e-8;   //��������� ��������� � 0
typedef double mytipe;  //��� ������, �������������� �� ���� ���������
mytipe **A; //�������, ��������� � �����
mytipe *B; //�������, ���������� �� �����


void vivod(const unsigned int DIM, const mytipe* const b, mytipe** const a); //����� ������� � �������
void vivod(const unsigned int DIM, mytipe** const a);  //����� �������
void vivod(const unsigned int DIM, const mytipe* const b);  //����� �������
void hod(const unsigned int DIM, mytipe* b, mytipe** a); //�������� ��� ������
bool findx(const unsigned int DIM, mytipe* b, mytipe** a);  //������ ��� ������ � ������� ������������� �������
void neviaz(const unsigned int DIM, const mytipe* const x);  //�������
void copy(const unsigned int DIM, mytipe* b, mytipe** a); //����������� ������� � �������
void copy(const unsigned int DIM, mytipe** a, mytipe** b); //����������� ������
void copy(const unsigned int DIM, mytipe* a, mytipe* b); //����������� ��������
void trans(const unsigned int DIM, mytipe** a);  // ���������������� �������
void multimatrix(const unsigned int DIM, mytipe** const a, mytipe** const b); // ������������ ������ ��� �������� ����������
mytipe norm(const unsigned int DIM, mytipe** const A, char flag); // 1 �� 4 ���� �������
mytipe norm(const unsigned int DIM, const mytipe* const b, char flag);  // 1 �� 3 ���� �������
void obusl(const unsigned int DIM, const mytipe* const x); //������ ����� ����� ���������������
void singlematr(const unsigned int DIM, mytipe** const E);  //���������� ������� ���������
void QR(const unsigned int DIM, mytipe** T, mytipe** R); //QR-��������
void matrvec(const unsigned int DIM, mytipe** A, mytipe* b);  //��������� ������� �� �������
int main() {
	setlocale(LC_ALL, "Russian"); //���������� ������� ����
	ifstream fin;    // ������� ���������� ��� ���������� �� �����
	fin.open("6(1).txt", ios_base::in | ios_base::app | ios_base::binary);  //������� ����
	unsigned int a;  // ������ ������������
	fin >> a;  //������� �� ����� ������ ������������
	cout << "Dim=" << a << endl;
	const unsigned int DIM1 = a;  //������ ������������� ������� �����������

	mytipe **ary;                      // |������� 
	ary = new mytipe *[DIM1];          // |������������ 
	for (int i = 0; i < DIM1; i++) {   // |������ 
		ary[i] = new mytipe[DIM1];     // |������� �
	}

	mytipe *b;                //|�������
	b = new mytipe[DIM1];     //|������ �����


	A = new mytipe *[DIM1];            // |��������  
	for (int i = 0; i < DIM1; i++) {   // |����� ��� 
		A[i] = new mytipe[DIM1];       // |������� �
	}

	B = new mytipe[DIM1];     //|�������� ����� ��� ������� b

	for (int i = 0; i < DIM1; i++)  //��������� �� ����� ������� � �������
	{
		for (int j = 0; j < DIM1; j++)
		{
			fin >> A[i][j];  //�������
		}
		fin >> B[i];  //�������
	}

	fin.close();  //��������� ����
	fin.clear();  //���������� ���� �� �������� ����������

	copy(DIM1, b, ary);  //����������� �������� ������� � �������

	vivod(DIM1, b, ary);  //������� ������� � �������


	bool flag = true;  //��� ����������� ������������� �������

	flag = findx(DIM1, b, ary);  //������ ��� ������ � ���������� ������������� �������

	if (flag)  //���� ������� �� ���������
	{
		vivod(DIM1, b, ary); //������� ����������������� �������


		hod(DIM1, b, ary);  // �������� ��� ������


		vivod(DIM1, b);  //������� ������� �����������

		neviaz(DIM1, b); //����� �������
	}


	//            !!!!!! QR-������ !!!!!!

	flag = false;

	mytipe **R;							// |��������
	R = new mytipe *[DIM1];             // |�����
	for (int i = 0; i < DIM1; i++) {    // |��� �������
		R[i] = new mytipe[DIM1];        // |R
	}


	copy(DIM1, b, R);   //����������� �������� ������� � �������

	vivod(DIM1, b, R);  //����� ������� � �������



	mytipe **T;                       // |������� ������������� ������� � �������� ������� = Q
	T = new mytipe *[DIM1];           // |
	for (int i = 0; i < DIM1; i++) {  // |
		T[i] = new mytipe[DIM1];      // |
	}	                               //|

	singlematr(DIM1, T);

	QR(DIM1, T, R);

	if (abs(R[DIM1 - 1][DIM1 - 1]) < EPSILON) { flag = true; cout << "������� ���������\n\n"; }
	else { matrvec(DIM1, T, b); }

	trans(DIM1, T); //�������� Q �� T, ������������

	cout << "Q=" << endl;

	vivod(DIM1, T);  //������� Q

	cout << "R=" << endl;

	vivod(DIM1, R); //������� R

	cout << "A=Q*R=" << endl;

	multimatrix(DIM1, T, R);
	if (!flag)  //���� ������� �����������
	{
		hod(DIM1, b, R);  //�������� ���(���������� ������� �)

		vivod(DIM1, b);  //����� ������� b


						 //!!!!!!!!!!!!!!������� � ������ ����� ����� ���������������!!!!!!!!!!

		obusl(DIM1, b);



		//!!!!!!!!!!!!!!����� ���������������!!!!!!!!!!!


		mytipe** inverseA;                  //�������� ������ ��� �������� �������
		inverseA = new mytipe *[DIM1];
		for (int i = 0; i < DIM1; i++) { inverseA[i] = new mytipe[DIM1]; }
		trans(DIM1, T);
		singlematr(DIM1, inverseA);

		for (int r = 0; r < DIM1; r++)
		{
			matrvec(DIM1, T, inverseA[r]);
			hod(DIM1, inverseA[r], R); //������� ������ ������� �������� �������
		}

		copy(DIM1, ary, A); //�������� ������ ������� � �������
		trans(DIM1, inverseA); //������������� �������� �������
		cout << "A^-1=";
		vivod(DIM1, inverseA);  //������� �������� �������
		cout << "A*A^-1=";
		multimatrix(DIM1, inverseA, ary); //����������� �������� � ������ �������
										  //       ����� ����� ��������������� �����:
		cout << "\ncond1 A = " << norm(DIM1, ary, '1')*norm(DIM1, inverseA, '1') << endl; //�������������� �����
		cout << "\ncond2 A = " << norm(DIM1, ary, '2')*norm(DIM1, inverseA, '2') << endl; //��������� �����
		cout << "\ncond(inf) A = " << norm(DIM1, ary, 'k')*norm(DIM1, inverseA, 'k') << endl; //���������� �����
		cout << "\ncond(max) A = " << norm(DIM1, ary, 'm')*norm(DIM1, inverseA, 'm') << endl; //max-�� �����

		for (int i = 0; i < DIM1; i++) //������� ������������ ������ �������
		{
			delete[] inverseA[i];
		}
		delete[] inverseA;//������� ������������ �������
	}
	else { cout << "\ncond A = infinity" << endl; }


	for (int i = 0; i < DIM1; i++) //������� ������������ ������ �������
	{
		delete[] ary[i];
	}
	delete[] ary;//������� ������������ �������
	for (int i = 0; i < DIM1; i++) //������� ������������ ������ �������
	{
		delete[] T[i];
	}
	delete[] T;//������� ������������ �������

	for (int i = 0; i < DIM1; i++) //������� ������������ ������ �������
	{
		delete[] R[i];
	}
	delete[] R;//������� ������������ �������
	delete[] b; //������� ������������ �������

	for (int i = 0; i < DIM1; i++) //������� ������������ ������ �������
	{
		delete[] A[i];
	}
	delete[] A;//������� ������������ �������
	delete[] B; //������� ������������ �������

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

bool findx(const unsigned int DIM, mytipe* b, mytipe** a) //������ ��� ������ � ��������� �� �������������
{
	mytipe eps;     //��������������� ���������� ��� ��������� ��������� (���������� max)
	mytipe* vspom;  //��������������� ������ ��� ������������ ����� �������
	mytipe vspom2;   //��������������� ���������� ��� ������������ ��������� �������
	bool flag = true;  //���������� ����������(�������� ������������� �������)


	for (int k = 0, h; k < (DIM); k++)
	{
		eps = abs(a[k][k]);  //������ ������������� �������� � ���������� ���������
		vspom = a[k];   //������ ������ �� ������ ������������� ��������(���� ������� � ����� �� ���� �� ����������)
		h = k;   //������ ������ ������ (�� �������, ��������� ����)

		if (abs(a[k][k]) < EPSILON) { flag = true; } //���� ������ ������� �������, �� ������������ ����
		else { flag = false; } //���� �� �������, �� �������������� ����

		for (int i = (k + 1); i < DIM; i++)
		{
			if (abs(a[i][k]) > eps && abs(a[i][k])>EPSILON) //���� ������� ������ ����������� max � ������� � �� ����
			{
				eps = a[i][k]; vspom = a[i];  h = i; flag = false;
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
void neviaz(const unsigned int DIM, const mytipe* const x) //�������
{
	mytipe* b1;          //��������� ������ ���
	b1 = new mytipe[DIM];  //������, ������� ��������� ��� ����������� ������� � �������


	for (int i = 0; i < DIM; i++)
	{
		b1[i] = 0;
		for (int j = 0; j < DIM; j++)
		{
			b1[i] += A[i][j] * x[j]; //�������, ������������ �������, ������� ������ �����

		}
		b1[i] -= B[i]; //�������� i-�� ���������� ������������� ������� � ���������
	}
	cout << "\n||b1- b||1 = " << norm(DIM, b1, '1') << endl;  //����� �������������� ����� (�������)
	cout << "\n||b1- b||2 = " << norm(DIM, b1, '2') << endl;  //����� ���������� ����� (�������)
	cout << "\n||b1- b||(inf) = " << norm(DIM, b1, 'k') << endl;  //����� ���������� ����� (�������)

	delete[] b1; //������������ ������ ������������� �������
}

void copy(const unsigned int DIM, mytipe* b, mytipe** a) //����������� ������� � �������
{

	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			a[i][j] = A[i][j];
		}
		b[i] = B[i];
	}
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

mytipe norm(const unsigned int DIM, const mytipe* const b, char flag) //3 ����� ������� � �������� �������� �����
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

void obusl(const unsigned int DIM, const mytipe* const x) //������ ����� ����� ���������������
{
	neviaz(DIM, x);  //�������

	mytipe* b;             //|��������� ������ ��� ��������� ���� ������� 
	b = new mytipe[DIM];   //|

	mytipe **ary;                     // |������� 
	ary = new mytipe *[DIM];          // |������������ 
	for (int i = 0; i < DIM; i++) {   // |������ 
		ary[i] = new mytipe[DIM];     // |������� �
	}

	copy(DIM, b, ary);  //���������� �������� ������� � �������

	mytipe **T;                       // |������� ������������� ������� � �������� ������� = Q
	T = new mytipe *[DIM];           // |
	for (int i = 0; i < DIM; i++) {  // |
		T[i] = new mytipe[DIM];      // |
	}	                               //|

	singlematr(DIM, T);

	QR(DIM, T, ary);

	bool flag;  //���� �� ������������� �������
	mytipe* deltaX;            //|�������� ������ ���
	deltaX = new mytipe[DIM];  //|������ �������� ������� ������� � � �����������
	const mytipe VOZM = 0.01;  //����������
	mytipe dx;   //������������� ����������� x
	mytipe db1 = VOZM / norm(DIM, b, '1');  //������������� ����������� b �� �������������� �����
	mytipe db2 = VOZM / norm(DIM, b, '2');  //������������� ����������� b �� ���������� �����
	mytipe dbk = VOZM / norm(DIM, b, 'k');  //������������� ����������� b �� ���������� �����
	mytipe norm1 = norm(DIM, x, '1');  //�������������� ����� ������� �
	mytipe norm2 = norm(DIM, x, '2');  //��������� ����� ������� �
	mytipe normk = norm(DIM, x, 'k');  //���������� ����� ������� �

	mytipe max1 = 0; //������������ ������ ����� ��������������� ����� �� �������������� �����
	mytipe max2 = 0; //������������ ������ ����� ��������������� ����� �� ���������� �����
	mytipe maxk = 0; //������������ ������ ����� ��������������� ����� �� ���������� �����

					 //������ ����������� ������� ������� ������
	for (int i = 0; i < 2 * DIM; i++)
	{
		if (i < DIM) { b[i] += VOZM; cout << i + 1 << ") b*="; vivod(DIM, b); } //���������� ����������
		else { b[i - DIM] -= VOZM; cout << i - DIM + 1 << ") b*="; vivod(DIM, b); } //�������� ����������

		flag = findx(DIM, b, ary); //�������� � ������������������ ���� � ������������ �������������

		matrvec(DIM, T, b);

		if (flag)  //���� ������� �����������
		{
			hod(DIM, b, ary);  //�������� ��� ������
		}
		else { cout << "������� ��������� => ������ ����� ��������������� ����������"; cin.get(); break; }

		for (int j = 0; j < DIM; j++) { deltaX[j] = abs(x[j] - b[j]); } //������ ��������

		cout << " x* = "; vivod(DIM, b); //������� ������������ �����8�� � �����������

		copy(DIM, b, B);  //���������� ������� b

		dx = norm(DIM, deltaX, '1') / norm1; //������� ������������� ����������� ������� ������� �������

		cout << "\n cond1 A >= " << dx / db1; //������� ������ ���� ��������������� �� ����� 1

		if (max1 < dx / db1) { max1 = dx / db1; } //��������� �������� 

		dx = norm(DIM, deltaX, '2') / norm2;
		cout << "; cond2 A >= " << dx / db2; //������� ������ ���� ��������������� �� ����� 2
		if (max2 < dx / db2) { max2 = dx / db2; }

		dx = norm(DIM, deltaX, 'k') / normk;
		cout << "; cond(inf) A >= " << dx / dbk << endl << endl; //������� ������ ���� ��������������� �� ���������� �����
		if (maxk < dx / dbk) { maxk = dx / dbk; }
	}
	if (max1 != 0) {
		cout << "������������ ������ ����� ��������������� �����:\n cond1 A >= " << max1 << " cond2 A >= " << max2 << " cond(inf) A >= " << maxk << endl << endl;
	}
	delete[] deltaX;
	delete[] b;
	for (int i = 0; i < DIM; i++) //������� ������������ ������ �������
	{
		delete[] ary[i];
	}
	delete[] ary;//������� ������������ �������

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

void QR(const unsigned int DIM, mytipe** T, mytipe** R) //QR-��������
{
	mytipe vspom2;   //��������������� ��������� ��� ����� ���������
	mytipe c;  //����������� ����������� �_ij
	mytipe s;  //����������� ����������� �_ij

	for (int k = 0; k < DIM; k++)  //���������� ������ T=Q^(-1) � R
								   //������������ ����� T_ij, ��������� �� ������� ���������
	{
		for (int i = k + 1; i < DIM; i++)  //������������ ����� T_ij, ������� ������ c ���������� �������
		{
			if (abs(R[i][k])>EPSILON)
			{
				c = R[k][k] / sqrt((R[k][k])*(R[k][k]) + (R[i][k])*(R[i][k]));//���������� c
				s = R[i][k] / sqrt((R[k][k])*(R[k][k]) + (R[i][k])*(R[i][k]));//���������� s
				for (int j = 0; j < DIM; j++)  //��������� T_ij �� A
				{
					vspom2 = T[k][j];
					T[k][j] = c * T[k][j] + s * T[i][j]; //��������� ���� ����� �
					T[i][j] = c * T[i][j] - s * vspom2;  //�� ������� ���������� �
				}
				for (int j = k; j < DIM; j++)
				{
					vspom2 = R[k][j];
					R[k][j] = c * R[k][j] + s * R[i][j]; //��������� ���� ����� �
					R[i][j] = c * R[i][j] - s * vspom2;    //�� ������� �
				}
			}

		}
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