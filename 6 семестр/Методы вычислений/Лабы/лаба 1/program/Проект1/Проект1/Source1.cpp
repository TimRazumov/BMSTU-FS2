#include "iostream"
#include <fstream>
#include <typeinfo>
#include <cmath>
#include <vector>

using namespace std;

const double EPSILON = 1e-14;   //��������� ��������� � 0

const double PI = 3.14159265;

typedef double mytipe;  //��� ������, �������������� �� ���� ���������

const mytipe Eps = 1e-6; //�������� ���������� ������

//typedef mytipe(*func)(const mytipe x);
//typedef mytipe*(*func_sys)(const mytipe* z0);




mytipe** Yakobi_Matr(int dim, char,mytipe* x, mytipe* v, mytipe h,mytipe tn); //���������� ������� �����
mytipe** Inverse(mytipe dim, mytipe** A);
void QR(const unsigned int DIM, mytipe** T, mytipe** R, mytipe** A); //QR-��������
void singlematr(const unsigned int DIM, mytipe** const E); //���������� ������� ���������
void matrvec(const unsigned int DIM, mytipe** A, mytipe* b); //��������� ������� �� �������
void trans(const unsigned int DIM, mytipe** a); //���������������� �������
void hod(const unsigned int DIM, mytipe* b, mytipe** a); //�������� ��� ������

mytipe error(int n, mytipe** yn, mytipe h);


void Newton(char, int dim,mytipe* x0, mytipe* v, mytipe* result, mytipe h,mytipe tn); //��� ������

mytipe** Euler_implicit(mytipe* u0, mytipe h, int n, int dim);
mytipe** simmetr(mytipe* u0, mytipe h, int n, int dim);
mytipe** Euler_explicit(mytipe* u0, mytipe h, int n, int dim); 
mytipe** Adams(mytipe* u0, mytipe h, int n, int dim);
mytipe** prediction_correction(mytipe* u0, mytipe h, int n, int dim);
mytipe** RungeKutta_2(mytipe* u0, mytipe h, int n, int dim);
mytipe** RungeKutta_4(mytipe* u0, mytipe h, int n, int dim); //�����-����� 4-��� �������

void copy(const unsigned int DIM, mytipe* b, const mytipe* a); //����������� ��������
void copy(const unsigned int DIM, mytipe** B, mytipe** const A); //����������� ������
void vivod(const unsigned int DIM, mytipe** const a); //��������� ������ �������
void vivod(const unsigned int DIM, const mytipe* const b); //��������� ������ �������
mytipe norm(const unsigned int DIM, const mytipe* const b, const char flag); //3 ����� ������� � �������� �������� �����
void mult_matrvect(const unsigned int DIM, mytipe** const A, const mytipe* const b, mytipe* x); //��������� ������� �� ������
void minus_vect(const unsigned int DIM, const mytipe* const b1, const mytipe* const b2, mytipe* x); //�������� ��������
void sum_vect(const unsigned int DIM, const mytipe* const b1, const mytipe* const b2, mytipe* x); //�������� ��������
mytipe* mult_vectnum(const unsigned int DIM, const mytipe* const b, const mytipe a); //��������� ������� �� �����

void copy_file(int dim, int n, mytipe** Y);


mytipe analit(mytipe tn)
{
	mytipe w = 1.;//sqrt(66.7);
	return tn*tn - tn + cos(w*tn) + sin(w*tn);//cos(w*tn);// �������
}

mytipe* f(const mytipe* x,const mytipe tn)
{

	//�������
	mytipe c = 1.;// 66.7;
	mytipe* z = new mytipe[2];
	z[0] = x[1]; 
	z[1] = -c*x[0] + tn*tn - tn + 2.;

	//���� 1
	/*mytipe* z = new mytipe[2];
	z[0] = 2 * x[0] + x[1] * x[1] - 1;
	z[1] = 6 * x[0] - x[1] * x[1] + 1;*/

	//���� 2
	/*mytipe* z = new mytipe[2];
	z[0] = 1- x[0]*x[0] - x[1] * x[1];
	z[1] = 2 * x[0];*/

	//������� 13 (��������� �������� � ������ ������)
	/*mytipe* z = new mytipe[2];
	mytipe delta=0.001;
	mytipe a=2.;
	mytipe b=0.3;
	mytipe F=0.05;
	mytipe v=sqrt(a)-0.2;
	z[0] = x[1];
	z[1] = -2 * delta*x[1] - a*x[0] - b*x[0] * x[0] * x[0] + F*cos(v*tn);*/

	//������� 20 (������-������)
	/*mytipe* z = new mytipe[2];
	mytipe r1 = 0.4;
	mytipe r2 = 0.1;
	mytipe b11 = 0.05;
	mytipe b12 = 0.1;
	mytipe b21 = 0.08; 
	mytipe b22 = 0.003;
	z[0] = r1*x[0] - b11*x[0] * x[0] - b12*x[0] * x[1];
	z[1] = -r2*x[1] - b22*x[1] * x[1] + b21*x[0] * x[1];*/


	//���� 3
    /*mytipe* z = new mytipe[3];
	mytipe sigma = 10;
	mytipe r = 28;
	mytipe b = 8 / 3;
	z[0] = sigma*(x[1] - x[0]); 
	z[1] = x[0] * (r - x[2]) - x[1];
	z[2] = x[0] * x[1] - b*x[2];*/

	return z;
}


mytipe* g(int dim,char F,const mytipe* s0, const mytipe* y0, mytipe h, const mytipe tn) //s0 --��������� ��������� ����������� ��� ����� ����, 
{                                                             //y0 -- ��� ���������� ������ u(ti)
	mytipe* z = new mytipe[dim];

	if(F=='e')//����� 
		for (int i = 0; i < dim; i++) z[i] = s0[i] - y0[i] - h  * f(s0, tn)[i];

	if (F == 's') //������������
		for (int i = 0; i < dim; i++) z[i] = s0[i] - y0[i] - h / 2 * (f(s0, tn)[i] + f(y0, tn)[i]);


	return z;
}

int main() {
	setlocale(LC_ALL, "Russian"); //���������� ������� ����
	mytipe h=0.2;
	mytipe t0=0.0;
	mytipe T=1.;
	int n = (T - t0) / h;
	int dim = 2;
	mytipe* u0 = new mytipe[dim]; 

	 u0[0] = 1; u0[1] = 0; //T = 1; �������
	
	// u0[0] = 0.0; u0[1] = -0.99996;//(h=0.0001; T=4;-0.99999;//-1.00001;) //���� 1 (2)
	//u0[0] = 7.0000; u0[1] = -1.00001;//���� 1 (1)

	//T = 4.0; u0[0] = 0.00; u0[1] = 1.14; //���� 2 (1)

	 //u0[0] = 0.1; u0[1] = 0.1; //T = 50.0; //������� 13

	 //u0[0] = 1.; u0[1] = 4.; //T = 50.0; //������� 20

	//T = 30.0; u0[0] = 2; u0[1] = 2; u0[2] = 2; //���� 3


	mytipe** Y = new mytipe*[n];
	for (int i = 0; i< n; i++) Y[i] = new mytipe[dim];

	//Y=Euler_implicit(u0, h, n, dim);
	//Y=Euler_explicit(u0, h, n, dim);
	//Y=simmetr(u0, h, n, dim);
	//Y = Adams(u0, h, n, dim);
	//Y=prediction_correction(u0, h, n, dim);
	//Y = RungeKutta_2(u0, h, n, dim);
	Y = RungeKutta_4(u0, h, n, dim);

	copy_file(dim,n,Y);


	cout <<" �����������: " <<error(n, Y, h);

	for (int i = 0; i < n; i++) delete[] Y[i];
	delete[] Y;

	cout << "\nEND:)";
	cin.get();
	return 0;
}

void copy_file(int dim, int n, mytipe** Y)
{
	ofstream fout;    // ������� ���������� ��� ������ � ����
	fout.open("iter.txt", ios_base::out | ios_base::trunc);
	for (unsigned int i = 0; i < n; i++) {
		for (unsigned int j = 0; j < dim; j++) {
			fout << Y[i][j] << " ";
		}
		fout << endl;
	}
	fout.close();
	fout.clear();
};

mytipe** Euler_explicit(mytipe* u0, mytipe h, int n, int dim)//�����
{
	mytipe** yn = new mytipe*[n];
	for (int i = 0; i < n; i++) yn[i] = new mytipe[dim];
	yn[0] = u0;

	for (int i = 0; i < n - 1; i++)
	{
		sum_vect(dim, yn[i], mult_vectnum(dim, f(yn[i], i*h), h), yn[i + 1]);
	}

	return yn;
}

mytipe** Euler_implicit(mytipe* u0, mytipe h, int n, int dim)//�������
{
	char flag = 'e';
	mytipe** yn = new mytipe* [n];
	for (int i = 0; i< n; i++) yn[i] = new mytipe[dim];
	
	yn[0] = u0;

	mytipe* vspom = new mytipe [dim];

	for (int i = 0; i < n-1; i++)
	{
		sum_vect(dim, yn[i], mult_vectnum(dim, f(yn[i],i*h),h), vspom);//��������� �����������
		Newton(flag,dim,vspom, yn[i],yn[i+1],h,(i+1)*h);
	}
	
return yn;
}

mytipe** simmetr (mytipe* u0, mytipe h, int n, int dim)
{
	char flag = 's';
	mytipe** yn = new mytipe*[n];
	for (int i = 0; i< n; i++) yn[i] = new mytipe[dim];

	yn[0] = u0;

	mytipe* vspom = new mytipe[dim];

	for (int i = 0; i < n - 1; i++)
	{
		sum_vect(dim, yn[i], mult_vectnum(dim, f(yn[i],i*h), h), vspom);//��������� �����������
		Newton(flag, dim, vspom, yn[i], yn[i + 1], h, (i + 1)*h);
	}
	
	delete[] vspom;
	return yn;
}



mytipe** Adams(mytipe* u0, mytipe h, int n, int dim)
{
	mytipe** yn = new mytipe*[n];
	for (int i = 0; i< n; i++) yn[i] = new mytipe[dim];
	yn[0] = u0;

	mytipe** y123 = new mytipe*[4];
	for (int i = 0; i< n; i++) yn[i] = new mytipe[dim];
	y123 = RungeKutta_4(u0, h, 4, dim);
	for (int i = 0; i < 4; i++) { copy(dim, yn[i], y123[i]); };

	mytipe* vspom = new mytipe [dim];
	for (int i = 3; i < n-1; i++) {
		for (int j = 0; j < dim; j++)
		sum_vect(dim, mult_vectnum(dim, f(yn[i],i*h), 55.), mult_vectnum(dim, f(yn[i-1],(i-1)*h), -59.),vspom);
		sum_vect(dim, vspom, mult_vectnum(dim, f(yn[i - 2],(i-2)*h), 37.), vspom);
		sum_vect(dim, vspom, mult_vectnum(dim, f(yn[i - 3], (i - 3)*h), -9.), vspom);
		sum_vect(dim, yn[i], mult_vectnum(dim, vspom, h/24.), yn[i + 1]); };

	delete[] vspom; 
	for (int i = 0; i < 4; i++) { delete[] y123[i]; }; delete[] y123;
	return yn;
}


mytipe** prediction_correction(mytipe* u0, mytipe h, int n, int dim)
{
	mytipe** yn = new mytipe*[n];
	for (int i = 0; i< n; i++) yn[i] = new mytipe[dim];
	yn[0] = u0;

	mytipe** y123 = new mytipe*[4];
	for (int i = 0; i< n; i++) yn[i] = new mytipe[dim];
	y123=RungeKutta_4(u0, h, 4, dim);
	for (int i = 0; i < 4; i++) { copy(dim,yn[i],y123[i]); };

	mytipe* vspom = new mytipe[dim];
	for (int i = 3; i < n - 1; i++) {
		for (int j = 0; j < dim; j++)
		sum_vect(dim, mult_vectnum(dim, f(yn[i], i*h), 55.), mult_vectnum(dim, f(yn[i - 1], (i - 1)*h), -59.), vspom);
		sum_vect(dim, vspom, mult_vectnum(dim, f(yn[i - 2], (i - 2)*h), 37.), vspom);
		sum_vect(dim, vspom, mult_vectnum(dim, f(yn[i - 3], (i - 3)*h), -9.), vspom);
		sum_vect(dim, yn[i], mult_vectnum(dim, vspom, h / 24.), yn[i + 1]);

		sum_vect(dim, mult_vectnum(dim, f(yn[i+1],(i +1)*h), 9.), mult_vectnum(dim, f(yn[i],i *h), 19.), vspom);
		sum_vect(dim, vspom, mult_vectnum(dim, f(yn[i - 1], (i - 1)*h), -5.), vspom);
		sum_vect(dim, vspom,f(yn[i - 2],(i - 2)*h), vspom);
		sum_vect(dim, yn[i], mult_vectnum(dim, vspom, h / 24.), yn[i + 1]);
	};

	delete[] vspom;
	for (int i = 0; i < 4; i++) { delete[] y123[i]; }; delete[] y123;
	return yn;
};

mytipe** RungeKutta_2(mytipe* u0, mytipe h, int n, int dim)
{
	mytipe* k1 = new mytipe[dim];
	mytipe* k2 = new mytipe[dim];
	mytipe* vspom = new mytipe[dim];
	mytipe** yn = new mytipe*[n];
	for (int i = 0; i< n; i++) yn[i] = new mytipe[dim];
	yn[0] = u0;
	for (int i = 0; i < n - 1; i++)
	{
		k1 = f(yn[i], i*h);
		sum_vect(dim, yn[i], mult_vectnum(dim, k1, h), vspom);
		k2 = f(vspom, (i + 1)*h);

		sum_vect(dim, k1, k2, vspom);
		sum_vect(dim, yn[i], mult_vectnum(dim, vspom, h / 2.), yn[i + 1]);
	}

	delete[] k1;
	delete[] k2;
	delete[] vspom;
	return yn;
}

mytipe** RungeKutta_4(mytipe* u0, mytipe h, int n, int dim)
{
	mytipe* k1 = new mytipe[dim];
	mytipe* k2 = new mytipe[dim];
	mytipe* k3 = new mytipe[dim];
	mytipe* k4 = new mytipe[dim];
	mytipe* vspom = new mytipe[dim];
	mytipe** yn = new mytipe*[n];
	for (int i = 0; i< n; i++) yn[i] = new mytipe[dim];
	yn[0] = u0;
	for (int i = 0; i < n - 1; i++)
	{
		k1 = f(yn[i],i*h);
		sum_vect(dim, yn[i], mult_vectnum(dim, k1, 0.5*h), vspom);
		k2 = f(vspom, (i + 0.5)*h);
		sum_vect(dim, yn[i], mult_vectnum(dim, k2, 0.5*h), vspom);
		k3 = f(vspom, (i + 0.5)*h);
		sum_vect(dim, yn[i], mult_vectnum(dim, k3, h), vspom);
		k4 = f(vspom, (i + 1.)*h);

		sum_vect(dim, k1, mult_vectnum(dim, k2, 2.), vspom);
		sum_vect(dim, vspom, mult_vectnum(dim, k3, 2.), vspom);
		sum_vect(dim, vspom, k4, vspom);
		sum_vect(dim, yn[i], mult_vectnum(dim, vspom, h / 6.), yn[i + 1]);
	}
		
		delete[] k1;
		delete[] k2;
		delete[] k3;
		delete[] k4;
		delete[] vspom;
   return yn;
}


void Newton(char flag, int dim,mytipe* x0,mytipe* v ,mytipe* result, mytipe h, mytipe tn) //n - ���-�� �������; x0 - ��� �������
{
	mytipe* xk = new mytipe[dim];
	mytipe* xk1 = new mytipe[dim];//��� �������
	mytipe** Yakob = new mytipe*[dim]; 
	for (int i = 0; i < dim; i++) Yakob[i] = new mytipe[dim];

	mytipe vspom;
	mytipe yakobian;
	int iter = 0;

	copy(dim, xk1, x0);

	do
	{
		iter++;
		copy(dim, xk, xk1); 

		copy(dim, Yakob, Yakobi_Matr(dim,flag,xk, v, h,tn));

		if (dim == 2){
		yakobian = Yakob[0][0] * Yakob[1][1] - Yakob[1][0] * Yakob[0][1];
		vspom = Yakob[0][0];
		Yakob[0][0] = Yakob[1][1] / yakobian;
		Yakob[1][1] = vspom / yakobian;
		Yakob[0][1] /= -yakobian;
		Yakob[1][0] /= -yakobian;}
		else Yakob = Inverse(dim, Yakob);
		

		mult_matrvect(dim, Yakob, g(dim,flag, xk,v,h,tn), xk1);
		minus_vect(dim, xk, xk1, xk1);
		minus_vect(dim, xk, xk1, xk);
	} while (norm(dim, xk, '2') > Eps);
	
	copy(dim, result, xk1);

	delete[] xk; delete[] xk1; 
	for (int i = 0; i < dim; i++) delete[] Yakob[i];
	delete[] Yakob;
}



mytipe error(int n, mytipe** yn,mytipe h)
{
	mytipe vspom; 
	mytipe error = fabs(yn[0][0] - analit(0));
	for (int i = 1; i < n; i++) 
	{ 
		vspom = fabs(yn[i][0] - analit(h*i));
		if (vspom > error) error = vspom;
	};
	return error;
}


mytipe** Yakobi_Matr(int dim,char flag,mytipe* x, mytipe* v, mytipe h, mytipe tn)//���������� ������� �����
{
	mytipe eps1 = 1e-5;

	mytipe** Yakob = new mytipe*[dim];
	for (int i = 0; i < dim; i++) Yakob[i] = new mytipe[dim];

	mytipe* fx = g(dim, flag, x, v, h, tn);
	
	mytipe*  vspom= new mytipe[dim];
	mytipe*  vspom1 = new mytipe[dim];

	for (int i = 0; i < dim; i++) {
		copy(dim, vspom, x);
		vspom[i] += eps1;
		copy(dim,vspom1,g(dim, flag, vspom, v, h, tn));
		for (int j = 0; j < dim; j++) Yakob[j][i] = (vspom1[j] - fx[j]) / eps1; 
	}
	
	delete[] vspom; 
	delete[] vspom1;
	delete[] fx;

	return Yakob;
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

void mult_matrvect(const unsigned int DIM, mytipe** const A, const mytipe* const b, mytipe* x) //��������� ������� �� ������
{
	for (int i = 0; i<DIM; i++) {
		x[i] = 0;
		for (int j = 0; j < DIM; j++) {
			x[i] += A[i][j] * b[j];
		}
	}
}

void minus_vect(const unsigned int DIM, const mytipe* const b1, const mytipe* const b2, mytipe* x) //�������� ��������
{
	for (int i = 0; i < DIM; i++) {
		x[i] = b1[i] - b2[i];
	}
}

void sum_vect(const unsigned int DIM, const mytipe* const b1, const mytipe* const b2, mytipe* x) //�������� ��������
{
	for (int i = 0; i < DIM; i++) {
		x[i] = b1[i] + b2[i];
	}
}


mytipe* mult_vectnum(const unsigned int DIM, const mytipe* const b, const mytipe a) //��������� ������� �� �����
{
	mytipe* x = new mytipe[DIM];
	for (int i = 0; i<DIM; i++) {
		x[i] = a * b[i];
	}
	return x;
}

void copy(const unsigned int DIM, mytipe* b, const mytipe* a) //����������� ��������
{
	for (int j = 0; j < DIM; j++)
	{
		b[j] = a[j];
	}
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

mytipe** Inverse(mytipe dim, mytipe** A)
{
	mytipe** inverseA = new mytipe *[dim]; mytipe** T = new mytipe*[dim]; mytipe** R = new mytipe*[dim];
		for (int i = 0; i < dim; i++) { T[i] = new mytipe[dim]; R[i] = new mytipe[dim]; inverseA[i] = new mytipe[dim];}
		singlematr(dim,inverseA);
		QR(dim, T, R, A);
		trans(dim,T);

		for (int r = 0; r < dim; r++)
		{
			matrvec(dim, T, inverseA[r]);
			hod(dim, inverseA[r], R); //������� ������ ������ �������� �������
		}
		trans(dim, inverseA); //�������� ����������������� �������� �������
		for (int i = 0; i < dim; i++) { delete[] T[i]; delete[] R[i];}
		delete[] R; delete[] T;
		return inverseA;
}

void QR(const unsigned int DIM, mytipe** T, mytipe** R, mytipe** A) //QR-��������
{
	mytipe vspom;   //��������������� ��������� ��� ����� ���������
	mytipe c;  //����������� ����������� �_ij
	mytipe s;  //����������� ����������� �_ij
	copy(DIM, R, A);
	singlematr(DIM, T);
	for (int k = 0; k < DIM; k++)  //���������� ������ T=Q^(-1) � R
								   //������������ ����� T_ij, ��������� �� ������� ���������
	{
		for (int i = k + 1; i < DIM; i++)  //������������ ����� T_ij, ������� ������ c ���������� �������
		{
			if (abs(R[i][k]) > EPSILON)
			{
				c = R[k][k] / sqrt((R[k][k])*(R[k][k]) + (R[i][k])*(R[i][k]));//���������� c
				s = R[i][k] / sqrt((R[k][k])*(R[k][k]) + (R[i][k])*(R[i][k]));//���������� s
				for (int j = 0; j < DIM; j++)  //��������� T_ij �� A
				{
					vspom = T[k][j];
					T[k][j] = c * T[k][j] + s * T[i][j]; //��������� ���� ����� �
					T[i][j] = c * T[i][j] - s * vspom;  //�� ������� ���������� �
				}
				for (int j = k; j < DIM; j++)
				{
					vspom = R[k][j];
					R[k][j] = c * R[k][j] + s * R[i][j]; //��������� ���� ����� �
					R[i][j] = c * R[i][j] - s * vspom;    //�� ������� �
				}
			}

		}
	}
	trans(DIM, T);//������� Q
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