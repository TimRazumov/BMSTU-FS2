#include "iostream"
#include <fstream>
#include <typeinfo>
#include <cmath>
#include "Matr.h"

using namespace std;

const double EPSILON = 1e-8;   //константа сравнения с 0
const double PI = 3.14159265;

typedef double mytipe;  //тип данных, использующийся во всей программе

typedef mytipe(*func)(const mytipe x);


void Uniform_grid(const unsigned int DIM, mytipe* x, const mytipe a, const mytipe b);
void Chebyshev_grid(const unsigned int DIM, mytipe* x, const mytipe a, const mytipe b);
void Values_in_nodes(const unsigned int DIM, const mytipe* const x, mytipe* fx, func F);
void Copy_to_file(const unsigned int n, const mytipe* const x, const mytipe* const Ln, const mytipe* const fx);

//форма Ньютона
mytipe* Polynom_Newton(const unsigned int DIM, const unsigned int DIM2, const mytipe* const xi, const mytipe* const fx, const mytipe* const x);
//форма Лагранжа
mytipe* Polynom_Lagrange(const unsigned int Nyzel, const unsigned int Nvalue, const mytipe* const xi, const mytipe* const yi, const mytipe* const x);
mytipe* Spline(int n, int N, mytipe* xi, mytipe* yi, mytipe* x);

mytipe* progon3d(int DIM, mytipe* a, mytipe* b, mytipe* c, mytipe* d);

mytipe error(int DIM, mytipe* a, mytipe* b);


mytipe Func(const mytipe x)
{
	return x*x;
	//return 1 / (1 + x*x);
	//return 1 / atan(1 + x*x);
	//return pow(4 * x*x*x + 2 * x*x - 4 * x + 2, sqrt(2)) + asin(1 / (5 + x - x*x)) - 5;
	//return x*x*x;
     //return 1.0;
	//return 1 / (1 + 25*x*x);
	//return pow(3, -x);
	//return exp(x);
}

int main() {
	setlocale(LC_ALL, "Russian"); //подключаем русский язык

	unsigned int n=5;//кол-во разбиений

	mytipe a = -1;//отрезок
	mytipe b = 1;

	const  unsigned int Nyzel = n + 1;//кол-во узлов

	mytipe* xi = new mytipe[Nyzel];         //узлы

    //Chebyshev_grid(n, xi, a, b);
    Uniform_grid(n, xi, a, b);

	mytipe* fxi = new mytipe[Nyzel];          //значения целевой функции в узлах

	Values_in_nodes(Nyzel, xi, fxi, Func);

	unsigned int N = 900;//разбиение отрезка для вычисления значений интерполянта 

	const  unsigned int Nvalue = N + 1;//кол-во вычисляемых знач интерп многочлена

	mytipe* x = new mytipe[Nvalue];          //узлы интерполянта

	Uniform_grid(N, x, xi[0], xi[n]);

	mytipe* pn = new mytipe[Nvalue];          //значения интерполянта 

    //pn = Spline(n, N, xi, fxi, x);
	pn=Polynom_Lagrange(Nyzel,Nvalue, xi, fxi, x);//форма Лагранжа
	//pn=Polynom_Newton(Nyzel, Nvalue, xi, fxi, x);//форма Ньютона

	 mytipe* f_true = new mytipe[Nvalue];
	 Values_in_nodes(Nvalue, x, f_true, Func);

	Copy_to_file(Nvalue, x, pn,f_true);//Копируем в файл для построения графиков в матлабе

	cout << "Интерполирование функции проведено успешно\n\n";

	cout << "Ошибка интерполяции: " << error(Nvalue, pn, f_true);

	delete[] x;
	delete[] fxi;
	delete[] xi;
	delete[] pn;
	delete[] f_true;

	cin.get();
	return 0;
}

mytipe* Polynom_Newton(const unsigned int DIM, const unsigned int DIM2, const mytipe* const xi, const mytipe* const fx, const mytipe* const x)
{
	mytipe* f_vspom = new mytipe[DIM];       //выделение памяти под
	mytipe* f_vspom2 = new mytipe[DIM];          //выделение памяти под
	copy(DIM, f_vspom2, fx);
	copy(DIM, f_vspom, fx);
	mytipe vspom1;
	int vspom2 = 0;
	mytipe* pn = new mytipe[DIM2];

	for (int j = 1; j < DIM; j++) {
		for (int k = j; k < DIM; k++) {
			f_vspom2[k] = (f_vspom[k] - f_vspom[k - 1]) / (xi[k] - xi[k - j]);
		}
		copy(DIM, f_vspom, f_vspom2);
	}

	for (int i = 0; i < DIM2; i++) {
		pn[i] = fx[0];
		vspom1 = x[i] - xi[0];
		pn[i] += vspom1*f_vspom[1];
		for (int j = 2; j < DIM; j++) {
			vspom1 *= x[i] - xi[j - 1];
			pn[i] += vspom1*f_vspom[j];
		}
	}

	delete[] f_vspom;
	delete[] f_vspom2;
	return pn;
}


mytipe* Polynom_Lagrange(const unsigned int Nyzel,const unsigned int Nvalue, const mytipe* const xi, const mytipe* const yi, const mytipe* const x)
{
	mytipe* Ln = new mytipe[Nvalue];
	mytipe* c= new mytipe[Nyzel];
	mytipe* Ci = new mytipe[Nyzel];

	for (int k = 0; k < Nyzel; k++) {
		c[k] = 1.;
		for (int j = 0; j < k; j++)  c[k] *= (xi[k] - xi[j]);
		for (int j = k+1; j < Nyzel; j++) c[k] *= (xi[k] - xi[j]);	
	}
	
	for (int i = 0; i < Nvalue; i++)
	{
		for (int k = 0; k < Nyzel; k++) {// вычисляем все Ck в точке i
			Ci[k] = 1.;
			for (int j = 0; j < k; j++)  Ci[k] *= (x[i] - xi[j]);
			for (int j = k + 1; j < Nyzel; j++) Ci[k] *= (x[i] - xi[j]);
			Ci[k] /= c[k];
		}

		Ln[i] = 0.;
		for (int k = 0; k < Nyzel; k++) Ln[i] += Ci[k] * yi[k];
	}

	delete[] c;
	delete[] Ci;
	return Ln;
}

void Copy_to_file(const unsigned int n, const mytipe* const x, const mytipe* const Ln, const mytipe* const fx)
{
	ofstream fout;    // создали переменную для записи в файл
	fout.open("polynom.txt", ios_base::out | ios_base::trunc);
	for (unsigned int i = 0; i < n; i++) {
		fout << x[i] <<" "<<Ln[i]<< " " << fx[i] << endl;
	}
	fout.close();
	fout.clear();
}

void Uniform_grid(const unsigned int DIM, mytipe* x, const mytipe a, const mytipe b)
{
	mytipe h = abs(b - a) / DIM;
	for (int i = 0; i <= DIM; i++) {
		x[i] = a + h*i;
	}
}

void Chebyshev_grid(const unsigned int DIM, mytipe* x, const mytipe a, const mytipe b)
{
	mytipe h = abs(b - a) / DIM;
	for (int i = 0, j = DIM; i <= DIM; i++, j--) {
		x[j] = (a + b + (b - a)*cos((2 * i + 1)*PI / 2 / (DIM+1))) / 2;
	}
}

void Values_in_nodes(const unsigned int DIM, const mytipe* const x, mytipe* fx, func F)
{
	for (int i = 0; i < DIM; i++) {
		fx[i] = F(x[i]);
	}
}

mytipe* Spline(int n, int N, mytipe* xi, mytipe* yi, mytipe* x)
{
	mytipe* a = new mytipe[n];
	mytipe* b = new mytipe[n];
	mytipe* c = new mytipe[n];
	mytipe* d = new mytipe[n];

	for (int i = 0; i < n; i++) { a[i] = yi[i]; }//определяем коэф a

	mytipe* h = new mytipe[n];
	mytipe* g = new mytipe[n];

	mytipe* diag3 = new mytipe[n - 1];
	mytipe* diag2 = new mytipe[n - 1];
	mytipe* diag1 = new mytipe[n - 1];
	mytipe* right = new mytipe[n - 1];

	for (int i = 0; i < n; i++) { h[i] = xi[i + 1] - xi[i];  g[i] = (yi[i + 1] - yi[i]) / h[i]; }
	for (int i = 0; i < n - 1; i++) {
		diag2[i] = 2 * (h[i + 1] + h[i]);
		diag1[i] = h[i];
		diag3[i] = h[i + 1];
		right[i] = 3 * (g[i + 1] - g[i]);
	}

	diag1[0] = 0.; diag3[n - 2] = 0.;

	right = progon3d(n - 1, diag1, diag2, diag3, right);

	c[0] = 0.0;//


	for (int i = 1; i < n; i++) { c[i] = right[i - 1]; };

	for (int i = 0; i < n - 1; i++) {
		b[i] = g[i] - (c[i + 1] + 2 * c[i])* h[i] / 3;
		d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
	}

	b[n - 1] = g[n - 1] - 2 * c[n - 1] * h[n - 1] / 3;
	d[n - 1] = -c[n - 1] / (3 * h[n - 1]);


	mytipe* S = new mytipe[N + 1];

	for (int j = 0, i = 0; j < n; j++)
		for (; (i < N + 1) && (x[i] <= xi[j + 1]); i++)
			S[i] = a[j] + b[j] * (x[i] - xi[j]) + c[j] * pow((x[i] - xi[j]), 2) + d[j] * pow((x[i] - xi[j]), 3);

	delete[] a;
	delete[] b;
	delete[] c;
	delete[] d;

	delete[] h;
	delete[] g;

	delete[] diag1;
	delete[] diag2;
	delete[] diag3;
	delete[] right;

	return S;
}


mytipe* progon3d(int DIM, mytipe* a, mytipe* b, mytipe* c, mytipe* d)

{
	mytipe* alfa = new mytipe[DIM];
	mytipe* betta = new mytipe[DIM];

	alfa[0] = -c[0] / b[0]; betta[0] = d[0] / b[0];

	for (int i = 1; i < DIM; i++) {
		alfa[i] = -c[i] / (b[i] + a[i] * alfa[i - 1]);
		betta[i] = (-a[i] * betta[i - 1] + d[i]) / (a[i] * alfa[i - 1] + b[i]);
	}

	for (int i = DIM - 2; i > -1; i--) {
		betta[i] += alfa[i] * betta[i + 1];
	}

	delete[] alfa;

	return  betta;
}

mytipe error(int DIM, mytipe* a, mytipe* b)
{
	mytipe MAXerror=0.0;
	for (int i = 0; i < DIM; i++) 
		if (fabs(a[i] - b[i]) > MAXerror) MAXerror = fabs(a[i] - b[i]);
		
	return MAXerror;
}





	


	

