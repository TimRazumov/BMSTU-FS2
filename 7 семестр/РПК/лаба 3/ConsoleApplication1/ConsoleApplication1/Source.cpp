//d/dx(k(x)du/dx)=f(x)

#include<iostream>
#include <fstream>

using namespace std;

const double Eps = 1e-8;//программный ноль

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
	int N = 100;//кол-во разбиений 
	int n = N + 1;//кол-во узлов
	double h = (b - a) / double(N);
	double Ua = 10.; double Ub = 20.;//граничные условия
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
	
	int** Elems = new int* [2];//таблица связности
	for (int i = 0; i < 2; i++) Elems[i] = new int [N];

	for (int j = 0; j < N; j++) { Elems[0][j] = j; Elems[1][j] = j + 1; };

	double M = 10e+10;//параметр штрафа
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

	ofstream fout;    // создали переменную для записи в файл
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


bool findx(const unsigned int DIM, double* b, double** a) //прямой ход Гаусса с проверкой на вырожденность
{
	double max;     //вспомогательная переменная для сравнения элементов (нахождения max)
	double* vspom;  //вспомогательная ссылка для перестановки строк матрицы
	double vspom2;   //вспомогательная переменная для перестановки элементов вектора
	bool flag = true;  //логическая переменная(проверка вырожденности матрицы)


	for (int k = 0, h; k < (DIM); k++)
	{
		max = fabs(a[k][k]);  //запись диагонального элемента в переменную сравнения
		vspom = a[k];   //запись ссылки на строку диагонального элемента(если условие в цикле ни разу не выполнится)
		h = k;   //запись номера строки (по причине, описанной выше)

		if (fabs(a[k][k]) < Eps) { flag = true; } //если первый элемент нулевой, то активировать флаг
		else { flag = false; } //если не нулевой, то деактивировать флаг

		for (int i = (k + 1); i < DIM; i++)
		{
			if (fabs(a[i][k]) > max && fabs(a[i][k])>Eps) //если элемент больше предыдущего max в столбце и не ноль
			{
				max = fabs(a[i][k]); vspom = a[i];  h = i; flag = false;
			}; //запомнить строку, изменить max и деактивировать флаг
		}

		if (!flag) //если столбец не нулевой
		{
			a[h] = a[k];  //обмен строками
			a[k] = vspom;
			vspom2 = b[h]; //обмен элементами столбца
			b[h] = b[k];
			b[k] = vspom2;
		}
		else { cout << "Матрица вырождена\n";   return false; }

		for (int i = k + 1; i < DIM; i++)
		{
			b[i] = b[i] - b[k] / a[k][k] * a[i][k]; //вычитание из элементов столбца
			for (int j = (DIM - 1); j >= 0; j--)
			{
				a[i][j] = a[i][j] - a[k][j] / a[k][k] * a[i][k]; //зануление всех элементов под диагональным
			}
		}

	}
	return true;
}

void hod(const unsigned int DIM, double* b, double** a) //обратный ход Гаусса
{

	for (int i = (DIM - 1); i >= 0; i--)
	{
		for (int j = i + 1; j < DIM; j++)
		{
			b[i] = b[i] - a[i][j] * b[j]; //вычитание из столбца правой части всех элементов матрицы этой же строки кроме элемента на глоавной диагонали
		}
		b[i] = b[i] / a[i][i]; //деление элемента столбца правой части на элемента главной диагонали
	}
}
