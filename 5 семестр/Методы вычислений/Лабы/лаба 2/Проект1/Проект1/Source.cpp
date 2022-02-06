#include "iostream"
#include <fstream>
#include <typeinfo>
#include <cmath>

using namespace std;

const double EPSILON = 1e-8;   //константа сравнения с 0

typedef double mytipe;  //тип данных, использующийся во всей программе

mytipe *x_real;


void vivod(const unsigned int DIM, const mytipe* const b, mytipe** const a); //вывод матрицы и столбца
void vivod(const unsigned int DIM, mytipe** const a);  //вывод матрицы
void vivod(const unsigned int DIM, const mytipe* const b);  //вывод вектора
mytipe neviaz(const unsigned int DIM, const mytipe* const x, mytipe** const A, const mytipe* const b, const char k); //невязка
void copy(const unsigned int DIM, mytipe** B, mytipe** const A); //копирование матриц
void copy(const unsigned int DIM, mytipe* b, const mytipe* a); //копирование векторов
void trans(const unsigned int DIM, mytipe** a);  // транспонирование матрицы
void multimatrix(const unsigned int DIM, mytipe** const a, mytipe** const b); // перемножение матриц без возврата результата
void mult_matr(const unsigned int DIM, mytipe** const A1, mytipe** const A2, mytipe** X); // перемножение матриц без возврата результата
void mult_matrvect(const unsigned int DIM, mytipe** const A, const mytipe* const b, mytipe* x); // перемножение матриц
void sum_vect(const unsigned int DIM, const mytipe* const b1, const mytipe* const b2, mytipe* x); //сложение векторов
void sum_matr(const unsigned int DIM, const mytipe** A1, const mytipe** A2, mytipe** X); //сложение матриц
void mult_matrnum(const unsigned int DIM, mytipe** const A, const mytipe a, mytipe** X); //умножене матрицы на число
void mult_vectnum(const unsigned int DIM, const mytipe* const b, const mytipe a, mytipe* x); //умножение вектора на число
mytipe norm(const unsigned int DIM, mytipe** const A, const char flag); // 1 из 4 норм матрицы
mytipe norm(const unsigned int DIM, const mytipe* const b, const char flag);  // 1 из 3 норм вектора
void singlematr(const unsigned int DIM, mytipe** const E);  //заполнение матрицы единичной
void matrvec(const unsigned int DIM, mytipe** A, mytipe* b);  //умножение матрицы на столбец

mytipe SimpleIter_C(const unsigned int DIM, mytipe** const A, const mytipe* const b, const mytipe TAU, const char); //нахождение матрицы С в методе простой итерации
mytipe Yakobi_C(const unsigned int DIM, mytipe** const A, const mytipe* const b, const char); //нахождение матрицы С в методе Якоби
mytipe Seidel_C(const unsigned int DIM, mytipe** const A, const mytipe* const b, const char); //нахождение матрицы С в методе Зейделя
mytipe Relax_C(const unsigned int DIM, mytipe** const A, const mytipe* const b, const mytipe W, const char); //нахождение матрицы С в методе релаксации

void SimpleIter(const unsigned int DIM, mytipe** const A, const mytipe* const b, const mytipe TAU, const char NORM, mytipe EPS, const mytipe normC, mytipe* xk); //метод простой итерации
void Yakobi(const unsigned int DIM, mytipe** const A, const mytipe* const b, const char NORM, mytipe EPS, const mytipe normC, mytipe* xk); //метод Якоби
void Seidel(const unsigned int DIM, mytipe** const A, const mytipe* const b, const char NORM, mytipe EPS, const mytipe normC, mytipe* xk); //метод Зейделя
void Relax(const unsigned int DIM, mytipe** const A, const mytipe* const b, const mytipe W, const char NORM, mytipe EPS, const mytipe normC, mytipe* xk); //метод релаксации

void C_Iter(const unsigned int DIM, mytipe** const A, const mytipe* const b, mytipe** const C, const mytipe* const y, const char NORM); //решение итерационных методов, используя С и y

void relax3d(const int var, const char NORM, const mytipe w, const mytipe EPS);

void C_l_d_u(mytipe** C, const int DIM, const char NORM);

int main() {
	setlocale(LC_ALL, "Russian"); //подключаем русский язык

	ifstream fin;    // создали переменную для считывания из файла
	fin.open("1.txt", ios_base::in | ios_base::app | ios_base::binary);  //открыли файл
	unsigned int a;  // размер пространства
	fin >> a;  //считали из файла размер пространства
	const unsigned int DIM1 = a;  //размер пространнства сделали константным


	mytipe **A;  //матрица, списанная с файла
	A = new mytipe *[DIM1];            // |выделили
	for (int i = 0; i < DIM1; i++) {   // |место под
		A[i] = new mytipe[DIM1];       // |матрицу А
	}

	mytipe *b;  //столбец, записанный из файла
	b = new mytipe[DIM1];     //|выделили место под столбец b


	for (int i = 0; i < DIM1; i++)  //считываем из файла матрицу и столбец
	{
		for (int j = 0; j < DIM1; j++)
		{
			fin >> A[i][j];  //матрицу
		}
		fin >> b[i];  //столбец
	}

	fin.close();  //закрываем файл
	fin.clear();  //отвязываем файл от файловой переменной

	cout << "Dim=" << a << endl;
	vivod(DIM1, b, A);    //выводим систему

	const char NORM = 'k';
	const mytipe EPS = 1e-3;

	mytipe *x;  //столбец начальных приближений
	x = new mytipe[DIM1];


	x_real = new mytipe[DIM1];
	x_real[0] = 5;
	x_real[1] = -7;
	x_real[2] = 12;
	x_real[3] = 4;

	//МЕТОД ПРОСТОЙ ИТЕРАЦИИ
	for (int i = 0; i < DIM1; i++) { x[i] = 0.0; }
	const mytipe TAU = 0.06;
	mytipe normC_SimplIter = SimpleIter_C(DIM1, A, b, TAU, NORM);
	SimpleIter(DIM1, A, b, TAU, NORM, EPS, normC_SimplIter, x);
	cout << "Вектор x = \n"; vivod(DIM1, x);
	cout << "норма невязки = " << neviaz(DIM1, x, A, b, NORM) << endl;



	// МЕТОД ЯКОБИ	
     for (int i = 0; i < DIM1; i++) { x[i] = 0.0; }
	mytipe normC_Yakobi = Yakobi_C(DIM1, A, b, NORM);
	Yakobi(DIM1, A, b, NORM, EPS, normC_Yakobi, x);
	cout << "Вектор x = \n"; vivod(DIM1, x);
	cout << "норма невязки = " << neviaz(DIM1, x, A, b, NORM) << endl;


	// МЕТОД ЗЕЙДЕЛЯ
	for (int i = 0; i < DIM1; i++) { x[i] = 0.0; }
	mytipe normC_Seidel = Seidel_C(DIM1, A, b, NORM);
	Seidel(DIM1, A, b, NORM, EPS, normC_Seidel, x);
	cout << "Вектор x = \n"; vivod(DIM1, x);
	cout << "норма невязки = " << neviaz(DIM1, x, A, b, NORM) << endl;


	//МЕТОД РЕЛАКСАЦИИ
	for (int i = 0; i < DIM1; i++) { x[i] = 0.0; }
	const mytipe W = 0.9;
	mytipe normC_Relax = Relax_C(DIM1, A, b, W, NORM);
	Relax(DIM1, A, b, W, NORM, EPS, normC_Relax, x);
	cout << "Вектор x = \n"; vivod(DIM1, x);
	cout << "норма невязки = " << neviaz(DIM1, x, A, b, NORM) << endl;


	//3-x диагональная матрица
	relax3d(14, NORM, 1.9, EPS);


	delete[] b; //удаляем динамический столбец

	for (int i = 0; i < DIM1; i++) //удаляем динамические строки матрицы
	{
		delete[] A[i];
	}
	delete[] A;//удаляем динамическую матрицу
	delete[]x; delete[]x_real;

	cin.get();
	return 0;
}

void vivod(const unsigned int DIM, const mytipe* const b, mytipe** const a) //процедура вывода матрица и столбца
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

void vivod(const unsigned int DIM, mytipe** const a) //процедура вывода матрицы
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

void vivod(const unsigned int DIM, const mytipe* const b) //процедура вывода столбца
{
	cout << '(';
	for (int i = 0; i < DIM; i++)
	{
		cout << b[i];
		if (i < DIM - 1) { cout << ','; } //красивая запись
	}
	cout << ")^т\n";
}




mytipe neviaz(const unsigned int DIM, const mytipe* const x, mytipe** const A, const mytipe* const b, const char k) //невязка
{
	mytipe* b1;          //выделение памяти под
	b1 = new mytipe[DIM];  //вектор, который получится при подстановке решения в систему
	mytipe nor;  //норма вектора разности

	for (int i = 0; i < DIM; i++)
	{
		b1[i] = 0;
		for (int j = 0; j < DIM; j++)
		{
			b1[i] += A[i][j] * x[j]; //подсчет, подстановкой решения, столбца правой части

		}
		b1[i] -= b[i]; //разность i-ой компоненты получившегося вектора и истенного
	}

	nor = norm(DIM, b1, k);
	delete[] b1; //освобождение памяти динамического вектора
	return nor;
}


void copy(const unsigned int DIM, mytipe** B, mytipe** const A) //копирование матриц
{

	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			B[i][j] = A[i][j];
		}
	}
}

void copy(const unsigned int DIM, mytipe* b, const mytipe* a) //копирование векторов
{
	for (int j = 0; j < DIM; j++)
	{
		b[j] = a[j];
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

void multimatrix(const unsigned int DIM, mytipe** const a, mytipe** const b) //перемножение матриц
{
	double t;
	mytipe** E;
	E = new mytipe *[DIM];              // выделяем место
	for (int i = 0; i < DIM; i++)       // под матрицу
	{                                   //  являющуюся результатом
		E[i] = new mytipe[DIM];         //  перемножения матриц a и b
	}

	for (int i = 0; i<DIM; i++) {
		for (int l = 0; l<DIM; l++) {
			t = 0;
			for (int j = 0; j<DIM; j++) {
				t += a[i][j] * b[j][l];  //перемножение строки матрицы a на столбец матрицы b
			}
			E[i][l] = t;   //запись результата перемножения в результирующую матрицу
		}
	}

	vivod(DIM, E);  //выводим результат

	for (int i = 0; i < DIM; i++) {   //удаляем динамические строки матрицы
		delete[] E[i];
	}
	delete[] E;             // удаляем динамическую матрицу
}

void mult_matr(const unsigned int DIM, mytipe** const A1, mytipe** const A2, mytipe** X) //перемножение матриц
{
	for (int i = 0; i<DIM; i++) {
		for (int l = 0; l<DIM; l++) {
			X[i][l] = 0;
			for (int j = 0; j<DIM; j++) {
				X[i][l] += A1[i][j] * A2[j][l];  //перемножение строки матрицы a на столбец матрицы b
			}
		}
	}
}

void mult_matrvect(const unsigned int DIM, mytipe** const A, const mytipe* const b, mytipe* x) //умножение матрицы на вектор
{
	for (int i = 0; i<DIM; i++) {
		x[i] = 0;
		for (int j = 0; j < DIM; j++) {
			x[i] += A[i][j] * b[j];
		}
	}
}

void mult_matrnum(const unsigned int DIM, mytipe** const A, const mytipe a, mytipe** X) //умножене матрицы на число
{
	for (int i = 0; i<DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			X[i][j] = a * A[i][j];
		}
	}
}

void mult_vectnum(const unsigned int DIM, const mytipe* const b, const mytipe a, mytipe* x) //умножение вектора на число
{
	for (int i = 0; i<DIM; i++) {
		x[i] = a * b[i];
	}
}

void sum_vect(const unsigned int DIM, const mytipe* const b1, const mytipe* const b2, mytipe* x) //сложение векторов
{
	for (int i = 0; i < DIM; i++) {
		x[i] = b1[i] + b2[i];
	}
}

void sum_matr(const unsigned int DIM, const mytipe** A1, const mytipe** A2, mytipe** X) //сложение матриц
{
	for (int i = 0; i<DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			X[i][j] = A1[i][j] + A2[i][j];
		}
	}
}

mytipe norm(const unsigned int DIM, mytipe** const A, const char flag) //4 нормы матрицы с запросом варианта нормы
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
				vspom += abs(A[j][i]); //суммируем элементы i-ой строки
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
				vspom += abs(A[j][i]);  //суммируем элементы j-ой строки
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

		for (int j = 0; j < DIM; j++)
		{
			for (int i = 0; i<DIM; i++)
			{
				if (abs(A[j][i]) > norm) { norm = abs(A[j][i]); } //находим максимальный элемент
			}
		}
		return DIM*norm;
	}
}

mytipe norm(const unsigned int DIM, const mytipe* const b, const char flag) //3 нормы вектора с запросом варианта нормы
{
	mytipe norm = 0;  //получаемая норма вектора
	if (flag == 'k')  //кубическая норма (max)
	{

		for (int i = 0; i < DIM; i++)
		{
			if (norm < abs(b[i])) { norm = abs(b[i]); }  //находим максимум
		}
		return norm;
	}

	if (flag == '1') //октаэдрическая норма (сумма элементов)
	{
		for (int j = 0; j < DIM; j++)
		{
			norm += abs(b[j]);  //производим сложение элементов
		}
		return norm;
	}

	if (flag == '2')  //шаровая норма  (Евклидова)
	{
		for (int i = 0; i<DIM; i++)
		{
			norm += b[i] * b[i];  //суммируем квадраты элементов
		}

		return sqrt(norm); //возвращаем корень квадратный из суммы квадратов
	}
}


void singlematr(const unsigned int DIM, mytipe** const E) //заполнение матрицы единичной
{
	for (int i = 0; i < DIM; i++)      //|делаем
	{									//|матрицу Т
		for (int j = 0; j < DIM; j++)  //|единичной
		{							    //|для
			E[i][j] = 0;				//|дальнейших
		}								//|преобразований
		E[i][i] = 1;					//|
	}
}


void matrvec(const unsigned int DIM, mytipe** A, mytipe* b) //умножение матрицы на столбец
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

mytipe SimpleIter_C(const unsigned int DIM, mytipe** const A, const mytipe* const b, const mytipe tau, const char NORM) //нахождение матрицы С в методе простой итерации
{
	mytipe **C;  //матрица на каждом шаге итерации
	C = new mytipe *[DIM];
	for (int i = 0; i < DIM; i++) {
		C[i] = new mytipe[DIM];
	}

	mytipe *y = new mytipe[DIM];           //вектор неизвестных на (к+1)-ом шаге

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

	cout << "\nматрица C=";
	vivod(DIM, C);  //выводим матрицу C
	cout << "вектор y=\n";
	vivod(DIM, y);
	mytipe norma = norm(DIM, C, NORM);
	cout << "\nнорма матрицы C=" << norma << endl;

	//C_Iter(DIM, A, b, C, y, NORM);

	for (int i = 0; i < DIM; i++) //удаляем динамические строки матрицы
	{
		delete[] C[i];
	}

	delete[] C;//удаляем динамическую матрицу
	delete[] y;
	return norma;
}

mytipe Yakobi_C(const unsigned int DIM, mytipe** const A, const mytipe* const b, const char NORM) //нахождение матрицы С в методе Якоби
{
	mytipe **C;  //матрица на каждом шаге итерации
	C = new mytipe *[DIM];
	for (int i = 0; i < DIM; i++) {
		C[i] = new mytipe[DIM];
	}

	mytipe *y = new mytipe[DIM];           //вектор неизвестных на (к+1)-ом шаге


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

	cout << "\nматрица C=";
	vivod(DIM, C);  //выводим матрицу C
	cout << "вектор y=\n";
	vivod(DIM, y);
	mytipe norma = norm(DIM, C, NORM);
	cout << "\nнорма матрицы C=" << norma << endl;

	//C_Iter(DIM, A, b, C, y, NORM);

	for (int i = 0; i < DIM; i++) //удаляем динамические строки матрицы
	{
		delete[] C[i];
	}

	delete[] C;//удаляем динамическую матрицу
	delete[] y;
	return norma;
}

mytipe Seidel_C(const unsigned int DIM, mytipe** const A, const mytipe* const b, const char NORM) //нахождение матрицы С в методе Зейделя
{
	mytipe** vspom = new mytipe *[DIM];;                  //выделяем память под обратную матрицу (A1+D)
	mytipe **C = new mytipe *[DIM];  //матрица на каждом шаге итерации
	for (int i = 0; i < DIM; i++) { vspom[i] = new mytipe[DIM]; C[i] = new mytipe[DIM]; }

	mytipe *y = new mytipe[DIM];           //вектор неизвестных на (к+1)-ом шаге

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
	cout << "\nматрица C=";
	vivod(DIM, C);  //выводим матрицу C
	cout << "вектор y=\n";
	vivod(DIM, y);
	mytipe norma = norm(DIM, C, NORM);
	cout << "\nнорма матрицы C=" << norma << endl;

	//C_l_d_u(C, DIM, NORM);
	//C_Iter(DIM, A, b, C, y, NORM);

	for (int i = 0; i < DIM; i++) //удаляем динамические строки матрицы
	{
		delete[] C[i]; delete[] vspom[i];
	}

	delete[] C;//удаляем динамическую матрицу
	delete[] y;
	return norma;
}

mytipe Relax_C(const unsigned int DIM, mytipe** const A, const mytipe* const b, const mytipe W, const char NORM) //нахождение матрицы С в методе релаксации
{
	mytipe** vspom = new mytipe *[DIM];;                  //выделяем память под обратную матрицу (A1+D)
	mytipe **C = new mytipe *[DIM];  //матрица на каждом шаге итерации
	for (int i = 0; i < DIM; i++) { vspom[i] = new mytipe[DIM]; C[i] = new mytipe[DIM]; }

	mytipe *y = new mytipe[DIM];           //вектор неизвестных на (к+1)-ом шаге
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

	cout << "\nматрица C=";
	vivod(DIM, C);  //выводим матрицу C
	cout << "вектор y=\n";
	vivod(DIM, y);
	mytipe norma = norm(DIM, C, NORM);
	cout << "\nнорма матрицы C=" << norma << endl;

	//C_l_d_u(C,DIM,NORM);
	//C_Iter(DIM, A, b, C, y, NORM);

	for (int i = 0; i < DIM; i++) //удаляем динамические строки матрицы
	{
		delete[] C[i]; delete[] vspom[i];
	}

	delete[] C;//удаляем динамическую матрицу
	delete[] y;
	return norma;
}

void SimpleIter(const unsigned int DIM, mytipe** const A, const mytipe* const b, const mytipe TAU, const char NORM, mytipe EPS, mytipe normC, mytipe* xk) //метод простой итерации
{
	cout << "\n\t\t\t----------------------\n\t\t\tМЕТОД ПРОСТОЙ ИТЕРАЦИИ\n\t\t\t----------------------\n ";
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
		cout << "норма погрешности на " << iter << "-ой итерации:  " << norm(DIM, V, NORM) << endl;
		for (int i = 0; i < DIM; i++) { V[i] = xk1[i] - xk[i]; }
		norm_new = norm(DIM, V, NORM);
		copy(DIM, xk, xk1);
	} while (norm_new > EPS);
	
	
	cout << "Метод простой итерации выполнен \n количество итераций = " << iter << endl; 
	delete[] xk1;
	delete[] V;
}

void Yakobi(const unsigned int DIM, mytipe** const A, const mytipe* const b, const char NORM, mytipe EPS, const mytipe normC, mytipe* xk) //метод Якоби
{
	cout << "\n\t\t\t-----------\n\t\t\tМЕТОД ЯКОБИ\n\t\t\t-----------\n ";
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
		cout << "норма погрешности на " << iter << "-ой итерации:  " << norm(DIM, V, NORM) << endl;

		for (int i = 0; i < DIM; i++) { V[i] = xk1[i] - xk[i]; }
		norm_new = norm(DIM, V, NORM);
		copy(DIM, xk, xk1);
	} while (norm_new > EPS);
	cout << "Метод Якоби выполнен \n количество итераций = " << iter << endl;
	delete[] xk1;
	delete[] V;
}

void Seidel(const unsigned int DIM, mytipe** const A, const mytipe* const b, const char NORM, mytipe EPS, const mytipe normC, mytipe* xk) //метод Зейделя
{
	cout << "\n\t\t\t-------------\n\t\t\tМЕТОД ЗЕЙДЕЛЯ\n\t\t\t-------------\n ";
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
		cout << "норма погрешности на " << iter << "-ой итерации:  " << norm(DIM, V, NORM) << endl;

		for (int i = 0; i < DIM; i++) { V[i] = xk1[i] - xk[i]; }
		norm_new = norm(DIM, V, NORM);
		copy(DIM, xk, xk1);
	} while (norm_new > EPS);

	cout << "Метод Зейделя выполнен \n количество итераций = " << iter << endl;
	delete[] xk1;
	delete[] V;
}

void Relax(const unsigned int DIM, mytipe** const A, const mytipe* const b, const mytipe W, const char NORM, mytipe EPS, const mytipe normC, mytipe* xk) //метод релаксации
{
	cout << "\n\t\t\t----------------\n\t\t\tМЕТОД РЕЛАКСАЦИИ\n\t\t\t----------------\n ";
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
		cout << "норма погрешности на " << iter << "-ой итерации:  " << norm(DIM, V, NORM) << endl;
		for (int i = 0; i < DIM; i++) { V[i] = xk1[i] - xk[i]; }
		norm_new = norm(DIM, V, NORM);
		copy(DIM, xk, xk1);
	} while (norm_new > EPS);

	cout << "Метод Релаксации выполнен \n количество итераций = " << iter << endl;
	delete[] xk1;
	delete[] V;
}

void C_Iter(const unsigned int DIM, mytipe** const A, const mytipe* const b, mytipe** const C, const mytipe* const y, const char NORM) //решение итерационных методов, используя С и y
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
	cout << n << " итераций\n";
	delete[] xk1;
}

void relax3d(const int var, const char NORM, const mytipe w, const mytipe EPS)
{
	const unsigned int DIM = 200 + var;  //размер пространнства сделали константным
	int iter = 0;
	mytipe *a;  //столбец, записанный из файла
	a = new mytipe[DIM];     //|выделили место под столбец b
	mytipe *b;  //столбец, записанный из файла
	b = new mytipe[DIM];     //|выделили место под столбец b
	mytipe *c;  //столбец, записанный из файла
	c = new mytipe[DIM];     //|выделили место под столбец b
	mytipe *d;  //столбец, записанный из файла
	d = new mytipe[DIM];     //|выделили место под столбец b
	mytipe *xk1;                //вектор неизвестных на (к+1)-ом шаге
	xk1 = new mytipe[DIM];
	mytipe *xk;                //вектор неизвестных на (к)-ом шаге
	xk = new mytipe[DIM];//начальное приближение

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

	mytipe *V = new mytipe[DIM];//вектор неизвестных на (к)-ом шаге

	do
	{
		copy(DIM, xk, xk1);
		xk1[0] = w*(d[0] - c[0] * xk[1]) / b[0] + (1.0 - w)*xk[0];
		for (int j = 1; j < DIM - 1; j++) { xk1[j] = w*(d[j] - a[j] * xk1[j - 1] - c[j] * xk[j + 1]) / b[j] + (1 - w)*xk[j]; }
		xk1[DIM - 1] = w*(d[DIM - 1] - a[DIM - 1] * xk1[DIM - 2]) / b[DIM - 1] + (1.0 - w)*xk[DIM - 1];

		for (int i = 0; i < DIM; i++) { V[i] = xk1[i] - xk[i]; }
		iter++;

	} while (norm(DIM, V, NORM) > EPS*(2.0 - w) / w);//(norm(DIM1, V, NORM) / (norm(DIM1, xk, NORM) + EPSILON) > EPS);

	cout << "\n\t\t\t------------------------\n\t\t\t3-х ДИАГОНАЛЬНАЯ МАТРИЦА\n\t\t\t------------------------\n\n ";
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
	cout << "\n\t\t\t\t" << iter << " итераций.\n";

	delete[] a; //удаляем динамический столбец
	delete[] b; //удаляем динамический столбец
	delete[] c; //удаляем динамический столбец
	delete[] d; //удаляем динамический столбец
	delete[] xk; //удаляем динамический столбец
	delete[] xk1; //удаляем динамический столбец
	delete[] V; //удаляем динамический столбец
}

void C_l_d_u(mytipe** C, const int DIM, const char NORM)
{
	mytipe **Cl = new mytipe *[DIM];   //матрица, списанная с файла
	mytipe **Cd = new mytipe *[DIM];   //матрица, списанная с файла
	mytipe **Cu = new mytipe *[DIM];   //матрица, списанная с файла
	for (int i = 0; i < DIM; i++) {   // |место под
		Cl[i] = new mytipe[DIM];       // |матрицу А
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
	cout << "\nматрица Cl=";
	vivod(DIM, Cl);  //выводим матрицу C

	cl = norm(DIM, Cl, NORM);
	cout << "норма матрицы Cl=" << cl << endl;

	cout << "\nматрица Cd=";
	vivod(DIM, Cd);  //выводим матрицу C

	cout << "норма матрицы Cd=" << norm(DIM, Cd, NORM) << endl;

	cout << "\nматрица Cu=";
	vivod(DIM, Cu);  //выводим матрицу C

	cu = norm(DIM, Cu, NORM);
	cout << "норма матрицы Cu=" << cu << "\n\n||Cl||+||Cu|| = " << cu + cl;

	for (int i = 0; i < DIM; i++) //удаляем динамические строки матрицы
	{
		delete[] Cl[i]; delete[] Cd[i]; delete[] Cu[i];
	}
	delete[] Cl; delete[] Cd; delete[] Cu;
}

