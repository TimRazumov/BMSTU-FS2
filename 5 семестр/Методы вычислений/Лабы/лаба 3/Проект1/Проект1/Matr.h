#pragma once//директива для того, чтобы файл подкл строго 1 раз

#include <iostream>

typedef double mytipe;


using namespace std;
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

void vivod(const unsigned int DIM, const mytipe* const b, mytipe** const a) //процедура вывода матрица и столбца
{
	cout << endl;
	for (int i = 0; i < DIM; i++) {
		for (int j = 0; j < DIM; j++) {
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

