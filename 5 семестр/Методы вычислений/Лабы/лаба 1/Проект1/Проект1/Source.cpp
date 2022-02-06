#include "iostream"
#include <fstream>
#include <typeinfo>
#include <cmath>

using namespace std;

const double EPSILON = 10e-8;   //константа сравнения с 0
typedef double mytipe;  //тип данных, использующийся во всей программе
mytipe **A; //матрица, списанная с файла
mytipe *B; //столбец, записанный из файла


void vivod(const unsigned int DIM, const mytipe* const b, mytipe** const a); //вывод матрицы и столбца
void vivod(const unsigned int DIM, mytipe** const a);  //вывод матрицы
void vivod(const unsigned int DIM, const mytipe* const b);  //вывод вектора
void hod(const unsigned int DIM, mytipe* b, mytipe** a); //обратный ход Гаусса
bool findx(const unsigned int DIM, mytipe* b, mytipe** a);  //прямой ход Гаусса с выходом выражденности матрицы
void neviaz(const unsigned int DIM, const mytipe* const x);  //невязка
void copy(const unsigned int DIM, mytipe* b, mytipe** a); //копирование матрицы и вектора
void copy(const unsigned int DIM, mytipe** a, mytipe** b); //копирование матриц
void copy(const unsigned int DIM, mytipe* a, mytipe* b); //копирование векторов
void trans(const unsigned int DIM, mytipe** a);  // транспонирование матрицы
void multimatrix(const unsigned int DIM, mytipe** const a, mytipe** const b); // перемножение матриц без возврата результата
mytipe norm(const unsigned int DIM, mytipe** const A, char flag); // 1 из 4 норм матрицы
mytipe norm(const unsigned int DIM, const mytipe* const b, char flag);  // 1 из 3 норм вектора
void obusl(const unsigned int DIM, const mytipe* const x); //оценка снизу числа обусловленности
void singlematr(const unsigned int DIM, mytipe** const E);  //заполнение матрицы единичной
void QR(const unsigned int DIM, mytipe** T, mytipe** R); //QR-алгоритм
void matrvec(const unsigned int DIM, mytipe** A, mytipe* b);  //умножение матрицы на столбец
int main() {
	setlocale(LC_ALL, "Russian"); //подключаем русский язык
	ifstream fin;    // создали переменную для считывания из файла
	fin.open("6(1).txt", ios_base::in | ios_base::app | ios_base::binary);  //открыли файл
	unsigned int a;  // размер пространства
	fin >> a;  //считали из файла размер пространства
	cout << "Dim=" << a << endl;
	const unsigned int DIM1 = a;  //размер пространнства сделали константным

	mytipe **ary;                      // |создали 
	ary = new mytipe *[DIM1];          // |динамический 
	for (int i = 0; i < DIM1; i++) {   // |массив 
		ary[i] = new mytipe[DIM1];     // |матрицы А
	}

	mytipe *b;                //|столбец
	b = new mytipe[DIM1];     //|правой части


	A = new mytipe *[DIM1];            // |выделили  
	for (int i = 0; i < DIM1; i++) {   // |место под 
		A[i] = new mytipe[DIM1];       // |матрицу А
	}

	B = new mytipe[DIM1];     //|выделили место под столбец b

	for (int i = 0; i < DIM1; i++)  //считываем из файла матрицу и столбец
	{
		for (int j = 0; j < DIM1; j++)
		{
			fin >> A[i][j];  //матрицу
		}
		fin >> B[i];  //столбец
	}

	fin.close();  //закрываем файл
	fin.clear();  //отвязываем файл от файловой переменной

	copy(DIM1, b, ary);  //копирование исходной матрицы и столбца

	vivod(DIM1, b, ary);  //выводим матрицу и столбец


	bool flag = true;  //для запоминания вырожденности системы

	flag = findx(DIM1, b, ary);  //прямой ход Гаусса и записываем вырожденность системы

	if (flag)  //если система не выраждена
	{
		vivod(DIM1, b, ary); //выводим верхнетреугольную матрицу


		hod(DIM1, b, ary);  // обратный ход Гаусса


		vivod(DIM1, b);  //выводим столбец неизвестных

		neviaz(DIM1, b); //норма невязки
	}


	//            !!!!!! QR-способ !!!!!!

	flag = false;

	mytipe **R;							// |выделили
	R = new mytipe *[DIM1];             // |место
	for (int i = 0; i < DIM1; i++) {    // |под матрицу
		R[i] = new mytipe[DIM1];        // |R
	}


	copy(DIM1, b, R);   //копирование исходной матрицы и столбца

	vivod(DIM1, b, R);  //вывод матрицы и столбца



	mytipe **T;                       // |создаем ортоганальную матрицу Т обратная которой = Q
	T = new mytipe *[DIM1];           // |
	for (int i = 0; i < DIM1; i++) {  // |
		T[i] = new mytipe[DIM1];      // |
	}	                               //|

	singlematr(DIM1, T);

	QR(DIM1, T, R);

	if (abs(R[DIM1 - 1][DIM1 - 1]) < EPSILON) { flag = true; cout << "Матрица вырождена\n\n"; }
	else { matrvec(DIM1, T, b); }

	trans(DIM1, T); //получаем Q из T, транспонируя

	cout << "Q=" << endl;

	vivod(DIM1, T);  //выводим Q

	cout << "R=" << endl;

	vivod(DIM1, R); //выводим R

	cout << "A=Q*R=" << endl;

	multimatrix(DIM1, T, R);
	if (!flag)  //если матрица невырождена
	{
		hod(DIM1, b, R);  //Обратный ход(нахождение вектора х)

		vivod(DIM1, b);  //вывод столбца b


						 //!!!!!!!!!!!!!!невязка и оценка снизу числа обусловленности!!!!!!!!!!

		obusl(DIM1, b);



		//!!!!!!!!!!!!!!Число обусловленности!!!!!!!!!!!


		mytipe** inverseA;                  //выделяем память под обратную матрицу
		inverseA = new mytipe *[DIM1];
		for (int i = 0; i < DIM1; i++) { inverseA[i] = new mytipe[DIM1]; }
		trans(DIM1, T);
		singlematr(DIM1, inverseA);

		for (int r = 0; r < DIM1; r++)
		{
			matrvec(DIM1, T, inverseA[r]);
			hod(DIM1, inverseA[r], R); //находим первый столбец обратной матрицы
		}

		copy(DIM1, ary, A); //копируем данную матрицу в рабочую
		trans(DIM1, inverseA); //транспонируем обратную матрицу
		cout << "A^-1=";
		vivod(DIM1, inverseA);  //выводим обратную матрицу
		cout << "A*A^-1=";
		multimatrix(DIM1, inverseA, ary); //перемножаем обратную и прямую матрицы
										  //       вывод числа обусловленности через:
		cout << "\ncond1 A = " << norm(DIM1, ary, '1')*norm(DIM1, inverseA, '1') << endl; //октоэдрическую норму
		cout << "\ncond2 A = " << norm(DIM1, ary, '2')*norm(DIM1, inverseA, '2') << endl; //евклидову норму
		cout << "\ncond(inf) A = " << norm(DIM1, ary, 'k')*norm(DIM1, inverseA, 'k') << endl; //кубическую норму
		cout << "\ncond(max) A = " << norm(DIM1, ary, 'm')*norm(DIM1, inverseA, 'm') << endl; //max-ую норму

		for (int i = 0; i < DIM1; i++) //удаляем динамические строки матрицы
		{
			delete[] inverseA[i];
		}
		delete[] inverseA;//удаляем динамическую матрицу
	}
	else { cout << "\ncond A = infinity" << endl; }


	for (int i = 0; i < DIM1; i++) //удаляем динамические строки матрицы
	{
		delete[] ary[i];
	}
	delete[] ary;//удаляем динамическую матрицу
	for (int i = 0; i < DIM1; i++) //удаляем динамические строки матрицы
	{
		delete[] T[i];
	}
	delete[] T;//удаляем динамическую матрицу

	for (int i = 0; i < DIM1; i++) //удаляем динамические строки матрицы
	{
		delete[] R[i];
	}
	delete[] R;//удаляем динамическую матрицу
	delete[] b; //удаляем динамический столбец

	for (int i = 0; i < DIM1; i++) //удаляем динамические строки матрицы
	{
		delete[] A[i];
	}
	delete[] A;//удаляем динамическую матрицу
	delete[] B; //удаляем динамический столбец

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

void hod(const unsigned int DIM, mytipe* b, mytipe** a) //обратный ход Гаусса
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

bool findx(const unsigned int DIM, mytipe* b, mytipe** a) //прямой ход Гаусса с проверкой на выражденность
{
	mytipe eps;     //вспомогательная переменная для сравнения элементов (нахождения max)
	mytipe* vspom;  //вспомогательная ссылка для перестановки строк матрицы
	mytipe vspom2;   //вспомогательная переменная для перестановки элементов вектора
	bool flag = true;  //логическая переменная(проверка вырожденности матрицы)


	for (int k = 0, h; k < (DIM); k++)
	{
		eps = abs(a[k][k]);  //запись диагонального элемента в переменную сравнения
		vspom = a[k];   //запись ссылки на строку диагонального элемента(если условие в цикле ни разу не выполнится)
		h = k;   //запись номера строки (по причине, описанной выше)

		if (abs(a[k][k]) < EPSILON) { flag = true; } //если первый элемент нулевой, то активировать флаг
		else { flag = false; } //если не нулевой, то деактивировать флаг

		for (int i = (k + 1); i < DIM; i++)
		{
			if (abs(a[i][k]) > eps && abs(a[i][k])>EPSILON) //если элемент больше предыдущего max в столбце и не ноль
			{
				eps = a[i][k]; vspom = a[i];  h = i; flag = false;
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
void neviaz(const unsigned int DIM, const mytipe* const x) //невязка
{
	mytipe* b1;          //выделение памяти под
	b1 = new mytipe[DIM];  //вектор, который получится при подстановке решения в систему


	for (int i = 0; i < DIM; i++)
	{
		b1[i] = 0;
		for (int j = 0; j < DIM; j++)
		{
			b1[i] += A[i][j] * x[j]; //подсчет, подстановкой решения, столбца правой части

		}
		b1[i] -= B[i]; //разность i-ой компоненты получившегося вектора и истенного
	}
	cout << "\n||b1- b||1 = " << norm(DIM, b1, '1') << endl;  //вывод октоэдрической нормы (невязка)
	cout << "\n||b1- b||2 = " << norm(DIM, b1, '2') << endl;  //вывод евклидовой нормы (невязка)
	cout << "\n||b1- b||(inf) = " << norm(DIM, b1, 'k') << endl;  //вывод кубической нормы (невязка)

	delete[] b1; //освобождение памяти динамического вектора
}

void copy(const unsigned int DIM, mytipe* b, mytipe** a) //копирование матрицы и вектора
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

void copy(const unsigned int DIM, mytipe** a, mytipe** b) //копирование матриц
{

	for (int i = 0; i < DIM; i++)
	{
		for (int j = 0; j < DIM; j++)
		{
			a[i][j] = b[i][j];
		}
	}
}

void copy(const unsigned int DIM, mytipe* a, mytipe* b) //копирование векторов
{
	for (int j = 0; j < DIM; j++)
	{
		a[j] = b[j];
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


mytipe norm(const unsigned int DIM, mytipe** const A, char flag) //4 нормы матрицы с запросом варианта нормы
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

mytipe norm(const unsigned int DIM, const mytipe* const b, char flag) //3 нормы вектора с запросом варианта нормы
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

void obusl(const unsigned int DIM, const mytipe* const x) //оценка снизу числа обусловленности
{
	neviaz(DIM, x);  //невязка

	mytipe* b;             //|выделение памяти под найденное нами решение 
	b = new mytipe[DIM];   //|

	mytipe **ary;                     // |создали 
	ary = new mytipe *[DIM];          // |динамический 
	for (int i = 0; i < DIM; i++) {   // |массив 
		ary[i] = new mytipe[DIM];     // |матрицы А
	}

	copy(DIM, b, ary);  //записываем исходные матрицу и столбец

	mytipe **T;                       // |создаем ортоганальную матрицу Т обратная которой = Q
	T = new mytipe *[DIM];           // |
	for (int i = 0; i < DIM; i++) {  // |
		T[i] = new mytipe[DIM];      // |
	}	                               //|

	singlematr(DIM, T);

	QR(DIM, T, ary);

	bool flag;  //флаг на вырожденность системы
	mytipe* deltaX;            //|выделяем память под
	deltaX = new mytipe[DIM];  //|вектор разности точного решения и с возмущением
	const mytipe VOZM = 0.01;  //возмущение
	mytipe dx;   //относительная погрешность x
	mytipe db1 = VOZM / norm(DIM, b, '1');  //относительная погрешность b по октоэдрической норме
	mytipe db2 = VOZM / norm(DIM, b, '2');  //относительная погрешность b по евклидовой норме
	mytipe dbk = VOZM / norm(DIM, b, 'k');  //относительная погрешность b по кубической норме
	mytipe norm1 = norm(DIM, x, '1');  //октоэдрическая норма вектора х
	mytipe norm2 = norm(DIM, x, '2');  //евклидова норма вектора х
	mytipe normk = norm(DIM, x, 'k');  //кубическая норма вектора х

	mytipe max1 = 0; //максимальная оценка числа обусловленности снизу по октоэдрической норме
	mytipe max2 = 0; //максимальная оценка числа обусловленности снизу по евклидовой норме
	mytipe maxk = 0; //максимальная оценка числа обусловленности снизу по кубической норме

					 //решаем возмущенную систему методом Гаусса
	for (int i = 0; i < 2 * DIM; i++)
	{
		if (i < DIM) { b[i] += VOZM; cout << i + 1 << ") b*="; vivod(DIM, b); } //прибавляем возмущение
		else { b[i - DIM] -= VOZM; cout << i - DIM + 1 << ") b*="; vivod(DIM, b); } //вычитаем возмущение

		flag = findx(DIM, b, ary); //приводим к верхнетреугольному виду с определением вырожденности

		matrvec(DIM, T, b);

		if (flag)  //если мтарица невырождена
		{
			hod(DIM, b, ary);  //обратный ход Гаусса
		}
		else { cout << "Матрица вырождена => оценка числа обусловленности невозможна"; cin.get(); break; }

		for (int j = 0; j < DIM; j++) { deltaX[j] = abs(x[j] - b[j]); } //вектор разности

		cout << " x* = "; vivod(DIM, b); //выводим получившееся решен8ие с возмущением

		copy(DIM, b, B);  //обновление столбца b

		dx = norm(DIM, deltaX, '1') / norm1; //считаем относительную погрешность вектора разнсти решений

		cout << "\n cond1 A >= " << dx / db1; //выводим оценку чила обусловленности по норме 1

		if (max1 < dx / db1) { max1 = dx / db1; } //проверяем максимум 

		dx = norm(DIM, deltaX, '2') / norm2;
		cout << "; cond2 A >= " << dx / db2; //выводим оценку чила обусловленности по норме 2
		if (max2 < dx / db2) { max2 = dx / db2; }

		dx = norm(DIM, deltaX, 'k') / normk;
		cout << "; cond(inf) A >= " << dx / dbk << endl << endl; //выводим оценку чила обусловленности по кубической норме
		if (maxk < dx / dbk) { maxk = dx / dbk; }
	}
	if (max1 != 0) {
		cout << "максимальная оценка числа обусловленности снизу:\n cond1 A >= " << max1 << " cond2 A >= " << max2 << " cond(inf) A >= " << maxk << endl << endl;
	}
	delete[] deltaX;
	delete[] b;
	for (int i = 0; i < DIM; i++) //удаляем динамические строки матрицы
	{
		delete[] ary[i];
	}
	delete[] ary;//удаляем динамическую матрицу

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

void QR(const unsigned int DIM, mytipe** T, mytipe** R) //QR-алгоритм
{
	mytipe vspom2;   //вспомогательный указатель для смены элементов
	mytipe c;  //вычисляемый коэффициент Т_ij
	mytipe s;  //вычисляемый коэффициент Т_ij

	for (int k = 0; k < DIM; k++)  //нахождение матриц T=Q^(-1) и R
								   //высчитывание новых T_ij, спускаясь по главной диагонали
	{
		for (int i = k + 1; i < DIM; i++)  //высчитывание новых T_ij, изменяя строки c ненулевого элемета
		{
			if (abs(R[i][k])>EPSILON)
			{
				c = R[k][k] / sqrt((R[k][k])*(R[k][k]) + (R[i][k])*(R[i][k]));//вычисление c
				s = R[i][k] / sqrt((R[k][k])*(R[k][k]) + (R[i][k])*(R[i][k]));//вычисление s
				for (int j = 0; j < DIM; j++)  //умножение T_ij на A
				{
					vspom2 = T[k][j];
					T[k][j] = c * T[k][j] + s * T[i][j]; //умножение двух строк Т
					T[i][j] = c * T[i][j] - s * vspom2;  //на столбец предыдущей Т
				}
				for (int j = k; j < DIM; j++)
				{
					vspom2 = R[k][j];
					R[k][j] = c * R[k][j] + s * R[i][j]; //умножение двух строк Т
					R[i][j] = c * R[i][j] - s * vspom2;    //на столбец А
				}
			}

		}
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