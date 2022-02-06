#pragma once//директива для того, чтобы файл подкл строго 1 раз

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


void copy(const unsigned int DIM, mytipe** B, mytipe** A) //копирование матриц
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

mytipe** multimatrix(const unsigned int DIM, mytipe** const a, mytipe** const b) //перемножение матриц
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
			t = 0.0;
			for (int j = 0; j<DIM; j++) {
				t += a[i][j] * b[j][l];  //перемножение строки матрицы a на столбец матрицы b
			}
			E[i][l] = t;   //запись результата перемножения в результирующую матрицу
		}
	}

	return E;
}

mytipe** mult_matrnum(const unsigned int DIM, mytipe** const A, const mytipe a) //умножене матрицы на число
{
	mytipe** E = new mytipe *[DIM];
	for (int i = 0; i < DIM; i++)       // под матрицу
	{                                   //  являющуюся результатом
		E[i] = new mytipe[DIM];         //  перемножения матриц a и b
	}

	for (int i = 0; i<DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			E[i][j] = a * A[i][j];
		}
	}
	return E;
}


mytipe* matrvec(const unsigned int DIM, mytipe** A, mytipe* b) //умножение матрицы на столбец
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
	return c;
}



mytipe* mult_vectnum(const unsigned int DIM, const mytipe* const b, const mytipe a) //умножение вектора на число
{
	mytipe* c;
	c = new mytipe[DIM];
	for (int i = 0; i<DIM; i++) {
		c[i] = a * b[i];
	}
	return c;
}

mytipe* sum_vect(const unsigned int DIM, const mytipe* const b1, const mytipe* const b2) //сложение векторов
{
	mytipe* c;
	c = new mytipe[DIM];
	for (int i = 0; i < DIM; i++) {
		c[i] = b1[i] + b2[i];
	}
	return c;
}

mytipe** sum_matr(const unsigned int DIM, mytipe** A1, mytipe** A2) //сложение матриц
{
	mytipe** E;
	E = new mytipe *[DIM];              // выделяем место
	for (int i = 0; i < DIM; i++)       // под матрицу
	{                                   //  являющуюся результатом
		E[i] = new mytipe[DIM];         //  перемножения матриц a и b
	}
	for (int i = 0; i<DIM; i++) {
		for (int j = 0; j < DIM; j++) {
			E[i][j] = A1[i][j] + A2[i][j];
		}
	}
	return E;
}



void singlematr(const unsigned int DIM, mytipe** E) //заполнение матрицы единичной
{
	for (int i = 0; i < DIM; i++)      //|делаем
	{									//|матрицу Т
		for (int j = 0; j < DIM; j++)  //|единичной
		{							    //|для
			E[i][j] = 0.;				//|дальнейших
		}								//|преобразований
		E[i][i] = 1.;					//|
	}
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
	mytipe EPS = 1e-20;

	for (int k = 0, h; k < (DIM); k++)
	{
		eps = abs(a[k][k]);  //запись диагонального элемента в переменную сравнения
		vspom = a[k];   //запись ссылки на строку диагонального элемента(если условие в цикле ни разу не выполнится)
		h = k;   //запись номера строки (по причине, описанной выше)

		if (abs(a[k][k]) < EPS) { flag = true; } //если первый элемент нулевой, то активировать флаг
		else { flag = false; } //если не нулевой, то деактивировать флаг

		for (int i = (k + 1); i < DIM; i++)
		{
			if (abs(a[i][k]) > eps && abs(a[i][k])>EPS) //если элемент больше предыдущего max в столбце и не ноль
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
		else { cout << "Матрица вырождена\n";   return flag; }

		for (int i = k + 1; i < DIM; i++)
		{
			b[i] = b[i] - b[k] / a[k][k] * a[i][k]; //вычитание из элементов столбца
			for (int j = (DIM - 1); j >= 0; j--)
			{
				a[i][j] = a[i][j] - a[k][j] / a[k][k] * a[i][k]; //зануление всех элементов под диагональным
			}
		}

	}
	return flag;
}

void QR(const unsigned int DIM, mytipe** T, mytipe** R, mytipe** A) //QR-алгоритм
{
	mytipe vspom;   //вспомогательный указатель для смены элементов
	mytipe c;  //вычисляемый коэффициент Т_ij
	mytipe s;  //вычисляемый коэффициент Т_ij
	copy(DIM, R, A);
	singlematr(DIM, T);
	for (int k = 0; k < DIM; k++)  //нахождение матриц T=Q^(-1) и R
								   //высчитывание новых T_ij, спускаясь по главной диагонали
	{
		for (int i = k + 1; i < DIM; i++)  //высчитывание новых T_ij, изменяя строки c ненулевого элемета
		{
			if (abs(R[i][k]) > EPSILON)
			{
				c = R[k][k] / sqrt((R[k][k])*(R[k][k]) + (R[i][k])*(R[i][k]));//вычисление c
				s = R[i][k] / sqrt((R[k][k])*(R[k][k]) + (R[i][k])*(R[i][k]));//вычисление s
				for (int j = 0; j < DIM; j++)  //умножение T_ij на A
				{
					vspom = T[k][j];
					T[k][j] = c * T[k][j] + s * T[i][j]; //умножение двух строк Т
					T[i][j] = c * T[i][j] - s * vspom;  //на столбец предыдущей Т
				}
				for (int j = k; j < DIM; j++)
				{
					vspom = R[k][j];
					R[k][j] = c * R[k][j] + s * R[i][j]; //умножение двух строк Т
					R[i][j] = c * R[i][j] - s * vspom;    //на столбец А
				}
			}

		}
	}
	trans(DIM, T);//МАТРИЦА Q
}

mytipe skalar(int N, mytipe* a, mytipe* b)
{
	mytipe r=0.0;
	for (int i = 0; i < N; i++) { r += a[i] * b[i]; };
	return r;
}

mytipe* proj(int N, mytipe* b, mytipe* a)
{
	return mult_vectnum(N, b, skalar(N, a, b) / skalar(N, b, b));
}