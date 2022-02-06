#include "iostream"
#include <fstream>
#include <typeinfo>
#include <cmath>


using namespace std;

const double EPSILON = 1e-14;   //константа сравнения с 0

const double PI = 3.14159265;

typedef double mytipe;  //тип данных, использующийся во всей программе

const mytipe Eps = 1e-6; //точность нахождения корней

typedef mytipe(*func)(const mytipe x);
typedef mytipe*(*func_sys)(const mytipe* x);


mytipe** Yakobi_Matr(mytipe* x, int m);//матрица Якоби(аналитическая)

void Uniform_grid(const unsigned int DIM, mytipe* x, const mytipe a, const mytipe b);
void Uniform_grid2(const mytipe h, mytipe** x, const mytipe a1, const mytipe b1, const mytipe a2, const mytipe b2);
void Func_value(const unsigned int DIM, mytipe* x, mytipe* f);
mytipe D(func f, mytipe x);
mytipe** Yakobi_Matr(func_sys f, mytipe* x); //Вычисление матрицы Якоби

mytipe* Local(mytipe a, mytipe b, int& N);
mytipe* Vilka(int N, mytipe* x);
mytipe* Newton(int N, mytipe* x);
mytipe* Hord(int N, mytipe* x);
mytipe Newton(func_sys f, mytipe* x0, mytipe* result); //Для систем

void Copy_to_file(const unsigned int n, mytipe** x, const int* const iter);
mytipe Newton_pic(func_sys f, mytipe* x0, mytipe* result, const mytipe a1, const mytipe b1, const mytipe a2, const mytipe b2); //n - кол-во решений; x0 - нач приближ

void copy(const unsigned int DIM, mytipe* b, const mytipe* a); //копирование векторов
void copy(const unsigned int DIM, mytipe** B, mytipe** const A); //копирование матриц
void vivod(const unsigned int DIM, mytipe** const a); //процедура вывода матрицы
void vivod(const unsigned int DIM, const mytipe* const b); //процедура вывода столбца
mytipe norm(const unsigned int DIM, const mytipe* const b, const char flag); //3 нормы вектора с запросом варианта нормы
void mult_matrvect(const unsigned int DIM, mytipe** const A, const mytipe* const b, mytipe* x); //умножение матрицы на вектор
void minus_vect(const unsigned int DIM, const mytipe* const b1, const mytipe* const b2, mytipe* x); //разность векторов


mytipe Func(const mytipe x)
{
	//return (x - 0.1)*(x - 0.22)*(x - 0.55)*(x - 0.7)*(x - 0.75);
	return sqrt(x + 1) - 1;
	//return 35 * x*x*x - 67 * x*x - 3 * x + 3;
	//return (x - 1)*(x - 1);
	//return x*x - 1;
	
}

mytipe FuncDot(const mytipe x)
{
	//return 5 * (0.024299 - 0.30238*x + 1.1907*x*x - 1.856*x*x*x + x*x*x*x);
	return 1/(2*sqrt(x+1));
	//return 3 * 35 * x*x - 67 * 2 * x - 3;
	
}

mytipe* func_system(const mytipe* x)
{
	mytipe* y = new mytipe[2];

	//y[0] = x[0] * x[0] - x[1] * x[1] - 15; //1
	//y[1] = x[0] * x[1] + 4;

	y[0] = x[0] * x[0] + x[1] * x[1] + x[0] + x[1] - 8; //2
	y[1] = x[0] * x[0] + x[1] * x[1] + x[0] * x[1] - 7;

	return y;
}

int main() {
	setlocale(LC_ALL, "Russian"); //подключаем русский язык


	mytipe a = -1;//отрезок
	mytipe b = 10;
	mytipe h;  //шаг разбиений
	cout << FuncDot(2) << " " << D(Func, 2);
	cout << "Введите шаг разбиения: "; cin >> h; cout << endl; cin.get();
	int N = round((b - a) / h) + 1;//кол-во узлов(кол-во разбиений N-1, шаг (b-a)/(N-1) )
	mytipe* x;
	mytipe* Solve;

	x = Local(a, b, N);
	//Solve = Vilka(N, x);
	Solve = Newton(N, x);
	//Solve = Hord(N, x);
	cout << "Корни: \n";
	vivod(N, Solve);

	//const int n = 2;// кол-во решений 
	//mytipe x1[2 * n]; x1[0] = -5; x1[1] = -3; x1[2] = 3; x1[3] = 5;//отрезки локализации (определены визуально)
	//mytipe x2[2 * n]; x2[0] = 0; x2[1] = 2; x2[2] = -2; x2[3] = 0;



	mytipe a1 = -10;//отрезок по x
	mytipe b1 = 10;
	mytipe a2 = -10;//отрезок по y
	mytipe b2 = 10;
	mytipe step;

	cout << "Введите шаг разбиения для системы уравнений: "; cin >> step; cout << endl; cin.get();

	int num1 = (b1 - a1) / step + 1;
	int num2 = (b2 - a2) / step + 1;


	mytipe** Sol = new mytipe*[num1*num2];
	for (int i = 0; i < num1*num2; i++) {
		Sol[i] = new mytipe[2];
	}
	mytipe** initial = new mytipe*[num1*num2];
	for (int i = 0; i < num1*num2; i++) {
		initial[i] = new mytipe[2];
	}

	int* count_iter = new int[num1*num2];

	Uniform_grid2(step, initial, a1, b1, a2, b2);

	for (int i = 0; i < num1*num2; i++) {
		count_iter[i] = Newton_pic(func_system, initial[i], Sol[i], a1, b1, a2, b2);
	}

	//cout << "Итерации: ";
	//for (int i = 0; i < num1*num2; i++){
	//	cout << "Количество итераций в узле " << i << " = " << count_iter[i] << endl;vivod(2, Sol[i]);
	//}

	Copy_to_file(num1*num2, initial, count_iter);
	//cout << "\nКорни нелинейного уравнения: \n";
	//for (int i = 0; i < n; i++) { cout << "X" << i + 1 << " = "; vivod(n, Sol[i]); };

	delete[] x;
	delete[] Solve;
	cout << "The End";
	cin.get();
	return 0;
}

void Copy_to_file(const unsigned int n, mytipe** x, const int* const iter)
{
	ofstream fout;    // создали переменную для записи в файл
	fout.open("iter.txt", ios_base::out | ios_base::trunc);
	for (unsigned int i = 0; i < n; i++) {
		fout << x[i][0] << " " << x[i][1] << " " << iter[i] << endl;
	}
	fout.close();
	fout.clear();
}

mytipe* Local(mytipe a, mytipe b, int& N)
{
	mytipe* x = new mytipe[N];         //ТАБЛИЦА ДЛЯ ЛОКАЛИЗАЦИИ КОРНЕЙ
	Uniform_grid(N - 1, x, a, b);

	mytipe* f = new mytipe[N];
	Func_value(N - 1, x, f);

	mytipe* vspom = new mytipe[2 * N];//максимально возможное кол-во корней
	int j = 0;

	for (int i = 1; i < N; i++)
	{
	
		 //if (fabs(f[i]) < EPSILON) 
		 //{
			 //vspom[j] = x[i]; vspom[j + 1] = x[i];   j += 2; i++;
		//if ((f[i] <EPSILON && f[i+1]>EPSILON) || (f[i] > EPSILON && f[i+1] < EPSILON)) {}
		//else {i++;}
		//}

		if ((f[i - 1] <EPSILON && f[i]>EPSILON) || (f[i - 1] > EPSILON && f[i] < EPSILON))
		{
			vspom[j] = x[i - 1]; vspom[j + 1] = x[i]; j += 2;
		}
	}
	if (fabs(f[N - 1]) < EPSILON) {
		vspom[j] = x[N - 1]; vspom[j + 1] = x[N - 1]; j += 2;
	}

	mytipe* Local = new mytipe[j];
	copy(j, Local, vspom);
	N = j / 2;
	cout << "Локализованно " << N << " корней\n\nОтрезки локализации:\n";
	for (int k = 0; k < j - 1; k += 2)  cout << "[" << Local[k] << ", " << Local[k + 1] << "]\n";

	delete[] x;
	delete[] f;
	delete[] vspom;
	return Local;
}

mytipe* Vilka(int N, mytipe* x)
{

	mytipe a;
	mytipe b;
	mytipe x0;
	mytipe* solve = new mytipe[N];
	int iter = 0;
	int SumIter = 0;

	for (int i = 0; i < N; i++)
	{
		a = x[2 * i];
		b = x[2 * i + 1];
		x0 = (a + b) / 2;
		while (b - a > 2 * Eps)
		{
			if (fabs(Func(a)) < EPSILON) { x0 = a; iter++; break; }
			else if (fabs(Func(b)) < EPSILON) { x0 = b; iter++; break; }
			else if (fabs(Func(x0)) < EPSILON) { iter++; break; }

			if (Func(a)*Func(x0) <= EPSILON)  b = x0;
			else  a = x0;
			x0 = (a + b) / 2;
			iter++;
		}
		solve[i] = x0;
		SumIter += iter;
		cout << iter << " + ";
		iter = 0;
	}

	cout << " = " << SumIter << endl;
	return solve;
}

mytipe* Newton(int N, mytipe* x)
{
	mytipe xk;
	mytipe xk1;
	mytipe b;
	mytipe a;
	//N = 1;
	mytipe* solve = new mytipe[N];
	int iter = 0;
	int SumIter = 0;
	
	for (int i = 0; i < N; i++)
	{
		//a = -1; b = 1;
		b = x[2 * i + 1];
		a = x[2 * i];
		xk1 = (Func(a)*b - Func(b)*a) / (Func(a) - Func(b));
		//a = 0; b = 1; xk1 = 1.5;
		do
		{
			xk = xk1;
			xk1 = xk - Func(xk) / D(Func, xk);
			//xk1 = xk - Func(xk) / FuncDot(xk);
			iter++;
			cout << xk1 << endl;
		} while (fabs(xk1 - xk) > Eps);
		solve[i] = xk1;
		SumIter += iter;
		cout << iter << " + ";
		iter = 0;
	}
	cout << " = " << SumIter << endl;
	return solve;
}

mytipe* Hord(int N, mytipe* x)
{
	mytipe xk;
	mytipe b;
	mytipe a;
	mytipe xk1;
	mytipe* solve = new mytipe[N];
	int iter = 0;
	int SumIter = 0;

	for (int i = 0; i < N; i++)
	{
		b = x[2 * i + 1];
		a = x[2 * i];
		xk1 = (Func(a)*b - Func(b)*a) / (Func(a) - Func(b));
		do
		{
			xk = xk1;
			xk1 = xk - Func(xk)*(b - xk) / (Func(b) - Func(xk));
			iter++;
		} while (fabs(xk1 - xk) > Eps);
		solve[i] = xk1;
		SumIter += iter;
		cout << iter << " + ";
		iter = 0;
	}
	cout << " = " << SumIter << endl;
	return solve;
}

mytipe Newton_pic(func_sys f, mytipe* x0, mytipe* result, const mytipe a1, const mytipe b1, const mytipe a2, const mytipe b2) //n - кол-во решений; x0 - нач приближ
{
	int dim = 2;//размерность системы

	mytipe* xk = new mytipe[dim];
	mytipe* xk1 = new mytipe[dim];//нач приближ
	mytipe** Yakob = new mytipe*[2]; Yakob[0] = new mytipe[2]; Yakob[1] = new mytipe[2];
	mytipe vspom;
	mytipe yakobian;
	int iter = 0;

	xk1[0] = x0[0]; xk1[1] = x0[1];
	do
	{
		iter++;
		copy(dim, xk, xk1);
		copy(dim, Yakob, Yakobi_Matr(f, xk));
		//copy(dim, Yakob, Yakobi_Matr(xk, 4));
		//vivod(dim, x0);
		//vivod(dim,Yakob);
		yakobian = Yakob[0][0] * Yakob[1][1] - Yakob[1][0] * Yakob[0][1];
		vspom = Yakob[0][0];
		Yakob[0][0] = Yakob[1][1] / yakobian;
		Yakob[1][1] = vspom / yakobian;
		Yakob[0][1] /= -yakobian;
		Yakob[1][0] /= -yakobian;

		mult_matrvect(dim, Yakob, f(xk), xk1);
		minus_vect(dim, xk, xk1, xk1);
		if (xk1[0]<a1 || xk1[0]>b1 || xk1[1]<a2 || xk1[1]>b2) { iter = 1000; break; }
		minus_vect(dim, xk, xk1, xk);
	} while (norm(dim, xk, '2') > Eps);
	result[0] = xk1[0];
	result[1] = xk1[1];

	return iter;
}

mytipe Newton(func_sys f, mytipe* x0, mytipe* result) //n - кол-во решений; x0 - нач приближ
{
	int dim = 2;//размерность системы

	mytipe* xk = new mytipe[dim];
	mytipe* xk1 = new mytipe[dim];//нач приближ
	mytipe** Yakob = new mytipe*[2]; Yakob[0] = new mytipe[2]; Yakob[1] = new mytipe[2];
	mytipe vspom;
	mytipe yakobian;
	int iter = 0;

	xk1[0] = x0[0]; xk1[1] = x0[1];
	do
	{
		iter++;
		copy(dim, xk, xk1);
		copy(dim, Yakob, Yakobi_Matr(f, xk));
		//copy(dim, Yakob, Yakobi_Matr(xk, 4));
		//vivod(dim, x0);
		//vivod(dim,Yakob);
		yakobian = Yakob[0][0] * Yakob[1][1] - Yakob[1][0] * Yakob[0][1];
		vspom = Yakob[0][0];
		Yakob[0][0] = Yakob[1][1] / yakobian;
		Yakob[1][1] = vspom / yakobian;
		Yakob[0][1] /= -yakobian;
		Yakob[1][0] /= -yakobian;

		mult_matrvect(dim, Yakob, f(xk), xk1);
		minus_vect(dim, xk, xk1, xk1);
		minus_vect(dim, xk, xk1, xk);
	} while (norm(dim, xk, '2') > Eps);
	result[0] = xk1[0];
	result[1] = xk1[1];

	return iter;
}


void Uniform_grid(const unsigned int DIM, mytipe* x, const mytipe a, const mytipe b)
{
	mytipe h = abs(b - a) / DIM;
	for (int i = 0; i <= DIM; i++) {
		x[i] = a + h*i;
	}
}

void Uniform_grid2(const mytipe h, mytipe** x, const mytipe a1, const mytipe b1, const mytipe a2, const mytipe b2)
{
	int i = 0;
	x[0][0] = a1;
	x[0][1] = a2;
	while (x[i][0] + h <= b1) {
		x[i + 1][0] = x[i][0] + h;
		x[i + 1][1] = x[i][1];
		i++;
	}
	while (x[i][1] + h <= b2) {
		x[i + 1][0] = a1;
		x[i + 1][1] = x[i][1] + h;
		i++;
		//cout << i << "  "; vivod(2, x[i]);
		while (x[i][0] + h <= b1) {
			x[i + 1][0] = x[i][0] + h;
			x[i + 1][1] = x[i][1];
			i++;
			//cout << i << "  "; vivod(2, x[i]);
		}

	}
}


void Func_value(const unsigned int DIM, mytipe* x, mytipe* f)
{

	for (int i = 0; i <= DIM; i++) f[i] = Func(x[i]);

}

mytipe D(func f, mytipe x)//численное вычисление производной в точке
{
	mytipe eps1 = 1e-10;
	return (f(x + eps1) - f(x)) / eps1;
}


mytipe** Yakobi_Matr(func_sys f, mytipe* x)//Вычисление матрицы Якоби
{
	mytipe eps1 = 1e-10;
	mytipe vspom1[2] = { x[0] + eps1, x[1] };
	mytipe vspom2[2] = { x[0], x[1] + eps1 };
	mytipe** Yakob = new mytipe*[2]; Yakob[0] = new mytipe[2]; Yakob[1] = new mytipe[2];
	mytipe* vspom;
	mytipe* fx;
	vspom = f(vspom1);
	fx = f(x);

	Yakob[0][0] = (vspom[0] - fx[0]) / eps1;
	Yakob[1][0] = (vspom[1] - fx[1]) / eps1;

	delete[] vspom;
	vspom = f(vspom2);

	Yakob[0][1] = (vspom[0] - fx[0]) / eps1;
	Yakob[1][1] = (vspom[1] - fx[1]) / eps1;

	delete[] vspom;
	delete[] fx;
	return Yakob;
}

mytipe** Yakobi_Matr(mytipe* x, int m)//матрица Якоби(аналитическая)
{
	mytipe** Y = new mytipe*[2];
	for (int i = 0; i< 2; i++) Y[i] = new mytipe[2];
	if (m == 4) {
		Y[0][0] = 2 * x[0];
		Y[0][1] = -2 * x[1];
		Y[1][0] = x[1];
		Y[1][1] = x[0];
	}
	if (m == 5) {
		Y[0][0] = 2 * x[0] + 1;
		Y[0][1] = 2 * x[1] + 1;
		Y[1][0] = 2 * x[0] + x[1];
		Y[1][1] = 2 * x[1] + x[0];
	}
	return Y;
};

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

void mult_matrvect(const unsigned int DIM, mytipe** const A, const mytipe* const b, mytipe* x) //умножение матрицы на вектор
{
	for (int i = 0; i<DIM; i++) {
		x[i] = 0;
		for (int j = 0; j < DIM; j++) {
			x[i] += A[i][j] * b[j];
		}
	}
}

void minus_vect(const unsigned int DIM, const mytipe* const b1, const mytipe* const b2, mytipe* x) //разность векторов
{
	for (int i = 0; i < DIM; i++) {
		x[i] = b1[i] - b2[i];
	}
}

void copy(const unsigned int DIM, mytipe* b, const mytipe* a) //копирование векторов
{
	for (int j = 0; j < DIM; j++)
	{
		b[j] = a[j];
	}
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