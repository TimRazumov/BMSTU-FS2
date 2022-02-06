#include "iostream"
#include <fstream>
#include <typeinfo>
#include <cmath>


using namespace std;

const double EPSILON = 1e-8;   //константа сравнения с 0
const double PI = 3.14159265;

typedef double mytipe;  //тип данных, использующийся во всей программе

typedef mytipe(*func)(const mytipe x, const mytipe u0, const mytipe L);
typedef mytipe(*Fanc)(const mytipe x, const mytipe t);
typedef mytipe(*Func)(const mytipe x);


void Copy_to_file(const unsigned int n, const mytipe h, mytipe* solve, char file_name[]); //n - кол-во узлов по пространству, h-шаг по пространству 
void Copy_to_file_anim(const unsigned int n, const unsigned int k, mytipe** solve, char file_name[]); //n - кол-во узлов по пространству, k-кол-во узлов по времени
void copy(int dim, mytipe* a, mytipe* b);
void vivod(const unsigned int DIM, const mytipe* const b);
mytipe norm(const unsigned int DIM, const mytipe* const b, const char flag); //3 нормы вектора с запросом варианта нормы
mytipe* minus_vect(const unsigned int DIM, const mytipe* const b1, const mytipe* const b2); //b1-b2

mytipe** Krest(int n, int k, mytipe h, mytipe tau, int test_name);
mytipe** Krest_a(int n, int k, mytipe h, mytipe tau, int test_name);

mytipe* progon3d(int DIM, mytipe* a, mytipe* b, mytipe* c, mytipe* d);

mytipe error(int k, int n, mytipe h, mytipe tao, mytipe** a, int test_name);


mytipe ux0(const mytipe x, int test_name)
{
	if (test_name == 1) { return sin(PI*x); }  //тест1
	else if (test_name == 2) { return x*(1 - x); } //тест2
	else if (test_name == 3) {   //тест3
		const mytipe Eps = 1e-7;  //точность аналитического решения
		const int k = 1 / 2 * (sqrt(2 / pow(PI, 2) / Eps) + 1);
		mytipe sum = 0.;
		for (int n = 0; n <= k; n++) {
			sum += 1 / pow(2 * n + 1, 3)*sin((2 * n + 1)*PI*x);
		}
		sum *= 8 / pow(PI, 3);
		return sum;
	}
	else if (test_name == 4) { 
		if ((x >= 0.0) && (x < 0.25)) { return 0; }
		else if ((x >= 0.25) && (x <= 0.75)) { return 1; }
		else if ((x > 0.75) && (x <= 1.0)) { return 0; }
		else { cout << "Interval error" << endl; }
	} //тест4
	else if (test_name == 13) { return 0.5*(pow(x, 2) + 1); } //тест13
	else if (test_name == 19) { return (x+0.5)*(x+1); } //тест19
	else { cout << "test_name isn't correct" << endl; }
}

mytipe dux0(const mytipe x, int test_name)
{
	if (test_name == 1) { return 0.; }  //тест1
	else if (test_name == 2) { return 0.; } //тест2
	else if (test_name == 3) { return 0.; } //тест3
	else if (test_name == 4) { return 0.; } //тест4
	else if (test_name == 13) { return x*sin(2 * x); } //тест13
	else if (test_name == 19) { return cos(x + 0.5); } //тест19
	else { cout << "test_name isn't correct" << endl; }
}

mytipe ddux0(const mytipe x, int test_name)
{
	if (test_name == 1) { return -PI*PI*sin(PI*x); }  //тест1
	else if (test_name == 2) { return -2; } //тест2
	else if (test_name == 3) { return -2; } //тест3
	else if (test_name == 4) { return 0; } //тест4
	else if (test_name == 13) { return 1; } //тест13
	else if (test_name == 19) { return 2; } //тест19
	else { cout << "test_name isn't correct" << endl; }
}

mytipe ddux0_a(const mytipe x1, const mytipe x2,const mytipe x3, const mytipe h, int test_name)
{
	if (test_name == 1) { return -PI*PI*sin(PI*x2); }  //тест1
	else if (test_name == 2) { return -2; } //тест2
	else if (test_name == 3) { return -2; } //тест3
	else if (test_name == 4) { return (x1-2*x2+x3)/h/h; } //тест4
	else if (test_name == 13) { return 1; } //тест13
	else if (test_name == 19) { return 2; } //тест19
	else { cout << "test_name isn't correct" << endl; }
}

mytipe u_0t(const mytipe t, int test_name)
{
	if (test_name == 1) { return 0.; }  //тест1
	else if (test_name == 2) { return 0.; } //тест2
	else if (test_name == 3) { return 0.; } //тест3
	else if (test_name == 4) { return 0.; } //тест4
	else if (test_name == 13) { return 0.5 + 3 * t; } //тест13
	else if (test_name == 19) { return 0.5; } //тест19
	else { cout << "test_name isn't correct" << endl; }
}

mytipe u_Lt(const mytipe t, int test_name)
{
	if (test_name == 1) { return 0.; }  //тест1
	else if (test_name == 2) { return 0.; } //тест2
	else if (test_name == 3) { return 0.; } //тест3
	else if (test_name == 4) { return 0.; } //тест4
	else if (test_name == 13) { return 1; } //тест13
	else if (test_name == 19) { return 3 - 2 * t; } //тест19
	else { cout << "test_name isn't correct" << endl; }
}

mytipe real_sol(const mytipe x, const mytipe t, int test_name)
{
	if (test_name == 1) { return sin(PI*x)*cos(PI*t); }  //тест1
	else if (test_name == 2) {   //тест2
		const mytipe Eps = 1e-7;  //точность аналитического решения
		const int k = 1 / 2 * (sqrt(2 / pow(PI, 2) / Eps) + 1);
		mytipe sum = 0.;
		for (int n = 0; n <= k; n++) {
			sum += 1 / pow(2 * n + 1, 3)*sin((2 * n + 1)*PI*x)*cos((2 * n + 1)*PI*t);
		}
		sum *= 8 / pow(PI, 3);
		return sum;
	}
	else if (test_name == 3) {   //тест3
		const mytipe Eps = 1e-7;  //точность аналитического решения
		const int k = 1 / 2 * (sqrt(2 / pow(PI, 2) / Eps) + 1);
		mytipe sum = 0.;
		for (int n = 0; n <= k; n++) {
			sum += 1 / pow(2 * n + 1, 3)*sin((2 * n + 1)*PI*x)*cos((2 * n + 1)*PI*t);
		}
		sum *= 8 / pow(PI, 3);
		return sum;
	}
	else if (test_name == 4) { cout << "for test_name " << test_name << " analitic solution doesn't exist" << endl; } //тест4
	else if (test_name == 13) { cout << "for test_name " << test_name << " analitic solution doesn't exist" << endl; } //тест13
	else if (test_name == 19) { cout << "for test_name " << test_name << " analitic solution doesn't exist" << endl; } //тест19
	else { cout << "test_name isn't" << endl; }
}

int main() {
	setlocale(LC_ALL, "Russian"); //подключаем русский язык
	mytipe L = 1.;//длина стержня
	mytipe T = 2.;//Время


	mytipe h = 0.005;  //шаг по пространству
	mytipe c = 1;
	mytipe courant = 1.;
	mytipe tau = courant*h / c; //шаг по времени
	int n = L / h + 1; //число узлов по пространству
	int k = T / tau + 1; //число узлов по времени

	int test_name = 4; //Номер теста

	mytipe** u;

	u = Krest_a(n, k, h, tau, test_name); //разностная схема типа крест

										//vivod(n, u[1]);
										//for (int i = 0; i < n; i++) { cout << sin(i*h)*exp(-tao) << " "; }

	//cout << "погрешность = " << error(k, n, h, tau, u, test_name);

	Copy_to_file_anim(n, k, u, "anim.dat");

	//char file_name1[] = "SOLVE1.txt";
	//Copy_to_file(n, h, u[0], file_name1);

	//char file_name2[] = "SOLVE2.txt";
	//Copy_to_file(n, h, u[1], file_name2);

	//char file_name3[] = "SOLVE3.txt";
	//Copy_to_file(n, h, u[2], file_name3);

	//char file_name4[] = "SOLVE4.txt";
	//Copy_to_file(n, h, u[8], file_name4);

	//char file_name5[] = "SOLVE5.txt";
	//Copy_to_file(n, h, u[10], file_name5);


	cout << "\nend:)";

	cin.get();
	return 0;
}

mytipe** Krest(int n, int k, mytipe h, mytipe tau, int test_name)
{
	const mytipe a = 1;

	mytipe** Y = new mytipe*[k];
	for (int j = 0; j < k; j++) { Y[j] = new mytipe[n]; };

	for (int i = 0; i < n; i++) { Y[0][i] = ux0(i*h, test_name); }

	for (int i = 0; i < n; i++) { Y[1][i] = Y[0][i] + tau*dux0(i*h, test_name) + pow(a*tau, 2) / 2 * ddux0(i*h, test_name); }

	for (int j = 2; j < k; j++) {
		Y[j][0] = u_0t(j*tau, test_name);
		Y[j][n - 1] = u_Lt(j*tau, test_name);
		for (int i = 1; i < n - 1; i++) {
			Y[j][i] = 2 * Y[j - 1][i] - Y[j - 2][i] + pow(a*tau / h, 2)*(Y[j - 1][i + 1] - 2 * Y[j - 1][i] + Y[j - 1][i - 1]);
		}
	}
	return Y;
}

mytipe** Krest_a(int n, int k, mytipe h, mytipe tau, int test_name)
{
	const mytipe a = 1;

	mytipe** Y = new mytipe*[k];
	for (int j = 0; j < k; j++) { Y[j] = new mytipe[n]; };

	for (int i = 0; i < n; i++) { Y[0][i] = ux0(i*h, test_name); }
	Y[1][0] = Y[0][0] + tau*dux0(0*h, test_name) + pow(a*tau, 2) / 2 * ddux0(0*h, test_name);
	for (int i = 1; i < n - 1; i++) { Y[1][i] = Y[0][i] + tau*dux0(i*h, test_name) + pow(a*tau, 2) / 2 * ddux0_a(ux0((i - 1)*h, test_name), ux0((i)*h, test_name), ux0((i + 1)*h, test_name), h, test_name); }
	Y[1][n-1] = Y[0][n-1] + tau*dux0((n-1)*h, test_name) + pow(a*tau, 2) / 2 * ddux0((n-1)*h, test_name);
	for (int j = 2; j < k; j++) {
		Y[j][0] = u_0t(j*tau, test_name);
		Y[j][n - 1] = u_Lt(j*tau, test_name);
		for (int i = 1; i < n - 1; i++) {
			Y[j][i] = 2 * Y[j - 1][i] - Y[j - 2][i] + pow(a*tau / h, 2)*(Y[j - 1][i + 1] - 2 * Y[j - 1][i] + Y[j - 1][i - 1]);
		}
	}
	return Y;
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

mytipe error(int k, int n, mytipe h, mytipe tau, mytipe** a, int test_name)
{
	mytipe MAXerror = 0.0;

	for (int j = 0; j < k; j++)
		for (int i = 0; i < n; i++)
			if (fabs(a[j][i] - real_sol(i*h, j*tau, test_name)) > MAXerror) {
				MAXerror = fabs(a[j][i] - real_sol(i*h, j*tau, test_name));
				cout << j << "  " << i << "  " << MAXerror << endl;
			};

	return MAXerror;
}

void Copy_to_file_anim(const unsigned int n, const unsigned int k, mytipe** solve, char file_name[]) //n - кол-во узлов по пространству, k-кол-во узлов по времени
{
	ofstream fout;    // создали переменную для записи в файл
	fout.open(file_name, ios_base::out | ios_base::trunc);
	for (unsigned int j = 0; j < k; j++) {
		for (unsigned int i = 0; i < n; i++) {
			fout << solve[j][i] << " ";
		};
		fout << endl;
	}
	fout.close();
	fout.clear();
}

void Copy_to_file(const unsigned int n, const mytipe h, mytipe* solve, char file_name[]) //n - кол-во узлов по пространству, h-шаг по пространству 
{
	ofstream fout;    // создали переменную для записи в файл
	fout.open(file_name, ios_base::out | ios_base::trunc);
	for (unsigned int i = 0; i < n; i++) {
		fout << " " << i*h << " " << solve[i] << " " << endl;
	};
	fout.close();
	fout.clear();
}


void copy(int dim, mytipe* a, mytipe* b)
{
	for (int i = 0; i < dim; i++) a[i] = b[i];
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

mytipe* minus_vect(const unsigned int DIM, const mytipe* const b1, const mytipe* const b2) //b1-b2
{
	mytipe* x = new mytipe[DIM];
	for (int i = 0; i < DIM; i++) {
		x[i] = b1[i] - b2[i];
	}
	return x;
}