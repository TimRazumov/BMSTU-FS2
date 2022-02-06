#include "iostream"
#include <fstream>
#include <typeinfo>
#include <cmath>

using namespace std;

const double EPSILON = 1e-8;   //константа сравнения с 0
const double PI = 3.14159265;

typedef double mytipe;  //тип данных, использующийся во всей программе

typedef mytipe(*Fanc)(const mytipe, const mytipe);
typedef mytipe(*Func)(const mytipe );

void vivod(const unsigned int DIM, const mytipe* const b);

mytipe* progon3d(int DIM, mytipe* a, mytipe* b, mytipe* c, mytipe* d);


mytipe ux0(const mytipe x,const mytipe y)
{  
	//return 100*cos(x)*cos(y);//свой тест 1
	//return 100*sin(x)*sin(y);//свой тест 2
	//return 100*sin(x)*cos(y);//свой тест 3
	//return 100*cos(x)*sin(y);//свой тест 4

	//для стац случаев (можно поставить что угодно)
	return 1.5;
}


mytipe M_u_G(const mytipe x)
{
	return 0.;
}

mytipe M_f(const mytipe x1, const mytipe x2)
{
	return 0.;
}

mytipe M1_solve(const mytipe x1, const mytipe x2,const mytipe t)
{
	return 100*cos(x1)*cos(x2)*exp(-2*t);
}


mytipe M2_solve(const mytipe x1, const mytipe x2, const mytipe t)
{
	return 100*sin(x1)*sin(x2)*exp(-2*t);
}

mytipe M3_solve(const mytipe x1, const mytipe x2, const mytipe t)
{
	return 100*sin(x1)*cos(x2)*exp(-2*t);
}


mytipe M4_solve(const mytipe x1, const mytipe x2, const mytipe t)
{
	return 100*cos(x1)*sin(x2)*exp(-2*t);
}

mytipe T1_u_G(const mytipe x)
{
	return 1.;
}
mytipe T1_f(const mytipe x1, const mytipe x2)
{
	return 0.;
}
mytipe T1_solve(const mytipe x1, const mytipe x2)
{
	return 1.;
}


mytipe T2_u_l(const mytipe x)
{
	return 1+x;
}
mytipe T2_u_r(const mytipe x)
{
	return 1+x;
}
mytipe T2_du_u(const mytipe x)
{
	return 1.;
}
mytipe T2_du_d(const mytipe x)
{
	return -1.;
}
mytipe T2_f(const mytipe x1, const mytipe x2)
{
	return 0.;
}
mytipe T2_solve(const mytipe x1, const mytipe x2)
{
	return 1+x2;
}


mytipe T3_du_l(const mytipe x)
{
	return 0.;
}
mytipe T3_du_r(const mytipe x)
{
	return 2.;
}
mytipe T3_u_u(const mytipe x)
{
	return 1+x*x;
}
mytipe T3_u_d(const mytipe x)
{
	return x*x;
}
mytipe T3_f(const mytipe x1, const mytipe x2)
{
	return -4.;
}
mytipe T3_solve(const mytipe x1, const mytipe x2)
{
	return(x1*x1 + x2*x2);
}

//нестационарный процесс
void fractional_step(int N1, int N2, int K,//Кол-во узлов
	mytipe h1, mytipe h2, mytipe tao,//Шаг
	Fanc ux0,//НУ 
	bool* cond,//0-первого рода, 1-второго рода
	Func u_left, Func u_right, Func u_down, Func u_up,//ГУ
	Fanc f
);

//стационарный процесс
void fractional_step(int N1, int N2,//Кол-во узлов
	mytipe h1, mytipe h2, mytipe tao,//Шаг
	Fanc ux0,//НУ 
	bool* cond,//0-первого рода, 1-второго рода
	Func u_left, Func u_right, Func u_down, Func u_up,//ГУ
	Fanc f,//правая часть
	mytipe eps// точность решения
);


int main() {
	setlocale(LC_ALL, "Russian"); //подключаем русский язык


	mytipe L = 1.;//длинна стержня
	mytipe T = 1.;//Время

	mytipe eps = 1e-2;


	int N1 = 200; //узлы
	int N2 = 200;
	int K = 400;
	
	mytipe h1 = L / (N1 - 1); 
	mytipe h2 = L / (N1 - 1);
	mytipe tao = T / (K - 1);


	bool* cond = new bool[4];



	//для теста 1
	cond[0] = 0;//левое 
	cond[1] = 0;//правое
	cond[2] = 0;//нижнее
	cond[3] = 0;//верхнее
	fractional_step(N1, N2, h1, h2, tao, ux0, cond, T1_u_G, T1_u_G, T1_u_G, T1_u_G, T1_f,eps);


	//для теста 2
	/*cond[0] = 0;//левое 
	cond[1] = 0;//правое
	cond[2] = 1;//нижнее
	cond[3] = 1;//верхнее
	fractional_step(N1, N2, h1, h2, tao, ux0, cond, T2_u_l, T2_u_r, T2_du_d, T2_du_u, T2_f,eps);*/

	//для теста 3
	/*cond[0] = 1;//левое 
	cond[1] = 1;//правое
	cond[2] = 0;//нижнее
	cond[3] = 0;//верхнее
	fractional_step(N1, N2, h1, h2, tao, ux0, cond, T3_du_l, T3_du_r, T3_u_d, T3_u_u, T3_f,eps);*/

	
	//мой тест 1
	/*cond[0] = 1;//левое
	cond[1] = 1;//правое
	cond[2] = 1;//нижнее
	cond[3] = 1;//верхнее
	fractional_step(N1, N2, K, h1, h2, tao, ux0, cond, M_u_G, M_u_G, M_u_G, M_u_G, M_f);*/

	//мой тест 2
	/*cond[0] = 0;//левое
	cond[1] = 0;//правое
	cond[2] = 0;//нижнее
	cond[3] = 0;//верхнее
	fractional_step(N1, N2, K, h1, h2, tao, ux0, cond, M_u_G, M_u_G, M_u_G, M_u_G, M_f);*/

	//мой тест 3
	/*cond[0] = 0;//левое
	cond[1] = 0;//правое
	cond[2] = 1;//нижнее
	cond[3] = 1;//верхнее
	fractional_step(N1, N2, K, h1, h2, tao, ux0, cond, M_u_G, M_u_G, M_u_G, M_u_G, M_f);*/


	//мой тест 4
	/*cond[0] = 1;//левое
	cond[1] = 1;//правое
	cond[2] = 0;//нижнее
	cond[3] = 0;//верхнее
	fractional_step(N1, N2, K, h1, h2, tao, ux0, cond, M_u_G, M_u_G, M_u_G, M_u_G, M_f);*/

	//mytipe*** u = fractional_step(N1, K, h1, tao, ux0, u_Gt);


	cout << "\nend:)";

	cin.get();
	return 0;
}

void fractional_step(int N1, int N2,//Кол-во узлов
	mytipe h1, mytipe h2, mytipe tao,//Шаг
	Fanc ux0,//НУ 
	bool* cond,//0-первого рода, 1-второго рода
	Func u_left, Func u_right, Func u_down, Func u_up,//ГУ
	Fanc f, mytipe eps
)
{

	ofstream fout;    // создали переменную для записи в файл
	fout.open("Nstat.dat", ios_base::out | ios_base::trunc);

	int iter = 0;
	mytipe error = 0.0;
	mytipe vspom1 = 0.0;

	const mytipe L1 = (N1 - 1)*h1;
	const mytipe L2 = (N2 - 1)*h2;

	mytipe** Y0 = new mytipe*[N1];//слой с прошлой итерации 
	mytipe** Y05 = new mytipe*[N1];//промежуточный слой
	mytipe** Y1 = new mytipe*[N1];//текущий слой

	for (int i = 0; i < N1; i++) {
		 Y0[i] = new mytipe[N2];
		 Y05[i] = new mytipe[N2];
	}; 

	//Коэффициенты для 1-го измерения
	mytipe* diag1_1=new mytipe[N1];
	mytipe* diag2_1 = new mytipe[N1];
	mytipe* diag3_1 = new mytipe[N1];
	mytipe* right1 = new mytipe[N1];

	//Коэффициенты для 2-го измерения
	mytipe* diag1_2 = new mytipe[N2];
	mytipe* diag2_2 = new mytipe[N2];
	mytipe* diag3_2 = new mytipe[N2];
	mytipe* right2 = new mytipe[N2];
	mytipe* vspom;

	diag1_1[0] = 0.;
	diag2_1[0] = (1 - cond[0]) - 2. * cond[0] * (1. / (h1*h1) + 1. / tao);
	diag3_1[0] = cond[0] * 2. / (h1*h1);

	for (int i = 1; i < N1-1; i++) { 
		diag1_1[i] = 1. / (h1*h1); 
		diag3_1[i] = diag1_1[i];
		diag2_1[i] = -2.*(diag1_1[i]+1/tao);
	}

	diag1_1[N1-1] = cond[1] * 2. / (h1*h1);
	diag2_1[N1 - 1] = (1 - cond[1]) - 2. * cond[1] * (1. / (h1*h1) + 1. / tao);
	diag3_1[N1 - 1] = 0.;

	diag1_2[0] = 0.;
	diag2_2[0] = (1 - cond[2]) - 2. * cond[2] * (1. / (h2*h2) + 1. / tao);
	diag3_2[0] = cond[2] * 2. / (h2*h2);

	for (int i = 1; i < N2-1; i++) {
		diag1_2[i] = 1. / (h2*h2);
		diag3_2[i] = diag1_2[i];
		diag2_2[i] = -2. * (diag1_2[i] + 1. / tao);
	}

	diag1_2[N2 - 1] = cond[3] * 2. / (h2*h2);
	diag2_2[N2 - 1] = (1 - cond[3]) - 2. * cond[3] * (1. / (h2*h2) + 1. / tao);
	diag3_2[N2 - 1] = 0.;

	//НУ
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++) Y0[i][j] = ux0(i*h1, j*h2);

	do 
	{ iter++;

		//ПЕРВОЕ ИЗМЕРЕНИЕ

		//j=0
		if (cond[2] == 0)//если нижнее ГУ 1-го рода
			for (int i = 0; i < N1; i++) Y05[i][0] = u_down(i*h1);
		else
		{
			right1[0] = (1 - cond[0])*u_left(0.) - cond[0] * (2 / tao*Y0[0][0] +
				2 / (h2*h2)*(Y0[0][1] - Y0[0][0]) + 2 / h2*u_down(0.) + f(0., 0.) + 2 / h1*u_left(0.));

			for (int i = 1; i < N1 - 1; i++)
				right1[i] = -(2 / tao*Y0[i][0] + 2 / (h2*h2)*(Y0[i][1] - Y0[i][0]) + 2 / h2*u_down(i*h1) + f(i*h1, 0.));

			right1[N1 - 1] = (1 - cond[1])*u_right(0.) - cond[1] * (2 / tao*Y0[N1 - 1][0] +
				2 / (h2*h2)*(Y0[N1 - 1][1] - Y0[N1 - 1][0]) + 2 / h2*u_down(L1) + f(L1, 0.) + 2 / h1*u_right(0.));


			vspom = progon3d(N1, diag1_1, diag2_1, diag3_1, right1);

			for (int i = 0; i < N1; i++)  Y05[i][0] = vspom[i];

			delete[] vspom;

		}


		for (int j = 1; j < N2 - 1; j++)
		{
			right1[0] = (1 - cond[0])*u_left(j*h2) - cond[0] * (2. / tao*Y0[0][j] + (Y0[0][j + 1] - 2. * Y0[0][j] + Y0[0][j - 1]) / (h2*h2) 
				+ f(0, j*h2) + 2. / h1*u_left(j*h2));

			for (int i = 1; i < N1 - 1; i++)
				right1[i] = -(2./tao *Y0[i][j] + (Y0[i][j + 1] - 2. * Y0[i][j] + Y0[i][j - 1]) / (h2*h2) + f(i*h1, j*h2));

			right1[N1 - 1] = (1 - cond[1])*u_right(j*h2) - cond[1] * (2. / tao*Y0[N1 - 1][j]
				+ (Y0[N1 - 1][j + 1] - 2. * Y0[N1 - 1][j] + Y0[N1 - 1][j - 1]) / (h2*h2) + f(L1, j*h2) + 2. / h1*u_right(j*h2));


			vspom = progon3d(N1, diag1_1, diag2_1, diag3_1, right1);

			for (int i = 0; i < N1; i++)  Y05[i][j] = vspom[i];

			delete[] vspom;
		}

		//j=N2
		if (cond[3] == 0)//если верхнее ГУ 1-го рода
			for (int i = 0; i < N1; i++) Y05[i][N2 - 1] = u_up(i*h1);
		else
		{
			right1[0] = (1 - cond[0])*u_left(L2) - cond[0] * (2 / tao*Y0[0][N2 - 1] +
				2 / (h2*h2)*(Y0[0][N2 - 2] - Y0[0][N2 - 1]) + 2 / h2*u_up(0.) + f(0., L2) + 2 / h1*u_left(L2));

			for (int i = 1; i < N1 - 1; i++)
				right1[i] = -(2 / tao*Y0[i][N2 - 1] + 2 / (h2*h2)*(Y0[i][N2 - 2] - Y0[i][N2 - 1]) + 2 / h2*u_up(i*h1) + f(i*h1, L2));

			right1[N1 - 1] = (1 - cond[1])*u_right(L2) - cond[1] * (2 / tao*Y0[N1 - 1][N2 - 1] +
				2 / (h2*h2)*(Y0[N1 - 1][N2 - 2] - Y0[N1 - 1][N2 - 1]) + 2 / h2*u_up(L1) + f(L1, L2) + 2 / h1*u_right(L2));


			vspom = progon3d(N1, diag1_1, diag2_1, diag3_1, right1);

			for (int i = 0; i < N1; i++)  Y05[i][N2 - 1] = vspom[i];

			delete[] vspom;
		}

		//ВТОРОЕ ИЗМЕРЕНИЕ

		//i=0
		if (cond[0] == 0)//если левое ГУ 1-го рода
		{
			Y1[0] = new mytipe[N2];
			for (int j = 0; j < N2; j++) Y1[0][j] = u_left(j*h2);
		}
		else
		{
			right2[0] = (1 - cond[2])*u_down(0.) - cond[2] * (2 / tao*Y05[0][0] +
				2 / (h1*h1)*(Y05[1][0] - Y05[0][0]) + 2 / h2*u_down(0.) + f(0., 0.) + 2 / h1*u_left(0.));

			for (int j = 1; j < N2 - 1; j++)
				right2[j] = -(2 / tao*Y05[0][j] + 2 / (h1*h1)*(Y05[1][j] - Y05[0][j]) + 2 / h1*u_left(j*h2) + f(0, j*h2));

			right2[N2 - 1] = (1 - cond[3])*u_up(0.) - cond[3] * (2 / tao*Y05[0][N2 - 1] +
				2 / (h1*h1)*(Y05[1][N2 - 1] - Y05[0][N2 - 1]) + 2 / h2*u_up(0.) + f(0,L2) + 2 / h1*u_left(L2));


			Y1[0] = progon3d(N2, diag1_2, diag2_2, diag3_2, right2);

		}


		for (int i = 1; i < N1 - 1; i++)
		{
			right2[0] = (1 - cond[2])*u_down(i*h1) - cond[2] * (2 / tao*Y05[i][0] +
				(Y05[i + 1][0] - 2 * Y05[i][0] + Y05[i - 1][0]) / (h1*h1) + f(i*h1, 0.) + 2 / h2*u_down(i*h1));

			for (int j = 1; j < N2 - 1; j++)
				right2[j] = -(2 / tao*Y05[i][j] + (Y05[i + 1][j] - 2 * Y05[i][j] + Y05[i - 1][j]) / (h1*h1) + f(i*h1, j*h2));

			right2[N2 - 1] = (1 - cond[3])*u_up(i*h1) - cond[3] * (2. / tao*Y05[i][N2 - 1]
				+ (Y05[i + 1][N2 - 1] - 2. * Y05[i][N2 - 1] + Y05[i - 1][N2 - 1]) / (h1*h1) + f(i*h1, L2) + 2. / h2*u_up(i*h1));

			Y1[i] = progon3d(N2, diag1_2, diag2_2, diag3_2, right2);
		}

		//i=N1
		if (cond[1] == 0)//если правое ГУ 1-го рода
		{
			Y1[N1-1] = new mytipe[N2];
			for (int j = 0; j < N2; j++) Y1[N1 - 1][j] = u_right(j*h2);
	
		}
		else
		{
			right2[0] = (1 - cond[2])*u_down(L1) - cond[2] * (2. / tao*Y05[N1-1][0] +
				2. / (h1*h1)*(Y05[N1-2][0] - Y05[N1-1][0]) + 2. / h2*u_down(L1) + f(L1, 0.) + 2. / h1*u_right(0.));

			for (int j = 1; j < N2 - 1; j++)
				right2[j] = -(2. / tao*Y05[N1-1][j] + 2. / (h1*h1)*(Y05[N1-2][j] - Y05[N1-1][j]) + 2. / h1*u_right(j*h2) + f(L1, j*h2));

			right2[N2 - 1] = (1 - cond[3])*u_up(L1) - cond[3] * (2. / tao*Y05[N1 - 1][N2-1] +
				2. / (h1*h1)*(Y05[N1-2][N2-1] - Y05[N1 - 1][N2-1]) + 2. / h2*u_up(L1) + f(L1,L2) + 2. / h1*u_right(L2));


			Y1[N1-1] = progon3d(N2, diag1_2, diag2_2, diag3_2, right2);

		}

		for (int i = 0; i < N1; i++) {

			for (int j = 0; j < N2; j++) Y0[i][j] = Y1[i][j];

			delete[] Y1[i];
		}

		
		error = 0.;
			for (int i = 0; i < N1; i++)
			for (int j = 0; j < N2; j++) {
			vspom1 = fabs(Y0[i][j] - Y05[i][j]);
			if (vspom1 > error) error = vspom1;
			};

	}while (error/tao>eps);
	

	for (int i = 0; i < N1; i++) {//запись в файл
		for (int j = 0; j < N2; j++)  fout << Y0[i][j] << " ";
		fout << endl;
	}

	cout  << iter << " итераций\n";

	delete[] right1;
	delete[] diag1_1;
	delete[] diag2_1;
	delete[] diag3_1;
	
	delete[] right2;
	delete[] diag1_2;
	delete[] diag2_2;
	delete[] diag3_2;

	for (int i = 0; i < N1; i++) {delete[] Y0[i]; delete[] Y05[i];}
	delete[] Y0; delete[] Y1; delete[] Y05;


	fout.close();
	fout.clear();
	
}




void fractional_step(int N1, int N2, int K,//Кол-во узлов
	mytipe h1, mytipe h2, mytipe tao,//Шаг
	Fanc ux0,//НУ 
	bool* cond,//0-первого рода, 1-второго рода
	Func u_left, Func u_right, Func u_down, Func u_up,//ГУ
	Fanc f
)
{

	ofstream fout;    // создали переменную для записи в файл
	fout.open("Nsolve.dat", ios_base::out | ios_base::trunc);

	mytipe error = 0.0;
	mytipe vspom1 = 0.0;
	
	const mytipe L1 = (N1 - 1)*h1;
	const mytipe L2 = (N2 - 1)*h2;

	mytipe** Y0 = new mytipe*[N1];//слой с прошлой итерации 
	mytipe** Y05 = new mytipe*[N1];//промежуточный слой
	mytipe** Y1 = new mytipe*[N1];//текущий слой

	for (int i = 0; i < N1; i++) {
		Y0[i] = new mytipe[N2];
		Y05[i] = new mytipe[N2];
	};

	//Коэффициенты для 1-го измерения
	mytipe* diag1_1 = new mytipe[N1];
	mytipe* diag2_1 = new mytipe[N1];
	mytipe* diag3_1 = new mytipe[N1];
	mytipe* right1 = new mytipe[N1];

	//Коэффициенты для 2-го измерения
	mytipe* diag1_2 = new mytipe[N2];
	mytipe* diag2_2 = new mytipe[N2];
	mytipe* diag3_2 = new mytipe[N2];
	mytipe* right2 = new mytipe[N2];
	mytipe* vspom;

	diag1_1[0] = 0.;
	diag2_1[0] = (1 - cond[0]) - 2. * cond[0] * (1. / (h1*h1) + 1. / tao);
	diag3_1[0] = cond[0] * 2. / (h1*h1);

	for (int i = 1; i < N1 - 1; i++) {
		diag1_1[i] = 1. / (h1*h1);
		diag3_1[i] = diag1_1[i];
		diag2_1[i] = -2.*(diag1_1[i] + 1 / tao);
	}

	diag1_1[N1 - 1] = cond[1] * 2. / (h1*h1);
	diag2_1[N1 - 1] = (1 - cond[1]) - 2. * cond[1] * (1. / (h1*h1) + 1. / tao);
	diag3_1[N1 - 1] = 0.;

	diag1_2[0] = 0.;
	diag2_2[0] = (1 - cond[2]) - 2. * cond[2] * (1. / (h2*h2) + 1. / tao);
	diag3_2[0] = cond[2] * 2. / (h2*h2);

	for (int i = 1; i < N2 - 1; i++) {
		diag1_2[i] = 1. / (h2*h2);
		diag3_2[i] = diag1_2[i];
		diag2_2[i] = -2. * (diag1_2[i] + 1. / tao);
	}

	diag1_2[N2 - 1] = cond[3] * 2. / (h2*h2);
	diag2_2[N2 - 1] = (1 - cond[3]) - 2. * cond[3] * (1. / (h2*h2) + 1. / tao);
	diag3_2[N2 - 1] = 0.;

	//НУ
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++) Y0[i][j] = ux0(i*h1, j*h2);

	for (int i = 0; i < N1; i++) {//запись в файл
		for (int j = 0; j < N2; j++)  fout << Y0[i][j] << " ";
		fout << endl;
	}

   for (int k = 1; k < K; k++)
	{
		//ПЕРВОЕ ИЗМЕРЕНИЕ

		//j=0
		if (cond[2] == 0)//если нижнее ГУ 1-го рода
			for (int i = 0; i < N1; i++) Y05[i][0] = u_down(i*h1);
		else
		{
			right1[0] = (1 - cond[0])*u_left(0.) - cond[0] * (2 / tao*Y0[0][0] +
				2 / (h2*h2)*(Y0[0][1] - Y0[0][0]) + 2 / h2*u_down(0.) + f(0., 0.) + 2 / h1*u_left(0.));

			for (int i = 1; i < N1 - 1; i++)
				right1[i] = -(2 / tao*Y0[i][0] + 2 / (h2*h2)*(Y0[i][1] - Y0[i][0]) + 2 / h2*u_down(i*h1) + f(i*h1, 0.));

			right1[N1 - 1] = (1 - cond[1])*u_right(0.) - cond[1] * (2 / tao*Y0[N1 - 1][0] +
				2 / (h2*h2)*(Y0[N1 - 1][1] - Y0[N1 - 1][0]) + 2 / h2*u_down(L1) + f(L1, 0.) + 2 / h1*u_right(0.));


			vspom = progon3d(N1, diag1_1, diag2_1, diag3_1, right1);

			for (int i = 0; i < N1; i++)  Y05[i][0] = vspom[i];

			delete[] vspom;

		}


		for (int j = 1; j < N2 - 1; j++)
		{
			right1[0] = (1 - cond[0])*u_left(j*h2) - cond[0] * (2. / tao*Y0[0][j] + (Y0[0][j + 1] - 2. * Y0[0][j] + Y0[0][j - 1]) / (h2*h2)
				+ f(0, j*h2) + 2. / h1*u_left(j*h2));

			for (int i = 1; i < N1 - 1; i++)
				right1[i] = -(2. / tao *Y0[i][j] + (Y0[i][j + 1] - 2. * Y0[i][j] + Y0[i][j - 1]) / (h2*h2) + f(i*h1, j*h2));

			right1[N1 - 1] = (1 - cond[1])*u_right(j*h2) - cond[1] * (2. / tao*Y0[N1 - 1][j]
				+ (Y0[N1 - 1][j + 1] - 2. * Y0[N1 - 1][j] + Y0[N1 - 1][j - 1]) / (h2*h2) + f(L1, j*h2) + 2. / h1*u_right(j*h2));


			vspom = progon3d(N1, diag1_1, diag2_1, diag3_1, right1);

			for (int i = 0; i < N1; i++)  Y05[i][j] = vspom[i];

			delete[] vspom;
		}

		//j=N2
		if (cond[3] == 0)//если верхнее ГУ 1-го рода
			for (int i = 0; i < N1; i++) Y05[i][N2 - 1] = u_up(i*h1);
		else
		{
			right1[0] = (1 - cond[0])*u_left(L2) - cond[0] * (2 / tao*Y0[0][N2 - 1] +
				2 / (h2*h2)*(Y0[0][N2 - 2] - Y0[0][N2 - 1]) + 2 / h2*u_up(0.) + f(0., L2) + 2 / h1*u_left(L2));

			for (int i = 1; i < N1 - 1; i++)
				right1[i] = -(2 / tao*Y0[i][N2 - 1] + 2 / (h2*h2)*(Y0[i][N2 - 2] - Y0[i][N2 - 1]) + 2 / h2*u_up(i*h1) + f(i*h1, L2));

			right1[N1 - 1] = (1 - cond[1])*u_right(L2) - cond[1] * (2 / tao*Y0[N1 - 1][N2 - 1] +
				2 / (h2*h2)*(Y0[N1 - 1][N2 - 2] - Y0[N1 - 1][N2 - 1]) + 2 / h2*u_up(L1) + f(L1, L2) + 2 / h1*u_right(L2));


			vspom = progon3d(N1, diag1_1, diag2_1, diag3_1, right1);

			for (int i = 0; i < N1; i++)  Y05[i][N2 - 1] = vspom[i];

			delete[] vspom;
		}

		//ВТОРОЕ ИЗМЕРЕНИЕ

		//i=0
		if (cond[0] == 0)//если левое ГУ 1-го рода
		{
			Y1[0] = new mytipe[N2];
			for (int j = 0; j < N2; j++) Y1[0][j] = u_left(j*h2);
		}
		else
		{
			right2[0] = (1 - cond[2])*u_down(0.) - cond[2] * (2 / tao*Y05[0][0] +
				2 / (h1*h1)*(Y05[1][0] - Y05[0][0]) + 2 / h2*u_down(0.) + f(0., 0.) + 2 / h1*u_left(0.));

			for (int j = 1; j < N2 - 1; j++)
				right2[j] = -(2 / tao*Y05[0][j] + 2 / (h1*h1)*(Y05[1][j] - Y05[0][j]) + 2 / h1*u_left(j*h2) + f(0, j*h2));

			right2[N2 - 1] = (1 - cond[3])*u_up(0.) - cond[3] * (2 / tao*Y05[0][N2 - 1] +
				2 / (h1*h1)*(Y05[1][N2 - 1] - Y05[0][N2 - 1]) + 2 / h2*u_up(0.) + f(0, L2) + 2 / h1*u_left(L2));


			Y1[0] = progon3d(N2, diag1_2, diag2_2, diag3_2, right2);

		}


		for (int i = 1; i < N1 - 1; i++)
		{
			right2[0] = (1 - cond[2])*u_down(i*h1) - cond[2] * (2 / tao*Y05[i][0] +
				(Y05[i + 1][0] - 2 * Y05[i][0] + Y05[i - 1][0]) / (h1*h1) + f(i*h1, 0.) + 2 / h2*u_down(i*h1));

			for (int j = 1; j < N2 - 1; j++)
				right2[j] = -(2 / tao*Y05[i][j] + (Y05[i + 1][j] - 2 * Y05[i][j] + Y05[i - 1][j]) / (h1*h1) + f(i*h1, j*h2));

			right2[N2 - 1] = (1 - cond[3])*u_up(i*h1) - cond[3] * (2. / tao*Y05[i][N2 - 1]
				+ (Y05[i + 1][N2 - 1] - 2. * Y05[i][N2 - 1] + Y05[i - 1][N2 - 1]) / (h1*h1) + f(i*h1, L2) + 2. / h2*u_up(i*h1));

			Y1[i] = progon3d(N2, diag1_2, diag2_2, diag3_2, right2);
		}

		//i=N1
		if (cond[1] == 0)//если правое ГУ 1-го рода
		{
			Y1[N1 - 1] = new mytipe[N2];
			for (int j = 0; j < N2; j++) Y1[N1 - 1][j] = u_right(j*h2);

		}
		else
		{
			right2[0] = (1 - cond[2])*u_down(L1) - cond[2] * (2. / tao*Y05[N1 - 1][0] +
				2. / (h1*h1)*(Y05[N1 - 2][0] - Y05[N1 - 1][0]) + 2. / h2*u_down(L1) + f(L1, 0.) + 2. / h1*u_right(0.));

			for (int j = 1; j < N2 - 1; j++)
				right2[j] = -(2. / tao*Y05[N1 - 1][j] + 2. / (h1*h1)*(Y05[N1 - 2][j] - Y05[N1 - 1][j]) + 2. / h1*u_right(j*h2) + f(L1, j*h2));

			right2[N2 - 1] = (1 - cond[3])*u_up(L1) - cond[3] * (2. / tao*Y05[N1 - 1][N2 - 1] +
				2. / (h1*h1)*(Y05[N1 - 2][N2 - 1] - Y05[N1 - 1][N2 - 1]) + 2. / h2*u_up(L1) + f(L1, L2) + 2. / h1*u_right(L2));


			Y1[N1 - 1] = progon3d(N2, diag1_2, diag2_2, diag3_2, right2);

		}

		for (int i = 0; i < N1; i++) {
			for (int j = 0; j < N2; j++)  Y0[i][j] = Y1[i][j];
				
			delete[] Y1[i];
		}

		for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++) {
		vspom1 = fabs(Y0[i][j] - M4_solve(i*h1, j*h2,k*tao));
		if (vspom1 > error) error = vspom1;
		};

		
		if (k % 4 == 0) 
			for (int i = 0; i < N1; i ++) {//запись в файл
				for (int j = 0; j < N2; j ++)  fout << Y0[i][j] << " ";
				fout << endl;
			}
		
	}



	cout << "Ошибка " << error;

	delete[] right1;
	delete[] diag1_1;
	delete[] diag2_1;
	delete[] diag3_1;

	delete[] right2;
	delete[] diag1_2;
	delete[] diag2_2;
	delete[] diag3_2;

	for (int i = 0; i < N1; i++) { delete[] Y0[i]; delete[] Y05[i]; }
	delete[] Y0; delete[] Y1; delete[] Y05;


	fout.close();
	fout.clear();

}




mytipe* progon3d(int DIM, mytipe* a, mytipe* b, mytipe* c, mytipe* d)//ax1+bx2+cx3=f
{
	mytipe* alfa = new mytipe[DIM];
	mytipe* betta = new mytipe[DIM];
	mytipe znam;

	alfa[0] = -c[0] / b[0]; betta[0] = d[0] / b[0];

	for (int i = 1; i < DIM; i++) {
		znam = b[i] + a[i] * alfa[i - 1];
		alfa[i] = -c[i] / znam;
		betta[i] = (d[i]-a[i] * betta[i - 1]) / znam;
	}

	for (int i = DIM - 2; i > -1; i--) betta[i] += alfa[i] * betta[i + 1];
	
	delete[] alfa;

	return  betta;
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


/*mytipe*** fractional_step(int n, int K,
	mytipe h, mytipe tao,
	Fanc ux0, Fanc u_Gt)
{
	const int N = n - 1;//колличество разбиений

	mytipe*** Y = new mytipe**[K];

	for (int i = 0; i < K; i++)
	{
		Y[i] = new mytipe*[n];
		for (int j = 0; j < n; j++) Y[i][j] = new mytipe[n];
	}

	for (int i = 0; i < n; i++) 
		for (int j = 0; j < n; j++) Y[0][i][j] = ux0(i*h,j*h); //нач условие

	mytipe* diag1 = new mytipe[N - 1];
	mytipe* diag2 = new mytipe[N - 1];
	mytipe* diag3 = new mytipe[N - 1];
	mytipe* right = new mytipe[N - 1];

	for (int i = 0; i < N - 1; i++) { diag1[i] = 1 / (h*h); diag3[i] = diag1[i]; diag2[i] = -2. * (diag1[i]+1/tao); }
	mytipe d1 = diag1[0];
	mytipe d3 = diag3[N-2];
	diag1[0] = 0.;
	diag3[N-2] = 0.;

	
	mytipe* vspom = new mytipe[N-1];
	mytipe** Y_up = new mytipe*[N - 1];
	mytipe error = 0.0;
	mytipe vspom1;
	
	for (int k = 1;  k < K; k++)
	{
		for (int i = 0; i < n; i++){ 
			Y[k][i][0] = u_Gt(tao*(k - 0.5),i*h); Y[k][i][N] = u_Gt(tao*(k - 0.5), i*h);
			Y[k][0][i] = u_Gt(tao*(k - 1), i*h); Y[k][N][i] = u_Gt(tao*(k - 1), i*h);
		}//граничные условия
			

		for (int j = 1; j <= N - 1; j++) { 
			for (int i = 1; i < N; i++) {
				right[i-1] = -(2. / tao*Y[k - 1][i][j] + (Y[k - 1][i][j + 1] - 2.*Y[k - 1][i][j] + Y[k - 1][i][j - 1]) / (h*h));
			}
			right[0] -= d1*Y[k][0][j]; right[N-2] -= d3*Y[k][N][j];
			
			vspom=progon3d(N-1,diag1, diag2, diag3, right); 
			
			for (int i = 1; i < N; i++) Y[k][i][j] = vspom[i - 1];

			//for (int i = 0; i < n; i++) cout<<Y[k][i][j]<<" ";

			delete[] vspom;
		}


		for (int i = 1; i <= N - 1; i++) {

			for (int j = 1; j < N; j++) right[j - 1] = -(2. / tao*Y[k][i][j] + (Y[k][i + 1][j] - 2.*Y[k][i][j] + Y[k][i - 1][j]) / (h*h));
		    right[0] -= d1*Y[k][i][0]; right[N - 2] -= d3*Y[k][i][N];
			Y_up[i-1] = progon3d(N - 1, diag1, diag2, diag3, right);
		
		}


		for (int i = 1; i < N; i++) {
			for (int j = 1; j < N; j++)  Y[k][i][j] = Y_up[i-1][j-1]; 
			delete[] Y_up[i-1];
		}

		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++) {
				vspom1 = fabs(Y[k][i][j]-100*sin(i*h)*sin(j*h)*exp(-2.*k*tao));
				if (vspom1 > error) error = vspom1;
			};
	}
	cout <<"Ошибка "<< error;


	delete[] right;
	delete[] diag1;
	delete[] diag2;
	delete[] diag3;
	delete[] Y_up;

	return Y;
	}*/