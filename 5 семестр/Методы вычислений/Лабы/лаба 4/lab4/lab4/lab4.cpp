#include "iostream"
#include <fstream>
#include <typeinfo>
#include <cmath>

using namespace std;

const double EPSILON = 1e-8;   //константа сравнения с 0
const double EPS_eigen = 1e-4;
const double PI = 3.14159265;

typedef double mytipe; //тип данных, использующийся во всей программе

#include "Matr.h"

mytipe* eigenvalues(const int N, mytipe** A);
void Hessenberg_matrix(const unsigned int DIM, mytipe** A);
mytipe** eigenvectors(const int N, mytipe** A, mytipe* lambda);
void Raley(int N, mytipe* lambda, mytipe** e,mytipe** A);
bool QThis_Eigenvalues_A(int N, mytipe* lambda, mytipe** A);
bool QThis_Eigenvectors_A(int N, mytipe* lambda, mytipe** e, mytipe** A);

int main() {
	setlocale(LC_ALL, "Russian");

	ifstream fin;    // создали переменную для считывания из файла
	fin.open("1.txt", ios_base::in | ios_base::app | ios_base::binary);  //открыли файл
	unsigned int a;  // размер пространства
	fin >> a;  //считали из файла размер пространства
	cout << "Dim=" << a << endl;
	const unsigned int DIM = a;  //размер пространнства сделали константным
	mytipe **A;                      // |создали 
	A = new mytipe *[DIM];          // |динамический 
	mytipe **e;                      // |создали 
	e = new mytipe *[DIM];          // |динамический
	for (int i = 0; i < DIM; i++) {   // |массив 
		A[i] = new mytipe[DIM];     // |матрицы А
		e[i] = new mytipe[DIM];
	}

	
   for (int i = 0; i < DIM; i++)  //считываем из файла матрицу и столбец
	{
		for (int j = 0; j < DIM; j++)
		{
			fin >> A[i][j];  //матрицу
		}
	}

	fin.close();  //закрываем файл
	fin.clear();  //отвязываем файл от файловой переменной

	cout << "\nA=";
	vivod(DIM, A);

	mytipe* lambda = new mytipe[DIM];

	cout << "\n\t\tQR-АЛГОРИТМ\n\n";
	lambda = eigenvalues(DIM, A);
	cout << endl;
	cout << "собственные значения матрицы A:"; vivod(DIM, lambda);
	//QThis_Eigenvalues_A(DIM,lambda,A);

	cout << "\n\t\tМЕТОД ОБРАТНЫХ ИТЕРАЦИЙ\n\n";
	cout << "собственные векторы матрицы A:\n\n";
	copy(DIM, e, eigenvectors(DIM, A, lambda));
	cout << endl;
	for (int i = 0; i < DIM; i++) { cout << "e" << i + 1 << " = "; vivod(DIM, e[i]); }
	//QThis_Eigenvectors_A(DIM, lambda, e, A);

	cout << "\n\t\tС ПОМОЩЬЮ СООТНОШЕНИЯ РЭЛЕЯ\n\n";
	Raley(DIM, lambda, e, A);
	cout << "\nсобственные значения матрицы A:";
	vivod(DIM, lambda);
	cout << "\nсобственные векторы матрицы A:\n";
	for (int i = 0; i < DIM; i++) { cout << "e" << i + 1 << " = "; vivod(DIM, e[i]); }

	//QThis_Eigenvalues_A(DIM, lambda, A);
	//QThis_Eigenvectors_A(DIM, lambda, e, A);

	delete[] lambda;
	for (int i = 0; i < DIM; i++) { delete[] A[i]; delete[] e[i];}
	delete[] A; delete[] e;

	cin.get();
	return 0;
}

mytipe* eigenvalues(const int N, mytipe** A)
{
	mytipe* lambda = new mytipe[N];
	mytipe sigma;
	mytipe**E=new mytipe*[N];
	mytipe**Ak = new mytipe*[N];
	mytipe**Qk = new mytipe*[N];
	mytipe**Rk = new mytipe*[N];
	mytipe* vspom = new mytipe [N-1];
	int iter = 0;
	int I = 0;
	for (int i = 0; i < N; i++) { E[i] = new mytipe[N]; Ak[i] = new mytipe[N]; Qk[i] = new mytipe[N]; Rk[i] = new mytipe[N];}; 
	singlematr(N, E);
	copy(N, Ak, A);
	Hessenberg_matrix(N, Ak);
	cout << "H = "; vivod(N, Ak);
	for (int n = N - 1,m=N; n > 0; n--,m--)
	{ 
		iter = 0;
		do
		{
				sigma = Ak[n][n];
				copy(m, Ak, sum_matr(m, Ak, mult_matrnum(m, E, -sigma)));//A-сигма*Е (сдвиг)
				QR(m, Qk, Rk, Ak);
				//copy(m, Ak, multimatrix(m, Rk, Qk));
				copy(m, Ak, sum_matr(m, multimatrix(m, Rk, Qk), mult_matrnum(m, E, sigma)));//RkQk+сигма*Е (обратный сдвиг)
				//if (iter % 5 == 0) { cout << "A="; vivod(N, Ak); }
				if (iter > 200) { cout << "ЗАЦИКЛИВАНИЕ!!!\n"; break; }
			    iter++;
			   //for (int k = 0; k < n; k++) { vspom[k] = Ak[k + 1][k]; };
		} while (fabs(Ak[n][n - 1])>EPS_eigen);
		lambda[n] = Ak[n][n];
		I += iter;
		cout << N-n << ") " << iter << " итераций\n";
	}
	cout << N << ") " << 0 << " итераций\nвсего ";
	cout << I << " итераций\n";
	lambda[0] = Ak[0][0];

	for (int i = 0; i < N; i++)
	{
	delete[] E[i];
	delete[] Ak[i];
	delete[] Qk[i];
	delete[] Rk[i];}
	delete[] E; delete[] Ak; delete[] Qk; delete[] Rk;
	return lambda;
}

void Hessenberg_matrix(const unsigned int DIM, mytipe** A)
{
	mytipe vspom;   //вспомогательный указатель для смены элементов
	mytipe c;  //вычисляемый коэффициент Т_ij
	mytipe s;  //вычисляемый коэффициент Т_ij

	for (int k = 0; k < DIM - 2; k++)  //нахождение матриц T=Q^(-1) и R
									   //высчитывание новых T_ij, перебирая столбцы
	{
		for (int l = k + 2; l < DIM; l++)  //высчитывание новых T_ij, перебирая строки под главной диагональю
		{

			c = A[k + 1][k] / sqrt((A[k + 1][k])*(A[k + 1][k]) + (A[l][k])*(A[l][k]));//вычисление c
			s = A[l][k] / sqrt((A[k + 1][k])*(A[k + 1][k]) + (A[l][k])*(A[l][k]));//вычисление s

			for (int j = k; j < DIM; j++)
			{
				vspom = A[k + 1][j];
				A[k + 1][j] = c * A[k + 1][j] + s * A[l][j]; //умножение двух строк Т
				A[l][j] = c * A[l][j] - s * vspom;    //на столбец А
			}
			for (int j = k; j < DIM; j++)
			{
				vspom = A[j][k + 1];
				A[j][k + 1] = c * A[j][k + 1] + s * A[j][l]; //умножение двух столбцов Т
				A[j][l] = c * A[j][l] - s * vspom;    //на строку А
			}
		}
	}
}

mytipe** eigenvectors(const int N, mytipe** A,mytipe* lambda)
{
	mytipe** e = new mytipe*[N];//матрица собственных векторов
	mytipe** Ak = new mytipe*[N];//Матрица А-лямбдаЕ
	mytipe** E = new mytipe*[N];//
	mytipe* e_old= new mytipe [N];
	int iter = 0;
	int I = 0;
	for (int i = 0; i < N; i++) { e[i] = new mytipe[N]; Ak[i] = new mytipe[N]; E[i] = new mytipe[N]; } singlematr(N, E); singlematr(N, e);
	for (int j = 0; j < N; j++)
	{   
		iter = 0;
	copy(N,Ak,sum_matr(N, A, mult_matrnum(N, E, -lambda[j])));
	do
	{
		copy(N, e_old, e[j]);
		findx(N, e[j], Ak);//решение СЛАУ
		hod(N, e[j], Ak);
		copy(N, e[j], mult_vectnum(N, e[j], 1 / norm(N, e[j], '2')));//нормировка	
		iter++;
	} while (fabs(fabs(skalar(N, e[j], e_old)) - 1)>EPS_eigen);//(norm(N, e_old, '2') > EPSILON); //(norm(N, matrvec(N, Ak, e[j]), '2')>EPS);
	cout << j + 1 << ") " << iter << " итераций\n";
	I += iter;
	}
	cout << "всего итераций:" << I<<endl;
	delete[] e_old;
	for (int i = 0; i < N; i++) { delete[] Ak[i];  delete[] E[i]; };
	delete[] E; delete[] Ak;

	return e;
}

void Raley(int N, mytipe* lambda, mytipe** e,mytipe** A)
{
	mytipe** Ak = new mytipe*[N];//Матрица А-лямбдаЕ
	mytipe** E = new mytipe*[N];//
	mytipe* vspom = new mytipe[N];
	int iter = 0;
	int I = 0;
	for (int i = 0; i < N; i++) { e[i] = new mytipe[N]; Ak[i] = new mytipe[N]; E[i] = new mytipe[N]; } 
	singlematr(N, E); singlematr(N, e);
	for (int j = 0; j < N; j++)
	{
		
		for (int i = 0; i < j; i++) {
			copy(N, vspom, mult_vectnum(N, proj(N, e[i], e[j]), -1));
			copy(N, e[j], sum_vect(N, e[j],vspom));	}//Ортогонализация	Грамма-Шмидта 

		copy(N, e[j], mult_vectnum(N, e[j], 1 / norm(N, e[j], '2')));//нормировка
		iter = 0;
		do
		{
			copy(N, vspom, e[j]);
			lambda[j] = skalar(N, matrvec(N, A, e[j]), e[j]);
			copy(N, Ak, sum_matr(N, A, mult_matrnum(N, E, -lambda[j])));
			findx(N, e[j], Ak);//решение СЛАУ
			hod(N, e[j], Ak);
			copy(N, e[j], mult_vectnum(N, e[j], 1 / norm(N, e[j], '2')));//нормировка	
			iter++;
		} while (fabs(fabs(skalar(N, e[j], vspom)) - 1)>EPS_eigen);//(norm(N, matrvec(N, Ak, e[j]), '2')>EPS);
		cout << j + 1 << ") " << iter << " итераций\n";
		I += iter;
	}

	cout << "всего итераций:" << I << endl;
	delete[] vspom;
	for (int i = 0; i < N; i++) { delete[] Ak[i];  delete[] E[i]; };
	delete[] E; delete[] Ak;
}

bool QThis_Eigenvalues_A(int N, mytipe* lambda, mytipe** A)
{
	mytipe** Ak = new mytipe*[N];//Матрица А-лямбдаЕ
	mytipe** E = new mytipe*[N];//
	mytipe sym = 0.0;
	for (int i = 0; i < N; i++) { Ak[i] = new mytipe[N]; E[i] = new mytipe[N]; };
	singlematr(N, E);
	for (int k = 0; k < N; k++) {
		for (int j = 0; j < N; j++) {
			sym = 0.0;
			copy(N, Ak, sum_matr(N, A, mult_matrnum(N, E, -lambda[k])));
			for (int i = 0; i < j; i++) { sym += fabs(Ak[j][i]); };
			for (int i = j + 1; i < N; i++) { sym += fabs(Ak[j][i]); };
			if (fabs(Ak[j][j]) - sym > EPS_eigen) {
				cout << endl << j + 1 << " значение не является собственным.\n" ;
				for (int i = 0; i < N; i++) { delete[] Ak[i];  delete[] E[i]; };
				delete[] E; delete[] Ak;
				return false;
			};
		};
	};
	for (int i = 0; i < N; i++) { delete[] Ak[i];  delete[] E[i]; };
	delete[] E; delete[] Ak;
	cout << "\nВсе значения являются собственными.\n";
	return true;
}

bool QThis_Eigenvectors_A(int N, mytipe* lambda, mytipe** e, mytipe** A)
{
	mytipe** Ak = new mytipe*[N];//Матрица А-лямбдаЕ
	mytipe** E = new mytipe*[N];//
	for (int i = 0; i < N; i++) { Ak[i] = new mytipe[N]; E[i] = new mytipe[N]; };
	singlematr(N, E);
	for (int j = 0; j < N; j++) {
		copy(N, Ak, sum_matr(N, A, mult_matrnum(N, E, -lambda[j])));
		if(norm(N, matrvec(N, Ak, e[j]), '2')>10*EPS_eigen){
			cout << endl << j + 1 << " вектор не является собственным.\n";
			for (int i = 0; i < N; i++) { delete[] Ak[i];  delete[] E[i]; };
			delete[] E; delete[] Ak;
			return false;
		};
	};

	for (int i = 0; i < N; i++) { delete[] Ak[i];  delete[] E[i]; };
	delete[] E; delete[] Ak;
	cout << "\nВсе вектора являются собственными.\n";
	return true;
}
