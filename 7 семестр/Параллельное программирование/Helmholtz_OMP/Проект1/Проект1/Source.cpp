#include "iostream"
#include <fstream>
#include <cmath>
#include <omp.h>

using namespace std;

typedef double mytipe;  //тип данных, использующийся во всей программе
const mytipe EPSILON = 1e-8;//константа сравнения с 0
const mytipe Pi = 3.14159265358979323846264338328;

mytipe f(mytipe x, mytipe y, mytipe k)
{
	return 2.*sin(Pi*y) + k*k*(1. - x)*x*sin(Pi*y) + Pi*Pi*(1. - x)*x*sin(Pi*y);
}

mytipe u(mytipe x, mytipe y)
{
	return (1. - x)*x*sin(Pi*y);
}

void print_file(int N, mytipe h, mytipe** A)
{
	ofstream fout;    // создали переменную для записи в файл
	fout.open("Solve.dat", ios_base::out | ios_base::trunc);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) fout << A[i][j] << " ";
		fout << endl;
	}
	fout.close();
	fout.clear();
}

int main() {
	setlocale(LC_ALL, "Russian"); //подключаем русский язык

	int N = 33;//Кол-во разбиений
	int n = N + 1;//Кол-во узлов
	mytipe k = 100.;
	mytipe h = 1. / N;
	mytipe eps = 1e-5;//точность

	mytipe** yk = new mytipe*[n];
	mytipe** yk1 = new mytipe*[n];

	for (int i = 0; i < n; i++) {
		yk[i] = new mytipe[n];
		yk1[i] = new mytipe[n];
	}


	for (int i = 0; i < n; i++) {//ГУ
		yk1[0][i] = 0.; yk1[N][i] = 0.;  yk1[i][0] = 0.; yk1[i][N] = 0.;
	}

	mytipe c = 1. / (k*k*h*h + 4.);
	mytipe d = h*h*c;
	mytipe* max = new mytipe[4];
	mytipe vspom = 0.;
	int iter = 0;
	mytipe t = 0.;



	for (int p = 2; p <= 2; p++)
	{
		omp_set_num_threads(p);//кол-во тредов, на которых ведется расчет

		if (p > 1)
			cout << "\n\n\n" << p << " треда\n";
		else
			cout << "1 тред\n";

		iter = 0;
		vspom = 0.;
		t = 0.;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)  yk[i][j] = 0.; //нач приближение

   //МЕТОД ЯКОБИ

		/*t -= omp_get_wtime();

		do {
			for (int i = 0; i < 4; i++)
				max[i] = 0.;

#pragma omp parallel for
			for (int i = 1; i < N; i++)
				for (int j = 1; j < N; j++)
				{
					yk1[i][j] = c*(yk[i + 1][j] + yk[i - 1][j] + yk[i][j + 1] + yk[i][j - 1]) + d*f(i*h, j*h, k);
					if (fabs(yk1[i][j] - yk[i][j]) > max[omp_get_thread_num()]) max[omp_get_thread_num()] = fabs(yk1[i][j] - yk[i][j]);
				}

			swap(yk, yk1);

			iter++;
			vspom = 0.;
			for (int i = 0; i < 4; i++)
				if (vspom < max[i]) vspom = max[i];

		} while(iter < 1500); //(vspom > eps && iter < 1500);

		t += omp_get_wtime();
		cout << "\n МЕТОД ЯКОБИ:\nt = " << t;

		vspom = 0.; for (int i = 0; i < 4; i++) max[i] = 0.;

#pragma omp parallel for
		for (int i = 1; i < N; i++)
			for (int j = 1; j < N; j++)
			{
				if (fabs(yk1[i][j] - u(i*h, j*h)) > max[omp_get_thread_num()]) max[omp_get_thread_num()] = fabs(yk1[i][j] - u(i*h, j*h));
			}

		for (int i = 0; i < 4; i++)
			if (vspom < max[i]) vspom = max[i];

		cout << "\nОшибка = " << vspom << "\nКолличество итераций = " << iter;

		*/

		//МЕТОД ЗЕЙДЕЛЯ

		iter = 0;
		vspom = 0.;
		t = 0.;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)  yk[i][j] = 0.; //нач приближение

		t -= omp_get_wtime();

		do {
			for (int i = 0; i < 4; i++)
				max[i] = 0.;

#pragma omp parallel 
			{
        #pragma omp for
				for (int i = 1; i < N; i++)
					for (int j = (i % 2 + 1); j < N; j += 2)
					{
						yk1[i][j] = c*(yk[i + 1][j] + yk[i - 1][j] + yk[i][j + 1] + yk[i][j - 1]) + d*f(i*h, j*h, k);
						if (fabs(yk1[i][j] - yk[i][j]) > max[omp_get_thread_num()]) max[omp_get_thread_num()] = fabs(yk1[i][j] - yk[i][j]);
					}
        #pragma omp barrier
        #pragma omp for
				for (int i = 1; i < N; i++)
					for (int j = ((i + 1) % 2 + 1); j < N; j += 2)
					{
						yk1[i][j] = c*(yk1[i + 1][j] + yk1[i - 1][j] + yk1[i][j + 1] + yk1[i][j - 1]) + d*f(i*h, j*h, k);
						if (fabs(yk1[i][j] - yk[i][j]) > max[omp_get_thread_num()]) max[omp_get_thread_num()] = fabs(yk1[i][j] - yk[i][j]);
					}
			}
			swap(yk, yk1);

			iter++;

			vspom = 0.;

			for (int i = 0; i < 4; i++)
				if (vspom < max[i]) vspom = max[i];

		} while (iter < 1); //(vspom > eps && iter < 1500);

		t += omp_get_wtime();
		cout << "\n\n МЕТОД ЗЕЙДЕЛЯ:\nt = " << t;

		vspom = 0.; for (int i = 0; i < 4; i++) max[i] = 0.;


#pragma omp parallel for
		for (int i = 1; i < N; i++)
			for (int j = 1; j < N; j++)
			{
				if (fabs(yk1[i][j] - u(i*h, j*h)) > max[omp_get_thread_num()]) max[omp_get_thread_num()] = fabs(yk1[i][j] - u(i*h, j*h));
			}

		for (int i = 0; i < 4; i++)
			if (vspom < max[i]) vspom = max[i];


		cout << "\nОшибка = " << vspom << "\nКолличество итераций = " << iter;

	}

	print_file(n,h,yk);

	for (int i = 0; i < n; i++) {
		delete[] yk[i];
		delete[] yk1[i];
	}

	delete[] max;
	delete[] yk;
	delete[] yk1;

	cout << "\n\nend";
	cin.get();
	return 0;
}