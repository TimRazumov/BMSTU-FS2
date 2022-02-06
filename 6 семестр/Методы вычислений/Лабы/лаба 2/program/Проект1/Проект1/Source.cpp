#include "iostream"
#include <fstream>
#include <typeinfo>
#include <cmath>


using namespace std;

const double EPSILON = 1e-8;   //константа сравнения с 0
const double PI = 3.14159265;

typedef double mytipe;  //тип данных, использующийся во всей программе

typedef mytipe(*func)(const mytipe x, const mytipe u0, const mytipe L);
typedef mytipe(*Fanc)(const mytipe x,const mytipe t);
typedef mytipe(*Func)(const mytipe x);


void Copy_to_file(const unsigned int n, const mytipe h, mytipe* solve, char file_name[]); //n - кол-во узлов по пространству, h-шаг по пространству 
void copy(int dim, mytipe* a, mytipe* b);
void vivod(const unsigned int DIM, const mytipe* const b);
mytipe norm(const unsigned int DIM, const mytipe* const b, const char flag); //3 нормы вектора с запросом варианта нормы
mytipe* minus_vect(const unsigned int DIM, const mytipe* const b1, const mytipe* const b2); //b1-b2

mytipe** Integro_interpolation(int n, int k, mytipe h, mytipe tao, Fanc ux0, Func u_0t, mytipe c, mytipe p, Func Kx, Func Pt);//слева равномерно, справа тепловой поток
mytipe** Integro_interpolation(int n, int k, mytipe h, mytipe tao, Fanc ux0, Func u_0t, Fanc u_Lt, mytipe c, mytipe p, Func Kx);//равномерно прогрето с обоих концов
mytipe** Integro_interpolation(int n, int k, mytipe h, mytipe tao, Fanc ux0, mytipe c, mytipe p, Func Kx, Func Pt1, Func Pt2);
mytipe** Kvazilin_explicit(int n, int k, mytipe h, mytipe tao, func ux0, func Ku); //n - кол-во узлов по пространству, k - по времени. 
mytipe** Kvazilin_implicit(const unsigned int n, const unsigned int k, mytipe h, mytipe tao, func ux0, func Ku, const unsigned int M); //n - кол-во узлов по пространству, k - по времени, M-число итерация в неявном цикле
mytipe** Kvazilin_implicit_Eps(const unsigned int n, const unsigned int k, mytipe h, mytipe tao, func ux0, func Ku, const mytipe Eps); //n - кол-во узлов по пространству, k - по времени. 
mytipe** Kvazilin_explicit2(int n, int k, mytipe h, mytipe tao, func ux0, func Ku);//n - кол-во узлов по пространству, k - по времени. 
mytipe* progon3d(int DIM, mytipe* a, mytipe* b, mytipe* c, mytipe* d);

mytipe error(int k, int n, mytipe h, mytipe tao, mytipe** a);


mytipe ux0(const mytipe x, const mytipe L)
{   //тест 2
	mytipe u0 = 0.04;
	//if (L-x<EPSILON) { return 0.; }
	//else if (x < EPSILON) { return 0.; }
	//else {return  u0 + x*(L - x); };

	//return u0;//равномерно прогрет
	return  u0 + x*(L - x);//тест 1

	//return sin(x);//свой тест
}

mytipe ux0_0(const mytipe x, const mytipe u0, const mytipe L)
{
	return 0.;
}

mytipe Kx(const mytipe x)
{
	//методичка 1,2 тест
	mytipe k1 = 0.1;
	mytipe k2 = 1.5;
	mytipe x1 = 0.25;
	mytipe x2 = 0.5;

	//if (x  <= x1) { return k1; }
	//else if (x < x2) { return k1*(x - x2) / (x1 - x2) + k2*(x - x1) / (x2 - x1); }
	//else return k2;

	return 1.;//свой тест

}

mytipe Ku(const mytipe u, const mytipe sigma, const mytipe kappa)
{
	return kappa*pow(u, sigma);
}

mytipe Pt(const mytipe t)//(тепловой поток)
{  // 13 вар 
	mytipe t0 = 0.5;
	mytipe Q = 10.;

	//if (t  < t0 && t>EPSILON) { return Q; }
	//else return 0.;

	return 0.;

}

mytipe u_0t(const mytipe t)
{	
	//13 вариант
	//mytipe u0 = 0.04;
	//return u0;
	
	return 0.;//свой тест
}

mytipe u_Lt(const mytipe L,const mytipe t)
{
	//13 вариант
	//mytipe u0 = 0.04;
	//return u0;
	
	return sin(L)*exp(-t);//свой тест
}

int main() {
	setlocale(LC_ALL, "Russian"); //подключаем русский язык
	mytipe L = 10.;//длинна стержня
	mytipe T = 3.02;//Время

	//13 вариант
	//mytipe c = 2.;
	//mytipe p = 0.5;

	//свой тест
	mytipe c = 1.;
	mytipe p = 1.;

	mytipe h = 0.2;
	mytipe tao = 0.0002;
	int n = L / h+1;
	int k = T / tao+1;


	mytipe** u = new mytipe*[k];
	for (int i = 0; i < k; i++) { u[i] = new mytipe[n]; };


	//u = Integro_interpolation(n, k, h, tao, ux0, u_0t, c, p, Kx, Pt);
	//u = Integro_interpolation(n, k, h, tao, ux0, u_0t,u_Lt, c, p, Kx);
	//u = Integro_interpolation(n, k, h, tao, ux0, c, p, Kx, Pt,Pt);
	//u = Kvazilin_explicit(n, k, h, tao, ux0_0, Ku);
	u = Kvazilin_explicit2(n, k, h, tao, ux0_0, Ku);
	//u = Kvazilin_implicit(n, k, h, tao, ux0_0, Ku, 3);
	//u = Kvazilin_implicit_Eps(n, k, h, tao, ux0_0, Ku, 1e-4);
	
	//vivod(n,u[1]);
	//for (int i = 0; i < n; i++) { cout<< sin(i*h)*exp(-tao)<<" "; }

	//cout <<"погрешность = " <<error(k, n, h, tao, u);

	char file_name1[] = "SOLVE1.txt";
	Copy_to_file(n, h, u[1000], file_name1);

	char file_name2[] = "SOLVE2.txt";
	Copy_to_file(n, h, u[10000], file_name2);

	char file_name3[] = "SOLVE3.txt";
	Copy_to_file(n, h, u[15000], file_name3);

	//char file_name4[] = "SOLVE4.txt";
	//Copy_to_file(n, h, u[8], file_name4);

	//char file_name5[] = "SOLVE5.txt";
	//Copy_to_file(n, h, u[10], file_name5);

	cout << "\nend:)";

	cin.get();
	return 0;
}


mytipe** Kvazilin_explicit(int n, int k, mytipe h, mytipe tao, func ux0, func Ku) //n - кол-во узлов по пространству, k - по времени. 
{
	const int N = n - 1;//колличество разбиений

	mytipe** Y = new mytipe*[k];
	Y[0] = new mytipe[n];

	for (int i = 0; i < n; i++) { Y[0][i] = ux0(i*h, 0, n*h); }

	mytipe sigma = 2;
	mytipe kappa = 0.5;
	mytipe c = 1;
	mytipe c1 = 5;
	mytipe u0 = sigma*pow(c1, 2) / kappa;

	mytipe u0t;
	mytipe ult;


	mytipe* a = new mytipe[2];

	mytipe* diag1 = new mytipe[N - 1];
	mytipe* diag2 = new mytipe[N - 1];
	mytipe* diag3 = new mytipe[N - 1];
	mytipe* right = new mytipe[N - 1];
	mytipe* vspom = new mytipe[N - 1];
	mytipe A0;
	mytipe BN;


	for (int j = 1; j < k; j++) {

		u0t = pow(u0*j*tao, 1 / sigma);
		if (n*h < c1*j*tao) { ult = pow(u0 / c1*(c1*j*tao - n*h), 1 / sigma); }
		else { ult = 0; };

		Y[j] = new mytipe[n];
		a[1] = 0.5*(Ku(Y[j - 1][1], sigma, kappa) + Ku(Y[j - 1][0], sigma, kappa));
		for (int i = 0; i < N - 1; i++) {
			a[0] = a[1];
			a[1] = 0.5*(Ku(Y[j - 1][i + 1], sigma, kappa) + Ku(Y[j - 1][i], sigma, kappa));
			diag1[i] = -a[0] / pow(h, 2);
			diag3[i] = -a[1] / pow(h, 2);
			diag2[i] = -(diag1[i] + diag3[i] - c / tao);
			right[i] = c / tao*Y[j - 1][i + 1];
		};
		A0 = diag1[0]; BN = diag3[N - 2];
		diag1[0] = 0.; diag3[N - 2] = 0.;

		right[0] -= A0*u0t; right[N - 2] -= BN*ult;

		copy(N - 1, vspom, progon3d(N - 1, diag1, diag2, diag3, right));

		for (int i = 1; i < N; i++) Y[j][i] = vspom[i - 1];

		Y[j][0] = u0t;
		Y[j][N] = ult;
	}

	delete[] diag1;
	delete[] diag2;
	delete[] diag3;
	delete[] right;

	delete[] a;
	delete[] vspom;

	return Y;
}

mytipe** Kvazilin_implicit(const unsigned int n, const unsigned int k, mytipe h, mytipe tao, func ux0, func Ku, const unsigned int M) //n - кол-во узлов по пространству, k - по времени. 
{
	const unsigned int N = n - 1;//колличество разбиений

	mytipe** Y = new mytipe*[k];
	Y[0] = new mytipe[n];

	for (int i = 0; i < n; i++) { Y[0][i] = ux0(i*h, 0, n*h); }

	mytipe sigma = 2;
	mytipe kappa = 0.5;
	mytipe c = 1;
	mytipe c1 = 5;
	mytipe u0 = sigma*pow(c1, 2) / kappa;

	mytipe u0t;
	mytipe ult;


	mytipe* a = new mytipe[2];

	mytipe* diag1 = new mytipe[N - 1];
	mytipe* diag2 = new mytipe[N - 1];
	mytipe* diag3 = new mytipe[N - 1];
	mytipe* right = new mytipe[N - 1];
	mytipe* vspom1 = new mytipe[N - 1];
	mytipe** vspom2 = new mytipe*[2];//для иттераций неявного решения
	vspom2[0] = new mytipe[n];
	vspom2[1] = new mytipe[n];
	mytipe A0;
	mytipe BN;


	for (int j = 1; j < k; j++) {

		u0t = pow(u0*j*tao, 1 / sigma);
		if (n*h < c1*j*tao) { ult = pow(u0 / c1*(c1*j*tao - n*h), 1 / sigma); }
		else { ult = 0; };

		Y[j] = new mytipe[n];
		copy(n, vspom2[1], Y[j - 1]);
		vspom2[1][0] = u0t;
		vspom2[1][N] = ult;
		for (int m = 1; m <= M; m++) {
			copy(n, vspom2[0], vspom2[1]);
			a[1] = 0.5*(Ku(vspom2[0][1], sigma, kappa) + Ku(vspom2[0][0], sigma, kappa));
			for (int i = 0; i < N - 1; i++) {
				a[0] = a[1];
				a[1] = 0.5*(Ku(vspom2[0][i + 1], sigma, kappa) + Ku(vspom2[0][i], sigma, kappa));
				diag1[i] = -a[0] / pow(h, 2);
				diag3[i] = -a[1] / pow(h, 2);
				diag2[i] = -(diag1[i] + diag3[i] - c / tao);
				right[i] = c / tao*vspom2[0][i + 1];
			};
			A0 = diag1[0]; BN = diag3[N - 2];
			diag1[0] = 0.; diag3[N - 2] = 0.;

			right[0] -= A0*u0t; right[N - 2] -= BN*ult;

			vspom1 = progon3d(N - 1, diag1, diag2, diag3, right);

			for (int i = 1; i < N; i++) vspom2[1][i] = vspom1[i - 1];

			delete[]vspom1;
		}
		copy(n, Y[j], vspom2[1]);
	}

	delete[] diag1;
	delete[] diag2;
	delete[] diag3;
	delete[] right;

	delete[] a;
	delete[] vspom2[0];
	delete[] vspom2[1];
	delete[] vspom2;

	return Y;
}

mytipe** Kvazilin_implicit_Eps(const unsigned int n, const unsigned int k, mytipe h, mytipe tao, func ux0, func Ku, const mytipe Eps) //n - кол-во узлов по пространству, k - по времени. 
{
	const unsigned int N = n - 1;//колличество разбиений

	mytipe** Y = new mytipe*[k];
	Y[0] = new mytipe[n];

	for (int i = 0; i < n; i++) { Y[0][i] = ux0(i*h, 0, n*h); }

	mytipe sigma = 2;
	mytipe kappa = 0.5;
	mytipe c = 1;
	mytipe c1 = 5;
	mytipe u0 = sigma*pow(c1, 2) / kappa;

	mytipe u0t;
	mytipe ult;


	mytipe* a = new mytipe[2];

	mytipe* diag1 = new mytipe[N - 1];
	mytipe* diag2 = new mytipe[N - 1];
	mytipe* diag3 = new mytipe[N - 1];
	mytipe* right = new mytipe[N - 1];
	mytipe* vspom1 = new mytipe[N - 1];
	mytipe** vspom2 = new mytipe*[2];//для иттераций неявного решения
	vspom2[0] = new mytipe[n];
	vspom2[1] = new mytipe[n];
	mytipe A0;
	mytipe BN;

	int iter;

	for (int j = 1; j < k; j++) {

		u0t = pow(u0*j*tao, 1 / sigma);
		if (n*h < c1*j*tao) { ult = pow(u0 / c1*(c1*j*tao - n*h), 1 / sigma); }
		else { ult = 0; };

		Y[j] = new mytipe[n];
		copy(n, vspom2[1], Y[j - 1]);
		vspom2[1][0] = u0t;
		vspom2[1][N] = ult;
		iter = 0;
		do {
			iter++;
			copy(n, vspom2[0], vspom2[1]);
			a[1] = 0.5*(Ku(vspom2[0][1], sigma, kappa) + Ku(vspom2[0][0], sigma, kappa));
			for (int i = 0; i < N - 1; i++) {
				a[0] = a[1];
				a[1] = 0.5*(Ku(vspom2[0][i + 1], sigma, kappa) + Ku(vspom2[0][i], sigma, kappa));
				diag1[i] = -a[0] / pow(h, 2);
				diag3[i] = -a[1] / pow(h, 2);
				diag2[i] = -(diag1[i] + diag3[i] - c / tao);
				right[i] = c / tao*vspom2[0][i + 1];
			};
			A0 = diag1[0]; BN = diag3[N - 2];
			diag1[0] = 0.; diag3[N - 2] = 0.;

			right[0] -= A0*u0t; right[N - 2] -= BN*ult;

			vspom1 = progon3d(N - 1, diag1, diag2, diag3, right);

			for (int i = 1; i < N; i++) vspom2[1][i] = vspom1[i - 1];

			delete[]vspom1;
		} while (norm(n, minus_vect(n, vspom2[1], vspom2[0]), 'k')>Eps);
		cout << "итераций на временном слое" << j << " = " << iter << endl;
		copy(n, Y[j], vspom2[1]);
	}

	delete[] diag1;
	delete[] diag2;
	delete[] diag3;
	delete[] right;

	delete[] a;
	delete[] vspom2[0];
	delete[] vspom2[1];
	delete[] vspom2;

	return Y;
}


mytipe** Kvazilin_explicit2(int n, int k, mytipe h, mytipe tao, func ux0, func Ku) //n - кол-во узлов по пространству, k - по времени. 
{
	const int N = n ;//колличество разбиений

	mytipe** Y = new mytipe*[k];
	Y[0] = new mytipe[n];

	for (int i = 0; i < n; i++) { Y[0][i] = ux0(i*h, 0, n*h); }

	mytipe sigma = 2;
	mytipe kappa = 0.5;
	mytipe c = 1;
	mytipe c1 = 5;
	mytipe u0 = sigma*pow(c1, 2) / kappa;

	mytipe u0t;
	mytipe ult;
	mytipe mu;


	mytipe* a = new mytipe[2];

	mytipe* diag1 = new mytipe[N-1];
	mytipe* diag2 = new mytipe[N-1];
	mytipe* diag3 = new mytipe[N-1];
	mytipe* right = new mytipe[N - 1];
	mytipe* vspom = new mytipe[N - 1];
	mytipe A0;
	mytipe BN;


	for (int j = 1; j < k; j++) {

		u0t = pow(u0*j*tao, 1 / sigma);
		if (n*h < c1*j*tao) {
			Y[j] = new mytipe[n];
			a[1] = 0.5*(Ku(Y[j - 1][1], sigma, kappa) + Ku(Y[j - 1][0], sigma, kappa));
			for (int i = 0; i < N - 1; i++) {
				a[0] = a[1];
				a[1] = 0.5*(Ku(Y[j - 1][i + 1], sigma, kappa) + Ku(Y[j - 1][i], sigma, kappa));
				diag1[i] = -a[0] / pow(h, 2);
				diag3[i] = -a[1] / pow(h, 2);
				diag2[i] = -(diag1[i] + diag3[i] - c / tao);
				right[i] = c / tao*Y[j - 1][i + 1];
			};
			A0 = diag1[0]; BN = diag3[N - 2];
			diag1[0] = 0.; diag3[N - 2] = 0.;
			diag2[N - 2] += 2*BN;
			diag1[N - 2] -= BN;

			right[0] -= A0*u0t; right[N - 2] -= 0.;

			copy(N - 1, vspom, progon3d(N - 1, diag1, diag2, diag3, right));

			for (int i = 1; i < N; i++) Y[j][i] = vspom[i - 1];

			Y[j][0] = u0t;
			//Y[j][N] = Y[j][N - 1] + Y[j - 1][N] - Y[j - 1][N - 1];
		}
		else {
			ult = 0;
			Y[j] = new mytipe[n];
			a[1] = 0.5*(Ku(Y[j - 1][1], sigma, kappa) + Ku(Y[j - 1][0], sigma, kappa));
			for (int i = 0; i < N - 1; i++) {
				a[0] = a[1];
				a[1] = 0.5*(Ku(Y[j - 1][i + 1], sigma, kappa) + Ku(Y[j - 1][i], sigma, kappa));
				diag1[i] = -a[0] / pow(h, 2);
				diag3[i] = -a[1] / pow(h, 2);
				diag2[i] = -(diag1[i] + diag3[i] - c / tao);
				right[i] = c / tao*Y[j - 1][i + 1];
			};
			A0 = diag1[0]; BN = diag3[N - 2];
			diag1[0] = 0.; diag3[N - 2] = 0.;

			right[0] -= A0*u0t; right[N - 2] -= BN*ult;

			copy(N - 1, vspom, progon3d(N - 1, diag1, diag2, diag3, right));

			for (int i = 1; i < N; i++) Y[j][i] = vspom[i - 1];

			Y[j][0] = u0t;
			Y[j][N] = ult;
		};
	}



	delete[] diag1;
	delete[] diag2;
	delete[] diag3;
	delete[] right;

	delete[] a;
	delete[] vspom;

	return Y;
}


mytipe** Integro_interpolation(int n, int k, mytipe h, mytipe tao, Fanc ux0, Func u_0t, mytipe c, mytipe p, Func Kx, Func Pt)
{
	const int N = n - 1;//колличество разбиений

	mytipe** Y = new mytipe*[k];
	Y[0] = new mytipe[n];

	for (int i = 0; i < n; i++) { Y[0][i] = ux0(i*h,(n - 1)*h); }

	mytipe sigma = 0.7;

	mytipe* a = new mytipe[N];



	mytipe* diag1 = new mytipe[N - 1];
	mytipe* diag2 = new mytipe[N - 1];
	mytipe* diag3 = new mytipe[N - 1];
	mytipe* right = new mytipe[N - 1];
	mytipe* vspom = new mytipe[N - 1];
	mytipe A0, BN; mytipe mu, kappa,vspom1;
	mytipe L = N*h; mytipe psi = 0.0;


	for (int i = 0; i < N; i++) { a[i] = Kx((i + 0.5)*h); };

	kappa = (sigma*a[N - 1] / h) / (c*p*h / (2 * tao) + sigma*a[N - 1] / h);

	for (int i = 0; i < N - 1; i++) {
		diag1[i] = sigma / h*a[i];
		diag3[i] = sigma / h*a[i + 1];
		diag2[i] = -(diag1[i] + diag3[i] + c*p*h / tao);
	}

	A0 = diag1[0]; BN = diag3[N - 2];

	diag1[0] = 0.; diag3[N - 2] = 0.;

	diag2[N - 2] += BN*kappa;

	for (int j = 1; j < k; j++) {

		Y[j] = new mytipe[n];
		Y[j][0] = u_0t(j*tao);

		for (int i = 0; i < N - 1; i++) right[i] = -(c*p*h / tao*Y[j - 1][i + 1] + (1 - sigma)*a[i] * (Y[j - 1][i + 2] - 2 * Y[j - 1][i + 1] + Y[j - 1][i]) / h);

		mu = (c*p*Y[j - 1][N] * h / (2 * tao) + sigma*Pt(tao*j) + (1 - sigma)*(Pt(tao*(j - 1)) - (Y[j - 1][N] - Y[j - 1][N - 1]) / h)) / (c*p*h / (2 * tao) + sigma*a[N - 1] / h);

		right[0] -= A0*u_0t(j*tao); right[N - 2] -= BN*mu;

		copy(N - 1, vspom, progon3d(N - 1, diag1, diag2, diag3, right));
		

		for (int i = 1; i < N; i++) Y[j][i] = vspom[i - 1];

		Y[j][N] = kappa*Y[j][N - 1] + mu;
	
	}
	
	delete[] diag1;
	delete[] diag2;
	delete[] diag3;
	delete[] right;

	delete[] a;
	delete[] vspom;

	return Y;
}

mytipe** Integro_interpolation(int n, int k, mytipe h, mytipe tao, Fanc ux0, mytipe c, mytipe p, Func Kx, Func Pt1,Func Pt2)
{
	const int N = n - 1;//колличество разбиений

	mytipe** Y = new mytipe*[k];
	Y[0] = new mytipe[n];

	for (int i = 0; i < n; i++) { Y[0][i] = ux0(i*h, (n - 1)*h); }

	mytipe sigma = 1.0;

	mytipe* a = new mytipe[N];



	mytipe* diag1 = new mytipe[N - 1];
	mytipe* diag2 = new mytipe[N - 1];
	mytipe* diag3 = new mytipe[N - 1];
	mytipe* right = new mytipe[N - 1];
	mytipe* vspom = new mytipe[N - 1];
	mytipe A0, BN; mytipe mu1, kappa1,mu2,kappa2, vspom1;
	mytipe L = N*h; mytipe I = 0.0;


	for (int i = 0; i < N; i++) { a[i] = Kx((i + 0.5)*h); };

	kappa1 = (sigma*a[0] / h) / (c*p*h / (2 * tao) + sigma*a[0] / h);
	kappa2 = (sigma*a[N - 1] / h) / (c*p*h / (2 * tao) + sigma*a[N - 1] / h);

	for (int i = 0; i < N - 1; i++) {
		diag1[i] = sigma / h*a[i];
		diag3[i] = sigma / h*a[i + 1];
		diag2[i] = -(diag1[i] + diag3[i] + c*p*h / tao);
	}

	A0 = diag1[0]; BN = diag3[N - 2];

	diag1[0] = 0.; diag3[N - 2] = 0.;

	diag2[0] += A0*kappa1;
	diag2[N - 2] += BN*kappa2;

	for (int j = 1; j < k; j++) {

		Y[j] = new mytipe[n];

		for (int i = 0; i < N - 1; i++) right[i] = -(c*p*h / tao*Y[j - 1][i + 1] + (1 - sigma)*a[i] * (Y[j - 1][i + 2] - 2 * Y[j - 1][i + 1] + Y[j - 1][i]) / h);

		mu1 = (c*p*Y[j - 1][0] * h / (2 * tao) + sigma*Pt1(tao*j) + (1 - sigma)*(Pt1(tao*(j - 1)) + (Y[j - 1][1] - Y[j - 1][0]) / h)) / (c*p*h / (2 * tao) + sigma*a[0] / h);
		mu2 = (c*p*Y[j - 1][N] * h / (2 * tao) + sigma*Pt2(tao*j) + (1 - sigma)*(Pt2(tao*(j - 1)) - (Y[j - 1][N] - Y[j - 1][N - 1]) / h)) / (c*p*h / (2 * tao) + sigma*a[N - 1] / h);

		right[0] -= A0*mu1; right[N - 2] -= BN*mu2;

		copy(N - 1, vspom, progon3d(N - 1, diag1, diag2, diag3, right));


		for (int i = 1; i < N; i++) Y[j][i] = vspom[i - 1];

		Y[j][0] = kappa1*Y[j][1] + mu1;
		Y[j][N] = kappa2*Y[j][N - 1] + mu2;

		I = 0.0;
		I += 0.5*Y[j][0] * h;
		for (int i = 1; i < N; i++) {
			I+=Y[j][i]*h;
		}
		I += 0.5*Y[j][N] * h;
		cout << "I_" << j-1 << " = " << I << "; ";
	}

	delete[] diag1;
	delete[] diag2;
	delete[] diag3;
	delete[] right;

	delete[] a;
	delete[] vspom;

	return Y;
}


mytipe** Integro_interpolation(int n, int k, mytipe h, mytipe tao, Fanc ux0, Func u_0t, Fanc u_Lt, mytipe c, mytipe p, Func Kx)
{
	const int N = n - 1;//колличество разбиений

	mytipe** Y = new mytipe*[k];
	Y[0] = new mytipe[n];

	for (int i = 0; i < n; i++) { Y[0][i] = ux0(i*h, (n - 1)*h); }

	//mytipe sigma = 0.2;
	mytipe sigma = 1.;

	mytipe* a = new mytipe[N];



	mytipe* diag1 = new mytipe[N - 1];
	mytipe* diag2 = new mytipe[N - 1];
	mytipe* diag3 = new mytipe[N - 1];
	mytipe* right = new mytipe[N - 1];
	mytipe* vspom = new mytipe[N - 1];
	mytipe A0, BN; mytipe mu, kappa, vspom1;
	mytipe L = N*h; mytipe psi = 0.0;

	
	for (int i = 0; i < N; i++) { a[i] = Kx((i + 0.5)*h); };


	for (int i = 0; i < N - 1; i++) {
		diag1[i] = sigma / h*a[i];
		diag3[i] = sigma / h*a[i + 1];
		diag2[i] = -(diag1[i] + diag3[i] + c*p*h / tao);
	}

	A0 = diag1[0]; BN = diag3[N - 2];

	diag1[0] = 0.; diag3[N - 2] = 0.;

	for (int j = 1; j < k; j++) {

		Y[j] = new mytipe[n];
		Y[j][0] = u_0t(j*tao);

		for (int i = 0; i < N - 1; i++) right[i] = -(c*p*h / tao*Y[j - 1][i + 1] + (1 - sigma)*a[i] * (Y[j - 1][i + 2] - 2 * Y[j - 1][i + 1] + Y[j - 1][i]) / h);


		right[0] -= A0*u_0t(j*tao);
		right[N - 2] -= BN* u_Lt(L,j*tao);

		copy(N - 1, vspom, progon3d(N - 1, diag1, diag2, diag3, right));

		for (int i = 1; i < N; i++) Y[j][i] = vspom[i - 1];

		Y[j][N] = u_Lt(L, j*tao);


		//for (int i = 1; i < N - 1; i++) {
		//vspom1 = (sin((i + 1)*h) - 2 * sin(i*h) + sin((i - 1)*h))*exp(-j*tao) / (h*h) - sin(i*h)*(exp(-(j + 1)*tao) - exp(-j*tao)) / tao;
		//if (fabs(vspom1) > fabs(psi))  psi = vspom1; 
		//};
	
	}

	//cout << "\nпогрешность аппроксимации = " << psi << endl;

	delete[] diag1;
	delete[] diag2;
	delete[] diag3;
	delete[] right;

	delete[] a;
	delete[] vspom;

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

mytipe error(int k, int n,mytipe h, mytipe tao, mytipe** a)
{
	mytipe MAXerror = 0.0;
	//mytipe psi = 0.;
	//mytipe vspom1;
	for (int j = 0; j < k; j++)
		for (int i = 0; i < n; i++)
		if (fabs(a[j][i] - sin(i*h)*exp(-j*tao)) > MAXerror) MAXerror = fabs(a[j][i] - sin(i*h)*exp(-j*tao));

	return MAXerror;
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