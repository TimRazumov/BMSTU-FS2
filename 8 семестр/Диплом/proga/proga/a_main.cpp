#include <iostream>
#include <fstream>
#include <typeinfo>
#include <cmath>
#include <vector>
#include <string>
#include <omp.h>
// #include <functional> 

#define NUM_TREADS 4
#define ZERO 1e-8 // константа сравнения с 0
#define PI 3.14159265358979324

using namespace std;
typedef double num_type;
typedef std::vector<num_type> Vector;
typedef std::vector<Vector> Matrix;
//using Func = std::function<num_type(const num_type, const num_type)>;
typedef num_type(*Func)(const num_type, const num_type);

typedef struct point
{
	num_type x1;
	num_type x2;
} Point;

typedef struct rect_mesh
{
	Point p_l;
	Point p_r;
	Point step;
} Rect_mesh;

typedef struct rect_mesh_time
{
	Rect_mesh rect_step;
	num_type time_end;
	num_type tau;
} Rect_mesh_t;

void write_file(const Matrix& S, const Matrix& I, const Matrix& R, int num_layer);
void read_file(Matrix& S, Matrix& I, Matrix& R, int num_layer);
Vector progon3d(const Vector& diag_l, const Vector& diag_c, const Vector& diag_r, const Vector& rhs);

num_type It0(const num_type x, const num_type y)
{
	return pow(0.5, 39)*(50 - x)*(50 - x)*x*x*(50 - y)*(50 - y)*y*y*exp(-pow(0.5, 3)*((x - 35)*(x - 35) + (y - 35)*(y - 35)));
	//return pow(0.49, 38)*(50 - x)*(50 - x)*x*x*(50 - y)*(50 - y)*y*y*exp(-((x - 25)*(x - 25) + (y - 25)*(y - 25)) / 10.);
	/*if (x < 3.)
		return 1.;//0.6*cos(PI*x / 6.)*cos(PI*x / 6.);
	return 0.;*/

}

num_type St0(const num_type x, const num_type y)
{
	return pow(0.5, 39)*(50 - x)*(50 - x)*x*x*(50 - y)*(50 - y)*y*y*exp(-pow(0.5, 8)*((x - 10)*(x - 10) + (y - 10)*(y - 10)));
	// return pow(0.5, 37)*(50 - x)*(50 - x)*x*x*(50 - y)*(50 - y)*y*y*exp(-((x - 25)*(x - 25) + (y - 25)*(y - 25)) / 90.);	
	/*if (x > 6.)
		return 1.;//cos(PI*x / 10.)*cos(PI*x / 10.);
	return 0.;*/
}

num_type Rt0(const num_type x, const num_type y)
{
	return fabs(0);
}

num_type Z(num_type xi, num_type yj, num_type xk, num_type ym)//интегральное ядро уравнения
{
	num_type Z0 = 1e-2, L0 = 1.;
	return Z0*exp(-sqrt((xi - xk)*(xi - xk) + (yj - ym)*(yj - ym)) / L0);
}

num_type D(num_type u)//коэффициент диффузии
{
	num_type m = 5.;
	num_type D0 = 1.;
	return  D0*exp(-1. / (m*m*u*u));
}

num_type Integrate(const Rect_mesh& rect, const Matrix& f); // двумерный метод центральных прямоугольников

															// нестационарный процесс
void fractional_step(
	const Rect_mesh_t& mesh,
	Func St0, Func It0, Func Rt0, // НУ
	num_type Q  // скорость вымирания
);

int main() {
	setlocale(LC_ALL, "Russian"); //подключаем русский язык
	num_type L = 50.;
	const Rect_mesh_t mesh{ { { 0, 0 },{ L, L },{ 1.25, 1.25 }/* 41 узел */ }, 1., 0.1 };

	num_type Q = 4*1e-3;

	fractional_step(mesh, St0, It0, Rt0, Q);

	std::cout << "\nend:)";

	cin.get();
	return 0;
}

void read_file(Matrix& S, Matrix& I, Matrix& R, int num_layer)
{
	ifstream f_S, f_I, f_R;    // создали переменную для записи в файл
	string file_name = "S_layer_" + to_string(num_layer) + ".dat";
	f_S.open(file_name, ios::out);
	file_name = "I_layer_" + to_string(num_layer) + ".dat";
	f_I.open(file_name, ios::out);
	file_name = "R_layer_" + to_string(num_layer) + ".dat";
	f_R.open(file_name, ios::out);

	for (int i = 0; i < S.size(); i++)
		for (int j = 0; j < S[i].size(); j++)
		{
			f_S >> S[i][j];
			f_I >> I[i][j];
			f_R >> R[i][j];
		}
	f_S.close(); f_S.clear();
	f_I.close(); f_I.clear();
	f_R.close(); f_R.clear();
}

void write_file(const Matrix& S, const Matrix& I, const Matrix& R, int num_layer)
{
	ofstream f_S, f_I, f_R;    // создали переменную для записи в файл
	string file_name = "S_layer_" + to_string(num_layer) + ".dat";
	f_S.open(file_name, ios::out);
	file_name = "I_layer_" + to_string(num_layer) + ".dat";
	f_I.open(file_name, ios::out);
	file_name = "R_layer_" + to_string(num_layer) + ".dat";
	f_R.open(file_name, ios::out);
	for (int i = 0, rows = S.size(); i < rows; i++)
	{
		for (int j = 0, cols = S[i].size(); j < cols; j++)
		{
			f_S << S[i][j] << " "; f_I << I[i][j] << " "; f_R << R[i][j] << " ";
		}
		f_S << endl; f_I << endl; f_R << endl;
	}
	f_S.close(); f_S.clear();
	f_I.close(); f_I.clear();
	f_R.close(); f_R.clear();
}


void fractional_step(
	const Rect_mesh_t& mesh,
	Func St0, Func It0, Func Rt0, // НУ 
	num_type Q  // скорость вымирания
)
{
	const int freq = 10, // частота записи слоев в файл
		iter_end = 3; // кол-во итераций во втором измерении
	int num_file = 1; // = 1; // номер файла
	num_type t = 0.; // время расчета
	omp_set_num_threads(NUM_TREADS); // кол-во тредов, на которых ведется расчет
	const num_type h1 = mesh.rect_step.step.x1, h2 = mesh.rect_step.step.x2, tau = mesh.tau;
	const int N1 = (mesh.rect_step.p_r.x1 - mesh.rect_step.p_l.x1) / h1 + 1,
		N2 = (mesh.rect_step.p_r.x2 - mesh.rect_step.p_l.x2) / h2 + 1;

	ofstream midd;
	midd.open("SIR_midd.dat", ios::out);
	Matrix S0(N1, Vector(N2)), S05(N1, Vector(N2)), S1(N1, Vector(N2)),
		I0(N1, Vector(N2)), I05(N1, Vector(N2)), I1(N1, Vector(N2)),
		R0(N1, Vector(N2)), R05(N1, Vector(N2)), R1(N1, Vector(N2));

	Matrix Int(N1, Vector(N2));

	// Коэффициенты для 1-го измерения
	Vector diag1_1(N1), diag2_1(N1), diag3_1(N1), right1_S(N1), right1_I(N1), right1_R(N1);

	// Коэффициенты для 2-го измерения
	Vector diag1_2(N2), diag2_2(N2), diag3_2(N2), right2_S(N2), right2_I(N2), right2_R(N2);

	// НУ
	if (num_file == 1)
	{
		for (int i = 0; i < N1; i++)
			for (int j = 0; j < N2; j++)
			{
				S0[i][j] = St0(i*h1, j*h2);
				I0[i][j] = It0(i*h1, j*h2);
				R0[i][j] = Rt0(i*h1, j*h2);
			}
		write_file(S0, I0, R0, 0);
		midd << 0 << " " << Integrate(mesh.rect_step, S0) << " " << Integrate(mesh.rect_step, I0)
			<< " " << Integrate(mesh.rect_step, R0) << std::endl;
	}
	else
	{
		read_file(S0, I0, R0, num_file - 1);
	}

	// проверка консервативности
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
		{
			S05[i][j] = S0[i][j] + I0[i][j] + R0[i][j];
		}

	num_type N0 = Integrate(mesh.rect_step, S05); // Общая численность
												  /*std::cout << "S0 = " << Integrate(mesh.rect_step, S0) / N0 << "   I0 = " << Integrate(mesh.rect_step, I0) / N0
												  << "    R0 = " << Integrate(mesh.rect_step, R0) / N0 << endl;*/

	std::cout << "N0_init = " << N0 << endl;

	for (int k = 1; k <= 3*4000/*4 * 5 * 600*/; k++)
	{
		t -= omp_get_wtime();

		// ПЕРВОЕ ИЗМЕРЕНИЕ
#pragma omp parallel for
		for (int i = 0; i < N1; i++) //  Z(i*h1, j*h2, n*h1, m*h2)
			for (int j = 0; j < N2; j++)
			{
				Int[i][j] = 0.;
				for (int n = 0; n < N1 - 1; n++)
				{
					for (int m = 0; m < N2 - 1; m++)
					{
						Int[i][j] += 0.25*h1*h2*(Z(i*h1, j*h2, n*h1, m*h2)*I0[n][m] + Z(i*h1, j*h2, n*h1, (m + 1)*h2)*I0[n][m + 1]
							+ Z(i*h1, j*h2, (n + 1)*h1, m*h2)*I0[n + 1][m] + Z(i*h1, j*h2, (n + 1)*h1, (m + 1)*h2)*I0[n + 1][m + 1]);
					}
				}
				Int[i][j] *= N0*S0[i][j];
			}

		// j = 0
		diag1_1[0] = 0.;
		diag3_1[0] = (D(R0[1][0]) + D(R0[0][0])) / (h1*h1);
		diag2_1[0] = -(diag3_1[0] + 2. / tau);

		for (int i = 0; i < N1; i++)
		{
			right1_S[i] = -(2. / tau*S0[i][0]
				+ (D(R0[i][1]) + D(R0[i][0])) / (h2*h2)*(S0[i][1] - S0[i][0])
				- Int[i][0]);
			right1_I[i] = -(2. / tau*I0[i][0]
				+ (D(R0[i][1]) + D(R0[i][0])) / (h2*h2)*(I0[i][1] - I0[i][0])
				+ Int[i][0] - Q*I0[i][0]);
			right1_R[i] = -(2. / tau*R0[i][0]
				+ (D(R0[i][1]) + D(R0[i][0])) / (h2*h2)*(R0[i][1] - R0[i][0])
				+ Q*I0[i][0]);
			if (i != N1 - 1 && i != 0)
			{
				diag1_1[i] = 0.5*(D(R0[i][0]) + D(R0[i - 1][0])) / (h1*h1);
				diag3_1[i] = 0.5*(D(R0[i + 1][0]) + D(R0[i][0])) / (h1*h1);
				diag2_1[i] = -(diag1_1[i] + diag3_1[i] + 2. / tau);
			}
		}

		diag1_1[N1 - 1] = (D(R0[N1 - 1][0]) + D(R0[N1 - 2][0])) / (h1*h1);
		diag2_1[N1 - 1] = -(diag1_1[N1 - 1] + 2. / tau);
		diag3_1[N1 - 1] = 0.;
		{
			Vector vspom = progon3d(diag1_1, diag2_1, diag3_1, right1_S);
			for (int i = 0; i < N1; i++)  S05[i][0] = vspom[i];
			vspom = progon3d(diag1_1, diag2_1, diag3_1, right1_I);
			for (int i = 0; i < N1; i++)  I05[i][0] = vspom[i];
			vspom = progon3d(diag1_1, diag2_1, diag3_1, right1_R);
			for (int i = 0; i < N1; i++)  R05[i][0] = vspom[i];
		}
		// j = 1..N2-1
		// #pragma omp parallel for private(right1_S, right1_I, right1_R, diag1_1, diag2_1, diag3_1)
		for (int j = 1; j < N2 - 1; j++)
		{

			diag1_1[0] = 0.;
			diag3_1[0] = (D(R0[1][j]) + D(R0[0][j])) / (h1*h1);
			diag2_1[0] = -(diag3_1[0] + 2. / tau);

			for (int i = 0; i < N1; i++)
			{
				right1_S[i] = -(2. / tau *S0[i][j]
					+ 0.5*(D(R0[i][j + 1]) + D(R0[i][j])) / (h2*h2)*(S0[i][j + 1] - S0[i][j])
					- 0.5*(D(R0[i][j]) + D(R0[i][j - 1])) / (h2*h2)*(S0[i][j] - S0[i][j - 1])
					- Int[i][j]);
				right1_I[i] = -(2. / tau *I0[i][j]
					+ 0.5*(D(R0[i][j + 1]) + D(R0[i][j])) / (h2*h2)*(I0[i][j + 1] - I0[i][j])
					- 0.5*(D(R0[i][j]) + D(R0[i][j - 1])) / (h2*h2)*(I0[i][j] - I0[i][j - 1])
					+ Int[i][j] - Q*I0[i][j]);
				right1_R[i] = -(2. / tau *R0[i][j]
					+ 0.5*(D(R0[i][j + 1]) + D(R0[i][j])) / (h2*h2)*(R0[i][j + 1] - R0[i][j])
					- 0.5*(D(R0[i][j]) + D(R0[i][j - 1])) / (h2*h2)*(R0[i][j] - R0[i][j - 1])
					+ Q*I0[i][j]);
				if (i != 0 && i != N1 - 1)
				{
					diag1_1[i] = 0.5*(D(R0[i][j]) + D(R0[i - 1][j])) / (h1*h1);
					diag3_1[i] = 0.5*(D(R0[i + 1][j]) + D(R0[i][j])) / (h1*h1);
					diag2_1[i] = -(diag1_1[i] + diag3_1[i] + 2. / tau);
				}
			}

			diag1_1[N1 - 1] = (D(R0[N1 - 1][j]) + D(R0[N1 - 2][j])) / (h1*h1);
			diag2_1[N1 - 1] = -(diag1_1[N1 - 1] + 2. / tau);
			diag3_1[N1 - 1] = 0.;
			{
				Vector vspom = progon3d(diag1_1, diag2_1, diag3_1, right1_S);
				for (int i = 0; i < N1; i++)  S05[i][j] = vspom[i];
				vspom = progon3d(diag1_1, diag2_1, diag3_1, right1_I);
				for (int i = 0; i < N1; i++)  I05[i][j] = vspom[i];
				vspom = progon3d(diag1_1, diag2_1, diag3_1, right1_R);
				for (int i = 0; i < N1; i++)  R05[i][j] = vspom[i];
			}
		}

		// j = N2
		diag1_1[0] = 0.;
		diag3_1[0] = (D(R0[N1 - 1][0]) + D(R0[N1 - 2][0])) / (h1*h1);
		diag2_1[0] = -(diag3_1[0] + 2. / tau);

		for (int i = 0; i < N1; i++)
		{
			right1_S[i] = -(2. / tau*S0[i][N2 - 1]
				+ (D(R0[i][N2 - 2]) + D(R0[i][N2 - 1])) / (h2*h2)*(S0[i][N2 - 2] - S0[i][N2 - 1])
				- Int[i][N2 - 1]);
			right1_I[i] = -(2. / tau*I0[i][N2 - 1]
				+ (D(R0[i][N2 - 2]) + D(R0[i][N2 - 1])) / (h2*h2)*(I0[i][N2 - 2] - I0[i][N2 - 1])
				+ Int[i][N2 - 1] - Q*I0[i][N2 - 1]);
			right1_R[i] = -(2. / tau*R0[i][N2 - 1]
				+ (D(R0[i][N2 - 2]) + D(R0[i][N2 - 1])) / (h2*h2)*(R0[i][N2 - 2] - R0[i][N2 - 1])
				+ Q*I0[i][N2 - 1]);
			if (i != 0 && i != N1 - 1)
			{
				diag1_1[i] = 0.5*(D(R0[i][N2 - 1]) + D(R0[i - 1][N2 - 2])) / (h1*h1);
				diag3_1[i] = 0.5*(D(R0[i + 1][N2 - 1]) + D(R0[i][N2 - 2])) / (h1*h1);
				diag2_1[i] = -(diag1_1[i] + diag3_1[i] + 2. / tau);
			}
		}

		diag1_1[N1 - 1] = (D(R0[N1 - 1][N2 - 1]) + D(R0[N1 - 2][N2 - 1])) / (h1*h1);
		diag2_1[N1 - 1] = -(diag1_1[N1 - 1] + 2. / tau);
		diag3_1[N1 - 1] = 0.;
		{
			Vector vspom = progon3d(diag1_1, diag2_1, diag3_1, right1_S);
			for (int i = 0; i < N1; i++)  S05[i][N2 - 1] = vspom[i];
			vspom = progon3d(diag1_1, diag2_1, diag3_1, right1_I);
			for (int i = 0; i < N1; i++)  I05[i][N2 - 1] = vspom[i];
			vspom = progon3d(diag1_1, diag2_1, diag3_1, right1_R);
			for (int i = 0; i < N1; i++)  R05[i][N2 - 1] = vspom[i];
		}

		/*==================================================================================================================================*/
		// ВТОРОЕ ИЗМЕРЕНИЕ
		for (int iter = 0; iter < iter_end; iter++)
		{
			if (iter != 0)
			{
#pragma omp parallel for
				for (int i = 0; i < N1; i++) // Z(i*h1, j*h2, n*h1, m*h2)
					for (int j = 0; j < N2; j++)
					{
						Int[i][j] = 0.;
						for (int n = 0; n < N1 - 1; n++)
						{
							for (int m = 0; m < N2 - 1; m++)
							{
								Int[i][j] += 0.25*h1*h2*(Z(i*h1, j*h2, n*h1, m*h2)*I0[n][m] + Z(i*h1, j*h2, n*h1, (m + 1)*h2)*I0[n][m + 1]
									+ Z(i*h1, j*h2, (n + 1)*h1, m*h2)*I0[n + 1][m] + Z(i*h1, j*h2, (n + 1)*h1, (m + 1)*h2)*I0[n + 1][m + 1]);
							}
						}
						Int[i][j] *= N0*S0[i][j];
					}
			}
			// i = 0
			diag1_2[0] = 0.;
			diag3_2[0] = (D(R0[0][1]) + D(R0[0][0])) / (h2*h2);
			diag2_2[0] = -(diag3_2[0] + 2. / tau);

			for (int j = 0; j < N2; j++)
			{
				right2_S[j] = -(2. / tau*S05[0][j]
					+ (D(R0[1][j]) + D(R0[0][j])) / (h1*h1)*(S05[1][j] - S05[0][j])
					- Int[0][j]);
				right2_I[j] = -(2. / tau*I05[0][j]
					+ (D(R0[1][j]) + D(R0[0][j])) / (h1*h1)*(I05[1][j] - I05[0][j])
					+ Int[0][j] - Q*I0[0][j]);
				right2_R[j] = -(2. / tau*R05[0][j]
					+ (D(R0[1][j]) + D(R0[0][j])) / (h1*h1)*(R05[1][j] - R05[0][j])
					+ Q*I0[0][j]);
				if (j != 0 && j != N2 - 1)
				{
					diag1_2[j] = 0.5*(D(R0[0][j]) + D(R0[0][j - 1])) / (h2*h2);
					diag3_2[j] = 0.5*(D(R0[0][j + 1]) + D(R0[0][j])) / (h2*h2);
					diag2_2[j] = -(diag1_2[j] + diag3_2[j] + 2. / tau);
				}
			}

			diag1_2[N2 - 1] = (D(R0[0][N2 - 2]) + D(R0[0][N2 - 1])) / (h2*h2);
			diag2_2[N2 - 1] = -(diag1_2[N2 - 1] + 2. / tau);
			diag3_2[N2 - 1] = 0.;

			S1[0] = progon3d(diag1_2, diag2_2, diag3_2, right2_S);
			I1[0] = progon3d(diag1_2, diag2_2, diag3_2, right2_I);
			R1[0] = progon3d(diag1_2, diag2_2, diag3_2, right2_R);

			// i = 1..N1-1
			// #pragma omp parallel for private(right2_S, right2_I, right2_R, diag1_2, diag2_2, diag3_2)
			for (int i = 1; i < N1 - 1; i++)
			{
				diag1_2[0] = 0.;
				diag3_2[0] = (D(R0[i][1]) + D(R0[i][0])) / (h2*h2);
				diag2_2[0] = -(diag3_2[0] + 2. / tau);

				for (int j = 0; j < N2; j++)
				{
					right2_S[j] = -(2. / tau*S05[i][j]
						+ 0.5*(D(R0[i + 1][j]) + D(R0[i][j]))*(S05[i + 1][j] - S05[i][j]) / (h1*h1)
						- 0.5*(D(R0[i][j]) + D(R0[i - 1][j]))*(S05[i][j] - S05[i - 1][j]) / (h1*h1)
						- Int[i][j]);
					right2_I[j] = -(2. / tau*I05[i][j]
						+ 0.5*(D(R0[i + 1][j]) + D(R0[i][j]))*(I05[i + 1][j] - I05[i][j]) / (h1*h1)
						- 0.5*(D(R0[i][j]) + D(R0[i - 1][j]))*(I05[i][j] - I05[i - 1][j]) / (h1*h1)
						+ Int[i][j] - Q*I0[i][j]);
					right2_R[j] = -(2. / tau*R05[i][j]
						+ 0.5*(D(R0[i + 1][j]) + D(R0[i][j]))*(R05[i + 1][j] - R05[i][j]) / (h1*h1)
						- 0.5*(D(R0[i][j]) + D(R0[i - 1][j]))*(R05[i][j] - R05[i - 1][j]) / (h1*h1)
						+ Q*I0[i][j]);
					if (j != 0 && j != N2 - 1)
					{
						diag1_2[j] = 0.5*(D(R0[i][j]) + D(R0[i][j - 1])) / (h2*h2);
						diag3_2[j] = 0.5*(D(R0[i][j + 1]) + D(R0[i][j])) / (h2*h2);
						diag2_2[j] = -(diag1_2[j] + diag3_2[j] + 2. / tau);
					}
				}

				diag1_2[N2 - 1] = (D(R0[i][N2 - 2]) + D(R0[i][N2 - 1])) / (h2*h2);
				diag2_2[N2 - 1] = -(diag1_2[N2 - 1] + 2. / tau);
				diag3_2[N2 - 1] = 0.;

				S1[i] = progon3d(diag1_2, diag2_2, diag3_2, right2_S);
				I1[i] = progon3d(diag1_2, diag2_2, diag3_2, right2_I);
				R1[i] = progon3d(diag1_2, diag2_2, diag3_2, right2_R);
			}

			// i = N1
			diag1_2[0] = 0.;
			diag3_2[0] = (D(R0[N1 - 1][1]) + D(R0[N1 - 1][0])) / (h2*h2);
			diag2_2[0] = -(diag3_2[0] + 2. / tau);

			for (int j = 0; j < N2; j++)
			{
				right2_S[j] = -(2. / tau*S05[N1 - 1][j]
					+ (D(R0[N1 - 2][j]) + D(R0[N1 - 1][j])) / (h1*h1)*(S05[N1 - 2][j] - S05[N1 - 1][j])
					- Int[N1 - 1][j]);
				right2_I[j] = -(2. / tau*I05[N1 - 1][j]
					+ (D(R0[N1 - 2][j]) + D(R0[N1 - 1][j])) / (h1*h1)*(I05[N1 - 2][j] - I05[N1 - 1][j])
					+ Int[N1 - 1][j] - Q*I0[N1 - 1][j]);
				right2_R[j] = -(2. / tau*R05[N1 - 1][j]
					+ (D(R0[N1 - 2][j]) + D(R0[N1 - 1][j])) / (h1*h1)*(R05[N1 - 2][j] - R05[N1 - 1][j])
					+ Q*I0[N1 - 1][j]);
				if (j != 0 && j != N2 - 1)
				{
					diag1_2[j] = 0.5*(D(R0[N1 - 1][j]) + D(R0[N1 - 1][j - 1])) / (h2*h2);
					diag3_2[j] = 0.5*(D(R0[N1 - 1][j + 1]) + D(R0[N1 - 1][j])) / (h2*h2);
					diag2_2[j] = -(diag1_2[j] + diag3_2[j] + 2. / tau);
				}
			}
			diag1_2[N2 - 1] = (D(R0[N1 - 1][N2 - 2]) + D(R0[N1 - 1][N2 - 1])) / (h2*h2);
			diag2_2[N2 - 1] = -(diag1_2[N2 - 1] + 2. / tau);
			diag3_2[N2 - 1] = 0.;

			S1[N1 - 1] = progon3d(diag1_2, diag2_2, diag3_2, right2_S);
			I1[N1 - 1] = progon3d(diag1_2, diag2_2, diag3_2, right2_I);
			R1[N1 - 1] = progon3d(diag1_2, diag2_2, diag3_2, right2_R);

			std::swap(S1, S0);
			std::swap(I1, I0);
			std::swap(R1, R0);
		}

		/*==================================================================================================================================*/

		if (k % freq == 0)
		{
			write_file(S0, I0, R0, num_file);
			num_file++;
			/*for (int i = 0; i < N1 - 1; i++)
			{
			std::cout << "x_i_0 = "<< i*h1 << " " <<(S0[i + 1][0] - S0[i][0]) / h1 << std::endl;
			std::cout << "x_i_L = " << i*h1 << " " << (S0[i + 1][N1 - 1] - S0[i][N1 - 1]) / h1 << std::endl;
			std::cout << "y_i_0 = " << i*h1 << " " << (S0[0][i + 1] - S0[0][i]) / h1 << std::endl;
			std::cout << "y_i_L = " << i*h1 << " " << (S0[N1 - 1][i + 1] - S0[N1 - 1][i]) / h1 << std::endl;
			}*/
			midd << k*tau << " " << Integrate(mesh.rect_step, S0) << " " << Integrate(mesh.rect_step, I0)
				<< " " << Integrate(mesh.rect_step, R0) << std::endl;
		}
		t += omp_get_wtime();
	} // end for T

	  // проверка консервативности
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
		{
			S05[i][j] = S0[i][j] + I0[i][j] + R0[i][j];
		}

	std::cout << "N0_end = " << Integrate(mesh.rect_step, S05) << endl << "\nt = " << t;
	midd.close(); midd.clear();
}

Vector progon3d(const Vector& diag_l, const Vector& diag_c, const Vector& diag_r, const Vector& rhs)
{ // d_l*x_{i-1} + d_c*x_i + d_r*x_{i+1} = rhs_i
	const int dim = diag_l.size();
	Vector alfa(dim), betta(dim);
	alfa[0] = -diag_r[0] / diag_c[0]; betta[0] = rhs[0] / diag_c[0];

	for (int i = 1; i < dim; i++)
	{
		num_type znam = diag_c[i] + diag_l[i] * alfa[i - 1];
		alfa[i] = -diag_r[i] / znam;
		betta[i] = (rhs[i] - diag_l[i] * betta[i - 1]) / znam;
	}

	for (int i = dim - 2; i > -1; i--) betta[i] += alfa[i] * betta[i + 1];

	return  betta;
}

num_type Integrate(const Rect_mesh& rect, const Matrix& f)
{
	num_type I = 0.;
	const num_type h1 = rect.step.x1, h2 = rect.step.x2;
	const int n1 = (rect.p_r.x1 - rect.p_l.x1) / h1, n2 = (rect.p_r.x2 - rect.p_l.x2) / h2; // n - кол-во разбиений
	for (int i = 0; i < n1; i++)
		for (int j = 0; j < n2; j++)
			I += 0.25*h1*h2*(f[i][j] + f[i][j + 1] + f[i + 1][j] + f[i + 1][j + 1]);

	return I;
}