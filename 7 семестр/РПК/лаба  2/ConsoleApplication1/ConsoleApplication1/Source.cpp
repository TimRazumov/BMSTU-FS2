#include <iostream>
#include <fstream>

using namespace std;

struct netStruct {
	int nx; //количество узлов по x
	int ny; //количество узлов по y
	double** uzel;
	double** triangle;
	double** boundary;
};

int main()
{
	double ax, ay, bx, by, h1, h2;
	int nx, ny;
	//	ofstream ofile("element.txt");
	//	ofile.precision(16);

	//cout << "Enter ax = "; cin >> ax;
	//cout << "\nEnter bx = "; cin >> bx;
	//cout << "\nEnter ay = "; cin >> ay;
	//cout << "\nEnter by = "; cin >> by;
	//cout << "\nEnter nx = "; cin >> nx;
	//cout << "\nEnter ny = "; cin >> ny;

	ax = 1;
	bx = 6;
	ay = 1;
	by = 6;
	nx = 5;
	ny = 5;

	h1 = (bx - ax) / (nx);
	h2 = (by - ay) / (ny);

	netStruct net;
	net.nx = nx + 1;
	net.ny = ny + 1;



	net.uzel = new double*[net.nx*net.ny];
	for (int i = 0; i < net.nx*net.ny; i++)
		net.uzel[i] = new double[2];

	for (int i = 0; i < net.nx; i++)
		for (int j = 0; j < net.ny; j++)
		{
			net.uzel[i*net.ny + j][0] = ax + h1*i;
			net.uzel[i*net.ny + j][1] = ay + h2*j;
		}



	net.triangle = new double*[2 * (net.nx - 1)*(net.ny - 1)];
	for (int i = 0; i < 2 * (net.nx - 1)*(net.ny - 1); i++)
		net.triangle[i] = new double[3];

	for (int i = 0; i < net.nx - 1; i++)
		for (int j = 0; j < net.ny - 1; j++)
		{
			net.triangle[2 * (i * (net.ny - 1) + j)][0] = i*(net.ny) + j + 1;
			net.triangle[2 * (i * (net.ny - 1) + j)][1] = (i + 1)*(net.ny) + j + 1;
			net.triangle[2 * (i * (net.ny - 1) + j)][2] = (i + 1)*(net.ny) + j + 1 + 1;
			net.triangle[2 * (i * (net.ny - 1) + j) + 1][0] = i*(net.ny) + j + 1;
			net.triangle[2 * (i * (net.ny - 1) + j) + 1][1] = (i + 1)*(net.ny) + j + 1 + 1;
			net.triangle[2 * (i * (net.ny - 1) + j) + 1][2] = i*(net.ny) + j + 1 + 1;
		}



	net.boundary = new double*[2 * (net.nx + net.ny - 2)];
	for (int i = 0; i < net.nx*net.ny; i++)
		net.boundary[i] = new double[2];

	for (int j = 0; j < net.ny - 1; j++)
	{
		net.boundary[j][0] = j * 2 + 2;
		net.boundary[j][1] = 3;
		net.boundary[net.ny + net.nx - 2 + j][0] = 2 * (net.nx - 1)*(net.ny - 1) - (j * 2 + 1);
		net.boundary[net.ny + net.nx - 2 + j][1] = 2;
	}
	for (int i = 0; i < net.nx - 1; i++)
	{
		net.boundary[net.ny - 1 + i][0] = (i + 1) * 2 * (net.ny - 1);
		net.boundary[net.ny - 1 + i][1] = 2;
		net.boundary[2 * net.ny + net.nx - 3 + i][0] = 2 * (net.nx - 1)*(net.ny - 1) - (i + 1) * 2 * (net.ny - 1) + 1;
		net.boundary[2 * net.ny + net.nx - 3 + i][1] = 1;
	}





	ofstream fileUzel("fileUzel.dat");
	fileUzel.precision(16);

	for (int i = 0; i < net.nx*net.ny; i++)
		fileUzel << net.uzel[i][0] << " " << net.uzel[i][1] << endl;
	fileUzel.close();


	ofstream fileTriangle("fileTriangle.dat");
	fileTriangle.precision(16);
	for (int i = 0; i < 2 * (net.nx - 1)*(net.ny - 1); i++)
		fileTriangle << net.triangle[i][0] << " " << net.triangle[i][1] << " " << net.triangle[i][2] << endl;
	fileTriangle.close();



	ofstream fileBoundary("fileBoundary.dat");
	fileBoundary.precision(16);
	for (int i = 0; i < 2 * (net.nx + net.ny - 2); i++)
		fileBoundary << net.boundary[i][0] << " " << net.boundary[i][1] << endl;
	fileBoundary.close();



	cout << "\nthe end\n";
	cin.get();

	return 0;
}
