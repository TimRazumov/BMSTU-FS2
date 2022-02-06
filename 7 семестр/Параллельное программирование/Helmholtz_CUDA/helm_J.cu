#include <iostream>
#include <fstream>
#include <cmath>


#define blocksizeX 32// blocksizeX*blocksizeY == n-2 т.к. не считаем в гу
#define blocksizeY 32
#define n 1026 //количество узлов
#define iter_end 100000//кол-во итераций
#define k 100.0
#define pi 3.1415926535897932385


using namespace std;

typedef double mytipe;  //тип данных, использующийся во всей программе

__constant__ mytipe c[2];//коэффиценты в расчётной формуле запишем в константную память


__global__ void KernelJacobe(mytipe* yk, mytipe*  yk1,mytipe* f)
{
	int row,col;// глобальное положение элемента в матрице
	col = blockIdx.x * blockDim.x + threadIdx.x + 1;//строка и столбец элемента текущей нити
	row = blockIdx.y * blockDim.y + threadIdx.y + 1;//пересчитываются только внутр эл-ты (без ГУ)

	int i  = row * n + col;//текущий индекс и индексы соседних элементов
	int up = (row + 1) * n + col;
	int down  = (row - 1) * n + col;
	int left = row * n + col - 1;
	int right = row * n + col + 1;

	//каждый поток вычисляет один раз
	yk1[i]=f[i]*c[0]+(yk[up]+yk[down]+yk[left]+yk[right])*c[1];
		
	return;
}


mytipe f(mytipe x, mytipe y)
{
	return 2.*sin(pi*y) + k*k*(1. - x)*x*sin(pi*y) + pi*pi*(1. - x)*x*sin(pi*y);
}

mytipe u(mytipe x, mytipe y)
{
	return (1. - x)*x*sin(pi*y);
}

mytipe error (mytipe* yk1, mytipe h)
{
	mytipe err = 0.0;
	mytipe maxerr = 0.0;
	
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			err = fabs(yk1[i*n + j] - u(i*h, j*h));
			if (err > maxerr)
				maxerr = err;
		}
	
	return maxerr;

}


int main() {
	
	mytipe h = 1. / mytipe(n - 1);//делим отрезок на кол-во разбиений

	mytipe* yk = new mytipe [n*n];
	mytipe* yk1 = new mytipe [n*n];
	mytipe* F = new mytipe [n*n];
	
	for(int i = 0; i < n; i++) 
		for(int j = 0; j < n; j++) {yk[i*n+j]=0.;yk1[i*n+j]=0.;F[i*n+j]=f(i*h,j*h);}

	int nbytes = n*n*sizeof(mytipe);
	
	
	cudaError_t SD;

	SD = cudaSetDevice(0);
	if (SD != cudaSuccess)//проверяем подключилась ли графическая карта
	{
		cout << "CUDA set device error" << endl;
		return 1;
	}
	
	//создаем указатели на графическом ядре (девайсе)
    mytipe* dev_yk = NULL;
	mytipe* dev_yk1 = NULL;
	mytipe* dev_F = NULL;

	//выделяем под них память
	cudaMalloc ((void **)&dev_yk, nbytes);
	cudaMalloc ((void **)&dev_yk1, nbytes);
	cudaMalloc ((void **)&dev_F, nbytes);

	
	dim3 threads(blocksizeX, blocksizeY);//кол-во тредов под один блок
	dim3 blocks((n-2)/blocksizeX, (n-2)/blocksizeY);//кол-во тредов под один блок. (n-2) т.к. на границе ничего не вычисляем
	
	cudaEvent_t start, stop;//счетчики время
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
		
	cudaEventRecord(start,0);
	cudaEventSynchronize(start);
	
	cudaMemcpy(dev_yk, yk, nbytes, cudaMemcpyHostToDevice);//копирование данных  с хоста на девайс 
	cudaMemcpy(dev_yk1, yk1, nbytes, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_F, F, nbytes, cudaMemcpyHostToDevice);
	mytipe host_c[2] =  {h*h / (4.0 + h*h*k*k),1.0 / (4.0 + h*h*k*k)};
	cudaMemcpyToSymbol(c,host_c,2*sizeof(mytipe),0, cudaMemcpyHostToDevice);//передача константной памяти

	
	for(int iter = 0; iter < iter_end; iter++)
	{
		KernelJacobe<<<blocks,threads>>>(dev_yk,dev_yk1,dev_F);
		swap(dev_yk, dev_yk1);
	}



	cudaMemcpy(yk1,dev_yk1,nbytes,cudaMemcpyDeviceToHost);
	//cudaMemcpy(yk,dev_yk,nbytes,cudaMemcpyDeviceToHost);
	cudaEventRecord(stop,0);	
	cudaEventSynchronize(stop);
	
	float dt;
	cudaEventElapsedTime(&dt,start,stop);
	cout << "dim "<< n << "x" << n<< endl;
	cout << "iter "<< iter_end<< endl;
	cout << "error  "<< error(yk1,h)<<endl;
	cout << "time " << dt/1000 << " seconds"<< endl;	
	cout<<"blocksize x: "<<blocksizeX<<" y: "<<blocksizeY<<endl;

	
	delete[] yk;
	delete[] yk1;
	delete[] F;

	
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	cudaFree(dev_yk);
	cudaFree(dev_yk1);
	cudaFree(dev_F);
	
	return 0;

}



/*__global__ void KernelJacobeLoc(double* yk, double*  yk1,double* f)
{
	int row,col;// глобальное положение элемента в матрице
	col = blockIdx.x * blockDim.x + threadIdx.x + 1;//строка и столбец элемента текущей нити
	row = blockIdx.y * blockDim.y + threadIdx.y + 1;//пересчитываются только внутр эл-ты (без ГУ)

	int i  = row * n + col;//текущий индекс и индексы соседних элементов
	int up = (row + 1) * n + col;
	int down  = (row - 1) * n + col;
	int left = row * n + col - 1;
	int right = row * n + col + 1;
	
	const int Nx_loc = blocksizeX+2;//(+2 для ГУ)
	const int Ny_loc = blocksizeY+2;

	__shared__ double  yk_loc[Nx_loc*Ny_loc];//лок матрица под блок 
	
	int col_loc = threadIdx.x + 1;//+1 т.к. в гу не считаем
	int row_loc = threadIdx.x + 1;
	
	//int i_loc  = row_loc * Nx_loc + col_loc;//текущий индекс и индексы соседних элементов в лок матр
	int up_loc = (row_loc + 1)* Nx_loc + col_loc;
	int down_loc  = (row_loc - 1) * Nx_loc + col_loc;
	int left_loc = row_loc * Nx_loc + col_loc - 1;
	int right_loc = row_loc * Nx_loc + col_loc + 1;
	
	yk_loc[up_loc] = yk[up];
	yk_loc[down_loc] = yk[down];
	yk_loc[left_loc] = yk[left];
	yk_loc[right_loc] = yk[right];
	
	//каждый поток вычисляет один раз
	yk1[i] = f[i]*c[0] + (yk_loc[up_loc] + yk_loc[down_loc] + yk_loc[left_loc] + yk_loc[right_loc])*c[1];
		
	return;
}*/