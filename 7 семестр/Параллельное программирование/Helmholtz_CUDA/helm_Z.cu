#include <iostream>
#include <fstream>
#include <cmath>


#define blocksizeX 512// blocksizeX*blocksizeY == (n-2)/2 т.к. не считаем в гу
#define blocksizeY 1
#define n 1026 //количество узлов
#define iter_end 100000//кол-во итераций
#define k 100.0
#define pi 3.1415926535897932385


using namespace std;

typedef double mytipe;  //тип данных, использующийся во всей программе

__constant__ mytipe c[2];//коэффиценты в расчётной формуле запишем в константную память


mytipe f(mytipe x, mytipe y)
{
	return 2.*sin(pi*y) + k*k*(1. - x)*x*sin(pi*y) + pi*pi*(1. - x)*x*sin(pi*y);
}


__global__ void KernelBlack(mytipe* y_red, mytipe*  y_black,mytipe* f_black)
{

	// глобальное положение элемента в матрице
	int col = blockIdx.x * blockDim.x + threadIdx.x + 1;//строка и столбец элемента текущей нити
	int row = blockIdx.y * blockDim.y + threadIdx.y + 1;//пересчитываются только внутр эл-ты (без ГУ)

	int i  = row * n/2 + col;//текущий индекс и индексы соседних элементов
	int up = (row + 1) * n/2 + col;
	int down  = (row - 1) * n/2 + col;
	int left = row * n/2 + col;
	int right = row * n/2 + col + 1;


	//каждый поток вычисляет один раз
        if ((row % 2)==0){
		y_black[i - 1] = f_black[i - 1]*c[0] + (y_red[up - 1] + y_red[down - 1] + y_red[left - 1] + y_red[right - 1])*c[1];
	} else {
        	y_black[i] = f_black[i]*c[0] + (y_red[up] + y_red[down] + y_red[left - 1] + y_red[right - 1])*c[1];
        };

       /* if ((row % 2)==0){
		y_black[i - 1] = f(row*h,(2*col - 1)*h);
	} else {
        	y_black[i] = f(row*h,2*col*h);
        };*/
	
	return;
}

__global__ void KernelRed(mytipe* y_red, mytipe*  y_black,mytipe* f_red)
{

	// глобальное положение элемента в матрице
	int col = blockIdx.x * blockDim.x + threadIdx.x + 1;//строка и столбец элемента текущей нити
	int row = blockIdx.y * blockDim.y + threadIdx.y + 1;//пересчитываются только внутр эл-ты (без ГУ)

	int i  = row * n/2 + col;//текущий индекс и индексы соседних элементов
	int up = (row + 1) * n/2 + col;
	int down  = (row - 1) * n/2 + col;
	int left = row * n/2 + col;
	int right = row * n/2 + col + 1;

	//каждый поток вычисляет один раз
         if ((row % 2)==0){
		y_red[i] = f_red[i]*c[0] + (y_black[up] + y_black[down] + y_black[left-1] + y_black[right-1])*c[1];
        } else {
        	y_red[i - 1] = f_red[i - 1]*c[0] + (y_black[up - 1] + y_black[down - 1] + y_black[left - 1] + y_black[right - 1])*c[1];
        };

/*if ((row % 2)==0){
		y_red[i] = f(row*h,2*col*h);
	} else {
        	y_red[i - 1]= f(row*h,(2*col - 1)*h);
        };*/

		
	return;
}


mytipe u(mytipe x, mytipe y)
{
	return (1. - x)*x*sin(pi*y);
}

mytipe error (mytipe* y_red, mytipe* y_black, mytipe h)
{
	mytipe err = 0.0;
	mytipe maxerr_r = 0.0;
	mytipe maxerr_b = 0.0;
	int m;
	
	for (int i = 0; i < n; i++){
		m =0;
	    for (int j = (i % 2); j < n; j += 2) {
			
			err = fabs(y_red[i*n/2 + m] - u(i*h, j*h));
			if (err > maxerr_r)
				maxerr_r = err;
		++m; 
		}
	}
	
	
	for (int i = 0; i < n; i++){ 
		m=0;
		for (int j = ((i + 1) % 2); j < n; j += 2) {
			err = fabs(y_black[i*n/2 + m] - u(i*h, j*h));
			if (err > maxerr_b)
				maxerr_b = err; 
		++m;
		}
	}
	
	if(maxerr_r > maxerr_b) return maxerr_r;
	else return maxerr_b;
	

}


int main() {
	
	mytipe h = 1. / mytipe(n - 1);//делим отрезок на кол-во разбиений

	mytipe* y_red = new mytipe [n*n/2];
	mytipe* y_black = new mytipe [n*n/2];
	mytipe* F_red = new mytipe [n*n/2];
	mytipe* F_black = new mytipe [n*n/2];
			
	
	for(int i = 0; i < n; i++) 
		for(int j = 0; j < n/2; j++) {y_black[i*n/2+j]=0.;y_red[i*n/2+j]=0.;}
	
	
	int m;
	
	for (int i = 0; i < n; i++){
		m =0;
	        for (int j = (i % 2); j < n; j += 2) { F_red[i*n/2+m] =f(i*h,j*h); ++m;}
	}
     
	for (int i = 0; i < n; i++){ 
		 m=0;
		for (int j = ((i + 1) % 2); j < n; j += 2) { F_black[i*n/2+m] =f(i*h,j*h); ++m;}
	}			
	
	
	
	cudaError_t SD;

	SD = cudaSetDevice(0);
	if (SD != cudaSuccess)//проверяем подключилась ли графическая карта
	{
		cout << "CUDA set device error" << endl;
		return 1;
	}
	
	//создаем указатели на графическом ядре (девайсе)
        mytipe* dev_y_red = NULL;
	mytipe* dev_y_black = NULL;
	mytipe* dev_F_red = NULL;
	mytipe* dev_F_black = NULL;
	
	int nbytes = n*n*sizeof(mytipe)/2;

	//выделяем под них память
	cudaMalloc ((void **)&dev_y_red, nbytes);
	cudaMalloc ((void **)&dev_y_black, nbytes);
	cudaMalloc ((void **)&dev_F_red, nbytes);
        cudaMalloc ((void **)&dev_F_black, nbytes);
	
	dim3 threads(blocksizeX, blocksizeY);//кол-во тредов под один блок
	dim3 blocks((n-2)/blocksizeX/2, (n-2)/blocksizeY);//кол-во блоков. (n-2) т.к. на границе ничего не вычисляем
	
	cudaEvent_t start, stop;//счетчики время
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
		
	cudaEventRecord(start,0);
	cudaEventSynchronize(start);
	
	cudaMemcpy(dev_y_red, y_red, nbytes, cudaMemcpyHostToDevice);//копирование данных  с хоста на девайс 
	cudaMemcpy(dev_y_black, y_black, nbytes, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_F_red, F_red, nbytes, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_F_black, F_black, nbytes, cudaMemcpyHostToDevice);
	mytipe host_c[2] =  {h*h / (4.0 + h*h*k*k),1.0 / (4.0 + h*h*k*k)};
	cudaMemcpyToSymbol(c,host_c,2*sizeof(mytipe),0, cudaMemcpyHostToDevice);//передача константной памяти

	
	for(int iter = 0; iter < iter_end; iter++)
	{            	
             KernelBlack<<<blocks,threads>>>(dev_y_red,dev_y_black,dev_F_black);
            // cudaDeviceSynchronize();
             KernelRed<<<blocks,threads>>>(dev_y_red,dev_y_black,dev_F_red);
	    // cudaDeviceSynchronize();
        }



	cudaMemcpy(y_red, dev_y_red, nbytes, cudaMemcpyDeviceToHost);
	cudaMemcpy(y_black, dev_y_black, nbytes,cudaMemcpyDeviceToHost);
	cudaEventRecord(stop,0);	
	cudaEventSynchronize(stop);



	float dt;
	cudaEventElapsedTime(&dt,start,stop);
	cout << "dim "<< n << "x" << n<< endl;
	cout << "iter "<< iter_end<< endl;
	cout << "error  "<< error(y_red, y_black, h)<<endl;
	cout << "time " << dt/1000 << " seconds"<< endl;	
	cout<<"blocksize x: "<<blocksizeX<<" y: "<<blocksizeY<<endl;

	
	delete[] y_red;
	delete[] y_black;
	delete[] F_red;
	delete[] F_black;

	
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	cudaFree(dev_y_red);
	cudaFree(dev_y_black);
	cudaFree(dev_F_red);
	cudaFree(dev_F_black);
	
	return 0;

}