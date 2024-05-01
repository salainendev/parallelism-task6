#include <iostream>
#include <boost/program_options.hpp>
#include <cmath>
#include <memory>
#include </opt/nvidia/hpc_sdk/Linux_x86_64/23.11/cuda/12.3/include/nvtx3/nvToolsExt.h>


template <class ctype> class Data
    {
      private:
        /// Length of the data array
        int len;
        /// Data array
    
      public:
        ctype *arr;
        /// Class constructor
        Data(int length)
        {
          len = length;
          arr = new ctype[len];
    #pragma acc enter data copyin(this)
    #pragma acc enter data create(arr[0:len])
        }

       
        

        /// Class destructor
        ~Data()
        {
    #pragma acc exit data delete(arr)
    #pragma acc exit data delete(this)
          delete arr;
          len = 0;
        }
    };







// собственно возвращает значение линейной интерполяции
double linearInterpolation(double x, double x1, double y1, double x2, double y2) {
    // делаем значение y(щначение клетки)используя формулу линейной интерполяции
    return y1 + ((x - x1) * (y2 - y1) / (x2 - x1));
}

// эта функция берёт матрицу, проходит по её элементам и обновляет их на основе среднего по 4-м соседям , заносит в новую матрицу считает разницу между новым значением
double oneIteraton(Data<double>& prevmatrix,Data<double>& curmatrix,int N){
    
    double error = 0.0;
      // Копируем данные на устройство
    #pragma acc enter data copyin(prevmatrix.arr[0:N*N], curmatrix.arr[0:N*N])

    // Обновляем данные на устройстве перед выполнением операций
    #pragma acc update device(prevmatrix.arr[0:N*N], curmatrix.arr[0:N*N])
    

    #pragma acc parallel loop reduction(max:error)
    for (size_t i = 1; i < N-1; i++)
    {
        #pragma acc loop
        for (size_t j = 1; j < N-1; j++)
        {
            curmatrix.arr[i*N+j]  = 0.25 * (prevmatrix.arr[i*N+j+1] + prevmatrix.arr[i*N+j-1] + prevmatrix.arr[(i-1)*N+j] + prevmatrix.arr[(i+1)*N+j]);
            error = fmax(error,fabs(curmatrix.arr[i*N+j] - prevmatrix.arr[i*N+j]));
        }
    }

    // Обновляем данные на хосте после выполнения операций
    #pragma acc update self(curmatrix.arr[0:N*N])

    // Освобождаем память на устройстве
    #pragma acc exit data delete(prevmatrix.arr[0:N*N], curmatrix.arr[0:N*N])


    return error;
}
void initMatrices(Data<double>& curmatrix,int N){
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            curmatrix.arr[i*N + j] = 0;
        }
        
    }
    curmatrix.arr[0] = 10.0;
    curmatrix.arr[N-1] = 20.0;
    curmatrix.arr[(N-1)*N + (N-1)] = 30.0;
    curmatrix.arr[(N-1)*N] = 20.0;
    // up граница
    for (size_t i = 1; i < N-1; i++)
    {
        curmatrix.arr[0*N+i] = linearInterpolation(i,0.0,curmatrix.arr[0],N-1,curmatrix.arr[N-1]);
    }
    // left граница 
    for (size_t i = 1; i < N-1; i++)
    {
        curmatrix.arr[i*N+0] = linearInterpolation(i,0.0,curmatrix.arr[0],N-1,curmatrix.arr[(N-1)*N]);
    }
    // right граница

    for (size_t i = 1; i < N-1; i++)
    {
        curmatrix.arr[i*N+(N-1)] = linearInterpolation(i,0.0,curmatrix.arr[N-1],N-1,curmatrix.arr[(N-1)*N + (N-1)]);
    }
    
    // нижняя граница 
    for (size_t i = 1; i < N-1; i++)
    {
        curmatrix.arr[(N-1)*N+i] = linearInterpolation(i,0.0,curmatrix.arr[(N-1)*N],N-1,curmatrix.arr[(N-1)*N + (N-1)]);
    }

}

void swap(Data<double>& curmatrix, Data<double>& prevmatrix, int N){
    
    #pragma acc enter data copyin(curmatrix.arr[0:N*N])

    #pragma acc update device(curmatrix.arr[0:N*N])

    double* curData = curmatrix.arr;
    double* prevData = prevmatrix.arr;

    
    std::copy(curData, curData + N*N, prevData);

    
       // Обновляем данные в prevmatrix на устройстве
    #pragma acc update device(prevmatrix.arr[0:N*N])

    // Обновляем данные в prevmatrix на хосте
    #pragma acc update self(prevmatrix.arr[0:N*N])

    // Освобождаем память, занятую для curmatrix.arr на устройстве
    #pragma acc exit data delete(curmatrix.arr)
    

}

int main(){
    int N = 16;
    int accuracy = 1e-6;
    int count_iter = 100;
    Data<double> curmatrix(N*N);
    initMatrices(std::ref(curmatrix),N);
    Data<double> prevmatrix(N*N);
    swap(std::ref(curmatrix),std::ref(prevmatrix),N);
   
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            /* code */
            std::cout << curmatrix.arr[i*N+j] << ' ';
            
        }
        std::cout << std::endl;
    }
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            /* code */
            std::cout << prevmatrix.arr[i*N+j] << ' ';
            
        }
        std::cout << std::endl;
    }
    

}