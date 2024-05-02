#include <iostream>
#include <boost/program_options.hpp>
#include <cmath>
#include <memory>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include </opt/nvidia/hpc_sdk/Linux_x86_64/23.11/cuda/12.3/include/nvtx3/nvToolsExt.h>
#include <omp.h>
namespace opt = boost::program_options;



// собственно возвращает значение линейной интерполяции
double linearInterpolation(double x, double x1, double y1, double x2, double y2) {
    // делаем значение y(щначение клетки)используя формулу линейной интерполяции
    return y1 + ((x - x1) * (y2 - y1) / (x2 - x1));
}

class Array
    {
      private:
        /// Length of the data array
        int len;
        /// Data array
    
      public:
        double *arr;
        /// Class constructor
        Array(int N) : len(N*N),arr(new double[N*N])
        {
        
        

          arr[0] = 10.0;
          arr[N-1] = 20.0;
          arr[(N-1)*N + (N-1)] = 30.0;
          arr[(N-1)*N] = 20.0;
              // инициализируем и потом сразу отправим на девайс
        for (size_t i = 1; i < N-1; i++)
        {
            arr[0*N+i] = linearInterpolation(i,0.0,arr[0],N-1,arr[N-1]);
            arr[i*N+0] = linearInterpolation(i,0.0,arr[0],N-1,arr[(N-1)*N]);
            arr[i*N+(N-1)] = linearInterpolation(i,0.0,arr[N-1],N-1,arr[(N-1)*N + (N-1)]);
            arr[(N-1)*N+i] = linearInterpolation(i,0.0,arr[(N-1)*N],N-1,arr[(N-1)*N + (N-1)]);
        }
        
        
    
    #pragma acc enter data create(this,arr[0:len])
        }

        /// Class destructor
    ~Array() {
        // Удаление данных на устройстве
        #pragma acc exit data delete(arr)
        
        delete[] arr;
    }
    };




// void saveMatrixToFile(const Data<double>& matrix, int N, const std::string& filename) {
//     std::ofstream outputFile(filename);
//     if (!outputFile.is_open()) {
//         std::cerr << "Unable to open file " << filename << " for writing." << std::endl;
//         return;
//     }

//     // Устанавливаем ширину вывода для каждого элемента
//     int fieldWidth = 10; // Ширина поля вывода, можно настроить по вашему усмотрению

//     // Записываем матрицу в файл с выравниванием столбцов
//     for (int i = 0; i < N; ++i) {
//         for (int j = 0; j < N; ++j) {
//             outputFile << std::setw(fieldWidth) << std::fixed << std::setprecision(4) << matrix.arr[i * N + j];
//         }
//         outputFile << std::endl;
//     }

//     outputFile.close();
// }



// double errorCalc(Data<double>& prevmatrix,Data<double>& curmatrix, int N){

//     double error = 0.0;
//         // Копируем данные на устройство
//     //#pragma acc enter data copyin(curmatrix.arr[0:N*N], prevmatrix.arr[0:N*N])

//     // Обновляем данные на устройстве
//     //#pragma acc update device(curmatrix.arr[0:N*N], prevmatrix.arr[0:N*N])
    

//     #pragma acc parallel loop reduction(max:error)
//     for (size_t i = 1; i < N-1; i++)
//     {
//         #pragma acc loop
//         for (size_t j = 1; j < N-1; j++)
//         {
            
//             error = fmax(error,fabs(curmatrix.arr[i*N+j] - prevmatrix.arr[i*N+j]));
//         }
//     }
    

    // Освобождаем память, занятую для массивов на устройстве
    //#pragma acc exit data delete(curmatrix.arr, prevmatrix.arr)

//     return error;

// }






// эта функция берёт матрицу, проходит по её элементам и обновляет их на основе среднего по 4-м соседям , заносит в новую матрицу считает разницу между новым значением
void oneIteraton(Array prevmatrix,Array curmatrix,int N){
    
    #pragma acc parallel loop present(prevmatrix.arr[0:N*N],curmatrix.arr[0:N*N])
    for (size_t i = 1; i < N-1; i++)
    {
        #pragma acc loop
        for (size_t j = 1; j < N-1; j++)
        {
            curmatrix.arr[i*N+j]  = 0.25 * (prevmatrix.arr[i*N+j+1] + prevmatrix.arr[i*N+j-1] + prevmatrix.arr[(i-1)*N+j] + prevmatrix.arr[(i+1)*N+j]);
            
        }
    }
}

// void swap(double* curmatrixarr, double* prevmatrixarr, int N){
    
//     //#pragma acc enter data copyin(curmatrix.arr[0:N*N])

    
//     #pragma acc kernels
//     {
//     double* curData = curmatrixarr;
//     double* prevData = prevmatrixarr;

//     std::copy(curData, curData + N*N, prevData);

//     }

    

    
       // Обновляем данные в prevmatrix на устройстве
    

    // Обновляем данные в prevmatrix на хосте
    // #pragma acc update self(prevmatrix.arr[0:N*N])

    // Освобождаем память, занятую для curmatrix.arr на устройстве
    //#pragma acc exit data delete(curmatrix.arr)
    

// }

void swap(Array &curmatrix,Array &prevmatrix,int N){
    #pragma acc parallel present(curmatrix.arr[0:N*N],prevmatrix.arr[0:N*N])
    {
            double* curData = curmatrix.arr;
            double* prevData = prevmatrix.arr;

            std::copy(curData, curData + N*N, prevData);
        
    }


}





int main(int argc, char const *argv[])
{
    int N = 13;
    int countIter= 100;
    double accuracy = 1e-6;
    double error = 999.0;
    int iter = 0;
    Array curmatrix(N);
    Array prevmatrix(N);
    #pragma acc enter data create(error)
    #pragma acc update device(curmatrix.arr[0:N*N],prevmatrix.arr[0:N*N],error) // инициализация матриц на гпу
    
    while (iter < 3){
            #pragma acc parallel loop present(prevmatrix.arr[0:N*N],curmatrix.arr[0:N*N])
            for (size_t i = 1; i < N-1; i++)
            {
                #pragma acc loop
                for (size_t j = 1; j < N-1; j++)
                {
                    curmatrix.arr[i*N+j]  = 0.25 * (prevmatrix.arr[i*N+j+1] + prevmatrix.arr[i*N+j-1] + prevmatrix.arr[(i-1)*N+j] + prevmatrix.arr[(i+1)*N+j]);
                    
                }
            }

                #pragma acc update self(prevmatrix.arr[0:N*N],curmatrix.arr[0:N*N])
                for (size_t i = 0; i < N; i++)
                {
                    for (size_t j = 0; j < N; j++)
                    {
                        /* code */
                        std::cout << curmatrix.arr[i*N+j] << ' ';
                        
                    }
                    std::cout << std::endl;
                }
std::cout << std::endl;

        swap(std::ref(curmatrix),std::ref(prevmatrix),N);
        iter++;


    }
    
    
    return 0;
}


// int main(int argc, char const *argv[]){

//     //парсим аргументы
//     opt::options_description desc("опции");
//     desc.add_options()
//         ("accuracy",opt::value<double>(),"точность")
//         ("cellsCount",opt::value<int>(),"размер матрицы")
//         ("iterCount",opt::value<int>(),"количество операций")
//         ("help","помощь")
//     ;

//     opt::variables_map vm;

//     opt::store(opt::parse_command_line(argc, argv, desc), vm);

//     opt::notify(vm);

//     if (vm.count("help")) {
//         std::cout << desc << "\n";
//         return 1;
//     }

    
//     // и это всё было только ради того чтобы спарсить аргументы.......



//     double start = omp_get_wtime();
//     int N = vm["cellsCount"].as<int>();
//     double accuracy = vm["accuracy"].as<double>();
//     int count_iter = vm["iterCount"].as<int>();
//     Data<double> curmatrix(N*N);
//     initMatrices(std::ref(curmatrix),N);
//     Data<double> prevmatrix(N*N);
//     swap(std::ref(curmatrix),std::ref(prevmatrix),N);
   
//     double error = 999.0;
//     int iter = 0;
//     while (count_iter>iter && error > accuracy)
//     {
//         oneIteraton(std::ref(prevmatrix),std::ref(curmatrix),N);
//         if ((iter+1)%100==0){
//         //#pragma acc update self(prevmatrix.arr[0:N*N],curmatrix.arr[0:N*N])
//         error = errorCalc(std::ref(prevmatrix),std::ref(curmatrix),N);
//         std::cout << "iteration: "<< iter+1 << ' ' << "error: " <<std::setprecision(8)<< error << std::endl;
//         }
//         swap(std::ref(curmatrix),std::ref(prevmatrix),N);
        
        
//         iter++;
//     }

//     if (N<=13){

//     for (size_t i = 0; i < N; i++)
//     {
//         for (size_t j = 0; j < N; j++)
//         {
//             /* code */
//             std::cout << curmatrix.arr[i*N+j] << ' ';
            
//         }
//         std::cout << std::endl;
//     }
//     }

//     double end = omp_get_wtime();
//     saveMatrixToFile(std::ref(curmatrix), N , "matrix.txt");
//     std::cout << end - start<<std::endl;
// }