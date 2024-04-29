#include <iostream>
#include <boost/program_options.hpp>
#include <cmath>
#include <memory>


// собственно возвращает значение линейной интерполяции
double linearInterpolation(double x, double x1, double y1, double x2, double y2) {
    // делаем значение y(щначение клетки)используя формулу линейной интерполяции
    return y1 + ((x - x1) * (y2 - y1) / (x2 - x1));
}

// эта функция берёт матрицу, проходит по её элементам и обновляет их на основе среднего по 4-м соседям , заносит в новую матрицу считает разницу между новым значением
double oneIteraton(std::unique_ptr<double[]>& prevmatrix,std::unique_ptr<double[]>& curmatrix,int N){
    
    double error = 0.0;

    
    #pragma acc parallel loop reduction(max:error) present(prevmatrix,curmatrix)
    for (size_t i = 1; i < N-1; i++)
    {
        for (size_t j = 1; j < N-1; j++)
        {
            curmatrix[i*N+j]  = 0.25 * (prevmatrix[i*N+j+1] + prevmatrix[i*N+j-1] + prevmatrix[(i-1)*N+j] + prevmatrix[(i+1)*N+j]);
            error = std::max(error,std::abs(curmatrix[i*N+j] - prevmatrix[i*N+j]));
        }
    }

    
    #pragma acc parallel loop present(prevmatrix,curmatrix)
    for (size_t i = 1; i < N-1; i++)
    {
        for (size_t j = 1; j < N-1; j++)
        {
            prevmatrix[i*N+j] = curmatrix[i*N+j];
        }
    }
    return error;
}






int main(int argc, char const *argv[])
{
    int N = 16;
    int count_iteration = 10;
    double precision  =  1e-6;
    std::unique_ptr<double[]> prevmatrix(new double[N*N]); // матрица температуры на предыдущем шаге
    std::unique_ptr<double[]> curmatrix(new double[N*N]); // матрица температуры на текущем шаге
    // заполнение матрицы нулями
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            curmatrix[i*N + j] = 0;
        }
        
    }
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            curmatrix[i*N + j] = 0;
        
        }
    }
    // заполнение углов матрицы по очереди ( как я понял по часовой стрелке)
    curmatrix[0] = 10.0;
    curmatrix[N-1] = 20.0;
    curmatrix[(N-1)*N + (N-1)] =30.0;
    curmatrix[(N-1)*N] = 20.0; 
    // линейная интерполяция между углами матриц 
    
    // верхняя граница 
    for (size_t i = 1; i < N-1; i++)
    {
        curmatrix[0*N+i] = linearInterpolation(i,0.0,curmatrix[0],N-1,curmatrix[N-1]);
    }
    // left граница 
    for (size_t i = 1; i < N-1; i++)
    {
        curmatrix[i*N+0] = linearInterpolation(i,0.0,curmatrix[0],N-1,curmatrix[(N-1)*N]);
    }
    // right граница

    for (size_t i = 1; i < N-1; i++)
    {
        curmatrix[i*N+(N-1)] = linearInterpolation(i,0.0,curmatrix[N-1],N-1,curmatrix[(N-1)*N + (N-1)]);
    }
    
    // нижняя граница 
    for (size_t i = 1; i < N-1; i++)
    {
        curmatrix[(N-1)*N+i] = linearInterpolation(i,0.0,curmatrix[(N-1)*N],N-1,curmatrix[(N-1)*N + (N-1)]);
    }
    

    // проинициализировали начальный уровень бытия матриц 
    // перед началом алгоритма предыдущая матрица равна текущей поэтому копируем её значения 

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            prevmatrix[i*N + j] = curmatrix[i*N + j];
        }
        
    }
    // сколько прошло итераций 
    int iter = 0;
    double error= 9999.0;

    while (iter < count_iteration && error>precision){
        error = oneIteraton(std::ref(prevmatrix),std::ref(curmatrix),N);
        iter++;
        std::cout << "iteration: " << iter << " " << "precision: " << error << std::endl;
    }
    // реализован продвинутый уровень бытия матриц 

    // это для вывода матрицы на экран (тест)
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            std::cout << curmatrix[i*N+j] << " ";
        }
       std::cout << std::endl;
    }
    
    

    return 0;
}
