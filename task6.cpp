#include <iostream>
#include <boost/program_options.hpp>
#include <cmath>



double linearInterpolation(double x, double x1, double y1, double x2, double y2) {
    // делаем значение y(щначение клетки)используя формулу линейной интерполяции
    return y1 + ((x - x1) * (y2 - y1) / (x2 - x1));
}

int main(int argc, char const *argv[])
{

    std::unique_ptr<double[]> matrix(new double[128*128]); // матрица температуры
    
    // заполнение матрицы нулями
    for (size_t i = 0; i < 128; i++)
    {
        for (size_t j = 0; j < 128; j++)
        {
            matrix[i*128 + j] = 0;
        }
        
    }
    // заполнение углов матрицы по очереди ( как я понял по часовой стрелке)
    matrix[0] = 10.0;
    matrix[127] = 20.0;
    matrix[127*128 + 127] =30.0;
    matrix[127*128] = 20.0; 
    // линейная интерполяция между углами матриц 

    for (size_t i = 0; i < 128; i++)
    {
        /* code */
    }
    

    return 0;
}
