# ABC-lab4
# Оптимизация доступа к памяти.

* На языке С/С++/C# реализовать функцию DGEMM BLAS последовательное умножение двух
квадратных матрицс элементами типа double. Обеспечить возможность задавать
размерностиматриц в качестве аргумента командной строки при запуске программы.
Инициализировать начальные значения матриц случайными числами.
 * Провести серию испытаний и построить график зависимости времени выполнения
программы от объёма входных данных. Например, для квадратных матриц с числом
строк/столбцов 1000, 2000, 3000, … 10000. 
 * Оценить предельные размеры матриц, которые можно перемножить на вашем
вычислительном устройстве. 
 * Реализовать дополнительную функцию DGEMM_opt_1, в которой выполняется
оптимизация доступа к памяти, за счет построчного перебора элементов обеих матриц.
 * Реализовать дополнительную функцию DGEMM_opt_2, в которой выполняется
оптимизация доступа к памяти, за счет блочного перебора элементов матриц. Обеспечить
возможность задавать блока, в качестве аргументафункции. 
 * Оценить ускорение умножения для матриц фиксированного размера, например,
1000х1000, 2000х2000, 5000х5000, 10000х10000. 
 * Для блочного умножения матриц определить размер блока, при котором достигается
максимальное ускорение.  
--------------------------------
* Предельные размеры матрицы: возьмем 80% от имеющихся 16 Gb RAM, получается свободных 12,8 Gb RAM для умножения матриц, 1 элемент типа double занимает 8 байт, поэтому делим 12,8 Gb на 8 байт. Получаем 1,6 млрд элементов, нам необходимо выделить память под 3 матрицы, делим 1,6 на 3 - получаем приблизительно 530 млн элементов для 1 матрицы. Этого хватит для умножения квадратных матриц размера приблизительно 23.000 х 23.000
