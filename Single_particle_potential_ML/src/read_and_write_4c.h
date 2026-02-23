//read_and_write_4c.h
#ifndef READ_AND_WRITE_4C_H
#define READ_AND_WRITE_4C_H

#include <iostream>
#include <cmath>
#include <string>

int check_for_data(int readMAX, int i, double *col_1, double *col_2);

int write_array_4c(int iMAX, double *array1, double *array2, double *array3, double *array4);

int append_array(int iMAX, double *array1, double *array2, int i);

int read_array(int iMAX, double *col_1, double *col_2);

#endif
