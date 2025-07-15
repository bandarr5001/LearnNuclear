//read_and_write.h
#ifndef READ_AND_WRITE_H
#define READ_AND_WRITE_H

#include <iostream>
#include <cmath>
#include <string>

int check_for_data(int readMAX, int i, double *col_1, double *col_2);

int write_array(int iMAX, double *array1, double *array2);

int append_array(int iMAX, double *array1, double *array2, int i);

int read_array(int iMAX, double *col_1, double *col_2);

#endif