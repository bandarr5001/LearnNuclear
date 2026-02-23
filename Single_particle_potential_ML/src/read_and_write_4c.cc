#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include "./read_and_write_4c.h"

int check_for_data(int readMAX, int i, double *col_1, double *col_2){
    std::string FileNameCheck;

    std::cout << "File to check: ";
    std::cin >> FileNameCheck; //user input
    std::cout << "Checking " << FileNameCheck << std::endl;
    std::ifstream checkFile(FileNameCheck);
    
    double dum1 = 0.0;
    double dum2 = 0.0;
    
    if (checkFile.is_open()){
        for(int i = 0; checkFile >> dum1 && i < readMAX; i += 1) {

            col_1[i] = dum1;
            col_2[i] = dum2;

            if(i == readMAX){ //Checks if the file is longer than 1000 lines. If it is it is too long for this code to read it.
                std::cout << "Exit Warning: File length exceeds readMAX temp.c line 8\n" << std::endl;
                checkFile.close();
                return i;
            }
        }
        checkFile.close();
        return i;
    }
    else {
        return 0;
    }
}

int append_array(int iMAX, double *array1, double *array2, int i) { //Writes the array into a file
    std::string FileNameAppend;

    std::cout << "File to append: ";
    std::cin >> FileNameAppend; //user input
    std::cout << "Appending to " << FileNameAppend << std::endl;

    std::ofstream outputFile(FileNameAppend, std::ios::app); // Creates or opens FileName.txt
    //cout << "iMAX:" << iMAX << endl;
    for(int count = 0; count < iMAX; count += 1) {
        //cout << "count before cycle:" << count << endl;
        //array1[count+i] = 2.7*(count+i);
        outputFile << array1[count+i] << " " << array2[count+i] << std::endl; //Writes count and array1 at position count to the file
        //cout << count+ii << " " << array1[count+ii] << endl;
        //cout << "count after cycle:" << count << endl;
    }

    outputFile.close(); //Close the file to prevent memory leaks

    return 1; //Returns 1 to show that the function ran and completed
}

int write_array_4c(int iMAX, double *array1, double *array2, double *array3, double *array4) { //Writes the array into a file
    std::string FileNameWrite;

    std::cout << "File to write: ";
    std::cin >> FileNameWrite; //user input
    std::cout << "Writing to " << FileNameWrite << std::endl;

    std::ofstream outputFile(FileNameWrite); // Creates FileName.txt
    //cout << "iMAX:" << iMAX << endl;
    for(int count = 0; count < iMAX; count += 1) {
        //cout << "count before cycle:" << count << endl;
        outputFile << array1[count] << " " << array2[count] << " " << array3[count] << " " << array4[count] << std::endl; //Writes array and array2 at position count to the file
        //cout << count+ii << " " << array1[count+ii] << endl;
        //cout << "count after cycle:" << count << endl;
    }

    outputFile.close(); //Close the file to prevent memory leaks

    return 1; //Returns 1 to show that the function ran and completed
}

int read_array(int readMAX, double *col_1, double *col_2) {
    std::string FileNameRead;

    std::cout << "File to read: ";
    std::cin >> FileNameRead; //user input
    std::cout << "Reading " << FileNameRead << std::endl;

    std::ifstream inputFile(FileNameRead); //Opens FileName.txt in order to read from it

    std::string line;

    //Have to be careful with readMAX, if it is too long it will read out a bunch of -nan values at the end
    //However, if the commented out condition is included if there are any -nan values the reading will stop
    for(int j = 0; std::getline(inputFile,line) && j < readMAX; j += 1) {
        std::istringstream iss(line);
        double dum1 = NAN;
        double dum2 = NAN;

        if (!(iss >> dum1)) dum1 = NAN;
        if (!(iss >> dum2)) dum2 = NAN;

        col_1[j] = dum1;
        col_2[j] = dum2;
        
        std::cout << col_1[j] << " " << col_2[j] << std::endl;
    }

    inputFile.close();

    return 2;
}
