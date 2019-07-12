#pragma once
/* Interface file containing function prototypes of general functions
Use of object to bind all function prototypes together*/

#include <iostream>
#include <vector>

struct Generalfunctions
{
public:
	// Function for printing integer arrays 
	// input_array[] = *input_array as arrays decay to pointers and only pointer is passed to functions. 
	void printIntarray_1D(int* input_array, int size_i)
	{
		// Start of array bracket 
		std::cout << "[";

		// March down i (cols)
		for (int i = 0; i < size_i; i++)
		{
			std::cout << input_array[i] << ",";
		}

		// End bracket 
		std::cout << "]" << std::endl;

	}


	// input_array[] = *input_array as arrays decay to pointers and only pointer is passed to functions. 
	void printIntarray_2D(int size_i, int size_j, int* input_array)
	{
		// Start of array bracket 
		std::cout << "[";



		// End bracket 
		std::cout << "]" << std::endl;

	}	// End function 


	void flux();

};