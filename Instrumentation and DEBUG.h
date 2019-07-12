#pragma once
/* Instrumentation and DEBUG */
// Contains macro for printing debug log system and timer

// External modules required
#include <iostream>
#include <chrono>

class Timer
{
public:
	const char* stackframeName;		// Name of stackframe for timer

	// Assign start and end variables to memory with weird type
	std::chrono::time_point<std::chrono::steady_clock> start, end;
	std::chrono::duration<float> duration;

	// Constructor (start timer here)
	Timer(const char* stackframe) : stackframeName(stackframe)
	{
		start = std::chrono::high_resolution_clock::now();
	}

	// Destructor (end timer here) 
	~Timer()
	{
		// end the timer
		end = std::chrono::high_resolution_clock::now();

		// Print time
		duration = end - start;
		std::cout << "Time taken for " << stackframeName << ": " << duration.count() << " s\n";
	}
};


// Print debug log system or not. 
#define DEBUG 0
#if		DEBUG==1
#define d_print(X) std::cout << X << std::endl; 
#elif	DEBUG==0
#define d_print(X) 
#endif

