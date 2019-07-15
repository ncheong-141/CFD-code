# CFD-code
CFD code adapted from a Matlab code. Focused on optimising code for execution time. 

Spatial discretion:		Unstructured 2D grid
Physics model:			  Euler method with Runga Kutta time stepping
Analysis:				      Aerofoil coefficient of lift/drag for the NACA0012

Optimization objectives:

  - Focus on optimising the code for execution time through various techniques such as:
    -> Efficient memory access and data structures. Keep data small and close together in memory. 
       -> Cacheline considerations (indirect write avoidance) 
       Keep memory which is use in a calculation next to each other in memory so it can be loaded 
       on the same cacheline. 
    -> Avoiding unnessecary work calculation/instruction wise. 
    -> Avoiding data dependencies for intruction level parallelism. 
    -> Strength reduction. Re-orginize code to allow for usage of fastest operations (reduce cycle time
       of a calculation.)
    -> Minimizing flow control (less need of the branch predictor and easier time for the prefetcher) 
    -> Temporal cache coherency, ensures calculations are in an order such that a variable it not accessed, 
	     discard, then immediatly accessed again after. (avoids loading cachelines unnessecarily.)  
    -> Prefer 32-bit ints to all other sizes. It is sweet spot for integers as ALU (Arithmetic logic unit)					
	     inside the CPU can do one 64-bit operations at a time and two 32-bit operations at time.  
    -> Prefer floats to doubles (if precision is not an issue) as more floats can be loaded onto the cacheline
       at once. 
    -> Avoid type conversions. (waste of time) 
    -> Use stack memory if data sizes do not change, accessing stack memory is much faster than dynamic usually. 
    -> Prefer inline functions than function calls. Stops CPU requiring to jump to a new stack frame. 

Results: 

  - 204 times the speed of Matlab code, where simply adapting the code from Matlab originally was 116 times faster. 
