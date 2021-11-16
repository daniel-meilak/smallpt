smallpt, a Path Tracer by Kevin Beason, with additional annotations by Dr. David Cline. Variables have also been given more readable names. Some C style code has been changed to C++.  

Make : `g++ -O3 -fopenmp smallpt.cpp -o smallpt`  
Remove `-fopenmp` for g++ version < 4.2  
Usage: `time ./smallpt 5000 && xv image.ppm`  