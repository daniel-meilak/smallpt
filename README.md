Playing around with smallpt Global Illumination

smallpt, a Path Tracer by Kevin Beason, 2008  
Make : `g++ -O3 -fopenmp smallpt.cpp -o smallpt`  
Remove `-fopenmp` for g++ version < 4.2  
Usage: `time ./smallpt 5000 && xv image.ppm`  