gnu: 
	g++ Bruteforce.cpp Forces.cpp -Ofast -fopenmp -march=native -mtune=native && time ./a.out

icc: 
	icpx -xhost Bruteforce.cpp Forces.cpp -qopenmp -qopt-report && time ./a.out

gnuBlock: 
	g++ BruteforceBlocked.cpp  Forces.cpp -Ofast -fopenmp -march=native -mtune=native && time ./a.out

iccBlock: 
	icpx -xhost BruteforceBlocked.cpp  Forces.cpp -qopenmp -qopt-report && time ./a.out


all: gnu icc