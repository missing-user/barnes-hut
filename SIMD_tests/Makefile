gnu: 
	g++ Bruteforce.cpp -Ofast -fopenmp -march=native -mtune=native && time ./a.out

icc: 
	icpx -xhost Bruteforce.cpp -qopenmp -qopt-report && time ./a.out

gnuBlock: 
	g++ BruteforceBlocked.cpp -Ofast -fopenmp -march=native -mtune=native && time ./a.out

iccBlock: 
	icpx -xhost BruteforceBlocked.cpp -qopenmp -qopt-report && time ./a.out


all: gnu icc