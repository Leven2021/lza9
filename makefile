seq:
	g++ ./src/sequential.cpp -o seq -O2 -std=c++11 -g
mpi:	
	mpic++ ./src/mpi.cpp -o mpi -std=c++11
pthread:
	g++ ./src/pthread.cpp -o pthread -lpthread -O2 -std=c++11
cuda:
	nvcc ./src/cuda.cu -o cuda -O2 --std=c++11
seqg:
	g++ ./src/sequential.cpp -o seqg -I/usr/include -L/usr/local/lib -L/usr/lib -lglut -lGLU -lGL -lm -DGUI -O2 -std=c++11
mpig:
	mpic++ ./src/mpi.cpp -o mpig -I/usr/include -L/usr/local/lib -L/usr/lib -lglut -lGLU -lGL -lm -DGUI -std=c++11
pthreadg:
	g++ ./src/pthread.cpp -o pthreadg -I/usr/include -L/usr/local/lib -L/usr/lib -lglut -lGLU -lGL -lm -lpthread -DGUI -O2 -std=c++11
/home/gpgpu-sim/CSC4005_2022Fall_Demo/cuda_simulator/GTX480_rundircudag:
	nvcc ./src/cuda.cu -o cudag -I/usr/include -L/usr/local/lib -L/usr/lib -lglut -lGLU -lGL -lm -O2 -DGUI --std=c++11
openmp:
	g++ ./src/openmp.cpp -o openmp -fopenmp -O2 -std=c++11
openmpg:
	g++ ./src/openmp.cpp -o openmpg -fopenmp -I/usr/include -L/usr/local/lib -L/usr/lib -lglut -lGLU -lGL -lm -O2 -DGUI -std=c++11
mpi_omp:
	mpic++ ./src/mpi_omp.cpp -o mpi_omp -fopenmp -std=c++11
mpi_ompg:
	mpic++ ./src/mpi_omp.cpp -o mpi_ompg -fopenmp -I/usr/include -L/usr/local/lib -L/usr/lib -lglut -lGLU -lGL -lm -O2 -DGUI -std=c++11
all:
	make seq
	make mpi
	make pthread
	make cuda
	make seqg
	make mpig
	make pthreadg
	make cudag
	make openmp
	make openmpg
	make mpi_omp
	make mpi_ompg
a4:
	make seq
	make mpi
	make pthread
	make cuda
	make openmp
	make mpi_omp
a4_gui:
	make seqg
	make mpig
	make pthreadg
	make openmpg
	make mpi_ompg
clean:
	rm -f seq mpi pthread seqg mpig pthreadg cuda cudag openmp openmpg mpi_omp mpi_ompg
