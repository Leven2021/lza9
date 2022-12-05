#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <mpi.h>
#include <omp.h>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"


int size; // problem size
int iteration;


int my_rank;
int world_size;

int n_omp_threads;


void initialize(float *data) {
    // intialize the temperature distribution
    int len = size * size;
    for (int i = 0; i < len; i++) {
        data[i] = wall_temp;
    }
}


void generate_fire_area(bool *fire_area){
    // generate the fire area
    int len = size * size;
    for (int i = 0; i < len; i++) {
        fire_area[i] = 0;
    }

    float fire1_r2 = fire_size * fire_size;
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            int a = i - size / 2;
            int b = j - size / 2;
            int r2 = 0.5 * a * a + 0.8 * b * b - 0.5 * a * b;
            if (r2 < fire1_r2) fire_area[i * size + j] = 1;
        }
    }

    float fire2_r2 = (fire_size / 2) * (fire_size / 2);
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            int a = i - 1 * size / 3;
            int b = j - 1 * size / 3;
            int r2 = a * a + b * b;
            if (r2 < fire2_r2) fire_area[i * size + j] = 1;
        }
    }
}


void update(float *data, float *new_data, int idx) {
    // TODO: update temperature for each point  (in parallelized way)
    if (idx < size * size) {
        int i = idx / size;
        if (i == 0 || i == size - 1) return;
        int j = idx % size;
        if (j == 0 || j == size - 1) return;

        float up = data[idx - size];
        float down = data[idx + size];
        float left = data[idx - 1];
        float right = data[idx + 1];
        float new_val = (up + down + left + right) / 4;
        new_data[idx] = new_val;
    }
}


void maintain_fire(float *data, bool *fire_area, int idx) {
    // TODO: maintain the temperature of the fire (in parallelized way)
    if (idx < size * size) {
        if (fire_area[idx]) data[idx] = fire_temp;
    }
}


void maintain_wall(float *data, int begin, int end) {
    // TODO: maintain the temperature of the wall
    for (int i = begin; i < end; i++) {
        data[i] = wall_temp;
        data[i + size * (size - 1)] = wall_temp;
    }
    for (int i = begin; i < end; i++) {
        data[size * i] = wall_temp;
        data[size * i + size - 1] = wall_temp;
    }
}


#ifdef GUI
void data2pixels(float *data, GLubyte* pixels, int begin, int end){
    // convert rawdata (large, size^2) to pixels (small, resolution^2) for faster rendering speed
    float factor_data_pixel = (float) size / resolution;
    float factor_temp_color = (float) 255 / fire_temp;

    omp_set_num_threads(n_omp_threads);
    #pragma omp parallel for
    for (int x = begin; x < end; x++){
        for (int y = 0; y < resolution; y++){
            int idx = x * resolution + y;
            int idx_pixel = idx * 3;
            int x_raw = x * factor_data_pixel;
            int y_raw = y * factor_data_pixel;
            int idx_raw = x_raw * size + y_raw;
            float temp = data[idx_raw];
            int color =  ((int) temp / 5 * 5) * factor_temp_color;
            pixels[idx_pixel] = color;
            pixels[idx_pixel + 1] = 255 - color;
            pixels[idx_pixel + 2] = 255 - color;
        }
    }
}

void plot(GLubyte* pixels){
    // visualize temprature distribution
    #ifdef GUI
    glClear(GL_COLOR_BUFFER_BIT);
    glDrawPixels(resolution, resolution, GL_RGB, GL_UNSIGNED_BYTE, pixels);
    glutSwapBuffers();
    #endif
}
#endif



void slave(){
    // TODO: MPI routine (one possible solution, you can use another partition method)
    int my_begin_row_id = size * my_rank / (world_size);
    int my_end_row_id = size * (my_rank + 1) / world_size;

    // TODO: Initialize a storage for temperature distribution of this area
    float* data_odd_local = new float[size * size];
    float* data_even_local = new float[size * size];
    bool* fire_area = new bool[size * size];
    // TODO: Receive initial temperature distribution of this area from master
    MPI_Recv(data_odd_local, size * size, MPI_FLOAT, 0, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(fire_area, size * size, MPI_CXX_BOOL, 0, 77, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // TODO: Initialize a storage for local pixels (pls refer to sequential version for initialization of GLubyte)
    #ifdef GUI
    GLubyte* local_pixels = new GLubyte[resolution * resolution * 3];
    #endif

    bool cont = true;
    int terminator;
    int count = 1;
    while (cont) {
        MPI_Recv(&terminator, 1, MPI_INT, 0, 777, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (terminator == 1) cont = false;
        // TODO: computation part
        if (count % 2 == 1) {
            omp_set_num_threads(n_omp_threads);
            #pragma omp parallel for
            for (int i = size * my_begin_row_id; i < size * my_end_row_id; i++) {
                update(data_odd_local, data_even_local, i);
            }

            omp_set_num_threads(n_omp_threads);
            #pragma omp parallel for
            for (int i = size * my_begin_row_id; i < size * my_end_row_id; i++) {
                maintain_fire(data_even_local, fire_area, i);
            }
            
            maintain_wall(data_even_local, my_begin_row_id, my_end_row_id);
        } else {
            omp_set_num_threads(n_omp_threads);
            #pragma omp parallel for
            for (int i = size * my_begin_row_id; i < size * my_end_row_id; i++) {
                update(data_even_local, data_odd_local, i);
            }

            omp_set_num_threads(n_omp_threads);
            #pragma omp parallel for
            for (int i = size * my_begin_row_id; i < size * my_end_row_id; i++) {
                maintain_fire(data_odd_local, fire_area, i);
            }
            
            maintain_wall(data_odd_local, my_begin_row_id, my_end_row_id);
        }
        
        // TODO: after computation, send border row to neighbours
        float* base = (count % 2 == 1) ? data_even_local : data_odd_local;
        if (my_rank % 2 == 1) {
            MPI_Recv(base + size * (my_begin_row_id - 1), size, MPI_FLOAT, my_rank - 1, 7777, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(base + size * my_begin_row_id, size, MPI_FLOAT, my_rank - 1, 77777, MPI_COMM_WORLD);
            if (my_rank != world_size - 1) {
                MPI_Send(base + size * (my_end_row_id - 1), size, MPI_FLOAT, my_rank + 1, 7777, MPI_COMM_WORLD);
                MPI_Recv(base + size * my_end_row_id, size, MPI_FLOAT, my_rank + 1, 77777, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        } else {
            if (my_rank != world_size - 1) {
                MPI_Send(base + size * (my_end_row_id - 1), size, MPI_FLOAT, my_rank + 1, 7777, MPI_COMM_WORLD);
                MPI_Recv(base + size * my_end_row_id, size, MPI_FLOAT, my_rank + 1, 77777, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            MPI_Recv(base + size * (my_begin_row_id - 1), size, MPI_FLOAT, my_rank - 1, 7777, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(base + size * my_begin_row_id, size, MPI_FLOAT, my_rank - 1, 77777, MPI_COMM_WORLD);
        }
        count++;

        #ifdef GUI
        // TODO: conver raw temperature to pixels (much smaller than raw data)
        int p_begin = resolution * my_rank / (world_size);
        int p_end = resolution * (my_rank + 1) / world_size;
        data2pixels(base, local_pixels, p_begin, p_end);
        // TODO: send pixels to master (you can use MPI_Byte to transfer anything to master, then you won't need to declare MPI Type :-) )
        MPI_Send(local_pixels + p_begin * resolution * 3, (p_end - p_begin) * resolution * 3, MPI_BYTE, 0, 777777, MPI_COMM_WORLD);
        #endif
    }

    // TODO: Remember to delete[] local_data and local_pixels.
    delete[] data_odd_local;
    delete[] data_even_local;
    delete[] fire_area;
    
    #ifdef GUI
    delete[] local_pixels;
    #endif
}


void master() {
    // TODO: MPI routine (one possible solution, you can use another partition method)
    int my_begin_row_id = 0;
    int my_end_row_id = size / world_size;

    float* data_odd = new float[size * size];
    float* data_even = new float[size * size];
    bool* fire_area = new bool[size * size];

    initialize(data_odd);
    generate_fire_area(fire_area);

    #ifdef GUI
    GLubyte* pixels;
    pixels = new GLubyte[resolution * resolution * 3];
    #endif

    int count = 1;
    double total_time = 0;

    // TODO: Send initial distribution to each slave process
    for (int i = 1; i < world_size; i++) {
        MPI_Send(data_odd, size * size, MPI_FLOAT, i, 7, MPI_COMM_WORLD);
        MPI_Send(fire_area, size * size, MPI_CXX_BOOL, i, 77, MPI_COMM_WORLD);
    }

    int terminator = 0;
    while (count <= iteration) {
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        if (count == iteration) terminator = 1;
        for (int i = 1; i < world_size; i++) {
            MPI_Send(&terminator, 1, MPI_INT, i, 777, MPI_COMM_WORLD);
        }
        // TODO: Computation of my part
        if (count % 2 == 1) {
            for (int i = size * my_begin_row_id; i < size * my_end_row_id; i++) {
                update(data_odd, data_even, i);
            }
            for (int i = size * my_begin_row_id; i < size * my_end_row_id; i++) {
                maintain_fire(data_even, fire_area, i);
            }
            maintain_wall(data_even, my_begin_row_id, my_end_row_id);
        } else {
            for (int i = size * my_begin_row_id; i < size * my_end_row_id; i++) {
                update(data_even, data_odd, i);
            }
            for (int i = size * my_begin_row_id; i < size * my_end_row_id; i++) {
                maintain_fire(data_odd, fire_area, i);
            }
            maintain_wall(data_odd, my_begin_row_id, my_end_row_id);
        }

        // TODO: Send border row to neighbours
        float* base = (count % 2 == 1) ? data_even : data_odd;
        MPI_Send(base + size * (my_end_row_id - 1), size, MPI_FLOAT, 1, 7777, MPI_COMM_WORLD);
        MPI_Recv(base + size * my_end_row_id, size, MPI_FLOAT, 1, 77777, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        double this_time = std::chrono::duration<double>(t2 - t1).count();
        total_time += this_time;
        printf("Iteration %d, elapsed time: %.6f\n", count, this_time);
        count++;

        #ifdef GUI
        int p_begin = 0;
        int p_end = resolution / world_size;
        data2pixels(base, pixels, p_begin, p_end);
        for (int i = 1; i < world_size; i++) {
            int begin_row = resolution * i / (world_size);
            int end_row = resolution * (i + 1) / world_size;
            int n = (end_row - begin_row) * resolution * 3;
            MPI_Recv(pixels + begin_row * resolution * 3, n, MPI_BYTE, i, 777777, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        plot(pixels);
        #endif
    }

    printf("Stop after %d iterations, elapsed time: %.6f, average computation time: %.6f\n", count-1, total_time, (double) total_time / (count-1));

    delete[] data_odd;
    delete[] data_even;
    delete[] fire_area;

    #ifdef GUI
    delete[] pixels;
    #endif
}


int main(int argc, char *argv[]) {
    size = atoi(argv[1]);
    iteration = atoi(argv[2]);
    n_omp_threads = atoi(argv[3]);

	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);


	if (my_rank == 0) {
        #ifdef GUI
        glutInit(&argc, argv);
        glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
        glutInitWindowPosition(0, 0);
        glutInitWindowSize(resolution, resolution);
        glutCreateWindow("Heat Distribution Simulation MPI & OpenMP Implementation");
        gluOrtho2D(0, resolution, 0, resolution);
        #endif
        master();
	} else {
        slave();
    }

	if (my_rank == 0){
		printf("Student ID: 119010211\n"); // replace it with your student id
		printf("Name: Ziang Liu\n"); // replace it with your name
		printf("Assignment 4: Heat Distribution Simulation MPI & OpenMP Implementation\n");
	}

	MPI_Finalize();

	return 0;
}

