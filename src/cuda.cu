#include <cuda.h>
#include <cuda_runtime.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// #include <chrono>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"


int block_size = 512; // cuda thread block size
int size; // problem size
int iteration;

__device__ int dsize = 10;


__global__ void initialize(float *data) {
    // TODO: intialize the temperature distribution (in parallelized way)
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    printf("dsize is %d\n", dsize);
    if (i < dsize * dsize) {
        data[i] = wall_temp;
    }
}


__global__ void generate_fire_area(bool *fire_area){
    // TODO: generate the fire area (in parallelized way)
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx < dsize * dsize) {
        int i = idx / dsize;
        int j = idx % dsize;
        fire_area[idx] = 0;

        float fire1_r2 = fire_size * fire_size;
        int a = i - dsize / 2;
        int b = j - dsize / 2;
        int r2 = 0.5 * a * a + 0.8 * b * b - 0.5 * a * b;
        if (r2 < fire1_r2) fire_area[i * dsize + j] = 1;

        float fire2_r2 = (fire_size / 2) * (fire_size / 2);
        a = i - 1 * dsize / 3;
        b = j - 1 * dsize / 3;
        r2 = a * a + b * b;
        if (r2 < fire2_r2) fire_area[i * dsize + j] = 1;
    }
}


__global__ void update(float *data, float *new_data) {
    // TODO: update temperature for each point  (in parallelized way)
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx < dsize * dsize) {
        int i = idx / dsize;
        if (i == 0 || i == dsize - 1) return;
        int j = idx % dsize;
        if (j == 0 || j == dsize - 1) return;

        float up = data[idx - dsize];
        float down = data[idx + dsize];
        float left = data[idx - 1];
        float right = data[idx + 1];
        float new_val = (up + down + left + right) / 4;
        new_data[idx] = new_val;
    }
}


__global__ void maintain_wall(float *data) {
    // TODO: maintain the temperature of the wall (sequential is enough)
    for (int i = 0; i < size; i++) {
        data[i] = wall_temp;
    }
    for (int i = size * (size - 1); i < size * size; i++) {
        data[i] = wall_temp;
    }
    for (int i = 0; i < size; i++) {
        data[size * i] = wall_temp;
        data[size * i + size - 1] = wall_temp;
    }
}


__global__ void maintain_fire(float *data, bool *fire_area) {
    // TODO: maintain the temperature of the fire (in parallelized way)
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < dsize * dsize) {
        if (fire_area[i]) data[i] = fire_temp;
    }
}


#ifdef GUI
__global__ void data2pixels(float *data, GLubyte* pixels){
    // TODO: convert rawdata (large, size^2) to pixels (small, resolution^2) for faster rendering speed (in parallelized way)
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx < resolution * resolution) {
        float factor_data_pixel = (float) dsize / resolution;
        float factor_temp_color = (float) 255 / fire_temp;

        int x = idx / resolution;
        int y = idx % resolution;

        int idx_pixel = idx * 3;
        int x_raw = x * factor_data_pixel;
        int y_raw = y * factor_data_pixel;
        int idx_raw = x_raw * dsize + y_raw;
        float temp = data[idx_raw];
        int color =  ((int) temp / 5 * 5) * factor_temp_color;
        pixels[idx_pixel] = color;
        pixels[idx_pixel + 1] = 255 - color;
        pixels[idx_pixel + 2] = 255 - color;
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


void master() {
    float *data_odd;
    float *data_even;
    bool *fire_area;

    cudaMalloc(&data_odd, size * size * sizeof(float));
    cudaMalloc(&data_even, size * size * sizeof(float));
    cudaMalloc(&fire_area, size * size * sizeof(bool));

    #ifdef GUI
    GLubyte *pixels;
    GLubyte *host_pixels;
    host_pixels = new GLubyte[resolution * resolution * 3];
    cudaMalloc(&pixels, resolution * resolution * 3 * sizeof(GLubyte));
    #endif

    int n_block_size = size * size / block_size + 1;
    int n_block_resolution = resolution * resolution / block_size + 1;

    initialize<<<n_block_size, block_size>>>(data_odd);
    generate_fire_area<<<n_block_size, block_size>>>(fire_area);
    
    int count = 1;
    double total_time = 0;

    while (count <= iteration){
        // std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        // TODO: modify the following lines to fit your need.
        if (count % 2 == 1) {
            update<<<n_block_size, block_size>>>(data_odd, data_even);
            maintain_fire<<<n_block_size, block_size>>>(data_even, fire_area);
            maintain_wall<<<1, 1>>>(data_even);
        } else {
            update<<<n_block_size, block_size>>>(data_even, data_odd);
            maintain_fire<<<n_block_size, block_size>>>(data_odd, fire_area);
            maintain_wall<<<1, 1>>>(data_odd);
        }

        // std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        // double this_time = std::chrono::duration<double>(t2 - t1).count();
        // total_time += this_time;
        // printf("Iteration %d, elapsed time: %.6f\n", count, this_time);
        count++;
        
        #ifdef GUI
        if (count % 2 == 1) {
            data2pixels<<<n_block_resolution, block_size>>>(data_even, pixels);
        } else {
            data2pixels<<<n_block_resolution, block_size>>>(data_odd, pixels);
        }
        cudaMemcpy(host_pixels, pixels, resolution * resolution * 3 * sizeof(GLubyte), cudaMemcpyDeviceToHost);
        plot(host_pixels);
        #endif

    }

    // printf("Stop after %d iterations, elapsed time: %.6f, average computation time: %.6f\n", count-1, total_time, (double) total_time / (count-1));


    cudaFree(data_odd);
    cudaFree(data_even);
    cudaFree(fire_area);

    #ifdef GUI
    cudaFree(pixels);
    delete[] host_pixels;
    #endif
    
}


int main(int argc, char *argv[]){
    
    size = atoi(argv[1]);
    iteration = atoi(argv[2]);

    int* temp_h = &size;
    int* temp_d;
    cudaMalloc(&temp_d, sizeof(int));
    cudaMemcpy(temp_d, temp_h, sizeof(int), cudaMemcpyHostToDevice);
    ini_size<<<1, 1>>>(temp_d);

    #ifdef GUI
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(resolution, resolution);
    glutCreateWindow("Heat Distribution Simulation CUDA Implementation");
    gluOrtho2D(0, resolution, 0, resolution);
    #endif

    master();

    printf("Student ID: 119010211\n"); // replace it with your student id
    printf("Name: Ziang Liu\n"); // replace it with your name
    printf("Assignment 4: Heat Distribution CUDA Implementation\n");

    return 0;

}


