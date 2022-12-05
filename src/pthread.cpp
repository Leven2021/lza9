#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <pthread.h>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"

int size;
int iteration;

int n_thd; // number of threads

typedef struct {
    //TODO: specify your arguments for threads
    int id;
    int count;
    float* data_odd;
    float* data_even;
    bool* fire_area;

    #ifdef GUI
    GLubyte* pix;
    #endif
    //TODO END
} Args;

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

void update(float *data, float *new_data, int begin, int end) {
    // TODO: update the temperature of each point, and store the result in `new_data` to avoid data racing
    if (begin == 0) begin = 1;
    if (end == size) end = size - 1;
    for (int i = begin; i < end; i++){
        for (int j = 1; j < size - 1; j++){
            int idx = i * size + j;
            float up = data[idx - size];
            float down = data[idx + size];
            float left = data[idx - 1];
            float right = data[idx + 1];
            float new_val = (up + down + left + right) / 4;
            new_data[idx] = new_val;
        }
    }
}

void maintain_fire(float *data, bool* fire_area, int begin, int end) {
    // TODO: maintain the temperature of fire
    for (int i = begin; i < end; i++){
        for (int j = 0; j < size; j++){
            int idx = i * size + j;
            if (fire_area[idx]) data[idx] = fire_temp;
        }
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


void* worker(void* args) {
    Args* my_arg = (Args*) args;
    int my_id = my_arg->id;
    int my_count = my_arg->count;

    int my_begin_row_id = size * my_id / (n_thd);
    int my_end_row_id = size * (my_id + 1) / n_thd;

    if (my_count % 2 == 1) {
        update(my_arg->data_odd, my_arg->data_even, my_begin_row_id, my_end_row_id);
        maintain_fire(my_arg->data_even, my_arg->fire_area, my_begin_row_id, my_end_row_id);
        maintain_wall(my_arg->data_even, my_begin_row_id, my_end_row_id);
    } else {
        update(my_arg->data_even, my_arg->data_odd, my_begin_row_id, my_end_row_id);
        maintain_fire(my_arg->data_odd, my_arg->fire_area, my_begin_row_id, my_end_row_id);
        maintain_wall(my_arg->data_odd, my_begin_row_id, my_end_row_id);
    }
}

#ifdef GUI
void* render(void* args) {
    Args* my_arg = (Args*) args;
    int my_id = my_arg->id;
    int count = my_arg->count;

    int begin_row = resolution * my_id / (n_thd);
    int end_row = resolution * (my_id + 1) / n_thd;

    if (count % 2 == 1) {
        data2pixels(my_arg->data_even, my_arg->pix, begin_row, end_row);
    } else {
        data2pixels(my_arg->data_odd, my_arg->pix, begin_row, end_row);
    }
}
#endif

void master() {
    float *data_odd;
    float *data_even;
    bool *fire_area;

    #ifdef GUI
    GLubyte* pixels = new GLubyte[resolution * resolution * 3];
    #endif
    
    data_odd = new float[size * size];
    data_even = new float[size * size];
    fire_area = new bool[size * size];

    generate_fire_area(fire_area);
    initialize(data_odd);

    int count = 1;
    double total_time = 0;

    while (count <= iteration) {
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        pthread_t thds[n_thd]; // thread pool
        Args args[n_thd]; // arguments for all threads
        for (int i = 0; i < n_thd; i++) {
            args[i].id = i;
            args[i].count = count;
            args[i].data_even = data_even;
            args[i].data_odd = data_odd;
            args[i].fire_area = fire_area;
        }
        for (int i = 0; i < n_thd; i++) pthread_create(&thds[i], NULL, worker, &args[i]);
        for (int i = 0; i < n_thd; i++) pthread_join(thds[i], NULL);

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        double this_time = std::chrono::duration<double>(t2 - t1).count();
        total_time += this_time;
        printf("Iteration %d, elapsed time: %.6f\n", count, this_time);
        count++;

        #ifdef GUI
        for (int i = 0; i < n_thd; i++) {
            args[i].id = i;
            args[i].count = count;
            args[i].data_even = data_even;
            args[i].data_odd = data_odd;
            args[i].fire_area = fire_area;
            args[i].pix = pixels;
        }
        for (int i = 0; i < n_thd; i++) pthread_create(&thds[i], NULL, render, &args[i]);
        for (int i = 0; i < n_thd; i++) pthread_join(thds[i], NULL);
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
    n_thd = atoi(argv[3]);


    #ifdef GUI
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(resolution, resolution);
    glutCreateWindow("Heat Distribution Simulation Pthread Implementation");
    gluOrtho2D(0, resolution, 0, resolution);
    #endif

    master();

    printf("Student ID: 119010211\n"); // replace it with your student id
    printf("Name: Ziang Liu\n"); // replace it with your name
    printf("Assignment 4: Heat Distribution Pthread Implementation\n");

	return 0;
}