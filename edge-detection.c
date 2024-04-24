#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include <float.h>
#include <stdbool.h>

// Global parameters/variables
#define t_grayscale unsigned char
#define GRAYSCALE_MAX UCHAR_MAX
#define GRAYSCALE_MIN 0
#define MPI_T_GRAYSCALE MPI_UNSIGNED_CHAR
#define N_GHOST_PER_SIZE 2  // (int) max(kernel_size1, kernel_size2, ...) / 2
#define ORIENTATION 0 // 0: horizontal, 1: vertical
#define HORIZONTAL 0
#define VERTICAL 1
#define M_PI 3.14159265358979323846 // Pi, the support of it is varies among compilers, so let's define it here.

// Gaussian filter parameters
#define GAUSSIAN_STD 1.5
#define GAUSSIAN_MEAN 0
#define GAUSSIAN_KERNEL_SIZE 5

// Sobel operator parameters
#define SOBEL_KERNEL_SIZE 3
typedef struct {
    t_grayscale magnitude;
    float angle;
} gradient_image; // the output of sobel

// Macros
#define I(N, P, p)  ((N+P-p-1)/P)
#define I_INVERSE(N, P, p, i) (p*(N/P) + MIN(N%P, p) + i)
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define LOCAL_CELL_COUNT(P, p, axis_main, axis_secondary) (I(axis_secondary, P, p) * axis_main) // Number of local cell (excluding ghost cells)
#define GHOST_CELL_COUNT(P, p, axis_main, axis_secondary) (N_GHOST_PER_SIZE * axis_main) // Number of ghost cell
#define ALL_CELL_COUNT(P, p, axis_main, axis_secondary) (LOCAL_CELL_COUNT(P, p, axis_main, axis_secondary) + 2 * GHOST_CELL_COUNT(P, p, axis_main, axis_secondary)) // Number of all cells
#define LOCAL_CELL_SIZE(P, p, axis_main, axis_secondary) (LOCAL_CELL_COUNT(P, p, axis_main, axis_secondary) * sizeof(t_grayscale)) // The size of the local cell
#define GHOST_CELL_SIZE(P, p, axis_main, axis_secondary) (GHOST_CELL_COUNT(P, p, axis_main, axis_secondary) * sizeof(t_grayscale)) // The size of the ghost cell
#define ALL_CELL_SIZE(P, p, axis_main, axis_secondary) (ALL_CELL_COUNT(P, p, axis_main, axis_secondary) * sizeof(t_grayscale)) // The size of all cells
#define LOCAL_CELL_OFFSET(P, p, axis_main, axis_secondary) (N_GHOST_PER_SIZE * axis_main) // The index of the first local cell
#define GHOST_CELL_HEAD_OFFSET(P, p, axis_main, axis_secondary) (0) // The index of the first ghost cell
#define GHOST_CELL_TAIL_OFFSET(P, p, axis_main, axis_secondary) (LOCAL_CELL_OFFSET(P, p, axis_main, axis_secondary) + LOCAL_CELL_COUNT(P, p, axis_main, axis_secondary)) // The index of the first ghost cell
#define IMAGE_IDX(i, j, axis_main) ((i)*(axis_main)+(j)) // The 


void error(const char *msg);
void communication(t_grayscale *ptr_cells, int P, int p, int axis_main, int axis_secondary);
double gaussian_distribution_2D(double x, double y, double mean, double std);
void gaussian_blur(t_grayscale *ptr_cells, int P, int p, int axis_main, int axis_secondary);
void sobel_operation(t_grayscale *ptr_cells, gradient_image *output, int P, int p, int axis_main, int axis_secondary);
void debug_struct_array(t_grayscale *output, gradient_image *input, int count);

int main(int argc, char *argv[]) {
    if (argc <= 2) {
        error("Too less argument provided. Usage: <input_file> <output_file>");
    }

    // MPI initialization
    int p, P, I;
    char *input_file = argv[1];
    char *output_file = argv[2];
    int axis_main, axis_secondary;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    
    // Read input file
    //** Read the width and height of the image
    MPI_File in_fh;
    MPI_Offset offset = 0;
    MPI_File_open(MPI_COMM_WORLD, input_file, MPI_MODE_RDONLY, MPI_INFO_NULL, &in_fh);
    if (in_fh == NULL) {
        error("Error opening input file.");
    }
    if (ORIENTATION == HORIZONTAL) {
        //** Read the width of the image.
        MPI_File_read_at_all(in_fh, offset, &axis_main, 1, MPI_INT32_T, MPI_STATUS_IGNORE);
        //** Read the height of the image.
        MPI_File_read_at_all(in_fh, offset += 4, &axis_secondary, 1, MPI_INT32_T, MPI_STATUS_IGNORE);
    } else {
        //** Read the height of the image.
        MPI_File_read_at_all(in_fh, offset, &axis_secondary, 1, MPI_INT32_T, MPI_STATUS_IGNORE);
        //** Read the width of the image.
        MPI_File_read_at_all(in_fh, offset += 4, &axis_main, 1, MPI_INT32_T, MPI_STATUS_IGNORE);
    }

    //** Read the image data
    //**** Represent the 2D array in 1D array
    t_grayscale *local_image = (t_grayscale *) malloc(ALL_CELL_SIZE(P, p, axis_main, axis_secondary));
    I = I(axis_secondary, P, p);
    // magnitude and angle struct after sobel operation
    gradient_image *local_image_gradient = (gradient_image *) malloc(ALL_CELL_COUNT(P, p, axis_main, axis_secondary) * sizeof(gradient_image));
    offset += I_INVERSE(axis_secondary, P, p, 0) * axis_main * sizeof(t_grayscale) + 4;
    MPI_File_read_at(in_fh, offset, local_image + axis_main * N_GHOST_PER_SIZE, LOCAL_CELL_COUNT(P, p, axis_main, axis_secondary), MPI_T_GRAYSCALE, MPI_STATUS_IGNORE);
    MPI_File_close(&in_fh);



    /*
    
    Canny Edge Detection Algorithm (Main Steps)
    
    */
    communication(local_image, P, p, axis_main, axis_secondary);
    gaussian_blur(local_image, P, p, axis_main, axis_secondary);
    
    communication(local_image, P, p, axis_main, axis_secondary);
    sobel_operation(local_image, local_image_gradient, P, p, axis_main, axis_secondary);
    
    t_grayscale *debug_out = (t_grayscale *) malloc(ALL_CELL_SIZE(P, p, axis_main, axis_secondary));
    debug_struct_array(debug_out, local_image_gradient, ALL_CELL_COUNT(P, p, axis_main, axis_secondary));

    // Write output file
    MPI_File out_fh;
    MPI_File_open(MPI_COMM_WORLD, output_file, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &out_fh);
    if(p == 0){
        if (ORIENTATION == HORIZONTAL) {
            MPI_File_write_at(out_fh, 0, &axis_main, 1, MPI_INT32_T, MPI_STATUS_IGNORE);
            MPI_File_write_at(out_fh, 4, &axis_secondary, 1, MPI_INT32_T, MPI_STATUS_IGNORE);
        } else {
            MPI_File_write_at(out_fh, 0, &axis_secondary, 1, MPI_INT32_T, MPI_STATUS_IGNORE);
            MPI_File_write_at(out_fh, 4, &axis_main, 1, MPI_INT32_T, MPI_STATUS_IGNORE);
        }
    }
    // MPI_File_seek(out_fh, offset, MPI_SEEK_SET);
    MPI_File_write_at(out_fh, offset, debug_out + LOCAL_CELL_OFFSET(P, p, axis_main, axis_secondary), LOCAL_CELL_COUNT(P, p, axis_main, axis_secondary), MPI_T_GRAYSCALE, MPI_STATUS_IGNORE);
    MPI_File_close(&out_fh);

    free(local_image);
    MPI_Finalize();
    return 0;
}

/**
 * @brief   Communication between ranks.
 *          The first rank will flip its local cells at head to its ghost cells at head.
 *          The last rank will flip its local cells at tail to its ghost cells at tail.
 *          
 * @param   ptr_cells   The pointer to the start of the local cells.
 * @param   P           Number of total ranks.
 * @param   p           Rank number.
 * @param   axis_main   The size of the main axis.
 * @param   axis_secondary The size of the secondary axis.
 */
void communication(t_grayscale *ptr_cells, int P, int p, int axis_main, int axis_secondary){
    
    t_grayscale *ptr_local_cell_to_send_head = ptr_cells + LOCAL_CELL_OFFSET(P, p, axis_main, axis_secondary);
    t_grayscale *ptr_local_cell_to_send_tail = ptr_cells + GHOST_CELL_TAIL_OFFSET(P, p, axis_main, axis_secondary) - GHOST_CELL_COUNT(P, p, axis_main, axis_secondary);
    t_grayscale *ptr_ghost_head = ptr_cells;
    t_grayscale *ptr_ghost_tail = ptr_cells + GHOST_CELL_TAIL_OFFSET(P, p, axis_main, axis_secondary);
    int ghost_cell_COUNT = GHOST_CELL_COUNT(P, p, axis_main, axis_secondary);
    // Rank 2n communicate with rank 2n+1
    if(p % 2 == 0){
        // Send then receive, to/from next rank
        if(p < P - 1){
            MPI_Send(ptr_local_cell_to_send_tail, ghost_cell_COUNT, MPI_T_GRAYSCALE, p + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(ptr_ghost_tail, ghost_cell_COUNT, MPI_T_GRAYSCALE, p + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else{
            // If last rank, flip alone the boarder locally
            for (int i = 0; i < N_GHOST_PER_SIZE; i++)
            {   
                memcpy(ptr_ghost_tail + i * axis_main, ptr_ghost_tail - (i + 2) * axis_main, axis_main);
            }
        }
    }
    else{
        // Receive then send, from/to previous rank (never be the first rank).
        MPI_Recv(ptr_ghost_head, ghost_cell_COUNT, MPI_T_GRAYSCALE, p - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(ptr_local_cell_to_send_head, ghost_cell_COUNT, MPI_T_GRAYSCALE, p - 1, 0, MPI_COMM_WORLD);
    }
    // Rank 2n communicate with rank 2n-1
    if(p % 2 == 0){
        // Send then receive, to/from previous rank
        if(p > 0){
            MPI_Send(ptr_local_cell_to_send_head, ghost_cell_COUNT, MPI_T_GRAYSCALE, p - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(ptr_ghost_head, ghost_cell_COUNT, MPI_T_GRAYSCALE, p - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else{
            // If first rank, flip alone the boarder locally
            for (int i = 0; i < N_GHOST_PER_SIZE; i++)
            {   
                memcpy(ptr_local_cell_to_send_head - (i + 1) * axis_main, ptr_local_cell_to_send_head + (i + 1) * axis_main, axis_main);
            }
        }
    }
    else{
        if(p < P - 1){
            MPI_Recv(ptr_ghost_tail, ghost_cell_COUNT, MPI_T_GRAYSCALE, p + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(ptr_local_cell_to_send_tail, ghost_cell_COUNT, MPI_T_GRAYSCALE, p + 1, 0, MPI_COMM_WORLD);
        }
        else{
            // If last rank, flip alone the boarder locally
            for (int i = 0; i < N_GHOST_PER_SIZE; i++)
            {   
                memcpy(ptr_ghost_tail + i * axis_main, ptr_ghost_tail - (i + 2) * axis_main, axis_main);
            }
        }
    }
}

/**
 * @brief   Return the G(x,y), where G is the PDF of 2D Gaussian distribution.
 * 
 * @param   x           Value of variable x.
 * @param   y           Value of variable y.
 * @param   mean        Mean of G.
 * @param   std         Standard deviation of G.
 * @return double       G(x,y).
 */
double gaussian_distribution_2D(double x, double y, double mean, double std){
    return (1 / (2 * std * std * M_PI)) * exp( - (x*x + y*y) / (2 * std * std) );
}

/**
 * @brief   Gaussian blur the image, in place.
 *          For the pixels on the edge, alone the secondary axis,
 *          they will be flipped to the other side of the edge when applying the filter.
 *          i.e. we have a 3*3 kernel and the image image:
 *          ________________
 *          | 11 12 13 ...
 *          | 21 22 23 ...
 *          | 31 32 33 ...
 *          | ...
 *          
 *          When apply kernel on the pixel 21, it will be:
 *                      |12 11 12|
 *          dot(kernel, |22 21 22|)
 *                      |32 31 32|
 * 
 * @param   ptr_cells   The pointer to the start of the local cells.
 * @param   P           Number of total ranks.
 * @param   p           Rank number.
 * @param   axis_main   The size of the main axis.
 * @param   axis_secondary The size of the secondary axis.
 */
void gaussian_blur(t_grayscale *ptr_cells, int P, int p, int axis_main, int axis_secondary){
    // Gaussian filter
    //** Create a Gaussian kernel
    double kernel[GAUSSIAN_KERNEL_SIZE][GAUSSIAN_KERNEL_SIZE];
    double norm = 0;
    t_grayscale *temp_image = (t_grayscale *) malloc(LOCAL_CELL_SIZE(P, p, axis_main, axis_secondary));
    for (int i = 0; i < GAUSSIAN_KERNEL_SIZE; i++)
    {
        for (int j = 0; j < GAUSSIAN_KERNEL_SIZE; j++)
        {
            kernel[i][j] = gaussian_distribution_2D(i - GAUSSIAN_KERNEL_SIZE / 2, j - GAUSSIAN_KERNEL_SIZE / 2, GAUSSIAN_MEAN, GAUSSIAN_STD);
            norm += kernel[i][j];
        }
    }

    for (int offset = 0; offset < LOCAL_CELL_COUNT(P, p, axis_main, axis_secondary); offset++)
    {   
        double sum = 0;
        for (int i =0; i < GAUSSIAN_KERNEL_SIZE; i++)
        {
            for (int j = 0; j < GAUSSIAN_KERNEL_SIZE; j++)
            {   
                int translated_idx = offset + (i - GAUSSIAN_KERNEL_SIZE / 2) + (j - GAUSSIAN_KERNEL_SIZE / 2) * axis_main;
                // Handel the edges
                if(offset % axis_main < GAUSSIAN_KERNEL_SIZE / 2){
                    translated_idx = offset + abs(i - GAUSSIAN_KERNEL_SIZE / 2) + (j - GAUSSIAN_KERNEL_SIZE / 2) * axis_main;
                }
                else if(offset % axis_main >= axis_main - GAUSSIAN_KERNEL_SIZE / 2){
                    translated_idx = offset - abs(i - GAUSSIAN_KERNEL_SIZE / 2) + (j - GAUSSIAN_KERNEL_SIZE / 2) * axis_main;
                }
                sum += kernel[i][j] * ptr_cells[LOCAL_CELL_OFFSET(P, p, axis_main, axis_secondary) + translated_idx];
            }
        }
        temp_image[offset] = (t_grayscale) (sum / norm);
    }
    memcpy(ptr_cells + LOCAL_CELL_OFFSET(P, p, axis_main, axis_secondary), temp_image, LOCAL_CELL_SIZE(P, p, axis_main, axis_secondary));
    free(temp_image);
}

/**
 * @brief   sobel operation
 * 
 */

void sobel_operation(t_grayscale *ptr_cells, gradient_image *output, int P, int p, int axis_main, int axis_secondary) {
    // sobel kernels
    int Gx[SOBEL_KERNEL_SIZE][SOBEL_KERNEL_SIZE] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
    int Gy[SOBEL_KERNEL_SIZE][SOBEL_KERNEL_SIZE] = {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};

    // the radius is only 1 but it doesn't matter if we follows 2.
    gradient_image *temp_gradients = (gradient_image *) malloc(LOCAL_CELL_COUNT(P, p, axis_main, axis_secondary) * sizeof(gradient_image));
    for (int offset = 0; offset < LOCAL_CELL_COUNT(P, p, axis_main, axis_secondary); offset++)
    {   
        double Gx_sum = 0, Gy_sum = 0;
        for (int i =0; i < SOBEL_KERNEL_SIZE; i++)
        {
            for (int j = 0; j < SOBEL_KERNEL_SIZE; j++)
            {   
                int translated_idx = offset + (i - SOBEL_KERNEL_SIZE / 2) + (j - SOBEL_KERNEL_SIZE / 2) * axis_main;
                // Handel the edges
                if(offset % axis_main < SOBEL_KERNEL_SIZE / 2){
                    translated_idx = offset + abs(i - SOBEL_KERNEL_SIZE / 2) + (j - SOBEL_KERNEL_SIZE / 2) * axis_main;
                }
                else if(offset % axis_main >= axis_main - SOBEL_KERNEL_SIZE / 2){
                    translated_idx = offset - abs(i - SOBEL_KERNEL_SIZE / 2) + (j - SOBEL_KERNEL_SIZE / 2) * axis_main;
                }
                Gx_sum += ptr_cells[LOCAL_CELL_OFFSET(P, p, axis_main, axis_secondary) + translated_idx] * Gx[i][j];
                Gy_sum += ptr_cells[LOCAL_CELL_OFFSET(P, p, axis_main, axis_secondary) + translated_idx] * Gy[i][j];
            }
        }
        // calculate the magnitude, with the sqrt method
        double G_magnitude = sqrt(Gx_sum * Gx_sum + Gy_sum * Gy_sum);
        temp_gradients[offset].magnitude = (t_grayscale)MIN(MAX(G_magnitude, GRAYSCALE_MIN), GRAYSCALE_MAX);
        // calculate the angle [rad]
        temp_gradients[offset].angle = atan2(Gy_sum, Gx_sum);
    }
    memcpy(output + LOCAL_CELL_OFFSET(P, p, axis_main, axis_secondary), temp_gradients, LOCAL_CELL_COUNT(P, p, axis_main, axis_secondary) * sizeof(gradient_image));
    free(temp_gradients);
}


/**
 * @brief   Exit the program with an error message.
 * 
 * @param   msg         The error message.
 */
void error(const char *msg) {
    printf("%s\n", msg);
    exit(1);
}

void debug_struct_array(t_grayscale *output, gradient_image *input, int count){
    for (int i = 0; i < count; i++)
    {
        output[i] = input[i].magnitude;
    }
}