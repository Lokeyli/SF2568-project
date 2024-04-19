#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define INPUT_FILE ".txt"
#define I(N, P, p)  ((N+P-p-1)/P)
#define I_INVERSE(N, P, p, i) (p*(N/P) + MIN(N%P, p) + i)
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define CHUNK_BEGIN(prt_chunk, n_ghost, width, t_size) (prt_chunk + n_ghost * width * t_size)
#define CHUNK_SIZE(P, p, width, height) (I(height, P, p) * width)
#define IMAGE_IDX(i, j, width) ((i)*(width)+(j))
#define GAUSSIAN_STD 1.5
#define GAUSSIAN_SIZE 5
#define t_grayscale unsigned char
#define MPI_T_GRAYSCALE MPI_UNSIGNED_CHAR

void error(char *msg);


int main(int argc, char *argv[]) {
    if (argc <= 2) {
        error("Too less argument provided. Usage: <input_file> <output_file>");
    }

    // MPI initialization
    int p, P, I;
    char *input_file = argv[1];
    char *output_file = argv[2];
    int width, height;
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
    MPI_File_read_at_all(in_fh, offset, &width, 1, MPI_INT32_T, MPI_STATUS_IGNORE);
    MPI_File_read_at_all(in_fh, offset += 4, &height, 1, MPI_INT32_T, MPI_STATUS_IGNORE);

    //** Read the image data
    //**** Represent the 2D array in 1D array
    t_grayscale *local_image = (t_grayscale *) malloc((2 * GAUSSIAN_SIZE + I(height, P, p)) * width * sizeof(t_grayscale));
    I = I(height, P, p);
    offset += I_INVERSE(height, P, p, 0) * width * sizeof(t_grayscale);
    MPI_File_read_at(in_fh, offset, CHUNK_BEGIN(local_image, 5, width, sizeof(t_grayscale)), CHUNK_SIZE(P, p, width, height), MPI_T_GRAYSCALE, MPI_STATUS_IGNORE);
    MPI_File_close(&in_fh);



    /*
    
    Canny Edge Detection Algorithm (Main Steps)
    
    */



    // Write output file
    MPI_File out_fh;
    MPI_File_open(MPI_COMM_WORLD, output_file, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &out_fh);
    if(p == 0){
        MPI_File_write(out_fh, &width, 1, MPI_INT32_T, MPI_STATUS_IGNORE);
        MPI_File_write(out_fh, &height, 1, MPI_INT32_T, MPI_STATUS_IGNORE);
    }
    MPI_File_seek(out_fh, offset, MPI_SEEK_SET);
    MPI_File_write_at(out_fh, offset, CHUNK_BEGIN(local_image, 5, width, sizeof(t_grayscale)), CHUNK_SIZE(P, p, width, height), MPI_T_GRAYSCALE, MPI_STATUS_IGNORE);
    MPI_File_close(&out_fh);

    free(local_image);
    MPI_Finalize();
    return 0;
}

void error(char *msg) {
    printf("%s\n", msg);
    exit(1);
}