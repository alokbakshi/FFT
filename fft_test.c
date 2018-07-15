#include <stdio.h>
#include <stdlib.h>
#include <time.h>


#include "fft_rec.h"

int main()
{
	
	const int input_size = 32768;

	srand(time(NULL));

	double* mem_buf = (double*) malloc( 6 * input_size * sizeof(double));

	double* input_re = mem_buf;
	double* input_img = mem_buf + input_size;

	double* output_re = mem_buf + 2 * input_size;
	double* output_img = mem_buf + 3 * input_size;

	double* rev_input_re = mem_buf + 4 * input_size;
	double* rev_input_img = mem_buf + 5 * input_size;

	
	for(int idx = 0; idx < input_size; ++idx)
	{
		input_re[idx] = (1.0 * rand()) / RAND_MAX;
		input_img[idx] = (1.0 * rand()) / RAND_MAX;
	}

	fft(output_re, output_img, input_re, input_img, input_size);
	inv_fft(rev_input_re, rev_input_img, output_re, output_img, input_size);


	for(int idx = 0; idx < input_size; ++idx)
	{

		double re_diff = input_re[idx] - rev_input_re[idx];
		double img_diff = input_img[idx] - rev_input_img[idx];

		printf("Re diff: %f\t Img diff: %f\n", re_diff, img_diff);
	}
	

	free(mem_buf);

	return 0;
}
		
