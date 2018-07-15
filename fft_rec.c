#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fft_rec.h"

#define PI 3.14159265358979323846


void fft_h(double* restrict, double* restrict, const double*, const double*, int, int);


/* size is assumed to be the power of 2 */
void fft(double* restrict output_re, double* restrict output_img,
	       const double* input_re,	const double* input_img, int size)
{
	fft_h(output_re, output_img, input_re, input_img, size, 1);
}

/* size is assumed to be the power of 2 */
void inv_fft(double* restrict output_re, double* restrict output_img,
	       const double* input_re,	const double* input_img, int size)
{
	fft_h(output_re, output_img, input_re, input_img, size, - 1);


	for(int idx = 0; idx < size; ++idx)
	{
		output_re[idx] /= size;
		output_img[idx] /= size;
	}
}
	

/* size is assumed to be the power of 2 */
void fft_h(double* restrict output_re, double* restrict output_img,
	       const double* input_re,	const double* input_img, int size, int sign)
{
	/* Deal with the boundary case first */
	if(size == 1)
	{
		output_re[0] = input_re[0];
		output_img[0] = input_img[0];
		return;
	}

	int h_size = (size / 2);
	
	double* buf_mem = (double*) malloc( 4 * size * sizeof(double) );
	
	double* input_ev_re = buf_mem;
	double* input_ev_img = buf_mem + h_size;

	double* input_odd_re = buf_mem + 2 * h_size;
	double* input_odd_img = buf_mem + 3 * h_size;

	double* output_ev_re = buf_mem + 4 * h_size;
	double* output_ev_img = buf_mem + 5 * h_size; 

	double* output_odd_re = buf_mem + 6 * h_size;
	double* output_odd_img = buf_mem + 7 * h_size;

	
	for(int idx = 0; idx < h_size; ++idx)
	{
		int d_idx = (idx << 1);
		
		input_ev_re[idx] = input_re[d_idx];
		input_ev_img[idx] = input_img[d_idx];

		input_odd_re[idx] = input_re[d_idx + 1];
		input_odd_img[idx] = input_img[d_idx + 1];
	}
	
	fft_h(output_ev_re, output_ev_img, input_ev_re, 
			input_ev_img, h_size, sign);

	fft_h(output_odd_re, output_odd_img, input_odd_re,
			input_odd_img, h_size, sign);

	double angle = (2 * PI) / size;
	
	double prim_root_re = cos(angle);
	double prim_root_img = sign * sin(angle);

	double factor_re = 1.0;
	double factor_img = 0.0;

	for(int idx = 0; idx < h_size; ++idx)
	{


		double s_re = factor_re * output_odd_re[idx] - factor_img * output_odd_img[idx];
	       	double s_img = factor_re * output_odd_img[idx] + factor_img * output_odd_re[idx];
	

		output_re[idx] = output_ev_re[idx] + s_re; 
		output_img[idx] = output_ev_img[idx] + s_img;
		

		output_re[h_size + idx] = output_ev_re[idx] - s_re;
		output_img[h_size + idx] = output_ev_img[idx] - s_img;

		
		double t0 = factor_re * prim_root_re - factor_img * prim_root_img;
		double t1 = factor_img * prim_root_re + factor_re * prim_root_img;
		
		factor_re = t0;
		factor_img = t1;
	}



	free(buf_mem);
}
