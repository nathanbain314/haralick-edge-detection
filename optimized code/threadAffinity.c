// Copyright (C) 2011-2012, Haldo Spontón <haldos@fing.edu.uy>
// Copyright (C) 2011-2012, Juan Cardelino <juanc@fing.edu.uy>

/*	
	This program is free software: you can use, modify and/or
	redistribute it under the terms of the GNU General Public
	License as published by the Free Software Foundation, either
	version 3 of the License, or (at your option) any later
	version. You should have received a copy of this license along
	this program. If not, see <http://www.gnu.org/licenses/>.
	
	This program is free software: you can use, modify and/or
	redistribute it under the terms of the simplified BSD
	License. You should have received a copy of this license along
	this program. If not, see
	<http://www.opensource.org/licenses/bsd-license.html>.
	
	This program is provided for research and education only: you can
	use and/or modify it for these purposes, but you are not allowed
	to redistribute this work or derivative works in source or
	executable form. A license must be obtained from the patent right
	holders for any other use.
*/

/** @file test_haralick.c
* \author Haldo Spontón <haldos@fing.edu.uy> & Juan Cardelino <juanc@fing.edu.uy>
* \date May, 2012
* \see ``Review of edge detectors´´ IPOL publication.
* \brief Implements the Haralick edge detection algorithm. Haralick's algorithm basically approaches the neighborhood of a pixel using a bicubic polynomial function. Then evaluate certain conditions on the parameters found for the model, equivalent to find zero crossings in the second derivative of the image in that neighborhood.
*/

//  Software Guide : BeginLatex
//  Haralick edge detectors, main C file.\\
//  
//  Parameters:
//	\begin{itemize}
//		\item \texttt{input\_image} - Input image.
//		\item \texttt{rhozero} - Threshold for the Haralick condition $|\frac{C_2}{2C_3}|\leq\rho_0$.
//		\item \texttt{padding\_method} - Padding method flag (in convolution): 0 means zero-padding, 1 means image boundary reflection.
//		\item \texttt{output\_image} - Output image (edges).
//	\end{itemize}
//
//	Includes:
//  Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
#define _GNU_SOURCE
#include <sched.h>
	#include "iio.c"
	#include "2dconvolution.c"
	#include <time.h>
  #include <pthread.h>
 
// Software Guide : EndCodeSnippet

/*!
 \fn int main_haralick(int argc, char *argv[])
 \brief Main function of the Haralick edge detection algorithms.
 @param input_image Input image filename.
 @param rho_zero Threshold for the Haralick edge condition.
 @param padding_method Padding method flag (in convolution): 0 means zero-padding, 1 means image boundary reflection.
 @param output Output image filename.
 \return None.
 \ingroup haralick
	\note The real name of this function is main. Function name was temporarily changed to the proper functioning of doxygen.
*/

#ifdef _POSIX_MONOTONIC_CLOCK
#ifdef CLOCK_MONOTONIC_RAW
static clockid_t CLOCKTYPE = CLOCK_MONOTONIC_RAW;
static const char* SCLOCKTYPE = "CLOCK_MONOTONIC_RAW";
#else
static clockid_t CLOCKTYPE = CLOCK_MONOTONIC 
static const char* SCLOCKTYPE = "CLOCK_MONOTONIC";
#endif /* CLOCK_MONOTONIC_RAW */
#else
#warning "CLOCK_MONOTONIC is unavailable, using CLOCK_REALTIME"
static clockid CLOCKTYPE = CLOCK_REALTIME;
static const char* SCLOCKTYPE = "CLOCK_REALTIME";
#endif /* _POSIX_MONOTONIC_CLOCK */

#define NANOSEC 1000000000LL

int num_threads;

struct thread_data
{
  float rhozero;
  float *edges;
  double **masks;
  double *aux;
  int w;
  int h;
  int thread_id;
	int *offsets;
};

void *haralick_thread(void *threadarg)
{
  struct thread_data *my_data;
  my_data = (struct thread_data *) threadarg;
  float rhozero, *edges;
  double **masks, *aux;
  int w, h, thread_id, *offsets;
  masks = (double**)my_data->masks;
  aux = my_data->aux;
  rhozero = my_data->rhozero;
  edges = my_data->edges;
  w = my_data->w;
  h = my_data->h;
  thread_id = my_data->thread_id;
	offsets = my_data->offsets;

  int i_zp, u, v, num_edges, size, i, start, end;
  size = w * h /(num_threads);
  start = thread_id * size;
  end = start + size;
  num_edges = 0;

    num_edges = 0;
    double k[10];
    double acum;
    double C2, C3, denom, sintheta, costheta, k1sqr, k2sqr;
		int w_1 = w - 1;
		i = start;
    i_zp=i+4*w+8*(i/w)+36;
		for( i ;i<end;i++, i_zp++){
        double *neighborhood = get_neighborhood(aux, i_zp, 5, offsets);
        // k1 to k10 (note: k1 (u=0) is not necessary)
        for(u=1;u<10;u++) k[u] = 0;
          acum = 0;
          for(v=0;v<25;v++){
            acum = neighborhood[v];
            k[1] += acum*masks[1][v];
            k[2] += acum*masks[2][v];
            k[3] += acum*masks[3][v];
            k[4] += acum*masks[4][v];
            k[5] += acum*masks[5][v];
            k[6] += acum*masks[6][v];
            k[7] += acum*masks[7][v];
            k[8] += acum*masks[8][v];
            k[9] += acum*masks[9][v];
          }
        //  k[u] = acum;

        // compute C2 and C3
        sintheta = k[1];
        costheta = k[2];
        k1sqr = sintheta * sintheta;
        k2sqr = costheta * costheta;
        denom = sqrt( k1sqr + k2sqr);
        C2 = k[3]*k1sqr + k[4]*sintheta*costheta + k[5]*k2sqr;
        C3 = k[6]*k1sqr*sintheta + k[7]*k1sqr*costheta +
           k[8]*sintheta*k2sqr + k[9]*costheta*k2sqr;
        //if ((fabs(C2 / (3*C3))<=rhozero)&&(C3<=0)) {
        if ((fabs(denom * C2 / (3*C3))<=rhozero)) {
          edges[i] = 255;
          num_edges += 1;
        }
        // free neighborhood
        free_neighborhood(neighborhood);
				if(i%w == w_1) i_zp += 8;
      }


	pthread_exit((void *)(int)num_edges);
}

//  Software Guide : BeginLatex
//	\vspace{0.5cm}
//	\Large{Main function} \\
//  Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
int main(int argc, char *argv[]) {
// Software Guide : EndCodeSnippet
	if (argc != 6) {
		printf("Usage: %s input_image rhozero padding_method output num_threads\n", argv[0]);
	} else {

		// Execution time:
		struct timespec begin, end;
	
		// Parameters
		float rhozero = atof(argv[2]);
		int i,j;
		int padding_method = atoi(argv[3]);
		num_threads = atoi(argv[5]) - 1;
	
		// Load input image (using iio)
//  Software Guide : BeginLatex
//	Load input image (using \textit{iio}): \\
//  Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
		int w, h, pixeldim;
		float *im_orig = iio_read_image_float_vec(argv[1], &w, &h, &pixeldim);
// Software Guide : EndCodeSnippet
		fprintf(stderr, "Input image loaded:\t %dx%d image with %d channel(s).\n", w, h, pixeldim);

		// Grayscale conversion (if necessary)
//  Software Guide : BeginLatex
//	Grayscale conversion (if necessary): explained in \ref{app:marr-hildreth}. \\ \\
//  Software Guide : EndLatex
		double *im = malloc(w*h*sizeof(double));
		if (im == NULL){
			fprintf(stderr, "Out of memory...\n");
			exit(EXIT_FAILURE);
		}
		int z;
		int zmax = w*h;
		if (pixeldim==3){
			for(z=0;z<zmax;z++){
				im[z] =  (double)(6968*im_orig[3*z] + 23434*im_orig[3*z + 1] + 2366*im_orig[3*z + 2])/32768;
			}
			fprintf(stderr, "images converted to grayscale\n");
		} else {
			for(z=0;z<zmax;z++){
				im[z] = (double)im_orig[z];
			}
			fprintf(stderr, "images are already in grayscale\n");
		}

		// Haralick's masks for computing k1 to k10
//		double masks[10][25] = { {-13,   2,   7,  2, -13,   2,  17,  22,  17,   2,  7,  22,  27, 22,  7,  2,  17, 22, 17,  2, -13,   2,  7,  2, -13},
//								 { 31,  -5, -17, -5,  31, -44, -62, -68, -62, -44,  0,   0,   0,  0,  0, 44,  62, 68, 62, 44, -31,   5, 17,  5, -31},
//								 { 31, -44,   0, 44, -31,  -5, -62,   0,  62,   5, -17, -68,  0, 68, 17, -5, -62,  0, 62,  5,  31, -44,  0, 44, -31},
// 								 {  2,   2,   2,  2,   2,  -1,  -1,  -1,  -1,  -1,  -2,  -2, -2, -2, -2, -1,  -1, -1, -1, -1,   2,   2,  2,  2,   2},
//								 {  4,   2,   0, -2,  -4,   2,   1,   0,  -1,  -2,   0,   0,  0,  0,  0, -2,  -1,  0,  1,  2,  -4,  -2,  0,  2,   4},
//								 {  2,  -1,  -2, -1,   2,   2,  -1,  -2,  -1,   2,   2,  -1, -2, -1,  2,  2,  -1, -2, -1,  2,   2,  -1, -2, -1,   2},
//								 { -1,  -1,  -1, -1,  -1,   2,   2,   2,   2,   2,   0,   0,  0,  0,  0, -2,  -2, -2, -2, -2,   1,   1,  1,  1,   1},
//								 { -4,  -2,   0,  2,   4,   2,   1,   0,  -1,  -2,   4,   2,  0, -2, -4,  2,   1,  0, -1, -2,  -4,  -2,  0,  2,   4},
//								 { -4,   2,   4,  2,  -4,  -2,   1,   2,   1,  -2,   0,   0,  0,  0,  0,  2,  -1, -2, -1,  2,   4,  -2, -4, -2,   4},
//								 { -1,   2,   0, -2,   1,  -1,   2,   0,  -2,   1,  -1,   2,  0, -2,  1, -1,   2,  0, -2,  1,  -1,   2,  0, -2,   1} };
//		double weights[10] = {175, 420, 420, 70, 100, 70, 60, 140, 140, 60};
//		for(i=0;i<10;i++){
//			for(j=0;j<25;j++){
//				masks[i][j] /= weights[i];
//			}
//		}
		// [update 23/02/2012]
		// New masks calculated by 2-d fitting using LS, with the function
		// f(x,y) = k1 + k2*x + k3*y + k4*x² + k5*xy + k6*y² + k7*x³ + k8*x²y + k9*xy² + k10*y³.
//  Software Guide : BeginLatex
//	Masks calculated by 2-d fitting (using LS) with the function: 
//	$$
//	f(x,y) = k_1 + k_2x + k_3y + k_4x^2 + k_5xy + k_6*y^2 + k_7x^3 + k_8x^2y + k_9xy^2 + k_{10}y^3 \\
//	$$
//  Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
		double masks[10][25] = { {   425,   275,  225,  275,  425,   
									 275,   125,   75,  125,  275,   
                                     225,    75,   25,   75,  225,   
									 275,   125,   75,  125,  275,   
									 425,   275,  225,  275,  425},
								 { -2260,  -620,    0,  620,  2260, 
								   -1660,  -320,    0,  320,  1660, 
								   -1460,  -220,    0,  220,  1460,
                                   -1660,  -320,    0,  320,  1660, 
                                   -2260,  -620,    0,  620, 2260},
								// matrix continues, too large to display in documentation.
// Software Guide : EndCodeSnippet
								 {  2260,  1660, 1460, 1660,  2260,   620,   320, 220,  320,  620,     0,    0,   0,   0,    0,  -620,  -320, -220,  -320,  -620, -2260, -1660, -1460, -1660, -2260},
								 {  1130,   620,  450,  620,  1130,   830,   320, 150,  320,  830,   730,  220,  50, 220,  730,   830,   320,  150,   320,   830,  1130,   620,   450,   620,  1130},
								 {  -400,  -200,    0,  200,   400,  -200,  -100,   0,  100,  200,     0,    0,   0,   0,    0,   200,   100,    0,  -100,  -200,   400,   200,     0,  -200,  -400},
								 {  1130,   830,  730,  830,  1130,   620,   320, 220,  320,  620,   450,  150,  50, 150,  450,   620,   320,  220,   320,   620,  1130,   830,   730,   830,  1130},
								 { -8260, -2180,    0, 2180,  8260, -6220, -1160,   0, 1160, 6220, -5540, -820,   0, 820, 5540, -6220, -1160,    0,  1160,  6220, -8260, -2180,     0,  2180,  8260},
								 {  5640,  3600, 2920, 3600,  5640,  1800,   780, 440,  780, 1800,     0,    0,   0,   0,    0, -1800,  -780, -440,  -780, -1800, -5640, -3600, -2920, -3600, -5640},
								 { -5640, -1800,    0, 1800,  5640, -3600,  -780,   0,  780, 3600, -2920, -440,   0, 440, 2920, -3600,  -780,    0,   780,  3600, -5640, -1800,     0,  1800,  5640},
								 {  8260,  6220, 5540, 6220,  8260,  2180,  1160, 820, 1160, 2180,     0,    0,   0,   0,    0, -2180, -1160, -820, -1160, -2180, -8260, -6220, -5540, -6220, -8260} };

		// Initialise edge image
//  Software Guide : BeginLatex
//	Initialise edge image using \texttt{calloc}:
//  Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
		float *edges = calloc(w*h,sizeof(float));
// Software Guide : EndCodeSnippet

		// Zero-padding
//  Software Guide : BeginLatex
//	Padding: a larger auxiliar image \texttt{aux} is required to compute the coefficients $k_1$ to $k_{10}$ in every pixel of the original image. \\
//	Two different methods are implemented: zero-padding (\texttt{padding\_method}$=0$) and reflection of original image (\texttt{padding\_method}$=1$).
//  Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
		int wx = (w+8);
		int hx = (h+8);
		double *aux = calloc(wx*hx,sizeof(double));
		int fila,col;
		int imax = wx*hx;
		if (padding_method == 0) {
			for(i=0;i<imax;i++){
				fila = (int)(i/wx);
				col = i-(wx*fila);	
				if ( (fila>=4)&&(col>=4)&&(fila<h+4)&&(col<w+4) ) {
					aux[i] = im[(col-4)+(w*(fila-4))];
				}
			}
		}
// Software Guide : EndCodeSnippet
// Software Guide : BeginCodeSnippet
		if (padding_method == 1) {
			int fila_refl, col_refl;
			for(i=0;i<imax;i++){
				fila = (int)(i/wx);
				col = i-(wx*fila);
				if (fila<4) {
					fila_refl = 7 - fila;
					if (col<4) { //zone1
						col_refl = 7 - col;
					} else if (col<w+4) {	//zone2
						col_refl = col;
					} else { //zone3
						col_refl = 2*w + 7 - col;
					}
				} else if (fila<h+4) {
					fila_refl = fila;
					if (col<4) { //zone4
						col_refl = 7 - col;
					} else if (col<w+4) { //image
						col_refl = col;
					} else { //zone5
						col_refl =  2*w + 7 - col;
					}
				} else {
					fila_refl = 2*h + 7 - fila;
					if (col<4) { //zone6
						col_refl =	7 - col;
					} else if (col<w+4) {	//zone7
						col_refl = col;
					} else { //zone8
						col_refl =  2*w + 7 - col;
					}
				}
				aux[i] = im[(col_refl-4)+(w*(fila_refl-4))];
			} //for
		}
// Software Guide : EndCodeSnippet
		
		// Haralick's algorithm
//  Software Guide : BeginLatex
//	Haralick algorithm: coefficients $k_1$ to $k_{10}$ are computed in every pixel of the original image 
//	(using the function \texttt{get\_neighbors\_offset} to get the index offsets of the neighbor pixels and 
//	the function \texttt{get\_neighborhood} to get the neighborhood of a pixel using those index offsets). 
//	Once the coefficients are calculated, are computed
//	$$
//	C_2 = \frac{k_2^2k_4 + k_2k_3k_5 + k_3^2k_6}{k_2^2 + k_3^2}
//	$$
//	and
//	$$
//	C_3 = \frac{k_2^3k_7 + k_2^2k_3k_8 + k_2k_3^2k_9 + k_3^3k_{10}}{(\sqrt{k_2^2 + k_3^2})^3},
//	$$
//	and then the edge condition is evaluated in every pixel; $|\frac{C_2}{2C_3}|\leq\rho_0$. \\
//  Software Guide : EndLatex
//  Software Guide : BeginCodeSnippet
		clock_gettime(CLOCKTYPE, &begin);

		pthread_t threads[num_threads];
		struct thread_data td[num_threads];
		pthread_attr_t attr;
		cpu_set_t cpus;
		pthread_attr_init(&attr);

		int u,v;

		int *offsets = get_neighbors_offset(wx, 5);

    for(i=0; i<num_threads; i++){
      td[i].rhozero = rhozero;
      td[i].edges = edges;
      td[i].masks = malloc(sizeof(double) * 10);
      for(u=1;u<10;u++)
      {   
      td[i].masks[u] = malloc(sizeof(double) * 25);
      for(v=0;v<25;v++)
      td[i].masks[u][v]=masks[u][v];
      }   
      td[i].aux = aux;
      td[i].w = w;
      td[i].h = h;
      td[i].thread_id = i;
			td[i].offsets = offsets;

			CPU_ZERO(&cpus);
			CPU_SET(i, &cpus);
			pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpus);

      pthread_create(&threads[i], &attr, haralick_thread, (void *) &td[i]);
    }  
		int i_zp, num_edges, start, stop;
		start = i * w * h / (num_threads+1);
		stop = w * h;
		num_edges = 0;
		
		double k[10];
		double acum;
		double C2, C3, denom, sintheta, costheta, k1sqr, k2sqr;
		int w_1 = w - 1;

    i = start;
    i_zp=i+4*w+8*(i/w)+36;

    for( i=start ;i<stop;i++, i_zp++){
				double *neighborhood = get_neighborhood(aux, i_zp, 5, offsets);
				// k1 to k10 (note: k1 (u=0) is not necessary)
				for(u=1;u<10;u++) k[u] = 0;
					acum = 0;
					for(v=0;v<25;v++){
						acum = neighborhood[v];
						k[1] += acum*masks[1][v];
						k[2] += acum*masks[2][v];
						k[3] += acum*masks[3][v];
						k[4] += acum*masks[4][v];
						k[5] += acum*masks[5][v];
						k[6] += acum*masks[6][v];
						k[7] += acum*masks[7][v];
						k[8] += acum*masks[8][v];
						k[9] += acum*masks[9][v];
					}
				//	k[u] = acum;
				
				// compute C2 and C3
				sintheta = k[1];
				costheta = k[2];
				k1sqr = sintheta * sintheta;
				k2sqr = costheta * costheta;
				denom = sqrt( k1sqr + k2sqr);
				C2 = k[3]*k1sqr + k[4]*sintheta*costheta + k[5]*k2sqr;
				C3 = k[6]*k1sqr*sintheta + k[7]*k1sqr*costheta +
					 k[8]*sintheta*k2sqr + k[9]*costheta*k2sqr;
				//if ((fabs(C2 / (3*C3))<=rhozero)&&(C3<=0)) {
				if ((fabs(denom * C2 / (3*C3))<=rhozero)) {
					edges[i] = 255;
					num_edges += 1;
				}
				// free neighborhood
				free_neighborhood(neighborhood);
				if(i%w == w_1) i_zp += 8;
		}

    void *returnValue;
    for(i=0;i<num_threads;i++)
    {
      pthread_join( threads[i], &returnValue);
      num_edges += (int)returnValue;
    }

//  Software Guide : EndCodeSnippet
		
		// Result
		fprintf(stderr, "%d edge points found...\n", num_edges);

		// Save output image
//	Software Guide : BeginLatex
//	Save output image (using \textit{iio}): \\
//	Software Guide : EndLatex
//	Software Guide : BeginCodeSnippet
		iio_save_image_float_vec(argv[4], edges, w, h, 1);
//	Software Guide : EndCodeSnippet
		fprintf(stderr, "Output Image saved in %s:\t %dx%d image with %d channel(s).\n", argv[4], w, h, pixeldim);
	
		// Free memory
		free(im_orig);
		free(im);
		free(aux);
		free(edges);

		fprintf(stderr, "haralick's edge detector computation done.\n");

		// Execution time:
		clock_gettime(CLOCKTYPE, &end);
		unsigned long long exectime = ((end).tv_sec - (begin).tv_sec) * NANOSEC + (end).tv_nsec - (begin).tv_nsec;
		fprintf(stderr, "\texecution time: %llu ns.\n", exectime);		

		return 0;
	
	} // else (argc)
}
