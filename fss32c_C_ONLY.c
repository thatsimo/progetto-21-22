#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <xmmintrin.h>

#define	type		float
#define	MATRIX		type*
#define	VECTOR		type*

typedef struct {
	MATRIX x; 							//posizione dei pesci
	VECTOR xh; 							//punto associato al minimo di f, soluzione del problema
	VECTOR c; 							//coefficienti della funzione
	VECTOR r; 							//numeri casuali
	int np; 							//numero di pesci, quadrato del parametro np
	int d; 								//numero di dimensioni del data set
	int iter; 							//numero di iterazioni
	type stepind; 						//parametro stepind
	type stepvol; 						//parametro stepvol
	type wscale; 						//parametro wscale
	int display;
	int silent;
} params;

typedef struct {
	type stepind_curr;
	type stepvol_curr;
	type w_sum;
	type w_old;
	type f_sum;
	VECTOR W;
	VECTOR delta_f;
	MATRIX delta_x;
	VECTOR x_new;
	VECTOR V;
	VECTOR f_curr;
	int r_i;
} support;

void* get_block(int size, int elements) { 
	return _mm_malloc(elements*size,16); 
}

void free_block(void* p) { 
	_mm_free(p);
}

VECTOR alloc_vector(int rows) {
	return (VECTOR) get_block(sizeof(type),rows);
}

MATRIX alloc_matrix(int rows, int cols) {
	return (MATRIX) get_block(sizeof(type),rows*cols);
}

void dealloc_matrix(MATRIX mat) {
	free_block(mat);
}

MATRIX load_data(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols, status, i;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}

	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);
	
	MATRIX data = alloc_matrix(rows,cols);
	status = fread(data, sizeof(type), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*k = cols;
	
	return data;
}

void save_data(char* filename, void* X, int n, int k) {
	FILE* fp;
	int i;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&k, 4, 1, fp);
		fwrite(&n, 4, 1, fp);
		for (i = 0; i < n; i++) {
			fwrite(X, sizeof(type), k, fp);
			//printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
			X += sizeof(type)*k;
		}
	}
	else{
		int x = 0;
		fwrite(&x, 4, 1, fp);
		fwrite(&x, 4, 1, fp);
	}
	fclose(fp);
}

support* initialize_support_data_stuct(params* input) {
	support* sup = malloc(sizeof(support));

	sup->stepind_curr = input->stepind;
	sup->stepvol_curr = input->stepvol;
	sup->w_sum = 0;
	sup->w_old = 0;
	sup->f_sum = 0;
	sup->W = alloc_vector(input->np);
	sup->delta_f = alloc_vector(input->np);
	sup->delta_x = alloc_matrix(input->np,input->d);
	sup->x_new = alloc_vector(input->d);
	sup->V = alloc_vector(input->d);
	sup->f_curr = alloc_vector(input->np);
	sup->r_i=0;

	return sup;
}

void find_and_assign_minimum(params* input, support* sup) {
	int k = 0;
	type min_f = sup->f_curr[0];

	for (int i = 0; i<input->np; i++) {
		if (sup->f_curr[i] < min_f) {
			min_f = sup->f_curr[i];
			k = i;
		}
	}

	printf("Valore ottimo : %f\n", min_f);

	input->xh = alloc_vector(input->d);

	for (int j  = 0; j<input->d; j++)
		input->xh[j] = input-> x[k*input->d+j];
}

void update_parameters(params* input, support* sup) {
	sup->stepind_curr = sup->stepind_curr - input->stepind / input->iter;
	sup->stepvol_curr = sup->stepvol_curr - input->stepvol / input->iter; 
}

void compute_weighted_avg(int np, int d, VECTOR acc, MATRIX num_1, VECTOR num_2, type den) {
    for (int i = 0; i < d; ++i)			
        acc[i] = 0;

    if (den == 0) return;

    for (int i = 0; i < np; i++)
        for (int j = 0; j < d; j++)
            acc[j] += num_1[i*d+j] * (num_2[i]/den);
}

void euclidian_distance(MATRIX x, int offset, VECTOR v,int d,type* dist) {
	type ret=0;

	for(int j=0;j<d;++j)
		ret+=(x[offset+j]-v[j])*(x[offset+j]-v[j]);

	*dist=sqrt(ret);
}

void volitive_movement(params* input, support* sup) {
	if( sup->w_sum == 0) {
		sup->w_sum = sup->w_old;
		return;
	}

    for (int i=0;i<input->d;++i) 
		sup->V[i]=0;

	if(sup->w_sum!=0)
		compute_weighted_avg(input->np,input->d,sup->V,input->x,sup->W,sup->w_sum);

	type sgn = (sup->w_old - sup->w_sum < 0 ) ? -1 : 1;

	for (int i = 0; i<input->np; i++) {
		type dist_x_i_B=0;

		euclidian_distance(input->x, i*input->d, sup->V, input->d,&dist_x_i_B);

		if (dist_x_i_B == 0) continue;

		for (int j = 0; j<input->d; j++)
			input->x[input->d*i+j] = input->x[input->d*i+j] + sgn * sup->stepvol_curr * input->r[sup->r_i] * ((input->x[input->d*i+j] - sup->V[j]) / dist_x_i_B);
		
		sup->r_i++;
	}
}

void instincitve_movement(params* input, support* sup) {
	if(sup->f_sum==0) return;

	int d=input->d,np=input->np;
	
	for(int i=0;i<input->d;++i)
		sup->V[i]=0;

	compute_weighted_avg(input->np,input->d,sup->V,sup->delta_x,sup->delta_f,sup->f_sum);

	for (int i = 0; i < np; i++)
		for(int j=0;j<d;++j)
			input->x[i*input->d+j]+=sup->V[j];
}

void alimentation_operator(params* input, support* sup) {
	type min_delta_f = sup->delta_f[0];
	
	for(int i=0;i<input->np;++i)
		if(sup->delta_f[i]<min_delta_f) min_delta_f=sup->delta_f[i];
	
    if (min_delta_f == 0) return;

    sup->w_old = sup->w_sum;
    sup->w_sum = 0;

    for (int i = 0; i< input->np-3; i+=4) {
        sup->W[i] += sup->delta_f[i]/min_delta_f;
        sup->w_sum += sup->W[i];

		sup->W[i+1] += sup->delta_f[i+1]/min_delta_f;
        sup->w_sum += sup->W[i+1];

		sup->W[i+2] += sup->delta_f[i+2]/min_delta_f;
        sup->w_sum += sup->W[i+2];

		sup->W[i+3] += sup->delta_f[i+3]/min_delta_f;
        sup->w_sum += sup->W[i+3];
    }

	for(int i=(input->np/4)*4;i<input->np;++i) {
		sup->W[i] += sup->delta_f[i]/min_delta_f;
        sup->w_sum += sup->W[i];
	}
}

type evaluate_f(MATRIX x, VECTOR c, int i, int d) {
	type quad=0;
	type scalar=0;
	
	for(int j=0;j<d;++j) {
		quad+=x[i*d+j]*x[i*d+j];
		scalar+=x[i*d+j]*c[j];
	}

	return expf(quad) + quad- scalar;
}


void individual_movement(int i, params* input, support* sup) {
	for (int j = 0; j < input->d; ++j) {
        sup->x_new[j] = input->x[i*input->d+j] + ((input->r[sup->r_i++])*2-1) * sup->stepind_curr;
    }	
				
	type f_new_i = evaluate_f(sup->x_new, input->c, 0, input->d);
	type delta_f_i = f_new_i - sup->f_curr[i];

	if (delta_f_i < 0) {
		sup->delta_f[i] = delta_f_i;
		sup->f_curr[i] = f_new_i;
		sup->f_sum = sup->f_sum + delta_f_i;

		for(int j = 0; j < input->d; ++j) {
			sup->delta_x[i*input->d+j] = sup->x_new[j] - input->x[i*input->d+j];				// aggiornamento delta_x
			input->x[i*input->d+j] = sup->x_new[j];										// aggiornamento posizione pesce	
		}
	} else {
		for (int j = 0; j < input->d; ++j) 
			sup->delta_x[i*input->d+j] = 0;

		sup->delta_f[i]=0;
	}
}

void fss(params* input) {
	support* sup = initialize_support_data_stuct(input);

	for (int i = 0; i < input->np-3; i+=4) {			
		sup->W[i] = input->wscale/2;
		sup->w_sum = sup->w_sum + sup->W[i];
		sup->f_curr[i] = evaluate_f(input->x, input->c, i, input->d);

		sup->W[i+1] = input->wscale/2;
		sup->w_sum = sup->w_sum + sup->W[i+1];
		sup->f_curr[i+1] = evaluate_f(input->x, input->c, i+1, input->d);

		sup->W[i+2] = input->wscale/2;
		sup->w_sum = sup->w_sum + sup->W[i+2];
		sup->f_curr[i+2] = evaluate_f(input->x, input->c, i+2, input->d);

		sup->W[i+3] = input->wscale/2;
		sup->w_sum = sup->w_sum + sup->W[i+3];
		sup->f_curr[i+3] = evaluate_f(input->x, input->c, i+3, input->d);
	}

	for (int i =(input->np/4)*4; i < input->np; i++) {
		sup->W[i] = input->wscale/2;
		sup->w_sum = sup->w_sum + sup->W[i];
		sup->f_curr[i] = evaluate_f(input->x, input->c, i, input->d);
	}

	for (int i = 0; i < input->iter; ++i) {
		sup->f_sum = 0;

		for (int j=0;j < input->np;++j)	
			individual_movement(j, input, sup);
		
		alimentation_operator(input, sup);
		instincitve_movement(input, sup);
		volitive_movement(input, sup);
		update_parameters(input, sup);
	}

	find_and_assign_minimum(input, sup);
}

int main(int argc, char** argv) {
	char fname[256];
	char* coefffilename = NULL;
	char* randfilename = NULL;
	char* xfilename = NULL;
	int i, j, k;
	double time;
	
	clock_t t;
	params* input = malloc(sizeof(params));

	input->x = NULL;
	input->xh = NULL;
	input->c = NULL;
	input->r = NULL;
	input->np = 25;
	input->d = 7;
	input->iter = 350;
	input->stepind = 1;
	input->stepvol = 0.1;
	input->wscale = 10;
	
	input->silent = 0;
	input->display = 1;

	if(argc <= 1){
		printf("%s -c <c> -r <r> -x <x> -np <np> -si <stepind> -sv <stepvol> -w <wscale> -it <itmax> [-s] [-d]\n", argv[0]);
		printf("\nParameters:\n");
		printf("\tc: il nome del file ds2 contenente i coefficienti\n");
		printf("\tr: il nome del file ds2 contenente i numeri casuali\n");
		printf("\tx: il nome del file ds2 contenente le posizioni iniziali dei pesci\n");
		printf("\tnp: il numero di pesci, default 25\n");
		printf("\tstepind: valore iniziale del parametro per il movimento individuale, default 1\n");
		printf("\tstepvol: valore iniziale del parametro per il movimento volitivo, default 0.1\n");
		printf("\twscale: valore iniziale del peso, default 10\n");
		printf("\titmax: numero di iterazioni, default 350\n");
		printf("\nOptions:\n");
		printf("\t-s: modo silenzioso, nessuna stampa, default 0 - false\n");
		printf("\t-d: stampa a video i risultati, default 0 - false\n");
		exit(0);
	}

	int par = 1;

	while (par < argc) {
		if (strcmp(argv[par],"-s") == 0) {
			input->silent = 1;
			par++;
		} else if (strcmp(argv[par],"-d") == 0) {
			input->display = 1;
			par++;
		} else if (strcmp(argv[par],"-c") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing coefficient file name!\n");
				exit(1);
			}
			coefffilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-r") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing random numbers file name!\n");
				exit(1);
			}
			randfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-x") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing initial fish position file name!\n");
				exit(1);
			}
			xfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-np") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing np value!\n");
				exit(1);
			}
			input->np = atoi(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-si") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing stepind value!\n");
				exit(1);
			}
			input->stepind = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-sv") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing stepvol value!\n");
				exit(1);
			}
			input->stepvol = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-w") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing wscale value!\n");
				exit(1);
			}
			input->wscale = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-it") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing iter value!\n");
				exit(1);
			}
			input->iter = atoi(argv[par]);
			par++;
		} else{
			printf("WARNING: unrecognized parameter '%s'!\n",argv[par]);
			par++;
		}
	}

	if(coefffilename == NULL || strlen(coefffilename) == 0){
		printf("Missing coefficient file name!\n");
		exit(1);
	}

	if(randfilename == NULL || strlen(randfilename) == 0){
		printf("Missing random numbers file name!\n");
		exit(1);
	}

	if(xfilename == NULL || strlen(xfilename) == 0){
		printf("Missing initial fish position file name!\n");
		exit(1);
	}

	int x,y;

	input->c = load_data(coefffilename, &input->d, &y);
	input->r = load_data(randfilename, &x, &y);
	input->x = load_data(xfilename, &x, &y);

	if(input->np < 0){
		printf("Invalid value of np parameter!\n");
		exit(1);
	}

	if(input->stepind < 0){
		printf("Invalid value of si parameter!\n");
		exit(1);
	}

	if(input->stepvol < 0){
		printf("Invalid value of sv parameter!\n");
		exit(1);
	}

	if(input->wscale < 0){
		printf("Invalid value of w parameter!\n");
		exit(1);
	}

	if(input->iter < 0){
		printf("Invalid value of it parameter!\n");
		exit(1);
	}

	if(!input->silent){
		printf("Coefficient file name: '%s'\n", coefffilename);
		printf("Random numbers file name: '%s'\n", randfilename);
		printf("Initial fish position file name: '%s'\n", xfilename);
		printf("Dimensions: %d\n", input->d);
		printf("Number of fishes [np]: %d\n", input->np);
		printf("Individual step [si]: %f\n", input->stepind);
		printf("Volitive step [sv]: %f\n", input->stepvol);
		printf("Weight scale [w]: %f\n", input->wscale);
		printf("Number of iterations [it]: %d\n", input->iter);
	}

	t = clock();
	fss(input);
	t = clock() - t;

	time = ((double)t)/CLOCKS_PER_SEC;

	if(!input->silent)
		printf("FSS time = %.3f secs\n", time);
	else
		printf("%.3f\n", time);

	sprintf(fname, "xh32_%d_%d_%d.ds2", input->d, input->np, input->iter);
	save_data(fname, input->xh, 1, input->d);

	if(input->display){
		if(input->xh == NULL)
			printf("xh: NULL\n");
		else{
			printf("xh: [");
			for(i=0; i<input->d-1; i++)
				printf("%f,", input->xh[i]);
			printf("%f]\n", input->xh[i]);
		}
	}

	if(!input->silent)
		printf("\nDone.\n");

	return 0;
}