/**************************************************************************************
* 
* CdL Magistrale in Ingegneria Informatica
* Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2020/21
* 
* Progetto dell'algoritmo Fish School Search 221 231 a
* in linguaggio assembly x86-64 + SSE
* 
* Fabrizio Angiulli, aprile 2019
* 
**************************************************************************************/

/*
* 
* Software necessario per l'esecuzione:
* 
*    NASM (www.nasm.us)
*    GCC (gcc.gnu.org)
* 
* entrambi sono disponibili come pacchetti software 
* installabili mediante il packaging tool del sistema 
* operativo; per esempio, su Ubuntu, mediante i comandi:
* 
*    sudo apt-get install nasm
*    sudo apt-get install gcc
* 
* potrebbe essere necessario installare le seguenti librerie:
* 
*    sudo apt-get install lib64gcc-4.8-dev (o altra versione)
*    sudo apt-get install libc6-dev-i386
* 
* Per generare il file eseguibile:
* 
* nasm -f elf64 fss64.nasm && gcc -m64 -msse -O0 -no-pie sseutils64.o fss64.o fss64c.c -o fss64c -lm && ./fss64c $pars
* 
* oppure
* 
* ./runfss64
* 
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

#define	type		double
#define	MATRIX		type*
#define	VECTOR		type*

typedef struct {
	MATRIX x; //posizione dei pesci
	VECTOR xh; //punto associato al minimo di f, soluzione del problema
	VECTOR c; //coefficienti della funzione
	VECTOR r; //numeri casuali
	int np; //numero di pesci, quadrato del parametro np
	int d; //numero di dimensioni del data set
	int iter; //numero di iterazioni
	type stepind; //parametro stepind
	type stepvol; //parametro stepvol
	type wscale; //parametro wscale
	int display;
	int silent;
} params;

/*
* 
*	Le funzioni sono state scritte assumento che le matrici siano memorizzate 
* 	mediante un array (double*), in modo da occupare un unico blocco
* 	di memoria, ma a scelta del candidato possono essere 
* 	memorizzate mediante array di array (double**).
* 
* 	In entrambi i casi il candidato dovr� inoltre scegliere se memorizzare le
* 	matrici per righe (row-major order) o per colonne (column major-order).
*
* 	L'assunzione corrente � che le matrici siano in row-major order.
* 
*/

void* get_block(int size, int elements) { 
	return _mm_malloc(elements*size,32); 
}

void free_block(void* p) { 
	_mm_free(p);
}

MATRIX alloc_matrix(int rows, int cols) {
	return (MATRIX) get_block(sizeof(type),rows*cols);
}

void dealloc_matrix(MATRIX mat) {
	free_block(mat);
}

/*
* 
* 	load_data
* 	=========
* 
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero
* 	successivi 4 byte: numero di colonne (M) --> numero intero
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri doubleing-point a precisione singola
* 
*****************************************************************************
*	Se lo si ritiene opportuno, � possibile cambiare la codifica in memoria
* 	della matrice. 
*****************************************************************************
* 
*/
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

/*
* 	save_data
* 	=========
* 
*	Salva su file un array lineare in row-major order
*	come matrice di N righe e M colonne
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 64 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 64 bit
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri interi o doubleing-point a precisione singola
*/
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

// PROCEDURE ASSEMBLY

extern void prova(params* input);

void fss(params* input){
	// -------------------------------------------------
	// Codificare qui l'algoritmo Fish Search School
	// -------------------------------------------------
}

int main(int argc, char** argv) {

	char fname[256];
	char* coefffilename = NULL;
	char* randfilename = NULL;
	char* xfilename = NULL;
	int i, j, k;
	clock_t t;
	double time;
	
	//
	// Imposta i valori di default dei parametri
	//

	params* input = malloc(sizeof(params));

	input->x = NULL;
	input->xh = NULL;
	input->c = NULL;
	input->r = NULL;
	input->np = 25;
	input->d = 2;
	input->iter = 350;
	input->stepind = 1;
	input->stepvol = 0.1;
	input->wscale = 10;
	
	input->silent = 0;
	input->display = 0;

	//
	// Visualizza la sintassi del passaggio dei parametri da riga comandi
	//

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

	//
	// Legge i valori dei parametri da riga comandi
	//

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

	//
	// Legge i dati e verifica la correttezza dei parametri
	//

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

	//
	// Visualizza il valore dei parametri
	//

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

	// COMMENTARE QUESTA RIGA!
	prova(input);
	//

	//
	// Fish School Search
	//

	t = clock();
	fss(input);
	t = clock() - t;
	time = ((double)t)/CLOCKS_PER_SEC;

	if(!input->silent)
		printf("FSS time = %.3f secs\n", time);
	else
		printf("%.3f\n", time);

	//
	// Salva il risultato di xh
	//
	sprintf(fname, "xh64_%d_%d_%d.ds2", input->d, input->np, input->iter);
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
