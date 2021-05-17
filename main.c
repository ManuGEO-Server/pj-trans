#include <stdlib.h>
#include <stdio.h>
#include <sndfile.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <unistd.h>

#define SWAP(a,b) ctmp=(a); (a)=(b); (b)=ctmp
#define ARRAY_LEN(x) ((int)(sizeof(x)/sizeof(x[0])))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define LOG2 10
#define BUFFER_LEN 1024
#define freqEchen 4096
#define HAUTEUR 40
#define LARGEUR 80
#define NBPOINT 20
#define bitrate 120

complex TW[BUFFER_LEN];
int freq_fondamentale;
SNDFILE *infile, *outfile;
SF_INFO sfinfo;

int bitrev(int inp, int numbits) {
	int rev=0;
	for(int i=0; i<numbits;i++) {
		rev=(rev << 1) | (inp & 1);
		inp >>=1;
	}

	return rev;
}

// fonction sfx_mix_mono_read_double
sf_count_t sfx_mix_mono_read_double(SNDFILE * file, double * data, sf_count_t datalen) {
	SF_INFO info;
	static double multi_data [2048];
	double mix=0.0;
	int k, ch, frames_read;
	sf_count_t dataout = 0;
	sf_command (file, SFC_GET_CURRENT_SF_INFO, &info, sizeof (info));
	if(info.channels == 1)
		return sf_read_double (file, data, datalen);
  
	while(dataout < datalen) {
		int this_read;
		this_read = MIN (ARRAY_LEN (multi_data) / info.channels, datalen);
		frames_read = sf_readf_double (file, multi_data, this_read);
		if(frames_read == 0)
			break;
      
		for(k = 0 ; k < frames_read ; k++) {
			mix = 0;
			for(ch = 0 ; ch < info.channels ; ch++)
				mix += multi_data [k * info.channels + ch];
          
			data [dataout + k] = mix / info.channels;
		}
      
		dataout += this_read;
	}

	return dataout;
}

// fft recurrente
void fftrec(complex *data, complex *result, unsigned int size, int log2n) {
	complex ypair [size], yimpair [size], Fimpair [size], Fpair [size];
	int n, N2=size/2;

	if(size > 1) {
		for(n=0;n<N2;n++) {
			ypair [n] = data [n]+ data [n+N2];
			yimpair [n] = (data [n] - data [n+N2]) * cexp(-2*I*M_PI*n/size);
		}
		
		fftrec(ypair, Fpair, N2, log2n);
		fftrec(yimpair, Fimpair, N2, log2n);
		
		for(n=0;n<N2;n++) {
			result[2*n] = Fpair[n];
			result[2*n+1] = Fimpair[n];
		}
	}
	
	else
		result[0] = data[0];
}

//fft interative
void twiddle(complex *TW, unsigned int size) {
	complex phi = cexp(-2*I*M_PI/size);
	TW[0] = 1;
	for(int i=1;i<size;i++)
		TW[i] = TW[i-1]*phi;
}

void fftiterTW(complex *data, unsigned int size, int log2n) {
	int i, j, N2, Bpair, Bimpair, Bp=1, N=size, b, n;
	complex impair, pair, ctmp;
	for(int k=0;k<log2n;k++) {
		N2=N/ 2;
		Bpair =0;
		for(b=0;b<Bp;b++) {
			Bimpair=Bpair+N2 ;
			for (n=0;n<N2;n++) {
				impair= data [Bpair+n] + data [Bimpair+n];
				pair = (data [Bpair+n] - data [Bimpair+n])*TW[n*size/N];
				data [Bpair+n] = pair;
				data [Bimpair+n] = impair;
			}
			
			Bpair = Bpair+N;
		}
		
		Bp = Bp*2;
		N = N/2;
	}
	
	for(i=0;i<size;i++) {
		j = bitrev(i,log2n);
		if(j>i)
			SWAP(data [j], data [i]);  
	}
	
	for(i=size-1;i>0;i--)
		data [i] = data [i-1];
  
	data [0] = ctmp;
}

//dft
void dft(complex *data, complex *result, unsigned int buffer) {
	int k,n;
	for(k=0;k<BUFFER_LEN;k++) {
		for(n=0;n<BUFFER_LEN;n++)
			result[k] = result[k] + data[n] * cexp(-2*I*M_PI*k*n/buffer);
	}
}

//conv double complex
complex convDoubleComplex(double val) {
	return val+I*0.0;
}

//stop
void stop(int t) {
	clock_t time = clock();
	while((clock()-time) <= (t*CLOCKS_PER_SEC/1000));
}

//maximum
double maximum(double tab[], int longueur) {
	double max = 0.0;
	for(int i=0; i<longueur; i++) {
		if(tab[i]>max) {
			max = tab[i];
			freq_fondamentale = i;
		}
	}
	
	return max;
}

//module complexe
double complexModule(complex c) {
	return sqrt(creal(c)*creal(c)+cimag(c)*cimag(c));
}

//afficher spectre
void afficherSpectre(double spectre[], double max) {
	system("clear");
	int c = max, l = HAUTEUR, i;
	while(l>=0) {
		while (c<l) {
			printf("\n");
			l--;
		}
		
		for(i=0; i<bitrate/2; i++) {
			if(spectre[i]>=c) {
				printf("#");
				spectre[i]=c-1;
			}
			
			else 
				printf(" ");
		}
		
		printf("\n");
		l--;
		c--;
	}
}

//main temps d'execution
//int main(int argc, char const *argv[]){
  //sf_count_t readcount;
  //double data[BUFFER_LEN];
  //double complex datac[BUFFER_LEN];
  //double complex spectre[BUFFER_LEN];
  //double complex resultDFT[BUFFER_LEN];
  //clock_t start, stop;
  //double tempsDFT, tempsFFTiterTW, tempsFFTrec;

  //infile = sf_open("lapur.wav", SFM_READ, &sfinfo);
  //if (infile==NULL){
      //printf("Fichier non lisible\n");
      //sf_perror(NULL);
      //return 1;
  //}
  //while((readcount = sfx_mix_mono_read_double(infile, data, BUFFER_LEN))>0){
      //for(int i=0;i<BUFFER_LEN;i++){
	  //datac[i] = convDoubleComplex(data[i]);
      //}
      //start = clock();
      //twiddle(TW,BUFFER_LEN);
      //fftiterTW(datac,BUFFER_LEN,LOG2);
      //stop = clock();
      //tempsFFTiterTW = (double)(stop-start)/(double)(CLOCKS_PER_SEC);
      //start = clock();
      //fftrec(datac,spectre,BUFFER_LEN,LOG2);
      //stop = clock();
      //tempsFFTrec = (double)(stop-start)/(double)(CLOCKS_PER_SEC);
      //start = clock();
      //dft(datac,resultDFT,BUFFER_LEN);
      //stop = clock();
      //tempsDFT = (double)(stop-start)/(double)(CLOCKS_PER_SEC);
  //}
  //sf_close(infile);
  //printf("FFT intérative + TW : %lf ms\n",tempsFFTiterTW*1000);
  //printf("FFT récursive : %lf ms\n",tempsFFTrec*1000);
  //printf("DFT : %lf ms\n",tempsDFT*1000);
  //return 0;
//}

// main affiche spectre
int main(int argc, char const *argv[]) {
	  
	if(argc != 2) {
		fprintf(stderr, "Erreur parametres\n");
		return 1;
	}
	  
	infile = sf_open(argv[1], SFM_READ, &sfinfo);
	
	if(infile==NULL) {
		fprintf(stderr, "Fichier introuvable\n");
		sf_perror(NULL);
		return 1;
	}
	  
	sf_count_t readcount;
	double spectre[bitrate], data[BUFFER_LEN], maximum, moyenne, min, max, samplerate = sfinfo.samplerate;
	double complex datac[BUFFER_LEN];
	int place, ech, i, j;
	  
	twiddle(TW,BUFFER_LEN);
		 
	while((readcount = sfx_mix_mono_read_double(infile, data, BUFFER_LEN))>0) {
		for(i=0; i<BUFFER_LEN; i++)
			datac[i] = convDoubleComplex(data[i]);
		  
		fftiterTW(datac,BUFFER_LEN,LOG2);
		ech = 0;
		moyenne = 0;
		place = 0;
		for(i=0; i<BUFFER_LEN; i++) {
			moyenne += complexModule(datac[i]);
			ech++;
			if(ech==BUFFER_LEN/bitrate) {
				spectre[place] = 20*log10(moyenne/ech);
				ech = 0;
				place++;
				moyenne = 0;
			}
				
			if(place == bitrate) {
				min = spectre[0];
				for(j=0; j < bitrate; j++)
					if(spectre[j] < min)
						min = spectre[j];
			}
		}
		  
		max = spectre[0];
		  
		for(j=0; j<bitrate/2; j++) {
			spectre[j] = spectre[j] - min;
			if(max < spectre[j]) {
				max = spectre[j];
				maximum = j;
			}
		}
			
		system("clear");
		afficherSpectre(spectre, max);
		printf("%f\n", (maximum*samplerate)/bitrate);
		stop(23);
	}
	  
	sf_close(infile);
	  
	return 0;
}
