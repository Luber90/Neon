/* arm-linux-gnueabihf-gcc main.c -static -mfpu=neon -ftree-vectorize -o main */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <arm_neon.h>

struct rgb
{
  char r;
  char g;
  char b;
};

struct image
{
  char p;
  int format;
  int width;
  int height;
  int intensity;
  //struct rgb **rgb;
  unsigned char *pixels;
};

double dot(double* a, double* b, int n){
  double result = 0;
  for(int i = 0; i < n; ++i){
    result += a[i]*b[i];
  }
}


void loadImg(uint8_t* name, unsigned char* img_pixels){
  struct image m;
  FILE *fp = fopen (name, "r");
  if (fp == 0)
    return;
  fscanf (fp, "%c%d\n", &m.p, &m.format);
  fscanf (fp, "%d %d\n", &m.width, &m.height);
  fscanf (fp, "%d\n", &m.intensity);
  //printf("%d %d\n", m.width, m.height);

  fread (img_pixels, sizeof (*img_pixels), m.width * m.height * 3, fp);
  fclose(fp);
}

void saveImg(uint8_t* name, unsigned char* img_pixels, int width, int height){
  FILE *fp = fopen (name, "w");
  if (fp == 0)
    return;
  fprintf (fp, "P5\n");
  fprintf (fp, "%d %d\n", width, height);
  fprintf (fp, "%d\n", 255);
  fwrite (img_pixels, sizeof (*img_pixels), width * height, fp);
  fclose (fp);
}

void reorder(uint8_t* bg_pivot, uint8_t width, uint8_t height, int nGauss,
             double* sigmas, double* omegas, double* mus,
             double bg_thresh){
  
  for(uint8_t i = 0; i < width; ++i){
    for(uint8_t j = 0; j < height; ++j){
      int changed = 0;
      double* ratios = (double*)malloc(nGauss*sizeof(double));
      for(int k = 0; k < nGauss; ++k){
        ratios[k] = omegas[(i*height*nGauss)+(j*nGauss)+k]/sqrt(sigmas[(i*height*nGauss)+(j*nGauss)+k]);
      }
      int* indices = (int*) malloc(nGauss*sizeof(int));
      for(int x = 0; x < nGauss; ++x) indices[x]=-1;
      int found = 0;
      while(found < nGauss){
        int max = -1;
        int index;
        for(int x = 1; x < nGauss; ++x){
          if(ratios[x]>max&&indices[x]==-1){
            max = ratios[x];
            index = x;
          }
        }
        indices[found] = index;
        found++;
      }
      for(int x = 0; x < nGauss; ++x){
        mus[(i*height*nGauss*3)+(j*nGauss*3)+(x*3)+0] = 
          mus[(i*height*nGauss*3)+(j*nGauss*3)+(indices[x]*3)+0];
        mus[(i*height*nGauss*3)+(j*nGauss*3)+(x*3)+1] = 
          mus[(i*height*nGauss*3)+(j*nGauss*3)+(indices[x]*3)+1];
        mus[(i*height*nGauss*3)+(j*nGauss*3)+(x*3)+2] = 
          mus[(i*height*nGauss*3)+(j*nGauss*3)+(indices[x]*3)+2];

        sigmas[(i*height*nGauss)+(j*nGauss)+x] =  
          sigmas[(i*height*nGauss)+(j*nGauss)+indices[x]];

        omegas[(i*height*nGauss)+(j*nGauss)+x] =  
          omegas[(i*height*nGauss)+(j*nGauss)+indices[x]];
      }
        double cummProb = 0.0;
        for(int x = 0; x < nGauss; ++x){
          cummProb+=omegas[(i*height*nGauss)+(j*nGauss)+x];
          if (cummProb>=bg_thresh){
              bg_pivot[i*height+j]=x;
              changed = 1;
              break;
          }
        }
        if(!changed){
          bg_pivot[i*height+j]=nGauss-2;
        }
    }
  }
}

void updateParams(uint8_t* labels, int width, int height, int nGauss,
                  double* sigmas, double* omegas, double* mus,
                  double bg_thresh, double lr, uint8_t* pivot_bg,
                  unsigned char* img){
  for(int i = 0; i < width*height; ++i) labels[i] = 0;
  for(int i = 0; i < width; ++i){
    for(int j = 0; j < height; ++j){
      unsigned char* pixel = (unsigned char*)malloc(3*sizeof(unsigned char));
      pixel[0] = img[i*height*3+j*3+0];
      pixel[1] = img[i*height*3+j*3+1];
      pixel[2] = img[i*height*3+j*3+2];
      int match = -1;
      for(int k = 0; k < nGauss; ++k){
        double Xmu[3] = {0.0};
        for(int x = 0; x < 3; ++x){
          Xmu[x] = pixel[x] - mus[(i*height*nGauss*3)+(j*nGauss*3)+(k*3)+x];
        }
        double tmp = 1/sigmas[(i*height*nGauss)+(j*nGauss)+k];
        double toDot[3] = {tmp*Xmu[0], tmp*Xmu[1], tmp*Xmu[2]};
        double dist = dot(Xmu, toDot, 3);

        //printf("%g %g\n", tmp, dot(toDot, Xmu,  3));
        /*if(i*j > ((width*height)/2)){
          printf("%d\n", i*j);
          printf("%d\n", width);
          printf("%d %d %d\n", pixel[0], pixel[1], pixel[2]);
          printf("%d %d %d\n", Xmu[0], Xmu[1], Xmu[2]);
          printf("%g\n", dist);
          while(1){}
        }*/
        //printf("%g\n", dist);
        if(dist < 6.25*sigmas[(i*height*nGauss)+(j*nGauss)+k]){
          match = k;
          break;
        }
      }
      if(match!=-1){
        omegas[(i*height*nGauss)+(j*nGauss)+0] = (1.0 - lr)*omegas[(i*height*nGauss)+(j*nGauss)+0];
        omegas[(i*height*nGauss)+(j*nGauss)+1] = (1.0 - lr)*omegas[(i*height*nGauss)+(j*nGauss)+1];
        omegas[(i*height*nGauss)+(j*nGauss)+2] = (1.0 - lr)*omegas[(i*height*nGauss)+(j*nGauss)+2];
        omegas[(i*height*nGauss)+(j*nGauss)+match] += lr;

        double rho = 0.0001;
        double tmp[3] = {pixel[0] - mus[(i*height*nGauss*3)+(j*nGauss*3)+(match*3)+0],
                        pixel[1] - mus[(i*height*nGauss*3)+(j*nGauss*3)+(match*3)+1],
                        pixel[2] - mus[(i*height*nGauss*3)+(j*nGauss*3)+(match*3)+2]};
        omegas[(i*height*nGauss)+(j*nGauss)+match] = rho*omegas[(i*height*nGauss)+(j*nGauss)+match]+
                                                       rho*dot(tmp, tmp, 3);
        mus[(i*height*nGauss*3)+(j*nGauss*3)+(match*3)+0] = rho*mus[(i*height*nGauss*3)+(j*nGauss*3)+(match*3)+0]+
                                                            rho*pixel[0];
        mus[(i*height*nGauss*3)+(j*nGauss*3)+(match*3)+1] = rho*mus[(i*height*nGauss*3)+(j*nGauss*3)+(match*3)+1]+
                                                            rho*pixel[1];
        mus[(i*height*nGauss*3)+(j*nGauss*3)+(match*3)+2] = rho*mus[(i*height*nGauss*3)+(j*nGauss*3)+(match*3)+2]+
                                                            rho*pixel[2];
        if(match > pivot_bg[i*height+j]){
          labels[i*height+j] = 255;
        }
      }
      else{
        mus[(i*height*nGauss*3)+(j*nGauss*3)+((nGauss-1)*3)+0] = pixel[0];
        mus[(i*height*nGauss*3)+(j*nGauss*3)+((nGauss-1)*3)+1] = pixel[1];
        mus[(i*height*nGauss*3)+(j*nGauss*3)+((nGauss-1)*3)+2] = pixel[2];
        labels[i*height+j] = 255;
      }

    }
  }

}

void MOG(int width, int height, int n, unsigned char* name){
  //setup
  int numOfGauss = 3;
  double BG_thresh = 0.6;
  double learning_rate = 0.01;
  double* mus;
  uint8_t* bg_pivots = (uint8_t*)malloc(width*height*sizeof(u_int8_t));
  uint8_t* labels = (uint8_t*)malloc(width*height*sizeof(u_int8_t));
  mus = (double*) malloc (width * height * numOfGauss*3*sizeof (double));
  
  double* sigma2;
  sigma2 = (double*) malloc (width * height * numOfGauss * sizeof (double));
  
  double* omegas;
  omegas = (double*) malloc (width * height * numOfGauss * sizeof (double));

  for(int i = 0; i < width; ++i){
    for(int j = 0; j < height; ++j){
      mus[(i*height*numOfGauss*3)+(j*numOfGauss*3)+(0*3)+0] = 122.0;
      mus[(i*height*numOfGauss*3)+(j*numOfGauss*3)+(0*3)+1] = 122.0;
      mus[(i*height*numOfGauss*3)+(j*numOfGauss*3)+(0*3)+2] = 122.0;
      mus[(i*height*numOfGauss*3)+(j*numOfGauss*3)+(1*3)+0] = 122.0;
      mus[(i*height*numOfGauss*3)+(j*numOfGauss*3)+(1*3)+1] = 122.0;
      mus[(i*height*numOfGauss*3)+(j*numOfGauss*3)+(1*3)+2] = 122.0;
      mus[(i*height*numOfGauss*3)+(j*numOfGauss*3)+(2*3)+0] = 122.0;
      mus[(i*height*numOfGauss*3)+(j*numOfGauss*3)+(2*3)+1] = 122.0;
      mus[(i*height*numOfGauss*3)+(j*numOfGauss*3)+(2*3)+2] = 122.0;

      sigma2[(i*height*numOfGauss)+(j*numOfGauss)+0] = 36.0;
      sigma2[(i*height*numOfGauss)+(j*numOfGauss)+1] = 36.0;
      sigma2[(i*height*numOfGauss)+(j*numOfGauss)+2] = 36.0;

      omegas[(i*height*numOfGauss)+(j*numOfGauss)+0] = 1.0/numOfGauss;
      omegas[(i*height*numOfGauss)+(j*numOfGauss)+1] = 1.0/numOfGauss;
      omegas[(i*height*numOfGauss)+(j*numOfGauss)+2] = 1.0/numOfGauss;
    }
  }
  for(int i = 0; i < n+1; i++){
    char *num;
    char buffer[20];
    char filename[50] = "moving/";
    strcat(filename, name);
    snprintf(buffer, 20, "%d", i);
    strcat(filename, buffer);
    strcat(filename, ".ppm");
    printf("%s\n", filename);
    unsigned char* inpixels = (unsigned char *) malloc (width * height * 3 *
			      sizeof (unsigned char));
    loadImg(filename, inpixels);
    reorder(bg_pivots, width, height, numOfGauss, sigma2, omegas, mus, BG_thresh);
    updateParams(labels, width, height, numOfGauss, sigma2, omegas,
                mus, BG_thresh, 0.01, bg_pivots, inpixels);
    char outname[50] = "out/out";
    strcat(outname, buffer);
    strcat(outname, ".ppm");
    saveImg(outname, labels, width, height);
  }


}

void neon_grayscale (uint8_t * dest, uint8_t * src, int num);
int main (int argc, char* argv[])
{
  int width = atoi(argv[1]);
  int height = atoi(argv[2]);
  int n = atoi(argv[3]);
  uint8_t* name = argv[4];

  unsigned char* inpixels,* outpixels;

  inpixels =
    (unsigned char *) malloc (width * height * 3 *
			      sizeof (unsigned char));
  loadImg(name, inpixels);

  outpixels =
    (unsigned char *) malloc (width * height * sizeof (unsigned char));


  neon_grayscale (outpixels, inpixels, width * height);

  saveImg("out.ppm", outpixels, width, height);
  MOG(width, height, n, name);
  return 0;
}

void neon_grayscale (uint8_t * dest, uint8_t * src, int num)
{
  int i;
  uint8x8_t r_ratio = vdup_n_u8 (77);
  uint8x8_t g_ratio = vdup_n_u8 (151);
  uint8x8_t b_ratio = vdup_n_u8 (28);
  num /= 8;			// NEON will work on 8 pixels a time
  for (i = 0; i < num; i++)
    {
      uint16x8_t temp;
      uint8x8x3_t rgb = vld3_u8 (src);
      uint8x8_t result;
      temp = vmull_u8 (rgb.val[0], r_ratio);
      temp = vmlal_u8 (temp, rgb.val[1], g_ratio);
      temp = vmlal_u8 (temp, rgb.val[2], b_ratio);
      result = vshrn_n_u16 (temp, 8);
      vst1_u8 (dest, result);
      src += 8 * 3;		// 3 x 8 pixels in RGB format
      dest += 8;		// One single 8 - bit value per pixel
    }
}
