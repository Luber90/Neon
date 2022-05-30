/* arm-linux-gnueabihf-gcc main.c -static -mfpu=neon -ftree-vectorize -o main */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
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

void neon_grayscale (uint8_t * dest, uint8_t * src, int num);
int
main (int argc, char* argv[])
{
  struct image m;
  unsigned char *outpixels;
  FILE *fp = fopen (argv[1], "r");	//fopen(...);
  if (fp == 0)
    return -1;
  FILE *fpout = fopen (argv[2], "w");
  if (fpout == 0)
    {
      fclose (fp);
      return -1;
    }
  fscanf (fp, "%c%d\n", &m.p, &m.format);
  fscanf (fp, "%d %d\n", &m.width, &m.height);
  fscanf (fp, "%d\n", &m.intensity);

  printf ("%c%d\n", m.p, m.format);
  printf ("%d %d\n", m.width, m.height);
  printf ("%d\n", m.intensity);

  //allocate array to hold the pixels
  m.pixels =
    (unsigned char *) malloc (m.width * m.height * 3 *
			      sizeof (unsigned char));
  fread (m.pixels, sizeof (*m.pixels), m.width * m.height * 3, fp);

  outpixels =
    (unsigned char *) malloc (m.width * m.height * sizeof (unsigned char));

  neon_grayscale (outpixels, m.pixels, m.width * m.height);

  fprintf (fpout, "P5\n");
  fprintf (fpout, "%d %d\n", m.width, m.height);
  fprintf (fpout, "%d\n", m.intensity);
  fwrite (outpixels, sizeof (*outpixels), m.width * m.height, fpout);
  fclose (fp);
  fclose (fpout);
  return 0;
}

void
neon_grayscale (uint8_t * dest, uint8_t * src, int num)
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
