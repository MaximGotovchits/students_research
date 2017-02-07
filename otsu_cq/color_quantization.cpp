#include <math.h>
#include <iostream>
#include <cstdlib>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

using namespace std;
using namespace cv;

#define MAXCOLORS  65536

#define RED       2
#define GREEN     1
#define BLUE      0

#define BOX_SIZE  33

struct box {
    int r0;                      /* min value, exclusive */
    int r1;                      /* max value, inclusive */
    int g0;
    int g1;
    int b0;
    int b1;
    int vol;
};

float           m2[BOX_SIZE][BOX_SIZE][BOX_SIZE];
long int        wt[BOX_SIZE][BOX_SIZE][BOX_SIZE];
long int        mr[BOX_SIZE][BOX_SIZE][BOX_SIZE];
long int        mg[BOX_SIZE][BOX_SIZE][BOX_SIZE];
long int        mb[BOX_SIZE][BOX_SIZE][BOX_SIZE];

void InitBoxes(void)
{
    int j;
    int m;
    int n;
    
    for (j = 0; j < BOX_SIZE; j++)
    {
        for (m = 0; m < BOX_SIZE; m++)
        {
            for (n = 0; n < BOX_SIZE; n++)
            {
                m2[j][m][n] = 0.0;
                wt[j][m][n] = 0;
                mr[j][m][n] = 0;
                mg[j][m][n] = 0;
                mb[j][m][n] = 0;
            }
        }
    }
}

void Hist3d(int *Ir, int *Ig, int *Ib, int num_pixels,
       long int *vwt, long int *vmr, long int *vmg, long int *vmb, float *m2,
       long int *Qadd)
{
    int ind;
    int r;
    int g;
    int b;
    int inr;
    int ing;
    int inb;
    int table[256];
    long int i;
    
    for(i=0; i<256; ++i)
    {
        table[i]=i * i;
    }
    
    for(i=0; i<num_pixels; ++i)
    {
        r = Ir[i];
        g = Ig[i];
        b = Ib[i];
        inr=(r>>3)+1;
        ing=(g>>3)+1;
        inb=(b>>3)+1;
        Qadd[i]=ind=(inr<<10)+(inr<<6)+inr+(ing<<5)+ing+inb;
        /* [inr][ing][inb] */
        ++vwt[ind];
        vmr[ind] += r;
        vmg[ind] += g;
        vmb[ind] += b;
        m2[ind] += (float)(table[r]+table[g]+table[b]);
    }
    
}

/* At conclusion of the histogram step, we can interpret
 *   wt[r][g][b] = sum over voxel of P(c)
 *   mr[r][g][b] = sum over voxel of r*P(c)  ,  similarly for mg, mb
 *   m2[r][g][b] = sum over voxel of c^2*P(c)
 * Actually each of these should be divided by 'NumPixels' to give the usual
 * interpretation of P() as ranging from 0 to 1, but we needn't do that here.
 */

/* We now convert histogram into moments so that we can rapidly calculate
 * the sums of the above quantities over any desired box.
 */

/* compute cumulative moments. */
void M3d(long int *vwt, long int *vmr, long int *vmg, long int *vmb, float *m2)
{
    long int ind1;
    long int ind2;
    int i;
    int r;
    int g;
    int b;
    long int line;
    long int line_r;
    long int line_g;
    long int line_b;
    long int area[BOX_SIZE];
    long int area_r[BOX_SIZE];
    long int area_g[BOX_SIZE];
    long int area_b[BOX_SIZE];
    float line2;
    float area2[BOX_SIZE];
    
    for(r=1; r<BOX_SIZE; ++r)
    {
        for(i=0; i<BOX_SIZE; ++i)
        {
            area2[i] = 0.0;
            area[i] = 0;
            area_r[i] = 0;
            area_g[i] = 0;
            area_b[i] = 0;
        }
        
        for(g=1; g<=32; ++g)
        {
            line2 = 0.0;
            line = 0;
            line_r = 0;
            line_g = 0;
            line_b = 0;
            for(b=1; b<=32; ++b)
            {
                ind1 = (r<<10) + (r<<6) + r + (g<<5) + g + b; /* [r][g][b] */
                line += vwt[ind1];
                line_r += vmr[ind1];
                line_g += vmg[ind1];
                line_b += vmb[ind1];
                line2 += m2[ind1];
                area[b] += line;
                area_r[b] += line_r;
                area_g[b] += line_g;
                area_b[b] += line_b;
                area2[b] += line2;
                ind2 = ind1 - 1089; /* [r-1][g][b] */
                vwt[ind1] = vwt[ind2] + area[b];
                vmr[ind1] = vmr[ind2] + area_r[b];
                vmg[ind1] = vmg[ind2] + area_g[b];
                vmb[ind1] = vmb[ind2] + area_b[b];
                m2[ind1] = m2[ind2] + area2[b];
            }
        }
    }
}

/* Compute sum over a box of any given statistic */
long int Vol(struct box *cube, long int mmt[BOX_SIZE][BOX_SIZE][BOX_SIZE])
{
    return( mmt[cube->r1][cube->g1][cube->b1]
           -mmt[cube->r1][cube->g1][cube->b0]
           -mmt[cube->r1][cube->g0][cube->b1]
           +mmt[cube->r1][cube->g0][cube->b0]
           -mmt[cube->r0][cube->g1][cube->b1]
           +mmt[cube->r0][cube->g1][cube->b0]
           +mmt[cube->r0][cube->g0][cube->b1]
           -mmt[cube->r0][cube->g0][cube->b0] );
}

/* The next two routines allow a slightly more efficient calculation
 * of Vol() for a proposed subbox of a given box.  The sum of Top()
 * and Bottom() is the Vol() of a subbox split in the given direction
 * and with the specified new upper bound.
 */

/* Compute part of Vol(cube, mmt) that doesn't depend on r1, g1, or b1 */
/* (depending on dir) */
long int Bottom(struct box *cube, int dir,
                long int mmt[BOX_SIZE][BOX_SIZE][BOX_SIZE])
{
    switch(dir)
    {
        case RED:
            return( -mmt[cube->r0][cube->g1][cube->b1]
                   +mmt[cube->r0][cube->g1][cube->b0]
                   +mmt[cube->r0][cube->g0][cube->b1]
                   -mmt[cube->r0][cube->g0][cube->b0] );
            break;
        case GREEN:
            return( -mmt[cube->r1][cube->g0][cube->b1]
                   +mmt[cube->r1][cube->g0][cube->b0]
                   +mmt[cube->r0][cube->g0][cube->b1]
                   -mmt[cube->r0][cube->g0][cube->b0] );
            break;
        case BLUE:
            return( -mmt[cube->r1][cube->g1][cube->b0]
                   +mmt[cube->r1][cube->g0][cube->b0]
                   +mmt[cube->r0][cube->g1][cube->b0]
                   -mmt[cube->r0][cube->g0][cube->b0] );
            break;
        default:
            cerr << "Internal error: unrecognized value for dir." << endl;
    }
    return -1;
}

/* Compute remainder of Vol(cube, mmt), substituting pos for */
/* r1, g1, or b1 (depending on dir) */
long int Top(struct box *cube, int dir, int pos,
             long int mmt[BOX_SIZE][BOX_SIZE][BOX_SIZE])
{
    switch(dir)
    {
        case RED:
            return( mmt[pos][cube->g1][cube->b1]
                   -mmt[pos][cube->g1][cube->b0]
                   -mmt[pos][cube->g0][cube->b1]
                   +mmt[pos][cube->g0][cube->b0] );
            break;
        case GREEN:
            return( mmt[cube->r1][pos][cube->b1]
                   -mmt[cube->r1][pos][cube->b0]
                   -mmt[cube->r0][pos][cube->b1]
                   +mmt[cube->r0][pos][cube->b0] );
            break;
        case BLUE:
            return( mmt[cube->r1][cube->g1][pos]
                   -mmt[cube->r1][cube->g0][pos]
                   -mmt[cube->r0][cube->g1][pos]
                   +mmt[cube->r0][cube->g0][pos] );
            break;
        default:
            cerr << "Internal error: unrecognized value for dir." << endl;
    }
    return -1;
}

/* Compute the weighted variance of a box */
/* NB: as with the raw statistics, this is really the variance * NumPixels */
float Var(struct box *cube)
{
    float dr;
    float dg;
    float db;
    float xx;
    float result;
    
    dr = (float) Vol(cube, mr);
    dg = (float) Vol(cube, mg);
    db = (float) Vol(cube, mb);
    xx = m2[cube->r1][cube->g1][cube->b1]
    -m2[cube->r1][cube->g1][cube->b0]
    -m2[cube->r1][cube->g0][cube->b1]
    +m2[cube->r1][cube->g0][cube->b0]
    -m2[cube->r0][cube->g1][cube->b1]
    +m2[cube->r0][cube->g1][cube->b0]
    +m2[cube->r0][cube->g0][cube->b1]
    -m2[cube->r0][cube->g0][cube->b0];
    
    result = xx - (dr*dr+dg*dg+db*db)/(float)Vol(cube,wt);
    return (float) fabs((float) result);
}

/* We want to minimize the sum of the variances of two subboxes.
 * The sum(c^2) terms can be ignored since their sum over both subboxes
 * is the same (the sum for the whole box) no matter where we split.
 * The remaining terms have a minus sign in the variance formula,
 * so we drop the minus sign and MAXIMIZE the sum of the two terms.
 */


float Maximize(struct box *cube, int dir, int first, int last, int *cut,
               long int whole_r, long int whole_g, long int whole_b,
               long int whole_w)
{
    long int half_r;
    long int half_g;
    long int half_b;
    long int half_w;
    long int base_r;
    long int base_g;
    long int base_b;
    long int base_w;
    int i;
    float temp;
    float max;
    
    base_r = Bottom(cube, dir, mr);
    base_g = Bottom(cube, dir, mg);
    base_b = Bottom(cube, dir, mb);
    base_w = Bottom(cube, dir, wt);
    max = 0.0;
    *cut = -1;
    for(i=first; i<last; ++i)
    {
        half_r = base_r + Top(cube, dir, i, mr);
        half_g = base_g + Top(cube, dir, i, mg);
        half_b = base_b + Top(cube, dir, i, mb);
        half_w = base_w + Top(cube, dir, i, wt);
        /* now half_x is sum over lower half of box, if split at i */
        if (half_w == 0)
        {
            /* subbox could be empty of pixels! */
            /* never split into an empty box */
            continue;
        }
        else
        {
            temp = ((float)half_r*half_r + (float)half_g*half_g +
                    (float)half_b*half_b)/half_w;
        }
        
        half_r = whole_r - half_r;
        half_g = whole_g - half_g;
        half_b = whole_b - half_b;
        half_w = whole_w - half_w;
        if (half_w == 0)
        {
            /* subbox could be empty of pixels! */
            /* never split into an empty box */
            continue;
        }
        else
        {
            temp += ((float)half_r*half_r + (float)half_g*half_g +
                     (float)half_b*half_b)/half_w;
        }
        
        if (temp > max)
        {
            max=temp;
            *cut=i;
        }
    }
    
    return(max);
}

int Cut(struct box *set1, struct box *set2)
{
    int dir;
    int cutr;
    int cutg;
    int cutb;
    float maxr;
    float maxg;
    float maxb;
    long int whole_r;
    long int whole_g;
    long int whole_b;
    long int whole_w;
    
    whole_r = Vol(set1, mr);
    whole_g = Vol(set1, mg);
    whole_b = Vol(set1, mb);
    whole_w = Vol(set1, wt);
    
    maxr = Maximize(set1, RED, set1->r0+1, set1->r1, &cutr,
                    whole_r, whole_g, whole_b, whole_w);
    maxg = Maximize(set1, GREEN, set1->g0+1, set1->g1, &cutg,
                    whole_r, whole_g, whole_b, whole_w);
    maxb = Maximize(set1, BLUE, set1->b0+1, set1->b1, &cutb,
                    whole_r, whole_g, whole_b, whole_w);
    
    if( (maxr>=maxg)&&(maxr>=maxb) )
    {
        dir = RED;
        if (cutr < 0)
        {
            return 0; /* can't split the box */
        }
    }
    else
    {
        if( (maxg>=maxr)&&(maxg>=maxb) ) {
            dir = GREEN;
        }
        else
        {
            dir = BLUE;
        }
    }
    
    set2->r1 = set1->r1;
    set2->g1 = set1->g1;
    set2->b1 = set1->b1;
    
    switch (dir)
    {
        case RED:
            set2->r0 = set1->r1 = cutr;
            set2->g0 = set1->g0;
            set2->b0 = set1->b0;
            break;
        case GREEN:
            set2->g0 = set1->g1 = cutg;
            set2->r0 = set1->r0;
            set2->b0 = set1->b0;
            break;
        case BLUE:
            set2->b0 = set1->b1 = cutb;
            set2->r0 = set1->r0;
            set2->g0 = set1->g0;
            break;
        default:
            cerr << "Internal error: unrecognized value for dir." << endl;
    }
    
    set1->vol=(set1->r1-set1->r0)*(set1->g1-set1->g0)*(set1->b1-set1->b0);
    set2->vol=(set2->r1-set2->r0)*(set2->g1-set2->g0)*(set2->b1-set2->b0);
    
    return 1;
}


void Mark(struct box *cube, int label, int *tag)
{
    int r;
    int g;
    int b;
    
    for(r=cube->r0+1; r<=cube->r1; ++r)
    {
        for(g=cube->g0+1; g<=cube->g1; ++g)
        {
            for(b=cube->b0+1; b<=cube->b1; ++b)
            {
                tag[(r<<10) + (r<<6) + r + (g<<5) + g + b] = label;
            }
        }
    }
}




int quantize_color(int *Ir, int *Ig, int *Ib, int num_pixels,
                   int *lut_r, int *lut_g, int *lut_b, int num_colors,
                   long int *Qadd, bool compute_output_image)
{
    struct box *cube;
    int *tag;
    int next;
    long int i;
    long int weight;
    int k;
    float *vv;
    float temp;
    
    cube = (struct box *) calloc(MAXCOLORS, sizeof(*cube));
    vv = (float *) calloc(MAXCOLORS, sizeof(*vv));
    
    InitBoxes();
    
    Hist3d(Ir, Ig, Ib, num_pixels, (long int *) &wt, (long int *) &mr,
           (long int *) &mg, (long int *) &mb, (float *) &m2, Qadd);
    
    M3d((long int *) &wt, (long int *) &mr, (long int *) &mg,
        (long int *) &mb, (float *) &m2);
    
    cube[0].r0 = cube[0].g0 = cube[0].b0 = 0;
    cube[0].r1 = cube[0].g1 = cube[0].b1 = 32;
    next = 0;
    for(i=1; i<num_colors; ++i)
    {
        if (Cut(&cube[next], &cube[i]))
        {
            /* volume test ensures we won't try to cut one-cell box */
            vv[next] = (cube[next].vol>1) ? Var(&cube[next]) : 0.0F;
            vv[i] = (cube[i].vol>1) ? Var(&cube[i]) : 0.0F;
        }
        else
        {
            vv[next] = 0.0;   /* don't try to split this box again */
            i--;              /* didn't create box i */
        }
        next = 0;
        temp = vv[0];
        for(k=1; k<=i; ++k)
        {
            if (vv[k] > temp) {
                temp = vv[k]; next = k;
            }
        }
        if (temp <= 0.0)
        {
            num_colors = i+1;
            /* Only got num_colors boxes */
            break;
        }
    }
    
    tag = (int *) malloc(BOX_SIZE*BOX_SIZE*BOX_SIZE * sizeof(*tag));
    
    for(k=0; k<num_colors; ++k)
    {
        Mark(&cube[k], k, tag);
        weight = Vol(&cube[k], wt);
        if (weight)
        {
            lut_r[k] = (int) (Vol(&cube[k], mr) / weight);
            lut_g[k] = (int) (Vol(&cube[k], mg) / weight);
            lut_b[k] = (int) (Vol(&cube[k], mb) / weight);
        }
        else
        {
            /* bogus box */
            lut_r[k] = lut_g[k] = lut_b[k] = 0;
        }
    }
    
    if (compute_output_image)
    {
        for (i=0; i<num_pixels; ++i)
        {
            Qadd[i] = tag[Qadd[i]];
        }
    }
    
    free(tag);
    free(vv);
    free(cube);
    
    return(num_colors);
}


Mat quantize(string image_name, int colors_number) {
    int num_pixels;
    int *image_bytes_red;
    int *image_bytes_green;
    int *image_bytes_blue;
    int *lut_red;
    int *lut_green;
    int *lut_blue;
    long int **out_image;
    long int *out_pr;
    int num_colors = colors_number;
    Mat input_image = imread(image_name);
    vector<Mat> channel(3);
    split(input_image, channel);
    num_pixels = input_image.rows * input_image.cols;
    image_bytes_red = (int*)malloc(sizeof(int) * num_pixels);
    image_bytes_green = (int*)malloc(sizeof(int) * num_pixels);
    image_bytes_blue = (int*)malloc(sizeof(int) * num_pixels);
    for (int j = 0; j < input_image.cols; ++j) {
        for (int i = 0; i < input_image.rows; ++i) {
            image_bytes_red[j + (i * input_image.cols)] = (int)channel[0].at<uchar>(i, j);
            image_bytes_green[j + (i * input_image.cols)] = (int)channel[1].at<uchar>(i, j);
            image_bytes_blue[j + (i * input_image.cols)] = (int)channel[2].at<uchar>(i, j);
        }
    }
    lut_red = (int *) malloc(num_colors * sizeof(*lut_red));
    lut_green = (int *) malloc(num_colors * sizeof(*lut_green));
    lut_blue = (int *) malloc(num_colors * sizeof(*lut_blue));
    out_image = (long int **)malloc(sizeof(long int *) * input_image.rows);
    for (int i = 0; i < input_image.rows; ++i) {
        out_image[i] = (long int*)malloc(sizeof(long int) * input_image.cols);
    }
    out_pr = (long int*)out_image[0];
    num_colors = quantize_color(image_bytes_red, image_bytes_green,
                                image_bytes_blue, num_pixels,
                                lut_red, lut_green, lut_blue, num_colors,
                                out_pr, true);
    Mat output_image;
    for (int i = 0; i < num_pixels; ++i) {
        int row = i / (num_pixels / input_image.rows);
        int col = i % (num_pixels / input_image.rows);   
        channel[0].at<uchar>(row, col) = lut_red[out_image[row][col]];
        channel[1].at<uchar>(row, col) = lut_green[out_image[row][col]];
        channel[2].at<uchar>(row, col) = lut_blue[out_image[row][col]];
    }
    merge(channel, output_image);
    imwrite("/Users/IIMaximII/Desktop/out1.png", output_image);
    free(lut_red);
    free(lut_green);
    free(lut_blue);
    return output_image;
}


int main() {
    string image_name = "/Users/IIMaximII/Desktop/imageQQ.png";
    quantize(image_name, 8);
    return 0;
}
