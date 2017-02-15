#include <cstdio>
#include <cstdlib>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>

using namespace std;
using namespace cv;

#define MAXCOLOR 256
#define	RED 2
#define	GREEN 1
#define BLUE 0

struct box {
    unsigned long r0;			 /* min value, exclusive */
    unsigned long r1;			 /* max value, inclusive */
    unsigned long g0;
    unsigned long g1;
    unsigned long b0;
    unsigned long b1;
    unsigned long vol;
};

/* Histogram is in elements 1..HISTSIZE along each axis,
 * element 0 is for base or marginal value
 * NB: these must start out 0!
 */

float m2[33][33][33];
long int wt[33][33][33], mr[33][33][33], mg[33][33][33], mb[33][33][33];
unsigned char *Ir, *Ig, *Ib;
int size; /*image size*/
int	K;    /*color look-up table size*/
unsigned short int *Qadd;

void Hist3d(long int* vwt, long int* vmr, long int* vmg, long int* vmb, float* m2)
/* build 3-D color histogram of counts, r/g/b, c^2 */
{
    int ind, r, g, b;
    int inr, ing, inb, table[256];
    long int i;
    
    for(i = 0; i < 256; ++i) {
        table[i]= i * i;
    }
    Qadd = (unsigned short int *)malloc(sizeof(short int)*size);
    if (Qadd==NULL) {
        printf("Not enough space\n");
        exit(1);
    }
    for(i=0; i<size; ++i) {
        r = Ir[i]; g = Ig[i]; b = Ib[i];
        inr=(r>>3)+1;
        ing=(g>>3)+1;
        inb=(b>>3)+1;
        Qadd[i]=ind=(inr<<10)+(inr<<6)+inr+(ing<<5)+ing+inb;
        /*[inr][ing][inb]*/
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
 * Actually each of these should be divided by 'size' to give the usual
 * interpretation of P() as ranging from 0 to 1, but we needn't do that here.
 */

/* We now convert histogram into moments so that we can rapidly calculate
 * the sums of the above quantities over any desired box.
 */


void M3d(long int* vwt, long int* vmr, long int* vmg, long int* vmb, float* m2) /* compute cumulative moments. */
{
    unsigned short int ind1, ind2;
    unsigned char i, r, g, b;
    long int line, line_r, line_g, line_b,
    area[33], area_r[33], area_g[33], area_b[33];
    float    line2, area2[33];
    
    for(r=1; r<=32; ++r){
        for(i=0; i<=32; ++i) {
            area2[i]=area[i]=area_r[i]=area_g[i]=area_b[i]=0;
        }
        for(g=1; g<=32; ++g){
            line2 = line = line_r = line_g = line_b = 0;
            for(b=1; b<=32; ++b){
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


long int Vol(struct box* cube, long int mmt[33][33][33])
/* Compute sum over a box of any given statistic */
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

long int Bottom(struct box* cube, unsigned char dir, long int mmt[33][33][33])
/* Compute part of Vol(cube, mmt) that doesn't depend on r1, g1, or b1 */
/* (depending on dir) */
{
    switch(dir){
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
    }
    printf("Buttom(...) error.");
    return -1;
}


long int Top(struct box* cube, unsigned char dir, int pos, long int mmt[33][33][33])
/* Compute remainder of Vol(cube, mmt), substituting pos for */
/* r1, g1, or b1 (depending on dir) */
{
    switch(dir){
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
    }
    printf("Top(...) error.");
    return -1;
}


float Var(struct box* cube)
/* Compute the weighted variance of a box */
/* NB: as with the raw statistics, this is really the variance * size */
{
    float dr, dg, db, xx;
    
    dr = Vol(cube, mr);
    dg = Vol(cube, mg);
    db = Vol(cube, mb);
    xx =  m2[cube->r1][cube->g1][cube->b1]
    -m2[cube->r1][cube->g1][cube->b0]
    -m2[cube->r1][cube->g0][cube->b1]
    +m2[cube->r1][cube->g0][cube->b0]
    -m2[cube->r0][cube->g1][cube->b1]
    +m2[cube->r0][cube->g1][cube->b0]
    +m2[cube->r0][cube->g0][cube->b1]
    -m2[cube->r0][cube->g0][cube->b0];
    double variance = xx - (dr*dr+dg*dg+db*db)/(float)Vol(cube,wt);
    cout << xx << " - " << (dr*dr+dg*dg+db*db)/(float)Vol(cube,wt) << endl;
    cout << "What f returns: " << variance << endl;
    cout << "Real variance: " << sqrt(variance / size) << endl;
    return xx - (dr*dr+dg*dg+db*db)/(float)Vol(cube,wt);
}

/* We want to minimize the sum of the variances of two subboxes.
 * The sum(c^2) terms can be ignored since their sum over both subboxes
 * is the same (the sum for the whole box) no matter where we split.
 * The remaining terms have a minus sign in the variance formula,
 * so we drop the minus sign and MAXIMIZE the sum of the two terms.
 */


float Maximize(struct box* cube, unsigned char dir, int first, int last, int* cut,
               long int whole_r, long int whole_g, long int whole_b, long int whole_w)
{
    long int half_r, half_g, half_b, half_w;
    long int base_r, base_g, base_b, base_w;
    int i;
    float temp, max;
    
    base_r = Bottom(cube, dir, mr);
    base_g = Bottom(cube, dir, mg);
    base_b = Bottom(cube, dir, mb);
    base_w = Bottom(cube, dir, wt);
    max = 0.0;
    *cut = -1;
    for(i=first; i<last; ++i){
        half_r = base_r + Top(cube, dir, i, mr);
        half_g = base_g + Top(cube, dir, i, mg);
        half_b = base_b + Top(cube, dir, i, mb);
        half_w = base_w + Top(cube, dir, i, wt);
        /* now half_x is sum over lower half of box, if split at i */
        if (half_w == 0) {      /* subbox could be empty of pixels! */
            continue;             /* never split into an empty box */
        } else
            temp = ((float)half_r*half_r + (float)half_g*half_g +
                    (float)half_b*half_b)/half_w;
            
            half_r = whole_r - half_r;
            half_g = whole_g - half_g;
            half_b = whole_b - half_b;
            half_w = whole_w - half_w;
            if (half_w == 0) {      /* subbox could be empty of pixels! */
                continue;             /* never split into an empty box */
            } else
                temp += ((float)half_r*half_r + (float)half_g*half_g +
                         (float)half_b*half_b)/half_w;
                
                if (temp > max) {
                    max=temp; *cut=i;
                }
    }
    return(max);
}

int Cut(struct box* set1, struct box* set2)
{
    unsigned char dir;
    int cutr, cutg, cutb;
    float maxr, maxg, maxb;
    long int whole_r, whole_g, whole_b, whole_w;
    
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
    
    if( (maxr>=maxg)&&(maxr>=maxb) ) {
        dir = RED;
        if (cutr < 0) {
            return 0; /* can't split the box */
        }
    }
    else {
        if( (maxg>=maxr)&&(maxg>=maxb) ) {
            dir = GREEN;
        }
        else {
            dir = BLUE;
        }
    }
    if (set1->r1 && set1->g1 && set1->b1) {
        set2->r1 = set1->r1;
        set2->g1 = set1->g1;
        set2->b1 = set1->b1;
    } else {
//        cerr << "Too many colors are chosen." << endl;
        return 0;
    }
    switch (dir) {
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
    }
    set1->vol=(set1->r1-set1->r0)*(set1->g1-set1->g0)*(set1->b1-set1->b0);
    set2->vol=(set2->r1-set2->r0)*(set2->g1-set2->g0)*(set2->b1-set2->b0);
    return 1;
}


void Mark(struct box* cube, int label, unsigned char* tag) {
    int r, g, b;
    
    for(r=cube->r0+1; r<=cube->r1; ++r) {
        for(g=cube->g0+1; g<=cube->g1; ++g) {
            for(b=cube->b0+1; b<=cube->b1; ++b) {
                tag[(r<<10) + (r<<6) + r + (g<<5) + g + b] = label;
            }
        }
    }
}


//vector<vector<vector<int>>> get_all_colors_number(Mat image) {
//    cout << 256*256*256 << endl;
//    vector<vector<vector<int>>> hist;
//    vector<int> colors_g(256, 0);
//    for (int i = 0; i < 256; ++i) {
//        vector<vector<int>> colors_r;
//        for (int j = 0; j < 256; ++j) {
//            colors_r.push_back(colors_g);
//        }
//        hist.push_back(colors_r);
//    }
//    int colors_number = 0;
//    vector<Mat> channel(3);
//    split(image, channel);
//    int rows = image.rows;
//    int cols = image.cols;
//    for (int i = 0; i < rows; ++i) {
//        for (int j = 0; j < cols; ++j) {
//            int r = (int)channel[0].at<uchar>(i, j);
//            int g = (int)channel[1].at<uchar>(i, j);
//            int b = (int)channel[2].at<uchar>(i, j);
//            if (hist[r][g][b] == 0) {
//                ++colors_number;
//            }
//            ++hist[r][g][b];
//        }
//    }
//    cout << "Number of colors: " << colors_number << endl;
//    return hist;
//}


Mat Quantize(string image_name, int colors_number) {
    Mat input_image = imread(image_name);
    vector<Mat> channel(3);
    split(input_image, channel);
    
//    auto g = get_all_colors_number(input_image);
    
    struct box cube[MAXCOLOR];
    unsigned char *tag;
    unsigned char lut_r[MAXCOLOR], lut_g[MAXCOLOR], lut_b[MAXCOLOR];
    int	next;
    long int i, weight;
    int	k;
    float vv[MAXCOLOR], temp;
    
    /* input R,G,B components into Ir, Ig, Ib;
     set size to width*height */
    
    
    size = channel[0].rows * channel[0].cols;
    Ir = (unsigned char*)malloc(sizeof(unsigned char) * size);
    Ig = (unsigned char*)malloc(sizeof(unsigned char) * size);
    Ib = (unsigned char*)malloc(sizeof(unsigned char) * size);
    
    for (unsigned i = 0; i < channel[0].rows; ++i) {
        for (unsigned j = 0; j < channel[0].cols; ++j) {
            Ir[j + (i * channel[0].cols)] = channel[0].at<uchar>(i, j);
            Ig[j + (i * channel[0].cols)] = channel[1].at<uchar>(i, j);
            Ib[j + (i * channel[0].cols)] = channel[2].at<uchar>(i, j);
        }
    }
    K = colors_number;
    
    clock_t begin = clock();
    
    Hist3d((long int *) &wt, (long int *) &mr, (long int *) &mg, (long int *) &mb, (float *) &m2);
    
//    printf("Histogram done\n");
    
    free(Ig); free(Ib); free(Ir);
    
    M3d((long int *) &wt, (long int *) &mr, (long int *) &mg,
        (long int *) &mb, (float *) &m2);
    
    cube[0].r0 = cube[0].g0 = cube[0].b0 = 0;
    cube[0].r1 = cube[0].g1 = cube[0].b1 = 32;
    next = 0;
    for(i=1; i<K; ++i){
        if (Cut(&cube[next], &cube[i])) {
            /* volume test ensures we won't try to cut one-cell box */
            vv[next] = (cube[next].vol>1) ? Var(&cube[next]) : 0.0;
            vv[i] = (cube[i].vol>1) ? Var(&cube[i]) : 0.0;
        } else {
            vv[next] = 0.0;   /* don't try to split this box again */
            i--;              /* didn't create box i */
        }
        next = 0; temp = vv[0];
        for(k=1; k<=i; ++k)
            if (vv[k] > temp) {
                temp = vv[k]; next = k;
            }
        if (temp <= 0.0) {
            K = i + 1;
            fprintf(stderr, "Only got %d boxes\n", K);
            break;
        }
    }
//    printf("Partition done\n");
    
    /* the space for array m2 can be freed now */
    
    tag = (unsigned char *)malloc(33*33*33);
    if (tag==NULL) {
        printf("Not enough space\n");
        exit(1);
    }
    for(k=0; k<K; ++k) {
        Mark(&cube[k], k, tag);
        weight = Vol(&cube[k], wt);
        if (weight) {
            lut_r[k] = Vol(&cube[k], mr) / weight;
            lut_g[k] = Vol(&cube[k], mg) / weight;
            lut_b[k] = Vol(&cube[k], mb) / weight;
        }
        else {
            fprintf(stderr, "bogus box %d\n", k);
            lut_r[k] = lut_g[k] = lut_b[k] = 0;
        }
    }
    
    for(i=0; i<size; ++i) {
        Qadd[i] = tag[Qadd[i]];
    }
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << elapsed_secs << endl;
    /* output lut_r, lut_g, lut_b as color look-up table contents,
     Qadd as the quantized image (array of table addresses). */
    
    
    Mat output_image;
    vector<Mat> output_channel;
    Mat out_r(channel[0].rows, channel[0].cols, CV_8UC1, Scalar(0));
    Mat out_g(channel[0].rows, channel[0].cols, CV_8UC1, Scalar(0));
    Mat out_b(channel[0].rows, channel[0].cols, CV_8UC1, Scalar(0));
    output_channel.push_back(out_r);
    output_channel.push_back(out_g);
    output_channel.push_back(out_b);
    unsigned row;
    unsigned col;
    for (unsigned i = 0; i < size; ++i) {
        row = i / channel[0].cols;
        col = i % channel[0].cols;
        output_channel[0].at<uchar>(row, col) = lut_r[Qadd[i]];
        output_channel[1].at<uchar>(row, col) = lut_g[Qadd[i]];
        output_channel[2].at<uchar>(row, col) = lut_b[Qadd[i]];
    }
//    cout << row << " " << col << endl;
    merge(output_channel, output_image);
    return output_image;
}


int main() {
//    for (int i = 1; i <= 5; ++i) {
//        Mat output_image = Quantize("/Users/IIMaximII/Desktop/mlmarerials/doc" + to_string(i) + ".jpg", 20);
//        imwrite("/Users/IIMaximII/Desktop/mlmarerials/out" + to_string(i) + ".jpg", output_image);
//    }
    Mat output_image = Quantize("/Users/IIMaximII/Desktop/mlmarerials/doc4.jpg", 4);
    imwrite("/Users/IIMaximII/Desktop/mlmarerials/out4.png", output_image);
    
}
