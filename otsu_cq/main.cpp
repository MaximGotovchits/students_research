#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <numeric>

using namespace cv;
using namespace std;

bool image_equals(const Mat& image1, const Mat& image2) {
    cout << countNonZero(image1 - image2) << endl;
    return ((countNonZero(image1 - image2) == 0) && (countNonZero(image2 - image1) == 0));
}

Mat binarize(Mat& image, int board) {
    for (int i = 0; i < image.rows; ++i) {
        for (int j = 0; j < image.cols; ++j) {
            if (image.at<uchar>(i, j) < board) {
                image.at<uchar>(i, j) = 0;
            } else {
                image.at<uchar>(i, j) = 255;
            }
        }
    }
    return image;
}

Mat otsu(Mat& image) {
    short board = 0;
    double first_class_mean = 0;
    double second_class_mean = 0;
    double first_class_prob;
    double second_class_prob;
    int channels[] = {0};
    float range[] = {0, 255};
    const float *ranges[] = {range};
    Mat histogram;
    int histSize = 256;
    calcHist(&image, 1, channels, Mat(), histogram, 1, &histSize, ranges, true, false);
    unsigned hist_sum = 0;
    hist_sum = accumulate(histogram.begin<uchar>(), histogram.end<uchar>(), 0);
    first_class_prob = (double) histogram.at<uchar>(0, 0) / hist_sum;
    second_class_prob = (double) 1 - first_class_prob;
    first_class_mean = 0;
    second_class_mean = (double) 1 - first_class_mean;
    double dispersion = (double) first_class_prob * second_class_prob * (first_class_mean - second_class_mean) * (first_class_mean - second_class_mean);
    double max_dispersion = dispersion;
    double first_class_sum = (double) histogram.at<uchar>(0, 0);
    double first_class_not_mean = 0;
    for (int j = 1; j < histogram.rows; j++) {
        first_class_sum += histogram.at<uchar>(0, j);
        first_class_prob = (double) first_class_sum / hist_sum;
        second_class_prob = (double) 1 - first_class_prob;
        first_class_not_mean += (j * histogram.at<uchar>(0, j));
        first_class_mean = (double) first_class_not_mean / first_class_prob;
        second_class_mean = 1 - first_class_mean;
        dispersion = (double) first_class_prob * second_class_prob * (first_class_mean - second_class_mean) * (first_class_mean - second_class_mean);
        if (dispersion > max_dispersion) {
            max_dispersion = dispersion;
            board = j;
        }
    }
    Mat result_image = binarize(image, board); // - светлее + темнее
    imwrite("new_image.jpg", result_image);
    return result_image;
}

int main(int argc, char** argv) {
    string image_name = "image.jpg";
    Mat image = imread(image_name, 0);
    Mat empty_image;
    if (!image.rows && !image.cols) {
        cerr << "Error: Image must be called \"image.jpg\"" << endl;
        return 1;
    }
    Mat opencv_otsu_result = imread("opencv_otsu_result.jpg");
    threshold(image, opencv_otsu_result, 0, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
    imwrite("opencv_otsu_result.jpg", opencv_otsu_result);
    Mat result = otsu(image);
    if (!image_equals(result, opencv_otsu_result)) {
        cerr << ":-(" << endl;
    }
    return 0;
}
