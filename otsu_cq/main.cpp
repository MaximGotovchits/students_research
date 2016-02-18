#include <opencv2/core/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

using namespace cv;
using namespace std;

bool image_equals(const Mat& image1, const Mat& image2) {
    if (image1.cols == image2.cols && image1.rows == image2.rows) {
        for (int i = 0; i < image1.rows; ++i) {
            for (int j = 0; j < image1.cols; ++j) {
                if (image1.at<uchar>(i,j) != image1.at<uchar>(i,j)) {
                    return 0;
                }
            }
        }
    } else {
        return 0;
    }
    return 1;
}

Mat binarize(Mat& image, const int board) {
    for (int i = 0; i < image.rows; ++i) {
        for (int j = 0; j < image.cols; ++j) {
            if (image.at<uchar>(i, j) <= board) {
                image.at<uchar>(i, j) = 0;
            } else {
                image.at<uchar>(i, j) = 255;
            }
        }
    }
    imwrite("new_image.jpg", image);
    return imread("new_image.jpg");
}

Mat otsu(Mat& image) {
    std::ofstream outfile("new_image.jpg");
    short board = 0;
    double first_class_mean = 0;
    double second_class_mean = 0;
    double first_class_prob;
    double second_class_prob;
    int max_scale = 0;
    for (int i = 0; i < image.rows; ++i) {
        for (int j = 0; j < image.cols; ++j) {
            if (image.at<uchar>(i,j) > max_scale) {
                max_scale = image.at<uchar>(i,j);
            }
        }
    }
    vector<long long> histogram(max_scale);
    for (int i = 0; i < image.rows; ++i) {
        for (int j = 0; j < image.cols; ++j) {
            ++histogram[image.at<uchar>(i, j)];
        }
    }
    unsigned hist_sum = 0;
    for (int i = 0; i < histogram.size(); ++i) {
        hist_sum += histogram[i];
    }
    first_class_prob = (double) histogram[0] / hist_sum;
    second_class_prob = (double) 1 - first_class_prob;
    first_class_mean = 0;
    second_class_mean = 1 - first_class_mean;
    double dispersion = (double) first_class_prob * second_class_prob * (first_class_mean - second_class_mean) * (first_class_mean - second_class_mean);
    double max_dispersion = dispersion;
    
    for (int i = 1; i < histogram.size(); ++i) {
        first_class_prob = 0;
        first_class_mean = 0;
        second_class_prob = 0;
        second_class_mean = 0;
        for (int j = 0; j <= i; ++j) {
            first_class_prob += histogram[j];
        }
        first_class_prob = (double) first_class_prob / hist_sum;
        second_class_prob = 1 - first_class_prob;
        for (int j = 0; j <= i; ++j) {
            first_class_mean = first_class_mean + j * histogram[j];
        }
        first_class_mean = (double) first_class_mean / first_class_prob;
        second_class_mean = 1 - first_class_mean;
        dispersion = (double) first_class_prob * second_class_prob * (first_class_mean - second_class_mean) * (first_class_mean - second_class_mean);
        if (dispersion > max_dispersion) {
            max_dispersion = dispersion;
            board = i;
        }
    }
    Mat result_image = binarize(image, board);
    outfile.close();
    return result_image;
}

int main(int argc, char** argv)
{
    string image_name = "image.jpg";
    Mat image = imread(image_name, 0);
    Mat empty_image;
    if (!image.rows && !image.cols) {
        cerr << "Error: Image must be called \"image.jpg\"" << endl;
        return 1;
    }
    std::ofstream outfile("opencv_otsu_result.jpg");
    Mat opencv_otsu_result = imread("opencv_otsu_result.jpg");
    threshold(image, opencv_otsu_result, 0, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
    imwrite("opencv_otsu_result.jpg", opencv_otsu_result);
    outfile.close();
    if (!image_equals(otsu(image), opencv_otsu_result)) {
        cerr << ":-(" << endl;
    }
    return 0;
}
