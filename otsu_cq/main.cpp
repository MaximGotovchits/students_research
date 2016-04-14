#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <deque>

using namespace cv;
using namespace std;

class pixels_class {
public:
    int first_board;
    int second_board;
    double dispersion;
    
    pixels_class() {}
    
    pixels_class(int _first_board, int _second_board, double _dispersion = 0) {
        first_board = _first_board;
        second_board = _second_board;
        dispersion = _dispersion;
    }
    
    bool operator < (const pixels_class& pc) const {
        return dispersion < pc.dispersion;
    }
};

bool image_equals(const Mat& image1, const Mat& image2) {
    return ((countNonZero(image1 - image2) == 0) && (countNonZero(image2 - image1) == 0));
}

void binarize(Mat image, Mat& result_image, int board, int front_color, int back_color) {
    if (front_color <= board && back_color >= board) {
        //cout << "OK" << endl;
    }
    for (int i = 0; i < image.rows; ++i) {
        for (int j = 0; j < image.cols; ++j) {
            if ((int)image.at<uchar>(i, j) <= board && (int)image.at<uchar>(i, j) >= front_color) {
                result_image.at<uchar>(i, j) = (uchar)front_color;
            } else if ((int)image.at<uchar>(i, j) > board && (int)image.at<uchar>(i, j) <= back_color) {
                result_image.at<uchar>(i, j) = (uchar)back_color;
            }
        }
    }
    //return result_image;
}

int get_board(Mat& image, int first_board, int second_board) {
    if (first_board == second_board) {
        return first_board;
    }
    double first_class_prob;
    double second_class_prob;
    vector<int> histogram(256);
    for (int j = 0; j < image.cols; ++j) {
        for (int i = 0; i < image.rows; ++i) {
            if ((int)image.at<uchar>(i, j) >= first_board && (int)image.at<uchar>(i, j) <= second_board) {
                ++histogram[(int)image.at<uchar>(i, j)];
            }
        }
    }
    int colors = 0;
    for (auto it = histogram.begin(); it != histogram.end(); ++it) {
        if (*it != 0) {
            ++colors;
        }
    }
    // cout << "COLORS: " << colors << endl;
    unsigned hist_sum = accumulate(histogram.begin(), histogram.end(), 0);
    double first_class_mean = 0;
    double second_class_mean = 0;
    first_class_prob = 0;
    second_class_prob = 0;
    double dispersion = 0;
    double max_dispersion = dispersion;
    double first_class_sum = 0;
    double first_class_not_mean = 0;
    double sum = 0.0;
    int j = 0;
    for (auto bin_value = histogram.begin(); bin_value != histogram.end(); ++bin_value) {
        sum += (j * (*bin_value));
        ++j;
    }
    int board = 0;
    for (int j = first_board; j <= second_board; ++j) {
        first_class_sum += histogram[j];
        if (first_class_sum == 0) {
            continue;
        }
        double second_class_sum = hist_sum - first_class_sum;
        if (second_class_sum == 0) {
            break;
        }
        first_class_not_mean += (float) (j * histogram[j]);
        first_class_mean = first_class_not_mean / first_class_sum;
        second_class_mean = (sum - first_class_not_mean) / second_class_sum;
        dispersion = (float) first_class_sum * second_class_sum * (first_class_mean - second_class_mean) * (first_class_mean - second_class_mean);
        if (dispersion > max_dispersion) {
            max_dispersion = dispersion;
            board = j;
        }
    }
    return board;
}

Mat binary_otsu(Mat& image) {
    int board = get_board(image, 0, 255);
    Mat result_image = image;
    binarize(image, result_image, board, 0, 255); // front_color(0) < back_color(255)
    return result_image;
}

Mat colorize(Mat& image, vector<int>& boards) {
    Mat result_image = image;
    imwrite("MUST_BE_GREY.png", image);
    sort(boards.begin(), boards.end());
    for (auto board = boards.begin(); board != boards.end(); ++board) {
        if (board == boards.begin()) {
            binarize(image, result_image, *board, 0, *(board + 1));
            continue;
        }
        if (board == boards.end() - 1) {
            binarize(image, result_image, *board, *(board - 1), 255);
            continue;
        }
        binarize(image, result_image, *board, *(board - 1), *(board + 1));
    }
    return result_image;
}

double get_dispersion(Mat& image, int first_board, int second_board) {
    double mean = 0.0;
    vector<int> class_pixels;
    for (int i = 0; i < image.rows; ++i) {
        for (int j = 0; j < image.cols; ++j) {
            if ((int)image.at<uchar>(i, j) >= first_board && (int)image.at<uchar>(i, j) < second_board) {
                mean += (int)image.at<uchar>(i, j);
                class_pixels.push_back((int)image.at<uchar>(i, j));
            }
        }
    }
    mean /= (double)class_pixels.size();
    double dispersion = 0.0;
    for (auto pixel = class_pixels.begin(); pixel != class_pixels.end(); ++pixel) {
        dispersion += (double) abs(mean - *pixel);
    }
    dispersion /= class_pixels.size();
    return dispersion;
}

void init_pixel_classes(Mat& image, deque<pixels_class>& pixel_classes) {
    pixel_classes[0].first_board = 0;
    pixel_classes[0].second_board = get_board(image, 0, 255);
    pixel_classes[0].dispersion = get_dispersion(image, pixel_classes[0].first_board, pixel_classes[0].second_board);
    pixel_classes[1].first_board = pixel_classes[0].second_board;
    pixel_classes[1].second_board = 255;
    pixel_classes[1].dispersion = get_dispersion(image, pixel_classes[1].first_board, pixel_classes[1].second_board);
}

Mat multithreshold_otsu(Mat& image, int color_amount) {
    // Incorrect cases:
    if (color_amount < 1 || color_amount > 256) {
        cerr << "Incorrect number of colors" << endl;
    }
    // Degenarate cases:
    if (color_amount == 256) {
        return image;
    }
    Mat result_image = image;
    if (color_amount == 1) {
        double mean = 0;
        for (int i = 0; i < image.rows; ++i) {
            for (int j = 0; j < image.cols; ++j) {
                mean += (int)image.at<uchar>(i, j);
            }
        }
        mean /= (image.rows * image.cols);
        int color;
        if (mean < 256 / 2) {
            color = 0;
        } else {
            color = 255;
        }
        for (int i = 0; i < image.rows; ++i) {
            for (int j = 0; j < image.cols; ++j) {
                result_image.at<uchar>(i, j) = (uchar)color;
            }
        }
        return result_image;
    }
    // Normal cases:
    if (color_amount == 2) {
        return binary_otsu(image);
    } else {
        deque<pixels_class> pixel_classes(2);
        init_pixel_classes(image, pixel_classes);
        int counter = 2;
        vector<int> boards;
        boards.push_back(pixel_classes[0].second_board);
        while(counter < color_amount) {
            auto max_subclass = max_element(pixel_classes.begin(), pixel_classes.end());
            int board = get_board(image, max_subclass -> first_board, max_subclass -> second_board);
            boards.push_back(board);
            pixel_classes.push_back(pixels_class(max_subclass -> first_board, board,
                                                 get_dispersion(image, max_subclass -> first_board, board)));
            pixel_classes.push_back(pixels_class(board, max_subclass -> second_board,
                                                 get_dispersion(image, board, max_subclass -> second_board)));
            pixel_classes.erase(max_subclass);
            ++counter;
        }
        return colorize(image, boards);
    }
    return result_image;
}

Mat multi_color_multhreshold_otsu(Mat image, vector<int> color_amount) {
    Mat result = image;
    Mat temp;
    cvtColor(image, temp, CV_BGR2GRAY);
    for (int k = 0; k < 3; ++k) {
        for (int j = 0; j < image.cols; ++j) {
            for (int i = 0; i < image.rows; ++i) {
                temp.at<uchar>(i, j) = image.at<Vec3b>(i, j)[k];
            }
        }
        Mat one_channel_result = multithreshold_otsu(temp, color_amount[k]);
        for (int j = 0; j < one_channel_result.cols; ++j) {
            for (int i = 0; i < one_channel_result.rows; ++i) {
                result.at<Vec3b>(i, j)[k] = one_channel_result.at<uchar>(i, j);
            }
        }
    }
    return result;
}

int main(int argc, char** argv) {
    string image_name = "image.png";
    Mat color_image = imread(image_name);
    Mat image = imread(image_name, 0);
    vector<int> color_amount = {2, 2, 2};
    Mat color_result = multi_color_multhreshold_otsu(color_image, color_amount);
    imwrite("colored_result.png", color_result);
    if (!image.rows && !image.cols) {
        cerr << "Error: Image must be called \"image.png\"" << endl;
        return 1;
    }
//    Mat opencv_otsu_result = imread("opencv_otsu_result.png");
//    threshold(image, opencv_otsu_result, 0, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
//    imwrite("opencv_otsu_result.png", opencv_otsu_result);
    Mat grey_result = multithreshold_otsu(image, 2);
    imwrite("grey_result.png", grey_result);
//    if (!image_equals(grey_result, opencv_otsu_result)) {
//        cerr << ":-(" << endl;
//    }
    return 0;
}
