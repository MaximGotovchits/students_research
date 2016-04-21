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

//class pixels_class_colored {
//public:
//    int first_board_r;
//    int second_board_r;
//    int first_board_g;
//    int second_board_g;
//    int first_board_b;
//    int second_board_b;
//
//    double dispersion;
//    
//    pixels_class_colored() {}
//    
//    pixels_class_colored(int _first_board_r, int _second_board_r,
//                         int _first_board_g, int _second_board_g,
//                         int _first_board_b, int _second_board_b,
//                         double _dispersion = 0) {
//        first_board_r = _first_board_r;
//        second_board_r = _second_board_r;
//        first_board_g = _first_board_g;
//        second_board_g = _second_board_g;
//        first_board_b = _first_board_b;
//        second_board_b = _second_board_b;
//
//        dispersion = _dispersion;
//    }
//    
//    bool operator < (const pixels_class_colored& pc) const {
//        return dispersion < pc.dispersion;
//    }
//
//};

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
    //result_image.setTo(front_color, image <= board & image >= front_color);
    //result_image.setTo(back_color, image > board & image <= back_color);
}

double& zero() {
    double zero = 0;
    return zero;
}

int get_board(Mat& image, int first_board, int second_board, double& dispersion = zero()) {
//    cout << dispersion;
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
    //cout << "COLORS: " << colors << endl;
    unsigned hist_sum = accumulate(histogram.begin(), histogram.end(), 0);
    double first_class_mean = 0;
    double second_class_mean = 0;
    
    first_class_prob = 0;
    second_class_prob = 0;
    
    dispersion = 0;
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

Mat get_channel(vector<Mat> channel) {
    vector<int> dispersions;
    dispersions.push_back(get_dispersion(channel[0], 0, 255));
    dispersions.push_back(get_dispersion(channel[1], 0, 255));
    dispersions.push_back(get_dispersion(channel[2], 0, 255));
    return channel[distance(dispersions.begin(), max_element(dispersions.begin(), dispersions.end()))];
}

class channel_board {
public:
    channel_board() {
        first_board = 0;
        second_board = 256;
    }
    
    int first_board;
    int second_board;
};

class volume_class {
public:
    volume_class() {
        for (int i = 0; i < 3; ++i) {
            channel_board default_channel;
            channel.push_back(default_channel);
            color.push_back(0);
        }
    }
    
    vector<channel_board> channel;
    vector<int> color;
};

double vol_get_dispersion(deque<volume_class>::iterator box, vector<Mat>& channel) {
    double dispersion = 0.0;
    for (int i = 0; i < 3; ++i) {
        dispersion += get_dispersion(channel[i], box -> channel[i].first_board, box -> channel[i].second_board);
    }
    return dispersion;
}

void add_box(deque<volume_class>& class_boxes, vector<Mat>& channel) {
    double dispersion = 0.0;
    double max_dispersion = 0.0;
    int channel_to_split = 0;
    int board = 0;
    int current_board = 0;
    auto box_to_split = class_boxes.begin();
    
    
    for (auto box = class_boxes.begin(); box != class_boxes.end(); ++box) {
        for (int j = 0; j < 3; ++j) {
            current_board = get_board(channel[j], box -> channel[j].first_board, box -> channel[j].second_board, dispersion);
            if (dispersion > max_dispersion) {
                max_dispersion = dispersion;
                channel_to_split = j;
                box_to_split = box;
                board = current_board;
            }
        }
    }
    vector<volume_class> pretenders;
    max_dispersion = vol_get_dispersion(box_to_split, channel);
    for (auto box = class_boxes.begin(); box != class_boxes.end(); ++box) {
        if((box -> channel[channel_to_split].first_board == box_to_split -> channel[channel_to_split].first_board) &&
           (box -> channel[channel_to_split].second_board == box_to_split -> channel[channel_to_split].second_board)) {
            dispersion = vol_get_dispersion(box, channel);
            if (dispersion > max_dispersion) {
                max_dispersion = dispersion;
                box_to_split = box;
            }
            
        }
    }
    
    volume_class first_box;
    volume_class second_box;
    for (int i = 0; i < 3; ++i) {
        if (i == channel_to_split) {
            first_box.channel[i].first_board = box_to_split -> channel[i].first_board;
            first_box.channel[i].second_board = board;
            second_box.channel[i].first_board = board;
            second_box.channel[i].second_board = box_to_split -> channel[i].second_board;
        } else {
            first_box.channel[i].first_board = box_to_split -> channel[i].first_board;
            first_box.channel[i].second_board = box_to_split -> channel[i].second_board;
            second_box.channel[i].first_board = box_to_split -> channel[i].first_board;
            second_box.channel[i].first_board = box_to_split -> channel[i].second_board;
        }
    }
    class_boxes.erase(box_to_split);
    class_boxes.push_back(first_box);
    class_boxes.push_back(second_box);
}

int get_channel_color(Mat& image, int first_board, int second_board) {
    double mean = 0.0;
    int pixel_amount = 0;
    for (int i = 0; i < image.rows; ++i) {
        for (int j = 0; j < image.cols; ++j) {
            if ((int)image.at<uchar>(i, j) >= first_board && (int)image.at<uchar>(i, j) <= second_board) {
                mean += (int)image.at<uchar>(i, j);
                ++pixel_amount;
            }
        }
    }
    mean /= pixel_amount;
    int color = mean;
    return color;
}

void get_box_color(deque<volume_class>::iterator box, vector<Mat>& channel) {
    for (int i = 0; i < 3; ++i) {
        box -> color[i] = get_channel_color(channel[i], box -> channel[i].first_board, box -> channel[i].second_board);
    }
}

void compute_colors(deque<volume_class>& class_boxes, vector<Mat>& channel) {
    for (auto box = class_boxes.begin(); box != class_boxes.end(); ++box) {
        get_box_color(box, channel);
    }
}

void box_colorize(deque<volume_class>::iterator box, vector<Mat>& channel) {
    for (int i = 0; i < 3; ++i) { // Choose channel.
        for (int j = 0; j < channel[i].rows; ++j) { // For each pixel in the channel.
            for (int k = 0; k < channel[i].cols; ++k) { // For each pixel in the channel.
                if ((int)channel[i].at<uchar>(j, k) >= box -> channel[i].first_board &&
                    (int)channel[i].at<uchar>(j, k) <= box -> channel[i].second_board) {
                    channel[i].at<uchar>(j, k) = (uchar)box -> color[i];
                }
            }
        }
    }
}

void vol_colorize(deque<volume_class>& class_boxes, vector<Mat>& channel) {
    for (auto box = class_boxes.begin(); box != class_boxes.end(); ++box) {
        box_colorize(box, channel);
    }
}

Mat multicolor_multhreshold_otsu(Mat image, int color_amount) {
    // Incorrect cases:
    if (color_amount < 1 || color_amount > 16777216) {
        cerr << "Incorrect number of colors" << endl;
    }
    // Degenarate cases:
    if (color_amount == 16777216) {
        return image;
    }
    if (color_amount == 1) {
        vector<Mat> channel(3);
        split(image, channel);
        for (int k = 0; k < 3; ++k) {
            double mean = 0;
            for (int i = 0; i < image.rows; ++i) {
                for (int j = 0; j < image.cols; ++j) {
                    mean += (int)image.at<uchar>(i, j);
                }
            }
            mean /= (channel[k].rows * channel[k].cols);
            int color;
            if (mean < 256 / 2) {
                color = 0;
            } else {
                color = 255;
            }
            for (int i = 0; i < channel[k].rows; ++i) {
                for (int j = 0; j < channel[k].cols; ++j) {
                    channel[k].at<uchar>(i, j) = (uchar)color;
                }
            }
        }
        Mat result;
        merge(channel, result);
        return result;
    }
    // Normal cases:
    if (color_amount == 2) {
        return binary_otsu(image);
    } else {
                /*
         1. Get channel (by max_dispersion ejection).
         3. Get area with max regural dispersion ().
         2. Get board for the chosen channel.
         3.
        */
        vector<Mat> channel(3);
        split(image, channel);
        int channel_to_split = 0;
        double max_dispersion = 0.0;
        double dispersion = 0.0;
        int board = 0;
        for (int i = 0; i < 3; ++i) {
            int current_board = get_board(channel[i], 0, 256);
            if (dispersion > max_dispersion) {
                max_dispersion = dispersion;
                channel_to_split = i;
                board = current_board;
            }
        }
        // init.
        volume_class first_box;
        first_box.channel[channel_to_split].first_board = 0;
        first_box.channel[channel_to_split].second_board = board;
        volume_class second_box;
        second_box.channel[channel_to_split].first_board = board;
        second_box.channel[channel_to_split].second_board = 255;
        deque<volume_class> class_boxes;
        class_boxes.push_back(first_box);
        // end init.
        for (int i = 1; i < color_amount; ++i) {
            add_box(class_boxes, channel);
        }
        compute_colors(class_boxes, channel);
        vol_colorize(class_boxes, channel);
        Mat result;
        merge(channel, result);
        return result;
    }
    
//    vector<Mat> channels(3);
//    split(image, channels);
//    for (int k = 0; k < 3; ++k) {
//        channels[k] = multithreshold_otsu(channels[k], color_amount[k]);
//    }
//    Mat result;
//    merge(channels, result);
//    return result;

    return Mat();
}

//void init_pixel_classes(Mat& image, deque<pixels_class>& pixel_classes) {
//    pixel_classes[0].first_board = 0;
//    pixel_classes[0].second_board = get_board(image, 0, 255);
//    pixel_classes[0].dispersion = get_dispersion(image, pixel_classes[0].first_board, pixel_classes[0].second_board);
//    pixel_classes[1].first_board = pixel_classes[0].second_board;
//    pixel_classes[1].second_board = 255;
//    pixel_classes[1].dispersion = get_dispersion(image, pixel_classes[1].first_board, pixel_classes[1].second_board);
//}

int main(int argc, char** argv) {
    string image_name = "image.png";
    Mat color_image = imread(image_name);
    Mat image = imread(image_name, 0);
//    for (int i = 2; i <= 64; i *= 2) {
        Mat color_result = multicolor_multhreshold_otsu(color_image, 256);
        string name = to_string(256) + ".png";
        imwrite(name, color_result);
  //  }
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
