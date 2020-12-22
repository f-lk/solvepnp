#ifndef VO1_FUNCTION_H
#define VO1_FUNCTION_H

#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <Core>
#include <Geometry>
#include <chrono>
#include <fstream>
using namespace std;
using namespace cv;

void calculateDis(Eigen::MatrixXd R, Eigen::MatrixXd t);
void Mat2Eigen(Mat M, Eigen::MatrixXd &E );
void calculateError(vector<Point3f> world3d, vector<Point2f> camera2d, Mat R, Mat t, string methodNmae_);
void calculateWorld(vector<Point2f> camera2d, vector<Point3f> world3d, Mat R, Mat t, string methodName_, string dataset_);


#endif //VO1_FUNCTION_H
