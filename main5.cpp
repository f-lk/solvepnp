#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <Core>
#include <Geometry>
#include <chrono>
#include <fstream>
#include "function.h"

using namespace std;
using namespace cv;

//enum SolvePnPMethod {
//    SOLVEPNP_ITERATIVE   = 0,
//    SOLVEPNP_EPNP        = 1, //!< EPnP: Efficient Perspective-n-Point Camera Pose Estimation @cite lepetit2009epnp
//    SOLVEPNP_P3P         = 2, //!< Complete Solution Classification for the Perspective-Three-Point Problem @cite gao2003complete
//    SOLVEPNP_DLS         = 3, //!< A Direct Least-Squares (DLS) Method for PnP  @cite hesch2011direct
//    SOLVEPNP_UPNP        = 4, //!< Exhaustive Linearization for Robust Camera Pose and Focal Length Estimation @cite penate2013exhaustive
//    SOLVEPNP_AP3P        = 5, //!< An Efficient Algebraic Solution to the Perspective-Three-Point Problem @cite Ke17
//    SOLVEPNP_IPPE        = 6, //!< Infinitesimal Plane-Based Pose Estimation @cite Collins14 \n
//    //!< Object points must be coplanar.
//    SOLVEPNP_IPPE_SQUARE = 7,
double eDisR, eDist;
int method=0; //   2/5/7
string methodName;
string txt=".txt";
string dataset="truth";
string inputFileDir="/home/flk/SlovePnP/pnpData/experiment_12.16/rowdata/distur2-5_rawdata/"+dataset+"/";
string testInputFileDir="/home/flk/SlovePnP/pnpData/experiment_12.17/";
string outputFileDir="/home/flk/SlovePnP/pnpData/experiment_12.17/";
float distur = 0.01;
int fittingNum =24;
vector<Point3f> worldP3Part;
vector<Point2f> cameraP2Part;
Mat K = ( Mat_<double> ( 3,3 ) << 4838.4, 0, 1280, 0, 3628.8, 720, 0, 0, 1 );
int maxy=98000;

int main ( int argc, char** argv )
{
    for(method;method<7;method++)
    {
        if(method==2||method==5||method==7||method==6)
            continue;
        else
        {
            switch (method) {
                case 0:
                    methodName="ITERATIVE";
                    break;
                case 1:
                    methodName="EPNP";
                    break;
                case 2:
                    methodName="P3P";
                    break;
                case 3:
                    methodName="DLS";
                    break;
                case 4:
                    methodName="UPNP";
                    break;
                case 5:
                    methodName="AP3P";
                    break;
                case 6:
                    methodName="IPPE";
                    break;
            }
            cout<<"the method is: "<<methodName<<endl;

            vector<Point3f> pts_3d;
            vector<Point2f> pts_2d;
            vector<Point3f> pts_3d_;
            vector<Point2f> pts_2d_;
            
            Point2f midleP;
            Point3f midleP3;
            float num1, num2, num3;

            ifstream fu(inputFileDir+"u_origin.txt");
            ifstream fv(inputFileDir+"v_origin.txt");
            ifstream fx(inputFileDir+"x_origin.txt");
            ifstream fy(inputFileDir+"y_origin.txt");
            ifstream fz(inputFileDir+"z_origin.txt");

            ifstream fTe_u(testInputFileDir+"test_u.txt");
            ifstream fTe_v(testInputFileDir+"test_v.txt");

            assert(fu.is_open()); //若失败,则输出错误消息,并终止程序运行

            vector<Point2f> testCamPoint;

            while (!fTe_u.eof()&&!fTe_v.eof())
            {
                fTe_u>>num1;
                fTe_v>>num2;
                midleP.x=num1;
                midleP.y=num2;
                testCamPoint.push_back(midleP);
            }

            while (!fu.eof()&&!fv.eof())
            {
                fu>>num1;
                fv>>num2;
                midleP.x=num1;
                midleP.y=num2;
                pts_2d.push_back(midleP);
            }
            cout<<"size of pts_2d"<<pts_2d.size()<<endl;
            while (!fx.eof()&&!fy.eof()&&!fz.eof())
            {
                fx>>num1;
                fy>>num2;
                fz>>num3;
                midleP3.x=num1;
                midleP3.y=num2;
                midleP3.z=num3;
                pts_3d.push_back(midleP3);
            }

            for(int i=0; i<pts_3d.size(); i++)
            {
                if(pts_3d[i].x==500||pts_3d[i].x==-3500||pts_3d[i].x==2000)
                {
                    pts_3d_.push_back(pts_3d[i]);
                    pts_2d_.push_back(pts_2d[i]);
                }
            }



            // add random distur
            srand( (unsigned)time( NULL ) );
            for(int i=0;i<pts_3d.size();i++)
            {
                pts_3d[i].x=pts_3d[i].x*(1+distur*((rand()%100/(double)101)-0.5));
                pts_3d[i].y=pts_3d[i].y*(1+distur*((rand()%100/(double)101)-0.5));
                pts_3d[i].z=pts_3d[i].z*(1+distur*((rand()%100/(double)101)-0.5));
                pts_2d[i].x=pts_2d[i].x*(1+distur*((rand()%100/(double)101)-0.5));
                pts_2d[i].y=pts_2d[i].y*(1+distur*((rand()%100/(double)101)-0.5));
            }

            Mat r, t;

            worldP3Part.clear();
            cameraP2Part.clear();

            for(int i=0;i<=fittingNum;i++)
            {
                worldP3Part.push_back(pts_3d[i]);
                cameraP2Part.push_back(pts_2d[i]);
            }

            // solve the pnp
            cout<<"method: "<<method<<endl;
            cout<<"the number of points used to solvepnp is : "<<worldP3Part.size()<<endl;
            solvePnP( worldP3Part, cameraP2Part, K, Mat(), r, t, false ,method); // 调用OpenCV 的 PnP 求解，可选择EPNP，DLS等方法
            Mat R;


            cv::Rodrigues ( r, R ); // r为旋转向量形式，用Rodrigues公式转换为矩阵

            calculateWorld(testCamPoint,R, t, methodName,"test1");

//            calculateError(pts_3d,pts_2d,R,t,methodName);

        }
    }
}










