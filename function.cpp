#include "function.h"
extern Mat K;
extern string inputFileDir;
extern string outputFileDir;
extern string dataset;
extern double eDisR, eDist;
extern float distur;
extern int fittingNum;
int countt=0;
float s=0;


void Mat2Eigen(Mat M, Eigen::MatrixXd &E )
{
    int r=M.rows;
    int c=M.cols;
    if(M.empty())
    {
        cout<<"the Mat is emputy"<<endl;
    }
    for(int i=0;i<r;i++)
    {
        for(int j=0;j<c;j++)
        {
            E(i,j)=M.at<double>(i,j);
        }
    }
}

void calculateError(vector<Point3f> world3d, vector<Point2f> camera2d, Mat R, Mat t, string methodName_)
{
    Eigen::MatrixXd trueR(3,3);
    Eigen::MatrixXd truet(3,1);

    float error=0,error_sum=0,error_ave=0,error_x=0,error_y=0,error_z=0, error_x_sum=0, error_y_sum=0,error_z_sum=0, error_x_ave=0, error_y_ave=0,error_z_ave=0;

    cv::Mat tempMat, tempMat2;
    Mat wcPoint;
    cv::Mat imagePoint = cv::Mat::ones(3, 1, cv::DataType<double>::type); //u,v,1

    trueR<<0.9961980177959597, 0.0196, -0.0849,
            0.0872, -0.2242505587701737, -0.9707,
            0, 0.9706297846156308, -0.2251093986976954;

    truet<<0,
            5261.946922898293,
            1208.48645468081;

    Eigen::MatrixXd RE(3,3);
    Eigen::MatrixXd tE(3,1);
    Eigen::MatrixXd KE(3,3);

    Mat2Eigen(R,RE);
    Mat2Eigen(t,tE);
    Mat2Eigen(K,KE);

    calculateDis(RE,tE);

    Eigen::MatrixXd p(3,1);
    Eigen::MatrixXd P(3,1);

    ofstream errorf(outputFileDir+methodName_+"_"+"error.txt");//+methodName_+"/"
    ofstream xErrorf(outputFileDir+methodName_+"-"+"xError.txt");
    ofstream yErrorf(outputFileDir+methodName_+"-"+"yError.txt");
    ofstream distancef(outputFileDir+methodName_+"_"+"distance.txt");

    errorf<<"the method is: "<<methodName_<<"\n"<<"the dataset is:"<<dataset<<"\n"<<"\n";
    errorf<<"the estimated R is: "<<R<<"\n"<<"the estimated t is: "<<t<<"\n";
    errorf<<"the euclidean distance for R is:"<<eDisR<<"\n"<<"the euclidean distance for t is:"<<eDist<<"\n"<<"\n";

    vector<float> verror_x;
    vector<float> verror_y;
    vector<float> verror_z;

    // compute the error
    for(int i=0; i<world3d.size(); i++)
    {
        countt=countt+1;

        imagePoint = cv::Mat::ones(3, 1, cv::DataType<double>::type); //u,v,1
        imagePoint.at<double>(0, 0) = camera2d[i].x;
        imagePoint.at<double>(1, 0) = camera2d[i].y;

        //compute the scale s
        tempMat=R.inv()*K.inv()*imagePoint;
        tempMat2=R.inv()*t;
        s=tempMat2.at<double>(2,0);
        s/=tempMat.at<double>(2,0);

        //compute the concidante in the world
        wcPoint = R.inv() * (s * K.inv() * imagePoint - t);
        Point3f worldPointLocal(wcPoint.at<double>(0, 0), wcPoint.at<double>(1, 0), wcPoint.at<double>(2, 0));

        error_x=abs(worldPointLocal.x-world3d[i].x);
        error_y=abs(worldPointLocal.y-world3d[i].y);
        error_z=abs(worldPointLocal.z-world3d[i].z);

        distancef<<world3d[i].y<<"\n";
        xErrorf<<error_x<<"\n";
        yErrorf<<error_y<<"\n";
        countt = 0;
    }

    errorf.close();
    xErrorf.close();
    yErrorf.close();
    distancef.close();
}

//calculate the euclidean distance
void calculateDis(Eigen::MatrixXd R, Eigen::MatrixXd t)
{
    double sum=0;
    Eigen::MatrixXd trueR(3,3);
    Eigen::MatrixXd truet(3,1);

    trueR<<0.9961980177959597, 0.08711778296237149, 3.508734154855642e-05,
            0.01964508855730578, -0.2242505587701737, -0.9743334939264073,
            -0.08487390550090182, 0.9706297846156308, -0.2251093986976954;

    truet<<-1.495050585197469,
            5261.946922898293,
            1208.48645468081;

    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            sum=sum+(R(i,j)-trueR(i,j))*(R(i,j)-trueR(i,j));
        }
    }
    eDisR=pow(sum,0.5);
    sum=0;
    for(int i=0;i<3;i++)
    {
        sum=sum+pow((t(i,0)-truet(i,0)),2);
    }
    eDist=pow(sum,0.5);
}

void calculateWorld(vector<Point2f> camera2d,vector<Point3f> world3d, Mat R, Mat t, string methodName_, string dataset_)
{
    ofstream worldPointf(outputFileDir+methodName_+"_"+"worldPoint.txt");

    Eigen::MatrixXd RE(3,3);
    Eigen::MatrixXd tE(3,1);
    Eigen::MatrixXd KE(3,3);

    Mat2Eigen(R,RE);
    Mat2Eigen(t,tE);
    Mat2Eigen(K,KE);
    calculateDis(RE,tE);

    worldPointf<<"the method is: "<<methodName_<<"\n"<<"the measurement equipment is gps"<<"\n";
//    worldPointf<<"the method is: "<<methodName_<<"\n"<<"the inputed camera plane dataset is:"<<dataset_<<"\n"<<"the distur is "<<distur<<"\n"<<"the number of fitting points is: "<<fittingNum+1<<"\n"<<"\n";
//    worldPointf<<"the estimated R is: "<<R<<"\n"<<"the estimated t is: "<<t<<"\n";
//    worldPointf<<"the euclidean distance for R is:"<<eDisR<<"\n"<<"the euclidean distance for t is:"<<eDist<<"\n"<<"\n";

    cv::Mat tempMat, tempMat2;
    Mat wcPoint;
    cv::Mat imagePoint = cv::Mat::ones(3, 1, cv::DataType<double>::type); //u,v,1

    for(int i=0; i<camera2d.size(); i++)
    {
        imagePoint = cv::Mat::ones(3, 1, cv::DataType<double>::type); //u,v,1
        imagePoint.at<double>(0, 0) = camera2d[i].x;
        imagePoint.at<double>(1, 0) = camera2d[i].y;

        //compute the scale s
        tempMat=R.inv()*K.inv()*imagePoint;
        tempMat2=R.inv()*t;
        s=tempMat2.at<double>(2,0);
        s/=tempMat.at<double>(2,0);

        //compute the concidante in the world
        wcPoint = R.inv() * (s * K.inv() * imagePoint - t);
        Point3f worldPointLocal(wcPoint.at<double>(0, 0), wcPoint.at<double>(1, 0), wcPoint.at<double>(2, 0));
        //worldPointf<<"the point in camrea plane is: "<<imagePoint.t()<<"   the point in world plane is: "<<worldPointLocal<<"\n";
        worldPointf<<worldPointLocal.x<<"   "<<worldPointLocal.y<<"      "<<world3d[i].x<<"   "<<world3d[i].y<<"      "<<abs(worldPointLocal.x-world3d[i].x)<<"   "<<abs(worldPointLocal.y-world3d[i].y)<<"\n";
    }
    worldPointf.close();
}


