#pragma once
#include <stdio.h>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <queue>
#include "math/czy_math_fit.hpp"

using namespace czy;

class LateralControl{
     public:
          //参数赋值
          double m,Cf,Cr,a,b,Iz,vx,vy,l;
          LateralControl(double speed_ms):vx(0.0){
               m=1500.0;
               Cf=-46085.0;
               Cr=-58915.0;
               a=1.016;
               b=1.562;
               Iz=2000.0;
               vx=speed_ms;
               vy=0.0;
               l=2.8;//轴距
          }

          //path路径信息
          double theta_dmin;
          double x_dmin;
          double y_dmin;
          double k_dmin;
          void path_info_extract(std::vector<double> &PosX, std::vector<double> &PosY, std::vector<double> &theta_r){
               x_dmin = vx*0.1;
               int index=0;
               for(int i=0;i<36;i++){
                    if(PosX[i]>=vx*1.0){
                         index=i;
                         break;
                    }
               }
               theta_dmin=(theta_r[index-1]+theta_r[index])/2;

               Fit F;
               F.polyfit(PosX,PosY,5,false);
               std::vector<double> f;
               F.getFactor(f);
               y_dmin= f[0]+f[1]*x_dmin+f[2]*pow(x_dmin,2)+f[3]*pow(x_dmin,3)+f[4]*pow(x_dmin,4)+f[5]*pow(x_dmin,5);
               double dyDdx = f[1]+2.0*f[2]*x_dmin+3.0*f[3]*pow(x_dmin,2)+4.0*f[4]*pow(x_dmin,3)+5.0*f[5]*pow(x_dmin,4);//一阶导数
               double dyDdx2 = 2.0*f[2]+6.0*f[3]*x_dmin+12.0*f[4]*pow(x_dmin,2)+20.0*f[5]*pow(x_dmin,3);//二阶导数
               k_dmin = dyDdx2/pow(1.0+dyDdx*dyDdx,1.5);
               if(isnan(k_dmin)) k_dmin = 0.0;
          }

          //计算err_k
          std::array<double,5> cal_err_k(){
               Eigen::Matrix<double,2,1> nor;
               nor << -sin(theta_dmin),cos(theta_dmin);
               Eigen::Matrix<double,2,1> tau;
               tau << cos(theta_dmin),sin(theta_dmin);
               Eigen::Matrix<double,2,1> d_err;
               double tp = 0.1;//预测时间
               double phi_dot = vx*k_dmin;
               double phi = 0.0;//当前横摆角
               double phi_pre = phi+phi_dot*tp;//预测横摆角
               double pre_pos_x=0.0+vx*tp*cos(phi)-vy*tp*sin(phi);
               double pre_pos_y=0.0+vy*tp*cos(phi)+vx*tp*sin(phi);//车身坐标系下，以车自身作为坐标原点，向前预测ts时间的距离
               d_err << pre_pos_x-x_dmin,pre_pos_y-y_dmin;
               double ed = nor.transpose()*d_err;
               double es = tau.transpose()*d_err;
               double theta_r = theta_dmin+k_dmin*es;
               double ed_dot = vy*cos(phi_pre-theta_r)+vx*sin(phi_pre-theta_r);
               double e_phi = phi_pre-theta_r;
               double s_dot = (vx*cos(phi_pre-theta_r)-vy*sin(phi_pre-theta_r))/(1-k_dmin*ed);
               double e_phi_dot = phi_dot-k_dmin*s_dot;
               double point_curvature = k_dmin;

               std::array<double,5> err_k;
               err_k[0] = ed;
               err_k[1] = ed_dot;
               err_k[2] = e_phi;
               err_k[3] = e_phi_dot;
               err_k[4] = point_curvature;
               return err_k;
          }

          //LQR横向控制方法
          Eigen::Matrix<double,1,4> dLQR(Eigen::Matrix4d &A,Eigen::Matrix<double,4,1> &B,
          Eigen::Matrix4d &Q,Eigen::Matrix<double,1,1> &R){
               double dt=0.01;//离散时间步长
               Eigen::Matrix4d eye;
               eye.setIdentity(4,4);//单位矩阵
               Eigen::Matrix4d Abar = (eye-dt*0.5*A).inverse()*(eye+dt*0.5*A);//离散A和B（中点欧拉法）
               Eigen::Matrix<double,4,1> Bbar = B*dt;
               Eigen::Matrix4d P_old = Q;//初值
               int numloop = 5000;//迭代次数限制
               double minThre = 10e-7;//目标阈值
               for (int i=0;i<numloop;i++){
                    Eigen::Matrix4d P_new = Q+Abar.transpose()*P_old*Abar-
                    Abar.transpose()*P_old*Bbar*(R+Bbar.transpose()*P_old*Bbar).inverse()
                    *Bbar.transpose()*P_old*Abar;
                    if (fabs((P_new-P_old).maxCoeff())<=minThre || i==5000){//达到收敛条件退出迭代-求最大系数
                              P_old = P_new;
                              break;
                    }
                    P_old = P_new;
               }
               Eigen::Matrix<double,1,4> k = (R+Bbar.transpose()*P_old*Bbar).inverse()*Bbar.transpose()*P_old*Abar;
               return k;
          }
          Eigen::Matrix<double,1,4> cal_k(){
               Eigen::Matrix4d A;
               A << 0,1,0,0,
                    0,(Cf+Cr)/m/vx,-(Cf+Cr)/m,(a*Cf-b*Cr)/m/vx,
                    0,0,0,1,
                    0,(a*Cf-b*Cr)/Iz/vx,-(a*Cf-b*Cr)/Iz,(a*a*Cf+b*b*Cr)/Iz/vx;
               Eigen::Matrix<double,4,1> B;
               B << 0, -Cf/m, 0, -a*Cf/Iz;
               Eigen::Matrix4d Q;
               Q(0,0)=3.0;//ed重视程度
               Q(1,1)=1;//ed_dot
               Q(2,2)=10.0;//ephi，应加大对ed和ephi的跟踪
               Q(3,3)=1.0;//ephi_dot
               Eigen::Matrix<double,1,1> R;
               R(0,0) = 30.0;//控制成本
               Eigen::Matrix<double,1,4> k = dLQR(A,B,Q,R);
               return k;
          }

          //前馈控制角度
          double cal_forward_angle(Eigen::Matrix<double,1,4> &k,double point_curvature){
               double k3 = k[2];
               double forward_angle = point_curvature*(a+b-b*k3-m*vx*vx/(a+b)*(b/Cf+a/Cr*k3-a/Cr));
               return forward_angle;
          }

          //计算前轮控制转角
          double cal_angle_frontwheel(Eigen::Matrix<double,1,4> &k,std::array<double,5> &err_k,double &forward_angle){
               Eigen::Matrix<double,4,1> err;
               err << err_k[0],err_k[1],err_k[2],err_k[3];
               double angle_cmd = -k*err+forward_angle;
          }

          //LQR横向控制器
          double LateralControlbyLQR(std::vector<double> &PosX, std::vector<double> &PosY,std::vector<double> &theta){
               path_info_extract(PosX, PosY, theta);
               std::array<double,5> err_k = cal_err_k();
               Eigen::Matrix<double,1,4> k = cal_k();
               double forward_angle = cal_forward_angle(k, err_k[4]);
               double angle_cmd = cal_angle_frontwheel(k,err_k,forward_angle);
               if (!isnan(angle_cmd)){
                    return angle_cmd>1.0?std::min(angle_cmd,1.0):std::max(angle_cmd,-1.0);
               }
               return 0.0;
          }
          //PID横向控制器
          double LateralControlbyPID(std::vector<double> &PosX, std::vector<double> &PosY,std::vector<double> &theta){
               int index_y=0,index_t=0;
               for(int i=0;i<36;i++){
                    if(PosX[i]>=0.2*vx&&index_y==0){
                         index_y=i;
                    }
                    if(PosX[i]>=1.6*vx){
                         index_t=i;
                         break;
                    }
               }
               double err_s=PosY[index_y],err_th=theta[index_t];
               static double errs_sum=0,errs_last=0,errth_sum=0,errth_last=0;
               double kp_s=0.15,ki_s=0.001,kd_s=0.06,kp_th=2.5,ki_th=0.001,kd_th=2;
               double angle_cmd = kp_s*err_s+ki_s*errs_sum+kd_s*(err_s-errs_last)
                                   +kp_th*err_th+ki_th*errth_sum+kd_th*(err_th-errth_last);
               errs_sum+=err_s;
               errs_last=err_s;
               errth_sum+=err_th;
               errth_last=err_th;
               return !isnan(angle_cmd)?angle_cmd:0.0;
          }
          //Stanley横向控制器
          double LateralControlbyStanley(std::vector<double> &PosX, std::vector<double> &PosY,std::vector<double> &theta){
               double lambda=2.5;
               int index_y=0,index_t=0;
               for(int i=0;i<36;i++){
                    if(PosX[i]>=0.2*vx&&index_y==0){
                         index_y=i;
                    }
                    if(PosX[i]>=1.6*vx){
                         index_t=i;
                         break;
                    }
               }
               double err=fabs(theta[index_y])<10e-3?
                    (PosY[index_y]/tan(theta[index_y])-(PosX[index_y]-a))*sin(theta[index_y]):PosY[index_y];//近似认为参考点的角度等于切线角度
               double delta_phi=theta[index_t];
               double angle_cmd = atan2(lambda*err,vx)+1.7*delta_phi;//角度追踪加权
               return !isnan(angle_cmd)?angle_cmd:0.0;
          }

          //均值滤波器
          std::vector<double> meanfilter(double x1,double x2,double x3){
               static std::queue<double> que1,que2,que3;
               static double sum1,sum2,sum3;
               std::vector<double> res={filter(que1,sum1,x1),filter(que2,sum2,x2),filter(que3,sum3,x3)};
               return res;
          }
          double filter(std::queue<double> &que,double &sum,double &x){
               while(que.size()>5){
                    que.pop();
               }
               if(que.size()<5){
                    que.push(x);
                    sum+=x;
               }
               else{
                    sum-=que.front();
                    sum+=x;
                    que.pop();
                    que.push(x);
               }
               return double(sum/que.size());
          }
          
};