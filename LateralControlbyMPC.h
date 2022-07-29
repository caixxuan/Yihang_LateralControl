#pragma once
#include <stdio.h>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <qpOASES.hpp>
#include "math/czy_math_fit.hpp"

using namespace std;
using namespace czy;
USING_NAMESPACE_QPOASES;

class MPC{
//MPC横向控制器-基于运动学的
public:
          //参数赋值
          double m,Cf,Cr,a,b,Iz,vx,vy,wheelbase_;
          MPC(double speed_ms):vx(0.0){
               m=1500.0;
               Cf=-46085.0;
               Cr=-58915.0;
               a=1.016;
               b=1.562;
               Iz=2000.0;
               vx=speed_ms;
               vy=0.0;
               wheelbase_=a+b;//轴距
          }
          Eigen::MatrixXd H,g;
          double LateralControlbyMPC(std::vector<double> &PosX, std::vector<double> &PosY,std::vector<double> &theta_r,double &yawrate){
               //特殊情况
               if(PosX[0]==PosX[35]) return 0.0;
               
               //目标二次型
               double x_dmin = vx*0.1;
               int index=0;
               for(int i=0;i<36;i++){
                    if(PosX[i]>=vx*1.0){
                         index=i;
                         break;
                    }
               }
               double theta_dmin=(theta_r[index-1]+theta_r[index])/2;
               Fit Ft;
               Ft.polyfit(PosX,PosY,5,false);
               std::vector<double> f;
               Ft.getFactor(f);
               double y_dmin= f[0]+f[1]*x_dmin+f[2]*pow(x_dmin,2)+f[3]*pow(x_dmin,3)+f[4]*pow(x_dmin,4)+f[5]*pow(x_dmin,5);
               double dyDdx = f[1]+2.0*f[2]*x_dmin+3.0*f[3]*pow(x_dmin,2)+4.0*f[4]*pow(x_dmin,3)+5.0*f[5]*pow(x_dmin,4);//一阶导数
               double dyDdx2 = 2.0*f[2]+6.0*f[3]*x_dmin+12.0*f[4]*pow(x_dmin,2)+20.0*f[5]*pow(x_dmin,3);//二阶导数
               double curvature_ = dyDdx2/pow(1.0+dyDdx*dyDdx,1.5);
               if(isnan(curvature_)) curvature_ = 0.0;//曲率
               //double delta_r = atan(wheelbase_ * curvature_);//参考方向角
               double delta_r = theta_dmin;
               int Nx=4,Nu=1,Np=5,Nc=3;//状态量数量，控制数量，预测长度，控制长度
               Eigen::MatrixXd A,B,C,Q,R,X0;
               A.resize(Nx,Nx);
               A << 0,1,0,0,
                    0,(Cf+Cr)/m/vx,-(Cf+Cr)/m,(a*Cf-b*Cr)/m/vx,
                    0,0,0,1,
                    0,(a*Cf-b*Cr)/Iz/vx,-(a*Cf-b*Cr)/Iz,(a*a*Cf+b*b*Cr)/Iz/vx;
               B.resize(Nx,Nu);
               B << 0, -Cf/m, 0, -a*Cf/Iz;
               C.resize(Nx,Nu);
               C << 0, (Cf*a-Cr*b)/m/vx-vx, 0, (Cr*b*b+Cf*a*a)/Iz/vx;
               double dt=0.1;//预测离散时间步长
               Eigen::MatrixXd eye;
               eye.setIdentity(4,4);//单位矩阵
               Eigen::MatrixXd Abar,Bbar,Cbar;
               Abar.resize(Nx,Nx);
               Abar = (eye-dt*0.5*A).inverse()*(eye+dt*0.5*A);//离散A和B（中点欧拉法）
               Bbar.resize(Nx,1);
               Bbar = B*dt;
               Cbar.resize(Nx,1);
               Cbar = dt*C*(yawrate);//vx/curvature_ 是横摆角速度 

               //计算初始状态
               X0=cal_err_k(x_dmin,y_dmin,curvature_,theta_dmin);
               Q.setZero(Nx,Nx);
               Q(0,0)=3.0;//状态重视度矩阵
               Q(1,1)=1.0;
               Q(2,2)=10.0;
               Q(3,3)=1.0;
               R.resize(Nu,Nu);
               R<<30.0;//控制成本
               double u_k;
               MPC_Matrices(Abar,Bbar,Cbar,Q,R,Np,X0);//计算H,g
               u_k=Prediction(Np,Nu,Nx);
               return isnan(u_k)?0.0:u_k;                                       
          }
          Eigen::MatrixXd cal_err_k(double &x_dmin, double &y_dmin, double &k_dmin, double &theta_dmin){
               Eigen::Matrix<double,2,1> nor;
               nor << -sin(theta_dmin),cos(theta_dmin);
               Eigen::Matrix<double,2,1> tau;
               tau << cos(theta_dmin),sin(theta_dmin);
               Eigen::Matrix<double,2,1> d_err;
               double tp = 0.0;//预测时间
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

               Eigen::MatrixXd err_k;
               err_k.resize(4,1);
               err_k<< ed,ed_dot,e_phi,e_phi_dot;
               return err_k;
          }

          void MPC_Matrices(Eigen::MatrixXd &A,Eigen::MatrixXd &B,Eigen::MatrixXd &C,Eigen::MatrixXd &Q,
          Eigen::MatrixXd &R,int &Np,Eigen::MatrixXd &X0){
               int Nx=A.rows(),Nu=B.cols();
               Eigen::MatrixXd aa, aa_l,aa_tmp;
               aa.resize(Np*Nx,1);
               aa_l.resize(Np*Nx,Nx);
               aa_tmp.setIdentity(Nx,Nx);
               for(int i=0;i<Np;i++){
                    aa_l.middleRows(i*Nx,Nx)=A*aa_tmp;
                    aa_tmp*=A;
               }
               aa=aa_l*X0;
               
               Eigen::MatrixXd kk,tmp,cc,cc_tmp;
               kk.setZero(Np*Nx,Np*Nu);
               tmp=A;
               kk.block(0,0,Nx,Nu)=B;
               for(int i=1;i<Np;i++){
                    Eigen::MatrixXd temp;
                    temp.resize(Nx,Np*Nu);
                    temp.block(0,0,Nx,Nu)=tmp*B;
                    temp.block(0,Nu,Nx,(Np-1)*Nu)=kk.middleRows((i-1)*Nx,Nx).leftCols((Np-1)*Nu);
                    kk.middleRows(i*Nx,Nx)=temp;
                    tmp=A*tmp;
               }
               
               cc.setZero(Np*Nx,1);
               cc_tmp=A;
               cc.topRows(Nx)=C;
               for(int i=1;i<Np;i++){
                    cc.middleRows(i*Nx,Nx)=cc_tmp*C+cc.middleRows((i-1)*Nx,Nx);
                    cc_tmp*=A;
               }
               
               Eigen::MatrixXd Q_bar,R_bar;
               Q_bar.resize(Np*Nx,Np*Nx);
               Eigen::MatrixXd eye;
               eye.setIdentity(Np,Np);
               Q_bar=Kroneck(eye,Q);
               R_bar.resize(Np*Nu,Np*Nu);
               R_bar=Kroneck(eye,R);
               H.resize(Np*Nu,Np*Nu);
               H=kk.transpose()*Q_bar*kk+R_bar;//Np*Nu
               g.resize(1,Np*Nu);
               g=2*(aa+cc).transpose()*Q_bar*kk;//1*(Np*Nu)
          }

          double Prediction(int &Np,int &Nu,int &Nx){
               int m_H=Np*Nu,n_H=Np*Nu;

               real_t H_[m_H*n_H],g_[m_H];
          
               for(int i=0;i<m_H;i++){
                    for(int j=0;j<n_H;j++){
                         H_[i*n_H+j]=H(i,j);
                    }
               }
               
               for(int i=0;i<m_H;i++){
                    g_[i]=g(0,i);
               }
               
               real_t A_[1*m_H]={1.0, 1.0, 1.0, 1.0, 1.0};
               real_t lbA[1]={-1000.0};
               real_t ubA[1]={1000.0};
               real_t lb[m_H]={-3.14/2, -3.14/2, -3.14/2, -3.14/2, -3.14/2};
               real_t ub[m_H]={3.14/2, 3.14/2, 3.14/2, 3.14/2, 3.14/2};
               
               // Setting up QProblem object.
	          QProblemB qp( m_H );

	          Options options;
	          qp.setOptions( options );
               options.setToMPC();
               //qp.setPrintLevel( PL_NONE );

	          // Solve first QP. 
	          int_t nWSR = 100;//迭代数量
               real_t cputime=0.5;//限制时间
	          qp.init( H_,g_,lb,ub, nWSR, &cputime );

               //return solution
               real_t xOpt[Nu];
               qp.getPrimalSolution( xOpt );
               double res=xOpt[0];
               
               return res;
          }
          Eigen::MatrixXd Kroneck(Eigen::MatrixXd &a, Eigen::MatrixXd &b){
               int i=a.rows(),j=a.cols(),m=b.rows(),n=b.cols();
               Eigen::MatrixXd c;
               c.resize(i*m,j*n);
               int index_x = c.rows();
               int index_y = c.cols();
               for (int x = 0;x<index_x;x++)
               {
                    for (int y = 0; y < index_y; y++)
                    {
                         int a_i = x / m;
                         int a_j = y / n;
                         int b_i = x % m;
                         int b_j = y % n;
                         c(x,y) = a(a_i,a_j)*b(b_i,b_j);
                    }
               }
               return c;
          }
};
