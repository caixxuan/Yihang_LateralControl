#include <stdio.h>
#include <qpOASES.hpp>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <bits/stdc++.h>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace std;
using namespace Eigen;

inline double norm(double x,double y) {return pow(x*x+y*y,0.5);}

void unit_process(vector<pair<double,double>> &nor){
        for(vector<pair<double,double>>::iterator it=nor.begin();it!=nor.end();it++){
                (*it).first=(*it).first/norm((*it).first,(*it).second);
                (*it).second=(*it).second/norm((*it).first,(*it).second);
        }
}

vector<double> kappa_req(vector<double> &x,vector<double> &y,int &n){
        vector<double> kappa(n);
        for(int i=1;i<n-1;i++){
                kappa[i]=norm(x[i+1]-2*x[i]+x[i-1],y[i+1]-2*y[i]+y[i-1])/norm(x[i+1]-x[i],y[i+1]-y[i]);
        }
        kappa[0]=kappa[1];
        kappa[n-1]=kappa[n-2];
        return kappa;
}

vector<vector<double>> read_csv(string &filename){
        ifstream inFile(filename);
        string lineStr;
        vector<vector<double>> strArray;//用来保存读取出来的数据，可以看成是一个二维数组，类型一般是string，其他类型可以转换
        // cout<<"the whole line is: "<<endl;
        while(getline(inFile,lineStr)) //这里的循环是每次读取一整行的数据,把结果保存在lineStr中，lineStr是用tab分割开的
        {
                //打印整行字符串
                // cout<<lineStr<<endl;
                //将结果保存为二维表结构
                stringstream ss(lineStr); //这里stringstream是一个字符串流类型，用lineStr来初始化变量 ss
                string str;
                vector<double> lineArray;
                //按照逗号进行分割
                while(getline(ss,str,'\t')){ //getline每次把按照逗号分割之后的每一个字符串都保存在str中
                        double num = stod(str);
                        lineArray.push_back(num); //这里将str保存在lineArray中
                }
                strArray.push_back(lineArray); //这里把lineArray保存在strArray。   这里的lineArray和lineArray有点类似于python中的list，只是固定了要保存的数据类型
        }
        // cout<<"--------------------------------------------"<<endl;
        // //打印一下我们保存的数据
        // cout<<"print what we have saved in the vector: "<<endl;
        // for(auto s : strArray)
        // {
        //         for(auto x : s){
        //                 cout<<x<< " ";
        //         }
        //         cout<<endl;
        // }
        return strArray;
}

int main(){
        //读取csv输入数据
        string filename_0 = "/home/caixuan/qptest/time/time18.5_results.csv";
        vector<vector<double>> RouOpt_data = read_csv(filename_0);//s, ub, lb, L
        int n=RouOpt_data.size();//DP粗优化轨迹点的个数
        vector<double> s(n),l(n),uub(n),llb(n);
        for(int i=0;i<n;i++){
                s[i]=RouOpt_data[i][0];
                uub[i]=RouOpt_data[i][1];
                llb[i]=RouOpt_data[i][2];
                l[i]=RouOpt_data[i][3];
        }
        vector<double> delta_s(n);
        for(int i=0;i<n;i++){
                delta_s[i]=s[i+1]-s[i];
        }
        delta_s[n-1]=delta_s[n-2];

        // // 绝对坐标系下的参考线（车道中心线）// //
        string filename_1 = "/home/caixuan/qptest/time/time18.5_ref_line_coeffs.csv";
        vector<vector<double>> RefLine_data = read_csv(filename_1);
        vector<double> c_x(4),c_y(4);//参考线的车身坐标系数
        for(int i=0;i<4;i++){
                c_x[i]=RefLine_data[1][i];
                c_y[i]=RefLine_data[1][i+4];
        }
        vector<double> x_ref(n),y_ref(n);//s对应的参考线绝对坐标
        for(int i=0;i<n;i++){
                x_ref[i]=c_x[0]+c_x[1]*s[i]+c_x[2]*s[i]*s[i]+c_x[3]*s[i]*s[i]*s[i];
                y_ref[i]=c_y[0]+c_y[1]*s[i]+c_y[2]*s[i]*s[i]+c_y[3]*s[i]*s[i]*s[i];
        }
        
        // //  测试输入 // // 
        /*
        int n=10;//需要优化的点
        //中心线(左负右正)
        vector<double> s_m={0,5,10,15,20,25,30,35,40,45},l_m(n,0);
        plt::named_plot("CENTER_LINE",l_m,s_m,"r--");
        //求单位法向量
        vector<pair<double,double>> tau(n-1),nor;//道路中心线的切向量和法向量
        for(int i=0;i<n-1;i++){
                tau[i]=make_pair(s_m[i+1]-s_m[i],l_m[i+1]-l_m[i]);
        }
        for(vector<pair<double,double>>::iterator it=tau.begin();it!=tau.end();it++){
                if(it==tau.begin()) nor.push_back(make_pair(-(*it).second,(*it).first));
                nor.push_back(make_pair(-(*it).second,(*it).first));
        }
        unit_process(nor);//单位化
        //求左右道路边界线坐标
        double Right[]={2,2,2,1,1,1,1,2,2,2},Left=2;
        vector<pair<double,double>> Nr(n),Nl(n);

        //假设粗优化的三阶多项式轨迹 s-y
        vector<double> s(n),y(n);//绝对坐标
        double dev[]={0.0,0.3,0.8,0.9,0.8,-0.1,-0.5,-1,-0.9};//假设粗优化点的横向偏移
        int i=0;
        vector<double> tmp_rx(n),tmp_ry(n),tmp_lx(n),tmp_ly(n);
        for(vector<pair<double,double>>::iterator itn=nor.begin(),itr=Nr.begin(),itl=Nl.begin();
        itn!=nor.end(),itr!=Nr.end(),itl!=Nl.end();itn++,itr++,itl++){
                (*itr).first=s_m[i]+(*itn).first*Right[i];
                (*itr).second=l_m[i]+(*itn).second*Right[i];
                (*itl).first=s_m[i]-(*itn).first*Left;
                (*itl).second=l_m[i]-(*itn).second*Left;
                tmp_rx[i]=(*itr).first;
                tmp_ry[i]=(*itr).second;
                tmp_lx[i]=(*itl).first;
                tmp_ly[i]=(*itl).second;
                s[i]=s_m[i]+(*itn).first*dev[i];
                y[i]=l_m[i]+(*itn).second*dev[i];
                i++;
        }
        
        //画图
        plt::plot(tmp_ry,tmp_rx,"g");
        plt::plot(tmp_ly,tmp_lx,"b");
        plt::named_plot("DP_LINE",y,s);
        
        //计算间距delta_s
        vector<double> delta_s(n);
        for(int i=1;i<n;i++){
                delta_s[i-1]=norm(s[i]-s[i-1],y[i]-y[i-1]);
        }
        delta_s[n-1]=delta_s[n-2];//最后一个近似等于

        // 完成准备工作 //
        */

        // 优化建模 //

        //需要优化的变量X={ n x l, n x l', n x l''} 一共n x 3个
        MatrixXd H;
        H.setZero(3*n,3*n);
        double wl=1,wl_p=50,wl_pp=100,wl_ppp=100;//平滑度权值
        double w_ref=1;//接近粗优化线权值
        for(int i=0;i<n;i++){
                H(i,i)=wl+w_ref;
                H(i+n,i+n)=wl_p;
                if(i==0 || i==(n-1)){
                        H(i+2*n,i+2*n)=wl_pp+wl_ppp/(delta_s[i]*delta_s[i]);
                }
                else H(i+2*n,i+2*n)=wl_pp+2*wl_ppp/(delta_s[i]*delta_s[i]);

                if((i+2*n+1) < 3*n){
                        H(i+2*n+1,i+2*n)=-1*wl_ppp/(delta_s[i]*delta_s[i]);
                        H(i+2*n,i+2*n+1)=-1*wl_ppp/(delta_s[i]*delta_s[i]);
                }
        }
        // H=2*H;

        vector<double> g(3*n,0);
        for(int i=0;i<n;i++){
                g[i]=-2*w_ref*l[i];
        }

        MatrixXd A;
        int M=3*(n-1);//不等式/等式限制数量
        A.setZero(M,3*n);
        
        for(int i=0,j=0;i<n-1 && j<n;i++,j++){
                A(i,2*n+j)=-1;
                A(i,2*n+j+1)=1;

                A(i+n-1,j)=-1;
                A(i+n-1,j+1)=1;

                A(i+2*(n-1),j+n)=-1;
                A(i+2*(n-1),j+n+1)=1;
                
                A(i+(n-1),j+n)=-delta_s[j];

                A(i+(n-1),j+2*n)=-1.0/3.0*delta_s[j]*delta_s[j];
                A(i+(n-1),j+2*n+1)=-1.0/6.0*delta_s[j]*delta_s[j];

                A(i+2*(n-1),j+2*n)=-0.5*delta_s[j];
                A(i+2*(n-1),j+2*n+1)=-0.5*delta_s[j];
        }
        vector<double> lbA(M,0),ubA(M,0);//优化变量联合等式约束
        double alpha_max=0.5,L=2.8,v=20;//最大转角，轴距，速度
        double jerk=alpha_max/L/v;//l'''的限制
        //运动学约束
        vector<double> kappa_s = kappa_req(x_ref,y_ref,n);//道路曲率
        for(int i=0;i<n-1;i++){
                lbA[i]=max(-jerk*delta_s[i],-(1/kappa_s[i]-abs(kappa_s[i])*L/tan(alpha_max)/kappa_s[i]));
                ubA[i]=min(jerk*delta_s[i],1/kappa_s[i]-abs(kappa_s[i])*L/tan(alpha_max)/kappa_s[i]);
        }
        vector<double> lb(3*n),ub(3*n);
        //lb = {lbs, ..., lowbound, ..., -kappa_bound+kappa_max, ...} kappa_bound道路曲率
        double lowbound_p=-2,upbound_p=2;//一阶导数 l' 上下限
        double kappa_max=tan(alpha_max)/L;//小汽车最小转弯半径一般为?
        
        for(int i=0;i<n;i++){
                lb[i]=llb[i];
                ub[i]=uub[i];
                lb[i+n]=lowbound_p;
                ub[i+n]=upbound_p;
                lb[i+2*n]=-kappa_s[i]-kappa_max;
                ub[i+2*n]=kappa_s[i]+kappa_max;
                cout<<kappa_s[i]<<' ';
        }
        lb[0]=l[0], ub[0]=l[0];//起终与终点施加硬约束
        lb[1]=l[1], ub[1]=l[1];
        // lb[n-1]=l[n-1],ub[n-1]=l[n-1];

        // QP求解 //
        USING_NAMESPACE_QPOASES;

        int N=3*n;
        real_t H_[N*N];
        for(int i=0;i<N;i++){
                for(int j=0;j<N;j++){
                        H_[i*N+j]=H(i,j);
                }
        }
	real_t g_[N];
        for(int i=0;i<N;i++){
                g_[i]=g[i];
        }
        real_t A_[M*N];
        for(int i=0;i<M;i++){
                for(int j=0;j<N;j++){
                        A_[i*N+j]=A(i,j);
                }
        }
        real_t lbA_[M];
        real_t ubA_[M];
        for(int i=0;i<M;i++){
                lbA_[i]=lbA[i];
                ubA_[i]=ubA[i];
        }
	real_t lb_[N];
        real_t ub_[N];
        for(int i=0;i<N;i++){
                lb_[i]=lb[i];
                ub_[i]=ub[i];
        }

	QProblem PathOpt( N,M );
        Options options;
	PathOpt.setOptions( options );

	int_t nWSR = 2000;
	PathOpt.init( H_,g_,A_,lb_,ub_,lbA_,ubA_, nWSR );

        real_t xOpt[N];
        PathOpt.getPrimalSolution( xOpt );

        vector<double> l_qp(n);
        for(int i=0;i<n;i++){
                l_qp[i]=xOpt[i];
                // cout<<l_qp[i]<<' ';
        }

        //frenet坐标下的DP(s,l)和QP(s,l_opt)优化轨迹扎转换到cartesian(x_dp,y_dp)&&(x_qp,y_qp)坐标系
        vector<double> theta_ref(n),x_dp(n),y_dp(n),x_qp(n),y_qp(n);
        for(int i=0;i<n-1;i++){
                theta_ref[i]=atan2(x_ref[i+1]-x_ref[i],y_ref[i+1]-y_ref[i]);
        }
        theta_ref[n-1]=theta_ref[n-2];
        for(int i=0;i<n;i++){
                x_dp[i]=x_ref[i]-l[i]*cos(theta_ref[i]);
                y_dp[i]=y_ref[i]+l[i]*sin(theta_ref[i]);
                x_qp[i]=x_ref[i]-l_qp[i]*cos(theta_ref[i]);
                y_qp[i]=y_ref[i]+l_qp[i]*sin(theta_ref[i]);
        }

        plt::named_plot("dp",l,s);
        plt::named_plot("qp",l_qp,s);
        plt::legend();
        plt::show();
        plt::axis("equal");

        plt::named_plot("ref",y_ref,x_ref,"--");
        plt::named_plot("dp",y_dp,x_dp);
        plt::named_plot("qp",y_qp,x_qp);
        plt::legend();
        plt::show();
        return 0;
}