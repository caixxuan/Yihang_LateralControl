#include <stdio.h>
#include <vector>
#include <random>
#include <unistd.h>
#include <time.h>
#include "matplotlibcpp.h"

namespace plt=matplotlibcpp;
using namespace std;

class Data{
    public:
        double x;
        double y;
        int Class;//每个点所属的类😫
        Data(double ,double ,int);
        Data operator+(Data u){
            Data* o=new Data(0.0,0.0,0);
            o->x=x+u.x;
            o->y=y+u.y;
            o->Class=u.Class;
            return *o;
        }
        Data operator=(Data u){
            x=u.x;
            y=u.y;
            Class=u.Class;
            return *this;
        }
};

Data::Data(double x_,double y_,int Class_):x(0.0),y(0.0),Class(0){x=x_,y=y_,Class=Class_;}


class kmeans{
    private:
        int n;//聚类数
        int numbers;//数据数量
        vector<Data> data;//待聚类数据
        vector<Data> center;//聚类中心
        double goodness;//评价得分，越小越好
    public:
        kmeans(int, int, vector<Data> );//构造函数初始化
        void kmeans_init();//初始化聚类中心
        void classify();//分类
        void centroid();//计算质心
        void euc_eval();//总体评价
        void kmeans_process();//kmeans主流程
};

kmeans::kmeans(int n_,int numbers_,vector<Data> data_):n(n_),numbers(numbers_),data(data_) {}

void kmeans::kmeans_init(){
    center.resize(n,{0.0,0.0,0});
    random_device rd;                   // random_device可以生成用来作为种子的随机的无符号整数值。
    mt19937 rd_gen;                     // mt19937是一种高效的随机数生成算法
    uniform_int_distribution<int> rd_dis(0,numbers-1);  //随机数源，随机数源调用随机数算法来生成随机数
    for(int j=0;j<n;j++){
        center[j]=data[rd_dis(rd_gen)];
    }
}

double euclidean(const Data &a,const Data &b){
    return sqrt(pow(a.x-b.x,2)+pow(a.y-b.y,2));
}
void kmeans::classify(){
    for(int i=0;i<numbers;i++){
        double dis=euclidean(data[i],center[0]);
        data[i].Class=0;
        for(int j=1;j<n;j++){
            if(euclidean(data[i],center[j])<dis) data[i].Class=j;
        }
    }
    // for(int i=0;i<numbers;i++){
    //     if(data[i].Class==0) plt::plot({data[i].x},{data[i].y},"r*");
    //     else if(data[i].Class==1) plt::plot({data[i].x},{data[i].y},"b*");
    //     else plt::plot({data[i].x},{data[i].y},"k*");
    // }
}

void kmeans::centroid(){
    vector<Data> tmp(n,{0,0,0});
    vector<int> sizep(n,0);//每一类的数量
    for(int i=0;i<numbers;i++){
        for(int j=0;j<n;j++){
            if(data[i].Class==j) {
                Data t=tmp[j];
                tmp[j]=t+data[i];
                sizep[j]++;
            }
        }
    }
    for(int j=0;j<n;j++){
        center[j].x=tmp[j].x/sizep[j];
        center[j].y=tmp[j].y/sizep[j];
    }
}

void kmeans::euc_eval(){
    goodness=0;
    for(int i=0;i<numbers;i++){
        for(int j=0;j<n;j++){
            if(data[i].Class==j) goodness+=euclidean(data[i],center[j]);
        }
    }
}

void kmeans::kmeans_process(){
    kmeans_init();
    double last_goodness=10e6;
    goodness=10e6;
    int step=0;
    vector<double> cost;
    while(step==0||abs(last_goodness-goodness)>2){
        last_goodness=goodness;
        classify();
        centroid();
        euc_eval();
        cost.emplace_back(goodness);
        ++step;
    }
    // 绘制图像
    for(int i=0;i<numbers;i++){
        if(data[i].Class==0) plt::plot({data[i].x},{data[i].y},"r*");
        else if(data[i].Class==1) plt::plot({data[i].x},{data[i].y},"b*");
        else if(data[i].Class==2) plt::plot({data[i].x},{data[i].y},"k*");
    }
    plt::show();

    plt::named_plot("cost",cost);
    plt::legend();
    plt::show();
}

int main(){
    int numbers=100;
    vector<Data> data;
    //构造数据，假设3类
    random_device rd2;                   // random_device可以生成用来作为种子的随机的无符号整数值。
    mt19937 rd_gen2(rd2());                     // mt19937是一种高效的随机数生成算法
    uniform_real_distribution<double> rd_dis0(0,20);  //随机数源，随机数源调用随机数算法来生成随机数
    uniform_real_distribution<double> rd_dis1(45,65);
    uniform_real_distribution<double> rd_dis2(80,95);
    for(int i=0;i<numbers;i++){
        if(i<35){
            data.push_back({rd_dis0(rd_gen2),rd_dis0(rd_gen2),0});
            // plt::plot({data[i].x},{data[i].y},"o");
        }
        else if(i<75){
            data.push_back({rd_dis1(rd_gen2),rd_dis1(rd_gen2),0});
            // plt::plot({data[i].x},{data[i].y},"^");
        }
        else{
            data.push_back({rd_dis2(rd_gen2),rd_dis2(rd_gen2),0});
            // plt::plot({data[i].x},{data[i].y},"*");
        }
        // cout<<data[i].x<<endl;
    }

    kmeans K(3,numbers,data);
    K.kmeans_process();
    return 0;
}