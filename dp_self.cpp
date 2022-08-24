#include <stdio.h>
#include <vector>
#include "matplotlibcpp.h"

namespace plt=matplotlibcpp;

using namespace std;

struct scatter{
        double s,l;
        int j_index;
        bool free;
        scatter(){}
        scatter(double s_,double l_,int j_index_,bool free_):s(s_),l(l_),j_index(j_index_),free(free_){}//1表示空，0表示有障碍物阻挡
        bool operator==(const scatter &dot){
                return s==dot.s && l==dot.l && free==dot.free;
        }
        void operator=(const scatter &dot){
                s=dot.s;
                l=dot.l;
                j_index=dot.j_index;
                free=dot.free;
        }
};

class RoadCons{
public:
        vector<vector<scatter>> DisDot;
        vector<double> s,l;
        void Sow(){
                for(double m=0;m<=20;m++){
                        s.push_back(m);
                }
                for(double n=-1.5;n<=1.5;n+=0.5){
                        l.push_back(n);
                }
                for(int i=0;i<s.size();i++){
                        vector<scatter> onestage;
                        for(int j=0;j<l.size();j++){
                                if((3<=i && i<6 && 0<=j && j<=4) || (12<=i && i<15 && 2<=j && j<l.size())){//障碍物
                                        onestage.push_back({s[i],l[j],j,false});
                                        plt::plot({l[j]},{s[i]},"x");
                                }
                                else{
                                        onestage.push_back({s[i],l[j],j,true});
                                        plt::plot({l[j]},{s[i]},"+");
                                } 
                        }
                        DisDot.push_back(onestage);
                }
        }
        scatter start_point,end_point;
        void StartPoint(int i, int j){
                start_point={s[i],l[j],j,true};
                plt::plot({start_point.l},{start_point.s},"r*");
        }
        void EndPoint(int i, int j){
                end_point={s[i],l[j],j,true};
                plt::plot({end_point.l},{end_point.s},"r*");
        }
};

class DP_solve:public RoadCons{
public:
        double cost_onestep(scatter s1,scatter s2){
                double w_s=1,w_l=2;
                return pow((s1.s-s2.s)*(s1.s-s2.s)+(s1.l-s2.l)*(s1.l-s2.l),0.5)+w_l*abs(s1.l-s2.l);//自定义成本函数 曼哈顿距离与欧氏距离 综合
                // return pow((s1.s-s2.s)*(s1.s-s2.s)+(s1.l-s2.l)*(s1.l-s2.l),0.5);//欧氏距离
        }
        void dp_cal(){
                int s_len=DisDot.size();
                int l_len=DisDot[0].size();
                vector<vector<int>> best_last_j(s_len,vector<int>(l_len,0));//记录每个节点的最佳前一转移的横向序列
                vector<vector<double>> cost(s_len,vector<double>(l_len,10e6));//记录每一节点的从起点来的最小累积成本
                for(int i=0;i<s_len;i++){
                        for(int j=0;j<l_len;j++){
                                //free==0 障碍物和是起点直接剪枝
                                if(DisDot[i][j].free==false || i==0) {
                                        if(DisDot[i][j]==start_point) cost[i][j]=0;
                                        cout<<cost[i][j]<<' ';
                                        continue;
                                }
                                //第二排，限制从起点出发
                                if(i==1){
                                        cost[i][j]=cost_onestep(start_point,DisDot[i][j]);
                                        best_last_j[i][j]=start_point.j_index;
                                        cout<<cost[i][j]<<' ';
                                        continue;
                                }
                                //free==1
                                double cost_candidate=10e6;
                                for(int j_last=0;j_last<l_len;j_last++){
                                        //终端成本
                                        if(i==s_len-1){
                                                //到达终点
                                                if(DisDot[i][j]==end_point){
                                                        if(cost_onestep(DisDot[i-1][j_last],DisDot[i][j])+cost[i-1][j_last] < cost_candidate){
                                                                cost_candidate=cost_onestep(DisDot[i-1][j_last],DisDot[i][j])+cost[i-1][j_last];
                                                                best_last_j[i][j]=DisDot[i-1][j_last].j_index;
                                                        }
                                                }
                                                //未到达终点
                                                else{
                                                        cost_candidate=cost_onestep(DisDot[i-1][j_last],DisDot[i][j])+cost[i-1][j_last]+10e6; 
                                                }
                                                continue;
                                        }
                                        //中间情况
                                        if(cost_onestep(DisDot[i-1][j_last],DisDot[i][j])+cost[i-1][j_last] < cost_candidate){
                                                cost_candidate=cost_onestep(DisDot[i-1][j_last],DisDot[i][j])+cost[i-1][j_last];
                                                best_last_j[i][j]=DisDot[i-1][j_last].j_index;
                                        }
                                }
                                cost[i][j]=cost_candidate;
                                cout<<cost[i][j]<<' ';
                        }
                        cout<<endl;
                }

                //找到每一排中每个节点的最佳前一转移点，就是决策点
                vector<double> s_opt(s_len),l_opt(s_len);
                int cur_best_j;
                cout<<"------"<<endl;
                for(int i=0;i<s_len;i++){
                        for(int j=0;j<l_len;j++){
                                cout<<best_last_j[i][j]<<' ';
                        }
                        cout<<endl;
                }

                for(int i=s_len-1;i>=0;i--){
                        s_opt[i]=DisDot[i][0].s;
                        if(i==s_len-1){
                                l_opt[i]=end_point.l;
                                cur_best_j=best_last_j[i][end_point.j_index];
                                continue;
                        }
                        
                        l_opt[i]=DisDot[i][cur_best_j].l;
                        cur_best_j=best_last_j[i][cur_best_j];
                        if(i==0){
                                l_opt[i]=start_point.l;
                        }
                }
                plt::plot(l_opt,s_opt);
                plt::show();
        }
                
};
int main(){
        DP_solve dp;//对象创建

        dp.Sow();//播种撒点-加障碍物0
        dp.StartPoint(0, 2);//起点位置
        dp.EndPoint(20, 5);//终点位置

        dp.dp_cal();//求解

        return 0;
}