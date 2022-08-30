#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <unistd.h>
#include <typeinfo>
#include <time.h>
#include "matplotlibcpp.h"

namespace plt=matplotlibcpp;
using namespace std;

class node {
private:
    double x, y;                // 节点坐标
    vector<double> pathX, pathY;// 路径 
    node* parent;               // 父节点-指针
    double cost;
public:
    node(double _x, double _y);
    double getX();
    double getY();
    void setParent(node*);
    node* getParent();
};

class RRT {
private:
    node* startNode, * goalNode;          // 起始节点和目标节点
    vector<vector<double>> obstacleList;  // 障碍物
    vector<node*> nodeList;               // 节点
    double stepSize;                       // 步长

    int goal_sample_rate;

    // 随机函数产生的是一种伪随机数，它实际是一种序列发生器，有固定的算法，只有当种子不同时，序列才不同，
    // 所以不应该把种子固定在程序中，应该用随机产生的数做种子，如程序运行时的时间等。
    random_device goal_rd;                // random_device可以生成用来作为种子的随机的无符号整数值。
    mt19937 goal_gen;                     // mt19937是一种高效的随机数生成算法
    uniform_int_distribution<int> goal_dis;  //随机数源，随机数源调用随机数算法来生成随机数

    random_device area_rd;
    mt19937 area_gen;
    uniform_real_distribution<double> area_dis;
public:
    RRT(node*, node*, const vector<vector<double>>&, double , int);//构造函数还可以这样写？
    node* getNearestNode(const vector<double>&);
    bool collisionCheck(node*);
    vector<node*> planning();
};


// node构造函数，初始化x,y,parent,cost
node::node(double _x, double _y) : x(_x), y(_y), parent(nullptr), cost(0) {}

double node::getX() { return x; }
double node::getY() { return y; }

void node::setParent(node* _parent) { parent = _parent; }
node* node::getParent() { return parent; }

// RRT类构造函数
RRT::RRT(node* _startNode, node* _goalNode, const vector<vector<double>>& _obstacleList,
         double _stepSize = 1.0, int _goal_sample_rate = 5)
    : startNode(_startNode), goalNode(_goalNode),
      obstacleList(_obstacleList),
      stepSize(_stepSize), goal_sample_rate(_goal_sample_rate),
      goal_gen(goal_rd()), goal_dis(uniform_int_distribution<int>(0, 5000)),
      area_gen(area_rd())
      {}

node* RRT::getNearestNode(const vector<double>& randomPosition)
{
    int minID = -1;
    double minDistance = numeric_limits<double>::max(); // 编译器允许的double类型的最大值

    cout << "nodeLists.size:" << nodeList.size() << endl;

    // 找到和随机位置距离最小的节点,通过遍历所有点
    for (int i = 0; i < nodeList.size(); i++)
    {
        // 在这里距离不需要开根号
        double distance = pow(nodeList[i]->getX()-randomPosition[0], 2)+pow(nodeList[i]->getY()-randomPosition[1], 2);
        if (distance < minDistance)
        {
            minDistance = distance;    // 更新最小距离，这里的距离应该也可以采用曼哈顿距离或者其他条件判断
            minID = i;                 // 通过minID去记录下distance最小时对应的节点ID
        }
    }

    // 返回找到的距离randomPosition最近的点
    return nodeList[minID];
}

bool isLineIntersectRectangle(double linePointX1,double linePointY1,double linePointX2,double linePointY2,
double rectangleLeftTopX,double rectangleLeftTopY,double rectangleRightBottomX,double rectangleRightBottomY)
    {
        double lineHeight = linePointY1 - linePointY2;
        double lineWidth = linePointX2 - linePointX1;  // 计算叉乘 
        double c = linePointX1 * linePointY2 - linePointX2 * linePointY1;
        if ((lineHeight * rectangleLeftTopX + lineWidth * rectangleLeftTopY + c >= 0 && lineHeight * rectangleRightBottomX + lineWidth * rectangleRightBottomY + c <= 0)
            || (lineHeight * rectangleLeftTopX + lineWidth * rectangleLeftTopY + c <= 0 && lineHeight * rectangleRightBottomX + lineWidth * rectangleRightBottomY + c >= 0)
            || (lineHeight * rectangleLeftTopX + lineWidth * rectangleRightBottomY + c >= 0 && lineHeight * rectangleRightBottomX + lineWidth * rectangleLeftTopY + c <= 0)
            || (lineHeight * rectangleLeftTopX + lineWidth * rectangleRightBottomY + c <= 0 && lineHeight * rectangleRightBottomX + lineWidth * rectangleLeftTopY + c >= 0))
        {
 
            if (rectangleLeftTopX > rectangleRightBottomX)
            {
                double temp = rectangleLeftTopX;
                rectangleLeftTopX = rectangleRightBottomX;
                rectangleRightBottomX = temp;
            }
            if (rectangleLeftTopY < rectangleRightBottomY)
            {
                double temp1 = rectangleLeftTopY;
                rectangleLeftTopY = rectangleRightBottomY;
                rectangleRightBottomY = temp1;
            }
            if ((linePointX1 < rectangleLeftTopX && linePointX2 < rectangleLeftTopX)
                || (linePointX1 > rectangleRightBottomX && linePointX2 > rectangleRightBottomX)
                || (linePointY1 > rectangleLeftTopY && linePointY2 > rectangleLeftTopY)
                || (linePointY1 < rectangleRightBottomY && linePointY2 < rectangleRightBottomY))
            {
                return false;
            }
            else
            {
                return true;
            }
        }
        else
        {
            return false;
        }
 
    }
// 检测new节点到父节点的连线是否collision free
bool RRT::collisionCheck(node* newNode) {
    //边界
    if(newNode->getX()<=0||newNode->getX()>=20||newNode->getY()<=-1.5||newNode->getY()>=1.5)
        return false;
    for (auto item : obstacleList){
        double x1=item[0],y1=item[1],x2=item[2],y2=item[3];
        double Ax=newNode->getParent()->getX(),Ay=newNode->getParent()->getY(),Bx=newNode->getX(),By=newNode->getY();
        if(isLineIntersectRectangle(Ax,Ay,Bx,By,x2,y1,x1,y2))
            return false;
    }
    return true;
}

vector<node*> RRT::planning() {
    int count = 0;

    // 画出起始位置和目标位置
    plt::plot({startNode->getY()},{startNode->getX()},"r*");
    plt::plot({goalNode->getY()},{goalNode->getX()},"r*");
    // 画出障碍物
    for (auto item : obstacleList){
        double x1=item[0],y1=item[1],x2=item[2],y2=item[3];
        plt::plot(vector<double>{y1,y1},vector<double>{x1,x2},"b-");
        plt::plot(vector<double>{y1,y2},vector<double>{x2,x2},"b-");
        plt::plot(vector<double>{y2,y2},vector<double>{x1,x2},"b-");
        plt::plot(vector<double>{y1,y2},vector<double>{x1,x1},"b-");
    }

    // RRT core code
    nodeList.push_back(startNode); // 每次开始都首先在节点列表中添加起点节点
    while(1)
    {
        // 生成一个随机位置(这个随机位置不是直接作为新节点去使用的，只是树的生长方向)
        vector<double> randomPosition;//random_node的x，y坐标

        if(goal_dis(goal_gen)>goal_sample_rate)   // 这里可以优化成直接用节点来表示
        {
            area_dis=uniform_real_distribution<double>(0, 20);
            double randX = area_dis(goal_gen);        // 在(0,20)之间随机产生一个值作为x坐标
            area_dis=uniform_real_distribution<double>(-1.5, 1.5);
            double randY = area_dis(goal_gen);        // 在(-1.5, 1.5)之间随机产生一个值作为y坐标
            randomPosition.push_back(randX);
            randomPosition.push_back(randY);
        }
        else{ // 找到了目标,将目标位置保存--有一定的几率把随机点往终点方向引
            randomPosition.push_back(goalNode->getX());
            randomPosition.push_back(goalNode->getY());
        }

        // 找到和新生成随机节点距离最近的节点
        node* nearestNode = getNearestNode(randomPosition);
        // 利用反正切计算角度,然后利用角度和步长计算新坐标
        double theta = atan2(randomPosition[1] - nearestNode->getY(), randomPosition[0] - nearestNode->getX());
        // 利用之前采样的位置，加上设定的步长，来得到一个new节点
        node* newNode = new node(nearestNode->getX()+stepSize*cos(theta),nearestNode->getY()+stepSize*sin(theta));
        newNode->setParent(nearestNode);

        if (!collisionCheck(newNode)) continue;
        nodeList.push_back(newNode);

        // 画出路径
        plt::plot(vector<double>{newNode->getY(),nearestNode->getY()},vector<double>{newNode->getX(),nearestNode->getX()},"g-");
        count++;

        if (sqrt(pow(newNode->getX()-goalNode->getX(),2)+pow(newNode->getY()-goalNode->getY(),2))<=stepSize)
        {
            cout << "The path has been found!" << endl;
            break;
        }
    }

    // 画出最终得到的路径
    vector<node*> path;
    path.push_back(goalNode);
    node* tmpNode = nodeList.back(); //返回节点列表的最后一个元素
    while (tmpNode->getParent() != nullptr)
    {
        plt::plot(vector<double>{tmpNode->getY(),tmpNode->getParent()->getY()},vector<double>{tmpNode->getX(),tmpNode->getParent()->getX()},"k");
        path.push_back(tmpNode);
        tmpNode = tmpNode->getParent();
    }
    plt::show();
    // 展示背景
    path.push_back(startNode);
    return path;
}

int main(int argc, char* argv[]) 
{
    // 障碍物,矩形，左下和右上坐标//{x1,y1,x2,y2}
    vector<vector<double>> obstacleList{
        {3, -1.5, 5, 0.5},
        {12, -0.5, 14, 1.5}
    };

    // 起始点和目标点
    node* startNode = new node(0, -0.5);
    node* goalNode = new node(20, 1.0);

    clock_t start,finish;
    double totaltime;
    start=clock();
 
    RRT rrt(startNode, goalNode, obstacleList, 0.5, 5);//(~,~,~,step,goal_sample_rate(终点方向几率))
    rrt.planning();

    finish=clock();
    totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
    cout<<"time"<<totaltime<<endl;
    return 0;
}
