/**
* This file is part of Fast-Planner.
*
* Copyright 2019 Boyu Zhou, Aerial Robotics Group, Hong Kong University of Science and Technology, <uav.ust.hk>
* Developed by Boyu Zhou <bzhouai at connect dot ust dot hk>, <uv dot boyuzhou at gmail dot com>
* for more information see <https://github.com/HKUST-Aerial-Robotics/Fast-Planner>.
* If you use this code, please cite the respective publications as
* listed on the above website.
*
* Fast-Planner is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* Fast-Planner is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with Fast-Planner. If not, see <http://www.gnu.org/licenses/>.
*/

#include <path_searching/kinodynamic_astar.h>
#include <sstream>
#include <plan_env/sdf_map.h>

using namespace std;
using namespace Eigen;

namespace fast_planner
{

/*
*@berif:  析构函数:删除所有路径中的点
*@param : void
*@return: void
*/
KinodynamicAstar::~KinodynamicAstar()
{
  for (int i = 0; i < allocate_num_; i++)
  {
    delete path_node_pool_[i];
  }
}

/*
*@berif:路径搜索在已知地图的前提下，给定起点和终点状态(位置和速度)实现混合A搜索
*       如果搜索成功返回一系列path_nodes_节点
*@param start_pt:     起始点的位置
*@param start_v:      起始点的速度
*@param start_a:      起始点的加速度
*@param end_pt:       终点的位置
*@param end_v:        终点的速度
*@param end_a:        终点的加速度
*@param init:         标志位,是否为第一次搜索
*@param dynamic:      标志位
*@param time_start:   搜索起始时间
*@return: 
*/
int KinodynamicAstar::search(Eigen::Vector3d start_pt, Eigen::Vector3d start_v, Eigen::Vector3d start_a,
                             Eigen::Vector3d end_pt, Eigen::Vector3d end_v, bool init, bool dynamic, double time_start)
{
  /*
   *主要是将起始点及目标点的三维位置转化至栅格地图的index, 并计算第一个扩展点的Heuristic cost
   */
  start_vel_ = start_v;                         // 起始点的速度
  start_acc_ = start_a;                         // 起始点的加速度

  PathNodePtr cur_node = path_node_pool_[0];    // 目前节点,在初始化中已经为搜索过程所用的节点分配好了内存
  cur_node->parent = NULL;                      // 起始点没有父节点
  cur_node->state.head(3) = start_pt;           // 当前节点的状态向量的前三个是起始点的位置
  cur_node->state.tail(3) = start_v;            // 当前节点的状态向量的后三个是起始点的速度
  cur_node->index = posToIndex(start_pt);       // 将起始点的位置转换为栅格地图的索引
  cur_node->g_score = 0.0;                      // 当前节点的g值

  Eigen::VectorXd end_state(6);                 // 终点的状态向量
  Eigen::Vector3i end_index;                    // 终点的位置对应的栅格地图的索引
  double time_to_goal;                          // 到底目标点的时间

  end_state.head(3) = end_pt;                   // 终点的状态向量,前三个为位置
  end_state.tail(3) = end_v;                    // 终点的状态向量,后三个为速度
  end_index = posToIndex(end_pt);               // 将终点的位置转换为栅格地图的索引
  // 求解f分数
  cur_node->f_score = lambda_heu_ * estimateHeuristic(cur_node->state, end_state, time_to_goal);
  cur_node->node_state = IN_OPEN_SET;           // 当前节点在开集中
  open_set_.push(cur_node);                     // 将当前节点放入开集
  use_node_num_ += 1;                           // 计算用过的节点数量

  // TODO:这是在干什么
  if (dynamic)
  {
    time_origin_ = time_start;
    cur_node->time = time_start;
    cur_node->time_idx = timeToIndex(time_start);
    expanded_nodes_.insert(cur_node->index, cur_node->time_idx, cur_node);
  }
  else
    expanded_nodes_.insert(cur_node->index, cur_node);      // 将起点加入expanded_nodes_

  PathNodePtr neighbor = NULL;                              // 邻近节点
  PathNodePtr terminate_node = NULL;                        // 终止节点
  bool init_search = init;                                  // 是否第一次搜索
  const int tolerance = ceil(1 / resolution_);              // 地图 resolution_ = 0.1, tolerance用来判断是否near_end,阈值为1m(10 index)

   /* ---------- search loop ---------- */
  while (!open_set_.empty())
  {
    /* ---------- get lowest f_score node ---------- */
    cur_node = open_set_.top();                             // 代价最小的是在开集最前面的节点,用于节点扩展

    /* ---------- determine termination ---------- */
    // TODO: Recording Horizon Local Planning  
    // 判断当前节点是否超出horizon,reach_horizon表示当前节点与起点的距离是否超出horizon范围
    bool reach_horizon = (cur_node->state.head(3) - start_pt).norm() >= horizon_;   // horizon_ = 9 m
    // 判断当前节点是否离终点较近
    bool near_end = abs(cur_node->index(0) - end_index(0)) <= tolerance &&
                    abs(cur_node->index(1) - end_index(1)) <= tolerance &&
                    abs(cur_node->index(2) - end_index(2)) <= tolerance;

    // 超出horizon范围或者在目标点附近
    if (reach_horizon || near_end)
    {
      terminate_node = cur_node;                // 当前节点为终止节点
      retrievePath(terminate_node);             // 路径存储在path_nodes_中

      if (near_end)
      {
        // Check whether shot traj exist
        estimateHeuristic(cur_node->state, end_state, time_to_goal);
        computeShotTraj(cur_node->state, end_state, time_to_goal);

        // TODO:说明目标点在起点附近, 不需进行搜索
        if (init_search)
          ROS_ERROR("Shot in first search loop!");
      }
    }
    
    // 超出horizon范围
    if (reach_horizon)
    {
      // is_shot_succ_ == false说明终点附近有障碍
      if (is_shot_succ_)
      {
        std::cout << "reach end" << std::endl;
        return REACH_END;
      }
      else
      {
        std::cout << "reach horizon" << std::endl;
        return REACH_HORIZON;
      }
    }

    // 在目标点附近
    if (near_end)
    {
      if (is_shot_succ_)
      {
        std::cout << "reach end" << std::endl;
        return REACH_END;
      }
      else if (cur_node->parent != NULL)
      {
        std::cout << "near end" << std::endl;
        return NEAR_END;
      }
      else
      {
        std::cout << "no path" << std::endl;
        return NO_PATH;
      }
    }

    /**************************** 节点扩张 ****************************/
    /* ---------- pop node and add to close set ---------- */
    open_set_.pop();                                      // 先top返回队首的值,然后再pop删除队首的值
    cur_node->node_state = IN_CLOSE_SET;                  // 将当前访问过的节点加入close_set
    iter_num_ += 1;                                       // 更新迭代次数

    // TODO: 时间是怎么进行处理的
    double res = 1 / 2.0, time_res = 1 / 1.0, time_res_init = 1 / 20.0;
    Eigen::Matrix<double, 6, 1> cur_state = cur_node->state;    // 当前状态
    Eigen::Matrix<double, 6, 1> pro_state;                      // 下一个时刻的状态                
    vector<PathNodePtr> tmp_expand_nodes;                       // tmp临时扩展节点
    Eigen::Vector3d um;                                         // 每一小段轨迹的输入
    double pro_t;                                               // 下一时刻
    vector<Eigen::Vector3d> inputs;                             // 输入:xyz轴的加速度ax,ay,az
    vector<double> durations;                                   // 轨迹前向积分的时间分辨率

    // 先按照之前的输入规划一段路径,若路径不满足条件,则划分加速度前向积分获得新的路径
    if (init_search)
    {
      inputs.push_back(start_acc_);                             // 将开始加速度存储储存到输入向量组
      for (double tau = time_res_init * init_max_tau_; 
           tau <= init_max_tau_ + 1e-3;
           tau += time_res_init * init_max_tau_)
        durations.push_back(tau);                               // tau: 每段轨迹前向积分的时间分辨率
      init_search = false;
    }
    // 按之前的加速度积分的轨迹不满足条件,细分加速度后重新规划
    else                                                  
    {
      for (double ax = -max_acc_; ax <= max_acc_ + 1e-3; ax += max_acc_ * res)
        for (double ay = -max_acc_; ay <= max_acc_ + 1e-3; ay += max_acc_ * res)
          for (double az = -max_acc_; az <= max_acc_ + 1e-3; az += max_acc_ * res)
          {
            um << ax, ay, az;                             // 每一小段轨迹的输入
            inputs.push_back(um);                         // 将加速度存储储存到输入向量组
          }
      for (double tau = time_res * max_tau_; tau <= max_tau_; tau += time_res * max_tau_)
        durations.push_back(tau);
    }

    for (int i = 0; i < inputs.size(); ++i)
      for (int j = 0; j < durations.size(); ++j)
      {
        um = inputs[i];
        double tau = durations[j];
        stateTransit(cur_state, pro_state, um, tau);
        pro_t = cur_node->time + tau;

        Eigen::Vector3d pro_pos = pro_state.head(3);        // 下一时刻的位置

        // Check if in close set 判断节点是否已经被扩展过
        Eigen::Vector3i pro_id = posToIndex(pro_pos);       // 下一个节点的位置id
        int pro_t_id = timeToIndex(pro_t);                  // 下一个节点的时间id
        PathNodePtr pro_node = dynamic ? expanded_nodes_.find(pro_id, pro_t_id) : expanded_nodes_.find(pro_id);
        if (pro_node != NULL && pro_node->node_state == IN_CLOSE_SET)
        {
          if (init_search)
            std::cout << "close" << std::endl;
          continue;
        }

        // Check maximal velocity 检查速度约束
        Eigen::Vector3d pro_v = pro_state.tail(3);          // 下一时刻的速度
        if (fabs(pro_v(0)) > max_vel_ || fabs(pro_v(1)) > max_vel_ || fabs(pro_v(2)) > max_vel_)
        {
          if (init_search)
            std::cout << "vel" << std::endl;
          continue;
        }

        // Check not in the same voxel 检查是否在相同栅格
        Eigen::Vector3i diff = pro_id - cur_node->index;
        int diff_time = pro_t_id - cur_node->time_idx;
        if (diff.norm() == 0 && ((!dynamic) || diff_time == 0))
        {
          if (init_search)
            std::cout << "same" << std::endl;
          continue;
        }

        // Check safety 碰撞检测
        Eigen::Vector3d pos;                            // 被检测时刻的位置
        Eigen::Matrix<double, 6, 1> xt;                 // 被检测时刻的状态
        bool is_occ = false;                            // 是否碰撞的标志位
        for (int k = 1; k <= check_num_; ++k)
        {
          double dt = tau * double(k) / double(check_num_);  
          stateTransit(cur_state, xt, um, dt);
          pos = xt.head(3);
          // 满足条件代表已经碰撞
          if (edt_environment_->sdf_map_->getInflateOccupancy(pos) == 1 )
          {
            is_occ = true;        
            break;
          }
        }
        if (is_occ)
        {
          if (init_search)
            std::cout << "safe" << std::endl;
          continue;
        }

        // 计算当前节点的g_score以及f_score
        double time_to_goal, tmp_g_score, tmp_f_score;
        tmp_g_score = (um.squaredNorm() + w_time_) * tau + cur_node->g_score;       // squaredNorm()向量的平方范数
        tmp_f_score = tmp_g_score + lambda_heu_ * estimateHeuristic(pro_state, end_state, time_to_goal);

        /****************************  节点剪枝 ****************************/
        // 判断当前临时扩展节点与current node的其他临时扩展节点是否在同一个voxel中,如果是的话，就要进行剪枝
        bool prune = false;                 // 标志位,是否进行剪枝
        for (int j = 0; j < tmp_expand_nodes.size(); ++j)
        {
          PathNodePtr expand_node = tmp_expand_nodes[j];        // 当前临时扩展节点
          if ((pro_id - expand_node->index).norm() == 0 && ((!dynamic) || pro_t_id == expand_node->time_idx))
          {
            prune = true;
            // 要判断当前临时扩展节点的fscore是否比同一个voxel的对比fscore
            // 如果是的话，则更新这一Voxel节点为当前临时扩展节点
            if (tmp_f_score < expand_node->f_score)
            {
              expand_node->f_score = tmp_f_score;
              expand_node->g_score = tmp_g_score;
              expand_node->state = pro_state;
              expand_node->input = um;
              expand_node->duration = tau;
              if (dynamic)
                expand_node->time = cur_node->time + tau;
            }
            break;
          }
        }

        // This node end up in a voxel different from others
        // 如果不剪枝的话，则首先判断当前临时扩展节点pro_node是否出现在open集中
        // 如果不是的话，则可以直接将pro_node加入open集中
        if (!prune)
        {
          if (pro_node == NULL)
          {
            pro_node = path_node_pool_[use_node_num_];
            pro_node->index = pro_id;
            pro_node->state = pro_state;
            pro_node->f_score = tmp_f_score;
            pro_node->g_score = tmp_g_score;
            pro_node->input = um;
            pro_node->duration = tau;
            pro_node->parent = cur_node;
            pro_node->node_state = IN_OPEN_SET;
            if (dynamic)
            {
              pro_node->time = cur_node->time + tau;
              pro_node->time_idx = timeToIndex(pro_node->time);
            }
            open_set_.push(pro_node);

            if (dynamic)
              expanded_nodes_.insert(pro_id, pro_node->time, pro_node);
            else
              expanded_nodes_.insert(pro_id, pro_node);

            tmp_expand_nodes.push_back(pro_node);

            use_node_num_ += 1;
            if (use_node_num_ == allocate_num_)
            {
              cout << "run out of memory." << endl;
              return NO_PATH;
            }
          }
          // 如果存在于open集但还未扩展的话，则比较当前临时扩展节点与对应VOXEL节点的fscore,若更小，则更新voxel中的节点。
          else if (pro_node->node_state == IN_OPEN_SET)
          {
            if (tmp_g_score < pro_node->g_score)
            {
              // pro_node->index = pro_id;
              pro_node->state = pro_state;
              pro_node->f_score = tmp_f_score;
              pro_node->g_score = tmp_g_score;
              pro_node->input = um;
              pro_node->duration = tau;
              pro_node->parent = cur_node;
              if (dynamic)
                pro_node->time = cur_node->time + tau;use_node_num_
            }
          }
          else
          {
            cout << "error type in searching: " << pro_node->node_state << endl;
          }
        }
      }
    // init_search = false;
  }

  cout << "open set empty, no path!" << endl;
  cout << "use node num: " << use_node_num_ << endl;
  cout << "iter num: " << iter_num_ << endl;
  return NO_PATH;
}

/*
*@berif:设置参数
*@param nh: ROS节点
*@return: void
*/
void KinodynamicAstar::setParam(ros::NodeHandle& nh)
{
  nh.param("search/max_tau", max_tau_, -1.0);                    // 如果考虑对时间维度进行划分才设置，这里未设置             
  nh.param("search/init_max_tau", init_max_tau_, -1.0);
  nh.param("search/max_vel", max_vel_, -1.0);                    // 速度限制
  nh.param("search/max_acc", max_acc_, -1.0);                    // 加速度限制
  nh.param("search/w_time", w_time_, -1.0);
  nh.param("search/horizon", horizon_, -1.0);                    // 限制全局规划的距离，保证实时性
  nh.param("search/resolution_astar", resolution_, -1.0);        // 空间分辨率
  nh.param("search/time_resolution", time_resolution_, -1.0);    // 时间维度分辨率
  nh.param("search/lambda_heu", lambda_heu_, -1.0);              // 启发函数权重
  nh.param("search/allocate_num", allocate_num_, -1);            // 最大节点数目
  nh.param("search/check_num", check_num_, -1);                  // 对中间状态安全检查
  nh.param("search/optimistic", optimistic_, true);
  tie_breaker_ = 1.0 + 1.0 / 10000;

  double vel_margin;
  nh.param("search/vel_margin", vel_margin, 0.0);
  max_vel_ += vel_margin;
}

/*
*@berif: 将起点到终点以及中间的节点添加到path_nodes_
*@param end_node: 终点
*@return: void
*/
void KinodynamicAstar::retrievePath(PathNodePtr end_node)
{
  // 将最后一个节点储存到path_nodes_,并作为path_nodes_第一个元素
  PathNodePtr cur_node = end_node;
  path_nodes_.push_back(cur_node);

  // 从end_node开始将父节点挨个储存到path_nodes_,直到起点
  while (cur_node->parent != NULL)
  {
    cur_node = cur_node->parent;
    path_nodes_.push_back(cur_node);
  }

  // 实现翻转路径数组,转为从起点开始到终点
  reverse(path_nodes_.begin(), path_nodes_.end());
}

/*
*@berif:利用庞特里亚金原理解决两点边值问题,得到最优解后用最优解的控制代价作为启发函数
*@param x1:
*@param x2:
*@param optimal_time: 优化时间的指针
*@return optimal_time: 最优控制时间
*/
double KinodynamicAstar::estimateHeuristic(Eigen::VectorXd x1, Eigen::VectorXd x2, double& optimal_time)
{
  const Vector3d dp = x2.head(3) - x1.head(3); 
  const Vector3d v0 = x1.segment(3, 3);                   // 取向量第3到第3+3个元素
  const Vector3d v1 = x2.segment(3, 3);                   // 取向量第3到第3+3个元素

  double c1 = -36 * dp.dot(dp);
  double c2 = 24 * (v0 + v1).dot(dp);
  double c3 = -4 * (v0.dot(v0) + v0.dot(v1) + v1.dot(v1));
  double c4 = 0;
  double c5 = w_time_;

  // 通过费拉里方法求解关于时间的一元四次方程
  std::vector<double> ts = quartic(c5, c4, c3, c2, c1);

  double v_max = max_vel_ * 0.5;
  double t_bar = (x1.head(3) - x2.head(3)).lpNorm<Infinity>() / v_max;
  ts.push_back(t_bar);

  double cost = 100000000;
  double t_d = t_bar;

  for (auto t : ts)
  {
    if (t < t_bar)
      continue;
    double c = -c1 / (3 * t * t * t) - c2 / (2 * t * t) - c3 / t + w_time_ * t;
    if (c < cost)
    {
      cost = c;
      t_d = t;
    }
  }

  optimal_time = t_d;

  return 1.0 * (1 + tie_breaker_) * cost;
}

/*
*@berif:利用庞特里亚金原理解一个两点边值问题
*@param:
*/
bool KinodynamicAstar::computeShotTraj(Eigen::VectorXd state1, Eigen::VectorXd state2, double time_to_goal)
{
  /* ---------- get coefficient ---------- */
  const Vector3d p0 = state1.head(3);
  const Vector3d dp = state2.head(3) - p0;
  const Vector3d v0 = state1.segment(3, 3);
  const Vector3d v1 = state2.segment(3, 3);
  const Vector3d dv = v1 - v0;
  double t_d = time_to_goal;
  MatrixXd coef(3, 4);
  end_vel_ = v1;

  Vector3d a = 1.0 / 6.0 * (-12.0 / (t_d * t_d * t_d) * (dp - v0 * t_d) + 6 / (t_d * t_d) * dv);
  Vector3d b = 0.5 * (6.0 / (t_d * t_d) * (dp - v0 * t_d) - 2 / t_d * dv);
  Vector3d c = v0;
  Vector3d d = p0;

  // 1/6 * alpha * t^3 + 1/2 * beta * t^2 + v0
  // a*t^3 + b*t^2 + v0*t + p0
  coef.col(3) = a, coef.col(2) = b, coef.col(1) = c, coef.col(0) = d;

  Vector3d coord, vel, acc;
  VectorXd poly1d, t, polyv, polya;
  Vector3i index;

  Eigen::MatrixXd Tm(4, 4);
  Tm << 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 3, 0, 0, 0, 0;

  /* ---------- forward checking of trajectory ---------- */
  double t_delta = t_d / 10;
  for (double time = t_delta; time <= t_d; time += t_delta)
  {
    t = VectorXd::Zero(4);
    for (int j = 0; j < 4; j++)
      t(j) = pow(time, j);

    for (int dim = 0; dim < 3; dim++)
    {
      poly1d = coef.row(dim);
      coord(dim) = poly1d.dot(t);
      vel(dim) = (Tm * poly1d).dot(t);
      acc(dim) = (Tm * Tm * poly1d).dot(t);

      // 速度、加速度不超限
      if (fabs(vel(dim)) > max_vel_ || fabs(acc(dim)) > max_acc_)
      {
        // cout << "vel:" << vel(dim) << ", acc:" << acc(dim) << endl;
        // return false;
      }
    }

    if (coord(0) < origin_(0) || coord(0) >= map_size_3d_(0) || coord(1) < origin_(1) || coord(1) >= map_size_3d_(1) ||
        coord(2) < origin_(2) || coord(2) >= map_size_3d_(2))
    {
      return false;
    }

    // if (edt_environment_->evaluateCoarseEDT(coord, -1.0) <= margin_) {
    //   return false;
    // }
    if (edt_environment_->sdf_map_->getInflateOccupancy(coord) == 1)
    {
      return false;
    }
  }
  coef_shot_ = coef;
  t_shot_ = t_d;
  is_shot_succ_ = true;
  return true;
}

/*
*@berif:求解一元三次方程
*@param a,b,c,d:一元三次方程的系数
*/
vector<double> KinodynamicAstar::cubic(double a, double b, double c, double d)
{
  vector<double> dts;

  double a2 = b / a;
  double a1 = c / a;
  double a0 = d / a;

  double Q = (3 * a1 - a2 * a2) / 9;
  double R = (9 * a1 * a2 - 27 * a0 - 2 * a2 * a2 * a2) / 54;
  double D = Q * Q * Q + R * R;
  if (D > 0)
  {
    double S = std::cbrt(R + sqrt(D));
    double T = std::cbrt(R - sqrt(D));
    dts.push_back(-a2 / 3 + (S + T));
    return dts;
  }
  else if (D == 0)
  {
    double S = std::cbrt(R);
    dts.push_back(-a2 / 3 + S + S);
    dts.push_back(-a2 / 3 - S);
    return dts;
  }
  else
  {
    double theta = acos(R / sqrt(-Q * Q * Q));
    dts.push_back(2 * sqrt(-Q) * cos(theta / 3) - a2 / 3);
    dts.push_back(2 * sqrt(-Q) * cos((theta + 2 * M_PI) / 3) - a2 / 3);
    dts.push_back(2 * sqrt(-Q) * cos((theta + 4 * M_PI) / 3) - a2 / 3);
    return dts;
  }
}

/*
*@berif:使用费拉里方法求解关于时间的一元四次方程
*@param a,b,c,d,e:一元四次方程的系数
*@return:
*/
vector<double> KinodynamicAstar::quartic(double a, double b, double c, double d, double e)
{
  vector<double> dts;

  double a3 = b / a;
  double a2 = c / a;
  double a1 = d / a;
  double a0 = e / a;

  // 一个元三次方程进行求解
  vector<double> ys = cubic(1, -a2, a1 * a3 - 4 * a0, 4 * a2 * a0 - a1 * a1 - a3 * a3 * a0);
  double y1 = ys.front();
  double r = a3 * a3 / 4 - a2 + y1;
  if (r < 0)
    return dts;

  double R = sqrt(r);
  double D, E;
  if (R != 0)
  {
    D = sqrt(0.75 * a3 * a3 - R * R - 2 * a2 + 0.25 * (4 * a3 * a2 - 8 * a1 - a3 * a3 * a3) / R);
    E = sqrt(0.75 * a3 * a3 - R * R - 2 * a2 - 0.25 * (4 * a3 * a2 - 8 * a1 - a3 * a3 * a3) / R);
  }
  else
  {
    D = sqrt(0.75 * a3 * a3 - 2 * a2 + 2 * sqrt(y1 * y1 - 4 * a0));
    E = sqrt(0.75 * a3 * a3 - 2 * a2 - 2 * sqrt(y1 * y1 - 4 * a0));
  }

  if (!std::isnan(D))
  {
    dts.push_back(-a3 / 4 + R / 2 + D / 2);
    dts.push_back(-a3 / 4 + R / 2 - D / 2);
  }
  if (!std::isnan(E))
  {
    dts.push_back(-a3 / 4 - R / 2 + E / 2);
    dts.push_back(-a3 / 4 - R / 2 - E / 2);
  }

  return dts;
}

void KinodynamicAstar::init()
{
  /* ---------- map params ---------- */
  this->inv_resolution_ = 1.0 / resolution_;
  inv_time_resolution_ = 1.0 / time_resolution_;
  edt_environment_->sdf_map_->getRegion(origin_, map_size_3d_);

  cout << "origin_: " << origin_.transpose() << endl;
  cout << "map size: " << map_size_3d_.transpose() << endl;

  /* ---------- pre-allocated node ---------- */
  path_node_pool_.resize(allocate_num_);
  for (int i = 0; i < allocate_num_; i++)
  {
    path_node_pool_[i] = new PathNode;
  }

  phi_ = Eigen::MatrixXd::Identity(6, 6);
  use_node_num_ = 0;
  iter_num_ = 0;
}

void KinodynamicAstar::setEnvironment(const EDTEnvironment::Ptr& env)
{
  this->edt_environment_ = env;
}

void KinodynamicAstar::reset()
{
  expanded_nodes_.clear();
  path_nodes_.clear();

  std::priority_queue<PathNodePtr, std::vector<PathNodePtr>, NodeComparator> empty_queue;
  open_set_.swap(empty_queue);

  for (int i = 0; i < use_node_num_; i++)
  {
    PathNodePtr node = path_node_pool_[i];
    node->parent = NULL;
    node->node_state = NOT_EXPAND;
  }

  use_node_num_ = 0;
  iter_num_ = 0;
  is_shot_succ_ = false;
  has_path_ = false;
}

/*
*@berif: 在完成路径搜索后按照预设的时间分辨率delta_t通过节点回溯和状态前向积分得到分辨率更高的路径点
*        如果最后的shot trajectory存在的话，则还要加上最后一段shot trajectory(即通过computeshottraj)算出来得到的
*@param delta_t:预设的时间分辨(原分辨率为duration,处理后的分辨率为delta_t)
*@return: 轨迹上采样点的pos
*/
std::vector<Eigen::Vector3d> KinodynamicAstar::getKinoTraj(double delta_t)
{
  vector<Vector3d> state_list;                  // 记录轨迹在采样点的pos

  /* ---------- get traj of searching ---------- */
  PathNodePtr node = path_nodes_.back();        // 从终点节点开始节点回溯,所以使用终点节点进行初始化
  Matrix<double, 6, 1> x0, xt;                  // pos_x, pos_y, pos_z, vel_x, vel_y, vel_z

  // 从终点开始向前推算记录轨迹在采样点的pos
  while (node->parent != NULL)
  {
    Vector3d ut = node->input;                  // 当前节点一小段轨迹上的input
    double duration = node->duration;           // 当前节点一小段轨迹上的duration
    x0 = node->parent->state;                   // 当前节点父节点的状态

    // 通过节点回溯和状态前向积分得到duration内,各个delta_t步长的状态
    for (double t = duration; t >= -1e-5; t -= delta_t)
    {
      stateTransit(x0, xt, ut, t);              // 推算上一个时刻的状态
      state_list.push_back(xt.head(3));         // 添加到列表进行记录
    }

    // 将节点设为原节点父节点并进行循环更新
    node = node->parent;
  }

  // 翻转state_list,使其从轨迹的起点开始
  reverse(state_list.begin(), state_list.end());

  /* ---------- get traj of one shot ---------- */
  /* ---------- 加上最后一段shot trajectory ----------*/
  if (is_shot_succ_)
  {
    Vector3d coord;
    VectorXd poly1d, time(4);

    for (double t = delta_t; t <= t_shot_; t += delta_t)
    {
      for (int j = 0; j < 4; j++)
        time(j) = pow(t, j);

      for (int dim = 0; dim < 3; dim++)
      {
        poly1d = coef_shot_.row(dim);
        coord(dim) = poly1d.dot(time);
      }
      state_list.push_back(coord);
    }
  }

  return state_list;
}

/*
*@berif:离散的获得一些轨迹点以及起始点的速度与加速度
*@param ts:
*@param point_set:
*@param start_end_derivatives:
*@return:
*/
void KinodynamicAstar::getSamples(double& ts, vector<Eigen::Vector3d>& point_set,
                                  vector<Eigen::Vector3d>& start_end_derivatives)
{
  /* ---------- path duration ---------- */
  double T_sum = 0.0;                       // 轨迹总时间
  if (is_shot_succ_)
    T_sum += t_shot_;
  PathNodePtr node = path_nodes_.back();
  while (node->parent != NULL)
  {
    T_sum += node->duration;
    node = node->parent;
  }
  // cout << "duration:" << T_sum << endl;

  // Calculate boundary vel and acc
  Eigen::Vector3d end_vel, end_acc;
  double t;
  if (is_shot_succ_)
  {
    t = t_shot_;
    end_vel = end_vel_;
    for (int dim = 0; dim < 3; ++dim)
    {
      Vector4d coe = coef_shot_.row(dim);
      end_acc(dim) = 2 * coe(2) + 6 * coe(3) * t_shot_;
    }
  }
  else
  {
    t = path_nodes_.back()->duration;
    end_vel = node->state.tail(3);
    end_acc = path_nodes_.back()->input;
  }

  // Get point samples
  int seg_num = floor(T_sum / ts);
  seg_num = max(8, seg_num);
  ts = T_sum / double(seg_num);
  bool sample_shot_traj = is_shot_succ_;
  node = path_nodes_.back();

  for (double ti = T_sum; ti > -1e-5; ti -= ts)
  {
    if (sample_shot_traj)
    {
      // samples on shot traj
      Vector3d coord;
      Vector4d poly1d, time;

      for (int j = 0; j < 4; j++)
        time(j) = pow(t, j);

      for (int dim = 0; dim < 3; dim++)
      {
        poly1d = coef_shot_.row(dim);
        coord(dim) = poly1d.dot(time);
      }

      point_set.push_back(coord);
      t -= ts;

      /* end of segment */
      if (t < -1e-5)
      {
        sample_shot_traj = false;
        if (node->parent != NULL)
          t += node->duration;
      }
    }
    else
    {
      // samples on searched traj
      Eigen::Matrix<double, 6, 1> x0 = node->parent->state;
      Eigen::Matrix<double, 6, 1> xt;
      Vector3d ut = node->input;

      stateTransit(x0, xt, ut, t);

      point_set.push_back(xt.head(3));
      t -= ts;

      // cout << "t: " << t << ", t acc: " << T_accumulate << endl;
      if (t < -1e-5 && node->parent->parent != NULL)
      {
        node = node->parent;
        t += node->duration;
      }
    }
  }
  reverse(point_set.begin(), point_set.end());

  // calculate start acc
  Eigen::Vector3d start_acc;
  if (path_nodes_.back()->parent == NULL)
  {
    // no searched traj, calculate by shot traj
    start_acc = 2 * coef_shot_.col(2);
  }
  else
  {
    // input of searched traj
    start_acc = node->input;
  }

  start_end_derivatives.push_back(start_vel_);
  start_end_derivatives.push_back(end_vel);
  start_end_derivatives.push_back(start_acc);
  start_end_derivatives.push_back(end_acc);
}

std::vector<PathNodePtr> KinodynamicAstar::getVisitedNodes()
{
  vector<PathNodePtr> visited;
  visited.assign(path_node_pool_.begin(), path_node_pool_.begin() + use_node_num_ - 1);
  return visited;
}

/*
 *@berif:将位置转换成索引
 *@param pt:当前位置
 *@return idx: 索引值
 */
Eigen::Vector3i KinodynamicAstar::posToIndex(Eigen::Vector3d pt)
{
  Vector3i idx = ((pt - origin_) * inv_resolution_).array().floor().cast<int>();
  return idx;
}

/*
 *@berif:将时间转换成索引
 *@param time:当前时间
 *@return idx: 索引值
 */
int KinodynamicAstar::timeToIndex(double time)
{
  int idx = floor((time - time_origin_) * inv_time_resolution_);      // time_origin_开始时间
  return idx;
}

/*
 *@berif:通过前向积分得到扩展节点的位置和速度
 *@param state0:转移前的状态变量
 *@param state1:转移后的状态变量
 *@param um:控制输入
 *@param tau:轨迹前向积分的时间分辨率
 *@return: void
 */
void KinodynamicAstar::stateTransit(Eigen::Matrix<double, 6, 1>& state0, Eigen::Matrix<double, 6, 1>& state1,
                                    Eigen::Vector3d um, double tau)
{
  for (int i = 0; i < 3; ++i)
    phi_(i, i + 3) = tau;

  Eigen::Matrix<double, 6, 1> integral;
  integral.head(3) = 0.5 * pow(tau, 2) * um;      // 加速度:a*t^2/2
  integral.tail(3) = tau * um;                    // 速度:a*t

  state1 = phi_ * state0 + integral;
}

}  // namespace fast_planner
