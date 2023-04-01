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

#include "bspline/non_uniform_bspline.h"
#include <ros/ros.h>

namespace fast_planner {

NonUniformBspline::NonUniformBspline(const Eigen::MatrixXd& points, const int& order,
                                     const double& interval) {
  setUniformBspline(points, order, interval);
}

NonUniformBspline::~NonUniformBspline() {}

/*
*@berif: 均匀B样条设置
*        获得控制点，轨迹次数，以及时间间隔的情况下，设置时间区间(Knot vector)
*@param points:   控制点
*@param order:    轨迹次数
*@param interval: 时间间隔
*@return: void
*/
void NonUniformBspline::setUniformBspline(const Eigen::MatrixXd& points, const int& order,
                                          const double& interval) {
  control_points_ = points;         // 控制点
  p_              = order;          // B-spline次数
  interval_       = interval;       // 时间间隔

  n_ = points.rows() - 1;           // 轨迹段数
  m_ = n_ + p_ + 1;                 // 涉及的节点数()

  /* ---------------  生成节点向量 --------------- */
  u_ = Eigen::VectorXd::Zero(m_ + 1);         // 定义节点向量(knot vector)
  for (int i = 0; i <= m_; ++i) {

    if (i <= p_) {                            // 前p+1个节点, 好像没什么用? 只要保证up = 0即可
      u_(i) = double(-p_ + i) * interval_;
    } else if (i > p_ && i <= m_ - p_) {      // 中间n-p个节点, up+1 = ts , up+2 = 2 * ts
      u_(i) = u_(i - 1) + interval_;
    } else if (i > m_ - p_) {                 // 后p+1个节点
      u_(i) = u_(i - 1) + interval_;
    }
  }
}

/*
*@berif:设置对应NonUniformBspline类的节点向量
*@param:设置NonUniformBspline类的节点向量
*@return:void
*/
void NonUniformBspline::setKnot(const Eigen::VectorXd& knot) { this->u_ = knot; }

/*
*@berif:获得对应NonUniformBspline类的节点向量
*@param:void
*@return:对应NonUniformBspline类的节点向量
*/
Eigen::VectorXd NonUniformBspline::getKnot() { return this->u_; }

/*
*@berif:B样条曲线的时间区间
*@param um:     时间区间下限
*@param um_p:   时间区间上限
*@return:       void
*/
void NonUniformBspline::getTimeSpan(double& um, double& um_p) {
  um   = u_(p_);
  um_p = u_(m_ - p_);
}

/*
TODO:
*@berif:
*@param:void
*@return:
*/
Eigen::MatrixXd NonUniformBspline::getControlPoint() { return control_points_; }

pair<Eigen::VectorXd, Eigen::VectorXd> NonUniformBspline::getHeadTailPts() {
  Eigen::VectorXd head = evaluateDeBoor(u_(p_));
  Eigen::VectorXd tail = evaluateDeBoor(u_(m_ - p_));
  return make_pair(head, tail);
}

/*
*@berif: B样条函数值计算
*        直接得到一个[t_{p} , t_{m−p}]作用域中的B样条函数值
*@param u: 时间间隔
*@return:  void
*/
Eigen::VectorXd NonUniformBspline::evaluateDeBoor(const double& u) {

  double ub = min(max(u_(p_), u), u_(m_ - p_));

  // determine which [ui,ui+1] lay in
  int k = p_;
  while (true) {
    if (u_(k + 1) >= ub) break;
    ++k;
  }

  /* deBoor's alg */
  vector<Eigen::VectorXd> d;
  for (int i = 0; i <= p_; ++i) {
    d.push_back(control_points_.row(k - p_ + i));
    // cout << d[i].transpose() << endl;
  }

  for (int r = 1; r <= p_; ++r) {
    for (int i = p_; i >= r; --i) {
      double alpha = (ub - u_[i + k - p_]) / (u_[i + 1 + k - r] - u_[i + k - p_]);
      // cout << "alpha: " << alpha << endl;
      d[i] = (1 - alpha) * d[i - 1] + alpha * d[i];
    }
  }

  return d[p_];
}

Eigen::VectorXd NonUniformBspline::evaluateDeBoorT(const double& t) {
  return evaluateDeBoor(t + u_(p_));
}

/*
*@berif:计算非均匀形式下的速度控制点及加速度控制点
*       用均匀B样条的控制点算是因为均匀B样条可以看做特殊形式的非均匀B样条
*@param: void
*@return:速度控制点
*/
Eigen::MatrixXd NonUniformBspline::getDerivativeControlPoints() {
  // B样条曲线的导数也是B样条曲线,其阶数变为了p_-1
  // 控制点 Qi = p_*(Pi+1-Pi)/(ui+p_+1-ui+1)
  Eigen::MatrixXd ctp = Eigen::MatrixXd::Zero(control_points_.rows() - 1, control_points_.cols());
  for (int i = 0; i < ctp.rows(); ++i) {
    ctp.row(i) =
        p_ * (control_points_.row(i + 1) - control_points_.row(i)) / (u_(i + p_ + 1) - u_(i + 1));
  }
  return ctp;
}

NonUniformBspline NonUniformBspline::getDerivative() {
  Eigen::MatrixXd   ctp = getDerivativeControlPoints();
  // 新定义一个NonUniformBspline对象,并将新获得一阶微分控制点,次数,Knot vector赋值给它
  NonUniformBspline derivative(ctp, p_ - 1, interval_);   // 初始化时将u设置为默认值

  /* 剪掉两个节点(u_--->knot) */
  // B样条的一阶微分是次数(p)-1，控制点数(n)-1的B样条曲线,因此相应的Knot vector-2
  // B-spline求导数后控制点个数(n)-1, 次数(p)-1, 因为区间范围是[U0, Un+p+1], 所以节点数-2
  // 由Vi = p_ * (Qi+1 - Qi) / (Ui+p+1 - ui+1)可知, V0 ~ Vn-1与收尾两个节点(U0, Un+p+1)无关
  Eigen::VectorXd knot(u_.rows() - 2);
  knot = u_.segment(1, u_.rows() - 2);  // 取向量第1个到第u_.rows()-2元素
  derivative.setKnot(knot);

  return derivative;
}

/*
*@berif:           设置节点间隔
*@param interval_: 节点间隔
*@return:           void
*/
double NonUniformBspline::getInterval() { return interval_; }

/*
*@berif:设置物理限制(速度限制、加速度限制、时间调整系数)
*@param vel: 速度限制
*@param acc: 加速度限制
*@return:    void
*/
void NonUniformBspline::setPhysicalLimits(const double& vel, const double& acc) {
  limit_vel_   = vel;
  limit_acc_   = acc;
  limit_ratio_ = 1.1;
}

/*
*@berif:根据公式计算每个控制点的速度和加速度是否超限,并求解最大速度和加速度
*@param show: 标志位:是否要显示检查结果
*@return:     标志位:是否可达
*/
bool NonUniformBspline::checkFeasibility(bool show) {
  bool fea = true;
  // SETY << "[Bspline]: total points size: " << control_points_.rows() << endl;

  Eigen::MatrixXd P         = control_points_;            // 定义控制点
  int             dimension = control_points_.cols();     // 控制点维度

  /* 检查速度可达性并插入点 */
  double max_vel = -1.0;                                  // 最大速度
  for (int i = 0; i < P.rows() - 1; ++i) {
    // 根据公式 Vi = p_ * (Qi+1 - Qi) / (Ui+p+1 - ui+1) 求节点速度
    Eigen::VectorXd vel = p_ * (P.row(i + 1) - P.row(i)) / (u_(i + p_ + 1) - u_(i + 1));

    // 如果速度超过速度限制,则不可达,并更新最大速度
    if (fabs(vel(0)) > limit_vel_ + 1e-4 || fabs(vel(1)) > limit_vel_ + 1e-4 ||
        fabs(vel(2)) > limit_vel_ + 1e-4) {

      if (show) cout << "[Check]: Infeasible vel " << i << " :" << vel.transpose() << endl;
      fea = false;          // 标志位:不满足速度要求,设置为不可达

      // 通过循环求解最大速度
      for (int j = 0; j < dimension; ++j) {
        max_vel = max(max_vel, fabs(vel(j)));
      }
    }
  }

  /* 检查加速度可达性并插入点 */
  double max_acc = -1.0;                                  // 最大加速度
  for (int i = 0; i < P.rows() - 2; ++i) {
    // 根据公式 (17) 求节点加速度
    Eigen::VectorXd acc = p_ * (p_ - 1) *
        ((P.row(i + 2) - P.row(i + 1)) / (u_(i + p_ + 2) - u_(i + 2)) -
         (P.row(i + 1) - P.row(i)) / (u_(i + p_ + 1) - u_(i + 1))) /
        (u_(i + p_ + 1) - u_(i + 2));

    // 如果加速度超过加速度限制,则不可达,并更新最大加速度
    if (fabs(acc(0)) > limit_acc_ + 1e-4 || fabs(acc(1)) > limit_acc_ + 1e-4 ||
        fabs(acc(2)) > limit_acc_ + 1e-4) {

      if (show) cout << "[Check]: Infeasible acc " << i << " :" << acc.transpose() << endl;
      fea = false;          // 标志位:不满足假速度要求,设置为不可达

      // // 通过循环求解最大加速度
      for (int j = 0; j < dimension; ++j) {
        max_acc = max(max_acc, fabs(acc(j)));
      }
    }
  }

  // 若vel, acc都不满足要求, ratio要保证调整后两者都符合要求
  double ratio = max(max_vel / limit_vel_, sqrt(fabs(max_acc) / limit_acc_));

  return fea;
}

/*
*@berif:根据最大速度和最大加速度调节
*@param:  void
*@return: 获得调整后的ratio
*/
double NonUniformBspline::checkRatio() {
  Eigen::MatrixXd P         = control_points_;              // 定义控制点
  int             dimension = control_points_.cols();       // 控制点维度

  // 求解最大速度
  double max_vel = -1.0;
  for (int i = 0; i < P.rows() - 1; ++i) {
    Eigen::VectorXd vel = p_ * (P.row(i + 1) - P.row(i)) / (u_(i + p_ + 1) - u_(i + 1));
    for (int j = 0; j < dimension; ++j) {
      max_vel = max(max_vel, fabs(vel(j)));
    }
  }
 // 求解最大加速度
  double max_acc = -1.0;
  for (int i = 0; i < P.rows() - 2; ++i) {
    Eigen::VectorXd acc = p_ * (p_ - 1) *
        ((P.row(i + 2) - P.row(i + 1)) / (u_(i + p_ + 2) - u_(i + 2)) -
         (P.row(i + 1) - P.row(i)) / (u_(i + p_ + 1) - u_(i + 1))) /
        (u_(i + p_ + 1) - u_(i + 2));
    for (int j = 0; j < dimension; ++j) {
      max_acc = max(max_acc, fabs(acc(j)));
    }
  }

  // 求解ratio
  double ratio = max(max_vel / limit_vel_, sqrt(fabs(max_acc) / limit_acc_));
  ROS_ERROR_COND(ratio > 2.0, "max vel: %lf, max acc: %lf.", max_vel, max_acc);

  return ratio;
}

/*
*@berif:  计算当前控制点是否超限，以及调整表比例
*@param show:  标志位:是否向命令行打印相关信息
*@return: 标志位:是否可达
*/
bool NonUniformBspline::reallocateTime(bool show) {
  // SETY << "[Bspline]: total points size: " << control_points_.rows() << endl;
  // cout << "origin knots:\n" << u_.transpose() << endl;
  bool fea = true;                                      // 定义标志位,判断是否可达

  Eigen::MatrixXd P         = control_points_;          // 控制点
  int             dimension = control_points_.cols();   // 控制点的维度

  double max_vel, max_acc;                              // 最大速度 最大加速度

  /* 检查速度可达性并插入点 */
  for (int i = 0; i < P.rows() - 1; ++i) {
    Eigen::VectorXd vel = p_ * (P.row(i + 1) - P.row(i)) / (u_(i + p_ + 1) - u_(i + 1));

    if (fabs(vel(0)) > limit_vel_ + 1e-4 || fabs(vel(1)) > limit_vel_ + 1e-4 ||
        fabs(vel(2)) > limit_vel_ + 1e-4) {

      fea = false;
      if (show) cout << "[Realloc]: Infeasible vel " << i << " :" << vel.transpose() << endl;

      max_vel = -1.0;
      for (int j = 0; j < dimension; ++j) {
        max_vel = max(max_vel, fabs(vel(j)));
      }
 
      // 防止μv和μa过大, 影响其它控制点, 这里引入limit_ratio_来做一个限制
      // 但这样不能保证一次调整后V,A一定满足要求, 因此需迭代调整.
      double ratio = max_vel / limit_vel_ + 1e-4;
      // 对ratio进行限制
      if (ratio > limit_ratio_) ratio = limit_ratio_;

      /* 
       * 对于当前控制点i有关的时间区间(t_{i},t_{i+pb+1})进行时间调整
       * 
       */
      double time_ori = u_(i + p_ + 1) - u_(i + 1);       // 原Ui+p+1 - Ui+1的时间
      double time_new = ratio * time_ori;                 // 调整后的时间
      double delta_t  = time_new - time_ori;              // 求解总增量
      double t_inc    = delta_t / double(p_);             // 时间调整共涉及p个节点区间, 这里将总增量均匀地加到各个区间上

      for (int j = i + 2; j <= i + p_ + 1; ++j) {         // 将增量加到 Ui+2 ~ Ui+p+1 上
        u_(j) += double(j - i - 1) * t_inc;
        if (j <= 5 && j >= 1) {
          // cout << "vel j: " << j << endl;
        }
      }

      // 将后续的节点区间依次后移, 保证整条时间轴不冲突
      for (int j = i + p_ + 2; j < u_.rows(); ++j) {      /
        u_(j) += delta_t;
      }
    }
  }

  /* 检查加速度可达性 */
  for (int i = 0; i < P.rows() - 2; ++i) {

    Eigen::VectorXd acc = p_ * (p_ - 1) *
        ((P.row(i + 2) - P.row(i + 1)) / (u_(i + p_ + 2) - u_(i + 2)) -
         (P.row(i + 1) - P.row(i)) / (u_(i + p_ + 1) - u_(i + 1))) /
        (u_(i + p_ + 1) - u_(i + 2));

    if (fabs(acc(0)) > limit_acc_ + 1e-4 || fabs(acc(1)) > limit_acc_ + 1e-4 ||
        fabs(acc(2)) > limit_acc_ + 1e-4) {

      fea = false;
      if (show) cout << "[Realloc]: Infeasible acc " << i << " :" << acc.transpose() << endl;

      max_acc = -1.0;
      for (int j = 0; j < dimension; ++j) {
        max_acc = max(max_acc, fabs(acc(j)));
      }

      double ratio = sqrt(max_acc / limit_acc_) + 1e-4;
      if (ratio > limit_ratio_) ratio = limit_ratio_;
      // cout << "ratio: " << ratio << endl;

      double time_ori = u_(i + p_ + 1) - u_(i + 2);
      double time_new = ratio * time_ori;
      double delta_t  = time_new - time_ori;
      double t_inc    = delta_t / double(p_ - 1);

      if (i == 1 || i == 2) {
        // cout << "acc i: " << i << endl;
        for (int j = 2; j <= 5; ++j) {
          u_(j) += double(j - 1) * t_inc;
        }

        for (int j = 6; j < u_.rows(); ++j) {
          u_(j) += 4.0 * t_inc;
        }

      } else {

        for (int j = i + 3; j <= i + p_ + 1; ++j) {
          u_(j) += double(j - i - 2) * t_inc;
          if (j <= 5 && j >= 1) {
            // cout << "acc j: " << j << endl;
          }
        }

        for (int j = i + p_ + 2; j < u_.rows(); ++j) {
          u_(j) += delta_t;
        }
      }
    }
  }

  return fea;
}

/*
*@berif: TODO
*@param ratio:  系数
*@return:       void 
*/
void NonUniformBspline::lengthenTime(const double& ratio) {
  int num1 = 5;
  int num2 = getKnot().rows() - 1 - 5;

  double delta_t = (ratio - 1.0) * (u_(num2) - u_(num1));
  double t_inc   = delta_t / double(num2 - num1);
  for (int i = num1 + 1; i <= num2; ++i) u_(i) += double(i - num1) * t_inc;
  for (int i = num2 + 1; i < u_.rows(); ++i) u_(i) += delta_t;
}

void NonUniformBspline::recomputeInit() {}

/*
*@berif: 通过对前端hybrid A*寻找到的初始路径进行拟合得到的
*@param ts:                   
*@param point_set:            轨迹点
*@param start_end_derivative:
*@param ctrl_pts:             控制节点
*@return:
*/
void NonUniformBspline::parameterizeToBspline(const double& ts, const vector<Eigen::Vector3d>& point_set,
                                              const vector<Eigen::Vector3d>& start_end_derivative,
                                              Eigen::MatrixXd&               ctrl_pts) {
  if (ts <= 0) {
    cout << "[B-spline]:time step error." << endl;
    return;
  }

  if (point_set.size() < 2) {
    cout << "[B-spline]:point set have only " << point_set.size() << " points." << endl;
    return;
  }

  if (start_end_derivative.size() != 4) {
    cout << "[B-spline]:derivatives error." << endl;
  }

  int K = point_set.size();          // 轨迹点个数

  // 求解A矩阵
  Eigen::Vector3d prow(3), vrow(3), arow(3);
  prow << 1, 4, 1;
  vrow << -1, 0, 1;
  arow << 1, -2, 1;

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(K + 4, K + 2);

  for (int i = 0; i < K; ++i) A.block(i, i, 1, 3) = (1 / 6.0) * prow.transpose();

  A.block(K, 0, 1, 3)         = (1 / 2.0 / ts) * vrow.transpose();
  A.block(K + 1, K - 1, 1, 3) = (1 / 2.0 / ts) * vrow.transpose();
  A.block(K + 2, 0, 1, 3)     = (1 / ts / ts) * arow.transpose();
  A.block(K + 3, K - 1, 1, 3) = (1 / ts / ts) * arow.transpose();

  // A.block(0, 0, K, K + 2) = (1 / 6.0) * A.block(0, 0, K, K + 2);
  // A.block(K, 0, 2, K + 2) = (1 / 2.0 / ts) * A.block(K, 0, 2, K + 2);
  // A.row(K + 4) = (1 / ts / ts) * A.row(K + 4);
  // A.row(K + 5) = (1 / ts / ts) * A.row(K + 5);

  // 求解b
  Eigen::VectorXd bx(K + 4), by(K + 4), bz(K + 4);
  for (int i = 0; i < K; ++i) {
    bx(i) = point_set[i](0);
    by(i) = point_set[i](1);
    bz(i) = point_set[i](2);
  }

  for (int i = 0; i < 4; ++i) {
    bx(K + i) = start_end_derivative[i](0);
    by(K + i) = start_end_derivative[i](1);
    bz(K + i) = start_end_derivative[i](2);
  }

  // 求解 Ax = b
  Eigen::VectorXd px = A.colPivHouseholderQr().solve(bx);
  Eigen::VectorXd py = A.colPivHouseholderQr().solve(by);
  Eigen::VectorXd pz = A.colPivHouseholderQr().solve(bz);

  // 转换成控制节点
  ctrl_pts.resize(K + 2, 3);
  ctrl_pts.col(0) = px;
  ctrl_pts.col(1) = py;
  ctrl_pts.col(2) = pz;
}

/*
*@berif:  获得B样条轨迹的时间区间
*@param:  void
*@return: B样条轨迹的时间区间
*/
double NonUniformBspline::getTimeSum() {
  double tm, tmp;
  getTimeSpan(tm, tmp);
  return tmp - tm;
}

/*
 *@berif: 获得节点的路径长度
 *@param res: 0.01增量
 *@return: 路径长度
 */
double NonUniformBspline::getLength(const double& res) {
  double          length = 0.0;
  double          dur    = getTimeSum();
  Eigen::VectorXd p_l    = evaluateDeBoorT(0.0), p_n;
  for (double t = res; t <= dur + 1e-4; t += res) {
    p_n = evaluateDeBoorT(t);
    length += (p_n - p_l).norm();
    p_l = p_n;
  }
  return length;
}

/*
 *@berif: 获得jerk值(位置的三阶导)
 *@param:void
 *@return: jerk值
 */
double NonUniformBspline::getJerk() {
  // 求三阶导
  NonUniformBspline jerk_traj = getDerivative().getDerivative().getDerivative();

  Eigen::VectorXd times     = jerk_traj.getKnot();
  Eigen::MatrixXd ctrl_pts  = jerk_traj.getControlPoint();
  int             dimension = ctrl_pts.cols();

  double jerk = 0.0;
  for (int i = 0; i < ctrl_pts.rows(); ++i) {
    for (int j = 0; j < dimension; ++j) {
      jerk += (times(i + 1) - times(i)) * ctrl_pts(i, j) * ctrl_pts(i, j);
    }
  }

  return jerk;
}

/*
 *@berif: 获得速度的最大值和速度的平均值
 *@param mean_v: 速度平均值
 *@param max_x:  速度最大值
 *@return:       void
 */
void NonUniformBspline::getMeanAndMaxVel(double& mean_v, double& max_v) {
  NonUniformBspline vel = getDerivative();            // 获得节点速度
  double            tm, tmp;                          //
  vel.getTimeSpan(tm, tmp);

  double max_vel = -1.0, mean_vel = 0.0;              // 速度最大值、速度平均值
  int    num = 0;                                     // 数量

  // 求速度的和以及速度最大值
  for (double t = tm; t <= tmp; t += 0.01) {
    Eigen::VectorXd vxd = vel.evaluateDeBoor(t);      // 节点速度
    double          vn  = vxd.norm();                 // 节点速度大小

    mean_vel += vn;                                   // 累加求速度和
    ++num;                                            // 计算速度的个数
    if (vn > max_vel) {
      max_vel = vn;
    }
  }

  mean_vel = mean_vel / double(num);                  // 求解平均值
  mean_v   = mean_vel;                                // 速度平均值
  max_v    = max_vel;                                 // 速度最大值
}

/*
 *@berif: 获得加速度的最大值和速度的平均值
 *@param mean_v: 加速度平均值
 *@param max_x:  加速度最大值
 *@return:       void
 */
void NonUniformBspline::getMeanAndMaxAcc(double& mean_a, double& max_a) {
  NonUniformBspline acc = getDerivative().getDerivative();
  double            tm, tmp;
  acc.getTimeSpan(tm, tmp);

  double max_acc = -1.0, mean_acc = 0.0;
  int    num = 0;
  for (double t = tm; t <= tmp; t += 0.01) {
    Eigen::VectorXd axd = acc.evaluateDeBoor(t);
    double          an  = axd.norm();

    mean_acc += an;
    ++num;
    if (an > max_acc) {
      max_acc = an;
    }
  }

  mean_acc = mean_acc / double(num);
  mean_a   = mean_acc;
  max_a    = max_acc;
}
}  // namespace fast_planner
