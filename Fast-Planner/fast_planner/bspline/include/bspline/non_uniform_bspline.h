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



#ifndef _NON_UNIFORM_BSPLINE_H_
#define _NON_UNIFORM_BSPLINE_H_

#include <Eigen/Eigen>
#include <algorithm>
#include <iostream>

using namespace std;

namespace fast_planner {

/*
 *@berif:不同维度的非均匀b样条的实现.由于均匀B样条是非均匀B样条的特殊情况,所以也可以表示均匀B样条
 */
class NonUniformBspline {
private:
  // 不同维度B样条曲线的控制点,每一行表示一个控制点,每一列表示一个维度(在三维空间的N个控制点Nx3)
  Eigen::MatrixXd control_points_;

  int             p_, n_, m_;  // p 阶数, n+1 控制点, m = n+p+1
  Eigen::VectorXd u_;          // 节点向量
  double          interval_;   // 节点间隔 \delta t

  Eigen::MatrixXd getDerivativeControlPoints();

  double limit_vel_, limit_acc_, limit_ratio_;  // 速度限制、加速度限制和时间调整比例限制

public:
  NonUniformBspline() {}
  NonUniformBspline(const Eigen::MatrixXd& points, const int& order, const double& interval);
  ~NonUniformBspline();

  // initialize as an uniform B-spline 使用标准B样条进行初始化
  void setUniformBspline(const Eigen::MatrixXd& points, const int& order, const double& interval);

  // get / set basic bspline info

  void                                   setKnot(const Eigen::VectorXd& knot);
  Eigen::VectorXd                        getKnot();
  Eigen::MatrixXd                        getControlPoint();
  double                                 getInterval();
  void                                   getTimeSpan(double& um, double& um_p);
  pair<Eigen::VectorXd, Eigen::VectorXd> getHeadTailPts();

  // compute position / derivative

  Eigen::VectorXd   evaluateDeBoor(const double& u);   // use u \in [up, u_mp]
  Eigen::VectorXd   evaluateDeBoorT(const double& t);  // use t \in [0, duration]
  NonUniformBspline getDerivative();

  // 3D B-spline interpolation of points in point_set, with boundary vel&acc
  // constraints
  // input : (K+2) points with boundary vel/acc; ts
  // output: (K+6) control_pts
  static void parameterizeToBspline(const double& ts, const vector<Eigen::Vector3d>& point_set,
                                    const vector<Eigen::Vector3d>& start_end_derivative,
                                    Eigen::MatrixXd&               ctrl_pts);

  /* check feasibility, adjust time */

  void   setPhysicalLimits(const double& vel, const double& acc);
  bool   checkFeasibility(bool show = false);
  double checkRatio();
  void   lengthenTime(const double& ratio);
  bool   reallocateTime(bool show = false);

  /* for performance evaluation */

  double getTimeSum();
  double getLength(const double& res = 0.01);
  double getJerk();
  void   getMeanAndMaxVel(double& mean_v, double& max_v);
  void   getMeanAndMaxAcc(double& mean_a, double& max_a);

  void recomputeInit();

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
}  // namespace fast_planner
#endif