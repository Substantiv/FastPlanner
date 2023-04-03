import sys
import pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent.parent))

import math
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np

from utils.angle import rot_mat_2d

# Covariance for EKF simulation
# EKF仿真的协方差
Q = np.diag([
    0.1,  # variance of location on x-axis                                      # 在x轴定位的方差
    0.1,  # variance of location on y-axis                                      # 在y轴定位的方差
    np.deg2rad(1.0),  # variance of yaw angle                                   # yaw的方差
    1.0  # variance of velocity                                                 # 速度的方差
]) ** 2  # predict state covariance                                             # 预测状态协方差
R = np.diag([1.0, 1.0]) ** 2  # Observation x,y position covariance             # 观测 x,y 位置的协方差

#  Simulation parameter 仿真参数
INPUT_NOISE = np.diag([1.0, np.deg2rad(30.0)]) ** 2                             # 输入噪声
GPS_NOISE = np.diag([0.5, 0.5]) ** 2                                            # GPS噪声

DT = 0.1  # time tick [s]                                                       # 时间间隔 [s]
SIM_TIME = 50.0  # simulation time [s]                                          # 仿真时间 [s]

show_animation = True                                                           # 展示动画


def calc_input():
    """
    Returns:[速度, 角速度]
    """
    v = 1.0  # [m/s]
    yawrate = 0.1  # [rad/s]
    u = np.array([[v], [yawrate]])
    return u


def observation(xTrue, xd, u):
    """
    Returns:
        xTrue:真实的状态变量,未加入噪声的状态变量(用于做对比)
        z:加入噪声后的GPS测量值
        xd:加入噪声后的状态变量
        ud:加入噪声后的IMU测量值
    """
    # 卡尔曼滤波器第一个公式
    xTrue = motion_model(xTrue, u)

    # add noise to gps x-y 给观测值增加噪声(GPS的噪声)
    z = observation_model(xTrue) + GPS_NOISE @ np.random.randn(2, 1)

    # add noise to input 给输入增加噪声(IMU的噪声)
    ud = u + INPUT_NOISE @ np.random.randn(2, 1)

    # 加入噪声后的观测
    xd = motion_model(xd, ud)

    return xTrue, z, xd, ud


def motion_model(x, u):
    """
        状态方程
    """
    F = np.array([[1.0, 0, 0, 0],
                  [0, 1.0, 0, 0],
                  [0, 0, 1.0, 0],
                  [0, 0, 0, 0]])

    B = np.array([[DT * math.cos(x[2, 0]), 0],
                  [DT * math.sin(x[2, 0]), 0],
                  [0.0, DT],
                  [1.0, 0.0]])

    x = F @ x + B @ u

    return x


def observation_model(x):
    """
        观测方程
    """
    H = np.array([
        [1, 0, 0, 0],
        [0, 1, 0, 0]
    ])

    z = H @ x

    return z


def jacob_f(x, u):
    """
    Jacobian of Motion Model

    motion model
    x_{t+1} = x_t+v*dt*cos(yaw)
    y_{t+1} = y_t+v*dt*sin(yaw)
    yaw_{t+1} = yaw_t+omega*dt
    v_{t+1} = v{t}
    so
    dx/dyaw = -v*dt*sin(yaw)
    dx/dv = dt*cos(yaw)
    dy/dyaw = v*dt*cos(yaw)
    dy/dv = dt*sin(yaw)
    """
    yaw = x[2, 0]
    v = u[0, 0]
    jF = np.array([
        [1.0, 0.0, -DT * v * math.sin(yaw), DT * math.cos(yaw)],
        [0.0, 1.0, DT * v * math.cos(yaw), DT * math.sin(yaw)],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0]])

    return jF


def jacob_h():
    # Jacobian of Observation Model
    jH = np.array([
        [1, 0, 0, 0],
        [0, 1, 0, 0]
    ])

    return jH


def ekf_estimation(xEst, PEst, z, u):
    """
    Return:
        xEst:更新后的状态变量
        PEst:更新后的协方差
    """
    #  Predict
    xPred = motion_model(xEst, u)                   # 预测计算状态变量的先验值
    jF = jacob_f(xEst, u)                           # 计算F矩阵的雅可比矩阵
    PPred = jF @ PEst @ jF.T + Q                    # 计算先验协方差

    #  Update
    jH = jacob_h()                                  # 计算H矩阵的雅可比矩阵
    zPred = observation_model(xPred)                
    y = z - zPred                                   # 校正部分的后验公式的一部分
    S = jH @ PPred @ jH.T + R                       # 用于计算卡尔曼增益公式的一部分
    K = PPred @ jH.T @ np.linalg.inv(S)             # 用于计算卡尔曼增益
    xEst = xPred + K @ y                            # 校正部分的后验公式
    PEst = (np.eye(len(xEst)) - K @ jH) @ PPred     # 校正部分的后验协方差
    return xEst, PEst


def plot_covariance_ellipse(xEst, PEst):  # pragma: no cover
    Pxy = PEst[0:2, 0:2]
    eigval, eigvec = np.linalg.eig(Pxy)

    if eigval[0] >= eigval[1]:
        bigind = 0
        smallind = 1
    else:
        bigind = 1
        smallind = 0

    t = np.arange(0, 2 * math.pi + 0.1, 0.1)
    a = math.sqrt(eigval[bigind])
    b = math.sqrt(eigval[smallind])
    x = [a * math.cos(it) for it in t]
    y = [b * math.sin(it) for it in t]
    angle = math.atan2(eigvec[1, bigind], eigvec[0, bigind])
    fx = rot_mat_2d(angle) @ (np.array([x, y]))
    px = np.array(fx[0, :] + xEst[0, 0]).flatten()
    py = np.array(fx[1, :] + xEst[1, 0]).flatten()
    plt.plot(px, py, "--r")


def main():
    print(__file__ + " start!!")

    time = 0.0

    # State Vector [x y yaw v]'
    xEst = np.zeros((4, 1))
    xTrue = np.zeros((4, 1))
    PEst = np.eye(4)

    xDR = np.zeros((4, 1))  # Dead reckoning

    # history
    hxEst = xEst
    hxTrue = xTrue
    hxDR = xTrue
    hz = np.zeros((2, 1))

    while SIM_TIME >= time:
        time += DT
        u = calc_input()

        xTrue, z, xDR, ud = observation(xTrue, xDR, u)

        xEst, PEst = ekf_estimation(xEst, PEst, z, ud)

        # store data history
        hxEst = np.hstack((hxEst, xEst))
        hxDR = np.hstack((hxDR, xDR))
        hxTrue = np.hstack((hxTrue, xTrue))
        hz = np.hstack((hz, z))

        if show_animation:
            plt.cla()
            # for stopping simulation with the esc key.
            plt.gcf().canvas.mpl_connect('key_release_event',
                    lambda event: [exit(0) if event.key == 'escape' else None])
            plt.plot(hz[0, :], hz[1, :], ".g")
            plt.plot(hxTrue[0, :].flatten(),
                     hxTrue[1, :].flatten(), "-b")
            plt.plot(hxDR[0, :].flatten(),
                     hxDR[1, :].flatten(), "-k")
            plt.plot(hxEst[0, :].flatten(),
                     hxEst[1, :].flatten(), "-r")
            plot_covariance_ellipse(xEst, PEst)
            plt.axis("equal")
            plt.grid(True)
            plt.pause(0.001)


if __name__ == '__main__':
    main()
