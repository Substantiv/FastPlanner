"""
    B样条曲线法实现车辆轨迹规划
"""
import numpy as np
import matplotlib.pyplot as plt
import copy
from celluloid import Camera  # 保存动图时用，


def BaseFunction(i=None, k=None, u=None, NodeVector=None):
    """第 i个k阶B样条基函数

    Args:
        i (_type_, optional): _description_. Defaults to None.
        k (_type_, optional): B样条阶数k. Defaults to None.
        u (_type_, optional): 自变量. Defaults to None.
        NodeVector (_type_, optional): 节点向量. array([u0,u1,u2,...,u_n+k],shape=[1,n+k+1].

    Returns:
        _type_: _description_
    """
    if k == 1:  # 0次B样条（1阶B样条）
        if u >= NodeVector[0, i] and u < NodeVector[0, i + 1]:
            Bik_u = 1
        else:
            Bik_u = 0
    else:
        # 公式中的两个分母
        denominator_1 = NodeVector[0, i + k - 1] - NodeVector[0, i]
        denominator_2 = NodeVector[0, i + k] - NodeVector[0, i + 1]
        # 如果遇到分母为 0的情况：
        # 1. 如果此时分子也为0，约定这一项整体为0；
        # 2. 如果此时分子不为0，则约定分母为1 。
        if denominator_1 == 0:
            denominator_1 = 1
        if denominator_2 == 0:
            denominator_2 = 1
        Bik_u = (u - NodeVector[0, i]) / denominator_1 * BaseFunction(i, k - 1, u, NodeVector) + \
                (NodeVector[0, i + k] - u) / denominator_2 * BaseFunction(i + 1, k - 1, u, NodeVector)

    return Bik_u


def U_quasi_uniform(n=None, k=None):
    """准均匀B样条的节点向量计算
    首末值定义为 0 和 1
    Args:
        n (_type_, optional): 控制点个数-1，控制点共n+1个. Defaults to None.
        k (_type_, optional): B样条阶数k， k阶B样条，k-1次曲线. Defaults to None.

    Returns:
        _type_: _description_
    """
    # 准均匀B样条的节点向量计算，共n+1个控制顶点，k-1次B样条，k阶
    NodeVector = np.zeros((1, n + k + 1))
    piecewise = n - k + 2  # B样条曲线的段数:控制点个数-次数

    if piecewise == 1:  # 只有一段曲线时，n = k-1
        NodeVector[0, n + 1:n + k + 1] = 1
    else:
        for i in range(n - k + 1):  # 中间段内节点均匀分布：两端共2k个节点，中间还剩(n+k+1-2k=n-k+1）个节点
            NodeVector[0, k + i] = NodeVector[0, k + i - 1] + 1 / piecewise

        NodeVector[0, n + 1:n + k + 1] = 1  # 末尾重复度k

    return NodeVector


def U_piecewise_B_Spline(n=None, k=None):
    """分段B样条的节点向量计算
    首末值定义为 0 和 1
    # 分段Bezier曲线的节点向量计算，共n+1个控制顶点，k阶B样条，k-1次曲线
    # 分段Bezier端节点重复度为k，内间节点重复度为k-1,且满足n/(k-1)为正整数
    Args:
        n (_type_, optional): 控制点个数-1，控制点共n+1个. Defaults to None.
        k (_type_, optional): B样条阶数k， k阶B样条，k-1次曲线. Defaults to None.

    Returns:
        _type_: _description_
    """

    NodeVector = np.zeros((1, n + k + 1))
    if n % (k - 1) == 0 and k - 1 > 0:  # 满足n是k-1的整数倍且k-1为正整数
        NodeVector[0, n + 1:n + k + 1] = 1  # 末尾n+1到n+k的数重复
        piecewise = n / (k - 1)  # 设定内节点的值
        if piecewise > 1:
            # for i in range(k-1): # 内节点重复k-1次
            NodeVector[0, k:n + 1] = 1 / piecewise  # 内节点重复度k-1
    else:
        print('error!需要满足n是k-1的整数倍且k-1为正整数')

    return NodeVector


if __name__ == '__main__':
    ## 数据定义
    k = 5  # k阶、k-1次B样条

    flag = 3  # 1,2,3分别绘制均匀B样条曲线、准均匀B样条曲线,分段B样条

    d = 3.5  # # 道路标准宽度
    P = np.array([
        [0, -d / 2],
        [10, -d / 2],
        [25, -d / 2 + 0.5],
        [25, d / 2 - 0.5],
        [40, d / 2],
        [50, d / 2]
    ])

    n = len(P) - 1  # 控制点个数-1

    ## 生成B样条曲线

    path = []  # 路径点数据存储
    Bik_u = np.zeros((n + 1, 1))

    if flag == 1:  # 均匀B样条很简单
        NodeVector = np.array([np.linspace(0, 1, n + k + 1)]
                              )  # 均匀B样条节点向量，首末值定义为 0 和 1
        # for u in np.arange(0,1,0.001):
        # u的范围为[u_{k-1},u_{n+2}],这样才是open的曲线，不然你可以使用[0,1]试试。
        for u in np.arange((k - 1) / (n + k + 1), (n + 2) / (n + k + 1) + 0.001, 0.001):
            for i in range(n + 1):
                Bik_u[i, 0] = BaseFunction(i, k, u, NodeVector)
            p_u = P.T @ Bik_u
            path.append(p_u)
    elif flag == 2:
        NodeVector = U_quasi_uniform(n, k)
        for u in np.arange(0, 1, 0.005):
            for i in range(n + 1):
                Bik_u[i, 0] = BaseFunction(i, k, u, NodeVector)
            p_u = P.T @ Bik_u
            path.append(p_u)
    elif flag == 3:
        NodeVector = U_piecewise_B_Spline(n, k)
        for u in np.arange(0, 1, 0.005):
            for i in range(n + 1):
                Bik_u[i, 0] = BaseFunction(i, k, u, NodeVector)
            p_u = P.T @ Bik_u
            path.append(p_u)
    path = np.array(path)

    ## 画图
    fig = plt.figure(1)
    # plt.ylim(-4, 4)
    # plt.axis([-10, 100, -15, 15])
    camera = Camera(fig)
    len_line = 50
    # 画灰色路面图
    GreyZone = np.array([[- 5, - d - 0.5], [- 5, d + 0.5],
                         [len_line, d + 0.5], [len_line, - d - 0.5]])
    for i in range(len(path)):
        # plt.cla()

        plt.fill(GreyZone[:, 0], GreyZone[:, 1], 'gray')
        # 画分界线
        plt.plot(np.array([- 5, len_line]), np.array([0, 0]), 'w--')

        plt.plot(np.array([- 5, len_line]), np.array([d, d]), 'w')

        plt.plot(np.array([- 5, len_line]), np.array([- d, - d]), 'w')

        plt.plot(P[:, 0], P[:, 1], 'ro')
        plt.plot(P[:, 0], P[:, 1], 'y')
        # 设置坐标轴显示范围
        # plt.axis('equal')
        plt.gca().set_aspect('equal')
        # 绘制路径

        plt.plot(path[0:i, 0], path[0:i, 1], 'g')  # 路径点
        plt.pause(0.001)
    #     camera.snap()
    # animation = camera.animate()
    # animation.save('trajectory.gif')

