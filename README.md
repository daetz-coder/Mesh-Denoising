网格去噪

> python可视化算法流程：[Mesh-Denoising/Non-Iterative, Feature-Preserving Mesh Smoothing.ipynb at main · daetz-coder/Mesh-Denoising (github.com)](https://github.com/daetz-coder/Mesh-Denoising/blob/main/Non-Iterative%2C Feature-Preserving Mesh Smoothing.ipynb)
>
> C++实现：[daetz-coder/Mesh-Denoising: Method for implementing grid denoising in “Non-Iterative, Feature-Preserving Mesh Smoothing“ (github.com)](https://github.com/daetz-coder/Mesh-Denoising)
>
> 相关的obj文件[Release v0.1 · daetz-coder/Mesh-Denoising (github.com)](https://github.com/daetz-coder/Mesh-Denoising/releases/tag/add-obj-files)

## 简介

本项目基于《Non-Iterative, Feature-Preserving Mesh Smoothing》论文，实现了一种非迭代、特征保留的网格去噪算法。项目通过Python环境下的Jupyter Notebook演示了算法原理与流程，包括利用高斯函数作为空间权重和影响权重，对目标顶点在局部三角面片上的投影进行加权平均，从而平滑噪声同时保留边缘特征。在C++实现部分，借助OpenMesh库完成了OBJ文件的读取、处理与保存，通过高斯权重计算和顶点投影等步骤对整个网格进行平滑处理，并记录每个顶点的移动日志。项目涵盖了算法演示、代码实现、实验环境配置和使用说明。

## 一、论文介绍

+ **[Non-Iterative, Feature-Preserving Mesh Smoothing](https://dl.acm.org/doi/abs/10.1145/1201775.882367)**

![image-20241224210505920](https://daetz-image.oss-cn-hangzhou.aliyuncs.com/img/202412242105005.png)

*We employ a spatial weight $f$ that depends on the distance $\|p - c_q  \|$between point $p$ and the centroid $c_q$ of region $q$. We also use an influence weight $g$ that depends on the distance $ \Pi_q(p) - p \|$between the prediction and the original position of $p$. Finally, we weight by the area $a_q$ of the triangles to account for variations in the sampling rate of the surface. The estimate $p'$ for a point on surface $S$ is then:*

$$
p' = \frac{1}{k(p)} \sum_{q \in S} \Pi_q(p) \, a_q \, f\left( \| c_q - p \| \right) \, g\left( \| \Pi_q(p) - p \| \right)
$$

*where*

$$
k(p) = \sum_{q \in S} a_q \, f\left( \| c_q - p \| \right) \, g\left( \| \Pi_q(p) - p \| \right)
$$

我们采用的空间权重 $f$ 取决于 与 $q$的质心 $c_q$ 之间的距离 $|| p-c_q||$。我们还使用影响权重 $g$ ，该影响权重取决于预测与$p$的原始位置之间的距离 $||\Pi_q(p)−p∥$ 。最后，我们用三角形的面积 $a_q$ 来衡量表面采样率的变化。表面$S$上一点的估计 $p'$是

论文中的核心公式是（文中式(3)）：

$$
p' = \frac{1}{k(p)} \sum_{q \in S} \Pi_q(p) \, a_q \, f\left( \| c_q - p \| \right) \, g\left( \| \Pi_q(p) - p \| \right)
$$

其中

$$
k(p) = \sum_{q \in S} a_q \, f\left( \| c_q - p \| \right) \, g\left( \| \Pi_q(p) - p \| \right)
$$

*Gaussians are used both for the spatial weight f and for the in uence weight g inthis paper. Other robust in uence weightscould also beused, but Gaussians haveperformed wellinour experiments, aswellasthe work ofothers [Smith and Brady 1997; Tomasi and Manduchi1998; Durand and Dorsey 2002].*

*本文中的空间权重 f 和影响权重 g 都使用高斯函数。也可以使用其他强大的影响权重，但高斯模型在我们的实验以及其他人的工作中表现良好*

作者提到本文中的空间权重 $f$ 和影响权重 $g $都使用**高斯函数**，所以我们在后续也使用高斯函数

+ $p$ 表示当前顶点 (需要被去噪的那个顶点)
+ $p'$ 表示该顶点去噪后（或平滑后）的新位置
+ $q\in S$表示网格上的所有三角面片(或者在实际实现时，往往取“局部邻域”即可)
+ $\Pi_{q}(p)$表示将顶点 $p$ 投影到三角面片 $q$ 上得到的“预测点”(prediction)
+ $a_q$ 表示三角面片 $q$ 的面积，用来补偿采样不均匀(如果一些区域三角面片更密集，就相当于对该区域要稍加“削弱”或“平衡”)
+ $c_q$ 表示三角面片 $q$ 的质心(即该三角形三个顶点的平均位置)
+ $f(\| c_q - p \|)$称为**空间权重(spatial weight)**，主要控制“离得远的三角形面片”对当前顶点的贡献要小
+ $g(\|\Pi_{q}(p) - p\|)$ 称为**影响权重(influence weight)**，主要控制“预测点离原始位置很远时，该预测应被视为异常/outlier，贡献要小”
+ $k(p)$ 是所有权重之和，用来归一化，保证计算出来的 $p'$ 是某种加权平均结果。

### 1、顶点移动方向

在每次滤波时，该公式会对顶点 $p$ 进行一次更新，计算出一个新的位置 $p'$。

直观来看，$p'$ 相当于所有面片给出的预测点 $\Pi_{q}$ 的加权平均，而每个预测点根据“离 $p$ 的空间距离(由 $f$ 控制)”以及“预测点和原始位置差异(由 $g$ 控制)”来分配权重。

这样就意味着：如果某个三角面片 $q$ 对当前顶点做出的预测和 $p$ 差异过大，可能表示在尖锐特征处跨过了棱线，那么 $g$ 就会让这个预测权重变小，相当于**自动排斥跨越尖锐特征的“错误”平滑**。

整个过程并不是“朝法向方向”或“朝平面”的简单移动，而是“由多个邻居面片对顶点进行投影，然后再综合加权”得到一个新的位置。最终效果是：对噪声点进行平滑的同时又能较好地保留边缘或特征。

### 2、除三角外不需连通

*Filtering amesh involvesevaluating Equation (3) for every vertex and then moving them asagroup totheir estimated positions. Note that no connectivityisrequired beyond triangles: wesimply use the Euclidean distance tothe centroid ofsurrounding triangles to nd the spatial neighborhoodof avertex. Awider spatial  lter includes alarger number ofneigh bors inthe estimate, and can therefore removea greater amountofnoise, orsmo oth larger features. The in uence weightdetermines when the predictions ofneigh bors are considered outliers (byaccording them less weight), and therebycon trols the size offeatures that are preservedinthe  ltering.*

过滤网格涉及对每个顶点评估公式（3），然后将它们作为一个组移动到其估计位置。请注意，除了三角形之外，不需要连通性：我们只需使用到周围三角形质心的欧几里得距离来查找顶点的空间邻域。更宽的空间滤波器在估计中包含更多数量的邻居，因此可以去除更大量的噪声或平滑其他更大的特征。影响权重决定了邻居的预测何时被视为异常值（通过赋予它们较小的权重），从而控制过滤中保留的特征的大小。

这是指在实际操作中，我们只需要知道当前顶点 $p$的“**空间邻居**”可以由与之相邻的三角面片(或其质心)来确定。也就是说：

+ **只关心顶点周围出现的三角面**来构造空间邻域(找出周围面片的质心距离较近者)，不需要再去关心更复杂的网格拓扑(比如网格上更高级别的区域分段之类)。
+ 所以只要能找到每个三角面片的几何中心 $c_q$，并且能把顶点 $p$ 投影到该三角面片上 $\Pi_q(p)$，就可以进行权重计算并做加权平均。
+ “不需要额外连通”可以理解为，我们不需要一个全局的网格数据结构来进行大规模的遍历，只需要基于局部三角面的几何信息就够了。换句话说，“只要知道局部的三角形信息和它们之间的空间关系(比如质心、面积等)”，就能完成一轮更新。

## 二、代码介绍

### 1、初始化$p$和$q$

+ $p$是需要处理的点
+ $q$泛指周围的三角
+ 我们在这里仅考虑周围的三个网格三角($q$) 对$p$进行网格去噪

```python
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
# 定义原始顶点 p
p = np.array([1, 2, 1])
print(f"Original vertex p: {p}\n")

# 定义三角面片 q1, q2, q3，每个三角形由三个顶点组成
triangles = {
    'q1': {
        'vertices': np.array([[0, 0, 3],
                              [3, 0, 0],
                              [0, 3, 0]]),
    },
    'q2': {
        'vertices': np.array([[2, 1, 1],
                              [4, 1, 0],
                              [3, 2, 2]]),
    },
    'q3': {
        'vertices': np.array([[1, 4, 1],
                              [3, 5, 1],
                              [2, 2, 3]]),
    }
}
```

![image-20241224195413357](https://daetz-image.oss-cn-hangzhou.aliyuncs.com/img/202412241954416.png)

### 2、计算$a_q$和$c_q$

+ $a_q$表示各个网格三角的面积
+ $c_q$表示网格三角的质心
+ **质心计算公式**

  对于一个二维平面上的三角形，质心是该三角形三个顶点的坐标的**算术平均值**。设三角形的三个顶点坐标为 $(x_1, y_1), (x_2, y_2), (x_3, y_3)$那么三角形的质心 $(x_c,y_c)$ 的坐标可以通过以下公式计算：

  $x_c = \frac{x_1 + x_2 + x_3}{3}$

  $y_c = \frac{y_1 + y_2 + y_3}{3}$

```python
# 计算每个三角形的平面方程、面积和质心
for key, tri in triangles.items():
    v1, v2, v3 = tri['vertices']
  
    # 计算平面法向量
    normal = np.cross(v2 - v1, v3 - v1)
    norm_length = np.linalg.norm(normal)
    normal_unit = normal / norm_length
    tri['normal'] = normal_unit
  
    # 计算平面方程常数项 c (n·v = c)
    c = np.dot(normal_unit, v1)
    tri['plane_eq'] = f"{normal_unit[0]:.2f}x + {normal_unit[1]:.2f}y + {normal_unit[2]:.2f}z = {c:.2f}"
  
    # 计算面积 (1/2 * |n|)
    area = 0.5 * norm_length
    tri['area'] = area
  
    # 计算质心
    centroid = (v1 + v2 + v3) / 3
    tri['centroid'] = centroid
  
    print(f"Triangle {key}:")
    print(f"  Vertices:\n{tri['vertices']}")
    print(f"  Plane equation: {tri['plane_eq']}")
    print(f"  Area: {tri['area']:.4f}")
    print(f"  Centroid: {tri['centroid']}\n")
  
```

```less
  Plane equation: 0.58x + 0.58y + 0.58z = 1.73
  Area: 7.7942
  Centroid: [1. 1. 1.]
  Plane equation: 0.27x + -0.80y + 0.53z = 0.27
  Area: 1.8708
  Centroid: [3.         1.33333333 1.        ]
  Plane equation: 0.30x + -0.60y + -0.75z = -2.83
  Area: 3.3541
  Centroid: [2.         3.66666667 1.66666667]
```

### 3、计算$p$的$\Pi_q(p)$

+ $\Pi_q(p)$ 用于表示原点$p$在各个面$q$上的投影点

```python
def project_point_to_plane(p, normal_unit, c):
    """
    投影点 p 到平面上。
    plane_eq: n·v = c
    normal_unit: 平面的单位法向量 n
    c: 平面方程常数项
    """
    distance = np.dot(normal_unit, p) - c
    projection = p - distance * normal_unit
    return projection
# 修正投影点计算
projections = {}
print("\nCalculating projections Πq(p):")
for key, tri in triangles.items():
    # 从平面方程字符串中提取常数项 c 并转换为浮点数
    c_str = tri['plane_eq'].split('=')[1].strip()
    try:
        c = float(c_str)
    except ValueError:
        print(f"Error: Cannot convert '{c_str}' to float for triangle {key}.")
        continue
    proj = project_point_to_plane(p, tri['normal'], c)
    projections[key] = proj
    print(f"Projection Π{key}(p): {proj}")

# 绘制图像
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 绘制三角面片
for key, tri in triangles.items():
    v1, v2, v3 = tri['vertices']
    verts = [v1, v2, v3]
    ax.add_collection3d(Poly3DCollection([verts], color='cyan', linewidths=1, edgecolors='r', alpha=0.5))

# 绘制原始顶点 p
ax.scatter(p[0], p[1], p[2], color='b', label='Original vertex p', s=100)

# 绘制投影点
for key, proj in projections.items():
    ax.scatter(proj[0], proj[1], proj[2], color='g', label=f'Projection Π{key}(p)', s=100)

# 设置坐标轴标签
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# 设置图像的显示范围
ax.set_xlim([0, 5])
ax.set_ylim([0, 5])
ax.set_zlim([0, 5])

# 显示图例
ax.legend()

# 显示图像
plt.show()
```

![image-20241224195857796](https://daetz-image.oss-cn-hangzhou.aliyuncs.com/img/202412241958860.png)

```less
Calculating projections Πq(p):
Projection Πq1(p): [0.66548263 1.66548263 0.66548263]
Projection Πq2(p): [1.28644625 1.14066125 1.5728925 ]
Projection Πq3(p): [0.64514591 2.70970819 1.88713524]
```

### 4、计算空间权重$ f(r)$

+ $ f(r)= f(\| c_q - p \|)$ 表示的是质心$c_q$ 与原点$q$的距离在高斯函数下得到权重

```python
# 定义高斯核函数
def gaussian_weight(r, sigma):
    return np.exp(- (r ** 2) / (2 * sigma ** 2))
# 参数
sigma_f = 1

# 计算空间权重 f(r)
f_weights = {}
print("\nCalculating spatial weights f(r):\n")
for key, tri in triangles.items():
    r = np.linalg.norm(tri['centroid'] - p)
    f = gaussian_weight(r, sigma_f)
    f_weights[key] = f
    print(f"f(||c{key} - p|| = {r:.3f}) = {f:.4f}")

```

```less
Calculating spatial weights f(r):

f(||cq1 - p|| = 1.000) = 0.6065
f(||cq2 - p|| = 2.108) = 0.1084
f(||cq3 - p|| = 2.055) = 0.1211
```

![image-20241224200106529](https://daetz-image.oss-cn-hangzhou.aliyuncs.com/img/202412242001592.png)

### 5、计算影响权重$g(Δ)$

+ $g(Δ)= g(\|\Pi_{q}(p) - p\|)$ 表示的是原点$p$对于$\Pi_{q}(p)$的距离使用高斯函数得到权重
+ 距离越近，权重越高，对$p'$的贡献越大
+ 当出现预测点离原始位置很远时，该预测应被视为异常/outlier，贡献要小

```bash
# 定义影响权重参数
sigma_g = 1
# 计算影响权重 g(Δ)
g_weights = {}
print("\nCalculating influence weights g(Δ):\n")
for key, proj in projections.items():
    delta = np.linalg.norm(proj - p)
    g = gaussian_weight(delta, sigma_g)
    g_weights[key] = g
    print(f"g(||Π{key}(p) - p|| = {delta:.3f}) = {g:.4f}")

```

```less
Calculating influence weights g(Δ):

g(||Πq1(p) - p|| = 0.579) = 0.8455
g(||Πq2(p) - p|| = 1.072) = 0.5631
g(||Πq3(p) - p|| = 1.190) = 0.4925
```

![image-20241224200216221](https://daetz-image.oss-cn-hangzhou.aliyuncs.com/img/202412242002295.png)

```python
import numpy as np
import matplotlib.pyplot as plt

def gaussian_weight(r, sigma):
    return np.exp(- (r ** 2) / (2 * sigma ** 2))

# Generate values of r from 0 to 10
r = np.linspace(0, 10, 100)

# Different values of sigma
sigmas = [0.5, 1, 1.5, 2]

# Plot Gaussian weights for each sigma
plt.figure(figsize=(8, 5))
for sigma in sigmas:
    plt.plot(r, gaussian_weight(r, sigma), label=f'sigma={sigma}')

plt.title('Gaussian Weight as a Function of Distance')
plt.xlabel('Distance (r)')
plt.ylabel('Weight')
plt.legend()
plt.grid(True)
plt.show()

```

![](https://daetz-image.oss-cn-hangzhou.aliyuncs.com/img/202412242050954.png)

### 6、计算$k(p)$

+ $k(p) = \sum_{q \in S} a_q \, f\left( \| c_q - p \| \right) \, g\left( \| \Pi_q(p) - p \| \right)$
+ $k(p_1)=a_{q_1} \times f(r_1) \times g(Δ_1)$
+ $k(p)=k(p_1)+k(p_2)+k(p_3)$
+ 也就表示三个平面的权重之和，用来归一化，保证计算出来的 $p'$ 是加权平均结果。

```python
# 计算综合权重 wq 和 k(p)
wq = {}
k_p = 0
print("\nCalculating combined weights wq and k(p):\n")
for key in triangles.keys():
    w = triangles[key]['area'] * f_weights[key] * g_weights[key]
    wq[key] = w
    k_p += w
    print(f"w{key} = a{key} * f(r) * g(Δ) = {triangles[key]['area']} * {f_weights[key]:.4f} * {g_weights[key]:.4f} ≈ {w:.4f}")

print(f"\nk(p) = wq1 + wq2 + wq3 ≈ {k_p:.4f}\n")

```

```less

Calculating combined weights wq and k(p):

wq1 = aq1 * f(r) * g(Δ) = 7.794228634059948 * 0.6065 * 0.8455 ≈ 3.9969
wq2 = aq2 * f(r) * g(Δ) = 1.8708286933869707 * 0.1084 * 0.5631 ≈ 0.1142
wq3 = aq3 * f(r) * g(Δ) = 3.3541019662496847 * 0.1211 * 0.4925 ≈ 0.2000

k(p) = wq1 + wq2 + wq3 ≈ 4.3111


```

### 7、计算$p'$

+ $p' = \frac{1}{k(p)} \sum_{q \in S} \Pi_q(p) \, a_q \, f\left( \| c_q - p \| \right) \, g\left( \| \Pi_q(p) - p \| \right)$区别与$k(p) = \sum_{q \in S} a_q \, f\left( \| c_q - p \| \right) \, g\left( \| \Pi_q(p) - p \| \right)$在于多了一个$\Pi_q(p)$也就是投影点

```bash
# 计算加权后的坐标
numerator = np.zeros(3)
print("Calculating weighted sum for new position p':\n")
for key in triangles.keys():
    weighted_proj = wq[key] * projections[key]
    numerator += weighted_proj
    print(f"w{key} * Π{key}(p) = {wq[key]:.4f} * {projections[key]} = {weighted_proj}")

p_prime = numerator / k_p
print(f"\nNew vertex position p' = {p_prime}\n")

```

```less
Calculating weighted sum for new position p':

wq1 * Πq1(p) = 3.9969 * [0.66666667 1.66666667 0.66666667] = [2.66463111 6.66157778 2.66463111]
wq2 * Πq2(p) = 0.1142 * [1.28571429 1.14285714 1.57142857] = [0.1467703  0.13046249 0.17938592]
wq3 * Πq3(p) = 0.2000 * [0.64444444 2.71111111 1.88888889] = [0.12891503 0.54233221 0.37785441]

New vertex position p' = [0.68202729 1.70125981 0.74733597]
```

```python
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# 绘制图像
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 绘制三角面片
for key, tri in triangles.items():
    verts = [list(tri['vertices'])]
    tri_collection = Poly3DCollection(verts, facecolors='cyan', linewidths=1, edgecolors='r', alpha=0.5)
    ax.add_collection3d(tri_collection)


# Plot original and new vertices
ax.scatter(*p, color='blue', s=100, label='Original vertex p')
ax.scatter(*p_prime, color='red', s=100, label='New vertex position p\'')

# Draw an arrow from p to p_prime using quiver
ax.quiver(p[0], p[1], p[2], p_prime[0]-p[0], p_prime[1]-p[1], p_prime[2]-p[2], color='green', length=1.0, arrow_length_ratio=0.3, normalize=True)

# Set labels and legend
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()

# Set viewing limits
ax.set_xlim([0, 5])
ax.set_ylim([0, 5])
ax.set_zlim([0, 5])

# Show plot
plt.show()
```

+ $p= [1, 2, 1]$
+ $ p' = [0.68202729,1.70125981,0.74733597]$
+ 其中蓝色代表原始点，红色表示经过网格去噪后移动的点，绿色箭头表示移动的方向，由于左下角这一块相较于原点的距离较近，故贡献较大，从而影响了移动方向

![image-20241224200505740](https://daetz-image.oss-cn-hangzhou.aliyuncs.com/img/202412242005816.png)

## 三、实验报告

写清楚实验内容、实验环境和实验方法，以及实验结果，可链接图像或录屏（不超过20M）

### 1、实验内容

本实验旨在重现论文《Non-Iterative, Feature-Preserving Mesh Smoothing》中提出的非迭代、特征保留的网格去噪算法。实验首先在Python环境下使用NumPy、SciPy、Trimesh和Matplotlib等库进行网格数据的读取与预处理，包括计算法向量和归一化处理，通过计算每个顶点的新位置，利用高斯函数作为空间权重和影响权重，对邻近三角面片的预测点进行加权平均，从而实现平滑去噪的同时保留网格的边缘特征。实验过程中设置了高斯函数的标准差参数σ_f和σ_g，以控制权重分配和特征保留的平衡。最终，通过代码实现顶点位置的更新，保存去噪后的网格并进行可视化对比，验证算法的有效性。整个实验流程模块化设计，确保了算法的可复现性和适应性，为进一步优化和应用提供了基础。并基于上述原理使用c++进行重现，展示去噪后的结果。

### 2、实验环境

```less
# Ubuntu 20.04.3 LTS
# Python 
Python 3.8.10
# requirements.txt
numpy
scipy
trimesh
matplotlib

# C++
C++17
CMake >= 3.10
GCC >= 7.0
OpenMesh
```

### 3、实验方法

#### 1)、main.cpp

```c++
#include "denoise_obj.hpp"           // 引入头文件 denoise_obj.hpp，该文件中声明了 denoise_obj 函数
#include <iostream>                  // 引入标准输入输出头文件，用于 std::cerr、std::cout 等

int main(int argc, char* argv[]) {   // main 函数入口，argc 表示命令行参数数量，argv 为参数列表
    if(argc != 3) {                  // 如果参数数量不等于 3（包括程序名在内），则说明参数不正确
        std::cerr << "Usage: denoise_obj <input.obj> <output.obj>\n"; // 向标准错误流输出提示信息
        return 1;                    // 返回非 0 值表示程序异常退出
    }

    std::string input_obj = argv[1]; // 将第二个参数（第一个用户参数）视为输入 OBJ 文件的路径
    std::string output_obj = argv[2];// 将第三个参数（第二个用户参数）视为输出 OBJ 文件的路径

    denoise_obj(input_obj, output_obj); // 调用去噪函数，传入输入和输出文件路径

    return 0;                       // 返回 0 表示程序成功执行
}

```

#### 2)、denoise_obj.hpp

```hpp
#ifndef DENOISE_OBJ_HPP             // 预处理指令：防止头文件被重复包含的保护宏（开始）
#define DENOISE_OBJ_HPP

#include <OpenMesh/Core/IO/MeshIO.hh>          // OpenMesh 库的 I/O 头文件，用于读写网格
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh> // OpenMesh 提供的三角网格类模版
#include <iostream>                            // 标准输入输出库
#include <vector>                              // std::vector 容器
#include <cmath>                               // 数学函数，如 std::sqrt、std::exp 等
#include <fstream>                             // 文件读写流
#include <sstream>                             // 字符串流
#include <ctime>                               // 时间相关操作，用于记录日志时间戳等

// 使用 OpenMesh 内置的三角网格类
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh; // 定义一个名为 MyMesh 的别名，表示三角网格类型

// 定义一个 3D 向量结构体
struct Vec3 {
    float x, y, z;                  // 三个浮点数，分别表示向量在 x、y、z 方向的分量

    Vec3();                         // 默认构造函数
    Vec3(float x, float y, float z);// 带参构造函数

    Vec3 operator+(const Vec3& other) const;   // 向量加法
    Vec3 operator-(const Vec3& other) const;   // 向量减法
    Vec3 operator*(float scalar) const;        // 向量与标量相乘
  
    float dot(const Vec3& other) const;        // 向量点乘
    Vec3 cross(const Vec3& other) const;       // 向量叉乘
  
    float length() const;                     // 计算向量长度
    Vec3 normalize() const;                   // 返回单位化后的向量

    void print() const;                       // 打印向量内容
};

// 计算高斯权重函数
float gaussian_weight(float r, float sigma);

// 计算面的质心
Vec3 compute_face_centroid(const MyMesh& mesh, const MyMesh::FaceHandle& face);

// 计算面的真实法线
Vec3 compute_face_normal(const MyMesh& mesh, const MyMesh::FaceHandle& face);

// 载入 OBJ 文件
bool load_obj(const std::string& file_path, MyMesh& mesh);

// 保存 OBJ 文件
bool save_obj(const std::string& file_path, const MyMesh& mesh);

// 顶点投影到切平面（此处带上基准点 base_point，用于计算向量）
Vec3 project_to_tangent_plane(const Vec3& vertex, const Vec3& normal, const Vec3& base_point);

// 对网格进行去噪，返回平滑后的顶点坐标
std::vector<Vec3> smooth_mesh(const MyMesh& mesh, float sigma_f = 1.0f, float sigma_g = 1.0f);

// 设置日志记录文件
void setup_logging(std::ofstream& log_file);

// 记录顶点变化到日志文件
void log_vertex_changes(std::ofstream& log_file, const Vec3& original_vertex, const Vec3& smoothed_vertex, size_t index);

// 主去噪函数，执行载入、平滑、日志记录、保存等完整流程
void denoise_obj(const std::string& input_obj, const std::string& output_obj);

#endif // DENOISE_OBJ_HPP            // 预处理指令：防止头文件被重复包含的保护宏（结束）

```

#### 3)、denoise_obj.cpp

```c++
#include "denoise_obj.hpp"   // 包含 denoise_obj.hpp，使用其声明的函数和结构体
#include <iostream>          // 标准输入输出库
#include <fstream>           // 文件读写流
#include <cmath>             // 数学函数库
#include <vector>            // std::vector 容器
#include <utility>           // std::pair 等实用工具

// Vec3 构造函数实现
Vec3::Vec3() : x(0), y(0), z(0) {} // 默认构造函数，初始化 x、y、z 均为 0

Vec3::Vec3(float x, float y, float z) : x(x), y(y), z(z) {} // 带参构造函数

Vec3 Vec3::operator+(const Vec3& other) const {
    return Vec3(x + other.x, y + other.y, z + other.z); // 向量加法
}

Vec3 Vec3::operator-(const Vec3& other) const {
    return Vec3(x - other.x, y - other.y, z - other.z); // 向量减法
}

Vec3 Vec3::operator*(float scalar) const {
    return Vec3(x * scalar, y * scalar, z * scalar);     // 向量与标量相乘
}

float Vec3::dot(const Vec3& other) const {
    return x * other.x + y * other.y + z * other.z;      // 点乘：x1*x2 + y1*y2 + z1*z2
}

Vec3 Vec3::cross(const Vec3& other) const {
    // 叉乘公式： (y1*z2 - z1*y2, z1*x2 - x1*z2, x1*y2 - y1*x2)
    return Vec3(
        y * other.z - z * other.y,
        z * other.x - x * other.z,
        x * other.y - y * other.x
    );
}

float Vec3::length() const {
    // 向量长度：sqrt(x^2 + y^2 + z^2)
    return std::sqrt(x * x + y * y + z * z);
}

Vec3 Vec3::normalize() const {
    float len = length();         // 计算向量的长度
    if(len == 0) return Vec3(0,0,0); // 如果长度为 0，返回一个零向量，避免除以零
    return Vec3(x / len, y / len, z / len); // 向量每个分量除以长度
}

void Vec3::print() const {
    // 打印向量内容到标准输出
    std::cout << "(" << x << ", " << y << ", " << z << ")";
}

// 计算高斯权重
float gaussian_weight(float r, float sigma) {
    // 高斯公式：exp( - (r^2) / (2*sigma^2) )
    return std::exp(- (r * r) / (2 * sigma * sigma));
}

// 计算面的质心
Vec3 compute_face_centroid(const MyMesh& mesh, const MyMesh::FaceHandle& face) {
    Vec3 centroid;            // 用于累加顶点坐标
    int count = 0;            // 统计面对应的顶点数
    // 遍历该面上的所有顶点
    for (auto v_it = mesh.cfv_begin(face); v_it != mesh.cfv_end(face); ++v_it) {
        const auto& point = mesh.point(*v_it); // 取得顶点坐标
        centroid = centroid + Vec3(point[0], point[1], point[2]); // 累加到质心向量
        count++;
    }
    if(count > 0)
        return centroid * (1.0f / count); // 取平均值得到质心
    else
        return Vec3(0, 0, 0);            // 若无顶点，返回零向量
}

// 计算面的真实法线
Vec3 compute_face_normal(const MyMesh& mesh, const MyMesh::FaceHandle& face) {
    // 用一个临时容器存储面的所有顶点坐标
    std::vector<Vec3> face_vertices;
    for (auto v_it = mesh.cfv_begin(face); v_it != mesh.cfv_end(face); ++v_it) {
        const auto& point = mesh.point(*v_it);
        face_vertices.emplace_back(point[0], point[1], point[2]);
    }

    if (face_vertices.size() < 3) {
        // 如果面顶点少于 3，无法计算法线，返回零向量
        return Vec3(0, 0, 0);
    }

    // 假设是三角面，取前三个顶点
    Vec3 A = face_vertices[0];
    Vec3 B = face_vertices[1];
    Vec3 C = face_vertices[2];

    Vec3 AB = B - A;    // 向量 AB
    Vec3 AC = C - A;    // 向量 AC
    Vec3 normal = AB.cross(AC).normalize(); // 叉乘后单位化，得到真实法线

    return normal;      // 返回法线
}

// 计算面的面积
float compute_face_area(const MyMesh& mesh, const MyMesh::FaceHandle& face) {
    // 用一个临时容器存储面的所有顶点坐标
    std::vector<Vec3> face_vertices;
    for (auto v_it = mesh.cfv_begin(face); v_it != mesh.cfv_end(face); ++v_it) {
        const auto& point = mesh.point(*v_it);
        face_vertices.emplace_back(point[0], point[1], point[2]);
    }

    if (face_vertices.size() < 3) {
        // 面顶点不足以计算面积，返回 0
        return 0.0f;
    }

    // 假设是三角形，取前三个顶点计算面积
    Vec3 A = face_vertices[0];
    Vec3 B = face_vertices[1];
    Vec3 C = face_vertices[2];

    Vec3 AB = B - A;
    Vec3 AC = C - A;
    Vec3 cross_product = AB.cross(AC);  // 叉乘结果向量
    float area = 0.5f * cross_product.length(); // 三角形面积 = 0.5 * |AB x AC|

    return area;
}

// 载入 OBJ 文件
bool load_obj(const std::string& file_path, MyMesh& mesh) {
    if (!OpenMesh::IO::read_mesh(mesh, file_path)) {
        // 如果读取失败，输出错误信息
        std::cerr << "Error loading OBJ file: " << file_path << std::endl;
        return false;
    }
    return true; // 读取成功，返回 true
}

// 保存 OBJ 文件
bool save_obj(const std::string& file_path, const MyMesh& mesh) {
    if (!OpenMesh::IO::write_mesh(mesh, file_path)) {
        // 如果写文件失败，输出错误信息
        std::cerr << "Error saving OBJ file: " << file_path << std::endl;
        return false;
    }
    return true; // 保存成功，返回 true
}

// 顶点投影到切平面（带基准点 base_point）
Vec3 project_to_tangent_plane(const Vec3& vertex, const Vec3& normal, const Vec3& base_point) {
    // 先计算 (vertex - base_point)，得到相对 base_point 的向量
    Vec3 vertex_to_base = vertex - base_point;
    // 计算在法线方向上的投影量，然后减去该投影量，实现投影到切平面
    return vertex - normal * normal.dot(vertex_to_base);
}

// 对网格进行去噪
std::vector<Vec3> smooth_mesh(const MyMesh& mesh, float sigma_f, float sigma_g) {
    std::vector<Vec3> smoothed_vertices;               // 存储去噪后的顶点坐标
    smoothed_vertices.reserve(mesh.n_vertices());      // 预留容量，避免动态扩张

    // 预先计算每个面的质心与法线，以及每个面的面积
    std::vector<std::pair<Vec3, Vec3>> face_centroids_normals; // (质心, 法线)
    std::vector<float> face_areas;                               // 面积数组
    for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
        Vec3 centroid = compute_face_centroid(mesh, *f_it);   // 计算质心
        Vec3 normal = compute_face_normal(mesh, *f_it);       // 计算真实法线
        face_centroids_normals.emplace_back(centroid, normal);

        float area = compute_face_area(mesh, *f_it);          // 计算面面积
        face_areas.push_back(area);
    }

    // 遍历每个顶点进行去噪
    size_t total_vertices = mesh.n_vertices();  // 顶点总数
    size_t current_vertex = 0;                  // 当前处理的顶点计数
    std::cout << "Smoothing vertices: " << std::endl; // 在控制台提示

    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it, ++current_vertex) {
        // 获取当前顶点坐标
        const auto& point = mesh.point(*v_it);
        Vec3 vertex(point[0], point[1], point[2]);

        // 用于存储对该顶点的各个面权重及投影结果
        std::vector<float> f_weights;      // 空间权重
        std::vector<float> g_weights;      // 影响权重
        std::vector<Vec3> projections;     // 存储投影点
        std::vector<float> triangle_areas; // 存储面面积

        f_weights.reserve(face_centroids_normals.size());
        g_weights.reserve(face_centroids_normals.size());
        projections.reserve(face_centroids_normals.size());
        triangle_areas.reserve(face_centroids_normals.size());

        // 针对每一个面，计算其权重和投影
        for (size_t i = 0; i < face_centroids_normals.size(); ++i) {
            const Vec3& centroid = face_centroids_normals[i].first;  // 面质心
            const Vec3& normal = face_centroids_normals[i].second;   // 面真实法线

            // 计算空间距离，用于 f_weights
            float distance = (vertex - centroid).length();
            float f_weight = gaussian_weight(distance, sigma_f);
            f_weights.push_back(f_weight);

            // 将顶点投影到该面的切平面
            Vec3 projection = project_to_tangent_plane(vertex, normal, centroid);
            projections.push_back(projection);

            // 计算投影距离，用于 g_weights
            float projection_dist = (projection - vertex).length();
            float g_weight = gaussian_weight(projection_dist, sigma_g);
            g_weights.push_back(g_weight);

            // 记录该面的面积
            float area = face_areas[i];
            triangle_areas.push_back(area);
        }

        // 计算加权和
        float total_weight = 0.0f;        // 总权重
        Vec3 weighted_sum(0, 0, 0);       // 投影后的加权累加
        for (size_t i = 0; i < projections.size(); ++i) {
            float weight = triangle_areas[i] * f_weights[i] * g_weights[i];
            weighted_sum = weighted_sum + (projections[i] * weight);
            total_weight += weight;
        }

        // 如果总权重不为零，则计算新的平滑顶点，否则保持原顶点
        if(total_weight != 0.0f){
            smoothed_vertices.push_back(weighted_sum * (1.0f / total_weight));
        }
        else{
            smoothed_vertices.push_back(vertex);
        }

        // 每处理 1000 个顶点，打印一次进度
        if(current_vertex % 1000 == 0){
            std::cout << "\rProcessed " << current_vertex << " / " 
                      << total_vertices << " vertices." << std::flush;
        }
    }
    // 最后输出一次完成信息
    std::cout << "\rProcessed " << total_vertices << " / " 
              << total_vertices << " vertices." << std::endl;

    return smoothed_vertices; // 返回所有平滑后的顶点数据
}

// 设置日志记录
void setup_logging(std::ofstream& log_file) {
    log_file.open("denoise_log.txt");            // 打开 denoise_log.txt 文件用于写日志
    if (!log_file.is_open()) {                   // 若文件无法打开，输出错误并退出
        std::cerr << "Error opening log file." << std::endl;
        exit(1);
    }
}

// 记录顶点变化
void log_vertex_changes(std::ofstream& log_file, const Vec3& original_vertex, const Vec3& smoothed_vertex, size_t index) {
    // 将原始顶点和平滑后顶点信息写入日志文件
    log_file << "Vertex " << index << " - Original: (" << original_vertex.x << ", " 
             << original_vertex.y << ", " << original_vertex.z << ") | Smoothed: ("
             << smoothed_vertex.x << ", " << smoothed_vertex.y << ", " 
             << smoothed_vertex.z << ")\n";
}

// 主去噪函数
void denoise_obj(const std::string& input_obj, const std::string& output_obj) {
    MyMesh mesh;  // 创建一个网格对象

    // 载入 OBJ 文件
    if (!load_obj(input_obj, mesh)) {
        return;   // 若载入失败，直接返回
    }

    // 调用去噪函数，得到平滑后的顶点
    std::vector<Vec3> smoothed_vertices = smooth_mesh(mesh, 1.0f, 1.0f);

    // 打开日志文件
    std::ofstream log_file;
    setup_logging(log_file);

    // 遍历所有顶点，将平滑前后的坐标记录到日志
    size_t index = 0;
    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it, ++index) {
        const auto& original = mesh.point(*v_it);
        Vec3 original_vertex(original[0], original[1], original[2]);
        Vec3 smoothed_vertex = smoothed_vertices[index];
        log_vertex_changes(log_file, original_vertex, smoothed_vertex, index);
    }
    log_file.close(); // 日志记录完毕后关闭文件

    // 使用平滑后的坐标更新网格
    index = 0;
    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it, ++index) {
        mesh.set_point(*v_it, OpenMesh::Vec3f(smoothed_vertices[index].x, 
                                              smoothed_vertices[index].y, 
                                              smoothed_vertices[index].z));
    }

    // 保存去噪后的 OBJ 文件
    if (!save_obj(output_obj, mesh)) {
        return;   // 若保存失败，直接返回
    }

    std::cout << "Denoised OBJ saved to " << output_obj << std::endl; // 输出保存成功的信息
}

```

### 4、实验结果

由于平滑的结果不是特别明显，这里除了直接比较对比图像，还将各个顶点的位置变化记录，输入日志文件，具体内容如下：

#### 1)、图像对比

![网格去噪_正面对比](https://daetz-image.oss-cn-hangzhou.aliyuncs.com/img/202412281944817.png)

![image-20241228193826166](https://daetz-image.oss-cn-hangzhou.aliyuncs.com/img/202412281938398.png)

![image-20250104114922951](https://daetz-image.oss-cn-hangzhou.aliyuncs.com/img/202501041149180.png)

![image-20250104115058843](https://daetz-image.oss-cn-hangzhou.aliyuncs.com/img/202501041150084.png)

![image-20250104115247255](https://daetz-image.oss-cn-hangzhou.aliyuncs.com/img/202501041152546.png)

![image-20250104115359582](https://daetz-image.oss-cn-hangzhou.aliyuncs.com/img/202501041153966.png)

#### 2)、顶点移动日志

```
Vertex 0 - Original: (-52.2326, 67.5457, -57.677) | Smoothed: (-52.114, 67.551, -57.7187)
Vertex 1 - Original: (-6.7165, 17.1235, 29.4826) | Smoothed: (-6.71261, 17.1312, 29.4632)
Vertex 2 - Original: (-53.6263, 66.6593, -57.0996) | Smoothed: (-53.5737, 66.5822, -57.1845)
Vertex 3 - Original: (-52.6775, 66.3879, -57.0771) | Smoothed: (-52.6105, 66.3819, -57.1432)
Vertex 4 - Original: (-53.3113, 67.1217, -57.4919) | Smoothed: (-53.2241, 67.076, -57.5762)
Vertex 5 - Original: (-51.7571, 67.1619, -57.0316) | Smoothed: (-51.7126, 67.1805, -57.1477)
Vertex 6 - Original: (26.8491, -9.25235, 19.0221) | Smoothed: (26.8152, -9.2391, 18.9756)
Vertex 7 - Original: (-51.9157, 68.0805, -57.2048) | Smoothed: (-51.8378, 68.0517, -57.3101)
Vertex 8 - Original: (-52.7414, 68.3093, -57.5367) | Smoothed: (-52.6312, 68.2986, -57.5991)
Vertex 9 - Original: (-54.4458, 67.1714, -56.7764) | Smoothed: (-54.4768, 67.1183, -56.8996)
"""
Vertex 21579 - Original: (-8.29693, 63.1762, -31.1635) | Smoothed: (-8.34546, 63.1384, -31.1378)
Vertex 21580 - Original: (-5.63182, 64.2246, -33.6254) | Smoothed: (-5.5907, 64.2076, -33.5633)
Vertex 21581 - Original: (-8.6238, 62.8402, -30.6838) | Smoothed: (-8.67226, 62.7952, -30.6705)

```

### 5、使用说明

本文的结构介绍如下

```less
├── CMakeLists.txt
├── include
│   └── denoise_obj.hpp
├── log
│   └── denoise_log.txt
├── Non-Iterative, Feature-Preserving Mesh Smoothing.ipynb
├── obj
│   ├── armadillo_denoised.obj
│   └── armadillo.obj
├── README.md
├── requirements.txt
├── result
│   ├── 网格去噪_正面对比.png
│   └── 网格去噪_背面对比.png
└── src
    ├── denoise_obj.cpp
    └── main.cpp
```

+ `CMakeLists.txt` 包含所需的 C++ 标准、源文件列表、依赖库等信息
+ `include/denoise_obj.hpp`声明了与网格去噪（Mesh Denoising）相关的类、函数和数据结构
+ `log/denoise_log.txt` 记录了网格去噪过程中所有顶点的移动情况
+ `Non-Iterative, Feature-Preserving Mesh Smoothing.ipynb`使用Jupyter Notebook 文件可视化网格去噪的原理
+ `obj/` 用于存放3D模型，其中 `armadillo_denoised.obj`表示已经去噪后的文件
+ `README.md` 表示项目的概述、安装指南、使用说明、功能介绍
+ `requirements.txt` 包含需要运行Jupyter Notebook 文件所需的依赖
+ `result/` 为了方便比较，我们把网格去噪后的结果正面和背面的对比图放入该文件下

#### 1)、python

```bash
pip install -r requirements.txt 
# start Jupyter Notebook 
```

#### 2)、c++

```bash
cd OpenMesh-11.0.0
mkdir build && cd build
cmake ..
make
sudo make install
```

![image-20241228200133841](https://daetz-image.oss-cn-hangzhou.aliyuncs.com/img/202412282001985.png)

```bash
# build
mkdir build 
cd build 
cmake ..
make

# remove 
make clean
cd ..
rm -rf build

# use params run
./denoise_obj ../obj/armadillo.obj ../obj/armadillo_denoised.obj
```
