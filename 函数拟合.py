##################################################
### Project: 空间统计基础函数拟合
### Script purpose: 《空间统计基础》疫情不均衡研究，实现函数拟合
### Date: 2020-07
### Author: Minghao Du
### email: 0108170318@csu.edu.cn / 695948191@qq.com
### Code: Python 3
##################################################

############## Load the packages ##############
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score

############## 读取数据 ##############
# 衰减
x_sars_de = np.arange(0,9,1)
y_sars_de = np.array([1601, 1083, 743, 550, 373, 338, 284, 236, 182])
y_sars_in_fit = np.array([1, 25, 85, 236])
x_mers_de = np.arange(0,4,1)
y_mers_de = np.array([448, 338, 291, 271])
x_h1n1_de = np.arange(0,9,1)
y_h1n1_de = np.array([2790, 2013, 1520, 1167, 969, 827, 721, 690, 438])
x_ebola_de = np.arange(0,5,1)
y_ebola_de = np.array([1779, 1133, 895, 652, 608])

# 增长
x_sars_in = np.arange(0,5,1)
y_sars_in = np.array([71, 1, 25, 85, 236])
x_sars_in_fit = np.arange(1,5,1)
y_sars_in_fit = np.array([1, 25, 85, 236])
x_mers_in = np.arange(0,32,1)
y_mers_in = np.array([4, 1 ,1, 1, 4, 11, 15, 27, 30, 29, 24, 26, 32, 18, 22, 31, 35,
                      26, 30, 34, 26, 37, 34, 44, 49, 21, 22, 38, 24, 39, 49, 56])
x_mers_in_fit = np.arange(3,32,1)
y_mers_in_fit = np.array([1, 4, 11, 15, 27, 30, 29, 24, 26, 32, 18, 22, 31, 35,
                      26, 30, 34, 26, 37, 34, 44, 49, 21, 22, 38, 24, 39, 49, 56])
x_h1n1_in = np.arange(0,18,1)
y_h1n1_in = np.array([38, 48, 158, 184, 210, 197, 241, 252, 297, 316, 410, 237, 303,
                      284, 255, 249, 261, 238])
x_ebola_in = np.arange(0,14,1)
y_ebola_in = np.array([6, 58, 14, 7, 18, 18, 16, 23, 59, 152, 243, 268, 261, 352])
x_ebola_in_fit = np.arange(3,14,1)
y_ebola_in_fit = np.array([7, 18, 18, 16, 23, 59, 152, 243, 268, 261, 352])
x_cov_in = np.arange(0,6,1)
y_cov_in = np.array([726, 547, 2082, 6786, 10457, 11424])
x_cov_in_fit = np.arange(1,6,1)
y_cov_in_fit = np.array([547, 2082, 6786, 10457, 11424])

############## 指数线性混合衰减函数的方程写出 ##############
# m=beta1, n=beta2, o=q0, p=a
def f(x, m, n, o, p):
    return np.piecewise(x, [x < p, x >= p], [lambda x:m*np.exp(-n*x)-o*x+o/(2*p)*(x**2),
                                               lambda x: m*np.exp(-n*x)-o*p/2])

############## 把单指数衰减函数的方程写出 ##############
def g(x, m, n):
    return m*np.exp(-n*x)

############## 把双指数衰减函数的方程写出 ##############
# m=N0, n=q0, p=p, o=a
def h(x, m, n, o, p):
    return np.piecewise(x, [x < o, x >= o], [lambda x:m*np.exp((x**2)*n/(2*o)-(n+p)*x),
                                               lambda x: m*np.exp(-p*x-n*o/2)])

############## 把逻辑斯蒂方程写出 ##############
# m=alpha1, n=alpha2, o=alpha3
def l(x, m, n, o):
    return m/(1+np.exp(n-o*x))

############## 函数拟合 ##############
# 衰减
p_sars_f , e_sars_f = optimize.curve_fit(f, x_sars_de, y_sars_de, p0=(1450, 0.1, 200, 2))
p_mers_f , e_mers_f = optimize.curve_fit(f, x_mers_de, y_mers_de, p0=(500, 0.1, 100, 1.9))
p_h1n1_f , e_h1n1_f = optimize.curve_fit(f, x_h1n1_de, y_h1n1_de, p0=(2600, 0.06, 600, 3))
p_ebola_f , e_ebola_f = optimize.curve_fit(f, x_ebola_de, y_ebola_de, p0=(1778, 0.1, 400, 1.5))
p_sars_g , e_sars_g = optimize.curve_fit(g, x_sars_de, y_sars_de, p0=(1550, 0.3))
p_mers_g , e_mers_g = optimize.curve_fit(g, x_mers_de, y_mers_de, p0=(432, 0.2))
p_h1n1_g , e_h1n1_g = optimize.curve_fit(g, x_h1n1_de, y_h1n1_de, p0=(2600, 0.2))
p_ebola_g , e_ebola_g = optimize.curve_fit(g, x_ebola_de, y_ebola_de, p0=(1700, 0.3))
p_sars_h , e_sars_h = optimize.curve_fit(h, x_sars_de, y_sars_de, p0=(1600, 1, 0.3, 0.28), maxfev=1000000)
p_mers_h , e_mers_h = optimize.curve_fit(h, x_mers_de, y_mers_de, p0=(448, 1, 0.2, 0.11), maxfev=1000000)
p_h1n1_h , e_h1n1_h = optimize.curve_fit(h, x_h1n1_de, y_h1n1_de, p0=(2700, 0.8, 0.4, 0.23), maxfev=1000000)
p_ebola_h , e_ebola_h = optimize.curve_fit(h, x_ebola_de, y_ebola_de, p0=(1700, 1, 0.1, 0.3), maxfev=1000000)

# 增长
p_sars_l , e_sars_l = optimize.curve_fit(l, x_sars_in_fit, y_sars_in_fit, p0=(495, 5.99, 1.47))
p_mers_l , e_mers_l = optimize.curve_fit(l, x_mers_in_fit, y_mers_in_fit, p0=(32.5, 3.4, 0.98), maxfev=10000)
p_h1n1_l , e_h1n1_l = optimize.curve_fit(l, x_h1n1_in, y_h1n1_in, p0=(282, 1.61, 0.66))
p_ebola_l , e_ebola_l = optimize.curve_fit(l, x_ebola_in_fit, y_ebola_in_fit, p0=(300, 10.4, 1))
p_cov_l , e_cov_l = optimize.curve_fit(l, x_cov_in_fit, y_cov_in_fit, p0=(11600, 5.1, 1.82))

############## 计算R² ##############
def r2(x, y, p1, p2, p3):
    y_pred_f = f(x, *p1)
    r2_score_f = r2_score(y, y_pred_f)
    print('r2_f:', r2_score_f)
    y_pred_g = g(x, *p2)
    r2_score_g = r2_score(y, y_pred_g)
    print('r2_g:', r2_score_g)    
    y_pred_h = h(x, *p3)
    r2_score_h = r2_score(y, y_pred_h)
    print('r2_h:', r2_score_h)
    return

r2(x_sars_de, y_sars_de, p_sars_f, p_sars_g, p_sars_h)
r2(x_mers_de, y_mers_de, p_mers_f, p_mers_g, p_mers_h)
r2(x_h1n1_de, y_h1n1_de, p_h1n1_f, p_h1n1_g, p_h1n1_h)
r2(x_ebola_de, y_ebola_de, p_ebola_f, p_ebola_g, p_ebola_h)

############## 绘图函数 ##############
# 衰减
def model_plot(x, y, LABELS, TITLE, p1, p2, p3, a):
    xd = np.linspace(0, a, 100)
    plt.scatter(x, y, c='blue')
    plt.plot(x, y, c='black')
    plt.plot(xd, f(xd, *p1), c='red')
    plt.plot(xd, g(xd, *p2), c='green')
    plt.plot(xd, h(xd, *p3), c='blue')
    plt.xticks(np.arange(0,a,1), labels = LABELS, rotation=45)
    plt.xlabel('Year',fontsize=16)
    plt.ylabel('Number',fontsize=16)
    plt.title(TITLE,fontsize=16, x=-0.1, y=0.95, fontweight='black')

# 增长
def model_plot_in(x, y, LABELS, TITLE, p, a):
    xd = np.linspace(0, a, 100)
    plt.scatter(x, y, c='blue')
    plt.plot(x, y, c='black')
    plt.plot(xd, l(xd, *p), c='red')
    plt.xticks(np.arange(0,a,1), labels = LABELS, rotation=45)
    plt.xlabel('Month',fontsize=16)
    plt.ylabel('Number',fontsize=16)
    plt.title(TITLE,fontsize=16, x=-0.1, y=0.95, fontweight='black')

############## 绘图 ##############
sars_label = ['2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010',
              '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018',
              '2019']
mers_label = ['2015', '2016', '2017', '2018', '2019']       
h1n1_label = ['2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018',
              '2019', '2020']
ebola_label = ['2015', '2016', '2017', '2018', '2019', '2020']

sars_label_in = ['2003/01', '2003/02', '2003/03', '2003/04', '2003/05']
mers_label_in = ['2013/01', '2013/02', '2013/03', '2013/04', '2013/05', '2013/06', '2013/07', 
                 '2013/08', '2013/09', '2013/10', '2013/11', '2013/12', '2014/01', '2014/02', 
                 '2014/03', '2014/04', '2014/05', '2014/06', '2014/07', '2014/08', '2014/09', 
                 '2014/10', '2014/11', '2014/12', '2015/01', '2015/02', '2015/03', '2015/04', 
                 '2015/05', '2015/06', '2015/07', '2015/08']
h1n1_label_in = ['2009/03', '2009/04', '2009/05', '2009/06', '2009/07', '2009/08', '2009/09', 
                 '2009/10', '2009/11', '2009/12', '2010/01', '2010/02', '2010/03', '2010/04', 
                 '2010/05', '2010/06', '2010/07', '2010/08']
ebola_label_in = ['2013/12', '2014/01', '2014/02', '2014/03', '2014/04', '2014/05', '2014/06', 
                  '2014/07', '2014/08', '2014/09', '2014/10', '2014/11', '2014/12', '2015/01']
cov_label_in = ['2020/01', '2020/02', '2020/03', '2020/04', '2020/05', '2020/06']

plt.figure(1, dpi=600, figsize=(35,20))

ax1 = plt.subplot(331)
model_plot_in(x_sars_in, y_sars_in, sars_label_in, 'a', p_sars_l, 5)
ax2 = plt.subplot(332)
model_plot_in(x_mers_in, y_mers_in, mers_label_in, 'b', p_mers_l, 32)
ax3 = plt.subplot(333)
model_plot_in(x_h1n1_in, y_h1n1_in, h1n1_label_in, 'c', p_h1n1_l, 18)
ax4 = plt.subplot(334)
model_plot_in(x_ebola_in, y_ebola_in, ebola_label_in, 'd', p_ebola_l, 14)
ax5 = plt.subplot(335)
model_plot_in(x_cov_in, y_cov_in, cov_label_in, 'e', p_cov_l, 6)

ax6 = plt.subplot(336)
model_plot(x_sars_de, y_sars_de, sars_label, 'f', p_sars_f, p_sars_g, p_sars_h, 10)
ax7 = plt.subplot(337)
model_plot(x_mers_de, y_mers_de, mers_label, 'g', p_mers_f, p_mers_g, p_mers_h, 5)
ax8 = plt.subplot(338)
model_plot(x_h1n1_de, y_h1n1_de, h1n1_label, 'h', p_h1n1_f, p_h1n1_g, p_h1n1_h, 10)
ax9 = plt.subplot(339)
model_plot(x_ebola_de, y_ebola_de, ebola_label, 'i', p_ebola_f, p_ebola_g, p_ebola_h, 6)

############################结束############################












