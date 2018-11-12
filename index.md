## An introduction to Kalman Filter Based tracking
我与卡尔曼滤波

卡尔曼滤波的核心公式

$X_{t}=&AX_{t-1} & +W_{t}\label{eq.state}$
Y_{t}=&CX_{t} &+ V_{t}\label{eq.measure}
\end{eqnarray}

状态向量$X \in R^{M}$, $W \in R^{M}$为高斯白噪声，符合正态分布$w \sim \mathcal{N}(0, \Sigma_{w}) $.
则给定前一时刻$X_{t-1}$状态时，$p(X_t | X_{t-1}) = \mathcal{N}(X_t | AX_{t-1}, \Sigma_{w}) $。
若假设$p(X_{t-1})=\mathcal{N}(X_{t-1} | \boldsymbol{\mu}_{t-1}, \Sigma_{t-1})$
根据状态转移矩阵的线性, $p(X_{t-1})$ 以及 $p(X_t | X_{t-1})$ 构成线性高斯模型 (Linear-Gaussian model),
推出$p(X_t)=\mathcal{N}(X_t | A\mu_{t-1}, A\Sigma_{t-1}A^T+\Sigma_w).$
可见$p(X_t)$为仍然同先验概率$p(X_{t-1})$一致，同为高斯分布。

$V \in R^{N}$为高斯白噪声，符合正态分布$v \sim \mathcal{N}(0, \Sigma_{v}) $.
且给定$X_t$时，$p(Y_t | X_{t}) = \mathcal{N}(Y_t | CX_{t}, \Sigma_{v}) $
又通过\eqref{eq_measure}，线性变换公式，又构成了对$Y_t$的线性高斯模型。
推出$p(Y_t)=\mathcal{N}(C\mu_{t}, C\Sigma_{X_t}C^T+\Sigma_v).$，这其中 $\mu_{t}=A\mu_{t-1}$, $\Sigma_{X_t}=A\Sigma_{t-1}A^T+\Sigma_w.$

注意根据$p(X_t)$及$p(Y_t | X_t)$可以推出联合概率密度$p(X_t, Y_t)$的表达式。

 $Z_t=\begin{bmatrix}X_{t}\\Y_{t}\end{bmatrix} $, 则$p(Z_t)=\mathcal{N}(\mu_z, \Sigma_z)$。
 $\mu_z = \begin{bmatrix}\mu_{x}\\ \mu_{y}\end{bmatrix} $，
 $\Sigma_z = \begin{bmatrix} &\Sigma_x & \Sigma_xC^T \\  &C\Sigma_x & \Sigma_{y|x}+C\Sigma_xC^T \end{bmatrix}$
 其Precision Matrix 精度矩阵为$\Lambda_z$, 则有 $\Sigma_z ={ \Lambda_z}^{-1}$, 根据分块矩阵求逆公式，可以得出：
 \beq
 \Lambda_z=\left( \begin{array}{rl} \Sigma_x^{-1}+C^T\Sigma_{y|x}C, & -C^T\Sigma_{y|x}  \\
                                               -C\Sigma_{y|x},                              & \Sigma_{y|x} \end{array} \right)
 \eeq
 已知联合分布$p(X_t, Y_t)$则可以将$Y_t$视为给定值，通过配方法推出后验分布$p(X_t | Y_t)$。
 \beq
 \label{eq.postor}
 \begin{split}
 p(X_t | Y_t)=\mathcal{N}(X_t |& \mu_x + \Sigma_xC^T(\Sigma_{y|x}+C\Sigma_xC^T)^{-1}(Y_t - \mu_y),\\
                                                 &\Sigma_x - \Sigma_xC^T(\Sigma_{y|x}+C\Sigma_xC^T)^{-1}C\Sigma_x)
 \end{split}
 \eeq
 
 令$K=\Sigma_xC^T(\Sigma_{y|x}+C\Sigma_xC^T)^{-1}$,则\eqref{eq.postor}变为：
 \beq
 \label{eq.postor.k}
 p(X_t | Y_t)=\mathcal{N}(X_t | \mu_x + K(Y_t - \mu_y),  (\mathit{I} - KC)\Sigma_x)
 \eeq
 
 可见公式\eqref{eq.postor.k}的均值和协方差矩阵与Kalman filter的递推公式紧密相关。
 若将$p(X_t)$看成是对预测的状态向量的先验概率分布，那么通过观测值$Y_t$，以及观测值与预测值的线性关系\eqref{eq.measure},
 可以推出似然函数概率密度$p(Y_t | X_t)$，用其乘以$p(X_t)$（修正先验）得到后验分布$p(X_t | Y_t)$，这就是贝叶斯方法在测量跟踪中的具体应用。
因为后验分布为高斯，那么最大化后验分布的$X_t $即为其均值$\mu_{x|y} =  \mu_x + K(Y_t - C\mu_x)$。将此\textcolor{red}{后验又作为$t+1$时刻的先验$p(X_{t})$}，就又可以递推出下一时刻的似然函数和后验分布。


值得注意的是联合分布$p(Z_t) = p(\begin{bmatrix}X_{t}\\Y_{t}\end{bmatrix} )$的关键参数协方差矩阵可以直接推导，也可以通过配方法先求精度矩阵，然后通过分块矩阵求逆得到。
法一：直接求解

\beq
\begin{split}
\Sigma_z &= E \left[ (Z - \mu_z) (Z - \mu_z)^T \right] \\
                &= E \left[ \left(\begin{array}{c} X-\mu_x \\ Y - \mu_y \end{array}\right)    \left( (X-\mu_x)^T, (Y - \mu_y)^T \right)\right] \\
                &= \left( \begin{array}{c} 
                             E \left[  (X - \mu_x)(X - \mu_x)^T \right],       E \left[  (X - \mu_x)(Y - \mu_y)^T \right] \\
                             E \left[  (Y - \mu_y)(X - \mu_x)^T \right],       E \left[  (Y - \mu_y)(Y - \mu_y)^T \right]
                             \end{array}
                             \right)      \\   
                &= \left( \begin{array}{cc} 
                             \Sigma_x     &  \Sigma_{xy} \\
                             \Sigma_{yx}  & \Sigma_{y}     
                             \end{array}
                             \right)     
\end{split}
\eeq
$\Sigma^T_{xy} = \Sigma_{yx}$
因为 $Y_t = CX_t + V_t$, 故 $\mu_y = E[Y_t] = CE[X_t] + E[V_t] = C\mu_x$
\beq
\begin{split}
&E \left[  (Y - \mu_y)(Y - \mu_y)^T \right] \\
&= E  \left[  (CX_t + V_t - C\mu_x)(CX_t + V_t - C\mu_x)^T \right] \\
&= E [C(X_t-\mu_x)(X_t-\mu_x)^TC^T] + E[V_t{V_t}^T] \\
&=C\Sigma_xC^T + \Sigma_v
\end{split}
\eeq

\beq
\begin{split}
&E \left[  (X - \mu_x)(Y - \mu_y)^T \right] \\
&= E  \left[  (X_t - \mu_x)(CX_t + V_t - C\mu_x)^T \right] \\
& = E [(X_t-\mu_x)(X_t-\mu_x)^TC^T] + E[X_t-\mu_x]E^T[{V_t}] \\
& =\Sigma_xC^T 
\end{split}
\eeq
因此
\beq
\Sigma_z = \left( \begin{array}{rl} 
                             \Sigma_x,   &   \Sigma_{x}C^T  \\
                            C\Sigma_{x}, &  C\Sigma_{x}C^T + \Sigma_v
                             \end{array}
                             \right)   
\eeq


\section{Iteration}
\begin{numcases}{}
  \bx_{t|t-1} &= $A\bx_{t-1} + W_t$ \label{eq_state} \\
% \notag  \\
 \by_{t}      &= $C\bx_{t|t-1} +  V_t$ \label{eq_measure}
\end{numcases} 

%\beq
%\begin{cases} \bx_{t|t-1} &= A\bx_{t-1} + W_t \\ 
%                       \by_{t}      &=C\bx_{t|t-1} +  V_t
% \end{cases}
% \eeq
如果将$p(\bx_{t-1})$高斯分布的均值写成$\hat{\bx}_{t-1}$，协方差写成$\Sigma_{t-1}$, 
将预测公式中的$\bx_{t}$改写成$\bx_{t|t-1}$,得到$\bx_{t|t-1} = A\bx_{t-1} + W_t$，此间的下标$t|t-1$是为了凸显该变量来自$t-1$时刻的预测先验，
与随后介绍$t$时刻的后验估计随机变量$\bx_{t}$做区别。
则根据线性高斯模型
\beq
\begin{cases} p(\bx_{t-1})& =\mathcal{N} (\hat{\bx}_{t-1}, \Sigma_{t-1}) \\ 
                       p(\bx_{t|t-1} | \bx_{t-1}) &= \mathcal{N} (A\bx_{t-1}, \Sigma_w)
 \end{cases}
 \eeq

推出预测先验概率
\beq
\label{eq.tprior}
\begin{split}
p(\bx_{t|t-1})&= \mathcal{N}(A\hat{\bx}_{t-1}, A\Sigma_{t-1}A^T+\Sigma_w) \\
                   &= \mathcal{N}(\hat{\bx}_{t|t-1}, \Sigma_{t|t-1})
\end{split}
\eeq
根据测量方程\eqref{eq_measure}, 可推在给定预测值$\bx_{t|t-1}$时$\by_t$的概率密度为
\beq
\label{eq.like}
p(\by_t|\bx_{t|t-1})=\mathcal{N}(C\bx_{t|t-1}, \Sigma_{v})
\eeq
联立等式\eqref{eq.tprior}和\eqref{eq.like}得到又一个线性高斯模型
\beq
\begin{cases} p(\bx_{t|t-1})          &= \mathcal{N}(\hat{\bx}_{t|t-1}, \Sigma_{t|t-1}) \\
                       p(\by_t|\bx_{t|t-1}) &=\mathcal{N}(C\bx_{t|t-1}, \Sigma_{v})
 \end{cases}
 \eeq
推出测量的概率分布和$t$时刻状态向量$\bx_t$的后验分布$p(\bx_t|\by_t)$
\beq
\label{eq.mdist}
p(\by_{t})= \mathcal{N}(C\hat{\bx}_{t|t-1}, C\Sigma_{t|t-1}C^T+\Sigma_v) 
\eeq

 \beq
 \label{eq.postor}
 \begin{split}
 p(\bx_t | \by_t)&=\mathcal{N}(\bx_t | ,  \hat{\bx}_{t} ,  \Sigma_{t} )\\
 \hat{\bx}_{t} &= \hat{\bx}_{t|t-1} + \Sigma_{t|t-1}C^T(C\Sigma_{t|t-1}C^T+\Sigma_v)^{-1}(\by_t - C\hat{\bx}_{t|t-1} )\\
 \Sigma_{t}     &= \Sigma_{t|t-1} - \Sigma_{t|t-1}C^T(C\Sigma_{t|t-1}C^T+\Sigma_v)^{-1}C\Sigma_{t|t-1}
 \notag
 \end{split}
 \eeq
 令 $K = \Sigma_{t|t-1}C^T(C\Sigma_{t|t-1}C^T+\Sigma_v)^{-1}$
则有
 \beq
 \label{eq.postor}
 \begin{split}
 \hat{\bx}_{t} &= \hat{\bx}_{t|t-1} + K(\by_t - C\hat{\bx}_{t|t-1} )\\
 \Sigma_{t}   &=(\mathit{I} - KC)\Sigma_{t|t-1}
 \notag
 \end{split}
 \eeq
 将$ p(\bx_t | \by_t)$作为新的先验$p(\bx_t)=\mathcal{N}(\bx_t | ,  \hat{\bx}_{t} ,  \Sigma_{t} )$，求解$t+1$的预测、似然和后验。
 
 算法：
 1 给定初始时刻的$p(\bx_0)$均值和方差，$(\hat{\bx}_{0}, \Sigma_0)$ \\
 2 for t = 1 : n \\
    将$(\hat{\bx}_{t-1}, \Sigma_{t-1})$代入公式\eqref{eq.tprior}，得到预测先验$p(\bx_{t|t-1})$的均值和方差$(\hat{\bx}_{t|t-1}), \Sigma_{t|t-1})$.
\begin{eqnarray}
         \hat{\bx}_{t|t-1} &=& A\hat{\bx}_{t-1}, \\
         \Sigma_{t|t-1}   &=&A\Sigma_{t-1}A^T+\Sigma_w
\end{eqnarray} 
    计算Kalman增益矩阵K：
    \beq
    K = \Sigma_{t|t-1}C^T(C\Sigma_{t|t-1}C^T+\Sigma_v)^{-1}
    \eeq 
    得出后验概率$p(\bx_t | \by_t)$的均值和协方差$(\hat{\bx}_{t}, \Sigma_{t})$，并将其均值作为$t$时刻的估计值。
\begin{eqnarray}
 \hat{\bx}_{t} &=& \hat{\bx}_{t|t-1} + K(\by_t - C\hat{\bx}_{t|t-1} )\\
 \Sigma_{t}   &=& (\mathit{I} - KC)\Sigma_{t|t-1}
\end{eqnarray} 
   将后验概率$p(\bx_t | \by_t)$的均值和协方差传递给滤波后的$p(\bx_t)=\mathcal{N}(\hat{\bx}_{t}, \Sigma_{t} )$ ，进入下一次迭代。
 此算法中的五个公式，即为Kalman Filter的五公式，K矩阵建立了预测和测量之间的联系。
 
 后验的高斯均值，相当于是建立在多个高斯分布的采样点上的最大值，暗含Monte carol采样的定理。但此采样必定收敛于最值（单模态）。
 这一句话，需要用采样点来反映随机特性，将其在图像上标示出来，让大家理解什么是加噪声的概念。
 
 其实Kalman滤波的核心思想是高斯随机过程，在目标动态运动，噪声随机的情况下，如何以最大的可能性跟踪目标，抵御干扰，那么就需要建立一个概率模型，以最大（后验）概率的参数估计来确定目标当前状态的向量。所以在实验的时候必须反映出目标位置的随机特性。
 


You can use the [editor on GitHub](https://github.com/joeyee/joeyee.github.io/edit/master/index.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/joeyee/joeyee.github.io/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://help.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.
