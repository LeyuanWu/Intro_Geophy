# 地磁轴向偶极子场：五个核心公式

> 轴向偶极子 $\mathbf{m} = m\hat{\mathbf{z}}$；场绕轴对称，**与经度 $\lambda$ 无关**；$\varphi$ 为纬度（北纬为正），$\hat{\boldsymbol{\varphi}}$ 指向北。

## 1. 偶极子场的纬度表达式

$$
\begin{cases}
   & \mathbf{m}\cdot\hat{\mathbf{r}} = m\sin\varphi \\
   & \hat{\mathbf{z}} = \sin\varphi\,\hat{\mathbf{r}} + \cos\varphi\,\hat{\boldsymbol{\varphi}}
\end{cases}
  \Longrightarrow
\mathbf{B}(r,\varphi) = \frac{\mu_0}{4\pi}\frac{m}{r^3}
\left(2\sin\varphi\,\hat{\mathbf{r}} - \cos\varphi\,\hat{\boldsymbol{\varphi}}\right) 
$$

$$
\Longrightarrow
B_r = \dfrac{\mu_0}{4\pi}\dfrac{2m\sin\varphi}{r^3}, \quad
B_{\varphi} = -\dfrac{\mu_0}{4\pi}\dfrac{m\cos\varphi}{r^3}
$$

## 2. 极点与赤道：2 倍关系

$$
B_{\text{极}} = \frac{\mu_0}{4\pi}\frac{2|m|}{r^3}\ (\varphi=\pm90^\circ,\ \text{纯径向}),
\quad
B_{\text{赤}} = \frac{\mu_0}{4\pi}\frac{|m|}{r^3}\ (\varphi=0^\circ,\ \text{纯水平})
\;\Longrightarrow\;
\boxed{B_{\text{极}} = 2\,B_{\text{赤}}}
$$

## 3. 数值（IGRF-1990，$g_1^0 = -29775$ nT）

$$
m = \frac{4\pi}{\mu_0}a^3 g_1^0
\;\Longrightarrow\;
B_{\text{赤}} = |g_1^0| \approx 29\,775\ \text{nT},
\qquad
B_{\text{极}} = 2|g_1^0| \approx 59\,550\ \text{nT}
$$

## 4. 场强随纬度的变化

$$
|\mathbf{B}(\varphi)| = \sqrt{B_r^2 + B_{\varphi}^2}
= \frac{\mu_0}{4\pi}\frac{|m|}{a^3}\sqrt{1 + 3\sin^2\varphi}
= B_{\text{赤}}\sqrt{1 + 3\sin^2\varphi},
\qquad B_{\text{赤}} \le |\mathbf{B}| \le 2\,B_{\text{赤}}
$$


$$
|\mathbf{B}(\varphi)| = B_{\text{赤}}\sqrt{1 + 3\sin^2\varphi}
\;\Longrightarrow\;
B_{\text{赤}}\ (\varphi=0^\circ) \;\le\; |\mathbf{B}| \;\le\; 2B_{\text{赤}}\ (\varphi=\pm90^\circ)
$$

## 5. 磁倾角与纬度：偶极子"纬度定理"

$$
\tan I = \frac{Z}{H} = \frac{-B_r}{B_{\varphi}} = 2\tan\varphi
\; \Longrightarrow \;
\varphi = \arctan\!\Big(\tfrac{1}{2}\tan I\Big)
$$

特例：

$$
\varphi = 0^\circ:\ I = 0^\circ; \qquad
\varphi = \pm 45^\circ:\ I \approx \pm 63.4^\circ; \qquad
\varphi = \pm 90^\circ:\ I = \pm 90^\circ
$$

该关系是**古地磁学重建古纬度**的基础：测得岩石中的磁倾角 $I$，即可反演古纬度

$$
\varphi = \arctan\!\left(\tfrac{1}{2}\tan I\right)
$$