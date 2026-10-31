# Derivation of the Magnetic Dipole Field


## 1. Analogy: Gravitational and Magnetic Monopoles

We begin by establishing the analogy between the gravitational field of a point mass and the hypothetical magnetic field of a point magnetic charge (monopole).

For a point mass $m$, the gravitational field $\mathbf{g}$ and its scalar potential $U$ (where $\mathbf{g} = \nabla U$) are:

$$
\text{Gravity:} \quad \mathbf{g} = -G\frac{m}{r^2}\mathbf{\hat{r}}, \quad \mathbf{g} = \nabla U \quad \Longleftrightarrow \quad U = G\frac{m}{r}
$$

By analogy, for a point magnetic charge $q$, the magnetic field $\mathbf{B}$ and its scalar potential $V$ (where $\mathbf{B} = -\nabla V$) are:

$$
\text{Magnetic:} \quad \mathbf{B} = \frac{\mu_0}{4\pi}\frac{q}{r^2}\mathbf{\hat{r}}, \quad \mathbf{B} = -\nabla V \quad \Longleftrightarrow \quad V = \frac{\mu_0}{4\pi}\frac{q}{r}
$$

## 2. Scalar Potential of a Magnetic Dipole

A magnetic dipole can be modeled as a pair of opposite magnetic charges, $+p$ and $-p$, separated by a small displacement vector $\mathbf{d}$. By the principle of superposition, the total scalar potential at a position $\mathbf{r}$ is:

$$
V = \frac{\mu_0 p}{4\pi}\left(\frac{1}{r_{+}} - \frac{1}{r_{-}}\right) = \frac{\mu_0 p}{4\pi}\left(\frac{r_{-} - r_{+}}{r_{+}r_{-}}\right)
$$

where $r_+$ is the distance to the positive charge and $r_-$ is the distance to the negative charge.

### Far-Field Approximation ($r \gg d$)

Let the charges be located at $\mp \frac{\mathbf{d}}{2}$. We can express the distances $r_-$ and $r_+$ as:

$$
\begin{cases}
r_{-} = \sqrt{\left(\mathbf{r}+\frac{\mathbf{d}}{2}\right) \cdot \left(\mathbf{r}+\frac{\mathbf{d}}{2}\right)} \approx r + \frac{1}{2}\hat{\mathbf{r}} \cdot \mathbf{d} \\[10pt]
r_{+} = \sqrt{\left(\mathbf{r}-\frac{\mathbf{d}}{2}\right) \cdot \left(\mathbf{r}-\frac{\mathbf{d}}{2}\right)} \approx r - \frac{1}{2}\hat{\mathbf{r}} \cdot \mathbf{d}
\end{cases}
$$

Using these approximations, we evaluate the terms in the potential equation:

1.  **Numerator:** $r_{-} - r_{+} \approx \left(r + \frac{1}{2}\hat{\mathbf{r}} \cdot \mathbf{d}\right) - \left(r - \frac{1}{2}\hat{\mathbf{r}} \cdot \mathbf{d}\right) = \mathbf{d} \cdot \hat{\mathbf{r}}$
2.  **Denominator:** $r_{+}r_{-} \approx r^2$ (ignoring higher-order terms of $d/r$)

Substituting these back into the potential expression yields:

$$
V = \frac{\mu_0 p}{4\pi}\frac{\mathbf{d} \cdot \hat{\mathbf{r}}}{r^2}
$$

Defining the **magnetic dipole moment** as $\mathbf{m} = p\mathbf{d}$, we obtain the standard dipole scalar potential:

$$
V = \frac{\mu_0}{4\pi}\frac{\mathbf{m} \cdot \hat{\mathbf{r}}}{r^2}
$$

## 3. Deriving the Magnetic Field $\mathbf{B}$

To find the magnetic field, we use the relation $\mathbf{B} = -\nabla V$. First, rewrite the potential using $\hat{\mathbf{r}} = \frac{\mathbf{r}}{r}$ to make differentiation easier:

$$
V = \frac{\mu_0}{4\pi} \frac{\mathbf{m} \cdot \mathbf{r}}{r^3}
$$

Now, apply the gradient operator:

$$
\mathbf{B} = -\frac{\mu_0}{4\pi} \nabla \left( \frac{\mathbf{m} \cdot \mathbf{r}}{r^3} \right)
$$

### Applying the Product Rule

Let $f = \mathbf{m} \cdot \mathbf{r}$ and $g = r^{-3}$. The product rule for gradients states $\nabla(fg) = g\,\nabla f + f\,\nabla g$.

**Term 1: Gradient of $\mathbf{m} \cdot \mathbf{r}$**
Since $\mathbf{m}$ is a constant vector:

$$
\nabla(\mathbf{m} \cdot \mathbf{r}) = \mathbf{m}
$$

**Term 2: Gradient of $r^{-3}$**
Using the chain rule and $\nabla r = \hat{\mathbf{r}} = \frac{\mathbf{r}}{r}$:

$$
\nabla(r^{-3}) = \frac{d}{dr}(r^{-3}) \nabla r = -3r^{-4} \left(\frac{\mathbf{r}}{r}\right) = -\frac{3\mathbf{r}}{r^5}
$$

### Combining the Terms

Substitute the two parts back into the product rule:

$$
\nabla \left( \frac{\mathbf{m} \cdot \mathbf{r}}{r^3} \right) = \frac{1}{r^3}\mathbf{m} + (\mathbf{m} \cdot \mathbf{r})\left(-\frac{3\mathbf{r}}{r^5}\right) = \frac{\mathbf{m}}{r^3} - \frac{3(\mathbf{m} \cdot \mathbf{r})\mathbf{r}}{r^5}
$$

Now, apply $\mathbf{B} = -\nabla V$:

$$
\mathbf{B} = -\frac{\mu_0}{4\pi} \left[ \frac{\mathbf{m}}{r^3} - \frac{3(\mathbf{m} \cdot \mathbf{r})\mathbf{r}}{r^5} \right] = \frac{\mu_0}{4\pi} \left[ \frac{3(\mathbf{m} \cdot \mathbf{r})\mathbf{r}}{r^5} - \frac{\mathbf{m}}{r^3} \right]
$$

## 4. Final Result

To express the result in terms of the unit vector $\hat{\mathbf{r}}$, we substitute $\mathbf{r} = r\hat{\mathbf{r}}$ and $\mathbf{m} \cdot \mathbf{r} = r(\mathbf{m} \cdot \hat{\mathbf{r}})$, then factor out $\frac{1}{r^3}$:

$$
\mathbf{B} = \frac{\mu_0}{4\pi} \left[ \frac{3(\mathbf{m} \cdot r\hat{\mathbf{r}})(r\hat{\mathbf{r}})}{r^5} - \frac{\mathbf{m}}{r^3} \right] = \frac{\mu_0}{4\pi} \left[ \frac{3(\mathbf{m} \cdot \hat{\mathbf{r}})\hat{\mathbf{r}}}{r^3} - \frac{\mathbf{m}}{r^3} \right]
$$

Factoring out $\frac{1}{r^3}$ gives the final, classic expression for the magnetic dipole field:

$$
\boxed{\mathbf{B} = \frac{\mu_0}{4\pi} \left[ \frac{3(\mathbf{m} \cdot \hat{\mathbf{r}})\hat{\mathbf{r}} - \mathbf{m}}{r^3} \right]}
$$