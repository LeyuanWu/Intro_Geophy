# Derivation of the Magnetic Dipole Field

## 1. Starting Point: The Biot-Savart Law
The magnetic field $\mathbf{B}$ produced by a steady current $I$ flowing through a closed loop is given by the Biot-Savart law:

$$
\mathbf{B} = \frac{\mu_0 I}{4\pi} \oint d\mathbf{l}' \times \frac{\mathbf{r} - \mathbf{r}'}{|\mathbf{r} - \mathbf{r}'|^3}
$$

where $\mathbf{r}$ is the observation point, $\mathbf{r}'$ is the source point on the loop, and $d\mathbf{l}'$ is the infinitesimal line element.

## 2. Far-Field Approximation ($r \gg r'$)
For an observation point far from the loop ($r \gg r'$), we expand the geometric term using a Taylor series. To first order in $r'/r$, the approximation is:

$$
\frac{\mathbf{r} - \mathbf{r}'}{|\mathbf{r} - \mathbf{r}'|^3} \approx \frac{\mathbf{r} - \mathbf{r}' + 3(\hat{\mathbf{r}} \cdot \mathbf{r}')\hat{\mathbf{r}}}{r^3}
$$

Substituting this back into the integral and distributing the cross product $d\mathbf{l}' \times (\dots)$, we can separate the expression into three distinct terms:

$$
\mathbf{B} \approx \frac{\mu_0 I}{4\pi r^3} \left[ \underbrace{\oint d\mathbf{l}' \times \mathbf{r}}_{\text{Term 1}} + \underbrace{\oint d\mathbf{l}' \times (-\mathbf{r}')}_{\text{Term 2}} + \underbrace{\oint 3(\hat{\mathbf{r}} \cdot \mathbf{r}')(d\mathbf{l}' \times \hat{\mathbf{r}})}_{\text{Term 3}} \right]
$$

## 3. Evaluating the Three Terms
We evaluate each term using standard vector calculus identities for closed loops, where $\mathbf{S}$ represents the vector area of the loop ($\mathbf{S} = \int d\mathbf{a}'$).

### Term 1: Monopole Term
Since $\mathbf{r}$ is constant with respect to the integration variable $\mathbf{r}'$, it can be factored out:

$$
\text{Term 1} = \left( \oint d\mathbf{l}' \right) \times \mathbf{r}
$$

For any closed loop, the integral of the line element is zero ($\oint d\mathbf{l}' = \mathbf{0}$). Thus:

$$
\text{Term 1} = \mathbf{0}
$$

*(This reflects the absence of magnetic monopoles.)*

### Term 2: Dipole Term (Part A)
Using the standard closed-loop identity relating the line integral to the vector area $\mathbf{S}$:

$$
\text{Term 2} = \oint d\mathbf{l}' \times (-\mathbf{r}') = 2\mathbf{S}
$$

### Term 3: Dipole Term (Part B)
Since $\hat{\mathbf{r}}$ is constant with respect to the integration, we can factor it out of the cross product and apply the closed-loop identity for the scalar-vector product:

$$
\text{Term 3} = 3 \left( \oint (\hat{\mathbf{r}} \cdot \mathbf{r}') \, d\mathbf{l}' \right) \times \hat{\mathbf{r}}
$$

Evaluating this using the vector area $\mathbf{S}$ yields:

$$
\text{Term 3} = 3(\mathbf{S} \cdot \hat{\mathbf{r}})\hat{\mathbf{r}} - 3\mathbf{S}
$$

## 4. Final Assembly
Substitute the evaluated terms back into the expanded Biot-Savart equation:

$$
\mathbf{B} = \frac{\mu_0 I}{4\pi r^3} \Big[ \mathbf{0} + 2\mathbf{S} + 3(\mathbf{S} \cdot \hat{\mathbf{r}})\hat{\mathbf{r}} - 3\mathbf{S} \Big]
$$

Combine the $\mathbf{S}$ terms ($2\mathbf{S} - 3\mathbf{S} = -\mathbf{S}$):

$$
\mathbf{B} = \frac{\mu_0 I}{4\pi r^3} \Big[ 3(\mathbf{S} \cdot \hat{\mathbf{r}})\hat{\mathbf{r}} - \mathbf{S} \Big]
$$

Finally, define the magnetic dipole moment as $\mathbf{m} = I\mathbf{S}$. Substituting $\mathbf{m}$ into the equation gives the standard form of the magnetic dipole field:

$$
\mathbf{B} = \frac{\mu_0}{4\pi r^3} \Big[ 3(\mathbf{m} \cdot \hat{\mathbf{r}})\hat{\mathbf{r}} - \mathbf{m} \Big]
$$

---

## Appendix: Detailed Derivation of Term 3

This appendix provides a rigorous proof of the identity used in evaluating Term 3:

$$
\oint (\hat{\mathbf{r}} \cdot \mathbf{r}') \, d\mathbf{l}' \times \hat{\mathbf{r}} = (\mathbf{S} \cdot \hat{\mathbf{r}})\hat{\mathbf{r}} - \mathbf{S}
$$

### Step A: Simplify Using Constant Vectors
Let $\mathbf{c} = \hat{\mathbf{r}}$. Since $\mathbf{c}$ is constant over the closed-loop integration, we can rearrange the cross product:

$$
\oint (\mathbf{c} \cdot \mathbf{r}') \, d\mathbf{l}' \times \mathbf{c} = \left( \oint (\mathbf{c} \cdot \mathbf{r}') \, d\mathbf{l}' \right) \times \mathbf{c}
$$

### Step B: Apply the Closed-Loop Vector Identity
We use the fundamental identity that relates a scalar-weighted line integral over a closed loop to the vector area $\mathbf{S}$:

$$
\oint f(\mathbf{r}') \, d\mathbf{l}' = \int_S d\mathbf{a}' \times \nabla' f(\mathbf{r}')
$$

Setting $f(\mathbf{r}') = \mathbf{c} \cdot \mathbf{r}'$, we compute the gradient:

$$
\nabla' (\mathbf{c} \cdot \mathbf{r}') = \mathbf{c}
$$

Therefore:

$$
\oint (\mathbf{c} \cdot \mathbf{r}') \, d\mathbf{l}' = \int_S d\mathbf{a}' \times \mathbf{c} = \left( \int_S d\mathbf{a}' \right) \times \mathbf{c} = \mathbf{S} \times \mathbf{c}
$$

### Step C: Evaluate the Outer Cross Product
Substitute the result from Step B back into the expression from Step A:

$$
\left( \oint (\mathbf{c} \cdot \mathbf{r}') \, d\mathbf{l}' \right) \times \mathbf{c} = (\mathbf{S} \times \mathbf{c}) \times \mathbf{c}
$$

### Step D: Apply the Vector Triple Product
Use the BAC-CAB rule $(\mathbf{A} \times \mathbf{B}) \times \mathbf{C} = (\mathbf{A} \cdot \mathbf{C})\mathbf{B} - (\mathbf{B} \cdot \mathbf{C})\mathbf{A}$:

$$
(\mathbf{S} \times \mathbf{c}) \times \mathbf{c} = (\mathbf{S} \cdot \mathbf{c})\mathbf{c} - (\mathbf{c} \cdot \mathbf{c})\mathbf{S}
$$

Since $\mathbf{c} = \hat{\mathbf{r}}$ is a unit vector, $\mathbf{c} \cdot \mathbf{c} = 1$:

$$
(\mathbf{S} \times \hat{\mathbf{r}}) \times \hat{\mathbf{r}} = (\mathbf{S} \cdot \hat{\mathbf{r}})\hat{\mathbf{r}} - \mathbf{S}
$$

### Step E: Include the Factor of 3
Multiplying by the prefactor of 3 from the original Term 3 expression confirms the result stated in Section 3:

$$
\oint 3(\hat{\mathbf{r}} \cdot \mathbf{r}')(d\mathbf{l}' \times \hat{\mathbf{r}}) = 3(\mathbf{S} \cdot \hat{\mathbf{r}})\hat{\mathbf{r}} - 3\mathbf{S}
$$