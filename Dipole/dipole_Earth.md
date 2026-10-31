## PPT-ready expressions

**Matching the dipole potential to the Gauss coefficients**

$$\underbrace{\frac{\mu_0}{4\pi}\frac{\mathbf{m}\cdot\hat{\mathbf{r}}}{r^2}}_{\text{dipole}} = \underbrace{\frac{a^3}{r^2}\left[ g_1^0\,P_1^0(\theta) + \left( g_1^1\cos\lambda + h_1^1\sin\lambda \right)P_1^1(\theta) \right]}_{\text{Gauss expansion}}$$

$$\frac{\mu_0}{4\pi r^2}\Big[\,m_z\cos\theta + m_x\sin\theta\cos\lambda + m_y\sin\theta\sin\lambda\,\Big] \;\Longleftrightarrow\; \frac{a^3}{r^2}\Big[\,g_1^0\cos\theta + g_1^1\sin\theta\cos\lambda + h_1^1\sin\theta\sin\lambda\,\Big]$$

$$\boxed{\;m_x = \frac{4\pi a^3}{\mu_0}\,g_1^1, \qquad m_y = \frac{4\pi a^3}{\mu_0}\,h_1^1, \qquad m_z = \frac{4\pi a^3}{\mu_0}\,g_1^0\;}$$

**Total dipole moment (magnitude)**

$$\;m = \frac{4\pi a^3}{\mu_0}\sqrt{\left(g_1^0\right)^2 + \left(g_1^1\right)^2 + \left(h_1^1\right)^2}\;$$

**Axial dipole moment** (component along rotation axis, $g_1^1=h_1^1=0$)

$$\;m_{\text{axial}} = \frac{4\pi a^3}{\mu_0}\,\vert g_1^0 \vert\;$$

**Dipole axis orientation** (geomagnetic pole position)

$$\tan\lambda_0 = \frac{h_1^1}{g_1^1}, \quad \tan\theta_0 = \frac{\sqrt{(g_1^1)^2+(h_1^1)^2}}{g_1^0}$$



### Appendix

#### 1. Spherical harmonic potential, internal sources

$$V = a \sum_{n=1}^{\infty} \left(\frac{a}{r}\right)^{n+1} \sum_{m=0}^{n} \left( g_n^m \cos m\lambda + h_n^m \sin m\lambda \right) P_n^m(\theta)$$

where $a$ is Earth's radius, $\theta$ colatitude, $\lambda$ longitude. The sum starts at $n=1$ because $g_0^0$ would describe a magnetic monopole and must vanish.

#### 2. Truncate at $n = 1$: the dipole contribution

$$V^{D} = \frac{a^3}{r^2} \left[ g_1^0 P_1^0(\theta) + \left( g_1^1 \cos\lambda + h_1^1 \sin\lambda \right) P_1^1(\theta) \right]$$

(Each harmonic coefficient is optimal in a least-squares sense regardless of how many terms are retained, so dropping $n>1$ does not bias $g_1^0, g_1^1, h_1^1$.)

#### 3. Substitute the degree-1 Legendre functions

With $P_1^0(\theta) = \cos\theta$ and $P_1^1(\theta) = \sin\theta$:

$$V^{D} = \frac{a^3}{r^2} \left[ g_1^0 \cos\theta + g_1^1 \cos\lambda \sin\theta + h_1^1 \sin\lambda \sin\theta \right]$$

#### 4. Convert to Cartesian coordinates

With origin at Earth's center, $z$ along the spin axis (north), $x$ on the Greenwich meridian:

$$x = r \sin\theta \cos\lambda, \qquad y = r \sin\theta \sin\lambda, \qquad z = r \cos\theta$$

so that

$$V^{D} = \frac{a^3}{r^2} \left[ g_1^0 \frac{z}{r} + g_1^1 \frac{x}{r} + h_1^1 \frac{y}{r} \right]$$

#### 5. Match against the centered-dipole potential

The potential of a dipole $\mathbf{m}$ at the origin is

$$V = C_{\mathrm{m}} \frac{\mathbf{m} \cdot \hat{\mathbf{r}}}{r^2} = \frac{C_{\mathrm{m}}}{r^2} \left[ m_x \frac{x}{r} + m_y \frac{y}{r} + m_z \frac{z}{r} \right]$$

with $C_{\mathrm{m}} = 1$ (dimensionless) in emu, or $C_{\mathrm{m}} = \mu_0/4\pi = 10^{-7}\ \mathrm{H\,m^{-1}}$ in SI.

**Conclusion:** the $n=1$ term of the expansion *is* the field of a geocentric dipole. Comparing coefficients of $x/r,\ y/r,\ z/r$ at $r = a$ in SI units:

$$m_x = \frac{4\pi}{\mu_0} a^3 g_1^1, \qquad m_y = \frac{4\pi}{\mu_0} a^3 h_1^1, \qquad m_z = \frac{4\pi}{\mu_0} a^3 g_1^0$$

#### 6. Numerical evaluation (IGRF 1990, Table 8.1)

With $a = 6.371\times10^{6}$ m, $g_1^0 = -29775$ nT, $g_1^1 = -1851$ nT, $h_1^1 = 5411$ nT:

$$m_x = -0.479\times10^{22}, \qquad m_y = 1.399\times10^{22}, \qquad m_z = -7.700\times10^{22} \quad (\mathrm{A\,m^2})$$

**Magnitude:**

$$m = \sqrt{m_x^2 + m_y^2 + m_z^2} = \frac{4\pi}{\mu_0} a^3 \sqrt{(g_1^0)^2 + (g_1^1)^2 + (h_1^1)^2} = 7.840\times10^{22}\ \mathrm{A\,m^2}$$

**Orientation** (colatitude and east longitude of the dipole axis):

$$\theta = \arccos\frac{m_z}{m} = 169.1^\circ, \qquad \lambda = \arccos\frac{m_x}{\sqrt{m_x^2 + m_y^2}} = 108.9^\circ$$

Hence the **geomagnetic poles**: the positive end of the dipole pierces the surface at $(\theta,\lambda) = (169.1^\circ, 108.9^\circ)$ → 79.1°S, 108.9°E; the negative end at $(10.9^\circ, 288.9^\circ)$ → 79.1°N, 71.1°W (Kane Basin). Because the field is not perfectly dipolar, these differ from the *magnetic* (dip) poles — per IGRF 1990 the north magnetic pole lies ~700 km from the geomagnetic pole. Relaxing the geocentric constraint gives the **eccentric dipole**, offset ~400 km along $+z$.

*(Check: substituting the Table 8.1 coefficients into the formulas above reproduces all of the book's reported values, e.g. $\tfrac{4\pi}{\mu_0}a^3 = 2.586\times10^{27}$ SI units gives $m_z \approx -7.70\times10^{22}$ A·m².)*

#### 7. The nondipole field

Removing the $n=1$ term from Eq. 8.16:

$$V^{N} = a \sum_{n=2}^{\infty} \left(\frac{a}{r}\right)^{n+1} \sum_{m=0}^{n} \left( g_n^m \cos m\lambda + h_n^m \sin m\lambda \right) P_n^m(\theta)$$

Differentiate the full expansion (8.16) with respect to $r$:

$$B_r = -\frac{\partial V}{\partial r} = \sum_{n=1}^{\infty} (n+1) \left(\frac{a}{r}\right)^{n+2} \sum_{m=0}^{n} \left( g_n^m \cos m\lambda + h_n^m \sin m\lambda \right) P_n^m(\theta)$$

At the **north geographic pole** ($r = a$, $\theta = 0^\circ$), since $P_n^m(0) = 1$ for $m = 0$ and $0$ otherwise:

$$B_r = \sum_{n=1}^{\infty} (n+1)\, g_n^0$$

Splitting dipole vs. nondipole parts, using Table 8.1:

$$B_r^{D} = 2 g_1^0 = -59{,}550\ \mathrm{nT}, \qquad B_r^{N} \approx \sum_{n=2}^{10} (n+1)\, g_n^0 \approx 3{,}518\ \mathrm{nT}$$

$$100 \times \frac{|B_r^{N}|}{|B_r^{D} + B_r^{N}|} = 6.3\ \text{percent}$$

Although computed for one component at one location, this is representative: the nondipole field contributes only ~10% of the main field, so the geomagnetic field is dipolar to a good first approximation. The text cautions, however, against attaching physical meaning to the dipole/nondipole split — the low-degree field originates in the fluid outer core, not at the geocenter, and dipole-plus-radial-dipole models (e.g., Alldredge and Hurwitz, cited therein) can fit the same field with a centered dipole of about twice the spherical-harmonic strength.