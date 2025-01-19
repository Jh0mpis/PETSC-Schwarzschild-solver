# Solving Schwarzschild using PETSC c library



## Schwarzschild Geodesics

The Schwarzschild line element in Boyer Lindquist is given by

$$ds^2 = g_{\mu\nu}dx^\mu dx^\nu = -f(r)dt^2+\dfrac{dr^2}{f(r)}+r^2d\Omega^2,$$

where $d\Omega^2 = d\theta^2+\sin^2\theta d\phi^2$, $g_{\mu\nu}$ is the Schwarzschild metric and $t\in (-\infty, \infty), r \in [0,\infty), \theta \in [0,\pi], \phi\in[0, 2\pi]$ are the manifold coordinates, moreover,

$$f(r) = 1-\dfrac{2M}{r},$$

where $M$ is the black Hole's mass and the event horizon is located at $r=2M$. The metric components come as follow

$$g_{tt}  = -f(r),$$
$$g_{rr}  = \dfrac{1}{f(r)},$$
$$g_{\theta\theta}  = r^3,$$ 
$$g_{\phi\phi}  = r^2\sin^2\theta.$$

The Langrangian for a free particle in the Schwarzschild space-time is

$$\mathcal{L}  = \dfrac{1}{2}g_{\mu\nu}\dot{x}^\mu\dot{x}^\nu$$
$$ = \dfrac{1}{2}\left(g_{tt}\dot{t}^2+g_{rr}\dot{r}^2+g_{\theta\theta}\dot{\theta}^2+g_{\phi\phi}\dot{\phi}^2\right),$$

where $\dot{x}^\mu = \dfrac{dx^\mu}{d\tau}$ for massive particles, where $\tau$ is the proper time (the time measured by the particle). For the massless particles we have $\dot{x}^\mu = \dfrac{dx^\mu}{d\lambda}$, where $\lambda$ is the affine parameter. Furthermore, the 4-velocity constrain can be expressed as follows

$$g_{\mu\nu}\dot{x}^\mu\dot{x}^\nu = -k,$$

where $k = 1$ for massive particles and $k = 0$ for massless particles (photons).

The Euler-Lagrange equations for the lagrangian can be written as follows

$$\dfrac{d}{\tau}\dfrac{\partial\mathcal{L}}{\partial\dot{x}^\mu} - \dfrac{\partial\mathcal{L}}{\partial x^\mu} = 0,  \text{(Massive),}$$
$$\dfrac{d}{\lambda}\dfrac{\partial\mathcal{L}}{\partial\dot{x}^\mu} - \dfrac{\partial\mathcal{L}}{\partial x^\mu} = 0,  \text{(Massless).}$$

Note that, $\dfrac{\partial\mathcal{L}}{\partial t} = \dfrac{\partial\mathcal{L}}{\partial\phi} = 0$, then $\dfrac{\partial\mathcal{L}}{\partial\dot{t}} = - E$ and $\dfrac{\partial\mathcal{L}}{\partial\dot{\phi}} = L_z$ are constants and are associated to the energy density of the particle and the angular moment, at the $z$, axis of the particle respectively, so we will obtain

$$\dfrac{\partial\mathcal{L}}{\partial\dot{t}} = - E   = g_{tt}\dot{t}$$

then,

$$\dot{t} = -\dfrac{E}{g_{tt}},$$

and

$$\dfrac{\partial\mathcal{L}}{\dot{\partial\phi}} = L_z  = g_{\phi\phi}\dot{\phi}$$

then,

$$\dot{\phi} = \dfrac{L_z}{g_{\phi\phi}}.$$

## Geodesics in the equatorial plane

Using the 4-velocity constrain we can see that

$$k = g_{tt}\dot{t}^2 + g_{rr}\dot{r}^2 + g_{\theta\theta}\dot{\theta}^2 + g_{\phi\phi}\dot{\phi}^2$$
$$ = \dfrac{E}{g_{tt}}+g_{rr}\dot{r}^2 + g_{\theta\theta}\dot{\theta}^2 + \dfrac{L_z^2}{g_{\theta\theta}},$$

then,

$$g_{rr}\dot{r}^2 + g_{\theta\theta}\dot{\theta}^2 = k - \dfrac{E}{g_{tt}}-\dfrac{L_z^2}{g_{\theta\theta}},$$

at the equatorial plane we have that $\theta = \dfrac{\pi}{2}$, $\dot{\theta} = 0$ and $\ddot{\theta} = 0$, so we fix the movement at the equatorial plane, using some algebra the energy equation becomes

$$\dfrac{1}{2}\dot{r}^2 = k\dfrac{1}{2g_{rr}} - \dfrac{E}{2g_{tt}g_{rr}}-\dfrac{L_z^2}{2g_{\theta\theta}g_{rr}},$$

recovering the values of each $g_{\mu\nu}$ and the $f(r)$ function we have

$$\dfrac{1}{2}\dot{r}^2 = k\dfrac{1}{2}f(r) + \dfrac{Ef(r)}{2f(r)}-\dfrac{L_z^2f(r)}{2r^2\sin^2(\pi/2)}$$
$$ = \dfrac{E+k}{2}-\dfrac{kM}{r}-\dfrac{L_z^2}{2r^2}+\dfrac{L_zM}{r^3}$$
$$ = \dfrac{E+k}{2} - U_{\text{eff}},$$

where, $U_{\text{eff}}$ is the effective potential measured by an Minkowski observer, is the potential that particle experiments produced by the spacetime curvature. So, the effective potential is

$$U_\text{eff} = \dfrac{kM}{r}+\dfrac{L_z^2}{2r^2}-\dfrac{L_zM}{r^3},$$

the first two terms are from the Newtonian gravity, the first one is the usual gravitational term (with $k=-1$ for massive particles) and the second one is a result of the repulsion force of a particle with an angular momentum. The last one term is the correction given by Einstein's General Relativity. By definition, we have that

$$\vec{f}  = -\vec{\nabla} U_\text{eff}$$
$$ f_r\hat{r} = -\dfrac{\partial U_\text{eff}}{\partial r}\hat{r},$$
$$\Rightarrow \ddot{r}  = \dfrac{kM}{r^2}+\dfrac{L_z^2}{r^3}-3\dfrac{L_zM}{r^4}$$

then, the dynamics of a massive particle in the equatorial plane at Schwarzschild metric are given as follows

$$\dot{t}  = -\dfrac{Er}{r-2M}, $$
$$\ddot{r}  = -\dfrac{M}{r^2}+\dfrac{L_z^2}{r^3}-3\dfrac{L_zM}{r^4},$$
$$\dot{\phi}  = \dfrac{L_z}{r^2}.$$

While for photons the dynamics are given by

$$\dot{t} = -\dfrac{Er}{r-2M}$$
$$\ddot{r} = \dfrac{L_z^2}{r^3}-3\dfrac{L_zM}{r^4},$$
$$\dot{\phi} = \dfrac{L_z}{r^2}.$$

