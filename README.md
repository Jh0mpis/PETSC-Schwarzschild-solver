# Solving Schwarzschild using PETSC c library



# Schwarzschild Geodesics

The Schwarzschild line element in Boyer Lindquist is given by

\begin{equation}
    ds^2 = g_{\mu\nu}dx^\mu dx^\nu = -f(r)dt^2+\dfrac{dr^2}{f(r)}+r^2d\Omega^2,
\end{equation}

where $d\Omega^2 = d\theta^2+\sin^2\theta d\phi^2$, $g_{\mu\nu}$ is the Schwarzschild metric and $t\in (-\infty, \infty), \,r \in [0,\infty),\, \theta \in [0,\pi], \, \phi\in[0, 2\pi]$ are the manifold coordinates, moreover,

\begin{gather}
    f(r) = 1-\dfrac{2M}{r},
    \label{eq:f}
\end{gather}

where $M$ is the black Hole's mass and the event horizon is located at $r=2M$. The metric components come as follow

\begin{equation}
    \begin{split}
        g_{tt} & = -f(r),\\
        g_{rr} & = \dfrac{1}{f(r)},\\
        g_{\theta\theta} & = r^2, \\        
        g_{\phi\phi} & = r^2\sin^2\theta. \\
    \end{split}
\end{equation}

The Langrangian for a free particle in the Schwarzschild space-time is

\begin{equation}
    \begin{split}
        \mathcal{L} & = \dfrac{1}{2}g_{\mu\nu}\dot{x}^\mu\dot{x}^\nu\\
        & = \dfrac{1}{2}\left(g_{tt}\dot{t}^2+g_{rr}\dot{r}^2+g_{\theta\theta}\dot{\theta}^2+g_{\phi\phi}\dot{\phi}^2\right),
    \end{split}
\end{equation}

where $\dot{x}^\mu = \dfrac{dx^\mu}{d\tau}$ for massive particles, where $\tau$ is the proper time (the time measured by the particle). For the massless particles we have $\dot{x}^\mu = \dfrac{dx^\mu}{d\lambda}$, where $\lambda$ is the affine parameter. Furthermore, the 4-velocity constrain can be expressed as follows

\begin{equation}
    g_{\mu\nu}\dot{x}^\mu\dot{x}^\nu = -k,
\end{equation}

where $k = 1$ for massive particles and $k = 0$ for massless particles (photons).\\

The Euler-Lagrange equations for the lagrangian can be written as follows

\begin{equation}
    \begin{split}
        \dfrac{d}{\tau}\partialD{\mathcal{L}}{\dot{x}^\mu} - \partialD{\mathcal{L}}{x^\mu} = 0, \, \text{(Massive),}\\
        \dfrac{d}{\lambda}\partialD{\mathcal{L}}{\dot{x}^\mu} - \partialD{\mathcal{L}}{x^\mu} = 0, \, \text{(Massless).}
    \end{split}
\end{equation}

Note that, $\partialD{\mathcal{L}}{t} = \partialD{\mathcal{L}}{\phi} = 0$, then $\partialD{\mathcal{L}}{\dot{t}} = - E$ and $\partialD{\mathcal{L}}{\dot{\phi}} = L_z$ are constants and are associated to the energy density of the particle and the angular moment, at the $z$, axis of the particle respectively, so we will obtain

\begin{equation}
    \begin{split}
        \partialD{\mathcal{L}}{\dot{t}}& = - E   = g_{tt}\dot{t}\\
        & \text{then,}\quad \dot{t} = -\dfrac{E}{g_{tt}},
    \end{split}
\end{equation}

and

\begin{equation}
    \begin{split}
        \partialD{\mathcal{L}}{\dot{\phi}}& = L_z   = g_{\phi\phi}\dot{\phi}\\
        & \text{then,}\quad \dot{\phi} = \dfrac{L_z}{g_{\phi\phi}}.
    \end{split}
\end{equation}
