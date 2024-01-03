# System Description and Analysis

This notebook models a simple spring loaded weapon system with two arm segments, one with the spring mounted to it and the other with the weapon mounted to it. The two arm segments have a fixed angle between them, and will have separate linear densities that, for the purposes of this model, are assumed to be constant. The spring arm has length $l_s$ and a linear density of $\lambda_s$. The weapon arm has a length of $l_w$ and a linear density of $\lambda_w$. The weapon itself has a mass of $m$ and a radius of $r$. The weapon is modeled as a disk with a constant density, a radius $r$ and a mass of $m$. Additionally, since the forces from the spring will be so high, the effects of gravity can be safely neglected.

In operation, the weapon system will be loaded by compressing a gas spring with a full length of $l_0$ and a final length of $l_f$. The gas spring will have a force on contact of $F_0$, a fully compressed force of $F_f$, and an initial pressure of $\rho_0$. For simplicity, the system will be arranged so that the spring is at rest at $\theta_0 = 0$ and fully loaded at $\theta_f = \pi$. After the weapon is loaded, the spring arm will be released, causing the spring to push the weapon forward and transforming the potential energy stored in the spring into kinetic energy. For the sake of this analysis, it will be assumed that the compression and release of the gas spring is a roughly adiabatic process. 

The goal of this model is to analyze the performance of various different commercially available gas springs, and find the overall system characteristics for each of them. System characteristics that will be analyzed include, but are not limited to, the maximum tip velocity, the maximum torque required to load the weapon, the maximum radial force experienced by the spring arm, the maximum angular force experienced by the spring arm, the efficiency of energy storage per unit of spring mass, and the efficiency of tip velocity per unit of spring mass.

## Calculating Cylinder Pressure and Volume

From the gas spring specifications, we know the initial pressure $\rho_0$ and the initial force $F_0$. From this information we can find the effective cross sectional area of the plunger perpendicular to the direction of movement $A$ as shown in Equation A.1.

\begin{equation}
\tag{A.1}
A = \frac{F_0}{\rho_0}
\end{equation}

Using this, we can find the pressure of the gas spring at the end of the stroke $\rho_f$ as shown in Equation A.2.

\begin{equation}
\tag{A.2}
\rho_f = \frac{F_f}{A}
\end{equation}

Next, we find the change in volume of the gas spring throughout the entire stroke $\Delta V$ as shown in Equation A.3.

\begin{equation}
\tag{A.3}
\Delta V = (l_f - l_0) A
\end{equation}

Since we know that $\rho_0 V_0 \approx \rho_f V_f$ and $V_f = V_0 + \Delta V$, we can find the initial volume of the gas spring $V_f$ as shown in Equation A.4.

\begin{equation}
\tag{A.4}
V_0 = \frac{\rho_f \Delta v}{\rho_0 - \rho_f}
\end{equation}

This simplifies to Equation A.4b.

\begin{equation}
\tag{A.4b}
V_o = \frac{F_0 F_f (l_f - l_0)}{\rho_0 (F_0 - F_f)}
\end{equation}

This allows us to find the volume at any given cylinder length $l$ as shown in Equation A.5.

\begin{equation}
\tag{A.5}
v(l) = v_0 + A (l - l_0)
\end{equation}

Similarly, we can find the pressure at any given cylinder length $l$ as shown in Equation A.6.

\begin{equation}
\tag{A.6}
\rho(l) = \frac{\rho_0 v_0}{v(l)}
\end{equation}

And finally, we can find the force at any given cylinder length $l$ as shown in Equation A.7.

## Calculating Force and Potential Energy

The force exerted by the plunger can be found by multiplying the pressure by the effective cross-sectional area of the plunger, as shown in Equation B.1.

\begin{equation}
\tag{B.1}
F(l) = \rho(l) A
\end{equation}

When the force exerted at any given angle is known, the potential energy stored in the spring can be found by integrating the force over the distance traveled by the spring piston. This is shown in Equation B.2.

\begin{equation}
\tag{B.2}
U(l) = \int_{l_0}^{l} F(x) dx
\end{equation}

Substituting in the value for $F(x)$ gives equation B.3.

\begin{equation}
\tag{B.3}
U(l) = \int_{l_0}^{l} \frac{A \rho_0 v_0}{v_0 + A (x - l_0)} dx
\end{equation}

This integral evaluates to Equation B.4.

\begin{equation}
\tag{B.4}
U(l) = \rho_0 v_0 \ln \left| \frac{v_0 + A(l - l_0)}{v_0} \right|
\end{equation}

## Calculating Torque

Since this is a relatively complicated system in a radial coordinate system, it will be modeled using Lagrangian mechanics. The general form of the Lagrangian is given by Equation C.1 below, where $L$ is the Lagrangian, $T$ is the kinetic energy, $U$ is the potential energy, and $\omega = \frac{d \theta}{dt}$ is the angular velocity of the spring arm.

\begin{equation}
\tag{C.1}
\mathcal{L}(\theta, \omega) = T(\theta, \omega) - U(\theta, \omega)
\end{equation}

In the case of this system, the kinetic energy is caused entirely by the rotational motion of the spring arm, weapon arm, and weapon disk and the potential energy is caused entirely by the elastic potential energy stored in the spring. However, in order to calculate the kinetic energy of the spring arm it is necessary to know the moment of inertia of the weapon system. This is given in Equation C.2.

\begin{equation}
\tag{C.2}
I = \frac{1}{3} \lambda_s l_s^3 + \frac{1}{3} \lambda_w l_w^3 + \frac{1}{2} m r^2 + m l_w^2
\end{equation}

The kinetic energy of the spring arm is given by Equation C.3.

\begin{equation}
\tag{C.3}
T(\theta, \omega) = \frac{1}{2} I \omega^2
\end{equation}

To calculate the potential energy of the spring, the length of the spring at any given angle $\theta$ must be calculated first. This is given by Equation C.4.

\begin{equation}
\tag{C.4}
L(\theta) = \sqrt{(l_0 - l_s)^2 + l_s^2 + 2 (l_0 - l_s) l_s \cos(\theta)}
\end{equation}

Once Equation C.4 is found, it is possible to find the potential energy stored in the spring using Equation B.4. This is shown in Equation C.5.

\begin{equation}
\tag{C.5}
U(\theta, \omega) = \rho_0 v_0 \ln \left| \frac{v_0 + A(L(\theta) - l_0)}{v_0} \right|
\end{equation}

Next, the partial derivative's of the Lagrangian with respect to $\theta$ and $\omega$ must be found. These are given by Equations C.6 and C.7.

\begin{equation}
\tag{C.6}
\frac{\partial \mathcal{L}}{\partial \theta} = -\frac{\rho_0 v_0 A l_s (l_0 - l_s) \sin(\theta)}{L(\theta)(v_0 + A (L(\theta) - l_0))}
\end{equation}

\begin{equation}
\tag{C.7}
\frac{\partial \mathcal{L}}{\partial \omega} = I \omega
\end{equation}

From these two equations, the full Euler-Lagrange equation can be found as shown in Equation C.8, where $\alpha = \frac{d \omega}{dt}$ is the angular acceleration of the spring arm.

\begin{equation}
\tag{C.8}
I \alpha = -\frac{\rho_0 v_0 A l_s (l_0 - l_s) \sin(\theta)}{L(\theta)(v_0 + A (L(\theta) - l_0))}
\end{equation}

And since torque is defined as the moment of inertia times the angular acceleration, Equation 8 can be rewritten as Equation C.9 to solve for the torque at any given angle $\theta$.

\begin{equation}
\tag{C.9}
\tau(\theta) = -\frac{\rho_0 v_0 A l_s (l_0 - l_s) \sin(\theta)}{L(\theta)(v_0 + A (L(\theta) - l_0))}
\end{equation}

## Calculating Radial and Angular Forces

Assuming that all the forces are acting in the same plane that the arms rest in, we can find the relation between the total force, the radial force, and the angular force acting on the arm as shown in Equation D.1.

\begin{equation}
\tag{D.1}
F_{total}^2 = F_\theta^2 + F_r^2
\end{equation}

The total force exerted will be the compression force of the gas spring given in Equation B.1, and the angular force will be the torque given in Equation C.9 divided by the length of the spring arm, $l_s$. This yields Equation D.2.

\begin{equation}
\tag{D.2}
\left( \frac{A\rho_0 v_0}{v_0 + A(L(\theta) - l_0)} \right)^2 = \left( -\frac{\rho_0 v_0 A l_s (l_0 - l_s) \sin(\theta)}{l_sL(\theta)(v_0 + A (L(\theta) - l_0))}\right)^2 + F_r^2
\end{equation}

Which simplifies to Equation D.3.

\begin{equation}
\tag{D.3}
F_r = \frac{F_0 v_0 \sqrt{L(\theta)^2 - (l_0 - l_s)^2 \sin^2(\theta)}}{L(\theta)(v_0 + A(L(\theta) - l_s))}
\end{equation}