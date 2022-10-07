# Simulated Tempering
## Theoretical background
### The expanded partition function
The expanded hamiltonian could be expressed as (1):
$$
H(\vec{p}, \vec{q}, T) = K(\vec{p}) + V(\vec{q}) + \lambda(T)
$$

For this system the expanded partition function reads (2):
$$
Z_{exp} \equiv \int_{T_1}^{T_2}\,dT \int_{\vec{p_1}}^{\vec{p_1}}\int_{\vec{q_1}}^{\vec{q_2}} e^{-\beta(T)H(\vec{p}, \vec{q}, T)}\,d\vec{q}\,d\vec{p}
$$
where $\beta(T) = \frac{1}{k_BT}$ (3), and we are supposing that the phase space is defined by the integral limits

Open the integrand (4):
$$
Z_{exp} = \int_{T_1}^{T_2} e^{-\beta(T)\lambda(T)}\,dT \int_{\vec{p_1}}^{\vec{p_1}}\int_{\vec{q_1}}^{\vec{q_2}} e^{-\beta(T)(K(\vec{p}) + V(\vec{q})}\,d\vec{q}\,d\vec{p}
$$
and because $H(\vec{p}, \vec{q}) = K(\vec{p}) + V(\vec{q})$ (5) and $Z(T) \equiv \int_{\vec{p_1}}^{\vec{p_1}}\int_{\vec{q_1}}^{\vec{q_2}} e^{-\beta(T)H(\vec{p}, \vec{q})}\,d\vec{q}\,d\vec{p}$ (6)

(7)
$$
Z_{exp} = \int_{T_1}^{T_2} e^{-\beta(T)\lambda(T)}Z(T)\,dT 
$$

(8) For the case that $T = \{T_m\}$ with $m\in [1,M] :$
$$
Z_{exp} = \sum_{m = 1}^{M} e^{-\beta(T_m)\lambda(T_m)}Z(T_m)
$$
If we says that the function $\lambda(T) \equiv -\frac{g(T)}{\beta(T)}$ (9), $g(T)$ is a general function of $T$, . For simplicity, let says that any function of temperature is $f_m \equiv {f(T_m)} and $. Then we could rewrite the last expression as:

(10)
$$
Z_{exp} = \sum_{m = 1}^{M} e^{g_m}Z_m
$$
We could interpret this expression as every temperature state is been weighted by a factor of $e^{g_m}$
### Probability of a temperature state
Let's stick with the expanded ensemble with a discrete set of temperatures and $\lambda(T) = -\frac{g(T)}{\beta(T)}$

We want that all "temperature-sate" have the same probability, for that we take any two state ($a$ and $b$) and impose the equality $P_a = P_b$ (is the same of $P(T_a) = P(T_b)$): Therefore:

(11)
$$
\frac{e^{g_a}Z_a}{Z_{exp}} = \frac{e^{g_b}Z_b}{Z_{exp}}
$$
Which gives (12):
$$
e^{g_a}Z_a = e^{g_b}Z_b
$$
And applying logarithms (13):
$$
\ln(Z_a)- \ln(Z_b)  = g_b - g_a 
$$
This is the general expression. Now, from statistic thermodynamics we know that the free energy could be expressed as (14):
$$
F(T) \equiv -k_BT\ln(Z(T)) = -\frac{\ln(Z(T))}{\beta(T)}
$$
Adding to the last expression (15):
$$
\beta_bF_b- \beta_aF_a  = g_b - g_a 
$$
So, the difference of the exponents of the weights is related with the difference on free energy between the two temperature states.
### Dealing with the left hand side

Now,let's try to get a more useful expression for the left hand side, because computing the free energy is not feasible in many biophysical systems.

We will start from the laws of thermodynamic:

First law (16):
$$
dU \equiv \delta{q} + \delta{w} = \delta{q} - PdV
$$

Second law (17):
$$
dS_s \geq \frac{\delta{q}}{T}
$$
Where $dS_s$ is the entropy of the system. To create an equality, we add to the right hand the term generation of entropy (18):
$$
dS_s = \frac{\delta{q}}{T} + \delta{S_i}:\delta{S_i} \geq 0
$$
(19) (18) in (16)

$$
dU = TdS_s - PdV - T\delta{S_i}
$$
(20) = (19) $ - d(TS_s)$
$$
dU -d(TS_s) = \cancel{TdS_s} - PdV - T\delta{S_i} - \cancel{TdS_s} - S_sdT
$$
$$
d(U-TS_s) = -S_sdT - PdV - T\delta{S_i}
$$
Now we defined $F \equiv U-TS_s$ (21) where $F \equiv f(T,V,n_i)$ (22)

(23) because of (21) in (20)
$$
d(F) = -S_sdT - PdV - T\delta{S_i}
$$

(24) Creating the differential form of (22)
$$
dF = \left( \frac{\partial{F}}{\partial{T}} \right)_{V,n}dT + \left( \frac{\partial{F}}{\partial{V}} \right)_{T,n}dV + \sum_{i = 1}^{N}\left( \frac{\partial{F}}{\partial{n_i}} \right)_{T,V,n \neq n_i}dn_i
$$
(25) Comparing (20) and (24)
$$
-S_s = \left( \frac{\partial{F}}{\partial{T}} \right)_{V,n}
$$
(26) = (25) in (21)
$$
F = U + T\left( \frac{\partial{F}}{\partial{T}} \right)_{V,n}
$$
(27) Dividing by $1/T^2$ and rearrange:

$$
\frac{F}{T^2} = \frac{U}{T^2} + \frac{1}{T}\left( \frac{\partial{F}}{\partial{T}} \right)_{V,n}
$$
$$
-\frac{U}{T^2} = \frac{1}{T}\left( \frac{\partial{F}}{\partial{T}} \right)_{V,n} - \frac{F}{T^2}
$$
$$
-\frac{U}{T^2} = \left( \frac{\partial{(F/T)}}{\partial{T}} \right)_{V,n}
$$
Notes that the right hand side is the derivative of a fraction.

(28) Integrating (27) from $T_a$ to $T_b$ and diving by $1/k_B$
$$
\int_{T_a}^{T_b}-\frac{U(T)}{T^2}\,dT = \int_{T_a}^{T_b}\,d\left( \frac{F(T)}{T} \right):V,n = cte
$$

$$
\int_{T_a}^{T_b}-\frac{U(T)}{T^2}\,dT = \frac{F(T_b)}{T_b} -\frac{F(T_a)}{T_a};1/k_B
$$
$$
-\frac{1}{k_B}\int_{T_a}^{T_b}\frac{U(T)}{T^2}\,dT = \beta_bF_b - \beta_aF_a
$$

(29) Comparing (28) with (15)
$$
g_b - g_a = -\frac{1}{k_B}\int_{T_a}^{T_b}\frac{U(T)}{T^2}\,dT
$$
With a similar analysis we could get (30):
$$
\beta_bG_b- \beta_aG_a  = -\frac{1}{k_B}\int_{T_a}^{T_b}\frac{H(T)}{T^2}\,dT
$$
So, now the problem fo calculate the Helmholtz energy in two states of temperature is translated to know how the internal energy change in between the temperature states $T_a$ and $T_b$.

### Estimating the integral $\int_{T_a}^{T_b}\frac{f(T)}{T^2}\,dT$
Here we will try to get a more convenient shape to the integral $\int_{T_a}^{T_b}\frac{f(T)}{T^2}\,dT$ (31)

Writhing $f(T)$ as a Taylor expansion (32):
$$
f(T) = \sum_{n = 0}^{\infty}\frac{f^{(n)}(T_0)}{n!}(T-T_0)^n
$$
Substituting (32) in (31) and taking the integral inside the serie -> (33)
$$
\int_{T_a}^{T_b}\sum_{n = 0}^{\infty}\frac{f^{(n)}(T_0)}{n!} \frac{(T-T_0)^n}{T^2}\,dT
$$
$$
\sum_{n = 0}^{\infty}\frac{f^{(n)}(T_0)}{n!}\int_{T_a}^{T_b} \frac{(T-T_0)^n}{T^2}\,dT
$$
Warning: I not completely sure if taking the integral inside the serie is completely correct. Without $1/T^2$ for sure because converge to $f(T)$ but when we divided by $1/T^2$ I am not sure.

With the Newton's binomial (34):
$$
(x+y)^n = \sum_{k=0}^n \binom{n}{k}x^{n-k}y^k
$$
We got (35)
$$
\sum_{n = 0}^{\infty}\frac{f^{(n)}(T_0)}{n!}\int_{T_a}^{T_b} \frac{\sum_{k=0}^n \binom{n}{k}T^{n-k}(-T_0)^k}{T^2}\,dT
$$
$$
\sum_{n = 0}^{\infty}\frac{f^{(n)}(T_0)}{n!}\int_{T_a}^{T_b} \sum_{k=0}^n \binom{n}{k}T^{n-k-2}(-T_0)^k\,dT
$$
$$
\sum_{n = 0}^{\infty} \left[ \frac{f^{(n)}(T_0)}{n!} \sum_{k=0}^n \left[ \binom{n}{k} (-T_0)^k \int_{T_a}^{T_b}T^{n-k-2}\,dT \right] \right]
$$
Then we got (36):
$$
\int_{T_a}^{T_b}\frac{f(T)}{T^2}\,dT = \sum_{n = 0}^{\infty} \frac{f^{(n)}(T_0)}{n!}C(n, T_0, T_a, T_b)
$$
Such as (37):
$$
C(n, T_0, T_a, T_b) = \sum_{k=0}^n \left[ \binom{n}{k} (-T_0)^k \int_{T_a}^{T_b}T^{n-k-2}\,dT \right]
$$

Finally, going back again to equation (30) we got (38):
$$
g_b - g_a = -\frac{1}{k_B} \sum_{n = 0}^{\infty} \frac{U^{(n)}(T_0)}{n!}C(n, T_0, T_a, T_b)
$$

### First and second order approximations
Let's calculate the first two constant of the serie (39)

$$
C(0, T_0, T_a, T_b) = \int_{T_a}^{T_b}T^{-2}\,dT = -\frac{1}{T}\Big|_{T_a}^{T_b} = \frac{1}{T_a} - \frac{1}{T_b}
$$
$$
C(1, T_0, T_a, T_b) = \int_{T_a}^{T_b}T^{-1}\,dT + (-T_0) \int_{T_a}^{T_b}T^{-2}\,dT = \ln{\frac{T_b}{T_a}} + T_0 \left[ \frac{1}{T_b} - \frac{1}{T_a}\right]
$$
Now, if $T_0 = \frac{T_a+T_b}{2}$:
#### First order
(40)
$$
g_b - g_a \approx -\frac{1}{k_B}U(T_0)C(0, T_0, T_a, T_b) = (\beta_b - \beta_a)U\left(\frac{T_a + T_b}{2}\right)
$$
This result is similar to the one got it by [10.1103/PhysRevE.76.016703](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.76.016703). They main differences are that they work with total energy, in this case is the internal energy; and they use $\frac{E_a + E_b}{2}$ rather than $U\left(\frac{T_a + T_b}{2}\right)$.

#### Second order
If we wold like go further in the proximation we could also included the second term of the serie (41):
$$
g_b - g_a \approx-\frac{1}{k_B}[ U(T_0)C(0, T_0, T_a, T_b) + U(T_0)'C(1, T_0, T_a, T_b)] 
$$
$$
= (\beta_b - \beta_a)U\left(\frac{T_a + T_b}{2}\right) +\left[ \frac{1}{k_B}\ln{\frac{T_a}{T_b}} + \frac{(T_a + T_b)(\beta_a-\beta_b)}{2}\right]U'\left(\frac{T_a + T_b}{2}\right)
$$
$$
= (\beta_b - \beta_a)U\left(\frac{T_a + T_b}{2}\right) +\left[ \frac{1}{k_B}\ln{\frac{T_a}{T_b}} + \frac{(T_a + T_b)(\beta_a-\beta_b)}{2}\right]C_v\left(\frac{T_a + T_b}{2}\right)
$$
Where $C_v$ is the heat capacity at $V = cte$

### Calculating the weights
Let suppose for simplicity that we only have five temperature states. The following results could be easily generalized.
We got that (42):
$$
g_b - g_a = I_{ab}
$$
The notation $I_{ab}$ is if equations (29) or (38) was used.

Now, we only need to solve the system of equations (43).
$$
\begin{pmatrix}
-1 & 1 & 0 &0 & 0 \\
0 & -1 & 1 & 0 & 0 \\
0 & 0 & -1 & 1 & 0 \\
0 & 0 & 0 & -1 & 1 \\
1 & 1 & 1 & 1 & 1 \\
\end{pmatrix}
\begin{pmatrix}
g_1\\
g_2\\
g_3\\
g_4\\
g_5\\   
\end{pmatrix}
=
\begin{pmatrix}
I_{12}\\
I_{23}\\
I_{34}\\
I_{45}\\
0\\   
\end{pmatrix}
$$ 
Where the last equation is a restriction: $\sum_{m = 1}^M g_m = 0$ 

## From the simulation perspective
In the end we need the difference between the exponents of the weights, it doesn't matter if we got them by equation (29) (exact solution) or (38) (still exact if we use infinite terms). If we know how exactly $U$ change with the temperature, we could go straight to equation (29), however if we just know the value of $U$ (and also $C_v$ for the second order approximation) in the discrete set of temperatures under investigation, we should go with equation (38) (more precisely 41)


## Compatibility between HMR and Simulated Tempering
The joint probability distribution function for the expanded ensemble $H(\vec{p}, \vec{q}, T) = K(\vec{p}) + V(\vec{q}) - \frac{g(T)}{\beta(T)}$ could be written as:

$$
pdf(\vec{p}, \vec{q}, T) = \frac {e^{-\beta(T)H(\vec{p}, \vec{q}, T)}}{Z_{exp}} = \frac {e^{-\beta(T)K(\vec{p})}e^{-\beta(T)V(\vec{q})}e^{g(T)}}{Z_{exp}}
$$
Where 
$$
Z_{exp} \equiv \sum_{m = 1}^{M} e^{g(T_m)}Z(T_m)
$$
The marginal probability mass function of $T$ (because $T$ is discrete)
$$
pmf(T_i) = \frac{e^{g(T_i)}Z(T_i)}{Z_{exp}}; \\
Z(T_i) =  \int_{\vec{q_1}}^{\vec{q_2}} e^{-\beta(T_i) V(\vec{q})}\,d\vec{q} \int_{\vec{p_1}}^{\vec{p_2}} e^{-\beta(T_i)K(\vec{p})}\,d\vec{p}
$$
For simulated tempering we want that all the temperature states have the same probability, if we have M states, then:
$$
pmf(T_i) = \frac{e^{g(T_i)}Z(T_i)}{Z_{exp}} = \frac{2}{M(M+1)}
$$
And this impositions still hold in the case that we change the masses despite that $Z_{exp}$,  $Z(T_i)$ and $g(T_i)$ will change. The numerator is the most important part to work in order to get equal probabilities, and because the changes on masses modify $Z(T_i)$, $g(T_i)$ will also be modified.  So, we are not able to use the weights obtained without HMR in a simulation with HMR. But our goal holds: "Have the same probability on each temperature state".



In the end we will **only** use the information from the target temperature ($T = T_a$). Therefore we will need the conditional probability distribution function $pdf(\vec{p}, \vec{q}|T = T_a)$.

To calculate such object, we need to find all the intersection between the ensemble ($\vec{p}, \vec{q}, T$) and the condition ($T = T_a$) and then dive it by all the elements that belongs to the condition.

In math, if we have $n$ continues random variables $\{R\}_n$ and we want to get the conditional probability over the random variable $R_x \in \{R\}_n$ , 
$$
pdf(\bold{R}:R_x \notin \bold{R}|R_x \Rightarrow C) = \frac{\int_{C} pdf(\bold{R}, R_x) \,dR_x}{\int \int_{C} pdf(\bold{R}, R_x) \,dR_x\,d\bold{R}}
$$
In our case $T$ is discrete, so is easier.

$$
pdf(\vec{p}, \vec{q}|T = T_a) = \frac{e^{-\beta(T_a)K(\vec{p})}e^{-\beta(T_a)V(\vec{q})}e^{g(T_a)}} {\int_{\vec{q_1}}^{\vec{q_2}} \int_{\vec{p_1}}^{\vec{p_2}} e^{-\beta(T_a)K(\vec{p})}e^{-\beta(T_a)V(\vec{q})}e^{g(T_a)} \,d\vec{p}\,d\vec{q}}
\\
pdf(\vec{p}, \vec{q}|T = T_a) = \frac{e^{-\beta(T_a)K(\vec{p})}e^{-\beta(T_a)V(\vec{q})}\cancel{e^{g(T_a)}}} {\cancel{e^{g(T_a)}} \int_{\vec{q_1}}^{\vec{q_2}} \int_{\vec{p_1}}^{\vec{p_2}} e^{-\beta(T_a)K(\vec{p})}e^{-\beta(T_a)V(\vec{q})} \,d\vec{p}\,d\vec{q}}
\\
pdf(\vec{p}, \vec{q}|T = T_a) = \frac{e^{-\beta(T_a)K(\vec{p})}e^{-\beta(T_a)V(\vec{q})}} {\int_{\vec{q_1}}^{\vec{q_2}} e^{-\beta(T_i) V(\vec{q})}\,d\vec{q} \int_{\vec{p_1}}^{\vec{p_2}} e^{-\beta(T_i)K(\vec{p})}\,d\vec{p}}
$$

Finally, we are interested on those properties that only depends on the coordinates, such as $A(\vec{q})$. If we want to get the expected value from the expanded ensemble given that $T = T_a$; then we need to compute the conditional expectation value:
$$
E[A(\vec{q})|T=T_a] = \int_{\vec{q_1}}^{\vec{q_2}} \int_{\vec{p_1}}^{\vec{p_2}} A(\vec{q})pdf(\vec{p}, \vec{q}|T = T_a) \,d\vec{p}\,d\vec{q}
\\
E[A(\vec{q})|T=T_a] = \frac{\int_{\vec{q_1}}^{\vec{q_2}} A(\vec{q})e^{-\beta(T_a)V(\vec{q})} \,d\vec{q} \cancel{\int_{\vec{p_1}}^{\vec{p_2}} e^{-\beta(T_a)K(\vec{p})}} \,d\vec{p}}{\int_{\vec{q_1}}^{\vec{q_2}} e^{-\beta(T_a) V(\vec{q})}\,d\vec{q} \cancel{\int_{\vec{p_1}}^{\vec{p_2}} e^{-\beta(T_a)K(\vec{p})}}\,d\vec{p}}
\\
E[A(\vec{q})|T=T_a] =\frac{\int_{\vec{q_1}}^{\vec{q_2}} A(\vec{q})e^{-\beta(T_a)V(\vec{q})} \,d\vec{q}} {\int_{\vec{q_1}}^{\vec{q_2}} e^{-\beta(T_a) V(\vec{q})}\,d\vec{q}}
$$
As we see in the last equation the expected value only depend on the coordinates.