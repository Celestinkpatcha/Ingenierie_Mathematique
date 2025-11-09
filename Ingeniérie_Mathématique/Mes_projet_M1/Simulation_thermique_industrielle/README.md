# 🧊 Simulation du refroidissement d’une pièce moulée en aluminium

## 1. Contexte industriel

Lors du moulage d’alliages d’aluminium, les pièces présentent souvent des **défauts internes** tels que des microfissures ou des contraintes résiduelles.  
Ces défauts proviennent d’un **refroidissement non uniforme** : certaines zones de la pièce se refroidissent trop vite par contact avec le moule, tandis que le cœur reste chaud.  
Ces différences de température engendrent des **gradients thermiques** qui provoquent des **contraintes mécaniques** susceptibles de fissurer la pièce.

**Objectif du projet :**  
Développer une simulation numérique du refroidissement d’une pièce en aluminium pour :
- visualiser la propagation de la température dans le temps,
- identifier les zones à refroidissement lent ou rapide,
- évaluer les gradients thermiques et les contraintes associées,
- tester plusieurs conditions aux limites (moule froid ou convection),
- proposer des pistes d’optimisation du procédé de moulage.

---

## 2. Modélisation physique




## 2. Modélisation physique

Le modèle repose sur l’équation de la chaleur transitoire en deux dimensions :

$$
\frac{\partial T}{\partial t}
= \alpha \left(
\frac{\partial^2 T}{\partial x^2}
+ \frac{\partial^2 T}{\partial y^2}
\right)
$$

avec :

- $T(x, y, t)$ : température (°C)
- $\alpha = \dfrac{k}{\rho c_p}$ : diffusivité thermique (m²/s)
- $k$ : conductivité thermique de l’aluminium (≈ 180 W·m⁻¹·K⁻¹)
- $\rho$ : masse volumique (≈ 2700 kg·m⁻³)
- $c_p$ : capacité thermique (≈ 900 J·kg⁻¹·K⁻¹)

Valeur typique utilisée :

$$
\alpha = 9 \times 10^{-5} \ \text{m}^2\text{/s}
$$

### 2.1. Condition initiale

La pièce sort du moule à haute température :

$$
T(x, y, 0) = T_{\text{cast}} = 700^\circ \text{C}
$$

### 2.2. Conditions aux limites

Deux types de conditions ont été étudiés.

**(1) Condition de Dirichlet (moule à température fixe)**

$$
T \big|_{\partial \Omega} = T_{\text{mold}} = 250^\circ \text{C}
$$

**(2) Condition de convection (loi de Newton)**

$$
- k \, \frac{\partial T}{\partial n} = h \, \bigl( T - T_\infty \bigr)
$$

où $h$ est le coefficient d’échange (W·m⁻²·K⁻¹) et $T_\infty$ la température du fluide/moule.

---

## 3. Discrétisation numérique

Le domaine rectangulaire étudié correspond à une section de la pièce moulée :

$$
L_x = 0{,}20 \ \text{m}, \qquad
L_y = 0{,}10 \ \text{m}
$$

Il est discrétisé en

- $N_x = 81$ points selon $x$
- $N_y = 41$ points selon $y$

d’où

$$
\Delta x = \frac{L_x}{N_x - 1},
\qquad
\Delta y = \frac{L_y}{N_y - 1}
$$

Dans la suite on prend $\Delta x = \Delta y$.

### 3.1. Schéma explicite (FTCS)

Le schéma aux différences finies utilisé est le schéma explicite “Forward Time, Central Space” :

$$
T_{i,j}^{n+1}
=
T_{i,j}^n
+
\alpha \, \Delta t \left[
\frac{T_{i+1,j}^n - 2 T_{i,j}^n + T_{i-1,j}^n}{\Delta x^2}
+
\frac{T_{i,j+1}^n - 2 T_{i,j}^n + T_{i,j-1}^n}{\Delta y^2}
\right]
$$

### 3.2. Condition de stabilité

Pour ce schéma explicite en 2D (avec $\Delta x = \Delta y$), le pas de temps doit vérifier :

$$
\Delta t \le \frac{\Delta x^2}{4 \, \alpha}
$$

Dans le code on choisit :

```python
dt_stable = dx**2 / (4 * alpha)
dt = 0.2 * dt_stable

## 🛠️ Technologies
- Python
- NumPy
- Matplotlib
- Spyder
- Jupyter Notebook

## 🚀 Installation
```bash
pip install -r requirements.txt

