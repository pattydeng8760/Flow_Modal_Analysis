# Modal Analysis Toolkit for Fluid Flows

This repository provides a small, modular toolkit for performing modal analysis of fluid-flow data, with an emphasis on methods commonly used in computational fluid dynamics (CFD) and experimental fluid mechanics:

- **Proper Orthogonal Decomposition (POD)**
- **Spectral Proper Orthogonal Decomposition (SPOD)**
- **Dynamic Mode Decomposition (DMD)**

The code is organized so that it can be applied to generic time-resolved fields (e.g., velocity, pressure, vorticity) extracted from simulations or measurements. The goal is to obtain energetically and dynamically relevant low-dimensional representations of complex flows, following the framework summarized in Taira *et al.* (2017), *AIAA Journal*.

---

## 1. Mathematical Background

This section summarizes the core mathematics of POD, SPOD, and DMD, using consistent notation. The description is oriented toward CFD / aeroacoustics applications where the input is a time series of flow snapshots.

### 1.1 Data Layout and Notation

Let $\mathbf{q}(\boldsymbol{\xi}, t)$ denote a scalar or vector flow quantity (e.g., one velocity component, pressure, vorticity) defined over a spatial grid with $n$ degrees of freedom and sampled at $m$ time instants $t_1, \dots, t_m$.

We form column vectors,

$$
\mathbf{x}_k \in \mathbb{R}^n, \quad \mathbf{x}_k 
= \mathbf{q}(\boldsymbol{\xi}, t_k),
$$

and collect them into a snapshot matrix,

$$
\mathbf{X} = 
\begin{bmatrix}
\mathbf{x}_1 & \mathbf{x}_2 & \cdots & \mathbf{x}_m
\end{bmatrix} 
\in \mathbb{R}^{n \times m}.
$$

In many flows, it is useful to separate the data into mean and fluctuations by means of Reynolds decomposition,

$$
\bar{\mathbf{x}} = \frac{1}{m} \sum_{k=1}^m \mathbf{x}_k, 
\qquad 
\mathbf{x}'_k = \mathbf{x}_k - \bar{\mathbf{x}},
$$

and define the fluctuation matrix

$$
\mathbf{X}' = 
\begin{bmatrix}
\mathbf{x}'_1 & \mathbf{x}'_2 & \cdots & \mathbf{x}'_m
\end{bmatrix}.
$$

All three methods (POD, SPOD, DMD) act on variants of this snapshot matrix or its frequency-domain counterpart.

---

### 1.2 Proper Orthogonal Decomposition (POD)

#### 1.2.1 Objective

POD seeks an orthonormal basis $\{ \boldsymbol{\phi}_j \}_{j=1}^r$ that optimally (in an $L^2$ or kinetic-energy sense) represents the fluctuations:

$$
\mathbf{x}'_k \approx \sum_{j=1}^r a_j(t_k)\, \boldsymbol{\phi}_j.
$$

The POD modes $\boldsymbol{\phi}_j$ are spatial structures; the coefficients $a_j(t_k)$ are temporal amplitudes. The modes are ordered such that mode 1 captures the largest fraction of energy, mode 2 the next largest, etc.

#### 1.2.2 Correlation and Eigenvalue Problem

The classical POD is defined via the spatial correlation (covariance) tensor

$$
\mathbf{C} = \frac{1}{m} \mathbf{X}' {\mathbf{X}'}^\top \in \mathbb{R}^{n \times n}.
$$

We seek eigenpairs,

$$
\mathbf{C} \boldsymbol{\phi}_j = \lambda_j \boldsymbol{\phi}_j, 
\qquad \lambda_1 \ge \lambda_2 \ge \dots \ge 0,
$$

where $\boldsymbol{\phi}_j$ are POD modes and $\lambda_j$ are the associated modal energies. Because $n$ is typically very large, one rarely forms $\mathbf{C}$ explicitly.

#### 1.2.3 Method of Snapshots / SVD Formulation

In practice, POD is computed via the method of snapshots or singular value decomposition (SVD). Define the temporal correlation matrix

$$
\mathbf{R} = \frac{1}{m} {\mathbf{X}'}^\top \mathbf{X}' \in \mathbb{R}^{m \times m}.
$$

Solve,

$$
\mathbf{R} \mathbf{a}_j = \lambda_j \mathbf{a}_j,
$$

where $\mathbf{a}_j \in \mathbb{R}^m$ are eigenvectors of $\mathbf{R}$.

The POD modes are reconstructed as,

$$
\boldsymbol{\phi}_j 
= \frac{1}{\sqrt{m \lambda_j}} \mathbf{X}' \mathbf{a}_j.
$$

Equivalently, perform an economy-sized SVD:

$$
\mathbf{X}' = \mathbf{U} \boldsymbol{\Sigma} \mathbf{V}^\top,
$$

with
- $\mathbf{U} \in \mathbb{R}^{n \times r}$: POD modes (left singular vectors),
- $\mathbf{V} \in \mathbb{R}^{m \times r}$: normalized temporal coefficients (right singular vectors),
- $\boldsymbol{\Sigma} = \text{diag}(\sigma_1, \dots, \sigma_r)$: singular values.

Then,

$$
\boldsymbol{\phi}_j = \mathbf{u}_j, \qquad 
\lambda_j = \frac{\sigma_j^2}{m},
$$

and the temporal coefficients are,

$$
a_j(t_k) = \sigma_j v_{kj},
$$

where $v_{kj}$ is the $k$-th component of the right singular vector $\mathbf{v}_j$.

#### 1.2.4 Reconstruction

The fluctuation field is approximated as

$$
\mathbf{x}'_k \approx \sum_{j=1}^r a_j(t_k) \boldsymbol{\phi}_j,
\quad
\mathbf{x}_k \approx \bar{\mathbf{x}} + \sum_{j=1}^r a_j(t_k) \boldsymbol{\phi}_j.
$$

By truncating to a small number of leading modes, one obtains a low-dimensional representation that captures a large fraction of the total energy.

#### 1.2.5 Benefits and Drawbacks of POD

**Benefits**

- **Energy optimality:** For a given number of modes, POD provides the best mean-square (energy) approximation of the data.
- **Orthogonality:** Modes are orthonormal, simplifying projection, reconstruction, and reduced-order modeling.
- **Noise separation:** Incoherent noise tends to appear in high-order modes and can be filtered by truncating the expansion.
- **Computationally efficient:** The snapshot/SVD formulation is well suited to large CFD datasets and can exploit optimized linear algebra libraries.

**Drawbacks**

- **Energy-based, not dynamical:** POD arranges modes by energy, not by dynamical importance. Energetically weak but dynamically crucial structures may appear in high-order modes.
- **Mixed frequency content:** Temporal coefficients $a_j(t)$ often contain multiple frequencies. Individual POD modes are not, in general, single-frequency structures.
- **Second-order statistics only:** POD is built on second-order correlations and does not directly capture higher-order statistics.
- **Truncation ambiguity:** It is not always obvious how many modes to retain; various ad hoc criteria (energy thresholds, spectral gaps) are used in practice.

---

### 1.3 Spectral Proper Orthogonal Decomposition (SPOD)

SPOD is a frequency-domain variant of POD tailored to statistically stationary flows. It yields modes that are simultaneously **orthogonal in space** and **monochromatic in time** (single frequency), making it particularly useful for separating coherent structures by frequency content.

#### 1.3.1 Data Segmentation and Fourier Transform

Assume we have a long time series of snapshots with uniform sampling period $\Delta t$. We divide the data into $n_b$ blocks (possibly overlapping), each containing $m_{\text{FFT}}$ snapshots:

$$
\mathbf{X}^{(\ell)} =
\begin{bmatrix}
\mathbf{x}_{1}^{(\ell)} & \cdots & \mathbf{x}_{m_{\text{FFT}}}^{(\ell)}
\end{bmatrix},
\quad \ell = 1, \dots, n_b.
$$

For each block $\ell$, compute a temporal discrete Fourier transform (e.g., via FFT) to obtain,

$$
\hat{\mathbf{X}}^{(\ell)} =
\begin{bmatrix}
\hat{\mathbf{x}}^{\omega_1, (\ell)} & 
\hat{\mathbf{x}}^{\omega_2, (\ell)} & 
\cdots & 
\hat{\mathbf{x}}^{\omega_{m_{\text{FFT}}}, (\ell)}
\end{bmatrix}
\in \mathbb{C}^{n \times m_{\text{FFT}}},
$$

where $\omega_k$ are discrete angular frequencies.

At a fixed frequency $\omega_k$, collect all realizations into,

$$
\hat{\mathbf{X}}^{\omega_k} =
\begin{bmatrix}
\hat{\mathbf{x}}^{\omega_k, (1)} &
\hat{\mathbf{x}}^{\omega_k, (2)} & 
\cdots &
\hat{\mathbf{x}}^{\omega_k, (n_b)}
\end{bmatrix}
\in \mathbb{C}^{n \times n_b}.
$$

#### 1.3.2 Cross-Spectral Density and Eigenvalue Problem

The cross-spectral density (CSD) tensor at frequency $\omega_k$ is,

$$
\mathbf{S}(\omega_k)
= \frac{1}{n_b} \hat{\mathbf{X}}^{\omega_k} 
\, \hat{\mathbf{X}}^{\omega_k\, *^\top}
\in \mathbb{C}^{n \times n},
$$

where $*$ denotes complex conjugation.

SPOD modes $\mathbf{\phi}_{\omega_k, j}$ and their energies $\lambda_{\omega_k, j}$ are defined as the eigenpairs of the CSD:

$$
\mathbf{S}(\omega_k) \, \mathbf{\phi}_{\omega_k, j}
= \lambda_{\omega_k, j} \, \mathbf{\phi}_{\omega_k, j},
\qquad
\lambda_{\omega_k, 1} \ge \lambda_{\omega_k, 2} \ge \dots \ge 0.
$$

Each mode $\mathbf{\phi}_{\omega_k, j}$ represents a spatial structure that oscillates at the single frequency $\omega_k$, with energy $\lambda_{\omega_k, j}$.

In practice, as with spatial POD, one uses an SVD of the reduced matrix at each frequency to compute SPOD modes efficiently.

#### 1.3.3 Benefits and Drawbacks of SPOD

**Benefits**

- **Frequency-resolved coherent structures:** Modes at each frequency are optimally energetic and orthogonal, providing a clean decomposition of broadband turbulent flows into frequency-dependent structures.
- **Connection to resolvent analysis:** For statistically stationary flows, leading SPOD modes often align with resolvent response modes, linking data-driven and operator-based perspectives.
- **Noise handling via Welch’s method:** Block averaging and windowing provide robust spectral estimates and reduce variance in the CSD.
- **Direct relevance to aeroacoustics:** SPOD naturally separates flow and acoustic structures by frequency, aiding interpretation of noise-generation mechanisms.

**Drawbacks**

- **Stationarity assumption:** SPOD relies on statistical stationarity. Strongly transient or non-stationary flows are not well described by a purely spectral approach.
- **Data requirements:** Requires long, high-quality time series to obtain converged spectral estimates across frequencies, which may be expensive to generate or store.
- **Parameter choices:** Results depend on choices such as block length, overlap, and windowing (Welch parameters). Poor choices can smear narrowband features or introduce spectral leakage.
- **Computational cost:** For very large datasets and many frequencies, repeatedly forming and decomposing CSD matrices can be demanding, though parallelization can mitigate this.

---

### 1.4 Dynamic Mode Decomposition (DMD)

DMD is a data-driven technique that approximates a best-fit linear operator governing the evolution of snapshots. It decomposes the flow into modes with single, well-defined growth/decay rates and frequencies. DMD can be interpreted as a finite-dimensional approximation of the Koopman operator for nonlinear systems.

#### 1.4.1 Snapshot Pairs and Linear Mapping

Assume uniformly spaced snapshots with time step $\Delta t$. Define two snapshot matrices:

$$
\mathbf{X} =
\begin{bmatrix}
\mathbf{x}_1 & \mathbf{x}_2 & \cdots & \mathbf{x}_{m-1}
\end{bmatrix},
\quad
\mathbf{X}' =
\begin{bmatrix}
\mathbf{x}_2 & \mathbf{x}_3 & \cdots & \mathbf{x}_m
\end{bmatrix}.
$$

DMD assumes the existence of a linear operator $\mathbf{A}$ such that,

$$
\mathbf{X}' \approx \mathbf{A}\, \mathbf{X},
$$

where $\mathbf{A} \in \mathbb{C}^{n \times n}$ advances the state by one time step.

In a least-squares sense, the best-fit operator is,

$$
\mathbf{A} = \mathbf{X}' \mathbf{X}^+,
$$

with $\mathbf{X}^+$ the Moore–Penrose pseudoinverse of $\mathbf{X}$. Because $n$ is large, we avoid forming $\mathbf{A}$ explicitly.

#### 1.4.2 Reduced DMD via SVD

Compute an economy SVD of $\mathbf{X}$:

$$
\mathbf{X} = \mathbf{U} \boldsymbol{\Sigma} \mathbf{V}^*,
$$

where
- $\mathbf{U} \in \mathbb{C}^{n \times r}$,
- $\boldsymbol{\Sigma} \in \mathbb{R}^{r \times r}$,
- $\mathbf{V} \in \mathbb{C}^{m-1 \times r}$,
and $^*$ denotes the conjugate transpose.

Project $\mathbf{A}$ onto the $r$-dimensional subspace spanned by $\mathbf{U}$:

$$
\tilde{\mathbf{A}} 
= \mathbf{U}^* \mathbf{A} \mathbf{U}
= \mathbf{U}^* \mathbf{X}' \mathbf{X}^+ \mathbf{U}
= \mathbf{U}^* \mathbf{X}' \mathbf{V} \boldsymbol{\Sigma}^{-1}
\in \mathbb{C}^{r \times r}.
$$

Solve the reduced eigenvalue problem,

$$
\tilde{\mathbf{A}} \mathbf{W} = \mathbf{W} \boldsymbol{\Lambda},
$$

where $\mathbf{\Lambda} = \text{diag}(\lambda_1, \dots, \lambda_r)$ contains the discrete-time DMD eigenvalues and columns of $\mathbf{W}$ are eigenvectors.

The **DMD modes** in the full state space are,

$$
\boldsymbol{\Phi} = \mathbf{X}' \mathbf{V} \boldsymbol{\Sigma}^{-1} \mathbf{W}
\in \mathbb{C}^{n \times r},
$$

with columns $\boldsymbol{\phi}_j$ giving spatial structures.

#### 1.4.3 Growth Rates and Frequencies

Each eigenvalue $\lambda_j$ encodes growth/decay and oscillation:

$$
\lambda_j = e^{(\mu_j + i \omega_j) \Delta t},
$$

where
- $\mu_j = \frac{1}{\Delta t} \Re(\log \lambda_j)$ is the growth/decay rate,
- $\omega_j = \frac{1}{\Delta t} \Im(\log \lambda_j)$ is the (angular) frequency.

Thus, the dynamics are approximated as

$$
\mathbf{x}_k \approx 
\sum_{j=1}^r b_j \boldsymbol{\phi}_j \lambda_j^{k-1},
$$

where coefficients $b_j$ are determined from the initial condition $\mathbf{x}_1$ by solving,

$$
\mathbf{x}_1 = \sum_{j=1}^r b_j \boldsymbol{\phi}_j.
$$

#### 1.4.4 Benefits and Drawbacks of DMD

**Benefits**

- **Single-frequency modes:** Each DMD mode has a single growth/decay rate and oscillation frequency, making it well suited for identifying coherent oscillatory structures (e.g., shedding modes, tonal noise features).
- **Data-driven:** Requires only time-resolved data, with no explicit knowledge of governing equations; applicable to simulations and experiments alike.
- **Koopman interpretation:** Under certain conditions, DMD provides a finite-dimensional approximation of Koopman eigenvalues and modes, giving a linear representation of nonlinear dynamics.
- **Flexible and extensible:** Numerous variants (sparsity-promoting DMD, compressed DMD, streaming DMD) exist to address noise, large datasets, or incomplete sampling.

**Drawbacks**

- **Time-resolved data requirement:** Standard DMD needs snapshot pairs with constant time spacing. Sparse or irregularly sampled data require specialized variants.
- **Nonlinear systems:** For strongly nonlinear dynamics, standard DMD may yield biased or noisy approximations unless enriched with nonlinear observables (extended DMD).
- **Mode ranking ambiguity:** Unlike POD, there is no unique energy-based ranking of DMD modes. Selecting “important” modes can be subjective and problem dependent.
- **Sensitivity to noise:** Noise and experimental uncertainties can significantly affect the spectrum, necessitating regularization or noise-robust formulations.

---

## 2. Code Structure

A typical organization for this repository might be:

- `__main__.py`  
  Entry point for running modal analyses from the command line (e.g., selecting method, field, time window).

- `__init__.py`  
  Makes the package importable as a Python module.

- `extract.py`  
  Utilities for loading and reshaping CFD / experimental data into snapshot matrices $\mathbf{X}$.

- `map_cut.py`  
  Mapping and cut-plane utilities (e.g., extracting 2D slices or specific spatial regions for analysis).

- `triple_decomposition.py`  
  Implements mean–coherent–stochastic decompositions of the flow field (e.g., $\mathbf{q} = \bar{\mathbf{q}} + \tilde{\mathbf{q}} + \mathbf{q}'$), which can be coupled with POD/DMD/SPOD.

- `dmd.py`  
  Core DMD implementation (SVD-based, with options for rank truncation, frequency selection, and reconstruction).
  
- `pod.py`  
  Core POD implementation (SVD-based, with options for rank truncation, frequency selection, and reconstruction).
  
- `ModalAnalysis.py`  
  High-level interface for performing POD, SPOD, and DMD on given datasets (e.g., method selection, parameter handling, and post-processing).

- `utils.py`  
  Shared numerical and plotting utilities (e.g., SVD wrappers, FFT helpers, energy spectra, mode visualization).

You can adapt this structure as needed; the mathematical framework above is agnostic to the specific file layout as long as the same data conventions are followed.

---

## 3. Typical Workflow

1. **Extract snapshots** Use `extract.py` and/or `map_cut.py` to build a snapshot matrix $\mathbf{X}$ (and optionally subtract the mean to obtain $\mathbf{X}'$) for the quantity of interest.

2. **Choose a modal method**
   - **POD** for energy-dominant structures and reduced-order modeling.
   - **SPOD** for statistically stationary flows and frequency-resolved coherent structures (e.g., jet noise, trailing-edge noise).
   - **DMD** for identifying oscillatory features with specific growth/decay rates and frequencies, and for linking to Koopman analysis.

3. **Compute modes**
   - Call the appropriate routine (e.g., `dmd.compute_dmd`, `modal.compute_pod`, `modal.compute_spod`) to obtain modes, eigenvalues/energies, and temporal coefficients.

4. **Analyze and visualize**
   - Plot spatial modes (e.g., contours over cut planes or surfaces).
   - Examine energy spectra or growth-rate vs frequency diagrams.
   - Reconstruct partial fields from selected modes to isolate mechanisms of interest (e.g., vortex shedding, coherent acoustic sources).

5. **Integrate with physics**
   - Interpret modes in terms of known flow phenomena (e.g., von Kármán shedding, Kelvin–Helmholtz rollers, tip vortices).
   - Use SPOD/DMD spectra to connect near-field structures to far-field acoustic measurements or predictions.

---

## 4. References

This toolkit and its mathematical description are closely aligned with the modal-analysis framework in:

- K. Taira *et al.*, “Modal Analysis of Fluid Flows: An Overview,” *AIAA Journal*, Vol. 55, No. 12, 2017, pp. 4013–4041.  
- P. Holmes, J. L. Lumley, G. Berkooz, and C. W. Rowley, *Turbulence, Coherent Structures, Dynamical Systems and Symmetry*, 2nd ed., Cambridge University Press, 2012.  
- J. N. Kutz, S. L. Brunton, B. W. Brunton, and J. L. Proctor, *Dynamic Mode Decomposition: Data-Driven Modeling of Complex Systems*, SIAM, 2016.

These works provide detailed derivations, numerical considerations, and illustrative examples for POD, SPOD, DMD, and related operator-based techniques.
