# JUelich Superconducting QUAntum Computer Emulator - JUSQUACE

## Table of contents

<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [JUelich Superconducting QUAntum Computer Emulator - JUSQUACE](#juelich-superconducting-quantum-computer-emulator-jusquace)
	- [Table of contents](#table-of-contents)
	- [Overview](#overview)
	- [How to clone and compile](#how-to-clone-and-compile)
	- [How to test the code](#how-to-test-the-code)
	- [Organization of the state vector](#organization-of-the-state-vector)
	- [Where to sample from the state vector](#where-to-sample-from-the-state-vector)
	- [The model Hamiltonian and control functions for single- and two-qubit operations](#the-model-hamiltonian-and-control-functions-for-single-and-two-qubit-operations)
	- [How to define a system and control pulses](#how-to-define-a-system-and-control-pulses)
	- [Examples](#examples)
		- [Spectroscopy](#spectroscopy)
		- [Eigenenergies](#eigenenergies)
	- [References](#references)

<!-- /TOC -->


## Overview

This repository contains the Julich Superconducting QUAntum Computer Emulator (JUSQUACE) which can be used to study the real-time dynamics of interacting time-dependent anharmonic oscillators (transmons) and time-independent harmonic oscillators (LC resonators). The interactions themselves are of the dipole-dipole type. The examples discussed in this README file can be reproduced by compiling and executing the files `time_evol.cpp` and `spectrum.cpp`. This README file also contains explanations on how to modify the device parameters and time-dependent functions to adapt the code to the user's needs. The simulation algorithm is based on the second-order product-formula algorithm to solve the time-dependent Schrödinger
equation (TDSE), see Refs. [[3-4]](#3),  for more details.  

## How to clone and compile

Clone the repository:

```bash
git clone https://jugit.fz-juelich.de/qip/jusquace
cd jusquace
```

Compile the programs with:

```bash
make time_evo
```
and
```bash
make spectrum
```

(tested with intel/19.0.2 on Linux)

## How to test the code

Verify that the programs runing `./main 0 &> res_map.dat &` (for `time_evol.cpp`) and `./main` (for `spectrum.cpp`) . You should find a file `.dat` files in your working directory. If so, you can use the scripts `.plt` to display the data. The corresponding figures should look like the ones you can find on the bottom of this page.

## Organization of the state vector
The state vector can be expressed as

```math
\begin{equation}
|\psi\rangle=\sum\limits_z c_{z} |z\rangle \label{eq1}\tag{1}
\end{equation}
```

and the corresponding basis states are tensor product states of the form

```math
\begin{equation}
|z\rangle = \otimes_{i} |n_{i}\rangle \otimes_{j} |n_{j}\rangle \label{eq2}\tag{2},
\end{equation}
```

where $\{|n_{j}\rangle\}$ denotes the eigenstates of transmon j and $\{|n_{i}\rangle\}$ represents the eigenstates of resonator j. Note that both $n_{j}$ and $n_{i}$ are elements of the index set $I_{3}=\{0,1,2,3\}$ since we model both, transmons and resonators with 4 states only. Please take into account that this setting is not valid for all physical scenarios of interest, see for example Section III A in Ref. [[1]](#1). The complex coefficients $c_{z}$ are stored as
```bash
Re(c_{z})=state[2*z+0]
Im(c_{z})=state[2*z+1]
```
in the main data structure `double *state` of size $= 2*4^{N}$, where $N=N_{T}+N_{R}$ and $N_{T}=$#Transmons and $N_{R}=$#Resonators.

Since the basis states are defined as tensor product states, we will discuss how to sample from the state vector. Assume that we simulate two qubits. Furthermore, we would like to determine the probabilities $`p(m_{1},m_{0})`$ for $`m_{1},m_{0} \in I_{3}`$ , this would mean we have to consider the state vector components

```bash
offset = 4
Re(n_{1},n_{0}) = state[2*(n_{0}+offset*n_{1})]
Im(n_{1},n_{0}) = state[2*(n_{0}+offset*n_{1})+1]
p(n_{1},n_{0})=Re(n_{1},n_{0})*Re(n_{1},n_{0})+Im(n_{1},n_{0})*Im(n_{1},n_{0}).
```
Note that we have to set the offset to four since we use the four lowest states to model the subsystems. More generally, we can use the function

```math
\begin{equation}
idx(|z\rangle) = \sum_{j=0}^{N_{T}-1} \mathrm{offset}^{j} n_{j} + \sum_{i=0}^{N_{R}-1} \mathrm{offset}^{N_{T}+i} n_{i} \label{eq3}\tag{3},
\end{equation}
```

to determine the index $idx(|z\rangle)$ for an arbitrary state vector component $c_{z}$.

## Where to sample from the state vector
We would recommend to sample from within the algorithm itself. The corresponding place is marked with comments.
```bash
pfa_time_evolution(double *state, const vector<uint64_t> tar_idx)
{
    for(unsigned int t=0;t<T;t++)
    {
        //main code
        => sample here from the state vector
    }
}
```
## The model Hamiltonian and control functions for single- and two-qubit operations
ALL INPUT PARAMETERS NEED TO BE IN UNITS OF [1/ns] !

The model Hamiltonian is defined as

```math
\begin{equation}
\hat{H}= \hat{H}_{\text{Res.}} + \hat{H}_{\text{Trans.}} + \hat{\mathcal{D}}_{\text{Flux.}} + \hat{\mathcal{D}}_{\text{Charge.}} + \hat{V}_{\text{Int.}}. \label{eq4}\tag{4}
\end{equation}
```

The first term

```math
\begin{equation}
\hat{H}_{\text{Res.}}= \sum_{i} \omega_{i}^{(R)} \hat{a}_{i}^{\dagger}\hat{a}_{i},\label{eq5}\tag{5}
\end{equation}
```

describes a collection of non-interacting time-independent harmonic oscillators (LC resonators). We simulate this Hamiltonian in a bare basis formed by the tensor product states of the harmonic oscillator $|n_{i} \rangle= | \psi_{i}^{(n)} \rangle$. The device parameter $\omega_{i}^{(R)}$ denotes the resonator frequency of the i-th resonator element. Similarly, the operators $\hat{a}_{i}^{\dagger}$ and $\hat{a}_{i}$ refer to the creation and annihilation operators of the i-th resonator element.

The second term

```math
\begin{equation}
\hat{H}_{\text{Trans.}}= \sum_{j} \left( \omega_{j}^{(q)}(t) \hat{b}_{j}^{\dagger}\hat{b}_{j} + \frac{\alpha_{j}^{(q)}(t)}{2}  \hat{b}_{j}^{\dagger}\hat{b}_{j} (\hat{b}_{j}^{\dagger}\hat{b}_{j} - \hat{I})\right)  \label{eq6}\tag{6}
\end{equation}
```

describes a collection of non-interacting, time-dependent, adiabatic anharmonic oscillators (transmons). In our model the adiabatic anharmonic oscillators represent the transmon elements or qubits. We simulate this Hamiltonian in a bare basis formed by the tensor product states of the time-dependent harmonic oscillator $|n_{j} \rangle= | \psi_{j}^{(n)}(t) \rangle$. The Hamiltonian is defined in terms of the creation and annihilation operators $\hat{b}_{j}^{\dagger}$ and $\hat{b}_{j}$. Both operators are defined with regard to the time-dependent harmonic oscillator eigenstates $|n_{j} \rangle$ in the same manner as the operators $\hat{a}_{i}^{\dagger}$ and $\hat{a}_{i}$. The functions $\omega_{j}^{(q)}(t)$ and $\alpha_{j}^{(q)}(t)$ are used to model the instantaneous energies $(E_{n}(t)- E_{0}(t))=n \alpha_{j}^{(q)}(t) + (\alpha_{j}^{(q)}(t)/2) n(n-1)$ of the j-th transmon for all times $t$. The simulation code allows us to choose between a low (high) accuracy first-order (high-order) approximation, i.e. we can choose between two different sets of functions $\omega_{j}^{(q)}(t)$ and $\alpha_{j}^{(q)}(t)$. Both approximations aim to approximate the energies of a specific circuit Hamiltonian, see Appendix B in Ref. [[1]](#1). Furthermore a derivation of Eq. ([6](#mjx-eqn-eq6)) can be found in Section II B of Ref. [[1]](#1).   

The third term

```math
\begin{equation}
\hat{\mathcal{D}}_{\text{Flux.}} = \sum_{j} \left( -i \sqrt{\frac{\xi_{j}(t)}{2}} \dot{\varphi}_{\text{eff.,j}}(t) (\hat{b}_{j}^{\dagger} -\hat{b}_{j}) + \frac{i}{4}\frac{\dot{\xi}_{j}(t)}{\xi_{j}(t)} (\hat{b}_{j}^{\dagger}\hat{b}_{j}^{\dagger} -\hat{b}_{j}\hat{b}_{j})\right)\label{eq7}\tag{7}
\end{equation}
```

is a non-linear flux driving term which results from the fact that we describe the transmons in a time-dependent basis with the basis states $| \psi_{j}^{(n)}(t) \rangle$. In Section II B of Ref. [[1]](#1) the origin of this term is discussed in more detail.

The fourth term

```math
\begin{equation}
\hat{\mathcal{D}}_{\text{Charge.}} = \sum_{j} \Omega(t) (\hat{b}_{j}^{\dagger} + \hat{b}_{j}) \label{eq8}\tag{8}
\end{equation}
```

is a linear charge driving term which describes a voltages sources capacitively coupled to the transmon qubits, see Appendix A and B in Ref. [[2]](#2).

The fifth term

```math
\begin{equation}
\hat{V}_{\text{Int.}} = \sum_{i,j} g_{i,j}(t)  \left((\hat{a}_{i}^{\dagger} + \hat{a}_{i}) \otimes(\hat{b}_{j}^{\dagger} + \hat{b}_{j}) \right)\label{eq9}\tag{9}
\end{equation}
```

describes time-dependent dipole-dipole interactions between the transmon and resonator elements. The time dependence results from the fact that we describe the transmons in a time-dependent basis, see Section II in Ref. [[1]](#1).   

## How to define a system and control pulses

The file `chip_library.cpp` contains a subroutine which can be used to define a device model for a system with the following parameters, see also Ref. [[2]](#2).

| $i$  | $\omega_{i}/2 \pi$ |$\omega_{i}/2 \pi$ | $\alpha_{i}/2 \pi$ | $G_{2,i}/2 \pi$  | $\phi_{0,i}/\pi$  |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| 0  | 0.0  | 4.2  | -0.320  | 0.300  | 0.0  |
| 1  | 0.0  | 5.2  | -0.295  | 0.300  | 0.0  |
| 2  | 45.0  | 0.0  | 0.0  | 0.000  | 0.0  |

## Examples

### Spectroscopy

We conjecture that the resonance frequency (for rabi oscillations between the computational state $`|00\rangle`$ and $`|01\rangle`$) is located around the frequency $`\omega^{(D)}=4.2`$ [GHz]. Therefore we would like to perform a spectroscopy of a small frequency band. If we initialize the system in the state $`|\psi\rangle=|00\rangle`$, we find the result:

<a href="/res_map.png">
    <p align="center">
    <img src="/res_map.eps" width="600">
    </p>
</a>

You should be able to reproduce this result with `./main 0 &> res_map.dat &`, after executing the command `make time_evo`.

### Eigenenergies

You should be able to reproduce this result with `./main.cpp example1`.

<a href="/res_map.png">
    <p align="center">
    <img src="/res_map.eps" width="600">
    </p>
</a>

You should be able to reproduce this result with `./main`, after executing the command `make spectrum`.

## References

<a name="1">[1]</a> H. Lagemann et al., "Numerical analysis of effective models for flux-tunable transmon systems", [Phys. Rev. A 106, 022615](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.106.022615) (2022)

<a name="2">[2]</a> H. Lagemann et al., "On the fragility of gate-error metrics in simulation models of flux-tunable transmon
quantum computers", [Arxiv](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.106.022615) (2022)

<a name="3">[3]</a> H. De Raedt, "Product formula algorithms for solving the time dependent Schrödinger equation", [*Comput. Phys. Rep.* **7**, 1](https://doi.org/10.1016/0167-7977(87)90002-5) (1987)

<a name="4">[4]</a> H. Lagemann , "Real-time simulations of transmon systems with
time-dependent Hamiltonian models", [Arxiv](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.106.022615) (2022)
