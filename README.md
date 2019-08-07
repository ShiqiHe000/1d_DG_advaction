# 1d_DG_advaction
The DG solver solves an one dimensional scaler advection problem. The domain is one element from -1 to 1. 

## Problem 
We approximate a one dimensional scalar advaction equation.

<a href="https://www.codecogs.com/eqnedit.php?latex=\left&space;\{\begin{matrix}&space;\varphi_{t}&plus;c&space;\varphi_{x}=0,&space;-1<x<1,&space;\\&space;\varphi(x,&space;0)=\varphi_{0}(x),&space;-1\leqslant&space;x\leqslant&space;1,&space;\\&space;\varphi(-1,&space;t)=g(t),&space;t>0.&space;\end{matrix}\right." target="_blank"><img src="https://latex.codecogs.com/gif.latex?\left&space;\{\begin{matrix}&space;\varphi_{t}&plus;c&space;\varphi_{x}=0,&space;-1<x<1,&space;\\&space;\varphi(x,&space;0)=\varphi_{0}(x),&space;-1\leqslant&space;x\leqslant&space;1,&space;\\&space;\varphi(-1,&space;t)=g(t),&space;t>0.&space;\end{matrix}\right." title="\left \{\begin{matrix} \varphi_{t}+c \varphi_{x}=0, -1<x<1, \\ \varphi(x, 0)=\varphi_{0}(x), -1\leqslant x\leqslant 1, \\ \varphi(-1, t)=g(t), t>0. \end{matrix}\right." /></a>

Wave speed : <a href="https://www.codecogs.com/eqnedit.php?latex=c" target="_blank"><img src="https://latex.codecogs.com/gif.latex?c" title="c" /></a>

If <a href="https://www.codecogs.com/eqnedit.php?latex=c>0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?c>0" title="c>0" /></a>, then the boundary condition is imposed at the left, and the solution is interpolated to the right. Vise a versa. 

## Implementation
The solver solves the problem with one element using polynomial with degree N. 

In order to obtain the time derivative, we first compute the boundary condition according to the sign of the wavespeed. To compute the boundary, it uses the user supplied function g(t). The boundary solution at the other boundary is evaluated by using a dot produce of the array solution vectors and the array of Lagrange interpolating polynomials. With the boundary values we can calculate the space derivative. The last step is to get the time derivative by multiplying the space derivative with the wavespeed. 

The next step is to integrate the time derivative in time. We choose third order Runge-Kutta to integrate in time. 

## Benchmark Solution
We solve the problem with initial condition <a href="https://www.codecogs.com/eqnedit.php?latex=\varphi_{0}(x)=e^{-ln(2)(x&plus;1)^2/\sigma&space;^2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\varphi_{0}(x)=e^{-ln(2)(x&plus;1)^2/\sigma&space;^2}" title="\varphi_{0}(x)=e^{-ln(2)(x+1)^2/\sigma ^2}" /></a>

We choose <a href="https://www.codecogs.com/eqnedit.php?latex=\triangle&space;t=1.5\times&space;10^{-4}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\triangle&space;t=1.5\times&space;10^{-4}" title="\triangle t=1.5\times 10^{-4}" /></a> to ensure the temporal errors were small relative to the spatial errors. 

### Error Perfomance
![alt text](https://github.com/ShiqiHe000/1d_DG_advaction/blob/master/advection_error.png)
