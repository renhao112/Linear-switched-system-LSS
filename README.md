# Linear-switched-system-LSS
The code for the second numerical example of the paper "On Krylov subspace implementation of ParaExp for linear switched systems" subjected to IEEE Control System Letters
## A large-scale linear switched system constructed from the 2D convection-diffusion equation
Consider the convection-diffusion equation

$$ u_t - \nu_{\sigma(t)}\Delta u + \boldsymbol{b} \nabla u=f_{\sigma(t)} \quad \text{in}\ [0,1]^2\times[0,T] $$

with zeros initial value condition and homogeneous Dirichlet bounday condition.
The diffusion parameter $\nu_{\sigma(t)}$ and the source term $f_{\sigma(t)}$ switch upon the time. Set $\sigma(t)=i, t\in[T_{i-1},T_i]$ with $T_0=0, T_5=T, T_i=2*i, i=1,\dots,5$. 
Let $\nu_1 = 0.06, \nu_2=0.04,\nu_3=0.02,\nu_4=0.08,\nu_5=0.1$, making it a convection-dominated problem, then central difference yields the LSS of the form

$$\dot{x} = A_ix(t)+B_i\tilde{u}(t), \quad t\in[T_{i-1},T_i],$$

where $A_i\in R^{n \times n}$ is the difference matrix with $n$ the number of the inner grid point. $B_i\in R^{n \times 4}$ comes from the discretization of $f_i$, however, it 
can be artificially given, seeing the code 'generate_lss.m', while the input $\tilde{u}(t)$ is given by 'ipt.m'.

For this example, we mainly interested in the algorithm's performance at different scales, the only input for generate_lss.m (namely generating the LSS) is the difference parameter
 $N_x$. The number of grid of each dimension is $N_x+2$, so the size of the LSS is $n=N_x \times N_x$.

 ## The description of each Matlab file
| file name | description |
| :-- | :--: |
| generate_lss | generate the LSS|
| ipt | input $\tilde{u}(t)$ |
| main_CD | the main function for the proposed ParaExp algorithm |
| EBK | the Krylov subspace method for inhomogeneous problem with zero initial value |
| SAI_appro | the Krylov subspace method for homogeneous problem |
| EBK_s | the Krylov subspace method for inhomogeneous problem with nonzero initial value|
| TR | the trapezoidal rule for inhomogeneous problem with zero initial value |
| TR_s | the trapezoidal rule for inhomogeneous problem with nonzero initial value |
| main_s | get the serial solution and the serial CPU time for the corresponding ParaExp algorithm |
 
 
