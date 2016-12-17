## Newton descent across the Einstein Field Equation constraints

solve the metric using the norm of G_ab = 8 pi T_ab
This is gonna be like my EinsteinFieldEquationSolution project, but using OpenCL
And now I'm going to do the legwork and finish the math to calculate the gradient descent directly, rather than numerically approximating it.
And I'm adding lots of metaprogramming, which means Lua markup, which is more difficult to write out in C++.
Which means I'm skipping the C++.  Hooray!

### How it works:

1. provide initial stress-energy conditions, in terms of primitives (not T_ab, because that depends on g_ab)
2. provide initial metric primitivies
3. solve the quadratic function G_ab = 8 pi T_ab, where G_ab is derived from the metric prims, and T_ab is derived from the stress-energy (and metric) prims

### Requires:

* luajit
* malkia's ufo
* cimgui
* lua-ext
* lua-vec
* lua-cl
* my ffi bindings (which include my ffi vector classes)
* template-lua
* symmath-lua

### Math work is here:

https://cdn.rawgit.com/thenumbernine/efesoln-cl-lua/master/efe.html

### Progress

For the spherical body problem it is converging slowly (1e-30) in flat space initial conditions,
but once it starts to develop surface features it converges very quickly.
I'm thinking I should use a line search to fix this.  Maybe Hessian as well, but that'd be complex, so I'll put it off until my automatic symbolic generation (which depends on me simplifying TensorRef * TensorRef in my symmath-lua)

Now I'm allowing update to alpha, beta, gamma separately, and testing it with alpha update only.  Goes much faster.
What happens is alpha varies greatly around the matter boundary (like it should), but then the large alphas get larger and small gets smaller, and the coordinates grow too much.
Maybe I should restrict alpha and det gamma to unit volumes?  Or only converge the linearized metric, leaving the background?

Now, using Newton method, I've messed with manually adjusting the line search.  Starting at flat space, wiggling the gradient back and forth, I can get it converging down towards 1e-6, but this is still nothing compared to initialization in stellar Schwarzschild, which is 1e-40.
It appears to be converging to a local minima formed by the discretized spatial differentiation calculations used by the curvature.
Starting at stellar Schwarzschild and converging to a minimia causes it to soon diverge as well.
Taking smaller steps shows that stepping in either direction along the gradient quickly walks out of the valley containing the minima.

This is probably due to the fact that my Newton descent is currently implemented based on the gradient of the continuous Einstein field equation, not the discretized.
The alpha values at the matter boundary are extremizing, but staying zero elsewhere where the matter gradient is zero.
This is probably because the gradient of the 1st kind Christoffel symbols is equal to zero for the continuous equation,
whereas the gradient of the discrete is equal to the difference of the left and right metrics.
If the discrete gradient was used, the alpha values might get propogated outward from the matter boundary. 

### TODO:

* find out why the conjugate residual solver is stopping after one iteration
* implement GMRES in OpenCL
* implement JFNK in OpenCL
* reduce allocations for solvers that aren't being used
* upon changing bodies, make sure the ffi.cdef code and headers are properly refreshing 

![](https://cdn.rawgit.com/thenumbernine/efesoln-cl-lua/master/images/pic.png)
