function ASMIntegrators(_, foreign) {
  'use asm';

  var f = foreign.f;

  function euler(y0, t0, h) {
    // Follow tangents on the curve over short distances.
    // You probably should use rk4 instead, which converges faster.
    // https://en.wikipedia.org/wiki/Numerical_methods_for_ordinary_differential_equations#Euler_method
    y0 = +y0;      // double
    t0 = +t0;      // double
    h = +h;        // double
    var yh = 0.0;  // -> double

    yh = y0 + h*+f(t0, y0);

    return +yh;
  }

  function rk4(y0, t0, h) {
    // Compute the next value by taking a weighted average of four samples of f
    // https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#The_Runge.E2.80.93Kutta_method
    y0 = +y0;      // double
    t0 = +t0;      // double
    h = +h;        // double
    var k1 = 0.0,  // double
        k2 = 0.0,  // double
        k3 = 0.0,  // double
        k4 = 0.0,  // double
        yh = 0.0;  // -> double

    k1 = +f(+(t0 )       , +(y0 ));
    k2 = +f(+(t0 + h/2.0), +(y0 + h/2.0*k1));
    k3 = +f(+(t0 + h/2.0), +(y0 + h/2.0*k2));
    k4 = +f(+(t0 + h)    , +(y0 + h *   k3));

    yh = y0 + h/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);

    return +yh;
  }

  function rk4general(y0, t0, h, Λ) {
    // Generalized RK4, where Λ can by 1, 3, 4, or 5. Λ=2 is the classical rk4.
    // the last part of:
    // https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#The_Runge.E2.80.93Kutta_method
    y0 = +y0;      // double
    t0 = +t0;      // double
    h = +h;        // double
    Λ = Λ|0;       // int -- should be in [1, 5], but this is not validated
    var λ  = 0.0,
        k1 = 0.0,  // double
        k2 = 0.0,  // double
        k3 = 0.0,  // double
        k4 = 0.0,  // double
        yh = 0.0;  // -> double

    λ = +~Λ;  // convert int to double

    k1 = +f(+(t0 )       , +(y0 ));
    k2 = +f(+(t0 + h/2.0), +(y0 + h/2.0*k1));
    k3 = +f(+(t0 + h/2.0), +(y0 + (0.5-1.0/λ)*k1*h + 1.0/λ*k2*h));
    k4 = +f(+(t0 + h)    , +(y0 + (1.0-λ/2.0)*k2*h + λ/2.0*k3*h));

    yh = y0 + h/6.0*(k1 + (4.0-λ)*k2 + λ*k3 + k4);

    return +yh;
  }


  return {
    euler: euler,
    rk4: rk4,
    rk4general: rk4general
  };
}


function Integrator(f) {
  // nicer wrapper around the asm.js constructor
  return ASMIntegrators(null, {f: f});
}


module.exports = {
  ASMIntegrators: ASMIntegrators,
  Integrator: Integrator
};
