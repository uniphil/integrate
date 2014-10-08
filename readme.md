Integrate
=========

Some numerical integrators for ordinary differential equations:

 * [Euler's Method](https://en.wikipedia.org/wiki/Euler_method)
 * [Fourth-order Runge-Kutta methods](https://en.wikipedia.org/wiki/Runge–Kutta_methods)

Want one that's not there? Open an issue, or better yet, a pull request adding it.


Fast
----

All of the methods are hand-written in [asm.js](http://asmjs.org), so they should be fast. Not that your integration method would ever be your bottleneck anyway, but hey... why not.


Install
-------

It's on NPM.

```
$ npm install integrate
```


Example
-----

There are two steps: constructing your integrator, and then calling it.

```javascript
var Integrate = require('integrate');

function myODE(t, y) {
  // ordinary differential function to integrate.
  return y;  // for an initial value y=0, this function is also known as 'e^x'
}

// step 1: get the integrators for your function
var integrate = Integrate.Integrator(myODE);

var y = 1,
    t = 0,
    step = 0.25;

while (true) {
  // step 2: integrate!
  console.log('t = ' + t + '   \t', y);
  if ((t += step) && t > 2) break;
  y = integrate.euler(y, t, step);
}

```

Should log this:

```
t = 0        1
t = 0.25     1.25
t = 0.5      1.5625
t = 0.75     1.953125
t = 1        2.44140625
t = 1.25     3.0517578125
t = 1.5      3.814697265625
t = 1.7      4.76837158203125
t = 2        5.9604644775390625
```


### How bad is euler?

Well for starters, that `5.960...` at the end of the example logs should actually be `7.389...`.

How far do we have to shrink our step size to get within `0.0001` of the exact solution? How much better is Runge-Kutta? The fourth-order Runge-Kutta integrator will evaulate your ODE four times for each step, so it has to converge at least 4x faster to be worthwhile...

Here's a quick and dirty test, computing how far we have to shrink the step size for the Euler method and for fourth-order Runge-Kutta, to get within `0.0001` of the exact solution for `t = 2`.

Our ODE, `y' = y` at `y(0) = 1` is actually just `e^x`, so we should be converging to `e^2`.


```javascript
var Integrate = require('integrate');

var integrate = Integrate.Integrator(function(_, y) { return y; });

var stopAt = 2,
    stopError = 0.0001;  // we are done when our error is <= this


function acceptableStep(f, stopAt, stopError) {
  var target = Math.pow(Math.E, stopAt),
      step = stopAt;

  while (true) {
    for (var t=0, y=1; t < stopAt; t+= step)
      y = f(y, t, step);
    if (Math.abs(target - y) <= stopError) break;
    step /= 2;
  }

  return step;
}

var acceptableEuler = acceptableStep(integrate.euler, stopAt, stopError);
var acceptableRk4 = acceptableStep(integrate.rk4, stopAt, stopError);
console.log('euler step with error < '+stopError+' at '+stopAt+': ', acceptableEuler);
console.log('rk4 step with error < '+stopError+' at t='+stopAt+': ', acceptableRk4);
```

Should log

```
euler step with error < 0.0001 at t=2:  0.00000762939453125
rk4 step with error < 0.0001 at t=2:  0.125
```

So, for a cost of 4x more evaluations per step, we get to run with a step size about 16,000x bigger with Runge-Kutta than with the Euler Method for similar accuracy. After our 4x evaluations per step penalty, we are still winning by about 4,000x the number of evaluations required in this example.

So, use rk4.


API
---

### Create an integrator from an ODE function

The nice way:

```javascript
var integrator = Integrate.Integrate(myODEFunction);
```

Shave off one wrapping function call:

```javascript
var integrator = Integrate.ASMIntegrators(null, {f: myODEFunction});
```

Both forms will return an identical object.

`myODEFunction` should accept two parameters: `t`, and `y`.


### Integrate with one of the numerical integration methods

* `euler` and `rk4` may both be called with three parameters: `y`, `t`, and `step`.

* `rk4general` takes a fourth parameter, `λ`, which should be an integer between 1 and 5, inclusive. See wikipedia: https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#The_Runge.E2.80.93Kutta_method

```
var yNext = integrator.euler(yLast, tNow, tStepSize);

var yNext = integrator.rk4(yLast, tNow, tStepSize);

var yNext = integrator.rk4general(yLast, tNow, 3);
```
