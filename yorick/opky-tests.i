/*
 * opky-tests.i --
 *
 * Tests for OPKY, the Yorick interface to OptimPack.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (c) 2014 Éric Thiébaut
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 *
 *-----------------------------------------------------------------------------
 */

require, "opky.i";

func opkt_rosenbrock_xinit(x0)
{
  x0(1::2) = -1.2;
  x0(2::2) =  1.0;
}

func opkt_rosenbrock_fg(x, gx)
{
  type = structof(x);
  c1 = type(1);
  c2 = type(2);
  c10 = type(10);
  c200 = type(200);
  r1 = 1::2;
  r2 = 2::2;
  x1 = x(r1);
  x2 = x(r2);
  t1 = c1 - x1;
  t2 = c10*(x2 - x1*x1);
  g2 = c200*(x2 - x1*x1);
  gx(r1) = -c2*(x1*g2 + t1);
  gx(r2) = g2;
  return sum(t1*t1) + sum(t2*t2);
}

func opkt_rosenbrock(n, single=, vmlm=) {
  type = (single ? float : double);
  x0 = array(type, n);
  opkt_rosenbrock_xinit, x0;
  x = opk_minimize(opkt_rosenbrock_fg, x0, vmlm=vmlm, verb=1n);
}

func opkt_nonnegative(m, n, mem=, single=, scale=, verb=, maxiter=, maxeval=)
{
  type = (single ? float : double);
  H = type(random_n(m, n));
  x0 = type(max(random_n(n), 0.0));
  std = 0.1; // noise level
  w = array(type(1.0/std^2), m);
  y = H(,+)*x0(+) + type(std*random_n(m));
  A = H(+,)*(w*H)(+,);
  b = H(+,)*(w*y)(+);
  mem = 5;
  opt = opk_vmlmc(n, mem, 1.0 /*0.1*sqrt(n)*/, single=single);
  task = opk_start(opt);
  x = array(type, n);
  gx = array(type, n);
  d = array(type, n);
  fx = type(0.0);
  for (;;) {
    if (task == OPK_TASK_PROJECT_X) {
      x = max(x, 0.0);
    } else if (task == OPK_TASK_COMPUTE_FG) {
      r = H(,+)*x(+) - y;
      wr = w*r;
      fx = 0.5*sum(wr*r);
      gx = H(+,)*wr(+);
      free_vars = double((gx < 0.0)|(x > 0.0));
    } else if (task == OPK_TASK_PROJECT_D) {
      d *= free_vars;
      dnorm = sqrt(sum(d*d));
    } else if (task == OPK_TASK_NEW_X || task == OPK_TASK_FINAL_X) {
      iter = opt.iterations;
      eval = opt.evaluations;
      if (task == OPK_TASK_FINAL_X) {
        last = 1;
      } else if (! is_void(maxiter) && iter >= maxiter) {
        last = 2;
      } else if (! is_void(maxeval) && eval >= maxeval) {
        last = 3;
      } else {
        last = 0;
      }
      if (verb) {
        if (iter == 0) {
          write, format="%s\n%s\n",
            "  ITER.  EVAL.  REST.  PROJ.          F(X)             ||PG(X)||",
            "----------------------------------------------------------------";
        }
        if (last || (iter%verb) == 0) {
          write, format="%6d %6d %6d %6d  %+24.15E  %10.3E  %10.3E\n",
            iter, eval, opt.restarts, opt.projections, fx,
            sqrt(sum(d*d)), dnorm;
        }
      }
      if (last != 0) {
        if (last == 2) {
          write,
            format="WARNING: maximum number of iterations exceeded (%d)\n",
            iter;
        } else if (last == 3) {
          write,
            format="WARNING: maximum number of evaluations exceeded (%d)\n",
            eval;
        }
        break;
      }
    } else if (task == OPK_TASK_WARNING) {
      write, format="WARNING: %s\n", "some warning occured";
      break;
    } else if (task == OPK_TASK_ERROR) {
      error, "some error occured";
    } else {
      error, "unexpected task";
    }
    task = opk_iterate(opt, x, fx, gx, d);
  }
  window, 20;
  fma;
  limits, square=0;
  limits;
  logxy, 0, 0;
  plp, x0, color="blue", symbol='#';
  plp, x, color="red", symbol='x';
  pause, 1;

  if (! is_func(op_mnb)) {
    // Try to load former version of OptimPack.
    include, "OptimPack1.i", 3;
  }
  if (is_func(op_mnb)) {
    // Use former version of OptimPack.
    ws = h_save(H, y, w);
    x1 = array(type, n);
    x1 = op_mnb(_opkt_nonnegative, x1, extra=ws, xmin=type(0),
                mem=mem, verb=1);
    plp, x1, color="cyan", symbol='+';
  }

}

func _opkt_nonnegative(x,&gx,ws)
{
  r = ws.H(,+)*x(+) - ws.y;
  wr = ws.w*r;
  gx = ws.H(+,)*wr(+);
  return 0.5*sum(wr*r);
}

/*
 * Local Variables:
 * mode: Yorick
 * tab-width: 8
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * fill-column: 79
 * coding: utf-8
 * ispell-local-dictionary: "american"
 * End:
 */
