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

func opkt_rosenbrock(n, single=) {
  type = (single ? float : double);
  x0 = array(type, n);
  opkt_rosenbrock_xinit, x0;
  x = opk_minimize(opkt_rosenbrock_fg, x0, verb=1n);
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
