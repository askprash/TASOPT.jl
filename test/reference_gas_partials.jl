# Array-based reference implementations of the pressure-ratio and enthalpy-change
# processes with their partials, gas_pratd and gas_delhd, as they were in
# src/engine/gascalc.jl before their removal (moved verbatim; nothing in
# src used them). They evaluate properties with TASOPT's current gassum,
# gassumd and gasfun, so they give the same numbers as the NamedTuple partials
# functions and the composition-direction helpers that replaced them in
# tfoper!, and unit_test_thermo_equivalence.jl checks that.
module ReferenceGasPartials
import TASOPT.engine: gassum, gassumd, gasfun

function gas_pratd(alpha, n, po, to, ho, so, cpo, ro, pratio, epol)


      itmax = 10
      ttol = 0.000001

      gexp = ro / (cpo * epol)
      t = to * pratio^gexp

      pile = log(pratio) / epol

      t_al = zeros(n)
      h_al = zeros(n)
      s_al = zeros(n)
      p_al = zeros(n)
      cp_al = zeros(n)
      r_al = zeros(n)

      dt = 0.0
      for iter = 1:itmax
            s, s_t, h, h_t, cp, r = gassum(alpha, n, t)

            #------ res( s(t,al), r[al], so, pratio )
            res = (s - so) / r - pile
            res_t = s_t / r

            dt = -res / res_t

            #c        write(*,*) iter, t, dt

            t = t + dt

            if (abs(dt) < ttol)

                  s, s_t, h, h_t, cp, cp_t, r = gassumd(alpha, n, t)
                  res_t = s_t / r

                  #---- res( so, pratio, al )
                  res_so = -1.0 / r
                  res_pi = -1.0 / (epol * pratio)
                  res_ep = pile / epol
                  #c    res_al =  si/r - (s-so)/r^2 * ri

                  t_so = -res_so / res_t
                  t_pi = -res_pi / res_t
                  t_ep = -res_ep / res_t
                  #c    t_al = -res_al/res_t 

                  #c          t(so,pi,ep,al)
                  #c       h( t(so,pi,ep,al),al )
                  #c       s( t(so,pi,ep,al),al )
                  #c      cp( t(so,pi,ep,al),al )
                  #c       r( al )
                  #c    

                  #---- p( po, so, pi )
                  p = po * pratio
                  p_so = 0.0
                  p_po = pratio
                  p_pi = po
                  p_ep = 0.0

                  #---- h( t(so,pi,ep,al),al )
                  h_so = h_t * t_so
                  h_pi = h_t * t_pi
                  h_ep = h_t * t_ep
                  #c    h_al = h_t*t_al + hi

                  #---- s( t(so,pi,ep,al),al )
                  s_so = s_t * t_so
                  s_pi = s_t * t_pi
                  s_ep = s_t * t_ep
                  #c    s_al = s_t*t_al + si

                  for i = 1:n
                        si, s_ti, hi, h_ti, cpi, ri = gasfun(i, t)
                        res_al = si / r - (s - so) / r^2 * ri
                        t_al[i] = -res_al / res_t
                        h_al[i] = h_t * (-res_al / res_t) + hi
                        s_al[i] = s_t * (-res_al / res_t) + si
                        p_al[i] = 0.0
                        cp_al[i] = cp_t * (-res_al / res_t) + cpi
                        r_al[i] = ri
                  end

                  return p, t, h, s, cp, r, p_po, p_so, t_so, h_so, s_so, p_pi, t_pi, h_pi, s_pi,
                  p_ep, t_ep, h_ep, s_ep,
                  p_al, t_al, h_al, s_al, cp_al, r_al
            end

      end
      println("gas_pratd: convergence failed.  dT =", dt)
      println("To pratio =", to, pratio)

end # gas_pratd

function gas_delhd(alpha, n, po, to, ho, so, cpo, ro, delh, epol)

      itmax = 15
      ttol = 0.000001


      t = to + delh / cpo
      t = max(t, 0.2 * to)

      t_al = zeros(n)
      h_al = zeros(n)
      s_al = zeros(n)
      p_al = zeros(n)
      cp_al = zeros(n)
      r_al = zeros(n)

      dt = 0.0
      for iter = 1:itmax

            s, s_t, h, h_t, cp, r = gassum(alpha, n, t)

            res = h - ho - delh
            res_t = h_t

            dt = -res / res_t

            rlx = 1.0
            if (rlx * dt < -0.4 * t)
                  rlx = -0.4 * t / dt
            end
            if (rlx * dt > 0.6 * t)
                  rlx = 0.6 * t / dt
            end

            t = t + rlx * dt

            if (abs(dt) < ttol)
                  s, s_t, h, h_t, cp, cp_t, r = gassumd(alpha, n, t)
                  res_t = h_t

                  res_ho = -1.0
                  res_dh = -1.0

                  t_ho = -res_ho / res_t
                  t_dh = -res_dh / res_t

                  p = po * exp(epol * (s - so) / r)
                  p_r = -p * epol * (s - so) / r^2
                  p_s = p * epol / r
                  p_so = -p * epol / r
                  p_po = exp(epol * (s - so) / r)
                  p_ep = p * (s - so) / r

                  p_t = p_s * s_t
                  p_ho = p_t * t_ho
                  p_dh = p_t * t_dh

                  h_ho = h_t * t_ho
                  h_dh = h_t * t_dh

                  s_ho = s_t * t_ho
                  s_dh = s_t * t_dh

                  for i = 1:n
                        si, s_ti, hi, h_ti, cpi, ri = gasfun(i, t)
                        res_al = hi
                        t_al[i] = -res_al / res_t
                        p_al[i] = p_t * (-res_al / res_t) + p_s * si + p_r * ri
                        h_al[i] = h_t * (-res_al / res_t) + hi
                        s_al[i] = s_t * (-res_al / res_t) + si
                        cp_al[i] = cp_t * (-res_al / res_t) + cpi
                        r_al[i] = ri
                  end

                  return p, t, h, s, cp, r, p_so, p_po, p_ep,
                  p_ho, t_ho, h_ho, s_ho, p_dh,
                  t_dh, h_dh, s_dh, p_al, t_al,
                  h_al, s_al, cp_al, r_al
            end

      end
      println("gas_delh: convergence failed.  dT =", dt)

end # gas_delhd

end # module
