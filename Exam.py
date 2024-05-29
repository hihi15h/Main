import math


def print_variables(**vars):
  for name, var in vars.items():
    if isinstance(var, float):
      print(f"{name}: {var:.3f}")  # format float with 2 decimal places
    else:
      print(f"{name}: {var}")


def print_precision(variable, nb, precision=3):
  print(f"{variable}: {nb:.{precision}f}")


def quadratic_roots(a, b, c):
  delta = b**2 - 4 * a * c
  if delta < 0:
    return "No real roots"
  elif delta == 0:
    return -b / (2 * a)
  else:
    return [(-b - math.sqrt(delta)) / (2 * a),
            (-b + math.sqrt(delta)) / (2 * a)]


def sin(x):
  return math.sin(math.radians(x))


def asin(x):
  return math.degrees(math.asin(x))


def tan(x):
  return math.tan(math.radians(x))


def atan(x):
  return math.degrees(math.atan(x))


def cos(x):
  return math.cos(math.radians(x))


def num1(ombre_bois, t, v_kmh, t_roche):
  global hauteur_chateau
  l_bois = 4.9 * t_roche**2
  ombre_chateau = v_kmh * t / 3.6
  hauteur_chateau = ombre_chateau * l_bois / ombre_bois
  print_precision("hauteur_chateau", hauteur_chateau)


def num2(f2, dm, hi2, di2, h_chateau):
  global do1
  do2 = 1 / (1 / f2 - 1 / di2)
  hi1 = -hi2 * do2 / di2
  di1 = dm - do2
  do1 = di1 * h_chateau / hi1
  print_precision("distance chateau", do1)


def num3(m_, vf, theta1, d_rampe, h_tour):
  global theta, m
  theta = theta1
  m = m_
  delta_h = h_tour / 2
  e_tot = m / 2 * vf**2 + 9.8 * m * delta_h
  k = 2 * (e_tot) / d_rampe**2
  print_precision("K", k)


def num4(m2, coef, A1, A2, theta_):
  F_tension = m2 * 9.8 * math.sin(math.radians(theta_)) + 9.8 * m2 * math.cos(
      math.radians(theta_)) * coef
  Fx = math.sin(math.radians(A2)) * F_tension
  Ft2 = Fx / math.sin(math.radians(A1))
  print_precision("Ft1", Ft2)
  print_precision("Ft2", F_tension)


def num5(type, R1, R2, Ver2, hi):
  R1 = R1 / 100
  R2 = R2 / 100
  if R2 > R1:
    inter = R1
    R1 = R2
    R2 = inter
  if type == "plan-convexe":
    R2 = math.inf
  elif type == "menisque convergent":
    R1 = -R1
    print("menisque")
  elif type == "biconcave":
    R1 = -R1
    R2 = -R2
  elif type == "plan-concave":
    R1 = -R1
    R2 = math.inf
  elif type == "menisque divergent":
    R2 = -R2
  f1 = 1 / ((1.52 - 1) * (1 / R1 + 1 / R2))
  ftot = 1 / (1 / f1 + Ver2)
  ho = -0.3
  Do = ftot * (1 - ho / hi)
  print_precision("Do", Do)


def num6(h_steeve, d_bord, prof):
  theta = asin((sin(atan(d_bord / h_steeve))) / 1.33)
  long = tan(theta) * prof + d_bord
  h = long / tan(theta)
  reponse = h - h_steeve - prof
  print_variables(theta=theta, long=long, h=h, reponse=reponse)


def num7(portee, l_pente, k, d_ressort, theta):
  global theta_8, vi
  vi = math.sqrt(9.8 * portee / sin(2 * theta))
  em = m * vi**2 / 2
  epr = k * d_ressort**2 / 2
  h = (em - epr) / m / 9.8
  theta_8 = theta
  angle = asin(h / l_pente)
  print_variables(vi=vi, em=em, epr=epr, h=h, angle=angle)


def num8(hi_):
  global dist_x, hi, t2
  hi = hi_
  a = -4.9
  b = vi * sin(theta_8)
  c = -hauteur_chateau / 2 + hi
  print_variables(a=a, b=b, c=c)
  t = quadratic_roots(a, b, c)
  t1 = 0
  t2 = 0
  for time in t:
    if time > t1:
      t2 = time
    else:
      t1 = time
  dist_x = vi * sin(90 - theta_8) * t2
  print_variables(t1=t1, t2=t2, dist_x=dist_x)


def num9():
  global dist_tot
  dist_tot = math.sqrt(dist_x**2 + (hauteur_chateau / 2 - hi)**2)
  t = dist_tot / 340 + t2
  print_variables(dist_tot=dist_tot, t=t)


def num10():
  vfy = math.sqrt((vi * sin(theta_8))**2 - 9.8 * 2 *
                  (hauteur_chateau / 2 - hi))
  vf = math.sqrt(vfy**2 + (vi * sin(90 - theta_8))**2)
  ek = m * vf**2 / 2
  print_variables(vfy=vfy, vf=vf, ek=ek)


def num11(r2):
  v = math.sqrt(6.67 * 10**-11 * 5.9736 * 10**24 / (6380 + r2) / 1000)
  print_variables(v=v)


def num12(cl0, cl1, cl2, cl3, cl4, cl5):
  cl_dict = {"r": 1, "b": 2, "n": 0}
  cl_list = ["vert", "rouge", "bleu", "blanc"]
  pantalons = cl_list[cl_dict[cl0] + cl_dict[cl3]]
  chemise = cl_list[cl_dict[cl1] + cl_dict[cl4]]
  pois = cl_list[cl_dict[cl2] + cl_dict[cl5]]
  print_variables(pantalons=pantalons, chemise=chemise, pois=pois)


def num13(vi, dx1, dt1, vf):
  a = 2 * (dx1 - vi * dt1) / dt1**2
  dx2 = (vf**2 - vi**2) / 2 / a - dx1
  dt = (vf - vi) / a
  dxp = math.sqrt(dx1**2 + dx2**2)
  v = dxp / dt
  angle = atan(dx2 / dx1)
  print_variables(a=a, dx1=dx1, dx2=dx2, dt=dt, dxp=dxp, v=v, angle=angle)


def num14(dx1, v1_, v2_, diff, dt):
  global v1
  v1 = v1_ / 3.6
  v2 = v2_ / 3.6
  vrel1 = v2 - v1
  dt2 = dt - (dx1 - diff) / vrel1
  vrel2 = diff / dt2
  rep = vrel2 - v2 + v1
  print_variables(v1=v1, v2=v2, vrel1=vrel1, dt2=dt2, rep=rep)


def num15(dx, h, lasso, m):
  a = -v1**2 / 2 / dx
  ft = -m * a / cos(asin(h / lasso))
  print_variables(a=a, ft=ft)


print("num 1\n")
num1(ombre_bois=0.6, t=357, v_kmh=4.51, t_roche=0.277)
print("\nnum 2\n")
num2(f2=8.7, dm=17.052, hi2=5, di2=12.18, h_chateau=hauteur_chateau)
print('\nnum 3\n')
num3(m_=4.26, vf=31.8, theta1=24.1, d_rampe=11.4, h_tour=hauteur_chateau)
print("\nnum 4\n")
num4(m2=2.53, coef=0.245, A1=19.7, A2=34.4, theta_=theta)
print("\nnum 5\n")
num5(type="menisque convergent", R1=10.4, R2=24.2, Ver2=-2.83608, hi=2.76)
print("\nnum 6\n")
num6(h_steeve=1.78, d_bord=27.7, prof=1.31)
print("\nnum 7\n")
num7(portee=156.75, l_pente=278.2, k=233, d_ressort=2.7, theta=53)
print("\nnum 8\n")
num8(hi_=634)
print("\nnum 9\n")
num9()
print("\nnum 10\n")
num10()
print("\nnum 11\n")
num11(5825)
print("\nnum 12\n")
num12("n", "n", "r", "b", "n", "b")
print("\nnum 13\n")
num13(7.07, 9242, 1423, 2.19)
print("\nnum 14\n")
num14(170, 9.3, 22.4, 25.3, 41.13)
print("\nnum 15\n")
num15(0.9, 2.31, 4.88, 57.9)
