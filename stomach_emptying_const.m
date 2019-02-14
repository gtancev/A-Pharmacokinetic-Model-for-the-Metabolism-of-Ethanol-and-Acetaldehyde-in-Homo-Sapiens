function k_s = stomach_emptying_const(k_s_max,D,a)

k_s = k_s_max/(1 + a*(D/1000)^2);

end