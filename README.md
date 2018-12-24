# Numerical Methods
Numerical Methods Project - CIN UFPE

Métodos:
* Euler
* Euler Inverso
* Euler Aprimorado
* Runge-Kutta 4
* Adams-Bashforth
* Adams-Moulton
* Fórmula Inversa de Diferenciação

# Bibliotecas:
* SymEnginehttps://github.com/symengine/SymEngine.jl

# entrada.txt:
Métodos de passo único:
* [euler|euler_inverso|euler_aprimorado|runge_kutta] y<sub>0</sub> t<sub>0</sub> h N f(t, y)
Múltiplos passos por lista de valores iniciais:
* adam_bashforth y<sub>0</sub>, y<sub>1</sub>, y<sub>2</sub>, [...], y<sub>K-2</sub> t<sub>0</sub> h N f(t, y) K
* [adam_multon|formula_inversa] y<sub>0</sub>, y<sub>1</sub>, y<sub>2</sub>, [...], y<sub>K-1</sub> t0 h N f(t, y) K

Múltiplos passos obtendo os valores iniciais por métodos anteriores:
* adam_bashforth_by_[euler|euler_inverso|euler_aprimorado|runge_kutta] y<sub>0</sub> t<sub>0</sub> h N f(t, y) K
* adam_multon_by_[euler|euler_inverso|euler_aprimorado|runge_kutta] y<sub>0</sub> t<sub>0</sub> h N f(t, y) K
* formula_inversa_by_[euler|euler_inverso|euler_aprimorado|runge_kutta] y<sub>0</sub> t<sub>0</sub> h N f(t, y) K


f(t, y): funcão diferencial

h: tamanho do passo

N: quantidade de pontos

t<sub>0</sub>: t inicial

y<sub>0</sub>: y inicial

K: ordem de métodos de múltiplos passos
# Execução
entrada.txt:
> runge_kutta 0 0 0.1 20 1-t+4*y

> adam_bashforth 0.0 0.1 0.23 0.402 0.6328 0 0.1 20 1-t+4*y 5

> adam_bashforth_by_euler 0 0 0.1 20 1-t+4*y 6

> adam_multon 0.0 0.1 0.23 0.402 0.6328 0 0.1 20 1-t+4*y 6

> adam_multon_by_euler 0 0 0.1 20 1-t+4*y 6

> formula_inversa 0.0 0.1 0.23 0.402 0.6328 0 0.1 20 1-t+4*y 6

> formula_inversa_by_euler 0 0 0.1 20 1-t+4*y 6

```sh
julia projeto.jl
```
