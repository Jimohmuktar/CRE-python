# Numerical methods analysis
def numerical_rules(param, h):
    n = len(param)
    f_0 = param[0]
    f_n = param[-1]
    sum_mid = sum(param[1:-1])
    sum_even = sum(param[1:-1:2])
    sum_odd = sum(param[2:-1:2])
    sum_mult3 = sum(param[2:-1:3])
    sum_not_mult3 = sum(f for i, f in enumerate(param[1:-1], 1) if i%3 != 0)
    
    if n % 2 != 0 and n % 3 != 0:               #Trapezoidal rule
        soln = (h/2)*(f_0 + 2*(sum_mid) - f_n)
        rule = "Trapezoidal"
     
    elif n % 2 == 0 and n % 3 != 0:             # Simpson's one-third rule
        soln = (h/3)*(f_0 + 4*sum_even + 2*sum_odd + f_n)
        rule = "Simpson one-third"
            
    elif n % 3 == 0:                            # Simpson's three-eight rule
        soln = 3*(h/8)*(f_0 + 3*sum_mult3 + 2*sum_not_mult3 + f_n)
        rule = "Simpson three-eight"
            
    else: print("Error")  
    print(rule, "was used")
    return soln
    