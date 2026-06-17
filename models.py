class ReactorModels:
    """ Compile ideal and parameter reactor models for conversion prediction for all orders reaction.
    i.e. tank in series model (tubular reactor modelled as a tank) with parameter number of tanks, 
    dispersion model (cstr with dispersion) with parameter dispersion coeficient, ideal plug flow reactor, 
    ideal continuous stirred tank reactor """

    def __init__ (self, mrt, variance, k):
        """ initializes an instance of a reactor """
        self.mrt = mrt
        self.variance = variance
        self.k = k

    # Tank in series model
    def tank_in_series (mrt, variance, k, conc, n):
        from scipy.optimize import fsolve
        import math
        nt = (mrt**2)/variance
        if n == 0:
            conversion = min(1.0, (k * mrt) / conc[0])
        elif n == 1:
            conversion = 1 - (1/((1+(mrt*k)/nt)**nt))
        else: 
            mrt_i = mrt/nt  # mrt per individual tank
            C_current = conc[0]
            # Iterate sequentially through each tank stage
            nt = math.ceil(nt)
            for _ in range(nt):
                def stage_residual(C_next):
                    return C_current - C_next - k * mrt_i * (C_next**n)
                C_current = fsolve(stage_residual, x0=C_current)[0]
            try:
                conversion = 1.0 - (C_current / conc[0])
            except ZeroDivisionError:
                print("zero concentration, can not divid by zero")
            else:
                pass
        conversion = round(conversion, 3)
        print("number of tanks: ", nt)
        print("conversion: ", conversion)
        return conversion

    # Dispersion model

    def dispersion(mrt, variance, k, conc, n):
        import numpy as np
        from scipy.optimize import fsolve
        from scipy.integrate import solve_bvp
        
        v_mrt = variance/mrt**2
        def variance_eq(pe):
            """Closed-closed varaince/tau^2 as function of pe"""
            if pe <= 0:
                return 1e6
            term = (2/pe) - (2/pe**2) * (1 - np.exp(-pe))
            return term - v_mrt
        
        # Initial guess. try a few because pe ca have 2 solution for high variance
        guesses = [0.5, 2, 5, 20, 100]
        solutions = set()

        for g in guesses:
            pe_sol = fsolve(variance_eq, g)[0]    
            # keep only positive, real solutions
            if pe_sol > 0 and abs(variance_eq(pe_sol)) < 1e-6:
                solutions.add(round(pe_sol,4))
        print(f"variance/mrt^2 = {v_mrt}")
        print(f"pe solutions: {sorted(solutions)}")

        # verify
        for pe in sorted(solutions):
            check = (2/pe) - (2/pe**2) * (1 - np.exp(-pe))
            print(f"pe = {pe:.4f} -> variance/mrt^2 calc = {check:.4f}")
        
        da = k * (conc[0]**(n - 1)) * mrt
        q = (1 + (4*da/pe))**(0.5)

        if n == 0:
            conversion = min(1.0, (k * mrt) / conc[0])
        elif n == 1:
            numerator = (4*q*np.exp(pe/2))
            denominator = ((1+q)**2 * np.exp(pe*q/2) - (1-q)**2 * np.exp(-pe*q/2))
            conversion = 1 - (numerator / denominator)
        else: 
            # Solve 2nd order non-linear ODE Boundary Value Problem (BVP)
            # Let y[0] = Gamma (dimensionless conc C/C_0), y[1] = dGamma/dz
            
            def ode_system(z, y):
                # dGamma = y[1]
                # d2Gamma/dz2 = Pe * (dGamma/dz + da * Gamma^n)
                return np.vstack((y[1], pe * (y[1] + da * (y[0]**n))))
            def boundary_conditions(ya, yb):
                # Inlet Danckwerts: ya[0] - (1/Pe)*ya[1] - 1 = 0
                # Outlet Danckwerts: yb[1] = 0
                return np.array([ya[0] - (1.0 / pe) * ya[1] - 1.0, yb[1]])
            
            # Grid initialization
            z_mesh = np.linspace(0, 1, 100)
            # Initial guess assuming linear drop across the length
            X_pfr = 1.0 - (1.0 + (n -1) * k * (conc[0]**(n - 1)) * mrt)**(1.0 / (1.0 - n))
            y_guess = np.vstack((np.linspace(1, 1 - X_pfr, 100), np.zeros(100)))

            sol = solve_bvp(ode_system, boundary_conditions, z_mesh, y_guess)
            Gamma_exit  = sol.y[0, -1]
            conversion = 1.0 - Gamma_exit

        conversion = round(conversion, 3)
        output = {
            "pe": pe.astype(float),
            "q": q.astype(float), 
            "da": da
        }
        print(output)
        return conversion

    def ideal_pfr_model(mrt, k, conc, n):
        """ predictive conversion for ideal plug flow reactor model"""
        import numpy as np
        if n == 0:
            conversion = min(1.0, (k * mrt) / conc[0])
        elif n == 1:
            conversion = 1 - np.exp(-k*mrt)
        elif n == 2:
            numerator =  k * conc[0] * mrt
            denominator = (1 + k * conc[0]*mrt)
            conversion = numerator / denominator
        elif n >= 3:
            conversion = 1.0 - (1.0 + (n -1) * k * (conc[0]**(n - 1)) * mrt)**(1.0 / (1.0 - n))
        else:
            print("Out of bound")

        conversion = round(conversion, 3)
        print(f"pfr conversion: {conversion:.3f}")
        return conversion

    def ideal_cstr_model(mrt, k, conc, n):
        """ predictive conversion for ideal continuous stirred tank reactor model """
        import numpy as np 
        from scipy.optimize import fsolve

        if n == 0: 
            conversion = min(1.0, (k * mrt) / conc[0])
        elif n == 1: 
            conversion = (k*mrt)/(1 + k*mrt)
        elif n == 2:
            numerator = (-1 + np.sqrt(1 + 4 * k * conc[0] * mrt))
            denominator = 2 * k * conc[0]
            conversion = numerator / denominator
        elif n >= 3: 
            # Solves for Algebraic equation
            def cstr_residual(X):
                return mrt - X / (k *(conc[0]**(n-1)) * (1 - X)**n)
            conversion = fsolve(cstr_residual, x0=0.5)[0]
        conversion = round(conversion, 3)
        print(f"cstr conversion: {conversion:.3f}")
        return conversion

