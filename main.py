from newton_rule import numerical_rules 
import numpy as np 
import pandas as pd 
import json

order = False
while order == False:
    n = int(input("What is the reaction order (0/1/2/3/4): "))

    if n >= 0:
        from models import ReactorModels as RM
        order = True
    else:
        print("Only order from 1 to 4 is accepted")
        order = False


try: 
    t_unit = input("Enter the time unit (s/min/hr): ")
    data_source = input("enter or load csv data? (e/l): ")
    if data_source =="l":
        filename = f"input/c_t_data{n}.csv"
        df = pd.read_csv(filename, header = None)
        time = df.iloc[0, 1:].astype(float).tolist()
        conc = df.iloc[1, 1:].astype(float).tolist()
        k = df.iloc[2, 1].astype(float)

        print("time:", time)
        print("conc:", conc)
        print("k_rate:", k)

except FileNotFoundError:
    print("file not found in the directory")

    # t_unit = input("Enter the time unit (s/min/hr): ")
    
    time = []
    while True:
        tm = float(input("Enter the time data (101010 to next)"))
        if tm == 101010:
            break
        time.append(tm)
    h = time[1] - time[0]
    print(time)

    # Concentration data entry
    conc = []
    for t in time:
        conc.append(float(input(f"Enter the concentration at {t} {t_unit}: ")))
    k = float(input("enter the value of the rate constant to continue, k: "))
    content = [
        f"time", time,
        f"conc", conc,
        f"k", k
    ]
    with open(f"input/c_t_data{n}.csv", 'w') as file_save:
        json.dump(content, file_save, indent=4)
else: pass

h = time[1]-time[0]
# Age time calculation 
param = conc
integral = numerical_rules(param, h)
age_time = [c/integral for c in conc]
print("age_time:", age_time)
print("concentration integral:", integral)

# time*age_time and time^2*age_time calculation
t_Et = [t*a for t, a in zip(time, age_time)]
t2_Et = [t*t*a for t, a in zip(time, age_time)]

# mean residence time
param = t_Et
mrt = numerical_rules(param, h)
mrt = round(mrt, 3)

param = t2_Et
variance = numerical_rules(param, h) - mrt**2
variance = round(variance, 3)
if k == None:
    k = float(input("enter the value of the rate constant to continue, k: "))
else: pass

# Tank in series model calculation
conversion_tis = RM.tank_in_series(mrt, variance, k, conc, n)
conversion_tis = round(conversion_tis, 3)
print("Tank in series model converstion: ", conversion_tis)

# Dispersion model calculation
conversion_dispersion = RM.dispersion(mrt, variance, k, conc, n)
conversion_dispersion = round(conversion_dispersion, 3)
print("Dispersion model converstion: ", conversion_dispersion)

# Ideal reactors model
conversion_pfr = RM.ideal_pfr_model(mrt, k, conc, n)
conversion_pfr = round(conversion_pfr, 3)

conversion_cstr = RM.ideal_cstr_model(mrt, k, conc, n)
conversion_cstr = round(conversion_cstr, 3)

print("time:", time)
print("conc:", conc)
print("age_time:", age_time)
print("t_Et:", t_Et)
print("t2_Et:", t2_Et)
print(f"the mean residence time is {mrt} and the variance is {variance}")

# Save the generated data table
output = {
    f"time ({t_unit})": time,
    f"conc (mol/{t_unit})": conc,
    f"E(t) ({t_unit})": age_time,
    f"t.E(t)": t_Et,
    f"t^2.E(t)": t2_Et
}
# Save the models conversion predictions
predictions = {
    f"conc integral": integral,
    f"mean residence time": mrt,
    f"Variance": variance,
    f"Tank in series model converstion": conversion_tis,
    f"Dispersion model converstion": conversion_dispersion,
    f"Ideal PFR model conversion": conversion_pfr,
    f"Ideal CSTR model convertion": conversion_cstr
}

with open(f"output/pred n={n}, k={k}.json", 'w') as f:
    json.dump(predictions, f, indent=4)

with open(f"output/output_table n={n}.json", 'w') as f2:
    json.dump(output, f2, indent=4)
