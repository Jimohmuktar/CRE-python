# Reaction Models: Chemical Reactor Conversion Prediction

## Overview

**Reaction Models** is a comprehensive Python package for predicting chemical conversion rates in various ideal and non-ideal reactor configurations. This project implements four distinct reactor models that simulate real-world chemical reactors under different operating conditions and reaction orders.

The package enables engineers and chemists to:
- Model ideal plug flow reactors (PFR) and continuous stirred tank reactors (CSTR)
- Account for non-ideal mixing using tank-in-series and dispersion models
- Handle reactions of various orders (zero-order through higher-order reactions)
- Calculate residence time distributions and reactor performance metrics
- Validate reactor performance through numerical analysis and integration techniques

---

## Key Features

- **Four Reactor Models:**
  - Ideal Plug Flow Reactor (PFR)
  - Ideal Continuous Stirred Tank Reactor (CSTR)
  - Tank-in-Series Model (accounts for non-ideal mixing)
  - Dispersion Model (Danckwerts boundary conditions)

- **Flexible Reaction Order Support:**
  - Zero-order reactions
  - First-order reactions
  - Second-order reactions
  - Higher-order reactions (n ≥ 3)

- **Advanced Numerical Methods:**
  - Trapezoidal integration rule
  - Simpson's one-third rule
  - Simpson's three-eighths rule
  - Adaptive rule selection based on data point count

- **Comprehensive Analysis:**
  - Mean residence time (MRT) calculation
  - Variance and residence time distribution analysis
  - Age time distribution (E(t)) computation
  - Boundary value problem solver for complex dispersion models

- **Data Handling:**
  - CSV input support for experimental data
  - JSON output for results and predictions
  - Flexible data entry via interactive console or file loading

---

## Theoretical Background

### Chemical Kinetics and Reaction Rates

The fundamental basis of this package is chemical kinetics, which describes how reaction rates depend on concentration, temperature, and reaction order.

#### Rate Law Expression

For a general reaction:
```
aA + bB → products
```

The rate of reaction is expressed as:

```
-dC_A/dt = k·C_A^n
```

Where:
- `C_A` = concentration of reactant A (mol/L)
- `t` = time (s, min, or hr)
- `k` = rate constant (units depend on reaction order)
- `n` = reaction order
- The negative sign indicates reactant depletion

#### Reaction Orders and Rate Laws

**Zero-Order Reactions (n = 0):**
```
-dC/dt = k
```
The reaction rate is independent of concentration. Integration gives:
```
C(t) = C_0 - k·t
```
Conversion at residence time τ:
```
X = min(1.0, k·τ / C_0)
```

**First-Order Reactions (n = 1):**
```
-dC/dt = k·C
```
Integration yields exponential decay:
```
C(t) = C_0·e^(-k·t)
```
Conversion at residence time τ:
```
X = 1 - e^(-k·τ)
```

**Second-Order Reactions (n = 2):**
```
-dC/dt = k·C^2
```
Integration gives:
```
1/C(t) = 1/C_0 + k·t
```
Conversion at residence time τ:
```
X = (k·C_0·τ) / (1 + k·C_0·τ)
```
For CSTR:
```
X = (-1 + √(1 + 4·k·C_0·τ)) / (2·k·C_0)
```

**Higher-Order Reactions (n ≥ 3):**
```
-dC/dt = k·C^n
```
For n > 1 (general case):
```
X = 1 - (1 + (n-1)·k·C_0^(n-1)·τ)^(1/(1-n))
```

### Reactor Models

#### 1. Ideal Plug Flow Reactor (PFR)

**Characteristics:**
- No mixing in the axial direction (flow direction)
- Complete radial mixing
- Residence time is constant for all fluid elements
- Plug flow assumption: δC/δz = 0 (perpendicular to flow)

**Residence Time in PFR:**
```
τ = V / Q = Length / Velocity
```

**Conversion Equations for PFR:**

For first-order reactions:
```
X_PFR = 1 - e^(-k·τ)
```

For second-order reactions:
```
X_PFR = (k·C_0·τ) / (1 + k·C_0·τ)
```

For higher-order reactions:
```
X_PFR = 1 - (1 + (n-1)·k·C_0^(n-1)·τ)^(1/(1-n))
```

**Advantages:**
- Maximum conversion for a given residence time
- Efficient for fast reactions
- Scalable design

---

#### 2. Ideal Continuous Stirred Tank Reactor (CSTR)

**Characteristics:**
- Complete mixing throughout the reactor
- Uniform concentration and temperature in the entire volume
- Exit concentration equals internal concentration
- Backmixing is complete

**Residence Time in CSTR:**
```
τ = V / Q = (1/q)
```
where q is the volumetric flow rate normalized by volume.

**Conversion Equations for CSTR:**

For zero-order reactions:
```
X_CSTR = min(1.0, k·τ / C_0)
```

For first-order reactions:
```
X_CSTR = (k·τ) / (1 + k·τ)
```

For second-order reactions (implicit equation solved numerically):
```
τ = X / (k·C_0^(n-1)·(1-X)^n)
```

**Advantages:**
- Simple to operate and control
- Excellent for exothermic reactions (easy heat removal)
- More stable against feed composition variations

**Disadvantages:**
- Lower conversion for same residence time compared to PFR
- Requires more intense mixing energy

---

#### 3. Tank-in-Series Model

**Characteristics:**
- Models a real reactor as N ideal CSTR units in series
- Bridges the gap between ideal CSTR (N=1) and ideal PFR (N→∞)
- Number of tanks calculated from variance of residence time distribution

**Number of Tanks Calculation:**

From residence time distribution statistics:
```
N = (τ_mean)² / σ²
```

Where:
- `τ_mean` = mean residence time
- `σ²` = variance

**Each Tank Equations:**

For the i-th tank with residence time `τ_i = τ_mean/N`:

For first-order reactions:
```
X_i = 1 - 1 / (1 + k·τ_i)^N
```

For higher-order reactions (solved iteratively):
```
C_i = C_(i-1) - k·τ_i·(C_i)^n
```

**Solution Method:**
- Sequential solution through each tank stage
- Nonlinear equation solver (fsolve) for implicit equations
- Iterative approach ensures convergence

**Physical Interpretation:**
- N=1: Single CSTR (complete backmixing)
- 2<N<10: Moderate mixing (real reactors)
- N>20: Approaches PFR behavior (minimal backmixing)

---

#### 4. Dispersion Model (Axial Dispersion)

**Characteristics:**
- Accounts for both convection and longitudinal dispersion
- Modeled as PFR with axial dispersion coefficient
- Uses Danckwerts boundary conditions at inlet/outlet
- More realistic for real reactors

**Dimensionless Parameters:**

**Peclet Number (Pe):**
```
Pe = U·L / D = τ·k_dispersion / D
```
Where:
- `U` = superficial velocity
- `L` = reactor length
- `D` = axial dispersion coefficient

**Damköhler Number (Da):**
```
Da = k·C_0^(n-1)·τ
```
Ratio of reaction rate to convection rate

**Variance Relationship:**
```
σ²/τ² = 2/Pe - 2/(Pe²)·(1 - e^(-Pe))
```

**Governing Differential Equations:**

Dimensionless dispersion model:
```
Pe·d²Γ/dz² - Pe·dΓ/dz - Da·Γ^n = 0
```

Where:
- `Γ = C/C_0` (dimensionless concentration)
- `z = Z/L` (dimensionless axial position)

**Boundary Conditions (Danckwerts):**
```
Inlet (z=0):   Γ_0 - (1/Pe)·dΓ/dz|_0 = 1
Outlet (z=1):  dΓ/dz|_1 = 0
```

**First-Order Reaction Solution:**
```
q = √(1 + 4·Da/Pe)
X = 1 - (4·q·e^(Pe/2)) / ((1+q)²·e^(Pe·q/2) - (1-q)²·e^(-Pe·q/2))
```

**Higher-Order Reactions:**
- Solved as boundary value problem (BVP)
- Using scipy.integrate.solve_bvp
- Requires numerical integration of 2nd-order ODE system

---

### Residence Time Distribution (RTD)

The age distribution E(t) is fundamental to all reactor models.

**Age Distribution Function:**
```
E(t) = C(t) / ∫₀^∞ C(t)dt
```

**Mean Residence Time:**
```
τ_mean = ∫₀^∞ t·E(t)dt
```

**Variance:**
```
σ² = ∫₀^∞ (t - τ_mean)²·E(t)dt = ∫₀^∞ t²·E(t)dt - τ_mean²
```

**RTD Moments:**
- 0th moment: ∫E(t)dt = 1 (normalization)
- 1st moment: ∫t·E(t)dt = τ_mean
- 2nd moment: ∫t²·E(t)dt = variance + τ_mean²

---

## Project Structure

```
Reaction Models/
├── README.md                 # This file
├── models.py                 # Core reactor model implementations
├── newton_rule.py            # Numerical integration methods
├── main.py                   # Interactive console application
├── quick_test.py             # Unit tests for all models
├── input/                    # Input data directory (CSV files)
│   ├── c_t_data0.csv        # Sample zero-order reaction data
│   ├── c_t_data1.csv        # Sample first-order reaction data
│   ├── c_t_data2.csv        # Sample second-order reaction data
│   └── c_t_data3.csv        # Sample higher-order reaction data
├── output/                   # Output directory (JSON results)
│   ├── pred n={n}, k={k}.json      # Conversion predictions
│   └── output_table n={n}.json     # Numerical analysis tables
└── __pycache__/              # Python bytecode cache
```

---

## Installation

### Prerequisites

- Python 3.7 or higher
- pip (Python package manager)

### Required Dependencies

```
numpy>=1.19.0           # Numerical computing
scipy>=1.5.0            # Scientific computing, optimization, integration
pandas>=1.1.0           # Data manipulation and analysis
```

### Setup Instructions

**Step 1: Clone or Download Repository**

```bash
git clone <repository-url>
cd "Reaction Models"
```

**Step 2: Create Virtual Environment (Recommended)**

```bash
# Windows
python -m venv venv
venv\Scripts\activate

# macOS/Linux
python3 -m venv venv
source venv/bin/activate
```

**Step 3: Install Dependencies**

```bash
pip install -r requirements.txt
```

Or manually install:

```bash
pip install numpy scipy pandas
```

**Step 4: Verify Installation**

Run the test suite:

```bash
python quick_test.py
```

Expected output:
```
....
----------------------------------------------------------------------
Ran 5 tests in X.XXXs

OK
```

---

## Usage

### Interactive Console Application

The main entry point is `main.py`, which provides an interactive interface for reactor analysis.

#### Running the Application

```bash
python main.py
```

#### Input Parameters

When running `main.py`, you'll be prompted for:

1. **Reaction Order (n):** 0, 1, 2, 3, or 4
2. **Time Unit:** seconds (s), minutes (min), or hours (hr)
3. **Data Source:** 
   - `e` = Enter data manually
   - `l` = Load from CSV file
4. **Concentration Data:** Time-concentration pairs
5. **Rate Constant (k):** Temperature-dependent kinetic parameter

#### Example Session

```
What is the reaction order (0/1/2/3/4): 1
Enter the time unit (s/min/hr): min
enter or load csv data? (e/l): l

Data loaded successfully:
time: [0, 5, 10, 15, 20, 25]
conc: [1.0, 0.82, 0.67, 0.55, 0.45, 0.37]
k_rate: 0.15

Processing...
mean residence time: 12.5 min
variance: 8.33
```

#### Manual Data Entry

```
What is the reaction order (0/1/2/3/4): 2
Enter the time unit (s/min/hr): s
enter or load csv data? (e/l): e

Enter the time data (101010 to stop): 0
Enter the time data (101010 to stop): 10
Enter the time data (101010 to stop): 20
Enter the time data (101010 to stop): 30
Enter the time data (101010 to stop): 101010

Enter the concentration at 0 s: 2.5
Enter the concentration at 10 s: 2.1
Enter the concentration at 20 s: 1.8
Enter the concentration at 30 s: 1.5

enter the value of the rate constant to continue, k: 0.05
```

---

## Module Documentation

### models.py

Core module containing all four reactor model implementations.

#### Class: ReactorModels

```python
class ReactorModels:
    """Compile ideal and parameter reactor models for conversion prediction"""
    
    def __init__(self, mrt, variance, k):
        """
        Initialize reactor instance
        
        Parameters:
        -----------
        mrt : float
            Mean residence time (in time units)
        variance : float
            Variance of residence time distribution
        k : float
            Rate constant (units depend on reaction order)
        """
```

#### Method: ideal_pfr_model()

```python
@staticmethod
def ideal_pfr_model(mrt, k, conc, n):
    """
    Predict conversion for ideal plug flow reactor
    
    Parameters:
    -----------
    mrt : float
        Mean residence time
    k : float
        Rate constant
    conc : list
        Concentration profile [C_0, C_1, ..., C_n]
    n : int
        Reaction order (0, 1, 2, or ≥3)
    
    Returns:
    --------
    conversion : float
        Fractional conversion (0.0 to 1.0)
    
    Equation:
    ---------
    For n=1: X = 1 - exp(-k*mrt)
    For n=2: X = (k*C_0*mrt)/(1 + k*C_0*mrt)
    For n≥3: X = 1 - (1 + (n-1)*k*C_0^(n-1)*mrt)^(1/(1-n))
    """
```

#### Method: ideal_cstr_model()

```python
@staticmethod
def ideal_cstr_model(mrt, k, conc, n):
    """
    Predict conversion for ideal continuous stirred tank reactor
    
    Parameters:
    -----------
    mrt : float
        Mean residence time
    k : float
        Rate constant
    conc : list
        Concentration profile [C_0, C_1, ..., C_n]
    n : int
        Reaction order (0, 1, 2, or ≥3)
    
    Returns:
    --------
    conversion : float
        Fractional conversion (0.0 to 1.0)
    
    Equation:
    ---------
    For n=1: X = (k*mrt)/(1 + k*mrt)
    For n=2: X = (-1 + √(1 + 4*k*C_0*mrt))/(2*k*C_0)
    For n≥3: Solved implicitly via fsolve
    """
```

#### Method: tank_in_series()

```python
@staticmethod
def tank_in_series(mrt, variance, k, conc, n):
    """
    Model real reactor as N ideal CSTRs in series
    
    Parameters:
    -----------
    mrt : float
        Mean residence time
    variance : float
        Variance of residence time distribution
    k : float
        Rate constant
    conc : list
        Concentration profile
    n : int
        Reaction order
    
    Returns:
    --------
    conversion : float
        Fractional conversion from tank-in-series model
    
    Algorithm:
    ----------
    1. Calculate number of tanks: N = (mrt)² / variance
    2. Residence time per tank: τ_i = mrt / N
    3. Solve sequentially through each stage
    4. For implicit equations, use scipy.optimize.fsolve
    """
```

#### Method: dispersion()

```python
@staticmethod
def dispersion(mrt, variance, k, conc, n):
    """
    Model reactor with axial dispersion (Danckwerts BCs)
    
    Parameters:
    -----------
    mrt : float
        Mean residence time
    variance : float
        Variance of residence time distribution
    k : float
        Rate constant
    conc : list
        Concentration profile
    n : int
        Reaction order
    
    Returns:
    --------
    conversion : float
        Fractional conversion from dispersion model
    
    Algorithm:
    ----------
    1. Calculate Peclet number from variance equation
    2. Solve for Pe (may have multiple solutions)
    3. Calculate Damköhler number
    4. For n≥3: Solve 2nd-order BVP using scipy.integrate.solve_bvp
    5. Extract exit conversion from solution
    
    Boundary Conditions:
    --------------------
    Inlet:  Γ_0 - (1/Pe)*dΓ/dz|_0 = 1
    Outlet: dΓ/dz|_1 = 0
    """
```

---

### newton_rule.py

Numerical integration methods for computing residence time distribution moments.

#### Function: numerical_rules()

```python
def numerical_rules(param, h):
    """
    Apply numerical integration rule based on data point count
    
    Parameters:
    -----------
    param : list
        Array of function values at equally-spaced points
    h : float
        Step size (spacing between points)
    
    Returns:
    --------
    integral : float
        Approximate integral value
    
    Method Selection:
    -----------------
    - If n is odd and n % 3 != 0:     Trapezoidal rule
    - If n is even and n % 3 != 0:    Simpson's one-third rule
    - If n % 3 == 0:                  Simpson's three-eighths rule
    """
```

#### Integration Rules

**Trapezoidal Rule:**
```
∫f dx ≈ (h/2) * [f_0 + 2∑f_i + f_n]
```
Error: O(h²)
Best for: Linear or nearly-linear functions

**Simpson's One-Third Rule:**
```
∫f dx ≈ (h/3) * [f_0 + 4∑f_even + 2∑f_odd + f_n]
```
Requires: Even number of intervals
Error: O(h⁴)
Best for: Polynomial functions (cubic and below)

**Simpson's Three-Eighths Rule:**
```
∫f dx ≈ (3h/8) * [f_0 + 3∑f_(3k+1) + 3∑f_(3k+2) + 2∑f_(3k+3) + f_n]
```
Requires: Number of intervals divisible by 3
Error: O(h⁴)
Best for: Higher precision with many points

---

## Examples

### Example 1: First-Order Reaction in PFR

**Problem:** A first-order irreversible reaction A → products occurs in a plug flow reactor. Given:
- Initial concentration: C₀ = 1.0 mol/L
- Rate constant: k = 0.15 min⁻¹
- Residence time: τ = 20 min

**Solution:**

```python
from models import ReactorModels as RM

# Parameters
mrt = 20.0
k = 0.15
conc = [1.0]
n = 1  # first-order

# Calculate conversion
conversion_pfr = RM.ideal_pfr_model(mrt, k, conc, n)
print(f"PFR Conversion: {conversion_pfr:.3f}")

# Expected: X = 1 - exp(-0.15*20) = 1 - exp(-3) = 0.950
```

**Output:**
```
pfr conversion: 0.950
```

---

### Example 2: Second-Order Reaction Comparison

**Problem:** Compare all four reactor models for a second-order reaction:
- Initial concentration: C₀ = 2.0 mol/L
- Rate constant: k = 0.1 L/(mol·min)
- Residence time: τ = 15 min
- Variance: σ² = 25 min²

**Solution:**

```python
from models import ReactorModels as RM

# Parameters
mrt = 15.0
variance = 25.0
k = 0.1
conc = [2.0]
n = 2  # second-order

# Calculate conversions with all models
conversion_pfr = RM.ideal_pfr_model(mrt, k, conc, n)
conversion_cstr = RM.ideal_cstr_model(mrt, k, conc, n)
conversion_tis = RM.tank_in_series(mrt, variance, k, conc, n)
conversion_disp = RM.dispersion(mrt, variance, k, conc, n)

print(f"PFR:  {conversion_pfr:.3f}")
print(f"CSTR: {conversion_cstr:.3f}")
print(f"Tank-in-Series: {conversion_tis:.3f}")
print(f"Dispersion: {conversion_disp:.3f}")
```

**Expected Results:**
```
PFR:  0.750          (best conversion)
CSTR: 0.545          (lowest conversion)
Tank-in-Series: 0.680 (between CSTR and PFR)
Dispersion: 0.670    (between CSTR and PFR)
```

**Interpretation:**
- PFR gives highest conversion due to no backmixing
- CSTR has lowest conversion due to complete mixing
- Tank-in-Series and Dispersion models fall between ideal reactors
- Higher variance increases backmixing effects (lower conversion)

---

### Example 3: Real Reactor Data Analysis

**Input CSV File (c_t_data1.csv):**
```
time,0,5,10,15,20,25
conc,1.0,0.82,0.67,0.55,0.45,0.37
k,0.15
```

**Running the Analysis:**

```bash
python main.py
# Select: order=1, time_unit='min', data_source='l'
```

**Output Files:**

`output/pred n=1, k=0.15.json`:
```json
{
    "conc integral": 4.125,
    "mean residence time": 12.5,
    "Variance": 8.33,
    "Tank in series model converstion": 0.953,
    "Dispersion model converstion": 0.96,
    "Ideal PFR model conversion": 0.988,
    "Ideal CSTR model convertion": 0.815
}
```

`output/output_table n=1.json`:
```json
{
    "time (min)": [0, 5, 10, 15, 20, 25],
    "conc (mol/min)": [1.0, 0.82, 0.67, 0.55, 0.45, 0.37],
    "E(t) (min)": [0.242, 0.199, 0.162, 0.133, 0.109, 0.090],
    "t.E(t)": [0.0, 0.994, 1.620, 1.995, 2.180, 2.250],
    "t^2.E(t)": [0.0, 4.970, 16.200, 29.925, 43.600, 56.250]
}
```

---

### Example 4: Using Numerical Integration

**Problem:** Calculate mean residence time from age distribution

```python
from newton_rule import numerical_rules

# Concentration measurements at equal time intervals (h=5)
conc = [1.0, 0.82, 0.67, 0.55, 0.45, 0.37]
h = 5.0  # time step

# Calculate integral (normalization factor)
integral = numerical_rules(conc, h)

# Age distribution
age_time = [c/integral for c in conc]

# Time-weighted age distribution
t_Et = [t*a for t, a in zip([0,5,10,15,20,25], age_time)]

# Calculate mean residence time
mrt = numerical_rules(t_Et, h)

print(f"Integral: {integral:.3f}")
print(f"Age Distribution: {[f'{a:.3f}' for a in age_time]}")
print(f"Mean Residence Time: {mrt:.3f} min")
```

**Output:**
```
Integral: 4.125
Age Distribution: ['0.242', '0.199', '0.162', '0.133', '0.109', '0.090']
Mean Residence Time: 12.500 min
Simpson one-third was used
```

---

## Testing

Comprehensive unit tests are provided in `quick_test.py`.

### Running Tests

```bash
python -m unittest quick_test.py -v
```

### Test Coverage

The test suite includes:

1. **Numerical Integration Test:**
   - Validates Trapezoidal/Simpson integration
   - Sample data: [0,3,5,5,4,2,1,0]
   - Expected result: 100.0

2. **Tank-in-Series Model Test:**
   - Tests first-order reaction (n=1)
   - Sample data: conc=[0,3,5,5,4,2,1,0], mrt=14.667, variance=48.212, k=0.3
   - Expected conversion: 0.953

3. **Dispersion Model Test:**
   - Tests first-order reaction
   - Same sample parameters as tank-in-series
   - Expected conversion: 0.96

4. **Ideal CSTR Model Test:**
   - Tests first-order reaction (n=1)
   - Expected conversion: 0.815

5. **Ideal PFR Model Test:**
   - Tests first-order reaction (n=1)
   - Expected conversion: 0.988

### Test Output Format

```
test_dispersion_model (quick_test.OverallCaseTest) ... ok
test_idealcstr_model (quick_test.OverallCaseTest) ... ok
test_idealpfr_model (quick_test.OverallCaseTest) ... ok
test_numerical_analysis (quick_test.OverallCaseTest) ... ok
test_tank_model (quick_test.OverallCaseTest) ... ok

----------------------------------------------------------------------
Ran 5 tests in 0.523s

OK
```

### Running Individual Tests

```bash
# Test only the PFR model
python -m unittest quick_test.OverallCaseTest.test_idealpfr_model -v

# Test only the CSTR model
python -m unittest quick_test.OverallCaseTest.test_idealcstr_model -v
```

---

## Output

The application generates two types of output files in the `output/` directory:

### 1. Predictions File: `pred n={n}, k={k}.json`

Contains all calculated reactor performance metrics:

```json
{
    "conc integral": 4.125,
    "mean residence time": 12.5,
    "Variance": 8.33,
    "Tank in series model converstion": 0.953,
    "Dispersion model converstion": 0.96,
    "Ideal PFR model conversion": 0.988,
    "Ideal CSTR model convertion": 0.815
}
```

**Fields:**
- `conc integral`: ∫C(t)dt - normalization of age distribution
- `mean residence time`: τ = ∫t·E(t)dt (minutes, seconds, or hours)
- `Variance`: σ² = ∫t²·E(t)dt - (τ)²
- `Tank in series model converstion`: X_TIS
- `Dispersion model converstion`: X_Dispersion
- `Ideal PFR model conversion`: X_PFR (maximum conversion)
- `Ideal CSTR model convertion`: X_CSTR (minimum conversion)

### 2. Data Table File: `output_table n={n}.json`

Contains numerical analysis and intermediate calculations:

```json
{
    "time (min)": [0, 5, 10, 15, 20, 25],
    "conc (mol/min)": [1.0, 0.82, 0.67, 0.55, 0.45, 0.37],
    "E(t) (min)": [0.242, 0.199, 0.162, 0.133, 0.109, 0.090],
    "t.E(t)": [0.0, 0.994, 1.620, 1.995, 2.180, 2.250],
    "t^2.E(t)": [0.0, 4.970, 16.200, 29.925, 43.600, 56.250]
}
```

**Columns:**
- `time`: Input time points
- `conc`: Input concentration measurements
- `E(t)`: Age distribution (age time) = C(t) / ∫C(t)dt
- `t.E(t)`: Time-weighted age distribution
- `t^2.E(t)`: Time²-weighted age distribution (for variance calculation)

---

## Requirements

### System Requirements

- **Python Version:** 3.7 or higher
- **Operating System:** Windows, macOS, or Linux
- **Memory:** Minimum 2GB RAM (recommended for large datasets)
- **Disk Space:** ~50MB for installation + dependencies

### Python Dependencies

```
numpy>=1.19.0           # Required for numerical operations
scipy>=1.5.0            # Required for optimization and integration
pandas>=1.1.0           # Required for data handling and CSV I/O
```

### Installation of Dependencies

Create a `requirements.txt` file:

```
numpy>=1.19.0
scipy>=1.5.0
pandas>=1.1.0
```

Install via pip:

```bash
pip install -r requirements.txt
```

---

## Troubleshooting

### Issue 1: ModuleNotFoundError: No module named 'scipy'

**Error Message:**
```
ModuleNotFoundError: No module named 'scipy'
```

**Solution:**
```bash
pip install scipy
# Or install all requirements
pip install -r requirements.txt
```

---

### Issue 2: File Not Found Error when Loading CSV

**Error Message:**
```
FileNotFoundError: file not found in the directory
```

**Solutions:**
1. Ensure CSV file exists in `input/` directory
2. Check filename format: `c_t_data{n}.csv` (where n is reaction order)
3. Verify path is correct: relative to current working directory

**CSV File Format:**
```
First row: time values (0, 5, 10, 15, 20, ...)
Second row: concentration values (1.0, 0.82, 0.67, ...)
Third row: rate constant k value
```

---

### Issue 3: ZeroDivisionError in Tank-in-Series Model

**Error Message:**
```
ZeroDivisionError: zero concentration, can not divid by zero
```

**Causes:**
- Reaction is so fast that all reactant is consumed
- Initial concentration is zero or near-zero
- Residence time is extremely long

**Solutions:**
1. Check initial concentration (conc[0] > 0)
2. Verify rate constant is reasonable
3. For very fast reactions, conversion is naturally 1.0 (100%)

---

### Issue 4: Convergence Problems in Dispersion Model

**Error Message:**
```
RuntimeWarning: Optimization failed. Could not find valid Peclet solution.
```

**Causes:**
- Very high variance/low Peclet number
- Unusual input parameters
- Multiple Pe solutions exist

**Solutions:**
1. Verify variance calculation is correct
2. Check that variance < mrt² (mathematical requirement)
3. Review concentration measurements for data quality

---

### Issue 5: Slow Execution or Memory Issues

**Symptoms:**
- Program runs very slowly
- Memory usage increases excessively
- BVP solver doesn't converge

**Solutions:**
1. Reduce number of data points (higher time step h)
2. Simplify grid in solve_bvp (fewer mesh points)
3. Check for infinite loops in manual data entry
4. Use file loading instead of manual entry for large datasets

---

### Issue 6: Unexpected Conversion Values

**Symptoms:**
- Conversion > 1.0 or < 0.0 (invalid range)
- Conversion is 1.0 for low rate constants
- Different models give very different results

**Solutions:**
1. Verify reaction order is correct
2. Check that rate constant units match time units
3. Ensure initial concentration is realistic
4. Validate input data (monotonically decreasing C(t))

---

## Contributing

### Reporting Bugs

If you encounter a bug, please:

1. Record the exact error message
2. Document input parameters (n, k, C₀, τ, σ²)
3. Include relevant output files (JSON)
4. Describe steps to reproduce

**Example Bug Report:**
```
Model: Tank-in-Series
Reaction Order: 2
Initial Concentration: 0.5 mol/L
Rate Constant: -0.1  # Invalid!
Error: ValueError in kinetics calculation
```

### Suggesting Improvements

Potential areas for enhancement:

1. **Reactor Models:**
   - Batch reactor model
   - Tubular reactor with plug flow
   - Packed bed reactor model
   - Fluidized bed model

2. **Reaction Kinetics:**
   - Temperature-dependent rate constants (Arrhenius)
   - Multiple parallel reactions
   - Series-parallel reaction networks
   - Catalyst deactivation models

3. **Numerical Methods:**
   - Runge-Kutta integration (RK4, RK45)
   - Adaptive step-size algorithms
   - Finite element method (FEM)
   - Collocation methods

4. **User Interface:**
   - Graphical user interface (GUI)
   - Web interface with visualization
   - Real-time plotting of results
   - 3D reactor visualization

5. **Data Features:**
   - Experimental data fitting
   - Parameter optimization
   - Sensitivity analysis
   - Monte Carlo uncertainty quantification

### Development Setup

```bash
# Clone repository
git clone <repository-url>
cd "Reaction Models"

# Create virtual environment
python -m venv venv
source venv/bin/activate  # or venv\Scripts\activate on Windows

# Install development dependencies
pip install -r requirements.txt
pip install pytest pytest-cov  # For testing

# Run tests
python -m pytest quick_test.py -v

# Code formatting (optional)
pip install black
black *.py
```

### Code Style Guidelines

- Follow PEP 8 style guide
- Use descriptive variable names
- Add docstrings to all functions
- Include type hints where possible
- Write unit tests for new features

---

## License

This project is provided for educational and research purposes.

**Usage Rights:**
- Free to use, modify, and distribute
- Cite this project in academic publications
- Include this LICENSE section in derivatives

---

## References

### Foundational Literature

1. **Chemical Reactor Theory**
   - Levenspiel, O. (1999). "Chemical Reaction Engineering" (3rd ed.). Wiley.
   - Froment, G. F., Bischoff, K. B., & De Wilde, J. (2011). "Chemical Reactor Analysis and Design" (3rd ed.). Wiley.

2. **Residence Time Distribution**
   - Zwietering, T. N. (1959). "The degree of mixing in continuous flow systems." Chemical Engineering Science, 11(1), 1-15.
   - Naumann, E. B., & Buffham, B. A. (1983). "Mixing in Continuous Flow Systems." Wiley-Interscience.

3. **Numerical Methods**
   - Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007). "Numerical Recipes: The Art of Scientific Computing" (3rd ed.).
   - Burden, R. L., & Faires, J. D. (2010). "Numerical Analysis" (9th ed.). Cengage Learning.

4. **Boundary Value Problems**
   - Ascher, U. M., Mattheij, R. B., & Russell, R. D. (1995). "Numerical Solution of Boundary Value Problems for Ordinary Differential Equations." SIAM.

### Relevant Standards

- **IUPAC Nomenclature:** International Union of Pure and Applied Chemistry
- **SI Units:** Bureau International des Poids et Mesures (BIPM)
- **Chemical Kinetics:** NIST Chemistry WebBook

### Related Software

- **MATLAB:** `ode45`, `bvp4c` for differential equation solving
- **Python Packages:**
  - SciPy: `scipy.integrate.solve_bvp`
  - NumPy: Array operations and numerical computing
  - SymPy: Symbolic mathematics
  - Matplotlib: Visualization (future enhancement)

---

## Author

**Jimoh Muktar** 
- Repository Maintainer
- Initial Implementation
- Documentation

---

## Appendix: Mathematical Derivations

### A1: PFR Conversion for Higher-Order Reactions

**Starting from rate law:**
```
-dC/dt = k·C^n,  for n > 1 and n ≠ 1
```

**Separation of variables:**
```
dC/C^n = -k·dt
∫[C₀ to C] C^(-n) dC = -k∫[0 to τ] dt
```

**Integration:**
```
[C^(1-n) / (1-n)] |[C₀ to C] = -k·τ
C^(1-n)/(1-n) - C₀^(1-n)/(1-n) = -k·τ
C^(1-n) - C₀^(1-n) = -k·(1-n)·τ
C^(1-n) = C₀^(1-n) - k·(1-n)·τ
C^(1-n) = C₀^(1-n)·[1 - (1-n)·k·C₀^(-(n-1))·τ]

Raising both sides to the power 1/(1-n):
C = C₀·[1 + (n-1)·k·C₀^(n-1)·τ]^(1/(1-n))
```

**Conversion:**
```
X = (C₀ - C)/C₀ = 1 - [1 + (n-1)·k·C₀^(n-1)·τ]^(1/(1-n))
```

### A2: Variance Relationship in Dispersion Model

**For a closed-closed reactor with dispersion:**

The variance (normalized by residence time squared) relates to Peclet number:

```
σ²/τ² = 2/Pe - 2/(Pe²)·(1 - e^(-Pe))
```

**For limiting cases:**

**1. Pe → ∞ (plug flow):**
```
σ²/τ² → 0
```
All fluid elements have same residence time.

**2. Pe → 0 (complete mixing):**
```
σ²/τ² → 1
```
Complete backmixing (CSTR behavior).

**3. Intermediate Pe:**
```
0 < σ²/τ² < 1
```
Partial backmixing (real reactor).

### A3: Numerical Integration Error Analysis

**Trapezoidal Rule Error:**
```
E_trap = -h³·f''(ξ)/12,  ξ ∈ [a,b]
|E_trap| = O(h²) per interval
|E_total| = O(h²)
```

**Simpson's One-Third Rule Error:**
```
E_Simpson = -h⁵·f⁴(ξ)/90,  ξ ∈ [a,b]
|E_Simpson| = O(h⁴) per interval
|E_total| = O(h⁴)
```

**Simpson's Three-Eighths Rule Error:**
```
E_3/8 = -3h⁵·f⁴(ξ)/80,  ξ ∈ [a,b]
|E_3/8| = O(h⁴) per interval
|E_total| = O(h⁴)
```

---

## Glossary

- **BVP:** Boundary Value Problem - ODE with boundary conditions at multiple points
- **CSTR:** Continuous Stirred Tank Reactor - perfectly mixed reactor
- **Da:** Damköhler number - ratio of reaction to convection timescales
- **E(t):** Age distribution function - probability density of residence times
- **mrt:** Mean Residence Time - average time fluid spends in reactor
- **n:** Reaction order - exponent in rate law
- **ODE:** Ordinary Differential Equation
- **Pe:** Peclet number - ratio of convection to dispersion
- **PFR:** Plug Flow Reactor - no backmixing
- **RTD:** Residence Time Distribution - statistical distribution of residence times
- **X:** Conversion - fraction of reactant consumed
- **τ:** Residence time - time spent in reactor
- **σ²:** Variance - measure of spread in residence time distribution

---
**Repository:** https://github.com/Jimohmuktar/CRE-python
**Last Updated:** 2026
**Version:** 1.0
