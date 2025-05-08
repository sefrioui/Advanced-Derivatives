# Advanced Derivatives – Implied Volatility Surface Construction

This repository contains the implementation of an implied volatility surface construction for the EURO STOXX 50 (STOXX50E) index using the **Andreasen-Huge** algorithm, as part of an assignment for the *Advanced Derivatives* course.

## Project Summary

We construct an implied volatility surface for the STOXX50E index using market data from March 2010. The approach follows the method introduced by Andreasen and Huge (2011), using a local volatility model and Dupire's forward PDE to interpolate and extrapolate implied volatilities.

The key assumptions are:
- Underlying index price is set to \( S_0 = 100 \)
- No interest rate or dividend yield (i.e., \( r = q = 0 \))
- Implied volatilities are available for selected strikes and maturities

##  Methodology

1. **Data Preprocessing**
   - Load market data from Excel (strike, maturity, implied vol)
   - Discretize the strike space with interval \( \Delta K \)

2. **Initial Option Prices**
   - Compute option prices at \( t = 0 \) using BSM and implied volatilities

3. **Dupire PDE and Finite Differences**
   - Use the Dupire forward PDE to project option prices across maturities
   - Apply finite difference methods to solve the PDE

4. **Calibration**
   - Minimize the squared error between market prices and Dupire model prices
   - Calibrate local volatility surface \( \tilde{\sigma}(S, t) \)

5. **Implied Vol Surface**
   - Invert BSM to get implied volatilities from computed prices
   - Plot implied volatility surface as a function of strike and time

6. **Interpolation of Missing Maturities**
   - Use tridiagonal system to interpolate intermediate maturities
   - Propagate prices using Andreasen-Huge method and back out volatilities

## 📈 Output

- Final result is a 3D implied volatility surface over strike and time
- See `Figure 1` for the plotted surface

## 🛠️ Tools & Libraries

- Python / MATLAB (depending on implementation)
- Excel for initial data loading
- Matplotlib or similar for visualization


## 📅 Date

November 2024

## 📄 Reference

- Andreasen, J., & Huge, B. (2011). *Volatility interpolation*. Risk Magazine.
