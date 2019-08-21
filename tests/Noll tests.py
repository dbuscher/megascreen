# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.2'
#       jupytext_version: 1.1.6
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
import numpy as np
import matplotlib.pyplot as plt
from joblib import Memory, Parallel, delayed
# Incantation to get latest version
import test_Noll
import importlib
importlib.reload(test_Noll)
from test_Noll import *


# %% [markdown]
# # Some wierdness in Table IV of Noll 1976
# Comparing values in the table with those computed using the formulae given in the paper shows errors several times the least significant decimal place printed in the table for the first few rows. Given that each line of the table is computed from the previous line, you would expect the errors to accumulate as you go down the table, but in general the reverse happens, and the error reduces as we go down the table.

# %%
def check_table_iv(const=0.023):
    resid=noll_piston_residual_analytic(const=const)
    for i in range(len(noll_table_iv)):
        print("{:8.5f} {:8.5f} {:8.2e}".format(noll_table_iv[i],resid,resid-noll_table_iv[i]))
        i=i+1
        n,m=jtonm(i+1)
        resid-=noll_covariance_analytic(n,n,const=const)
check_table_iv()

# %% [markdown]
# We can reduce the error to acceptable values by setting the constant in front of the power spectrum to 0.02284 instead of the value 0.023 given in Noll's paper. 

# %%
check_table_iv(0.02284) 

# %% [markdown]
# Using a higher-precision value of this constant (whether using Fried's constant of 6.88 as the only low-precision number, or computing it to high precision by evaluating Fried's formula) yields an increased discrepancy. 

# %%
fried=2*((24/5)*gamma(6/5))**(5/6)
print("fried=",fried)
const=fried*gamma(11/6)*(5/6)/(2*np.pi**(8/3)*gamma(1/6))
print("const =",const)
check_table_iv(const)

# %% [markdown]
# So did why did Noll use a value for this constant which is neither the value quoted in the his own paper nor a more precise value?

# %% [markdown]
# # Test the simulator against Winker 1991

# %%
# Use the joblib Memory decorator to avoid having to recompute long-running calculations
memory = Memory("cache", verbose=0)
@memory.cache
def MemoWinker(
    diameter=32,
    L0Min=16,
    L0Max=8000,
    numL0=20,
    numIter=100,
    maxRadial=2,
    nfftOuter=256,
    nfftInner=256,
    randomSeed=12345,
):
    """Memoised wrapper for the Winker function"""
    np.random.seed(randomSeed)
    return Winker(
        diameter, L0Min, L0Max, numL0, numIter, maxRadial, nfftOuter, nfftInner
    )


t = MemoWinker(numIter=1000, diameter=100, randomSeed=5432112)


# %%
def PlotWinker(t, diameter):
    """Reproduce Fig 2 from Winker 1991"""
    # Plot theoretical values
    L0=np.logspace(0,3,50)
    resid=winker_residual(L0,2)
    for n in range(len(resid)):
        plt.loglog(1/(L0),resid[n],ls="dotted")
    # Plot simulation values
    L0 = 2 * t["L0"] / diameter
    for z in [0, 2, 5]:
        plt.loglog(1 / L0, t["Z" + str(z)], label="Z_max=" + str(z+1))
    plt.legend()

def winker_residual(L0,max_n):
    resid=np.empty((max_n+1,len(L0)))
    for i in range(len(L0)):
        L=L0[i]
        resid[0,i]=winker_piston_residual(L0=L)
        for n in range(1,max_n+1):
            resid[n,i]=resid[n-1,i]-winker_variance_quad(n,R=1,r0=2,L0=L)*(n+1)
    return resid

plt.figure(figsize=(12,10))
plt.xlabel(r"$R/L_0$")
plt.ylabel("$Z_j$")
plt.title("Simulation of Winker Fig 2")
PlotWinker(t, 100)

# %%
