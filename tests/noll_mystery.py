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

# %% [markdown]
# # Notes on Table IV of Noll 1976
#
# This is only tangentially related to the MegaScreen program, but came up when testing the accuracy of the phase screen generator (actually when testing the accuracy of the evaluation of the integrals in Winker 1991, used to compare with the simulations). These notes are here mostly to remind me of what I found.
#
# Comparing values in the first few rows of Table IV of the Noll 1976 paper with those computed using the formulae given in the paper shows errors several times the least significant decimal place printed in the table. Given that each line of the table is computed from the previous line, you would expect the errors to accumulate as you go down the table, but in general the reverse happens, and the error reduces as we go down the table.

# %%
import numpy as np
from scipy.special import gamma
from zernike import jtonm

pi=np.pi

def noll_covariance_analytic(ni, nj, m=None, const=0.0229):
    """Evaluate covariance of Zernikes using equation A2 from Noll 1976

    Returns covariance of radial order ni and nj  and azimuthal
    order m.  Note that if m_j != m_i then the correlation is zero. 
    Does not check if m is allowed."""
    if m == None:
        m = ni % 2  # Lowest valid m
    return (
        const
        * 2 ** (-5 / 3)
        * (-1) ** ((ni + nj - 2 * m) / 2)
        * np.sqrt((ni + 1) * (nj + 1))
        * pi ** (8.0 / 3)
        * gamma(14.0 / 3.0)
        * gamma((ni + nj - 5.0 / 3) / 2)
        / (
            gamma((ni - nj + 17 / 3.0) / 2)
            * gamma((nj - ni + 17.0 / 3) / 2)
            * gamma((ni + nj + 23 / 3.0) / 2)
        )
    )


def noll_piston_residual_analytic(p=8 / 3, const=0.0229):
    """Evaluates the analytic solution for the piston-removed
    residual atmospheric aberrations using the last equation in 
    the Appendix of Noll 1976

    This should give :math:`\Delta_1` from Table IV of that paper
    but the values deviate by more than 1 significant figure.
    """
    return (
        pi
        * gamma(p + 2)
        / (
            2 ** p
            * (gamma((p + 3) / 2)) ** 2
            * gamma((p + 5) / 2)
            * gamma((1 + p) / 2)
            * np.sin((pi / 2) * (p - 1))
        )
        * (2 * pi) ** (5 / 3)
        * (const * 2 * pi)
        * 2 ** (-5 / 3)
    )


# Table IV from Noll 1976
noll_table_iv = [
    1.0299,
    0.582,
    0.134,
    0.111,
    0.0880,
    0.0648,
    0.0587,
    0.0525,
    0.0463,
    0.0401,
    0.0377,
    0.0352,
    0.0328,
    0.0304,
    0.0279,
    0.0267,
    0.0255,
    0.0243,
    0.0232,
    0.0220,
    0.0208,
]


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
# Using a higher-precision value of this constant (whether using Fried's constant of 6.88 as the only low-precision number, or computing it to high precision by evaluating Fried's formula) yields almost as large a discrepancy as using the value of 0.23. 

# %%
fried=2*((24/5)*gamma(6/5))**(5/6)
print("fried=",fried)
const=fried*gamma(11/6)*(5/6)/(2*np.pi**(8/3)*gamma(1/6))
print("const =",const)
check_table_iv(const)

# %% [markdown]
# So did why did Noll use a value for this constant which is neither the value quoted in the his own paper nor a more precise value?
