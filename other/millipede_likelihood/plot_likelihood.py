import numpy as np
import matplotlib.pyplot as plt
import scipy.special
from scipy.stats import poisson

# the Poisson PDF
def poisson_pdf(k, lam):
    return np.power(lam, k) * np.exp(-lam) / scipy.special.factorial(k)

# the Poisson CDF
def poisson_cdf(k, lam):
    return scipy.special.gammaincc(k+1, lam) # already regularized (thanks Python)

x_vals = np.linspace(0,20,15000)
lam = 12

pdf = poisson_pdf(x_vals, lam)
cdf = poisson_cdf(x_vals, lam)

pdf_nll = -np.log(pdf)
pdf_nll_prime = np.gradient(pdf_nll, x_vals[1] - x_vals[0])
pdf_nll_primeprime = np.gradient(pdf_nll_prime, x_vals[1] - x_vals[0])

cdf_nll = -np.log(cdf)
cdf_nll_prime = np.gradient(cdf_nll, x_vals[1] - x_vals[0])
cdf_nll_primeprime = np.gradient(cdf_nll_prime, x_vals[1] - x_vals[0])

min_index = np.argmin(-np.log(pdf))
# print("Min Index {}, X Value There is {}".format(min_index, x_vals[min_index]))
# max_index = np.argmax(pdf)
# print("Max PDF Index {}, X Value There is {}".format(max_index, x_vals[max_index]))


fig, ax = plt.subplots(1,3,figsize=(12,4))

# likelihood
ax[0].plot(x_vals, pdf, label=r'PDF')
ax[0].plot(x_vals, cdf, linestyle='--', label=r'CDF')
ax[0].set_xlabel('k')
ax[0].set_ylabel(r'Likelihood: $\mathcal{L}$', fontsize=12)
ax[0].legend()

# negative log likelihood
ax[1].plot(x_vals, pdf_nll, label='PDF')
ax[1].plot(x_vals, cdf_nll, linestyle='--', label='CDF')
ax[1].set_xlabel('k')
ax[1].set_ylabel(r'Negative Log Likelihood: ($-\log\mathcal{L}$)', fontsize=12)


# second derivative of NLL
ax[2].plot(x_vals, pdf_nll_prime, label='PDF')
ax[2].plot(x_vals, cdf_nll_prime, linestyle='--', label='CDF')
ax[2].set_xlabel('k')
ax[2].set_ylabel(r"Second Derivative of NLL: ($-\log\mathcal{L}$)''", fontsize=12)


plt.tight_layout()
plt.subplots_adjust(wspace=0.3)
fig.suptitle("Poisson, $\lambda$ = {}".format(lam))
fig.subplots_adjust(top=0.9)

fig.savefig('likelihood.png')
del fig, ax