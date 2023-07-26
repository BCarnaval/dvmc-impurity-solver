"""Plot averages.
"""
import pandas as pd
import matplotlib.pyplot as plt


def main() -> None:
    # Read data from text file using pandas
    data = pd.read_csv("./cdmft.tsv", sep='\t', header=0)
    x, y = data['mu'], data['ave_mu']

    # Plot average density as a function of chemical potential
    plt.xlabel(r"$\mu$")
    plt.ylabel(r"$n$")
    plt.plot(x, y, 'o--', label="dVMC", alpha=0.6)
    plt.legend()
    plt.show()

    return


if __name__ == "__main__":
    main()
