"""Plot averages values.
"""
import pandas as pd
import matplotlib.pyplot as plt


def main() -> None:
    # Read data from text file using pandas
    data_exact = pd.read_csv("./expected/ref/exact_Lieb-Wu.tsv", sep="\t", header=0)
    data_dvmc = pd.read_csv("./cdmft.tsv", sep='\t', header=0)

    # Plot average density as a function of chemical potential
    plt.xlabel(r"$\mu$")
    plt.ylabel(r"$n$")
    plt.plot(data_exact['mu'], data_exact['ave_mu'], '--', lw=0.75, label="Lieb-Wu")
    plt.plot(data_dvmc['mu'], data_dvmc['ave_mu'],
             'o', color="C5", alpha=0.8, label="cdmft-dvmc")
    plt.legend()
    plt.savefig("./ave_mu.pdf")
    plt.show()

    return


if __name__ == "__main__":
    main()
