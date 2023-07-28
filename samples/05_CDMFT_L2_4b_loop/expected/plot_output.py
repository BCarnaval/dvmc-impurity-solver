"""Plot averages values.
"""
import pandas as pd
import matplotlib.pyplot as plt


def main() -> None:
    # Read data from text file using pandas
    data_ed = pd.read_csv("./cdmft_ED.tsv", sep="\t", header=0)
    data_dvmc = pd.read_csv("./cdmft.tsv", sep='\t', header=0)

    # Plot average density as a function of chemical potential
    plt.xlabel(r"$\mu$")
    plt.ylabel(r"$n$")
    plt.plot(data_ed['mu'], data_ed['ave_mu'], '--', lw=0.75, label="ED")
    plt.plot(data_dvmc['mu'], data_dvmc['ave_mu'],
             'o', label="dVMC", color="C5", alpha=0.8)
    plt.legend()
    plt.savefig("./ave_mu.pdf", dpi=800)
    plt.show()

    return


if __name__ == "__main__":
    main()
