#!/usr/bin/env python3
import sys


def main():
    Np = Ns * Ns * 2
    f = open('greenone.def', 'w')
    f.write(
        "===============================\n" +
        "NCisAjs         "+str(Np)+"\n" +
        "===============================\n" +
        "======== Green functions ======\n" +
        "===============================\n"
    )
    for ii in range(Ns):
        for ss in range(2):
            for jj in range(Ns):
                f.write(
                    '%5s %6s %6s %6s \n' % (str(ii), str(ss), str(jj), str(ss)))
    f.close()

    Np = Ns * Ns * Ns * Ns * 2 * 2
    f = open('greentwo.def', 'w')
    f.write(
        "===============================\n" +
        "NCisAjsCktAltDC     "+str(Np)+"\n" +
        "===============================\n" +
        "=Green functions for Sq AND Nq=\n" +
        "===============================\n"
    )
    for ii in range(Ns):
        for ss in range(2):
            for jj in range(Ns):
                for ii2 in range(Ns):
                    for ss2 in range(2):
                        for jj2 in range(Ns):
                            # for ss2 in range(2):
                            f.write('%5s %6s %6s %6s %6s %6s %6s %6s\n' % (str(ii), str(
                                ss), str(jj), str(ss),  str(ii2), str(ss2), str(jj2), str(ss2)))
    f.close()

    return


if __name__ == "__main__":
    if len(sys.argv) == 2:
        Ns = int(sys.argv[1])
    else:
        print("example:\n$ makeFullGreen.py 8")
        sys.exit()

    main()
