"""
Last changed: 19.10.2019 XX:YY by Erlend Tiberg North
"""
# Imports
import os
import numpy as np
import matplotlib.pyplot as plt

def options():
    # This function displays the cases available for evaluation
    # and asks the user for a choice of file
    print("What do you wish to evaluate?")
    print("1:   lambda.txt")
    print("2:   integrationpoints.txt")
    print("3:   montecarlo.txt")
    print("4:   Timings")
    print()
    print("-1:   Quit Program")

    while True:
        try:
            filenum = int(input("\nWhat filename do you want to evaluate? "))
            break
        except ValueError:
            print("\nCould not convert your input to int")
            print("Please specify file using corresponding number")
    return filenum

def eval_lambda(filename):

    n = []
    la = []
    error_legendre = []

    with open(filename) as infile:
        infile.readline()
        lines = infile.readlines()
        for line in lines:
            words = line.split(" , ")
            n.append(int(words[0]))
            la.append(float(words[1]))
            error_legendre.append(float(words[2]))

    n = np.array(n)
    la = np.array(la)
    error_legendre = np.array(error_legendre)

    plt.plot(la,error_legendre)
    plt.title("Plot of error ($n=%d$) with varying $\lambda$" % n[0])
    plt.xlabel("$\lambda$")
    plt.ylabel("Error")
    plt.grid()
    plt.savefig("../images/error-" + filename[:-4] + ".png")
    plt.show()

def eval_integrationpoints(filename):

    n = []
    la = []
    error_legendre = []
    error_laguerre = []

    with open(filename) as infile:
        infile.readline()
        lines = infile.readlines()
        for line in lines:
            words = line.split(" , ")
            n.append(int(words[0]))
            la.append(float(words[1]))
            error_legendre.append(float(words[2]))
            error_laguerre.append(float(words[3]))

    plt.plot(n, error_legendre, label="Gauss-Legendre")
    plt.plot(n, error_laguerre, label="Gauss-Legendre")
    plt.title("Plot of error ($\lambda=%d$) with varying n" % la[0])
    plt.xlabel("n")
    plt.ylabel("Error")
    plt.grid()
    plt.savefig("../images/error-" + filename[:-4] + ".png")
    plt.show()

def eval_montecarlo(filename):

    n = []
    la = []
    error_BMC = []
    error_SMC = []
    error_PSMC = []
    time_span_BMC = []
    time_span_SMC = []
    time_span_PSMC = []

    with open(filename) as infile:
        infile.readline()
        lines = infile.readlines()
        for line in lines:
            words = line.split(" , ")
            n.append(int(words[0]))
            la.append(float(words[1]))
            error_BMC.append(float(words[2]))
            error_SMC.append(float(words[3]))
            error_PSMC.append(float(words[4]))
            time_span_BMC.append(float(words[5]))
            time_span_SMC.append(float(words[6]))
            time_span_PSMC.append(float(words[7]))

        la = np.array(la)
        error_BMC = np.array(error_BMC)
        error_SMC = np.array(error_SMC)
        error_PSMC = np.array(error_PSMC)
        time_span_BMC = np.array(time_span_BMC)
        time_span_SMC = np.array(time_span_SMC)
        time_span_PSMC = np.array(time_span_PSMC)

        plt.subplot(2, 1, 1)
        plt.plot(n, error_BMC, label="Brute Force Monte Carlo")
        plt.plot(n, error_SMC, label="Spherical Monte Carlo")
        plt.plot(n, error_PSMC, label="Parallized Spherical Monte Carlo")
        plt.title("Plot of error ($\lambda=%d$) with varying n" % la[0])
        plt.xlabel("n")
        plt.ylabel("Error")
        plt.xlim(-0.5, 10500)
        plt.ylim(-0.01, 0.2)
        plt.grid()
        plt.legend()

        plt.subplot(2, 1, 2)
        plt.plot(n, error_BMC, label="Brute Force Monte Carlo")
        plt.plot(n, error_SMC, label="Spherical Monte Carlo")
        plt.plot(n, error_PSMC, label="Parallized Spherical Monte Carlo")
        plt.xlabel("n")
        plt.ylabel("Error")
        plt.xlim(10500, 100000000)
        plt.ylim(-0.0001, 0.001)
        plt.grid()
        plt.legend()
        plt.savefig("../images/error-montecarlo.png")
        plt.show()

def eval_timings(filename1, filename2):

    la = []
    error_legendre = []
    error_laguerre = []
    time_span_gauss_legendre = []
    time_span_gauss_laguerre = []
    error_BMC = []
    error_SMC = []
    error_PSMC = []
    time_span_BMC = []
    time_span_SMC = []
    time_span_PSMC = []

    with open(filename1) as infile1:
        infile1.readline()
        lines = infile1.readlines()
        for line in lines:
            words = line.split(" , ")
            la.append(float(words[1]))
            error_legendre.append(float(words[2]))
            error_laguerre.append(float(words[3]))
            time_span_gauss_legendre.append(float(words[4]))
            time_span_gauss_laguerre.append(float(words[5]))

    with open(filename2) as infile2:
        infile2.readline()
        lines = infile2.readlines()
        for line in lines:
            words = line.split(" , ")
            error_BMC.append(float(words[2]))
            error_SMC.append(float(words[3]))
            error_PSMC.append(float(words[4]))
            time_span_BMC.append(float(words[5]))
            time_span_SMC.append(float(words[6]))
            time_span_PSMC.append(float(words[7]))

    la = np.array(la)
    error_legendre = np.array(error_legendre)
    error_laguerre = np.array(error_laguerre)
    time_span_gauss_legendre = np.array(time_span_gauss_legendre)
    time_span_gauss_laguerre = np.array(time_span_gauss_laguerre)
    error_BMC = np.array(error_BMC)
    error_SMC = np.array(error_SMC)
    error_PSMC = np.array(error_PSMC)
    time_span_BMC = np.array(time_span_BMC)
    time_span_SMC = np.array(time_span_SMC)
    time_span_PSMC = np.array(time_span_PSMC)

    plt.subplot(2, 1, 1)
    plt.plot(time_span_gauss_legendre, error_legendre, label="Gauss-Legendre")
    plt.plot(time_span_gauss_laguerre, error_laguerre, label="Gauss-Laguerre")
    plt.title("Plot of error ($\lambda=%d$) as function of time" % la[0])
    plt.xlabel("Time (s)")
    plt.ylabel("Error")
    plt.grid()
    plt.legend()
    plt.xlim(-0.5, 10)
    plt.ylim(-0.001, 0.05)


    plt.subplot(2, 1, 2)
    plt.plot(time_span_BMC, error_BMC, label="Brute Force Monte Carlo")
    plt.plot(time_span_SMC, error_SMC, label="Spherical Monte Carlo")
    plt.plot(time_span_PSMC, error_PSMC, label="Parallized Spherical Monte Carlo")
    plt.xlabel("Time (s)")
    plt.ylabel("Error")
    plt.grid()
    plt.legend()
    plt.xlim(-0.5, 10)
    plt.ylim(-0.001, 0.05)
    plt.savefig("../images/method-timings-small.png")
    plt.show()

    plt.subplot(2, 1, 1)
    plt.plot(time_span_gauss_legendre, error_legendre, label="Gauss-Legendre")
    plt.plot(time_span_gauss_laguerre, error_laguerre, label="Gauss-Laguerre")
    plt.title("Plot of error ($\lambda=%d$) as function of time" % la[0])
    plt.xlabel("Time (s)")
    plt.ylabel("Error")
    plt.grid()
    plt.legend()
    plt.xlim(10, 800)
    plt.ylim(-0.001, 0.008)

    plt.subplot(2, 1, 2)
    plt.plot(time_span_BMC, error_BMC, label="Brute Force Monte Carlo")
    plt.plot(time_span_SMC, error_SMC, label="Spherical Monte Carlo")
    plt.plot(time_span_PSMC, error_PSMC, label="Parallized Spherical Monte Carlo")
    plt.xlabel("Time (s)")
    plt.ylabel("Error")
    plt.grid()
    plt.legend()
    plt.xlim(10, 800)
    plt.ylim(-0.001, 0.008)
    plt.savefig("../images/method-timings-large.png")
    plt.show()

def main():
    # Variables for setup
    Running = True
    cases = {
        1: "lambda.txt",
        2: "integrationpoints.txt",
        3: "montecarlo.txt",
        4: "timings",
        -1: "quit"
    }

    # Welcome
    print("\nHello, and welcome to data evaluation\n")


    while Running: # Main loop for program

        #Choose file
        while True:
            case = cases.get(options())
            print("Case: " + case)
            if case == "quit":
                break
            elif case == "timings":
                if not os.path.isfile("integrationpoints.txt"):
                    print("integrationpoints.txt was not found. \
                    Please verify it's existence.")
                    print("NB! This option requires the files \
                    integrationpoints.txt and montecarlo.txt.")
                elif not os.path.isfile("montecarlo.txt"):
                    print("montecarlo.txt was not found. \
                    Please verify it's existence.")
                    print("NB! This option requires the files \
                    integrationpoints.txt and montecarlo.txt.")
                else:
                    print("You have chosen %s!" % case)
                    break
            elif not os.path.isfile(case):
                print("The choosen filename was not found. Please verify \
                it's existence.")
            else:
                print("You have chosen %s!" % case)
                break

        if case == cases.get(1):
            eval_lambda(case)
        elif case == cases.get(2):
            eval_integrationpoints(case)
        elif case == cases.get(3):
            # EVALUATE FOR MONTECARLO
            eval_montecarlo(case)
        elif case == cases.get(4):
            eval_timings("integrationpoints.txt", "montecarlo.txt")
        elif case == "quit":
            print("Quitting program..")
            Running = False
        print("Evaluation complete. Returning to main menu..\n")

if __name__ == "__main__":
    main()
