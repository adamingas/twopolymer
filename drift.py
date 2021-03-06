from __future__ import division

import sys
import os
import math
import linecache
import config
import argparse
import numpy as np
import timeit
import matplotlib.pyplot as plt

"""
Version 3.0.0Alpha
Changed to two polymer chains. Simulate run for a number of constants and find maximum separation squared.
"""

def walk(k,chi_element):
    """
    Performs the random walk.
    There are 4 possibilities of a walk, which are combinations of the binary states of hydrodynamic
    interactions and thermal noise.
    k is the weissenberg number
    """


    if str.upper(noise) in ("YES","Y","NOISE"):
        for j in xrange(runs):
            out = open("cnf{}_Wi{}_chi{}".format(j, k,chi_element), "w")
            radius = np.random.uniform()
            costh = np.random.uniform(-1,1,size=2)
            fi = np.random.uniform(size=2)*2*math.pi
            sinth = np.sqrt(1-costh**2)


            particles = [np.array(initial_positions[0:3]),np.array(initial_positions[3:6])]

            com1 = (particles[1] + particles[0])/2
            radius2 = np.linalg.norm(particles[1] - particles[0])

            particles.append(com1 + np.array([radius*math.cos(fi[0])*sinth[0],radius*math.sin(fi[0])*sinth[0],radius*costh[0]]))
            particles.append(particles[2] + np.array([radius2*math.cos(fi[1])*sinth[1],radius2*math.sin(fi[1])*sinth[1],radius2*costh[1]]))
            particles_interm = [None]*4
            polyvec1 = particles[1] - particles[0]
            polyvec2 = particles[3] - particles[2]
            # out.write("{} {} {} {} {} {} {} {} {}\n".format(time_step * (1), polyvec1[0], polyvec1[1], polyvec1[2],
            # np.linalg.norm(polyvec1),polyvec2[0], polyvec2[1], polyvec2[2],np.linalg.norm(polyvec2)))
            out.write("{} {} {} {} {}\n".format(0, str(particles[0]).strip("array([,])"),
                        str(particles[1]).strip("array([,])"), str(particles[2]).strip("array([,])"),
                        str(particles[3]).strip("array([,])")))

            for i in xrange(steps):


                # Constructing the oseen tensor (stokeslet matrix) and a list of the separation vectors between particles
                # The list of vectors goes like: pointing from particle ()- to particle ()
                # 1-2, 1-3, 2-3, 1-4, 2-4, 3-4
                oseen_tensor, sep_lst = oseen(particles)
                # Chi_element is divided by two because the time step is half in the intermediate step
                noise_vector = noise_producer(oseen_tensor, chi_element/2)
                # Updating new position of particles
                particles_interm[0] = particles[0] + time_step*function(0,1,2,3,sep_lst[0],sep_lst[5],particles,oseen_tensor,k)/2 + noise_vector[Part1]
                particles_interm[1] = particles[1] + time_step*function(1,0,2,3,-sep_lst[0],sep_lst[5],particles,oseen_tensor,k)/2 + noise_vector[Part2]
                particles_interm[2] = particles[2] + time_step*function(2,3,0,1,sep_lst[5],sep_lst[0],particles,oseen_tensor,k)/2 + noise_vector[Part3]
                particles_interm[3] = particles[3] + time_step*function(3,2,0,1,-sep_lst[5],sep_lst[0],particles,oseen_tensor,k)/2 + noise_vector[Part4]

                oseen_tensor, sep_lst = oseen(particles_interm)
                noise_vector = noise_producer(oseen_tensor, chi_element)

                particles[0] = particles[0] + time_step*function(0,1,2,3,sep_lst[0],sep_lst[5],particles_interm,oseen_tensor,k) + noise_vector[Part1]
                particles[1] = particles[1] + time_step*function(1,0,2,3,-sep_lst[0],sep_lst[5],particles_interm,oseen_tensor,k) + noise_vector[Part2]
                particles[2] = particles[2] + time_step*function(2,3,0,1,sep_lst[5],sep_lst[0],particles_interm,oseen_tensor,k) + noise_vector[Part3]
                particles[3] = particles[3] + time_step*function(3,2,0,1,-sep_lst[5],sep_lst[0],particles_interm,oseen_tensor,k) + noise_vector[Part4]


                # polyvec1 = particles[1] - particles[0]
                # polyvec2 = particles[3] - particles[2]
                # out.write("{} {} {} {} {} {} {} {} {}\n".format(time_step * (i + 1),polyvec1[0],polyvec1[1],polyvec1[2],np.linalg.norm(polyvec1),
                #                                     polyvec2[0], polyvec2[1], polyvec2[2], np.linalg.norm(polyvec2)))
                out.write("{} {} {} {} {}\n".format(time_step * (i+1), str(particles[0]).strip("array([,])"),
                str(particles[1]).strip("array([,])"),str(particles[2]).strip("array([,])"),str(particles[3]).strip("array([,])")))
            update_progress(j / (runs))
            out.close()

def function(ref, refpair, p3, p4, sepvec, sepvec2, particles, oseen_tensor, k):

    """

    :param ref: The particle on which the force acts on
    :param refpair: The particle on the same chain of the reference particle
    :param p3: The 3rd particle
    :param p4: The 4th particle

    :param sepvec: Separation vector from reference particle to its pair on the same chain
    :param sepvec2: Separation vector between the other two particles from smallest particle number
                    to largest (1 to 2, or, 3 to 4)
    :param particles:
    :return: The velocity of the reference particle at the current time step
    """
    magvec = np.linalg.norm(sepvec)**2
    magvec2 = np.linalg.norm(sepvec2)**2

    hydro_forces = ((sepvec) / (2 * (1 - np.linalg.norm(sepvec) ** 2))) + (
        -oseen_tensor[refpair*3:(refpair+1)*3, ref*3:(ref+1)*3].dot(sepvec)/(1-magvec)) + (oseen_tensor[3*p3:3*(p3+1), ref*3:(ref+1)*3].dot(sepvec2)
        - oseen_tensor[3*p4:3*(p4+1), ref*3:(ref+1)*3].dot(sepvec2))/(1-magvec2)

    vector = np.array([k * particles[int(ref)][1], 0, 0]) + hydro_forces



    return vector
def oseen(particles):
    # Construct list of separations between particles
    sep_lst = []
    oseen_tensor = np.zeros((12,12))
    for i in range(4):
        for j in range(i):
            # separation vector
            svec = particles[i] - particles[j]
            # separation magnitude
            smagn = np.linalg.norm(svec)**2
            sep_lst.append(svec)
            # The units might not be correct
            oseen_tensor[i*3:(i+1)*3,j*3:(j+1)*3] = (3*ar_ratio/8)*(np.identity(3)*(smagn + 2*epsilon_squared)+ (np.outer(svec,svec)))/(smagn +epsilon_squared)**(1.5)

    # Filling in the diagonal elements of the oseen tensor and returning the symmetric matrix

    np.fill_diagonal(oseen_tensor,1/2)
    return np.maximum(oseen_tensor,oseen_tensor.T),sep_lst
def noise_producer(matrix,chi_element):
    """
    This method calculates the noise exerted on the 2 spheres and returns a 3 dimensional vector.
    :param matrix: Stokeslet matrix (oseen tensor)
    :return: Random noise vector
    """
    # Creating the sigma matrix
    sigma_m = np.zeros((12,12))
    sigma_m[0,0] = math.sqrt(matrix[0,0])
    for i in range(12):
        for j in range(i):
            sigma_m[i,j] = matrix[i,j]
            k = 0
            while k < j:
                sigma_m[i,j] -= sigma_m[i,k]*sigma_m[j,k]
                k += 1

            sigma_m[i,j] = sigma_m[i,j]/sigma_m[j,j]
        sigma_m[i,i] = math.sqrt(matrix[i,i] - np.sum(sigma_m[i,0:i]**2))


    return sigma_m.dot(np.random.normal(0,math.sqrt(chi_element* time_step),12))
def analyse():
    """
    This function goes to the folder of the previously run simulation, and averages over
    all the runs (which, because of thermal noise, are different). It writes the mean
    squared separation and mean separation of each constant simulated in a file. It also
    creates a file with the maximum separation of every constant and at which time it occurred
    """

    global results
    os.chdir("Drift_cnf:{}_Wi:{}-{}_chi:{}-{}_steps:{}_ts:{}_ra{}_noise{}".format(runs,constant[0],constant[-1],
        chi[0],chi[-1],steps,time_step,ar_ratio,noise))

    # max_file = open("AMax_File_con"
    #                 ":{}-{}_numc:{}_h{}_s{}_ts{}_ra{}_n{}".format(
    #     constant[0], constant[-1], len(constant), hydro, steps, time_step, ar_ratio, noise), "w")
    # max_file.write("#Constant, Maxseparation over initial separation, time of max separation\n")


    for x in chi:
        for v in constant:
            #Rfile = open(os.getcwd() + "/Results.out", "a")
            # results is the file where all the results will be printed out


            #The following way of reading files is because of a limitation
            #in the number of files a computer can have open at the same time
            thousands_of_runs = int(math.ceil(runs / 1000))

            # Reads every thousand runs of a simulation
            for k in range(thousands_of_runs):
                # Opens the first 1000 runs in a dictionary, then opens the next 1000 and so on.
                filedata = {i: open("cnf{}_Wi{}_chi{}".format(i,v,x), "r") for i in xrange(k * 1000, min(runs, (k + 1) * 1000))}
                # Mean separation and Mean square separation lists that contain temporary files
                # with the respective values for every thousand runs. They are deleted afterwards
                # Averages_list contains temporary files which hold the average values of various parameters for every
                # thousand runs of the simulation


                # Adding squared separation and separation together
                # to average noise

                num = 0
                for file in filedata.values():
                    num +=1

                    results = open("RES_{}_{}_{}_{}_{}_{}_{}.out".format(num, v, x, steps, time_step, ar_ratio, noise),"w")

                    results.write("# Time-step, Mean Separation1, Mean Squared Separation1, Mean Separation2, Mean Squared Separation2, Mean COM separation, Angle 1, Angle 2 \n")

                    for lines in xrange(steps + 1):
                        token = str.split(file.readline())
                        t = float(token[0])
                        particles = []
                        for i in range(4):
                            particles.append(np.array([float(token[3*i +j +1]) for j in range(3)]))
                        # Polymer separation vector and centre of mass vecotr
                        polyvec1 = particles[1] - particles[0]
                        polymagn1 = np.linalg.norm(polyvec1)
                        polyvec2 = particles[3] - particles[2]
                        polymagn2 = np.linalg.norm(polyvec2)
                        comvec = (particles[3] + particles[2] - particles[1] - particles[0])/2
                        commag = np.linalg.norm(comvec)
                        # Angle with shear axis
                        angle1 = np.degrees(np.arccos(np.clip(np.dot(polyvec1,np.array([1,0,0]))/polymagn1,-1.0,1.0)))
                        angle2 = np.degrees(np.arccos(np.clip(np.dot(polyvec2, np.array([1, 0, 0])) / polymagn2, -1.0, 1.0)))


                        sepsq1 = polymagn1 ** 2
                        sep1 = polymagn1
                        sepsq2 = polymagn2 ** 2
                        sep2 = polymagn2
                        comsep = commag
                        totangle1 = angle1
                        totangle2 = angle2
                        results.write(
                        "{} {} {} {} {} {} {} {}\n".format(t, sep1, sepsq1, sep2, sepsq2,
                                                       comsep,angle1, angle2))

                    results.close()
                update_progress(num / runs)
    os.chdir("..")

def plot():
    """
    This method calculates the maximum separation between the polymers and the maximum separation
    between particles of the same polymer
    :return:
    """
    os.chdir("Drift_cnf:{}_Wi:{}-{}_chi:{}-{}_steps:{}_ts:{}_ra{}_noise{}".format(runs,constant[0],constant[-1],
        chi[0],chi[-1],steps,time_step,ar_ratio,noise))

    for j,chi_element in enumerate(chi):

        for i,wi in enumerate(constant):
            plt.figure()
            color = iter(plt.cm.rainbow(np.linspace(0, 1,runs)))
            for cnf in range(runs):
                results_array = np.loadtxt(
                    "RES_{}_{}_{}_{}_{}_{}_{}.out".format(cnf+1,wi, chi_element, steps, time_step, ar_ratio, noise), skiprows=1)

                c = next(color)
                plt.plot(np.arange(steps+1) * time_step,results_array[:,5] , c=c)
            plt.legend()
            plt.ylabel("Centres of Mass separation")
            plt.xlabel("Time elapsed")
            plt.title("C.o.M. Separation vs Time")
            plt.savefig("Drift_ts:{}_steps:{}_conf:{}_Wi:{}.png".format(time_step,steps,runs,wi))




            print("wi number {} done, {} left".format(wi, len(constant) - i - 1))
        print("chi number {} done, {} left".format(chi_element, len(chi) -j-1))
    os.chdir("..")
def average_sepparation():
    print ("Average Sepparation calculation")
    os.chdir("Shear_Wi:{}-{}_chi:{}-{}_steps:{}_ts:{}_ra{}_noise{}".format(constant[0], constant[-1],
        chi[0], chi[-1], steps, time_step, ar_ratio,noise))


    for j,chi_element in enumerate(chi):
        avsep_file = open(("Av.Sep_hydro:{}_chi:{}_Wi:{}-{}_numWi:{}").format(hydro, chi_element,
                                                                constant[0],constant[-1],len(constant)),"w")
        for i,wi in enumerate(constant):
            wisep_array = np.loadtxt("RES_{}_{}_{}_{}_{}_{}.out".format(wi, chi_element, steps, time_step, ar_ratio, noise))


            average1 = np.mean(wisep_array[500:,1])
            average2 = np.mean(wisep_array[500:, 3])
            mean = (average1 + average2)/2
            avsep_file.write("{} {}\n".format(wi,mean))
            print("wi number {} done, {} left".format(wi, len(constant) - i - 1))
        print("chi number {} done, {} left".format(chi_element, len(chi) -j-1))
    os.chdir("..")

def angle():
    print ("Angle Calculation")
    os.chdir("Shear_Wi:{}-{}_chi:{}-{}_steps:{}_ts:{}_ra{}_noise{}".format(constant[0], constant[-1],
        chi[0], chi[-1], steps, time_step, ar_ratio,noise))


    for j,chi_element in enumerate(chi):
        avsep_file = open(("Angle_hydro:{}_chi:{}_Wi:{}-{}_numWi:{}").format(hydro, chi_element,
                                                                constant[0],constant[-1],len(constant)),"w")
        for i,wi in enumerate(constant):
            file_array = np.loadtxt("RES_{}_{}_{}_{}_{}_{}.out".format(wi, chi_element, steps, time_step, ar_ratio, noise))

            angle_array1 = np.mean(file_array[500:,6])
            angle_array2 = np.mean(file_array[500:,7])
            avsep_file.write("{} {}\n".format(wi,(angle_array1+angle_array2)/2))
            print("wi number {} done, {} left".format(wi, len(constant) - i - 1))
        print("chi number {} done, {} left".format(chi_element, len(chi) -j-1))
    os.chdir("..")

def simulate():
    """
    For version 1 the creation of folder is inside the for loop.
    For version 2 all simulations are created in a single folder but separate files.
    :return:
    """
    #This code simulates the walk of the polymer and stores the max separation in a file

    try:
        os.mkdir("Drift_cnf:{}_Wi:{}-{}_chi:{}-{}_steps:{}_ts:{}_ra{}_noise{}".format(runs,constant[0],constant[-1],
        chi[0],chi[-1],steps,time_step,ar_ratio,noise))
    except OSError as e:
        if e.errno != 17:
            raise

    os.chdir("Drift_cnf:{}_Wi:{}-{}_chi:{}-{}_steps:{}_ts:{}_ra{}_noise{}".format(runs,constant[0],constant[-1],
        chi[0],chi[-1],steps,time_step,ar_ratio,noise))



    for o,k in enumerate(chi):
        for n,j in enumerate(constant):

            walk(j,k)
            print("wi number {} done, {} left".format(j,len(constant)-n -1))
        print("chi number {} done, {} left".format(k, len(chi) -o-1))
    os.chdir("..")

def update_progress(progress):
    barLength = 10  # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength * progress))
    text = "\rPercent: [{0}] {1}% {2}".format("#" * block + "-" * (barLength - block), progress * 100, status)
    sys.stdout.write(text)
    sys.stdout.flush()

if __name__ == "__main__":
    """
    Program gets variables from config file. 
    To run the program, arguments have to be passed while initialising it in the form of flags.
    The -walk argument runs the random walk simulation using the parameters of the config file
    The -analyse argument finds the folder that matches the parameters of the config file
    and subsequently runs through the files to produce a mean separation squared and max separation
    files for every constant.
    The -max argument is used only when there is no noise in the simulation so as to avoid running
    the computationally expensive -analyse method 
    """
    #Constants used in the simulation
    steps = int(config.steps)
    runs = int(config.runs)
    constant = config.constant #Weissenberg numbers in a numpy array
    time_step = float(config.time_step)
    initial_positions = config.initial_positions
    # init_separation = math.sqrt(xinitial**2 + yinitial**2 + zinitial**2)
    hydro = str(config.hydro)
    ar_ratio = float(config.ar_ratio)
    noise = str(config.noise)
    chi = config.chi
    epsilon_squared = 4*ar_ratio*ar_ratio
    # Lists that signify particle position in the Oseen tensor to make code more readable
    Part1 = range(0,3)
    Part2 = range(3,6)
    Part3 = range(6,9)
    Part4 = range(9,12)
    parser = argparse.ArgumentParser(description="Program version 2"
                                                 "The program simulates the motion of a polymer in shear flow.\n "
                                                 "The model is of a finite extensibility non-linear elastic spring"
                                                 "(FENE).\n"
                                                 "Parameters of the simulation can be found in the config.py file.\n"
                                                 "\nVersion 1: Has no thermal fluctuations so analyser doesn't do anything\n"
                                     "Version 2: There are hydrodynamic interactions between the two ends of the polymer chain.\n"
                                     "Version 3: There is thermal noise acting on the polymer")
    # parser.add_argument("echo", help="echo the string you use here")
    parser.add_argument("-a", "--analyse", help="Run analyser. For version 1 this does nothing (v1) since there "
                                                "are no thermal fluctuations.",
                        action="store_true")
    parser.add_argument("-w", "--walk", help="Simulate walks with parameters of config.py file", action="store_true")

    parser.add_argument("-p", "--plot", help="Plot configurations", action="store_true")

    parser.add_argument("-av", "--average", help="Averages separation for each weisenberg number", action="store_true")

    parser.add_argument("-an", "--angle", help="Averages angle between separation vector and x-axis for each weisenberg number", action="store_true")

    args = parser.parse_args()
    # print args.echo
    if args.walk:
        simulate()
    if args.analyse:
        print timeit.timeit(analyse, number=1)
    if args.average:
        average_sepparation()
    if args.angle:
        angle()
    if args.plot:
        plot()