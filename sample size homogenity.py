# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 11:19:24 2022
@author: BRGUBO

Homogenity of variance:  power / Type II error calculation 
Based on Levene´s test

The goal of this app is to analyse the homogenity of variance of a continueous variable:
    
    
    Generete theorical histogram of both populations - visual approach.
    Calculate the error rate and stat power when a defined sample size were taken from the population
    Estimate the sample size needed to achieve a desired goal of ERROR / POWER
    
Version: 1.0


Version improvements:
    
    - all variables needs to be clear before running
    
Version: 1.1

    Capability for the user to define the P level (Type I error)
    
    
suggested future improvements:
    
    - all variables needs to be clear before running
    - exit() function with a bug 
    - histogram with a bell shape line
    


"""
version = "1.1"


print(chr(27) + "[2J")
print("Comparing Homogenity of Variance between 2 independent groups - Sample Size Analyser\n")
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from art import *
import winsound
import os
import sys

# reset all variables
from IPython import get_ipython;   
get_ipython().magic('reset -sf')



# clear screen
print('\033[H\033[J', end='')
# build title - ASC II title
art = text2art(f"Homogeinity of CV ", font='smslant')
print(art) # print title
print("Sample size analyser & estimation\n")





def error():
    # sound variables
    frequency = 2500  # Set Frequency To 2500 Hertz
    duration = 1000  # Set Duration To 1000 ms == 1 second

    art = text2art("ERROR", font='broadway')
    winsound.Beep(frequency, duration)
    print('\n')
    print(art)
    print('\nReview your input data - execution aborted!\n')
    sys.exit()
    
    

# manage p level input
def value_p(p_input):
    if int(p_input) == 1:
        p_opt = 0.01
        return 0.01
    elif int(p_input) == 2:
        p_opt = 0.05
        return p_opt
    elif int(p_input) == 3:
        p_opt = 0.1
        return p_opt
    else:
        print("ERROR - input only 1, 2 or 3 ")
        error()
        #p_opt = 0.05 # assume p=0.05 as std value
        



# input data ---------------------------------------------------

average_input = float(input("Input the population average: "))
sd_dev_input = float(input("Input the standard deviation: "))
input_diff = (float(input("Expected improvement of homogeinity with Probiotic (%): ")))
num_rep = int(input("Number of repetitions / group: "))
str_stat_power_level = input("Desired Stat Power (usually 80%): ")

                              
if str_stat_power_level == "":
   stat_power_level = 80
else:
    stat_power_level = float(str_stat_power_level)
                              
p_option = float(input('\nP value (type I error):\n[ 1 ] -> p=0,01\n[ 2 ] -> p=0,05\n[ 3 ] -> p=0,1\n'))
P_level = value_p(p_option)
#P_level = 0.05 # define the P level (Type I error)

# calculation of expected std deviation 
expect_sd_dev = sd_dev_input * (1 - input_diff/100)


# calculation of coeficient of variation 
CV_control = round(sd_dev_input/average_input * 100,2)
CV_treatment = round(expect_sd_dev/average_input * 100,2)


rep=3000  # number of simulations 
a = 0
x = 0


# ------------------------------------------------------------

for i in range(rep):

    simul_control = np.random.normal(average_input,sd_dev_input,num_rep) # generate control group
    simul_treatment = np.random.normal(average_input,expect_sd_dev,num_rep) # generate control group

    stat, p = stats.levene(simul_control,simul_treatment )

    if p<= P_level:
        result = "*"
        a = a + 1
    else:
        result = " "
    
    # print the first 20 trials
    if x < 21:
        print(f"trial nº: {x} --> Stat = {round(stat,1)} , P = {round(p,3)} {result} ")
        
    x = x +1
    
percentage_positive = round(a / rep * 100,1)
false_neg = 100 - percentage_positive

print(f"\n ***** Summary*****")
print(f"Average of CONTROL group: {average_input} +/- {round(sd_dev_input,2)}  --- CV= {CV_control}")
print(f"Average of TREATED group: {average_input} +/- {round(expect_sd_dev,2)} --- CV= {CV_treatment}")
print(f"\nNumber of planned repetitions / group: {num_rep}")
print(f"Defined P value = {P_level} (Type I Error)")
print("-"*80)  
print(f"\nThe False Negative rate (Type II error) of the trial is estimated by: {round(false_neg,1)} %")
print(f"The stat power rate is estimated by: {round(100-false_neg,1)} %")
print("The standard Type II error rate is 20% or 80% stat power\n")


# -------------------------------------------------------------------
# create Histogram 

n_hist = 2000000 # num of replicate for histogram

hist_control = np.random.normal(average_input,sd_dev_input,n_hist) # generate control group
hist_treatment = np.random.normal(average_input,expect_sd_dev,n_hist) # generate control group

plt.hist([hist_control,hist_treatment],
         bins=50,
         label=["Control", "Probiotic"],
         )
plt.legend()
plt.title("Theorical Distribution of Populations")
plt.show()

#------------------------------------------------------------------------
"""
Estimation of sample size needed by montecarlo methodology

"""
def n_optimization (n_rep, actual_power, power, treat_effect):

    if treat_effect <= 10:    # consider high impact of treatment on homogenity redution
    
        if actual_power < power*.5:
            n_sample = n_rep + 300
        elif actual_power < power*.8:
            n_sample = n_rep + 100
        elif actual_power < power*.9:
            n_sample = n_rep + 50
        elif actual_power < power*.96:
            n_sample = n_rep + 20
        elif actual_power < power*.98:
            n_sample = n_rep + 5
        else:
            n_sample = n_rep + 1
    
    elif treat_effect <= 15: # consider medium impact of treatment on homogenity redution
        if actual_power < power*.5:
            n_sample = n_rep + 100
        elif actual_power < power*.8:
            n_sample = n_rep + 50
        elif actual_power < power*.9:
            n_sample = n_rep + 20
        elif actual_power < power*.96:
            n_sample = n_rep + 5
        elif actual_power < power*.98:
            n_sample = n_rep + 2
        else:
            n_sample = n_rep + 1    
    elif treat_effect <= 30: # consider low impact of treatment on homogenity redution
        if actual_power < power*.5:
            n_sample = n_rep + 5
        elif actual_power < power*.8:
            n_sample = n_rep + 2
        elif actual_power < power*.9:
            n_sample = n_rep + 20
        else:
            n_sample = n_rep + 1  
    else:
        n_sample = n_rep + 1
    return n_sample





percentage_positive = 0

# ------------------------------------------------------------

rep = 5  # min sample size 
trials_simul = 3000 # numbers of trials simulations - if big number = slower 

print("\n --------RUNNING SAMPLE SIZE ESTIMATION -----------")
print("\n .... WAIT - it can take some time....\n")

while percentage_positive <= stat_power_level:
    a = 0
    for i in range(trials_simul):

        simul_control = np.random.normal(average_input,sd_dev_input,rep) # generate control group
        simul_treatment = np.random.normal(average_input,expect_sd_dev,rep) # generate control group
    
        stat, p = stats.levene(simul_control,simul_treatment )
    
        if p<= P_level:
            a = a + 1
                 
        
        #print(f"trial nº: {x} --> Stat = {round(stat,1)} , P = {round(p,3)}")
        
    percentage_positive = round(a / trials_simul * 100,1) # calc of stat power
    
    rep= n_optimization (rep, percentage_positive, stat_power_level, input_diff ) # function - algoritm for n optimization
    #print(percentage_positive, rep)
    

print("-"*80)         
print(f"Aprox. sample size / group needed is:  {rep}  --> Estimated stat Power: {percentage_positive}%, P_level <= {P_level} ")        
print(f"Planned Stat Power: {stat_power_level} %")
print("-"*80)         
print("\nFor this simulation Levene's test was applied to compare homogeneity of variances among groups.")
print(f"*** Code version: {version} ***")








# simul_data_control = np.random.normal(control_mean, control_std_desv, n_samples_per_treatment)