"""
Description:
------------
This script simulates a stochastic reaction system using a Stochastic Simulation Algorithm (SSA).
The system represented is a model of concentration-time profile of a pharmaceutic drug. Further details
can be found in the PDF located in this directory. 

The script reads initial conditions and reaction parameters from text files and then:
  - Computes reaction rates based on the current state.
  - Determines the time until the next reaction and which reaction occurs.
  - Runs the simulation until termination conditions are met.
  - Extracts infection times from the trajectory and saves/plots the results.
  
Additionally, the code estimates the uncertainty in vaccine efficacy (for three vaccines: Moderna, BioNTech, and AstraZeneca)
using a bootstrap method and plots histograms for the sampled incidence rate parameters.

Example:
  If you have a file 'Input.txt' containing initial parameters and reaction rates in 'V_Bio_kinf.txt', running
  this script will simulate the process, output the trajectory to a file, and display several plots that illustrate
  the infection dynamics and the bootstrap distributions for vaccine efficacies.
"""

import numpy as np
import matplotlib.pyplot as plt

# --------------------------
# A: Load fixed input quantities
# --------------------------
# Read initial state and parameters from file
input_data = np.loadtxt('Input.txt')

seed_num = input_data[0]
num_traj = int(input_data[1])
np.random.seed(seed=int(seed_num))

# --------------------------
# B: Define functions for the SSA
# --------------------------

def ReactionRates(k, X):
    """
    Compute the reaction rates (propensities) based on the current state.
    
    Parameters:
      k : array_like
          Reaction rate constants (e.g., [k_drop, k_inf]).
      X : array_like
          Current state vector. For example, X[0] may be the number of participants.
    
    Returns:
      R : ndarray
          A 2x1 array with computed reaction rates.
    
    Example:
      If X[0] = 100 and k = [0.05, 0.01],
      then R[0] = 0.05*100 = 5 and R[1] = 0.01*100 = 1.
    """
    R = np.zeros((2, 1))
    R[0] = k[0] * X[0]  # Reaction rate for the first reaction
    R[1] = k[1] * X[0]  # Reaction rate for the second reaction
    return R

def Time_to_next_reaction(total_rate):
    """
    Calculate the time until the next reaction occurs.
    
    Uses the formula: tau = (1/total_rate) * log(1/r)
    where r is a random number between 0 and 1.
    
    Parameters:
      total_rate : float
          Sum of all reaction propensities.
    
    Returns:
      tau : float
          Time increment until the next reaction.
          
    Example:
      If total_rate = 10 and a random r=0.5, then tau = (1/10)*log(2).
    """
    r = np.random.rand()
    # Avoid division by zero in logarithm by ensuring r is never 0
    while r == 0:
        r = np.random.rand()
    return (1.0 / total_rate) * np.log(1.0 / r)

def FindReactionIndex(propensities):
    """
    Determine which reaction will occur next based on propensities.
    
    This function uses a random number to select the reaction index.
    
    Parameters:
      propensities : array_like
          An array of reaction propensities.
    
    Returns:
      index : int
          The index of the chosen reaction.
          
    Example:
      If propensities = [5, 1] and a random number r is generated,
      the function will choose index 0 if the cumulative sum condition is met.
    """
    r = np.random.rand()
    # Ensure r is not 0 to avoid issues in scaling
    while r == 0:
        r = np.random.rand()
    total = np.sum(propensities)
    cumulative = np.cumsum(propensities)
    # Count how many elements in cumulative are less than r * total
    index = int(np.sum(cumulative < r * total))
    return index

def SSA_to0(Stochiometry, X0, Reaction_vec):
    """
    Perform the Stochastic Simulation Algorithm (SSA) until termination.
    
    The simulation stops when the first component of the state X is <= 0 or if the
    total propensity becomes zero.
    
    Parameters:
      Stochiometry : ndarray
          Stoichiometric matrix representing state changes for each reaction.
      X0 : array_like
          Initial state vector.
      Reaction_vec : array_like
          Reaction rate constants.
    
    Returns:
      X_store : ndarray
          Array of state vectors over time.
      T_store : ndarray
          Array of corresponding time points.
          
    Example:
      Given a starting population X0 = [100, 0] and a stoichiometric matrix S,
      the function simulates how the state evolves until the process stops.
    """
    # Storage for simulation history
    T_store = []
    X_store = []
    
    # Initialize simulation variables
    t = 0.0
    x = X0.copy()  # Use copy to avoid modifying the input array
    X_store.append(x)
    T_store.append(t)
    
    while True:
        # Compute propensities based on the current state
        propensities = ReactionRates(Reaction_vec, x)
        total_rate = np.sum(propensities)  # Renamed variable to avoid conflict with built-in "sum"
        
        # If total_rate is nonzero, compute time to next reaction
        if total_rate != 0:
            tau = Time_to_next_reaction(total_rate)
        
        # Termination condition: if there are no more reactants or no reactions can occur
        if (x[0] <= 0) or (total_rate == 0):
            return np.array(X_store), np.array(T_store)
        else:
            t += tau  # Advance time
            # Determine which reaction occurs next
            j = FindReactionIndex(propensities)
            # Update state based on the stoichiometric change for the reaction j
            x = x + Stochiometry[:, j]
            # Store the updated state and time
            X_store.append(x)
            T_store.append(t)

# --------------------------
# C: Set up simulation parameters and perform a single simulation run
# --------------------------

n_sim = 100
n_cases = 8

# Load reaction rate parameters for BioNTech from file
k_inf_full = np.loadtxt('V_Bio_kinf.txt')

# Initial conditions: number of participants and initial infected count
participants = 17411
av_sur_time = 0.127160990178623

# Use the first k_inf value for the initial simulation
k_inf_1 = k_inf_full[0]

# Compute the "drop" rate from the average survival time (example: drop = 1/av_sur_time - k_inf)
k_drop = []
for x in k_inf_full:
    k_drop.append((1.0 / av_sur_time) - x)

# Set initial state vector: [number of participants, initial infected count]
X = np.array([participants, 0])
# Compute drop rate for first parameter
k_drop_1 = (1.0 / av_sur_time) - k_inf_1
# Combine drop and infection rates into an array
k0 = np.array([k_drop_1, k_inf_1])

# Stoichiometric matrix: columns represent state changes for each reaction
S = np.array([[-1, -1],
              [ 0,  1]])

# Run a single SSA simulation (from Homework Task 1a)
states, times = SSA_to0(S, X, k0)

# Extract times when the number of infected increases
infection_indices = np.where(np.diff(states[:, 1]) > 0)[0] + 1
infection_times = times[infection_indices]
infection_times = np.insert(infection_times, 0, 0)  # Include initial time

# Extract corresponding infected state values
infection_states = states[infection_indices, 1]
infection_states = np.insert(infection_states, 0, 0)  # Include initial state

# Save trajectory (time, state) to file
Output = np.concatenate((np.array(infection_times, ndmin=2),
                         np.array(infection_states, ndmin=2)), axis=0)
filename = 'V_BionTech_Task2Infected.txt'
np.savetxt(filename, Output, delimiter=',', fmt='%1.3f')
print(f"Saved trajectory 1 to {filename}")

# Plot infection dynamics using a step plot
plt.step(infection_times, infection_states)
plt.xlabel('Time')
plt.ylabel('Amount of Infected')
plt.title('Infection Dynamics Over Time')
plt.show()

# --------------------------
# D: Run multiple simulation trajectories to gather statistics
# --------------------------
final_infected = []
final_times = []
num_simulations = 100

# Loop over a set number of simulations using different reaction rate parameters
for i in range(num_simulations):
    # For each simulation, select the i-th parameters from the loaded arrays
    k_i = k_inf_full[i]
    k_do = k_drop[i]
    # Combine the parameters into an array (note: here we use the new k values)
    k = np.array([k_do, k_i])
    
    # Run the simulation using the current parameters (use k instead of k0)
    states, times = SSA_to0(S, X, k)
    
    # Extract infection times and states as before
    infection_indices = np.where(np.diff(states[:, 1]) > 0)[0] + 1
    infection_times = times[infection_indices]
    infection_times = np.insert(infection_times, 0, 0)
    
    infection_states = states[infection_indices, 1]
    infection_states = np.insert(infection_states, 0, 0)
    
    # Save the final infected count and corresponding time from this simulation
    final_infected.append(infection_states[-1])
    final_times.append(infection_times[-1])

# Print statistics of the final infected counts across simulations
print("Mean infected: ", np.mean(final_infected))
print("Standard deviation: ", np.std(final_infected))

# Plot histogram of final infected counts
plt.hist(final_infected, density=False, alpha=0.75, edgecolor='black', linewidth=1.2)
plt.ylabel('Amount of Runs')
plt.xlabel('Amount of Infected')
plt.axvline(np.mean(final_infected), color='red', label="Mean Infected")
plt.axvline(n_cases, color='black', label="Case Threshold")
plt.title('Distribution of Final Infected Counts')
plt.legend()
plt.show()

# --------------------------
# E: Bootstrap estimates for vaccine efficacy (Task 2c)
# --------------------------
# The efficacy is defined as: efficacy = 1 - (k_inf_V / k_inf_P)
# where k_inf_V and k_inf_P are the incidence rates for the vaccine and placebo groups, respectively.

# --- Moderna ---
k_inf_V_MOD = np.loadtxt('V_Moderna_kinf.txt')
k_inf_P_MOD = np.loadtxt('P_Moderna_kinf.txt')

efficacy = 1 - k_inf_V_MOD[0] / k_inf_P_MOD[0]
mean_efficacy = 1 - np.mean(k_inf_V_MOD) / np.mean(k_inf_P_MOD)

print("Moderna")
print("Efficacy: ", efficacy)
print("Mean Efficacy: ", mean_efficacy)

n_bootstrap = 10000
n_samples = 1000
quotients = np.zeros(n_bootstrap)

for i in range(n_bootstrap):
    indices = np.random.randint(0, n_samples, size=n_samples)
    num_sample = k_inf_V_MOD[indices]
    den_sample = k_inf_P_MOD[indices]
    quotients[i] = 1 - np.mean(num_sample) / np.mean(den_sample)

mean_quotient = np.mean(quotients)
std_quotient = np.std(quotients)
confidence_interval = np.percentile(quotients, [2.5, 97.5])

print('Moderna bootstrap estimates:')
print('mean', mean_quotient)
print('std', std_quotient)
print('ci_lower', confidence_interval[0])
print('ci_upper', confidence_interval[1])

plt.hist(quotients, density=False, alpha=0.75, edgecolor='black', linewidth=1.2)
plt.ylabel('Runs')
plt.xlabel('Efficacy')
plt.title("Vaccine Efficacy Simulation: Moderna")
plt.axvline(mean_quotient, color='black', label="Mean")
plt.axvline(confidence_interval[0], color='red', label="Lower CI")
plt.axvline(confidence_interval[1], color='red', label="Upper CI")
plt.legend()
plt.show()

# --- BioNTech ---
k_inf_V_BIO = np.loadtxt('V_Bio_kinf.txt')
k_inf_P_BIO = np.loadtxt('P_Bio_kinf.txt')

efficacy = 1 - k_inf_V_BIO[0] / k_inf_P_BIO[0]
mean_efficacy = 1 - np.mean(k_inf_V_BIO) / np.mean(k_inf_P_BIO)

print("BioNTech")
print("Efficacy: ", efficacy)
print("Mean Efficacy: ", mean_efficacy)

quotients = np.zeros(n_bootstrap)
for i in range(n_bootstrap):
    indices = np.random.randint(0, n_samples, size=n_samples)
    num_sample = k_inf_V_BIO[indices]
    den_sample = k_inf_P_BIO[indices]
    quotients[i] = 1 - np.mean(num_sample) / np.mean(den_sample)

mean_quotient = np.mean(quotients)
std_quotient = np.std(quotients)
confidence_interval = np.percentile(quotients, [2.5, 97.5])

print('BioNTech bootstrap estimates:')
print('mean', mean_quotient)
print('std', std_quotient)
print('ci_lower', confidence_interval[0])
print('ci_upper', confidence_interval[1])

plt.hist(quotients, density=False, alpha=0.75, edgecolor='black', linewidth=1.2)
plt.ylabel('Runs')
plt.xlabel('Efficacy')
plt.title("Vaccine Efficacy Simulation: BioNTech")
plt.axvline(mean_quotient, color='black', label="Mean")
plt.axvline(confidence_interval[0], color='red', label="Lower CI")
plt.axvline(confidence_interval[1], color='red', label="Upper CI")
plt.legend()
plt.show()

# --- AstraZeneca ---
k_inf_V_A = np.loadtxt('V_Astra_kinf.txt')
k_inf_P_A = np.loadtxt('P_Astra_kinf.txt')

efficacy = 1 - k_inf_V_A[0] / k_inf_P_A[0]
mean_efficacy = 1 - np.mean(k_inf_V_A) / np.mean(k_inf_P_A)

print("AstraZeneca")
print("Efficacy: ", efficacy)
print("Mean Efficacy: ", mean_efficacy)

quotients = np.zeros(n_bootstrap)
for i in range(n_bootstrap):
    indices = np.random.randint(0, n_samples, size=n_samples)
    num_sample = k_inf_V_A[indices]
    den_sample = k_inf_P_A[indices]
    quotients[i] = 1 - np.mean(num_sample) / np.mean(den_sample)

mean_quotient = np.mean(quotients)
std_quotient = np.std(quotients)
confidence_interval = np.percentile(quotients, [2.5, 97.5])

print('AstraZeneca bootstrap estimates:')
print('mean', mean_quotient)
print('std', std_quotient)
print('ci_lower', confidence_interval[0])
print('ci_upper', confidence_interval[1])

plt.hist(quotients, density=False, alpha=0.75, edgecolor='black', linewidth=1.2)
plt.ylabel('Runs')
plt.xlabel('Efficacy')
plt.title("Vaccine Efficacy Simulation: AstraZeneca")
plt.axvline(mean_quotient, color='black', label="Mean")
plt.axvline(confidence_interval[0], color='red', label="Lower CI")
plt.axvline(confidence_interval[1], color='red', label="Upper CI")
plt.legend()
plt.show()

# --------------------------
# F: Plot histograms for incidence rate parameters (Task 2a)
# --------------------------

# Example plots for the incidence rate parameters for each vaccine and placebo group

# BioNTech Placebo
plt.hist(k_inf_P_BIO, density=True, alpha=0.75, edgecolor='black', linewidth=1.2)
plt.ylabel('Samples')
plt.xlabel('k_inf')
plt.title("Inversely Sampled Incidence Rates: BioNTech Placebo")
plt.show()

# BioNTech Vaccine
plt.hist(k_inf_V_BIO, density=True, alpha=0.75, edgecolor='black', linewidth=1.2)
plt.ylabel('Samples')
plt.xlabel('k_inf')
plt.title("Inversely Sampled Incidence Rates: BioNTech Vaccine")
plt.show()

# Moderna Placebo
plt.hist(k_inf_P_MOD, density=True, alpha=0.75, edgecolor='black', linewidth=1.2)
plt.ylabel('Samples')
plt.xlabel('k_inf')
plt.title("Inversely Sampled Incidence Rates: Moderna Placebo")
plt.show()

# Moderna Vaccine
plt.hist(k_inf_V_MOD, density=True, alpha=0.75, edgecolor='black', linewidth=1.2)
plt.ylabel('Samples')
plt.xlabel('k_inf')
plt.title("Inversely Sampled Incidence Rates: Moderna Vaccine")
plt.show()

# AstraZeneca Placebo
plt.hist(k_inf_P_A, density=True, alpha=0.75, edgecolor='black', linewidth=1.2)
plt.ylabel('Samples')
plt.xlabel('k_inf')
plt.title("Inversely Sampled Incidence Rates: AstraZeneca Placebo")
plt.show()

# AstraZeneca Vaccine
plt.hist(k_inf_V_A, density=True, alpha=0.75, edgecolor='black', linewidth=1.2)
plt.ylabel('Samples')
plt.xlabel('k_inf')
plt.title("Inversely Sampled Incidence Rates: AstraZeneca Vaccine")
plt.show()