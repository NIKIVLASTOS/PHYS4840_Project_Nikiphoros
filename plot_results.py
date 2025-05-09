#===============                                            
# File: plot_results.py                                                # Python to visualize 1D Schrodinger results
#===================
import numpy as np                                                     
import matplotlib.pyplot as plt                                        
import os                                                              

# ##### Read input.txt#####
def read_clean_input(path):                                            # This Function to rea and clean input.txt
    with open(path, "r") as f:                                         # Open file for reading
        return [line.strip() for line in f if line.strip() and not line.strip().startswith("!")]  # return non-empty, non-comment lines

lines = read_clean_input("input.txt")                                  # read  lines
xmin, xmax, Nx_check = map(float, lines[0].split())                    # parse first line: grid limits and Nx
potential_type = int(lines[1])                                         # parse second line: potential type
nstates = int(lines[2])                                                # parse third line: number of eigenstates
V0 = float(lines[3]) if potential_type == 3 else 0.0                   # parse V0 if potential type is 3 (finite square well)

# #### Load data ##
x = np.loadtxt("xgrid.txt")                                            # load spatial grid
eigenvalues = np.loadtxt("eigenvalues.txt")                            # load computed eigenvalues
wavefunctions = np.loadtxt("wavefunctions.txt", ndmin=2)               # load wavefunctions (2D array forced)

# ## Ensure wavefunctions are shape (nstates, Nx) ##
if wavefunctions.shape[0] != nstates:                                  # if shape is (Nx, nstates), transpose
    wavefunctions = wavefunctions.T                                    # transpose to get (nstates, Nx)

neigs, Nx = wavefunctions.shape                                        # get dimensions of wavefunction array
assert len(x) == Nx                                                    # make sure x-grid length matches Nx

# #### Define potential ######
def V(x, potential_type, xmin, xmax, V0):                              # Function to evaluate potential at x
    if potential_type == 1:                                            # Harmonic oscillator
        return 0.5 * 100 * x**2                                       # Scaled harmonic oscillator (match scaling to Fortran laid out scaling facotr which is 100 right now)
    elif potential_type == 2:                              # Infinite square well from xmin+1 to xmax-1
        wall_left = xmin + 1.0
        wall_right = xmax - 1.0
        V_inf_plot = 100                          # I use 100 instead of 1e12 to make the walls visible in the plots, This can be change by User (Dr. Joyce or grader)
        return np.where((x >= wall_left) & (x <= wall_right), 0.0, V_inf_plot)
    elif potential_type == 3:                                          # Finite square well        
        wall_left = xmin + 2.0
        wall_right = xmax - 2.0
        return np.where((x >= wall_left) & (x <= wall_right), 0.0, V0)
    elif potential_type == 4:
        return np.where(x < 0, 0, 5.0)                                 # Step up to 5.0 at x = 0
    elif potential_type == 5:                                        # step then infinite wall
        return np.piecewise(x,
            [x < 0, (x >= 0) & (x < 3), x >= 3],
            [0.0, 5.0, 1e12])
    elif potential_type == 6:                                        # Barrier over finite domain
        return np.where((x >= 0) & (x <= 5), 5.0, 0.0)
    else:
        raise ValueError("Invalid potential type.")                    # rrror for invalid potential code (e.g. if you enter >6)

Vx = V(x, potential_type, xmin, xmax, V0)                              # evaluate full potential array

# ### Create output directory #
os.makedirs("plots", exist_ok=True)                                    # make directory to store plots (no error if exists)


potential_titles = {                                                   # names for the TITles for each plot
    1: "Harmonic Potential",
    2: "Infinite Square Well",
    3: "Finite Square Well",
    4: "Step Potential",
    5: "Stepped Trap Potential",
    6: "Step Barrier"
}
pot_name = potential_titles.get(potential_type, "Unknown Potential")

# # Diagnostic plot: raw wavefunctions over grid ###
plt.figure(figsize=(10, 6))                                            # Create figure for raw wavefunctions
for n in range(nstates):                                               # Loop over states
    plt.plot(x, wavefunctions[n], label=f"$\\psi_{n}(x)$")             # plot each raw wavefunction
plt.title("Raw wavefunctions (no energy offset or scaling)")           # Title
plt.xlabel("x")                                                        # X-axis label
plt.grid(True)                                                         # Add grid
plt.legend()                                                           # Add legend
plt.tight_layout()                                                     # Tight layout
plt.savefig("plots/raw_wavefunctions.png")                             # save raw wavefunction plot

#####  Plot settings #####
wf_scale = 5.0                                                        # Vertical scaling for visibility
ymin = -5                                                              # Lower y-limit
ymax = max(eigenvalues[:nstates]) + 25                                 # Upper y-limit (based on max eigenvalue)
xlim_range = (x[0] - 2, x[-1]+ 2)                                             # X-axis limits from grid

if potential_type in [2, 3,4,5,6]:                                      # Change energy scale so you can see different wave functions better, not needed for harmonic oscillator
    energy_scale = 15.0
else:
    energy_scale = 1.0



##### Plot wavefunctions  ####
plt.figure(figsize=(10, 6))                                            # create figure
plt.plot(x, Vx, 'r--', label="V(x)")                                   # Plot potential
for n in range(nstates):                                               # Plot wavefunctions shifted by energy
    norm_factor = np.max(np.abs(wavefunctions[n]))
    plt.plot(x, energy_scale * eigenvalues[n] + wf_scale * wavefunctions[n] / norm_factor, label=rf"$\psi_{{{n}}}(x)$")
for n in range(nstates):                                               # Draw horizontal lines at each eigenvalue
    plt.axhline(eigenvalues[n], color='gray', linestyle='-', linewidth=0.5)
plt.axhline(0, color='black', linestyle='--', linewidth=0.5)           # Line at E = 0


plt.title(rf"Wavefunctions $\psi_n(x)$ in {pot_name}")
plt.xlabel("x")                                                        # Xaxis label
plt.ylabel("Amplitude (shifted to $E_n$)")                             # Yaxis label
plt.ylim(ymin, ymax)                                                   # Set yaxis limits
plt.xlim(*xlim_range)                                                  # set xaxis limits
plt.grid(True)                                                         # Add grid
plt.legend()                                                           # Add legend
plt.tight_layout()                                                     # tight layout
plt.savefig("plots/wavefunctions.png")                                 # save wavefunctions plot


#### Plot probability densities ##
plt.figure(figsize=(10, 6))                                            # figure for probability densities
plt.plot(x, Vx, 'r--', label="V(x)")                                   # Plot potential
for n in range(nstates):                                               # Plot shifted probability densities
    norm_factor = np.max(wavefunctions[n]**2)
    plt.plot(x, energy_scale * eigenvalues[n] + wf_scale * wavefunctions[n]**2 / norm_factor, label=rf"$|\psi_{{{n}}}(x)|^2$")
for n in range(nstates):                                               # Energy level lines
    plt.axhline(eigenvalues[n], color='gray', linestyle='-', linewidth=0.5)
plt.axhline(0, color='black', linestyle='--', linewidth=0.5)           # Line at E = 0

plt.title(rf"Probability Densities $|\psi_n(x)|^2$ in {pot_name}")  # Title
plt.xlabel("x")                                                        # Xaxis label
plt.ylabel("Probability (shifted to $E_n$)")                           # Yaxis label
plt.ylim(ymin, ymax)                                                   # Yaxis range
plt.xlim(*xlim_range)                                                  # Xaxis range
plt.grid(True)                                                         # Grid
plt.legend()                                                           # Legend
plt.tight_layout()                                                     # Layout
plt.savefig("plots/probability_densities.png")                         # Save plot

print("Plots saved in 'plots/' directory.")                            # Print completion message

