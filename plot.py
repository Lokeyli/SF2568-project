## Read the given txt files from command line arguments and plot the data
## Usage: python plot.py <file1> <file2> <file3> ...
## Example: python plot.py data1.txt data2.txt data3.txt

import sys
import matplotlib.pyplot as plt
import numpy as np

# Read the given txt files from command line arguments
for i in range(1, len(sys.argv)):
    with open(sys.argv[i], "r") as f:
        data = f.readlines()
        x = [float(line.split()[0]) for line in data]
        y = [float(line.split()[2]) for line in data]
        plt.scatter(x, y, label=sys.argv[i])

plt.xlabel("Number of Core")
plt.ylabel("Time Taken/s")
plt.title("Time Taken against Number of Core")
plt.legend()
plt.savefig("docs/figures/time-taken.png")

# Plot another graph, the speedup
# Using P = 1 as the baseline

# Clean the plot
plt.clf()


speedups = []
X = []
filename_postfix = sys.argv[-1][-1]
for i in range(1, len(sys.argv)):
    with open(sys.argv[i], "r") as f:
        data = f.readlines()
        x = [float(line.split()[0]) for line in data]
        y = [float(line.split()[2]) for line in data]
        # Sort y and x base on x, in ascending order
        x, y = zip(*sorted(zip(x, y)))
        # Calculate speedup
        speedup = [y[0] / time for time in y]
        plt.scatter(x, speedup, label=sys.argv[i])
        speedups.append(speedup)
        X = x
# Calculate the average speedup
average_speedup = np.mean(speedups, axis=0)
plt.plot(X, average_speedup, label="Average")
# Calculate the average efficiency, e = speedup / P
average_efficiency = average_speedup / X

plt.xlabel("Number of Core")
plt.ylabel("Speedup")
plt.title("Speedup against Number of Core")
plt.legend()
plt.savefig(f"docs/figures/speedup-compare-blocked-nonblocked-{filename_postfix}.png")

# Plot another graph, the efficiency
# Using P = 1 as the baseline

# Clean the plot
plt.clf()

plt.plot(X, average_efficiency, label="Average")
plt.xlabel("Number of Core")
plt.ylabel("Efficiency")
plt.title("Efficiency against Number of Core")
plt.legend()
plt.savefig("docs/figures/efficiency.png")
