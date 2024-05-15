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

for i in range(1, len(sys.argv)):
    with open(sys.argv[i], "r") as f:
        data = f.readlines()
        x = [float(line.split()[0]) for line in data]
        y = [1 / float(line.split()[2]) for line in data]
        plt.scatter(x, y, label=sys.argv[i])

plt.xlabel("Number of Core")
plt.ylabel("Speedup")
plt.title("Speedup against Number of Core")
plt.legend()
plt.savefig("docs/figures/speedup.png")
