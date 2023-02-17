import subprocess
import os
import time
import numpy as np

def run_test(n=8000, duration = 1, timestep = 0.1, brute_force = False, cores = 1):
  command = f"./build/src/barnes-hut -n{n} -d{duration} -t{timestep} { ' --brute_force' if brute_force else ''}"
  env = os.environ.copy()
  env["OMP_NUM_THREADS"] = str(cores)

  # measure execution time of the command
  start = time.time()
  subprocess.run(command, shell=True, env=env)
  end = time.time()

  return end-start

avg = np.zeros((16,10,2))
for i in range(10):
  # create a 16x10x2 numpy array to store the results
  results = np.zeros((16,10,2))
  # vary cores, brute force and n
  for cores in range(1,17):
    for brute_force in range(2):
      for n in range(10):
        print(f"cores: {cores}, brute force: {brute_force}, n: {n*1000+1000}")
        results[cores-1, n, brute_force] = run_test(n=n*1000+1000, brute_force=brute_force, cores=cores)
  avg += results
  # print the results
  print(results)

# save the results to a file
np.save("results.npy", avg)
