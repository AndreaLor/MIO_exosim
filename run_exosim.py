from exosim_n.run_files import run_exosim
import sys
import numpy as np

def go():
    input_file = sys.argv[1]
    run_exosim.run(input_file)

if __name__ == '__main__':
    go()
