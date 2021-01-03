# What is TumGenomSim?

TumGenomSim is a simulator that generates multi-regional tumor sequencing data. It has been tested on Ubuntu and RedHat linux, but should run without modifications on all other Unix-based platforms, including OS X.

## Dependencies

* [Boost](https://www.boost.org/users/history/version_1_68_0.html) (1.68)
* [Cmake](https://cmake.org/) (>=3.10.2)

## Install

Given the required dependencies, the installation should only require a few steps:

   git clone https://github.com/hdetering/tumgenomsim.git
   mkdir tumgenomsim/build && cd tumgenomsim/build
   cmake ..
   make

The executable can be found at `tumgenomsim/build/src/tgs`.

## Command line
Since most of the user parameters are provided in a config file, the command line arguments to run TumGenomSim are kept to a minimum. Where arguments overlap with the config file, the command line arguments override those in the file.

    $ ./tgs -h
    
    TumGenomSim 1.2.0
    
    Available options:
      -v [ --version ]          print version string
      -h [ --help ]             print help message
      -c [ --config ] arg       config file
      -t [ --tree ] arg         file containing user defined clone tree (Newick 
                                format)
      -o [ --out-dir ] arg      output directory
      -p [ --threads ] arg (=1) number of parallel threads
      -s [ --seed ] arg         random seed
    
The command line parameters control the following aspects of a TumGenomSim run:
`-v/--version`: just print the program version and exit
`-h/--help`: 	  print usage info
`-c/--config`: 	path to config file
`-t/--tree`: 	  clone tree to use in simulation
`-o/--output`: 	path to output directory (will be created if not exists)
`-p/--threads`:	how many threads to use in parallel execution
`-s/--seed`:	  seed to initialize random number generator

For a complete list of user parameters see the [annotated example config](https://github.com/hdetering/tumgenomsim/blob/master/config.yml).
