#	Lattice Anonymous Credentials
This repository contains the code by Sven Argo, Corentin Jeudy, Georg Land for the paper **Practical Post-Quantum Signatures for Privacy** by Sven Argo, Tim Güneysu, Corentin Jeudy, Georg Land, Adeline Roux-Langlois, Olivier Sanders.

The repository contains
- the parameter selection/estimation scripts, and
- the implementation

## Parameter Scripts
In the `scripts` subdirectory, we provide parameter estimation scripts which depends on the lattice-estimator.
The lattice-estimator is included as submodule, so please do a recursive pull.

## Code Running Instructions
The code can be built and run with docker, which is explained after the subsequent explanation of building from scratch without docker.

Requirements:
- FLINT, which depends on GMP and MPFR
- AES-NI instructions
- cmake

FLINT is included as submodule to this repository, so please do a recursive pull.
Once FLINT is installed, the following commands allow building the tests and benchmarks:
```shell
cd code
mkdir build
cd build
cmake ..
make -j
```

Then, running `./test` or `./bench` runs the tests or the benchmarks, respectively.

## Building with Docker
For information on how to install docker on your system visit the [docker docs](https://docs.docker.com/).
Depending on your setup, you may need to prefix the following commands with `sudo`.

### 1. Build the docker image
Navigate to the top-level directory which contains this README file. Then run (notice the trailing dot!)
```
$ docker build -t lac-docker .
```
This takes some time as docker needs to pull the base-image for Ubuntu and then builds all dependencies
(i.e. GMP, MPFR and FLINT) as well as the test and benchmark executables from scratch.
In particular the `make check` instructions are *very* time consuming.

The libraries are built with standard options and fine-tuning is possible by adjusting the corresponding
`make` and `configure` commands in the Dockerfile. However, installation paths *must not* be changed.

Note that this step needs to be repeated for any changes in the code base or Dockerfile but not
for subsequent runs of the code.

### 2. Run the docker image
Issue the following command to start the docker container, where `XXX` is either `./test` or `./bench` to select
the tests or benchmarks respectively.
```
$ docker run -t lac-docker XXX
```
You should now see the output of either the test or the benchmarking executables.

### 3. Stopping the running docker image
To get a list of running docker containers use `$ docker ps`. Then you can use `$ docker stop <container_id>` to stop
the execution. For more information see the [docker docs](https://docs.docker.com/engine/reference/builder/).

### 4. Cleaning up
To remove all stopped docker containers use `$docker container prune`.
To remove docker images use `$docker images` to get a list of all images and then `$docker rm <image_id>` to remove
the image in question.
