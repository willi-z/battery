# Battery Simulation with Fenics

# With Docker
## 0. First Time

````
docker run -ti -v ${PWD}:/home/fenics/shared --name battery dolfinx/lab:stable
````

## 1. Starting
````
docker start battery
````

## 2. Stopping
````
docker stop new-project
````


## Reference
[Fenics Documentation](https://docs.fenicsproject.org/)
[Fenics Python Documentation](https://docs.fenicsproject.org/dolfinx/v0.5.1/python/)
[Fenics Docker](https://fenics.readthedocs.io/projects/containers/en/latest/index.html)


# Bare Python

## Prerequirements

```sh
sudo apt install -y python3-pip
sudo apt install -y build-essential libssl-dev libffi-dev python3-dev
sudo apt install -y python3-venv
```

## Setup

```sh
python3 -m venv env
```

##  Workflow

```sh
source env/bin/activate
```