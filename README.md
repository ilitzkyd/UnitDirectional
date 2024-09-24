# UnitDirectional Project

This project involves simulations of uniaxial stress fibers in the 
y-direction, using **FEniCS** for solving boundary value problems in a 
3D mesh. The project is structured to define and solve fiber contraction 
problems within different material domains, including gel, cytoplasm, 
and nucleus.

## Project Structure

UnitDirectional/

├── 06092019_G1/   

├── dolfin/  

│   └── python/

│       └── setup.py   

├── UnitDir.py  

└── README.md          



## Dependencies

The project uses the following key dependencies:

- **FEniCS** (for solving PDEs)
- **meshio** (for mesh file handling)
- **mshr** (for mesh generation)

Make sure these are installed before running the project. 
Additionally, you can use the `setup.py` file in the `dolfin/python/` 
directory to install the required packages for **Texar**.

## Running the Project

1. Ensure all dependencies are installed by running:
```bash
pip install -r requirements.txt
```
2. Run the sumulation:
```bash
python UnitDir.py
```

3. Results will be saved as `.pvd` files in the `result/` directory.

## Mesh Files

The folder `06092019_G1/` contains pre-processed mesh files in `.xdmf`, 
`.h5`, and `.msh` formats. These are used for setting up the simulation 
in **UnitDir.py**.

## Texar Installation

If you're working with Texar, ensure Python 3.6+ is installed and run 
the `setup.py` script in the `dolfin/python/` folder.


## Contributions

Feel free to contribute by opening an issue or creating a pull request. 
Make sure to follow the guidelines provided in the [Texar repository](https://github.com/asyml/texar).


