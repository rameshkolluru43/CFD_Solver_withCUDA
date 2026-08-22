from setuptools import setup, Extension
import pybind11

ext_modules = [
    Extension(
        "CFD_Solver",
        ["Euler_Python_Wrapper.cpp"],
        include_dirs=[
            pybind11.get_include(),  # Ensure pybind11 is found
            "../",  # Base path
            "../Grid_Readers",
            "../Basic_Files",
            "../Fluxes",
            "../Numerical_Methods",
            "../Boundary_Conditions",
            "../Test_Cases",
            "../Euler_Solver",
            "/opt/homebrew/include",  # VTK & system includes
            "/opt/homebrew/Cellar/jsoncpp/1.9.6/include",  # jsoncpp headers
            "/opt/homebrew/include/nlohmann",  # Ensure JSON headers are included
        ],
        library_dirs=[
            "/opt/homebrew/lib",  # Ensure system libraries are found
            "/opt/homebrew/Cellar/jsoncpp/1.9.6/lib",  # Ensure jsoncpp is found
        ],
        libraries=[
            "jsoncpp",  # Link jsoncpp properly
            "vtkCommonCore-9.4",
            "vtkCommonDataModel-9.4",
            "vtkIOXML-9.4",
            "vtkCommonExecutionModel-9.4",
            "vtkFiltersCore-9.4",
            "vtkCommonMath-9.4",
            "vtkCommonMisc-9.4",
            "vtkCommonSystem-9.4",
            "vtkCommonTransforms-9.4",
            "vtkCommonColor-9.4",
            "vtkCommonComputationalGeometry-9.4",
            "vtksys-9.4",
        ],
        language="c++",
        extra_compile_args=["-std=c++17"],
    )
]

setup(
    name="CFD_Solver",
    ext_modules=ext_modules,
)
