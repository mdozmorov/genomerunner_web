from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import glob

setup(
        name="grtk",
        version="0.0.1-HEAD",
        author="Cory Giles",
        author_email="mail@corygil.es",
        description="Utilities for manipulating genomic region data.",

        packages=["grtk"],
        scripts=glob.glob("script/*"),
        cmdclass= {"build_ext": build_ext},
        ext_modules = [
            Extension(
                "grtk.bbi",
                include_dirs=["src"],
                sources=["grtk/bbi.pyx", "src/bbi.cpp", "src/mm.cpp"],
                extra_compile_args=["-std=c++0x"],
                libraries=["z"],
                language="c++")
            ]
)
