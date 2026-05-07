# Building and Releasing

## Building from source
First, create a Python `venv`, like so:
```sh
cd python_bindings
python3 -m venv .env
source .env/bin/activate
```

We use `maturin` and `pyo3` to develop, build and publish. Make sure that the `maturin` executable is available, for example by running 
```sh
pip install maturin
```
or by using the package manager of your choice, e.g. `dnf install maturin` or `apt install python3-maturin`.

For developing, run
```sh
maturin develop
```
which will install the crate as a python module in the current virtual environment. You can then e.g. use it from a REPL, like so:
```sh
python3
>>> import tsalign
>>> alignment = tsalign.align("ACGT", "AGT")
>>> alignment.cigar()
'1M1D2M'
>>> ...
```

## Release
Release by simply adjusting the version in `Cargo.toml`, create and push a tag that is called `py_lib_tsalign-v<version>`. A GitHub action should™ publish to the [PyPI repository](https://pypi.org/project/tsalign/).
