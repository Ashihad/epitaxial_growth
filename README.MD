# Master thesis

# Requirements
- C++17 compiler, available in `$PATH`
- CMake, version >=3.22
- Modern Python interpreter (tested on version 3.13.7)

# Compilation steps
## Release build
```
mkdir release/
cd release/
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
```

## Debug build (-O0 -ggdb3)
```
mkdir debug/
cd debug/
cmake .. -DCMAKE_BUILD_TYPE=Debug
cmake --build .
```

## Running program
```
release/app/main
or
debug/app/main
```

### Python setup
Using `venv` package is recommended, here is a simple setup:
```
mkdir .venv
python -m venv .venv
source .venv/bin/activate
``` 
Next, install all required Python packages in `venv`: 
```
pip install -r requirements.txt
```
To disable `venv` run `deactivate` in shell