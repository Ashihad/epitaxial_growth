all: build

install:
	sudo apt-get install gcc g++ cmake make doxygen git llvm pkg-config curl zip unzip tar python3-dev clang-format clang-tidy

prepare:
	rm -rf build
	mkdir build

build: prepare
	@cd build 
	make 
	@cd -

clang-tidy:
	cd build && $(MAKE) objects_clangtidy
	cd build && $(MAKE) epitaxal_mc_clangtidy