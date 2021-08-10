project:=bandersnatch
all_projects:=$(project)

all: example

example: 
	cargo build --release
	cp ../target/release/libbanderpy.so python/libbanderpy.so
	python3 python/example.py

clean: 
	cargo clean
	rm python/*.so
