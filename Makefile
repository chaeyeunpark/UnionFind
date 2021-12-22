
clean:
	rm -rf ./build/* ./UnionFindPy.egg-info 
	$(MAKE) -C ./docs clean
	$(MAKE) -C ./union_find/cpp clean

format:
	@./build_utils/format_cpp.sh

doc:
	$(MAKE) -C ./docs html
