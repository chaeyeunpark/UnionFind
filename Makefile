

clean:
	rm -rf ./build ./externals/* ./UnionFindPy.egg-info

format:
	@./build_utils/format_cpp.sh
