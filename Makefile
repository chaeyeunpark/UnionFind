.PHONY: clean format doc

help:
	@echo "Please use \`make <target>\` where <target> is one of"
	@echo "  clean"
	@echo "  format-cpp"
	@echo "  doc"

clean:
	rm -rf ./build/* ./UnionFindPy.egg-info UnionFindPy/_union_find_py.cpython-*
	$(MAKE) -C ./docs clean
	$(MAKE) -C ./UnionFindPy/cpp clean

format-cpp:
	$(MAKE) -C ./UnionFindPy/cpp format

doc:
	$(MAKE) -C ./docs html
