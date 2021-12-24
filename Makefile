.PHONY: clean format doc

help:
	@echo "Please use \`make <target>\` where <target> is one of"
	@echo "  clean"
	@echo "  format-cpp"
	@echo "  doc"

clean:
	rm -rf ./build/* ./UnionFindPy.egg-info union_find/_union_find_py.cpython-*
	$(MAKE) -C ./docs clean
	$(MAKE) -C ./union_find/cpp clean

format-cpp:
	$(MAKE) -C ./UnionFindPy/cpp format

doc:
	$(MAKE) -C ./docs html
