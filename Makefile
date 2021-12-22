.PHONY: clean format doc

help:
	@echo "Please use \`make <target>\` where <target> is one of"
	@echo "  clean"
	@echo "  format-cpp"
	@echo "  doc"

clean:
	rm -rf ./build/* ./UnionFindPy.egg-info 
	$(MAKE) -C ./docs clean
	$(MAKE) -C ./union_find/cpp clean

format-cpp:
	@./build_utils/format_cpp.sh

doc:
	$(MAKE) -C ./docs html
