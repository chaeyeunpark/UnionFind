.PHONY: help clean format

help:
	@echo "Please use \`make <target>\` where <target> is one of"
	@echo "  clean"
	@echo "  format"

clean:
	rm -rf ./build/* externals/catch.hpp externals/nlohmann

format:
	@./build_utils/format_cpp.sh
