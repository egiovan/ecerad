install:
	fpm install --profile release --prefix ./local
	cd py_ecerad && $(MAKE)

clean:
	fpm clean --all
	cd py_ecerad && $(MAKE) clean