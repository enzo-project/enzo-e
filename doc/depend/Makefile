.PHONY: depend depend-eps depend-png
depend: depend-eps depend-png

depend-eps:
	dot -Teps depend-libraries.dot > depend-libraries.eps
	dot -Teps depend-includes.dot > depend-includes.eps
depend-png:
	dot -Tpng depend-libraries.dot > depend-libraries.png
	dot -Tpng depend-includes.dot > depend-includes.png
