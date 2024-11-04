###############################################################################
#
#   make check            : Check syntax of python files using flake8.
#   make clean            : Return the folder to its original state.
#
###############################################################################

FLAKE8 = flake8

FLAKE8_FILES = ctplanet examples

.PHONY: clean check

help:
	@echo "Commands:"
	@echo ""
	@echo "  check                Check syntax of python files using flake8"
	@echo "  clean                Return the folder to its original state"
	@echo ""

clean:
	@-rm -rf examples/InSight
	@-rm -rf examples/figs/*.dat
	@-rm -rf examples/figs/*.png
	@-rm -rf examples/figs/*.sh

check:
	@$(FLAKE8) --extend-ignore=E741,W605 --exclude=ctplanet/_version.py $(FLAKE8_FILES)
