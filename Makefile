PACKAGES := clinica test
POETRY ?= poetry
CONDA ?= conda
CONDA_ENV ?= "./env"

.PHONY: help
help: Makefile
	@echo "Commands:"
	@sed -n 's/^##//p' $<

.PHONY: check.lock
check.lock:
	@$(POETRY) lock --check

.PHONY: clean.doc
clean.doc:
	@$(RM) -rf site

.PHONY: config.pypi
config.pypi:
ifdef PYPI_TOKEN
	@$(POETRY) config pypi-token.pypi "${PYPI_TOKEN}"
else
	$(error "Missing API token for PyPI repository")
endif

## doc			: Build the documentation.
.PHONY: doc
doc: clean.doc install.doc
	@$(POETRY) run mkdocs build

## env			: Bootstap an environment.
.PHONY: env
env: env.dev

.PHONY: env.conda
env.conda:
	@$(CONDA) env create -p $(CONDA_ENV)

.PHONY: env.dev
env.dev: install

## format			: Format the codebase.
.PHONY: format
format: install.dev format.black format.isort

.PHONY: format.black
format.black:
	$(info Formatting code with black)
	@$(POETRY) run black --quiet $(PACKAGES)

.PHONY: format.isort
format.isort:
	$(info Formatting code with isort)
	@$(POETRY) run isort --quiet $(PACKAGES)

## install		: Install the project.
.PHONY: install
install: check.lock
	@$(POETRY) install

.PHONY: install.dev
install.dev: check.lock
	@$(POETRY) install --only dev

.PHONY: install.doc
install.doc: check.lock
	@$(POETRY) install --only docs

## lint			: Lint the codebase.
.PHONY: lint
lint: install.dev lint.black lint.isort

.PHONY: lint.black
lint.black:
	$(info Linting code with black)
	@$(POETRY) run black --check --diff $(PACKAGES)

.PHONY: lint.isort
lint.isort:
	$(info Linting code with isort)
	@$(POETRY) run isort --check --diff $(PACKAGES)

## lock 		: Refresh locked dependencies.
.PHONY: lock
lock:
	@$(POETRY) lock --no-update

## publish		: Publish the package to PyPI.
.PHONY: publish
publish: publish.pypi

.PHONY: publish.pypi
publish.pypi: config.pypi
	@$(POETRY) publish --build

.PHONY: spellcheck
spellcheck: install.dev
	@$(POETRY) run codespell

.PHONY: test
test: install
	@$(POETRY) run python -m pytest -v test/unittests
