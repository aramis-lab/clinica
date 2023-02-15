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

# clean			: Clean project tree.
.PHONY: clean
clean: clean.doc clean.py clean.test

.PHONY: clean.doc
clean.doc:
	@$(RM) -rf site/

.PHONY: clean.py
clean.py:
	@find . -name __pycache__ -exec $(RM) -r {} +

.PHONY: clean.test
clean.test:
	@$(RM) -r .pytest_cache/

## doc			: Build the documentation.
.PHONY: doc
doc: clean.doc install.doc
	@$(POETRY) run mkdocs build

## env			: Bootstrap an environment.
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

.PHONY: spellcheck
spellcheck: install.dev
	@$(POETRY) run codespell

.PHONY: test
test: install
	@$(POETRY) run python -m pytest -v test/unittests
