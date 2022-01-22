PACKAGES=clinica test

.PHONY: build
build:
	poetry build

.PHONY: clean.doc
clean.doc:
	$(RM) -rf site

.PHONY: config.testpypi
config.testpypi:
	poetry config repositories.testpypi https://test.pypi.org/legacy

.PHONY: doc
doc: clean.doc env.doc
	poetry run mkdocs build

.PHONY: env
env: env.dev env.doc

.PHONY: env.dev
env.dev:
	poetry install

.PHONY: env.doc
env.doc:
	poetry install --extras docs

.PHONY: format
format: format.black format.isort

.PHONY: format.black
format.black: env.dev
	poetry run black --quiet $(PACKAGES)

.PHONY: format.isort
format.isort: env.dev
	poetry run isort --quiet $(PACKAGES)

.PHONY: lint
lint: lint.black lint.isort

.PHONY: lint.black
lint.black: env.dev
	poetry run black --check --diff $(PACKAGES)

.PHONY: lint.isort
lint.isort: env.dev
	poetry run isort --check --diff $(PACKAGES)

.PHONY: publish
publish: publish.pypi

.PHONY: publish.pypi
publish.pypi: build
	poetry publish

.PHONY: publish.testpypi
publish.testpypi: build config.testpypi
	poetry publish --repository testpypi
