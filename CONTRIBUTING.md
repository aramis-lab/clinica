# Contributing to Clinica

There are many ways in which you can contribute to the ongoing development of this project. For example, you can:

* [Report a bug](https://github.com/aramis-lab/clinica/issues/new?assignees=&labels=bug&template=bug_report.md&title=).
* [Discuss the current state of the code](https://github.com/aramis-lab/clinica/issues/new?assignees=&labels=&template=discussion.md&title=).
* [Submit a bugfix](https://github.com/aramis-lab/clinica/compare).
* [Propose a new feature](https://github.com/aramis-lab/clinica/issues/new?assignees=&labels=enhancement&template=feature_request.md&title=).
* [Become a maintainer](mailto:clinica-user@googlegroups.com).


## Asking questions, reporting a bug, or discussing the state of the code

You have two ways for asking questions, reporting a bug, or discussing the current state of the code (feedbacks are more than welcome).

- [Open an issue on GitHub](https://github.com/aramis-lab/clinica/issues).
- If you are less familiar with GitHub, you can still contact us using the [Google Group dedicated to Clinica](https://groups.google.com/g/clinica-user).

## Making direct contributions

### Open an issue on GitHub

The very first thing to do is to open an issue on GitHub in which you describe the bug, or the feature you want to work on.

By doing this, you will get immediate feedback from the developer team on the contribution you are planning to make (is it already being addressed by someone else, is it useful and within the scope of the project...).

This preliminary discussion will save you a lot of time as the developer team will also be able to point you to the right place to modify code for instance.

### Setup your local branch

Once the issue has been discussed and you are ready to implement your contribution, you need to get a copy of the source code on which you can make edits.

We use a fairly standard process in the open source world, where the contributors fork the project and work on their fork prior to making a Pull Request to the upstream repository.

If this is your first time contributing to an open source project, please read this very clear [tutorial](https://github.com/firstcontributions/first-contributions).

You will find how to **fork** this project, **clone** your fork, make a **branch**, add modifications and create a **pull request** using the command line or [GUI tools](https://github.com/firstcontributions/first-contributions#tutorials-using-other-tools).

To sum up these steps:

- Go to the [repository](https://github.com/aramis-lab/clinica) and click the "Fork" button to create your own copy of the project.
- Clone the project to your local computer:

```{.sourceCode .bash}
git clone git@github.com:<username>/clinica.git
```

where `<username>` is your GitHub username.

- Create a branch for the feature you want to work on. Since the branch name will appear on the pull request, use a meaningful name:

```{.sourceCode .bash}
git switch -c <branch_name>
```

(e.g. `git switch -c fix-something-important`)

- Commit locally as you progress (`git add` and `git commit`)
- To push your changes to your fork on GitHub, type:

```{.sourceCode .bash}
git push origin <branch_name>
```

In order to normalize code style across Clinica, we use a combination of different tools (formatter, linter...).

The easiest way to use these tools with the same set of rules as the one we use in Clinica is to install a pre-commit hook, which will trigger them automatically when doing a commit locally.

Please see [instructions here](https://pre-commit.com/) about how to install the pre-commit hook on your local machine (shortly: run `pre-commit install` inside the cloned folder).

### Submit a Pull Request

At this point, your are happy with your contribution and you want to merge it with the official source code of Clinica.

A first thing to verify is whether your changes broke some unit tests. You can verify this by running:

```bash
pytest -v test/unittests
```

This will run the entire suite of unit tests (should take a few minutes) and tell you how many passed and how many failed.

If all tests passed, it is a good sign, and you can proceed with opening a Pull Request. If some tests failed, you can inspect which one failed and see whether this is expected or not.

Sometime, a fix or a new feature will break existing unit tests and these will need to be updated, so broken unit tests don't mean that we will reject the PR, but these will need to be fixed during the review process.

To submit a Pull Request, go to your fork repository on GitHub. The new branch will show up with a green Pull Request button. Simply click it.

This will open a new page on which you need to provide a title and a description, and both are important.

The title should start with a special tag indicating the nature of the PR. Chose the appropriate tag carefully among this list:

- `[CI]`: The Pull Request is a contribution to the continuous integration.
- `[DOC]`: The Pull Request is a documentation contribution.
- `[ENH]`: The Pull Request is an enhancement, usually a new feature or a functional improvement to an existing one.
- `[FIX]`: The Pull Request is a bug fix.
- `[MAINT]`: The Pull Request is doing maintenance related things like upgrading dependencies.
- `[REF]`: The Pull Request is doing a code refactoring without adding functionalities.
- `[TEST]`: The Pull Request is adding or modifying tests.

The description depends a bit on the context but should clearly starts with `Closes #XXXX` where XXXX is the number of the issue you opened initially.

If the issue describes both the context and what needs to be done, then this could be enough. If there are additional information you can provide on your implementation, then this is the right place to write them.

### The review process

Your pull request will then be under review. A friendly discussion with the reviewer(s) will be done in order to improve the quality of your contribution if needed.

If so, update your pull request until your changes are pushed up. The pull request will update automatically and will finally be merged by the reviewer(s).

### Keeping your fork up to date

If you plan to contribute regularly to Clinica, you have to synchronize your fork regularly with the main repository. This can be done by following these steps:

1. Add the main Clinica repository to your remote:

   ```{.sourceCode .bash}
   git remote add upstream https://github.com/aramis-lab/clinica.git
   ```

   `upstream` is an arbitrary name, you can call it whatever you want, but upstream is a convention.

2. Get the data associated with the remote (it will only download it and do nothing else):

   ```{.sourceCode .bash}
   git fetch upstream
   ```

3. Let's say you want the latest version of the `dev` branch. Go to your own `dev` branch, then pull from the `upstream` remote:

   ```{.sourceCode .bash}
   git checkout dev
   git pull upstream dev
   ```

Now the `dev` branch of your own fork is up to date with `aramislab/clinica` dev.

### Solving conflicts in pull request

It may happen that your pull request has conflicts (your feature branch tries to modify the same portion of code as a recently made commit on the `dev` branch).

You must integrate the latest commit of `dev` into your feature branch. One way of doing it would be to merge `upstream/dev` directly into your feature branch.

We will rather rebase it, so that our git history in Clinica can keep a linear history.

1. Go to your feature branch, and perform an interactive rebasing:

    ```{.sourceCode .bash}
    git checkout newfeature
    git rebase -i upstream/dev
    ```

    Keep only the commits that must be applied, corresponding to `newfeature`, at the top of `upstream/dev` (more information [here](https://thoughtbot.com/blog/git-interactive-rebase-squash-amend-rewriting-history)).

    The conflict will appear once you validate the rebasing.

2. Solve the conflicts, following git instructions.

3. Once it is done, you can force push (because you have rewritten history):

    ```{.sourceCode .bash}
    git push --force
    ```

    or

    ```{.sourceCode .bash}
    git push --force origin newfeature
    ```

Now your Pull Request does not have conflicts anymore!

### Integrating a PR into the main code

Pull Requests should be "squashed and merged" and the name used should reflect the convention described above.

Using this scheme enables Clinica to display a simple history on `dev` and `main` branches where the nature of the changes are easily understandable.

## Coding style conventions

For Python files, we use the [numpydoc](https://numpydoc.readthedocs.io/en/latest/format.html) style for docstrings.

The [PEP 8](https://www.python.org/dev/peps/pep-0008/) convention is used but some rules are ignored (they are listed in the [`.pep8speaks.yml` configuration file](https://github.com/aramis-lab/clinica/blob/dev/.pep8speaks.yml).

For Markdown files, currently we don't have a consistent style (but try to keep line length under 80 characters).

## Documentation

This project uses the [MkDocs](https://www.mkdocs.org/) tool with the [Material theme](https://squidfunk.github.io/mkdocs-material/) and extra plugins to generate the documentation.

A folder named `docs` contains all the files to build the documentation. These files are written in Markdown format.

The documentation is built automatically for each Pull Request. You can preview the output for your PR by replacing and copying this URL in your browser:

```
https://aramislab.paris.inria.fr/clinica/docs/public/<your-PR-ID>
```

## License

By contributing to Clinica, you agree that your contributions will be licensed under its [**License**](https://github.com/aramis-lab/clinica/blob/dev/LICENSE.txt).
