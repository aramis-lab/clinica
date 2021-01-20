# Contributing to this repo

There are many ways in which you can contribute to the ongoing development of
this project. For example, you can:

* [Report a bug](https://github.com/aramis-lab/clinica/issues/new?assignees=&labels=bug&template=bug_report.md&title=).
* [Discuss of the current state of the code](https://github.com/aramis-lab/clinica/issues/new?assignees=&labels=&template=discussion.md&title=).
* [Submit a bugfix](https://github.com/aramis-lab/clinica/compare).
* [Propose a new feature](https://github.com/aramis-lab/clinica/issues/new?assignees=&labels=enhancement&template=feature_request.md&title=).
* [Become a maintainer](mailto:clinica-user@googlegroups.com).


## Asking questions or reporting a bug or discussing the state of the code

You have two ways for asking questions or reporting a bug or discussing about
he current state of the code (feedback is more than welcome).

- [Open an issue on GitHub](https://github.com/aramis-lab/clinica/issues).

- If you are less familiar with GitHub, you can still write us using [Google
  Groups dedicated to Clinica](https://groups.google.com/g/clinica-user).


## Making direct contributions (e.g. submitting a fix)

If this is your first time contributing to an open source project, please read
this very clear
[tutorial](https://github.com/firstcontributions/first-contributions). You will
find how to **fork** this project, **clone** your fork, make a **branch**, add
modifications and create a **pull request** using command line or [GUI
tools](https://github.com/firstcontributions/first-contributions#tutorials-using-other-tools).

To sum up these steps:

- Go to the [repository](https://github.com/aramis-lab/clinica) and click the
  "Fork" button to create your own copy of the project.

- Clone the project to your local computer:

```Shell
$ git clone git@github.com:<username>/clinica.git
```
where `<username>` is your GitHub username.

- Create a branch for the feature you want to work on. Since the branch name
  will appear on the pull request, use a meaningful name:

```Shell
$ git checkout -b <branch_name>
```
(e.g. `git checkout -b fix_something_important`)

- Commit locally as you progress (`git add` and `git commit`)

- To push your changes to your fork on GitHub, type:

```Shell
$ git push origin <branch_name>
```

- Go to GitHub. The new branch will show up with a green Pull Request button.
  Simply click it.

- You pull request will be under review. A friendly discussion with the
  reviewer(s) will be done in order to improve the quality of your contribution
  if needed. If so, update your pull request until your changes are pushed up.
  The pull request will update automatically and will finally be merged by the
  reviewer(s).

### Keep your fork up to date

If you plan to contribute regularly to Clinica, you have to synchronise your
fork regularly with the main repository. This can be done by following these
steps:

1. Add the main Clinica repository to your remote:

   ```
   git remote add upstream https://github.com/aramis-lab/clinica.git
   ```

   `upstream` is an arbitrary name, you can call it whatever you want, but
   upstream is a convention.

2. Get the data associated with the remote (it will only download it and do
   nothing else):

   ```
   git fetch upstream
   ```

3. Let's say you want the latest version of the `dev` branch. Go to your own
   dev branch, then pull from the `upstream` remote:

   ```
   git checkout dev
   git pull upstream dev
   ```
Now your `dev` branch of your own fork is up to date with `aramislab/clinica`
dev.

### Solving conflicts in pull request 

It may happen that your pull request has conflicts (your feature branch tries
to modify the same portion of code as a recently made commit on the `dev`
branch). You must integrate the latest commit of dev into your feature branch.
One way of doing it would be to merge `upstream/dev` directly into your feature
branch. We will rather rebase it, so that our git history in Clinica can keep a
linear history.

1. Go to your feature branch, and perform an interactive rebasing: 
   ```
   git checkout newfeature
   git rebase -i upstream/dev
   ```

   Keep only the commits that must be applied at the top of `upstream/dev` that
   correspond to your newfeature (more information
   [here](https://thoughtbot.com/blog/git-interactive-rebase-squash-amend-rewriting-history)

   The conflict will appear once you validate the rebasing.

2. Solve the conflicts, following git instructions.

3. Once it is done, you can force push (because you have rewritten history):
   ```
   git push --force
   ```

   or 
   ```
   git push --force origin newfeature
   ```

Now your Pull Request does not have conflicts anymore!

### Integrating PR into main code

There is no a fully-defined policy for integrating PR into the main code, i.e.,
a PR can be "merged", "squashed and merged" or "rebased and merged", depending
on the nature of the PR and the organisation of the commits in your fork. For
example, a PR that fixes a typo will be integrated by "squash and merge" while
a PR proposing a new pipeline will be probably integrated by "rebase and
merge". In this case we would like to keep intermediate commit messages (if we
consider them important for the project). 

## Coding/style conventions

For Python files, we use [Google Python Style
Guide](https://google.github.io/styleguide/pyguide.html) for docstrings. [PEP
8](https://www.python.org/dev/peps/pep-0008/) convention is used but some rules
are ignored (they are listed on the [`.pep8speaks.yml` configuration
file](https://github.com/aramis-lab/clinica/blob/dev/.pep8speaks.yml).

For Markdown files, currently we don't have a consistent style (but try to keep
line length under 80 characters :)).

## Testing locally

### Source code
> Clinica uses private data for both converters and pipelines.

> However, Jenkins logs are public so you can see the behaviour of Clinica
> tests on PR.


### Documentation

This project uses [MkDocs](https://www.mkdocs.org/) tool with [Material
theme](https://squidfunk.github.io/mkdocs-material/) and extra plugins to
generate the documentation.

A folder named `docs` contains all the files to build the documentation. These
files are written in Markdown format.  Documentation is built automatically for
each PR. You can preview the output for your PR by remplacing and copying this
URL in your browser:

```
https://aramislab.paris.inria.fr/clinica/docs/public/<your-PR-ID>
```

## License

By contributing, you agree that your contributions will be licensed under its
[**License**](https://github.com/aramis-lab/clinica/blob/dev/LICENSE.txt).
