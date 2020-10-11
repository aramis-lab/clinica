# Contributing to this repo

There are many ways in which to contribute to the ongoing improvement of this list of tutorials and other resources useful for open science and neuroimaging, EEG and MEG.



## Asking questions or reporting an issue

You have two ways for asking questions or reporting an issue:

- [Open an issue on GitHub](https://github.com/aramis-lab/clinica/issues).

- If you are less familiar with GitHub, you can still write us using [Google Groups dedicated to Clinica](https://groups.google.com/g/clinica-user).



## Making direct contributions

If this is your first time contributing to an open source project, please read this very clear [tutorial](https://github.com/firstcontributions/first-contributions). You will find how to **fork** this project, **clone** your fork, make a **branch**, add modifications and create a **pull request** using command line or [GUI tools](https://github.com/firstcontributions/first-contributions#tutorials-using-other-tools).

To sum up these steps:

- Go to the [repository](https://github.com/aramis-lab/clinica) and click the "Fork" button to create your own copy of the project.

- Clone the project to your local computer:
```Shell
$ git clone git@github.com:<username>/clinica.git
```
where `<username>` is your GitHub username.

- Create a branch for the feature you want to work on. Since the branch name will appear on the pull request, use a meaningful name:
```Shell
$ git checkout -b <branch_name>
```
(e.g. `git checkout -b fix_something_important`)

- Commit locally as you progress (`git add` and `git commit`)

- To push your changes to your fork on GitHub, type:
```Shell
$ git push origin <branch_name>
```

- Go to GitHub. The new branch will show up with a green Pull Request button. Simply click it.

- You pull request will be under review. A friendly discussion with the reviewer(s) will be done in order to improve the quality of your contribution if needed. If so, update your pull request until your changes are pushed up. The pull request will update automatically and will finally be merged by the reviewer(s).



## Coding/style conventions

Currently, this project does not have a consistent style for Markdown files.

For Python files, we use [Google Python Style Guide](https://google.github.io/styleguide/pyguide.html) for docstrings. [PEP 8](https://www.python.org/dev/peps/pep-0008/) convention is used but some rules are ignored (they are listed on the [`.pep8speaks.yml` configuration file](https://github.com/aramis-lab/clinica/blob/dev/.pep8speaks.yml).



## Testing locally

### Source code
> Currently, Clinica is using private data for both converters and pipelines.

> However, Jenkins logs are public so you can see the behaviour of Clinica tests on PR.


### Documentation

This project uses [MkDocs](https://www.mkdocs.org/) tool with [Material theme](https://squidfunk.github.io/mkdocs-material/) and extra plugins to generate the documentation.

- To test locally, you will need to install Clinica using [developer installation](http://www.clinica.run/doc/Installation/#install-clinica)
(If you are working on your *fork*, simply replace `git@github.com:aramis-lab/clinica.git` by `git@github.com:<username>/clinica.git` where `<username>` is your GitHub username)

- Once done, you need to run MkDocs. Simply type:
```Shell
$ mkdocs serve
```

- Finally, open up [`http://127.0.0.1:8000/`](http://127.0.0.1:8000/) in your browser, and you should see the default home page being displayed.

Please note that any modifications on the Markdown files or the configuration file (`mdkocs.yml`) will automatically update the [localhost website](http://127.0.0.1:8000).



## License

By contributing, you agree that your contributions will be licensed under its
[**License**](https://github.com/aramis-lab/clinica/blob/dev/LICENSE.txt).
